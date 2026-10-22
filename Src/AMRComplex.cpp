#include "AMRComplex.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>
#include <map>
#include <array>
#include <algorithm>
#include <cstring>
#include <vector>

using namespace amrex;

namespace amrcomplex {

namespace {

// Marks field data that no fill or propagation has touched. A VALID node
// still carrying it means the value plumbing has a hole the identity
// plumbing does not, which is exactly the drift this design guards against.
constexpr Real FIELD_SENTINEL = 1.0e30;

// A 64-bit mix reduced into a Mersenne prime, so that a plain sum over items
// is order-independent and cannot overflow for any realistic item count.
constexpr long CHECKSUM_MOD = 2147483647L;   // 2^31 - 1

AMREX_FORCE_INLINE
void mixIn (unsigned long long& h, unsigned long long v) noexcept
{
    h ^= v;
    h *= 1099511628211ULL;          // FNV-1a prime
    h ^= (h >> 29);
}

AMREX_FORCE_INLINE
long foldChecksum (unsigned long long h) noexcept
{
    return long(h % (unsigned long long)(CHECKSUM_MOD));
}

// Lexicographic index of iv within bx, k-j-i order. A pure function of
// (iv,bx), so every rank agrees without communicating.
AMREX_FORCE_INLINE
long lexIndex (const IntVect& iv, const Box& bx) noexcept
{
    const IntVect lo  = bx.smallEnd();
    const IntVect len = bx.size();
#if (AMREX_SPACEDIM == 2)
    return long(iv[1]-lo[1])*long(len[0]) + long(iv[0]-lo[0]);
#else
    return (long(iv[2]-lo[2])*long(len[1]) + long(iv[1]-lo[1]))*long(len[0])
         +  long(iv[0]-lo[0]);
#endif
}

// Faces of the dual cell in VTK_HEXAHEDRON corner numbering.
#if (AMREX_SPACEDIM == 2)
using FaceKey = std::array<long,4>;
const int NFACE    = 4;
const int NFCORNER = 2;                       // a 2D face is an edge
const int FACE[4][4] = {{0,1,0,1},{1,2,1,2},{2,3,2,3},{3,0,3,0}};
#else
using FaceKey = std::array<long,4>;
const int NFACE    = 6;
const int NFCORNER = 4;
const int FACE[6][4] = {{0,1,2,3},{4,5,6,7},{0,1,5,4},
                        {1,2,6,5},{2,3,7,6},{3,0,4,7}};
#endif

#if (AMREX_SPACEDIM == 3)
// Six tets sharing the main diagonal 0-6; they tile the unit cube exactly.
// Any decomposition gives the same volume for a hexahedron with planar faces.
// Stitch cells need not have planar faces, so a volume residual can mean
// either a genuine gap or a face-diagonal inconsistency between neighbours --
// reportVerify() says so rather than pretending otherwise.
const int TET[6][4] = {{0,1,2,6},{0,2,3,6},{0,3,7,6},
                       {0,7,4,6},{0,4,5,6},{0,5,1,6}};
#endif

// The outer facets of that split -- each simplex's facets not shared with
// another simplex of the same cell -- oriented outward, with the cell face
// each lies on. The wall-flux reference integrates over exactly these, so
// it splits every face the way the volume does.
#if (AMREX_SPACEDIM == 2)
// Triangles {0,1,2},{0,2,3}, counter-clockwise: their outer edges run
// 0->1, 1->2, 2->3, 3->0, one per face.
const int NBFACET = 4;
const int BFACET[4][2]     = {{0,1},{1,2},{2,3},{3,0}};
const int BFACET_FACE[4]   = {0, 1, 2, 3};
#else
// Tet (0,x,y,6) is positively oriented, so its outward facets are (x,y,6)
// and (0,y,x) -- the two not containing the diagonal 0-6.
const int NBFACET = 12;
const int BFACET[12][3]    = {{1,2,6},{0,2,1}, {2,3,6},{0,3,2}, {3,7,6},{0,7,3},
                              {7,4,6},{0,4,7}, {4,5,6},{0,5,4}, {5,1,6},{0,1,5}};
const int BFACET_FACE[12]  = {3,0, 4,0, 4,5, 1,5, 1,2, 3,2};
#endif

} // anonymous namespace


AMRComplex::AMRComplex (const Params& params)
    : m_params(params)
{
    PlotFileData pf(m_params.infile);
    m_time = pf.time();

    int finest = pf.finestLevel();
    if (m_params.finestLevel >= 0) {
        AMREX_ALWAYS_ASSERT(m_params.finestLevel <= pf.finestLevel());
        finest = m_params.finestLevel;
    }
    m_nlev = finest + 1;

    m_grids.resize(m_nlev);
    m_dmap.resize(m_nlev);
    m_geom.resize(m_nlev);
    m_nodes.resize(m_nlev);
    m_ncomp = int(m_params.varnames.size());
    m_field.resize(m_nlev);
    if (m_nlev > 1) { m_refRatio.resize(m_nlev-1); }

    RealBox rb(pf.probLo().data(), pf.probHi().data());
    const int coord = pf.coordSys();

    for (int lev = 0; lev < m_nlev; ++lev) {
        m_grids[lev] = pf.boxArray(lev);
        m_dmap[lev]  = pf.DistributionMap(lev);
        m_geom[lev]  = Geometry(pf.probDomain(lev), &rb, coord,
                                m_params.is_per.data());
        if (lev < m_nlev-1) {
            const int r = pf.refRatio(lev);
            m_refRatio[lev] = IntVect(AMREX_D_DECL(r,r,r));
        }
    }

    numberNodes();    // pure BoxArray bookkeeping; must precede the fill
    classifyNodes(pf);
}

//
// Global numbering.
//
// A BoxArray is globally replicated: ba[i] is the same box on every rank for
// any NProcs. Only the DistributionMapping depends on rank count. So the
// BoxArray's own index order is already a deterministic global ordering, and
// offsets built from box sizes need no communication at all.
//
// gid counts every cell of the hierarchy, not only VALID ones. That costs
// nothing -- no array is allocated per gid -- and buys a great deal: gid
// becomes a pure function of (level, iv), so it can be assigned at fill time
// and then rides FillBoundary and ParallelCopy along with the identity. There
// is no second pass to propagate ids into ghosts. Ids are therefore sparse
// over VALID nodes; formats needing contiguous numbering can compact later.
//
void AMRComplex::numberNodes ()
{
    m_levelOffset.resize(m_nlev+1, 0);
    m_boxOffset.resize(m_nlev);

    for (int lev = 0; lev < m_nlev; ++lev) {
        const BoxArray& ba = m_grids[lev];
        m_boxOffset[lev].resize(ba.size()+1, 0);
        long running = 0;
        for (int i = 0; i < ba.size(); ++i) {
            m_boxOffset[lev][i] = running;
            running += long(ba[i].numPts());
        }
        m_boxOffset[lev][ba.size()] = running;
        m_levelOffset[lev+1] = m_levelOffset[lev] + running;
    }

    m_localRanges.clear();
    m_nLocal = 0;
    const int me = ParallelDescriptor::MyProc();
    for (int lev = 0; lev < m_nlev; ++lev) {
        for (int i = 0; i < m_grids[lev].size(); ++i) {
            if (m_dmap[lev][i] != me) { continue; }
            const long len = m_boxOffset[lev][i+1] - m_boxOffset[lev][i];
            m_localRanges.push_back({m_levelOffset[lev] + m_boxOffset[lev][i], m_nLocal, len});
            m_nLocal += len;
        }
    }
}

long AMRComplex::localIndex (long gid) const noexcept
{
    auto it = std::upper_bound(m_localRanges.begin(), m_localRanges.end(), gid,
                               [] (long g, const LocalRange& r) { return g < r.gidLo; });
    if (it == m_localRanges.begin()) { return -1; }
    --it;
    return (gid < it->gidLo + it->len) ? it->localLo + (gid - it->gidLo) : -1;
}

long AMRComplex::localGid (long i) const noexcept
{
    auto it = std::upper_bound(m_localRanges.begin(), m_localRanges.end(), i,
                               [] (long x, const LocalRange& r) { return x < r.localLo; });
    --it;
    return it->gidLo + (i - it->localLo);
}

int AMRComplex::gidLevel (long gid) const noexcept
{
    const auto it = std::upper_bound(m_levelOffset.begin(), m_levelOffset.end(), gid);
    return int(it - m_levelOffset.begin()) - 1;
}

int AMRComplex::gidOwner (long gid) const noexcept
{
    const int lev = gidLevel(gid);
    const auto& bo = m_boxOffset[lev];
    const long  g  = gid - m_levelOffset[lev];
    const int   i  = int(std::upper_bound(bo.begin(), bo.end(), g) - bo.begin()) - 1;
    return m_dmap[lev][i];
}

void AMRComplex::classifyNodes (PlotFileData& pf)
{
    for (int lev = 0; lev < m_nlev; ++lev) {

        m_nodes[lev].define(m_grids[lev], m_dmap[lev], 1, 1);

        // Field values live alongside the identity, on the same grids with the
        // same ghost width, and are propagated by the same calls. The sentinel
        // makes an unpropagated value detectable rather than plausible.
        if (m_ncomp > 0) {
            m_field[lev].define(m_grids[lev], m_dmap[lev], m_ncomp, 1);
            m_field[lev].setVal(FIELD_SENTINEL);
            for (int n = 0; n < m_ncomp; ++n) {
                const MultiFab src = pf.get(lev, m_params.varnames[n]);
                m_field[lev].ParallelCopy(src, 0, n, 1);
            }
        }

        // Everything, ghosts included, starts INIT naming its own cell. A
        // ghost still INIT after both propagation steps lies outside the
        // hierarchy -- beyond a non-periodic wall -- and elements touching it
        // are simply not generated.
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box& gbx = mfi.fabbox();
            auto const& a  = m_nodes[lev].array(mfi);
            AMREX_LOOP_3D(gbx, i, j, k,
            {
                a(i,j,k) = makeNode(lev, IntVect(AMREX_D_DECL(i,j,k)), Node::INIT);
            });
        }

        // Valid cells, with their gid.
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box& vbx  = mfi.validbox();
            const long base = m_levelOffset[lev] + m_boxOffset[lev][mfi.index()];
            auto const& a   = m_nodes[lev].array(mfi);
            AMREX_LOOP_3D(vbx, i, j, k,
            {
                const IntVect iv(AMREX_D_DECL(i,j,k));
                Node n   = makeNode(lev, iv, Node::VALID);
                n.gid    = base + lexIndex(iv, vbx);
                a(i,j,k) = n;
            });
        }

        // Cells covered by the next finer level are not part of the complex.
        // Only valid cells need marking; ghosts inherit their owner's type
        // through the propagation below.
        if (lev < m_nlev-1) {
            iMultiFab mask = makeFineMask(m_grids[lev], m_dmap[lev],
                                          m_grids[lev+1], m_refRatio[lev],
                                          0, 1);
            for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
                const Box& vbx = mfi.validbox();
                auto const& a  = m_nodes[lev].array(mfi);
                auto const& m  = mask.const_array(mfi);
                AMREX_LOOP_3D(vbx, i, j, k,
                {
                    if (m(i,j,k) == 1) { a(i,j,k).type = Node::COVERED; }
                });
            }
        }

        // Same-level ghosts, including the periodic wrap. FillBoundary copies
        // the Node value and cannot shift an IntVect stored inside it, so the
        // ghost beyond the high face ends up holding the low-face node itself.
        // That copy *is* the periodic identification: no arithmetic, and no
        // opportunity for the coordinate-based defect this design removes.
        m_nodes[lev].FillBoundary(m_geom[lev].periodicity());
        if (m_ncomp > 0) { m_field[lev].FillBoundary(m_geom[lev].periodicity()); }

        // Whatever ghosts remain are covered by the coarser level, if any.
        if (lev > 0) { stampCoarseParents(lev); }

        // Proper nesting is not a nicety here. A level's ghost ring must find
        // VALID parents one level down; if the parent level is itself covered
        // there by this level, the stitch has nothing to attach to and the
        // interface silently leaks volume rather than failing. Catch it at
        // creation instead of explaining a 3% residual later.
        if (lev > 0) {
            // Ghosts that must have found a parent: those inside the domain,
            // and those across a periodic face. Only a ghost beyond a wall may
            // legitimately stay INIT. Checking the domain alone exempted every
            // wrapped ghost, so a level reaching a periodic face with no
            // coarser level on the far side left its seam elements INIT and
            // dropped them -- a gap one stitch layer thick, silently.
            Box perDom = m_geom[lev].Domain();
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (m_geom[lev].isPeriodic(d)) { perDom.grow(d, 1); }
            }
            Long badGhosts = 0;
            for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
                const Box  gbx = mfi.fabbox() & perDom;
                const Box& vbx = mfi.validbox();
                auto const& a  = m_nodes[lev].const_array(mfi);
                AMREX_LOOP_3D(gbx, i, j, k,
                {
                    if (vbx.contains(IntVect(AMREX_D_DECL(i,j,k)))) { continue; }
                    // Only ghosts *stamped from the coarser level* matter here.
                    // A ghost carrying a node of this level is a neighbouring
                    // box's cell, and it is perfectly legitimate for that to be
                    // COVERED by the next level up.
                    const Node& n = a(i,j,k);
                    if (n.type == Node::INIT ||
                        (n.level == lev-1 && n.type == Node::COVERED)) { ++badGhosts; }
                });
            }
            ParallelDescriptor::ReduceLongSum(badGhosts);
            if (badGhosts > 0) {
                Abort("AMRComplex: level " + std::to_string(lev) + " has "
                      + std::to_string(badGhosts) + " ghost cells whose parent at "
                      "level " + std::to_string(lev-1) + " is COVERED or absent. The hierarchy "
                      "is not properly nested: level " + std::to_string(lev)
                      + " reaches the edge of level " + std::to_string(lev-1)
                      + ", so the coarse-fine stitch has no valid coarse node to "
                      "attach to and the interface would leak volume silently.");
            }
        }
    }
}

//
// Fill fine ghosts lying outside this level with their coarse parent's node.
// A block of ref^DIM fine positions resolves to the one coarse node, so
// elements straddling the interface carry repeated gids and are degenerate
// hexahedra. That degeneracy is the stitching: no case table, and no
// restriction on refinement ratio.
//
void AMRComplex::stampCoarseParents (int lev)
{
    // Walk up the levels, not just one. The stitch was written as if a fine
    // ghost always has a parent at lev-1, which is what proper nesting
    // guarantees INSIDE the domain -- and only inside it. Across a periodic
    // wrap the neighbour is whatever happens to be on the far face, and there
    // is no nesting rule relating the two. A doubly periodic flamesheet
    // refines to level 3 at the flame and leaves the far end of the third
    // direction at level 1 or 0, so a level-2 ghost wrapping onto that face
    // finds no level-1 cell at all. There is nothing wrong with that
    // hierarchy; the stitch was simply too narrow. Falling through to the next
    // coarser level costs one extra pass, only over ghosts still unresolved,
    // and the element it produces is a degenerate hexahedron exactly like any
    // other -- a block of ratio^DIM fine positions resolving to one coarse
    // node is the whole mechanism, and it does not care whether ratio is 2 or
    // 4 or 8.
    IntVect rr(AMREX_D_DECL(1,1,1));

    for (int src = lev-1; src >= 0; --src) {
        rr *= m_refRatio[src];

        // Only ghosts that nothing has claimed yet. Counted globally, because
        // the copies below are collective and every rank must agree on whether
        // this pass happens at all.
        Long unresolved = 0;
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box  gbx = mfi.fabbox();
            const Box& vbx = mfi.validbox();
            auto const& f  = m_nodes[lev].const_array(mfi);
            AMREX_LOOP_3D(gbx, i, j, k,
            {
                if (vbx.contains(IntVect(AMREX_D_DECL(i,j,k)))) { continue; }
                if (f(i,j,k).type == Node::INIT) { ++unresolved; }
            });
        }
        ParallelDescriptor::ReduceLongSum(unresolved);
        if (unresolved == 0) { break; }
        if (m_params.singleLevelStamp && src < lev-1) { break; }

        // Parents of the fine ghost ring: coarsen the level's boxes to the
        // source level, then grow by one source cell, which certainly covers
        // them -- coarsen(V.grow(1)) is contained in coarsen(V).grow(1) at any
        // ratio.
        BoxArray cgba = BoxArray(m_grids[lev]).coarsen(rr);
        cgba.grow(1);

        FabArray<BaseFab<Node>> cnodes(cgba, m_dmap[lev], 1, 0);
        MultiFab cfield;
        if (m_ncomp > 0) {
            cfield.define(cgba, m_dmap[lev], m_ncomp, 0);
            cfield.setVal(FIELD_SENTINEL);
            cfield.ParallelCopy(m_field[src], 0, 0, m_ncomp, IntVect(0), IntVect(0),
                                m_geom[src].periodicity());
        }

        for (MFIter mfi(cnodes); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& a = cnodes.array(mfi);
            AMREX_LOOP_3D(bx, i, j, k,
            {
                a(i,j,k) = makeNode(src, IntVect(AMREX_D_DECL(i,j,k)), Node::INIT);
            });
        }

        cnodes.ParallelCopy(m_nodes[src], 0, 0, 1, IntVect(0), IntVect(0),
                            m_geom[src].periodicity());

        // Identity and value are stamped in the same pass, under the same test.
        // Two separate loops would be two chances for a stitch node's value and
        // its identity to disagree, and nothing downstream would notice.
        //
        // A source cell that is itself INIT is not there at all -- outside this
        // level's BoxArray -- so the ghost is left unresolved for the next pass
        // up. A COVERED source IS stamped, deliberately: it means the source
        // level is refined where this ghost needs a parent, which is a nesting
        // error rather than a missing level, and stamping it is what lets the
        // check in classifyNodes name it precisely instead of reporting the
        // generic "no owner anywhere" further down.
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box& gbx = mfi.fabbox();
            auto const& f  = m_nodes[lev].array(mfi);
            auto const& c  = cnodes.const_array(mfi);
            const int nc   = m_ncomp;
            Array4<Real> ff, cf;
            if (nc > 0) { ff = m_field[lev].array(mfi); cf = cfield.array(mfi); }
            AMREX_LOOP_3D(gbx, i, j, k,
            {
                if (f(i,j,k).type == Node::INIT) {
                    const IntVect civ = amrex::coarsen(IntVect(AMREX_D_DECL(i,j,k)), rr);
                    if (c(civ).type == Node::INIT) { continue; }
                    f(i,j,k) = c(civ);
                    for (int n = 0; n < nc; ++n) {
                        ff(i,j,k,n) = cf(AMRCX_IJK(civ),n);
                    }
                }
            });
        }
    }
}

//
// Ownership of elements, by base-cell region.
//
// The obvious rule -- "the box whose valid region contains the base cell" --
// is wrong, and wrong in a way that silently loses exactly the elements this
// design exists to produce. The element bridging a coarse-fine interface on
// the *low* side has its base in a fine ghost stamped from the coarser level:
// no fine box holds it as valid, and the coarse box that would hold it has a
// COVERED corner and is rejected. Under the naive rule nobody generates it,
// leaving a gap one stitch-layer thick around every refined region.
//
// So the base range is the valid box grown by one on the low side, and a base
// cell in that shell belongs to the lowest-indexed box that would visit it.
// BoxArray is globally replicated, so this is deterministic and needs no
// communication: the same element is produced by exactly one rank, at any
// rank count.
//
// This is answered once per box as box algebra, not once per cell as a query.
// The shell is at most AMREX_SPACEDIM boxes, so the whole ownership question
// costs a handful of BoxArray::intersections calls per grid rather than one
// or two per ghost cell.
//
namespace {

// bl minus everything in cut.
BoxList subtractFrom (const BoxList& bl, const BoxList& cut)
{
    if (cut.isEmpty()) { return bl; }
    BoxList out;
    for (const Box& b : bl) { out.join(amrex::complementIn(b, cut)); }
    return out;
}

} // anonymous namespace

BoxList AMRComplex::ownedBaseRegions (int lev, int boxIdx) const
{
    const BoxArray& ba = m_grids[lev];
    const Box       V  = ba[boxIdx];

    BoxList owned;
    owned.push_back(V);          // valid boxes partition the level: all mine

    Box B = V;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { B.growLo(d,1); }
    BoxList shell = amrex::boxDiff(B, V);
    if (shell.isEmpty()) { return owned; }

    std::vector<std::pair<int,Box>> isects;

    // Drop shell cells that are valid at this level, whether directly or
    // through the periodic wrap. Either way the identical element -- same
    // node set, since identity is topological -- is produced where the cell
    // is valid, and generating it here too would double-count the low face.
    Vector<IntVect> shifts;
    const Box& dom = m_geom[lev].Domain();

    BoxList remaining;
    for (const Box& P : shell) {
        BoxList acc; acc.push_back(P);

        Vector<IntVect> allShifts(1, IntVect::TheZeroVector());
        if (m_geom[lev].isAnyPeriodic() && !dom.contains(P)) {
            m_geom[lev].periodicShift(dom, P, shifts);
            for (const auto& sv : shifts) { allShifts.push_back(sv); }
        }

        for (const auto& sv : allShifts) {
            ba.intersections(Box(P).shift(sv), isects);
            if (isects.empty()) { continue; }
            BoxList cut;
            for (const auto& pr : isects) { cut.push_back(Box(pr.second).shift(-sv)); }
            acc = subtractFrom(acc, cut);
            if (acc.isEmpty()) { break; }
        }
        remaining.join(acc);
    }

    // Of the genuine stitch bases that survive, a cell is visited by box A
    // exactly when it lies in A's low-grown valid box. The lowest such index
    // owns it; anything claimed by a lower index is not ours.
    //
    // That comparison is made on ONE name per base cell: its image inside the
    // domain. Across a periodic seam the same stitch base is a ghost of two
    // boxes under two raw indices -- (-1,j) seen from the box at the low end,
    // (n-1,j) from the box at the high end -- and comparing raw indices let
    // each box own it under its own name, so the element was generated twice.
    // The volume overshot by exactly those elements, and only on the low side
    // of a refined region that reaches a seam: the flame-sheet configuration.
    // This is convention 1 -- identity is topological, never positional --
    // applied to base cells as well as nodes.
    //
    // Ownership is decided on the canonical name; the owned region is handed
    // back in this box's own frame, since that is where its ghost cells, and
    // so the element's corners, are read.
    // NEGATIVE CONTROL. Params::legacyOwnership reinstates the defect the
    // canonical-name logic below exists to remove: decide ownership by
    // comparing RAW box indices. Across a periodic seam the same stitch base
    // is a ghost of two boxes under two different raw indices, so each owns it
    // under its own name and the element is generated twice. Confined to the
    // low side of a refined region that reaches a seam, which is why every
    // centred-box synthetic case missed it. Tools/amrComplexShapeTest.sh sets
    // this and requires the volume gate to FAIL; if it ever passes here, the
    // gate has stopped measuring anything.
    if (m_params.legacyOwnership) {
        for (const Box& P : remaining) {
            BoxList claimed;
            Box q = P;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { q.growHi(d,1); }
            ba.intersections(q, isects);
            for (const auto& pr : isects) {
                if (pr.first >= boxIdx) { continue; }
                Box lo = ba[pr.first];
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { lo.growLo(d,1); }
                const Box ov = lo & P;
                if (ov.ok()) { claimed.push_back(ov); }
            }
            owned.join(amrex::complementIn(P, claimed));
        }
        return owned;
    }

    const Geometry& geom = m_geom[lev];
    Vector<IntVect> imageShifts(1, IntVect::TheZeroVector());
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!geom.isPeriodic(d)) { continue; }
        const int n = dom.length(d);
        const int m = imageShifts.size();
        for (int i = 0; i < m; ++i) {
            IntVect lo = imageShifts[i], hi = imageShifts[i];
            lo[d] -= n; hi[d] += n;
            imageShifts.push_back(lo);
            imageShifts.push_back(hi);
        }
    }

    // The domain, unbounded in non-periodic directions: a shell cell beyond a
    // wall has no image to move to, and needs none -- its corners are never
    // VALID, so it generates nothing either way.
    Box canonDom = dom;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!geom.isPeriodic(d)) {
            canonDom.setSmall(d, std::numeric_limits<int>::min()/2);
            canonDom.setBig  (d, std::numeric_limits<int>::max()/2);
        }
    }

    for (const Box& P : remaining) {
        for (const IntVect& s : imageShifts) {
            // C: the part of P whose canonical name is P shifted by s.
            const Box C = Box(P).shift(s) & canonDom;
            if (!C.ok()) { continue; }

            // Box A claims a canonical cell if any image of it lies in A's
            // low-grown box.
            BoxList claimed;
            for (const IntVect& t : imageShifts) {
                const Box Ct = Box(C).shift(t);
                Box q = Ct;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { q.growHi(d,1); }
                ba.intersections(q, isects);
                for (const auto& pr : isects) {
                    if (pr.first >= boxIdx) { continue; }
                    Box lo = ba[pr.first];
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { lo.growLo(d,1); }
                    const Box ov = lo & Ct;
                    if (ov.ok()) { claimed.push_back(Box(ov).shift(-t)); }
                }
            }

            BoxList mine = amrex::complementIn(C, claimed);
            for (Box& b : mine) { b.shift(-s); }
            owned.join(mine);
        }
    }

    return owned;
}

//
// Volume of the domain. This is the reference only when every direction is
// periodic; with a wall, verify() measures the reference from the wall faces
// instead -- see there.
//
Real AMRComplex::expectedVolume () const
{
    const auto problo = m_geom[0].ProbLo();
    const auto probhi = m_geom[0].ProbHi();
    Real v = 1.0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { v *= probhi[d] - problo[d]; }
    return v;
}


RealVect AMRComplex::nodePosition (const Node& n) const noexcept
{
    const auto plo = m_geom[n.level].ProbLo();
    const auto dx  = m_geom[n.level].CellSize();
    return RealVect(AMREX_D_DECL(plo[0] + (n.iv[0]+Real(0.5))*dx[0],
                                 plo[1] + (n.iv[1]+Real(0.5))*dx[1],
                                 plo[2] + (n.iv[2]+Real(0.5))*dx[2]));
}

AMRComplex::VerifyReport AMRComplex::verify (bool faceCensus) const
{
    VerifyReport r;
    r.censusTaken = faceCensus;

    const auto problo = m_geom[0].ProbLo();
    const auto probhi = m_geom[0].ProbHi();
    RealVect L;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { L[d] = probhi[d] - problo[d]; }
    // The reference volume.
    //
    // Fully periodic, it is the domain. With a wall it is not simply the
    // domain less half a cell at the wall: the complex stops half a cell from
    // the wall at whatever level covers the wall there, and slopes between
    // levels where a refined patch on the wall ends. Assuming level 0
    // everywhere was exact for every synthetic case -- no refined level ever
    // reached a wall -- and wrong for the first real one, a flame sheet
    // refined down to its inflow wall.
    //
    // So the reference is measured, by the divergence theorem with
    // F = x_w e_w for a non-periodic direction w: the volume of a region is
    // the flux of F out through its boundary, and the boundary of a correct
    // complex consists only of faces lying against a wall -- every corner in
    // the first cell layer at a wall, at its own level. Summing the flux over
    // exactly those faces of every element gives the volume the complex must
    // have. Faces anywhere else are interior and must cancel in pairs; a gap,
    // an overlap, a fold or a face-diagonal disagreement leaves flux that does
    // not cancel, and the residual shows it. Each face is split into the same
    // triangles the volume's simplex split puts on it, so for a correct
    // complex the two agree to rounding however the faces bend.
    //
    // Blind spot, stated rather than hidden: an element whose every corner
    // lies against one wall -- a thin wedge where a refined patch on the wall
    // ends -- has all its faces counted as wall faces, so its absence would
    // move both sums equally. The face census still sees it.
    int wdir = -1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!m_params.is_per[d]) { wdir = d; break; }
    }
    Real wallSum = 0.0, wallComp = 0.0;
    auto kahan = [] (Real& sum, Real& comp, Real term) {
        const Real ky = term - comp;
        const Real ks = sum + ky;
        comp = (ks - sum) - ky;
        sum  = ks;
    };

    // Compensated (Kahan) summation. Naive summation over a few million
    // elements accumulates ~N*eps of relative error -- about 5e-10 at 2e6
    // elements -- which would force a tolerance loose enough to hide a real
    // gap of the same size. The gate is only useful if it stays sharp as the
    // problem grows.
    Real volSum = 0.0, volComp = 0.0;
    Long nElt = 0, nNod = 0;
    Long degenFaces = 0, nUnfilled = 0;
    Long nodeSum = 0, eltSum = 0;
    unsigned long long facePar = 0ULL, faceParIn = 0ULL;
    std::map<FaceKey,int> faceCount;
    // Element identity is its corner node set, sorted so that a duplicate
    // generated in another box's frame hashes the same. Gated with the census
    // for cost: one map insertion per element.
    std::map<std::array<long,Element::NNODE>,long> eltCount;
    Vector<Long> shape(Element::NNODE+1, 0);

    for (int lev = 0; lev < m_nlev; ++lev) {
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            auto const& a  = m_nodes[lev].const_array(mfi);

            const int nc = m_ncomp;
            Array4<const Real> fv;
            if (nc > 0) { fv = m_field[lev].const_array(mfi); }
            AMREX_LOOP_3D(vbx, i, j, k,
            {
                if (a(i,j,k).type != Node::VALID) { continue; }
                ++nNod;
                for (int n = 0; n < nc; ++n) {
                    if (fv(i,j,k,n) >= FIELD_SENTINEL) { ++nUnfilled; break; }
                }

                // Node hash: identity, global id, and the exact bits of every
                // propagated field value. A rank-dependent difference in the
                // field plumbing shows up here and nowhere else.
                const Node& nd = a(i,j,k);
                unsigned long long h = 1469598103934665603ULL;
                mixIn(h, (unsigned long long)(nd.level));
                AMREX_D_TERM(mixIn(h, (unsigned long long)(nd.iv[0]));,
                             mixIn(h, (unsigned long long)(nd.iv[1]));,
                             mixIn(h, (unsigned long long)(nd.iv[2])););
                mixIn(h, (unsigned long long)(nd.gid));
                for (int n = 0; n < nc; ++n) {
                    unsigned long long bits;
                    const double dv = double(fv(i,j,k,n));
                    std::memcpy(&bits, &dv, sizeof(bits));
                    mixIn(h, bits);
                }
                nodeSum = (nodeSum + foldChecksum(h)) % CHECKSUM_MOD;
            });
        }
    }

    // The element pass goes through the canonical traversal, so verify() and
    // every other consumer agree by construction about which elements exist.
    forEachElement([&] (const Node (&cn)[Element::NNODE], const Real*)
    {
        if (faceCensus) {
            std::array<long,Element::NNODE> ek;
            for (int c = 0; c < Element::NNODE; ++c) { ek[c] = cn[c].gid; }
            std::sort(ek.begin(), ek.end());
            ++eltCount[ek];
        }
        ++nElt;

        int distinct = 0;
        for (int c = 0; c < Element::NNODE; ++c) {
            bool seen = false;
            for (int q = 0; q < c; ++q) {
                if (cn[q].gid == cn[c].gid) { seen = true; break; }
            }
            if (!seen) { ++distinct; }
        }
        ++shape[distinct];

        {
            long g[Element::NNODE];
            for (int c = 0; c < Element::NNODE; ++c) { g[c] = cn[c].gid; }
            std::sort(g, g + Element::NNODE);
            unsigned long long h = 1469598103934665603ULL;
            for (int c = 0; c < Element::NNODE; ++c) {
                mixIn(h, (unsigned long long)(g[c]));
            }
            eltSum = (eltSum + foldChecksum(h)) % CHECKSUM_MOD;
        }

        for (int f = 0; f < NFACE; ++f) {
            FaceKey key;
            key[0] = key[1] = key[2] = key[3] = -1;
            int nd = 0;
            for (int c = 0; c < NFCORNER; ++c) {
                const long g = cn[FACE[f][c]].gid;
                bool seen = false;
                for (int q = 0; q < nd; ++q) { if (key[q]==g) { seen=true; break; } }
                if (!seen) { key[nd++] = g; }
            }
            // A face is degenerate once it has collapsed below the dimension
            // of a face: fewer than 2 distinct nodes in 2D, where a face is an
            // edge, and fewer than 3 in 3D. Testing against 3 unconditionally
            // marked every 2D face degenerate and left the whole check
            // silently inert -- the same defect as a check that never fires.
            if (nd < AMREX_SPACEDIM) { ++degenFaces; continue; }
            std::sort(key.begin(), key.begin()+nd);
            if (faceCensus) { ++faceCount[key]; }
            unsigned long long fh = 1469598103934665603ULL;
            for (int q = 0; q < nd; ++q) { mixIn(fh, (unsigned long long)(key[q])); }
            facePar ^= fh;

            // A face on a non-periodic wall: every corner sits at the domain
            // face of its own level. Coarse and fine corners of a stitch face
            // both satisfy that at their own resolution, so this is exact at a
            // refined wall and needs no geometry. Such a face belongs to one
            // element only and can never cancel, so it is excluded from the
            // interior parity -- the one that answers "is this watertight".
            bool wallFace = false;
            for (int d = 0; d < AMREX_SPACEDIM && !wallFace; ++d) {
                if (m_params.is_per[d]) { continue; }
                bool allLo = true, allHi = true;
                for (int c = 0; c < NFCORNER; ++c) {
                    const Node& nn = cn[FACE[f][c]];
                    const Box&  dd = m_geom[nn.level].Domain();
                    if (nn.iv[d] != dd.smallEnd(d)) { allLo = false; }
                    if (nn.iv[d] != dd.bigEnd(d))   { allHi = false; }
                }
                if (allLo || allHi) { wallFace = true; }
            }
            if (!wallFace) { faceParIn ^= fh; }
        }

        // Corner offsets from corner 0, minimum image in periodic directions.
        const RealVect p0 = nodePosition(cn[0]);
        RealVect q[Element::NNODE];
        for (int c = 0; c < Element::NNODE; ++c) {
            const RealVect p = nodePosition(cn[c]);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real sdd = p[d] - p0[d];
                if (m_params.is_per[d]) {
                    if (sdd >  0.5*L[d]) { sdd -= L[d]; }
                    if (sdd < -0.5*L[d]) { sdd += L[d]; }
                }
                q[c][d] = sdd;
            }
        }

        if (wdir >= 0) {
            // Which walls each corner lies against: bit 2d for the low wall
            // in d, 2d+1 for the high one. A face is a wall face when all its
            // corners share a bit.
            unsigned mask[Element::NNODE];
            for (int c = 0; c < Element::NNODE; ++c) {
                const Box& dl = m_geom[cn[c].level].Domain();
                mask[c] = 0u;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (m_params.is_per[d]) { continue; }
                    if (cn[c].iv[d] == dl.smallEnd(d)) { mask[c] |= 1u << (2*d);   }
                    if (cn[c].iv[d] == dl.bigEnd(d))   { mask[c] |= 1u << (2*d+1); }
                }
            }
            bool onWall[NFACE];
            for (int f = 0; f < NFACE; ++f) {
                unsigned m = ~0u;
                for (int c = 0; c < NFCORNER; ++c) { m &= mask[FACE[f][c]]; }
                onWall[f] = (m != 0u);
            }
            const Real x0 = p0[wdir];
            for (int b = 0; b < NBFACET; ++b) {
                if (!onWall[BFACET_FACE[b]]) { continue; }
#if (AMREX_SPACEDIM == 2)
                // Edge p->q of a counter-clockwise triangle: outward normal
                // times length is (dy, -dx).
                const RealVect& pa = q[BFACET[b][0]];
                const RealVect& pb = q[BFACET[b][1]];
                const Real n[2] = { pb[1]-pa[1], -(pb[0]-pa[0]) };
                const Real xm = x0 + Real(0.5)*(pa[wdir] + pb[wdir]);
                kahan(wallSum, wallComp, xm * n[wdir]);
#else
                // Outward-oriented triangle: area vector (b-a)x(c-a)/2, and
                // x_w is linear over it, so its mean is the vertex mean.
                const RealVect& pa = q[BFACET[b][0]];
                const RealVect& pb = q[BFACET[b][1]];
                const RealVect& pc = q[BFACET[b][2]];
                const RealVect e1 = pb - pa, e2 = pc - pa;
                const Real n[3] = { Real(0.5)*(e1[1]*e2[2]-e1[2]*e2[1]),
                                    Real(0.5)*(e1[2]*e2[0]-e1[0]*e2[2]),
                                    Real(0.5)*(e1[0]*e2[1]-e1[1]*e2[0]) };
                const Real xm = x0 + (pa[wdir] + pb[wdir] + pc[wdir])/Real(3.0);
                kahan(wallSum, wallComp, xm * n[wdir]);
#endif
            }
        }

#if (AMREX_SPACEDIM == 2)
        {
            const int TRI[2][3] = {{0,1,2},{0,2,3}};
            for (int t = 0; t < 2; ++t) {
                const RealVect& u = q[TRI[t][1]];
                const RealVect& v = q[TRI[t][2]];
                const Real term  = Real(0.5)*std::abs(u[0]*v[1] - u[1]*v[0]);
                const Real ky    = term - volComp;
                const Real ksum  = volSum + ky;
                volComp = (ksum - volSum) - ky;
                volSum  = ksum;
            }
        }
#else
        for (int t = 0; t < 6; ++t) {
            const RealVect& u = q[TET[t][1]];
            const RealVect& v = q[TET[t][2]];
            const RealVect& w = q[TET[t][3]];
            const Real det = u[0]*(v[1]*w[2]-v[2]*w[1])
                           - u[1]*(v[0]*w[2]-v[2]*w[0])
                           + u[2]*(v[0]*w[1]-v[1]*w[0]);
            const Real term  = std::abs(det)/Real(6.0);
            const Real ky    = term - volComp;
            const Real ksum  = volSum + ky;
            volComp = (ksum - volSum) - ky;
            volSum  = ksum;
        }
#endif
    });


    ParallelDescriptor::ReduceRealSum(volSum);
    ParallelDescriptor::ReduceRealSum(wallSum);
    const Real dvol = (wdir < 0) ? expectedVolume() : wallSum;
    r.domainVolume = dvol;
    ParallelDescriptor::ReduceLongSum(nElt);
    ParallelDescriptor::ReduceLongSum(nNod);
    ParallelDescriptor::ReduceLongSum(nUnfilled);
    ParallelDescriptor::ReduceLongSum(nodeSum);
    ParallelDescriptor::ReduceLongSum(eltSum);
    ParallelDescriptor::ReduceLongSum(shape.dataPtr(), shape.size());

    for (const auto& fc : faceCount) {
        if      (fc.second == 1) { ++r.facesBoundary; }
        else if (fc.second == 2) { ++r.facesInterior; }
        else                     { ++r.facesBad; }
    }
    r.facesDegenerate = degenFaces;
    if (faceCensus) {
        long dup = 0;
        for (const auto& e : eltCount) { if (e.second > 1) { dup += e.second - 1; } }
        ParallelDescriptor::ReduceLongSum(dup);
        r.duplicateElements = dup;
    }

    // XOR is associative and commutative, so this is independent of both the
    // visit order and the distribution -- unlike the incidence census above,
    // which is irreducibly local.
#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1) {
        unsigned long long tmp[2] = {facePar, faceParIn};
        unsigned long long red[2] = {0ULL, 0ULL};
        MPI_Allreduce(tmp, red, 2, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                      ParallelDescriptor::Communicator());
        facePar   = red[0];
        faceParIn = red[1];
    }
#endif
    r.faceParity = facePar;
    r.faceParityInterior = faceParIn;
    r.nodesUnfilled   = nUnfilled;
    r.nodeChecksum    = nodeSum % CHECKSUM_MOD;
    r.eltChecksum     = eltSum  % CHECKSUM_MOD;

    r.elementVolume = volSum;
    r.relError      = (dvol > 0 ? (volSum - dvol)/dvol : Real(0.0));
    r.nElements     = nElt;
    r.nNodes        = nNod;
    for (int s = 0; s <= Element::NNODE; ++s) { r.shape[s] = shape[s]; }

    return r;
}

void AMRComplex::reportVerify (const VerifyReport& r) const
{
    Print() << "\n=== AMRComplex verify ===\n";
    // Lines prefixed "~" are diagnostics, not invariants: their last bits
    // legitimately depend on how the work was partitioned. The volume sum is
    // compensated per rank but the cross-rank reduction is a plain sum, so
    // the residual moves by an ulp or so with the rank count -- real, benign,
    // and not something to compare for equality. The tolerance test in the
    // RESULT line is the comparable form of the same statement.
    Print() << "  levels         : " << m_nlev << "\n";
    Print() << "  VALID nodes    : " << r.nNodes << "\n";
    Print() << "  elements       : " << r.nElements << "\n";
    Print() << "~ expected volume: " << r.domainVolume << "\n";
    Print() << "~ element volume : " << r.elementVolume << "\n";
    Print() << "~ relative error : " << r.relError << "\n";

    Print() << "  shape histogram (distinct corner nodes):\n";
    for (int s = 0; s <= Element::NNODE; ++s) {
        if (r.shape[s] == 0) { continue; }
        const char* nm = "";
#if (AMREX_SPACEDIM == 2)
        if      (s == 4) { nm = " quad";     }
        else if (s == 3) { nm = " triangle"; }
#else
        if      (s == 8) { nm = " hex";     }
        else if (s == 6) { nm = " prism";   }
        else if (s == 5) { nm = " pyramid"; }
        else if (s == 4) { nm = " tet";     }
#endif
        Print() << "    " << s << nm << " : " << r.shape[s] << "\n";
    }

    if (m_ncomp > 0) {
        Print() << "  fields         : " << m_ncomp << " component(s), "
                << r.nodesUnfilled << " VALID node(s) unfilled"
                << (r.fieldsFilled() ? " (ok)" : "  <-- PROBLEM") << "\n";
    }
    // Printed in the format Tools/amrComplexShapeTest.sh parses; it reads this
    // line to tell "the volume is wrong because elements are duplicated" from
    // "the volume is wrong for some other reason", which is the difference
    // between a diagnostic control and a red light.
    if (r.duplicateElements >= 0) {
        Print() << "  duplicate elts : " << r.duplicateElements << "\n";
    } else {
        Print() << "  duplicate elts : not counted (face_census=1 to count)\n";
    }
    Print() << "  checksum nodes : " << r.nodeChecksum << "\n";
    Print() << "  checksum elts  : " << r.eltChecksum  << "\n";
    Print() << "  checksum faces : " << r.faceParity << "\n";

    bool allPeriodic = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!m_params.is_per[d]) { allPeriodic = false; }
    }
    // faceParityInterior discounts the domain walls, so the closure test is
    // meaningful whatever the periodicity. With every direction periodic
    // there are no walls to discount and the two parities coincide.
    Print() << "  closed complex : " << (r.faceClosed() ? "yes" : "NO")
            << "  (every interior face shared by exactly two elements)\n";
    if (!allPeriodic) {
        Print() << "  checksum intr  : " << r.faceParityInterior << "\n";
    }

    if (r.censusTaken) {
        Print() << "~ face census (rank-local, exact only at one rank):\n";
        Print() << "~   interior (incidence 2) : " << r.facesInterior << "\n";
        Print() << "~   boundary (incidence 1) : " << r.facesBoundary << "\n";
        Print() << "~   non-manifold (>2)      : " << r.facesBad << "\n";
        Print() << "~   degenerate (<3 nodes)  : " << r.facesDegenerate << "\n";
        Print() << "~   face-conforming        : "
                << (r.faceConforming() ? "yes" : "NO") << "\n";
    } else {
        Print() << "~ face census    : not taken\n";
    }

    if (!r.passed()) {
        Print() << "  RESULT: FAIL -- the complex does not tile the domain.\n"
                << "    A non-zero residual is either a genuine gap or overlap at a\n"
                << "    coarse-fine interface, or a face-diagonal inconsistency between\n"
                << "    neighbouring stitch cells with non-planar faces. Separate the\n"
                << "    two by re-running with finestLevel=0, where no stitch cells\n"
                << "    exist: a residual there means the fault lies elsewhere.\n";
    } else {
        Print() << "  RESULT: PASS -- volume closes.\n";
    }
    Print() << std::endl;
}

} // namespace amrcomplex
