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

} // anonymous namespace


AMRComplex::AMRComplex (const Params& params)
    : m_params(params)
{
    PlotFileData pf(m_params.infile);

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
    buildElements();
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
            Long badGhosts = 0, orphans = 0;
            const Box& dm = m_geom[lev].Domain();
            IntVect  sampleIV(AMREX_D_DECL(0,0,0));
            bool     haveSample = false;
            for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
                const Box  gbx = mfi.fabbox();
                const Box& vbx = mfi.validbox();
                auto const& a  = m_nodes[lev].const_array(mfi);
                AMREX_LOOP_3D(gbx, i, j, k,
                {
                    const IntVect iv(AMREX_D_DECL(i,j,k));
                    if (vbx.contains(iv)) { continue; }
                    const Node& n = a(i,j,k);

                    // Only ghosts *stamped from the coarser level* matter here.
                    // A ghost carrying a node of this level is a neighbouring
                    // box's cell, and it is perfectly legitimate for that to be
                    // COVERED by the next level up.
                    if (n.level < lev && n.type == Node::COVERED) { ++badGhosts; }

                    // A ghost still INIT after both propagation steps has no
                    // owner anywhere in the hierarchy. That is legitimate only
                    // outside a NON-periodic wall, where the complex simply
                    // stops half a cell short. Anywhere within the
                    // periodically-extended domain it is a hole: every element
                    // that would bridge it is silently not generated, and the
                    // interface leaks volume with nothing at all to indicate
                    // it. The case that produces it is a level touching a
                    // periodic boundary without wrapping -- its ghost then
                    // coarsens to a cell of the next level down that does not
                    // exist there, because the level across the wrap is two or
                    // more levels coarser and this stitch only ever spans one.
                    if (n.type == Node::INIT) {
                        bool inHier = true;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            if ((iv[d] < dm.smallEnd(d) || iv[d] > dm.bigEnd(d))
                                && !m_params.is_per[d]) { inHier = false; }
                        }
                        if (inHier) {
                            ++orphans;
                            if (!haveSample) { sampleIV = iv; haveSample = true; }
                        }
                    }
                });
            }
            ParallelDescriptor::ReduceLongSum(badGhosts);
            ParallelDescriptor::ReduceLongSum(orphans);
            if (badGhosts > 0) {
                Abort("AMRComplex: level " + std::to_string(lev) + " has "
                      + std::to_string(badGhosts) + " ghost cells whose parent at "
                      "a coarser level is COVERED. The hierarchy "
                      "is not properly nested: level " + std::to_string(lev)
                      + " reaches the edge of level " + std::to_string(lev-1)
                      + ", so the coarse-fine stitch has no valid coarse node to "
                      "attach to and the interface would leak volume silently.");
            }
            if (orphans > 0) {
                if (haveSample) {
                    AllPrint() << "AMRComplex: level " << lev
                               << " orphan ghost at index "
                               << AMREX_D_TERM(sampleIV[0], << "," << sampleIV[1],
                                                            << "," << sampleIV[2])
                               << " (domain " << dm << ")\n";
                }
                Abort("AMRComplex: level " + std::to_string(lev) + " has "
                      + std::to_string(orphans) + " ghost cells inside the "
                      "periodically-extended domain with no owner at any level. "
                      "The stitch walks up every level to level 0, which covers "
                      "the whole domain, so this should be unreachable: a ghost "
                      "inside the domain always has a level-0 parent. Reaching it "
                      "means the plotfile's level 0 does not cover its own domain, "
                      "or the levels are not a refinement of one another. The "
                      "elements that would bridge these ghosts are simply not "
                      "generated, which is a silent gap, so this aborts rather "
                      "than warns."
                      + std::string(m_params.singleLevelStamp
                          ? " single_level_stamp is set, which stops the stitch "
                            "after one level: if that was deliberate, this abort "
                            "IS the expected result and the negative control has "
                            "done its job."
                          : ""));
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
// So the base set for a level is its footprint U grown by one cell on the low
// side, and the whole question is how to partition that set among the boxes so
// that every base cell is visited exactly once, by a rule that is a pure
// function of (cell, BoxArray) and therefore needs no communication.
//
// The rule is corner order. Number the 2^DIM corner offsets as bitmasks, x in
// bit 0; a base cell c then belongs to the box holding the corner c + delta*,
// where delta* is the SMALLEST mask with c + delta* in U, taken modulo the
// domain in periodic directions. For c in U that is delta* = 0, so a valid cell
// belongs to its own box exactly as before. For a stitch base outside U it
// names one corner, hence one box, without reference to any other box's index.
//
// The earlier rule was "the lowest-indexed box that would visit it". That is
// also deterministic, but it is not a function of the cell alone: it compares
// index ranges literally, and on a torus one cell has up to 2^DIM literal
// indices. A level spanning two periodic directions touches both sides of both
// wraps, and a base cell on a wrap *edge* -- index (i,-1,-1) to one box and
// (i,-1,N-1) to another -- was claimed by every box that could see any one of
// its representations. Those elements were generated two and three times over.
// Because the duplicated set is only the wrap edges and corners, it surfaced as
// a small POSITIVE volume residual, on real doubly periodic data and on nothing
// else: every synthetic case had a refined region touching at most one side of
// any wrap, where a single representation exists and the old rule is right.
// Corner order has no representation to disagree about -- the test is always
// "is this corner in the footprint", and the footprint is a set on the torus.
//
// This is answered once per box as box algebra, not once per cell as a query.
// The dominant cost is the delta' = 0 removal, which discards at a stroke the
// whole of V shifted onto itself and leaves only thin slivers for the rest.
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

AMREX_FORCE_INLINE
IntVect deltaOf (int d) noexcept
{
    return IntVect(AMREX_D_DECL((d>>0)&1, (d>>1)&1, (d>>2)&1));
}

// bl minus every cell c whose image c + off is valid at this level, taken
// modulo the domain in periodic directions. This is the only place in the
// ownership rule where periodicity enters, and it enters as box algebra: a
// Node's stored IntVect is never shifted, here or anywhere else.
BoxList removeShiftedValid (const BoxList& bl, const IntVect& off,
                            const BoxArray& ba, const Geometry& geom,
                            bool useWrap)
{
    const Box& dom = geom.Domain();
    Vector<IntVect> shifts;
    std::vector<std::pair<int,Box>> isects;

    BoxList out;
    for (const Box& A : bl) {
        BoxList acc;
        acc.push_back(A);
        const Box T = Box(A).shift(off);

        Vector<IntVect> allShifts(1, IntVect::TheZeroVector());
        if (useWrap && geom.isAnyPeriodic() && !dom.contains(T)) {
            geom.periodicShift(dom, T, shifts);
            for (const auto& sv : shifts) { allShifts.push_back(sv); }
        }

        for (const auto& sv : allShifts) {
            ba.intersections(Box(T).shift(sv), isects);
            if (isects.empty()) { continue; }
            BoxList cut;
            for (const auto& pr : isects) {
                cut.push_back(Box(pr.second).shift(-sv-off));
            }
            acc = subtractFrom(acc, cut);
            if (acc.isEmpty()) { break; }
        }
        out.join(acc);
    }
    return out;
}

} // anonymous namespace

BoxList AMRComplex::ownedBaseRegions (int lev, int boxIdx) const
{
    const BoxArray& ba   = m_grids[lev];
    const Geometry& geom = m_geom[lev];
    const Box       V    = ba[boxIdx];

    BoxList owned;
    owned.push_back(V);            // delta* = 0: valid cells, always ours

    for (int d = 1; d < (1 << AMREX_SPACEDIM); ++d) {
        const IntVect off = deltaOf(d);
        BoxList acc;
        acc.push_back(Box(V).shift(-off));   // cells whose corner d lands in V
        for (int e = 0; e < d && !acc.isEmpty(); ++e) {
            // A smaller mask already places this cell in some box -- this one
            // or another. Either way it is not ours under mask d.
            acc = removeShiftedValid(acc, deltaOf(e), ba, geom,
                                     m_params.legacyOwnership == 0);
        }
        owned.join(acc);
    }

    return owned;
}

//
// Volume the complex should occupy: the full domain in periodic directions,
// and the domain less one level-0 cell in non-periodic ones, where
// cell-centred data leaves the complex inset by dx/2 at each wall.
//
Real AMRComplex::expectedVolume (bool* exact) const
{
    const auto problo = m_geom[0].ProbLo();
    const auto probhi = m_geom[0].ProbHi();

    // Which levels own a cell adjacent to each non-periodic wall. The complex
    // stops half a cell short of the wall, and that half cell belongs to
    // whichever level is VALID there -- not necessarily level 0. A regridded
    // hierarchy refines right up to an inflow or outflow boundary, and
    // assuming level 0 there understates the reference and turns a correct
    // complex into a positive residual, which is the opposite of every other
    // failure mode and so reads as an overlap.
    Vector<Long> wall(2*AMREX_SPACEDIM*std::max(m_nlev,1), 0L);
    for (int lev = 0; lev < m_nlev; ++lev) {
        const Box& dom = m_geom[lev].Domain();
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            auto const& a  = m_nodes[lev].const_array(mfi);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (m_params.is_per[d]) { continue; }
                for (int side = 0; side < 2; ++side) {
                    const int w = (side == 0 ? dom.smallEnd(d) : dom.bigEnd(d));
                    if (vbx.smallEnd(d) > w || vbx.bigEnd(d) < w) { continue; }
                    Box sl = vbx;
                    sl.setSmall(d,w);
                    sl.setBig(d,w);
                    Long cnt = 0;
                    AMREX_LOOP_3D(sl, i, j, k,
                    {
                        if (a(i,j,k).type == Node::VALID) { ++cnt; }
                    });
                    wall[(d*2+side)*m_nlev+lev] += cnt;
                }
            }
        }
    }
    if (wall.size() > 0) {
        ParallelDescriptor::ReduceLongSum(wall.dataPtr(), int(wall.size()));
    }

    bool ex = true;
    Real v = 1.0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real L = probhi[d] - problo[d];
        if (m_params.is_per[d]) { v *= L; continue; }
        Real t = 0.0;
        for (int side = 0; side < 2; ++side) {
            int nl = 0, which = 0;
            for (int lev = 0; lev < m_nlev; ++lev) {
                if (wall[(d*2+side)*m_nlev+lev] > 0) { ++nl; which = lev; }
            }
            if (nl > 1) { ex = false; }
            t += Real(0.5)*m_geom[which].CellSize(d);
        }
        v *= (L - t);
    }
    if (exact) { *exact = ex; }
    return v;
}


//
// One element per valid base cell whose 2^DIM corners are all VALID. Valid
// boxes partition a level, so every element is produced exactly once by
// exactly one rank: no global element set, no deduplication, no gather.
//
void AMRComplex::buildElements ()
{
    m_elements.clear();
    forEachElement([&] (const Node (&cn)[Element::NNODE], const Real*)
    {
        Element e;
        for (int c = 0; c < Element::NNODE; ++c) { e.n[c] = cn[c].gid; }
        m_elements.push_back(e);
    });

    if (m_params.verbose) {
        Long ne = Long(m_elements.size());
        ParallelDescriptor::ReduceLongSum(ne);
        Print() << "AMRComplex: " << ne << " elements" << std::endl;
    }
}

RealVect AMRComplex::nodePosition (const Node& n) const noexcept
{
    const auto plo = m_geom[n.level].ProbLo();
    const auto dx  = m_geom[n.level].CellSize();
    return RealVect(AMREX_D_DECL(plo[0] + (n.iv[0]+Real(0.5))*dx[0],
                                 plo[1] + (n.iv[1]+Real(0.5))*dx[1],
                                 plo[2] + (n.iv[2]+Real(0.5))*dx[2]));
}

AMRComplex::VerifyReport AMRComplex::verify () const
{
    VerifyReport r;

    const auto problo = m_geom[0].ProbLo();
    const auto probhi = m_geom[0].ProbHi();
    RealVect L;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { L[d] = probhi[d] - problo[d]; }
    bool refExact = true;
    const Real dvol = expectedVolume(&refExact);
    r.domainVolume    = dvol;
    r.referenceExact  = refExact;

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
    const bool doCensus = (m_params.faceCensus != 0);
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
            if (doCensus) { ++faceCount[key]; }
            unsigned long long fh = 1469598103934665603ULL;
            for (int q = 0; q < nd; ++q) { mixIn(fh, (unsigned long long)(key[q])); }
            facePar ^= fh;

            // A face on a non-periodic wall: every corner sits at the domain
            // face of its own level. Coarse and fine corners of a stitch face
            // both satisfy that at their own resolution, so this is exact at a
            // refined wall and needs no geometry.
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
    ParallelDescriptor::ReduceLongSum(nElt);
    ParallelDescriptor::ReduceLongSum(nNod);
    ParallelDescriptor::ReduceLongSum(nUnfilled);
    ParallelDescriptor::ReduceLongSum(nodeSum);
    ParallelDescriptor::ReduceLongSum(eltSum);
    ParallelDescriptor::ReduceLongSum(shape.dataPtr(), shape.size());

    if (doCensus) {
        for (const auto& fc : faceCount) {
            if      (fc.second == 1) { ++r.facesBoundary; }
            else if (fc.second == 2) { ++r.facesInterior; }
            else                     { ++r.facesBad; }
        }
        r.censusRun = true;
    }
    r.facesDegenerate = degenFaces;

    // Exact duplicate-element census. A diagnostic, not the invariant: it costs
    // a map over every element, and faceParityInterior states the same thing
    // for nothing. It is here because when it does fire it says precisely what
    // is wrong, and prints where.
    if (doCensus) {
        std::map<std::array<long,Element::NNODE>,int> ecount;
        forEachElement([&] (const Node (&cn)[Element::NNODE], const amrex::Real*)
        {
            std::array<long,Element::NNODE> k;
            for (int c = 0; c < Element::NNODE; ++c) { k[c] = cn[c].gid; }
            std::sort(k.begin(), k.end());
            ++ecount[k];
        });
        for (const auto& e : ecount) {
            if (e.second > 1) { r.duplicateElements += e.second - 1; }
        }
        int shown = 0;
        forEachElement([&] (const Node (&cn)[Element::NNODE], const amrex::Real*)
        {
            if (shown >= m_params.reportBadFaces) { return; }
            std::array<long,Element::NNODE> k;
            for (int c = 0; c < Element::NNODE; ++c) { k[c] = cn[c].gid; }
            std::sort(k.begin(), k.end());
            if (ecount[k] < 2) { return; }
            ++shown;
            AllPrint() << "~ DUP x" << ecount[k] << " :";
            for (int c = 0; c < Element::NNODE; ++c) {
                AllPrint() << " L" << cn[c].level << ":"
                           << AMREX_D_TERM(cn[c].iv[0], << "," << cn[c].iv[1],
                                                        << "," << cn[c].iv[2]);
            }
            AllPrint() << "\n";
        });
    }

    // Localisation. A second pass rather than bookkeeping in the first: the
    // incidence of a face is not known until every element has been seen, so
    // there is nothing useful to record the first time through.
    if (doCensus && m_params.reportBadFaces > 0) {
        int shown = 0;
        forEachElement([&] (const Node (&cn)[Element::NNODE], const amrex::Real*)
        {
            if (shown >= m_params.reportBadFaces) { return; }
            for (int f = 0; f < NFACE; ++f) {
                if (shown >= m_params.reportBadFaces) { return; }
                FaceKey key;
                key[0] = key[1] = key[2] = key[3] = -1;
                int nd = 0;
                for (int c = 0; c < NFCORNER; ++c) {
                    const long g = cn[FACE[f][c]].gid;
                    bool seen = false;
                    for (int q = 0; q < nd; ++q) { if (key[q]==g) { seen=true; break; } }
                    if (!seen) { key[nd++] = g; }
                }
                if (nd < AMREX_SPACEDIM) { continue; }
                std::sort(key.begin(), key.begin()+nd);
                const auto it = faceCount.find(key);
                if (it == faceCount.end() || it->second == 2) { continue; }
                ++shown;
                AllPrint() << "~ face incidence " << it->second << " (face " << f
                           << " of element, corners L:level/index)\n~   face   :";
                for (int c = 0; c < NFCORNER; ++c) {
                    const Node& n = cn[FACE[f][c]];
                    AllPrint() << " L" << n.level << ":"
                               << AMREX_D_TERM(n.iv[0], << "," << n.iv[1], << "," << n.iv[2]);
                }
                AllPrint() << "\n~   element:";
                for (int c = 0; c < Element::NNODE; ++c) {
                    const Node& n = cn[c];
                    AllPrint() << " L" << n.level << ":"
                               << AMREX_D_TERM(n.iv[0], << "," << n.iv[1], << "," << n.iv[2]);
                }
                AllPrint() << "\n";
            }
        });
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
    ParallelDescriptor::ReduceLongSum(r.duplicateElements);
    r.faceParityInterior = faceParIn;
    r.faceParity = facePar;
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
    Print() << "  checksum nodes : " << r.nodeChecksum << "\n";
    Print() << "  checksum elts  : " << r.eltChecksum  << "\n";
    Print() << "  checksum faces : " << r.faceParity << "\n";
    Print() << "  checksum intfac: " << r.faceParityInterior << "\n";

    // Not gated on periodicity any more. Faces on a non-periodic wall are
    // excluded from the parity, so this is a statement about the interior of
    // the complex at any periodicity -- and it is the only watertightness test
    // that stays sharp where the scalar volume reference does not.
    Print() << "  closed complex : " << (r.faceClosed() ? "yes" : "NO")
            << "  (every interior face shared by an even number of elements)\n";
    if (r.censusRun) {
        Print() << "  duplicate elts : " << r.duplicateElements
                << (r.duplicateElements == 0 ? " (ok)" : "  <-- OVERLAP") << "\n";
    }

    if (r.censusRun) {
        Print() << "~ face census (rank-local, exact only at one rank):\n";
        Print() << "~   interior (incidence 2) : " << r.facesInterior << "\n";
        Print() << "~   boundary (incidence 1) : " << r.facesBoundary << "\n";
        Print() << "~   non-manifold (>2)      : " << r.facesBad << "\n";
        Print() << "~   degenerate (<3 nodes)  : " << r.facesDegenerate << "\n";
        Print() << "~   face-conforming        : "
                << (r.faceConforming() ? "yes" : "NO") << "\n";
    } else {
        Print() << "~ face census        : not run (face_census=1 to enable; it is\n"
                << "~                       an ordered-map pass over every face of every\n"
                << "~                       element, and rank-local rather than an\n"
                << "~                       invariant -- checksum faces above is the\n"
                << "~                       rank-independent form and always runs)\n";
        Print() << "~   degenerate (<3 nodes)  : " << r.facesDegenerate << "\n";
    }

    if (!r.referenceExact && !r.passed()) {
        Print() << "  RESULT: NOT APPLICABLE -- a non-periodic wall is adjacent to\n"
                << "    cells of more than one level, so the complex stops a different\n"
                << "    distance short of the wall in different places, its union is\n"
                << "    not a box, and the scalar reference above is a bound rather\n"
                << "    than an identity. Re-run with is_per=1 in every direction: the\n"
                << "    reference is then exactly the domain volume and the residual is\n"
                << "    a sharp gate. The complex itself is built topologically and does\n"
                << "    not change with is_per except in the wrap identification.\n";
    } else if (!r.passed()) {
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
