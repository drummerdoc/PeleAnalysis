#include "AMRComplex.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>
#include <vector>

using namespace amrex;

namespace amrcomplex {

namespace {

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

// Corner offsets of the dual cell, VTK_QUAD / VTK_HEXAHEDRON order.
AMREX_FORCE_INLINE
IntVect cornerOffset (int c) noexcept
{
#if (AMREX_SPACEDIM == 2)
    const int ox[4] = {0,1,1,0};
    const int oy[4] = {0,0,1,1};
    return IntVect(ox[c], oy[c]);
#else
    const int ox[8] = {0,1,1,0,0,1,1,0};
    const int oy[8] = {0,0,1,1,0,0,1,1};
    const int oz[8] = {0,0,0,0,1,1,1,1};
    return IntVect(ox[c], oy[c], oz[c]);
#endif
}

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
    classifyNodes();
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

void AMRComplex::classifyNodes ()
{
    for (int lev = 0; lev < m_nlev; ++lev) {

        m_nodes[lev].define(m_grids[lev], m_dmap[lev], 1, 1);

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

        // Whatever ghosts remain are covered by the coarser level, if any.
        if (lev > 0) { stampCoarseParents(lev); }

        // Proper nesting is not a nicety here. A level's ghost ring must find
        // VALID parents one level down; if the parent level is itself covered
        // there by this level, the stitch has nothing to attach to and the
        // interface silently leaks volume rather than failing. Catch it at
        // creation instead of explaining a 3% residual later.
        if (lev > 0) {
            Long badGhosts = 0;
            for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
                const Box  gbx = mfi.fabbox() & m_geom[lev].Domain();
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
                    if (n.level == lev-1 && n.type == Node::COVERED) { ++badGhosts; }
                });
            }
            ParallelDescriptor::ReduceLongSum(badGhosts);
            if (badGhosts > 0) {
                Abort("AMRComplex: level " + std::to_string(lev) + " has "
                      + std::to_string(badGhosts) + " ghost cells whose parent at "
                      "level " + std::to_string(lev-1) + " is COVERED. The hierarchy "
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
    const IntVect rr = m_refRatio[lev-1];

    // Parents of the fine ghost ring: coarsen the level's boxes, then grow by
    // one coarse cell, which certainly covers them.
    BoxArray cgba = BoxArray(m_grids[lev]).coarsen(rr);
    cgba.grow(1);

    FabArray<BaseFab<Node>> cnodes(cgba, m_dmap[lev], 1, 0);

    for (MFIter mfi(cnodes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& a = cnodes.array(mfi);
        AMREX_LOOP_3D(bx, i, j, k,
        {
            a(i,j,k) = makeNode(lev-1, IntVect(AMREX_D_DECL(i,j,k)), Node::INIT);
        });
    }

    cnodes.ParallelCopy(m_nodes[lev-1], 0, 0, 1, IntVect(0), IntVect(0),
                        m_geom[lev-1].periodicity());

    for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
        const Box& gbx = mfi.fabbox();
        auto const& f  = m_nodes[lev].array(mfi);
        auto const& c  = cnodes.const_array(mfi);
        AMREX_LOOP_3D(gbx, i, j, k,
        {
            if (f(i,j,k).type == Node::INIT) {
                const IntVect civ = amrex::coarsen(IntVect(AMREX_D_DECL(i,j,k)), rr);
                f(i,j,k) = c(civ);
            }
        });
    }
}

//
// Ownership of an element, by its base cell.
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
bool AMRComplex::ownsBaseCell (int lev, const IntVect& iv, int boxIdx) const
{
    std::vector<std::pair<int,Box>> isects;

    // Valid cells belong to the box holding them; valid boxes partition a level.
    m_grids[lev].intersections(Box(iv,iv), isects);
    if (!isects.empty()) { return isects[0].first == boxIdx; }

    // A base cell outside the level's domain in a periodic direction is the
    // wrap image of a cell inside it. If that image is valid, the identical
    // element -- same node pair, since identity is topological -- is already
    // produced there, and generating it here too would double-count the whole
    // low face. Only if the image is not valid is this a genuine stitch base
    // reached through the wrap.
    const Box& dom = m_geom[lev].Domain();
    if (!dom.contains(iv)) {
        IntVect jv = iv;
        bool wrapped = false;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!m_params.is_per[d]) { continue; }
            const int len = dom.length(d);
            while (jv[d] < dom.smallEnd(d)) { jv[d] += len; wrapped = true; }
            while (jv[d] > dom.bigEnd(d))   { jv[d] -= len; wrapped = true; }
        }
        if (wrapped) {
            m_grids[lev].intersections(Box(jv,jv), isects);
            if (!isects.empty()) { return false; }
        }
    }

    // A stamped ghost. Candidates are the boxes meeting this element's corner
    // footprint; of those, the ones that actually visit iv are the ones whose
    // low-grown valid box contains it.
    const Box footprint(iv, iv + IntVect::TheUnitVector());
    m_grids[lev].intersections(footprint, isects);

    int owner = std::numeric_limits<int>::max();
    for (const auto& pr : isects) {
        Box lo = m_grids[lev][pr.first];
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { lo.growLo(d,1); }
        if (lo.contains(iv)) { owner = std::min(owner, pr.first); }
    }
    return owner == boxIdx;
}

//
// Volume the complex should occupy: the full domain in periodic directions,
// and the domain less one level-0 cell in non-periodic ones, where
// cell-centred data leaves the complex inset by dx/2 at each wall.
//
Real AMRComplex::expectedVolume () const
{
    const auto problo = m_geom[0].ProbLo();
    const auto probhi = m_geom[0].ProbHi();
    const auto dx0    = m_geom[0].CellSize();
    Real v = 1.0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real L = probhi[d] - problo[d];
        v *= (m_params.is_per[d] ? L : L - dx0[d]);
    }
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

    for (int lev = 0; lev < m_nlev; ++lev) {
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            Box bbx = mfi.validbox();
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { bbx.growLo(d,1); }
            const int myBox = mfi.index();
            auto const& a  = m_nodes[lev].const_array(mfi);
            AMREX_LOOP_3D(bbx, i, j, k,
            {
                const IntVect iv(AMREX_D_DECL(i,j,k));
                if (!ownsBaseCell(lev, iv, myBox)) { continue; }
                Element e;
                bool ok = true;
                for (int c = 0; c < Element::NNODE; ++c) {
                    const IntVect jv = iv + cornerOffset(c);
                    const Node& n = a(AMREX_D_DECL(jv[0],jv[1],jv[2]));
                    if (n.type != Node::VALID) { ok = false; break; }
                    e.n[c] = n.gid;
                }
                if (ok) { m_elements.push_back(e); }
            });
        }
    }

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
    const Real dvol = expectedVolume();
    r.domainVolume = dvol;

    // Compensated (Kahan) summation. Naive summation over a few million
    // elements accumulates ~N*eps of relative error -- about 5e-10 at 2e6
    // elements -- which would force a tolerance loose enough to hide a real
    // gap of the same size. The gate is only useful if it stays sharp as the
    // problem grows.
    Real volSum = 0.0, volComp = 0.0;
    Long nElt = 0, nNod = 0;
    Vector<Long> shape(Element::NNODE+1, 0);

    for (int lev = 0; lev < m_nlev; ++lev) {
        for (MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            auto const& a  = m_nodes[lev].const_array(mfi);

            AMREX_LOOP_3D(vbx, i, j, k,
            {
                if (a(i,j,k).type == Node::VALID) { ++nNod; }
            });

            Box bbx = mfi.validbox();
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { bbx.growLo(d,1); }
            const int myBox = mfi.index();

            AMREX_LOOP_3D(bbx, i, j, k,
            {
                const IntVect iv(AMREX_D_DECL(i,j,k));
                if (!ownsBaseCell(lev, iv, myBox)) { continue; }
                Node cn[Element::NNODE];
                bool ok = true;
                for (int c = 0; c < Element::NNODE; ++c) {
                    const IntVect jv = iv + cornerOffset(c);
                    const Node& n = a(AMREX_D_DECL(jv[0],jv[1],jv[2]));
                    if (n.type != Node::VALID) { ok = false; break; }
                    cn[c] = n;
                }
                // AMREX_LOOP_3D expands to plain nested loops, so this
                // continues the innermost loop -- which is what is wanted.
                if (!ok) { continue; }

                ++nElt;

                int distinct = 0;
                for (int c = 0; c < Element::NNODE; ++c) {
                    bool seen = false;
                    for (int p = 0; p < c; ++p) {
                        if (cn[p].gid == cn[c].gid) { seen = true; break; }
                    }
                    if (!seen) { ++distinct; }
                }
                ++shape[distinct];

                // Corner offsets from corner 0, minimum image in periodic
                // directions. Identity is topological, so a wrap element's
                // corners are a domain length apart in raw coordinates; using
                // those directly is precisely the defect this design removes.
                const RealVect p0 = nodePosition(cn[0]);
                RealVect q[Element::NNODE];
                for (int c = 0; c < Element::NNODE; ++c) {
                    const RealVect p = nodePosition(cn[c]);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        Real s = p[d] - p0[d];
                        if (m_params.is_per[d]) {
                            if (s >  0.5*L[d]) { s -= L[d]; }
                            if (s < -0.5*L[d]) { s += L[d]; }
                        }
                        q[c][d] = s;
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
        }
    }

    ParallelDescriptor::ReduceRealSum(volSum);
    ParallelDescriptor::ReduceLongSum(nElt);
    ParallelDescriptor::ReduceLongSum(nNod);
    ParallelDescriptor::ReduceLongSum(shape.dataPtr(), shape.size());

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
    Print() << "  levels         : " << m_nlev << "\n";
    Print() << "  VALID nodes    : " << r.nNodes << "\n";
    Print() << "  elements       : " << r.nElements << "\n";
    Print() << "  expected volume: " << r.domainVolume << "\n";
    Print() << "  element volume : " << r.elementVolume << "\n";
    Print() << "  relative error : " << r.relError << "\n";

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
