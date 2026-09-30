//
// Synthetic multi-level plotfile generator, for exercising AMRComplex without
// production data. Levels are properly nested regions, refined by
// ref_ratio each time; the field is a smooth analytic function so that any
// later isosurface has a known answer too.
//
//   make -f GNUmakefile.AMRComplex EBASE=makeTestPlt
//   ./makeTestPlt3d.gnu.ex n_cell=32 nlevels=2 outfile=plt_test
//
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <algorithm>
#include <string>
#include <cmath>

using namespace amrex;

namespace {

// Carve a refined region out of `allowed`, on a lattice of blk-sized blocks
// aligned to the index origin, so that every block boundary is also a coarse
// cell boundary for any ref_ratio dividing blk.
//
// A single centred convex box -- the only shape this generator produced until
// now -- exercises none of the geometry a regridded hierarchy actually has:
// no step in the coarse-fine interface, no re-entrant corner where a base
// cell draws some of its corners from a same-level neighbour and the rest
// from the coarser level, and no contact with the domain wall. A real
// plotfile is made almost entirely of those configurations.
amrex::BoxList carveRegion (const amrex::BoxList& allowed, int blk,
                            amrex::Real frac, const std::string& shape,
                            const std::string& anchor)
{
    const Box ab = allowed.minimalBox();

    IntVect tlo, thi, nblk;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int base   = ((ab.smallEnd(d) + blk - 1)/blk)*blk;
        const int navail = (ab.bigEnd(d) + 1 - base)/blk;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(navail >= 1,
            "carveRegion: allowed region is narrower than one block; "
            "raise n_cell or lower refine_blk");
        int nk = int(frac*Real(navail));
        nk = std::max(1, std::min(navail, nk));
        int off = 0;
        if      (anchor == "lo")     { off = 0;             }
        else if (anchor == "hi")     { off = navail - nk;   }
        else if (anchor == "center") { off = (navail-nk)/2; }
        else { amrex::Abort("makeTestPlt: unknown refine_anchor " + anchor); }
        tlo[d]  = base + off*blk;
        thi[d]  = tlo[d] + nk*blk - 1;
        nblk[d] = nk;
    }
    amrex::ignore_unused(thi);

#if (AMREX_SPACEDIM == 3)
    const int nkb = nblk[2];
#else
    const int nkb = 1;
#endif

    BoxList blocks;
    for (int kb = 0; kb < nkb; ++kb) {
    for (int jb = 0; jb < nblk[1]; ++jb) {
    for (int ib = 0; ib < nblk[0]; ++ib) {
        amrex::ignore_unused(kb);
        bool keep = true;
        if (shape == "box") {
            keep = true;
        } else if (shape == "lshape") {
            // Drop the upper corner octant: one re-entrant corner, and in 3D
            // three re-entrant edges meeting at it.
            bool up = (ib >= (nblk[0]+1)/2) && (jb >= (nblk[1]+1)/2);
#if (AMREX_SPACEDIM == 3)
            up = up && (kb >= (nblk[2]+1)/2);
#endif
            keep = !up;
        } else if (shape == "checker") {
            // Punch out every other block in the first two directions. The
            // region stays connected, and each punched block contributes four
            // re-entrant edges -- hundreds of them from a handful of blocks.
            keep = !((ib % 2 == 1) && (jb % 2 == 1));
        } else {
            amrex::Abort("makeTestPlt: unknown refine_shape " + shape);
        }
        if (!keep) { continue; }
        const IntVect blo(AMREX_D_DECL(tlo[0]+ib*blk, tlo[1]+jb*blk, tlo[2]+kb*blk));
        blocks.push_back(Box(blo, blo + IntVect(blk-1)));
    }}}

    // Keep only what the parent allows. `allowed` may itself be non-convex --
    // it is, as soon as any level above used a non-convex shape -- and the
    // lattice above was laid out on its bounding box, so this is not a
    // formality.
    BoxList out;
    for (const Box& b : blocks) {
        for (const Box& a : allowed) {
            const Box ov = a & b;
            if (ov.ok()) { out.push_back(ov); }
        }
    }
    out.simplify();
    return out;
}

} // anonymous namespace


int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        // n_cell and prob_hi may be given per direction, so that cells can be
        // anisotropic (dx != dy) and the domain non-cubic. Every case up to
        // now has been a unit cube with one n_cell, which exercises none of
        // the per-direction arithmetic in the complex -- and elongated
        // domains are the normal case for a pipe or a channel.
        int n_cell = 32;              pp.query("n_cell", n_cell);
        Vector<int> ncv(AMREX_SPACEDIM, n_cell);
        if (const int nn = pp.countval("n_cell_vec")) {
            AMREX_ALWAYS_ASSERT(nn == AMREX_SPACEDIM);
            pp.getarr("n_cell_vec", ncv, 0, AMREX_SPACEDIM);
        }
        Vector<Real> phi(AMREX_SPACEDIM, 1.0);
        if (const int np = pp.countval("prob_hi")) {
            AMREX_ALWAYS_ASSERT(np == AMREX_SPACEDIM);
            pp.getarr("prob_hi", phi, 0, AMREX_SPACEDIM);
        }
        int max_grid_size = 16;       pp.query("max_grid_size", max_grid_size);
        int nlevels = 2;              pp.query("nlevels", nlevels);
        int ref_ratio = 2;            pp.query("ref_ratio", ref_ratio);
        std::string outfile("plt_test"); pp.query("outfile", outfile);

        // Fraction of the domain refined at each successive level.
        Real refine_frac = 0.5;       pp.query("refine_frac", refine_frac);
        int  cell_k = 4;              pp.query("cell_k", cell_k);

        // Shape and placement of each refined level. The defaults reproduce
        // the original behaviour exactly: one centred convex box per level.
        std::string refine_shape("box");
        pp.query("refine_shape", refine_shape);
        std::string refine_anchor("center");
        pp.query("refine_anchor", refine_anchor);
        int refine_blk = max_grid_size;
        pp.query("refine_blk", refine_blk);
        AMREX_ALWAYS_ASSERT(refine_blk >= ref_ratio && refine_blk % ref_ratio == 0);

        AMREX_ALWAYS_ASSERT(nlevels >= 1);
        AMREX_ALWAYS_ASSERT(ref_ratio == 2 || ref_ratio == 4);

        RealBox rb({AMREX_D_DECL(Real(0.0),Real(0.0),Real(0.0))},
                   {AMREX_D_DECL(phi[0],phi[1],phi[2])});
        // Periodicity is not recorded in a plotfile; AMRComplex takes it as an
        // input. Writing with is_per=0 here is therefore not a statement about
        // how the file may later be interpreted.
        Array<int,AMREX_SPACEDIM> isper{{AMREX_D_DECL(0,0,0)}};

        Vector<Geometry>            geom(nlevels);
        Vector<BoxArray>            grids(nlevels);
        Vector<DistributionMapping> dmap(nlevels);
        Vector<MultiFab>            data(nlevels);
        Vector<IntVect>             rr(std::max(nlevels-1,1), IntVect(AMREX_D_DECL(ref_ratio,ref_ratio,ref_ratio)));

        BoxList core;
        Box dom(IntVect(AMREX_D_DECL(0,0,0)),
                IntVect(AMREX_D_DECL(ncv[0]-1,ncv[1]-1,ncv[2]-1)));

        // The sphere stays a sphere in physical space whatever the domain
        // shape, so 4*pi*r^2 remains the reference. Radius is set from the
        // shortest side so it always fits.
        Real Lmin = phi[0];
        for (int d = 1; d < AMREX_SPACEDIM; ++d) { Lmin = std::min(Lmin, phi[d]); }
        const Real radius = Real(0.3)*Lmin;

        for (int lev = 0; lev < nlevels; ++lev) {
            if (lev > 0) { dom.refine(ref_ratio); }
            geom[lev].define(dom, &rb, 0, isper.data());

            if (lev == 0) {
                grids[lev].define(dom);
            } else {
                // `core` is the parent's nesting-safe region: parent cells at
                // least one parent cell from the parent's complement, with the
                // domain boundary counted as safe. Anything carved inside its
                // refinement is properly nested by construction, and is free
                // to touch the domain wall -- which a regridded hierarchy does
                // and a centred box never does.
                BoxList allowed;
                for (const Box& b : core) {
                    allowed.push_back(Box(b).refine(ref_ratio));
                }
                BoxList bl = carveRegion(allowed, refine_blk, refine_frac,
                                         refine_shape, refine_anchor);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!bl.isEmpty(),
                    "refined region carved away to nothing; raise n_cell, "
                    "raise refine_frac, or lower refine_blk");
                grids[lev].define(bl);
            }
            grids[lev].maxSize(max_grid_size);

            // Nesting-safe region for the next level, by box algebra rather
            // than by construction, so that it stays right for a non-convex
            // footprint. The domain boundary is not part of the complement, so
            // a level touching the wall does not forbid its child from
            // touching it as well.
            {
                BoxList comp = grids[lev].complementIn(dom);
                core.clear();
                if (comp.isEmpty()) {
                    core.push_back(dom);
                } else {
                    BoxList grown;
                    for (const Box& b : comp) { grown.push_back(amrex::grow(b,1)); }
                    core = amrex::complementIn(dom, grown);
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!core.isEmpty() || lev == nlevels-1,
                    "level footprint leaves no properly nested room for the "
                    "next level; lower refine_frac or use refine_shape=box");
            }
            dmap[lev].define(grids[lev]);
            // dist, plus one linear field per direction. A linear function
            // interpolated linearly along an edge is exact, so a component
            // carrying "x" must come back on the surface equal to the
            // vertex's own x coordinate -- an exact identity to check the
            // component sampling against, not an approximate one.
            // ... plus a cellular field, cos(2 pi k x/L) summed over
            // directions. On a periodic domain it has exactly k^DIM maxima, so
            // just below its maximum the number of connected components above
            // threshold is known exactly -- and the maxima sit ON the domain
            // faces, so finding k^DIM of them is also a test that the wrap is
            // in the adjacency.
            data[lev].define(grids[lev], dmap[lev], 2+AMREX_SPACEDIM, 0);

            const auto dx  = geom[lev].CellSizeArray();
            const auto plo = geom[lev].ProbLoArray();
            GpuArray<Real,AMREX_SPACEDIM> phiA;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { phiA[d] = phi[d]; }
            const Real rad = radius;
            const Real kk  = Real(cell_k);

            for (MFIter mfi(data[lev]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const& a = data[lev].array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    AMREX_D_TERM(const Real x = plo[0] + (i+Real(0.5))*dx[0];,
                                 const Real y = plo[1] + (j+Real(0.5))*dx[1];,
                                 const Real z = plo[2] + (k+Real(0.5))*dx[2];);
                    // Distance to a sphere of radius 0.3 centred in the domain:
                    // an isosurface at 0 has area 4*pi*r^2 exactly.
                    Real r2 = 0;
                    AMREX_D_TERM(r2 += (x-Real(0.5)*phiA[0])*(x-Real(0.5)*phiA[0]);,
                                 r2 += (y-Real(0.5)*phiA[1])*(y-Real(0.5)*phiA[1]);,
                                 r2 += (z-Real(0.5)*phiA[2])*(z-Real(0.5)*phiA[2]););
                    a(i,j,k,0) = std::sqrt(r2) - rad;
                    AMREX_D_TERM(a(i,j,k,1) = x;,
                                 a(i,j,k,2) = y;,
                                 a(i,j,k,3) = z;);
                    Real cc = 0;
                    AMREX_D_TERM(cc += std::cos(Real(2.0)*Real(M_PI)*kk*x/phiA[0]);,
                                 cc += std::cos(Real(2.0)*Real(M_PI)*kk*y/phiA[1]);,
                                 cc += std::cos(Real(2.0)*Real(M_PI)*kk*z/phiA[2]););
                    a(i,j,k,1+AMREX_SPACEDIM) = cc;
                });
            }
        }

        Vector<const MultiFab*> mfs(nlevels);
        Vector<int> steps(nlevels, 0);
        for (int lev = 0; lev < nlevels; ++lev) { mfs[lev] = &data[lev]; }

        Vector<std::string> vnames = {"dist", AMREX_D_DECL("xco","yco","zco"), "cells"};
        WriteMultiLevelPlotfile(outfile, nlevels, mfs, vnames, geom,
                                Real(0.0), steps, rr);

        Print() << "wrote " << outfile << ": " << nlevels << " level(s), base";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << ncv[d]; }
        Print() << ", domain";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << phi[d]; }
        Print() << ", sphere r=" << radius << ", ref_ratio " << ref_ratio << "\n";
        Print() << "   refine_shape " << refine_shape << ", refine_anchor "
                << refine_anchor << ", refine_blk " << refine_blk << "\n";
        for (int lev = 0; lev < nlevels; ++lev) {
            const Box& d0 = geom[lev].Domain();
            const Box  bb = grids[lev].minimalBox();
            Print() << "   level " << lev << ": " << grids[lev].size()
                    << " boxes, " << grids[lev].numPts() << " cells, fill "
                    << Real(grids[lev].numPts())/Real(d0.numPts())
                    << ", touches";
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Print() << " " << (bb.smallEnd(d) == d0.smallEnd(d) ? "lo" : "--")
                        << (bb.bigEnd(d)   == d0.bigEnd(d)   ? "hi" : "--");
            }
            Print() << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}
