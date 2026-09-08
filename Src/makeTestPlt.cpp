//
// Synthetic multi-level plotfile generator, for exercising AMRComplex without
// production data. Levels are nested boxes centred in the domain, refined by
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

using namespace amrex;

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

        Box prev;
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
                prev = dom;
            } else {
                // Refine the previous level's region and shrink it, so each
                // level is properly nested inside its parent with a real
                // buffer of coarse cells. Without that buffer a level can
                // exactly cover its parent, leaving the next level down with
                // nothing but COVERED parents to stitch to -- which produces a
                // silent gap rather than an error.
                Box fbx = Box(prev).refine(ref_ratio);
                // Shrink proportionally in each direction, so an elongated
                // domain gets an elongated refined region rather than a cube
                // that may not fit.
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    int w = int(refine_frac*fbx.length(d));
                    w = std::max(8, (w/8)*8);
                    const int shrink = (fbx.length(d) - w)/2;
                    fbx.grow(d, -((shrink/8)*8));
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                        fbx.length(d) >= 8, "refined region shrank away; raise n_cell");
                }
                grids[lev].define(fbx);
                prev = fbx;
            }
            grids[lev].maxSize(max_grid_size);
            dmap[lev].define(grids[lev]);
            // dist, plus one linear field per direction. A linear function
            // interpolated linearly along an edge is exact, so a component
            // carrying "x" must come back on the surface equal to the
            // vertex's own x coordinate -- an exact identity to check the
            // component sampling against, not an approximate one.
            data[lev].define(grids[lev], dmap[lev], 1+AMREX_SPACEDIM, 0);

            const auto dx  = geom[lev].CellSizeArray();
            const auto plo = geom[lev].ProbLoArray();
            GpuArray<Real,AMREX_SPACEDIM> phiA;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { phiA[d] = phi[d]; }
            const Real rad = radius;

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
                });
            }
        }

        Vector<const MultiFab*> mfs(nlevels);
        Vector<int> steps(nlevels, 0);
        for (int lev = 0; lev < nlevels; ++lev) { mfs[lev] = &data[lev]; }

        Vector<std::string> vnames = {"dist", AMREX_D_DECL("xco","yco","zco")};
        WriteMultiLevelPlotfile(outfile, nlevels, mfs, vnames, geom,
                                Real(0.0), steps, rr);

        Print() << "wrote " << outfile << ": " << nlevels << " level(s), base";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << ncv[d]; }
        Print() << ", domain";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << phi[d]; }
        Print() << ", sphere r=" << radius << ", ref_ratio " << ref_ratio << "\n";
        for (int lev = 0; lev < nlevels; ++lev) {
            Print() << "   level " << lev << ": " << grids[lev].size()
                    << " boxes, " << grids[lev].numPts() << " cells\n";
        }
    }
    amrex::Finalize();
    return 0;
}
