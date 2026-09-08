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

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        int n_cell = 32;              pp.query("n_cell", n_cell);
        int max_grid_size = 16;       pp.query("max_grid_size", max_grid_size);
        int nlevels = 2;              pp.query("nlevels", nlevels);
        int ref_ratio = 2;            pp.query("ref_ratio", ref_ratio);
        std::string outfile("plt_test"); pp.query("outfile", outfile);

        // Fraction of the domain refined at each successive level.
        Real refine_frac = 0.5;       pp.query("refine_frac", refine_frac);

        AMREX_ALWAYS_ASSERT(nlevels >= 1);
        AMREX_ALWAYS_ASSERT(ref_ratio == 2 || ref_ratio == 4);

        RealBox rb({AMREX_D_DECL(Real(0.0),Real(0.0),Real(0.0))},
                   {AMREX_D_DECL(Real(1.0),Real(1.0),Real(1.0))});
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
                IntVect(AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)));

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
                int w = int(refine_frac*fbx.length(0));
                w = std::max(8, (w/8)*8);
                const int shrink = (fbx.length(0) - w)/2;
                fbx.grow(-((shrink/8)*8));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    fbx.length(0) >= 8, "refined region shrank away; raise n_cell");
                grids[lev].define(fbx);
                prev = fbx;
            }
            grids[lev].maxSize(max_grid_size);
            dmap[lev].define(grids[lev]);
            data[lev].define(grids[lev], dmap[lev], 1, 0);

            const auto dx  = geom[lev].CellSizeArray();
            const auto plo = geom[lev].ProbLoArray();

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
                    AMREX_D_TERM(r2 += (x-Real(0.5))*(x-Real(0.5));,
                                 r2 += (y-Real(0.5))*(y-Real(0.5));,
                                 r2 += (z-Real(0.5))*(z-Real(0.5)););
                    a(i,j,k) = std::sqrt(r2) - Real(0.3);
                });
            }
        }

        Vector<const MultiFab*> mfs(nlevels);
        Vector<int> steps(nlevels, 0);
        for (int lev = 0; lev < nlevels; ++lev) { mfs[lev] = &data[lev]; }

        WriteMultiLevelPlotfile(outfile, nlevels, mfs, {"dist"}, geom,
                                Real(0.0), steps, rr);

        Print() << "wrote " << outfile << ": " << nlevels << " level(s), base "
                << n_cell << "^" << AMREX_SPACEDIM << ", ref_ratio " << ref_ratio << "\n";
        for (int lev = 0; lev < nlevels; ++lev) {
            Print() << "   level " << lev << ": " << grids[lev].size()
                    << " boxes, " << grids[lev].numPts() << " cells\n";
        }
    }
    amrex::Finalize();
    return 0;
}
