//
// Checks the curvature fields amrGradStats writes against an exact answer.
//
// The test surface is the one makeTestPlt already provides: "dist" is the
// signed distance to a sphere, so every level set of it is a sphere of radius
// r about the same centre, and at a cell whose centre is a distance r from
// that point the exact answers are
//
//   |grad dist| = 1        n = rhat
//   div(n)      = 2/r      (the sum of the principal curvatures)
//   K           = 1/r^2    k1 = k2 = 1/r
//
// at every r, not only on the zero level set. So one plotfile tests the whole
// curvature chain over a continuum of curvatures, with no reference run and
// nothing to tune. A distance function is also the one field for which the
// unit normal needs no normalisation, which separates errors in the curvature
// from errors in the regularised normal that feeds it.
//
// Three things this is meant to catch, each of which would otherwise pass
// quietly:
//
//   a Gaussian curvature assembled with the wrong invariant. K and (div n)^2
//   differ by a factor on a sphere and by nothing at all in a flat region, so
//   a mistake here shows up only where the surface is curved -- which is
//   where the number is used.
//
//   the principal pair swapped, or the discriminant clamped the wrong way. On
//   a sphere k1 and k2 are equal, so the sphere alone cannot see it; the
//   ordering is checked instead by requiring k1 >= k2 everywhere in the
//   measured region, which a sign error violates immediately.
//
//   a filter that removes the structure along with the noise. The level sets
//   here are smooth and well resolved, so smoothing must barely move them.
//   Run the same check on a smoothed field and the error must stay small:
//   that is the claim the pass-count design rests on, stated as a test.
//
//   amrCurvVerify3d.gnu.ex infile=plt_sph_g field=dist rlo=0.15 rhi=0.35
//
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_Reduce.H>

#include <iomanip>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        ParmParse pp;
        std::string infile;      pp.get("infile", infile);
        std::string field("dist"); pp.query("field", field);
        Real rlo = 0.15;         pp.query("rlo", rlo);
        Real rhi = 0.35;         pp.query("rhi", rhi);
        // Relative tolerance on the worst cell in the annulus. The default is
        // loose on purpose: this is a second derivative of a field that is
        // only piecewise smooth at the centre, and the point of the check is
        // that the answer is the right quantity, not that it is converged.
        Real tol = 0.05;         pp.query("tol", tol);
        int finestLevel = -1;    pp.query("finestLevel", finestLevel);

        PlotFileData pf(infile);
        int finest = pf.finestLevel();
        if (finestLevel >= 0) { finest = std::min(finest, finestLevel); }
        const int nlev = finest + 1;

        // The centre of the sphere makeTestPlt builds is the centre of the
        // domain; allow it to be given for any other case.
        GpuArray<Real,AMREX_SPACEDIM> ctr;
        {
            Vector<Real> c(AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                c[d] = Real(0.5)*(pf.probLo()[d] + pf.probHi()[d]);
            }
            pp.queryarr("centre", c, 0, AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { ctr[d] = c[d]; }
        }

        const std::string nMean  = "meanCurv_"  + field;
        const std::string nGauss = "gaussCurv_" + field;
        const std::string nK1    = "curv1_"     + field;
        const std::string nK2    = "curv2_"     + field;

        const Vector<std::string>& have = pf.varNames();
        auto present = [&](const std::string& s) {
            return std::find(have.begin(), have.end(), s) != have.end();
        };
        if (!present(nMean)) {
            Print() << "ERROR: " << infile << " has no " << nMean
                    << "; run amrGradStats with curvOf=" << field << "\n";
            amrex::Finalize();
            return 2;
        }
        const bool haveGauss = present(nGauss) && present(nK1) && present(nK2);
        if (!haveGauss) {
            Print() << "note: no " << nGauss << " in this file; checking the"
                    << " mean curvature only (gaussCurv=1 adds the rest)\n";
        }

        // Worst relative error in each quantity, and the worst violation of
        // the ordering, over uncovered cells in the annulus.
        ReduceOps<ReduceOpMax,ReduceOpMax,ReduceOpMax,ReduceOpSum> rop;
        ReduceData<Real,Real,Real,Real> rdat(rop);
        using RT = typename ReduceData<Real,Real,Real,Real>::Type;

        for (int lev = 0; lev < nlev; ++lev) {
            const BoxArray& ba            = pf.boxArray(lev);
            const DistributionMapping& dm = pf.DistributionMap(lev);
            GpuArray<Real,AMREX_SPACEDIM> dx, plo;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                dx[d]  = pf.cellSize(lev)[d];
                plo[d] = pf.probLo()[d];
            }

            iMultiFab cov(ba, dm, 1, 0);
            cov.setVal(0);
            if (lev < nlev-1) {
                const int r = pf.refRatio(lev);
                const iMultiFab fm = makeFineMask(ba, dm, IntVect(0),
                                                  pf.boxArray(lev+1),
                                                  IntVect(AMREX_D_DECL(r,r,r)),
                                                  Periodicity::NonPeriodic(), 0, 1);
                iMultiFab::Copy(cov, fm, 0, 0, 1, 0);
            }

            const MultiFab mc = pf.get(lev, nMean);
            MultiFab gk, k1f, k2f;
            if (haveGauss) {
                gk  = pf.get(lev, nGauss);
                k1f = pf.get(lev, nK1);
                k2f = pf.get(lev, nK2);
            } else {
                gk.define(ba, dm, 1, 0);  gk.setVal(0.0);
                k1f.define(ba, dm, 1, 0); k1f.setVal(0.0);
                k2f.define(ba, dm, 1, 0); k2f.setVal(0.0);
            }

            const bool hg = haveGauss;
            for (MFIter mfi(mc); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const& m = mc.const_array(mfi);
                auto const& g = gk.const_array(mfi);
                auto const& p = k1f.const_array(mfi);
                auto const& q = k2f.const_array(mfi);
                auto const& c = cov.const_array(mfi);
                rop.eval(bx, rdat,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> RT
                {
                    if (c(i,j,k)) { return {Real(0.0),Real(0.0),Real(0.0),Real(0.0)}; }
                    const int iv[3] = {i,j,k};
                    Real r2 = 0;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        const Real x = plo[d] + (iv[d]+Real(0.5))*dx[d] - ctr[d];
                        r2 += x*x;
                    }
                    const Real r = std::sqrt(r2);
                    if (r < rlo || r > rhi) {
                        return {Real(0.0),Real(0.0),Real(0.0),Real(0.0)};
                    }
                    // div(n) is the SUM of the principal curvatures, and a
                    // level set of a distance function in DIM dimensions is
                    // a sphere with DIM-1 of them, each 1/r. So 2/r in 3D
                    // and 1/r in 2D -- the 3D value used unconditionally
                    // would report every correct 2D curvature as 100% wrong.
                    const Real exact = Real(AMREX_SPACEDIM - 1)/r;
                    const Real eM = std::abs(m(i,j,k) - exact)/exact;
                    Real eG = 0, eP = 0, bad = 0;
                    if (hg) {
                        eG = std::abs(g(i,j,k) - Real(1.0)/(r*r))*(r*r);
                        eP = amrex::max(std::abs(p(i,j,k) - Real(1.0)/r),
                                        std::abs(q(i,j,k) - Real(1.0)/r))*r;
                        bad = (q(i,j,k) > p(i,j,k)) ? Real(1.0) : Real(0.0);
                    }
                    return {eM, eG, eP, bad};
                });
            }
        }

        RT hv = rdat.value(rop);
        Real eM = amrex::get<0>(hv);
        Real eG = amrex::get<1>(hv);
        Real eP = amrex::get<2>(hv);
        Real nb = amrex::get<3>(hv);
        ParallelDescriptor::ReduceRealMax(eM);
        ParallelDescriptor::ReduceRealMax(eG);
        ParallelDescriptor::ReduceRealMax(eP);
        ParallelDescriptor::ReduceRealSum(nb);

        Print() << "=== amrCurvVerify ===\n"
                << "  file           : " << infile << "\n"
                << "  field          : " << field << "\n"
                << "  annulus        : " << rlo << " <= r <= " << rhi << "\n"
                << "  tolerance      : " << tol << " (relative, worst cell)\n"
                << std::scientific << std::setprecision(6)
                << "  div(n) vs " << (AMREX_SPACEDIM-1) << "/r  : " << eM << "\n";
        if (haveGauss) {
            Print() << "  K vs 1/r^2     : " << eG << "\n"
                    << "  k1,k2 vs 1/r   : " << eP << "\n"
                    << std::defaultfloat
                    << "  cells with k2 > k1 : " << Long(nb) << "\n";
        }
        Print() << std::defaultfloat;

        const bool ok = (eM < tol)
                     && (!haveGauss || (eG < tol && eP < tol && nb == 0));
        if (ok) {
            Print() << "  RESULT: PASS\n";
        } else {
            Print() << "  RESULT: FAIL -- the curvature fields do not agree"
                    << " with the exact sphere.\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
