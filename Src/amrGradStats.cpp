//
// Coarse-fine-aware gradient fields for a multi-level plotfile.
//
// Writes a plotfile carrying the requested input variables unchanged plus, for
// each of them, the cell-centred gradient and its magnitude; for each pair, the
// alignment of their gradients; and optionally the unit normal and mean
// curvature of one named field. The output is a plotfile on the same hierarchy,
// so amrJointFit consumes it with no special handling -- derived quantities are
// produced on the block-structured grid with the solver's own operators, and
// the topology and statistics tools downstream never grow a stencil.
//
// Everything here runs where the data is. MLMG's apply and getFluxes are already
// GPU-resident, average_face_to_cellcenter is a ParallelFor inside AMReX, and the
// pointwise work below -- magnitudes, the alignment pair, the unit normal -- is
// ParallelFor too, so no component ever comes back to the host. The only host
// traffic is the plotfile read and write. The per-variable composite pass costs
// one apply and one getFluxes over the hierarchy; at 3e8 cells that is the
// dominant term and it is the same kernel PeleLMeX runs.
//
// The alignment pair is the point of this for the flame work. For two progress
// variables a and b,
//
//   cos = grad a . grad b / (|grad a| |grad b|)     sign: do they agree
//   sin = |grad a x grad b| / (|grad a| |grad b|)   degree of transversality
//
// sin is identically zero wherever b is locally a function of a -- that is,
// wherever a single progress variable suffices, for ANY monotone relation
// between them. So it is a pointwise, threshold-free, reference-free measure of
// where the single-progress-variable picture fails. The same idiom as the
// Takeno flame index, applied to a different pair.
//
// Curvature, and what the names mean. meanCurv_<f> is div(n), which is the
// SUM of the principal curvatures and so twice the mean curvature in the
// differential-geometry sense. The name predates this note and downstream
// files refer to it, so it stays; curv1_ and curv2_ are the principal pair
// themselves and gaussCurv_ their product, so anyone needing the strict mean
// curvature can halve the first or average the other two.
//
// gaussCurv=1 adds the Gaussian and principal curvatures, the shape index and
// the curvedness, at the cost of carrying the normal's Jacobian through the
// curvature step. 3D only. smooth_passes=N low-passes the curvOf field N
// times before anything is differentiated, and renames every quantity derived
// from it so a filtered column cannot be confused with an unfiltered one.
// Together they are how a curvature field is interrogated: sweep N and watch
// the reported realizability numbers.
//
// What those numbers mean, calibrated against makeTestPlt's exact sphere and
// its "noisy" field, which is the same sphere perturbed cell to cell by a
// stated fraction of a cell:
//
//   mean severity   surface                     error in div(n)
//   3e-4            exact, resolved             0.4%
//   2e-3            perturbed by 0.2 cells      1000%
//   4e-3            perturbed by 1 cell         4700%
//   8e-3            perturbed by 4 cells        6000%
//
// So the meter has a low dynamic range and errs towards saying nothing is
// wrong: an order of magnitude above the clean baseline already means the
// curvature magnitudes are not to be trusted at all. It is a necessary
// condition, not a sufficient one -- a violation proves noise, and its
// absence proves little.
//
// The sobering half of the same calibration: filtering recovers an order of
// magnitude and no more. At a fifth of a cell of surface noise, eight passes
// take div(n) from ten times the true value to about twice it. Conditional
// curvature statistics survive filtering; a curvature extremum does not
// become trustworthy, and no amount of smoothing makes it so.
//
//   amrGradStats3d.gnu.ex infile=plt is_per="1 1 0" \
//       vars="temp Y(H2)" curvOf=temp carry="HeatRelease I_R(H2)" outfile=plt_g
//
//   amrGradStats3d.gnu.ex infile=plt is_per="1 1 0" vars="temp Y(H2)" \
//       curvOf=temp gaussCurv=1 smooth_passes=2 outfile=plt_g_s2
//
#include "AMRGradient.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Reduce.H>

#include <iomanip>
#include <limits>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        std::string infile;  pp.get("infile", infile);
        std::string outfile; pp.query("outfile", outfile);
        int finestLevel = -1; pp.query("finestLevel", finestLevel);
        std::string curvOf;  pp.query("curvOf", curvOf);
        Real nfloor = 0;     pp.query("normal_floor", nfloor);
        // Negative control only; see AMRGradient.H.
        int perLevel = 0;    pp.query("per_level", perLevel);
        // A COUNT of explicit low-pass passes over the curvOf field, not a
        // smoothing time. Applied before anything is differentiated, so every
        // derived quantity is consistent with one field. See AMRGradient.H.
        int nsmooth = 0;     pp.query("smooth_passes", nsmooth);
        // Gaussian and principal curvatures alongside the mean. 3D only: a
        // curve in the plane has one principal curvature and no Gaussian one.
        int doGauss = 0;     pp.query("gaussCurv", doGauss);

        const int nv = pp.countval("vars");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nv >= 1, "vars is required");
        Vector<std::string> vars;  pp.getarr("vars", vars, 0, nv);
        Vector<std::string> carry;
        if (const int n = pp.countval("carry")) { pp.getarr("carry", carry, 0, n); }

        // Names as they appear in the input plotfile. vars is what the OUTPUT
        // components are named after, and smoothing makes the two differ.
        Vector<std::string> varSrc = vars;

        int vc = -1;
        for (int v = 0; v < nv; ++v) { if (vars[v] == curvOf) { vc = v; } }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(curvOf.empty() || vc >= 0,
                                         "curvOf must name one of vars");
        if (nsmooth > 0) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vc >= 0,
                "smooth_passes needs curvOf to say which field to filter");
            // Everything derived from the filtered field is named for it, so a
            // filtered diagnostic can never be mistaken for an unfiltered one
            // in a downstream joint fit -- which matters, because the point of
            // the pass count is to compare files across values of it. The raw
            // field rides along untouched so the filter stays auditable.
            carry.push_back(curvOf);
            vars[vc] = "smooth" + std::to_string(nsmooth) + "_" + curvOf;
            curvOf   = vars[vc];
        }
        if (doGauss && AMREX_SPACEDIM != 3) {
            Print() << "gaussCurv ignored: there is no Gaussian curvature in "
                    << AMREX_SPACEDIM << "D\n";
            doGauss = 0;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!doGauss || !curvOf.empty(),
                                         "gaussCurv needs curvOf");

        Array<int,AMREX_SPACEDIM> isper{{AMREX_D_DECL(0,0,0)}};
        Array<int,AMREX_SPACEDIM> symd {{AMREX_D_DECL(0,0,0)}};
        {
            Vector<int> t(AMREX_SPACEDIM,0);
            pp.queryarr("is_per", t, 0, AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { isper[d] = t[d]; }
            std::fill(t.begin(), t.end(), 0);
            pp.queryarr("sym_dir", t, 0, AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { symd[d] = t[d]; }
        }

        PlotFileData pf(infile);
        int finest = pf.finestLevel();
        if (finestLevel >= 0) { finest = std::min(finest, finestLevel); }
        const int nlev = finest + 1;

        RealBox rb(pf.probLo().data(), pf.probHi().data());
        Vector<Geometry> geom(nlev);
        Vector<BoxArray> grids(nlev);
        Vector<DistributionMapping> dmap(nlev);
        Vector<IntVect> rr(std::max(nlev-1,1), IntVect(AMREX_D_DECL(2,2,2)));
        for (int lev = 0; lev < nlev; ++lev) {
            grids[lev] = pf.boxArray(lev);
            dmap[lev]  = pf.DistributionMap(lev);
            geom[lev]  = Geometry(pf.probDomain(lev), &rb, pf.coordSys(), isper.data());
            if (lev < nlev-1) {
                const int r = pf.refRatio(lev);
                rr[lev] = IntVect(AMREX_D_DECL(r,r,r));
            }
        }

        const bool doCurv = !curvOf.empty();
        // layout: [vars][carry][per var: grad(DIM)+mag][per pair: cos,sin]
        //         [if curv: normal(DIM), meanCurv]
        const int nPair  = nv*(nv-1)/2;
        const int iVar   = 0;
        const int iCarry = iVar + nv;
        const int iGrad  = iCarry + int(carry.size());
        const int iPair  = iGrad + nv*(AMREX_SPACEDIM+1);
        const int iNorm  = iPair + 2*nPair;
        const int iShape = iNorm + (doCurv ? AMREX_SPACEDIM+1 : 0);
        const int nShape = 5;   // K, k1, k2, shape index, curvedness
        const int ncomp  = iShape + (doGauss ? nShape : 0);

        Vector<MultiFab> S(nlev);
        for (int lev = 0; lev < nlev; ++lev) {
            S[lev].define(grids[lev], dmap[lev], ncomp, 1);
            S[lev].setVal(0.0);
            for (int v = 0; v < nv; ++v) {
                const MultiFab src = pf.get(lev, varSrc[v]);
                MultiFab::Copy(S[lev], src, 0, iVar+v, 1, 0);
            }
            for (int c = 0; c < int(carry.size()); ++c) {
                const MultiFab src = pf.get(lev, carry[c]);
                MultiFab::Copy(S[lev], src, 0, iCarry+c, 1, 0);
            }
            S[lev].FillBoundary(geom[lev].periodicity());
        }

        // Filtering, before anything is differentiated. Curvature is a second
        // derivative, so it amplifies the cell-to-cell component of the field
        // by 1/dx^2 and reports the result as structure. The honest way to
        // find out how much of an extremum is that component is to remove it
        // and see what is left, which is a sweep over the pass count, not a
        // single run with a tolerance in it.
        if (nsmooth > 0) {
            Vector<MultiFab> phi(nlev);
            Vector<MultiFab*> pphi(nlev);
            for (int lev = 0; lev < nlev; ++lev) {
                phi[lev].define(grids[lev], dmap[lev], 1, 1);
                phi[lev].setVal(0.0);
                MultiFab::Copy(phi[lev], S[lev], iVar+vc, 0, 1, 0);
                phi[lev].FillBoundary(geom[lev].periodicity());
                pphi[lev] = &phi[lev];
            }
            amrcomplex::AMRGradient::smooth(pphi, geom, grids, dmap, isper,
                                            symd, nsmooth, 4);
            Real dmax = 0;
            for (int lev = 0; lev < nlev; ++lev) {
                MultiFab::Subtract(phi[lev], S[lev], iVar+vc, 0, 1, 0);
                dmax = std::max(dmax, phi[lev].norm0(0, 0));
                MultiFab::Add(phi[lev], S[lev], iVar+vc, 0, 1, 0);
                MultiFab::Copy(S[lev], phi[lev], 0, iVar+vc, 1, 0);
            }
            ParallelDescriptor::ReduceRealMax(dmax);
            for (int lev = 0; lev < nlev; ++lev) {
                S[lev].FillBoundary(iVar+vc, 1, geom[lev].periodicity());
            }
            // What the filter took out, in the units of the field. If this is
            // comparable to the field's own range the pass count is too high
            // and resolved structure is going with the noise.
            Print() << nsmooth << " smoothing pass(es) on " << varSrc[vc]
                    << " : largest change " << std::setprecision(6)
                    << std::scientific << dmax << std::defaultfloat << "\n";
        }

        // Gradients, one composite pass per variable.
        for (int v = 0; v < nv; ++v) {
            Vector<MultiFab> phi(nlev), g(nlev);
            Vector<MultiFab*> pphi(nlev), pg(nlev);
            for (int lev = 0; lev < nlev; ++lev) {
                phi[lev].define(grids[lev], dmap[lev], 1, 1);
                phi[lev].setVal(0.0);
                MultiFab::Copy(phi[lev], S[lev], iVar+v, 0, 1, 0);
                phi[lev].FillBoundary(geom[lev].periodicity());
                g[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
                g[lev].setVal(0.0);
                pphi[lev] = &phi[lev]; pg[lev] = &g[lev];
            }
            amrcomplex::AMRGradient::compute(pphi, pg, geom, grids, dmap, isper,
                                             symd, 4, perLevel);
            for (int lev = 0; lev < nlev; ++lev) {
                const int base = iGrad + v*(AMREX_SPACEDIM+1);
                MultiFab::Copy(S[lev], g[lev], 0, base, AMREX_SPACEDIM, 0);
                for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    auto const& a = S[lev].array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real s = 0;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            s += a(i,j,k,base+d)*a(i,j,k,base+d);
                        }
                        a(i,j,k,base+AMREX_SPACEDIM) = std::sqrt(s);
                    });
                }
            }
        }

        // Pairwise alignment.
        int ip = 0;
        for (int v = 0; v < nv; ++v) {
        for (int w = v+1; w < nv; ++w, ++ip) {
            const int bv = iGrad + v*(AMREX_SPACEDIM+1);
            const int bw = iGrad + w*(AMREX_SPACEDIM+1);
            const int oc = iPair + 2*ip;
            for (int lev = 0; lev < nlev; ++lev) {
                for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    auto const& a = S[lev].array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const Real ma = a(i,j,k,bv+AMREX_SPACEDIM);
                        const Real mb = a(i,j,k,bw+AMREX_SPACEDIM);
                        const Real den = ma*mb;
                        if (den <= Real(0.0)) { a(i,j,k,oc)=0; a(i,j,k,oc+1)=0; return; }
                        Real dot = 0;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            dot += a(i,j,k,bv+d)*a(i,j,k,bw+d);
                        }
#if (AMREX_SPACEDIM == 2)
                        const Real cx = a(i,j,k,bv+0)*a(i,j,k,bw+1)
                                      - a(i,j,k,bv+1)*a(i,j,k,bw+0);
                        const Real cr = std::abs(cx);
#else
                        const Real c0 = a(i,j,k,bv+1)*a(i,j,k,bw+2)-a(i,j,k,bv+2)*a(i,j,k,bw+1);
                        const Real c1 = a(i,j,k,bv+2)*a(i,j,k,bw+0)-a(i,j,k,bv+0)*a(i,j,k,bw+2);
                        const Real c2 = a(i,j,k,bv+0)*a(i,j,k,bw+1)-a(i,j,k,bv+1)*a(i,j,k,bw+0);
                        const Real cr = std::sqrt(c0*c0+c1*c1+c2*c2);
#endif
                        a(i,j,k,oc)   = dot/den;
                        a(i,j,k,oc+1) = cr/den;
                    });
                }
            }
        }}

        // Unit normal and mean curvature of one field, each divergence term a
        // composite pass of its own -- the same treatment as the gradient, which
        // is the whole point.
        if (doCurv) {
            const int bv = iGrad + vc*(AMREX_SPACEDIM+1);

            Real gmax = 0;
            for (int lev = 0; lev < nlev; ++lev) {
                gmax = std::max(gmax, S[lev].max(bv+AMREX_SPACEDIM, 0));
            }
            ParallelDescriptor::ReduceRealMax(gmax);
            // Default 1e-2 of the peak gradient. The first version of this used
            // 1e-8, which is no floor at all: in the uniform preheat and product
            // regions grad T is numerical noise rather than zero, the "unit"
            // normal built from it points in an arbitrary direction that changes
            // cell to cell, and its divergence is 1/dx-scale garbage. On a lean
            // H2 flame that produced curvatures of 4e5 1/m -- a 2.5 micron radius
            // on a 20 micron mesh, which is not a structure, it is a symptom.
            const Real floor = (nfloor > 0 ? nfloor : Real(1.0e-2)*gmax);

            for (int lev = 0; lev < nlev; ++lev) {
                for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    auto const& a = S[lev].array(mfi);
                    const Real fl = floor;
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        // Divide by max(|grad|, floor) rather than zeroing
                        // below the floor. Zeroing makes the normal field jump
                        // from 0 to unit length across the floor contour, and
                        // the divergence of that jump is a spurious curvature
                        // sitting on a level set of |grad| -- an artefact that
                        // looks exactly like a flame structure. This form is
                        // continuous everywhere: it is the unit normal wherever
                        // the gradient is strong, and decays smoothly to zero
                        // where there is nothing to orient.
                        const Real m = amrex::max(a(i,j,k,bv+AMREX_SPACEDIM), fl);
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            a(i,j,k,iNorm+d) = a(i,j,k,bv+d)/m;
                        }
                    });
                }
                S[lev].FillBoundary(iNorm, AMREX_SPACEDIM, geom[lev].periodicity());
                S[lev].setVal(0.0, iNorm+AMREX_SPACEDIM, 1, 0);
            }
            // How much of the domain the floor is actually governing. If this
            // is a large fraction then meanCurv is meaningful on only a thin
            // shell and should be read conditioned on gradMag, not alone.
            {
                // Volume-weighted over uncovered cells, as a device reduction:
                // the host loop this replaces was the one place in this driver
                // that pulled a component back across the bus, and it did it
                // over the whole hierarchy to produce two scalars.
                ReduceOps<ReduceOpSum,ReduceOpSum> rop;
                ReduceData<Real,Real> rdat(rop);
                using RT = typename ReduceData<Real,Real>::Type;
                for (int lev = 0; lev < nlev; ++lev) {
                    Real dV = 1;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { dV *= geom[lev].CellSize(d); }
                    iMultiFab cov(grids[lev], dmap[lev], 1, 0);
                    cov.setVal(0);
                    if (lev < nlev-1) {
                        const iMultiFab fm = makeFineMask(grids[lev], dmap[lev],
                                                          IntVect(0), grids[lev+1],
                                                          rr[lev], geom[lev].periodicity(),
                                                          0, 1);
                        iMultiFab::Copy(cov, fm, 0, 0, 1, 0);
                    }
                    const int mag = bv+AMREX_SPACEDIM;
                    const Real fl = floor;
                    for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        auto const& a = S[lev].const_array(mfi);
                        auto const& c = cov.const_array(mfi);
                        rop.eval(bx, rdat,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> RT
                        {
                            if (c(i,j,k)) { return {Real(0.0), Real(0.0)}; }
                            return {dV, (a(i,j,k,mag) < fl ? dV : Real(0.0))};
                        });
                    }
                }
                RT hv = rdat.value(rop);
                Real vall = amrex::get<0>(hv);
                Real vlo  = amrex::get<1>(hv);
                ParallelDescriptor::ReduceRealSum(vlo);
                ParallelDescriptor::ReduceRealSum(vall);
                Print() << "normal floor for " << curvOf << " : " << std::setprecision(6)
                        << std::scientific << floor
                        << "  (peak |grad| " << gmax << ")\n"
                        << "   uncovered volume below it : "
                        << std::setprecision(4) << std::fixed
                        << (vall > 0 ? Real(100.0)*vlo/vall : Real(0)) << " %\n"
                        << std::defaultfloat;
            }

            // The DIM composite passes below each produce one ROW of the
            // Jacobian of the unit normal, J_dj = d(n_d)/d(x_j). Only its
            // trace has been kept until now, and the trace is the mean
            // curvature. The rest of the matrix is the rest of the shape
            // operator and costs nothing extra to retain, because the passes
            // that produce it are already being run.
            //
            // n has unit length, so n is a null vector of J and the third
            // eigenvalue is zero: the eigenvalues are {k1, k2, 0}. Hence
            //
            //   trace(J)                  = k1 + k2
            //   (trace(J)^2 - trace(J^2))/2 = k1 k2 = K
            //
            // and the principal pair follows from the quadratic. No second
            // invariant has to be assembled by hand and no Hessian of the
            // field is formed -- the normal already carries the
            // normalisation that turns a Hessian into a shape operator.
            Vector<MultiFab> Jac;
            if (doGauss) {
                Jac.resize(nlev);
                for (int lev = 0; lev < nlev; ++lev) {
                    Jac[lev].define(grids[lev], dmap[lev],
                                    AMREX_SPACEDIM*AMREX_SPACEDIM, 0);
                    Jac[lev].setVal(0.0);
                }
            }

            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Vector<MultiFab> phi(nlev), g(nlev);
                Vector<MultiFab*> pphi(nlev), pg(nlev);
                for (int lev = 0; lev < nlev; ++lev) {
                    phi[lev].define(grids[lev], dmap[lev], 1, 1);
                    phi[lev].setVal(0.0);
                    MultiFab::Copy(phi[lev], S[lev], iNorm+d, 0, 1, 0);
                    phi[lev].FillBoundary(geom[lev].periodicity());
                    g[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
                    g[lev].setVal(0.0);
                    pphi[lev] = &phi[lev]; pg[lev] = &g[lev];
                }
                amrcomplex::AMRGradient::compute(pphi, pg, geom, grids, dmap, isper,
                                                 symd, 4, perLevel);
                for (int lev = 0; lev < nlev; ++lev) {
                    MultiFab::Add(S[lev], g[lev], d, iNorm+AMREX_SPACEDIM, 1, 0);
                    if (doGauss) {
                        MultiFab::Copy(Jac[lev], g[lev], 0, d*AMREX_SPACEDIM,
                                       AMREX_SPACEDIM, 0);
                    }
                }
            }

            if (doGauss) {
                for (int lev = 0; lev < nlev; ++lev) {
                    for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        auto const& a = S[lev].array(mfi);
                        auto const& J = Jac[lev].const_array(mfi);
                        const int iK = iShape;
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            const Real tr = a(i,j,k,iNorm+AMREX_SPACEDIM);
                            Real tr2 = 0;
                            for (int p = 0; p < AMREX_SPACEDIM; ++p) {
                                for (int q = 0; q < AMREX_SPACEDIM; ++q) {
                                    tr2 += J(i,j,k,p*AMREX_SPACEDIM+q)
                                         * J(i,j,k,q*AMREX_SPACEDIM+p);
                                }
                            }
                            const Real K    = Real(0.5)*(tr*tr - tr2);
                            const Real disc = tr*tr - Real(4.0)*K;
                            // disc >= 0 for any real surface. Where it is
                            // negative the field is not describing a surface
                            // at this scale; the root is clamped so the
                            // outputs stay finite, and the share of volume
                            // that needed clamping is reported below as the
                            // measure of how much of this field is noise.
                            const Real r = std::sqrt(amrex::max(disc, Real(0.0)));
                            const Real k1 = Real(0.5)*(tr + r);
                            const Real k2 = Real(0.5)*(tr - r);
                            a(i,j,k,iK)   = K;
                            a(i,j,k,iK+1) = k1;
                            a(i,j,k,iK+2) = k2;
                            // Koenderink's pair. The shape index says WHICH
                            // shape on a scale from -1 (cup) through 0
                            // (saddle) to +1 (cap), independent of how
                            // sharply curved it is; the curvedness carries
                            // the scale. Splitting them lets a crest and a
                            // cusp of different sharpness be compared as
                            // shapes, which a curvature magnitude cannot do.
                            const Real den = k1 - k2;
                            // At an umbilic point k1 == k2, every direction
                            // is principal and the arctan has no argument.
                            // The limit is not zero: approaching umbilicity
                            // with both curvatures positive the index tends
                            // to +1, a cap, and to -1 for a cup. A spherical
                            // crest is very nearly umbilic, so returning the
                            // saddle value 0 there -- as the obvious guard
                            // does -- would classify exactly the feature
                            // this is for as its opposite.
                            a(i,j,k,iK+3) = (den > Real(0.0))
                                ? Real(2.0/3.14159265358979323846)
                                      *std::atan((k1 + k2)/den)
                                : (tr > Real(0.0) ? Real(1.0)
                                                  : (tr < Real(0.0) ? Real(-1.0)
                                                                    : Real(0.0)));
                            a(i,j,k,iK+4) = std::sqrt(Real(0.5)*(k1*k1 + k2*k2));
                        });
                    }
                }

                // The noise meter. Measured only where the normal is a unit
                // vector to begin with: below the floor it is deliberately
                // not, so the identity the realizability test rests on does
                // not hold there and counting those cells would just restate
                // the floor fraction.
                ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum> rop;
                ReduceData<Real,Real,Real> rdat(rop);
                using RT = typename ReduceData<Real,Real,Real>::Type;
                // The violation is measured RELATIVE to the size of the two
                // terms, and what is reported is how far past the boundary a
                // cell is, not merely that it is past it.
                //
                // The reason is umbilic surface. A spherical cap has k1 == k2
                // and therefore sits exactly on kappa^2 = 4K, so the ordinary
                // truncation error of a second derivative puts about half of
                // any such region on the forbidden side. Measured on an exact
                // distance-to-a-sphere field -- umbilic everywhere, and whose
                // curvatures this build reproduces to 0.7% -- a strict test
                // calls 83% of the field unrealizable. The severity there is
                // 3e-4, which is the useful statement: the violations are the
                // size of the discretisation error, not the size of the
                // answer. Noise makes K and kappa^2 unrelated rather than
                // nearly equal, so it shows up as severity of order one.
                //
                // Hence a fraction counted only above a severity cut well
                // clear of truncation, and the mean severity beside it.
                const Real sevcut = Real(0.1);
                for (int lev = 0; lev < nlev; ++lev) {
                    Real dV = 1;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) { dV *= geom[lev].CellSize(d); }
                    iMultiFab cov(grids[lev], dmap[lev], 1, 0);
                    cov.setVal(0);
                    if (lev < nlev-1) {
                        const iMultiFab fm = makeFineMask(grids[lev], dmap[lev],
                                                          IntVect(0), grids[lev+1],
                                                          rr[lev], geom[lev].periodicity(),
                                                          0, 1);
                        iMultiFab::Copy(cov, fm, 0, 0, 1, 0);
                    }
                    const int mag = bv+AMREX_SPACEDIM;
                    const int iK  = iShape;
                    const Real fl = floor;
                    for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        auto const& a = S[lev].const_array(mfi);
                        auto const& c = cov.const_array(mfi);
                        rop.eval(bx, rdat,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> RT
                        {
                            if (c(i,j,k) || a(i,j,k,mag) < fl) {
                                return {Real(0.0), Real(0.0), Real(0.0)};
                            }
                            const Real tr = a(i,j,k,iNorm+AMREX_SPACEDIM);
                            const Real KK = a(i,j,k,iK);
                            const Real viol = Real(4.0)*KK - tr*tr;
                            const Real scale = tr*tr + Real(4.0)*std::abs(KK);
                            if (scale <= Real(0.0) || viol <= Real(0.0)) {
                                return {dV, Real(0.0), Real(0.0)};
                            }
                            // Severity in [0,1]: how far past the boundary,
                            // as a share of what the two terms are worth.
                            const Real sev = amrex::min(viol/scale, Real(1.0));
                            return {dV, (sev > sevcut ? dV : Real(0.0)), dV*sev};
                        });
                    }
                }
                RT hv = rdat.value(rop);
                Real vin  = amrex::get<0>(hv);
                Real vbad = amrex::get<1>(hv);
                Real vsev = amrex::get<2>(hv);
                ParallelDescriptor::ReduceRealSum(vin);
                ParallelDescriptor::ReduceRealSum(vbad);
                ParallelDescriptor::ReduceRealSum(vsev);
                Print() << "unrealizable curvature above the floor : "
                        << std::setprecision(4) << std::fixed
                        << (vin > 0 ? Real(100.0)*vbad/vin : Real(0))
                        << " % of volume past severity " << sevcut
                        << ",  mean severity "
                        << (vin > 0 ? vsev/vin : Real(0))
                        << std::defaultfloat << "\n"
                        << "   (kappa^2 < 4K cannot happen on a surface, so"
                        << " this is how much of the\n"
                        << "    curvature field is point-to-point noise, and"
                        << " how badly. A resolved\n"
                        << "    surface gives mean severity of order 1e-4."
                        << " Raise smooth_passes until\n"
                        << "    these stop falling, and no further.)\n";
            }
        }

        Vector<std::string> nm(ncomp);
        for (int v = 0; v < nv; ++v) { nm[iVar+v] = vars[v]; }
        for (int c = 0; c < int(carry.size()); ++c) { nm[iCarry+c] = carry[c]; }
        const char* dn[3] = {"x","y","z"};
        for (int v = 0; v < nv; ++v) {
            const int b = iGrad + v*(AMREX_SPACEDIM+1);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { nm[b+d] = "grad"+std::string(dn[d])+"_"+vars[v]; }
            nm[b+AMREX_SPACEDIM] = "gradMag_"+vars[v];
        }
        ip = 0;
        for (int v = 0; v < nv; ++v) {
        for (int w = v+1; w < nv; ++w, ++ip) {
            nm[iPair+2*ip]   = "cosTheta_"+vars[v]+"_"+vars[w];
            nm[iPair+2*ip+1] = "sinTheta_"+vars[v]+"_"+vars[w];
        }}
        if (doCurv) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { nm[iNorm+d] = "normal"+std::string(dn[d])+"_"+curvOf; }
            nm[iNorm+AMREX_SPACEDIM] = "meanCurv_"+curvOf;
        }
        if (doGauss) {
            nm[iShape+0] = "gaussCurv_"+curvOf;
            nm[iShape+1] = "curv1_"+curvOf;
            nm[iShape+2] = "curv2_"+curvOf;
            nm[iShape+3] = "shapeIndex_"+curvOf;
            nm[iShape+4] = "curvedness_"+curvOf;
        }

        if (outfile.empty()) { outfile = infile + "_g"; }
        Vector<const MultiFab*> mfp(nlev);
        Vector<int> steps(nlev, 0);
        for (int lev = 0; lev < nlev; ++lev) { mfp[lev] = &S[lev]; }
        WriteMultiLevelPlotfile(outfile, nlev, mfp, nm, geom, Real(0.0), steps, rr);
        Print() << "wrote " << outfile << " with " << ncomp << " components\n";

        // Report each component's range, because the next tool in the chain
        // needs it. amrJointFit takes clo/chi per coordinate, and for a derived
        // quantity like meanCurv_temp there is no way to guess those: pick them
        // too narrow and most of the volume falls out of range, too wide and
        // the bins are wasted on empty space. The range is over valid cells on
        // every level, covered ones included, which can only widen it.
        Print() << "\n"
                << std::setw(34) << std::left  << "   component"
                << std::setw(16) << std::right << "min"
                << std::setw(16) << std::right << "max" << "\n";
        for (int i = 0; i < ncomp; ++i) {
            Real lo =  std::numeric_limits<Real>::max();
            Real hi = -std::numeric_limits<Real>::max();
            for (int lev = 0; lev < nlev; ++lev) {
                lo = std::min(lo, S[lev].min(i, 0));
                hi = std::max(hi, S[lev].max(i, 0));
            }
            ParallelDescriptor::ReduceRealMin(lo);
            ParallelDescriptor::ReduceRealMax(hi);
            Print() << std::setw(34) << std::left  << ("   " + nm[i])
                    << std::setw(16) << std::right << std::setprecision(6)
                    << std::scientific << lo
                    << std::setw(16) << std::right << hi << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}
