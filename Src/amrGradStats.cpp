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
//   amrGradStats3d.gnu.ex infile=plt is_per="1 1 0" \
//       vars="temp Y(H2)" curvOf=temp carry="HeatRelease I_R(H2)" outfile=plt_g
//
#include "AMRGradient.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>

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

        const int nv = pp.countval("vars");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nv >= 1, "vars is required");
        Vector<std::string> vars;  pp.getarr("vars", vars, 0, nv);
        Vector<std::string> carry;
        if (const int n = pp.countval("carry")) { pp.getarr("carry", carry, 0, n); }

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
        const int ncomp  = iNorm + (doCurv ? AMREX_SPACEDIM+1 : 0);

        Vector<MultiFab> S(nlev);
        for (int lev = 0; lev < nlev; ++lev) {
            S[lev].define(grids[lev], dmap[lev], ncomp, 1);
            S[lev].setVal(0.0);
            for (int v = 0; v < nv; ++v) {
                const MultiFab src = pf.get(lev, vars[v]);
                MultiFab::Copy(S[lev], src, 0, iVar+v, 1, 0);
            }
            for (int c = 0; c < int(carry.size()); ++c) {
                const MultiFab src = pf.get(lev, carry[c]);
                MultiFab::Copy(S[lev], src, 0, iCarry+c, 1, 0);
            }
            S[lev].FillBoundary(geom[lev].periodicity());
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
            int vc = -1;
            for (int v = 0; v < nv; ++v) { if (vars[v] == curvOf) { vc = v; } }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vc >= 0, "curvOf must name one of vars");
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
                Real vlo = 0, vall = 0;
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
                    for (MFIter mfi(S[lev]); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        auto const& a = S[lev].const_array(mfi);
                        auto const& c = cov.const_array(mfi);
                        AMREX_LOOP_3D(bx, i, j, k,
                        {
                            if (c(i,j,k)) { continue; }
                            vall += dV;
                            if (a(i,j,k,bv+AMREX_SPACEDIM) < floor) { vlo += dV; }
                        });
                    }
                }
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
                }
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
