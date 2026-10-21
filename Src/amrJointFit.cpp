//
// Joint statistics of a response against a set of progress variables.
//
// Does heat release behave like a function of (T, Y_H2)? The honest form of
// that question is: among all functions of those two variables, how much of
// the variance of heat release does the best one explain? The best one is the
// conditional mean, so there is no fit to argue about -- only a number, and a
// table of how that number moves as you add variables and refine bins.
//
// The driver is amrJointFit.cpp, not amrJointStats.cpp: on a case-insensitive
// filesystem the latter is the SAME FILE as AMRJointStats.cpp, so whichever was
// written second wins and the other translation unit vanishes. The link then
// fails on every symbol of the library with "referenced from _main", which is
// the fingerprint of this mistake and not of a missing source in the makefile.
//
//   amrJointFit3d.gnu.ex infile=plt is_per="1 1 0" \
//       coords="temp Y(H2)" clo="300 0" chi="2000 0.03" cbins="64 64" \
//       resp="HeatRelease I_R(H2)" outfile=jt.dat outplt=plt_resid
//
#include "AMRJointStats.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        amrcomplex::AMRComplex::Params cp;
        pp.get("infile", cp.infile);
        pp.query("finestLevel", cp.finestLevel);
        pp.query("verbose", cp.verbose);

        Vector<std::string> coords, resp;
        const int nc = pp.countval("coords");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nc >= 1, "coords is required");
        pp.getarr("coords", coords, 0, nc);
        const int nr = pp.countval("resp");
        if (nr > 0) { pp.getarr("resp", resp, 0, nr); }

        // forEachNode carries field values in a fixed stack buffer, so the
        // total has a hard ceiling. Fail now rather than silently truncate.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nc + nr <= 16,
            "coords plus resp must be at most 16 components");

        cp.varnames.clear();
        for (const auto& s : coords) { cp.varnames.push_back(s); }
        for (const auto& s : resp)   { cp.varnames.push_back(s); }

        Vector<int> isper(AMREX_SPACEDIM, 0);
        pp.queryarr("is_per", isper, 0, AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { cp.is_per[d] = isper[d]; }

        amrcomplex::AMRJointStats::Params jp;
        jp.nCoord = nc;
        jp.nResp  = nr;
        jp.verbose = cp.verbose;
        AMREX_ALWAYS_ASSERT(pp.countval("clo")   == nc);
        AMREX_ALWAYS_ASSERT(pp.countval("chi")   == nc);
        pp.getarr("clo", jp.clo, 0, nc);
        pp.getarr("chi", jp.chi, 0, nc);
        jp.cbins.resize(nc, 64);
        if (const int n = pp.countval("cbins")) {
            AMREX_ALWAYS_ASSERT(n == nc);
            pp.getarr("cbins", jp.cbins, 0, nc);
        }

        std::string out;    pp.query("outfile", out);
        std::string outplt; pp.query("outplt", outplt);

        Print() << "AMRJointStats: " << cp.infile << ", periodicity";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << cp.is_per[d]; }
        // Every line that is not a table row starts with '#'. That is not
        // cosmetic: a harness that selects data rows by their leading
        // whitespace will otherwise pick up these too, and they carry a
        // different number of fields, which silently shifts every column in a
        // field-by-field comparison. That cost a debugging round trip once.
        Print() << "\n# coords:";
        for (int d = 0; d < nc; ++d) {
            Print() << " " << coords[d] << "[" << jp.clo[d] << "," << jp.chi[d]
                    << "]/" << jp.cbins[d];
        }
        Print() << "\n# resp  :";
        for (int r = 0; r < nr; ++r) { Print() << " " << resp[r]; }
        Print() << std::endl;

        amrcomplex::AMRComplex cx(cp);
        const auto cr = cx.verify();
        // Unconditional gate: expectedVolume measures the wall contribution by
        // the divergence theorem, so the reference is exact for any level
        // pattern and a non-zero residual is always a defect.
        if (!cr.passed()) {
            cx.reportVerify(cr);
            Abort("complex does not verify");
        }
        if (!cr.faceClosed()) {
            cx.reportVerify(cr);
            Abort("complex is not watertight; the volume measure is not "
                  "trustworthy and neither is anything conditioned on it");
        }

        amrcomplex::AMRJointStats js(cx, jp);

        // The volume in range is the gate on the binning itself: a large
        // dropped fraction means clo/chi are wrong, and every R^2 below is
        // then conditional on an unstated subset of the domain.
        const Real vtot = js.volumeTotal() + js.volumeDropped();
        Print() << "\n# volume in range : " << js.volumeTotal()
                << "  (" << 100.0*js.volumeTotal()/vtot << "% of "
                << vtot << ")\n";
        Print() << "# volume dropped  : " << js.volumeDropped()
                << "  (" << 100.0*js.volumeDropped()/vtot << "%)\n";

        // Explained variance, over every coordinate subset and every bin width
        // the one accumulated table supports. R^2 rises monotonically as bins
        // shrink and reaches 1 trivially when each bin holds one cell, so the
        // comparison that means anything is across subsets at MATCHED bin
        // counts -- which is why nOcc is printed beside it.
        for (int r = 0; r < nr; ++r) {
            Print() << "\n# response " << resp[r] << "\n";
            // The volume of a bin is an integer count per level times that
            // level's cell volume, so it is bit-identical at any rank count.
            // R^2 is not and cannot be: it is a ratio of floating-point sums
            // whose cross-rank reduction is a plain MPI sum, so it agrees to a
            // few ulp rather than exactly, and where R^2 is itself near zero
            // those ulp are its leading digits.
            Print() << "# R^2 is a floating-point reduction: equal across rank"
                    << " counts to a few ulp, not bit for bit.\n";
            Print() << "#" << std::string(72,'-') << "\n";
            Print() << "# " << std::setw(26) << std::left << "coords" << std::right
                    << std::setw(8)  << "coarsen"
                    << std::setw(12) << "nOcc"
                    << std::setw(16) << "R^2"
                    << std::setw(16) << "1-R^2" << "\n";
            for (unsigned mask = 1; mask < (1u<<unsigned(nc)); ++mask) {
                std::ostringstream nm;
                bool first = true;
                int  minb  = 1<<30;
                for (int d = 0; d < nc; ++d) {
                    if (mask & (1u<<unsigned(d))) {
                        if (!first) { nm << ","; }
                        nm << coords[d];
                        first = false;
                        minb = std::min(minb, jp.cbins[d]);
                    }
                }
                for (int co = 1; co <= minb; co *= 2) {
                    if (minb/co < 2 && co > 1) { break; }
                    const auto f = js.fit(r, mask, co);
                    Print() << "  " << std::setw(26) << std::left << nm.str()
                            << std::right << std::setprecision(10)
                            << std::setw(8)  << co
                            << std::setw(12) << f.nOccupied
                            << std::setw(16) << f.r2
                            << std::setw(16) << (Real(1.0)-f.r2) << "\n";
                }
            }
        }
        Print() << std::endl;

        // The bin table, at the finest binning and the full coordinate set.
        if (!out.empty() && ParallelDescriptor::IOProcessor()) {
            std::ofstream ofs(out.c_str());
            ofs << "#";
            for (int d = 0; d < nc; ++d) { ofs << " " << coords[d]; }
            ofs << " volume cells";
            for (int r = 0; r < nr; ++r) {
                ofs << " " << resp[r] << "_mean " << resp[r] << "_sd";
            }
            ofs << "\n";
            ofs << std::setprecision(10);
            for (long b = 0; b < js.nBins(); ++b) {
                const Real vb = js.volumeOfBin(b);
                if (vb <= Real(0.0)) { continue; }
                for (int d = 0; d < nc; ++d) { ofs << js.binCentre(b,d) << " "; }
                ofs << vb << " " << js.cellsInBin(b);
                for (int r = 0; r < nr; ++r) {
                    ofs << " " << js.conditionalMean(b,r)
                        << " " << std::sqrt(std::max(Real(0.0), js.conditionalVar(b,r)));
                }
                ofs << "\n";
            }
            ofs.close();
            Print() << "wrote " << out << "\n";
        }

        // The residual field, as a plotfile on the same hierarchy. psi is a
        // function on range space, so it is evaluated at every cell including
        // the covered ones -- those are excluded from the statistics but there
        // is no reason to leave holes in a field anyone will look at.
        if (!outplt.empty() && nr > 0) {
            const int ncomp = 1 + 3*nr;
            Vector<std::string> vn;
            vn.push_back("inbin");
            for (int r = 0; r < nr; ++r) {
                vn.push_back(resp[r]);
                vn.push_back(resp[r]+"_psi");
                vn.push_back(resp[r]+"_resid");
            }
            const int nlev = cx.nLevels();
            Vector<MultiFab> mf(nlev);
            Vector<const MultiFab*> mfp(nlev);
            Vector<Geometry> gm(nlev);
            Vector<int> steps(nlev, 0);
            Vector<IntVect> rr(std::max(nlev-1,1), IntVect(AMREX_D_DECL(2,2,2)));
            for (int lev = 0; lev < nlev; ++lev) {
                gm[lev] = cx.geom(lev);
                if (lev < nlev-1) { rr[lev] = cx.refRatio(lev); }
                mf[lev].define(cx.boxArray(lev), cx.distributionMap(lev), ncomp, 0);
                const MultiFab& src = cx.field(lev);
                for (MFIter mfi(mf[lev]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    auto const& a = mf[lev].array(mfi);
                    auto const& f = src.const_array(mfi);
                    AMREX_LOOP_3D(bx, i, j, k,
                    {
                        Real v[16];
                        for (int c = 0; c < nc+nr; ++c) { v[c] = f(i,j,k,c); }
                        const long b = js.binOf(v);
                        a(i,j,k,0) = (b >= 0) ? Real(1.0) : Real(0.0);
                        for (int r = 0; r < nr; ++r) {
                            const Real q   = v[nc+r];
                            const Real psi = (b >= 0) ? js.conditionalMean(b,r) : q;
                            a(i,j,k,1+3*r+0) = q;
                            a(i,j,k,1+3*r+1) = psi;
                            a(i,j,k,1+3*r+2) = q - psi;
                        }
                    });
                }
                mfp[lev] = &mf[lev];
            }
            WriteMultiLevelPlotfile(outplt, nlev, mfp, vn, gm, Real(0.0), steps, rr);
            Print() << "wrote " << outplt << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}
