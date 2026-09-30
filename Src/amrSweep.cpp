//
// Threshold sweep: the smallest thing that produces a scientific result.
//
// For a candidate field, walk the threshold and report at each value how many
// connected regions the threshold picks out, how their sizes are distributed,
// and how much isosurface there is. Run it on several candidate surrogates and
// the tables are directly comparable, which is the point -- what a "cell" is,
// is exactly what is not yet known.
//
// What to look for: a plateau in N(v), or a gap in the size distribution.
// Either is a threshold-independent scale and so a defensible cell definition.
// A featureless continuum is also an answer, and a more interesting one.
//
// This is a coarse sampling of a merge tree. When the merge tree is built it
// must reproduce this N(v) exactly at every sampled threshold.
//
//   amrSweep3d.gnu.ex infile=plt isoCompName=HeatRelease \
//       vlo=1e8 vhi=1e10 nval=40 log=1 is_per="1 1 0"
//
#include "AMRSegment.H"
#include "AMRContour.H"

#include <AMReX_ParmParse.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        amrcomplex::AMRComplex::Params p;
        pp.get("infile", p.infile);
        pp.query("finestLevel", p.finestLevel);
        pp.query("verbose", p.verbose);

        std::string var;
        pp.get("isoCompName", var);
        p.varnames = {var};

        Vector<int> isper(AMREX_SPACEDIM, 0);
        pp.queryarr("is_per", isper, 0, AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { p.is_per[d] = isper[d]; }

        Real vlo = 0, vhi = 1;
        pp.get("vlo", vlo);
        pp.get("vhi", vhi);
        int nval = 20;      pp.query("nval", nval);
        int logspace = 0;   pp.query("log", logspace);
        int above = 1;      pp.query("above", above);
        int doArea = 1;     pp.query("area", doArea);
        // Components below this volume are counted separately rather than
        // dropped: at a low threshold most of them are single cells, and how
        // many there are is itself informative about the noise floor.
        Real vmin = 0;      pp.query("min_volume", vmin);

        std::string out;    pp.query("outfile", out);

        amrcomplex::AMRComplex cx(p);
        const auto cr = cx.verify();
        // The volume residual is a gate only where the reference is exact --
        // see AMRComplex::expectedVolume. On a hierarchy that refines up to a
        // non-periodic wall the reference is a bound, not an identity, so
        // refusing to run on a non-zero residual would refuse every real
        // plotfile. Report it, say what to do about it, and carry on.
        if (!cr.passed()) {
            cx.reportVerify(cr);
            if (cr.referenceExact) { Abort("complex does not verify"); }
            Print() << "WARNING: the volume reference is a bound, not an "
                    << "identity, because a non-periodic wall carries more than\n"
                    << "         one level. Continuing. Re-run amrComplexVerify "
                    << "with is_per=1 in every\n         direction to gate this "
                    << "hierarchy sharply.\n";
        }

        Print() << "\n# sweep of " << var << " in " << p.infile
                << "   (" << (above ? "above" : "below") << " threshold)\n";
        Print() << "#" << std::string(94,'-') << "\n";
        Print() << "#" << "     threshold"
                << "   nComp" << "  nKept"
                << "        volume" << "     largestVol" << "      medianVol"
                << "          area\n";

        std::ofstream ofs;
        if (!out.empty() && ParallelDescriptor::IOProcessor()) {
            ofs.open(out.c_str());
            ofs << "# threshold nComp nKept volume largestVol medianVol area\n";
        }

        for (int i = 0; i < nval; ++i) {
            const Real f = (nval == 1) ? Real(0.0) : Real(i)/Real(nval-1);
            const Real v = logspace
                ? std::exp(std::log(vlo) + f*(std::log(vhi)-std::log(vlo)))
                : vlo + f*(vhi-vlo);

            amrcomplex::AMRSegment seg(cx, 0, v, above != 0);

            std::vector<Real> vols;
            vols.reserve(seg.components().size());
            for (const auto& c : seg.components()) { vols.push_back(c.volume); }
            std::sort(vols.begin(), vols.end());

            long nKept = 0;
            for (Real x : vols) { if (x >= vmin) { ++nKept; } }
            const Real largest = vols.empty() ? Real(0) : vols.back();
            const Real median  = vols.empty() ? Real(0) : vols[vols.size()/2];

            Real area = 0;
            if (doArea) {
                amrcomplex::AMRContour ct(cx, 0, v);
                area = ct.verify().area;
            }

            Print() << " " << std::setw(14) << v
                    << " " << std::setw(7)  << seg.nComponents()
                    << " " << std::setw(6)  << nKept
                    << " " << std::setw(13) << seg.totalVolume()
                    << " " << std::setw(14) << largest
                    << " " << std::setw(14) << median
                    << " " << std::setw(13) << area << "\n";
            if (ofs.is_open()) {
                ofs << v << " " << seg.nComponents() << " " << nKept << " "
                    << seg.totalVolume() << " " << largest << " " << median
                    << " " << area << "\n";
            }
        }
        Print() << std::endl;
        if (ofs.is_open()) { ofs.close(); }
    }
    amrex::Finalize();
    return 0;
}
