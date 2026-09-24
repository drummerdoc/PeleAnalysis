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
// vlo and vhi default to the field's range over the hierarchy, which is
// printed either way. Timings are max over ranks; memory is peak resident
// set, max per rank and summed over ranks -- the numbers that decide whether
// a production snapshot fits.
//
#include "AMRSegment.H"
#include "AMRContour.H"

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace amrex;

namespace {

// Peak resident set of this process in GB, from /proc; 0 where unavailable.
Real peakRssGB ()
{
    std::ifstream f("/proc/self/status");
    std::string key;
    while (f >> key) {
        if (key == "VmHWM:") { long kb = 0; f >> kb; return Real(kb) / (1024.0*1024.0); }
        std::string rest; std::getline(f, rest);
    }
    return 0;
}

void reportMemory (const char* when)
{
    Real mx = peakRssGB(), sum = mx;
    ParallelDescriptor::ReduceRealMax(mx);
    ParallelDescriptor::ReduceRealSum(sum);
    Print() << "# peak memory " << when << ": " << mx << " GB max/rank, "
            << sum << " GB total\n";
}

Real elapsedMax (Real t0)
{
    Real t = ParallelDescriptor::second() - t0;
    ParallelDescriptor::ReduceRealMax(t);
    return t;
}

} // namespace

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

        int nval = 20;      pp.query("nval", nval);
        int logspace = 0;   pp.query("log", logspace);
        int above = 1;      pp.query("above", above);
        int doArea = 1;     pp.query("area", doArea);
        // Components below this volume are counted separately rather than
        // dropped: at a low threshold most of them are single cells, and how
        // many there are is itself informative about the noise floor.
        Real vmin = 0;      pp.query("min_volume", vmin);

        std::string out;    pp.query("outfile", out);

        Real t0 = ParallelDescriptor::second();
        amrcomplex::AMRComplex cx(p);
        const Real tBuild = elapsedMax(t0);
        t0 = ParallelDescriptor::second();
        const auto cr = cx.verify(false);   // census off: a map over every face
        if (!cr.passed()) { cx.reportVerify(cr); Abort("complex does not verify"); }
        const Real tVerify = elapsedMax(t0);

        Print() << "# complex: " << cx.nLevels() << " levels, " << cr.nNodes
                << " nodes, " << cr.nElements << " elements\n";
        Print() << "# time build " << tBuild << " s, verify " << tVerify << " s\n";
        reportMemory("after build");

        // Range over every level's nodes. COVERED nodes carry the averaged-down
        // value, which lies inside the range of the fine data, so this bounds
        // the valid range without a mask.
        Real fmin = cx.field(0).min(0), fmax = cx.field(0).max(0);
        for (int lev = 1; lev < cx.nLevels(); ++lev) {
            fmin = std::min(fmin, cx.field(lev).min(0));
            fmax = std::max(fmax, cx.field(lev).max(0));
        }
        Print() << "# " << var << " range [" << fmin << ", " << fmax << "]\n";

        Real vlo = fmin, vhi = fmax;
        pp.query("vlo", vlo);
        pp.query("vhi", vhi);
        if (logspace && !(vlo > 0 && vhi > 0)) {
            Abort("log=1 needs vlo and vhi positive; set them explicitly");
        }

        Print() << "\n# sweep of " << var << " in " << p.infile
                << "   (" << (above ? "above" : "below") << " threshold)\n";
        Print() << "#" << std::string(94,'-') << "\n";
        Print() << "#" << "     threshold"
                << "   nComp" << "  nKept"
                << "        volume" << "     largestVol" << "      medianVol"
                << "          area" << "    exchanged" << "   tSeg" << "  tArea\n";

        std::ofstream ofs;
        if (!out.empty() && ParallelDescriptor::IOProcessor()) {
            ofs.open(out.c_str());
            ofs << "# " << var << " in " << p.infile << ", "
                << (above ? "above" : "below") << " threshold, range ["
                << fmin << ", " << fmax << "]\n";
            ofs << "# threshold nComp nKept volume largestVol medianVol area\n";
            ofs << std::setprecision(10);
        }

        for (int i = 0; i < nval; ++i) {
            const Real f = (nval == 1) ? Real(0.0) : Real(i)/Real(nval-1);
            const Real v = logspace
                ? std::exp(std::log(vlo) + f*(std::log(vhi)-std::log(vlo)))
                : vlo + f*(vhi-vlo);

            t0 = ParallelDescriptor::second();
            amrcomplex::AMRSegment seg(cx, 0, v, above != 0);
            const Real tSeg = elapsedMax(t0);

            std::vector<Real> vols;
            vols.reserve(seg.components().size());
            for (const auto& c : seg.components()) { vols.push_back(c.volume); }
            std::sort(vols.begin(), vols.end());

            long nKept = 0;
            for (Real x : vols) { if (x >= vmin) { ++nKept; } }
            const Real largest = vols.empty() ? Real(0) : vols.back();
            const Real median  = vols.empty() ? Real(0) : vols[vols.size()/2];

            Real area = 0;
            t0 = ParallelDescriptor::second();
            if (doArea) {
                amrcomplex::AMRContour ct(cx, 0, v);
                area = ct.verify().area;
            }
            const Real tArea = elapsedMax(t0);

            Print() << " " << std::setw(14) << v
                    << " " << std::setw(7)  << seg.nComponents()
                    << " " << std::setw(6)  << nKept
                    << " " << std::setw(13) << seg.totalVolume()
                    << " " << std::setw(14) << largest
                    << " " << std::setw(14) << median
                    << " " << std::setw(13) << area
                    << " " << std::setw(12) << seg.nExchanged()
                    << " " << std::setw(6) << std::setprecision(3) << tSeg
                    << " " << std::setw(6) << tArea << std::setprecision(6) << "\n";
            if (ofs.is_open()) {
                ofs << v << " " << seg.nComponents() << " " << nKept << " "
                    << seg.totalVolume() << " " << largest << " " << median
                    << " " << area << "\n";
            }
        }
        Print() << std::endl;
        reportMemory("after sweep");
        if (ofs.is_open()) { ofs.close(); }
    }
    amrex::Finalize();
    return 0;
}
