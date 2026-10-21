//
// Driver: build the complex from a plotfile and run the verification gate.
//
#include "AMRComplex.H"
#include <AMReX_ParmParse.H>

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
        pp.query("legacy_ownership", p.legacyOwnership);
        pp.query("single_level_stamp", p.singleLevelStamp);

        // The face census is an argument to verify(), not a Params field: it
        // costs a std::map over every face instance on the rank and is off by
        // default. amrComplexVerify is the tool whose job is to run it, so it
        // exposes the switch. report_bad_faces, which printed the corners of
        // the first N non-manifold faces, has no counterpart in this
        // implementation and is accepted-and-ignored rather than silently
        // dropped, so an old command line says what happened.
        int faceCensus = 0;   pp.query("face_census", faceCensus);
        int reportBad  = 0;
        if (pp.query("report_bad_faces", reportBad) && reportBad > 0) {
            faceCensus = 1;
            amrex::Print() << "NOTE report_bad_faces is not available in this "
                           << "build; enabling face_census instead.\n";
        }

        // Field components to carry on the nodes, by plotfile variable name.
        if (const int nv = pp.countval("varnames")) {
            p.varnames.resize(nv);
            pp.getarr("varnames", p.varnames, 0, nv);
        }

        Vector<int> isper(AMREX_SPACEDIM, 0);
        pp.queryarr("is_per", isper, 0, AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { p.is_per[d] = isper[d]; }

        Print() << "AMRComplex: " << p.infile << ", periodicity";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << p.is_per[d]; }
        if (!p.varnames.empty()) {
            Print() << "AMRComplex: fields";
            for (const auto& v : p.varnames) { Print() << " " << v; }
            Print() << std::endl;
        }
        Print() << std::endl;

        amrcomplex::AMRComplex cx(p);
        const auto rep = cx.verify(faceCensus != 0);
        cx.reportVerify(rep);
    }
    amrex::Finalize();
    return 0;
}
