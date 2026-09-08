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
        const auto rep = cx.verify();
        cx.reportVerify(rep);
    }
    amrex::Finalize();
    return 0;
}
