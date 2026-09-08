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

        Vector<int> isper(AMREX_SPACEDIM, 0);
        pp.queryarr("is_per", isper, 0, AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { p.is_per[d] = isper[d]; }

        Print() << "AMRComplex: " << p.infile << ", periodicity";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << p.is_per[d]; }
        Print() << std::endl;

        amrcomplex::AMRComplex cx(p);
        const auto rep = cx.verify();
        cx.reportVerify(rep);
    }
    amrex::Finalize();
    return 0;
}
