//
// Driver: build the complex, contour it, and check the surface.
//
#include "AMRContour.H"
#include <AMReX_ParmParse.H>
#include <cmath>

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

        std::string var("dist");
        pp.query("isoCompName", var);
        p.varnames = {var};

        Real isoVal = 0.0;
        pp.query("isoVal", isoVal);

        Vector<int> isper(AMREX_SPACEDIM, 0);
        pp.queryarr("is_per", isper, 0, AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { p.is_per[d] = isper[d]; }

        amrcomplex::AMRComplex cx(p);
        const auto cr = cx.verify();
        if (!cr.passed()) {
            cx.reportVerify(cr);
            Abort("complex does not verify; not contouring it");
        }

        // diagonal_rule=0 selects the deliberately inconsistent subdivision,
        // as a negative control on the closed-surface check.
        int diagRule = 1;
        pp.query("diagonal_rule", diagRule);
        if (diagRule == 0) {
#if (AMREX_SPACEDIM == 2)
            Print() << "NOTE: diagonal_rule has no effect in 2D. Two cells "
                    << "share an EDGE, which has no diagonal, so the "
                    << "consistency question this control exists to test "
                    << "does not arise below 3D. The control passing here is "
                    << "vacuous, not evidence.\n";
#else
            Print() << "WARNING: using the INCONSISTENT diagonal rule "
                    << "(negative control)\n";
#endif
        }
        amrcomplex::AMRContour ct(cx, 0, isoVal, diagRule != 0);
        const auto sr = ct.verify();
        ct.reportVerify(sr);

        // Analytic reference. makeTestPlt writes a signed distance to a sphere
        // of radius 0.3, so the zero isosurface has area 4*pi*r^2 exactly. The
        // extracted area is piecewise linear on tets and so converges from
        // below; what matters is that the error falls with resolution, not
        // that it is zero.
        Real radius = -1.0;
        pp.query("sphereRadius", radius);
        if (radius > 0) {
#if (AMREX_SPACEDIM == 2)
            const Real exact = Real(2.0)*Real(M_PI)*radius;      // circumference
#else
            const Real exact = Real(4.0)*Real(M_PI)*radius*radius;
#endif
            Print() << "  exact measure  : " << exact << "\n";
            Print() << "  relative error : " << (sr.area-exact)/exact << "\n";
        }
        Print() << std::endl;
    }
    amrex::Finalize();
    return 0;
}
