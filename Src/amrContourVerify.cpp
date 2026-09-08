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

        // The iso component comes first; any further names are carried on the
        // surface as sampled data, exactly as isosurface.cpp does.
        std::string var("dist");
        pp.query("isoCompName", var);
        p.varnames = {var};
        if (const int nv = pp.countval("varnames")) {
            Vector<std::string> extra(nv);
            pp.getarr("varnames", extra, 0, nv);
            for (const auto& e : extra) { if (e != var) { p.varnames.push_back(e); } }
        }

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
        // A serial test for a parallel defect. Within one rank the vertex
        // cache hides any element-dependence in the computation completely;
        // it surfaces only as two ranks producing copies that differ in the
        // last bits, which is what broke the MEF byte-identity claim.
        int checkDet = 1;
        pp.query("check_vertex_determinism", checkDet);

        amrcomplex::AMRContour ct(cx, 0, isoVal, diagRule != 0, checkDet != 0);
        if (checkDet) {
            const long mm = ct.determinismMismatches();
            Print() << "  vertex determinism: " << mm
                    << " cached vertices differ on recompute"
                    << (mm == 0 ? "  PASS\n" : "  FAIL\n");
        }

        std::string outfile;
        pp.query("outfile", outfile);
        if (!outfile.empty()) {
            char buf[72]; std::snprintf(buf, sizeof(buf), "%g", isoVal);
            ct.writeMEF(outfile, std::string(buf), p.varnames);
        }
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
        // Exactness check on the carried components. makeTestPlt writes a
        // linear field per direction; linear interpolation of a linear
        // function along an edge is exact, so each sampled value must equal
        // the vertex's own coordinate to roundoff. Anything larger means the
        // sampling is using a different parameter from the position, or the
        // wrong endpoint values.
        bool checkLinear = false;
        pp.query("check_linear_components", checkLinear);
        if (checkLinear && cx.nComp() >= 1+AMREX_SPACEDIM) {
            const auto& V = ct.vertices();
            const auto& D = ct.vertexData();
            const int nc = cx.nComp();
            Real worst = 0;
            for (std::size_t i = 0; i < V.size(); ++i) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const Real got = D[std::size_t(i)*std::size_t(nc)+std::size_t(1+d)];
                    worst = std::max(worst, std::abs(got - V[i][d]));
                }
            }
            ParallelDescriptor::ReduceRealMax(worst);
            Print() << "  carried-component exactness: worst |sampled - coord| = "
                    << worst << (worst < Real(1.0e-12) ? "  PASS\n" : "  FAIL\n");
        }
        Print() << std::endl;
    }
    amrex::Finalize();
    return 0;
}
