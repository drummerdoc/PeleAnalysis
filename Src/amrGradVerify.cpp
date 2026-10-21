//
// Verify a gradient plotfile against constants that hold exactly, and report
// the error separately on the cells where the coarse-fine treatment is the only
// thing that can be wrong.
//
// The whole claim AMRGradient makes is about coarse-fine interfaces. A global
// norm hides that: the interface is a surface, so a ring of badly wrong values
// around every refined patch is a vanishing fraction of the volume and barely
// moves an L2 over the domain. So this reports three regions per check:
//
//   all   every uncovered valid cell
//   int   those not within 2 cells of a non-periodic domain face (where the
//         Neumann condition makes a linear field inexact for honest reasons)
//   cf    those in "int" that are also within 2 cells of a coarse-fine
//         interface -- on a fine level, near the edge of that level's
//         footprint; on a coarse level, near the part covered by finer cells
//
// "cf" is the measurement. For a linear field, MLMG's interface interpolation
// is exact, so the expected value there is the same roundoff as anywhere else.
// Build the operator one level at a time instead (per_level=1 in amrGradStats,
// which is what Src/curvature.cpp does) and "cf" is the only region that moves.
//
//   amrGradVerify3d.gnu.ex infile=plt_g is_per="0 0 0" \
//       check="gradx_xco grady_xco gradMag_dist" expect="1 0 1"
//
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <iomanip>

using namespace amrex;

namespace {

constexpr int RAD = 2;

// 1 on cells to keep, 0 elsewhere. Covered cells are dropped, cells near a
// non-periodic domain face are dropped from "int" and below.
void
buildMasks (const BoxArray& ba, const DistributionMapping& dm,
            const Geometry& geom, const BoxArray* finer, const IntVect& rr,
            const Array<int,AMREX_SPACEDIM>& is_per,
            Real exr, const Array<Real,AMREX_SPACEDIM>& exc,
            iMultiFab& keepAll, iMultiFab& keepInt, iMultiFab& keepCF)
{
    const Box& dom = geom.Domain();

    // uncovered
    iMultiFab cov(ba, dm, 1, RAD);
    cov.setVal(0);
    if (finer) {
        const iMultiFab fm = makeFineMask(ba, dm, IntVect(RAD), *finer, rr,
                                          geom.periodicity(), 0, 1);
        iMultiFab::Copy(cov, fm, 0, 0, 1, RAD);
    }

    // inside this level's footprint, including through periodic wraps: valid
    // cells are 1, everything else (including outside the domain) is 0.
    iMultiFab own(ba, dm, 1, RAD);
    own.setVal(0, 0, 1, RAD);
    own.setVal(1, 0, 1, 0);
    own.FillBoundary(geom.periodicity());

    keepAll.define(ba, dm, 1, 0);
    keepInt.define(ba, dm, 1, 0);
    keepCF .define(ba, dm, 1, 0);

    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> exca;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { exca[d] = exc[d]; }
    const Real exr2 = exr*exr;

    for (MFIter mfi(keepAll); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& ka = keepAll.array(mfi);
        auto const& ki = keepInt.array(mfi);
        auto const& kc = keepCF .array(mfi);
        auto const& cv = cov.const_array(mfi);
        auto const& ow = own.const_array(mfi);
        AMREX_LOOP_3D(bx, i, j, k,
        {
            const IntVect iv(AMREX_D_DECL(i,j,k));
            if (cv(i,j,k)) { ka(i,j,k)=0; ki(i,j,k)=0; kc(i,j,k)=0; continue; }
            if (exr2 > Real(0.0)) {
                Real r2 = 0;
                AMREX_D_TERM(
                    { const Real p = plo[0]+(i+Real(0.5))*dx[0]-exca[0]; r2 += p*p; },
                    { const Real p = plo[1]+(j+Real(0.5))*dx[1]-exca[1]; r2 += p*p; },
                    { const Real p = plo[2]+(k+Real(0.5))*dx[2]-exca[2]; r2 += p*p; });
                if (r2 < exr2) { ka(i,j,k)=0; ki(i,j,k)=0; kc(i,j,k)=0; continue; }
            }
            ka(i,j,k) = 1;

            bool nearDomainFace = false;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (is_per[d]) { continue; }
                if (iv[d] - dom.smallEnd(d) < RAD) { nearDomainFace = true; }
                if (dom.bigEnd(d) - iv[d] < RAD)   { nearDomainFace = true; }
            }
            ki(i,j,k) = nearDomainFace ? 0 : 1;
            if (!ki(i,j,k)) { kc(i,j,k)=0; continue; }

            // CF-adjacent: a neighbour within RAD is either outside this
            // level's footprint (fine side) or covered by finer cells
            // (coarse side).
            bool cf = false;
            for (int kk = -RAD; kk <= RAD && !cf; ++kk) {
            for (int jj = -RAD; jj <= RAD && !cf; ++jj) {
            for (int ii = -RAD; ii <= RAD && !cf; ++ii) {
                const int ni = i+ii;
                const int nj = j+AMREX_D_PICK(0,jj,jj);
                const int nk = k+AMREX_D_PICK(0,0,kk);
                if (!ow(ni,nj,nk) || cv(ni,nj,nk)) { cf = true; }
            }}}
            kc(i,j,k) = cf ? 1 : 0;
        });
    }
}

struct Acc { Real linf = 0; Real l2 = 0; Real vol = 0; Long n = 0; };

void
accumulate (const MultiFab& q, const iMultiFab& keep, Real target, Real dV,
            Acc& a)
{
    for (MFIter mfi(q); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& v = q.const_array(mfi);
        auto const& m = keep.const_array(mfi);
        AMREX_LOOP_3D(bx, i, j, k,
        {
            if (!m(i,j,k)) { continue; }
            const Real e = std::abs(v(i,j,k) - target);
            a.linf = std::max(a.linf, e);
            a.l2  += e*e*dV;
            a.vol += dV;
            a.n   += 1;
        });
    }
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        std::string infile;  pp.get("infile", infile);
        int finestLevel = -1; pp.query("finestLevel", finestLevel);

        const int nchk = pp.countval("check");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nchk >= 1, "check is required");
        Vector<std::string> chk; pp.getarr("check", chk, 0, nchk);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pp.countval("expect") == nchk,
            "expect must have one value per check entry");
        Vector<Real> expect; pp.getarr("expect", expect, 0, nchk);

        Array<int,AMREX_SPACEDIM> isper{{AMREX_D_DECL(0,0,0)}};
        {
            Vector<int> t(AMREX_SPACEDIM,0);
            pp.queryarr("is_per", t, 0, AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { isper[d] = t[d]; }
        }

        // Optional: drop cells within exclude_r of exclude_c. A signed-distance
        // field has a cone tip at the centre of its zero set, where it is not
        // differentiable at all; the O(1) error there sits in an O(dx^DIM)
        // volume, so it contributes dx^(DIM/2) to the L2 and caps the measured
        // convergence rate at 1.5 however good the scheme is. Excluding it is
        // what makes the rate a statement about the discretisation.
        Real exr = 0;
        Array<Real,AMREX_SPACEDIM> exc{{AMREX_D_DECL(0,0,0)}};
        pp.query("exclude_r", exr);
        {
            Vector<Real> t(AMREX_SPACEDIM, 0);
            if (pp.queryarr("exclude_c", t, 0, AMREX_SPACEDIM)) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { exc[d] = t[d]; }
            }
        }

        PlotFileData pf(infile);
        int finest = pf.finestLevel();
        if (finestLevel >= 0) { finest = std::min(finest, finestLevel); }
        const int nlev = finest + 1;

        RealBox rb(pf.probLo().data(), pf.probHi().data());
        Vector<Geometry> geom(nlev);
        Vector<BoxArray> grids(nlev);
        Vector<DistributionMapping> dmap(nlev);
        for (int lev = 0; lev < nlev; ++lev) {
            grids[lev] = pf.boxArray(lev);
            dmap[lev]  = pf.DistributionMap(lev);
            geom[lev]  = Geometry(pf.probDomain(lev), &rb, pf.coordSys(), isper.data());
        }

        Vector<iMultiFab> kAll(nlev), kInt(nlev), kCF(nlev);
        for (int lev = 0; lev < nlev; ++lev) {
            const BoxArray* finer = (lev < nlev-1) ? &grids[lev+1] : nullptr;
            IntVect rr(AMREX_D_DECL(2,2,2));
            if (lev < nlev-1) {
                const int r = pf.refRatio(lev);
                rr = IntVect(AMREX_D_DECL(r,r,r));
            }
            buildMasks(grids[lev], dmap[lev], geom[lev], finer, rr, isper,
                       exr, exc, kAll[lev], kInt[lev], kCF[lev]);
        }

        Print() << "# amrGradVerify " << infile << "  levels 0.." << finest << "\n";
        Print() << "# CHECK component region Linf L2 volume ncells\n";

        int bad = 0;
        for (int c = 0; c < nchk; ++c) {
            const std::string& comp = chk[c];
            const Real target = expect[c];

            Acc a[3];
            for (int lev = 0; lev < nlev; ++lev) {
                const MultiFab q = pf.get(lev, comp);
                Real dV = 1;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { dV *= pf.cellSize(lev)[d]; }
                accumulate(q, kAll[lev], target, dV, a[0]);
                accumulate(q, kInt[lev], target, dV, a[1]);
                accumulate(q, kCF [lev], target, dV, a[2]);
            }
            const char* rn[3] = {"all","int","cf"};
            for (int r = 0; r < 3; ++r) {
                ParallelDescriptor::ReduceRealMax(a[r].linf);
                ParallelDescriptor::ReduceRealSum(a[r].l2);
                ParallelDescriptor::ReduceRealSum(a[r].vol);
                ParallelDescriptor::ReduceLongSum(a[r].n);
                const Real l2 = (a[r].vol > 0) ? std::sqrt(a[r].l2/a[r].vol) : Real(0);
                Print() << "CHECK " << comp << " " << rn[r]
                        << " " << std::setprecision(12) << std::scientific << a[r].linf
                        << " " << l2
                        << " " << a[r].vol
                        << " " << a[r].n << "\n";
                if (a[r].n == 0) { ++bad; }
            }
        }
        if (bad) { Print() << "# WARNING " << bad << " empty region(s)\n"; }
    }
    amrex::Finalize();
    return 0;
}
