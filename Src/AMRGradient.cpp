#include "AMRGradient.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

namespace amrcomplex {

void
AMRGradient::compute (const Vector<MultiFab*>& phi,
                      const Vector<MultiFab*>& grad,
                      const Vector<Geometry>& geom,
                      const Vector<BoxArray>& grids,
                      const Vector<DistributionMapping>& dmap,
                      const Array<int,AMREX_SPACEDIM>& is_per,
                      const Array<int,AMREX_SPACEDIM>& sym_dir,
                      int maxOrder,
                      int perLevel)
{
    const int nlev = int(geom.size());
    AMREX_ALWAYS_ASSERT(int(phi.size())  == nlev);
    AMREX_ALWAYS_ASSERT(int(grad.size()) == nlev);

    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    // No coarsening: the operator is never solved with, only applied.
    info.setMaxCoarseningLevel(0);

    std::array<LinOpBCType,AMREX_SPACEDIM> lo_bc, hi_bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (is_per[d]) {
            lo_bc[d] = hi_bc[d] = LinOpBCType::Periodic;
        } else if (sym_dir[d]) {
            lo_bc[d] = hi_bc[d] = LinOpBCType::reflect_odd;
        } else {
            lo_bc[d] = hi_bc[d] = LinOpBCType::Neumann;
        }
    }

    // The group of levels handled by one operator: the whole hierarchy for the
    // composite path, one level at a time for the negative control.
    const int ngroup = perLevel ? nlev : 1;
    for (int g = 0; g < ngroup; ++g) {
        const int lev0 = perLevel ? g : 0;
        const int nl   = perLevel ? 1 : nlev;

        Vector<Geometry>            lgeom (geom .begin()+lev0, geom .begin()+lev0+nl);
        Vector<BoxArray>            lgrids(grids.begin()+lev0, grids.begin()+lev0+nl);
        Vector<DistributionMapping> ldmap (dmap .begin()+lev0, dmap .begin()+lev0+nl);
        Vector<MultiFab*>           lphi  (phi  .begin()+lev0, phi  .begin()+lev0+nl);

        MLPoisson poisson(lgeom, lgrids, ldmap, info);
        poisson.setMaxOrder(maxOrder);
        poisson.setDomainBC(lo_bc, hi_bc);
        for (int l = 0; l < nl; ++l) { poisson.setLevelBC(l, lphi[l]); }

        Vector<Array<MultiFab,AMREX_SPACEDIM>> face(nl);
        Vector<MultiFab> lap(nl);
        for (int l = 0; l < nl; ++l) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                face[l][d].define(amrex::convert(lgrids[l],
                                                 IntVect::TheDimensionVector(d)),
                                  ldmap[l], 1, 0);
            }
            lap[l].define(lgrids[l], ldmap[l], 1, 1);
        }

        MLMG mlmg(poisson);
        // The apply is what establishes the coarse-fine correction. getFluxes
        // without it returns fluxes that have not seen the interface.
        mlmg.apply(GetVecOfPtrs(lap), lphi);
        mlmg.getFluxes(GetVecOfArrOfPtrs(face), lphi, MLMG::Location::FaceCenter);

        for (int l = 0; l < nl; ++l) {
            average_face_to_cellcenter(*grad[lev0+l], 0,
                                       GetArrOfConstPtrs(face[l]));
            // MLPoisson's flux is -grad.
            grad[lev0+l]->mult(Real(-1.0), 0, AMREX_SPACEDIM, 0);
        }
    }
}

void
AMRGradient::smooth (const Vector<MultiFab*>& phi,
                     const Vector<Geometry>& geom,
                     const Vector<BoxArray>& grids,
                     const Vector<DistributionMapping>& dmap,
                     const Array<int,AMREX_SPACEDIM>& is_per,
                     const Array<int,AMREX_SPACEDIM>& sym_dir,
                     int npass,
                     int maxOrder)
{
    if (npass <= 0) { return; }

    const int nlev = int(geom.size());
    AMREX_ALWAYS_ASSERT(int(phi.size()) == nlev);

    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);

    std::array<LinOpBCType,AMREX_SPACEDIM> lo_bc, hi_bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (is_per[d]) {
            lo_bc[d] = hi_bc[d] = LinOpBCType::Periodic;
        } else if (sym_dir[d]) {
            lo_bc[d] = hi_bc[d] = LinOpBCType::reflect_odd;
        } else {
            lo_bc[d] = hi_bc[d] = LinOpBCType::Neumann;
        }
    }

    // ONE amplitude for the whole hierarchy, set by the finest spacing
    // present. This is the detail that makes the filter usable at all.
    //
    // The obvious choice is a per-level amplitude -- each level's own
    // Nyquist-killer, c*dx_l^2 -- so that every level has its cell-to-cell
    // component removed exactly. That choice is wrong, and wrong in the
    // specific way this whole file exists to avoid: the correction applied on
    // the coarse side of an interface is then four times the correction on
    // the fine side, which puts a STEP in the field along every coarse-fine
    // interface, and the curvature of a step is a ring artefact around every
    // refined patch. It was measured: on an exact distance-to-a-sphere field
    // whose curvatures are otherwise reproduced to 0.4%, one per-level pass
    // took the worst-cell error to 4% and three passes to 10%, all of it on
    // the interface, while the exact answer is unchanged by smoothing a
    // radially symmetric field at all.
    //
    // A single amplitude is a single physical filter width, so nothing is
    // discontinuous anywhere. Taking it from the finest spacing makes the
    // filter exactly the Nyquist-killer there -- which is where the mesh is
    // and where the noise worth removing lives -- and leaves coarser levels
    // damped by (dx_fine/dx_l)^2 of that, which is weaker but stable at every
    // level and, on a properly nested hierarchy, applied to regions nobody is
    // extracting curvature from.
    //
    // The smallest spacing within a level is the one that counts when cells
    // are not cubic: with dx_min in the numerator no direction's contribution
    // exceeds its isotropic value, so a pass can never amplify.
    Real hmin = geom[0].CellSize(0);
    for (int l = 0; l < nlev; ++l) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            hmin = std::min(hmin, geom[l].CellSize(d));
        }
    }
    const Real fac = hmin*hmin/Real(4*AMREX_SPACEDIM);

    Vector<MultiFab> lap(nlev);
    for (int l = 0; l < nlev; ++l) {
        lap[l].define(grids[l], dmap[l], 1, 1);
    }

    for (int p = 0; p < npass; ++p) {
        // Rebuilt every pass: setLevelBC captures the field as it stands, and
        // the field is what the pass changes.
        MLPoisson poisson(geom, grids, dmap, info);
        poisson.setMaxOrder(maxOrder);
        poisson.setDomainBC(lo_bc, hi_bc);
        for (int l = 0; l < nlev; ++l) { poisson.setLevelBC(l, phi[l]); }

        MLMG mlmg(poisson);
        mlmg.apply(GetVecOfPtrs(lap), phi);

        for (int l = 0; l < nlev; ++l) {
            MultiFab::Saxpy(*phi[l], fac, lap[l], 0, 0, 1, 0);
            phi[l]->FillBoundary(geom[l].periodicity());
        }
    }
}

} // namespace amrcomplex
