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

} // namespace amrcomplex
