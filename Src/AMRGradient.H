#ifndef AMRGRADIENT_H_
#define AMRGRADIENT_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

namespace amrcomplex {

//
// Cell-centred gradients on an AMR hierarchy, with the coarse-fine interface
// treated the way the solver treats it.
//
// Why not a difference stencil. A centred difference on a level's own ghost
// data is wrong at a coarse-fine interface in a specific way: the ghost cell is
// an interpolation of coarse data, so the stencil silently mixes two
// resolutions and the result is neither the coarse nor the fine gradient. The
// error sits exactly on the interface, which is a surface, so it does not
// converge away -- it just gets thinner. For a quantity like |grad T| that is
// then used to classify the flame, a ring of wrong values around every refined
// patch is not an acceptable artefact.
//
// What this does instead is what Src/grad.cpp does, which is the right answer
// and which this is lifted from: build a COMPOSITE MLPoisson over the whole
// hierarchy, apply it, and take the face-centred fluxes the operator computes.
// MLMG's composite apply is precisely the thing that reconciles the two sides
// of a coarse-fine interface -- it is the same reconciliation the solver uses,
// so the gradient is consistent with the discretisation the data came from.
// The fluxes are then averaged to cell centres.
//
// Two details that are easy to get wrong and are both load-bearing:
//
//   mlmg.apply() must be called BEFORE getFluxes. The coarse-fine correction is
//   established by the apply; getFluxes on its own returns fluxes that have not
//   seen it. Src/grad.cpp carries a one-line comment to this effect and it is
//   the whole reason that call is there.
//
//   MLPoisson's flux is -grad, not grad. The sign is flipped here once, so no
//   caller has to remember.
//
// Note on Src/curvature.cpp: it computes all the derived flame quantities this
// is meant to feed -- curvature, normals, strain -- but it builds its operators
// one level at a time, inside a loop over levels, which is a single-level
// operator per level and therefore does NOT get the composite treatment. Its
// derived fields carry the interface artefact described above. That is the gap
// this fills.
//
class AMRGradient
{
public:

    // phi: one component per level, with at least one ghost cell, already
    // filled on the valid region. is_per is the periodicity; a non-periodic
    // direction gets a homogeneous Neumann condition unless sym_dir marks it
    // as a symmetry plane, in which case it is reflect_odd -- the same choice
    // Src/grad.cpp makes, so results are comparable with it.
    //
    // grad is filled with AMREX_SPACEDIM components per level, cell-centred.
    //
    // perLevel is a NEGATIVE CONTROL, not an option to use: set to 1 it builds
    // one single-level operator per level, which is what Src/curvature.cpp
    // does. It exists so a test can demonstrate that the coarse-fine gate it
    // claims to check actually fires -- a gate nothing can fail is not a gate.
    static void compute (const amrex::Vector<amrex::MultiFab*>& phi,
                         const amrex::Vector<amrex::MultiFab*>& grad,
                         const amrex::Vector<amrex::Geometry>& geom,
                         const amrex::Vector<amrex::BoxArray>& grids,
                         const amrex::Vector<amrex::DistributionMapping>& dmap,
                         const amrex::Array<int,AMREX_SPACEDIM>& is_per,
                         const amrex::Array<int,AMREX_SPACEDIM>& sym_dir,
                         int maxOrder = 4,
                         int perLevel = 0);

    // A counted number of explicit low-pass passes over phi, in place.
    //
    // This is not a diffusion solve run to a time. The time-based form asks
    // the user for a smoothing time, whose effect depends on the mesh
    // spacing: the same number filters a level-3 field differently from a
    // level-1 one, and no value of it means "remove the point-to-point
    // component and nothing else". A pass count does mean that, because a
    // pass is the same operation however fine the mesh is.
    //
    // One pass is  phi <- phi + c h^2 lap(phi)  with c = 1/(4*DIM) and h the
    // FINEST spacing in the hierarchy -- one amplitude everywhere, not one
    // per level. c is the value that annihilates the cell-to-cell
    // oscillation, the mode alternating sign in every direction at once,
    // exactly. A wave resolved by n cells survives each pass with amplitude
    // 1 - (pi/n)^2/DIM, so at ten cells per wavelength a pass costs 3% in 3D
    // while the Nyquist component costs everything. That is the design
    // intent: aim the filter at what the discretisation cannot represent,
    // not at the structure it can.
    //
    // A per-level amplitude is the tempting choice and it is wrong; the
    // reasoning, and the measurement that settles it, are in the .cpp.
    //
    // The Laplacian is the composite one, so the filter treats the
    // coarse-fine interface the way the gradient does and does not deposit a
    // ring of its own around every refined patch.
    //
    // What it cannot do: recover a surface from cell-scale wrinkling.
    // Measured on a distance function perturbed by a fifth of a cell, the
    // curvature starts out an order of magnitude too large and eight passes
    // bring it back only to about twice the true value. Filtering makes
    // conditional curvature statistics usable; it does not make an extremum
    // trustworthy, and nothing does.
    //
    // phi needs one ghost cell; it is left filled on exit.
    static void smooth (const amrex::Vector<amrex::MultiFab*>& phi,
                        const amrex::Vector<amrex::Geometry>& geom,
                        const amrex::Vector<amrex::BoxArray>& grids,
                        const amrex::Vector<amrex::DistributionMapping>& dmap,
                        const amrex::Array<int,AMREX_SPACEDIM>& is_per,
                        const amrex::Array<int,AMREX_SPACEDIM>& sym_dir,
                        int npass,
                        int maxOrder = 4);
};

} // namespace amrcomplex

#endif
