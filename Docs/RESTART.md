# Restarting this work in a new session

One document. It used to be two, written by two parallel lines of work on two
machines and concatenated at the merge; each half had its own "paste this to
start", its own conventions list, and its own stale "next" section, and the
first half sent the reader to `Docs/reeb-graph-assessment.md`, which no longer
exists. This replaces both.

Origin session: `https://claude.ai/code/session_01D4f1Detq1YmcvLsYuCywQA`, also
in the `Claude-Session:` trailer of every commit on `amr-complex`. Nothing of
value lives in that session. Everything is committed. A new session needs this
file, the repo, and whichever of the two notes under `Docs/` the task touches.

Technical documents in this repo are written natively in LaTeX, never generated
from Markdown. `Docs/design/` is the design note for the complex, contour and
segmentation machinery; `Docs/ThermodiffusiveLens.tex` is the science note for
the lens analysis. Both build with `make` in their directory. This file is the
only Markdown that is maintained.

---

## Paste this to start

> I'm continuing work in `~/src/PeleAnalysis`, branch `amr-complex`. Read
> `Docs/RESTART.md` first. Then, depending on the task: for the complex,
> contour or segmentation machinery read the design note sources under
> `Docs/design/`, starting with `sec-problem.tex` and `sec-roadmap.tex`; for
> the gradient, curvature, joint-statistics or lens analysis read
> `Docs/ThermodiffusiveLens.tex`, sections 8 and 9 for what the tools do and
> how they are used. Do not re-derive the numbers in those notes; they are
> settled and several of them were corrected in place more than once. Then
> tell me what you think the next step is before writing any code.

---

## Where things are

| | |
|---|---|
| git working copy | the Mac, `~/src/PeleAnalysis`, branch `amr-complex`, remote `origin` = `drummerdoc/PeleAnalysis` |
| Mac toolchain | `COMP=llvm USE_MPI=TRUE`; executables are `*3d.llvm.MPI.ex` |
| compute and data | `peart`, Linux, `~/src/PeleAnalysis`, default `gnu` toolchain |
| AMReX | `../../amrex` relative to the repo, which is the default and correct on both |
| production plotfile | `peart:/media/datastore/marcusd/P8atm/plt29791` |

A tool picks its executable as the NEWEST match, so a stale serial build left
in `Src/` will be preferred over a maintained MPI one. Build in `Src/`, or not
in `Src/` at all.

## The case

Lean premixed hydrogen, 8 atm, Li and Dryer mechanism. Nine species — H2, O2,
H2O, H, O, OH, HO2, H2O2, N2 — with the production rate of each, all in the
plotfile. Unburnt state `T_u = 693 K`, `Y_u(H2) = 1.325792e-2`, which is
roughly `phi = 0.46`. Periodicity `is_per="1 1 0"`: walls in z. `plt29791` has
four levels, 373 M cells, 98 GB, with 5213 boxes on level 1 filling 64% of it
and levels 1 and 2 hard against the walls.

## What this is

Post-processing for AMReX plotfiles, aimed at characterising the cellular
structures of the thermodiffusive instability in lean premixed hydrogen
flames. The open question — what a flame *cell* actually is — is why the
tooling compares candidate definitions rather than computing one, and why the
measures are built to be robust to how a cell is defined.
`Docs/design/sec-problem.tex` has the full statement.

## What is built and verified

All suites pass in 2D and 3D at 1 to 8 ranks. The complex is verified exactly
on all four levels of real PeleLMeX flame sheets.

| file | what |
|---|---|
| `Src/AMRComplex.{H,cpp}` | conforming multilevel cell complex; node identity is `(level, IntVect)` |
| `Src/AMRContour.{H,cpp}` | marching simplices, oriented, carried components, MEF output |
| `Src/AMRSegment.{H,cpp}` | connected components of a threshold set |
| `Src/AMRJointStats.{H,cpp}` | volume-weighted joint statistics in range space |
| `Src/AMRParallel.{H,cpp}` | the only place MPI calls appear |
| `Src/AMRGradient.{H,cpp}` | composite coarse-fine-aware gradient, and the counted low-pass filter |
| `Src/amrComplexVerify.cpp` | volume closure, face parity, checksums |
| `Src/amrContourVerify.cpp` | surface area and closure against analytic references |
| `Src/amrSweep.cpp`, `amrLabel.cpp` | component counts versus threshold; labelled output |
| `Src/amrJointFit.cpp` | conditional means and explained variance against coordinate subsets |
| `Src/amrGradStats.cpp` | gradients, alignment, curvatures, filtering; writes a plotfile |
| `Src/amrGradVerify.cpp` | gradient error by region, including within two cells of a coarse-fine interface |
| `Src/amrCurvVerify.cpp` | curvature against the exact sphere |
| `Src/makeTestPlt.cpp` | synthetic hierarchies with analytic references, including a noisy field |
| `Src/GNUmakefile.AMRComplex` | standalone, AMReX only, `DIM ?= 3` |
| `Tools/amrBuildAndTest.sh` | builds every target, runs every suite |
| `Tools/runLensPipeline.sh` | the lens analysis end to end |
| `Tools/lensAnalysis.py` | every number in `Docs/ThermodiffusiveLens.tex` |

## Build and test

    cd ~/src/PeleAnalysis
    MPIRUN="mpirun --oversubscribe" Tools/amrBuildAndTest.sh

Drop `--oversubscribe` off OpenMPI. `DIM=2`, `RANKS=`, `COMP=`, `SKIP_BUILD=1`
all honoured. On the Mac add `COMP=llvm`. Six suites: complex and contour,
shapes, MEF, joint statistics, gradient statistics, and curvature. The
curvature suite runs in 3D only, deliberately: a plane curve has one principal
curvature and no Gaussian one, so four of its five sections would report a
vacuous pass in 2D.

## Running the lens analysis

    Tools/runLensPipeline.sh --infile <plotfile> --case <name> \
        --tu 693 --yu 1.325792e-2 --finest 3 --nproc 32

`--tu` and `--yu` are case properties and must be the prescribed inflow
values, not defaults: the Lewis estimates are very sensitive to them. Three
stages, selectable with `--stage grad|stats|analyse|all`. Stage one is the
expensive one and is what a GPU would help; stages two and three are
traversals and reductions.

Stage one writes `stage1.meta` beside its output with the temperature range
and the peak `|grad T|`, and later stages read it. Do not go back to scraping
the log for those: `--stage stats` used to erase the log it was about to read
and fall back to a compiled-in default, which silently binned two runs of the
same field differently.

Resource guide from `plt29791`, roughly 2.8 GB of output per component:

| finestLevel | cells | stage-1 output | peak RSS |
|---|---|---|---|
| 1 | 25 M | ~4 GB | ~6 GB |
| 2 | 121 M | ~17 GB | ~30 GB |
| 3 | 373 M | ~53 GB | ~95 GB |

**`finestLevel` is not a resolution knob.** It chooses how much of one
plotfile's hierarchy the analysis reads, and the levels above it have already
been averaged down into the levels below, so a level-1 analysis is a coarsened
view of the same solution. The sequence measures what a diagnostic loses on a
coarsened view and says nothing about whether the simulation resolves
anything. It is not a grid-convergence study and must not be reported as one.

What that costs, measured on `plt29791` at level 1 against level 3: the lens
does not care (same thresholds to four figures), the Lewis estimate does (21%
in the coldest band), and curvature admits more extremum at every level.

## What the analysis has established

On `plt29791` at full depth. These are settled; do not re-derive them.

**The lens.** `{T > t2} AND {Y > t1}` — hot AND rich, the off-band quadrant.
Getting this backwards gives the burnt corner and an enrichment below one. At
a 0.870% volume budget it captures 38.6% of fuel consumption against the
41.3% of the optimal region of that volume, an efficiency of 0.935 and an
enrichment of 44.6 times, at `t2 = 1445.6 K` and `t1 = 1.3465e-3`. A single
fixed `t2 = 1410 K` retains 0.916 of the two-parameter optimum across budgets.
The optimality condition is that the mean consumption density on the two faces
of the quadrant be equal.

**Lewis number.** The flux-ratio estimator gives `Le = 0.374 +/- 0.014`, flat
across 713 to 1400 K. Flatness is the result: a transport property should not
vary through the preheat zone, and on a coarsened view it appears to, by 28%
of its mean. The band-curvature estimator is biased about 45% high because
`Y_L(T)` is an average over flame elements of differing enrichment; prefer the
flux ratio.

**The asymmetry, which is the robust statement.** Split by the sign of the
excess-enthalpy scalar `S = theta + y - 1`: apparent Le is 0.409 on the rich
side, the lens, and 0.307 on the lean side, the cusps. The ratio 1.33 is
unchanged from level 2 on, while both absolute values are still moving. Same
mixture, same T, same pressure on both sides, so the difference is
differential diffusion and nothing else.

**The radical pool says the same thing, independently.** With
`Phi = (Y(HO2)+Y(H2O2))/(Y(H)+Y(O)+Y(OH))`, the lens side is about half as
termination-weighted as the cusps at the same temperature — ratio 0.52, flat
from 713 K through the preheat — rising monotonically through 1 at 1466 K.
Net `I_R(H)` is strongly negative through the preheat, peaking at 1122 K,
where the lens destroys 866 against the cusps' 302, a factor of 2.9. H is the
most diffusive species present; it streams upstream out of the reaction zone
and is destroyed in the preheat, and the lens receives nearly three times as
much.

**Not the crossover.** `Phi = 1` at 1283 K on the lens side and 1330 K on the
cusp side, well above the formal second-limit crossover of about 1149 K for
the unburnt composition at 8 atm, which is itself well above the 1088 K where
the curvature sensitivity of consumption turns. The turning point is not the
pool changing character; that hypothesis was tested and failed. A better one,
untested: the turn coincides with the peak of net H destruction at 1122 K,
which is one bin away. The pool lagging the rate crossover is expected — one
is a rate equality, the other an accumulation.

**Curvature predicts little that `(T, Y)` does not.** Explained variance of
`I_R(H2)`: temperature alone 0.283, `Y(H2)` alone 0.609, the pair 0.927, and
adding curvature 0.935 while multiplying the bin count elevenfold. The reading
is not that curvature is irrelevant but that it is what causes the `(T, Y)`
departure, so conditioning on the pair has already conditioned on it. The
direct test — explained variance of curvature given `(T, Y)` — has not been
run and is cheap.

## Conventions that are load-bearing

Each was learned from a defect. Breaking one reintroduces it.

1. **Identity is `(level, IntVect)`, never a coordinate.** Periodic
   identification then falls out of `FillBoundary` with no arithmetic and
   coarse-fine stitching needs no case table. This includes base cells:
   element ownership is decided on the base cell's image inside the domain.
2. **Anything entering a checksum is named by a gid or a tuple of gids, never
   an array index.** Two separate MPI defects were exactly this.
3. **Any ordering written to a file is canonicalised explicitly**, never
   inherited from traversal order. Nodes and elements both.
4. **Guarded regions are untested regions.** Serial builds never compile
   `#ifdef AMREX_USE_MPI`. Keep the guard around the MPI calls only — that is
   what `AMRParallel` is for.
5. **Declare then resize.** `Vector<T> v(std::size_t(n))` declares a function.
6. **Scientific outputs are integer-derived where possible.** Component
   volumes are counts per level times that level's cell volume, so they are
   bit-identical at any rank count rather than equal to a tolerance.
7. **Gradients and divergences are taken on the block-structured grid with the
   solver's own operators** — a composite `MLPoisson` applied, then its face
   fluxes averaged to centres — never with a difference stencil on a level's
   own ghost data, which mixes two resolutions in a ring around every refined
   patch. Never on a derived or dual structure. `apply()` must precede
   `getFluxes`, and `MLPoisson`'s flux is minus the gradient.
8. **A filter amplitude is uniform over the hierarchy.** Per-level amplitude
   puts a step in the field at every coarse-fine interface and the curvature
   of a step is a ring artefact. See the measurement in `AMRGradient.cpp`.
9. **We are in the GPU era.** Everything pointwise is a `ParallelFor` and
   every diagnostic reduction is a `ReduceOps`. There are no host loops over
   field data left in `amrGradStats`; do not add one.
10. **Two things named for one field must say which field.** Anything derived
    from a filtered field carries the pass count in its name, so a filtered
    column cannot be read as an unfiltered one in a later joint fit.

## Negative controls — do not delete these

Every one exists because the obvious check was caught passing on a broken
mesh. They are the most valuable thing in the tree.

| control | breaks | what still passes |
|---|---|---|
| `diagonal_rule=0` | face-diagonal consistency | area, and better than the correct rule |
| `orient_elements=0` | element winding | area and closure, both identical |
| `check_vertex_determinism=1` | element-dependent vertex data | everything, in serial |
| improper nesting | stitch with no valid coarse parent | nothing — it aborts |
| `single_level_stamp=1` | the coarse-parent walk past one level | nothing — it aborts on ghosts it cannot resolve |
| `legacy_ownership=1` | ownership without the periodic wrap | the face census; the volume gate must fail it |
| `per_level=1` in `amrGradStats` | the composite coarse-fine treatment | nothing — the gradient is wrong by 1e14 |
| `makeTestPlt noisy` | the surface, at cell scale | the field itself, and every first derivative of it |

**Area is not a topology test. Closure is not an orientation test. A smooth
analytic field cannot test a filter**, because there is nothing for it to
remove, so the test passes however little the filter does.

Synthetic geometry that caught real defects: `refine_shape`/`refine_anchor`
driving 13 shape cases, and `region=L|slab|sheet`. A centred box exercises
none of them. Section 6 of the shape suite is the wall-bounded closure case,
which the suite lacked until a merge dropped the wall-aware face parity and
every production reduction aborted while the complex was watertight.

## Calibrations you should not re-derive

The realizability meter in `amrGradStats` — the share of volume where
`kappa^2 < 4K`, which cannot happen on a real surface — against
`makeTestPlt`'s exact sphere and its noisy variant:

| mean severity | surface | error in `div(n)` |
|---|---|---|
| 3e-4 | exact, resolved | 0.4% |
| 2e-3 | perturbed 0.2 cells | 1000% |
| 4e-3 | perturbed 1 cell | 4700% |
| 8e-3 | perturbed 4 cells | 6000% |

So the meter has little dynamic range and errs towards saying nothing is
wrong. An order of magnitude above the clean baseline already means the
curvature magnitudes are worthless. A violation proves noise; its absence
proves little.

Filtering buys an order of magnitude and no more: at a fifth of a cell of
surface noise, eight passes take `div(n)` from ten times the true value to
about twice it. **Conditional curvature statistics on a filtered field are
usable. A curvature extremum is not, and no pass count makes it so** —
consistent with `max|kappa|` tracking `1/dx` at every level of the real case.

## Next

1. **Species and curvature in one reduction.** The only way to test whether
   the 1088 K turn in the curvature sensitivity is the 1122 K peak of net H
   destruction. Needs stage one re-run with the radicals carried, about 75 GB
   at level 3, so run level 2 first at about 26 GB.
2. **Explained variance of curvature given `(T, Y)`.** One reduction. Settles
   whether the lens is a curvature diagnostic that needs no second
   derivatives, which would dodge the whole noise problem.
3. **`Y_L^0(T)` from a 1-D laminar flame.** The weakest link in the chain: it
   defines `S`, so it decides which side of the band every bin falls on, and
   the asymmetries above are measured against a centreline taken from the 3-D
   field itself.
4. **The lens boundary as a complex.** Two isosurface pieces — part of
   `{T = t2}` inside `{Y > t1}`, part of `{Y = t1}` inside `{T > t2}` — joined
   on the seam curve where both thresholds are met. The area fraction carried
   by the Y piece is zero exactly when Le is one, so it is an order parameter
   for thermodiffusive departure, per connected component. Validation comes
   free: the Gaussian curvature integrated over a closed component is
   `2*pi*chi`, and `chi` is already available from element counts, so two
   independent routes to one integer.
5. Curvature statistics on the lens surface: area-weighted shape index and
   curvedness, conditioned on which piece and on distance from the seam.
6. Merge tree, gated on reproducing `amrSweep`'s `N(v)` exactly at every
   sampled threshold; then overlap matching, for both time tracking and
   cross-surrogate comparison.
7. 3D figures for `Docs/ThermodiffusiveLens.tex` from `amrLabel` output.

Open `[DECIDE]`: whether a cell is thresholded from above or below may differ
per surrogate, so it stays a runtime choice.

## Known limits and loose ends

- `Docs/reeb-additions-from-newmachine.md` is 191 lines from the other machine
  that have not been ported into `Docs/design/`. Until they are, it is the
  only place some of that reasoning is written down.
- `makeTestPlt` carries two option vocabularies for one idea —
  `refine_shape`/`refine_anchor` and `region` — because both are load-bearing
  for different suites. The TODO is in the source. Collapsing them means
  rewriting case lists in two suites.
- `legacy_ownership` and `single_level_stamp` are hand-verified and wired into
  section 4 and 5 of the shape suite, but not into its parallel case arrays.
- The full `plt01800` sweep peaks at 83 GB of 122; about 36 GB of that is
  added during the sweep and not attributed.
- Reading one field from a 180 GB plotfile dominates wall time: 446 s of the
  full-hierarchy build.
- `writeMEF` gathers to the IO processor, the same wall `isosurface.cpp` has.
  `.pvtu` output would retire it.
- The walled volume reference cannot see an element whose every corner lies on
  one wall; the face census still does.
- Bare MPI executables hang outside `mpirun` on the Mac; launch with
  `mpirun -n 1` even for serial runs.
- `curvature.cpp` builds one single-level operator per level and so carries the
  coarse-fine artefact `AMRGradient` exists to avoid. It is the older tool;
  prefer `amrGradStats`. Its `do_smooth` is a diffusion solve to a time, which
  is the form `smooth_passes` replaced.
