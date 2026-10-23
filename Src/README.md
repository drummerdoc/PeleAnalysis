# Tools in this directory

Each entry is one driver, what it does, and a line that actually runs. The
options are real — `ParmParse` ignores a name it does not recognise without
saying so, so an invented option produces a run that looks fine and answers a
different question. Prefix with `mpiexec -n <N>` for the MPI builds.

Build any of them standalone, AMReX only:

    make -f GNUmakefile.AMRComplex USE_MPI=TRUE DIM=3 EBASE=<name> -j8

Add `COMP=llvm` on the Mac. `AMREX_HOME` defaults to `../../amrex`.
`Tools/amrBuildAndTest.sh` builds every target below and runs every suite;
start there rather than building one tool by hand.

Fuller treatments: `Docs/design/` for the complex, contour and segmentation
machinery, `Docs/ThermodiffusiveLens.tex` sections 8 and 9 for the gradient,
curvature and statistics tools and how the lens analysis uses them, and
`Docs/RESTART.md` for the state of the work and the conventions that are
load-bearing.

---

## Added in this session

### `amrCurvVerify.cpp`

Checks the curvature fields `amrGradStats` writes against an exact answer.
`makeTestPlt`'s `dist` field is the signed distance to a sphere, so every
level set of it is a sphere and at radius `r` the answers are `div(n)` =
`(DIM-1)/r`, `K` = `1/r^2` and `k1` = `k2` = `1/r` — at every `r`, so one
plotfile tests the whole chain over a continuum of curvatures with nothing to
tune. Also asserts `k1 >= k2`, which a sphere cannot check by itself because
it is umbilic.

    amrCurvVerify3d.gnu.MPI.ex infile=plt_sph_g field=dist rlo=0.15 rhi=0.35

`tol` is the relative tolerance on the worst cell, default 0.05; `centre`
overrides the sphere centre, which defaults to the middle of the domain.

---

## Modified in this session

### `amrGradStats.cpp`

Writes a plotfile carrying the input variables plus, for each, the
cell-centred gradient and its magnitude; for each pair, the alignment of
their gradients; and for one named field the unit normal and curvatures. All
gradients come from a composite `MLPoisson` applied over the whole hierarchy,
so the coarse-fine interface is treated the way the solver treats it rather
than by a difference stencil on interpolated ghost data. The output is a
plotfile on the same hierarchy, so `amrJointFit` consumes it unchanged.

    amrGradStats3d.gnu.MPI.ex infile=plt29791 outfile=plt_g finestLevel=3 \
        is_per="1 1 0" vars="temp Y(H2)" curvOf=temp \
        carry="HeatRelease I_R(H2)"

New this session: `gaussCurv=1` adds the Gaussian curvature, the principal
pair, the shape index and the curvedness, in 3D, for one extra field of
memory and no extra solves. `smooth_passes=N` low-passes the `curvOf` field
N times before anything is differentiated — a count, not a smoothing time —
and renames everything derived from it, so a filtered column cannot be read
as an unfiltered one later.

    amrGradStats3d.gnu.MPI.ex infile=plt29791 outfile=plt_g_s2 finestLevel=3 \
        is_per="1 1 0" vars="temp Y(H2)" curvOf=temp gaussCurv=1 \
        smooth_passes=2

It reports the normal floor and what fraction of volume sits below it, and
with `gaussCurv` the share of volume where `kappa^2 < 4K`, which cannot
happen on a real surface. Read that number against the calibration table in
`Docs/RESTART.md` before trusting any curvature magnitude.

### `AMRGradient.{H,cpp}`

The library behind the above: `compute` for the composite gradient, and new
this session `smooth` for the counted low-pass. One amplitude for the whole
hierarchy, taken from the finest spacing — a per-level amplitude puts a step
in the field at every coarse-fine interface, and the curvature of a step is a
ring artefact around every refined patch. The measurement that settles it is
in the `.cpp`. Not a driver; no invocation of its own.

### `makeTestPlt.cpp`

Synthetic multi-level hierarchies with analytic references: properly nested
levels of a requested shape, carrying a distance function, one linear field
per direction, a cellular field with a known number of maxima, and new this
session `noisy` — the same sphere displaced cell to cell by a stated fraction
of a cell, from a hash of the index so it is identical on any machine at any
rank count.

    makeTestPlt3d.gnu.MPI.ex n_cell=64 nlevels=3 refine_shape="random random" \
        refine_anchor="lo center" refine_blk=8 refine_frac=0.5 \
        noise_amp=0.0015625 outfile=plt_test

`noisy` exists because a smooth analytic field cannot test a filter: there is
nothing for a low-pass to remove, so the test passes however little the
filter does. `region=L|slab|sheet` is a second, older vocabulary for the
refinement pattern, kept because two suites depend on it.

### `AMRComplex.{H,cpp}`

The conforming multilevel cell complex everything else is built on; node
identity is `(level, IntVect)`, never a coordinate. Restored this session:
`faceParityInterior`, the face invariant taken over interior faces only. A
face on a domain wall belongs to one element and cannot cancel, so the
all-faces parity is nonzero by construction wherever there is a wall, and a
merge that dropped the interior one left every production reduction aborting
with "complex is not watertight" on a complex that was watertight. Not a
driver; exercised through `amrComplexVerify`.

---

## Used for the analysis, unchanged

### `amrComplexVerify.cpp`

Builds the complex and checks it: volume against the domain by the divergence
theorem, face closure, unfilled nodes, and order-independent checksums over
the node and element sets. The checksums are the actual MPI invariant;
matching element counts is weak evidence.

    amrComplexVerify3d.gnu.MPI.ex infile=plt29791 is_per="1 1 0" \
        varnames=temp face_census=1

`face_census=1` adds a rank-local incidence census and the duplicate-element
count, which costs a map holding every face instance — most of the verify's
memory at production size. `legacy_ownership=1` and `single_level_stamp=1`
are negative controls and must fail.

### `amrJointFit.cpp`

Volume-weighted joint statistics in range space: for a set of coordinates it
accumulates the conditional mean and standard deviation of each response per
bin, and reports how much of each response's variance the best function of
every coordinate subset explains. There is no fit to argue about — the best
function is the conditional mean.

    amrJointFit3d.gnu.MPI.ex infile=plt_g is_per="1 1 0" \
        coords="temp Y(H2)" clo="693 0" chi="2210 1.325792e-2" \
        cbins="128 192" resp="HeatRelease I_R(H2)" outfile=ty_fine.dat

It reads a raw plotfile as happily as a derived one, so a reduction that
needs no gradients does not need stage one at all:

    amrJointFit3d.gnu.MPI.ex infile=plt29791 is_per="1 1 0" \
        coords="temp Y(H2)" clo="693 0" chi="2212 1.325792e-2" \
        cbins="128 192" \
        resp="Y(H) Y(O) Y(OH) Y(HO2) Y(H2O2) I_R(H)" \
        outfile=kinetics_TY.dat

Pick `clo`/`chi` from the component ranges `amrGradStats` prints: too narrow
and most of the volume falls out of range, too wide and the bins are spent on
empty space.

### `amrGradVerify.cpp`

Measures gradient error against fields whose gradient is known, separately
over the whole domain, the interior, and within two cells of a coarse-fine
interface. The three regions are the point: the interface error is a surface,
so it does not converge away, and an all-domain norm hides it.

    amrGradVerify3d.gnu.MPI.ex infile=plt_g is_per="0 0 0" \
        check="gradx_xco grady_xco gradMag_dist" expect="1 0 1"

Parallel `check`/`expect` arrays rather than `name=value` pairs because
`ParmParse` chokes on `=` inside a quoted value.

### `amrSweep.cpp`

Connected-component count against threshold for one field — the `N(v)` curve
— with optional areas and a minimum volume to suppress single-cell
components.

    amrSweep3d.gnu.MPI.ex infile=plt29791 isoCompName=HeatRelease \
        vlo=1e8 vhi=1e10 nval=40 log=1 is_per="1 1 0"

`vlo`/`vhi` default to the field's range over the hierarchy. For `I_R(H2)`,
which is negative where fuel is consumed, sweep from below with `above=0`.

### `amrLabel.cpp`

Segmentation written out: an integer component label per cell at one or more
thresholds, so any field can be conditioned on cells defined by another.
`surfaces=1` also writes the component surfaces.

    amrLabel3d.gnu.MPI.ex infile=plt29791 isoCompName=HeatRelease \
        is_per="1 1 0" thresholds="4.35e10 4.27e11 6.88e11" \
        above="0 0 1" outfile=plt_seg

The option is `surfaces`, not `surf`; `surf_dir`, `surf_suffix`, `surf_max`
and `surf_skip` control where the surface files go and how many are written.

### `amrContourVerify.cpp`

Isosurface area and closure against analytic references, with the contour
negative controls: `diagonal_rule=0`, `orient_elements=0` and
`check_vertex_determinism=1` each break something real while leaving the
obvious check passing.

    amrContourVerify3d.gnu.MPI.ex infile=plt_test is_per="1 1 1" isoVal=0

### `surfMEFtoDAT.cpp`

Converts MEF surface output to Tecplot-readable ASCII.

    surfMEFtoDAT3d.gnu.MPI.ex infile=surf.mef outfile=surf.dat

---

## Scripts that drive these

In `Tools/`, documented in their own headers:

| script | what |
|---|---|
| `amrBuildAndTest.sh` | builds every target above, runs all six suites |
| `runLensPipeline.sh` | the lens analysis end to end, three stages |
| `lensAnalysis.py` | six analysis modes over `amrJointFit` tables |
| `amrCurvatureTest.sh` | twelve checks on the curvature chain and the filter |
| `amrGradStatsTest.sh` | nineteen checks per dimension on the gradient |
| `amrComplexShapeTest.sh` | non-convex and wall-adjacent refinement |
| `amrComplexMpiTest.sh` | rank-count invariance |
| `amrContourMefTest.sh` | surface output |
| `amrJointStatsTest.sh` | the reduction |
