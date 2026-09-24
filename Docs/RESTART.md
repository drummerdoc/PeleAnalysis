# Restarting this work in a new session

Origin session: `https://claude.ai/code/session_01D4f1Detq1YmcvLsYuCywQA`
(also in the `Claude-Session:` trailer of every commit on `amr-complex`).

Nothing of value lives in that session. Everything is committed. A new session
needs this file, `Docs/reeb-graph-assessment.md`, and the repo.

---

## Paste this to start

> I'm continuing work on AMRComplex/AMRContour in ~/src/PeleAnalysis, branch
> `amr-complex`. Read `Docs/RESTART.md` first, then `Docs/reeb-graph-assessment.md`
> sections 14 to 18 — the verified state, the target problem, the sweep, what
> the first real snapshot exposed, and the memory work that made it fit. Do not re-derive the earlier sections;
> they are history and several of their conclusions were later corrected in
> place. Then tell me what you think the next step is before writing any code.

---

## What this is

Post-processing for AMReX plotfiles, aimed at characterising the cellular
structures of the thermodiffusive instability in lean premixed hydrogen
flames. The open scientific question — what a flame *cell* actually is — is
why the tooling is built around comparing candidate definitions rather than
computing one. Section 15 of the assessment has the full statement.

## What is built and verified

On `amr-complex`, all suites passing in 2D and 3D at 1–8 ranks, and the
complex verified exactly on all four levels of a real PeleLMeX flame sheet
(`~/DATA/plt01800`, 631M elements), which sweeps in 83 GB at 64 ranks
(assessment §17, §18):

| file | what |
|---|---|
| `Src/AMRComplex.{H,cpp}` | conforming multilevel cell complex; node identity is `(level, IntVect)` |
| `Src/AMRContour.{H,cpp}` | marching simplices, oriented, carried components, MEF output |
| `Src/AMRSegment.{H,cpp}` | connected components of a threshold set |
| `Src/AMRParallel.{H,cpp}` | the only place MPI calls appear |
| `Src/amrComplexVerify.cpp`, `amrContourVerify.cpp`, `amrSweep.cpp` | drivers |
| `Src/makeTestPlt.cpp` | synthetic hierarchies with analytic references |
| `Src/GNUmakefile.AMRComplex` | standalone, AMReX only, `DIM ?= 3` |
| `Tools/amrBuildAndTest.sh` | builds everything, runs every suite |

## Build and test

    cd ~/src/PeleAnalysis
    MPIRUN="mpirun --oversubscribe" Tools/amrBuildAndTest.sh

Drop `--oversubscribe` off OpenMPI. `DIM=2`, `RANKS=`, `COMP=`, `SKIP_BUILD=1`
all honoured. `AMREX_HOME` defaults to `../../amrex`, correct for this layout.

## Conventions that are load-bearing

These were each learned from a defect. Breaking one reintroduces it.

1. **Identity is `(level, IntVect)`, never a coordinate.** Periodic
   identification then falls out of `FillBoundary` with no arithmetic, and
   coarse-fine stitching needs no case table. **This includes base cells:**
   element ownership is decided on the base cell's image inside the domain.
   A seam base had two raw names and was owned once under each (§17.1).
2. **Anything entering a checksum is named by a gid or a tuple of gids, never
   by an array index.** Two separate MPI defects were exactly this.
3. **Any ordering written to a file is canonicalised explicitly**, never
   inherited from traversal order. Nodes *and* elements.
4. **Guarded regions are untested regions.** Serial builds never compile
   `#ifdef AMREX_USE_MPI`, and no MPI is available to the assistant. Keep the
   guard around the MPI calls only — that is what `AMRParallel` is for.
5. **Declare then resize.** `Vector<T> v(std::size_t(n))` declares a function.
   Three most-vexing-parses so far.
6. **Scientific outputs are integer-derived where possible.** Component
   volumes are counts per level times that level's cell volume, so they are
   bit-identical at any rank count rather than comparable to a tolerance.

## Negative controls — do not delete these

Every one exists because the obvious check was caught passing on a broken
mesh. They are the most valuable thing in the tree.

| control | breaks | what still passes |
|---|---|---|
| `diagonal_rule=0` | face-diagonal consistency | area, and *better* than the correct rule |
| `orient_elements=0` | element winding | area *and* closure, both identical |
| `check_vertex_determinism=1` | element-dependent vertex data | everything, in serial |
| improper nesting | stitch with no valid coarse parent | nothing — it aborts |
| `mpi_SH` forced periodic in its wall direction | nesting across a periodic face | nothing — it aborts (it used to drop a layer silently) |
| seam-ownership mutant (scratch only, §17.3) | one stitch element per seam, owned twice | the walled volume check must fail it; it does |

**Area is not a topology test. Closure is not an orientation test.** And a
reference measured from the mesh must be shown able to fail: the wall
reference was mutation-tested for exactly that reason.

Synthetic geometry that caught real defects: `makeTestPlt region=L|slab|sheet`
(re-entrant corners, refined levels meeting periodic seams, a refined level
on a wall). A centred box exercises none of them.

## Next

1. **Representative snapshots.** The four-marker sweep on the full
   hierarchy of `plt01800`: `I_R(H2)`, `temp`, `HeatRelease`, `Y(H2O2)`,
   `is_per="1 1 0"`, `N(v)` curves side by side. `I_R(H2)` is negative where
   fuel is consumed: sweep it from below (`above=0`), `vhi=0`, `vlo` its
   minimum.
2. Segmentation written as a plotfile (integer label per cell).
3. **Merge tree** — gated on reproducing `amrSweep`'s `N(v)` exactly at every
   sampled threshold.
4. Overlap matching, serving both time tracking and cross-surrogate comparison.
5. Boundary invariants for the b1 events (a crescent closing into a ring).

Open `[DECIDE]`: whether a cell is thresholded from above or below — may
differ per surrogate, so it stays a runtime choice.

Deferred deliberately: wall extrapolation (assessment §7.2); `.pvtu` output to
retire the gather in `writeMEF`; `curvature.cpp` on tet-derived normals in the
stitch band, which needs real geometry.

## Known limits

- The full `plt01800` sweep peaks at 83 GB of 122; about 36 GB of that is
  added during the sweep and not yet attributed (§18.2).
- Reading one field from the 180 GB plotfile dominates wall time: 446 s of
  the full-hierarchy build.
- `writeMEF` gathers to the IO processor. Same wall `isosurface.cpp` has.
- The walled volume reference cannot see an element whose every corner lies
  on one wall; the face census still does (§17.3).
- Bare MPI executables hang outside `mpirun` on the development machine;
  always launch with `mpirun -n 1` even for serial runs.
