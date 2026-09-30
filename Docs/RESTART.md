# Restarting this work in a new session

Origin session: `https://claude.ai/code/session_01D4f1Detq1YmcvLsYuCywQA`
(also in the `Claude-Session:` trailer of every commit on `amr-complex`).

Nothing of value lives in that session. Everything is committed. A new session
needs this file, the design note under `Docs/design/`, and the repo.

The design note is LaTeX (`cd Docs/design && make`). It was rewritten natively
from an earlier Markdown draft, which is removed — that draft grew incrementally
with corrections patched in place, so reading it top to bottom absorbed the
mistakes along with the fixes. The LaTeX document is organised as a document and
states the settled position once. Earlier drafts remain in git history.

---

## Paste this to start

> I'm continuing work on AMRComplex/AMRContour in ~/src/PeleAnalysis, branch
> `amr-complex`. Read `Docs/RESTART.md`, then the design note sources under
> `Docs/design/` — start with `sec-problem.tex` and `sec-roadmap.tex`, which give
> the target and the current state, then whichever of `sec-complex`,
> `sec-contour`, `sec-segment`, `sec-verification`, `sec-topology` or `sec-gpu`
> the task touches. Then tell me what you think the next step is before writing
> any code. Note: technical documents in this repo are written natively in
> LaTeX, not generated from Markdown.

---

## What this is

Post-processing for AMReX plotfiles, aimed at characterising the cellular
structures of the thermodiffusive instability in lean premixed hydrogen
flames. The open scientific question — what a flame *cell* actually is — is
why the tooling is built around comparing candidate definitions rather than
computing one. `Docs/design/sec-problem.tex` has the full statement.

## What is built and verified

On `amr-complex`, all suites passing in 2D and 3D at 1–8 ranks:

| file | what |
|---|---|
| `Src/AMRComplex.{H,cpp}` | conforming multilevel cell complex; node identity is `(level, IntVect)` |
| `Src/AMRContour.{H,cpp}` | marching simplices, oriented, carried components, MEF output |
| `Src/AMRSegment.{H,cpp}` | connected components of a threshold set |
| `Src/AMRParallel.{H,cpp}` | the only place MPI calls appear |
| `Src/amrComplexVerify.cpp`, `amrContourVerify.cpp`, `amrSweep.cpp` | drivers |
| `Src/makeTestPlt.cpp` | synthetic hierarchies with analytic references |
| `Src/GNUmakefile.AMRComplex` | standalone, AMReX only, `DIM ?= 3` |
| `Src/AMRSegment.{H,cpp}` | connected components of a threshold set |
| `Src/AMRParallel.{H,cpp}` | the only place MPI calls appear |
| `Tools/amrBuildAndTest.sh` | builds everything, runs every suite |
| `Docs/design/` | the design note, LaTeX |

## Build and test

    cd ~/src/PeleAnalysis
    MPIRUN="mpirun --oversubscribe" Tools/amrBuildAndTest.sh

Drop `--oversubscribe` off OpenMPI. `DIM=2`, `RANKS=`, `COMP=`, `SKIP_BUILD=1`
all honoured. `AMREX_HOME` defaults to `../../amrex`, correct for this layout.

## Conventions that are load-bearing

These were each learned from a defect. Breaking one reintroduces it.

1. **Identity is `(level, IntVect)`, never a coordinate.** Periodic
   identification then falls out of `FillBoundary` with no arithmetic, and
   coarse-fine stitching needs no case table.
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

**Area is not a topology test. Closure is not an orientation test.**

## Next

1. **Representative snapshots.** One plotfile with two or three candidate
   surrogates — heat release rate, temperature, an intermediate — swept
   individually with `amrSweep`, and the `N(v)` curves laid side by side.
   Answers whether a plateau exists, gives real timings, and shows whether the
   surrogates agree, before anything is built on an assumption about any of
   the three.
2. Segmentation written as a plotfile (integer label per cell), so any field
   can be conditioned on cells defined by any other field.
3. **Merge tree** — gated on reproducing `amrSweep`'s `N(v)` exactly at every
   sampled threshold.
4. Overlap matching, serving both time tracking and cross-surrogate comparison.
5. Boundary invariants for the b1 events (a crescent closing into a ring).

Open `[DECIDE]`: whether a cell is thresholded from above or below — may
differ per surrogate, so it stays a runtime choice.

Deferred deliberately: wall extrapolation; `.pvtu` output to retire the gather
in `writeMEF`; `curvature.cpp` on tet-derived normals in the stitch band, which
needs real geometry. See `Docs/design/sec-roadmap.tex`.

Open decision recorded in `Docs/design/sec-gpu.tex`: the implementation is
host-only, and whether the merge tree is written for host or device should be
settled before it is written, not after.

## Known limits

- `writeMEF` gathers to the IO processor. Same wall `isosurface.cpp` has.
- `ownedBaseRegions` and the segmentation maps have not been run at
  production size.
- Everything so far is synthetic. Real regridded hierarchies will exercise the
  proper-nesting assertion in a way synthetic ones cannot; if it aborts, that
  is the check working.
