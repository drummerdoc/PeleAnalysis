# Periodic-aware Reeb graphs from AMReX plotfiles: assessment

Question: can `Src/isosurface.cpp` be the basis for computing periodic-aware
Reeb graphs of AMR data, and if not, what respects coarse-fine interfaces
properly?

Short answer: the **front end** of `isosurface.cpp` is a good basis and should be
reused. The **engine** is not, and cannot be extended into one — not because the
marching-cubes code is poor, but because a Reeb graph is a property of the scalar
field over the whole domain, and marching cubes produces a triangle soup at one
isovalue with no domain connectivity attached. The object every Reeb/contour-tree
algorithm consumes is a *cell complex with vertex scalars and explicit
adjacency*; nothing in `isosurface.cpp` builds one.

The recommendation is to build that complex — an AMR dual mesh with stitch cells
at coarse-fine interfaces — and derive both the topology and the isosurface from
it. That also removes, rather than patches, the periodic defect already
documented in `isosurface.cpp`'s CHECK 1.

---

## 1. What is reusable from `isosurface.cpp`

Roughly the first 250 lines of `main()` are directly reusable, and are the part
that is tedious to get right:

- `PlotFileData` open, component selection by name, `finestLevel` override.
- Per-level `Geometry` construction with `is_per` from `ParmParse`, and
  `growPeriodicDomain`.
- `FillPatchSingleLevel` / `FillPatchTwoLevels` to give each level valid ghost
  data from its coarse parent.
- The fine-covered mask: coarsen each `grids[lev+1]` box by `refRatio`, mark the
  intersection, *including the periodic-shift pass*. In the dual-mesh
  formulation this same criterion selects exactly which coarse cells contribute
  a dual vertex, so the code carries over almost verbatim.
- The MEF/XDMF writers, for dumping the derived isosurface once the complex
  exists.

Separately, and more importantly: **`Src/amrToFE.cpp` already builds the
conforming multilevel complex** that §4A calls for. It is the real starting
point; see there.

## 2. What does not carry over, and why

### 2.1 One isovalue is not a Reeb graph

Carr-Snoeyink-Axen (contour tree), Pascucci et al. (online Reeb), Tierny et al.
(loop surgery), Gueunet et al. (FTR) all sweep *sorted vertices* of a mesh with
union-find over *cell adjacency*. `Polygonise` throws that away: it emits
`Triangle`s of `PMapIt` into a `PMap` keyed by `Edge`, and the global `eltSet`
stores three node indices per element. There is no path from that to critical
points and their pairing.

There is a cheap approximation — extract isosurfaces at many levels and track
components between adjacent levels — but it misses every critical point between
sampled levels and never gives exact saddles. Worth ~50 lines as a sanity toy;
not worth building on.

### 2.2 The coarse-fine seam is approximately closed, not watertight

Current behaviour: `Polygonise` bails if any of the eight corners has
`mask < 0`, so coarse dual cells touching a covered cell emit nothing; the fine
level covers them via its ghost band. Coarse and fine nodes on the interface
plane then land at *different positions* and are welded only if they compare
bitwise-close in `Node::operator<`. The seam therefore carries cracks and
T-junctions.

For an area integral that is a percent-level error. For topology it is fatal: a
crack splits one component into two, inventing a Reeb arc and a spurious
minimum. Widening the weld tolerance does not fix it — see 2.3.

### 2.3 `Node::operator<` is not a strict weak ordering

```cpp
if (std::sqrt(sum) < epsilon_DEF) return false;   // then lexicographic
```

With `a ~ b` and `b ~ c` but `a !~ c`, transitivity fails, so `std::set<Node>`
behaviour becomes insertion-order dependent. At `epsilon_DEF = 1e-15` this is
essentially exact equality and mostly benign. **Any attempt to raise the
tolerance to weld coarse nodes to fine nodes turns this into non-deterministic
topology** — different rank counts, different graph. If node welding survives
into the new tool at all, it needs a snap-to-grid key (quantise to a fixed
lattice, hash the integer triple), not an epsilon comparison.

### 2.4 `PCInterp` in the ghost band

The coarse-to-fine fill is piecewise constant. The field a fine box sees in its
ghost region is a staircase of coarse values. Two consequences for topology:
(a) it is not the same function the coarse level sees, so the two sides of the
seam genuinely disagree; (b) it manufactures exactly-flat plateaus, which are
degenerate for any sweep algorithm and create spurious critical points at every
interface. Conservative-linear or node-bilinear interpolation is the minimum,
but the real fix is a single function definition on a single complex.

### 2.5 Everything collates to rank 0

`Collate()` gathers all nodes and elements to `IOProcessor`, which then rebuilds
`nodeSet`/`eltSet` serially. For "large snapshots" this is already the binding
constraint for isosurfaces; a merge tree over the full field is a much larger
object. A distributed design is required from the start, not retrofitted.

### 2.6 The known periodic defect is structural

CHECK 1's comment records it: across a periodic face the crossing is interpolated
from raw coordinates rather than the minimum image, so the node lands mid-domain
(observed inflating arc length by 45%, worst plotfile 109%). This happens because
connectivity is inferred *from coordinates*. In a dual mesh, periodicity is a
property of the connectivity — the cells on the high face simply reference the
low-face vertices — and coordinates are never used to decide who touches whom.
The bug class disappears rather than being tested for.

---

## 3. The structural point: merge tree vs contour tree vs Reeb graph

This decides how much work the project is, and should be settled first.

| | domain | shape | cost |
|---|---|---|---|
| Join / split (merge) tree | any | **always a tree, even on a 3-torus** | union-find sweep, distributed implementations exist |
| Contour tree | join ⋈ split | equals the Reeb graph **only if the domain is simply connected** | cheap, but *invalid on a periodic box* |
| Reeb graph | any | carries loops | needs loop surgery or a genuine Reeb algorithm |

A triply periodic box is not simply connected, so the contour tree route is not
available: this is exactly what "periodic-aware" costs. The standard practical
construction is **loop surgery** (Tierny, Gyulassy, Pascucci, TVCG 2009): locate
the loop-carrying edges, cut the mesh along them so it becomes simply connected,
compute the contour tree, then glue the cuts back. FTR (Gueunet et al.) is the
current task-parallel implementation in TTK.

### What a Reeb loop actually is

A cycle in the Reeb graph is a one-parameter family of contours that leaves a
level, travels, and returns to where it started — you can walk around it. Cycles
come from non-trivial 1-cycles in the **domain**, not from the shape of any one
contour. On a triply periodic box, `b1(T^3) = 3`, so the Reeb graph can carry up
to three independent loops, and each one means a structure that closes on itself
through the periodic wrap: a percolating sheet, a connected pocket network that
spans the box and rejoins itself.

That is the whole content of "periodic-aware" at the graph level. On a
non-periodic box `b1 = 0`, the Reeb graph is a tree, and it is the contour tree.

### What the Reeb graph does *not* tell you

It does **not** give the genus of an individual level-set component. A saddle
that adds a handle to a contour without merging or splitting it is a *degree-2*
node in the Reeb graph, and degree-2 nodes are conventionally dropped — they sit
in the middle of an arc and are invisible in the output. So handle counts,
tunnels, and genus evolution are not readable from a Reeb graph, and building one
to get them would be the wrong tool.

The right tool for that is much cheaper: extract the isosurface, split it into
connected components, and compute the Euler characteristic per component
directly, `chi = V - E + F`, then `g = (2 - chi)/2` for a closed orientable
component. Exact, O(elements), no sweep at all. Its one prerequisite is that the
surface be **watertight** — which is precisely what §2.2 says the current
coarse-fine seam is not, and precisely what the dual mesh of §4A delivers. So
the genus question is an argument for the dual mesh, and an argument *against*
the Reeb graph.

### [DECIDE] Which question is being asked?

| Question | Tool | Effort |
|---|---|---|
| How many disconnected pockets/kernels, when do they merge or split? | **periodic merge tree** — union-find with periodic adjacency, nothing more | smallest |
| Genus / handles / tunnels of a given isosurface component | **Euler characteristic per component** on a watertight surface | trivial, *given* §4A |
| Does a structure percolate, wrap the box, close on itself through the wrap? Families of contours that rejoin? | **Reeb graph with loops** — loop surgery | largest |

All three need the dual mesh of §4A. Only the third needs loop surgery. Note the
first two are periodic-*correct* without being periodic-*aware* in the graph
sense: periodicity enters as adjacency in the union-find, or as connectivity in
the surface, and never produces a cycle in the output.

---

## 4. Respecting coarse-fine interfaces: three routes

### Route A — AMR dual mesh, conforming across levels (recommended)

**This already exists in the tree: `Src/amrToFE.cpp`.** It is the construction
this section was going to specify, and it is the right basis — not
`isosurface.cpp`, and not `flattenAMRFile.cpp`.

#### What it does

- `Node` identity is **`(level, IntVect)`** — topological, never coordinates:

  ```cpp
  struct Node {
    enum typeEnum{INIT=0, COVERED=1, VALID=2};
    bool operator< (const Node& rhs) const {
      if (level < rhs.level) return true;
      if ((level == rhs.level) && (iv<rhs.iv)) return true;
      return false; }
    int level; IntVect iv; int grid; typeEnum type; };
  ```

  `IntVect::operator<` is lexicographic, so this is a valid strict weak ordering
  — unlike `isosurface.cpp`'s epsilon comparison (§2.3). Two elements share a
  vertex because they name the same `(level,iv)`, not because two floats came out
  close. This is the whole reason it is watertight.

- With `connect_cc=1` (the default) vertices sit at **cell centres** and elements
  join neighbouring centres — the dual mesh. (`connect_cc=0` puts nodes at cell
  corners and copies the cell value out; that is a piecewise-constant voxel mesh
  and is useless for topology. Keep the default.)

- **The coarse-fine transition.** For `lev>0` it takes the ring of cells within
  `ref` of the fine grids but outside them (`GetBndryCells`) and **stamps the
  coarse parent's Node into the fine ghost positions**: for each coarse cell
  `iv`, all `ref^d` fine positions `ref*iv + ivt` receive
  `Node(iv, lev-1, VALID)`. Coarse cells under a finer grid are marked `COVERED`.
  An element is emitted only when all `2^d` corners are `VALID`.

  So the transition elements are emitted at the *fine* level, and their
  coarse-side corners **are the coarse level's own nodes, by identity**. Because
  `ref^d` fine positions collapse onto one coarse Node, those elements carry
  repeated node indices — they are **degenerate hexahedra**, i.e. pyramids,
  prisms and tets written in brick form. That is a legitimate alternative to
  Weber et al.'s explicit stitch-cell case tables, and it obtains conformity
  without any case table at all. It also means `ref_ratio` is not restricted to
  2, which retires that `[DECIDE]`.

#### What is missing for this purpose

1. **Periodicity is duplication, not identification.** `nGrowPer` grows
   `subboxArray` past the domain face and appends periodically-shifted copies of
   the boxes. The stamped IntVects there are *outside* the domain, so they are
   **different Nodes** from their partners on the far side. The output is an
   extended slab with duplicated material — right for a seamless picture, wrong
   for topology: the complex is not a torus, so it has no 1-cycles and its Reeb
   graph cannot have loops.

   The fix needs no arithmetic at all — see §7.3. Because Node identity is
   *stored as data*, a periodic `FillBoundary` on the `MultiNodeFab` copies the
   far-side Node verbatim into the near-side ghost, IntVect included and
   unshifted. The ghost cell at `hi+1` then literally contains
   `Node(lo, lev, VALID)`, which *is* the identification. Elements close through
   the face and the complex becomes a genuine 3-torus.

   Consequence to design around: on a torus the coordinates cannot be
   single-valued — a wrap element has corners whose positions differ by `L`.
   That is fine for the graph (coordinates are decoration there) and wrong for
   viz. **[DECIDE]** carry both: a topological output with identified nodes for
   the sweep, and the existing duplicated-slab output for pictures.

2. **Output format.** It writes Tecplot (`tec`) or `flt`. For the topology it
   needs to emit tets (or a VTU that TTK reads). Degenerate hexes must be
   subdivided — Freudenthal, with the diagonal chosen by a parity rule on the
   base index so neighbouring cells agree across shared faces, and degenerate
   tets dropped. The parity rule needs care where levels meet, since the two
   sides index differently.

3. **Verification that it is actually gap-free.** The construction should tile
   the interface exactly, but that should be measured, not assumed. **Sum the
   element volumes and compare with the domain volume** — the volumetric
   analogue of the edge check below. Any gap or overlap at the transition shows
   up immediately, per level pair.

4. It predates the modern AMReX API (`DataServices`/`AmrData`, `BL_SPACEDIM`)
   and is single-`std::map`-on-one-rank in the same way `isosurface.cpp` is
   (§2.5). Fine for development, a wall at production size.

### Route B — flatten to a uniform periodic grid (reference implementation)

`Src/flattenAMRFile.cpp` already exists. Flatten to the finest level (or a chosen
analysis level `L`), then the domain is a uniform grid and TTK handles it
directly: `ttk::PeriodicImplicitTriangulation` plus `FTRGraph` / `ttkMergeTree` /
`ttkReebGraph`. Zero new topology code.

**It does remove the coarse-fine crack**, because the output has no coarse-fine
interfaces at all — `WriteSingleLevelPlotfile` of one `MultiFab` over
`geom.Domain()`. That class of defect cannot occur. Four caveats, though, before
treating it as a clean reference:

1. **The default throws the refinement away.** `output_level` defaults to 0. In
   `PltFileManager::fillPatchFromPlt` that hits the `lev0rr == 1` branch — a
   plain `FillPatchSingleLevel` from plt level 0 — and every fine level is
   ignored. It must be run with `output_level=<finestLevel>`. Above level 0 the
   loop over `pltlev` progressively overlays each plt level, `break`ing when the
   resolution matches, so each region ends up sourced from the finest plt level
   at or below `output_level` that covers it. That part is right.

2. **It does not fix the periodic crack.** The periodic node-placement defect in
   CHECK 1 is a property of `isosurface.cpp`'s coordinate handling at the wrap,
   not of the AMR hierarchy. It is level-count-independent and survives
   flattening untouched. Flattening buys the coarse-fine crack only.

3. **It is not the AMR function.** Where the hierarchy held only coarse data, the
   flattened field is an `mf_cell_cons_interp` upsampling (`interp_type`
   defaults to 1). That interpolation is slope-limited, so it does not
   manufacture new extrema inside a refined block — good — but it does put r^d
   distinct values on a linear ramp where the simulation only ever computed one.
   Fine-scale critical points found in coarsely-resolved regions are artifacts of
   the interpolation, not features of the field. They should be low-persistence
   and die under simplification, but they inflate both the graph and the cost.
   **[DECIDE]** whether the reference answer is defined on the interpolated field
   or on the AMR field — these are different functions, and Route A cannot
   reproduce Route B exactly if it is the latter. Agreement should be stated as
   "same graph after persistence simplification at threshold T", not "identical".

4. **Check the periodicity flag.** Plotfiles do not record periodicity — that is
   why `isosurface.cpp` takes `is_per` from `ParmParse`. Inside
   `fillPatchFromPlt` the `dummyBCRec` is `int_dir` only where
   `a_level_geom.isPeriodic(idim)` is true, and `foextrap` otherwise. If the
   Geometry from `pltData.getGeom()` does not carry periodicity, the
   interpolation at the domain faces first-order-extrapolates instead of
   wrapping, perturbing a boundary layer one coarse cell thick — exactly where
   the periodic topology lives. Worth verifying on an actual run before trusting
   any periodic reference answer.

Cheap improvement while in there: `nvars = plt_vars.size()` flattens *every*
component in the plotfile. A Pele plotfile can carry 50-100. Adding the
`comps`/`sComp`/`nComp` subsetting that `isosurface.cpp` already has may be the
difference between feasible and not, since memory grows by the full AMR
compression factor (often 10-1000x) on top of that.

With those addressed it is the right **verification gate** — the same pattern as
the 1-D suite for CAMR — including on a deliberately over-refined hierarchy,
where the AMR answer must not change when refinement is added to an
already-resolved field.

### An aside worth having either way: measure watertightness directly

"Watertight" is a property of a surface, not of a dataset, and it is directly
testable rather than reasoned about. Build a map from each sorted vertex-index
pair to the number of incident triangles; on a closed surface **every edge has
exactly two**. Count the exceptions and report them.

That is ~20 lines, it is exact, and it is a far stronger check than CHECK 1's
edge-length heuristic — which detects mis-placed nodes but says nothing about
cracks. Run it on the current AMR output and it will quantify the coarse-fine
seam immediately; run it on the flattened output and it will settle the question
in this section empirically instead of by argument.

It also covers the one thing flattening does *not* guarantee: marching-cubes
face ambiguity. The standard `edgeTable`/`triTable` pair resolves ambiguous
faces by a fixed convention rather than by the asymptotic decider (Chernyaev's
MC33), so even on a perfectly uniform grid the extracted surface can be closed
but topologically wrong — a spurious tunnel or a merged pair of sheets. For an
area integral that does not matter. For a Reeb graph it is the answer. If the
edge check passes but the graph looks wrong, this is the next suspect.

### Route C — hierarchical / local-global merge trees

Compute local merge trees per box, then merge hierarchically. This is what the
LBL group (Bremer, Weber, Pascucci, Day, Bell) did for LMC combustion AMR data —
the direct ancestor of PeleLMeX — for tracking burning structures. Distributed
and scalable. It gives a merge tree, not a Reeb graph, so it is Route C only if
the [DECIDE] in §3 comes out "merge tree".

### Route D — build the mesh here, do the topology in TTK

Worth taking seriously as the smallest project that reaches a correct answer:
implement only the AMR dual mesh with stitch cells (§A steps 1-4), write it as a
VTK unstructured grid with the periodic vertex identification baked into the
connectivity, and let TTK compute the merge tree / Reeb graph. TTK's periodic
support is limited to *uniform* grids, but a periodic **unstructured** mesh needs
no special support at all — the wrap is already in the connectivity. This is the
single highest-leverage observation in this note.

---

## 5. Suggested staging

Each stage is independently useful and independently checkable.

1. **Route B end to end.** Flatten a small case, run TTK, get a merge tree and a
   Reeb graph. Establishes reference answers and file formats, subject to the
   four caveats in §4B. No new C++.
2. **Run `amrToFE.cpp` as it stands** on a multilevel case with `connect_cc=1`,
   and verify the volume sum against the domain volume, per level pair. This
   answers "is it watertight" by measurement. Nothing is built on it until this
   passes.
3. **Periodic node identity.** Canonicalise the IntVect modulo the domain
   (§4A.1). Verify: on a single-level periodic case the element count and volume
   sum must be unchanged, and the node count must drop by exactly the duplicated
   slab.
4. **Tet output.** Freudenthal subdivision with a consistent parity rule,
   degenerate tets dropped, written as VTU. Verify volume again, and feed to TTK.
   On a single-level case the graph must match stage 1 exactly.
5. **Two levels, then N.** Verify the graph does not change when refinement is
   added to an already-resolved field — modulo persistence simplification, per
   §4B.3.
6. **Distributed.** Only if the target size demands it (§6). Global Morton IDs,
   no single-rank map.
7. **Watertight isosurface as a by-product.** Extract from the same complex.
   `isosurface.cpp` then becomes the comparison gate rather than the engine, and
   its CHECK 2 — currently inactive because the periodic clip was never
   implemented — should finally have something to check.

## 6. Open decisions

- **[DECIDE]** Merge tree, per-component Euler characteristic, or full Reeb graph
  with loop surgery (§3). Sets the scope of everything. Only the third needs the
  loops; the first two need the dual mesh but no sweep over the domain topology.
- ~~Require `ref_ratio == 2`, or build the 4:1 stitch tables.~~ Retired:
  `amrToFE.cpp`'s degenerate-hex stitching needs no case table and no ratio
  restriction.
- **[DECIDE]** Emit one output or two — identified-node topological mesh for the
  sweep, duplicated-slab mesh for viz (§4A.1).
- **[DECIDE]** TTK as a library dependency, or topology implemented in-tree.
- **[DECIDE]** Which directions are periodic in practice, and whether any case
  mixes periodic with physical boundaries (the dual mesh must terminate cleanly
  at a non-periodic face — half-cells at the boundary, or clip).
- ~~Cell-centred vs node-centred dual vertices.~~ Settled: cell-centred,
  `amrToFE.cpp`'s `connect_cc=1`.
- **[DECIDE]** Target problem size. Above ~10^9 cells the collate-to-rank-0
  pattern is out and stage 4 becomes mandatory rather than optional.

---

## References

- Weber, Kreylos, Ligocki, Shalf, Hagen, Hamann, Joy, *Extraction of Crack-free
  Isosurfaces from Adaptive Mesh Refinement Data*, VisSym 2001.
  https://escholarship.org/uc/item/89k672nm
- Weber, Childs et al., *Efficient parallel extraction of crack-free isosurfaces
  from AMR data*, LDAV 2012. https://ieeexplore.ieee.org/document/6378973/
- Tierny, Gyulassy, Pascucci et al., *Loop surgery for volumetric meshes: Reeb
  graphs reduced to contour trees*, IEEE TVCG 2009.
  https://hal.science/hal-01211176
- Gueunet, Fortin, Jomier, Tierny, *Contour Forests: fast multi-threaded
  augmented contour trees*, LDAV 2016.
  https://julien-tierny.github.io/stuff/papers/gueunet_ldav16.pdf
- Bremer, Weber, Pascucci, Day, Bell, *Analyzing and tracking burning structures
  in lean premixed hydrogen flames*, IEEE TVCG 16(2), 2010.
  https://dav.lbl.gov/archive/Publications/2010/LBNL-2276E.pdf
- Tierny et al., *The Topology ToolKit*, IEEE TVCG 2018.
  https://topology-tool-kit.github.io/ — `ttk::PeriodicImplicitTriangulation`,
  `ttk::ftr::FTRGraph`, `ttkMergeTree`, `ttkReebGraph`.

---

## 7. Side questions

### 7.1 Rewriting `isosurface.cpp` on top of `amrToFE.cpp`

Mostly a deletion exercise, with one real new piece.

**What goes away.** The fine-covered mask and the `mask<0` bail; `g1box` /
`rm_external_elements`; the geometric node welding and its broken `operator<`
(§2.3); the `PCInterp` ghost-fill inconsistency (§2.4); the periodic
node-placement defect (§2.6); and `Collate`'s geometric uniquify, which becomes
a plain index remap because connectivity is already global by Node key. Call it
~700 lines removed.

**What replaces `Polygonise`.** Contour the complex, not the grid. Subdivide each
(possibly degenerate) hex into tets and run **marching tetrahedra**: 16 cases,
three of them distinct — empty, one vertex separated (a triangle), two-two (a
quad, split into two triangles). No case tables worth the name, ~100 lines
against `Polygonise`'s ~400 of `edgeTable`/`triTable`.

It also **retires the ambiguity caveat in §4B**: marching tets has no ambiguous
faces, so the extracted surface is topologically consistent by construction.
There is no MC33 question to answer.

**The structural win.** Key the vertex cache on the **unordered pair of Node
keys** — `{(level,iv), (level,iv)}` canonicalised by `Node::operator<` — rather
than on `Edge{IV_l, IV_r}` plus coordinate welding. A crossing on a shared edge
is then computed once and found by identity from both sides, at any level and
across the periodic wrap. That single change is what makes the surface watertight,
and the same table is what the Reeb graph's arcs will reference, so the two tools
share it.

**Real work items, honestly:**

1. **Tet subdivision with consistent face diagonals** — Freudenthal with a parity
   rule on the base index, degenerate tets dropped, and care where levels meet
   because the two sides index differently. This is the only hard part, and it is
   the *same* work item the Reeb graph needs. Do it once.
2. Field values on nodes — largely present already (`fileData`/`FillVar`).
3. The single-rank `std::map` (§4A.4) is inherited, not fixed.

**Risks to weigh, not ignore:**

- **Sliver tets near coarse-fine interfaces.** Degenerate hexes subdivide into
  poorly-shaped tets, so triangle aspect ratios in the transition band will be
  worse than today's. Areas are unaffected (exact for the linear interpolant),
  but **normals and curvature are not** — and `curvature.cpp` consumes this
  output. Check that downstream before committing.
- **The surface changes.** Piecewise-linear-on-tets is a different interpolant
  from marching cubes on a dual cell, and the result depends on the diagonal
  choice. Expect percent-level area differences. Validate against an analytic
  case (a sphere, `4*pi*r^2`) rather than against the current code.

### 7.2 Extending the surface from the half-cell inset to the walls

First, a scoping point: **periodic directions need no extension.** The half-cell
inset there is not a clipping choice, it is the missing wrap. Once §7.3 is in,
the dual cell between the last and first cell centres exists and the contour
closes through it. There is no wall to reach. CHECK 2's "the periodic clip was
never implemented" is satisfied by identification, not by clipping.

So this question is only about **genuine non-periodic walls**, where cell-centred
data really does leave the complex spanning `[dx/2, L-dx/2]`.

In `amrToFE.cpp`'s index-space framework the extension is clean: grow the box by
one cell in each non-periodic direction, stamp a Node at the ghost index with a
value taken by `foextrap` from its interior neighbour, and **place that node's
coordinate at the wall** rather than at the ghost cell centre — otherwise the
complex overhangs the domain by `dx/2`. That is a special case in coordinate
emission only, keyed on the IntVect lying outside the domain; connectivity,
element generation and contouring are all untouched.

With constant extrapolation the field has no wall-normal variation in that layer,
so `grad f` is tangential there and **the surface meets the wall at 90 degrees**.
That is a defensible placeholder, but it is a choice with visible consequences —
see §7.2.1. It should not be the permanent answer without measuring what it costs
on real cases.

Two consequences to plan for:

- The fiddly part is not the faces but the **edges and corners of the domain**,
  where two or three extrapolations meet. Index space makes this bookkeeping
  rather than geometry, but it still has to be written and tested.
- **Reported measures will change**, legitimately, and by more than rounding on
  wall-dominated cases. Re-validate against an analytic case rather than against
  the current numbers. And note for §3: a surface with boundary has
  `chi = 2 - 2g - b` for `b` boundary loops, so the per-component genus
  bookkeeping needs the boundary count once walls are in play.

Estimate: ~150-250 lines, plus the extrapolation decision.

#### 7.2.1 What zero-gradient extrapolation actually does to the surface

Worth being concrete, because the artifacts are specific and two of them are
measurable before committing to anything.

**The surface gets a prismatic collar.** With `f` constant in the wall-normal
direction across the boundary layer, the surface there is a straight extrusion of
the contour on the inset plane, out to the wall. Consequences:

- **Area gains exactly `P * dx/2`**, where `P` is the length of the contour's
  intersection with the inset plane. The bias is one-signed and exactly
  computable, which is useful: it can be reported rather than guessed at.
- **A curvature ridge at a fixed `dx/2` from every wall.** The join is C0 but not
  C1, so the extrusion meets the interior surface in a kink. `curvature.cpp` will
  see a spurious high-curvature band standing off each wall at the same distance.
  Of all the artifacts this is the one most likely to be mistaken for physics.
- **Exact ties in the sweep.** Every boundary node equals its interior parent, so
  the Reeb/merge sweep sees `N_wall` tied pairs. With a consistent tie-break
  (the `(level,iv)` key) these are zero-persistence pairs — they inflate the
  sweep and the graph but do not change it after simplification. The exception
  worth naming: an extremum adjacent to a wall becomes a **plateau extremum** of
  two or more vertices, which needs flat-region handling rather than a tie-break.

**The physics failure is bigger than any of those.** Zero normal gradient is
*exactly right* for a symmetry, slip, or adiabatic wall, and *exactly wrong* for a
Dirichlet wall. On an isothermal wall a temperature isosurface should terminate on
the wall wherever `T_wall` crosses the isovalue; zero-gradient will simply
suppress that intersection. No amount of care in the extrapolation fixes a
boundary condition it was never told about.

#### 7.2.2 How to decide, cheaply

The error is bounded by how obliquely the surface meets the wall. Where the true
surface is already near-perpendicular, zero-gradient is nearly exact; where it is
near-tangent, the collar is badly wrong in both area and in where the
wall-intersection curve sits.

So: **measure the distribution of incidence angles on the inset plane.** The
contour and the wall normal are both already in hand, so this is a histogram over
the boundary elements and costs nothing. If the mass sits near 90 degrees,
zero-gradient is fine and the question is closed. If there is a tail toward
tangency, it is not, and the boundary condition has to be supplied.

Report alongside it the collar's area contribution `P * dx/2` as a fraction of
total area. Sub-percent means stop worrying; double digits means the choice is
load-bearing.

**Analytic test that isolates the wall treatment:** a sphere centred *on* a wall.
The exact answer is `2*pi*r^2`, the incidence angle is 90 degrees everywhere, and
zero-gradient should be near-exact — so this validates the machinery. Then
re-run with the sphere centre offset outside the domain so the surface meets the
wall obliquely; the error there is the error the scheme will make on real data.

#### 7.2.3 Options, in increasing order of honesty and cost

- **Zero-gradient (`foextrap`).** Right for symmetry/adiabatic/slip. The
  placeholder. Artifacts as above.
- **Limited linear extrapolation.** Use the two interior cells, then clamp the
  result to their range. Removes the curvature ridge and the ties, and the
  limiter prevents an overshoot from manufacturing a crossing the data does not
  support — the same slope-limiting idea as `cell_cons_interp`. Degenerates to
  zero-gradient exactly where it is limited, so it is never worse.
- **Supplied boundary conditions.** The only option that is *right* at a
  Dirichlet wall. Plotfiles carry no BC metadata, so this has to be a
  user-supplied table, e.g.
  `wall_bc.<comp>.<face> = neumann | dirichlet <value>`.
- **Do not extrapolate; clip and record.** Keep the complex at the inset plane,
  accept an open surface, and report the boundary explicitly: area, plus the
  wall-intersection perimeter as a separate number, plus `b` for the
  `chi = 2 - 2g - b` bookkeeping of §3. This invents nothing. It understates area
  by exactly the collar term, and it looks wrong in pictures.

**[DECIDE]** which of these is the default. The recommendation is: limited linear
as the default, zero-gradient as an option, a BC table when a case needs it, and
the §7.2.2 diagnostics printed either way so the choice is never silent.

### 7.3 Yes — use `FillBoundary(geom.periodicity())` on the `MultiNodeFab`

This is the right mechanism and it is better than the modulo-canonicalisation I
proposed earlier. The reason is the same property that makes the whole tool work:
**identity is data, not position.**

`FillBoundary` copies Node *values* into ghost cells. It does not, and cannot,
shift an `IntVect` that is stored inside the copied object. So after a periodic
fill the ghost at `domain.bigEnd()+1` contains `Node(domain.smallEnd(), lev,
VALID)` — the far-side node, by identity, with no arithmetic performed. Elements
built across the face reference the far-side node automatically. The
identification is a side effect of the copy semantics.

The machinery is already in the file: `pData.FillBoundary(geom[lev]->periodicity())`
is used for the *field* data in the `nGrowPer` path. What is missing is the same
call on the `MultiNodeFab`. `Node` is trivially copyable (`int`, `IntVect`,
`int`, enum — no pointers, unlike `isosurface.cpp`'s), and the file already
instantiates `FabArray<BaseFab<Node>>`, so this should compile as-is.

Four things to get right:

1. **Ordering against the coarse-parent stamping.** The sequence has to be:
   fill valid nodes, periodic `FillBoundary`, *then* stamp coarse parents into
   fine ghost regions, then mark `COVERED`. Where a fine grid abuts the periodic
   face *and* has a coarse neighbour, the stamp and the wrap both target the same
   ghost cells. `GetBndryCells` already does its own `periodicShift`, so the two
   paths overlap and the precedence needs to be decided deliberately rather than
   inherited from statement order.
2. **Make periodicity explicit.** `geom[lev].reset(new Geometry(ProbDomain, &rb))`
   passes no `is_per`, so it picks up `Geometry`'s ParmParse statics
   (`geometry.is_periodic`, as the help text documents). Take an explicit
   `is_per` array the way `isosurface.cpp` does — this is exactly the flag that
   decides whether the answer is a torus or a slab, and it should not arrive
   implicitly.
3. **`nGrow = 1` suffices** for the wrap: the dual cell straddling the face needs
   only the far-side cell centre.
4. **Check for double-generated wrap elements.** Elements are enumerated over the
   *grown* box and deduplicated by `std::set<Element>`, whose `operator<` compares
   node slots in fixed order without canonicalising rotation. Across the wrap the
   two views should produce identical slot orders — so the set should dedup
   correctly — but that is worth confirming rather than assuming. Restricting
   element generation to bases inside the **valid** box makes it generated-once
   by construction and removes the question.

Keep `nGrowPer` as the visualization path: on a torus the coordinates cannot be
single-valued (a wrap element's corners differ by `L`), so the duplicated slab
remains the right thing for pictures while the identified mesh is the right thing
for the sweep. That is the [DECIDE] already recorded in §4A.1.

---

## 8. One tool or two?

Neither. **One shared library, two thin drivers.**

### Why not one tool

The two algorithms share the complex and nothing else. Contouring at one isovalue
touches only the cells straddling it: local, `O(surface)` memory, no
communication beyond a final gather, seconds of runtime. The Reeb sweep visits
every vertex in global sorted order: `O(N)` memory, a total order across ranks,
union-find, tree combination, communication throughout, minutes to hours.

Folding them into one binary would put a global sweep behind every isosurface
request. That is a serious regression for the existing chain — `curvature.cpp`,
`binMEF`, `stream.cpp` and the whole `*MEF` family in `Src/` consume this output
routinely, and the MEF format is load-bearing across a dozen tools.

There is one real sense in which the Reeb tool *should* emit contours: given the
graph, "extract the contour at value `v` belonging to arc `A`" yields an
**arc-labelled** isosurface, where every connected component carries its
identity in the graph. That is strictly more than `isosurface.cpp` produces and
it is exactly what feature tracking needs. But it is a different product from
"give me the surface at `v`", and it should not be the only way to get one.

### Why not two independent tools

The complex, the tet subdivision and the edge-keyed vertex cache are the hard,
correctness-critical parts. Two copies will diverge — and divergence between two
paths that ought to agree is precisely the failure already in the tree twice
over: the coarse-fine seam (§2.2) and the periodic node placement (§2.6).

### The shape

- **`AMRComplex`** — plotfile to conforming multilevel complex: node identity,
  periodic `FillBoundary` identification (§7.3), coarse-parent stamping, tet
  subdivision, node values, global IDs, wall treatment (§7.2). The `amrToFE.cpp`
  core, modernised. This is where the difficulty lives.
- **`Contour`** — marching tets on the complex, vertex cache keyed on the
  unordered Node pair. ~200 lines.
- **driver `isosurface2`** — complex + contour + measures + MEF/XDMF. No sweep.
- **driver `reebGraph`** — complex + sweep, optionally calling `Contour` for
  arc-labelled output.

### Sequencing, and the reason for it

**Build `AMRComplex` first and prove it with the cheap consumer.** The isosurface
driver is the validator for the complex, not merely another client of it: it runs
in seconds, it has analytic references (sphere area; hemisphere on a wall per
§7.2.2), it has the existing CHECK 1 and CHECK 2 — CHECK 2 finally being live
once identification is in — and it has a downstream chain that will complain
loudly and specifically if the surface is wrong. Get that green before writing a
line of sweep code, for the same reason the 1-D suite runs before the 2-D pipe
break.

Retiring `isosurface.cpp` then becomes a consequence rather than a goal: it goes
when `isosurface2` reproduces it on the analytic cases, passes both checks, and
`curvature.cpp` is happy with the tet-derived normals (§7.1, the sliver risk).
Two features have to be carried across rather than dropped: the MEF format
contract with the `*MEF` toolchain, and `build_distance_function` (SDFGen /
`make_level_set3`), which is independent of all of this and simply needs a home.

---

## 9. What is actually needed before starting `AMRComplex`

Very little more planning. Two measurements, one decision, one page — call it half
a day — and then write code. The risk of planning further is that everything past
this point would be planning against an unverified assumption.

### 9.1 Two things to find out first, because they are cheap and they fork the plan

**(a) Does `amrToFE.cpp` still build and run?** (~1 hour.) It uses `DataServices`
/ `AmrData` and `BL_SPACEDIM`, which are legacy, and `Src/GNUmakefile` currently
has other targets commented out. Determine which code is live before designing
around it. The fork it settles: port in place, versus lift the `Node` identity and
coarse-parent stamping into new code and leave the old tool alone.

**(b) Is the construction actually gap-free?** (~2 hours, and the single
highest-information action available.) Everything in §4A, §7 and §8 rests on the
degenerate-hex stitching tiling the coarse-fine interface exactly. That has been
read, not measured. On a synthetic two- and three-level hierarchy:

- **Sum the element volumes and compare with the domain volume**, per level pair.
  A gap or an overlap at the transition shows up immediately and unambiguously.
- **Histogram the elements by distinct-node-count** — 8 hex, 6 prism, 5 pyramid,
  4 tet, and anything else is a surprise worth understanding. This is not just a
  check: it enumerates exactly which degenerate shapes the transition band
  actually produces, which is the input to the tet-subdivision work in §7.1. It
  turns "write a Freudenthal subdivision with a parity rule" from a guess into a
  finite list.

If (b) fails, the plan changes materially — Weber's explicit stitch-cell tables
come back, and §4A's "no case table at all" claim is wrong. Better to know now
than after `AMRComplex` has an interface built on it.

### 9.2 Settled: the complex must be distributed

Requirement stated: large files, parallel throughout. That closes the earlier
[DECIDE], and it has to be honoured **up front rather than retrofitted**, because
the global-ID scheme is the contract both drivers share. See §10 — the API port
is small, but the parallelism is not a port at all.

### 9.3 One page, not a design document

A header sketch of `AMRComplex`: what it owns, what it hands out, what each
driver calls. Enough to keep `isosurface2` and `reebGraph` honest about sharing
it, and no more. Specifically worth pinning down: whether it hands out an explicit
tet list plus node arrays, or an abstract vertices/edges/incidence interface. Both
consumers need explicit tets, so explicit is the simpler answer unless there is a
reason to keep the sweep cell-type agnostic.

### 9.4 What does *not* need planning now

Deferred and genuinely separable — the complex is correct without any of them,
and each can be added against a working, validated base:

- Wall extrapolation (§7.2). Set aside deliberately.
- The Reeb sweep itself, loop surgery, persistence simplification.
- Output formats beyond whatever the volume check needs.
- Tet subdivision parity-rule details — write it, test it on the two-level case,
  and let §9.1(b)'s histogram tell you the case list.

---

## 10. Porting `amrToFE.cpp`, and the part that is not a port

### 10.1 The API port is small; agreed

`DataServices`/`AmrData` to `PlotFileData`, `BL_SPACEDIM` to `AMREX_SPACEDIM`,
raw loops to `MFIter`. Models in the tree:

- **`isosurface.cpp`** is the front-end to copy: `PlotFileData`, component
  selection by name, per-level `Geometry` built with an explicit `is_per` from
  `ParmParse` (§7.3.2 wants exactly this), `FillPatchTwoLevels` for level-to-level
  data, and the fine-covered mask with its periodic-shift pass.
- **`avgPlotfiles.cpp`** for the modern multi-level idiom — `Vector<Geometry>`,
  `refRatios`, `WriteMultiLevelPlotfile`.
- **`conditionalMean.cpp`** is *not* a porting exemplar: it is still on
  `DataServices`/`AmrData` itself. Its parallel pattern — `ReduceRealSum` into a
  fixed-size bin array — is the easy case and does not transfer to unstructured
  output.

### 10.2 The parallelism is not a port, because it is not there

`amrToFE.cpp` is **single-rank-correct only** today, and the API has nothing to do
with it:

- `nodeMap[ifab(iv,0)] = cnt++` builds a **rank-local** `std::map<Node,int>` with
  a **rank-local counter**. Two ranks both start at 0, so node IDs collide.
- The `Collate(tmpData,connData,MYLEN)` call at line 816 is **commented out**.
  Nothing gathers, nothing renumbers.

So the serialness is in the data structure, not the interface, and porting the API
will not touch it. Worth knowing before assuming the effort is small: the port is
small, the parallel rewrite is the actual work.

Note also that **no tool in `Src/` currently writes a distributed unstructured
mesh.** MEF, Tecplot, and the XDMF path in `isosurface.cpp` all gather to
`IOProcessor` first. This is new code, not a port.

### 10.3 The design falls out of the node identity

The good news is that the hard problem is already solved by the existing
representation: **`(level, IntVect)` is a rank-independent global key.** No map is
needed to *identify* a node — only to *number* it. So:

1. Add an `int gid` field to `Node`. It stays trivially copyable, so it still
   rides `FillBoundary` and `ParallelCopy` (§7.3).
2. Fill `VALID`/`COVERED` per level as now.
3. Periodic `FillBoundary` on the `MultiNodeFab` — the wrap identification of
   §7.3, free.
4. Coarse-parent stamping into fine ghosts. A fine box's coarse parent may live on
   another rank, so this is a `ParallelCopy` from a coarsened-fine `BoxArray` —
   the standard AMReX two-level pattern, and the idiom `isosurface.cpp` already
   uses via `FillPatchTwoLevels`.
5. **Number from the BoxArray index** (see §10.5) — contiguous and
   rank-count-independent:

   ```
   gid = levelOffset[lev] + boxOffset[lev][i] + localIndex(iv within box i)
   ```

   where `i` is the box's index in `grids[lev]` and `localIndex` walks the valid
   box lexicographically, skipping `COVERED`.
6. **`FillBoundary` + coarse-to-fine `ParallelCopy` again** to propagate `gid`
   into ghost and stamped positions. Now every corner an element can reference
   carries a globally valid ID, locally.
7. **Generate elements with the base cell restricted to the valid box.** Each
   element is then produced exactly once, by exactly one rank — no
   `std::set<Element>`, no dedup, no gather. This also retires the open question
   in §7.3.4.

Net effect: no global map, no global element set, `O(local)` memory, and one
small all-reduce of communication. The `Collate` machinery is deleted rather than
repaired.

### 10.4 Parallel output is the other half

Gathering to `IOProcessor` fails at the target size regardless of how well the
complex is built. Least-resistance path: each rank writes its own `.vtu` piece
plus one small `.pvtu` master. ParaView and TTK read it directly, it needs no
collective I/O layer, and it matches the element-generated-once ownership from
§10.3.7 with no extra bookkeeping. HDF5/XDMF is the alternative if a single file
is wanted — `convert2hdf5.cpp` and `isosurface.cpp`'s XDMF writer are starting
points, but both currently assume gathered data.

### 10.5 Numbering: the BoxArray index is already rank-independent

No sort is needed, and no Morton code either. **A `BoxArray` is a globally
replicated object**: `ba[i]` is the same box on every rank for any `NProcs`. Only
the `DistributionMapping` depends on the rank count. So the BoxArray's own index
order is already a deterministic global ordering of boxes, and the numbering in
§10.3.5 follows directly.

The only piece needing communication is `boxOffset[lev][i]`, the count of `VALID`
nodes in boxes `0..i-1` at that level:

- each rank fills its own boxes' counts into a `Vector<int>` of size `ba.size()`,
  zeros elsewhere;
- one `ParallelDescriptor::ReduceIntSum` over that vector;
- a local prefix sum.

Every rank then holds an identical offset table. Cost is one all-reduce of
`ba.size()` ints — negligible even at `1e5` boxes per level. The per-box `VALID`
count is itself rank-independent, since `VALID`/`COVERED` depends only on the grid
hierarchy and the counting is over valid boxes, which partition the level.

This is strictly better than either option previously recorded here: contiguous
*and* reproducible, so output files **are** bitwise comparable between a 1-rank
and a 64-rank run. The earlier [DECIDE] between an `Exscan` offset
(rank-dependent) and a Morton code (sparse, needing sort-and-compaction) is
retired — neither is necessary.

**Answering the literal question about a sort routine:** `DistributionMapping`
does expose one —

```cpp
static std::vector<std::vector<int> >
makeSFC (const BoxArray& ba, bool use_box_vol=true,
         int nprocs=ParallelContext::NProcsSub());
```

It returns box indices in space-filling-curve order, bucketed by rank; passing
`nprocs=1` should yield a single bucket holding the pure SFC permutation, which is
deterministic and independent of the run's actual rank count. That is not needed
for correctness — plain BoxArray index order already is — but it is worth
considering as an *optional* permutation applied in place of the raw index order,
because SFC ordering gives spatial locality in the node numbering, which helps
cache behaviour and output compression. Verify the `nprocs=1` semantics before
relying on them.

Regression gating should still prefer invariants — volume sum, element-shape
histogram, surface area, graph structure — over file diffs, since those catch
real defects rather than incidental reordering. But with this scheme a file diff
across rank counts is available as a cheap additional check, which it would not
have been otherwise.

---

## 11. Can `AMRComplex` be distributed, and written plotfile-style?

Distributed: yes, by construction. Plotfile-style: **yes for the half that
matters, and the other half should not be stored at all.**

### 11.0 A note on `Collate`

The §10.5 scheme is not quite what `Collate` does. `Collate` is `MPI_Gatherv` in
**rank order**, so the assembled ordering is rank-count dependent — which is
precisely why `isosurface.cpp`'s output is not reproducible across rank counts
today. The BoxArray-index scheme replaces that ordering rather than reproducing
it. The *shape* of the two is similar (per-owner counts, then offsets), which is
probably the resemblance; the difference is that the offsets are keyed on a global
box index rather than on a rank.

### 11.1 The complex is a view, not data

Worth stating before designing a file format: **the complex is a deterministic
function of the grid hierarchy, `is_per`, and a small parameter set** (wall
treatment, subdivision rule). The hierarchy is already in the plotfile. So in the
general case the mesh does not need to be stored at all — it is regenerated on
read, exactly, from inputs that are already on disk.

Persisting it is therefore an optimisation (avoiding regeneration across many
analyses of one snapshot) or an interchange requirement (handing it to TTK), never
a correctness requirement. That should shape the default: **write the parameters,
not the mesh.**

### 11.2 The node-centred half maps onto a plotfile exactly

A node *is* a cell — `(level, IntVect)`. So every per-node quantity is naturally
cell-centred multi-level data, and `WriteMultiLevelPlotfile` writes it distributed
with no new I/O code:

- node type (`VALID` / `COVERED`),
- the analysis fields themselves,
- and, later, **the topology results**: merge-tree branch ID, Reeb arc ID,
  segmentation label, persistence.

That last point is the important one. **The segmentation — which cell belongs to
which arc or branch — is exactly a plotfile**: integer labels, cell-centred, on
the same hierarchy the data came from. It round-trips through `PlotFileData`, it
opens in Amrvis/ParaView/yt with no conversion, and it feeds straight back into
the existing tools (`conditionalMean` conditioned on segment ID, and so on). It is
also how the LBL burning-structures work presented its results (§4C).

Two details:

- **Do not store `gid`.** It is reconstructible from the box offsets, and
  plotfiles store `Real` — in a single-precision build, integer IDs above `2^24`
  stop being exact. Labels have small range and are safe; `gid` is neither needed
  nor safe.
- **`COVERED` cells** still occupy storage at their level. That is normal
  plotfile behaviour and viz tools already resolve the hierarchy correctly; the
  type field just makes the complex's own view explicit.

### 11.3 The connectivity does not map onto a plotfile, and need not

Transition elements have corners on two levels, so they belong to no single
level's `BoxArray`. They are irreducibly unstructured.

But §10.3.7 already gives every element a unique owner — the rank holding the
valid box containing its base cell — so the element list is partitioned by
`(level, box index)` with no extra bookkeeping. When it does need to go to disk,
that maps directly onto per-rank `.vtu` pieces with a `.pvtu` master (§10.4), or
onto an HDF5 dataset with per-box offsets if a single file is wanted.

### 11.4 The resulting two-file story

- **`<name>_complex`** — an AMReX multi-level plotfile: node type, fields, and
  topology labels. Native, distributed, visualisable, and the natural output of
  the whole pipeline.
- **`<name>_mesh.pvtu`** (+ pieces) — the unstructured view, for TTK and
  ParaView. Optional, regenerable, written only when asked for.

One gap to note against §7.2: wall boundary-layer nodes lie *outside* the domain
and so have no cell in the hierarchy. They would need either a grown domain in the
written file or separate storage. Deferred with the rest of the wall work.

---

## 12. Result of the gate (section 9.1(b))

Answered, on branch `amr-complex`. **The construction is gap-free**, but not
under the ownership rule §10.3.7 proposed, and the difference is the whole
finding.

### 12.1 The naive ownership rule loses the stitch layer

"An element belongs to the box whose valid region contains its base cell" is
wrong. The element bridging a coarse-fine interface on the **low** side has its
base in a fine ghost stamped from the coarser level: no fine box holds it as
valid, and the coarse box that would hold it has a `COVERED` corner and is
rejected. Nobody generates it. Measured cost on a two-level 32^3 case:
**1.93% of domain volume**, a gap one stitch-layer thick around every refined
region — large enough to matter, small enough to be mistaken for something else
if it had not been measured.

The rule that works: base range is the valid box grown by one **on the low side**,
with the lowest-indexed visiting box owning that shell, and periodic wrap images
excluded so the low face is not counted twice. Still generated-once, still
deterministic, still no communication.

Getting this wrong was worth the cost of finding it: it is exactly the class of
defect that produces a plausible-looking mesh and a quietly wrong answer.

### 12.2 Verified

15 combinations — 1, 2 and 3 levels; `ref_ratio` 2 and 4; fully periodic,
non-periodic, and mixed — all with **relative volume error exactly zero**.

Expected volume is the domain in periodic directions and the domain less one
level-0 cell in non-periodic ones, since cell-centred data leaves the complex
inset by `dx/2` at each wall (§7.2).

Two incidental results worth keeping:

- **Compensated summation is necessary, not fastidious.** Naive summation costs
  ~5e-10 relative at 2e6 elements. A tolerance loose enough to absorb that would
  hide a real gap of the same size, and the gate has to stay sharp as problems
  grow. With Kahan summation the residual is exactly zero at every size tested.
- **`ref_ratio` is unrestricted.** Ratio 4 verifies exactly as ratio 2 does,
  confirming §4A's claim that degenerate-hex stitching needs no case table. The
  `[DECIDE]` on requiring 2:1 stays retired.

### 12.3 The shape histogram, i.e. the tet-subdivision case list

Across every case, elements are only ever **hexes, prisms and pyramids**. No
tets, and nothing unexpected. On the three-level 64^3 periodic case: 702656
hexes, 23808 prisms, 13056 pyramids.

So the Freudenthal subdivision of §7.1 has three cases to cover, not an open-
ended family. That was the point of asking for the histogram rather than only
the volume.

### 12.4 Improper nesting is now an error, not a silent leak

The first three-level attempt failed with a ~3% residual. Cause: the test
generator produced a level that exactly covered its parent, so the next level's
ghosts had only `COVERED` parents to stitch to and the interface had nothing to
attach to.

That was a defect in the generator, but it exposed a real robustness gap:
nothing detected it. `classifyNodes()` now aborts at construction when a ghost
stamped from level `lev-1` comes back `COVERED`, naming the levels involved.
Prevent the bad state at creation rather than explaining the residual later.

### 12.5 What was built

- `Src/AMRComplex.H`, `Src/AMRComplex.cpp` — the complex.
- `Src/amrComplexVerify.cpp` — driver running the gate.
- `Src/makeTestPlt.cpp` — synthetic properly-nested hierarchies with an
  analytic field (sphere SDF, so a later isosurface has a known area too).
- `Src/GNUmakefile.AMRComplex` — standalone, AMReX-only. Deliberately
  independent of `Make.Analysis` and therefore of PelePhysics, which is not
  checked out in this tree and which `AMRComplex` does not need.

### 12.6 Not yet done

- **MPI is unmeasured.** No MPI in the environment used. The rank-independence
  of the numbering (§10.5) and of the ownership rule (§12.1) is argued but not
  demonstrated. First real test: identical element count, volume and histogram
  at 1 and N ranks, and — since numbering is rank-independent by construction —
  identical output files.
- `ownsBaseCell` calls `BoxArray::intersections` per ghost base cell. Correct
  and hash-backed, but it is the obvious hot spot; restrict it to the ghost
  shell or precompute an ownership mask when it starts to matter.
- Field values on nodes; tet subdivision; `.vtu`/`.pvtu` output; the plotfile
  writer of §11.2.

---

## 13. State of play: the complex and the contour are done and verified

Branch `amr-complex`. Sections 9 through 12 are answered; what follows is what
was actually built and what it now guarantees.

### 13.1 What exists

| file | what it is |
|---|---|
| `Src/AMRComplex.{H,cpp}` | the conforming multilevel complex |
| `Src/AMRContour.{H,cpp}` | marching tetrahedra on it |
| `Src/amrComplexVerify.cpp` | the complex gate |
| `Src/amrContourVerify.cpp` | the surface gate |
| `Src/makeTestPlt.cpp` | synthetic properly-nested hierarchies, sphere SDF |
| `Src/GNUmakefile.AMRComplex` | standalone, AMReX-only, no PelePhysics |
| `Tools/amrComplexMpiTest.sh` | the four-section MPI harness |

### 13.2 What is verified

- **The complex tiles exactly.** Relative volume error zero across 1-3 levels,
  `ref_ratio` 2 and 4, periodic, non-periodic and mixed.
- **It is face-conforming**, which the volume gate cannot see: volume is
  additive over any decomposition, so a T-junctioned complex closes exactly.
  Measured by face incidence in serial and by an XOR parity globally.
- **The surface is closed with Euler characteristic 2** in every case, and the
  area converges at **second order**: 32^3 `-2.799e-3`, 64^3 `-6.985e-4`,
  128^3 `-1.745e-4`, ratios 4.007 and 4.003.
- **All of it is rank-invariant.** 91 of 91 checks, 1 to 8 ranks including odd
  counts, across SFC/KNAPSACK/ROUNDROBIN distributions and varying
  `max_grid_size`. Byte-identical, not merely close.

### 13.3 The subdivision rule, settled

Cone each cell from its **minimum-gid corner**. That puts the min-gid diagonal
on all six faces automatically -- for the three containing the apex it follows
from a corner minimal over the cell being minimal over any face containing it;
for the three opposite it is imposed. Two cells sharing a face see the same
four gids and split it identically, across levels included.

The parity rule §7.1 worried about is not needed and should not be used:
parity is keyed on indices, which is precisely what differs between levels.

### 13.4 Two lessons that cost real time

**Area is not a topology test.** Forcing an inconsistent diagonal produces
10401 boundary edges, 445 non-manifold edges and `chi = -2976`, while
reporting an area error of `-0.0012` -- *better* than the correct rule's
`-0.0016`. A surface riddled with cracks can score better on area than a sound
one. The negative control is checked in (`diagonal_rule=0`) so the gate can be
re-falsified whenever the subdivision changes.

**Never let an array index into a checksum.** The same bug appeared twice: the
face census summed rank-local counts, and the edge parity hashed rank-local
vertex indices. Both times the structure was correct and only the
instrumentation was wrong, and both times the symptom was alarming -- an open
surface reported on a closed one. Standing rule: anything entering a checksum
is named by a gid or a tuple of gids, never by a position in an array.

### 13.5 Next

- The Reeb sweep. It needs the tet edges this already builds, a total order on
  nodes (gid serves, and is rank-independent), and union-find. §3's
  `[DECIDE]` -- merge tree, per-component Euler characteristic, or full Reeb
  graph with loop surgery -- is still open and still sets the scope.
- `isosurface2`: MEF/XDMF output from `AMRContour`, then the comparison that
  retires `isosurface.cpp` (§8).
- **Unmeasured, and needs real data rather than another synthetic test:**
  whether `curvature.cpp` is happy with normals derived from tets in the
  stitch band. Degenerate hexes subdivide into poorly-shaped tets, and area is
  insensitive to that while curvature is not (§7.1).
- Wall extrapolation (§7.2), still deliberately deferred.

---

## 14. Verified state (supersedes 13.5-13.6)

Branch `amr-complex`. All suites pass in 2D and 3D, at 1-8 ranks, under every
DistributionMapping strategy: `Tools/amrBuildAndTest.sh`.

### 14.1 What is built

| file | what |
|---|---|
| `Src/AMRComplex.{H,cpp}` | conforming multilevel complex, 2D and 3D |
| `Src/AMRContour.{H,cpp}` | marching simplices, oriented, with carried components and MEF output |
| `Src/amrComplexVerify.cpp`, `Src/amrContourVerify.cpp` | the gates |
| `Src/makeTestPlt.cpp` | synthetic hierarchies: anisotropic cells, non-cubic domains, linear probe fields |
| `Src/GNUmakefile.AMRComplex` | standalone, AMReX only, `DIM ?= 3` |
| `Tools/amrComplexMpiTest.sh`, `Tools/amrContourMefTest.sh`, `Tools/amrBuildAndTest.sh` | the suites |

### 14.2 What is guaranteed

- **The complex tiles exactly and is face-conforming** -- 1-3 levels,
  `ref_ratio` 2 and 4, periodic/non-periodic/mixed, isotropic and anisotropic
  cells, cubic and elongated domains.
- **The surface is closed, has the right Euler characteristic, and is
  consistently oriented** with normals toward increasing f. Enclosed volume by
  the divergence theorem confirms outward, not merely consistent.
- **Second-order convergence**: 3D areas at 32/64/128 give ratios 4.007 and
  4.003; 2D lengths at 32/64/128/256 give 4.01, 4.00, 4.00.
- **Carried components are exact**: a linear field comes back on the surface
  equal to the vertex coordinate, to zero.
- **Everything is rank-invariant byte-for-byte**, including the MEF file.

### 14.3 The four negative controls

Every one of them was added because the obvious check had already been shown
to pass on a broken mesh. They are the most valuable thing in the tree.

| control | what it breaks | what still passes |
|---|---|---|
| `diagonal_rule=0` | face-diagonal consistency | area -- and *better* than the correct rule |
| `orient_elements=0` | element winding | area *and* closure, both identical |
| `check_vertex_determinism` | element-dependent vertex data | everything, in serial |
| improper nesting abort | stitch with no valid coarse parent | nothing -- it aborts |

The recurring lesson: **area is not a topology test, and closure is not an
orientation test.** Each is blind to the defect the next one catches.

### 14.4 Every MPI defect found was instrumentation

Five in a row -- rank-local face counts, rank-local vertex indices, a
whitespace-split periodicity list, a renamed sed terminator, an unsorted
element list -- plus one real one, the element-frame vertex position. The
structure was right nearly every time; the things measuring it were not.
Standing rule that came out of it: **anything entering a checksum is named by
a gid or a tuple of gids, never by a position in an array**, and any ordering
written to a file is canonicalised explicitly rather than inherited from
traversal.

### 14.5 Next

- **`[DECIDE]` (§3), unchanged and still gating the sweep:** merge tree,
  per-component Euler characteristic, or full Reeb graph with loop surgery.
  The second is nearly free now -- the contour already computes chi globally
  and orientation is settled, so it is per-component labelling on top of what
  exists.
- Run on a real plotfile. Everything so far is synthetic. Watch the nesting
  assertion against hierarchies from actual regridding, and the cost of
  `ownedBaseRegions` and the contour's `std::map` at production size.
- `curvature.cpp` on tet-derived normals in the stitch band -- the sliver risk
  of §7.1, now testable since orientation is consistent, but it needs real
  geometry.
- `.pvtu` output, to retire the gather-to-IOProc in `writeMEF`.
- Wall extrapolation (§7.2), still deliberately deferred.

---

## 15. The actual target, and what it settles

### 15.1 The problem

Characterising the effects and time-dependent dynamics of the thermodiffusive
instability in lean premixed hydrogen flames -- the cellular structures these
flames spontaneously form, whose sizes and separations change dramatically
with mixture and with the turbulence of the approach flow. Configurations
include a doubly periodic flamesheet with inflow/outflow in the third
direction, and growing flame kernels, in 2D and 3D, under different transport
models.

**What is not known, and what the tool exists to explore: what a cell IS.**
There are chemistry arguments about when a local flame is self-sustaining,
notions of back-supported structure, and a set of observational surrogates
that have traditionally stood in for the presence of a flame -- fuel
consumption rate, temperature above a threshold, a flame intermediate. Which
to use, and at what threshold, is open.

A complication specific to these flames: because they are not unity-Lewis,
the progress variables do not align. Contours of fuel, temperature and
intermediates cross one another. Each may tell a different story in detail
while giving a similar global picture -- and whether that is so is itself one
of the questions.

### 15.2 What that settles

**The spec is not "compute the topology of field X". It is "compare candidate
definitions within one frame."** So the pipeline must be field-agnostic and
cheap to re-run, and its outputs must be directly comparable across
definitions. Adding a candidate surrogate should be a command-line argument.

**[DECIDE] resolved (§3): merge tree, not Reeb graph.** Cells are connected
components of a super-level set, and their appearance, growth and merging is
exactly a merge-tree story. Two consequences:

- **Persistence is the robustness measure that was asked for.** Rather than
  choosing a threshold, compute all thresholds at once; each cell then carries
  a birth value, a merge value, and their difference -- how far the knob must
  move to make it vanish. The reportable quantity is `N(p)`, the number of
  cells with persistence above `p`. A plateau means a defensible scale; a
  featureless continuum means no threshold is defensible, which is itself the
  finding.
- **Periodicity costs nothing.** A join tree is always a tree, on a torus as
  much as anywhere: components merge monotonically and cannot form a cycle.
  Loop surgery is needed for the contour tree and the Reeb graph, neither of
  which is required here. The doubly and triply periodic cases, with
  arbitrary boundary planes, come free -- periodicity enters only as adjacency
  in the union-find, which `AMRComplex` already supplies.

**Two additions the merge tree does not cover.**

- The crescent whose tips join to leave one contour inside another is a *b1*
  event: the region stays one connected component and the merge tree is blind
  to it. It needs per-component boundary invariants -- count the boundary
  curves of each labelled region in 2D, Euler characteristic per surface
  component in 3D. Small, and it sits on the segmentation.
- Euler characteristic is the wrong invariant in 2D, where every closed curve
  has `chi = 0`. The 2D statement is the number of boundary components and
  their homology class on the torus; a curve that wraps has no inside, so
  "contained within" is conditional on being contractible.

### 15.3 One engine serves both comparisons

Matching segmentations between consecutive snapshots by spatial overlap is the
same operation as matching a fuel-consumption segmentation against a
temperature segmentation on the same snapshot. Built once, it gives both:
whether a cell persists in time, and whether two surrogates identify the same
cells. The second measures the open question directly -- where surrogates
agree the structure is definition-independent, and where they diverge that is
differential diffusion expressed as topology.

Segmentation as a plotfile is what makes this work in practice: an integer
label per cell of the hierarchy feeds straight into `conditionalMean`, so
statistics of any field inside cells defined by any other field is a
one-liner. That is the common context the surrogates need to be compared in.

**A hypothesis worth testing, not asserting:** a merge-tree branch that never
becomes an independent maximum -- a saddle-separated shoulder of a larger
structure -- is structurally distinct from one that is its own high-persistence
branch. If back-supported structures are the former and self-sustaining ones
the latter, persistence measures the degree of independence. Physics judgement
required.

### 15.4 Build order, and why the sweep comes first

`N(v)` and `A(v)` -- the count of super-level-set components and the isosurface
area as functions of threshold -- can be had now, by looping the existing
contour and adding a union-find over the complex. That is the smallest thing
that produces a scientific result, it runs on real snapshots before either the
merge tree or the tracker exists, and it answers whether plateaus or scale
gaps are present before anything is designed around the assumption that they
are.

It is also a coarse sampling of the merge tree, so nothing is wasted: **when
the merge tree is built it must reproduce the sweep's `N(v)` exactly at every
sampled threshold.** That is a free and sharp regression gate, of the kind
that has caught four real defects in this work already.

Staging:

1. **Sweep** -- region components and sizes, surface area, per threshold.
   Distributed union-find over the complex; component labels are the minimum
   gid in the component, hence rank-independent by construction.
2. **Segmentation output** as a plotfile, and contour component counts.
3. **Merge tree**, gated on reproducing the sweep.
4. **Overlap matching**, serving time tracking and cross-surrogate comparison.
5. Boundary invariants for the b1 events.

**Open [DECIDE]:** whether a cell is a region thresholded from above (high
reaction rate, high temperature) or from below. The sweep supports both, and
the answer may differ per surrogate, so it stays a runtime choice rather than
a design commitment.

---

## 16. Sweep: built and verified

Stage 1 of §15.4 is done. `Src/AMRSegment.{H,cpp}` and `Src/amrSweep.cpp`, on
branch `amr-complex`; all suites pass in 2D and 3D at 1-8 ranks.

    amrSweep3d.<comp>.MPI.ex infile=<plt> isoCompName=<field> \
        vlo=<lo> vhi=<hi> nval=<n> [log=1] [above=0] [min_volume=<v>] \
        is_per="1 1 0" outfile=<table>

Per threshold: number of connected regions, total volume, largest and median
component volume, isosurface area.

### 16.1 Verified against closed forms

`makeTestPlt` gained a cellular field, `cos(2 pi k x/L)` summed over
directions, which on a periodic domain has exactly `k^DIM` maxima and saddles
between neighbours at a value that can be written down:

| case | above the saddle | below | analytic saddle |
|---|---|---|---|
| 3D, k=3 | 27 components | 1 | exactly 1 |
| 2D, k=3 | 9 components | 1 | exactly 0 |

The maxima sit *on* the domain faces, so finding exactly `k^DIM` is also a
test that the wrap is in the adjacency rather than assumed. The sphere case
checks volume and area against closed forms to 0.06% and 0.04%. Counts and
volumes are unchanged across `max_grid_size` 8/16/32/64 and across a
three-level hierarchy.

At a threshold sitting exactly on a saddle value the count is genuinely
ill-defined -- the 3D run reads 20 at `v = 1.000`. That is the degeneracy
persistence exists to quantify, not a defect.

### 16.2 Two things the tests forced

**Volumes are integer-derived.** Cell volumes take one value per level, so the
tally counts integers per (component, level) and multiplies at the end. Volume
is therefore an exact function of integers and bit-identical at any rank
count. Accumulating reals would have made the primary scientific output move
in its last bits with the decomposition -- comparable only to a tolerance, and
avoidable.

**Only shared nodes are exchanged.** The first version published a
`(gid, root)` pair for every participating node: `O(N)` all-to-all, which
passes at four ranks on a 64^3 box and fails on real snapshots. Now only
ghosts, owned nodes on a box boundary, and owned nodes adjacent to a COVERED
cell -- the last because they may be stamped as coarse parents into a finer
box's ghost region, and omitting them would fail to join components across a
coarse-fine interface, silently and only under MPI.

### 16.3 A process note worth keeping

`#ifdef AMREX_USE_MPI` regions are never compiled by a serial build, and no
MPI is available in the environment used to develop this. A most-vexing-parse
sat undetected in a gather for exactly that reason. The response was
structural: `AMRParallel` reduces the conditional region to the MPI calls
themselves, with every ordinary statement outside the guard and compiled by
every build. **Guarded regions are untested regions -- keep them minimal.**

### 16.4 Next

Representative snapshots. One file carrying two or three candidate surrogates
-- heat release rate, temperature, an intermediate -- swept individually, and
the `N(v)` curves laid side by side. That single run answers whether a plateau
exists, gives real timings for `ownedBaseRegions` and the segmentation maps at
production size, and shows whether the surrogates agree, before the merge tree
is built on an assumption about any of the three.

Real regridded hierarchies will exercise the proper-nesting assertion in a way
synthetic ones cannot. If it aborts, that is the check working.

---

## 17. First real snapshot: three defects, and the scaling wall

`~/DATA/plt01800`: a 3D flame sheet from PeleLMeX, periodic in x and y, walls
in z. Four levels, `ref_ratio` 2, about 708M cells (411M on level 3), 180 GB
on disk. Every refined level spans the full x-y period and reaches down to
the z = 0 inflow wall. Blocking factor 8, boxes of 16 and 32.

Every synthetic case refined one centred box, which never reaches a periodic
seam or a wall. The real hierarchy does both at every level, and the volume
gate failed at two levels. Three defects, each reproduced synthetically and
closed exactly; commit `01b398a`.

### 17.1 A stitch element owned twice at a periodic seam

On the low side of a refined region a stitch element's base is a ghost cell,
so no box holds it as valid. `ownedBaseRegions` assigns it with two tests:
drop it if it is valid through the wrap, else give it to the lowest-indexed
box that visits it. The first test considered periodic images; the second
compared raw indices. At a seam the same base is `(-1,j)` to the box at the
low end and `(n-1,j)` to the box at the high end, so each owned it under its
own name.

On a 16x16 2D slab the overshoot is exactly one element, area `9/4096`. On
the snapshot (levels 0-1, all periodic) it was 4048 unmatched faces, every
one on the x or y seam, and 1122 duplicated elements.

Fix: ownership is decided on one name per base cell, its image inside the
domain, and a claim is checked against every periodic image of the other
boxes' regions. The owned region is handed back in the box's own frame,
because that is where its ghost cells -- and so the element's corners -- are
read. **Convention 1 now applies to base cells as well as nodes.** Nodes had
been canonicalised under the wrap; base cells had not, so one element could
carry two names and ownership was decided per name.

The first diagnosis was wrong and is worth recording: 7488 elements with 7
distinct corners suggested re-entrant corners of the refined region. An
L-shaped region closes exactly in 2D and 3D and contains 7-corner elements
too. They are a legitimate stitch shape at a concave edge.

### 17.2 The nesting assertion exempted ghosts across a periodic face

It examined only ghosts inside the domain. A level that reaches a periodic
face with no coarser level on the far side leaves its wrapped ghosts INIT,
and their elements are dropped: a gap one stitch layer thick, silently. The
snapshot read as periodic in z showed exactly this -- level 2 sits on z = 0
and wraps onto the top of the domain, where there is no level 1 -- as a
residual of exactly `5/3072`, the layer between a level-2 and a level-0 node
row across the full x-y area. Any ghost inside the domain or across a
periodic face that finds no parent now aborts. That hierarchy is not nested
for periodic z, which it never was physically, and the check now says so.

### 17.3 The volume reference assumed level 0 at every wall

`expectedVolume` subtracted half a level-0 cell at each wall. With a refined
level on a wall the complex stops half a *fine* cell from it, and slopes
between levels where a refined patch on the wall ends. The leading term
accounted for 99% of the snapshot's residual; the slope is the rest.

The reference is now measured, by the divergence theorem with `F = x_w e_w`
for a non-periodic direction `w`: the volume of the complex is the flux of
`F` through its boundary, and a correct complex has no boundary except faces
lying against a wall. So the reference is that flux summed over exactly those
faces of every element, each split into the same triangles the volume's
simplex split puts on it. Faces away from the walls must cancel in pairs, so
a gap, an overlap, a fold or a face-diagonal disagreement still leaves a
residual. Fully periodic, the reference is still the domain.

- Reproduces `(15/16)^DIM` on one level with walls, the old closed form, and
  `1 - 1/256 - 1/64` on a wall-anchored sheet: half a level-2 cell at one
  wall, half a level-0 cell at the other.
- **Mutation test.** With 17.1 reinstated in a scratch copy, the walled cases
  fail by 2.3e-3 (2D) and 5.1e-3 (3D). The first mutant tried did not
  reproduce the defect at all -- clipping to the domain without the image
  shifts drops the out-of-domain name outright -- and its PASS meant nothing.
  A reference measured from the mesh must be shown able to fail.
- **Blind spot, documented at the code:** an element whose every corner lies
  against one wall has all its faces counted as wall faces, so its absence
  moves both sums equally. The face census still sees it.

### 17.4 New synthetic geometry and controls

`makeTestPlt region=`:

| region | shape | exercises |
|---|---|---|
| `L` | finest level with the high quadrant removed | re-entrant corner (2D), edge (3D) |
| `slab` | spans all but the last direction | coarse-fine interfaces meeting the periodic seams |
| `sheet` | slab on the low wall of the last direction | a refined level on a wall |

`mpi_RE` and `mpi_SL` run in every section of the MPI suite, including the
fully walled setting. `mpi_SH` has its own section: rank invariance under its
physical periodicities, and a negative control that must abort on nesting
when the wall direction is forced periodic.

### 17.5 Verified on the snapshot

x and y periodic, z walls, 64 ranks:

| levels | elements | relative error | verify wall time |
|---|---|---|---|
| 0-1 | 47,871,576 | 0 | 3.9 s |
| 0-2 | 262,139,256 | 0 | 12.5 s |

Both without a carried field. Reading one field adds I/O: the levels 0-2
sweep spent 189 s building, almost all of it reading `HeatRelease` from the
180 GB plotfile.

### 17.6 The scaling wall

A levels 0-2 sweep on `HeatRelease` peaked at **89 GB total** after build and
verify at 64 ranks, about 300 bytes per cell. The full hierarchy is 708M cells
on a 122 GB machine, so it does not fit. Contributors, by construction:

- `m_elements`, an explicit `Element` per cell: 64 bytes.
- `Node` (32 bytes) and the field (8 bytes per component) per cell, both with
  a ghost layer.
- `verify()`'s face census, a `std::map` over every face instance on the rank.
- In `AMRSegment`, hash maps and sets per in-set node, and `m_ownedSet` over
  every owned node.

And the merge: `mergeAcrossRanks` all-gathers every shared node's
`(gid, root)` to every rank. On level 0 that was 24% of all nodes at the
lowest threshold; with level-3 boxes averaging about 24^3 it is a similar
fraction of hundreds of millions, replicated on every rank.

### 17.7 Next

1. **Memory.** Drop the stored element list; make the face census optional
   outside tests; segmentation containers sized to the in-set nodes only.
2. **Merge.** Reduce to edges between local components before the global
   exchange: their number is set by component pieces per box, not by nodes.
3. **The four-marker sweep** on the full hierarchy: `I_R(H2)`, `temp`,
   `HeatRelease`, `Y(H2O2)`.
