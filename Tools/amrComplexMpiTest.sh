#!/usr/bin/env bash
#
# MPI verification for AMRComplex.
#
# What this is testing, and why these tests and not others.
#
# AMRComplex claims that neither the numbering nor the element ownership
# depends on how the data is distributed. Both claims are structural rather
# than empirical -- a BoxArray is globally replicated, so its index order is
# the same on every rank, and gid is a pure function of (level, IntVect) --
# but they have never been run under MPI. These tests are what turns the
# argument into a measurement.
#
# Comparing element *counts* across rank counts would be weak: two different
# meshes can have the same count. So the exe emits two order-independent
# checksums -- a sum of per-item hashes, invariant to visit order and to who
# holds what, but sensitive to any difference in membership, connectivity, or
# the exact bits of the propagated field values. Those are the real invariant.
# The node checksum in particular is the only thing that can see a
# rank-dependent difference in the coarse-fine field propagation.
#
# Varying the rank count alone is a weaker test than it looks, because it
# conflates two things: how many ranks there are, and which box lands where.
# So the second section holds the rank count fixed and varies the
# DistributionMapping strategy instead, which changes the assignment without
# changing anything else.
#
# ---------------------------------------------------------------------------
# Running it
#
# Build both executables with MPI. AMREX_HOME defaults to ../../amrex relative
# to Src, which is already right for a tree laid out as src/PeleAnalysis and
# src/amrex side by side; COMP must match whatever you normally build with.
#
#   cd Src
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=makeTestPlt -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=amrContourVerify -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=amrSweep -j8
#
# Then, from any writable directory:
#
#   /path/to/PeleAnalysis/Tools/amrComplexMpiTest.sh
#
# Everything is written under ./amrComplex_mpi_test; set WORK to put it
# elsewhere. The script exits nonzero with the number of failures, so it can
# be used as a gate.
#
#   Tools/amrComplexMpiTest.sh 1 2 4 6      # choose the rank counts
#   DIM=2 Tools/amrComplexMpiTest.sh        # test the 2D build instead
#   WORK=/scratch/acx Tools/amrComplexMpiTest.sh
#
# Launcher. The script invokes "$MPIRUN -n <N> <exe>", which suits mpiexec,
# mpirun and srun alike:
#
#   OpenMPI on a laptop  MPIRUN="mpirun --oversubscribe" Tools/amrComplexMpiTest.sh
#                        -- required whenever a requested rank count exceeds
#                           the physical core count, which 7 and 8 easily do.
#                           MPICH and Intel MPI do not need it.
#   Slurm                MPIRUN=srun Tools/amrComplexMpiTest.sh
#                        -- from inside an salloc or a batch script sized for
#                           the largest rank count requested.
#
# The cases are small on purpose: the largest is a three-level 64^3 hierarchy,
# about 740k elements, a few seconds per run. This is a correctness gate, not
# a scaling study, and it is meant to be cheap enough to run on every change.
# ---------------------------------------------------------------------------
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
RANKS=${*:-"1 2 3 4 7 8"}
MPIRUN=${MPIRUN:-mpiexec}

# Both dimensions are supported. DIM selects which set of executables to test;
# the periodicity settings and the analytic reference follow from it.
DIM=${DIM:-3}
DTAG="${DIM}d"

shopt -s nullglob
VERIFY=("$SRC"/amrComplexVerify${DTAG}.*MPI.ex)
CONTOUR=("$SRC"/amrContourVerify${DTAG}.*MPI.ex)
SWEEP=("$SRC"/amrSweep${DTAG}.*MPI.ex)
MAKEPLT=("$SRC"/makeTestPlt${DTAG}.*MPI.ex)
shopt -u nullglob

if [ ${#VERIFY[@]} -eq 0 ] || [ ${#MAKEPLT[@]} -eq 0 ]; then
    echo "ERROR: no MPI executables found in $SRC."
    echo "Build them with USE_MPI=TRUE (see the header of this script)."
    exit 1
fi
VERIFY=${VERIFY[0]}
MAKEPLT=${MAKEPLT[0]}
CONTOUR=${CONTOUR[0]:-}
SWEEP=${SWEEP[0]:-}

WORK=${WORK:-$(pwd)/amrComplex_mpi_test}
mkdir -p "$WORK" && cd "$WORK" || exit 1

echo "verify : $VERIFY"
echo "workdir: $WORK"
echo "dim    : $DIM"
echo "ranks  : $RANKS"
echo

# The invariant block.
#
# verify() marks with a leading "~" every line whose last bits legitimately
# depend on how the work was partitioned, and those are dropped here:
#
#   - the volume sum is compensated per rank, but the cross-rank reduction is
#     a plain sum, so the residual moves by an ulp or so with the rank count.
#     Real, benign, and not something to compare for equality -- the tolerance
#     test in the RESULT line is the comparable form of that statement.
#   - the face incidence census is irreducibly rank-local: each rank sees only
#     its own elements, so a face shared across a rank boundary is counted
#     once on each side and looks like two boundary faces. The rank-independent
#     face invariant is "checksum faces", an XOR over every face instance,
#     which cancels shared faces wherever their elements live.
#
# The census is now off by default in the exe, because on a 333M-element
# hierarchy it is an ordered-map pass over 2e9 faces whose destructor alone
# overflows the stack. Every case here is under a million elements, where it
# costs nothing, so the suite asks for it explicitly and keeps the coverage.
#
# What remains must be byte-identical.
invariants () {
    sed -n '/=== AMRComplex verify ===/,/RESULT/p' | grep -v '^$' | grep -v '^~'
}

FAIL=0
pass () { echo "  PASS  $1"; }
fail () { echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

# --------------------------------------------------------------------------
# Test data. Generated once, at one rank, so the plotfile's own BoxArray is
# fixed for the whole run: a different BoxArray would change the gids
# legitimately and there would be nothing to compare.
# --------------------------------------------------------------------------
echo "generating test hierarchies..."
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=1              outfile=mpi_L1 > /dev/null 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=2              outfile=mpi_L2 > /dev/null 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=64 nlevels=3              outfile=mpi_L3 > /dev/null 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=2 ref_ratio=4  outfile=mpi_R4 > /dev/null 2>&1
# Geometry a centred box never produces. mpi_RE: the finest level is L-shaped,
# so it has a re-entrant corner (2D) or edge (3D). mpi_SL: every refined level
# spans the domain in all but the last direction -- the flame-sheet layout --
# so coarse-fine interfaces meet the periodic seams. The first real snapshot
# failed on exactly this: a stitch base at the seam had two raw names and was
# owned, and generated, once under each.
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=2 region=L     outfile=mpi_RE > /dev/null 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=3 region=slab  outfile=mpi_SL > /dev/null 2>&1
# mpi_SH: a slab sitting on the low wall of the last direction, as a flame
# sheet refined down to its inflow is. Every refined level reaches that wall,
# so the complex stops half a FINE cell from it -- which the old level-0
# reference got wrong. Not in CASES: with that direction periodic it is not
# properly nested (the finest level would wrap onto the coarse far side), so
# the all-periodic sections cannot run it. Section 1b covers it instead.
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=3 region=sheet outfile=mpi_SH > /dev/null 2>&1
echo

CASES="mpi_L1 mpi_L2 mpi_L3 mpi_R4 mpi_RE mpi_SL"

# Parallel arrays, not a single string: each periodicity setting contains
# spaces, so a space-separated list would be split into nonsense.
if [ "$DIM" = "2" ]; then
    PER_TAG=(per11 per00 per10)
    PER_VAL=("1 1" "0 0" "1 0")
    ALLPER="1 1"
else
    PER_TAG=(per111 per000 per110)
    PER_VAL=("1 1 1" "0 0 0" "1 1 0")
    ALLPER="1 1 1"
fi

# --------------------------------------------------------------------------
# 1. Rank invariance.
#
# Every reported quantity must be byte-identical to the single-rank run:
# element and node counts, volume residual, shape histogram, face census,
# unfilled-node count, and both checksums. Not "close" -- identical. The
# numbering is rank-independent by construction, so there is no rounding
# story that would excuse a difference.
# --------------------------------------------------------------------------
echo "=== 1. rank invariance ==============================================="
for c in $CASES; do
  for idx in "${!PER_TAG[@]}"; do
    tag=${PER_TAG[$idx]}; per=${PER_VAL[$idx]}
    $MPIRUN -n 1 "$VERIFY" infile=$c is_per="$per" varnames=dist face_census=1 2>/dev/null \
        | invariants > "ref_${c}_${tag}.txt"
    if ! grep -q "PASS" "ref_${c}_${tag}.txt"; then
        fail "$c per=$per : the 1-rank baseline itself does not pass"
        continue
    fi
    for n in $RANKS; do
      [ "$n" = "1" ] && continue
      $MPIRUN -n "$n" "$VERIFY" infile=$c is_per="$per" varnames=dist face_census=1 2>/dev/null \
          | invariants > "out_${c}_${tag}_n${n}.txt"
      if diff -q "ref_${c}_${tag}.txt" "out_${c}_${tag}_n${n}.txt" > /dev/null; then
          pass "$c per=$per  n=$n"
      else
          fail "$c per=$per  n=$n"
          diff "ref_${c}_${tag}.txt" "out_${c}_${tag}_n${n}.txt" | sed 's/^/        /'
      fi
    done
  done
done
echo

# --------------------------------------------------------------------------
# 1b. A refined level on a wall.
#
# Rank invariance under the periodicities the sheet actually has -- spanned
# directions periodic or not, wall direction a wall -- where passing at all
# requires the measured wall reference. And the negative control: forcing the
# wall direction periodic must ABORT on nesting. Before the nesting check
# looked across periodic faces, that case quietly dropped a whole layer of
# seam elements instead.
# --------------------------------------------------------------------------
echo "=== 1b. refined level on a wall ======================================"
c=mpi_SH
for idx in 1 2; do
  tag=${PER_TAG[$idx]}; per=${PER_VAL[$idx]}
  $MPIRUN -n 1 "$VERIFY" infile=$c is_per="$per" varnames=dist 2>/dev/null \
      | invariants > "ref_${c}_${tag}.txt"
  if ! grep -q "PASS" "ref_${c}_${tag}.txt"; then
      fail "$c per=$per : the 1-rank baseline itself does not pass"
      continue
  fi
  for n in $RANKS; do
    [ "$n" = "1" ] && continue
    $MPIRUN -n "$n" "$VERIFY" infile=$c is_per="$per" varnames=dist 2>/dev/null \
        | invariants > "out_${c}_${tag}_n${n}.txt"
    if diff -q "ref_${c}_${tag}.txt" "out_${c}_${tag}_n${n}.txt" > /dev/null; then
        pass "$c per=$per  n=$n"
    else
        fail "$c per=$per  n=$n"
        diff "ref_${c}_${tag}.txt" "out_${c}_${tag}_n${n}.txt" | sed 's/^/        /'
    fi
  done
done
if $MPIRUN -n 1 "$VERIFY" infile=$c is_per="$ALLPER" varnames=dist 2>&1 \
       | grep -q "not properly nested"; then
    pass "$c per=$ALLPER  aborts on nesting (negative control)"
else
    fail "$c per=$ALLPER  should abort on nesting across the periodic face"
fi
echo

# --------------------------------------------------------------------------
# 2. Distribution invariance at fixed rank count.
#
# Sharper than varying the rank count, because it changes only which box
# lands on which rank. A bug that happens to be insensitive to rank count --
# say, one that depends on a box's neighbours being local -- shows up here.
# --------------------------------------------------------------------------
echo "=== 2. distribution invariance (4 ranks, varying strategy) ==========="
NP=4
for c in $CASES; do
  $MPIRUN -n $NP "$VERIFY" infile=$c is_per="$ALLPER" varnames=dist face_census=1 \
      DistributionMapping.strategy=SFC 2>/dev/null | invariants > "str_${c}_sfc.txt"
  for strat in KNAPSACK ROUNDROBIN; do
    $MPIRUN -n $NP "$VERIFY" infile=$c is_per="$ALLPER" varnames=dist face_census=1 \
        DistributionMapping.strategy=$strat 2>/dev/null | invariants > "str_${c}_${strat}.txt"
    if diff -q "str_${c}_sfc.txt" "str_${c}_${strat}.txt" > /dev/null; then
        pass "$c  SFC vs $strat"
    else
        fail "$c  SFC vs $strat"
        diff "str_${c}_sfc.txt" "str_${c}_${strat}.txt" | sed 's/^/        /'
    fi
  done
done
echo

# --------------------------------------------------------------------------
# 3. Decomposition independence of the answer, not of the numbering.
#
# A different max_grid_size gives a different BoxArray, so the gids genuinely
# differ and the checksums are *expected* to differ -- comparing them would be
# wrong. What must not change is the geometry: the complex still tiles the
# domain exactly, is still face-conforming, and still has no unfilled nodes.
# --------------------------------------------------------------------------
echo "=== 3. decomposition independence (varying max_grid_size) ============"
for mgs in 8 16 32; do
  $MPIRUN -n 1 "$MAKEPLT" n_cell=64 nlevels=3 max_grid_size=$mgs \
      outfile=mpi_mgs$mgs > /dev/null 2>&1
  out=$($MPIRUN -n 4 "$VERIFY" infile=mpi_mgs$mgs is_per="$ALLPER" varnames=dist face_census=1 2>/dev/null)
  res=$( echo "$out" | awk '/RESULT/{print $2}')
  closed=$(echo "$out" | awk '/closed complex/{print $4}')
  unf=$(  echo "$out" | awk '/component\(s\)/{print $5}')
  # The rank-local census is not usable here; "closed complex" is the global
  # statement that every face is shared by exactly two elements.
  if [ "$res" = "PASS" ] && [ "$closed" = "yes" ] && [ "$unf" = "0" ]; then
      pass "max_grid_size=$mgs  (volume=$res closed=$closed unfilled=$unf)"
  else
      fail "max_grid_size=$mgs  (volume=$res closed=$closed unfilled=$unf)"
  fi
done
echo

# --------------------------------------------------------------------------
# 4. The contour, same treatment.
#
# The surface is where floating point finally enters: vertices are
# interpolated along edges, so nothing here is bitwise-exact the way the
# complex is. The invariants are the integer ones -- vertex and triangle
# counts, and the edge parity, an XOR over every triangle-edge instance that
# is zero exactly when the surface is closed.
#
# Area is deliberately not compared for equality, for the same reason the
# volume is not.
# --------------------------------------------------------------------------
if [ -n "$CONTOUR" ]; then
echo "=== 4. contour invariance ============================================"
# End the range on the "closed" line, not on a longer phrase inside it: this
# extractor was silently broken once already when that line was reworded, and
# an unterminated sed range runs to end of file and sweeps in whatever follows.
cinv () { sed -n '/=== AMRContour verify ===/,/^  closed /p' | grep -v '^$' | grep -v '^~'; }
for c in $CASES; do
  $MPIRUN -n 1 "$CONTOUR" infile=$c is_per="$ALLPER" isoVal=0 2>/dev/null | cinv > "ct_${c}_ref.txt"
  for n in $RANKS; do
    [ "$n" = "1" ] && continue
    $MPIRUN -n "$n" "$CONTOUR" infile=$c is_per="$ALLPER" isoVal=0 2>/dev/null | cinv > "ct_${c}_n${n}.txt"
    if diff -q "ct_${c}_ref.txt" "ct_${c}_n${n}.txt" > /dev/null; then
        pass "contour $c  n=$n"
    else
        fail "contour $c  n=$n"
        diff "ct_${c}_ref.txt" "ct_${c}_n${n}.txt" | sed 's/^/        /'
    fi
  done
done
echo
else
echo "=== 4. contour invariance: SKIPPED (no MPI amrContourVerify built) ==="
echo
fi

# --------------------------------------------------------------------------
# 5. The threshold sweep.
#
# The scientific output, so the bar is the same as everywhere else: identical,
# not close. Component counts are integers, and component volumes are formed
# from integer counts per level times that level's cell volume, so they are an
# exact function of integers rather than a sum of reals in partition order.
# Nothing here has any business changing with the rank count.
#
# Run with area=0: surface area is a genuine floating-point sum and is already
# covered, in its own terms, by section 4.
# --------------------------------------------------------------------------
if [ -n "$SWEEP" ]; then
echo "=== 5. sweep invariance =============================================="
# Columns 1-7 only: threshold through area. The exchanged-pair count and the
# timings after them depend on the decomposition by design.
sinv () { grep -E '^ ' | grep -v '^#' | awk '{print $1,$2,$3,$4,$5,$6,$7}'; }
for c in $CASES; do
  for fld in cells dist; do
    if [ "$fld" = "cells" ]; then LO=2.5; HI=-1.0; else LO=-0.2; HI=0.15; fi
    $MPIRUN -n 1 "$SWEEP" infile=$c isoCompName=$fld is_per="$ALLPER" \
        vlo=$LO vhi=$HI nval=12 area=0 2>/dev/null | sinv > "sw_${c}_${fld}_ref.txt"
    if [ ! -s "sw_${c}_${fld}_ref.txt" ]; then
        fail "sweep $c $fld : 1-rank run produced no table"; continue
    fi
    for n in $RANKS; do
      [ "$n" = "1" ] && continue
      $MPIRUN -n "$n" "$SWEEP" infile=$c isoCompName=$fld is_per="$ALLPER" \
          vlo=$LO vhi=$HI nval=12 area=0 2>/dev/null | sinv > "sw_${c}_${fld}_n${n}.txt"
      if diff -q "sw_${c}_${fld}_ref.txt" "sw_${c}_${fld}_n${n}.txt" > /dev/null; then
          pass "sweep $c $fld  n=$n"
      else
          fail "sweep $c $fld  n=$n"
          diff "sw_${c}_${fld}_ref.txt" "sw_${c}_${fld}_n${n}.txt" | sed 's/^/        /' | head -8
      fi
    done
  done
done
echo
else
echo "=== 5. sweep invariance: SKIPPED (no MPI amrSweep built) ============="
echo
fi

# --------------------------------------------------------------------------
echo "======================================================================"
if [ $FAIL -eq 0 ]; then
    echo "ALL PASS"
else
    echo "$FAIL FAILURE(S) -- see the diffs above."
    echo
    echo "Reading a failure:"
    echo "  checksum elts differs, counts agree  -> ownership: some element is"
    echo "     produced by a different rank, twice, or not at all. Look at"
    echo "     ownedBaseRegions and the lower-index claim."
    echo "  checksum nodes differs, elts agree   -> field propagation: the"
    echo "     coarse-fine stamp or the periodic FillBoundary is giving"
    echo "     different values depending on locality. Identity is fine;"
    echo "     values are not."
    echo "  face census differs                  -> connectivity differs across"
    echo "     ranks, which should be impossible if gid is rank-independent."
    echo "     Suspect the numbering before anything else."
    echo "  volume differs in the last digits    -> compensated summation is"
    echo "     per-rank; the cross-rank reduction is a plain sum. Only"
    echo "     plausible at very large rank counts, and it is the one"
    echo "     difference here that would be benign."
fi
exit $FAIL
