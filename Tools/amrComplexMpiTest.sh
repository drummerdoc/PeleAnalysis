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

shopt -s nullglob
VERIFY=("$SRC"/amrComplexVerify3d.*MPI.ex)
MAKEPLT=("$SRC"/makeTestPlt3d.*MPI.ex)
shopt -u nullglob

if [ ${#VERIFY[@]} -eq 0 ] || [ ${#MAKEPLT[@]} -eq 0 ]; then
    echo "ERROR: no MPI executables found in $SRC."
    echo "Build them with USE_MPI=TRUE (see the header of this script)."
    exit 1
fi
VERIFY=${VERIFY[0]}
MAKEPLT=${MAKEPLT[0]}

WORK=${WORK:-$(pwd)/amrComplex_mpi_test}
mkdir -p "$WORK" && cd "$WORK" || exit 1

echo "verify : $VERIFY"
echo "workdir: $WORK"
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
echo

CASES="mpi_L1 mpi_L2 mpi_L3 mpi_R4"

# Parallel arrays, not a single string: each periodicity setting contains
# spaces, so a space-separated list would be split into nonsense.
PER_TAG=(per111 per000 per110)
PER_VAL=("1 1 1" "0 0 0" "1 1 0")

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
done
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
  $MPIRUN -n $NP "$VERIFY" infile=$c is_per="1 1 1" varnames=dist \
      DistributionMapping.strategy=SFC 2>/dev/null | invariants > "str_${c}_sfc.txt"
  for strat in KNAPSACK ROUNDROBIN; do
    $MPIRUN -n $NP "$VERIFY" infile=$c is_per="1 1 1" varnames=dist \
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
  out=$($MPIRUN -n 4 "$VERIFY" infile=mpi_mgs$mgs is_per="1 1 1" varnames=dist 2>/dev/null)
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
