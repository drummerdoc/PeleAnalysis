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
# Build first, with MPI:
#
#   cd Src
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE AMREX_HOME=<path> -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE AMREX_HOME=<path> EBASE=makeTestPlt -j8
#
# Then:  Tools/amrComplexMpiTest.sh [rank counts...]     (default: 1 2 3 4 7 8)
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

# The invariant block: everything verify() reports except timing. Blank lines
# and the banner are dropped so the comparison is on content only.
invariants () { sed -n '/=== AMRComplex verify ===/,/RESULT/p' | grep -v '^$'; }

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
PERS="1|1 1 1 0|0 0 0 m|1 1 0"

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
  for pp in $PERS; do
    tag=${pp%%|*}; per=${pp#*|}
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
  relerr=$(echo "$out" | awk '/relative error/{print $4}')
  conf=$(  echo "$out" | awk '/face-conforming/{print $3}')
  unf=$(   echo "$out" | awk '/component\(s\)/{print $5}')
  bad=$(   echo "$out" | awk '/non-manifold/{print $4}')
  if [ "$relerr" = "0" ] && [ "$conf" = "yes" ] && [ "$unf" = "0" ] && [ "$bad" = "0" ]; then
      pass "max_grid_size=$mgs  (relerr=$relerr conforming=$conf unfilled=$unf nonmanifold=$bad)"
  else
      fail "max_grid_size=$mgs  (relerr=$relerr conforming=$conf unfilled=$unf nonmanifold=$bad)"
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
