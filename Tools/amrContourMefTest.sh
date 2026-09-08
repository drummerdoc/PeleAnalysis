#!/usr/bin/env bash
#
# MEF write test under MPI.
#
# writeMEF is the one path in AMRContour that gathers to the IO processor, and
# gathering is where every MPI defect in this work has lived so far. The claim
# it makes is unusually strong and therefore unusually easy to falsify:
# duplicate surface vertices are merged by their node-edge KEY rather than by
# comparing coordinates, and the global node order is a sort on those keys, so
# the output file should be BYTE-IDENTICAL at any rank count. isosurface.cpp's
# MEF files never were, because their node order came from the gather order.
#
# So the test is just cmp. No tolerance, no field-by-field comparison: either
# the bytes match or the claim is wrong.
#
# It also round-trips each file through surfMEFtoDAT, because a file can be
# self-consistent and still not be readable by the toolchain that matters.
#
# Build first (COMP to match your usual):
#
#   cd Src
#   for e in makeTestPlt amrContourVerify surfMEFtoDAT; do
#     make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=$e -j8
#   done
#
# Then:
#
#   Tools/amrContourMefTest.sh                 # 3D, ranks 1 2 3 4
#   DIM=2 Tools/amrContourMefTest.sh           # 2D
#   MPIRUN="mpirun --oversubscribe" Tools/amrContourMefTest.sh
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
DIM=${DIM:-3}
DTAG="${DIM}d"
RANKS=${*:-"1 2 3 4"}
MPIRUN=${MPIRUN:-mpiexec}

if [ "$DIM" = "2" ]; then PER="1 1"; VARS="xco yco"; else PER="1 1 1"; VARS="xco yco zco"; fi

shopt -s nullglob
CONTOUR=("$SRC"/amrContourVerify${DTAG}.*MPI.ex)
MAKEPLT=("$SRC"/makeTestPlt${DTAG}.*MPI.ex)
TODAT=("$SRC"/surfMEFtoDAT${DTAG}.*MPI.ex)
shopt -u nullglob

if [ ${#CONTOUR[@]} -eq 0 ] || [ ${#MAKEPLT[@]} -eq 0 ]; then
    echo "ERROR: no MPI executables for DIM=$DIM in $SRC (see header)."
    exit 1
fi
CONTOUR=${CONTOUR[0]}; MAKEPLT=${MAKEPLT[0]}; TODAT=${TODAT[0]:-}

WORK=${WORK:-$(pwd)/amrContour_mef_test}
mkdir -p "$WORK" && cd "$WORK" || exit 1

echo "contour: $CONTOUR"
echo "dim    : $DIM      ranks: $RANKS"
echo "workdir: $WORK"
echo

FAIL=0
pass () { echo "  PASS  $1"; }
fail () { echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

# Generated once at one rank so the plotfile's BoxArray is fixed: a different
# BoxArray would change the gids legitimately and there would be nothing to
# compare.
echo "generating test hierarchies..."
$MPIRUN -n 1 "$MAKEPLT" n_cell=64 nlevels=1 outfile=mef_L1 > /dev/null 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=64 nlevels=3 outfile=mef_L3 > /dev/null 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=32 nlevels=2 ref_ratio=4 outfile=mef_R4 > /dev/null 2>&1
echo

echo "=== 1. MEF files are byte-identical across rank counts ==============="
for c in mef_L1 mef_L3 mef_R4; do
  for n in $RANKS; do
    $MPIRUN -n "$n" "$CONTOUR" infile=$c is_per="$PER" isoVal=0 \
        varnames="$VARS" check_linear_components=1 \
        outfile=${c}_n${n}.mef > ${c}_n${n}.out 2>&1
    if [ ! -s "${c}_n${n}.mef" ]; then fail "$c n=$n produced no file"; continue; fi
    # The carried-component identity must hold at every rank count too.
    if ! grep -q "exactness.*PASS" "${c}_n${n}.out"; then
        fail "$c n=$n carried-component exactness"
    fi
  done
  for n in $RANKS; do
    [ "$n" = "1" ] && continue
    if cmp -s "${c}_n1.mef" "${c}_n${n}.mef"; then
        pass "$c  1 vs $n ranks  (byte-identical)"
    else
        fail "$c  1 vs $n ranks"
        echo "        $(cmp "${c}_n1.mef" "${c}_n${n}.mef" 2>&1 | head -1)"
        echo "        sizes: $(stat -f%z "${c}_n1.mef" 2>/dev/null || stat -c%s "${c}_n1.mef") vs $(stat -f%z "${c}_n${n}.mef" 2>/dev/null || stat -c%s "${c}_n${n}.mef")"
    fi
  done
done
echo

echo "=== 2. every file is readable by surfMEFtoDAT ========================"
if [ -n "$TODAT" ]; then
  for c in mef_L1 mef_L3 mef_R4; do
    for n in $RANKS; do
      rm -f ${c}_n${n}.dat
      $MPIRUN -n 1 "$TODAT" infile=${c}_n${n}.mef > /dev/null 2>&1
      if [ -s "${c}_n${n}.dat" ]; then
          hdr=$(sed -n '2p' ${c}_n${n}.dat)
          case "$hdr" in
            *"N="*"E="*) pass "$c n=$n  -> ${hdr}" ;;
            *) fail "$c n=$n  unexpected DAT header: $hdr" ;;
          esac
      else
          fail "$c n=$n  surfMEFtoDAT produced nothing"
      fi
    done
  done
else
  echo "  SKIPPED (no MPI surfMEFtoDAT built)"
fi
echo

echo "======================================================================"
if [ $FAIL -eq 0 ]; then
    echo "ALL PASS"
else
    echo "$FAIL FAILURE(S)."
    echo
    echo "Reading a failure:"
    echo "  files differ, same size      -> the node ORDER differs, so the sort"
    echo "     on node-edge keys is not producing the same order everywhere."
    echo "     Suspect the key comparison, not the gather."
    echo "  files differ, sizes differ   -> the vertex MERGE differs: some"
    echo "     vertex is duplicated or dropped. Suspect the dedup map, and"
    echo "     check whether two ranks disagree about a node-edge key."
    echo "  exactness fails only at n>1  -> the carried components are not"
    echo "     surviving the gather, i.e. the per-vertex data is being packed"
    echo "     or unpacked with the wrong stride."
fi
exit $FAIL
