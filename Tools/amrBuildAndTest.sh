#!/usr/bin/env bash
#
# Build every AMRComplex/AMRContour executable in both dimensions and run every
# test suite. One command, so there is no sequence to get wrong.
#
#   Tools/amrBuildAndTest.sh
#
# Environment:
#   COMP=llvm            compiler (default: llvm on macOS, gnu elsewhere)
#   MPIRUN=...           launcher; OpenMPI on a laptop needs
#                        MPIRUN="mpirun --oversubscribe"
#   RANKS="1 2 3 4"      rank counts for the MPI suites
#   DIMS="2 3"           dimensions to build and test
#   JOBS=8               make -j
#   SKIP_BUILD=1         run the tests against what is already built
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
# The default follows the platform. It used to be llvm unconditionally, which
# is right on the Mac this tree was developed on and wrong on a Linux machine,
# where it built llvm executables beside a stale gnu set and the suites that
# picked alphabetically then tested the stale ones. COMP is exported so every
# suite tests exactly what was built here, whatever else is in Src.
case $(uname -s) in
    Darwin) DEFCOMP=llvm ;;
    *)      DEFCOMP=gnu  ;;
esac
COMP=${COMP:-$DEFCOMP}
export COMP
JOBS=${JOBS:-8}
DIMS=${DIMS:-"2 3"}
RANKS=${RANKS:-"1 2 3 4"}
export MPIRUN=${MPIRUN:-mpiexec}
SKIP_BUILD=${SKIP_BUILD:-0}

# amrGradStats and amrGradVerify are last on purpose: they are the only targets
# that pull in Boundary, AmrCore and LinearSolvers, so they are also the only
# ones whose first build spends minutes compiling AMReX. Putting them last means
# a mistake in the cheap targets is reported in seconds rather than after that.
TARGETS="amrComplexVerify amrContourVerify amrSweep amrLabel amrJointFit makeTestPlt surfMEFtoDAT"
TARGETS="$TARGETS amrGradStats amrGradVerify amrCurvVerify"
LIBSRC="AMRComplex.cpp AMRContour.cpp AMRSegment.cpp AMRParallel.cpp AMRJointStats.cpp"
LIBSRC="$LIBSRC AMRGradient.cpp"
LOGDIR=${LOGDIR:-$(pwd)/amr_test_logs}
mkdir -p "$LOGDIR"

echo "=================================================================="
echo " build : COMP=$COMP  DIMS=$DIMS  jobs=$JOBS"
echo " tests : ranks=$RANKS  launcher=$MPIRUN"
echo " logs  : $LOGDIR"
echo "=================================================================="
echo

# ------------------------------------------------------------------------
# Step zero: the one build failure this tree keeps repeating.
#
# A driver whose name is the case-variant of a library source is the SAME FILE
# on a case-insensitive filesystem. Writing both leaves one file, one object is
# built instead of two, and the link fails on every symbol of the library with
# "referenced from _main" -- which looks like a makefile omission and is not.
# It has happened twice: amrComplex.cpp against AMRComplex.cpp, and
# amrJointStats.cpp against AMRJointStats.cpp. Both are cheap to rule out, so
# they are ruled out before anything is compiled.
#
#   a library source must not define main
#   a target must have a source that does
#   no target may be the case-variant of a library source
# ------------------------------------------------------------------------
# The return type may sit on its own line, as it does in this tree's older
# drivers, so the test is on the call signature rather than on "int main".
# The leading character class keeps it from matching domain( or _main.
has_main () { grep -qE '(^|[^_[:alnum:]])main[[:space:]]*\(' "$1"; }

printf "  checking source/driver names ... "
NAMEFAIL=""
for f in $LIBSRC; do
    if [ ! -f "$SRC/$f" ]; then
        NAMEFAIL="$NAMEFAIL\n    missing library source: $f"
    elif has_main "$SRC/$f"; then
        NAMEFAIL="$NAMEFAIL\n    $f defines main, so it is a driver, not a library source."
        NAMEFAIL="$NAMEFAIL\n      On a case-insensitive filesystem this is what a driver"
        NAMEFAIL="$NAMEFAIL\n      named as its case-variant looks like after the overwrite."
    fi
    lf=$(printf '%s' "${f%.cpp}" | tr 'A-Z' 'a-z')
    for t in $TARGETS; do
        lt=$(printf '%s' "$t" | tr 'A-Z' 'a-z')
        if [ "$lf" = "$lt" ]; then
            NAMEFAIL="$NAMEFAIL\n    target $t collides with library source $f ignoring case"
        fi
    done
done
for t in $TARGETS; do
    if [ ! -f "$SRC/$t.cpp" ]; then
        NAMEFAIL="$NAMEFAIL\n    missing driver source: $t.cpp"
    elif ! has_main "$SRC/$t.cpp"; then
        NAMEFAIL="$NAMEFAIL\n    $t.cpp does not define main"
    fi
done
if [ -n "$NAMEFAIL" ]; then
    echo "FAILED"
    printf "%b\n" "$NAMEFAIL"
    exit 1
fi
echo "ok"
echo

if [ "$SKIP_BUILD" != "1" ]; then
  cd "$SRC" || exit 1
  for d in $DIMS; do
    for t in $TARGETS; do
      printf "  building %-18s DIM=%s ... " "$t" "$d"
      if make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP="$COMP" DIM="$d" \
              EBASE="$t" -j"$JOBS" > "$LOGDIR/build_${t}_${d}d.log" 2>&1; then
          # Trust the timestamp, not make's exit status: a stale executable
          # from an interrupted link has bitten this tree once already.
          exe=$(ls -t "$SRC"/${t}${d}d.*MPI.ex 2>/dev/null | head -1)
          if [ -n "$exe" ] && [ -x "$exe" ]; then
              echo "ok  ($(basename "$exe"))"
          else
              echo "FAILED (no executable produced)"; exit 1
          fi
      else
          echo "FAILED -- see $LOGDIR/build_${t}_${d}d.log"
          grep -m5 "error:" "$LOGDIR/build_${t}_${d}d.log" | sed 's/^/      /'
          exit 1
      fi
    done
  done
  echo
fi

RC=0
run_suite () {   # name  script  dim  extra-args...
  local name=$1 script=$2 dim=$3; shift 3
  local log="$LOGDIR/${name}_${dim}d.log"
  echo "------------------------------------------------------------------"
  echo " $name  (DIM=$dim)"
  echo "------------------------------------------------------------------"
  # RANKS goes in the environment as well as positionally: the older suites
  # read it from "$@", the newer ones from the environment.
  if DIM="$dim" RANKS="$RANKS" WORK="$LOGDIR/work_${name}_${dim}d" \
        bash "$HERE/$script" $RANKS > "$log" 2>&1; then
      tail -1 "$log"
  else
      RC=1
      echo "FAILURES -- full log at $log"
      # Show both shapes of failure. A suite that ran and failed reports
      # lines beginning "  FAIL"; a suite that never started -- no usable
      # executable, a missing input -- writes an ERROR and nothing else, and
      # greping only for FAIL printed a single trailing line with the
      # diagnosis scrolled off above it.
      grep -E "^ERROR|^  ERROR" "$log" | head -12
      grep -A4 "  FAIL" "$log" | head -40
      tail -1 "$log"
  fi
  echo
}

for d in $DIMS; do
  run_suite complex_and_contour amrComplexMpiTest.sh    "$d"
  run_suite shapes              amrComplexShapeTest.sh  "$d"
  run_suite mef                 amrContourMefTest.sh    "$d"
  run_suite joint_stats         amrJointStatsTest.sh    "$d"
  run_suite grad_stats          amrGradStatsTest.sh     "$d"
  # 3D only: there is no Gaussian curvature of a plane curve, so sections 1
  # to 4 of that suite have nothing to check in 2D. Leaving it out of the 2D
  # pass is deliberate; running it there would report a vacuous pass.
  if [ "$d" = "3" ]; then
    run_suite curvature         amrCurvatureTest.sh     "$d"
  fi
done

echo "=================================================================="
if [ $RC -eq 0 ]; then echo " ALL SUITES PASS"; else echo " SOME SUITES FAILED -- see $LOGDIR"; fi
echo "=================================================================="
exit $RC
