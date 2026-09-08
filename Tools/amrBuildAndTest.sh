#!/usr/bin/env bash
#
# Build every AMRComplex/AMRContour executable in both dimensions and run all
# three test suites. One command, so there is no sequence to get wrong.
#
#   Tools/amrBuildAndTest.sh
#
# Environment:
#   COMP=llvm            compiler (default llvm, to match this tree's builds)
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
COMP=${COMP:-llvm}
JOBS=${JOBS:-8}
DIMS=${DIMS:-"2 3"}
RANKS=${RANKS:-"1 2 3 4"}
export MPIRUN=${MPIRUN:-mpiexec}
SKIP_BUILD=${SKIP_BUILD:-0}

TARGETS="amrComplexVerify amrContourVerify makeTestPlt surfMEFtoDAT"
LOGDIR=${LOGDIR:-$(pwd)/amr_test_logs}
mkdir -p "$LOGDIR"

echo "=================================================================="
echo " build : COMP=$COMP  DIMS=$DIMS  jobs=$JOBS"
echo " tests : ranks=$RANKS  launcher=$MPIRUN"
echo " logs  : $LOGDIR"
echo "=================================================================="
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
  if DIM="$dim" WORK="$LOGDIR/work_${name}_${dim}d" \
        bash "$HERE/$script" $RANKS > "$log" 2>&1; then
      tail -1 "$log"
  else
      RC=1
      echo "FAILURES -- full log at $log"
      grep -A4 "  FAIL" "$log" | head -40
      tail -1 "$log"
  fi
  echo
}

for d in $DIMS; do
  run_suite complex_and_contour amrComplexMpiTest.sh  "$d"
  run_suite mef                 amrContourMefTest.sh  "$d"
done

echo "=================================================================="
if [ $RC -eq 0 ]; then echo " ALL SUITES PASS"; else echo " SOME SUITES FAILED -- see $LOGDIR"; fi
echo "=================================================================="
exit $RC
