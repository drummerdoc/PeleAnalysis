#!/usr/bin/env bash
#
# End-to-end thermodiffusive-lens analysis of a PeleLMeX plotfile.
#
# Two stages. The first derives gradient, alignment and curvature fields on the
# block-structured grid with coarse-fine-aware composite operators and writes a
# new plotfile; the second reduces that to volume-weighted joint statistics in
# range space and analyses them. Stage one is the expensive one and is the part
# that benefits from a GPU; stage two is a traversal and is MPI/OpenMP only.
#
# The whole point of keeping this as one script is that the thresholds and bin
# ranges are case properties, not defaults: T_u and Y_u in particular must be
# the prescribed inflow values, because the Lewis estimates are very sensitive
# to them. Set them per case, do not inherit them.
#
#   Tools/runLensPipeline.sh --infile plt29791 --case lean-h2 --tu 693 --yu 1.325792e-2
#
# Resource guide, from a 373M-cell four-level case (98 GB, 35 components):
# roughly 2.8 GB of output per component at full depth, and peak resident
# memory of about 20 component-equivalents.
#
#   finestLevel   cells     stage-1 output   stage-1 peak RSS
#        1         25 M         ~4 GB             ~6 GB
#        2        121 M        ~17 GB            ~30 GB
#        3        373 M        ~53 GB            ~95 GB
#
# Curvature needs the finest level available: max|kappa| scales as 1/dx at
# every level tested, so the extremum is always mesh noise, and the thermal
# front spans only 2.7 cells at level 1 against 10.7 at level 3. The lens and
# Lewis results, which involve no second derivatives, are converged at level 1.
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
DIM=${DIM:-3}
DTAG="${DIM}d"

INFILE=""; CASE=""; TU=""; YU=""
FINEST=${FINEST:-1}
ISPER=${ISPER:-"1 1 0"}
NPROC=${NPROC:-8}
MPIRUN=${MPIRUN:-mpiexec}
OUTDIR=""
STAGE=${STAGE:-all}
TMAX=""; KMAX=${KMAX:-100000}

usage () {
    sed -n '2,40p' "$0" | sed 's/^#\s\?//'
    cat <<'EOF'

Required:
  --infile PATH      plotfile to analyse
  --tu K             unburnt temperature (the prescribed inflow value)
  --yu FRAC          unburnt fuel mass fraction

Optional:
  --case NAME        label for outputs (default: basename of infile)
  --outdir DIR       where results go (default: ./lens_<case>_L<finest>)
  --finest N         finest AMR level to use (default 1)
  --isper "a b c"    periodicity (default "1 1 0")
  --nproc N          MPI ranks (default 8)
  --tmax K           upper temperature bin edge (default: from the data)
  --kmax VAL         curvature bin half-range (default 1e5)
  --stage grad|stats|analyse|all
EOF
    exit 1
}

while [ $# -gt 0 ]; do
    case $1 in
        --infile) INFILE=$2; shift 2 ;;
        --case)   CASE=$2;   shift 2 ;;
        --tu)     TU=$2;     shift 2 ;;
        --yu)     YU=$2;     shift 2 ;;
        --outdir) OUTDIR=$2; shift 2 ;;
        --finest) FINEST=$2; shift 2 ;;
        --isper)  ISPER=$2;  shift 2 ;;
        --nproc)  NPROC=$2;  shift 2 ;;
        --tmax)   TMAX=$2;   shift 2 ;;
        --kmax)   KMAX=$2;   shift 2 ;;
        --stage)  STAGE=$2;  shift 2 ;;
        -h|--help) usage ;;
        *) echo "unknown option: $1"; usage ;;
    esac
done
[ -n "$INFILE" ] && [ -n "$TU" ] && [ -n "$YU" ] || usage
[ -d "$INFILE" ] || { echo "no such plotfile: $INFILE"; exit 2; }
[ -n "$CASE" ] || CASE=$(basename "$INFILE")
[ -n "$OUTDIR" ] || OUTDIR="$(pwd)/lens_${CASE}_L${FINEST}"

pick () { ls -t "$SRC/${1}${DTAG}".*ex 2>/dev/null | head -1; }
launcher () {
    case $1 in *MPI*) command -v "$MPIRUN" >/dev/null 2>&1 \
        && printf '%s -n %s' "$MPIRUN" "$NPROC" ;;
    esac
}
GRD=$(pick amrGradStats); FIT=$(pick amrJointFit)
for pair in "amrGradStats:$GRD" "amrJointFit:$FIT"; do
    nm=${pair%%:*}; px=${pair#*:}
    if [ -z "$px" ] || [ ! -x "$px" ]; then
        echo "ERROR: no ${DTAG} build of $nm in $SRC"
        echo "  make -f GNUmakefile.AMRComplex USE_MPI=TRUE DIM=$DIM EBASE=$nm -j8"
        exit 2
    fi
done
GRDRUN="$(launcher "$GRD") $GRD"
FITRUN="$(launcher "$FIT") $FIT"
PY=${PYTHON:-python3}

mkdir -p "$OUTDIR" || exit 2
LOG="$OUTDIR/pipeline.log"
: > "$LOG"
G="$OUTDIR/${CASE}_g"

say  () { echo "$@" | tee -a "$LOG"; }
step () {
    local tag=$1 rc
    shift
    { echo; echo "### $tag"; echo "+ $*"; } >> "$LOG"
    "$@" >> "$LOG" 2>&1
    rc=$?
    [ $rc -eq 0 ] || { echo "    (exit $rc)" >> "$LOG"; say "  FAILED: $tag (see $LOG)"; }
    return $rc
}

say "case        : $CASE"
say "plotfile    : $INFILE"
say "finestLevel : $FINEST     periodicity: $ISPER     ranks: $NPROC"
say "T_u         : $TU K       Y_u: $YU"
say "outputs     : $OUTDIR"
say ""

# ---------------------------------------------------------------- stage one
if [ "$STAGE" = all ] || [ "$STAGE" = grad ]; then
    say "stage 1: composite gradients, alignment and curvature"
    step "amrGradStats" $GRDRUN infile="$INFILE" outfile="$G" \
         finestLevel="$FINEST" is_per="$ISPER" vars="temp Y(H2)" \
         curvOf=temp carry="HeatRelease I_R(H2)" || exit 1
    say "  wrote $G"
    say "  component ranges and the normal floor are in $LOG"
    grep -A3 "normal floor" "$LOG" | sed 's/^/    /' | tee -a /dev/null
fi
[ -d "$G" ] || { say "stage 1 output missing: $G"; exit 1; }

# Upper temperature edge from the data unless given.
if [ -z "$TMAX" ]; then
    TMAX=$(awk '/^   temp /{printf "%.0f", $3+10}' "$LOG" | tail -1)
    [ -n "$TMAX" ] || TMAX=2210
fi
say "temperature bins span $TU to $TMAX K"

# ---------------------------------------------------------------- stage two
if [ "$STAGE" = all ] || [ "$STAGE" = stats ]; then
    say ""
    say "stage 2: joint statistics"
    step "joint TY (fine)" $FITRUN infile="$G" is_per="$ISPER" \
         coords="temp Y(H2)" clo="$TU 0" chi="$TMAX $YU" cbins="128 192" \
         resp="HeatRelease I_R(H2)" outfile="$OUTDIR/ty_fine.dat"
    step "joint TY + gradients" $FITRUN infile="$G" is_per="$ISPER" \
         coords="temp Y(H2)" clo="$TU 0" chi="$TMAX $YU" cbins="128 192" \
         resp="gradMag_temp gradMag_Y(H2)" outfile="$OUTDIR/lewis2d.dat"
    step "joint T x curvature" $FITRUN infile="$G" is_per="$ISPER" \
         coords="temp meanCurv_temp" clo="$TU -$KMAX" chi="$TMAX $KMAX" \
         cbins="48 32" resp="Y(H2) gradMag_Y(H2) I_R(H2)" \
         outfile="$OUTDIR/postulate.dat"
    step "joint TY x curvature, flame-gated" $FITRUN infile="$G" is_per="$ISPER" \
         coords="temp Y(H2) meanCurv_temp gradMag_temp" \
         clo="$TU 0 -$KMAX 1e7" chi="$TMAX $YU $KMAX 1e9" cbins="20 20 16 4" \
         resp="I_R(H2) HeatRelease" outfile="$OUTDIR/ty_curv_gated.dat"
fi

# -------------------------------------------------------------- stage three
if [ "$STAGE" = all ] || [ "$STAGE" = analyse ]; then
    say ""
    say "stage 3: analysis"
    for m in frontier band; do
        [ -f "$OUTDIR/ty_fine.dat" ] || continue
        $PY "$HERE/lensAnalysis.py" $m "$OUTDIR/ty_fine.dat" \
            --Tu "$TU" --Yu "$YU" > "$OUTDIR/${m}.txt" 2>&1 \
            && say "  $OUTDIR/${m}.txt"
    done
    if [ -f "$OUTDIR/lewis2d.dat" ]; then
        $PY "$HERE/lensAnalysis.py" lewis "$OUTDIR/lewis2d.dat" \
            --Tu "$TU" --Yu "$YU" > "$OUTDIR/lewis.txt" 2>&1 \
            && say "  $OUTDIR/lewis.txt"
    fi
    if [ -f "$OUTDIR/postulate.dat" ]; then
        $PY "$HERE/lensAnalysis.py" postulate "$OUTDIR/postulate.dat" \
            --Tu "$TU" > "$OUTDIR/postulate.txt" 2>&1 \
            && say "  $OUTDIR/postulate.txt"
    fi
    if [ -f "$OUTDIR/ty_curv_gated.dat" ]; then
        $PY "$HERE/lensAnalysis.py" r2 "$OUTDIR/ty_curv_gated.dat" \
            > "$OUTDIR/r2_gated.txt" 2>&1 \
            && say "  $OUTDIR/r2_gated.txt"
    fi
fi

say ""
say "done. headline numbers:"
[ -f "$OUTDIR/frontier.txt" ] && grep -E "^\s+0\.870%|fixed t2" "$OUTDIR/frontier.txt" | sed 's/^/  /' | tee -a "$LOG"
[ -f "$OUTDIR/lewis.txt" ]    && grep -E "713-|Le_apparent" "$OUTDIR/lewis.txt" | sed 's/^/  /' | tee -a "$LOG"
exit 0
