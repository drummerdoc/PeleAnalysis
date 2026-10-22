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
# finestLevel does not change the resolution of the simulation. It chooses how
# much of this plotfile's hierarchy the analysis reads. The levels above it are
# discarded, and their data has already been averaged down into the levels
# below, so a level-1 analysis is a coarsened view of the same solution. The
# sequence therefore measures how much a diagnostic loses when it is computed
# from a coarsened view. It says nothing about whether the simulation itself
# resolves anything, and it is not a grid-convergence study.
#
# Measured on plt29791, level 1 against level 3:
#
#   the lens is insensitive. lens/optimal 0.929 against 0.935, enrichment
#   42.9x against 44.6x, and the same two thresholds to four figures. The lens
#   can be located on a coarsened view and the answer carries over.
#
#   the Lewis estimate is not. The flux ratio rises 21% in the coldest band,
#   and its variation across the preheat zone collapses from 28% of the mean
#   to 7%. Averaging down flattens |grad Y| more than |grad T|, because the
#   fuel layer is the thinner of the two, so a coarsened view biases the ratio
#   low and bends what should be a flat profile. Run it at full depth.
#
#   curvature admits more extremum with every level added. Whether what is
#   being admitted is structure or point-to-point noise is a separate
#   question, and a smoothing sweep answers it, not a level sweep.
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
META="$OUTDIR/stage1.meta"
G="$OUTDIR/${CASE}_g"

# The log is appended to, never truncated. A --stage stats or --stage analyse
# run used to erase the log it was about to read its bin ranges out of, and
# then silently fall back to the compiled-in defaults. That is how an L3 run
# ended up binning temperature to 2210 K while the L2 run that wrote the same
# field used 2212, and the two sets of thresholds were then not comparable.
# Anything a later stage needs from stage one is written to stage1.meta as
# well, so a later stage never depends on scraping prose.
{
    echo
    echo "======================================================================"
    echo "run $(date '+%Y-%m-%d %H:%M:%S')  stage=$STAGE  finest=$FINEST  ranks=$NPROC"
    echo "======================================================================"
} >> "$LOG"

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

    # The two numbers later stages need: the upper temperature edge and the
    # peak |grad T| that sets the flame gate. Taken from the last occurrence,
    # so a log carrying several runs yields this run's values.
    m_tmax=$(awk '/^   temp /{v=$3} END{if (v != "") printf "%.0f", v+10}' "$LOG")
    m_gmax=$(awk '/normal floor for temp/{v=$9} END{gsub(/[^0-9eE+.-]/,"",v); print v}' "$LOG")
    m_floor=$(awk '/normal floor for temp/{v=$6} END{print v}' "$LOG")
    {
        echo "TMAX_META=${m_tmax}"
        echo "GMAX_META=${m_gmax}"
        echo "FLOOR_META=${m_floor}"
    } > "$META"
fi
[ -d "$G" ] || { say "stage 1 output missing: $G"; exit 1; }

TMAX_META=""; GMAX_META=""; FLOOR_META=""
[ -f "$META" ] && . "$META"

# Precedence for every derived bin edge: what was asked for on the command
# line, then what stage one measured, then a compiled-in default that is only
# ever right by luck and says so.
if [ -z "$TMAX" ]; then
    if [ -n "$TMAX_META" ]; then
        TMAX=$TMAX_META
    else
        TMAX=2210
        say "WARNING: no stage1.meta and no --tmax; using the default upper edge $TMAX K."
        say "         Bin edges will not match a run that had the measured value."
    fi
fi
say "temperature bins span $TU to $TMAX K"

# The flame gate. This is a window on |grad T| that keeps the flame and drops
# the quiescent gas on either side of it. Its edges have to follow the data:
# the gate was once written as a fixed 1e7 to 1e9, and at a peak |grad T| of
# 9.65e7 every cell in range fell into the first of its four bins, so the
# coordinate was inert and the gated reduction was quietly an ungated one.
# Scaling off the measured peak keeps the four bins populated at any case.
if [ -n "$GMAX_META" ]; then
    GLO=$(awk -v g="$GMAX_META" 'BEGIN{printf "%.6e", g/10.0}')
    GHI=$(awk -v g="$GMAX_META" 'BEGIN{printf "%.6e", g*1.05}')
    say "flame gate  : |grad T| in $GLO to $GHI  (peak $GMAX_META)"
else
    GLO=1e7; GHI=1e8
    say "WARNING: no measured peak |grad T|; gating on the default $GLO to $GHI."
    say "         Check that gradMag_temp has more than one occupied bin in r2_gated.txt."
fi

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
         clo="$TU 0 -$KMAX $GLO" chi="$TMAX $YU $KMAX $GHI" cbins="20 20 16 4" \
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
# Pick the rows by what they are, not by the value they happened to have on
# the case this script was first written for. A grep for the literal 0.870%
# matched nothing the moment a budget list changed.
if [ -f "$OUTDIR/frontier.txt" ]; then
    awk '/^ *[0-9]/{n++; if (n==3) print} /fixed t2/{print}' \
        "$OUTDIR/frontier.txt" | sed 's/^/  /' | tee -a "$LOG"
fi
if [ -f "$OUTDIR/lewis.txt" ]; then
    awk '/^ *[0-9]+-[0-9]+/{n++; if (n==1) print} /Le_apparent|sign of S/{print}' \
        "$OUTDIR/lewis.txt" | sed 's/^/  /' | tee -a "$LOG"
fi
if [ -f "$OUTDIR/r2_gated.txt" ]; then
    ng=$(awk '$1=="gradMag_temp"{print $2; exit}' "$OUTDIR/r2_gated.txt")
    if [ -n "$ng" ] && [ "$ng" -le 1 ]; then
        say "  WARNING: the flame gate has $ng occupied bin, so the gated reduction is ungated."
    fi
fi
exit 0
