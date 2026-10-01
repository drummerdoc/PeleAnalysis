#!/usr/bin/env bash
#
# The standard joint-fit ladder for a premixed flame plotfile.
#
# One question, asked six ways: how much of the variance of heat release (and of
# fuel consumption) is explained by a pair of progress variables, and how does
# that answer move with (a) which second variable you pick and (b) how much of
# the AMR hierarchy you include.
#
#   The signal pair is (T, Y_fuel). Fuel is strongly sub-unity Lewis, so it is
#   transported ahead of the front and its contours do not follow the isotherms.
#
#   The control pair is (T, Y_product). Water is made where reaction happens and
#   diffuses at roughly the thermal rate, so these two should be nearly
#   parallel. Whatever R^2 gain the control shows is numerics, discretisation
#   and ordinary non-unity-Lewis background. The gain the signal pair shows ABOVE
#   the control is the thermodiffusive departure, with the common mode cancelled.
#   Running only the signal pair measures something; running both measures
#   something minus nothing, which is the part worth quoting.
#
#   The finestLevel ladder is the resolution question. If R^2 and the shape of
#   the conditional mean are the same at 2, 3 and 4 levels, the structure is
#   resolved. If they move, it is a grid feature. Bin counts are held FIXED
#   across the ladder on purpose: R^2 rises with bin count for free, so varying
#   both at once would confound resolution with binning.
#
# Usage:
#   Tools/amrJointFitLadder.sh <plotfile> [outprefix]
#
# Environment (defaults are PeleLMeX names and a lean hydrogen flame):
#   NP=8                 ranks
#   TVAR=temp            FUEL="Y(H2)"      PROD="Y(H2O)"
#   RESP="HeatRelease I_R(H2)"
#   TLO=690  THI=2210    FLO=0 FHI=0.0133  PLO=0 PHI=0.154
#   CBINS="128 128"      ISPER="1 1 0"     LEVELS="1 2 -1"   (-1 means all)
#   DIM=3
#
# Check the ranges against your own file before trusting the numbers: a large
# "volume dropped" in the output means the bins do not span the data and every
# R^2 below it is then conditional on an unstated subset of the domain.
# ---------------------------------------------------------------------------
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
PLT=${1:-}
OUT=${2:-jt}
if [ -z "$PLT" ] || [ ! -d "$PLT" ]; then
    echo "usage: $0 <plotfile> [outprefix]"; exit 1
fi

NP=${NP:-8}
DIM=${DIM:-3}
TVAR=${TVAR:-temp}
FUEL=${FUEL:-"Y(H2)"}
PROD=${PROD:-"Y(H2O)"}
RESP=${RESP:-"HeatRelease I_R(H2)"}
TLO=${TLO:-690};  THI=${THI:-2210}
FLO=${FLO:-0};    FHI=${FHI:-0.0133}
PLO=${PLO:-0};    PHI=${PHI:-0.154}
CBINS=${CBINS:-"128 128"}
ISPER=${ISPER:-"1 1 0"}
LEVELS=${LEVELS:-"1 2 -1"}
MPIRUN=${MPIRUN:-mpirun}

shopt -s nullglob
EXE=("$SRC"/amrJointFit${DIM}d.*MPI.ex)
shopt -u nullglob
[ ${#EXE[@]} -gt 0 ] || { echo "no amrJointFit${DIM}d MPI executable in $SRC"; exit 1; }
EXE=${EXE[0]}

echo "exe      : $EXE"
echo "plotfile : $PLT"
echo "ranks    : $NP      periodicity: $ISPER      bins: $CBINS"
echo "signal   : $TVAR x $FUEL        control: $TVAR x $PROD"
echo "responses: $RESP"
echo

run () {   # tag  coords  clo  chi  extra...
    local tag=$1 co=$2 lo=$3 hi=$4; shift 4
    local log="${OUT}_${tag}.log"
    printf "  %-22s " "$tag"
    if $MPIRUN -n "$NP" "$EXE" infile="$PLT" is_per="$ISPER" \
         coords="$co" clo="$lo" chi="$hi" cbins="${CB:-$CBINS}" resp="$RESP" \
         outfile="${OUT}_${tag}.dat" "$@" > "$log" 2>&1; then
        local drop
        drop=$(sed -n 's/^# volume dropped  :.*(\(.*\)%)/\1/p' "$log")
        echo "ok   (dropped ${drop:-?}%)"
    else
        echo "FAILED -- see $log"
        tail -3 "$log" | sed 's/^/      /'
    fi
}

for L in $LEVELS; do
    lab=$L; [ "$L" = "-1" ] && lab=all
    # The residual field rides along with the full-depth signal run rather than
    # costing a second pass over the whole hierarchy for the same numbers.
    extra=""
    [ "$L" = "-1" ] && extra="outplt=${OUT}_resid"
    run "sig_L${lab}" "$TVAR $FUEL" "$TLO $FLO" "$THI $FHI" finestLevel=$L $extra
    run "ctl_L${lab}" "$TVAR $PROD" "$TLO $PLO" "$THI $PHI" finestLevel=$L
done

# All three coordinates, to see whether the triple closes what the pair leaves
# open.
# A prefix assignment on a shell function call persists after it in bash, so
# the override is explicit and undone rather than relied on to be scoped.
CB="48 48 48"
run "triple" "$TVAR $FUEL $PROD" "$TLO $FLO $PLO" "$THI $FHI $PHI" finestLevel=-1
unset CB
echo

# ---------------------------------------------------------------------------
# Summary. For each run and response: R^2 for each single coordinate and for
# the pair, all at the finest binning, plus the pair re-read at the coarsening
# whose occupancy is closest to the richest single's. That last row is the
# comparison that means something -- R^2 rises with bin count for free, so the
# single and the pair have to be read at a matched number of occupied bins
# before the difference between them can be called information.
# ---------------------------------------------------------------------------
echo "=================================================================="
echo " explained variance"
echo "=================================================================="
for f in ${OUT}_*.log; do
    [ -f "$f" ] || continue
    tag=${f#${OUT}_}; tag=${tag%.log}
    awk -v tag="$tag" '
      /^# response/ { resp=$3; next }
      /^  [a-z]/ && NF==5 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
          key=$1; co=$2+0; occ=$3+0; r2=$4+0
          n=split(key, parts, ",")
          if (co==1) {
              if (n==1) { if (occ>bestocc[resp]) bestocc[resp]=occ
                          print " " tag "  " resp "  " key "  occ=" occ "  R2=" r2 > "/dev/stdout" }
              else      { pairkey[resp]=key; pairocc[resp]=occ; pairr2[resp]=r2 }
          }
          # Distance in log occupancy, not in occupancy. The coarsening ladder
          # is geometric, so an absolute difference always prefers the coarser
          # side -- 16 looks closer to 32 than 64 does, though both are a
          # factor of two away.
          if (n>1 && occ>0 && bestocc[resp]>0) {
              d = log(occ/bestocc[resp]); if (d<0) d=-d
              if (!(resp in bd) || d<bd[resp]) { bd[resp]=d; bk[resp]=key
                                                 bo[resp]=occ; br[resp]=r2; bc[resp]=co } }
      }
      END {
          for (r in pairkey)
              printf " %s  %s  %s  occ=%d  R2=%.6g      [vs single occ=%d: coarsen %d, occ=%d, R2=%.6g]\n",
                     tag, r, pairkey[r], pairocc[r], pairr2[r],
                     bestocc[r], bc[r], bo[r], br[r]
      }' "$f"
done
echo
echo "bin tables : ${OUT}_*.dat"
echo "residual   : ${OUT}_resid   (inbin, <resp>, <resp>_psi, <resp>_resid)"
echo "plot       : Tools/plotJointFit.py ${OUT}_sig_Lall.dat --log-volume"
