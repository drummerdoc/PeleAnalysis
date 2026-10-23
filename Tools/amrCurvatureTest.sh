#!/usr/bin/env bash
#
# Curvature and filtering tests for amrGradStats.
#
# Why this exists.
#
# The curvature fields had no test at all. meanCurv_ was being read off a
# production plotfile and reported, and the only check on it was whether the
# number looked plausible -- which for a second derivative of a noisy field is
# no check, because an implausible number and a wrong one are the same thing
# and a plausible wrong one is worse. Adding Gaussian and principal curvature
# made that untenable: those three are built from invariants of the same
# matrix, so a mistake in assembling them is invisible in a flat region, shows
# up only where the surface is curved, and the curved regions are the entire
# point.
#
# What makes a real test possible is that makeTestPlt's "dist" field is the
# signed distance to a sphere. Every one of its level sets is a sphere about
# the same centre, so at radius r the exact answers are div(n) = 2/r,
# K = 1/r^2 and k1 = k2 = 1/r -- at every r, not just on one isosurface. One
# plotfile therefore tests the whole chain over a continuum of curvatures with
# nothing to tune and no reference run to drift against.
#
# The filter needs the opposite kind of field, and a smooth analytic one
# cannot provide it: there is nothing for a low-pass filter to remove, so a
# test on it passes whatever the filter does, including nothing. Hence
# makeTestPlt's "noisy" field -- the same sphere displaced cell to cell by a
# stated fraction of a cell. It is the positive control: without filtering the
# curvature of it is wrong by an order of magnitude, and the test is that
# filtering recovers most of that while leaving the clean case alone.
#
# The five sections, and what each would catch:
#
#   1. the exact sphere. A wrong invariant, a swapped principal pair, a
#      discriminant clamped the wrong way, or a sign error in the normal.
#
#   2. the filter on a clean field. A filter that eats resolved structure.
#      The exact answer does not depend on the radial profile of a radially
#      symmetric field, so smoothing cannot change the right answer here --
#      any drift is the filter damaging what it should leave. This is the
#      section that caught a per-level filter amplitude, which puts a step in
#      the field at every coarse-fine interface: it took the error from 0.4%
#      to 4% on one pass and 10% on three, all of it on the interface.
#
#   3. the filter on a noisy field. A filter that does nothing. The error
#      must fall substantially and monotonically with the pass count.
#
#   4. the realizability meter. A meter that cannot distinguish the two
#      fields above is not measuring anything.
#
#   5. the kinetics analysis mode, against a table whose answer was planted.
#
# Build first:
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE DIM=3 EBASE=makeTestPlt -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE DIM=3 EBASE=amrGradStats -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE DIM=3 EBASE=amrCurvVerify -j8
#
# 3D only: there is no Gaussian curvature of a plane curve, so there is
# nothing for sections 1 to 4 to check in 2D.
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
MPIRUN=${MPIRUN:-mpiexec}
PY=${PYTHON:-python3}
NC=${NC:-64}

pick () { ls -t "$SRC/${1}3d".*ex 2>/dev/null | head -1; }
launcher () {
    # Test the FIRST WORD of MPIRUN, not the whole string. The documented way
    # to run these suites is MPIRUN="mpirun --oversubscribe", and command -v
    # on that whole string finds nothing, so this returned empty and the MPI
    # executable was launched bare -- which does not fail, it HANGS, because a
    # bare MPI binary waits for a launcher that is never coming. The older
    # suites have no such test and simply always use $MPIRUN, which is why
    # they were unaffected.
    case $1 in *MPI*) command -v "${MPIRUN%% *}" >/dev/null 2>&1 \
        && printf '%s -n 1' "$MPIRUN" ;;
    esac
}

MAKEPLT=$(pick makeTestPlt)
GRAD=$(pick amrGradStats)
CURV=$(pick amrCurvVerify)
for pair in "makeTestPlt:$MAKEPLT" "amrGradStats:$GRAD" "amrCurvVerify:$CURV"; do
    nm=${pair%%:*}; px=${pair#*:}
    if [ -z "$px" ] || [ ! -x "$px" ]; then
        echo "ERROR: no 3d build of $nm in $SRC"
        echo "  make -f GNUmakefile.AMRComplex DIM=3 EBASE=$nm -j8"
        exit 1
    fi
done
MK="$(launcher "$MAKEPLT") $MAKEPLT"
GR="$(launcher "$GRAD") $GRAD"
CV="$(launcher "$CURV") $CURV"

WORK=${WORK:-$(pwd)/amrCurvature_test}
mkdir -p "$WORK" && cd "$WORK" || exit 1

echo "makeTestPlt : $MAKEPLT"
echo "amrGradStats: $GRAD"
echo "amrCurvVerify: $CURV"
echo "workdir     : $WORK   base cells: $NC"
echo

FAIL=0
pass () { echo "  PASS  $1"; }
fail () { echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

# One cell of the finest level, for expressing the noise amplitude. Two levels
# at ref_ratio 2, so the finest spacing is 1/(2*NC) on a unit domain.
DXF=$($PY -c "print(1.0/(2*$NC))")
AMP=$($PY -c "print(0.2/(2*$NC))")

err_of () { awk '/div\(n\) vs/{print $NF}' "$1"; }
sev_of () { awk '/unrealizable/{print $NF}' "$1"; }

echo "generating the two hierarchies (noise amplitude $AMP = 0.2 dx_fine)..."
rm -rf plt_clean plt_noisy
$MK n_cell=$NC nlevels=2 refine_shape=box refine_anchor=center \
    refine_frac=0.5 noise_amp=$AMP outfile=plt_clean > gen.log 2>&1
if [ ! -d plt_clean ]; then
    fail "generator produced nothing, see $WORK/gen.log"
    echo "$FAIL failure(s)"; exit "$FAIL"
fi
# One plotfile carries both fields: "dist" is clean, "noisy" is the same
# sphere perturbed. Nothing else differs between them, which is what makes
# sections 2 and 3 a controlled comparison rather than two separate runs.
echo

# --------------------------------------------------------------------------
echo "=== 1. the exact sphere ==============================================="
rm -rf g_clean_0
$GR infile=plt_clean vars="dist" curvOf=dist gaussCurv=1 is_per="0 0 0" \
    outfile=g_clean_0 > s_clean_0.log 2>&1
$CV infile=g_clean_0 field=dist rlo=0.15 rhi=0.35 > v_clean_0.log 2>&1
if grep -q "RESULT: PASS" v_clean_0.log; then
    pass "div(n), K, k1, k2 against the exact sphere: $(err_of v_clean_0.log) worst-cell"
else
    fail "exact sphere: $(grep -E 'div\(n\)|K vs|k1,k2|k2 > k1' v_clean_0.log | tr -s ' ')"
fi
echo

# --------------------------------------------------------------------------
echo "=== 2. the filter leaves a resolved surface alone ====================="
E0=$(err_of v_clean_0.log)
for n in 1 3 5; do
    rm -rf g_clean_$n
    $GR infile=plt_clean vars="dist" curvOf=dist gaussCurv=1 smooth_passes=$n \
        is_per="0 0 0" outfile=g_clean_$n > s_clean_$n.log 2>&1
    $CV infile=g_clean_$n field=smooth${n}_dist rlo=0.15 rhi=0.35 \
        > v_clean_$n.log 2>&1
    EN=$(err_of v_clean_$n.log)
    # The exact answer does not depend on the radial profile, so the error
    # must not grow. A factor of two of headroom, against 10x for the defect
    # this is here to catch.
    if ! grep -q "RESULT: PASS" v_clean_$n.log; then
        fail "$n pass(es) on the clean sphere: error $EN, was $E0"
    elif [ "$($PY -c "print(1 if $EN > 2.0*$E0 else 0)")" = "1" ]; then
        fail "$n pass(es) on the clean sphere: error grew $E0 -> $EN"
    else
        pass "$n pass(es): error $EN, unsmoothed $E0"
    fi
done
echo

# --------------------------------------------------------------------------
echo "=== 3. the filter removes cell-scale noise ============================"
rm -rf g_noisy_0
$GR infile=plt_clean vars="noisy" curvOf=noisy gaussCurv=1 is_per="0 0 0" \
    outfile=g_noisy_0 > s_noisy_0.log 2>&1
$CV infile=g_noisy_0 field=noisy rlo=0.15 rhi=0.35 tol=1e9 > v_noisy_0.log 2>&1
N0=$(err_of v_noisy_0.log)
# Without filtering the curvature of the perturbed surface must be badly
# wrong. If it is not, the perturbation is not reaching the curvature and
# everything below is vacuous.
if [ "$($PY -c "print(1 if $N0 > 1.0 else 0)")" = "1" ]; then
    pass "unfiltered noisy field is wrong by $N0 (relative), as it must be"
else
    fail "unfiltered noisy field errs only $N0; the control is not controlling"
fi
PREV=$N0
for n in 1 3 8; do
    rm -rf g_noisy_$n
    $GR infile=plt_clean vars="noisy" curvOf=noisy gaussCurv=1 smooth_passes=$n \
        is_per="0 0 0" outfile=g_noisy_$n > s_noisy_$n.log 2>&1
    $CV infile=g_noisy_$n field=smooth${n}_noisy rlo=0.15 rhi=0.35 tol=1e9 \
        > v_noisy_$n.log 2>&1
    NN=$(err_of v_noisy_$n.log)
    if [ "$($PY -c "print(1 if $NN < $PREV else 0)")" = "1" ]; then
        pass "$n pass(es) on the noisy field: error $PREV -> $NN"
    else
        fail "$n pass(es) on the noisy field: error did not fall, $PREV -> $NN"
    fi
    PREV=$NN
done
RED=$($PY -c "print(f'{$N0/$PREV:.1f}')")
if [ "$($PY -c "print(1 if $N0/$PREV > 3.0 else 0)")" = "1" ]; then
    pass "eight passes recover a factor of $RED"
else
    fail "eight passes recover only a factor of $RED; the filter is not working"
fi
echo

# --------------------------------------------------------------------------
echo "=== 4. the realizability meter separates the two fields ==============="
SC=$(sev_of s_clean_0.log)
SN=$(sev_of s_noisy_0.log)
if [ -z "$SC" ] || [ -z "$SN" ]; then
    fail "meter: no severity reported; gaussCurv may not have run"
elif [ "$($PY -c "print(1 if $SC > 1e-3 else 0)")" = "1" ]; then
    fail "meter: clean sphere reports severity $SC, which is not a clean baseline"
elif [ "$($PY -c "print(1 if $SN < 3.0*$SC else 0)")" = "1" ]; then
    fail "meter: noisy $SN against clean $SC; cannot tell them apart"
else
    pass "meter: clean $SC, noisy $SN, ratio $($PY -c "print(f'{$SN/$SC:.1f}')")"
fi
echo

# --------------------------------------------------------------------------
echo "=== 5. the kinetics mode, against a planted answer ===================="
# A table built so that Phi = 1 at a known temperature on each side of the
# band, with a known rich-side lift. Checks the pool arithmetic, the side
# split and the log interpolation together; a linear interpolation of a
# quantity that falls by orders of magnitude lands several tens of K low,
# which is the size of the effect being looked for in the real data.
cat > mkkin.py <<'PYEOF'
import sys
import numpy as np
dst = sys.argv[1]
TC, WID, LIFT = 1100.0, 120.0, 1.2
Ts = np.arange(700.0, 1800.0, 10.0)
rows = []
for T in Ts:
    for side, dy in ((+1, +1.0e-4), (-1, -1.0e-4)):
        phi = np.exp((TC - T)/WID)*(LIFT if side > 0 else 1.0)
        br = 1.0e-5*np.exp((T - 1600.0)/200.0)
        tm = br*phi
        rows.append([T, 5.0e-3 + dy, 1.0, 100.0,
                     0.5*br, 0.05*br, 0.2*br, 0.02*br, 0.3*br, 0.03*br,
                     0.8*tm, 0.08*tm, 0.2*tm, 0.02*tm,
                     1.0e-3*(T - TC), 1.0e-4])
hdr = ("# temp Y(H2) volume cells Y(H)_mean Y(H)_sd Y(O)_mean Y(O)_sd "
       "Y(OH)_mean Y(OH)_sd Y(HO2)_mean Y(HO2)_sd Y(H2O2)_mean Y(H2O2)_sd "
       "I_R(H)_mean I_R(H)_sd")
with open(dst, 'w') as fh:
    fh.write(hdr + "\n")
    np.savetxt(fh, np.array(rows), fmt='%.10g')
print(f"{TC:.0f} {TC + WID*np.log(LIFT):.1f} {LIFT:.2f}")
PYEOF
PLANT=$($PY -I mkkin.py kin_synth.dat 2>/dev/null)
WANT_LEAN=$(echo "$PLANT" | awk '{print $1}')
WANT_RICH=$(echo "$PLANT" | awk '{print $2}')
WANT_LIFT=$(echo "$PLANT" | awk '{print $3}')
$PY "$HERE/lensAnalysis.py" kinetics kin_synth.dat --Tu 693 \
    --bands 700,800,900,1000,1100,1200,1300,1400,1500,1600,1700 \
    > kin.log 2>&1
GOT_LEAN=$(awk '/lean side/{print $(NF-1)}' kin.log)
GOT_RICH=$(awk '/rich side/{print $(NF-1)}' kin.log)
# exit on the first match: the band label appears again in the I_R(H) table
# further down, and taking both turns the comparison into a syntax error.
GOT_LIFT=$(awk '/^ +1000-1100/{print $NF; exit}' kin.log)
if [ -z "$GOT_LEAN" ] || [ -z "$GOT_RICH" ]; then
    fail "kinetics: no crossing reported, see $WORK/kin.log"
else
    ok=1
    for p in "$GOT_LEAN:$WANT_LEAN:lean" "$GOT_RICH:$WANT_RICH:rich"; do
        g=${p%%:*}; rest=${p#*:}; w=${rest%%:*}; lab=${rest#*:}
        if [ "$($PY -c "print(1 if abs($g-$w) > 15.0 else 0)")" = "1" ]; then
            fail "kinetics: $lab crossing $g K, planted $w K"
            ok=0
        fi
    done
    if [ "$($PY -c "print(1 if abs($GOT_LIFT-$WANT_LIFT) > 0.03 else 0)")" = "1" ]; then
        fail "kinetics: rich/lean ratio $GOT_LIFT, planted $WANT_LIFT"
        ok=0
    fi
    [ "$ok" = "1" ] && pass "kinetics: crossings $GOT_LEAN/$GOT_RICH K against planted $WANT_LEAN/$WANT_RICH, ratio $GOT_LIFT"
fi
if grep -q "changes sign between" kin.log; then
    pass "kinetics: net I_R(H) sign change located"
else
    fail "kinetics: net I_R(H) sign change not reported"
fi
echo

echo "======================================================================="
if [ "$FAIL" -eq 0 ]; then
    echo "all curvature tests pass"
else
    echo "$FAIL failure(s)"
fi
exit "$FAIL"
