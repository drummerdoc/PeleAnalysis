#!/usr/bin/env bash
#
# Non-convex and wall-adjacent refinement for AMRComplex.
#
# Why this exists.
#
# Every synthetic hierarchy the suite built until now was one centred, convex
# box per level, strictly inset from its parent and from the domain. A real
# regridded hierarchy is nothing like that. A representative lean-hydrogen
# flamesheet plotfile has 5213 boxes at level 1 filling 64% of the level, 5012
# at level 2 filling 36%, 17818 at level 3 filling 12%, and levels 1 and 2 sit
# hard against the domain walls. So the suite was passing on the one geometry
# the code finds easy and had never once been shown the geometry it is for.
#
# Two properties are new here, and each breaks a different assumption:
#
#   a non-convex footprint. A base cell in the low-side stitch shell can then
#   draw some of its corners from a same-level neighbour box and the rest from
#   the coarser level. That produces element shapes -- 7 distinct corners out
#   of 8 in 3D -- which a single convex box cannot produce at all, and it is
#   where a step in the coarse-fine interface either tiles or leaks.
#
#   Contact with the domain wall. A refined level reaching a non-periodic wall
#   makes the complex stop dx_fine/2 short of it rather than dx_coarse/2, and
#   through the periodic wrap it makes a coarse-fine stitch whose two sides are
#   at opposite ends of the domain.
#
#   A mismatched periodic wrap. A level spanning the domain touches both sides
#   of a periodic wall, and its coverage on the two sides is uncorrelated, so
#   the wrap is a mixture of same-level and coarse-fine stitches at different
#   transverse positions. refine_shape=random with refine_frac=1 makes that on
#   purpose; it is the regime a regridded doubly-periodic hierarchy is mostly
#   made of.
#
# The gate is run with every direction periodic. That is not a convenience: on
# a fully periodic domain the union of the complex is the whole domain, so the
# reference volume is exactly Lx*Ly*Lz with no inset term and no dependence on
# which level happens to touch which wall. The residual is then a sharp
# identity rather than a bound, and any nonzero value is a real gap or overlap.
# The complex is built topologically, so is_per changes only the wrap
# identification -- nothing about the stitch geometry being tested here.
#
# Build (see Tools/amrComplexMpiTest.sh for the full set):
#
#   cd Src
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=makeTestPlt -j8
#
# Then, from any writable directory:
#
#   /path/to/PeleAnalysis/Tools/amrComplexShapeTest.sh
#   DIM=2 Tools/amrComplexShapeTest.sh
#   RANKS="1 4" MPIRUN="mpirun --oversubscribe" Tools/amrComplexShapeTest.sh
#
# On failure the script re-runs the failing case at one rank with
# report_bad_faces=8, which prints the level and index of the corners of the
# first few faces whose incidence is not two. At one rank on a fully periodic
# domain every such face is a defect, so that output localises the gap.
# ---------------------------------------------------------------------------
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
RANKS=${RANKS:-"1 4"}
MPIRUN=${MPIRUN:-mpiexec}
DIM=${DIM:-3}
DTAG="${DIM}d"

shopt -s nullglob
VERIFY=("$SRC"/amrComplexVerify${DTAG}.*MPI.ex)
MAKEPLT=("$SRC"/makeTestPlt${DTAG}.*MPI.ex)
shopt -u nullglob

if [ ${#VERIFY[@]} -eq 0 ] || [ ${#MAKEPLT[@]} -eq 0 ]; then
    echo "ERROR: no MPI executables found in $SRC."
    echo "Build them with USE_MPI=TRUE (see the header of this script)."
    exit 1
fi
VERIFY=${VERIFY[0]}
MAKEPLT=${MAKEPLT[0]}

WORK=${WORK:-$(pwd)/amrComplex_shape_test}
mkdir -p "$WORK" && cd "$WORK" || exit 1

if [ "$DIM" = "2" ]; then ALLPER="1 1"; else ALLPER="1 1 1"; fi

echo "verify : $VERIFY"
echo "workdir: $WORK"
echo "dim    : $DIM   is_per: $ALLPER   ranks: $RANKS"
echo

invariants () {
    sed -n '/=== AMRComplex verify ===/,/RESULT/p' | grep -v '^$' | grep -v '^~'
}

FAIL=0
pass () { echo "  PASS  $1"; }
fail () { echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

# Parallel arrays: tag, shape, anchor, levels, base cells, extent, fill.
# refine_blk is held at 8 so that a 32-cell base has enough blocks for the
# patterns to punch anything out at all. shape and anchor may be multi-word,
# one entry per refined level, and are deliberately left unquoted so ParmParse
# sees them as lists.
#
# rnd_span is the case that matters most: refine_frac=1 makes level 1 span the
# whole domain, so it touches both sides of every periodic wall while its
# coverage on the two sides is uncorrelated. The wrap is then a mixture of
# same-level and coarse-fine stitches at different transverse positions, which
# is what a regridded doubly-periodic hierarchy actually looks like and what no
# symmetric patch produces.
#
# The three-level cases keep level 2 interior (anchor "... center"). A level
# that touches a periodic boundary while its parent does not reach the far side
# of that boundary has no parent across the wrap at all, and that is a hard
# error rather than a shape to test -- section 3 covers it as a negative
# control.
TAG=(box_lo box_hi lsh_ctr lsh_lo  chk_ctr chk_lo  rnd_ctr rnd_span box_lo_L3  rnd_L3          chk_ctr_L3)
SHP=(box    box    lshape  lshape  checker checker random  random   "box box"  "random random" "checker checker")
ANC=(lo     hi     center  lo      center  lo      center  lo       "lo center" "lo center"    "center center")
NLV=(2      2      2       2       2       2       2       2        3          3               3)
NCL=(32     32     32      32      32      32      32      32       48         48              48)
FRC=(0.5    0.5    0.5     0.5     0.5     0.5     0.5     1.0      0.5        0.5             0.5)
FIL=(0.6    0.6    0.6     0.6     0.6     0.6     0.6     0.6      0.6        0.6             0.6)

echo "generating hierarchies..."
for i in "${!TAG[@]}"; do
    t=${TAG[$i]}
    $MPIRUN -n 1 "$MAKEPLT" n_cell=${NCL[$i]} nlevels=${NLV[$i]} \
        refine_shape=${SHP[$i]} refine_anchor=${ANC[$i]} refine_blk=8 \
        refine_frac=${FRC[$i]} refine_fill=${FIL[$i]} outfile=shp_$t \
        > gen_$t.log 2>&1
    if [ ! -d "shp_$t" ]; then
        fail "$t : generator failed, see gen_$t.log"
        continue
    fi
    grep -E "^   (level|refine_shape)" gen_$t.log | sed "s/^/  [$t]/"
done
echo

# --------------------------------------------------------------------------
# 1. The volume gate, on every shape.
# --------------------------------------------------------------------------
echo "=== 1. volume closes on non-convex and wall-adjacent levels ==========="
for i in "${!TAG[@]}"; do
    t=${TAG[$i]}
    [ -d "shp_$t" ] || continue
    $MPIRUN -n 1 "$VERIFY" infile=shp_$t is_per="$ALLPER" varnames=dist \
        face_census=1 > ver_$t.log 2>&1
    if grep -q "RESULT: PASS" ver_$t.log; then
        nm=$(sed -n 's/^~ *non-manifold (>2) *: *//p' ver_$t.log)
        du=$(sed -n 's/^  duplicate elts *: *//p' ver_$t.log | awk '{print $1}')
        if [ "${nm:-0}" != "0" ]; then
            fail "$t : volume closes but $nm non-manifold face(s)"
        elif [ "${du:-0}" != "0" ]; then
            fail "$t : volume closes but $du element(s) generated more than once"
        elif ! grep -q "closed complex : yes" ver_$t.log; then
            fail "$t : volume closes but the complex is not face-closed"
        else
            sh=$(sed -n '/shape histogram/,/^  fields/p' ver_$t.log \
                 | sed -n 's/^    //p' | tr '\n' ' ')
            pass "$t : $sh"
        fi
    else
        fail "$t : $(sed -n 's/^~ relative error *: *//p' ver_$t.log)"
        echo "        localising, first 8 faces with incidence != 2:"
        $MPIRUN -n 1 "$VERIFY" infile=shp_$t is_per="$ALLPER" varnames=dist \
            report_bad_faces=8 2>&1 | grep '^~ *\(face\|element\)' \
            | sed 's/^/        /'
    fi
done
echo

# --------------------------------------------------------------------------
# 2. Rank invariance, on the same shapes.
#
# The ownership rule for the stitch shell -- lowest-indexed visiting box wins
# -- is the part most likely to go wrong on a non-convex footprint, because a
# shell cell can now be visited by several boxes on several ranks. If it does
# go wrong, the element set differs between rank counts and the checksums
# diverge even where the volume still happens to close.
# --------------------------------------------------------------------------
echo "=== 2. rank invariance on the same shapes ============================="
for i in "${!TAG[@]}"; do
    t=${TAG[$i]}
    [ -d "shp_$t" ] || continue
    ref=""
    for n in $RANKS; do
        $MPIRUN -n "$n" "$VERIFY" infile=shp_$t is_per="$ALLPER" varnames=dist \
            2>/dev/null | invariants > "inv_${t}_n$n.txt"
        if [ -z "$ref" ]; then ref="inv_${t}_n$n.txt"; continue; fi
        if cmp -s "$ref" "inv_${t}_n$n.txt"; then
            pass "$t : $n ranks identical to $(basename "$ref")"
        else
            fail "$t : $n ranks differ"
            diff "$ref" "inv_${t}_n$n.txt" | sed 's/^/        /'
        fi
    done
done
echo

# --------------------------------------------------------------------------
# 3. Negative control: the orphan-ghost abort.
#
# A level that touches a periodic boundary while its parent does not reach the
# far side of that boundary has, across the wrap, a neighbour two or more levels
# coarser. The stitch spans one level, so that ghost finds no parent at all and
# every element that would bridge it is silently not generated -- a gap with
# nothing to indicate it. Until this control existed the code produced exactly
# that, and the only symptom was a volume residual on real data.
#
# So this case must ABORT, and it must abort naming the cause. A run that
# completes here is the failure.
# --------------------------------------------------------------------------
echo "=== 3. negative control: level 2 at a wall its parent does not wrap ===="
$MPIRUN -n 1 "$MAKEPLT" n_cell=48 nlevels=3 refine_shape=box box \
    refine_anchor=lo lo refine_blk=8 refine_frac=0.5 outfile=shp_orphan \
    > gen_orphan.log 2>&1
if [ ! -d "shp_orphan" ]; then
    fail "orphan control : generator failed, see gen_orphan.log"
else
    $MPIRUN -n 1 "$VERIFY" infile=shp_orphan is_per="$ALLPER" varnames=dist \
        > ver_orphan.log 2>&1
    if grep -q "no owner at any level" ver_orphan.log; then
        pass "orphan control : aborted, and named the cause"
    elif grep -q "RESULT:" ver_orphan.log; then
        fail "orphan control : ran to completion instead of aborting"
        sed -n 's/^~ relative error/        relative error/p' ver_orphan.log
    else
        fail "orphan control : failed for some other reason, see ver_orphan.log"
        tail -3 ver_orphan.log | sed 's/^/        /'
    fi
fi
echo

# --------------------------------------------------------------------------
# 4. Negative control: the ownership rule without the periodic wrap.
#
# Base-cell ownership is corner order -- a base cell belongs to the box holding
# the lowest-numbered corner offset that lands in the level's footprint, taken
# modulo the domain. legacy_ownership=1 drops the modulo, which is the defect
# that rule replaced: on a torus one base cell has up to 2^DIM literal indices,
# and a level spanning two periodic directions then has each wrap-edge element
# claimed by every box that can see any one of them.
#
# So rnd_span must FAIL with legacy_ownership=1, and fail in all three ways at
# once: duplicated elements, a non-zero interior face parity, and a POSITIVE
# volume residual. A check that cannot be made to fail is not a check.
# --------------------------------------------------------------------------
echo "=== 4. negative control: ownership without the periodic wrap ==========="
if [ -d "shp_rnd_span" ]; then
    $MPIRUN -n 1 "$VERIFY" infile=shp_rnd_span is_per="$ALLPER" varnames=dist \
        face_census=1 legacy_ownership=1 > ver_legacy.log 2>&1
    du=$(sed -n 's/^  duplicate elts *: *//p' ver_legacy.log | awk '{print $1}')
    if grep -q "RESULT: PASS" ver_legacy.log; then
        fail "legacy control : passed, so the volume gate cannot see the defect"
    elif [ "${du:-0}" = "0" ]; then
        fail "legacy control : failed but reported no duplicate elements"
    elif grep -q "closed complex : yes" ver_legacy.log; then
        fail "legacy control : failed but reported the complex face-closed"
    else
        pass "legacy control : $du duplicates, not face-closed, residual $(sed -n 's/^~ relative error *: *//p' ver_legacy.log)"
    fi
else
    fail "legacy control : shp_rnd_span was never generated"
fi
echo

echo "======================================================================="
if [ "$FAIL" -eq 0 ]; then
    echo "all shape tests pass"
else
    echo "$FAIL failure(s)"
fi
exit "$FAIL"
