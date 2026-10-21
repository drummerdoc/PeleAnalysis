#!/usr/bin/env bash
#
# Verification for AMRJointStats.
#
# What is being claimed, and what each test measures.
#
# The tool reports R^2 = 1 - E[var(q | coords)] / var(q), the share of the
# variance of a response that the best possible function of the chosen
# coordinates explains. Four things can be wrong with that: the volume measure,
# the variance arithmetic, the marginalisation that produces coordinate
# subsets, and the coarsening that produces bin-width dependence. Each has an
# exact or analytic reference here, so none of them is taken on trust.
#
#   Volume. The traversal weights each VALID cell by the volume of its own
#   cell, and uncovered cells partition the domain, so the volume in range must
#   equal the domain volume exactly when the bins span the data. This is the
#   same integer identity AMRComplex::verify gates on, re-measured through a
#   different code path.
#
#   Negative control. makeTestPlt writes xco, yco, zco holding the coordinates
#   exactly. The volume-weighted mean of yco over any union of whole cells is
#   the midpoint rule for a linear function, which is exact on any box
#   partition -- so conditioning yco on xco explains nothing and R^2 must be
#   zero, not small. With bins aligned to the coarsest grid no cell straddles a
#   bin boundary and the result is exact at every coarsening. A check that
#   cannot return zero on independent data is not measuring dependence.
#
#   Convergence. dist is exactly a function of (xco,yco,zco), so R^2 -> 1, and
#   the residual is the within-bin variation of a smooth field: 1 - R^2 = O(h^2)
#   in the bin width. Successive coarsenings must therefore show a ratio of 4.
#   Two coarsenings are excluded at each end on purpose: the finest, where the
#   bin is no wider than the coarsest cell and the asymptotics do not apply, and
#   the coarsest, where 1 - R^2 is no longer small.
#
#   Marginalisation and coarsening. Summing the table over an axis must give
#   exactly the lower-dimensional subset, and summing adjacent bins must give
#   exactly the coarser binning. Both are checked against runs made directly at
#   that subset and that bin count, to the last printed digit. These are the
#   law-of-total-variance steps, and getting them wrong would quietly inflate
#   every R^2.
#
#   Symmetry. The sphere is isotropic but the hierarchy is not -- level 1 is a
#   hash-selected spanning region anchored at the low corner. So the three
#   single-coordinate R^2 values must agree, and the three pairs must agree,
#   to within the asymmetry of the refinement pattern and nothing more. This
#   catches an axis-ordering or stride error that no other test here would see.
#
# Build:
#   cd Src
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=amrJointFit -j8
#   make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm EBASE=makeTestPlt -j8
#
# Then from any writable directory:
#   /path/to/PeleAnalysis/Tools/amrJointStatsTest.sh
#   DIM=2 Tools/amrJointStatsTest.sh
#   RANKS="1 4" MPIRUN="mpirun --oversubscribe" Tools/amrJointStatsTest.sh
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
JS=("$SRC"/amrJointFit${DTAG}.*MPI.ex)
MAKEPLT=("$SRC"/makeTestPlt${DTAG}.*MPI.ex)
shopt -u nullglob
if [ ${#JS[@]} -eq 0 ] || [ ${#MAKEPLT[@]} -eq 0 ]; then
    echo "ERROR: no MPI executables found in $SRC."
    echo "Build them with USE_MPI=TRUE (see the header of this script)."
    exit 1
fi
JS=${JS[0]}; MAKEPLT=${MAKEPLT[0]}

WORK=${WORK:-$(pwd)/amrJointStats_test}
mkdir -p "$WORK" && cd "$WORK" || exit 1

# Per-dimension, because cbins takes one value per coordinate and a stray
# third entry in 2D aborts the run rather than being ignored.
if [ "$DIM" = "2" ]; then
    ALLPER="1 1"; CO="xco yco"; CLO="0 0"; CHI="1 1"
    B8="8 8"; B16="16 16"; B24="24 24"
else
    ALLPER="1 1 1"; CO="xco yco zco"; CLO="0 0 0"; CHI="1 1 1"
    B8="8 8 8"; B16="16 16 16"; B24="24 24 24"
fi

echo "exe    : $JS"
echo "workdir: $WORK"
echo "dim    : $DIM   ranks: $RANKS"
echo

FAIL=0
pass () { echo "  PASS  $1"; }
fail () { echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

# R^2 for a named coordinate set at a given coarsening.
r2 () { awk -v k="$1" -v c="$2" '$1==k && $2==c {print $3}' ; }

# Table data rows only. The structural test is the point: a leading-whitespace
# match alone also catches the run's own header lines, which carry a different
# number of fields and so shift every column in a field-by-field comparison.
# A data row is exactly "name coarsen nOcc R^2 1-R^2" with the middle two
# integral, and the exe now prefixes everything else with '#', so this is belt
# and braces.
datarows () {
    awk '/^  [a-z]/ && NF==5 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/'
}
# the row fields we care about: name coarsen r2
rows  () { datarows | awk '{print $1, $2, $4}' ; }
# with the occupancy, for the rank-invariance comparison
rows4 () { datarows | awk '{print $1, $2, $3, $4}' ; }

echo "generating hierarchies..."
$MPIRUN -n 1 "$MAKEPLT" n_cell=96 nlevels=2 refine_shape=random refine_anchor=lo \
    refine_blk=8 refine_frac=1.0 refine_fill=0.6 outfile=js_big > gen_big.log 2>&1
$MPIRUN -n 1 "$MAKEPLT" n_cell=48 nlevels=3 refine_shape=random box \
    refine_anchor=lo lo refine_blk=8 refine_frac=1.0 refine_fill=0.6 \
    outfile=js_jump > gen_jump.log 2>&1
for d in js_big js_jump; do
    [ -d "$d" ] || { fail "generator failed for $d"; }
done
grep -E "^   level" gen_big.log | sed 's/^/  [big]/'
echo

# --------------------------------------------------------------------------
echo "=== 1. the volume measure =============================================="
# --------------------------------------------------------------------------
for d in js_big js_jump; do
    [ -d "$d" ] || continue
    out=$($MPIRUN -n 1 "$JS" infile=$d is_per="$ALLPER" coords="$CO" \
            clo="$CLO" chi="$CHI" cbins="$B8" resp=dist 2>&1)
    vin=$(echo "$out"  | sed -n 's/^# volume in range : *\([^ ]*\).*/\1/p')
    vout=$(echo "$out" | sed -n 's/^# volume dropped  : *\([^ ]*\).*/\1/p')
    if [ "$vin" = "1" ] && [ "$vout" = "0" ]; then
        pass "$d : volume in range exactly 1, nothing dropped"
    else
        fail "$d : volume in range $vin, dropped $vout (want 1 and 0)"
    fi
done
echo

# --------------------------------------------------------------------------
echo "=== 2. negative control: an independent response ======================="
# --------------------------------------------------------------------------
out=$($MPIRUN -n 1 "$JS" infile=js_big is_per="$ALLPER" coords=xco clo=0 chi=1 \
        cbins=96 resp=yco 2>&1)
bad=$(echo "$out" | rows | awk '$1=="xco" && ($3+0) != 0 {print}' | wc -l | tr -d ' ')
ncoar=$(echo "$out" | rows | awk '$1=="xco"' | wc -l | tr -d ' ')
if [ "$ncoar" -lt 3 ]; then
    fail "negative control : only $ncoar coarsenings reported"
elif [ "$bad" = "0" ]; then
    pass "negative control : R^2 exactly 0 at all $ncoar coarsenings"
else
    fail "negative control : $bad of $ncoar coarsenings gave nonzero R^2"
    echo "$out" | rows | awk '$1=="xco"' | sed 's/^/        /'
fi
echo

# --------------------------------------------------------------------------
echo "=== 3. second-order convergence on an exact function ==================="
# --------------------------------------------------------------------------
KEY=$(echo "$CO" | tr ' ' ',')
out=$($MPIRUN -n 1 "$JS" infile=js_big is_per="$ALLPER" coords="$CO" \
        clo="$CLO" chi="$CHI" cbins="$B24" resp=dist 2>&1)
echo "$out" | rows | awk -v k="$KEY" '$1==k {print $2, 1.0-$3}' > conv.txt
nconv=$(wc -l < conv.txt | tr -d ' ')
if [ "$nconv" -lt 3 ]; then
    fail "convergence : only $nconv points"
else
    # The O(h^2) statement is an expansion in a small residual, so the ratio is
    # only tested where the residual still is small. RESMAX is that cut, not a
    # tolerance on the answer: above it the leading term no longer dominates
    # and a ratio below 4 is the correct behaviour, not a failure. Pairs above
    # the cut are reported and not judged.
    RESMAX=0.2
    ok=1
    nchk=0
    prev=""
    msg=""
    while read -r co e; do
        if [ -n "$prev" ]; then
            ratio=$(awk -v a="$e" -v b="$prev" 'BEGIN{printf "%.3f", a/b}')
            inwin=$(awk -v e="$e" -v m="$RESMAX" 'BEGIN{print (e<m) ? 1 : 0}')
            good=$(awk -v r="$ratio" 'BEGIN{print (r>3.7 && r<4.3) ? 1 : 0}')
            if [ "$inwin" = "1" ]; then
                nchk=$((nchk+1))
                msg="$msg $ratio"
                if [ "$good" != "1" ]; then ok=0; fi
            else
                msg="$msg ($ratio)"
            fi
        fi
        prev=$e
    done < conv.txt
    if [ "$nchk" -lt 2 ]; then ok=0; fi
    if [ "$ok" = "1" ]; then
        pass "convergence : 1-R^2 ratios$msg, $nchk tested against 4.0 (bracketed ones are above the small-residual cut)"
    else
        fail "convergence : 1-R^2 ratios$msg, $nchk tested, want 4.0 +/- 0.3"
    fi
fi
echo

# --------------------------------------------------------------------------
echo "=== 3b. resonant binning: a second, harder null ======================="
# --------------------------------------------------------------------------
# cells is cos(2 pi k x/L) summed over directions, k = 4 by default. At 8 bins
# per side the sample points in every bin are centred on a zero crossing of the
# cosine, so each bin mean vanishes and R^2 must be zero -- at 8 bins and every
# coarsening of it.
#
# This is a different test from the yco control, not a repeat of it. There the
# response was linear and the conditional mean exactly constant, so zero came
# out of a cancellation of identical terms. Here the response swings through its
# full range inside every bin and the within-bin sum of squares is large: zero
# has to emerge as a small difference between two big numbers, which is exactly
# the path a one-pass S2 - S1^2/V formulation gets wrong.
#
# Gated on a SINGLE-level hierarchy, where it is exact. On a multi-level one it
# is not, and the reason is worth recording rather than tolerancing away: the
# volume-weighted mean uses cell-centre values, which is the midpoint rule, and
# the midpoint rule is exact for a linear function but carries an O(dx^2) error
# for a cosine. Within one level that error follows the cosine and cancels over
# a half period; across a coarse-fine interface the two levels carry different
# dx and it does not. The residual below is therefore a measure of the
# hierarchy's discretisation, not of this tool, and it falls as the grid is
# refined. Measured here: exactly 0 at one level, and 1e-11 to 3e-8 at two.
# --------------------------------------------------------------------------
$MPIRUN -n 1 "$MAKEPLT" n_cell=96 nlevels=1 outfile=js_one > gen_one.log 2>&1
if [ ! -d js_one ]; then
    fail "resonant null : generator failed for js_one"
else
    out3b=$($MPIRUN -n 1 "$JS" infile=js_one is_per="$ALLPER" coords="$CO" \
            clo="$CLO" chi="$CHI" cbins="$B8" resp=cells 2>&1)
    nrow=$(echo "$out3b" | rows | grep -c .)
    nz=$(echo "$out3b" | rows | awk '($3+0) != 0' | grep -c .)
    if [ "$nrow" -lt 2 ]; then
        fail "resonant null : only $nrow rows"
    elif [ "$nz" = "0" ]; then
        pass "resonant null : R^2 exactly 0 on all $nrow rows, one level"
    else
        fail "resonant null : $nz of $nrow rows nonzero at one level"
        echo "$out3b" | rows | sed 's/^/        /'
    fi
fi

out3c=$($MPIRUN -n 1 "$JS" infile=js_big is_per="$ALLPER" coords="$CO" \
        clo="$CLO" chi="$CHI" cbins="$B8" resp=cells 2>&1)
worst=$(echo "$out3c" | rows | awk '{v=($3+0); if(v<0)v=-v; if(v>m)m=v} END{printf "%.3e", m}')
okz=$(awk -v w="$worst" 'BEGIN{print (w<1e-6)?1:0}')
if [ "$okz" = "1" ]; then
    pass "resonant floor : $worst at two levels, the mixed-dx midpoint error"
else
    fail "resonant floor : $worst at two levels, want < 1e-6"
fi
echo

# --------------------------------------------------------------------------
echo "=== 4. marginalisation and coarsening ================================="
# --------------------------------------------------------------------------
echo "$out" | rows | awk '$1=="xco" {print $2, $3}' > marg_from_full.txt
$MPIRUN -n 1 "$JS" infile=js_big is_per="$ALLPER" coords=xco clo=0 chi=1 \
        cbins=24 resp=dist 2>&1 | rows | awk '$1=="xco" {print $2, $3}' > marg_direct.txt
if cmp -s marg_from_full.txt marg_direct.txt; then
    pass "marginalisation : {xco} from the full table is identical to a direct run"
else
    fail "marginalisation : differs from a direct run"
    diff marg_from_full.txt marg_direct.txt | sed 's/^/        /'
fi

a=$(echo "$out" | rows | awk '$1=="xco" && $2==2 {print $3}')
b=$($MPIRUN -n 1 "$JS" infile=js_big is_per="$ALLPER" coords=xco clo=0 chi=1 \
        cbins=12 resp=dist 2>&1 | rows | awk '$1=="xco" && $2==1 {print $3}')
if [ "$a" = "$b" ] && [ -n "$a" ]; then
    pass "coarsening : 24 bins at coarsen 2 is identical to 12 bins ($a)"
else
    fail "coarsening : 24/coarsen2 gave $a, 12 bins gave $b"
fi
echo

# --------------------------------------------------------------------------
echo "=== 5. permutation symmetry of an isotropic field ====================="
# --------------------------------------------------------------------------
sing=$(echo "$out" | rows | awk '$2==1 && $1 !~ /,/ {print $3}')
pair=$(echo "$out" | rows | awk '$2==1 && gsub(/,/,",")==1 {print $3}')
spread () { awk '{v[NR]=$1} END{mn=v[1];mx=v[1];for(i=1;i<=NR;i++){if(v[i]<mn)mn=v[i];if(v[i]>mx)mx=v[i]}; printf "%.3e", mx-mn}'; }
ds=$(echo "$sing" | spread)
dp=$(echo "$pair" | spread)
ns=$(echo "$sing" | grep -c .)
np=$(echo "$pair" | grep -c .)
oks=$(awk -v d="$ds" 'BEGIN{print (d<1e-4)?1:0}')
okp=$(awk -v d="$dp" 'BEGIN{print (d<1e-4)?1:0}')
if [ "$ns" -lt "$DIM" ]; then
    fail "symmetry : only $ns single-coordinate rows, expected $DIM"
elif [ "$oks" = "1" ] && [ "$okp" = "1" ]; then
    pass "symmetry : singles spread $ds, pairs spread $dp"
else
    fail "symmetry : singles spread $ds, pairs spread $dp (want < 1e-4)"
    echo "$sing" | sed 's/^/        single /'
    echo "$pair" | sed 's/^/        pair   /'
fi
echo

# --------------------------------------------------------------------------
echo "=== 6. rank invariance ================================================"
# --------------------------------------------------------------------------
# Two different claims, gated differently, because only one of them can be
# exact.
#
# The volume of a bin is an integer count per level times that level's cell
# volume, reduced as integers and multiplied in a fixed order. That IS bit
# identical at any rank count, and the bin table's coordinate, volume and cell
# columns are compared byte for byte. If the partition ever changed which cells
# land in which bin, this is where it would show, with no tolerance to hide in.
#
# R^2 is a ratio of floating-point sums whose cross-rank reduction is a plain
# MPI sum, so it agrees to a few ulp and not exactly. Within a rank the sums are
# compensated, so the floor is set by the cross-rank step: a few times 1e-15 on
# R^2. The gate is 1e-12 absolute, three orders above that floor and far below
# anything a real distribution dependence could produce -- such a defect moves
# leading digits, not the fifteenth. The worst difference actually seen is
# printed on success rather than assumed, so the floor is measured every run.
#
# This matters most exactly where it looks worst. Where R^2 is itself near zero,
# as it is for an independent response at the coarsest aggregation, those few
# ulp ARE its leading digits: 1.2130630e-11 against 1.2130297e-11 is two
# spellings of nothing, differing by 3e-16.
# --------------------------------------------------------------------------
ref=""
refbin=""
for n in $RANKS; do
    $MPIRUN -n "$n" "$JS" infile=js_jump is_per="$ALLPER" coords="$CO" \
        clo="$CLO" chi="$CHI" cbins="$B16" resp="dist cells" \
        outfile=bin_n$n.dat 2>&1 | rows4 > inv_n$n.txt
    cut -d' ' -f1-$((DIM+2)) bin_n$n.dat > binexact_n$n.txt
    if [ -z "$ref" ]; then ref="inv_n$n.txt"; refbin="binexact_n$n.txt"; continue; fi

    if cmp -s "$refbin" "binexact_n$n.txt"; then
        pass "$n ranks : bin volumes and cell counts byte-identical"
    else
        fail "$n ranks : bin volumes or cell counts differ"
        diff "$refbin" "binexact_n$n.txt" | head -8 | sed 's/^/        /'
    fi

    paste -d' ' "$ref" "inv_n$n.txt" > pair.txt
    read -r worst where bad <<EOF2
$(awk '{
        if ($1!=$5 || $2!=$6 || $3!=$7) { bad=bad+1 }
        d=$4-$8; if (d<0) d=-d;
        if (d>m) { m=d; w=$1"/coarsen"$2 }
      } END { printf "%.3e %s %d\n", m, (w==""?"-":w), bad+0 }' pair.txt)
EOF2
    if [ "${bad:-0}" != "0" ]; then
        fail "$n ranks : $bad row(s) differ in key or occupancy, not just value"
        # Print them. The first version of this branch reported the count and
        # nothing else, and the rows in question turned out to be identical --
        # the mismatch was a column shift in this comparison, not in the data.
        # A test that withholds its evidence costs more than it saves.
        awk '$1!=$5 || $2!=$6 || $3!=$7 {print "        " $0}' pair.txt | head -6
    else
        okd=$(awk -v d="$worst" 'BEGIN{print (d<1e-12)?1:0}')
        if [ "$okd" = "1" ]; then
            pass "$n ranks : occupancy identical, worst dR^2 $worst at $where"
        else
            fail "$n ranks : worst dR^2 $worst at $where, want < 1e-12"
            paste -d' ' "$ref" "inv_n$n.txt" \
              | awk '{d=$4-$8; if(d<0)d=-d; if(d>1e-12) print}' | head -6 \
              | sed 's/^/        /'
        fi
    fi
done
echo

echo "======================================================================="
if [ "$FAIL" -eq 0 ]; then
    echo "all joint-stats tests pass"
else
    echo "$FAIL failure(s)"
fi
exit "$FAIL"
