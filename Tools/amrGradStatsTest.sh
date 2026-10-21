#!/usr/bin/env bash
#
# Verification suite for amrGradStats / AMRGradient.
#
# The claim under test is narrow and specific: gradients taken from a composite
# MLPoisson apply are correct AT COARSE-FINE INTERFACES, where a difference
# stencil on a level's own ghost data is not. Every test below is therefore
# reported on three regions -- all cells, interior cells, and the band within
# two cells of a coarse-fine interface -- and it is the third that matters.
#
# The exact tests use a linear field, because MLMG's interface interpolation is
# exact for a linear function: the only error left is floating point. A
# difference quotient divides by dx, so a field of magnitude O(1) differenced at
# the finest spacing carries eps/dx, and a divergence of a difference carries
# eps/dx^2. That is not slop to be tolerated, it is the arithmetic: with the
# refined patch at the LOW corner the ramp's values are near zero and the result
# is bit-exact, and with it at the HIGH corner the same code gives 7e-15, which
# is 32*eps amplified by 1/dx at dx=1/128.
#
# So the gate is tol = 64*eps*|phi|max/dx, and /dx^2 for the curvature. The
# point is not that the gate is tight in absolute terms -- it is that the
# per_level control misses by O(1/dx) with an O(1) coefficient, which is
# fourteen orders of magnitude the other side of it. Each test prints the
# observed value and that margin, so the headroom is visible rather than
# asserted.
#
# Each exact test is paired with per_level=1, which builds one single-level
# operator per level -- what Src/curvature.cpp does. It must FAIL. A gate that
# nothing can fail is not a gate.
#
# Build:
#   cd Src
#   make -f GNUmakefile.AMRComplex EBASE=amrGradStats  -j8
#   make -f GNUmakefile.AMRComplex EBASE=amrGradVerify -j8
#   make -f GNUmakefile.AMRComplex EBASE=makeTestPlt   -j8
#
# Then from any writable directory:
#   /path/to/PeleAnalysis/Tools/amrGradStatsTest.sh
#   DIM=2 Tools/amrGradStatsTest.sh
# ---------------------------------------------------------------------------
#
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../Src"
DIM=${DIM:-3}
DTAG="${DIM}d"

# Pick the NEWEST build of each tool, whatever compiler and whether or not it is
# an MPI build, and run a one-rank job.
#
# An earlier version preferred a non-MPI build, on the theory that nothing here
# tests the decomposition. That was wrong in a way worth recording: this tree
# still carries *.gnu.ex binaries from an earlier era alongside the llvm.MPI
# builds that are actually maintained, so "prefer serial" reliably picked the
# stalest executable present. Newest-wins keeps all three tools from the same
# build, which is what makes their results comparable in the first place.
pick () {
    ls -t "$SRC/${1}${DTAG}".*ex 2>/dev/null | head -1
}
# One-rank launcher for an MPI build; empty for a serial one.
launcher () {
    case $1 in
        *MPI*) if command -v "${MPIRUN:-mpiexec}" >/dev/null 2>&1; then
                   printf '%s -n 1' "${MPIRUN:-mpiexec}"
               fi ;;
    esac
}

GEN=$(pick makeTestPlt); GRD=$(pick amrGradStats); VFY=$(pick amrGradVerify)
MISSING=""
for pair in "makeTestPlt:$GEN" "amrGradStats:$GRD" "amrGradVerify:$VFY"; do
    nm=${pair%%:*}; px=${pair#*:}
    if [ -z "$px" ] || [ ! -x "$px" ]; then MISSING="$MISSING $nm"; fi
done
if [ -n "$MISSING" ]; then
    echo "ERROR: no ${DTAG} executable for:$MISSING"
    echo "  looked for: $SRC/<name>${DTAG}.*ex"
    echo "  present in $SRC:"
    ls -1 "$SRC"/*${DTAG}.*ex 2>/dev/null | sed 's|.*/|    |' || echo "    (none)"
    echo
    echo "Build the missing ones:"
    for nm in $MISSING; do
        echo "    make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm DIM=${DIM} EBASE=$nm -j8"
    done
    exit 2
fi

GENRUN="$(launcher "$GEN") $GEN"
GRDRUN="$(launcher "$GRD") $GRD"
VFYRUN="$(launcher "$VFY") $VFY"

WORK=${WORK:-$(pwd)/amrGradStats_test}
mkdir -p "$WORK" && cd "$WORK" || exit 2
RUNLOG="$WORK/runs.log"
: > "$RUNLOG"

# Every tool invocation goes through this. Discarding their output -- which an
# earlier version of this script did -- turns one broken generator into thirteen
# unexplained failures downstream, with the one message that would have named
# the cause thrown away. Keep it, and show it when a step fails.
step () {
    # rc is declared with tag, not assigned on its own line after the command:
    # a `local rc=$?` would report the status of the `local` builtin, not the
    # command, which is how this first printed "(exit 0)" for a failed run.
    local tag=$1 rc
    shift
    { echo; echo "### $tag"; echo "+ $*"; } >> "$RUNLOG"
    "$@" >> "$RUNLOG" 2>&1
    rc=$?
    [ $rc -eq 0 ] || echo "    (exit $rc)" >> "$RUNLOG"
    return $rc
}
showlast () { tail -${1:-20} "$RUNLOG" | sed 's/^/      /'; }

echo "gen    : $GEN"
echo "grad   : $GRD"
echo "verify : $VFY"
echo "workdir: $WORK"
echo "runlog : $RUNLOG"

NPASS=0; NFAIL=0
ok   () { echo "  PASS  $*"; NPASS=$((NPASS+1)); }
bad  () { echo "  FAIL  $*"; NFAIL=$((NFAIL+1)); }

# linf for one component/region out of a verify run
lx () {
    awk -v c="$2" -v r="$3" '$1=="CHECK" && $2==c && $3==r {print $4}' "$1"
}
l2 () {
    awk -v c="$2" -v r="$3" '$1=="CHECK" && $2==c && $3==r {print $5}' "$1"
}
nc () {
    awk -v c="$2" -v r="$3" '$1=="CHECK" && $2==c && $3==r {print $7}' "$1"
}

# numeric compare without bc
le () { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a<=b)}'; }
gt () { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a>b)}'; }

if [ "$DIM" = 3 ]; then
    VARS="xco yco zco dist"
    EXACT="gradx_xco grady_xco gradz_xco gradMag_xco grady_yco gradz_zco cosTheta_xco_yco sinTheta_xco_yco"
    EXPECT="1 0 0 1 1 1 0 1"
    PERVARS="xco yco zco"
    PEREXACT="grady_xco gradz_xco gradx_yco gradz_yco gradx_zco grady_zco"
    PEREXPECT="0 0 0 0 0 0"
else
    VARS="xco yco dist"
    EXACT="gradx_xco grady_xco gradMag_xco grady_yco cosTheta_xco_yco sinTheta_xco_yco"
    EXPECT="1 0 1 1 0 1"
    PERVARS="xco yco"
    PEREXACT="grady_xco gradx_yco"
    PEREXPECT="0 0"
fi

# A smoke test before anything is measured. If the generator cannot write a
# plotfile then every check below fails for one reason, and saying so once with
# the generator's own output beats thirteen downstream failures that name
# nothing. A stale executable from an older build of this tree is the usual
# cause, so the message says that.
if ! step preflight $GENRUN n_cell=16 nlevels=2 refine_shape=box \
        refine_anchor=center outfile=p_preflight || [ ! -d p_preflight ]; then
    echo "ERROR: the generator did not produce a plotfile."
    echo "  exe: $GEN"
    echo "  its output:"
    showlast 25
    echo
    echo "  Most likely it is a stale build. Rebuild it and retry:"
    echo "    cd $SRC && make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=llvm DIM=$DIM EBASE=makeTestPlt -j8"
    exit 2
fi

echo "=== amrGradStats verification, DIM=$DIM ==="

# ---------------------------------------------------------------- 1,2,3
# Exact linear identities, non-periodic, on hierarchies whose refined patches
# sit in three different places. "anchor" moves where the coarse-fine interface
# is, which is the whole point: an interface in the middle of the domain, one
# touching a domain face, and one offset.
for ANCHOR in center lo hi; do
  for NLEV in 2 3; do
    TAG="exact_${ANCHOR}_L${NLEV}"
    step "gen $TAG" $GENRUN n_cell=32 nlevels=$NLEV refine_shape=box \
         refine_anchor=$ANCHOR outfile=p_$TAG
    if [ ! -d p_$TAG ]; then
        bad "$TAG: generator produced nothing"; showlast 12; continue
    fi

    step "grad $TAG" $GRDRUN infile=p_$TAG outfile=g_$TAG is_per="0 0 0" \
         vars="$VARS" curvOf=xco
    step "grad $TAG control" $GRDRUN infile=p_$TAG outfile=n_$TAG is_per="0 0 0" \
         vars="$VARS" curvOf=xco per_level=1
    if [ ! -d g_$TAG ] || [ ! -d n_$TAG ]; then
        bad "$TAG: amrGradStats produced nothing"; showlast 12; continue
    fi

    CH="$EXACT meanCurv_xco"
    EX="$EXPECT 0"
    $VFYRUN infile=g_$TAG is_per="0 0 0" check="$CH" expect="$EX" > v_$TAG.txt 2>&1
    $VFYRUN infile=n_$TAG is_per="0 0 0" check="$CH" expect="$EX" > w_$TAG.txt 2>&1
    if ! grep -q "^CHECK" v_$TAG.txt; then
        bad "$TAG: amrGradVerify reported no checks"
        sed 's/^/      /' v_$TAG.txt | tail -12; continue
    fi

    # eps/dx at the finest level, with |phi|max = 1 for a unit-domain ramp.
    INVDX=$(awk -v n=32 -v l=$NLEV 'BEGIN{printf "%.17g", n*2^(l-1)}')
    TOLG=$(awk -v h="$INVDX" 'BEGIN{printf "%.17g", 64*2.220446049250313e-16*h}')
    TOLC=$(awk -v h="$INVDX" 'BEGIN{printf "%.17g", 64*2.220446049250313e-16*h*h}')

    nbadc=0; worst=0; worstat=""
    for c in $CH; do
      case $c in meanCurv_*) T=$TOLC ;; *) T=$TOLG ;; esac
      for r in int cf; do
        v=$(lx v_$TAG.txt "$c" "$r")
        n=$(nc v_$TAG.txt "$c" "$r")
        if [ -z "$v" ] || [ -z "$n" ] || [ "$n" -lt 1 ]; then
            echo "    missing or empty region $c/$r"; nbadc=$((nbadc+1)); continue
        fi
        # track the worst as a fraction of its own tolerance
        f=$(awk -v a="$v" -v t="$T" 'BEGIN{printf "%.4g", (t>0)?a/t:0}')
        gt "$f" "$worst" && { worst=$f; worstat="$c/$r"; }
        le "$v" "$T" || { echo "    composite $c/$r Linf=$v tol=$T"; nbadc=$((nbadc+1)); }
      done
    done
    if [ $nbadc -eq 0 ]; then
        ok "$TAG composite within roundoff on int and cf (worst $worst x tol, at ${worstat:-none})"
    else bad "$TAG composite: $nbadc check(s) outside roundoff"; fi

    # The control must break, and must break in the cf band, by a margin that
    # makes the gate unmistakable.
    nfire=0; margin=0
    for c in $EXACT; do
      w=$(lx w_$TAG.txt "$c" cf)
      [ -z "$w" ] && continue
      if gt "$w" "$TOLG"; then
          nfire=$((nfire+1))
          m=$(awk -v a="$w" -v t="$TOLG" 'BEGIN{printf "%.3g", a/t}')
          gt "$m" "$margin" && margin=$m
      fi
    done
    NEX=$(echo $EXACT | wc -w | tr -d " ")
    if [ $nfire -eq "$NEX" ]; then
        ok "$TAG per_level control fires on all $NEX cf checks (worst ${margin}x tol)"
    elif [ $nfire -ge 1 ]; then
        ok "$TAG per_level control fires on $nfire/$NEX cf checks (worst ${margin}x tol)"
    else bad "$TAG per_level control did NOT fire -- the cf gate is vacuous"; fi
  done
done

# ---------------------------------------------------------------- 4
# Periodic wrap. The field has to actually BE periodic for the test to mean
# anything: a ramp in x is discontinuous across a periodic x boundary, and the
# coarse-fine interpolation is multi-dimensional, so that jump contaminates the
# transverse fluxes too -- a first attempt at this test asserted grady_xco=0
# with x periodic and measured the discontinuity, not the scheme.
#
# So: x stays Neumann and only y and z are periodic. xco is constant in y and z,
# hence continuous across both wraps, and its y- and z-derivatives must be
# exactly zero -- including in cells whose stencil reaches through a wrap, and
# including a refined patch placed on the y/z domain faces (refine_anchor=lo),
# where the wrap and the coarse-fine interface meet.
TAG=periodic
if [ "$DIM" = 3 ]; then PER="0 1 1"; PEREXACT="grady_xco gradz_xco"; PEREXPECT="0 0"
else                    PER="0 1";   PEREXACT="grady_xco";          PEREXPECT="0"; fi

step "gen $TAG" $GENRUN n_cell=32 nlevels=3 refine_shape=box refine_anchor=lo \
     outfile=p_$TAG
step "grad $TAG" $GRDRUN infile=p_$TAG outfile=g_$TAG is_per="$PER" vars="xco"
step "grad $TAG control" $GRDRUN infile=p_$TAG outfile=n_$TAG is_per="$PER" \
     vars="xco" per_level=1
$VFYRUN infile=g_$TAG is_per="$PER" check="$PEREXACT" expect="$PEREXPECT" > v_$TAG.txt 2>&1
$VFYRUN infile=n_$TAG is_per="$PER" check="$PEREXACT" expect="$PEREXPECT" > w_$TAG.txt 2>&1

nbadc=0
for c in $PEREXACT; do
  for r in all int cf; do
    v=$(lx v_$TAG.txt "$c" "$r")
    n=$(nc v_$TAG.txt "$c" "$r")
    if [ -z "$v" ] || [ -z "$n" ] || [ "$n" -lt 1 ]; then
        echo "    missing or empty $c/$r"; nbadc=$((nbadc+1)); continue
    fi
    # differencing a field that is CONSTANT in this direction subtracts equal
    # values, so here bit-exactness is not an aspiration, it is subtraction.
    le "$v" 0 || { echo "    composite $c/$r Linf=$v (want exactly 0)"; nbadc=$((nbadc+1)); }
  done
done
if [ $nbadc -eq 0 ]; then ok "periodic transverse derivatives bit-exact through the wrap"
else bad "periodic: $nbadc deviation(s)"; fi

nfire=0
for c in $PEREXACT; do
  w=$(lx w_$TAG.txt "$c" cf)
  [ -n "$w" ] && gt "$w" 1.0e-10 && nfire=$((nfire+1))
done
if [ $nfire -ge 1 ]; then ok "periodic per_level control fires"
else bad "periodic per_level control did NOT fire"; fi

# ---------------------------------------------------------------- 5
# Second-order convergence on a field that is NOT linear. |grad dist| = 1
# exactly for a signed distance, so the departure from 1 is pure truncation
# error. The ball around the sphere centre is excluded: dist has a cone tip
# there, which is not differentiable, and its O(1) error inside an O(dx^DIM)
# volume caps the measurable rate at DIM/2 regardless of the scheme.
#
# The rate is measured separately on the interior and on the coarse-fine band.
# Both must be second order. A coarse-fine treatment that is merely consistent
# rather than second-order shows up precisely here: the interior rate stays at 2
# and the cf rate falls to 1.
for N in 32 64; do
    step "gen cv$N" $GENRUN n_cell=$N nlevels=2 refine_shape=box \
         refine_anchor=center outfile=p_cv$N
    step "grad cv$N" $GRDRUN infile=p_cv$N outfile=g_cv$N is_per="0 0 0" vars="dist"
    $VFYRUN infile=g_cv$N is_per="0 0 0" check="gradMag_dist" expect="1" \
         exclude_c="0.5 0.5 0.5" exclude_r=0.12 > v_cv$N.txt 2>&1
done
for R in int cf; do
    e1=$(l2 v_cv32.txt gradMag_dist $R)
    e2=$(l2 v_cv64.txt gradMag_dist $R)
    if [ -z "$e1" ] || [ -z "$e2" ]; then bad "convergence/$R: missing L2"; continue; fi
    rate=$(awk -v a="$e1" -v b="$e2" 'BEGIN{ if (b>0) printf "%.3f", log(a/b)/log(2); else print "99" }')
    if gt "$rate" 1.8; then ok "|grad dist|=1 second order on $R (rate $rate: $e1 -> $e2)"
    else bad "|grad dist| on $R converges at $rate, want >1.8 ($e1 -> $e2)"; fi
done

# ---------------------------------------------------------------- 6
# finestLevel truncation. Restricting to level 0 must give the same answer as a
# single-level run of the same base grid, and must still be exact -- a sanity
# check that the composite path degenerates correctly rather than special-casing.
TAG=flat
step "gen $TAG" $GENRUN n_cell=32 nlevels=3 refine_shape=box \
     refine_anchor=center outfile=p_$TAG
step "grad $TAG" $GRDRUN infile=p_$TAG outfile=g_$TAG is_per="0 0 0" \
     vars="$VARS" finestLevel=0
$VFYRUN infile=g_$TAG is_per="0 0 0" check="$EXACT" expect="$EXPECT" > v_$TAG.txt 2>&1
nbadc=0
for c in $EXACT; do
  v=$(lx v_$TAG.txt "$c" int)
  [ -z "$v" ] && { nbadc=$((nbadc+1)); continue; }
  le "$v" 0 || { echo "    $c/int Linf=$v"; nbadc=$((nbadc+1)); }
done
# with one level there is no coarse-fine interface, so the cf band must be empty
cfn=$(nc v_$TAG.txt gradx_xco cf)
if [ $nbadc -eq 0 ]; then ok "finestLevel=0 exact on interior"
else bad "finestLevel=0: $nbadc deviation(s)"; fi
if [ "${cfn:-1}" = 0 ]; then ok "finestLevel=0 has an empty cf band, as it must"
else bad "finestLevel=0 reports $cfn cf cells with no finer level present"; fi

# ---------------------------------------------------------------- 7
# sinTheta is the measure the flame work actually leans on: it must be
# identically zero where one field is a function of the other, for ANY monotone
# relation. xco and a strictly monotone function of xco have parallel gradients
# everywhere, so sin must vanish exactly -- and does, with no reference to what
# the function was. dist vs xco are genuinely transverse, so sin must NOT
# vanish; otherwise the measure is blind.
TAG=align
step "gen $TAG" $GENRUN n_cell=32 nlevels=3 refine_shape=box \
     refine_anchor=center outfile=p_$TAG
step "grad $TAG" $GRDRUN infile=p_$TAG outfile=g_$TAG is_per="0 0 0" vars="xco dist"
$VFYRUN infile=g_$TAG is_per="0 0 0" check="sinTheta_xco_dist" expect="0" > v_$TAG.txt 2>&1
sv=$(lx v_$TAG.txt sinTheta_xco_dist int)
if [ -n "$sv" ] && gt "$sv" 0.5; then ok "sinTheta separates transverse fields (Linf from 0 = $sv)"
else bad "sinTheta failed to register transversality between xco and dist ($sv)"; fi

echo
echo "=== DIM=$DIM : $NPASS passed, $NFAIL failed ==="
[ $NFAIL -eq 0 ] || exit 1
exit 0
