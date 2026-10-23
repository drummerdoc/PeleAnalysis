# Shared executable selection for the test suites and the pipeline.
# Sourced, not run:  . "$HERE/amrPickExe.sh"
#
# Why this exists.
#
# The suites each grew their own way of finding the binary to test, and the
# ways disagreed. Four of them took the first match of a glob, which is
# ALPHABETICAL order, so with both a gnu and an llvm build present in Src the
# gnu one always won however old it was. Three took the newest by mtime. On a
# machine where the harness builds llvm and a two-day-old gnu build is still
# lying in Src, the first group silently tested the old binary: the shape
# suite reported three wall-closure failures against an executable that
# predated the check, while the executable the harness had just built passed.
#
# Nothing about that is visible in the failure. A stale binary gives wrong
# answers confidently, and the suite blames the mesh.
#
# So: one rule, stated once, and it refuses rather than guesses.
#
#   1. If COMP is set and $SRC/<base><dtag>.$COMP.*ex exists, that is the one.
#      The harness exports COMP, so a suite it launches tests exactly what it
#      just built, whatever else is in the directory.
#   2. Otherwise the newest match by mtime, which is the better default when a
#      person is building one thing by hand.
#   3. Either way, if the chosen executable is older than the newest source
#      file in Src, that is reported and treated as an error. Every defect of
#      this kind is a stale binary against updated sources, and that is
#      exactly the comparison that catches it.
#
# pick_exe <base> [dtag] [mpionly]
#   base     driver name, e.g. amrComplexVerify
#   dtag     defaults to $DTAG, else ${DIM}d
#   mpionly  non-empty to require an MPI build
# Prints the path. Returns 1 and prints nothing on failure.

pick_exe () {
    local base=$1
    local dtag=${2:-${DTAG:-${DIM:-3}d}}
    local mpionly=${3:-}
    local suffix="*ex"
    [ -n "$mpionly" ] && suffix="*MPI.ex"

    local chosen=""
    if [ -n "${COMP:-}" ]; then
        chosen=$(ls -t "$SRC/${base}${dtag}.${COMP}."$suffix 2>/dev/null | head -1)
    fi
    [ -n "$chosen" ] || chosen=$(ls -t "$SRC/${base}${dtag}."$suffix 2>/dev/null | head -1)
    if [ -z "$chosen" ] || [ ! -x "$chosen" ]; then
        # Say which it is. "No usable build" covers two quite different
        # situations -- nothing was ever built, or what was built is out of
        # date -- and they are fixed the same way but diagnosed differently,
        # so the message should not make the reader guess.
        echo "ERROR: no ${dtag} executable matching ${base}${dtag}.${COMP:+${COMP}.}${suffix} in $SRC" >&2
        if ls "$SRC/${base}${dtag}."*ex >/dev/null 2>&1; then
            echo "  present but not matching:" >&2
            ls -t "$SRC/${base}${dtag}."*ex 2>/dev/null | sed 's|.*/|    |' >&2
        else
            echo "  nothing of that name at all; the build did not run or did not finish" >&2
        fi
        return 1
    fi

    # Staleness, reported and not enforced.
    #
    # This started as a refusal and has now produced two false stops and no
    # true ones. First it compared against the newest of all forty-odd
    # sources in Src, so a new driver nothing else links stopped three
    # suites. Narrowing it to the makefile's own list left AMRGradient in
    # every tool's dependencies when only amrGradStats links it, and stopped
    # the same suites again on executables that were perfectly current.
    #
    # The information is worth having; the refusal is not. A check that stops
    # a run it should not understand stops being read. So it warns, loudly,
    # with both timestamps, and lets the suite proceed -- and if the binary
    # really is stale the suite's own failures now come with this line above
    # them, which is all the refusal was ever for.
    #
    # The list mirrors GNUmakefile.AMRComplex: five library sources always,
    # AMRGradient only for the target that links it, and the driver itself.
    # Keep them in step.
    local deps="AMRComplex AMRContour AMRSegment AMRParallel AMRJointStats"
    case $base in amrGradStats) deps="$deps AMRGradient" ;; esac
    local files="" d
    for d in $deps; do
        [ -f "$SRC/$d.cpp" ] && files="$files $SRC/$d.cpp"
        [ -f "$SRC/$d.H" ]   && files="$files $SRC/$d.H"
    done
    [ -f "$SRC/${base}.cpp" ] && files="$files $SRC/${base}.cpp"

    local newest
    newest=$(ls -t $files 2>/dev/null | head -1)
    if [ -n "$newest" ] && [ "$newest" -nt "$chosen" ]; then
        echo "WARNING: $(basename "$chosen") is older than $(basename "$newest"), which it is built from." >&2
        echo "           executable: $(date -r "$chosen" '+%Y-%m-%d %H:%M' 2>/dev/null)" >&2
        echo "           source    : $(date -r "$newest" '+%Y-%m-%d %H:%M' 2>/dev/null)" >&2
        echo "         Testing it anyway. If what follows fails, rebuild first:" >&2
        echo "         make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=${COMP:-<comp>} DIM=${dtag%d} EBASE=${base} -j8" >&2
    fi
    printf '%s' "$chosen"
}

# require_exe <var> <base> [mpionly] -- pick or exit with a usable message.
require_exe () {
    local var=$1 base=$2 mpionly=${3:-}
    local p
    if ! p=$(pick_exe "$base" "" "$mpionly"); then
        echo "ERROR: no usable ${DTAG:-${DIM:-3}d} build of $base in $SRC" >&2
        exit 1
    fi
    eval "$var=\$p"
}
