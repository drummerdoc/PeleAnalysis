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
        return 1
    fi

    # The staleness comparison. Only sources that are linked into these
    # drivers live in Src, so the newest of them is the right reference.
    local newest
    newest=$(ls -t "$SRC"/*.cpp "$SRC"/*.H 2>/dev/null | head -1)
    if [ -n "$newest" ] && [ "$newest" -nt "$chosen" ]; then
        echo "ERROR: $(basename "$chosen") is older than $(basename "$newest")." >&2
        echo "  Rebuild before testing. A stale binary answers confidently and" >&2
        echo "  the failure then looks like a defect in the mesh or the data." >&2
        echo "  make -f GNUmakefile.AMRComplex USE_MPI=TRUE COMP=${COMP:-<comp>} DIM=${dtag%d} EBASE=${base} -j8" >&2
        return 1
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
