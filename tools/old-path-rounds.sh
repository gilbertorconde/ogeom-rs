#!/usr/bin/env bash
# Count the round requests the blend engine still hands to the current
# fillet path (ADR 0001): that path can be deleted only once none do.
#
# Runs every test of ogeom-fillet with the `old-path-log` feature, which
# logs each request leaving the engine with the asking test's name, then
# prints the number of requests, of distinct tests, and the tests by file.
# A test may ask the current path on purpose, as an oracle; such tests are
# ported to the engine before the path is deleted.
#
# Usage: tools/old-path-rounds.sh [extra cargo test arguments]
set -euo pipefail
cd "$(dirname "$0")/.."
log="$(mktemp)"
output="$(mktemp)"
trap 'rm -f "$log" "$output"' EXIT
if OGEOM_OLD_PATH_LOG="$log" cargo test -p ogeom-fillet --features old-path-log \
    --no-fail-fast "$@" > "$output" 2>&1; then
    :
else
    status=$?
    cat "$output" >&2
    echo "round count failed: cargo test exited with status $status" >&2
    exit "$status"
fi
requests=$(wc -l < "$log")
tests=$(awk '{print $2}' "$log" | sort -u)
count=$(printf '%s\n' "$tests" | grep -c . || true)
echo "round requests to the current path: $requests, from $count tests"
printf '%s\n' "$tests" | while read -r t; do
    [ -n "$t" ] || continue
    name="${t##*::}"
    file=$(grep -rl --include='*.rs' "fn ${name}(" crates/ogeom-fillet/tests crates/ogeom-fillet/src \
        | head -1 || true)
    echo "$(basename "${file:-unknown}") $name"
done | sort | awk '{n[$1]++; t[$1]=t[$1] " " $2} END {for (f in n) print n[f], f ":" t[f]}' | sort -rn
