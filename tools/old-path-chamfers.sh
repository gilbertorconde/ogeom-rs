#!/usr/bin/env bash
# Count the chamfer requests the blend engine still hands to the current
# chamfer path, the wedges (ADR 0001): that path can be deleted only once
# none do. The chamfer side's twin of tools/old-path-rounds.sh.
#
# Runs every test of ogeom-fillet and of ogeom with the `old-path-log`
# feature, which logs each request leaving the engine with the asking
# test's name, then prints the number of chamfer requests, of distinct
# tests, and the tests by file.
#
# Usage: tools/old-path-chamfers.sh [extra cargo test arguments]
set -euo pipefail
cd "$(dirname "$0")/.."
log="$(mktemp)"
trap 'rm -f "$log"' EXIT
OGEOM_OLD_PATH_LOG="$log" cargo test -p ogeom-fillet --features old-path-log \
    --no-fail-fast "$@" > /dev/null 2>&1 || true
OGEOM_OLD_PATH_LOG="$log" cargo test -p ogeom --features ogeom-fillet/old-path-log \
    --no-fail-fast "$@" > /dev/null 2>&1 || true
chamfers="$(mktemp)"
trap 'rm -f "$log" "$chamfers"' EXIT
grep '^chamfer' "$log" > "$chamfers" || true
requests=$(wc -l < "$chamfers")
tests=$(awk '{print $2}' "$chamfers" | sort -u)
count=$(printf '%s\n' "$tests" | grep -c . || true)
echo "chamfer requests to the current path: $requests, from $count tests"
printf '%s\n' "$tests" | while read -r t; do
    [ -n "$t" ] || continue
    name="${t##*::}"
    file=$(grep -rl --include='*.rs' "fn ${name}(" crates/ogeom-fillet/tests crates/ogeom/tests \
        crates/ogeom-fillet/src | head -1 || true)
    echo "$(basename "${file:-unknown}") $name"
done | sort | awk '{n[$1]++; t[$1]=t[$1] " " $2} END {for (f in n) print n[f], f ":" t[f]}' | sort -rn
