#!/usr/bin/env bash
# Count old-path entry attempts across both blend and facade tests.
# Requests include refusals and nested entries, not just successful builds.
# Usage: tools/old-path-all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
if (( $# )); then
    echo "usage: tools/old-path-all.sh (runs both complete suites, including ignored tests)" >&2
    exit 2
fi
logs="$(mktemp -d)"
trap 'rm -f "$logs"/*; rmdir "$logs"' EXIT
for package in ogeom-fillet ogeom; do
    log="$logs/$package.log"
    output="$logs/$package.output"
    : > "$log"
    echo "Counting $package, including ignored tests..." >&2
    if OGEOM_OLD_PATH_LOG="$log" cargo test -p "$package" --profile iterate \
        --features ogeom-fillet/old-path-log --no-fail-fast \
        -- --include-ignored --test-threads=1 > "$output" 2>&1; then
        :
    else
        status=$?
        cat "$output" >&2
        echo "old-path count failed: cargo test -p $package exited with status $status" >&2
        exit "$status"
    fi
done
python3 - "$logs" <<'PY'
import collections
import pathlib
import re
import sys

logs = pathlib.Path(sys.argv[1])
categories = ("round", "chamfer", "vertex", "variable", "marched-backedge")
round_dispatches = {
    "fillet_edge", "fillet_edges", "fillet_edges_with", "fillet_edge_distances", "blend_edges",
}


def classify(entry):
    if ":" in entry:
        category, name = entry.split(":", 1)
        if category not in categories or not name:
            raise ValueError(f"unknown old-path entry: {entry}")
        return category, name
    if entry == "fillet_edges_variable":
        return "variable", entry
    if entry.startswith("chamfer_"):
        return "chamfer", entry
    if entry in round_dispatches:
        return "round", entry
    raise ValueError(f"uncategorized old-path entry: {entry}")


entries = {
    entry
    for source in pathlib.Path("crates/ogeom-fillet/src").rglob("*.rs")
    for entry in re.findall(r'\bto_current_path\("([^"\n]+)"\)', source.read_text())
}
counts = {classify(entry): collections.Counter() for entry in entries}
for package in ("ogeom-fillet", "ogeom"):
    output = (logs / f"{package}.output").read_text()
    results = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;", output)
    if not results or sum(int(passed) for passed, _, _ in results) == 0:
        raise ValueError(f"{package}: no successful test execution summaries")
    if any(int(failed) or int(ignored) for _, failed, ignored in results):
        raise ValueError(f"{package}: failed or ignored tests in an all-tests count")
    print(f"{package}: {sum(int(passed) for passed, _, _ in results)} tests passed; cargo exit 0")
    for line in (logs / f"{package}.log").read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(f"{package}: malformed old-path record: {line!r}")
        entry, test, executable = fields
        if entry not in entries or not test or test == "?" or not executable:
            raise ValueError(f"{package}: unidentified old-path record: {line!r}")
        counts[classify(entry)][(package, executable, test)] += 1

print("\nRequests are entry invocations, including refusals and nested calls.")
print("Distinct tests are keyed by package, test executable and test name; totals are unions.")
print("marched-backedge covers engine-reachable old constructors, including their old-path callers.")
print("\n| Category | Entry point | Requests | Distinct tests |")
print("|---|---|---:|---:|")
all_tests = collections.Counter()
for category in categories:
    category_tests = collections.Counter()
    for (kind, entry), tests in sorted(counts.items()):
        if kind == category:
            category_tests.update(tests)
            print(f"| {category} | `{entry}` | {sum(tests.values())} | {len(tests)} |")
    all_tests.update(category_tests)
    print(f"| **{category}** | **subtotal** | **{sum(category_tests.values())}** | **{len(category_tests)}** |")
print(f"| **all** | **total** | **{sum(all_tests.values())}** | **{len(all_tests)}** |")
PY
