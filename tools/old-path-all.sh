#!/usr/bin/env bash
# Count old-path entry attempts across both blend and facade tests.
# Requests include refusals and nested entries, not just successful builds.
# Usage: tools/old-path-all.sh [--evidence-dir NEW_DIRECTORY]
set -euo pipefail
cd "$(dirname "$0")/.."
if (( $# == 0 )); then
    logs="$(mktemp -d)"
    trap 'rm -f "$logs"/*; rmdir "$logs"' EXIT
elif (( $# == 2 )) && [[ "$1" == --evidence-dir ]]; then
    mkdir -- "$2"
    logs="$(cd "$2" && pwd)"
    echo "Retaining count evidence in $logs" >&2
else
    echo "usage: tools/old-path-all.sh [--evidence-dir NEW_DIRECTORY]" >&2
    exit 2
fi
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
import json
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
public_requests = {key: set() for key in counts}
nested = collections.Counter()
request_owners = {}
records = []
unscoped = 0
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
        if len(fields) != 6:
            raise ValueError(f"{package}: malformed old-path record: {line!r}")
        entry, test, executable, request, api, depth = fields
        if entry not in entries or not test or test == "?" or not executable:
            raise ValueError(f"{package}: unidentified old-path record: {line!r}")
        if not request.isdecimal() or not depth.isdecimal() or (int(request) == 0) != (api == "-"):
            raise ValueError(f"{package}: invalid request identity: {line!r}")
        key = classify(entry)
        target = re.sub(r"-[0-9a-f]{16}$", "", pathlib.Path(executable).stem)
        family = f"{package}/{target}"
        if target == package.replace("-", "_") and "::" in test:
            family += "/" + test.rsplit("::", 1)[0]
        records.append({
            "category": key[0], "entry": key[1], "family": family, "test": test,
            "package": package, "executable": executable,
            "public_request": int(request) or None, "public_api": api,
            "nested": int(depth) > 0,
        })
        counts[key][(package, executable, test)] += 1
        nested[key] += int(depth) > 0
        if int(request):
            identity = (package, executable, int(request))
            owner = (test, api)
            if request_owners.setdefault(identity, owner) != owner:
                raise ValueError(f"{package}: reused request identity: {line!r}")
            public_requests[key].add(identity)
        else:
            unscoped += 1

(logs / "records.json").write_text(json.dumps(records, indent=2) + "\n")
print("\nFallback records and constructor invocations include refusals and declined attempts.")
print("Nested invocations are the subset inside another logged entry's call stack.")
print("Distinct public requests deduplicate the outermost public API call across all its entries.")
print(f"Private/direct invocations without a public API scope: {unscoped}; these add no public request.")
print("Distinct tests use package, executable and test name; test and request totals are unions.")
print("marched-backedge covers engine-reachable old constructors, including their old-path callers.")
print("\n| Category | Entry point | Public fallback records | Constructor invocations | Nested invocations | Distinct tests | Distinct public requests |")
print("|---|---|---:|---:|---:|---:|---:|")
all_tests = collections.Counter()
all_requests = set()
all_dispatches = all_constructors = all_nested = 0
for category in categories:
    category_tests = collections.Counter()
    category_requests = set()
    dispatches = constructors = category_nested = 0
    for key, tests in sorted(counts.items()):
        kind, entry = key
        if kind == category:
            category_tests.update(tests)
            category_requests.update(public_requests[key])
            is_dispatch = entry in entries
            n = sum(tests.values())
            d, c = (n, 0) if is_dispatch else (0, n)
            dispatches += d
            constructors += c
            category_nested += nested[key]
            print(f"| {category} | `{entry}` | {d} | {c} | {nested[key]} | {len(tests)} | {len(public_requests[key])} |")
    all_tests.update(category_tests)
    all_requests.update(category_requests)
    all_dispatches += dispatches
    all_constructors += constructors
    all_nested += category_nested
    print(f"| **{category}** | **subtotal** | **{dispatches}** | **{constructors}** | **{category_nested}** | **{len(category_tests)}** | **{len(category_requests)}** |")
print(f"| **all** | **total** | **{all_dispatches}** | **{all_constructors}** | **{all_nested}** | **{len(all_tests)}** | **{len(all_requests)}** |")
PY
