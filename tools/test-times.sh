#!/usr/bin/env bash
# Time every test of the default tier and hold it to its budget.
#
# The default tier is what `cargo t` runs: every test not marked
# `#[ignore = "heavy"]`. Each of its tests is meant to finish within
# BUDGET seconds in the `iterate` profile; a test over it belongs in the
# heavy tier, which `tools/check.sh` and CI run with `--include-ignored`.
#
# Prints the slowest tests and exits non-zero when one runs past LIMIT
# seconds. LIMIT sits well above BUDGET because the machine may be shared:
# a test just over the budget is reported, a test far past it fails.
#
# Per-test times come from libtest's `--report-time`, an unstable flag;
# RUSTC_BOOTSTRAP lets the stable toolchain's test harness accept it. Only
# the harness reads it here: the code under test is built as usual.
#
# Usage: tools/test-times.sh [cargo test arguments, e.g. -p ogeom]
set -euo pipefail
cd "$(dirname "$0")/.."

budget="${BUDGET:-5}"
limit="${LIMIT:-15}"
shown="${SHOWN:-20}"
log="$(mktemp)"
trap 'rm -f "$log"' EXIT

if [ "$#" -eq 0 ]; then
    set -- --workspace
fi

RUSTC_BOOTSTRAP=1 cargo test --profile iterate --no-fail-fast "$@" -- \
    -Z unstable-options --report-time --format pretty >"$log" 2>&1 || {
    grep -E "panicked|FAILED|^error" "$log" | head -40
    echo "test-times: the tests did not pass" >&2
    exit 1
}

python3 - "$log" "$budget" "$limit" "$shown" <<'EOF'
import re
import sys

log, budget, limit, shown = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
suite = None
rows = []
for line in open(log):
    m = re.search(r"Running (\S+)", line)
    if m:
        suite = m.group(1)
    m = re.match(r"test (\S+) \.\.\. ok <([\d.]+)s>", line)
    if m:
        rows.append((float(m.group(2)), suite, m.group(1)))
rows.sort(reverse=True)
total = sum(r[0] for r in rows)
over = [r for r in rows if r[0] > budget]
past = [r for r in rows if r[0] > limit]
print(f"{len(rows)} tests, {total:.0f} s in all; {len(over)} over the {budget:g} s budget")
for t, s, n in rows[:shown]:
    mark = "  past the limit" if t > limit else ("  over budget" if t > budget else "")
    print(f"{t:7.1f} s  {s}  {n}{mark}")
if past:
    print(f"{len(past)} tests past {limit:g} s: mark them heavy or make them cheaper")
    sys.exit(1)
EOF
