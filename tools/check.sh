#!/usr/bin/env bash
# Verify the workspace: format, lints, tests. Exits non-zero on any failure.
#
# Exists because grepping cargo's output for "ok" is a trap: a suite that fails
# still prints "ok" lines for the suites that passed, so a failing property test
# can hide behind a green-looking summary. This checks exit codes.
set -euo pipefail
cd "$(dirname "$0")/.."

# Note for anyone invoking this: `./tools/check.sh | tail -3` reports tail's exit
# status, not this script's, so a failing run looks like a passing one. Run it
# bare, or keep the pipeline under `set -o pipefail`.

echo "== fmt =="
cargo fmt --all -- --check

echo "== clippy =="
cargo clippy --workspace --all-targets --all-features -- -D warnings

echo "== test =="
# Property tests draw fresh cases each run, so a single green run proves less
# than it looks. Repeat enough to shake out a boundary case that only some
# seeds reach.
runs="${OG_TEST_RUNS:-3}"
for i in $(seq 1 "$runs"); do
    echo "-- run $i of $runs --"
    cargo test --workspace --all-features --no-fail-fast
done

echo "== docs =="
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --quiet

passing=$(cargo test --workspace --no-fail-fast 2>&1 \
    | grep -E '^test result:' \
    | awk -F'ok\\. ' '{split($2,a," "); s+=a[1]} END {print s}')
echo
echo "OK — $passing tests passing"
