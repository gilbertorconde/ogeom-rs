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
# One full pass over everything, captured so the count at the end comes from
# this run rather than from running the world again.
log="$(mktemp)"
trap 'rm -f "$log"' EXIT
echo "-- full run --"
cargo test --workspace --all-features --no-fail-fast 2>&1 | tee "$log"

# Property tests draw fresh cases each run, so a single green run proves less
# than it looks — but only the *unit* suites hold property tests, and the
# corpus integration suites are deterministic and heavy. Repeat what benefits
# from repetition and leave the rest at one honest pass.
runs="${OGEOM_TEST_RUNS:-2}"
for i in $(seq 1 "$runs"); do
    echo "-- unit re-run $i of $runs --"
    cargo test --workspace --all-features --lib --no-fail-fast
done

echo "== bench (informational) =="
# Watched, not gated: a loaded box would turn a performance gate into a
# coin flip. The ratios are calibrated against a fixed arithmetic spin, so
# they mean the same thing across machines.
cargo run --quiet --release -p ogeom-bench -- --check tools/ogeom-bench/baseline.json || true

echo "== docs =="
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --quiet

echo "== book =="
# The guide's code blocks are included by anchor from crates/ogeom/tests/book.rs,
# which the test pass above just ran — so a book that builds is a book whose
# examples passed. Install the builder with `cargo install mdbook --locked`;
# an absent mdbook fails here on purpose, because a skipped gate is a silent one.
mdbook build docs/book

echo "== parity =="
# The audit gate: every verdict in docs/parity/parity.toml cites evidence that
# still exists — symbols against the rustdoc just built, tests against the
# tree — and docs/PARITY.md matches the committed index. Runs entirely against
# committed files; no reference checkout is consulted here or in CI.
python3 tools/parity.py check

passing=$(grep -E '^test result:' "$log" \
    | awk -F'ok\\. ' '{split($2,a," "); s+=a[1]} END {print s}')
echo
echo "OK — $passing tests passing"
