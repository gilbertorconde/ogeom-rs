#!/usr/bin/env bash
# The quick check for work on the mesh converter: formatting, lints and the
# tests of the crates it lives in and reads through, built with the iterate
# profile (release speed, incremental rebuilds). It is not the gate:
# ./tools/check.sh still runs before a change leaves this area.
#
# Set OGEOM_CONVERTER_BENCH to a command to run it last, for a benchmark kept
# outside the repository (meshes that may not be committed).
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== fmt =="
cargo fmt --all -- --check
echo "== comment rot =="
node tools/lint-comment-rot.mjs --all
echo "== clippy =="
cargo clippy -p ogeom-algo -p ogeom-io --all-targets -- -D warnings
cargo clippy -p ogeom --test mesh_solid --test mesh_steps --test mesh_corpus --test boolean_scale -- -D warnings
echo "== test =="
cargo test --profile iterate -q -p ogeom-algo --lib
cargo test --profile iterate -q -p ogeom-io --lib
# The converter's heavy tests too: this is the area's own check.
cargo test --profile iterate -q -p ogeom --test mesh_solid --test mesh_steps --test mesh_corpus --test boolean_scale -- --include-ignored
if [ -n "${OGEOM_CONVERTER_BENCH:-}" ]; then
    echo "== bench =="
    eval "$OGEOM_CONVERTER_BENCH"
fi
echo "OK"
