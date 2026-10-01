---
name: performance-optimizer
description: "Measure a performance problem in the kernel, profile it, and keep only changes that beat a predeclared threshold in paired runs. Use for a slow STEP read, a slow boolean or mesh, a benchmark drift reported by ogeom-bench, or a stress scenario that takes too long. Not for wrong results (use systematic-debugging) and not for a review (use review)."
---

# Performance optimizer

Optimize only a measured problem. Keep a change only when comparable
evidence shows it improves the agreed metric with no broken constraint.

## 1. Define

Before any edit, write down:

- the symptom and the workload that shows it: a corpus file, a stress
  scenario, a benchmark name;
- one primary metric of the same kind as the symptom:

| Symptom | Metric and source |
|---|---|
| A benchmark drifted | the calibrated ratio from `cargo run --release -p ogeom-bench -- --check tools/ogeom-bench/baseline.json` |
| A STEP or IGES file reads slowly | wall time of the read in a release build of a repro program, or the matching benchmark |
| A boolean, blend or mesh is slow | wall time of the operation on the named inputs, median of repeated runs |
| A stress scenario is slow | the per-case seconds from `ogeom-stress --times --scenario <name>` |
| A test job is slow | `cargo test -p <crate> -- -Z unstable-options --report-time` on nightly, or the wall time of the test binary |

- the target, and the minimum gain that keeps an experiment (beyond the
  noise you measured);
- the hard constraints: every result identical to the last ulp, or the
  difference stated and bounded; tolerances unchanged; no `unsafe`; no new
  dependency without the `deny.toml` check; the stress baseline and the
  parity ledger unchanged.

Stop with `BLOCKED` when the problem does not reproduce.

## 2. Baseline

- Run the relevant correctness tests first, so a baseline failure is not
  blamed on an experiment.
- Fix the build mode (`--release`), the inputs, and the warm-up. Cover the
  reported case and the boundary that could reverse the conclusion: the
  small shape as well as the large one, the planar as well as the spline.
- Record raw results, the median, the spread and the machine. Milliseconds
  move with the machine; the bench harness reports every number as a ratio
  to a fixed arithmetic spin so that runs compare across machines. Use the
  ratio when comparing against a recorded baseline.
- Interleave baseline and candidate runs (A, B, A, B) rather than all
  before then all after; a background build or a thermal throttle moves
  everything run in one block.
- Report an inconclusive result as inconclusive. Never rerun until a gain
  appears.

## 3. Profile and hypothesize

Profile the whole path before one function: `cargo flamegraph` on the repro
program, or `perf record` on the test binary, in release with debug
symbols. Separate root cost from first-call cost, debug builds and
instrumentation. Write a short ordered hypothesis list: expected gain,
mechanism, files, risk, how to verify. Drop a hypothesis with no measurable
mechanism or that needs speculative scale.

The usual mechanisms in this kernel: a surface evaluated three times where
one jet would do, a span located by linear search, an allocation per
evaluation, a bounding box rebuilt per query, a section marched at a step
finer than its tolerance needs, work repeated per face that depends only on
the surface.

## 4. One experiment at a time

- Apply the smallest change that tests one mechanism.
- For caching and memoization, protect invalidation: a cached evaluation
  keyed on a parameter must be invalidated when the entity's geometry or
  tolerance changes. For parallelism, keep the result order deterministic
  and the arena append-only.
- Run the focused tests, then repeat the exact baseline measurement.
- `KEEP` only when the gain beats the predeclared threshold and every
  constraint passes. Otherwise `DISCARD` and revert only that experiment.
  Never lower the threshold after seeing the result.
- After a kept change, take a new baseline before the next hypothesis.

## 5. Finish

Run `tools/check.sh`. When a kept change moves a benchmark by design, say
so in the commit message; the bench baseline is informational and is
rewritten only on purpose. A kept change that alters any result is a
correctness change and goes through review.

## 6. Report

Metric, workload and machine; the baseline and final numbers with spread;
the hypothesis ledger with each `KEEP` or `DISCARD` and its evidence; the
checks run; and the verdict: `IMPROVED`, `NO_CHANGE` or `BLOCKED`. A
benchmark gain is not a proven gain on the user's file; say which it is.
