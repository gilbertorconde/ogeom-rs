---
name: review
description: "Read-only review of a branch diff or PR for wrong geometry, broken invariants and release risk, with an optional independent panel (at most three lenses, at most two rounds). Use before opening a PR, before merge, or when asked for an adversarial or multi-agent review of a change. Not for a whole-kernel audit (that is docs/REVIEW.md), a test-portfolio audit (use test-audit) or a red CI run (use ci-iterate)."
---

# Review

One lead owns scope, evidence and verdict. The review is read-only: no
edits, commits, pushes, PR comments or issue writes unless the user asks.

## 1. Scope and tier

1. Resolve base, head and the full diff: `git diff origin/main...HEAD`. If
   the output is truncated, read each changed file until every changed line
   is seen. Require a clean tree, or name the uncommitted changes you
   reviewed.
2. State the protected outcome, the acceptance criteria and the non-goals
   from the PR body, the linked issue, the parity row or the request. Mark a
   guess `UNKNOWN`.
3. Drop diff noise (`Cargo.lock`, a regenerated `docs/PARITY.md`, corpus
   files) but keep `docs/parity/parity.toml`, `deny.toml`, the baselines
   under `tools/`, and anything under `.github/`.
4. Pick the tier:

| Tier | When | Depth |
|---|---|---|
| trivial | 10 changed lines or fewer, no sensitive path | direct read, no panel |
| lite | 100 lines or fewer, no sensitive path | direct read, one lens at most |
| full | larger, or any sensitive path | direct read plus the panel as needed |

Sensitive paths always get `full`: `ogeom-core` (tolerances, predicates,
arenas), `ogeom-topo`, `ogeom-bool`, `ogeom-intersect`, the STEP and IGES
readers and writers in `ogeom-io`, `docs/DATA_MODEL.md`, `docs/SCOPE.md`,
`tools/check.sh` and the workflows.

## 2. Rules

Each finding cites the rule it breaks. The slugs come from
`CONTRIBUTING.md`, `docs/PLAN.md` and `docs/DATA_MODEL.md`:

| Slug | Rule |
|---|---|
| `measured-not-asserted` | A claim about geometry is backed by a test against a closed form, a published value or an independent computation |
| `approximation-says-so` | Sampling, fitting and tolerances are stated; no number hides its own error |
| `refused-by-name` | What the kernel cannot do returns an error that says which thing and why, and a test pins it; never a plausible wrong shape |
| `no-stubs` | Implemented, or refused by name and tracked in `docs/PLAN.md` |
| `tolerance-justified` | No tolerance loosened to pass a test without the reason in the same commit; tolerances in tests are explicit and justified |
| `invariants-not-in-a-pr` | Nothing in `docs/DATA_MODEL.md` changes without being argued as a design change |
| `history-emitted` | Every operation reports what became of its inputs |
| `tolerances-only-grow` | Per-entity tolerances widen, never shrink, and containment holds across boundary relationships |
| `present-tense-comments` | Comments describe present behaviour: no change narration, no issue or commit as the reason, no mention of a file the author tested with |
| `independence` | Nothing names, links, vendors or mirrors another kernel; the field's vocabulary only |
| `scope` | `docs/SCOPE.md` decides; usage data orders the work and never widens the scope |
| `no-unwrap-in-library` | `unwrap`, `expect`, lossy casts and `unsafe` in library code need an `allow` with a reason and a documented `# Panics` |

## 3. Direct review

For each changed file, map what can go quietly wrong: inputs and their
domains, tolerance use, orientation and handedness, parameter domains
(periodic, seams, poles, degenerate edges), location chains, history, and
the errors returned. Then check every item for every file:

- a plausible wrong result that passes `check`: the wrong side kept, an
  orientation flipped, a void lost, a volume off by one face;
- handedness: mirrors, negative scales, reversed frames and axes;
- tolerance: a comparison at the wrong scale, `confusion` where
  `intersection` or `approximation` is meant, a tolerance read from one
  entity and applied to another, a test that compares at a tolerance the
  code does not state;
- parameterization: periodic domains, a seam walked twice, poles and cone
  apexes, unclamped knots, a surface window mistaken for a trim;
- locations: a result built on a node without the input's location, a
  placement baked where it must not be or left unbaked where it must;
- convergence: a Newton step or a march that returns a value without having
  converged, a sample count that bounds correctness;
- degenerate placements: tangent, coincident, through a vertex, zero
  length, near-coplanar inside the weld distance;
- units: a scale assumed to be millimetres;
- errors: a failure turned into an empty shape or a flag; a refusal without
  a name, or without a test that pins its text;
- superseded code, aliases and dual paths removed; the parity row and
  `CHANGELOG.md` updated when behaviour changes.

Trace each critical scenario from input through the algorithm to the
measured result. Read unchanged code only to prove an affected path. Build a
counterexample where you can: a repro program in a scratch crate outside the
repository that depends on the workspace by path, a replay with
`ogeom-stress --case <scenario>/<part>/<n>`, or `ogeom-cli`.

## 4. Independent panel

Use no panel when direct evidence settles the verdict. Otherwise read
[references/independent-review.md](references/independent-review.md) before
you launch anyone. Limits: at most three lenses, at most two rounds. Every
reviewer gets the same frozen packet and never sees a sibling's findings.
After round two the lead checks the remaining evidence directly.

## 5. Materiality gate

Accept a finding only when all hold:

- the diff introduced, exposed or worsened it;
- it has `file:line`, a defect class, a concrete failure scenario and the
  smallest fix;
- it is not already handled elsewhere (search the code and the tests).

Reject taste, theory, generic practice, hypothetical scale, and a design
that is only different. Deduplicate by root cause. Every correctness claim
needs a counterexample: a shape, a placement, a parameter, a number.

## 6. Verdict

Severity: `P0` a wrong answer that passes `check`, `P1` blocks release,
`P2` important, `P3` minor.

- `FAIL`: an open P0 or P1, unmet acceptance, a required gate the change
  broke, or a rule slug broken on a sensitive path.
- `CONCERNS`: only non-blocking risk remains.
- `PASS`: every acceptance criterion has evidence and no finding survives.
- `BLOCKED`: required evidence (a gate, a corpus file, a closed form) is
  missing and has no credible substitute.

Report in this order: verdict, scope (base, head, tier, exclusions),
findings (severity, `file:line`, slug, failure scenario, fix), checks run
with results (`tools/check.sh` is the full gate; name which parts ran), and
what stays unproven. "No findings" is a valid result.
