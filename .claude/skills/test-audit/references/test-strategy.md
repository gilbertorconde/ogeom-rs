# Test strategy

Use this before adding tests for a feature, and during an audit when
deciding what a test is worth. It turns risks into one decision per risk and
per affected test.

## 1. Risk map

1. State the protected outcome and the acceptance criteria. Stop with
   `BLOCKED` when there is no concrete behaviour to protect.
2. Trace each critical path from input through the algorithm to a
   measurable result: a volume, a face count, a distance, a round trip, a
   refusal.
3. Mark the behaviour that matters most here: a wrong shape that passes
   `check`, orientation and handedness, tolerance containment, location
   chains, history, the format read back, the refusal by name.
4. List plausible defect classes: wrong success, rejected valid input,
   accepted invalid input, boundary error (seam, pole, apex, domain end),
   degenerate placement (tangent, coincident, through a vertex), scale
   (unit, tiny, huge), mirror, partial result, non-convergence reported as
   success.
5. Drop behaviour a dependency already guarantees and no result depends on.
   Rank the rest qualitatively; do not invent numbers.

## 2. One action per risk and per affected test

| Action | Use when |
|---|---|
| `KEEP` | Trusted, unique proof is still valid |
| `ADD` | A material risk has no proof |
| `UPDATE` | The intent is valuable but the boundary, setup or oracle changed |
| `MERGE` | Proof can be consolidated without losing a scenario or failure localization |
| `DELETE` | The proof is obsolete, duplicate, trivial or untrustworthy |
| `NO_TEST` | Another control or an accepted residual risk covers it; name which |

The action is separate from the run state (`PASS`, `FAIL`, `BLOCKED`,
`UNPROVEN`). Zero selected tests prove nothing.

## 3. Deletion and merge need proof

- Delete or merge only when the basis is obsolete, or other evidence covers
  every still-required behaviour and failure mode with equal or better
  trust.
- A regression guard for a fixed wrong shape stays. So does the only proof
  of a rare critical edge.
- Coverage shows execution, not proof. A slow corpus suite that uniquely
  proves a reader is not low value; a fast test that proves the standard
  library is.
- Test count, pass rate and coverage are not quality targets.

## 4. Level, oracle and gate

- One owner per contract, at the strongest boundary that is cheap and
  deterministic enough. Another level needs its own distinct risk.
- The oracle is independent of the code under test: a closed form, a
  published value, a round trip, a second construction, an invariant, a law
  over random inputs. Never recompute the expected value with the
  implementation's own logic.
- Control randomness: property tests use a seeded strategy and record their
  failing cases; the stress harness is seeded.
- Say which gate owns the test (`tools/check.sh`, the nightly stress run,
  or local only) and its measured run cost when it adds more than a few
  seconds.

## 5. Output

For each decision: basis, protected outcome, risk, existing test or gap,
action, level, oracle, gate, and run state. Verdict: `READY`,
`INCONCLUSIVE` (name the next evidence step) or `BLOCKED`.
