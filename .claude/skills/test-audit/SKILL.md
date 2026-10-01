---
name: test-audit
description: "Invoke whenever writing, changing, reviewing or sweeping tests. Authoring gate for new tests (what it protects, what regression fails it, why existing proof does not, what oracle it measures against) plus an audit workflow for low-value, implementation-coupled or duplicative tests. Not a correctness review of a change (use review) and not a debugging process (use systematic-debugging)."
---

# Test audit

Three modes, one value bar. Authoring mode gates every new or changed test
at write time. Audit mode runs focused sweeps of tests that re-assert the
source, duplicate stronger proof, or couple behaviour to implementation.
Campaign mode prunes one crate's whole test surface; before starting one,
read [CAMPAIGN.md](CAMPAIGN.md).

Two references carry the method.
[references/test-strategy.md](references/test-strategy.md) maps risks to one
portfolio action each and states the proof a deletion needs.
[references/kernel-methodology.md](references/kernel-methodology.md) names
the owning layer and oracle for each kind of contract in this repository.

## Authoring gate

Before adding any test, answer five questions; a missing answer means do not
add it yet:

1. What observable behaviour, invariant or law does it protect?
2. What credible regression makes it fail? A wrong side kept, a tolerance
   violated, a void lost, a refusal that stops being raised.
3. Why does existing coverage not already catch that failure? Each contract
   has one primary owner at the strongest boundary; another layer needs its
   own distinct risk. Prefer adding a case to an existing table or a
   property's strategy over a near-duplicate test.
4. What is the oracle, and is it independent of the code under test? A
   closed form, a published value, a round trip, a second construction of
   the same result, an invariant (volume, area, validity, containment), a
   property over random inputs. Never the implementation's own output
   recorded as the expectation, unless the test says it is a snapshot and
   why that is enough.
5. Is every tolerance in the assertion explicit and justified in a comment
   next to it?

Then check the test against every [junk pattern](#junk-patterns); a match
fails the gate unless the [retention bar](#retention-bar) names the contract
it independently guards. A test that would break under behaviour-preserving
refactoring is asserting implementation, not behaviour; rewrite it at the
owning boundary before landing it.

Bug regression tests must fail on the pre-fix code for the intended reason
and pass after the owner-boundary repair. A regression test that never
demonstrably failed proves nothing. One regression at the owner boundary
covers the bug; do not replay the same scenario at every layer it crosses.
Name the test for the geometric condition, never for the file or the user
that showed it.

## Junk patterns

The shared checklist for both modes: the authoring gate rejects a new test
that matches one, and audits hunt for existing tests that do.

- assertion-free runs: a construction whose only check is that it did not
  panic, where `check` or a measurement was available;
- self-comparisons: a value computed by the kernel compared with the same
  value computed the same way;
- expectations recorded from the implementation's own output without a
  stated reason, so the test pins today's number rather than the truth;
- a tolerance chosen to make the assertion pass, with no reason;
- `is_ok()` or `is_err()` alone where the result's content or the refusal's
  text is the contract;
- a refusal test that accepts any error, so a panic turned into a different
  error still passes;
- private predicate or call-shape tests duplicated at real boundaries;
- duplicate invocations of the same contract across crates: the same box
  volume asserted in four suites;
- copied inventories: a list of symbols, formats or variants re-asserted
  from the source;
- tests whose only purpose is keeping a test-only `pub` item, feature flag
  or helper alive;
- dead production code whose only callers are tests;
- a "looks right" test: a mesh or a drawing asserted only by triangle count
  or by not being empty;
- a test named for a condition stronger than the input exercises, such as
  "a tangent drill" that places the drill a millimetre clear.

## Value bar

Tests justify their maintenance cost by protecting behaviour, a credible
regression, or an independently meaningful contract. In an audit, an
existing test that must change for behaviour-preserving source
reorganization is suspect, not automatically deletable; the authoring gate
still rejects new ones.

Before judging a candidate, read the complete test and production owner, its
entry point, callers, callees, sibling implementations, overlapping tests
and relevant history. Read `CONTRIBUTING.md` and `docs/PLAN.md` first. When
the test claims a published value, find the publication.

## Discovery

Keep discovery read-only and report evidence before editing. For broad
scope, run parallel discovery lanes when available, one per crate group:

- `ogeom-core`, `ogeom-math`, `ogeom-geom`: laws and properties;
- `ogeom-topo`, `ogeom-algo`, `ogeom-mesh`: construction and measurement;
- `ogeom-intersect`, `ogeom-bool`, `ogeom-fillet`, `ogeom-offset`,
  `ogeom-heal`, `ogeom-hlr`: the algorithms;
- `ogeom-io`, `ogeom-doc`, and the corpus suites under `crates/ogeom/tests`;
- a cross-cutting pattern sweep.

Outside campaign mode, prefer a few high-confidence candidates over a large
speculative inventory.

## Retention bar

Keep a test when it independently enforces a public API, a format
(STEP, IGES, the native format, the mesh formats), a law, a tolerance
contract, a refusal by name, a parity verdict's evidence, a book example,
or a documented default. Also keep:

- a regression with a credible failure mode, especially one that reproduced
  a wrong shape passing `check`;
- a test cited as evidence in `docs/parity/parity.toml`; deleting it fails
  the parity gate, and the verdict must be re-evidenced first;
- a test included by anchor into the book (`crates/ogeom/tests/book.rs`);
- a retained test that fails on the baseline: treat it as a possible kernel
  bug, reproduce it, and repair the owner rather than deleting it.

Static or slow is not a deletion reason. A corpus suite that takes seconds
and uniquely proves a reader on outside-authored input is not low value.

## Candidate evidence

Record every field below before editing. A missing field means the
candidate is not ready for deletion:

- exact test name and location;
- what failure it can actually detect;
- non-test callers of the covered production or support item;
- stronger remaining owner-boundary proof, or why no proof is needed;
- the parity verdicts and book anchors that cite it;
- relevant history and the reason the test exists;
- production or test-support deletion unlocked;
- risk and the focused validation command.

## Edit shape

Choose one coherent owner-boundary batch. Delete obsolete test-only items
and dead production paths instead of preserving them. Move retained
regressions to their canonical owners. Consolidate repeated assertions into
one property or one table.

Prefer net-negative production lines. Do not add replacement tests that
restate the same implementation, and do not convert uncertain candidates
into cleanup to increase deletion counts.

## Validation

Never edit source or tests while a test runner is running in the checkout.

1. Run the smallest owner and sibling tests: `cargo test -p <crate> --test
   <file>` or `cargo test -p <crate> <filter>`.
2. Run `cargo fmt --all -- --check` and `git diff --check`.
3. Run `python3 tools/parity.py check` when a cited test moved or went.
4. Run `tools/check.sh` before handing off: it repeats the unit suites so
   a property that fails only on some seeds is seen.
5. Inspect `git diff --numstat`; report production and tooling separately
   from tests and test support.
6. After final audit edits, run the review skill over the diff.

## Landing

Commit only when authorized, one scope per commit in the repo convention.
Never push; hand off. Land one coherent batch at a time; after landing,
refresh from `main` and rerun read-only discovery for the next batch.

## Handoff

Report:

- root cause and removed low-value categories;
- production owner simplifications;
- retained false positives and why they remain valuable;
- focused and full proof actually run;
- production versus test lines;
- named follow-ups.
