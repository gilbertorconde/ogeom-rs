# Contributing to ogeom

## Independence

ogeom is an independent implementation. It does not depend on, vendor,
link against, bundle, or commit any existing CAD kernel, and nothing in this
repository will pull one in. `cargo build` needs a Rust toolchain and nothing
else — no C compiler, no system libraries, no submodules. Keep it that way.

**Never take another kernel's code into this one.** Not as a dependency, not
as a vendored subtree, not as a file, not as a fragment pasted into a
function. Nothing in `crates/` may be a piece of somebody else's build.

**Do not build a carbon copy of one.** This is the rule the others serve. A
kernel that mirrors another's class hierarchy, its decomposition, its call
graph and its file layout is a translation wearing different identifiers, and
that is what is being avoided — not the words the field uses. What crosses
over from a reference implementation is *understanding*: what an algorithm
must handle, which cases exist, what a format's records mean. The design here
is arrived at here.

**Use the field's vocabulary.** A boundary representation is a b-rep, a point
in the plane is a `Point2`, a curve in a surface's parameters is a pcurve, a
blend is a fillet. That is what the field calls these things, and inventing
private jargon would make this kernel harder to read for no gain at all.
`docs/PLAN.md` and `docs/DATA_MODEL.md` use conventional names deliberately.
Naming a *concept* is not naming a dependency.

**Do not name another kernel in anything committed.** The vocabulary is the
field's; the product names are not ours to carry. Where a format needs its own
magic bytes to be readable at all, those bytes are data and go in as data.

Be aware of what all this means legally. Most existing kernels are copyleft;
ogeom is MIT OR Apache-2.0. Renaming identifiers does not change whether a
work is derived from another — copyright follows the expression, not the names
— which is exactly why the rule above is about *structure* and not about
words. Formats are the safe case: a file format is not copyrightable, and
implementing one from its published description is interoperation, which is
how both the STEP and the `.brep` support were built.

Prefer, in order:

- published algorithm specifications and papers — Shewchuk on robust predicates,
  Piegl & Tiller on NURBS, the marching-intersection and surface-surface
  literature, the published specifications for boolean pipelines;
- format standards — ISO 10303 for STEP, and the rest;
- first principles, and your own tests.

A reference checkout kept locally goes under `vendor/` or `reference/`, both
of which are gitignored. It never becomes a build or test dependency, and it
never appears by name in anything committed.

If you have contributed to another CAD kernel's source, say so in your pull
request so we can be careful about which areas you touch.

## Scope

**Parity with the reference kernel's modelling modules — FoundationClasses,
ModelingData, ModelingAlgorithms, DataExchange — and nothing else.**
Visualization, the application framework and the test harness are out.
`docs/SCOPE.md` is normative: it states the rule, how a case is decided
mechanically from the reference tree's own module and toolkit files, and what
the rule is not.

Two things follow that are worth stating here, because both have been got wrong:

Parity is a claim about **capability, not structure**. It is not a licence to
mirror another kernel's class hierarchy or decomposition — see *Independence*
above, which still holds in full. Where we do a job differently on purpose, the
parity record says `divergent` and gives the reasoning; that is an answer, not a
gap.

Usage data is a **sequencing** input, never a scope input. `docs/api_surface.json`
profiles how one application exercises the reference; it says what to get right
first. A capability inside the four modules is in scope whether or not that
application ever calls it.

## Invariants are not negotiable in a pull request

`docs/DATA_MODEL.md` is normative. A change that breaks one of its invariants —
flattening the location chain, giving an edge a single curve, adding an operation
that does not emit history, conflating same/equal/partner — is a design change and
has to be argued as one, not slipped in.

The reason is not purity. Each is cheap now and effectively impossible to retrofit
across a kernel's worth of algorithms later, and each has a concrete failure mode
written next to it.

## Correctness

Geometry code fails quietly. A boolean that produces a plausible-looking but
wrong solid does not throw; it corrupts a document six operations later. So:

- **State the property, then test it.** Round-trips, composition laws, tolerance
  containment, orientation consistency, antisymmetry of predicates. Property
  tests over laws are worth more than a pile of examples — see
  `crates/ogeom-core/tests/properties.rs`.
- **Validate against ground truth you can compute independently.** Analytic
  results for analytic inputs; closed-form volumes and areas; known benchmark
  datasets. Not "it looks right in the viewer".
- **There is no external oracle**, and comparing against another kernel is not
  available to us — it would mean vendoring one. `docs/PLAN.md`'s
  [Verification](docs/PLAN.md#verification) section sets out the five
  directions that stand in for one, and which kinds of defect each is there to
  catch. Reach for the one that fits what you are adding; the cheapest useful
  one is usually asking whether the same result built two ways agrees.
- **Never loosen a tolerance to make a test pass** without saying why in the same
  commit. Numerical tolerances in tests are explicit and justified in a comment.
- **Failures are values.** An algorithm that did not converge returns that fact.
  It does not return an empty shape and set a flag.

## Practical

- **Run `./tools/check.sh`** before review — format, lints, tests and docs, with
  the test suite repeated so a property test that only fails on some seeds does
  not slip through. Do not verify by grepping cargo's output for "ok": a run
  with a failing suite still prints "ok" for every suite that passed, and a real
  failure hides behind a green-looking summary.
- Workspace lints forbid `unsafe` and warn on `unwrap`/`expect` and lossy numeric
  casts in library code — a kernel is arithmetic end to end, and those are how
  wrong answers get shipped. A deliberate exception carries an
  `#[allow(..., reason = "...")]` and a documented `# Panics` section.
- New dependencies need a permissive license (`deny.toml` enforces this), and must
  not introduce a C toolchain requirement or break WASM. Raise it explicitly if
  you think an exception is warranted. Check
  [Dependencies considered](docs/PLAN.md#appendix-dependencies-considered)
  first — several plausible-looking crates have already been turned down for
  reasons that have not changed, and one is scheduled to be re-examined at a
  specific point rather than whenever it next comes up. Add your decision there
  either way: the record is what stops the question being re-litigated.
- Public items are documented. `missing_docs` is a warning and CI runs with
  `-D warnings`.
