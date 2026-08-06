# Contributing to ogeom

## Independence

ogeom is an independent implementation. It does not depend on, vendor,
link against, bundle, or commit any existing CAD kernel, and nothing in this
repository will pull one in. `cargo build` needs a Rust toolchain and nothing
else — no C compiler, no system libraries, no submodules. Keep it that way.

**Never take another kernel's code into this one.** Not as a dependency, not
as a vendored subtree, not as a file, not as a fragment pasted into a
function. Nothing in `crates/` may be a piece of somebody else's build.

**Never name another kernel in this repository's source or documentation.**
Not in a comment, not in an error message, not in a commit message, not in a
design note. Concepts get our own names throughout; where a format or an
algorithm has a conventional name, use the neutral one — a file format is
named by its own extension, an algorithm by what it does.

Where a reference implementation informs the work, what crosses over is
*understanding*: what an algorithm must handle, what cases exist, what a
format's records mean. What is written here is written here, in our own
structure and our own vocabulary.

Be aware of what that means legally. Most existing kernels are copyleft;
ogeom is MIT OR Apache-2.0. Renaming identifiers does not change whether a
work is derived from another — copyright follows the expression, not the
names — so the further an implementation stays from somebody else's
structure, the sounder its footing. Formats are the safe case: a file format
is not copyrightable, and implementing one from its published description is
interoperation, which is exactly how the STEP support was built.

Prefer, in order:

- published algorithm specifications and papers — Shewchuk on robust predicates,
  Piegl & Tiller on NURBS, the marching-intersection and surface-surface
  literature, the published specifications for boolean pipelines;
- format standards — ISO 10303 for STEP, and the rest;
- first principles, and your own tests.

A reference checkout kept locally goes under `vendor/` or `reference/`, both
of which are gitignored. It never becomes a build or test dependency, and it
never appears by name in anything committed.

Naming and vocabulary are a different matter. `docs/PLAN.md` and
`docs/DATA_MODEL.md` cite conventional names for concepts, because that is how
the field talks about itself and inventing private jargon would help nobody.
That is a glossary, not a dependency.

If you have contributed to another CAD kernel's source, say so in your pull
request so we can be careful about which areas you touch.

## Scope

Scope is set by what a CAD kernel has to do, not by what any particular
application asks for. "Application X implements this itself" and "application X
rarely calls this" are not arguments for leaving something out — see
`docs/PLAN.md`, which explains why an earlier draft of that document got this
wrong and what it cost.

Usage data is a **sequencing** input: it tells us what to make fast and get right
first. It is not a scope input.

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
