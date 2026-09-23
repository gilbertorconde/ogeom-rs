# Contributing to ogeom

## Independence

ogeom is an independent implementation. It does not depend on, vendor, link
against, bundle or commit any existing CAD kernel, and nothing in this
repository may pull one in. `cargo build` needs a Rust toolchain and nothing
else: no C compiler, no system libraries, no submodules. Keep it that way.

**Never take another kernel's code into this one.** Not as a dependency, a
vendored subtree, a file, or a fragment pasted into a function. Nothing in
`crates/` may be part of someone else's build.

**Do not build a copy of another kernel.** This is the rule the others serve.
A kernel that mirrors another's class hierarchy, decomposition, call graph and
file layout is a translation with different identifiers. That structure is
what we avoid, not the field's vocabulary. What you may take from a reference
implementation is *understanding*: what an algorithm must handle, which cases
exist, what a format's records mean. The design itself is worked out here.

**Use the field's vocabulary.** A boundary representation is a b-rep, a point
in the plane is a `Point2`, a curve in a surface's parameters is a pcurve, a
blend is a fillet. Private jargon would make the kernel harder to read for no
gain. `docs/PLAN.md` and `docs/DATA_MODEL.md` use conventional names on
purpose. Naming a *concept* is not naming a dependency.

**Do not name another kernel in anything committed.** The vocabulary belongs to
the field; product names do not belong to us. Where a format needs its own
magic bytes to be readable, those bytes are data and are committed as data.

### The legal side

- Most existing kernels are copyleft. ogeom is MIT OR Apache-2.0.
- Renaming identifiers does not stop a work being derived from another.
  Copyright follows the expression, not the names. That is why the rule above is
  about *structure*, not words.
- File formats are the safe case. A file format is not copyrightable, and
  implementing one from its published description is interoperation. The STEP
  and `.brep` support were both built this way.

### Sources, in order of preference

1. Published algorithm specifications and papers: Shewchuk on robust
   predicates, Piegl & Tiller on NURBS, the marching-intersection and
   surface-surface literature, the published specifications for boolean
   pipelines.
2. Format standards: ISO 10303 for STEP, and the rest.
3. First principles, and your own tests.

A local reference checkout goes under `vendor/` or `reference/`. Both are
gitignored. A checkout is never a build or test dependency, and it is never
named in anything committed.

If you have contributed to another CAD kernel's source, say so in your pull
request, so we can be careful about which areas you work on.

## Scope

**Parity with the reference kernel's modelling modules (FoundationClasses,
ModelingData, ModelingAlgorithms, DataExchange) and nothing else.**
Visualization, the application framework and the test harness are out.

`docs/SCOPE.md` is normative. It states the rule, how to decide a case
mechanically from the reference tree's own module and toolkit files, and what
the rule does not mean.

Two points are repeated here because both have been misread before:

- **Parity is about capability, not structure.** It does not allow mirroring
  another kernel's class hierarchy or decomposition; *Independence* above still
  holds in full. Where we deliberately do a job differently, the parity record
  says `divergent` and gives the reasoning. That is an answer, not a gap.
- **Usage data sets the order of work, never the scope.**
  `docs/api_surface.json` profiles how one application uses the reference
  kernel, which tells us what to get right first. A capability inside the four
  modules is in scope whether or not that application ever calls it.

## Invariants cannot be changed in a pull request

`docs/DATA_MODEL.md` is normative. A change that breaks one of its invariants
is a design change and must be argued as one, not slipped in. Examples:
flattening the location chain, giving an edge a single curve, adding an
operation that does not emit history, conflating same/equal/partner.

The reason is practical. Each invariant is cheap to keep now and effectively
impossible to add back later across all of a kernel's algorithms. Each one has
its concrete failure mode written next to it.

## Correctness

Geometry code fails quietly. A boolean that returns a plausible but wrong solid
does not throw; it corrupts a document six operations later. So:

- **State the property, then test it.** Round-trips, composition laws, tolerance
  containment, orientation consistency, antisymmetry of predicates. Property
  tests over laws are worth more than many examples. See
  `crates/ogeom-core/tests/properties.rs`.
- **Validate against ground truth you can compute independently.** Analytic
  results for analytic inputs, closed-form volumes and areas, known benchmark
  datasets. "It looks right in the viewer" does not count.
- **There is no external oracle.** Comparing against another kernel is not
  possible, because it would mean vendoring one. Instead, check against closed
  forms, round trips, invariants (volume, area, validity), and properties over
  random inputs. The cheapest useful check is usually: build the same result
  two ways and confirm they agree.
- **Never loosen a tolerance to make a test pass** without explaining why in the
  same commit. Numerical tolerances in tests are explicit and justified in a
  comment.
- **Failures are values.** An algorithm that did not converge returns that fact.
  It does not return an empty shape and set a flag.

## Practical rules

- **Comments describe current behaviour**, never the change that produced it.
  - No "used to", "formerly", "now returns", "since X landed".
  - No issue, PR or commit references as the reason for a behaviour. Those
    belong in the commit message.
  - `node tools/lint-comment-rot.mjs --all` enforces this in `tools/check.sh`.
    `--pedantic` adds an advisory tier. Put `lint-comment-rot: ignore` on a
    line to exempt it where a reference is genuinely needed.
  - The tracked pre-commit hook runs the lint over the lines a commit adds.
    Enable it once per clone with `git config core.hooksPath .githooks`.
- **Run `./tools/check.sh`** before review. It runs formatting, lints, tests and
  docs, and repeats the test suite so a property test that fails only on some
  seeds does not slip through. Do not verify by grepping cargo's output for
  "ok": a run with a failing suite still prints "ok" for every suite that
  passed, so a real failure can hide behind a green-looking summary.
- **Workspace lints** forbid `unsafe`, and warn on `unwrap`/`expect` and lossy
  numeric casts in library code. A kernel is arithmetic from end to end, and
  these are how wrong answers get shipped. A deliberate exception needs an
  `#[allow(..., reason = "...")]` and a documented `# Panics` section.
- **New dependencies** need a permissive license (`deny.toml` enforces this).
  They must not require a C toolchain or break WASM. If you think an exception
  is warranted, raise it explicitly.
  - Record your decision there either way. The record stops the question being
    argued again.
- **Public items are documented.** `missing_docs` is a warning, and CI runs with
  `-D warnings`.

## Releasing

The fifteen library crates under `crates/` are published together, at one
version, from `[workspace.package]`. The three crates under `tools/` have
`publish = false`; they exist only for this repository.

```sh
cargo publish --dry-run --workspace   # packages, verifies and orders, uploads nothing
cargo publish --workspace             # the same, for real
```

Cargo works out the order from the dependency graph, so one command publishes
all fifteen. Know these three points before changing the arrangement:

- **Versions move together.** Every crate inherits `version.workspace`, and the
  internal requirements in `[workspace.dependencies]` are pinned to the same
  number. Bump both, or the workspace stops resolving. Also bump the
  requirements in `outside/Cargo.toml`, because that workspace depends on these
  crates by version.
- **`ogeom-mesh` dev-depends on `ogeom-algo` by path only.** `ogeom-algo`
  depends on `ogeom-mesh`, so a version requirement here would be a cycle no
  registry can resolve. Cargo drops a path-only dev-dependency from the
  published manifest, so the cycle stays inside this repository. Do not "tidy"
  it to `.workspace = true`.
- **The corpus is not published.** `tests/corpus/` sits at the repository root,
  outside every package, and each crate's `exclude` lists the suites that read
  it. A published tarball contains only tests it can run. The repository and
  `tools/check.sh` still run all of them.
