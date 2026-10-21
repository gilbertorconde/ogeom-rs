# Contributing to ogeom

## Independence

ogeom is an independent implementation. It does not depend on, vendor, link
against, bundle or commit any existing CAD kernel, and nothing in this
repository may pull one in. `cargo build` needs a Rust toolchain and nothing
else: no C compiler, no system libraries, no submodules. Keep it that way.

**Never take another kernel's code into this one.** Not as a dependency, a
vendored subtree, a file, or a fragment pasted into a function. Nothing in
`crates/` may be part of someone else's build.

**Do not read another kernel's implementation source to work on this one.**
Not its `.cxx`, `.cpp`, `.c`, `.lxx` or `.gxx` files, and not the bodies of
its headers. Do not keep a checkout of one next to this repository. The same
rule binds tools: never give an AI assistant, a script or a search index
another kernel's implementation source as context for work here, and never
ask one to explain or summarise how another kernel's code does something.

**Do not build a copy of another kernel.** A kernel that mirrors another's
class hierarchy, decomposition, call graph and file layout is a translation
with different identifiers. The design is worked out here, from the sources
below.

**Use the field's vocabulary.** A boundary representation is a b-rep, a point
in the plane is a `Point2`, a curve in a surface's parameters is a pcurve, a
blend is a fillet. Private jargon would make the kernel harder to read for no
gain. Naming a *concept* is not naming a dependency.

**Do not name another kernel in anything committed.** The vocabulary belongs to
the field; product names do not belong to us. Two exceptions: `PROVENANCE.md`,
which records where the kernel's design came from, and data a format needs
to be readable (magic bytes, a writer's name in a file header), committed as
data.

### The legal side

- Most existing kernels are copyleft. ogeom is MIT OR Apache-2.0.
- Renaming identifiers does not stop a work being derived from another.
  Copyright follows the expression, not the names. That is why the rules above
  cover reading, not only copying.
- File formats are the safe case. A file format is not copyrightable, and
  implementing one from its published description is interoperation. The STEP
  and `.brep` support were both built this way.

### Sources

1. Published algorithm specifications and papers: Shewchuk on robust
   predicates, Piegl & Tiller on NURBS, the marching-intersection and
   surface-surface literature, the published descriptions of boolean
   pipelines.
2. Format standards: ISO 10303 for STEP, and the rest.
3. A kernel's public user documentation, for what an operation is expected to
   do. Not its source.
4. First principles, and your own tests.

If you have contributed to another CAD kernel's source, say so in your pull
request, so we can be careful about which areas you work on.

`PROVENANCE.md` states how these rules have been applied over the project's
history.

## Scope

**A complete solid-modelling kernel and its data exchange, and nothing else.**
Rendering, application document frameworks and scripting harnesses are out.

`docs/SCOPE.md` is normative. It states the rule, how to decide a case, and
what the rule does not mean. `docs/CAPABILITIES.md` records, capability by
capability, where the kernel stands.

Where we deliberately do a job differently from the usual way, the capability
record says `divergent` and gives the reasoning. That is an answer, not a gap.

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
- **Tests have a budget.** Every test in the default tier finishes within
  5 s in the `iterate` profile. A test that needs longer (a sweep, a corpus
  part, a fine mesh) is marked `#[ignore = "heavy"]`; prefer keeping a
  small representative case in the default tier beside it.
  - `cargo t` (an alias for `cargo test --profile iterate`) is the working
    loop: the crates you touched, plus the heavy tests of that area by name
    with `-- --include-ignored`.
  - Do not test with `--release`: it links every test binary with link-time
    optimisation and one codegen unit, ten seconds or more each, for the same
    answers.
  - `tools/check.sh` and CI run both tiers. `tools/test-times.sh` prints the
    slowest tests of the default tier and fails past 15 s.
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

## Skills for coding agents

`.claude/skills/<name>/SKILL.md` holds the working methods an agent follows
for a kind of task. They cite the rules above; they never weaken one. Read
the skill before the task it covers.

| Skill | Use when |
|---|---|
| `review` | reviewing a branch or PR for wrong geometry, broken invariants and release risk |
| `systematic-debugging` | any bug, failure, wrong shape or refusal, before proposing a fix |
| `test-audit` | writing, changing, reviewing or pruning tests |
| `performance-optimizer` | a measured slowness or a benchmark drift |
| `ci-iterate` | a CI run is red, or a PR has review feedback to address |
| `comment-cleanup` | before committing, when the comment-rot lint fails, or to trim or unify comments |

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
