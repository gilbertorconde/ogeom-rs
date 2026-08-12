# Introduction

**ogeom** is a boundary-representation CAD kernel, written from scratch in
Rust. It is not a binding, a wrapper, or a port: it depends on no existing
CAD kernel, vendors none, and links against none.

What it covers: parametric curves and surfaces, shared B-rep topology with
per-entity tolerances and operation history, booleans, blends and chamfers,
offsets and shells, sweeps and lofts, healing, tessellation, product
structure, 2D drawings with hidden line removal, and data exchange — STEP in
both directions with assemblies, colours, semantic PMI and saved views, IGES
in both directions, the native document format, and the mesh formats.

Two things distinguish how this kernel is built, and both show up in how
this book is written:

**Claims are measured, not asserted.** Every numeric statement in the test
suite is held to a closed form or a round trip. The same rule applies here:
every code block in this guide is included from
[a real test file](https://github.com/gilbertorconde/ogeom-rs/blob/main/crates/ogeom/tests/book.rs)
that runs in CI, so an example that stops compiling — or stops being true —
fails the build instead of lying in the documentation.

**Completeness is audited, not felt.** The kernel's target is written down
in [Scope](kernel/scope.md): parity with the reference kernel's four
modelling modules, and nothing else. Where it stands against that target is
[a committed, machine-checked ledger](kernel/parity-ledger.md) — every
public header of those modules accounted for by a named capability, every
verdict citing symbols and tests the build verifies.

## Where things are

- **This guide** — concepts and workflows, in reading order.
- **[The API reference](https://gilbertorconde.github.io/ogeom-rs/api/ogeom/)** —
  rustdoc for the `ogeom` umbrella crate. Depend on that crate; the crates
  under it are an implementation detail and will move.
- **[The repository](https://github.com/gilbertorconde/ogeom-rs)** — source,
  issues, and the governing documents this book includes verbatim.
