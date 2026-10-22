# Introduction

**ogeom** is a boundary-representation (B-rep) CAD kernel in Rust. It does
not wrap, vendor or link any existing CAD kernel, and it is licensed
LGPL-2.1 as a derivative work of Open CASCADE Technology, its reference
kernel ([provenance](kernel/provenance.md)).

## What it covers

- Parametric curves and surfaces.
- Shared B-rep topology with per-entity tolerances and operation history.
- Booleans, blends and chamfers, offsets and shells, sweeps and lofts.
- Healing and tessellation.
- Product structure (assemblies) and 2D drawings with hidden-line removal.
- Data exchange:
  - STEP, read and write, with assemblies, colours, semantic PMI and saved
    views.
  - IGES, read and write.
  - The native document format.
  - Mesh formats. Meshes can be turned back into solids, with their curved
    regions recognized.

## How this book is kept honest

- **Examples are tested.** Every code block in this guide comes from
  [a test file](https://github.com/gilbertorconde/ogeom-rs/blob/main/crates/ogeom/tests/book.rs)
  that runs in CI, so a broken or false example fails the build. The test
  suite checks every numeric claim against a closed form or a round trip.
- **Completeness is recorded.** [Scope](kernel/scope.md) sets the target:
  a complete solid-modelling kernel with its data exchange, plus meshes
  into solids. [The capability ledger](kernel/capabilities.md) is committed
  and machine-checked: each verdict cites symbols and tests the build
  verifies.

## Where things are

- **This guide:** concepts and workflows, in reading order.
- **[The API reference](https://gilbertorconde.github.io/ogeom-rs/api/ogeom/):**
  rustdoc for the `ogeom` umbrella crate. Depend on that crate. The crates
  under it are an implementation detail and will change.
- **[The repository](https://github.com/gilbertorconde/ogeom-rs):** source,
  issues, and the governing documents this book includes verbatim.
