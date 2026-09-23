# Introduction

**ogeom** is a boundary-representation (B-rep) CAD kernel written from
scratch in Rust. It does not wrap, vendor or link any existing CAD kernel.

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
- **Completeness is audited.** [Scope](kernel/scope.md) sets the target:
  parity with the reference kernel's four modelling modules, plus meshes
  into solids. [The parity ledger](kernel/parity-ledger.md) is committed
  and machine-checked: it maps every public header of those modules to a
  named capability, and each verdict cites symbols and tests the build
  verifies.

## Where things are

- **This guide:** concepts and workflows, in reading order.
- **[The API reference](https://gilbertorconde.github.io/ogeom-rs/api/ogeom/):**
  rustdoc for the `ogeom` umbrella crate. Depend on that crate. The crates
  under it are an implementation detail and will change.
- **[The repository](https://github.com/gilbertorconde/ogeom-rs):** source,
  issues, and the governing documents this book includes verbatim.
