# Getting started

The kernel is not on crates.io yet. Depend on it by git:

```toml
[dependencies]
ogeom = { git = "https://github.com/gilbertorconde/ogeom-rs" }
```

Use the `ogeom` umbrella crate. It re-exports the whole API as modules
(`ogeom::algo`, `ogeom::boolean`, `ogeom::topo`, `ogeom::io`, ...). The
`ogeom-*` crates underneath are an implementation detail and their
boundaries will change.

## A first solid

A block with a hole through it: two primitives, one boolean, and a volume
check.

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:first_solid}}
```

The same patterns apply to the whole API:

- **`Model` owns all data.** Geometry, topology, tolerances and history
  live in one [`Model`]. Operations take `&mut model`. A [`Shape`] is a
  cheap handle into the model; copying it copies no geometry. See
  [the data model](data-model.md).
- **Every operation takes a `Tolerances`.** There is no global epsilon.
  `Tolerances::millimetres()` is the preset for models in millimetres. See
  [Tolerances](tolerances.md).
- **Every operation returns a `Built`.** `built.shape` is the result. The
  rest of `Built` is the history: which input entities generated or were
  modified into which outputs. Every operation fills it in, so parametric
  applications can rely on it.
- **Errors are values.** Everything returns `Result`. When the kernel
  cannot produce a correct result, it returns a
  [named refusal](refusals.md) instead of bad geometry.

[`Model`]: https://gilbertorconde.github.io/ogeom-rs/api/ogeom/topo/struct.Model.html
[`Shape`]: https://gilbertorconde.github.io/ogeom-rs/api/ogeom/topo/struct.Shape.html
