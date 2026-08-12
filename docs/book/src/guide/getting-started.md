# Getting started

The kernel is not on crates.io yet; depend on it by git:

```toml
[dependencies]
ogeom = { git = "https://github.com/gilbertorconde/ogeom-rs" }
```

Depend on `ogeom`, the umbrella crate. It re-exports the whole API as
modules — `ogeom::algo`, `ogeom::boolean`, `ogeom::topo`, `ogeom::io` and
so on. The `ogeom-*` crates underneath are an implementation detail; their
boundaries will move, and the umbrella is what stays put.

## A first solid

A block with a hole through it: two primitives, one boolean, and a
measurement — because a result you have not measured is a result you are
assuming.

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:first_solid}}
```

Everything in that example generalises:

- **`Model` owns everything.** Geometry, topology, tolerances, history —
  all entities live in one [`Model`], and operations take `&mut model`.
  A [`Shape`] is a cheap handle into it, not the data itself; copying one
  copies nothing. The [data model](data-model.md) chapter is the full
  story.
- **Every operation takes a `Tolerances`.** There is no global epsilon.
  `Tolerances::millimetres()` is the preset for models measured in
  millimetres; [Tolerances](tolerances.md) explains what the values mean
  and the rules they obey.
- **Every operation returns a `Built`.** `built.shape` is the result;
  the rest of `Built` is the operation's history — which input entities
  generated or were modified into which outputs — and it is populated by
  every operation, always. Parametric applications are built on that
  promise.
- **Failures are values.** Everything returns `Result`. When the kernel
  cannot do something honestly, it [refuses by name](refusals.md) rather
  than producing garbage geometry.

[`Model`]: https://gilbertorconde.github.io/ogeom-rs/api/ogeom/topo/struct.Model.html
[`Shape`]: https://gilbertorconde.github.io/ogeom-rs/api/ogeom/topo/struct.Shape.html
