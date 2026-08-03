# openGeometry

A boundary-representation CAD kernel, written from scratch in Rust.

**Status: the foundation, the geometry, the B-rep model, construction and
tessellation are built and tested. Intersection, booleans and everything that
depends on them are not started.** See
[Milestones](docs/SCOPE.md#milestones) for what each of those means and how it
is decided.

## What this is

A complete, independent solid modeling kernel: parametric curves and surfaces,
shared B-rep topology, booleans, blending, offsetting, tessellation, sketching,
assemblies, drawing generation and data exchange.

It is not a binding, a wrapper, or a port. It depends on no existing CAD kernel,
vendors none, and links against none.

The scope is set by what a CAD kernel has to do — see [`docs/SCOPE.md`](docs/SCOPE.md).
Applications are consumers of the kernel, not definitions of it. If some capability
happens to be one that a given application implements for itself, that is not a
reason to leave it out.

## The data model

[`docs/DATA_MODEL.md`](docs/DATA_MODEL.md) is normative — everything in `crates/`
implements it. Twelve invariants, each cheap to honour now and effectively
impossible to retrofit across a kernel's worth of algorithms later:

- a shape is a `(topology node, location, orientation)` triple, cheap to copy;
- location is a **chain** of transforms, not a flat matrix, so instancing works;
- orientation composes multiplicatively on descent;
- identity is a trichotomy — same / equal / partner — each with a matching hasher;
- tolerances are **per entity** and only ever widen;
- an edge carries a *list* of representations, including one pcurve per adjacent face;
- every operation emits history.

Two places where the conventional design is wrong and this one diverges:

- **Stable entity identity.** Kernels that identify topology by pointer produce new
  entities on every operation, so downstream references die — that *is* the
  topological naming problem, and applications have spent years working around it.
  Here, provenance is recorded at creation.
- **Predicates behind a trait**, so the robustness strategy can change without
  rewriting algorithms.

## Layout

`•` has working code; `·` is a crate that exists so its dependents can name it,
and holds nothing yet.

```
• crates/og-core        arenas, identity, errors, tolerances, predicates
• crates/og-math        primitives, transforms, B-spline basis, solvers
• crates/og-geom        curves and surfaces behind adaptor traits
• crates/og-topo        the B-rep data model
• crates/og-algo        construction, primitives, measurement, classification
• crates/og-mesh        tessellation, and mesh → B-rep
• crates/og-io          native format, STEP, IGES, glTF, STL, DXF, and the rest
• crates/og             the public umbrella API
· crates/og-intersect   curve/curve, curve/surface, surface/surface
· crates/og-bool        general fuse and the boolean filters
· crates/og-heal        healing, simplification, canonical recognition
· crates/og-fillet      blends and chamfers
· crates/og-offset      offset, shell, sweep, loft, form features
· crates/og-hlr         hidden line removal for drawings
· crates/og-sketch      2D geometric constraint solver
· crates/og-doc         assemblies, appearance, PMI, undo/redo, persistence
· crates/og-select      BVH picking and selection
• tools/ogcli           command-line front end
• tools/ogview          software renderer, so a wrong result is visible
• tools/apisurf         API-usage profiler (analysis only; see docs/SCOPE.md)
```

None of the eight with code is *finished* — each implements the part of
[`docs/SCOPE.md`](docs/SCOPE.md) its milestone called for and no more. What is
owed and where it lands is in that document's deferred table, so a gap is a
scheduled debt rather than a surprise found by whoever hits it.

## Building

```sh
cargo build --workspace
cargo test  --workspace
```

No C toolchain, no system libraries, no submodules.

## Prior art, honestly

Building a B-rep kernel is hard in a way that is easy to underestimate.
[Fornjot](https://github.com/hannobraun/fornjot) — a funded, full-time, greenfield
Rust attempt — was archived in June 2026 having never reached a working boolean.
[truck](https://github.com/ricosjp/truck) has real NURBS and topology after years
of corporate backing, but its boolean surface is two functions and its STEP support
is export-only. The one demonstrably successful greenfield Rust kernel is
closed-source and venture-funded.

Surface/surface intersection, and the tolerance-managed topology rebuild that
follows it, are where these efforts die. This project is staged so every milestone
is independently useful, and so the intersection work is gated on a measured
benchmark before the boolean pipeline is committed to. If it does not clear the
bar, this ships as a geometry library instead of consuming years.

## License

MIT OR Apache-2.0, at your option.
