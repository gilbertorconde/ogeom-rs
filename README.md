# ogeom

A boundary-representation CAD kernel, written from scratch in Rust.

**Status: geometry, topology, construction, tessellation, intersection,
booleans, STEP import with healing, blends and chamfers, offsets, shells,
sweeps and lofts are built and tested — milestones M1 through M4. Drawings,
sketching, assemblies and the application-facing layers are not started.** See
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
• crates/ogeom-core        arenas, identity, errors, tolerances, predicates
• crates/ogeom-math        primitives, transforms, B-spline basis, solvers
• crates/ogeom-geom        curves and surfaces behind adaptor traits
• crates/ogeom-topo        the B-rep data model
• crates/ogeom-algo        construction, primitives, measurement, classification
• crates/ogeom-mesh        tessellation, and mesh → B-rep
• crates/ogeom-io          native format, STEP, IGES, glTF, STL, DXF, and the rest
• crates/ogeom-intersect   curve/curve, curve/surface, surface/surface
• crates/ogeom-bool        general fuse and the boolean filters
• crates/ogeom-heal        healing, simplification, canonical recognition
• crates/ogeom-fillet      blends and chamfers
• crates/ogeom-offset      offset, shell, sweep, loft, form features
• crates/ogeom             the public umbrella API
· crates/ogeom-hlr         hidden line removal for drawings
· crates/ogeom-sketch      2D geometric constraint solver
· crates/ogeom-doc         assemblies, appearance, PMI, undo/redo, persistence
· crates/ogeom-select      BVH picking and selection
• tools/ogeom-cli              command-line front end
• tools/ogeom-view             software renderer, so a wrong result is visible
• tools/apisurf            API-usage profiler (analysis only; see docs/SCOPE.md)
```

None of the crates with code is *finished* — each implements the part of
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
follows it, are where these efforts die. This project staged that work behind a
measured benchmark before committing to the boolean pipeline, and the bar was
cleared: booleans run on native and imported geometry, with history, and every
result either measures or refuses with the reason. The staging discipline stays —
each milestone is independently useful, and every known gap is a written entry in
the deferred table rather than a silent wrong answer.

## License

MIT OR Apache-2.0, at your option.
