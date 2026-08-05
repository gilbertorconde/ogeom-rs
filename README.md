# ogeom

A boundary-representation CAD kernel, written from scratch in Rust.

**Status: geometry, topology, construction, tessellation, intersection,
booleans, STEP import with healing, blends and chamfers, offsets, shells,
sweeps and lofts are built and tested. Drawings, sketching, assemblies and the
application-facing layers are not started.**

## What this is

A complete, independent solid modeling kernel: parametric curves and surfaces,
shared B-rep topology, booleans, blending, offsetting, tessellation, sketching,
assemblies, drawing generation and data exchange.

It is not a binding, a wrapper, or a port. It depends on no existing CAD kernel,
vendors none, and links against none.

## The data model

Every crate implements one shared B-rep model. Geometry (curves, surfaces) and
topology (vertices, edges, faces, shells, solids) live in shared arenas; a
shape is a cheap handle into them, and the same node can appear many times —
placed, mirrored, instanced — without being copied. The model rests on a small
set of invariants, each cheap to honour from the start and effectively
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
• tools/ogeom-cli          command-line front end
• tools/ogeom-view         software renderer, so a wrong result is visible
• tools/apisurf            API-usage profiler (analysis only)
```

## Building

```sh
cargo build --workspace
cargo test  --workspace
```

## License

MIT OR Apache-2.0, at your option.
