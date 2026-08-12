# ogeom

A boundary-representation CAD kernel, written from scratch in Rust.

**Status:** geometry, topology, construction, tessellation, intersection,
booleans, blends and chamfers, offsets, shells, sweeps and lofts, healing,
product structure, 2D drawings with hidden line removal and sections, and
exchange — STEP and IGES in both directions, STEP carrying assemblies,
colours, semantic PMI and saved views, plus the native format, STL, DXF,
glTF, OBJ, PLY, VRML and 3MF — are built and tested.

**How complete that is against a real kernel is measured, and gated.** The
target is `docs/SCOPE.md`: parity with the reference kernel's modelling
modules, and nothing else. `docs/PARITY.md` audits every public header of
those modules against 97 named capabilities — every verdict citing symbols
and tests the build verifies — and `tools/check.sh` fails if the audit and
the code drift apart. No capability is absent; the `partial` rows say exactly
what restriction each one carries.

**Documentation:** [the guide](https://gilbertorconde.github.io/ogeom-rs/) and
[the API reference](https://gilbertorconde.github.io/ogeom-rs/api/ogeom/).
Every code block in the guide is included from a test this repository runs, so
an example that stops being true fails the build.

## What this is

An independent solid modeling kernel: parametric curves and surfaces, shared
B-rep topology, booleans, blending, offsetting, tessellation, assemblies,
drawing generation and data exchange. A kernel and nothing else — what that
excludes, and why, is `docs/SCOPE.md`.

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

Every crate listed has working code in it.

```
• crates/ogeom-core        arenas, identity, errors, tolerances, predicates
• crates/ogeom-math        primitives, transforms, B-spline basis, solvers
• crates/ogeom-geom        curves and surfaces behind adaptor traits
• crates/ogeom-topo        the B-rep data model
• crates/ogeom-algo        construction, primitives, measurement, classification
• crates/ogeom-mesh        tessellation
• crates/ogeom-io          native format, STEP, IGES, glTF, STL, DXF, OBJ, PLY, …
• crates/ogeom-intersect   curve/curve, curve/surface, surface/surface
• crates/ogeom-bool        general fuse and the boolean filters
• crates/ogeom-heal        healing, simplification, upgrading
• crates/ogeom-fillet      blends and chamfers
• crates/ogeom-offset      offset, shell, sweep, loft, form features
• crates/ogeom-hlr         hidden line removal, sections, 2D drawings
• crates/ogeom-doc         product structure, appearance, semantic PMI
• crates/ogeom             the public umbrella API
• tools/ogeom-cli          command-line front end
• tools/ogeom-view         software renderer, so a wrong result is visible
• tools/ogeom-bench        benchmarks, watched by the gate but not gating
• tools/apisurf            reference-index and API-usage tooling (analysis only)
• tools/parity.py          generates and checks the parity audit
```

Working code that is not kernel lives in `outside/`, which is a separate
workspace the kernel excludes. `docs/SCOPE.md` says what the boundary is and
`outside/README.md` says why each crate is on the far side of it.

```
• outside/crates/ogeom-recognize   feature recognition and process planning
• outside/crates/ogeom-reverse     mesh → B-rep, canonical surface recognition
• outside/crates/ogeom-select      BVH picking, selection, draft and thickness
• outside/crates/ogeom-sketch      2D geometric constraint solver
```

## Building

```sh
cargo build --workspace
cargo test  --workspace
```

The full gate — format, lints, repeated tests, docs, the book and the parity
audit — is one script, and it is what CI runs:

```sh
bash tools/check.sh   # needs mdbook: cargo install mdbook --locked
```

`outside/` is excluded from that workspace and builds on its own:

```sh
cd outside && cargo test --workspace
```

## License

MIT OR Apache-2.0, at your option.
