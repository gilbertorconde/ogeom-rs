# ogeom

A B-rep solid modelling kernel in Rust. It does not wrap or link any other
CAD kernel, and it is licensed LGPL-2.1 as a derivative work of Open
CASCADE Technology, its reference kernel: see
[`PROVENANCE.md`](PROVENANCE.md).

- Guide: <https://gilbertorconde.github.io/ogeom-rs/>
- API reference: <https://gilbertorconde.github.io/ogeom-rs/api/ogeom/>
- Changes: [`CHANGELOG.md`](CHANGELOG.md)

## What it does

| Area | Included |
|---|---|
| Modelling | primitives, extrusions, revolutions, sweeps, lofts, booleans, fillets and chamfers, offsets and shells, drafts, face removal |
| Geometry | lines, conics, B-splines and NURBS, elementary and freeform surfaces, curve and surface intersection |
| Analysis | volume, area, mass properties, distances, validity checks, healing |
| Meshing | tessellation at a stated deflection, per-face meshes that agree along shared edges |
| Drawings | hidden line removal, sections, DXF output |
| Exchange | STEP (assemblies, colours, PMI, views) and IGES in both directions; STL, OBJ, PLY, glTF, 3MF, VRML, DXF, and a native format |
| Meshes to solids | `solid_from_mesh` turns a triangle mesh into a solid, merging flat regions into planar faces and rebuilding cylinders, cones, spheres, tori, extrusions, surfaces of revolution and fitted B-spline patches |

## Scope and completeness

The target is a complete solid-modelling kernel with its data exchange, plus
turning meshes into solids. [`docs/SCOPE.md`](docs/SCOPE.md) defines the
boundary.

[`docs/CAPABILITIES.md`](docs/CAPABILITIES.md) records 90 named capabilities.
Each verdict cites the code and tests that back it, and the build fails if the
two drift apart. Rows marked `partial` state their limits.

[`PROVENANCE.md`](PROVENANCE.md) records where the code and design came from.

When the kernel cannot do something, it returns an error that says so. It does
not return a guessed result.

## Crates

| Crate | Contents |
|---|---|
| `ogeom` | the public API, re-exporting everything below |
| `ogeom-core` | arenas, identity, errors, tolerances, predicates |
| `ogeom-math` | points, vectors, transforms, B-spline basis, solvers |
| `ogeom-geom` | curves and surfaces |
| `ogeom-topo` | the B-rep data model |
| `ogeom-algo` | construction, measurement, classification, mesh to solid |
| `ogeom-mesh` | tessellation |
| `ogeom-intersect` | curve and surface intersection |
| `ogeom-bool` | booleans and face removal |
| `ogeom-heal` | healing and simplification |
| `ogeom-fillet` | fillets and chamfers |
| `ogeom-offset` | offsets, shells, sweeps, lofts, drafts |
| `ogeom-hlr` | hidden line removal and drawings |
| `ogeom-io` | file formats |
| `ogeom-doc` | assemblies, colours, names, PMI |

`tools/` holds a command-line tool, a software renderer, benchmarks and the
capability-ledger scripts. `outside/` holds working code that is out of scope (constraint
solving, feature recognition, picking) in a separate workspace.

## Design

- A shape is a handle: a topology node, a chain of placements and an
  orientation. The same node can be placed many times without copying.
- Every entity keeps its own tolerance, which only grows.
- Every entity records where it came from, so references survive a rebuild.
- Every operation returns a history of what became of its inputs.

[`docs/DATA_MODEL.md`](docs/DATA_MODEL.md) has the details.

## Building

```sh
cargo build --workspace
cargo test --workspace
```

The full check that CI runs (formatting, lints, tests, docs, the guide and the
capability ledger):

```sh
bash tools/check.sh   # needs mdbook: cargo install mdbook --locked
```

### In your own crate

The workspace's release profile does not reach a crate that depends on
ogeom, and neither does a target CPU. Two settings in the dependent crate
make the kernel faster with results the same to the bit:

```toml
# Cargo.toml of the crate that depends on ogeom
[profile.release]
lto = "thin"
codegen-units = 1
```

```toml
# .cargo/config.toml, where every machine it runs on has AVX2 and FMA
# (x86-64-v3: Intel Haswell, AMD Excavator and later)
[build]
rustflags = ["-C", "target-cpu=x86-64-v3"]
```

Without a target CPU each of the kernel's fused multiply-adds is a call
into the C library rather than one instruction. Both are rounded once, so
the answers do not change. Measured converting and refining a 29,003-face
solid on one thread, against a dependent crate's default release build:

| Settings | Instructions | CPU time |
| --- | --- | --- |
| default | 9.69e11 | 60.6 s |
| thin LTO, one codegen unit | 9.39e11 (-3%) | 57.1 s (-6%) |
| `x86-64-v3` | 7.44e11 (-23%) | 46.9 s (-23%) |
| both | 7.22e11 (-26%) | 45.7 s (-25%) |

The refined solid, and each of `mesh_corpus`'s fifteen parts converted
back, came out the same to the bit under all four.

## License

The GNU Lesser General Public License, version 2.1, and no other version:
`LGPL-2.1-only`. The licence text is [`COPYING.LESSER`](COPYING.LESSER) and
the GPL it refers to is [`COPYING`](COPYING).

ogeom is a derivative work of Open CASCADE Technology in the licensing
sense: it was built with OCCT as its reference kernel, and it follows
OCCT's public API vocabulary, data-model conventions and default
tolerances. [`PROVENANCE.md`](PROVENANCE.md) records the full provenance
and the audit behind that statement.
