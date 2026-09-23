# Exchange

Everything here is in `ogeom::io`.

- Exact formats (STEP, IGES, native) carry whole documents: the model plus
  [product structure, PMI and views](documents.md).
- Mesh formats carry tessellations at the deflection you chose.

## STEP

Read and write, at document level:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:step_roundtrip}}
```

Round trips preserve assemblies with instancing, names, colours, semantic
and presentation PMI, datum systems and saved views.

Bodies come back as:

- `solids`;
- `shells`, for parts exported as faces instead of a solid. A surface model
  stays a shell, under its product like any other body.

`read_step` returns a `StepImport`. Its `report` lists by name every
entity the reader met but did not translate, so nothing is dropped
silently.

### Boundary curves off their surface

Real exports often have boundary curves that sit off the surfaces they
trim.

| Offset | What the reader does |
|---|---|
| Under 1 mm | Heals it: fits the trim, widens the edge's tolerance to the measured offset, and emits a warning. |
| Over 1 mm | Treats the boundary as not describing that surface. The face is read untrimmed and refuses to mesh. It is listed in `report.untrimmed_faces` with its file id and face shape, so you can highlight it or pass it to the healer. `check` also reports these faces as broken. |

`report.summary` groups the per-edge warnings (often thousands) into one
entry per kind: count, worst measured value and an example id. Use it for
a status bar. `warnings` keeps the full text.

### Progress and cancellation

Scope a `Watch` around the call to receive each stage. The readers report
solids as `(done, total)`, so a progress bar can be determinate. The
`Watch`'s canceller stops the work at the next checkpoint:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:watching_an_import}}
```

## IGES

`read_iges` and `write_iges` work at document level like STEP. They cover
the core entity set real files use:

- curve and surface entities, including conic arcs of every kind, ruled
  surfaces, and offset curves and surfaces;
- trimmed surfaces;
- transforms;
- colour;
- the manifold solid B-rep.

`IgesReport` names anything outside that set. IGES round trips are tested
by volume like STEP, including periodic cases (spheres, tori) where seam
handling is error-prone.

## Native format and `.brep`

- `native::write_document` and `native::read_document` round-trip the whole
  document (exact geometry, tolerances, structure, PMI, views, notes) with
  no loss. Use it between ogeom sessions.
- `brep::write` and `brep::read` store a single shape as text, for
  model-level interchange.

## Mesh and drawing formats

| Format | Read | Write |
|---|---|---|
| STL (ascii and binary) | yes | yes |
| glTF / GLB | yes | GLB |
| OBJ | yes | yes |
| PLY | yes | yes |
| VRML | no | yes |
| 3MF (deflated or stored, multi-part) | yes | yes |
| DXF (2D drawings) | yes | yes |

`read_3mf` returns one placed mesh per build item. It flattens components,
follows multi-part packages from the production extension, scales the
model's unit to millimetres, keeps a uniform object colour when the file
has one, and warns about anything read with a caveat.

The mesh writers take the tessellation you built, so the error is the
deflection you chose. DXF is the output for
[HLR drawings](meshing.md#drawings): visible and hidden polylines.

### Meshes to solids

`algo::solid_from_mesh` turns a mesh from any of these formats into a solid
you can model on:

- It builds topology from the mesh's own connectivity.
- It merges coplanar triangles into planar faces. An STL cube comes back as
  six faces.
- It rebuilds regions lying on a cylinder, cone, sphere or torus as that
  surface. A meshed bore becomes a cylinder again, and a meshed ball one
  spherical face.

Its report says where an open mesh is open, and which curved regions could
not be rebuilt exactly and stayed faceted.
