# Exchange

Everything here is `ogeom::io`. The exact formats carry documents — model
plus [product structure, PMI, views](documents.md) — because that is what
the files actually contain; the mesh formats carry tessellations at the
deflection you chose.

## STEP

Both directions, document-level:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:step_roundtrip}}
```

Assemblies with instancing, names, colours, semantic and presentation
PMI, datum systems and saved views all survive the trip. `read_step`
returns a `StepImport` whose `report` lists, by name, every entity the
reader met and did not translate — the file's inventory of what was and
was not understood, instead of a silent partial import.

Real exports carry slop: boundary curves that sit off the surfaces they
trim. Under a millimetre the reader heals it — the trim is fitted, the
edge's tolerance widens to the measured offset, and a warning says so —
because that is the file's own error, honestly carried. Beyond a
millimetre the boundary is not describing that surface at all, and a
fitted trim would be invented geometry drawn with confidence: the face
reads untrimmed, refuses to mesh, and `report.untrimmed_faces` carries it
— file id and face shape both — so a consumer can mark the exact gap or
hand the face straight to the healer. `check` reports the
same faces as broken from the model side.

A large import is worth watching. A `Watch` scoped around the call hears
each stage — the readers announce their solids as `(done, total)`, so a
progress bar can be determinate — and its canceller stops the work at
the next checkpoint:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:watching_an_import}}
```

## IGES

`read_iges` / `write_iges`, same document-level shape, covering the core
entity set a real importer meets: the curve and surface entities, trimmed
surfaces, transforms, colour, and the manifold solid B-rep. The
`IgesReport` names what fell outside that set. The suite holds IGES to
the same standard as STEP: round trips measured by volume, including the
periodic cases (spheres, tori) where seam handling is where importers
usually break.

## The native format and `.brep`

`native::write_document` / `native::read_document` round-trip the entire
document — exact geometry, tolerances, structure, PMI, views, notes —
with no translation loss; it is the format to use between ogeom sessions.
`brep::write` / `brep::read` carry a single shape in a text form, for
interchange at the model level.

## Mesh and drawing formats

| Format | Read | Write |
|---|---|---|
| STL (ascii and binary) | yes | yes |
| glTF / GLB | yes | GLB |
| OBJ | — | yes |
| PLY | — | yes |
| VRML | — | yes |
| 3MF | — | yes |
| DXF (2D drawings) | yes | yes |

The mesh writers take the tessellation you built at your chosen
deflection — the error budget is yours, stated once. DXF is the outlet
for [HLR drawings](meshing.md#drawings): visible and hidden polylines,
ready for a title block.
