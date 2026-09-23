# Tessellation and drawings

## Meshing

`ogeom::mesh::tessellate` triangulates a shape at a stated `Deflection` —
the maximum distance the mesh may stand off the exact geometry, with an
angular bound alongside. The triangulation attaches to the model
(`triangulation_of`, `polyline_of` read it back), so a shape carries its
mesh the way it carries its geometry. `simplify` decimates an existing
mesh toward a `Target`; `hatch_face` cross-hatches a face for section
fills.

`triangulate` returns one welded mesh for a whole shape. A viewer that
keeps one mesh per face, to pick and colour faces, meshes them one at a
time — and must still agree with each face's neighbours along every
shared edge, or the pieces crack apart where a narrow face drew its edges
finer than asked. `edge_chords_for` answers that agreement once per shape,
and `triangulate_face_with` draws each face to it; `tessellate` stores its
faces the same way.

Deflection is the honesty parameter throughout: everything downstream
that consumes the mesh — [mass properties](measurement.md), the mesh
exchange formats — inherits exactly the error you chose here, no more
and no less.

## Drawings

`ogeom::hlr` turns solids into 2D drawings the classical way: exact
hidden-line removal, not a rendered picture.

- **`project`** takes shapes and a view direction and returns a
  `Drawing` of `DrawnCurve`s, each tagged with its `Visibility` (visible
  or hidden) and its `Source` — which model edge, silhouette or outline
  produced it. Silhouettes of curved faces are marched on the exact
  surfaces.
- **`section`** cuts a shape with a plane and returns the `SectionView` —
  the cut face outlines ready for hatching; `broken_section` is the
  partial-depth variant.

The `Source` tag on every drawn curve is what makes drawings live: a
dimension attached to a drawn line can find the model edge it measures
after a rebuild, through the same history machinery everything else uses.
