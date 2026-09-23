# Tessellation and drawings

## Meshing

`ogeom::mesh::tessellate` triangulates a shape at a given `Deflection`:
the maximum distance between the mesh and the exact geometry, plus an
angular bound.

- The triangulation is stored on the model. Read it back with
  `triangulation_of` and `polyline_of`.
- `triangulate` returns one welded mesh for a whole shape.
- `simplify` decimates an existing mesh toward a `Target`.
- `hatch_face` cross-hatches a face, for section fills.

### One mesh per face

A viewer that picks and colours individual faces needs one mesh per face.
Adjacent face meshes must use the same points along every shared edge, or
cracks appear (for example where a narrow face samples its edges more
finely than requested).

1. Call `edge_chords_for` once per shape to fix the shared edge points.
2. Call `triangulate_face_with` for each face using those points.

`tessellate` stores its faces the same way.

### Deflection controls accuracy downstream

Everything that consumes the mesh ([mass properties](measurement.md), the
mesh exchange formats) carries exactly the error you chose here.

## Drawings

`ogeom::hlr` produces 2D drawings by exact hidden-line removal (not a
rendered image).

- **`project`** takes shapes and a view direction and returns a `Drawing`
  of `DrawnCurve`s. Each curve has:
  - a `Visibility` (visible or hidden);
  - a `Source`: the model edge, silhouette or outline that produced it.

  Silhouettes of curved faces are traced on the exact surfaces.
- **`section`** cuts a shape with a plane and returns a `SectionView` with
  the cut face outlines, ready for hatching.
- **`broken_section`** is the partial-depth variant.

Because each drawn curve has a `Source`, a dimension attached to a drawn
line can find the model edge it measures after a rebuild, using the same
history as every other operation.
