# Measurement and checking

Everything here is in `ogeom::algo`.

## Mass properties

`volume_properties`, `surface_properties` and `linear_properties` return
`MassProperties`: the measure (`mass`), the centroid, and the inertia
tensor about the centroid.

- Area and volume are integrated on the exact surfaces, trim and all,
  whenever every face has pcurves. The result then reports a `deflection`
  of zero and does not depend on the one passed in.
- Otherwise they are computed on a tessellation at the given `Deflection`,
  which the result reports. The answer converges as the deflection gets
  smaller. Lengths are always measured this way.
- An inside-out shell would measure a negative volume. It is reported as an
  error, not returned as a negative number.

## Distances and projections

| Function | Returns |
|---|---|
| `distance_between_shapes` | Minimum distance and the `ClosestPair` that realises it. |
| `project_on_curve`, `project_on_surface`, `project_on_planar_curve` | Parameter and distance of a point's projection onto the geometry. |
| `project_edge_onto_plane` | An edge projected onto a plane as an exact 2D curve: point, line, circle, ellipse, parabola, hyperbola or B-spline (clamped); any other curve fitted. |
| `curve_length` | Length along a curve. |
| `parameter_at_length` | Parameter at a given arc length. |
| `points_by_count`, `points_by_spacing` | Points spread along a curve. |

## Bounds and classification

- `shape_bounds`, `curve_bounds`, `surface_bounds`, `vertex_bounds`:
  axis-aligned bounding boxes.
- `oriented_bounds`: a tight oriented box (`Obb`).
- `classify_in_solid`: is a point inside, outside or on a solid. The
  `_exact` variants work on the exact geometry.
- `classify_on_face`: the same test for a point on a face.

## Validity

`check` verifies the rules the builders promise: edges within their faces'
tolerance, vertices within their edges', wires closed where they claim to
be. It returns a `Diagnosis` of named `Problem`s, each with a `Severity`.
It accepts everything the kernel legitimately builds, so what it flags is
really broken.

`check` also reports a solid whose boundary passes through itself: a face
folded through itself, or two faces crossing away from the edges they
share. The faces' meshes find where they might cross, and each crossing is
confirmed on the surfaces before it is named.

`check_self_intersection` and `check_tessellation` are separate because
they are more expensive. `check_self_intersection` measures every pair of
faces that share no topology for touching, not only for crossing. Call
them when needed.
