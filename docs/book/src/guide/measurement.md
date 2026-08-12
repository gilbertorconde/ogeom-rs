# Measurement and checking

Everything here is `ogeom::algo`.

## Mass properties

`volume_properties`, `surface_properties` and `linear_properties` return
`MassProperties` — the measure itself (`mass`), the centroid, and the
inertia tensor about it. They are computed on a tessellation at a stated
`Deflection`, and that is worth knowing: the answer converges as the
deflection tightens, and the suite's own volume assertions carry bounds
derived from exactly that. An inside-out shell measures negative and is
reported as the error it is, not as a negative volume.

## Distances and projections

`distance_between_shapes` returns the minimum distance and the
`ClosestPair` realising it. `project_on_curve`, `project_on_surface` and
`project_on_planar_curve` drop a point onto geometry and report the
parameter and distance. `curve_length` and `parameter_at_length` measure
along curves; `points_by_count` / `points_by_spacing` discretise them.

## Extent and classification

`shape_bounds`, `curve_bounds`, `surface_bounds` and `vertex_bounds`
return axis-aligned boxes; `oriented_bounds` fits the tight `Obb`.
`classify_in_solid` answers in / out / on for a point (with `_exact`
variants that work on the exact geometry), `classify_on_face` the same
on a surface.

## Validity

`check` runs the validity analysis — the same containment rules the
builders promise, verified: edges within their faces' tolerance, vertices
within their edges', wires closed where they claim to be. It returns a
`Diagnosis` of named `Problem`s with `Severity`, and it is calibrated to
the builders — what the kernel legitimately builds, the checker accepts;
what it flags is genuinely broken. `check_self_intersection` and
`check_tessellation` are the deeper, more expensive questions, separate
because you should get to choose when to pay for them.
