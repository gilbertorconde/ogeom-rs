# Healing

Imported geometry is often imperfect: gaps between faces, edges whose 3D
curve and surface curves (pcurves) disagree, one surface split into many
patches. `ogeom::heal` repairs what it can within tolerance and reports
what it did. Each repair measures what it achieved. A shape that cannot be
repaired within the stated tolerance comes back with a named diagnosis.

| Function | What it does |
|---|---|
| `sew` (in `ogeom::algo`) | Stitches faces into shells by merging boundaries that coincide within tolerance. Uses each entity's own tolerance: a vertex merges with what lies inside its tolerance, and the survivor widens to cover what it absorbed. |
| `repair_same_parameter` | Re-fits an edge's pcurves until they agree with its 3D curve within tolerance. The report gives the achieved deviation per edge. |
| `unify_same_domain` | Merges adjacent faces on the same surface and adjacent edges on the same curve. Undoes the fragmentation that exchange formats cause. |
| `merge_edges` | Joins chains of edges into single edges where the geometry allows. |
| `reduce_tolerances` | Narrows entity tolerances to what the geometry actually measures. This is the only operation that narrows tolerances, and it does so by re-measuring. |
| `canonical_simplify` | Replaces NURBS geometry that is exactly analytic with the analytic form (for example a plane stored as a bicubic patch, or a circle stored as a rational spline). Each match is verified at every sample, and the report includes the worst deviation. A surface that is only almost a cylinder stays a spline. |
| `reanchor_periodic_rings` | Moves the seam of periodic faces so later algorithms see a consistent parametrisation. |
| `fix_shape` | One-call repair for a shape of unknown quality (see below). Returns a `FixReport`. |
| `Reshape` | The primitive underneath: a recorded substitution of entities that rebuilds everything referencing them. |

## `fix_shape` steps

1. Diagnose.
2. Put each wire's edges end to end.
3. Collapse edges shorter than their own vertices' tolerances.
4. Fit missing pcurves.
5. Sew loose faces.
6. Tighten tolerances.
7. Widen any vertex tighter than its edges (and any edge tighter than its
   faces) to restore the [containment rule](tolerances.md#3-the-containment-rule).
8. Diagnose again.

The `FixReport` lists what it did and what the checker still finds. Small
faces and small solids are left in place, because removing them is
defeaturing (see `remove_faces` in [Booleans](booleans.md)).
