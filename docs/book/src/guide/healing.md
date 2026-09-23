# Healing

Imported geometry arrives imperfect: gaps between faces, edges whose 3D
curve and surface curves disagree, one surface split into a dozen
patches. `ogeom::heal` repairs what can be repaired honestly and reports
what it did.

- **`sew`** (in `ogeom::algo`) stitches faces into shells, merging
  coincident-within-tolerance boundaries. It honours each entity's stated
  tolerance — a vertex whose tolerance says "I am this uncertain" merges
  with what falls inside that reach, and the survivor widens to cover
  what it absorbed.
- **`repair_same_parameter`** re-fits an edge's surface curves until 3D
  curve and pcurves agree within tolerance, and its report says what
  deviation was achieved — measured, per edge.
- **`unify_same_domain`** merges adjacent faces that lie on the same
  surface and adjacent edges on the same curve — undoing the
  fragmentation exchange formats inflict.
- **`merge_edges`** joins chains of edges into single edges where the
  geometry allows it.
- **`reduce_tolerances`** narrows entity tolerances to what the geometry
  actually measures — the only honest direction-reversal in the tolerance
  model, because it re-derives the values instead of asserting them.
- **`canonical_simplify`** recognises exact analytic geometry hiding in
  NURBS clothing — the plane that came through exchange as a bicubic
  patch, the circle written as a rational spline — and substitutes the
  canonical form. Recognition is verified at every sample and the report
  carries a worst-deviation certificate; a surface that is *almost* a
  cylinder stays a spline.
- **`reanchor_periodic_rings`** re-anchors the seam of periodic faces so
  downstream algorithms see a consistent parametrisation.
- **`fix_shape`** is the one-call pass over a shape nobody promised was
  well-formed: it diagnoses, puts a wire's edges end to end, collapses
  edges shorter than their own vertices' tolerances, fits missing
  pcurves, sews loose faces, tightens tolerances, widens any vertex left
  tighter than its edges (and any edge than its faces) to restore the
  containment rule, and diagnoses again.
  Its `FixReport` says what it did and what the checker still sees;
  small faces and solids are left in place, since removing them is
  defeaturing.
- **`Reshape`** is the primitive underneath: a recorded substitution of
  entities that rebuilds everything referencing them.

The philosophy is the kernel's general one: each repair measures what it
achieved rather than declaring success, and a shape that cannot be
repaired within stated tolerance comes back with a named diagnosis, not a
silent pass.
