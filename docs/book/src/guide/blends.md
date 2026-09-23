# Blends and chamfers

Everything here is in `ogeom::fillet`.

## Rounding one edge

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:fillet_an_edge}}
```

| Function | Result |
|---|---|
| `fillet_edge` | Constant-radius round. |
| `fillet_edge_variable` | Round with a radius law along the edge. |
| `chamfer_edge` | Symmetric flat bevel. |
| `chamfer_edge_distances` | Asymmetric bevel (two distances). |
| `chamfer_edge_angle` | Bevel from a distance and an angle. |

Fillets are not limited to planes. Between cylinders, cones, spheres, tori
and fitted patches, the rolling ball is marched numerically instead of
solved in closed form. The blend surface is fitted through its cross-section
arcs and merged into the solid the same way.

## Several edges at once

Edges filleted or chamfered in separate calls each stop flush against
whatever they end on. Passed in one call, they join:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:edges_asked_together}}
```

- **`fillet_edges`** trims neighbouring bands against each other, joins
  tangent chains without a seam, and fills every vertex where three or more
  filleted edges meet with a rolling-ball corner patch.
- **`round_vertex`** is that corner tool on its own:
  - where one ball can touch all the faces, it produces a sphere patch;
  - where no single ball can (for example the apex of a rectangular
    pyramid), it produces the exact envelope of the rolling ball: spheres
    joined by cylinders.
- **`chamfer_edges`** and **`chamfer_edges_with`** bevel a set of edges in
  one operation. Every wedge is built on the solid as it was before the
  call, so the bevels mitre where they meet.

## Blends between faces without a shared edge

- **`blend_faces`** rolls a constant-radius ball between two faces that
  need not share an edge.
- **`march_blend`** is the underlying marcher. It traces the contact circle
  and reports why it stopped (`BlendStop`). Use it to blend up to an
  obstruction on purpose.

## Tangency checks

A blend must end tangent to the faces it joins. `analyse_blend` measures
the achieved contact and returns it as `BlendContact`. Fillets report their
own tangency deviation, and a blend that cannot reach tangency within
tolerance is refused.
