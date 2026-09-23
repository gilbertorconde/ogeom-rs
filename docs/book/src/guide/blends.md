# Blends and chamfers

Everything here is `ogeom::fillet`.

## Rounding an edge

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:fillet_an_edge}}
```

`fillet_edge` rounds at constant radius; `fillet_edge_variable` takes a
radius law along the edge. `chamfer_edge` cuts a symmetric flat,
`chamfer_edge_distances` an asymmetric one, and `chamfer_edge_angle` a
distance-and-angle one.

A fillet is not limited to planes. Between cylinders, cones, spheres, tori
and fitted patches the ball is marched rather than solved in closed form,
the blend is fitted through its own arcs, and the result melts into the
solid the same way.

## Edges asked together

Edges rounded or bevelled one call at a time each stop flush against
whatever they end on. Asked together, they meet:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:edges_asked_together}}
```

`fillet_edges` trims neighbouring bands against each other, joins a
tangent chain without a seam, and closes any vertex where three or more of
its edges meet with the rolling ball's own patch. `round_vertex` is that
corner tool on its own: at a vertex one ball touches it leaves a sphere
patch, and at one no single ball touches — a rectangular pyramid's apex —
the exact envelope of the rolling ball, spheres joined by cylinders.
`chamfer_edges` and `chamfer_edges_with` bevel a set of edges as one
operation, every wedge built on the solid as it stands before the call, so
the bevels mitre where they meet.

## Blends without a shared edge

`blend_faces` rolls a constant-radius ball between two faces that need
not share an edge at all — the general face–face blend. `march_blend` is
the machinery underneath, exposed: it marches the contact circle and
reports how it stopped (`BlendStop`), which callers can use to blend up
to an obstruction deliberately.

## Honesty at the tangent line

A blend's job is to end tangent to the faces it joins, and near-tangency
is where blend algorithms traditionally lie. `analyse_blend` measures the
achieved contact (`BlendContact`) instead of asserting it: the fillet
reports its own tangency deviation, and a blend that cannot achieve
tangency within tolerance is refused rather than delivered looking
smooth.
