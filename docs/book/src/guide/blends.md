# Blends and chamfers

Everything here is in `ogeom::fillet`. Edge blends go through one call,
`BlendRequest`: add the edges and what each sweeps, set any options, then
build.

## Rounding one edge

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:fillet_an_edge}}
```

What an edge sweeps is a `BlendSection`:

| Request | Result |
|---|---|
| `.round(edges, r)`, or `.edge(e, BlendSection::Round(r))` | Constant-radius round. |
| `.chain(c, BlendSection::RoundLaw(law))` | Round with a radius law along the edge or chain. |
| `.edge(e, BlendSection::RoundDistances { .. })` | Round meeting each face at its own distance. |
| `.chamfer(edges, d)`, or `.edge(e, Chamfer::Symmetric(d))` | Symmetric flat bevel. |
| `.edge(e, Chamfer::Distances { .. })` | Asymmetric bevel (two distances). |
| `.edge(e, Chamfer::Angle { .. })` | Bevel from a distance and an angle. |

Fillets are not limited to planes. Between cylinders, cones, spheres, tori
and fitted patches, the rolling ball is marched numerically instead of
solved in closed form. The blend surface is fitted through its cross-section
arcs and merged into the solid the same way.

## Several edges at once

Edges blended in separate requests each stop flush against whatever they
end on. Added to one request, they join:

```rust
{{#include ../../../../crates/ogeom/tests/book.rs:edges_asked_together}}
```

- **Rounds** in one request trim neighbouring bands against each other, join
  tangent chains without a seam, and fill every vertex where three or more
  filleted edges meet with a rolling-ball corner patch
  (`.round_corner(..)` chooses the corner).
- **`round_vertex`** is that corner tool on its own:
  - where one ball can touch all the faces, it produces a sphere patch;
  - where no single ball can (for example the apex of a rectangular
    pyramid), it produces the exact envelope of the rolling ball: spheres
    joined by cylinders.
- **Bevels** in one request, from `.chamfer(edges, d)` or one
  `.edge(e, spec)` per edge, are built in one operation. Every wedge is
  built on the solid as it was before the request, so the bevels mitre
  where they meet (`.chamfer_corner(..)` chooses the corner).

The same request also reads, without building: `.limit`, `.largest` and
`.seat` give the largest size up to a bound at which it builds, varying
what `Vary` names (`Vary::Scale` scales every length). `.build_checked`
names every edge it refuses, and `.build_partial` builds the largest subset
that does build.

## Blends between faces without a shared edge

- **`blend_faces`** rolls a constant-radius ball between two faces of one
  solid that need not share an edge. Planes, and curved faces whose
  surfaces share a direction or an axis, get an exact cylinder or torus.
  Curved faces that meet along edges of the solid are blended along those
  edges, marched where the seat has no closed form. Curved faces that share
  no edge and no direction or axis (a B-spline face, a cylinder at a slant
  to a plane) are marched round where their surfaces cross. Where the ball
  leaves one face, the round ends in a flat cap in the ball's section.
- **`march_blend`** is the underlying marcher. It traces the contact circle
  and reports why it stopped (`BlendStop`). Use it to blend up to an
  obstruction on purpose.

## Tangency checks

A blend must end tangent to the faces it joins. `analyse_blend` measures
the achieved contact and returns it as `BlendContact`. Fillets report their
own tangency deviation, and a blend that cannot reach tangency within
tolerance is refused.

Each `BlendContact` carries the gap (G0), the angle between the normals
(G1) and the step in normal curvature square to the edge (G2), all sampled
at stations along the shared edge. `face_curvature_samples` gives the
principal curvatures on a grid over one face, for curvature maps and zebra
stripes.

## Blends: one request

The option-named edge-blend calls are deprecated since 0.22.0 and will be
removed in a later release. Each deprecation note names its exact
replacement. The common paths:

| Before | On the request |
|---|---|
| `fillet_edge(m, s, e, r, t)` | `BlendRequest::new().round([e], r).build(m, s, t)` |
| `fillet_edges(m, s, &e, r, t)` | `BlendRequest::new().round(e, r).build(m, s, t)` |
| `fillet_edges_variable(m, s, &c, &law, t)` | `.chain(c, BlendSection::RoundLaw(law)).build(..)` |
| `fillet_edges_corner(.., corner, ..)` | `.round(e, r).round_corner(corner).build(..)` |
| `fillet_edges_corner_checked(..)` | `.round(e, r).round_corner(corner).build_checked(..)` (see below) |
| `fillet_edge_distances(m, s, e, f, d1, d2, p, t)` | `.edge(e, BlendSection::RoundDistances { face: f, on_face: d1, on_other: d2, profile: p }).build(..)` |
| `blend_edges(m, s, &pairs, corner, t)` | one `.edge(e, blend)` per pair (`Blend` converts), then `.chamfer_corner(corner).build(..)` |
| `chamfer_edge(m, s, e, d, t)` | `BlendRequest::new().chamfer([e], d).build(m, s, t)` |
| `chamfer_edges(m, s, &e, d, t)` | `BlendRequest::new().chamfer(e, d).build(m, s, t)` |
| `chamfer_edges_with(m, s, &specs, t)` | one `.edge(e, spec)` per spec, then `.build(..)` |
| `chamfer_edges_cornered(.., corner, ..)` | the same, then `.chamfer_corner(corner).build(..)` |
| `chamfer_edges_checked(..)` | the same, then `.build_checked(..)` (see below) |
| `largest_fillet_radius_all(m, s, &e, up, t)` | `.round(e, 1.0).largest(m, s, Vary::Scale, up, t)` |
| `largest_fillet_distances(m, s, e, f, d1, d2, up, t)` | `.edge(e, RoundDistances { .. }).largest(.., Vary::Scale, ..)` |
| `largest_chamfer_distance(m, s, e, spec, up, t)` | `.edge(e, spec at measure 1.0).largest(.., Vary::Scale, ..)` |
| `Blend` | `BlendSection` (`From<Blend>`) |

- **The checked builds name edges differently.** The request names every
  edge it refuses, in request order, so the request without them builds.
  The old checked calls named the first culprit and its neighbours.
- **The one-edge calls built alone.** `fillet_edge`, `chamfer_edge` and
  their kin had a one-edge construction; the request builds a single edge
  the way it builds many.
- **Unchanged:** `FilletCorner`, `ChamferCorner`, `RadiusLaw`,
  `RoundProfile`, `Chamfer`, `BlendRefusal`, `fillet_faces`,
  `fillet_sheet_edges`, `blend_faces`, `analyse_blend` and
  `make_blend_surface`.
