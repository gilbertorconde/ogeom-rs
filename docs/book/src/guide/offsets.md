# Offsets, shells and features

Everything here is in `ogeom::offset`.

## Offsetting

| Function | What it does |
|---|---|
| `offset_shape` | Moves a solid's boundary along its normals (outward, or inward with a negative distance) and rebuilds the intersections where offset surfaces collide. |
| `offset_faces` | Offsets only the named faces, and the faces around them follow: pull a wall out, grow a bore, raise a boss. The topology stays the same, and the history maps each face to its result. |
| `move_faces` | Moves the named faces by a translation or rotation, the faces around them following, as `offset_faces` does. |
| `make_thick_solid` | Shelling: removes the named faces, offsets the rest inward and joins them. Gives a hollow part with an opening. |
| `make_thick_solid_sided` | Shelling with the wall on either side of the faces or both (`Walls::both`, `Walls::split(behind, ahead)`). Each opening leaves one flat rim on the removed face's surface. |
| `offset_sheet` | Moves a face or a shell along its normals, trim and joins kept. A face with no parallel surface in its own family (a B-spline) is fitted, and its tolerance holds the measured deviation. |
| `make_thick_sheet` | Thickens a face or an open shell into a solid: the sheet and its offset (or half the thickness each way), closed by side faces along the free edges. Faces meeting at a crease join on the mitre. |
| `make_thick_sheet_sided` | `make_thick_sheet` with its own depth to each side of the sheet. |
| `offset_wire` | 2D offset of a planar wire. The caller picks the `Join` style (arcs or intersections). Useful for tool-path-like outlines. |
| `thicken_open_wire` | The closed outline of an open planar chain of lines and arcs thickened between two signed offsets to its left, ends closed by a `Cap`. |
| `apply_draft` | Tilts faces by a draft angle about a neutral plane, for moulded parts. |

### Where a wire offset collapses

An inward offset that runs past a narrow part of the profile crosses
itself. `offset_wire` follows the common CAD rule for this: split the raw
offset where it crosses itself and drop the loops that are invalid. The
result is the boundary of the region at the offset distance or more from
the wire. In detail:

- A piece that comes closer to the wire than the offset is dropped.
- Two surviving pieces on the same line or circle that run the same way
  are one stretch of the boundary, so it is kept once. This happens when
  a neck closes and the walls on either side of it line up.
- Two surviving pieces on the same span that run opposite ways enclose
  nothing, so both are dropped. This is a neck exactly twice the offset
  wide.
- When the collapse cuts the region apart, the result is a compound with
  one wire per island.

### Sharp corners on curved sides

`Join::Intersection` carries both pieces on along their own line or
circle to their meeting nearest the corner. A line and a circle, or two
circles, may not meet ahead of the corner (each piece extended by less
than half a turn); then the join is refused by name and the arc join is
the way out.

### Thickening an open chain

`thicken_open_wire(model, chain, plane, left, right, cap, join, tol)`
builds the band between two offsets of an open chain. Both offsets are
signed to the chain's left about the plane's normal, so `(2, 0)` keeps
the chain as the right side and `(1, -1)` straddles it. Each side is the
chain's own offset: the parallel line, the concentric arc. A side on the
outside of a turn is joined across the gap by `join`; a side on the
inside is trimmed to where its pieces cross. The ends take a `Cap`:

- `Flat`: a line square to the chain through its end point.
- `Round`: a half circle across the band's width, past the end.
- `Extended`: both sides carried on past the end by half the band's
  width, then joined square across.

Refused by name: a closed wire, a band of no width, a side that reaches
past an arc's centre (an annular sector needs an inner radius), a chain
that turns straight back on itself, and an outline that crosses itself
(a thickness too great for the chain's turns).

## Sweeps

The general sweeps share this module's machinery:

- `make_pipe` and `make_pipe_skinned`: sweep along a single spine edge.
- `make_loft` and `make_loft_skinned`: loft through profile sections.
- `make_evolved`: sweep a profile along a planar spine.
- `make_filling`: fill the loop four edges bound with a fitted patch face.
- `make_filling_n`: fill a hole bounded by any number of edges with one
  face bounded by those same edges, meeting each neighbouring face at G0,
  G1 or G2 and passing through points and curves inside the hole. The
  gaps, angles and curvature differences it reached come back measured
  per side.

## Features

Feature operations combine a sketch with a solid in one step:

- `feature_prism`: bosses and pockets.
- `feature_revol`: revolved bosses and grooves.
- `feature_rib` and `feature_slot`.
- `make_rib`: a rib from an open chain, grown until the solid stops it.

`make_rib(model, solid, chain, plane, thickness, along_normal, tol)`
follows the usual rib feature:

1. Each end of the chain is carried straight on along its tangent into
   the wall it points at. An end that points at nothing stays put.
2. The web is thickened across the plane (`along_normal` false:
   `thickness.0` along the normal, `thickness.1` against it) and grows in
   the plane square to the chain's chord. Or it is thickened in the plane
   (`along_normal` true: `thickness.0` to the chain's left, `thickness.1`
   to its right, flat ends, sharp corners) and grows along the normal.
3. The growth is tried both ways. A way the solid does not close (the
   web runs out past the solid) meets no wall. Where both ways close,
   the smaller web is the rib. Where neither does, the rib is refused by
   name.
4. The web's sides are the lines through the chain's (carried) ends along
   the growth, so a web next to a wall shorter than the chain ends square
   there instead of running round the wall's end.

Each is a constrained boolean internally and records history like every
other operation.

`normal_projection` projects a wire onto a shape along the shape's normals
(for engraving). It returns the `Projected` curves on the target faces.

`split_face` (in `ogeom::heal`) cuts one face along curves that lie on it
(`Projection::OnFace`) or are dropped onto it along its normals
(`AlongNormals`) or along a direction (`Along`). Each curve must cross the
face from boundary to boundary or close inside it. The pieces share the new
edges, and the history maps the face to its pieces.

A stretch of a curve along the face's boundary, or along an earlier cut,
divides nothing and is passed over. An entry of `edges` that only runs
outside the face or along its boundary misses it, and is refused. A curve
that crosses the face anywhere does not miss it, even where part of it runs
outside.

`split_face_with` takes `SplitOptions`. With `skip_misses`, a curve that
misses the face (its projection never lands, or its trace lies outside the
face or along its boundary) is left out and its index returned in
`Split::missed`, instead of refusing the split. Each
`PieceOrigin` in `Split::pieces` lists the entries of `edges` that close into
a loop round the piece, so a sketch region can be matched to its piece
without sampling it.
