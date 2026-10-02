# Offsets, shells and features

Everything here is in `ogeom::offset`.

## Offsetting

| Function | What it does |
|---|---|
| `offset_shape` | Moves a solid's boundary along its normals (outward, or inward with a negative distance) and rebuilds the intersections where offset surfaces collide. |
| `offset_faces` | Offsets only the named faces, and the faces around them follow: pull a wall out, grow a bore, raise a boss. The topology stays the same, and the history maps each face to its result. |
| `move_faces` | Moves the named faces by a translation or rotation, the faces around them following, as `offset_faces` does. |
| `make_thick_solid` | Shelling: removes the named faces, offsets the rest inward and joins them. Gives a hollow part with an opening. |
| `offset_sheet` | Moves a face or a shell along its normals, trim and joins kept. A face with no parallel surface in its own family (a B-spline) is fitted, and its tolerance holds the measured deviation. |
| `make_thick_sheet` | Thickens a face or an open shell into a solid: the sheet and its offset (or half the thickness each way), closed by side faces along the free edges. |
| `offset_wire` | 2D offset of a planar wire. The caller picks the `Join` style (arcs or intersections). Useful for tool-path-like outlines. |
| `apply_draft` | Tilts faces by a draft angle about a neutral plane, for moulded parts. |

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

Each is a constrained boolean internally and records history like every
other operation.

`normal_projection` projects a wire onto a shape along the shape's normals
(for engraving). It returns the `Projected` curves on the target faces.

`split_face` (in `ogeom::heal`) cuts one face along curves that lie on it
(`Projection::OnFace`) or are dropped onto it along its normals
(`AlongNormals`) or along a direction (`Along`). Each curve must cross the
face from boundary to boundary or close inside it. The pieces share the new
edges, and the history maps the face to its pieces.
