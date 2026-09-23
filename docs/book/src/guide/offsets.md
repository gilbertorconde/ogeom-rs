# Offsets, shells and features

Everything here is in `ogeom::offset`.

## Offsetting

| Function | What it does |
|---|---|
| `offset_shape` | Moves a solid's boundary along its normals (outward, or inward with a negative distance) and rebuilds the intersections where offset surfaces collide. |
| `make_thick_solid` | Shelling: removes the named faces, offsets the rest inward and joins them. Gives a hollow part with an opening. |
| `offset_wire` | 2D offset of a planar wire. The caller picks the `Join` style (arcs or intersections). Useful for tool-path-like outlines. |
| `apply_draft` | Tilts faces by a draft angle about a neutral plane, for moulded parts. |

## Sweeps

The general sweeps share this module's machinery:

- `make_pipe` and `make_pipe_skinned`: sweep along a single spine edge.
- `make_loft` and `make_loft_skinned`: loft through profile sections.
- `make_evolved`: sweep a profile along a planar spine.
- `make_filling`: build an N-sided patch face that fills a boundary.

## Features

Feature operations combine a sketch with a solid in one step:

- `feature_prism`: bosses and pockets.
- `feature_revol`: revolved bosses and grooves.
- `feature_rib` and `feature_slot`.

Each is a constrained boolean internally and records history like every
other operation.

`normal_projection` projects a wire onto a shape along the shape's normals
(for engraving). It returns the `Projected` curves on the target faces.
