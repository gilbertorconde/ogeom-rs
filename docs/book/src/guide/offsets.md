# Offsets, shells and features

Everything here is `ogeom::offset`.

## Offsetting

- **`offset_shape`** moves a solid's boundary outward (or inward, with a
  negative distance) along its normals, rebuilding the intersections
  where offset surfaces collide.
- **`make_thick_solid`** is shelling: remove the named faces, offset the
  rest inward, and join — the hollowed casting with an opening.
- **`offset_wire`** is the 2D counterpart on a planar wire, with the
  `Join` style (arcs or intersections) chosen by the caller — the basis
  of tool-path-like outlines.
- **`apply_draft`** tilts faces by a draft angle about a neutral plane,
  the moulding operation.

## Sweeps that live here

The general sweeps share this crate's machinery: `make_pipe` (and
`make_pipe_skinned`) along a spine wire, `make_loft` (and
`make_loft_skinned`) through profile sections, `make_evolved` sweeping a
profile along a planar spine. `make_filling` builds the face that fills a
boundary — the N-sided patch.

## Features

The feature operations combine a sketch with a solid in one step:
`feature_prism` (bosses and pockets), `feature_revol` (revolved bosses
and grooves), `feature_rib` and `feature_slot`. Each is a constrained
boolean under the hood and emits the same history every operation does.

`normal_projection` projects a wire onto a shape along the shape's own
normals — the engraving projection — returning the `Projected` curves on
the target faces.
