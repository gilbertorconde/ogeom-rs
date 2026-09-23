# Making shapes

Everything here is in `ogeom::algo` unless stated otherwise. Every
operation returns a `Built` with the result and its history.

## Primitives

`make_box`, `make_cylinder`, `make_cone`, `make_sphere`, `make_torus`,
`make_wedge` and `make_half_space` each take a `Frame` (origin and
orientation) plus their dimensions. Degenerate dimensions are refused by
name, not clamped.

## Bottom-up construction

Build other shapes in B-rep order, from vertices up:

| Level | Functions |
|---|---|
| Vertex | `make_vertex` |
| Edge | `make_edge` (on a curve), `make_edge_between` (between points) |
| Wire | `make_wire` (ordered edges), `make_wire_unordered` (any order), `make_polygon` (straight sides) |
| Face | `make_face`, `make_face_on` (over a surface), `make_face_with_pcurves` (when the boundary's surface parametrisation matters) |
| Shell, solid, compound | `make_shell`, `make_solid`, `make_compound` |

Helpers:

- `sew` stitches faces that share boundaries within tolerance into shells.
- `is_wire_closed` and `is_shell_closed` check closure before you build the
  next level.

## Sweeps and fitting

- `make_prism` extrudes. `make_revolution` revolves.
- The general sweeps live in `ogeom::offset`, which holds their shared
  machinery: `make_pipe` (along a wire), `make_loft` (through sections),
  `make_evolved` (along a planar profile). See
  [Offsets, shells and features](offsets.md).
- `interpolate` fits a curve through points.
- `approximate` and `approximate_within` fit a curve near points to a
  stated tolerance. The fit reports the deviation it achieved.
- `make_text` renders text as wires, for engraving.

## History

Every constructor and every operation in later chapters records history:
which inputs `generated` which outputs, which were `modified`, and which
`is_deleted`. Stable references into a rebuilt model (for example "fillet
that edge") are resolved through this history. §7 of
[the data model](data-model.md) defines the contract.
