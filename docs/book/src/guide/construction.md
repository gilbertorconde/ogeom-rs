# Making shapes

Everything here lives in `ogeom::algo` unless said otherwise, and every
operation returns a `Built` carrying the result and its history.

## Primitives

`make_box`, `make_cylinder`, `make_cone`, `make_sphere`, `make_torus`,
`make_wedge` and `make_half_space` each take a `Frame` — an origin and an
orientation — and their dimensions. Degenerate dimensions are refused by
name, not clamped.

## Bottom-up

When a shape is not a primitive, it is built the way the B-rep is
structured: `make_vertex`, `make_edge` (on a curve) and `make_edge_between`
(between points), `make_wire` from ordered edges — or
`make_wire_unordered` when the order is not known — `make_polygon` for the
straight-sided case, `make_face` and `make_face_on` over a surface,
`make_face_with_pcurves` when the boundary's surface parametrisation
matters, then `make_shell`, `make_solid`, `make_compound`. `sew` stitches
faces that share boundaries within tolerance into shells; `is_wire_closed`
and `is_shell_closed` answer the question their names ask before you
promote anything.

## Sweeps and fitting

`make_prism` extrudes and `make_revolution` revolves. The general sweeps —
`make_pipe` along a wire, `make_loft` through sections, `make_evolved`
along a planar profile — live in `ogeom::offset`, which owns the sweep
machinery they share.

Curves through data come from `interpolate` (through points) and
`approximate` / `approximate_within` (near points, to a stated tolerance —
the fit reports what it achieved). `make_text` renders text as wires for
engraving.

## History

Every constructor and every operation in the chapters that follow emits
history: which inputs `generated` which outputs, which were `modified`,
which `is_deleted`. This is not optional bookkeeping — stable references
into a rebuilt model ("fillet *that* edge") are resolved by walking it,
and §7 of [the data model](data-model.md) is the contract.
