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

Where blends meet, the kernel handles the meeting: filleting the three
edges at a box corner produces the closed corner patch where they collide,
and the suite measures that solid against its closed-form volume.

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
