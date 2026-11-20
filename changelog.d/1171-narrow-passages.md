- **Added: narrow passages for fillets** (#1171).
  - `Overflow` gains `narrow_blend_made`, `narrow_whole_face` and
    `narrow_partial_face`, on by default. A narrow passage is a concave
    zone tighter than the radius: two walls of a floor no farther apart
    than twice the radius, or a corner face tighter than the radius. Each round there
    keeps its constant section and runs on through the zone, and the
    rounds meet where they cross (ADR 0001, amendment "overflow handling",
    "Narrow passages").
  - Built for a planar floor's rim with planar walls square to it: a
    through slot or a stadium pocket no wider than twice the radius, and a
    pocket whose corner faces are tighter than the radius (tangent, or at
    kinks with whole face). These were refused before.
  - Each strategy takes what the cheaper ones take. With all three off a
    passage is refused by name. Partial face is not needed by any built
    case yet.
