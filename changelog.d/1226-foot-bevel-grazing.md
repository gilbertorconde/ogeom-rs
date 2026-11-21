- **Fixed: a boss foot bevel whose leg touches an edge round's tangent line
  is refused by a name that says why** (#1226).
  - On a block with a rounded edge and a fused boss, the foot bevel was
    read as building up to 1.0, refused just above it and building again
    from 2.71. Only the limit reading jumped from 1.0001 to 2.71: every
    tenth from 1.1 to 2.7 builds, valid, with its fill within 7e-9 of a
    sampled reference.
  - At 1.0 the bevel's leg on the floor touches the round's tangent line.
    Just past it, the leg crosses that line at a grazing angle, and the
    floor and the round are left wedges too sharp for the booleans. Where
    neither blend into round nor trim through builds such a bevel (the
    leg meeting the line at under a hundredth of a radian),
    `chamfer_edges` and `chamfer_edges_overflow` now say so by name
    instead of giving the engine's generic refusal.
