- **Added: overflow options on rolling corners, variable radii and two distances (#1222).**
  - `fillet_edges_corner_overflow`, `fillet_edges_variable_overflow`,
    `fillet_edge_distances_overflow` and `blend_edges_overflow` take the
    `Overflow` set; `Overflow::default()` is each existing call.
  - A rolling corner passes its rounds to `fillet_edges_overflow`, so blend
    into round, hold edge and trim through switch off as for setback corners.
  - A law holding one radius, two equal elliptical distances and
    `blend_edges` rounds of one radius now blend a boss's foot into the edge
    round it reaches, where they were refused before.
  - A law that varies and two different distances take trim through alone;
    with trim through off, a result changing another face is refused by name.
