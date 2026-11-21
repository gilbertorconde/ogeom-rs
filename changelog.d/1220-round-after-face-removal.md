- **Fixed: a round whose far feature is taken off the solid no longer panics
  on an edge the removal rebuilt** (#1220).
  - `remove_faces` now traces each kept face to the face the result holds,
    through the sewing that closes the wound. Before, it traced it to a face
    from before the sewing.
  - The engine's rounds on a solid with a feature taken off find each asked
    edge again between its two faces' images when the removal rebuilt it.
    They no longer hand on the old edge, which panicked in `fillet_edges`,
    `fillet_edges_overflow` and `fillet_edges_corner` with "no entry found
    for key". The answer no longer depends on how the model numbers its
    nodes, so it is the same after a refused request as on the solid alone.
  - An edge asked that is not an edge of the solid is refused by name
    (`OgeomError::Dangling`) where the curved engine reads it.
