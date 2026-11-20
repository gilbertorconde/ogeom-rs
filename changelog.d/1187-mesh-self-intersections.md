- **Added: self-intersections of a mesh, found and resolved** (#1187).
  - `mesh::self_intersections` lists every pair of triangles that meet
    anywhere but at the corners they share, decided exactly by the mesh
    boolean's predicates, lower number first, sorted. Any triangulation is
    taken, open or not; a hierarchy of boxes keeps a large scan fast.
  - `mesh::resolve_self_intersections` cuts a closed, 2-manifold mesh
    along every crossing and keeps the boundary of the points it winds
    about at least once, by exact winding numbers: overlapping shells
    merge, a nested shell is dropped, an inside-out one is kept as a
    cavity, coplanar faces facing the same way are kept once. `Resolved`
    carries the mesh, the count of crossing pairs, a `cut` flag per
    triangle and the rounding bound, as `MeshBoolean` does. A mesh that
    does not cross itself comes back with the same triangles.
  - Where rounding the new vertices tilts a sliver past its neighbour,
    the result is resolved again, at most three more times; a flat
    triangle rounding leaves is removed by flipping it into its
    neighbour, which moves no vertex.
