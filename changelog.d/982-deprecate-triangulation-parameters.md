- **Deprecated: `Triangulation::parameters`** (#982).
  - A vertex at an apex holds only one of its `u` there. The field is
    still filled for this release and removed in a later one.
  - Migration: replace `mesh.parameters[t[k] as usize]` with
    `mesh.triangle_parameters(i).unwrap()[k]` (or index
    `corner_parameters[i][k]`), `mesh.parameters.len() ==
    mesh.positions.len()` with `mesh.has_parameters()`, and a write of the
    field with `set_vertex_parameters` once the triangles are final.
  - Breaking: `Triangulation` has two new public fields, so a struct
    literal needs `..Triangulation::default()`.
