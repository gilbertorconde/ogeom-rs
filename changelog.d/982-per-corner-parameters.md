- **Added: per-corner surface parameters and normals on `Triangulation`**
  (#982).
  - `Triangulation::corner_parameters` holds each triangle's three `(u, v)`
    in the order of its indices, and `corner_normals` its three normals. A
    point with more than one `u` (an apex, a pole, later a seam) keeps the
    right one in each triangle, so the mean of a triangle's corner
    parameters is a point on its own patch. At an apex the corner holds the
    apex's `v` and the middle of the other two corners' `u`.
  - `triangle_parameters(i)` and `triangle_normals(i)` read a triangle's
    corners, falling back to the per-vertex arrays for a mesh built with
    those alone. `has_parameters`, `vertex_parameters` and
    `set_vertex_parameters` cover the per-vertex cases.
  - The native format writes a triangle whose corners differ from its
    vertices as an `fc` line, in version 4. A file with no such triangle is
    written at the version it was before.
