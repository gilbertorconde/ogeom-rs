- **Added: remeshing and subdivision of an `Editable`** (#1186).
  - `Editable::remesh` takes a `Remesh`: a target edge length, passes,
    an optional crease angle, an optional region and an optional face
    limit. Each pass splits edges longer than 4/3 of the target (longest
    first), collapses those shorter than 4/5 where the mesh stays a
    2-manifold and no face turns over, flips toward valence 6 (4 on a
    rim), and relaxes each vertex across its tangent plane before
    projecting it back onto the input. Crease and rim vertices slide
    along their edges; corners stay. Faces with a corner outside the
    region stay as they were. `RemeshReport` counts the splits, collapses
    and flips and gives the largest and RMS distance from the new surface
    to the old.
  - `Editable::subdivide` cuts every face into four, `levels` times:
    `Subdivision::Flat` at edge midpoints, `Subdivision::Smooth` by Loop's
    scheme with creases and rims kept sharp. It returns a `SubdivideReport`.
  - Both carry source vertices as `simplify` does, and refuse a result past
    `max_faces` by name, leaving the mesh as it was.
