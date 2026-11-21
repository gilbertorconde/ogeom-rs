- **Fixed: a face that runs to an apex or a pole meshed a triangle of no
  area** (#982).
  - A cone's apex or a sphere's pole is now one vertex of the face's mesh,
    with a fan of real triangles round it. Before, the chart ring ran along
    the collapsed row and left two copies of the apex a whole turn apart,
    joined by a triangle with two corners at the apex.
  - On the step-conformance parts `cone_point`, `cone_point_vertex` and
    `dome`, at chords 0.01/3, 0.1 and 1.0, no triangle has two corners at
    one point, each apex is one vertex per face, the chord is met and the
    welded bodies close.
