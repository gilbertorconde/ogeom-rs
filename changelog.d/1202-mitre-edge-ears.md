- **Fixed: two curved faces could mesh the same triangle across the edge
  they share** (#1202).
  - Where an edge's polyline bulges out of a curved face's straight chart
    chords, the face could close off a point inside the edge with its two
    neighbours along it. That triangle's corners are the edge's alone. The
    face across, its border bulging the same way, could draw the very same
    triangle, and welded, the chord under the two was four triangles' side.
    A curved face now turns that chord to the triangle across it, so the
    point reaches into the face. A plane keeps such a triangle, which lies
    in it.
  - Two quarter cylinders mitred on a quarter ellipse (a block's two
    rounded edges meeting in a corner) drew such a triangle at chord 0.1;
    now neither face does at any chord from 0.0125 to 0.5, and the two
    weld with no side of more than two triangles. The step-conformance
    parts mesh exactly as before.
