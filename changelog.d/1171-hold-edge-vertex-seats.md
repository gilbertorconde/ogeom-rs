- **Added: hold edge rests the ball on a vertex, and holds more rounds**
  (#1171).
  - Where a held edge stops and the next goes on at a convex corner of the
    face, the ball rests on the vertex between them: its centre stays the
    radius from it. The pieces on the edges and the vertex meet G1, and
    both edges are kept (ADR 0001, amendment "overflow handling", "Hold
    edge, more seats").
  - Several held stretches on one band, a chain of edges along one line
    over faces split in the same two planes, and end faces leaning up to
    60 degrees from square are now held. They went to trim through before.
  - A ball resting on two edges at once (a pocket's corner) is still not
    held and goes to trim through, as before. So does the boss underside
    of amendment "M5c", whose edge is concave.
  - A B-spline surface's own boundary iso-curve now gets the chart's
    straight side as its exact pcurve, also on a patch closing to a point.
