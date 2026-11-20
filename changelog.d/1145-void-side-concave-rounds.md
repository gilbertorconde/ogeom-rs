- **Added: concave rounds the engine declined are rounded from the void's
  side** (#1145). Where every edge of a `fillet_edges` request is concave
  and no route builds it, the void beside the edges (a box about them less
  the solid) is rounded convex and the difference fused on as the fill.
  It is refused where a fill reaches the box or the result is not one sound
  solid. Three concave rounds meeting on a pocket's cylinder (a D pocket's
  floor corners with its uprights) now close on the ball's sphere under
  both `Setback` and `Rolling`, and a D pocket rounded all over rolls
  (ADR 0001, amendment "rolling corners on arcs, and concave corners from
  the void's side").
