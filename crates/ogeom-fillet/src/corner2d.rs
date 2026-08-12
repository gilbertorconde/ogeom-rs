//! 2D fillets and chamfers: rounding and beveling a wire's corner.
//!
//! The sketch-plane cousins of the edge blends. A corner where two straight
//! edges of a wire meet is replaced by a tangent arc (the fillet) or a
//! straight cut at set distances (the chamfer); the two edges are trimmed
//! back on their own curves, and the wire is rebuilt with the connector in
//! the corner's place. The tangent construction for corners with curved
//! sides — a line meeting an arc, two arcs — is the 2D tangency problem
//! proper — docs/PARITY.md, fillet.corners-2d — rather than approximated here.

use crate::support::edge_curve;
use ogeom_algo::{Built, History, edge_vertices, make_edge_between, make_vertex, make_wire};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{CircleCurve, Curve};
use ogeom_math::{Circle, Direction, Frame, Point, Vector};
use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore};

/// Round a corner of a wire with an arc tangent to both of its edges.
///
/// `vertex` names the corner; the two edges meeting there must be straight.
/// The result is a new wire with the two edges trimmed to the tangency points
/// and the arc between them: the corner vertex is deleted, the edges are
/// modified into their trimmed selves, and the arc is generated from the
/// vertex.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the vertex is
/// not a corner of the wire between two straight edges, the edges are
/// collinear, `radius` is not a usable length, or the tangency points fall
/// off either edge.
pub fn fillet_corner_2d(
    model: &mut Model,
    wire: &Shape,
    vertex: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    let corner = corner_of(model, wire, vertex, tol)?;
    corner.opening(tol)?;
    // The tangent circle's centre lies on each side's offset locus — the
    // parallel line for a straight side, the concentric circle for an arc —
    // and where two loci cross is a candidate. GccAna's question, answered
    // the same way: enumerate the loci, intersect in closed form, and keep
    // the qualified candidate nearest the corner: centre on the corner's
    // inner side, both tangency feet on the edges themselves.
    let plane_z = {
        let n = corner.sides[0].away.cross(corner.sides[1].away);
        Direction::new(n, tol)?
    };
    let frame = Frame::new(
        corner.point,
        plane_z,
        Direction::new(corner.sides[0].away, tol)?,
        tol,
    )?;
    let flat = |p: Point| {
        let l = frame.to_local(p);
        ogeom_math::Point2::new(l.x, l.y)
    };
    let lift = |q: ogeom_math::Point2| {
        frame.origin() + frame.x().vector() * q.x + frame.y().vector() * q.y
    };
    let loci = |side: &Side| -> OgeomResult<Vec<ogeom_geom::PlanarCurve>> {
        Ok(match &side.curve {
            Curve::Line(l) => {
                let o = flat(l.axis().location);
                let d = flat(l.axis().location + l.axis().direction.vector()) - o;
                let d = d / d.magnitude();
                let n = ogeom_math::Vector2::new(-d.y, d.x);
                let axis = |shift: f64| {
                    ogeom_math::Direction2::new(d, tol)
                        .map(|dir| ogeom_math::Axis2::new(o + n * shift, dir))
                };
                vec![
                    ogeom_geom::Line2d::over(axis(radius)?, -1e6, 1e6)?.into(),
                    ogeom_geom::Line2d::over(axis(-radius)?, -1e6, 1e6)?.into(),
                ]
            }
            Curve::Circle(c) => {
                let centre = flat(c.circle().centre());
                let big = c.circle().radius();
                let mut out: Vec<ogeom_geom::PlanarCurve> = vec![
                    ogeom_geom::Circle2d::new(ogeom_math::Circle2::new(
                        ogeom_math::Frame2::new(centre, ogeom_math::Direction2::X),
                        big + radius,
                        tol,
                    )?)
                    .into(),
                ];
                if big - radius > tol.confusion() {
                    out.push(
                        ogeom_geom::Circle2d::new(ogeom_math::Circle2::new(
                            ogeom_math::Frame2::new(centre, ogeom_math::Direction2::X),
                            big - radius,
                            tol,
                        )?)
                        .into(),
                    );
                }
                out
            }
            _ => ogeom_bail!(Construction, "a corner side is neither line nor arc"),
        })
    };
    // Where does a candidate centre touch a side, and how far along the edge
    // is that from the corner in the side's own parameter?
    let foot = |side: &Side, c2: ogeom_math::Point2| -> Option<(ogeom_math::Point2, f64)> {
        match &side.curve {
            Curve::Line(l) => {
                let o = flat(l.axis().location);
                let d = flat(l.axis().location + l.axis().direction.vector()) - o;
                let d = d / d.magnitude();
                let t = o + d * ((c2 - o).dot(d));
                let corner2 = flat(side.curve.point_at(side.at, tol).ok()?);
                let away2 = {
                    let a = flat(lift(corner2) + side.away) - corner2;
                    a / a.magnitude()
                };
                Some((t, (t - corner2).dot(away2)))
            }
            Curve::Circle(circle) => {
                let centre = flat(circle.circle().centre());
                let big = circle.circle().radius();
                let v = c2 - centre;
                let m = v.magnitude();
                if m <= tol.confusion() {
                    return None;
                }
                let t = centre + v * (big / m);
                let corner2 = flat(side.curve.point_at(side.at, tol).ok()?);
                let a = (corner2 - centre).y.atan2((corner2 - centre).x);
                let b = (t - centre).y.atan2((t - centre).x);
                let mut delta = b - a;
                let tau = core::f64::consts::TAU;
                while delta > core::f64::consts::PI {
                    delta -= tau;
                }
                while delta < -core::f64::consts::PI {
                    delta += tau;
                }
                // Positive when the foot lies along the away direction: the
                // circle's own winding in the chart says which sign that is.
                let winding = {
                    let d1 = flat(lift(corner2) + side.away * 1.0) - corner2;
                    let radial = corner2 - centre;
                    (radial.x * d1.y - radial.y * d1.x).signum()
                };
                Some((t, delta * winding))
            }
            _ => None,
        }
    };
    let corner2 = flat(corner.point);
    let bis2 = {
        let a = flat(corner.point + corner.sides[0].away) - corner2;
        let b = flat(corner.point + corner.sides[1].away) - corner2;
        let u = a + b;
        u / u.magnitude()
    };
    let mut best: Option<(ogeom_math::Point2, [ogeom_math::Point2; 2], [f64; 2], f64)> = None;
    for la in loci(&corner.sides[0])? {
        for lb in loci(&corner.sides[1])? {
            let found = ogeom_intersect::intersect_curves_2d(
                &la,
                &lb,
                ogeom_intersect::CurveCurveOptions::default(),
                tol,
            )?;
            for crossing in &found.crossings {
                let c2 = crossing.point;
                if (c2 - corner2).dot(bis2) <= tol.confusion() {
                    continue;
                }
                let (Some((t0, d0)), Some((t1, d1))) =
                    (foot(&corner.sides[0], c2), foot(&corner.sides[1], c2))
                else {
                    continue;
                };
                // Both feet forward of the corner, within their edges.
                let room0 = corner.sides[0].room;
                let room1 = corner.sides[1].room;
                if d0 <= tol.parametric()
                    || d1 <= tol.parametric()
                    || d0 >= room0 - tol.parametric()
                    || d1 >= room1 - tol.parametric()
                {
                    continue;
                }
                let score = c2.distance(corner2);
                if best.as_ref().is_none_or(|(_, _, _, held)| score < *held) {
                    best = Some((c2, [t0, t1], [d0, d1], score));
                }
            }
        }
    }
    let Some((centre2, feet, deltas, _)) = best else {
        ogeom_bail!(
            Construction,
            "no circle of radius {radius} is tangent to both sides within \
             their own extents"
        );
    };
    let centre = lift(centre2);
    let contacts = [lift(feet[0]), lift(feet[1])];
    let trim_deltas = deltas;
    let arc = |model: &mut Model, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let w1 = contacts[0] - centre;
        let w2 = contacts[1] - centre;
        let z = Direction::new(w1.cross(w2), tol)?;
        let x = Direction::new(w1, tol)?;
        let frame = Frame::new(centre, z, x, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let sweep = w1.cross(w2).magnitude().atan2(w1.dot(w2));
        Ok(make_edge_between(
            model,
            Curve::Circle(CircleCurve::new(circle)),
            (0.0, sweep),
            from,
            to,
            tol,
        )?
        .shape)
    };
    rebuild(model, wire, vertex, &corner, trim_deltas, arc, tol)
}

/// Cut a corner of a wire, trimming `first` back along the earlier edge and
/// `second` along the later, joined by a straight segment.
///
/// "Earlier" and "later" follow the wire's own traversal order through the
/// corner.
///
/// # Errors
///
/// As [`fillet_corner_2d`].
pub fn chamfer_corner_2d(
    model: &mut Model,
    wire: &Shape,
    vertex: &Shape,
    first: f64,
    second: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for distance in [first, second] {
        if !distance.is_finite() || distance <= tol.confusion() {
            ogeom_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    let corner = corner_of(model, wire, vertex, tol)?;
    corner.opening(tol)?;
    // Distances are arc lengths along each side; on an arc that is an angle.
    let delta = |side: &Side, d: f64| -> f64 {
        match &side.curve {
            Curve::Circle(c) => d / c.circle().radius(),
            _ => d,
        }
    };
    let deltas = [
        delta(&corner.sides[0], first),
        delta(&corner.sides[1], second),
    ];
    let contacts = [
        corner.sides[0].curve.point_at(
            deltas[0].mul_add(corner.sides[0].sense, corner.sides[0].at),
            tol,
        )?,
        corner.sides[1].curve.point_at(
            deltas[1].mul_add(corner.sides[1].sense, corner.sides[1].at),
            tol,
        )?,
    ];
    let cut = |model: &mut Model, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let line = ogeom_geom::LineCurve::segment(contacts[0], contacts[1], tol)?;
        let curve = Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        Ok(make_edge_between(model, curve, domain, from, to, tol)?.shape)
    };
    rebuild(model, wire, vertex, &corner, deltas, cut, tol)
}

/// One side of a corner: the wire edge running into or out of it.
struct Side {
    /// Position in the wire's ordered edge list.
    index: usize,
    /// The occurrence as the wire uses it, orientation included.
    used: Shape,
    /// The side's own curve.
    curve: Curve,
    /// Unit tangent at the corner, pointing along the edge.
    away: Vector,
    /// The corner's parameter on the edge's own curve.
    at: f64,
    /// `+1` when walking away from the corner increases the parameter.
    sense: f64,
    /// How much parameter the edge has to give before its far end.
    room: f64,
    /// The far end's vertex.
    far: Shape,
}

/// A corner of a wire: the shared point and its two sides, in traversal
/// order — `sides[0]` runs into the corner, `sides[1]` out of it.
struct Corner {
    point: Point,
    edges: Vec<Shape>,
    sides: [Side; 2],
}

impl Corner {
    /// The opening angle between the two sides, strictly inside `(0, π)`.
    fn opening(&self, tol: Tolerances) -> OgeomResult<f64> {
        let (a, b) = (self.sides[0].away, self.sides[1].away);
        let angle = a.cross(b).magnitude().atan2(a.dot(b));
        if angle <= tol.angular() || angle >= core::f64::consts::PI - tol.angular() {
            ogeom_bail!(
                Construction,
                "the corner's edges are collinear; there is no corner to blend"
            );
        }
        Ok(angle)
    }
}

/// Find the corner `vertex` makes in `wire`: the two adjacent straight edges
/// and the geometry the trims run on.
fn corner_of(model: &Model, wire: &Shape, vertex: &Shape, tol: Tolerances) -> OgeomResult<Corner> {
    if model.kind_of(wire)? != ShapeType::Wire {
        ogeom_bail!(Construction, "a 2D blend rounds a corner of a wire");
    }
    let edges = explore(model, wire, Filter::OfType(ShapeType::Edge))?;
    if edges.len() < 2 {
        ogeom_bail!(Construction, "a corner needs at least two edges");
    }
    let n = edges.len();
    let mut found: Option<(usize, usize)> = None;
    for i in 0..n {
        let j = (i + 1) % n;
        let Some((_, end)) = edge_vertices(model, &edges[i])? else {
            continue;
        };
        let Some((start, _)) = edge_vertices(model, &edges[j])? else {
            continue;
        };
        if end.node() == vertex.node() && start.node() == vertex.node() {
            found = Some((i, j));
            break;
        }
    }
    let Some((i, j)) = found else {
        ogeom_bail!(
            Construction,
            "the vertex is not a corner between two consecutive edges of the \
             wire"
        );
    };

    let point = {
        let Some(node) = model.node(vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            ogeom_bail!(Construction, "vertex node holds no point");
        };
        vertex.transform(model.datums())?.apply(data.point)
    };

    let side = |model: &Model, index: usize, corner_at_end: bool| -> OgeomResult<Side> {
        let used = edges[index].clone();
        let (curve, range) = edge_curve(model, &used, tol)?;
        if !matches!(curve, Curve::Line(_) | Curve::Circle(_)) {
            ogeom_bail!(
                Construction,
                "a corner blend against a free-form side has no closed-form \
                 tangency; lines and arcs do"
            );
        }
        // The corner sits at the traversal end (or start), which for a
        // reversed occurrence is the stored range's other bound.
        let reversed = used.orientation() == Orientation::Reversed;
        let at_high = corner_at_end != reversed;
        let (at, sense, room) = if at_high {
            (range.1, -1.0, range.1 - range.0)
        } else {
            (range.0, 1.0, range.1 - range.0)
        };
        let away = {
            let d = curve.d1_at(at, tol)? * sense;
            let m = d.magnitude();
            if m <= tol.confusion() {
                ogeom_bail!(Construction, "a corner edge is degenerate at its corner");
            }
            d / m
        };
        let Some((start, end)) = edge_vertices(model, &used)? else {
            ogeom_bail!(Construction, "a corner edge has no bounding vertices");
        };
        let far = if corner_at_end { start } else { end };
        Ok(Side {
            index,
            used,
            curve,
            away,
            at,
            sense,
            room,
            far,
        })
    };
    let sides = [side(model, i, true)?, side(model, j, false)?];
    Ok(Corner {
        point,
        edges,
        sides,
    })
}

/// Trim both sides, build the connector between the new vertices, and
/// reassemble the wire with history.
fn rebuild(
    model: &mut Model,
    wire: &Shape,
    vertex: &Shape,
    corner: &Corner,
    trims: [f64; 2],
    connector: impl FnOnce(&mut Model, &Shape, &Shape) -> OgeomResult<Shape>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut history = History::new();
    let mut trimmed: Vec<Shape> = Vec::with_capacity(2);
    let mut joints: Vec<Shape> = Vec::with_capacity(2);
    for (side, trim) in corner.sides.iter().zip(trims) {
        if trim >= side.room - tol.parametric() {
            ogeom_bail!(
                Construction,
                "a trim of {trim} consumes the whole edge; the blend reaches \
                 past the corner's neighbours"
            );
        }
        let (curve, range) = edge_curve(model, &side.used, tol)?;
        let contact = trim.mul_add(side.sense, side.at);
        let new_range = if side.sense > 0.0 {
            (contact, range.1)
        } else {
            (range.0, contact)
        };
        let joint = make_vertex(model, curve.point_at(contact, tol)?).shape;
        // The trimmed edge runs between the far vertex and the new tangency
        // vertex, on the same curve, with the same orientation flag the wire
        // used before.
        let (from, to) = if side.sense > 0.0 {
            (joint.clone(), side.far.clone())
        } else {
            (side.far.clone(), joint.clone())
        };
        let mut new_edge = make_edge_between(model, curve, new_range, &from, &to, tol)?.shape;
        if side.used.orientation() == Orientation::Reversed {
            new_edge = new_edge.reversed();
        }
        history.modify(&side.used, new_edge.clone());
        trimmed.push(new_edge);
        joints.push(joint);
    }

    let joined = connector(model, &joints[0], &joints[1])?;
    history.generate(vertex, joined.clone());
    history.delete(vertex);

    let mut edges: Vec<Shape> = Vec::with_capacity(corner.edges.len() + 1);
    for (k, e) in corner.edges.iter().enumerate() {
        if k == corner.sides[0].index {
            edges.push(trimmed[0].clone());
            edges.push(joined.clone());
        } else if k == corner.sides[1].index {
            edges.push(trimmed[1].clone());
        } else {
            edges.push(e.clone());
        }
    }
    // The corner pair may wrap the list's end; rotate so the two halves stay
    // adjacent in traversal order.
    if corner.sides[1].index < corner.sides[0].index {
        edges.rotate_left(corner.sides[1].index + 1);
    }
    let built = make_wire(model, &edges, tol)?;
    history.modify(wire, built.shape.clone());
    Ok(Built::new(built.shape, history))
}
