//! 2D fillets and chamfers: rounding and beveling a wire's corner.
//!
//! The sketch-plane cousins of the edge blends. A corner where two straight
//! edges of a wire meet is replaced by a tangent arc (the fillet) or a
//! straight cut at set distances (the chamfer); the two edges are trimmed
//! back on their own curves, and the wire is rebuilt with the connector in
//! the corner's place. The tangent construction for corners with curved
//! sides — a line meeting an arc, two arcs — is the 2D tangency problem
//! proper, recorded in the deferred table rather than approximated here.

use crate::support::edge_curve;
use og_algo::{Built, History, edge_vertices, make_edge_between, make_vertex, make_wire};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve3d as _;
use og_geom::{CircleCurve, Curve};
use og_math::{Circle, Direction, Frame, Point, Vector};
use og_topo::{Filter, Model, Orientation, Shape, ShapeType, explore};

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
/// [`OgError::Construction`](og_core::OgError::Construction) if the vertex is
/// not a corner of the wire between two straight edges, the edges are
/// collinear, `radius` is not a usable length, or the tangency points fall
/// off either edge.
pub fn fillet_corner_2d(
    model: &mut Model,
    wire: &Shape,
    vertex: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        og_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    let corner = corner_of(model, wire, vertex, tol)?;
    let opening = corner.opening(tol)?;
    // The tangency points sit where a circle of this radius touches both
    // sides, and the centre sits on the bisector at the matching height.
    let trim = radius / (opening / 2.0).tan();
    let centre = {
        let bisector = {
            let u = corner.sides[0].away + corner.sides[1].away;
            u / u.magnitude()
        };
        corner.point + bisector * (radius / (opening / 2.0).sin())
    };
    let contacts = [
        corner.point + corner.sides[0].away * trim,
        corner.point + corner.sides[1].away * trim,
    ];
    let arc = |model: &mut Model, from: &Shape, to: &Shape| -> OgResult<Shape> {
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
    rebuild(model, wire, vertex, &corner, [trim, trim], arc, tol)
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
) -> OgResult<Built> {
    for distance in [first, second] {
        if !distance.is_finite() || distance <= tol.confusion() {
            og_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    let corner = corner_of(model, wire, vertex, tol)?;
    corner.opening(tol)?;
    let contacts = [
        corner.point + corner.sides[0].away * first,
        corner.point + corner.sides[1].away * second,
    ];
    let cut = |model: &mut Model, from: &Shape, to: &Shape| -> OgResult<Shape> {
        let line = og_geom::LineCurve::segment(contacts[0], contacts[1], tol)?;
        let curve = Curve::Line(line);
        let domain = og_geom::Curve3d::domain(&curve);
        Ok(make_edge_between(model, curve, domain, from, to, tol)?.shape)
    };
    rebuild(model, wire, vertex, &corner, [first, second], cut, tol)
}

/// One side of a corner: the wire edge running into or out of it.
struct Side {
    /// Position in the wire's ordered edge list.
    index: usize,
    /// The occurrence as the wire uses it, orientation included.
    used: Shape,
    /// Unit direction from the corner along this edge.
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
    fn opening(&self, tol: Tolerances) -> OgResult<f64> {
        let (a, b) = (self.sides[0].away, self.sides[1].away);
        let angle = a.cross(b).magnitude().atan2(a.dot(b));
        if angle <= tol.angular() || angle >= core::f64::consts::PI - tol.angular() {
            og_bail!(
                Construction,
                "the corner's edges are collinear; there is no corner to blend"
            );
        }
        Ok(angle)
    }
}

/// Find the corner `vertex` makes in `wire`: the two adjacent straight edges
/// and the geometry the trims run on.
fn corner_of(model: &Model, wire: &Shape, vertex: &Shape, tol: Tolerances) -> OgResult<Corner> {
    if model.kind_of(wire)? != ShapeType::Wire {
        og_bail!(Construction, "a 2D blend rounds a corner of a wire");
    }
    let edges = explore(model, wire, Filter::OfType(ShapeType::Edge))?;
    if edges.len() < 2 {
        og_bail!(Construction, "a corner needs at least two edges");
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
        og_bail!(
            Construction,
            "the vertex is not a corner between two consecutive edges of the \
             wire"
        );
    };

    let point = {
        let Some(node) = model.node(vertex) else {
            og_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            og_bail!(Construction, "vertex node holds no point");
        };
        vertex.transform(model.datums())?.apply(data.point)
    };

    let side = |model: &Model, index: usize, corner_at_end: bool| -> OgResult<Side> {
        let used = edges[index].clone();
        let (curve, range) = edge_curve(model, &used)?;
        let Curve::Line(_) = &curve else {
            og_bail!(
                Construction,
                "a corner blend with a curved side is the 2D tangency \
                 problem, recorded in the deferred table"
            );
        };
        // The corner sits at the traversal end (or start), which for a
        // reversed occurrence is the stored range's other bound.
        let reversed = used.orientation() == Orientation::Reversed;
        let at_high = corner_at_end != reversed;
        let (at, sense, room) = if at_high {
            (range.1, -1.0, range.1 - range.0)
        } else {
            (range.0, 1.0, range.1 - range.0)
        };
        let far_param = if at_high { range.0 } else { range.1 };
        let away = {
            let far_point = curve.point_at(far_param, tol)?;
            (far_point - point) / point.distance(far_point)
        };
        let Some((start, end)) = edge_vertices(model, &used)? else {
            og_bail!(Construction, "a corner edge has no bounding vertices");
        };
        let far = if corner_at_end { start } else { end };
        Ok(Side {
            index,
            used,
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
    connector: impl FnOnce(&mut Model, &Shape, &Shape) -> OgResult<Shape>,
    tol: Tolerances,
) -> OgResult<Built> {
    let mut history = History::new();
    let mut trimmed: Vec<Shape> = Vec::with_capacity(2);
    let mut joints: Vec<Shape> = Vec::with_capacity(2);
    for (side, trim) in corner.sides.iter().zip(trims) {
        if trim >= side.room - tol.parametric() {
            og_bail!(
                Construction,
                "a trim of {trim} consumes the whole edge; the blend reaches \
                 past the corner's neighbours"
            );
        }
        let (curve, range) = edge_curve(model, &side.used)?;
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
