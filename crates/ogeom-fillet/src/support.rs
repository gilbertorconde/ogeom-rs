//! Shared scaffolding for the subtractive blends.
//!
//! A chamfer and a constant-radius fillet on a straight edge between planar
//! faces differ only in the face that replaces the edge: a bevel plane for
//! one, a tangent cylinder for the other. Everything around that face (the
//! seat on the solid, the legs running along the adjacent faces, faces
//! assembled from explicit curves with exact pcurves) is one piece of
//! scaffolding, kept here so the two operations cannot drift apart.

use ogeom_algo::{make_edge, make_edge_between, make_revolution_band, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{CircleCurve, Curve, CylinderSurface, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, explore};

/// Where a blend sits on a solid: a straight edge and the two planar faces
/// meeting there, reduced to the numbers the wedge construction runs on.
pub(crate) struct Seat {
    /// The edge's start, at its lower parameter.
    pub start: Point,
    /// The edge's end.
    pub end: Point,
    /// Unit direction from start to end.
    pub along: Vector,
    /// Outward unit normals of the two faces, in discovery order.
    pub normals: [Vector; 2],
    /// The two faces themselves, in the same order, so a caller naming a
    /// face can find which leg it owns.
    pub faces: [Shape; 2],
    /// Whether the edge is convex: material inside the dihedral, so a blend
    /// subtracts. A concave edge's blend adds, with every sign mirrored.
    pub convex: bool,
}

impl Seat {
    /// On face `i`, the unit direction perpendicular to the edge that walks
    /// away from the *other* face, into the material the blend cuts back
    /// along.
    pub fn leg(&self, i: usize, tol: Tolerances) -> OgeomResult<Vector> {
        let own = self.normals[i];
        let other = self.normals[1 - i];
        let mut t = own.cross(self.along);
        if t.dot(other) > 0.0 {
            t = -t;
        }
        let m = t.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "a face is tangent to its own edge");
        }
        Ok(t / m)
    }
}

/// Whether `candidate` is the same edge occurrence as `edge`: the same node,
/// placed the same way in the world. Two references through different
/// location paths that compose to one placement are the same occurrence; the
/// same node re-instanced elsewhere (a prism's profile edge at both caps)
/// is not.
pub(crate) fn same_occurrence(
    model: &Model,
    candidate: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> bool {
    if candidate.node() != edge.node() {
        return false;
    }
    if candidate.is_same(edge) {
        return true;
    }
    let (Ok(a), Ok(b)) = (
        candidate.transform(model.datums()),
        edge.transform(model.datums()),
    ) else {
        return false;
    };
    use ogeom_math::Point;
    [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, 0.0, 1.0),
    ]
    .into_iter()
    .all(|p| a.apply(p).distance(b.apply(p)) <= tol.confusion())
}

/// An edge's 3D curve and range, placed in the world.
pub(crate) fn edge_curve(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64))> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "expected an edge");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "the edge has no curve to blend along");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    // Baked into the world: the edge lives wherever its placement puts it,
    // and every seat is measured there. A rigid placement preserves the
    // parameterization, so the stored range carries over unchanged.
    use ogeom_geom::Transformable as _;
    let placed = geometry
        .clone()
        .transformed(&edge.transform(model.datums())?, tol)?;
    Ok((placed, *range))
}

/// A straight edge's line, placed in the world, as a point on it and its
/// unit direction; `None` for an edge with no curve or a curved one.
fn edge_line(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<(Point, Vector)>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    if data.curve3d().is_none() {
        return Ok(None);
    }
    let (Curve::Line(line), _) = edge_curve(model, edge, tol)? else {
        return Ok(None);
    };
    let axis = line.axis();
    Ok(Some((axis.location, axis.direction.vector())))
}

/// Whether every face of `solid` bounding `edge` is planar: the straight
/// seat's own forms speak only planes, and a straight edge on a curved
/// face (a ruling of a drum where it meets a wall) is the march's.
pub(crate) fn hosts_planar(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let planar = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
        if !planar {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `end`, a straight edge's end, slid back along `inward` to where its
/// vertex `(at, reach)` stands: the vertex's foot on the line, or where the
/// line meets one of `meeting` (the lines, as `(point, unit direction)`, of
/// the other straight edges of its faces leaving the vertex), whichever
/// lies furthest back. A meeting counts where it stands within `reach` of
/// the vertex and of `end`; the end moves only where the vertex stands
/// within `reach` of it and the line runs past the place found. `end`
/// itself otherwise.
fn slid_back(
    end: Point,
    inward: Vector,
    (at, reach): (Point, f64),
    meeting: &[(Point, Vector)],
    tol: Tolerances,
) -> Point {
    if end.distance(at) > reach {
        return end;
    }
    let mut overrun = (at - end).dot(inward);
    for &(through, along) in meeting {
        // Where the line comes nearest the other: on the face both lie
        // in, where they cross.
        let normal = inward.cross(along);
        let across = normal.dot(normal);
        if across <= tol.angular() * tol.angular() {
            continue;
        }
        let back = (through - end).cross(along).dot(normal) / across;
        if back.abs() <= reach && (end + inward * back).distance(at) <= reach {
            overrun = overrun.max(back);
        }
    }
    if overrun > 0.0 {
        end + inward * overrun
    } else {
        end
    }
}

/// Find the seat of a blend: the straight edge's ends and direction, and the
/// outward normals of the exactly two planar faces of `solid` meeting there.
///
/// Refuses tangent edges, which have no corner. Convexity is read from the
/// solid and carried in the seat: a convex edge's wedge subtracts, a
/// concave edge's adds, with every sign mirrored.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// not straight, is not shared by exactly two planar faces of `solid`, or the
/// edge is tangent.
pub(crate) fn planar_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Seat> {
    let (curve, range) = edge_curve(model, edge, tol)?;
    let Curve::Line(_) = &curve else {
        ogeom_bail!(
            Construction,
            "blending a curved edge needs the marching blend machinery; this \
             is the straight-edge form"
        );
    };
    let mut start = curve.point_at(range.0, tol)?;
    let mut end = curve.point_at(range.1, tol)?;
    let along = (end - start) / start.distance(end);
    let vertices = ogeom_algo::edge_vertices(model, edge)?;
    let mut normals: Vec<Vector> = Vec::new();
    let mut faces: Vec<Shape> = Vec::new();
    // The lines of the other straight edges of the two faces leaving each
    // of the edge's vertices: where each face's boundary turns there.
    let mut meeting: [Vec<(Point, Vector)>; 2] = [Vec::new(), Vec::new()];
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let bounds = explore(model, &face, Filter::OfType(ShapeType::Edge))?;
        if !bounds.iter().any(|e| same_occurrence(model, e, edge, tol)) {
            continue;
        }
        if let Some((a, b)) = &vertices {
            for other in &bounds {
                if same_occurrence(model, other, edge, tol) {
                    continue;
                }
                let Some((p, q)) = ogeom_algo::edge_vertices(model, other)? else {
                    continue;
                };
                let leaves = [a, b].map(|v| {
                    same_occurrence(model, &p, v, tol) || same_occurrence(model, &q, v, tol)
                });
                if leaves == [false, false] {
                    continue;
                }
                let Some(line) = edge_line(model, other, tol)? else {
                    continue;
                };
                for (lines, leaves) in meeting.iter_mut().zip(leaves) {
                    if leaves {
                        lines.push(line);
                    }
                }
            }
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
            ogeom_bail!(
                Construction,
                "blending an edge of a curved face needs the marching blend \
                 machinery; this is the planar form"
            );
        };
        let placement = face.transform(model.datums())?;
        let mut normal = placement.apply_vector(plane.plane().normal().vector());
        if face.orientation() == ogeom_topo::Orientation::Reversed {
            normal = -normal;
        }
        normals.push(normal);
        faces.push(face);
    }
    if normals.len() != 2 {
        ogeom_bail!(
            Construction,
            "a blend needs an edge shared by exactly two faces, found {}",
            normals.len()
        );
    }
    // The edge ends where its vertices stand, which its line may run past
    // by up to their tolerance, and a vertex stands only within its
    // tolerance of the corner where each of the edge's faces turns: where
    // the line meets the next edge of that face's boundary. A blend's or a
    // chamfer's cap stands square to the edge at its end, and a cap past
    // the vertex's foot or past such a corner cuts into a neighbour's that
    // ends there, where the two meet flush. So an end is slid back along
    // the line by exactly the line's overrun of the furthest back of those
    // places, each within the vertex's tolerance; it leaves at most that
    // tolerance of material, rather than cutting into a neighbour. Each end
    // is paired with its own vertex, and where the two vertices'
    // tolerances span the whole edge, which vertex an end belongs to is not
    // known and neither end moves; no end moves by more than its vertex's
    // tolerance, so the two slides never meet.
    if let Some((a, b)) = &vertices {
        let stand = |v: &Shape| -> OgeomResult<(Point, f64)> {
            let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Dangling, "an edge's vertex is not in this model");
            };
            let at = v.transform(model.datums())?.apply(data.point);
            Ok((at, data.tolerance.get().max(tol.confusion())))
        };
        let (first, second) = (stand(a)?, stand(b)?);
        if !a.is_same(b) && start.distance(end) > first.1 + second.1 {
            let [at_first, at_second] = &meeting;
            let ((start_at, start_meets), (end_at, end_meets)) = if start.distance(first.0)
                + end.distance(second.0)
                <= start.distance(second.0) + end.distance(first.0)
            {
                ((first, at_first), (second, at_second))
            } else {
                ((second, at_second), (first, at_first))
            };
            start = slid_back(start, along, start_at, start_meets, tol);
            end = slid_back(end, -along, end_at, end_meets, tol);
        }
    }

    let mut seat = Seat {
        start,
        end,
        along,
        normals: [normals[0], normals[1]],
        faces: [faces[0].clone(), faces[1].clone()],
        convex: true,
    };
    // Convexity is read from the face itself, not derived from the normals:
    // the leg construction cannot answer it, because it *chooses* its side.
    // Which way the first face actually extends from the edge (sampled
    // against its own trim) leans behind the other face's plane on a convex
    // edge and in front of it on a concave one.
    let raw = {
        let t = normals[0].cross(along);
        let m = t.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "a face is tangent to its own edge");
        }
        t / m
    };
    let mid = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    let mut face_side: Option<Vector> = None;
    'scales: for scale in [1e-3, 1e-2, 5e-2] {
        let eps = start.distance(end) * scale;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        for dir in [raw, -raw] {
            if crate::support::on_face_side(
                model,
                &seat.faces[0],
                mid + dir * eps,
                deflection,
                tol,
            )? {
                face_side = Some(dir);
                break 'scales;
            }
        }
    }
    let Some(extends) = face_side else {
        ogeom_bail!(
            Construction,
            "cannot read which way the edge's face extends; the face is \
             thinner than the probe can resolve"
        );
    };
    let lean = extends.dot(seat.normals[1]);
    if lean.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the edge's faces are tangent; there is no corner to blend"
        );
    }
    seat.convex = lean < 0.0;
    Ok(seat)
}

/// Where a revolved blend sits: a circular rim shared by one perpendicular
/// planar cap and one coaxial cylindrical wall, reduced to the numbers the
/// revolved wedge runs on.
///
/// Four seats, one parameterization. With `sigma` the wall's outward radial
/// sign and `tau` telling whether the wall extends away from the cap's
/// outward side, `tau` alone decides whether the wedge subtracts (the
/// external rim and the hole's rim, both convex) or fuses (the boss base and
/// the blind hole's floor, both concave).
pub(crate) struct RevolvedSeat {
    /// The rim's centre.
    pub centre: Point,
    /// The rim's radius, also the wall's.
    pub radius: f64,
    /// The cap's outward unit normal.
    pub up: Vector,
    /// The rim frame's `x`, so every ring built on the seat shares a
    /// parameter origin.
    pub x_ref: Direction,
    /// The wall's outward radial sign: `+1` material inside, `-1` a bore.
    pub sigma: f64,
    /// `+1` when the wall extends away from the cap's outward side (the rim
    /// configurations) and `-1` alongside it, the concave seats.
    pub tau: f64,
    /// The planar cap.
    pub cap_face: Shape,
    /// The cylindrical wall.
    pub wall_face: Shape,
    /// Whether the rims on the rim's circle (edges bounding walls on its
    /// cylinder) run its whole turn, in one piece or in several.
    pub closes: bool,
}

impl RevolvedSeat {
    /// Whether the wedge fuses rather than subtracts.
    pub const fn additive(&self) -> bool {
        self.tau < 0.0
    }

    /// A frame at `origin` sharing the seat's axis and parameter origin.
    pub fn frame_at(&self, origin: Point, tol: Tolerances) -> OgeomResult<Frame> {
        Frame::new(origin, Direction::new(self.up, tol)?, self.x_ref, tol)
    }

    /// A full ring on the seat's axis: radius `r` in the plane through
    /// `origin`.
    pub fn ring(
        &self,
        model: &mut Model,
        origin: Point,
        r: f64,
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let circle = Circle::new(self.frame_at(origin, tol)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        let domain = ogeom_geom::Curve3d::domain(&curve);
        Ok(make_edge(model, curve, domain, tol)?.shape)
    }
}

/// Find a revolved blend's seat: the exactly one perpendicular planar cap and
/// one coaxial cylindrical wall meeting at a circular rim.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// rim's faces are not that pair, or the wall has no far ring to read its
/// side from.
pub(crate) fn revolved_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    tol: Tolerances,
) -> OgeomResult<RevolvedSeat> {
    let rim_circle = rim.circle();
    let rim_centre = rim_circle.centre();
    let rim_radius = rim_circle.radius();
    let rim_axis = rim_circle.frame().z().vector();

    // The two faces at the rim: exactly one plane and one coaxial cylinder.
    let mut cap: Option<(Shape, Vector)> = None;
    let mut wall: Option<(Shape, f64)> = None;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let placement = face.transform(model.datums())?;
        let reversed = face.orientation() == Orientation::Reversed;
        match model.geometry().surface(data.surface) {
            Some(SurfaceGeometry::Plane(p)) => {
                let mut normal = placement.apply_vector(p.plane().normal().vector());
                if reversed {
                    normal = -normal;
                }
                if normal.cross(rim_axis).magnitude() > tol.angular()
                    || p.plane()
                        .distance_to(placement.inverse()?.apply(rim_centre))
                        > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's cap is not the perpendicular plane through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                cap = Some((face.clone(), normal));
            }
            Some(SurfaceGeometry::Cylinder(c)) => {
                let cyl = c.cylinder();
                let axis_point = placement.apply(cyl.frame().origin());
                let axis_z = placement.apply_vector(cyl.frame().z().vector());
                let off_axis = {
                    let to_rim = rim_centre - axis_point;
                    (to_rim - axis_z * to_rim.dot(axis_z)).magnitude()
                };
                if axis_z.cross(rim_axis).magnitude() > tol.angular()
                    || off_axis > tol.confusion()
                    || (cyl.radius() - rim_radius).abs() > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's wall is not the coaxial cylinder through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                let sigma = if reversed { -1.0 } else { 1.0 };
                wall = Some((face.clone(), sigma));
            }
            Some(_) => ogeom_bail!(
                Construction,
                "the rim meets a face that is neither plane nor cylinder; \
                 that seat needs the marching blend machinery"
            ),
            None => ogeom_bail!(Dangling, "face refers to a surface not in this model"),
        }
    }
    let (Some((cap_face, up_raw)), Some((wall_face, sigma))) = (cap, wall) else {
        ogeom_bail!(
            Construction,
            "a revolved blend needs the edge shared by one planar cap and one \
             cylindrical wall"
        );
    };
    let up = up_raw / up_raw.magnitude();

    // Which side of the cap the wall extends: read at the edge itself, just
    // above and just below the cap's plane, where the wall stands on one
    // side only; a wall that runs on past the cap elsewhere (a boss standing
    // over a block's edge, its wall going down beside the block) stands on
    // both sides of the plane somewhere, so its far rings alone do not say.
    // Where the probes do not decide, the wall band's other ring does.
    // `tau` positive means away from the cap's outward side (the rim
    // configurations) and negative means alongside it, the concave seats.
    let local_tau = {
        let (curve, range) = edge_curve(model, edge, tol)?;
        let at = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
        let eps = rim_radius * 1e-2;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        let stands = |p: Point| -> OgeomResult<bool> {
            Ok(
                ogeom_algo::classify_on_face(model, &wall_face, p, deflection, tol)?
                    == ogeom_algo::Containment::In,
            )
        };
        match (stands(at + up * eps)?, stands(at - up * eps)?) {
            (true, false) => Some(-1.0),
            (false, true) => Some(1.0),
            _ => None,
        }
    };
    let tau = if let Some(tau) = local_tau {
        tau
    } else {
        let mut side = None;
        for e in explore(model, &wall_face, Filter::OfType(ShapeType::Edge))? {
            if same_occurrence(model, &e, edge, tol) {
                continue;
            }
            // Any circular edge of the wall at another height says which
            // side the wall extends: whether the ring survives as one
            // closed edge or as the arcs a boolean rebuilt it into.
            let Ok((curve, _)) = edge_curve(model, &e, tol) else {
                continue;
            };
            let Curve::Circle(c) = curve else {
                continue;
            };
            let lean = (c.circle().centre() - rim_centre).dot(up);
            if lean.abs() > tol.confusion() * 10.0 {
                side = Some(lean);
                break;
            }
        }
        let Some(lean) = side else {
            ogeom_bail!(
                Construction,
                "the rim's wall has no far ring to read its side from; the \
                 partial seat needs the marching blend machinery"
            );
        };
        if lean < 0.0 { 1.0 } else { -1.0 }
    };

    let frame = Frame::new(
        rim_centre,
        Direction::new(up, tol)?,
        rim_circle.frame().x(),
        tol,
    )?;
    let closes = circle_covered(model, solid, &rim_circle, &frame, tol)?;
    Ok(RevolvedSeat {
        centre: rim_centre,
        radius: rim_radius,
        up,
        x_ref: rim_circle.frame().x(),
        sigma,
        tau,
        cap_face,
        wall_face,
        closes,
    })
}

/// Whether the rims on `circle` run its whole turn between them: the
/// edges of `solid` on the circle that bound a wall on the circle's own
/// cylinder, each edge's arc taken as angles about `frame`, whose origin
/// is the circle's centre, and no gap left between them round the turn.
/// An edge on the circle between two faces of the cap's plane (where a
/// drum's top meets a neck fused onto it) is no rim.
fn circle_covered(
    model: &Model,
    solid: &Shape,
    circle: &Circle,
    frame: &Frame,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use core::f64::consts::TAU;
    let reach = tol.confusion() * 100.0;
    let angle_of = |p: Point| {
        let local = frame.to_local(p);
        local.y.atan2(local.x)
    };
    // Each arc as its start and its sweep, turning positively about the
    // frame's axis.
    let mut arcs: Vec<(f64, f64)> = Vec::new();
    let axis = circle.frame().z().vector();
    let mut on_walls: Vec<Shape> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(SurfaceGeometry::Cylinder(c)) = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
        else {
            continue;
        };
        let placement = face.transform(model.datums())?;
        let cylinder = c.cylinder();
        let axis_point = placement.apply(cylinder.frame().origin());
        let axis_z = placement.apply_vector(cylinder.frame().z().vector());
        let to_centre = circle.centre() - axis_point;
        let off_axis = (to_centre - axis_z * to_centre.dot(axis_z)).magnitude();
        if axis_z.cross(axis).magnitude() <= tol.angular()
            && off_axis <= reach
            && (cylinder.radius() - circle.radius()).abs() <= reach
        {
            on_walls.extend(ogeom_topo::explore_unique(model, &face, ShapeType::Edge)?);
        }
    }
    for edge in on_walls {
        let Ok((Curve::Circle(c), range)) = edge_curve(model, &edge, tol) else {
            continue;
        };
        let c = c.circle();
        if c.centre().distance(circle.centre()) > reach
            || (c.radius() - circle.radius()).abs() > reach
            || c.frame()
                .z()
                .vector()
                .cross(circle.frame().z().vector())
                .magnitude()
                > tol.angular()
        {
            continue;
        }
        if ogeom_algo::edge_vertices(model, &edge)?.is_some_and(|(a, b)| a.is_same(&b)) {
            return Ok(true);
        }
        let curve = Curve::Circle(CircleCurve::new(c));
        let [a0, am, a1] = [range.0, f64::midpoint(range.0, range.1), range.1]
            .map(|t| curve.point_at(t, tol).map(angle_of));
        let (a0, am, a1) = (a0?, am?, a1?);
        let (to_mid, to_end) = ((am - a0).rem_euclid(TAU), (a1 - a0).rem_euclid(TAU));
        arcs.push(if to_mid <= to_end {
            (a0.rem_euclid(TAU), to_end)
        } else {
            (a1.rem_euclid(TAU), (a0 - a1).rem_euclid(TAU))
        });
    }
    if arcs.is_empty() {
        return Ok(false);
    }
    arcs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let slack = tol.confusion() * 1e3 / circle.radius();
    let first = arcs[0].0;
    let mut covered = first;
    for (start, sweep) in arcs {
        if start > covered + slack {
            return Ok(false);
        }
        covered = covered.max(start + sweep);
    }
    Ok(covered >= first + TAU - slack)
}

/// Whether a blend along an open arc of a rim must run the whole turn.
///
/// An open arc is two different seats: a piece a boolean split off a full
/// rim, whose blend must run the whole turn, or a rim that genuinely stops,
/// a stadium's rounded end or a rounded corner. The solid's edges on the
/// rim's circle answer first: a rim that stops leaves a gap between them
/// somewhere round the turn. The two faces answer then: the rim goes on
/// round the arc's complement only where the wall stands just below it and
/// the cap lies just beside it, on the cap's own side of the rim. A wall
/// that carries on past where the cap ends (a boss standing over the edge
/// of a block) holds no rim there.
pub(crate) fn runs_whole_turn(
    model: &Model,
    seat: &RevolvedSeat,
    arc: &CircleCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Curve3d as _;
    if !seat.closes {
        return Ok(false);
    }
    let frame = seat.frame_at(seat.centre, tol)?;
    let angle_of = |t: f64| -> OgeomResult<f64> {
        let local = frame.to_local(arc.point_at(t, tol)?);
        Ok(local.y.atan2(local.x))
    };
    // The arc's own middle: which way round it runs about the seat's axis
    // depends on whether the cap faces along the circle's normal or against
    // it, and the parameter's middle is the arc's middle either way.
    let mid = angle_of(f64::midpoint(range.0, range.1))?;
    let complement = mid + core::f64::consts::PI;
    let dir = frame.x().vector() * complement.cos() + frame.y().vector() * complement.sin();
    let eps = seat.radius * 1e-2;
    let probe = seat.centre + dir * seat.radius - seat.up * (seat.tau * eps);
    let deflection = ogeom_mesh::Deflection {
        chord: eps * 0.1,
        ..ogeom_mesh::Deflection::default()
    };
    if ogeom_algo::classify_on_face(model, &seat.wall_face, probe, deflection, tol)?
        != ogeom_algo::Containment::In
    {
        return Ok(false);
    }
    // The cap lies inside the rim where the material is inside the wall and
    // the wall runs away from the cap's outward side, or where both are
    // turned; outside it otherwise.
    let outward = -seat.sigma * seat.tau;
    let beside = seat.centre + dir * (seat.radius + outward * eps);
    Ok(
        ogeom_algo::classify_on_face(model, &seat.cap_face, beside, deflection, tol)?
            == ogeom_algo::Containment::In,
    )
}

/// The flanks every revolved wedge shares: the band of the wall down to the
/// tangency ring, the annulus of the cap out to its own, and the three rings
/// bounding them. Only the face between the two tangency rings differs:
/// a quarter-tube for the fillet, a cone for the chamfer.
pub(crate) struct RevolvedFlanks {
    /// The band of the wall between the rim and `wall_ring`.
    pub wall_band: Shape,
    /// The annulus of the cap between the rim and `cap_ring`.
    pub annulus: Shape,
    /// The tangency ring on the wall, `wall_depth` from the rim.
    pub wall_ring: Shape,
    /// The tangency ring on the cap, at radius `cap_rho`.
    pub cap_ring: Shape,
}

/// Build the flanks: legs `wall_depth` down the wall and in to `cap_rho`
/// along the cap, each coincident with the solid's own face, aligned when
/// subtracting and opposed when fusing, which is what the melt needs.
pub(crate) fn revolved_flanks(
    model: &mut Model,
    seat: &RevolvedSeat,
    wall_depth: f64,
    cap_rho: f64,
    tol: Tolerances,
) -> OgeomResult<RevolvedFlanks> {
    let wall_level = seat.centre - seat.up * (seat.tau * wall_depth);
    let apex_ring = seat.ring(model, seat.centre, seat.radius, tol)?;
    let wall_ring = seat.ring(model, wall_level, seat.radius, tol)?;
    let cap_ring = seat.ring(model, seat.centre, cap_rho, tol)?;

    let wall_band = {
        let origin = seat.centre - seat.up * (wall_depth + 1.0);
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(seat.frame_at(origin, tol)?, seat.radius, tol)?,
            (0.0, 2.0 * wall_depth + 2.0),
        )?
        .into();
        let band = make_revolution_band(model, &surface, &wall_ring, &apex_ring, tol)?;
        if seat.sigma * seat.tau < 0.0 {
            band.reversed()
        } else {
            band
        }
    };

    let annulus = {
        let plane = Plane::through(seat.centre, Direction::new(seat.up * seat.tau, tol)?);
        let reach = (seat.radius + wall_depth + cap_rho) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        // Both rings turn about `up`. The wider is the outer boundary and
        // turns positively about the face's normal, the narrower the hole
        // and turns against it, so the face keeps itself on their left.
        let (wide, narrow) = if cap_rho > seat.radius {
            (&cap_ring, &apex_ring)
        } else {
            (&apex_ring, &cap_ring)
        };
        let with_normal = seat.tau > 0.0;
        let turned = |ring: &Shape, positive: bool| {
            if positive == with_normal {
                ring.clone()
            } else {
                ring.reversed()
            }
        };
        let outer = ogeom_algo::make_wire(model, &[turned(wide, true)], tol)?.shape;
        let inner = ogeom_algo::make_wire(model, &[turned(narrow, false)], tol)?.shape;
        let face = ogeom_algo::make_face(model, surface.clone(), &[outer, inner], tol)?.shape;
        let surface_id = {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "the face just built is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "the face holds no face data");
            };
            data.surface
        };
        for pedge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            let (curve, prange) = edge_curve(model, &pedge, tol)?;
            let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) else {
                ogeom_bail!(
                    Construction,
                    "an annulus edge has no closed-form pcurve on its plane"
                );
            };
            ogeom_algo::attach_pcurve(
                model,
                &pedge,
                pcurve,
                surface_id,
                ogeom_topo::Location::identity(),
                prange,
            )?;
        }
        face
    };

    Ok(RevolvedFlanks {
        wall_band,
        annulus,
        wall_ring,
        cap_ring,
    })
}

/// Whether `probe` sits strictly inside the face's trim.
pub(crate) fn on_face_side(
    model: &Model,
    face: &Shape,
    probe: ogeom_math::Point,
    deflection: ogeom_mesh::Deflection,
    tol: Tolerances,
) -> OgeomResult<bool> {
    Ok(
        ogeom_algo::classify_on_face(model, face, probe, deflection, tol)?
            == ogeom_algo::Containment::In,
    )
}

/// An edge along the segment from `from` to `to`, parameterized by arc
/// length, joining two existing vertices.
///
/// A wire chains through shared vertex *objects*, not through coincident
/// coordinates, which is why this takes the vertices and not just the
/// points.
pub(crate) fn segment_between(
    model: &mut Model,
    from: (&Shape, Point),
    to: (&Shape, Point),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let line = LineCurve::segment(from.1, to.1, tol)?;
    let curve = Curve::Line(line);
    let domain = curve.domain();
    Ok(make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape)
}

/// A face on `surface` bounded by `edges` in traversal order, with an exact
/// same-parameter pcurve attached to every edge.
///
/// [`ogeom_algo::make_face_with_pcurves`] with one wire: the blend keeps this
/// thin name because every wedge face is a single loop.
pub(crate) fn face_from_edges(
    model: &mut Model,
    surface: SurfaceGeometry,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &[edges.to_vec()], tol)?.shape)
}

/// Record every face of `built`'s result that no face of `solid` reaches
/// through its history as generated from `edge`: a blend or bevel face
/// comes from the tool the boolean took it from, and the edge it replaces
/// is where a caller naming faces through the operation finds it.
pub(crate) fn credit_new_faces(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    built: &mut ogeom_algo::Built,
) -> OgeomResult<()> {
    // A wedge set aside for later leaves the solid as it was: nothing new.
    if built.shape.is_same(solid) {
        return Ok(());
    }
    let inputs = ogeom_topo::explore_unique(model, solid, ShapeType::Face)?;
    for face in ogeom_topo::explore_unique(model, &built.shape, ShapeType::Face)? {
        let reached = inputs.iter().any(|g| {
            built
                .history
                .modified(g)
                .iter()
                .chain(built.history.generated(g))
                .any(|x| x.is_same(&face))
        });
        if !reached {
            built.history.generate(edge, face);
        }
    }
    Ok(())
}

/// The solid's boundary, prepared to say whether a rolling ball's centre
/// still stands in the material, as a seat running out past its edge's end
/// asks at a few points near that end.
///
/// The boundary answers against the true surfaces, and a face's trimming
/// rings are drawn only when a probe or its ray comes near the face, so
/// the cost follows the faces near the probe rather than the whole solid.
pub(crate) fn material_probe(
    model: &Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::SolidBoundary> {
    ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)
}

/// [`material_probe`] built the first time a seat's run-out asks, and
/// kept for the rest of that seat. A shape that bounds no inside (an open
/// shell) has no probe, and a seat on it does not run out.
#[derive(Default)]
pub(crate) struct LazyProbe(Option<Option<ogeom_algo::SolidBoundary>>);

impl LazyProbe {
    /// The probe of `solid`, or `None` where it bounds no inside.
    pub(crate) fn get(
        &mut self,
        model: &Model,
        solid: &Shape,
        tol: Tolerances,
    ) -> OgeomResult<Option<&ogeom_algo::SolidBoundary>> {
        if self.0.is_none() {
            self.0 = Some(match material_probe(model, solid, tol) {
                Ok(probe) => Some(probe),
                Err(ogeom_core::OgeomError::Construction(_)) => None,
                Err(e) => return Err(e),
            });
        }
        Ok(self.0.as_ref().and_then(Option::as_ref))
    }
}

/// The faces of `solid` following the crease of `edge` at its end `at`:
/// the faces other than `hosts` meeting the edge's vertex there.
fn following_faces(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    hosts: [&Shape; 2],
    at: Point,
) -> OgeomResult<Option<Vec<Shape>>> {
    let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(None);
    };
    let point_of = |v: &Shape| -> OgeomResult<Point> {
        let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Dangling, "an edge's vertex is not in this model");
        };
        Ok(v.transform(model.datums())?.apply(data.point))
    };
    let vertex = if point_of(&a)?.distance(at) <= point_of(&b)?.distance(at) {
        a
    } else {
        b
    };
    let mut following: Vec<Shape> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        if hosts.iter().any(|h| h.is_same(&face)) {
            continue;
        }
        let meets = explore(model, &face, Filter::OfType(ShapeType::Vertex))?
            .iter()
            .any(|v| v.is_same(&vertex));
        if meets {
            following.push(face);
        }
    }
    Ok(Some(following))
}

/// The faces of `solid` following the crease of `edge` at its end `at`:
/// the faces other than `hosts` meeting the edge's vertex there, each a
/// plane as a point on it and its outward unit normal. `None` where one of
/// them is not a plane, or none follows.
pub(crate) fn following_planes(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    hosts: [&Shape; 2],
    at: Point,
) -> OgeomResult<Option<Vec<(Shape, Point, Vector)>>> {
    let Some(faces) = following_faces(model, solid, edge, hosts, at)? else {
        return Ok(None);
    };
    let mut following: Vec<(Shape, Point, Vector)> = Vec::new();
    for face in faces {
        let Some(SurfaceGeometry::Plane(plane)) = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
        else {
            return Ok(None);
        };
        let placement = face.transform(model.datums())?;
        let origin = placement.apply(plane.plane().origin());
        let mut normal = placement.apply_vector(plane.plane().normal().vector());
        if face.orientation() == Orientation::Reversed {
            normal = -normal;
        }
        following.push((face, origin, normal));
    }
    Ok((!following.is_empty()).then_some(following))
}

/// A face a straight band runs out across, as its run-out reads it.
#[derive(Clone, Copy)]
enum Following {
    /// A plane: a point on it and its outward unit normal.
    Plane(Point, Vector),
    /// A cylinder holding the material within it: a point of its axis,
    /// the axis's unit direction and its radius.
    Boss(Point, Vector, f64),
    /// A cone holding the material within its one nappe: its apex, the
    /// axis's unit direction into that nappe, and `1 + tan^2` of its half
    /// angle. A point `x` stands within it where `q(x - apex) <= 0` and
    /// `(x - apex).axis >= 0`, `q` the quadratic form
    /// `q(v) = v.v - k (v.axis)^2` of [`Following::cone_form`].
    Taper(Point, Vector, f64),
}

impl Following {
    /// How far along the unit `outward` the line from `p` crosses the
    /// face's surface leaving the material: `None` where it runs along it,
    /// or misses a boss.
    fn exit(self, p: Point, outward: Vector, tol: Tolerances) -> Option<f64> {
        match self {
            Self::Plane(origin, normal) => {
                let rate = outward.dot(normal);
                (rate.abs() > tol.angular()).then(|| (origin - p).dot(normal) / rate)
            }
            Self::Boss(foot, axis, big) => {
                let flat = |v: Vector| v - axis * v.dot(axis);
                let (q, u) = (flat(p - foot), flat(outward));
                let a = u.dot(u);
                if a <= tol.angular() * tol.angular() {
                    return None;
                }
                let b = q.dot(u);
                let d = b.mul_add(b, -a * big.mul_add(-big, q.dot(q)));
                (d >= 0.0).then(|| (d.sqrt() - b) / a)
            }
            Self::Taper(apex, axis, k) => {
                // Along the line the form is `a s^2 + 2 b s + c`; the line
                // leaves the nappe where the form rises through zero, on the
                // nappe's side of the apex. A line whose stretch within the
                // nappe runs on without end never leaves it.
                let form = Self::cone_form(axis, k);
                let w = p - apex;
                let (a, b, c) = (form(outward, outward), form(outward, w), form(w, w));
                let rising = |s: f64| a.mul_add(s, b) >= 0.0 && (w + outward * s).dot(axis) >= 0.0;
                if a.abs() <= tol.angular() * tol.angular() {
                    if b.abs() <= tol.angular() {
                        return None;
                    }
                    let s = -c / (2.0 * b);
                    return rising(s).then_some(s);
                }
                let d = b.mul_add(b, -a * c);
                if d < 0.0 {
                    return None;
                }
                let root = d.sqrt();
                [(-b - root) / a, (-b + root) / a]
                    .into_iter()
                    .filter(|&s| rising(s))
                    .reduce(f64::max)
            }
        }
    }

    /// The symmetric form `(v, w) -> v.w - k (v.axis)(w.axis)` whose
    /// quadratic `q(v) = v.v - k (v.axis)^2` is negative within a cone about
    /// the unit `axis` through the origin, `k` being `1 + tan^2` of its half
    /// angle: `v.v - (v.axis)^2` is the squared distance off the axis, and
    /// `tan^2 (v.axis)^2` the squared radius there.
    fn cone_form(axis: Vector, k: f64) -> impl Fn(Vector, Vector) -> f64 {
        move |v: Vector, w: Vector| v.dot(w) - k * v.dot(axis) * w.dot(axis)
    }

    /// The farthest along the unit `outward` any point of the convex hull
    /// of `points` leaves a boss: seen along the axis, a point `q` off it
    /// leaves at `(sqrt(big^2 - w^2) - q.u) / |u|`, `u` the shadow of
    /// `outward` and `w` the point's offset square to it. That is a concave
    /// function of the point with no stationary point inside the hull (its
    /// gradient along `u` never vanishes), so its greatest value lies on an
    /// edge of the hull, and every edge is a segment between two of the
    /// points. Along a segment it is concave in the offset, greatest where
    /// its derivative vanishes or at an end. Points whose line misses the
    /// boss never stand in it. A cone's is [`Following::taper_hull_exit`].
    /// `None` for a plane, or for a boss the hull's lines all miss.
    fn hull_exit(self, points: &[Point], outward: Vector, tol: Tolerances) -> Option<f64> {
        if let Self::Taper(apex, axis, k) = self {
            return Self::taper_hull_exit(points, outward, (apex, axis, k), tol);
        }
        let Self::Boss(foot, axis, big) = self else {
            return None;
        };
        let flat = |v: Vector| v - axis * v.dot(axis);
        let u = flat(outward);
        let m = u.magnitude();
        if m <= tol.angular() {
            return None;
        }
        let along = u / m;
        let across = axis.cross(along);
        // Each point as (its offset across, its standing along).
        let seen: Vec<(f64, f64)> = points
            .iter()
            .map(|&p| {
                let q = flat(p - foot);
                (q.dot(across), q.dot(along))
            })
            .collect();
        let leave = |(w, l): (f64, f64)| -> Option<f64> {
            let h = big.mul_add(big, -w * w);
            (h >= 0.0).then(|| (h.sqrt() - l) / m)
        };
        let mut best: Option<f64> = None;
        let mut take = |s: Option<f64>| {
            if let Some(s) = s {
                best = Some(best.map_or(s, |b: f64| b.max(s)));
            }
        };
        for (i, &a) in seen.iter().enumerate() {
            take(leave(a));
            for &b in &seen[i + 1..] {
                let (dw, dl) = (b.0 - a.0, b.1 - a.1);
                if dw.abs() <= tol.confusion() {
                    continue;
                }
                // d/dt [sqrt(big^2 - w^2) - l] = -w dw / sqrt(big^2 - w^2)
                // - dl vanishes where w / sqrt(big^2 - w^2) = -dl / dw.
                let k = -dl / dw;
                let w = k * big / k.mul_add(k, 1.0).sqrt();
                let t = (w - a.0) / dw;
                if (0.0..=1.0).contains(&t) {
                    take(leave((w, t.mul_add(dl, a.1))));
                }
            }
        }
        best
    }

    /// [`Following::hull_exit`] for a cone. The material within a cone's
    /// nappe is convex, so the exit along `outward` is a concave function
    /// of the point among those whose lines meet the nappe: two such
    /// points' exits, mixed, place the mixed point's line within the nappe
    /// that far along. It falls by one per unit moved along `outward`, so
    /// its greatest value over the hull lies on an edge of it, a segment
    /// between two of the points. Along a segment `a + t e` the exit `s`
    /// is the far root of `q(a - apex + t e + s u) = 0`, and where it is
    /// greatest its derivative in `t` vanishes: `form(e, x) = 0` at the
    /// exit `x`, linear in `t` and `s`. Solved for `t` and put back, the
    /// form's root is a quadratic in `s`, whose roots give the turns of the
    /// segment the greatest exit may stand at, each then measured as it
    /// stands. A segment along which `form(e, e)` vanishes has its exit
    /// monotone between its ends.
    fn taper_hull_exit(
        points: &[Point],
        outward: Vector,
        (apex, axis, k): (Point, Vector, f64),
        tol: Tolerances,
    ) -> Option<f64> {
        let form = Self::cone_form(axis, k);
        let cone = Self::Taper(apex, axis, k);
        let mut best: Option<f64> = None;
        let mut take = |s: Option<f64>| {
            if let Some(s) = s {
                best = Some(best.map_or(s, |b: f64| b.max(s)));
            }
        };
        for (i, &a) in points.iter().enumerate() {
            take(cone.exit(a, outward, tol));
            for &b in &points[i + 1..] {
                let e = b - a;
                let ee = form(e, e);
                if ee.abs() <= tol.angular() * e.dot(e) {
                    continue;
                }
                // The turn `t = -(form(e, u) s + form(e, w)) / ee` at which
                // `form(e, x)` vanishes, as `x = w' + s u'`.
                let w = a - apex;
                let w_ = w - e * (form(e, w) / ee);
                let u_ = outward - e * (form(e, outward) / ee);
                let (qa, qb, qc) = (form(u_, u_), form(w_, u_), form(w_, w_));
                let roots: Vec<f64> = if qa.abs() <= tol.angular() * tol.angular() {
                    if qb.abs() <= tol.angular() {
                        Vec::new()
                    } else {
                        vec![-qc / (2.0 * qb)]
                    }
                } else {
                    let d = qb.mul_add(qb, -qa * qc);
                    if d < 0.0 {
                        Vec::new()
                    } else {
                        vec![(-qb - d.sqrt()) / qa, (-qb + d.sqrt()) / qa]
                    }
                };
                for s in roots {
                    let t = -form(e, outward).mul_add(s, form(e, w)) / ee;
                    if (0.0..=1.0).contains(&t) {
                        take(cone.exit(a + e * t, outward, tol));
                    }
                }
            }
        }
        best
    }
}

/// The faces a straight band runs out across past the end `at` of `edge`,
/// leaving it along `outward`: the planes following the crease there, or,
/// where `boss` allows, the one face following it where that is a cylinder
/// or a cone holding the material within it that the crease leaves the
/// material through, cutting the corner away obliquely. `None` elsewhere: a crease
/// running on into a boss (a pin standing at its end) runs into the
/// material.
fn run_out_faces(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    hosts: [&Shape; 2],
    (at, outward): (Point, Vector),
    boss: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<(Shape, Following)>>> {
    use ogeom_geom::Transformable as _;
    if let Some(planes) = following_planes(model, solid, edge, hosts, at)? {
        return Ok(Some(
            planes
                .into_iter()
                .map(|(face, origin, normal)| (face, Following::Plane(origin, normal)))
                .collect(),
        ));
    }
    if !boss {
        return Ok(None);
    }
    let Some(faces) = following_faces(model, solid, edge, hosts, at)? else {
        return Ok(None);
    };
    let [face] = faces.as_slice() else {
        return Ok(None);
    };
    let Some(surface) = model
        .node(face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .cloned()
    else {
        return Ok(None);
    };
    // The wall's axis, and what it follows as: a cylinder by its radius, a
    // cone by its apex and the axis into the nappe the vertex stands on.
    let (foot, axis, following) =
        match surface.transformed(&face.transform(model.datums())?, tol)? {
            SurfaceGeometry::Cylinder(placed) => {
                let cylinder = placed.cylinder();
                let foot = cylinder.axis().location;
                let axis = cylinder.axis().direction.vector();
                (foot, axis, Following::Boss(foot, axis, cylinder.radius()))
            }
            SurfaceGeometry::Cone(placed) => {
                let cone = placed.cone();
                let foot = cone.axis().location;
                let mut axis = cone.axis().direction.vector();
                let slope = cone.half_angle().tan();
                let apex = foot - axis * (cone.reference_radius() / slope);
                if (at - apex).dot(axis) < 0.0 {
                    axis = -axis;
                }
                (
                    foot,
                    axis,
                    Following::Taper(apex, axis, slope.mul_add(slope, 1.0)),
                )
            }
            _ => return Ok(None),
        };
    let off = at - foot;
    let radial = off - axis * off.dot(axis);
    let held = model
        .tolerance_of(face)?
        .map_or(tol.confusion(), ogeom_core::Tolerance::get)
        .max(tol.confusion() * 1e3);
    let Some(normal) = crate::marched::face_normal_within(model, face, at, held, tol)? else {
        return Ok(None);
    };
    let sign = if face.orientation() == Orientation::Reversed {
        -1.0
    } else {
        1.0
    };
    // The material lies within the cylinder or cone where the face's
    // outward normal points off the axis, and the crease leaves it there
    // where it runs along that normal.
    let normal = normal * sign;
    if normal.dot(radial) <= 0.0 || normal.dot(outward) <= tol.angular() {
        return Ok(None);
    }
    Ok(Some(vec![(face.clone(), following)]))
}

/// The plane of the face a straight crease runs into at its end `at`,
/// leaving it along `outward`: the one face following the crease there,
/// planar, and turned against `outward`, so its material overhangs the
/// crease's corner (a rib running up under a slanted roof). A band along
/// the crease is trimmed by that plane, or it cuts into what the face
/// bounds. The plane comes as a point on it and its outward unit normal.
pub(crate) fn overhang_at(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    hosts: [&Shape; 2],
    at: Point,
    outward: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Vector)>> {
    let Some(following) = following_planes(model, solid, edge, hosts, at)? else {
        return Ok(None);
    };
    let [(_, origin, normal)] = following.as_slice() else {
        return Ok(None);
    };
    Ok((outward.dot(*normal) < -tol.angular()).then_some((*origin, *normal)))
}

/// How far a convex straight seat's band runs out past the end `at` of
/// `edge`, along `outward`, where the faces meeting the crease there cut
/// the corner away obliquely: a floor giving way to a ramp, the wall the
/// ball rolls on carrying on above it, or a boss's cylinder or cone the
/// crease meets. The band carries on along the host
/// it stays tangent to until its whole section has left the material, and
/// the cut trims it against the face that follows.
///
/// `section` holds points of the band's section square to the crease at
/// `at`: its corners (the crease and the two contacts) and points between.
/// Each is carried past `at` and asked whether it stands in the material,
/// first a short way past the cap and then a stride at a time over four
/// radii. Where a ramp is the only face that follows, the material in
/// the section shrinks across a plane; the section lies within the
/// triangle of its corners, and a plane leaves a triangle last at a
/// corner, so the corners' answers are exact. Where a boss's cylinder or
/// cone is the one face that follows, the triangle leaves it last along one
/// of its sides, found in closed form ([`Following::hull_exit`]). A point
/// leaving through a face may cross it on the face's own edge, which counts
/// as on it.
///
/// `None` where the cap stands flush: nothing of the section stands in the
/// material just past `at` (a wall square to the crease, or one leaning
/// back over it), or the section does not leave the material within reach,
/// or a point of it leaves the material, comes back into it, or leaves
/// other than through a face meeting the crease's end there. That
/// last is a step or a ledge past the end rather than a ramp across it:
/// material no ball on this edge touches.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
pub(crate) fn run_out_reach(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    hosts: [&Shape; 2],
    (at, outward): (Point, Vector),
    section: &[Point],
    (radius, boss): (f64, bool),
    probe: &mut LazyProbe,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let Some(following) = run_out_faces(model, solid, edge, hosts, (at, outward), boss, tol)?
    else {
        return Ok(None);
    };
    if following.is_empty() {
        return Ok(None);
    }
    let Some(probe) = probe.get(model, solid, tol)? else {
        return Ok(None);
    };
    let held = |p: Point| -> OgeomResult<bool> {
        Ok(probe.holds(model, p, tol)? != ogeom_algo::Containment::Out)
    };
    let stride = radius / 8.0;
    // A short way past the cap, but well clear of the band the boundary
    // reads a point as on it within.
    let first = (stride / 8.0).max(tol.confusion() * 1e5).min(stride);
    let margin = radius * 0.25;
    let mut stations = vec![first];
    stations.extend((1..=32).map(|k| stride * f64::from(k)));
    // Per point, the last station it stood in the material at and the
    // first it had left it by.
    let mut last_held: Vec<Option<f64>> = vec![None; section.len()];
    let mut left_by: Vec<Option<f64>> = vec![None; section.len()];
    let mut clear_from: Option<f64> = None;
    let mut reached = false;
    for (k, &s) in stations.iter().enumerate() {
        for (i, &p) in section.iter().enumerate() {
            if held(p + outward * s)? {
                if left_by[i].is_some() {
                    return Ok(None);
                }
                last_held[i] = Some(s);
            } else if left_by[i].is_none() {
                left_by[i] = Some(s);
            }
        }
        if k == 0 && last_held.iter().all(Option::is_none) {
            return Ok(None);
        }
        if last_held
            .iter()
            .zip(&left_by)
            .all(|(h, l)| h.is_none() || l.is_some())
        {
            let from = *clear_from.get_or_insert(s);
            if s >= from + margin {
                reached = true;
                break;
            }
        }
    }
    if !reached {
        return Ok(None);
    }
    // Each point that stood in the material past the cap leaves it through
    // a face following the crease: its line crosses that face's plane
    // between the last station it stood in and the first it had left by,
    // give or take the first station's distance (the band the boundary
    // reads a point as on it within), and the crossing lies on the face.
    let deflection = ogeom_mesh::Deflection {
        chord: (radius * 1e-3).max(tol.confusion() * 1e2),
        ..ogeom_mesh::Deflection::default()
    };
    let mut leaves = 0.0_f64;
    for (i, &p) in section.iter().enumerate() {
        let (Some(held_to), Some(left)) = (last_held[i], left_by[i]) else {
            continue;
        };
        let mut through: Option<f64> = None;
        for (face, surface) in &following {
            let Some(s) = surface.exit(p, outward, tol) else {
                continue;
            };
            if s < held_to - first || s > left + first {
                continue;
            }
            let crossing = p + outward * s;
            if ogeom_algo::classify_on_face(model, face, crossing, deflection, tol)?
                != ogeom_algo::Containment::Out
                || on_face_edge(model, face, crossing, tol)?
            {
                through = Some(through.map_or(s, |t: f64| t.min(s)));
            }
        }
        let Some(s) = through else {
            return Ok(None);
        };
        leaves = leaves.max(s);
    }
    // A boss's cylinder or cone leaves the section last not at a corner of
    // it but where its hull leaves the boss last.
    for (_, surface) in &following {
        if let Some(s) = surface.hull_exit(section, outward, tol) {
            leaves = leaves.max(s);
        }
    }
    let reach = leaves + margin;
    for &p in section {
        if held(p + outward * reach)? {
            return Ok(None);
        }
    }
    Ok(Some(reach))
}

/// Whether `point` lies on an edge of `face`, within the edge's own
/// tolerance of its curve over its range. The face's trim, drawn as a
/// polyline in its chart from points of each edge's curve, stands off an
/// edge whose image in the chart bends more than the edge does in space,
/// so a point on the edge itself can read as outside the trim.
fn on_face_edge(model: &Model, face: &Shape, point: Point, tol: Tolerances) -> OgeomResult<bool> {
    for edge in ogeom_topo::explore_unique(model, face, ShapeType::Edge)? {
        let Ok((curve, (a, b))) = edge_curve(model, &edge, tol) else {
            continue;
        };
        let within = model
            .tolerance_of(&edge)?
            .map_or(tol.confusion(), ogeom_core::Tolerance::get)
            .max(tol.confusion());
        let gap = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(point)) };
        // The nearest of a few stations brackets the nearest point, and a
        // golden search closes on it. Every distance taken is to a point of
        // the edge, so a point read as on it is on it; one missed keeps the
        // band capped.
        let steps = 64;
        let at = |k: usize| {
            #[allow(clippy::cast_precision_loss, reason = "a station index")]
            let share = k as f64 / steps as f64;
            (b - a).mul_add(share, a)
        };
        let mut best = (f64::INFINITY, 0);
        for k in 0..=steps {
            let d = gap(at(k))?;
            if d < best.0 {
                best = (d, k);
            }
        }
        let (mut lo, mut hi) = (at(best.1.saturating_sub(1)), at((best.1 + 1).min(steps)));
        let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
        let mut nearest = best.0;
        for _ in 0..80 {
            let (x, y) = (hi - (hi - lo) * ratio, lo + (hi - lo) * ratio);
            let (dx, dy) = (gap(x)?, gap(y)?);
            nearest = nearest.min(dx).min(dy);
            if dx < dy {
                hi = y;
            } else {
                lo = x;
            }
        }
        if nearest <= within {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A blend's wedge set aside to be applied with others in one boolean.
pub(crate) struct Wedge {
    /// The wedge as a solid.
    pub(crate) solid: Shape,
    /// Fused on (a concave edge) rather than cut away.
    pub(crate) additive: bool,
    /// The edges it rounds: one, or two whose bands it joins.
    pub(crate) edges: Vec<Shape>,
    /// What the blend is refused as when the boolean cannot apply the
    /// wedge, where its builder knows a contact the melt does not resolve.
    pub(crate) melt_refusal: Option<&'static str>,
    /// Faces of the wedge credited to an edge other than the one it is set
    /// aside with: a wedge joining two bands credits each band's faces to
    /// its own edge.
    pub(crate) credit: Vec<(Shape, Shape)>,
}

/// Name the refusal the wedge last set aside is refused as, should the
/// boolean fail to apply it. Nothing when no caller is collecting.
pub(crate) fn note_melt_refusal(refusal: &'static str) {
    WEDGES.with(|held| {
        if let Some(last) = held.borrow_mut().as_mut().and_then(|w| w.last_mut()) {
            last.melt_refusal = Some(refusal);
        }
    });
}

thread_local! {
    /// Where [`apply_wedge`] sets its wedge aside, while a caller collects.
    static WEDGES: std::cell::RefCell<Option<Vec<Wedge>>> = const { std::cell::RefCell::new(None) };
}

/// Run `f` with every wedge [`apply_wedge`] builds set aside rather than
/// applied, and hand the wedges back with its result.
pub(crate) fn collecting_wedges<T>(f: impl FnOnce() -> T) -> (T, Vec<Wedge>) {
    struct Restore(Option<Vec<Wedge>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let outer = self.0.take();
            WEDGES.with(|held| *held.borrow_mut() = outer);
        }
    }
    let outer = WEDGES.with(|held| held.borrow_mut().replace(Vec::new()));
    let restore = Restore(outer);
    let out = f();
    let wedges = WEDGES
        .with(|held| held.borrow_mut().take())
        .unwrap_or_default();
    drop(restore);
    (out, wedges)
}

/// Sew the wedge's faces, demand a closed shell, and apply it to the solid:
/// subtracted on a convex edge, fused on a concave one. Either way the
/// history reads the same truth: the edge the blend replaces is gone.
pub(crate) fn apply_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: Option<&Shape>,
    faces: &[Shape],
    additive: bool,
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::Built> {
    let wedge = wedge_solid(model, faces, true, tol)?;
    apply_wedge_solid(model, solid, edge, &wedge, additive, tol)
}

/// Sew the wedge's faces into a solid, demanding a closed shell.
///
/// `fitted` says whether any of the faces' pcurves may be fitted, which
/// stands off its edge between its samples: every edge of such a wedge
/// states what its pcurves measure along it before the boolean reads it.
/// A wedge whose every pcurve is its edge's closed-form image stands
/// within rounding of its edges, under any tolerance an edge states, and
/// has nothing to measure.
pub(crate) fn wedge_solid(
    model: &mut Model,
    faces: &[Shape],
    fitted: bool,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let sewn = ogeom_algo::sew(model, faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the blend wedge did not close");
    }
    if fitted {
        ogeom_algo::state_pcurve_gaps(model, &sewn.shells[0], tol)?;
    }
    // Put together raw: a wedge is oriented by construction, and asks no
    // mesh to say which way it faces.
    model.add_solid(std::slice::from_ref(&sewn.shells[0]))
}

/// The wedge with what lies behind a plane cut away: the plane as a point
/// on it and its outward unit normal, the part of the wedge on the
/// normal's back side within `span` of `origin` gone.
pub(crate) fn trimmed_wedge(
    model: &mut Model,
    wedge: &Shape,
    (origin, normal): (Point, Vector),
    span: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let inward = Direction::new(-normal, tol)?;
    let across = {
        let x = if normal.x.abs() < 0.9 {
            Vector::new(1.0, 0.0, 0.0)
        } else {
            Vector::new(0.0, 1.0, 0.0)
        };
        Direction::new(x - normal * x.dot(normal), tol)?
    };
    let frame = Frame::new(origin, inward, across, tol)?;
    let corner = origin - (frame.x().vector() + frame.y().vector()) * span;
    let behind = ogeom_algo::make_box(
        model,
        Frame::new(corner, inward, across, tol)?,
        (span * 2.0, span * 2.0, span * 2.0),
        tol,
    )?
    .shape;
    let cut = ogeom_bool::cut(model, wedge, &behind, tol)?.shape;
    let solids = explore(model, &cut, Filter::OfType(ShapeType::Solid))?;
    let [solid] = solids.as_slice() else {
        ogeom_bail!(Construction, "the blend wedge trimmed by a face fell apart");
    };
    Ok(solid.clone())
}

/// Apply a wedge to the solid: subtracted on a convex edge, fused on a
/// concave one, or set aside where a caller collects wedges.
pub(crate) fn apply_wedge_solid(
    model: &mut Model,
    solid: &Shape,
    edge: Option<&Shape>,
    wedge: &Shape,
    additive: bool,
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::Built> {
    let wedge = wedge.clone();
    let set_aside = WEDGES.with(|held| {
        held.borrow_mut().as_mut().map(|wedges| {
            wedges.push(Wedge {
                solid: wedge.clone(),
                additive,
                edges: edge.into_iter().cloned().collect(),
                melt_refusal: None,
                credit: Vec::new(),
            });
        })
    });
    if set_aside.is_some() {
        let mut history = ogeom_algo::History::new();
        if let Some(edge) = edge {
            history.delete(edge);
        }
        return Ok(ogeom_algo::Built::new(solid.clone(), history));
    }
    let mut result = if additive {
        ogeom_bool::fuse(model, solid, &wedge, tol)?
    } else {
        ogeom_bool::cut(model, solid, &wedge, tol)?
    };
    if let Some(edge) = edge {
        result.history.delete(edge);
    }
    Ok(result)
}

/// A planar face over explicit corners, with `outward` as its plane normal.
///
/// The corners must be coplanar and `outward` perpendicular to them: the
/// callers know both by construction, which is why this takes the normal
/// instead of rediscovering it.
pub(crate) fn planar_face(
    model: &mut Model,
    corners: &[Point],
    outward: Vector,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let normal = Direction::new(outward, tol)?;
    let plane = Plane::through(corners[0], normal);
    let mut reach = 1.0_f64;
    for p in corners {
        reach = reach.max(p.distance(corners[0]) * 2.0);
    }
    let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
    // The ring keeps the face on its left about `outward`: corners given
    // turning the other way about it are walked in reverse. The turn is
    // summed about the first corner: the same sum as about any point, but
    // about a far one it is a difference of huge products.
    let turn = (0..corners.len()).fold(Vector::ZERO, |sum, i| {
        let (p, q) = (corners[i], corners[(i + 1) % corners.len()]);
        sum + (p - corners[0]).cross(q - corners[0])
    });
    let mut corners = corners.to_vec();
    if turn.dot(outward) < 0.0 {
        corners.reverse();
    }
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|p| make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::with_capacity(corners.len());
    for i in 0..corners.len() {
        let j = (i + 1) % corners.len();
        edges.push(segment_between(
            model,
            (&vertices[i], corners[i]),
            (&vertices[j], corners[j]),
            tol,
        )?);
    }
    face_from_edges(model, surface.into(), &edges, tol)
}

/// How far `face` reaches from the line through `from` along `leg`: the
/// greatest distance any of its edges stands back in that direction.
pub(crate) fn face_reach(
    model: &Model,
    face: &Shape,
    from: ogeom_math::Point,
    leg: ogeom_math::Vector,
    tol: Tolerances,
) -> OgeomResult<f64> {
    use ogeom_geom::Curve3d as _;
    let mut reach = 0.0_f64;
    for edge in ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let placement = edge.transform(model.datums())?;
        for k in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 16.0;
            let p = placement.apply(geometry.point_at(t, tol)?);
            reach = reach.max((p - from).dot(leg));
        }
    }
    Ok(reach)
}

/// Whether a rolling ball of `radius` sits on both faces along `edge`: at
/// stations along it, each face's contact is set back into the face by
/// `radius · tan(φ/2)` (φ the turn between the faces' outward normals,
/// taken on the surfaces where the edge runs), found on the face's surface
/// and classified against the face. A contact outside a face is a band
/// that runs past it; the refusal names the edge and the distance. Edges
/// that are not shared by exactly two faces, and faces meeting tangent,
/// are left to the blend itself. The contacts found come back, each with
/// its face and its setback, for a chain to check its bands against each
/// other.
///
/// `anywhere` asks only that the ball sit somewhere along the edge, for a
/// band no other band of the request is checked against. A face
/// narrowing along the edge (a triangle's side running down to its tip, a
/// pyramid's side toward its base) holds the ball along part of the edge
/// and lets the contact run off it elsewhere. Where the contact then
/// stands outside the solid, the face it ran off ended at the solid's
/// boundary, and the band is cut off there by the face that follows, as
/// a band running up under a roof is. Where it stands in the material,
/// the face ended against more material (a boss on the wall), and a band
/// there would cut under it, touching no face: that station refuses. A
/// ball that fits at no station is refused too.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where a contact falls outside its face (for `anywhere`, in the
/// material, or at every station).
pub(crate) fn ball_fits(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    anywhere: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Contact>> {
    use ogeom_geom::Surface as _;
    let mut contacts = Vec::new();
    let Ok((curve, range)) = edge_curve(model, edge, tol) else {
        return Ok(contacts);
    };
    // Each face holding the edge, with the edge's sense in that face's walk.
    let mut sides: Vec<(Shape, bool)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        for e in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if same_occurrence(model, &e, edge, tol) {
                // The walk as the face presents it: exploring under the
                // face composes its orientation into the edge's already.
                let reversed = e.orientation() == Orientation::Reversed;
                sides.push((face.clone(), reversed));
                break;
            }
        }
    }
    let [(face0, rev0), (face1, rev1)] = sides.as_slice() else {
        return Ok(contacts);
    };
    let normal_at = |face: &Shape, p: Point| -> OgeomResult<Option<(SurfaceGeometry, Vector)>> {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            return Ok(None);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(None);
        };
        use ogeom_geom::Transformable as _;
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let found = ogeom_algo::project_on_surface(&placed, p, 16, tol)?;
        let (u, v) = found.parameters;
        let Ok(n) = placed.normal_at(u, v, tol) else {
            return Ok(None);
        };
        let n = if face.orientation() == Orientation::Reversed {
            -n.vector()
        } else {
            n.vector()
        };
        Ok(Some((placed, n)))
    };
    // The planes a straight crease between planes runs into at its ends,
    // which trim its band: a contact behind one is cut off with the band's
    // end, not a band running past its face.
    let mut overhangs: Vec<(Point, Vector)> = Vec::new();
    let planar = |face: &Shape| {
        model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)))
    };
    if matches!(curve, Curve::Line(_)) && planar(face0) && planar(face1) {
        let tangent = {
            let d = curve.d1_at(range.0, tol)?;
            d / d.magnitude()
        };
        for (t, outward) in [(range.0, -tangent), (range.1, tangent)] {
            let at = curve.point_at(t, tol)?;
            if let Some(plane) = overhang_at(model, solid, edge, [face0, face1], at, outward, tol)?
            {
                overhangs.push(plane);
            }
        }
    }
    let trimmed_off = |p: Point| {
        overhangs
            .iter()
            .any(|(origin, normal)| (p - *origin).dot(*normal) < -tol.confusion())
    };
    // The edge's curve comes placed. Stations stay clear of its ends,
    // where a neighbouring face in a chain carries the band on.
    let mut missed: Option<String> = None;
    let mut held = false;
    let mut probe = LazyProbe::default();
    'stations: for fraction in [0.25, 0.5, 0.75] {
        let t = range.0 + (range.1 - range.0) * fraction;
        let p = curve.point_at(t, tol)?;
        let d = curve.d1_at(t, tol)?;
        let m = d.magnitude();
        if m <= tol.confusion() {
            continue;
        }
        let tangent = d / m;
        let (Some((s0, n0)), Some((s1, n1))) = (normal_at(face0, p)?, normal_at(face1, p)?) else {
            continue;
        };
        let turn = n0.dot(n1).clamp(-1.0, 1.0).acos();
        if turn <= tol.angular() * 1e3 {
            continue;
        }
        let setback = radius * (turn * 0.5).tan();
        // Into each face: across the edge, square to the normal, on the
        // side a hair's step from the edge lands inside the face. The walk
        // the orientation flags imply is the first guess; the face itself
        // settles it.
        let into = |face: &Shape,
                    surface: &SurfaceGeometry,
                    n: Vector,
                    reversed: bool|
         -> OgeomResult<Option<Vector>> {
            let walk = if reversed { -tangent } else { tangent };
            let d = n.cross(walk);
            let l = d.magnitude();
            if l <= 0.0 {
                return Ok(None);
            }
            let d = d / l;
            // A hair's step, and a hair whatever the radius: a step scaled
            // by a huge setback lands off the face on both sides.
            let hair = (setback * 1e-3).clamp(tol.confusion() * 10.0, tol.confusion() * 1e5);
            let lands = |step: Vector| -> OgeomResult<bool> {
                let (u, v) = ogeom_algo::project_on_surface(surface, p + step, 16, tol)?.parameters;
                chart_holds(model, face, ogeom_math::Point2::new(u, v), tol)
            };
            let near = [lands(d * hair)?, lands(-d * hair)?];
            // A hair's step lands in the face on both sides where the
            // face's trim turns back close to the edge, inside the chord
            // its chart is drawn with. The whole setback, landing in the
            // face on one side only, settles it there.
            if near == [true, true] {
                let far = [lands(d * setback)?, lands(-d * setback)?];
                if far == [false, true] {
                    return Ok(Some(-d));
                }
            }
            Ok(match near {
                [true, _] => Some(d),
                [false, true] => Some(-d),
                [false, false] => None,
            })
        };
        let (Some(d0), Some(d1)) = (into(face0, &s0, n0, *rev0)?, into(face1, &s1, n1, *rev1)?)
        else {
            continue;
        };
        for (face, surface, step) in [(face0, &s0, d0), (face1, &s1, d1)] {
            let guess = p + step * setback;
            let found = ogeom_algo::project_on_surface(surface, guess, 16, tol)?;
            let (u, v) = found.parameters;
            // On its own face, or on any face of the solid it has run on
            // to (a coplanar piece the face was split from, the next face
            // of a chain): anywhere on the boundary the ball can sit.
            if !chart_holds(model, face, ogeom_math::Point2::new(u, v), tol)?
                && !on_boundary(model, solid, found.point, tol)?
            {
                // A convex seat's band is cut off where its contact stands
                // behind a face the crease runs into: walking into the
                // first face from a convex edge goes behind the second.
                if d0.dot(n1) < 0.0 && trimmed_off(found.point) {
                    continue;
                }
                let refusal = format!(
                    "a blend of radius {radius} on the edge through {p:?} sets back {setback} \
                     across a face that does not reach so far"
                );
                // A solid that bounds no inside says nothing of where the
                // contact stands, and the station is let run off.
                let cut_off = match probe.get(model, solid, tol)? {
                    Some(probe) => {
                        probe.holds(model, found.point, tol)? == ogeom_algo::Containment::Out
                    }
                    None => true,
                };
                if !anywhere || !cut_off {
                    ogeom_bail!(Construction, "{refusal}");
                }
                missed = Some(refusal);
                continue 'stations;
            }
            contacts.push(Contact {
                face: face.node(),
                at: found.point,
                setback,
            });
        }
        held = true;
    }
    if let Some(refusal) = missed
        && !held
    {
        ogeom_bail!(Construction, "{refusal}");
    }
    Ok(contacts)
}

/// Where a rolling ball touches one of an edge's faces.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Contact {
    /// The face touched.
    pub face: ogeom_topo::TShapeId,
    /// The point of contact.
    pub at: Point,
    /// How far back from the edge it stands.
    pub setback: f64,
}

/// Refuse a chain whose bands overlap across a face they share.
///
/// Each edge's band covers the strip of a face between the edge and its
/// contact line. Two edges of one face that do not meet at a vertex (the
/// two long edges of a narrow strip, a wall's two rims) each need their own
/// strip; where one's contact stands inside the other's strip, the two
/// balls overlap and no pair of blends is the rounding. Edges meeting at a
/// vertex are corner mates, whose bands meet there by design.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where two bands overlap.
pub(crate) fn bands_clear(
    model: &Model,
    edges: &[Shape],
    contacts: &[Vec<Contact>],
    tol: Tolerances,
) -> OgeomResult<()> {
    let vertices = |edge: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
        Ok(ogeom_algo::edge_vertices(model, edge)?
            .map(|(a, b)| vec![a.node(), b.node()])
            .unwrap_or_default())
    };
    // Each edge's box: a contact further from it than a strip is further
    // from the edge too, and the scan along the edge is not needed.
    let mut boxes = Vec::with_capacity(edges.len());
    for edge in edges {
        boxes.push(ogeom_algo::shape_bounds(model, edge, tol)?);
    }
    for (i, first) in edges.iter().enumerate() {
        for (j, second) in edges.iter().enumerate() {
            if i == j {
                continue;
            }
            let (ends, others) = (vertices(first)?, vertices(second)?);
            if ends.iter().any(|v| others.contains(v)) {
                continue;
            }
            let Ok((curve, range)) = edge_curve(model, second, tol) else {
                continue;
            };
            for contact in &contacts[i] {
                // The second edge's own setback on this face, if it has one.
                let Some(strip) = contacts[j]
                    .iter()
                    .filter(|c| c.face == contact.face)
                    .map(|c| c.setback)
                    .reduce(f64::max)
                else {
                    continue;
                };
                if boxes[j].distance_to(contact.at) >= strip {
                    continue;
                }
                let away = distance_to_edge(&curve, range, contact.at, tol)?;
                if away < strip - tol.confusion() * 1e3 {
                    ogeom_bail!(
                        Construction,
                        "two blends of the chain overlap across a face they share: a \
                         contact stands {away} from the other edge, inside its band of {strip}"
                    );
                }
            }
        }
    }
    Ok(())
}

/// The distance from `p` to a curve over `range`: the nearest of a scan,
/// refined in the interval around it.
fn distance_to_edge(
    curve: &Curve,
    range: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    use ogeom_geom::Curve3d as _;
    const SCAN: u32 = 64;
    let at = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(p)) };
    let step = (range.1 - range.0) / f64::from(SCAN);
    let mut best = (range.0, at(range.0)?);
    for k in 1..=SCAN {
        let t = range.0 + step * f64::from(k);
        let d = at(t)?;
        if d < best.1 {
            best = (t, d);
        }
    }
    let (mut lo, mut hi) = ((best.0 - step).max(range.0), (best.0 + step).min(range.1));
    for _ in 0..60 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if at(a)? < at(b)? {
            hi = b;
        } else {
            lo = a;
        }
    }
    Ok(at(f64::midpoint(lo, hi))?.min(best.1))
}

thread_local! {
    /// Face meshes kept while a caller checks many edges against one
    /// unchanging solid: a face with a hundred holes is otherwise meshed
    /// once per edge on it.
    static FACE_MESHES: std::cell::RefCell<Option<ogeom_core::FastMap<ogeom_topo::SameKey, std::rc::Rc<ogeom_topo::Triangulation>>>> =
        const { std::cell::RefCell::new(None) };
}

/// Run `f` with face meshes kept from one call of [`chart_holds`] to the
/// next. The solid must not change while it runs.
pub(crate) fn keeping_face_meshes<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(
        Option<ogeom_core::FastMap<ogeom_topo::SameKey, std::rc::Rc<ogeom_topo::Triangulation>>>,
    );
    impl Drop for Restore {
        fn drop(&mut self) {
            let outer = self.0.take();
            FACE_MESHES.with(|held| *held.borrow_mut() = outer);
        }
    }
    let outer = FACE_MESHES.with(|held| held.borrow_mut().replace(ogeom_core::FastMap::default()));
    let _restore = Restore(outer);
    f()
}

/// A face's mesh at the default deflection, kept where a caller asked.
fn face_mesh(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<std::rc::Rc<ogeom_topo::Triangulation>> {
    let key = ogeom_topo::SameKey(face.clone());
    if let Some(mesh) =
        FACE_MESHES.with(|held| held.borrow().as_ref().and_then(|m| m.get(&key).cloned()))
    {
        return Ok(mesh);
    }
    let mesh = std::rc::Rc::new(ogeom_mesh::triangulate_face(
        model,
        face,
        ogeom_mesh::Deflection::default(),
        tol,
    )?);
    FACE_MESHES.with(|held| {
        if let Some(meshes) = held.borrow_mut().as_mut() {
            meshes.insert(key, std::rc::Rc::clone(&mesh));
        }
    });
    Ok(mesh)
}

/// Whether a chart point lies in a face's region: inside one of the
/// triangles its own mesh lays over its chart, or on one's side.
fn chart_holds(
    model: &Model,
    face: &Shape,
    at: ogeom_math::Point2,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mesh = face_mesh(model, face, tol)?;
    let chart = |i: u32| {
        let (u, v) = mesh.parameters[i as usize];
        ogeom_math::Point2::new(u, v)
    };
    let span = mesh
        .parameters
        .iter()
        .fold(0.0_f64, |m, (u, v)| m.max(u.abs()).max(v.abs()))
        .max(1.0);
    let eps = span * 1e-9;
    Ok(mesh.triangles.iter().any(|t| {
        let [a, b, c] = t.map(chart);
        let cross = |p: ogeom_math::Point2, q: ogeom_math::Point2| {
            (q.x - p.x) * (at.y - p.y) - (q.y - p.y) * (at.x - p.x)
        };
        let (d1, d2, d3) = (cross(a, b), cross(b, c), cross(c, a));
        let neg = d1 < -eps || d2 < -eps || d3 < -eps;
        let pos = d1 > eps || d2 > eps || d3 > eps;
        !(neg && pos)
    }))
}

/// Whether a point lies on some face of `solid`: on its surface, within
/// the model's reach, and inside its chart region.
fn on_boundary(model: &Model, solid: &Shape, point: Point, tol: Tolerances) -> OgeomResult<bool> {
    use ogeom_geom::Transformable as _;
    let reach = tol.confusion() * 100.0;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let found = ogeom_algo::project_on_surface(&placed, point, 16, tol)?;
        if found.distance > reach {
            continue;
        }
        let (u, v) = found.parameters;
        if chart_holds(model, &face, ogeom_math::Point2::new(u, v), tol)? {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use ogeom_core::FastMap;

    use ogeom_core::Tolerances;
    use ogeom_math::{Frame, Point};
    use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    /// A triangle's farthest exit from a cone along a direction square to
    /// its axis, against the greatest exit of points spread over the
    /// triangle. The triangle straddles the plane through the axis along
    /// that direction, where a point's exit is greatest, so the greatest
    /// stands inside a side, past every corner's: the corners alone fall
    /// short of it.
    #[test]
    fn a_cone_s_farthest_exit_over_a_triangle_is_its_hull_exit() {
        use ogeom_math::Vector;
        // Apex at the origin, the nappe up `z`, a tenth's slope: radius
        // `z / 10`.
        let slope = 0.1_f64;
        let cone = super::Following::Taper(Point::ORIGIN, Vector::Z, slope.mul_add(slope, 1.0));
        let outward = Vector::X;
        let corners = [
            Point::new(2.3, -0.3, 25.0),
            Point::new(2.3, 0.4, 25.0),
            Point::new(2.35, 0.1, 24.6),
        ];
        let got = cone.hull_exit(&corners, outward, T).unwrap();
        // Seen along the axis, the point leaves where its line meets the
        // circle of the cone's radius at its height.
        let leave = |p: Point| (p.z * slope).mul_add(p.z * slope, -p.y * p.y).sqrt() - p.x;
        let n = 600;
        let mut best = f64::NEG_INFINITY;
        for i in 0..=n {
            for j in 0..=(n - i) {
                let (a, b) = (f64::from(i) / f64::from(n), f64::from(j) / f64::from(n));
                let p = corners[0] + (corners[1] - corners[0]) * a + (corners[2] - corners[0]) * b;
                best = best.max(leave(p));
            }
        }
        let at_corners = corners
            .iter()
            .map(|&p| leave(p))
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            best - at_corners > 1e-3,
            "{best} against the corners' {at_corners}"
        );
        // The spread lies in the triangle, so its best is no greater than
        // the hull's. The greatest stands on the first side, whose points
        // the spread takes 1.2e-3 apart, and the exit's second derivative
        // across the axis there is about one over the radius, 0.4: the
        // spread misses it by under 1e-7.
        assert!(
            got >= best - 1e-12 && got - best < 1e-5,
            "{got} against {best}"
        );
    }

    /// How many edges two faces of `shape` walk the same way.
    fn walked_one_way(model: &Model, shape: &Shape) -> usize {
        let mut walks: FastMap<_, Vec<Orientation>> = FastMap::default();
        for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
            for wire in model.children_of(&face).unwrap() {
                for edge in model.children_of(&wire).unwrap() {
                    if !model
                        .node(&edge)
                        .unwrap()
                        .data()
                        .as_edge()
                        .unwrap()
                        .degenerate
                    {
                        walks
                            .entry(edge.node())
                            .or_default()
                            .push(edge.orientation());
                    }
                }
            }
        }
        walks
            .values()
            .filter(|w| w.len() == 2 && w[0] == w[1])
            .count()
    }

    /// The edge of `shape` whose middle stands nearest `at`.
    fn edge_at(model: &Model, shape: &Shape, at: Point) -> Shape {
        let off = |e: &Shape| {
            let b = ogeom_algo::shape_bounds(model, e, T).unwrap();
            b.low().unwrap().midpoint(b.high().unwrap()).distance(at)
        };
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .min_by(|a, b| off(a).total_cmp(&off(b)))
            .unwrap()
    }

    /// The wedges a straight edge's fillet and chamfer and a drum rim's
    /// fillet take away walk each edge between two of their faces once each
    /// way, and hold the volume the blend removes: a quarter round's
    /// section, or a half square's, along the edge or round the rim at the
    /// section's centroid.
    #[test]
    fn blend_wedges_walk_each_edge_once_each_way() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
            .unwrap()
            .shape;
        let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 4.0, T)
            .unwrap()
            .shape;
        let ridge = edge_at(&model, &block, Point::new(5.0, 0.0, 4.0));
        let rim = edge_at(&model, &drum, Point::new(0.0, 0.0, 4.0));
        let pi = core::f64::consts::PI;
        let section = 1.0 - pi / 4.0;
        let centroid = (5.0 / 6.0 - pi / 4.0) / section;
        let cases = [
            (&block, &ridge, false, section * 10.0),
            (&block, &ridge, true, 0.5 * 10.0),
            (&drum, &rim, false, 2.0 * pi * (3.0 - centroid) * section),
        ];
        for (solid, edge, bevel, removed) in cases {
            let (made, wedges) = super::collecting_wedges(|| {
                if bevel {
                    crate::chamfer_edge(&mut model, solid, edge, 1.0, T)
                } else {
                    crate::fillet_edge(&mut model, solid, edge, 1.0, T)
                }
            });
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("one wedge, not {}", wedges.len());
            };
            assert_eq!(walked_one_way(&model, &wedge.solid), 0, "bevel {bevel}");
            let volume = ogeom_algo::volume_properties(
                &model,
                &wedge.solid,
                ogeom_mesh::Deflection::default(),
                T,
            )
            .unwrap()
            .mass;
            assert!(
                (volume - removed).abs() < removed * 1e-6,
                "bevel {bevel}: {volume} against {removed}"
            );
        }
    }

    /// A straight seat's wedge, on a convex edge and a concave one, has
    /// every pcurve its edge's closed-form image: each edge's pcurves stand
    /// within the tolerance it states, so the wedge is built without
    /// measuring them.
    #[test]
    fn straight_seat_wedges_hold_their_pcurves_within_their_edges() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
            .unwrap()
            .shape;
        let step = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 3.0, 8.0), T)
            .unwrap()
            .shape;
        let ell = ogeom_bool::fuse(&mut model, &block, &step, T)
            .unwrap()
            .shape;
        let ridge = edge_at(&model, &block, Point::new(5.0, 0.0, 4.0));
        let crease = edge_at(&model, &ell, Point::new(5.0, 3.0, 4.0));
        for (solid, edge) in [(&block, &ridge), (&ell, &crease)] {
            let (made, wedges) =
                super::collecting_wedges(|| crate::fillet_edge(&mut model, solid, edge, 1.0, T));
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("one wedge, not {}", wedges.len());
            };
            for e in explore_unique(&model, &wedge.solid, ShapeType::Edge).unwrap() {
                let stated = model
                    .node(&e)
                    .unwrap()
                    .data()
                    .as_edge()
                    .unwrap()
                    .tolerance
                    .get();
                let gap = ogeom_algo::edge_pcurve_gap(&model, &e, T).unwrap().unwrap();
                assert!(
                    gap <= stated,
                    "a pcurve stands {gap} off an edge stating {stated}"
                );
            }
        }
    }

    /// The volume, validity and face, edge and vertex counts of a solid.
    fn summary(model: &Model, shape: &Shape) -> (f64, bool, [usize; 3]) {
        let volume =
            ogeom_algo::volume_properties(model, shape, ogeom_mesh::Deflection::default(), T)
                .unwrap()
                .mass;
        let valid = ogeom_algo::check(model, shape, T).unwrap().is_valid();
        let count = |kind| explore_unique(model, shape, kind).unwrap().len();
        (
            volume,
            valid,
            [
                count(ShapeType::Face),
                count(ShapeType::Edge),
                count(ShapeType::Vertex),
            ],
        )
    }

    /// A straight band's wedge cuts the block it was built for by trimming
    /// the faces it sits on. The same wedge on the block with a blind hole
    /// reaching up into it from below, whose wall its blend would cross, is
    /// left to the general fuse. Either way the solid is the general fuse's.
    #[test]
    fn a_seated_wedge_is_cut_by_trimming_only_where_nothing_reaches_into_it() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
            .unwrap()
            .shape;
        let ridge = edge_at(&model, &block, Point::new(5.0, 0.0, 4.0));
        let (made, wedges) =
            super::collecting_wedges(|| crate::fillet_edge(&mut model, &block, &ridge, 1.0, T));
        made.unwrap();
        let [wedge] = &wedges[..] else {
            panic!("one wedge, not {}", wedges.len());
        };
        // Up to a tenth below the ridge, a tenth in from the front: inside
        // the corner the band rounds away.
        let at = Frame::new(
            Point::new(5.0, 0.1, -1.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let drill = ogeom_algo::make_cylinder(&mut model, at, 0.05, 4.9, T)
            .unwrap()
            .shape;
        let drilled = ogeom_bool::cut(&mut model, &block, &drill, T)
            .unwrap()
            .shape;
        for (solid, seated) in [(&block, 1), (&drilled, 0)] {
            let mut general = model.clone();
            let want = ogeom_bool::cut(&mut general, solid, &wedge.solid, T).unwrap();
            let (got, work) = ogeom_bool::work::counted(|| {
                ogeom_bool::cut_seated(&mut model, solid, &wedge.solid, T).unwrap()
            });
            assert_eq!(work.seated_cuts, seated);
            let (want, got) = (summary(&general, &want.shape), summary(&model, &got.shape));
            assert!(got.1, "the cut is valid");
            assert_eq!((got.1, got.2), (want.1, want.2));
            assert!(
                (got.0 - want.0).abs() <= want.0 * 1e-9,
                "{} against {}",
                got.0,
                want.0
            );
        }
    }

    /// A stadium prism from `z0` up `height`: two runs of `length` along
    /// x joined by semicircles of radius `r`, centred on (0, 0) and
    /// (`length`, 0).
    fn stadium(model: &mut Model, length: f64, r: f64, z0: f64, height: f64) -> Shape {
        use ogeom_math::{Circle, Direction, Plane, Vector};
        let corner = |x: f64, y: f64| Point::new(x, y, z0);
        let [tl, tr, br, bl] = [(0.0, r), (length, r), (length, -r), (0.0, -r)]
            .map(|(x, y)| ogeom_algo::make_vertex(model, corner(x, y)).shape);
        let arc = |model: &mut Model, cx: f64, x: Direction, from: &Shape, to: &Shape| {
            let frame = Frame::new(corner(cx, 0.0), Direction::Z, x, T).unwrap();
            let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(
                Circle::new(frame, r, T).unwrap(),
            ));
            ogeom_algo::make_edge_between(model, curve, (0.0, core::f64::consts::PI), from, to, T)
                .unwrap()
                .shape
        };
        let top = super::segment_between(model, (&tl, corner(0.0, r)), (&tr, corner(length, r)), T)
            .unwrap();
        let right = arc(model, length, -Direction::Y, &br, &tr);
        let bottom =
            super::segment_between(model, (&br, corner(length, -r)), (&bl, corner(0.0, -r)), T)
                .unwrap();
        let left = arc(model, 0.0, Direction::Y, &tl, &bl);
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(corner(0.0, 0.0), Direction::Z),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face_with_pcurves(
            model,
            plane.into(),
            &[vec![left, bottom.reversed(), right, top.reversed()]],
            T,
        )
        .unwrap()
        .shape;
        ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
            .unwrap()
            .shape
    }

    /// The wedges an arc of a rim takes away or fills walk each edge between
    /// two of their faces once each way, whichever side of the rim the cap
    /// and the material stand on, and hold the volume the blend removes or
    /// adds: a quarter round's section swept round the half turn at its
    /// centroid's radius, inside the rim where the wall faces out and
    /// outside it where the wall faces the axis or the rim is concave.
    #[test]
    fn rim_arc_wedges_walk_each_edge_once_each_way() {
        let mut model = Model::new();
        let (length, r) = (10.0, 5.0);
        let post = stadium(&mut model, length, r, 0.0, 4.0);
        let block = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-10.0, -10.0, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (30.0, 20.0, 4.0),
            T,
        )
        .unwrap()
        .shape;
        let bore = stadium(&mut model, length, r, -1.0, 6.0);
        let slot = ogeom_bool::cut(&mut model, &block, &bore, T).unwrap().shape;
        let plate = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-10.0, -10.0, -4.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (30.0, 20.0, 4.0),
            T,
        )
        .unwrap()
        .shape;
        let boss = ogeom_bool::fuse(&mut model, &plate, &post, T)
            .unwrap()
            .shape;
        let pi = core::f64::consts::PI;
        let section = 1.0 - pi / 4.0;
        let centroid = (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi);
        let inside = pi * (r - centroid) * section;
        let outside = pi * (r + centroid) * section;
        let at = |z: f64| Point::new(length + r / 2.0, 0.0, z);
        let cases = [
            ("post top", &post, at(4.0), inside),
            ("post foot", &post, at(0.0), inside),
            ("slot top", &slot, at(4.0), outside),
            ("slot foot", &slot, at(0.0), outside),
            ("boss foot", &boss, at(0.0), outside),
        ];
        for (name, solid, near, expected) in cases {
            let rim = edge_at(&model, solid, near);
            let (made, wedges) =
                super::collecting_wedges(|| crate::fillet_edge(&mut model, solid, &rim, 1.0, T));
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("{name}: one wedge, not {}", wedges.len());
            };
            assert_eq!(walked_one_way(&model, &wedge.solid), 0, "{name}");
            let volume = ogeom_algo::volume_properties(
                &model,
                &wedge.solid,
                ogeom_mesh::Deflection::default(),
                T,
            )
            .unwrap()
            .mass;
            assert!(
                (volume - expected).abs() < expected * 1e-6,
                "{name}: {volume} against {expected}"
            );
        }
    }

    /// The edge of `shape` whose curve is the first `kind` accepts.
    fn edge_on(model: &Model, shape: &Shape, kind: impl Fn(&ogeom_geom::Curve) -> bool) -> Shape {
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| super::edge_curve(model, e, T).is_ok_and(|(curve, _)| kind(&curve)))
            .unwrap()
    }

    /// The volume of a solid, through the mesh where a face is fitted.
    fn volume(model: &Model, solid: &Shape) -> f64 {
        ogeom_algo::volume_properties(model, solid, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    }

    /// A bore down a ball beside its pole, rounded where it leaves the far
    /// side: the ball rolling round the rim touches the sphere along a rail
    /// that passes over the pole, so the sphere's leg is the band from the
    /// rail to the pole with the rim cut from it. Its seam to the pole runs
    /// clear of the rim however the bore is turned about the ball's axis:
    /// the wedge meshes closed to what it measures, the leg's mesh covers
    /// its area, and the blend is valid and measures the same at each turn.
    #[test]
    fn a_pole_leg_meshes_to_its_area_as_the_bore_turns() {
        let reference = pole_leg_blend(0.0);
        for turn in [0.9, 3.2, 5.6] {
            let whole = pole_leg_blend(turn);
            assert!(
                (whole - reference).abs() < reference * 1e-6,
                "turn {turn}: {whole} against {reference}"
            );
        }
    }

    /// The bore of [`a_pole_leg_meshes_to_its_area_as_the_bore_turns`]
    /// turned so the rim's middle is its point nearest the pole: the side
    /// probe there stands a tenth of a millimetre either side of the rim,
    /// half a millimetre from the pole, and the ball's touch on the sphere
    /// lies across the pole from it. The blend seats at another point of
    /// the rim and measures what it does at the other turns.
    #[test]
    fn a_pole_leg_blends_with_the_rim_s_middle_nearest_the_pole() {
        let reference = pole_leg_blend(0.0);
        let whole = pole_leg_blend(4.7);
        assert!(
            (whole - reference).abs() < reference * 1e-6,
            "{whole} against {reference}"
        );
    }

    /// The far rim of a bore of radius 1.5 down a ball of radius 10, its
    /// axis 2 from the ball's and `turn` round it, rounded at radius 1:
    /// checks the wedge's mesh, the sphere leg's mesh and the blend's
    /// validity, and gives the blend's volume.
    fn pole_leg_blend(turn: f64) -> f64 {
        use ogeom_math::Direction;
        let chord = 1e-2;
        let deflection = ogeom_mesh::Deflection::with_chord(chord).unwrap();
        let (c, s) = (turn.cos(), turn.sin());
        let mut model = Model::new();
        let ball = ogeom_algo::make_sphere(&mut model, Frame::WORLD, 10.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(
            Point::new(2.0 * c, 2.0 * s, -20.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let bore = ogeom_algo::make_cylinder(&mut model, frame, 1.5, 40.0, T)
            .unwrap()
            .shape;
        let drilled = ogeom_bool::cut(&mut model, &ball, &bore, T).unwrap().shape;
        // The rim on the ball's top, on the bore's side away from the axis.
        let top = (100.0_f64 - 3.5 * 3.5).sqrt();
        let rim = edge_at(&model, &drilled, Point::new(3.5 * c, 3.5 * s, top));
        let (made, wedges) =
            super::collecting_wedges(|| crate::fillet_edge(&mut model, &drilled, &rim, 1.0, T));
        made.unwrap();
        let [wedge] = &wedges[..] else {
            panic!("turn {turn}: one wedge, not {}", wedges.len());
        };
        let exact = volume(&model, &wedge.solid);
        let mesh = ogeom_mesh::triangulate(&model, &wedge.solid, deflection, T).unwrap();
        assert!(mesh.is_closed(), "turn {turn}");
        assert!(
            (mesh.volume() - exact).abs() < exact * 0.05,
            "turn {turn}: wedge mesh {} against {exact}",
            mesh.volume()
        );
        let leg = explore_unique(&model, &wedge.solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|face| {
                model
                    .node(face)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface))
                    .is_some_and(|g| matches!(g, ogeom_geom::SurfaceGeometry::Sphere(_)))
            })
            .unwrap();
        let area = ogeom_algo::surface_properties(&model, &leg, deflection, T)
            .unwrap()
            .mass;
        let leg_mesh = ogeom_mesh::triangulate_face(&model, &leg, deflection, T).unwrap();
        let meshed: f64 = leg_mesh
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| leg_mesh.positions[i as usize]);
                (b - a).cross(c - a).magnitude() * 0.5
            })
            .sum();
        // Flat facets sagging by at most the chord under a sphere of radius
        // 10 fall short of its area by under twice chord / 10.
        assert!(
            (meshed - area).abs() < area * 2.0 * chord / 10.0,
            "turn {turn}: leg mesh {meshed} against {area}"
        );
        let blended = crate::fillet_edge(&mut model, &drilled, &rim, 1.0, T)
            .unwrap()
            .shape;
        assert!(
            ogeom_algo::check(&model, &blended, T).unwrap().is_valid(),
            "turn {turn}"
        );
        volume(&model, &blended)
    }

    /// Where two bores cross inside a block, the wedge of the band rounding
    /// each edge where they meet meshes to what it measures. Round a closed
    /// loop each leg is the patch of one bore's wall between two closed
    /// rings, both on one branch of the wall's chart, so the mesh covers
    /// the patch and not the wall round between them.
    #[test]
    fn a_band_round_crossing_bores_meshes_to_its_volume() {
        use ogeom_math::Direction;
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 20.0), T)
            .unwrap()
            .shape;
        let bore = |model: &mut Model, at: Point, axis: Direction, x: Direction, r: f64| {
            let frame = Frame::new(at, axis, x, T).unwrap();
            ogeom_algo::make_cylinder(model, frame, r, 22.0, T)
                .unwrap()
                .shape
        };
        let a = bore(
            &mut model,
            Point::new(-1.0, 10.0, 10.0),
            Direction::X,
            Direction::Y,
            4.0,
        );
        let b = bore(
            &mut model,
            Point::new(10.0, -1.0, 12.0),
            Direction::Y,
            Direction::Z,
            3.0,
        );
        let once = ogeom_bool::cut(&mut model, &block, &a, T).unwrap().shape;
        let bored = ogeom_bool::cut(&mut model, &once, &b, T).unwrap().shape;
        let meeting: Vec<Shape> = explore_unique(&model, &bored, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                super::edge_curve(&model, e, T)
                    .is_ok_and(|(curve, _)| matches!(curve, ogeom_geom::Curve::BSpline(_)))
            })
            .collect();
        assert_eq!(meeting.len(), 3);
        for edge in &meeting {
            let (made, wedges) =
                super::collecting_wedges(|| crate::fillet_edge(&mut model, &bored, edge, 1.0, T));
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("one wedge, not {}", wedges.len());
            };
            let exact = volume(&model, &wedge.solid);
            let mesh =
                ogeom_mesh::triangulate(&model, &wedge.solid, ogeom_mesh::Deflection::default(), T)
                    .unwrap();
            assert!(mesh.is_closed());
            assert!(
                (mesh.volume() - exact).abs() < exact * 0.25,
                "mesh {} against {exact}",
                mesh.volume()
            );
        }
    }

    /// The wedges a marched band takes away or fills walk each edge between
    /// two of their faces once each way and pass the checker, as do the
    /// blended solids, and each holds
    /// what the blend adds to or removes from the solid: a closed band round
    /// a branch's foot, a band running out through a block's walls at both
    /// ends, and a band pinched where two equal drums turn tangent.
    #[test]
    #[ignore = "heavy"]
    fn marched_wedges_walk_each_edge_once_each_way() {
        use ogeom_geom::Curve;
        use ogeom_math::{Direction, Vector};
        let mut model = Model::new();
        let across = |x: f64, y: f64, z: f64| {
            Frame::new(Point::new(x, y, z), Direction::X, Direction::Y, T).unwrap()
        };
        let main = ogeom_algo::make_cylinder(&mut model, across(-20.0, 0.0, 0.0), 10.0, 40.0, T)
            .unwrap()
            .shape;
        let branch = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
            .unwrap()
            .shape;
        let branched = ogeom_bool::fuse(&mut model, &main, &branch, T)
            .unwrap()
            .shape;
        let foot = edge_on(&model, &branched, |c| {
            use ogeom_geom::Curve3d as _;
            let (lo, hi) = c.domain();
            matches!(c, Curve::BSpline(_))
                && c.point_at(lo, T)
                    .and_then(|p| c.point_at(hi, T).map(|q| p.distance(q)))
                    .is_ok_and(|d| d <= 1e-6)
        });

        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let tilt = 0.35_f64;
        let axis = Direction::new(Vector::new(0.0, tilt.cos(), -tilt.sin()), T).unwrap();
        let drum = ogeom_algo::make_cylinder(
            &mut model,
            Frame::new(Point::new(10.0, -5.0, 11.5), axis, Direction::X, T).unwrap(),
            4.0,
            30.0,
            T,
        )
        .unwrap()
        .shape;
        let grooved = ogeom_bool::cut(&mut model, &block, &drum, T).unwrap().shape;
        let groove = edge_at(&model, &grooved, Point::new(10.0, 14.84, 0.0));

        let upright = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
            .unwrap()
            .shape;
        let lying = ogeom_algo::make_cylinder(&mut model, across(-20.0, 0.0, 10.0), 5.0, 40.0, T)
            .unwrap()
            .shape;
        let crossed = ogeom_bool::fuse(&mut model, &upright, &lying, T)
            .unwrap()
            .shape;
        let pinch = edge_on(&model, &crossed, |c| matches!(c, Curve::Ellipse(_)));

        let cases = [
            ("closed band", &branched, &foot, 1.5),
            ("run-out band", &grooved, &groove, 1.0),
            ("pinched band", &crossed, &pinch, 1.0),
        ];
        for (name, solid, edge, radius) in cases {
            let (made, wedges) =
                super::collecting_wedges(|| crate::fillet_edge(&mut model, solid, edge, radius, T));
            made.unwrap();
            let [wedge] = &wedges[..] else {
                panic!("{name}: one wedge, not {}", wedges.len());
            };
            assert_eq!(walked_one_way(&model, &wedge.solid), 0, "{name}");
            // Each edge stands within its own tolerance of every pcurve it
            // carries, measured along its whole length.
            for e in ogeom_topo::explore_unique(&model, &wedge.solid, ShapeType::Edge).unwrap() {
                let stated = model
                    .node(&e)
                    .unwrap()
                    .data()
                    .as_edge()
                    .unwrap()
                    .tolerance
                    .get();
                if let Some(gap) = ogeom_algo::edge_pcurve_gap(&model, &e, T).unwrap() {
                    assert!(
                        gap <= stated,
                        "{name}: a pcurve stands {gap:e} off an edge stating {stated:e}"
                    );
                }
            }
            // The checker asks the walk too, and that every vertex is held
            // as loosely as the edges meeting it.
            let found = ogeom_algo::check(&model, &wedge.solid, T).unwrap();
            assert!(found.is_usable(), "{name}: {found}");
            // What the blend adds is the wedge; what it removes is the part
            // of the wedge inside the solid, the band running on past the
            // walls where it runs out.
            let blended = crate::fillet_edge(&mut model, solid, edge, radius, T)
                .unwrap()
                .shape;
            let found = ogeom_algo::check(&model, &blended, T).unwrap();
            assert!(found.is_usable(), "{name}: the blended solid: {found}");
            let moved = (volume(&model, &blended) - volume(&model, solid)).abs();
            let held = if wedge.additive {
                volume(&model, &wedge.solid)
            } else {
                let inside = ogeom_bool::common(&mut model, &wedge.solid, solid, T)
                    .unwrap()
                    .shape;
                volume(&model, &inside)
            };
            assert!(
                (held - moved).abs() < moved * 1e-3,
                "{name}: the wedge holds {held}, the blend moves {moved}"
            );
        }
    }
}
