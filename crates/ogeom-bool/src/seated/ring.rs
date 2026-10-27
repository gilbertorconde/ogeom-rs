//! A cut or a fill by revolved blend wedges sitting on a rim, made by
//! trimming the two faces each sits on rather than by the general fuse.
//!
//! A rim's blend wedge is three faces, each a whole turn about one axis: a
//! band of a cylinder (the wall seat), an annulus of a plane square to the
//! axis (the cap seat), and the quarter of a torus's tube between them (the
//! blend). Through the axis, the wedge is the square between the tube's
//! centre, the wall and the cap, less the quarter disc of the tube. Where
//! the wall seat lies on a cylindrical face of the solid and the cap seat
//! on a planar one, each facing the same way, and nothing else of the solid
//! comes into the ring that square sweeps, the cut is known without
//! intersecting anything: the cap face loses the annulus, the wall face the
//! band, and the blend closes what is left along the two contact circles.
//! A fill, which rounds a concave rim (a boss's base, a blind hole's
//! floor), is the same trim with the two faces facing against the seats.
//!
//! Each condition is checked, and where one fails the cut is not made here:
//!
//! - The wedge is a torus face, a cylinder face and a plane face at no
//!   placement, about one axis: the cylinder's radius is the tube's distance
//!   from the axis and the tube's radius apart, the plane stands the tube's
//!   radius from the tube's level, square to the axis.
//! - Each face is one loop but the annulus, which is two. The rim (cap and
//!   wall), the wall's contact circle at the tube's level (wall and blend)
//!   and the cap's at the tube's distance (cap and blend) are each one
//!   closed circle on the axis. The wall's seam is a line from its contact
//!   to the rim, and the blend's a meridian of the tube, a quarter turn
//!   about the tube's centre between the two contacts, through the square.
//! - The cap seat lies on a planar face of the solid facing the same way (the
//!   other way for a fill), one of whose loops is the rim: every edge on the
//!   rim's circle, a whole turn. The one other face walking those edges lies
//!   on the wall seat's cylinder, facing as the cap face does against its
//!   seat, and walks the rim between the two walks of one seam line, which
//!   runs on past the tube's level.
//! - Every other edge of the wall's face stands past the tube's level, on the
//!   side away from the cap, and every other edge of the cap's face past the
//!   cap's contact circle, on the side away from the wall, each by more than
//!   the boolean's weld reach and the tolerances on both sides. The two seats
//!   are then the parts of their faces between the rim and the contacts.
//! - Every other face of the solid stands, by its conservative box and within
//!   the tolerance it states, out of the ring the square sweeps: past the
//!   cap's plane or the tube's level, or past the wall's or the contact's
//!   radius.
//! - The wall's seam meets the wall's contact either at the meridian's end,
//!   within both tolerances, or further from it than the boolean welds.
//! - Lumps of a compound tool stand apart, box from box. For a fill, the
//!   general fuse sets none of the tool's faces aside.
//! - The faces made, with the solid's other faces, walk every edge twice,
//!   once each way.
//!
//! The ring's inside then meets no face of the solid: the two the seats lie
//! on bound it, on the cap's plane and the wall's cylinder, and every other
//! stands out of it. It is one connected region, beside each seat on the
//! solid's material side for a cut, so it lies inside the solid, and the
//! wedge with it. The solid less the wedge is bounded by the solid's faces
//! with the seats taken out, and by the blend turned to face out of the
//! material. For a fill each seat has the solid's material on its far side,
//! so the ring lies outside the solid, the wedge with it, and the union is
//! bounded by the solid's faces with the seats taken out and by the blend
//! as the wedge faces it.
//!
//! The result is put together as the general fuse puts its own together:
//! the faces it would set aside kept, every other face made again, a closed
//! edge holding one vertex split at its middle as the fuse splits it, and
//! each contact circle split where the seam and the meridian meet it, at its
//! middle where they are one.

use super::{Copies, Solid, Walked, curve_of, lumps_apart, read_solid, walk_face};
use ogeom_algo::{Built, History, make_wire};
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerance, Tolerances};
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{
    EdgeData, EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId,
    TShapeId, VertexData,
};

/// Whether every lump of a tool holds a torus face: a tool for the ring
/// trim, not the straight one.
pub(super) fn is_ringed(model: &Model, b: &Shape) -> OgeomResult<bool> {
    let lumps = match model.kind_of(b)? {
        ShapeType::Solid => vec![b.clone()],
        ShapeType::Compound => model.children_of(b)?,
        _ => return Ok(false),
    };
    if lumps.is_empty() {
        return Ok(false);
    }
    for lump in &lumps {
        let mut torus = false;
        for face in ogeom_topo::explore_unique(model, lump, ShapeType::Face)? {
            torus |= matches!(
                surface_of(model, &face),
                Some((_, SurfaceGeometry::Torus(_)))
            );
        }
        if !torus {
            return Ok(false);
        }
    }
    Ok(true)
}

/// A rim's axis: a point of it at the tube's level and its unit direction
/// toward the cap.
#[derive(Clone, Copy)]
struct Axis {
    origin: Point,
    up: Vector,
}

impl Axis {
    /// How far along the axis `p` stands from the tube's level.
    fn along(&self, p: Point) -> f64 {
        (p - self.origin).dot(self.up)
    }

    /// How far from the axis `p` stands.
    fn off(&self, p: Point) -> f64 {
        let d = p - self.origin;
        (d - self.up * d.dot(self.up)).magnitude()
    }

    /// The least and greatest distance from the axis over a box: the
    /// greatest at a corner, the least that of the box's shadow across the
    /// axis, the hull of its corners' shadows, from the axis's own.
    fn box_off(&self, corners: &[Point]) -> (f64, f64) {
        let seed = if self.up.x.abs() < 0.9 {
            Vector::new(1.0, 0.0, 0.0)
        } else {
            Vector::new(0.0, 1.0, 0.0)
        };
        let u = self.up.cross(seed);
        let u = u / u.magnitude();
        let v = self.up.cross(u);
        let shadow: Vec<(f64, f64)> = corners
            .iter()
            .map(|c| {
                let d = *c - self.origin;
                (d.dot(u), d.dot(v))
            })
            .collect();
        let most = shadow.iter().map(|(x, y)| x.hypot(*y)).fold(0.0, f64::max);
        (hull_distance(&shadow), most)
    }

    /// The least and greatest of [`Axis::along`] over a box's corners.
    fn box_along(&self, corners: &[Point]) -> (f64, f64) {
        corners
            .iter()
            .map(|c| self.along(*c))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), a| {
                (lo.min(a), hi.max(a))
            })
    }
}

/// How far along the axis and how far from it an edge of the solid reaches:
/// the least and greatest of each. A segment and a circle square to the axis
/// are read exactly, a circle over its whole turn, and any other curve by
/// its conservative box.
fn edge_reach(
    axis: &Axis,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<((f64, f64), (f64, f64))> {
    match curve {
        Curve::Line(_) => {
            let (s, e) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
            let (a0, a1) = (axis.along(s), axis.along(e));
            // The nearest the segment comes to the axis line.
            let d = e - s;
            let across = d - axis.up * d.dot(axis.up);
            let from = s - axis.origin;
            let from = from - axis.up * from.dot(axis.up);
            let length2 = across.dot(across);
            let t = if length2 > 0.0 {
                (-from.dot(across) / length2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let least = (from + across * t).magnitude();
            Ok((
                (a0.min(a1), a0.max(a1)),
                (least, axis.off(s).max(axis.off(e))),
            ))
        }
        Curve::Circle(c) if parallel(c.circle().frame().z().vector(), axis.up, tol) => {
            let circle = c.circle();
            let (a, d, r) = (
                axis.along(circle.centre()),
                axis.off(circle.centre()),
                circle.radius(),
            );
            Ok(((a, a), ((r - d).abs(), r + d)))
        }
        _ => {
            let corners = ogeom_algo::measure::curve_bounds_over(curve, range, tol)?.corners();
            Ok((axis.box_along(&corners), axis.box_off(&corners)))
        }
    }
}

/// How far along the axis and how far from it a face of the solid reaches,
/// the least and greatest of each: by its conservative box, narrowed to the
/// plane's own height for a plane square to the axis, and to the radii a
/// cylinder or a torus about the axis reaches over its whole turn.
fn face_reach(
    model: &Model,
    axis: &Axis,
    face: &Shape,
    bounds: &ogeom_math::Aabb,
    reach: f64,
    tol: Tolerances,
) -> ((f64, f64), (f64, f64)) {
    let corners = bounds.corners();
    let (mut along, mut off) = (axis.box_along(&corners), axis.box_off(&corners));
    let narrow = |range: &mut (f64, f64), lo: f64, hi: f64| {
        if lo <= range.1 && hi >= range.0 {
            *range = (range.0.max(lo), range.1.min(hi));
        }
    };
    match surface_of(model, face) {
        Some((_, SurfaceGeometry::Plane(p))) => {
            let plane = p.plane();
            if parallel(plane.normal().vector(), axis.up, tol) {
                let a = axis.along(plane.frame().origin());
                narrow(&mut along, a, a);
            }
        }
        Some((_, SurfaceGeometry::Cylinder(c))) => {
            let cylinder = c.cylinder();
            if parallel(cylinder.frame().z().vector(), axis.up, tol)
                && axis.off(cylinder.frame().origin()) <= reach
            {
                narrow(&mut off, cylinder.radius(), cylinder.radius());
            }
        }
        Some((_, SurfaceGeometry::Torus(t))) => {
            let torus = t.torus();
            if parallel(torus.frame().z().vector(), axis.up, tol)
                && axis.off(torus.frame().origin()) <= reach
            {
                let (major, minor) = (torus.major_radius(), torus.minor_radius());
                let a = axis.along(torus.frame().origin());
                narrow(&mut off, (major - minor).max(0.0), major + minor);
                narrow(&mut along, a - minor, a + minor);
            }
        }
        _ => {}
    }
    (along, off)
}

/// The distance from the origin to the convex hull of some points in a
/// plane: none where the origin lies inside it.
fn hull_distance(points: &[(f64, f64)]) -> f64 {
    let mut pts = points.to_vec();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    pts.dedup();
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    // The hull counterclockwise, by the monotone chain.
    let mut hull: Vec<(f64, f64)> = Vec::with_capacity(pts.len() * 2);
    for pass in 0..2 {
        let start = hull.len();
        let ordered: Box<dyn Iterator<Item = &(f64, f64)>> = if pass == 0 {
            Box::new(pts.iter())
        } else {
            Box::new(pts.iter().rev())
        };
        for &p in ordered {
            while hull.len() >= start + 2
                && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
    }
    let segment = |a: (f64, f64), b: (f64, f64)| {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let length2 = dx * dx + dy * dy;
        let t = if length2 > 0.0 {
            (-(a.0 * dx + a.1 * dy) / length2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (a.0 + dx * t).hypot(a.1 + dy * t)
    };
    match hull.len() {
        0 => pts.first().map_or(0.0, |p| p.0.hypot(p.1)),
        1 => hull[0].0.hypot(hull[0].1),
        n => {
            let inside =
                n >= 3 && (0..n).all(|i| cross(hull[i], hull[(i + 1) % n], (0.0, 0.0)) >= 0.0);
            if inside {
                return 0.0;
            }
            (0..n)
                .map(|i| segment(hull[i], hull[(i + 1) % n]))
                .fold(f64::INFINITY, f64::min)
        }
    }
}

/// A rim wedge as the trim reads it.
struct Ring {
    axis: Axis,
    /// The tube's distance from the axis, where the cap's contact runs.
    tube: f64,
    /// The wall's radius.
    wall_radius: f64,
    /// The tube's radius, the cap's height over the tube's level.
    height: f64,
    /// The cap seat.
    cap: Shape,
    /// The wall seat.
    wall: Shape,
    /// `+1` where the wall seat faces away from the axis, `-1` toward it.
    wall_sense: f64,
    /// The blend.
    blend: Shape,
    /// The wall's and the cap's contacts.
    wall_ring: Shape,
    cap_ring: Shape,
    /// The blend's meridian, its seam.
    meridian: Shape,
    /// The widest tolerance the wedge states.
    tolerance: f64,
}

/// A face's surface id and geometry, where the face stands at no placement.
fn surface_of(model: &Model, face: &Shape) -> Option<(SurfaceId, SurfaceGeometry)> {
    let NodeData::Face(data) = model.node(face)?.data() else {
        return None;
    };
    if !face.location().is_identity() || !data.location.is_identity() {
        return None;
    }
    Some((
        data.surface,
        model.geometry().surface(data.surface)?.clone(),
    ))
}

/// The sign a face's orientation and its frame's handedness give the
/// surface's normal: a chart's natural normal is `du x dv`, which a
/// left-handed frame turns against its `z`.
fn facing(face: &Shape, frame: &ogeom_math::Frame) -> f64 {
    let turned = face.orientation() == Orientation::Reversed;
    let mirrored = frame.handedness() == ogeom_math::Handedness::Left;
    if turned == mirrored { 1.0 } else { -1.0 }
}

/// A planar face's outward unit normal and a point of its plane, where it
/// stands at no placement.
fn plane_out(model: &Model, face: &Shape) -> Option<(Point, Vector)> {
    let Some((_, SurfaceGeometry::Plane(p))) = surface_of(model, face) else {
        return None;
    };
    let frame = p.plane().frame();
    Some((frame.origin(), frame.z().vector() * facing(face, &frame)))
}

/// Whether two vectors run within the angular tolerance of one line.
fn parallel(a: Vector, b: Vector, tol: Tolerances) -> bool {
    a.cross(b).magnitude() <= tol.angular() * a.magnitude() * b.magnitude()
}

/// A closed circle on the axis: its radius and how far along the axis it
/// stands, where its curve is a circle whose centre lies on the axis and
/// whose plane is square to it, both within `reach`, and the edge runs its
/// whole turn on one vertex.
fn ring_on(
    model: &Model,
    axis: &Axis,
    edge: &Shape,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    let Some((Curve::Circle(c), range, _)) = curve_of(model, edge) else {
        return Ok(None);
    };
    let Some((from, to)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(None);
    };
    let circle = c.circle();
    if from.node() != to.node()
        || ((range.1 - range.0) - core::f64::consts::TAU).abs() > tol.angular()
        || axis.off(circle.centre()) > reach
        || !parallel(circle.frame().z().vector(), axis.up, tol)
    {
        return Ok(None);
    }
    Ok(Some((circle.radius(), axis.along(circle.centre()))))
}

/// The edge a face's loop walks twice, and the edges it walks once.
fn seam_and_rest(ring: &[Walked]) -> Option<(Shape, Vec<Shape>)> {
    let mut seam: Option<Shape> = None;
    let mut rest = Vec::new();
    for (k, walked) in ring.iter().enumerate() {
        let twice = ring
            .iter()
            .enumerate()
            .any(|(j, w)| j != k && w.edge.node() == walked.edge.node());
        if twice {
            if seam
                .as_ref()
                .is_some_and(|s| s.node() != walked.edge.node())
            {
                return None;
            }
            seam = Some(walked.edge.oriented(Orientation::Forward));
        } else {
            rest.push(walked.edge.oriented(Orientation::Forward));
        }
    }
    Some((seam?, rest))
}

/// The wedge read as a [`Ring`], `None` where it is not one as the module
/// asks.
fn read_ring(model: &Model, wedge: &Shape, tol: Tolerances) -> OgeomResult<Option<Ring>> {
    let shells = model.children_of(wedge)?;
    let [shell] = shells.as_slice() else {
        return Ok(None);
    };
    let faces = model.children_of(shell)?;
    if faces.len() != 3 {
        return Ok(None);
    }
    let mut cap: Option<(Shape, Vec<Vec<Walked>>, Point, Vector)> = None;
    let mut wall: Option<(Shape, Vec<Walked>, ogeom_math::Cylinder)> = None;
    let mut blend: Option<(Shape, Vec<Walked>, ogeom_math::Torus)> = None;
    for face in &faces {
        let Some(mut rings) = walk_face(model, face)? else {
            return Ok(None);
        };
        match surface_of(model, face) {
            Some((_, SurfaceGeometry::Plane(_))) if cap.is_none() && rings.len() == 2 => {
                let Some((origin, outward)) = plane_out(model, face) else {
                    return Ok(None);
                };
                cap = Some((face.clone(), rings, origin, outward));
            }
            Some((_, SurfaceGeometry::Cylinder(c))) if wall.is_none() && rings.len() == 1 => {
                wall = Some((face.clone(), rings.remove(0), c.cylinder()));
            }
            Some((_, SurfaceGeometry::Torus(t))) if blend.is_none() && rings.len() == 1 => {
                blend = Some((face.clone(), rings.remove(0), t.torus()));
            }
            _ => return Ok(None),
        }
    }
    let (Some(cap), Some(wall), Some(blend)) = (cap, wall, blend) else {
        return Ok(None);
    };
    let mut tolerance = tol.confusion();
    for face in &faces {
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in ogeom_topo::explore_unique(model, face, kind)? {
                tolerance = tolerance.max(model.tolerance_of(&sub)?.map_or(0.0, |t| t.get()));
            }
        }
    }
    let reach = tolerance + tol.confusion();
    let torus = blend.2;
    let frame = torus.frame();
    let (tube, radius) = (torus.major_radius(), torus.minor_radius());
    let cylinder = wall.2;
    let wall_radius = cylinder.radius();
    // The cap's side of the tube's level.
    let mut axis = Axis {
        origin: frame.origin(),
        up: frame.z().vector(),
    };
    let height = axis.along(cap.2);
    if height < 0.0 {
        axis.up = -axis.up;
    }
    let height = height.abs();
    if tube <= reach
        || (height - radius).abs() > reach
        || ((wall_radius - tube).abs() - radius).abs() > reach
        || !parallel(cap.3, axis.up, tol)
        || cap.3.dot(axis.up) <= 0.0
        || !parallel(cylinder.frame().z().vector(), axis.up, tol)
        || axis.off(cylinder.frame().origin()) > reach
    {
        return Ok(None);
    }
    // The edges each two faces share.
    let shared = |x: &[Walked], y: &[Walked]| -> Vec<Shape> {
        let mut out: Vec<Shape> = Vec::new();
        for w in x {
            if y.iter().any(|v| v.edge.node() == w.edge.node())
                && !out.iter().any(|e| e.node() == w.edge.node())
            {
                out.push(w.edge.oriented(Orientation::Forward));
            }
        }
        out
    };
    let cap_edges: Vec<Walked> = cap.1.iter().flatten().cloned().collect();
    let (rim, cap_ring, wall_ring) = match (
        shared(&cap_edges, &wall.1).as_slice(),
        shared(&cap_edges, &blend.1).as_slice(),
        shared(&wall.1, &blend.1).as_slice(),
    ) {
        ([rim], [cap_ring], [wall_ring]) => (rim.clone(), cap_ring.clone(), wall_ring.clone()),
        _ => return Ok(None),
    };
    // Each loop: the annulus's two rings, the wall's and the blend's two
    // rings and a seam each.
    if cap.1.iter().any(|ring| ring.len() != 1) {
        return Ok(None);
    }
    let (Some((wall_seam, wall_rings)), Some((meridian, blend_rings))) =
        (seam_and_rest(&wall.1), seam_and_rest(&blend.1))
    else {
        return Ok(None);
    };
    let pair = |rings: &[Shape], a: &Shape, b: &Shape| {
        rings.len() == 2
            && rings.iter().any(|e| e.node() == a.node())
            && rings.iter().any(|e| e.node() == b.node())
    };
    if wall.1.len() != 4
        || blend.1.len() != 4
        || !pair(&wall_rings, &rim, &wall_ring)
        || !pair(&blend_rings, &cap_ring, &wall_ring)
    {
        return Ok(None);
    }
    // The rings where the square's corners put them.
    let at = |edge: &Shape, r: f64, a: f64| -> OgeomResult<bool> {
        Ok(ring_on(model, &axis, edge, reach, tol)?
            .is_some_and(|(own, along)| (own - r).abs() <= reach && (along - a).abs() <= reach))
    };
    if !at(&rim, wall_radius, height)?
        || !at(&wall_ring, wall_radius, 0.0)?
        || !at(&cap_ring, tube, height)?
    {
        return Ok(None);
    }
    // The wall's seam is a line, and the meridian a quarter turn through
    // the square, its middle standing inside it.
    let Some((Curve::Line(_), _, _)) = curve_of(model, &wall_seam) else {
        return Ok(None);
    };
    let Some((curve @ Curve::Circle(_), range, _)) = curve_of(model, &meridian) else {
        return Ok(None);
    };
    let sweep = (range.1 - range.0).abs();
    let middle = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    let (lo, hi) = (tube.min(wall_radius), tube.max(wall_radius));
    let (off, along) = (axis.off(middle), axis.along(middle));
    if (sweep - core::f64::consts::FRAC_PI_2).abs() > tol.angular()
        || off <= lo + reach
        || off >= hi - reach
        || along <= reach
        || along >= height - reach
    {
        return Ok(None);
    }
    // The wall seat faces away from the tube: out from the axis where the
    // wall stands beyond the tube.
    let wall_sense = if wall_radius > tube { 1.0 } else { -1.0 };
    if facing(&wall.0, &cylinder.frame()) != wall_sense {
        return Ok(None);
    }
    Ok(Some(Ring {
        axis,
        tube,
        wall_radius,
        height,
        cap: cap.0,
        wall: wall.0,
        wall_sense,
        blend: blend.0,
        wall_ring,
        cap_ring,
        meridian,
        tolerance,
    }))
}

/// A ring placed on the solid: the faces its seats lie on, the cap face's
/// loop that is the rim, the rim's edges and the wall face's seam.
struct Placed {
    ring: Ring,
    cap_host: usize,
    cap_loop: usize,
    wall_host: usize,
    rims: Vec<TShapeId>,
    seam: Shape,
}

/// Find the faces a ring's seats lie on and check the module's conditions
/// on them; `None` where one fails. The faces face the seats' way for a cut
/// and against it for a fill.
fn place(
    model: &Model,
    solid: &mut Solid,
    ring: Ring,
    fill: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Placed>> {
    let reach = ring.tolerance + tol.confusion();
    let axis = ring.axis;
    let host_sense = if fill { -1.0 } else { 1.0 };
    // The cap face: a plane facing the way the cap seat faces, or against it
    // for a fill, through the rim, one of whose loops runs on the rim's
    // circle alone.
    let mut found: Option<(usize, usize)> = None;
    for (i, face) in solid.faces.iter().enumerate() {
        let Some((origin, outward)) = plane_out(model, face) else {
            continue;
        };
        if !parallel(outward, axis.up, tol)
            || outward.dot(axis.up) * host_sense <= 0.0
            || (axis.along(origin) - ring.height).abs() > reach
        {
            continue;
        }
        for (l, walked) in solid.loops[i].iter().enumerate() {
            let mut on_rim = true;
            let mut turn = 0.0;
            for w in walked {
                let Some((Curve::Circle(circle), range, own)) = curve_of(model, &w.edge) else {
                    on_rim = false;
                    break;
                };
                let circle = circle.circle();
                let near = own + reach;
                if axis.off(circle.centre()) > near
                    || (axis.along(circle.centre()) - ring.height).abs() > near
                    || (circle.radius() - ring.wall_radius).abs() > near
                    || !parallel(circle.frame().z().vector(), axis.up, tol)
                {
                    on_rim = false;
                    break;
                }
                turn += (range.1 - range.0).abs();
            }
            if on_rim && (turn - core::f64::consts::TAU).abs() <= tol.angular() {
                if found.is_some() {
                    return Ok(None);
                }
                found = Some((i, l));
            }
        }
    }
    let Some((cap_host, cap_loop)) = found else {
        return Ok(None);
    };
    let rims: Vec<TShapeId> = solid.loops[cap_host][cap_loop]
        .iter()
        .map(|w| w.edge.node())
        .collect();
    // The wall face: the other face walking every rim edge.
    let mut walls: FastSet<usize> = FastSet::default();
    for rim in &rims {
        let Some(by) = solid.walkers.get(rim) else {
            return Ok(None);
        };
        if by.len() != 2 || !by.contains(&cap_host) {
            return Ok(None);
        }
        walls.extend(by.iter().copied().filter(|&f| f != cap_host));
    }
    let wall_host = match walls.iter().copied().collect::<Vec<_>>().as_slice() {
        [one] => *one,
        _ => return Ok(None),
    };
    let wall_face = &solid.faces[wall_host];
    let Some((_, SurfaceGeometry::Cylinder(c))) = surface_of(model, wall_face) else {
        return Ok(None);
    };
    let cylinder = c.cylinder();
    if facing(wall_face, &cylinder.frame()) != ring.wall_sense * host_sense
        || (cylinder.radius() - ring.wall_radius).abs() > reach
        || !parallel(cylinder.frame().z().vector(), axis.up, tol)
        || axis.off(cylinder.frame().origin()) > reach
    {
        return Ok(None);
    }
    // In the wall face's loop the rim is walked whole, between the two walks
    // of one seam line.
    let mut seam: Option<Shape> = None;
    let mut in_loop = None;
    for (l, walked) in solid.loops[wall_host].iter().enumerate() {
        let n = walked.len();
        let on = |k: usize| rims.contains(&walked[k % n].edge.node());
        if !(0..n).any(on) {
            continue;
        }
        if in_loop.is_some() {
            return Ok(None);
        }
        in_loop = Some(l);
        let Some(start) = (0..n).find(|&k| on(k) && !on(k + n - 1)) else {
            return Ok(None);
        };
        let len = (0..n).take_while(|&k| on(start + k)).count();
        if len != rims.len() || (len..n).any(|k| on(start + k)) {
            return Ok(None);
        }
        let (before, after) = (&walked[(start + n - 1) % n], &walked[(start + len) % n]);
        if before.edge.node() != after.edge.node() || before.to != walked[start].from {
            return Ok(None);
        }
        seam = Some(before.edge.oriented(Orientation::Forward));
    }
    let Some(seam) = seam else {
        return Ok(None);
    };
    // The seam a line along the axis, running from the rim past the tube's
    // level.
    let Some((line @ Curve::Line(_), range, seam_tolerance)) = curve_of(model, &seam) else {
        return Ok(None);
    };
    if !parallel(line.d1_at(range.0, tol)?, axis.up, tol) {
        return Ok(None);
    }
    let near = seam_tolerance + reach + crate::sub_shapes::reach(tol);
    let Some((first, last)) = ogeom_algo::edge_vertices(model, &seam)? else {
        return Ok(None);
    };
    let ends: Vec<f64> = [first, last]
        .iter()
        .filter_map(|v| super::vertex_at(model, v).map(|(p, _)| axis.along(p)))
        .collect();
    let [e0, e1] = ends.as_slice() else {
        return Ok(None);
    };
    if e0.min(*e1) >= -near || (e0.max(*e1) - ring.height).abs() > near {
        return Ok(None);
    }
    // Every other edge of the wall face stands past the tube's level, away
    // from the cap.
    for walked in solid.loops[wall_host].iter().flatten() {
        let id = walked.edge.node();
        if rims.contains(&id) || id == seam.node() {
            continue;
        }
        let Some((curve, range, own)) = curve_of(model, &walked.edge) else {
            return Ok(None);
        };
        let ((_, high), _) = edge_reach(&axis, &curve, range, tol)?;
        if high >= -(own + near) {
            return Ok(None);
        }
    }
    // Every other edge of the cap face stands past the cap's contact, away
    // from the wall.
    for (l, walked) in solid.loops[cap_host].iter().enumerate() {
        if l == cap_loop {
            continue;
        }
        for w in walked {
            let Some((curve, range, own)) = curve_of(model, &w.edge) else {
                return Ok(None);
            };
            let (_, (least, most)) = edge_reach(&axis, &curve, range, tol)?;
            let clear = if ring.wall_radius > ring.tube {
                most <= ring.tube - (own + near)
            } else {
                least >= ring.tube + (own + near)
            };
            if !clear {
                return Ok(None);
            }
        }
    }
    // Every other face of the solid stands out of the ring.
    let (lo, hi) = (
        ring.tube.min(ring.wall_radius),
        ring.tube.max(ring.wall_radius),
    );
    for i in 0..solid.faces.len() {
        if i == cap_host || i == wall_host {
            continue;
        }
        let (bounds, stated) = solid.face_box(model, i, tol)?;
        let margin = stated + tol.confusion();
        let ((low, high), (least, most)) =
            face_reach(model, &axis, &solid.faces[i], &bounds, reach, tol);
        if low >= ring.height - margin
            || high <= margin
            || most <= lo + margin
            || least >= hi - margin
        {
            continue;
        }
        return Ok(None);
    }
    Ok(Some(Placed {
        ring,
        cap_host,
        cap_loop,
        wall_host,
        rims,
        seam,
    }))
}

/// A closed curve's parameter at a point standing on it, within a whole
/// turn after `start`: read off the curve at two parameters a quarter turn
/// apart, whichever way it runs.
fn turn_at(curve: &Curve, start: f64, p: Point, tol: Tolerances) -> OgeomResult<Option<f64>> {
    let Curve::Circle(c) = curve else {
        return Ok(None);
    };
    let centre = c.circle().centre();
    let x = curve.point_at(0.0, tol)? - centre;
    let y = curve.point_at(core::f64::consts::FRAC_PI_2, tol)? - centre;
    let d = p - centre;
    let t = d.dot(y).atan2(d.dot(x));
    Ok(Some(
        t + ((start - t) / core::f64::consts::TAU).ceil() * core::f64::consts::TAU,
    ))
}

/// The edge over part of its range: every description `keep` holds, each
/// running over the curve's own range, over `over` instead; `None` where
/// one runs over another range.
fn piece(
    data: &EdgeData,
    range: (f64, f64),
    over: (f64, f64),
    keep: &dyn Fn(&EdgeRepr) -> bool,
) -> Option<EdgeData> {
    let mut data = data.clone();
    data.representations.retain(|r| keep(r));
    for r in &mut data.representations {
        let (EdgeRepr::Curve3d { range: own, .. }
        | EdgeRepr::PCurve { range: own, .. }
        | EdgeRepr::Seam { range: own, .. }) = r
        else {
            return None;
        };
        if *own != range {
            return None;
        }
        *own = over;
    }
    Some(data)
}

/// Whether a description is the 3D curve or an image on `surface`.
fn on_surface(r: &EdgeRepr, surface: SurfaceId) -> bool {
    match r {
        EdgeRepr::Curve3d { .. } => true,
        EdgeRepr::PCurve {
            surface: s,
            location,
            ..
        }
        | EdgeRepr::Seam {
            surface: s,
            location,
            ..
        } => *s == surface && location.is_identity(),
        _ => false,
    }
}

/// A closed edge's pieces between the parameters `cuts` (ascending, its
/// range's start first and end last), the vertex at each inner cut given:
/// each piece the edge over its part of the range, holding what `keep`
/// holds.
fn cut_ring(
    model: &mut Model,
    edge: &Shape,
    cuts: &[f64],
    at: &[Shape],
    keep: &dyn Fn(&EdgeRepr) -> bool,
) -> OgeomResult<Option<Vec<Shape>>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()).cloned() else {
        return Ok(None);
    };
    let Some(EdgeRepr::Curve3d { range, .. }) = data.curve3d().cloned() else {
        return Ok(None);
    };
    let mut out = Vec::with_capacity(cuts.len() - 1);
    for k in 0..cuts.len() - 1 {
        let Some(part) = piece(&data, range, (cuts[k], cuts[k + 1]), keep) else {
            return Ok(None);
        };
        let from = at[k % at.len()].clone();
        let to = at[(k + 1) % at.len()].clone();
        out.push(model.add_edge(part, &[from, to])?);
    }
    Ok(Some(out))
}

/// A new vertex standing on `curve` at `t`, as uncertain as `tolerance`.
fn vertex_on(
    model: &mut Model,
    curve: &Curve,
    t: f64,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    Ok(model.add_vertex(VertexData {
        point: curve.point_at(t, tol)?,
        tolerance: Tolerance::new(tolerance.max(tol.confusion()))?,
    }))
}

/// Which way a walked circle turns about the axis: `+1` counterclockwise
/// seen from the cap's side.
fn turning(model: &Model, axis: &Axis, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<f64>> {
    let Some((curve, range, _)) = curve_of(model, &edge.oriented(Orientation::Forward)) else {
        return Ok(None);
    };
    let t = f64::midpoint(range.0, range.1);
    let p = curve.point_at(t, tol)?;
    let mut d = curve.d1_at(t, tol)?;
    if edge.orientation() == Orientation::Reversed {
        d = -d;
    }
    let out = p - axis.origin;
    let out = out - axis.up * out.dot(axis.up);
    let spin = axis.up.dot(out.cross(d));
    Ok((spin.abs() > tol.confusion() * d.magnitude()).then_some(spin.signum()))
}

/// A ring's contact circle in pieces: the edges made, in the order the
/// circle's parameter runs from the ring's own vertex, and in that order
/// from the vertex the face of the solid walks it from.
struct Cut {
    pieces: Vec<Shape>,
    around: Vec<Shape>,
}

impl Cut {
    /// The pieces walked once round from the vertex the face of the solid
    /// walks them from, turning `sense` about the axis, where the curve's
    /// parameter turns `own`.
    fn walked(&self, sense: f64, own: f64) -> Vec<Shape> {
        if sense == own {
            self.around.clone()
        } else {
            self.around.iter().rev().map(Shape::reversed).collect()
        }
    }

    /// The pieces in place of the whole ring walked as `walked` walks it.
    fn instead_of(&self, walked: &Shape) -> Vec<Shape> {
        if walked.orientation() == Orientation::Reversed {
            self.pieces.iter().rev().map(Shape::reversed).collect()
        } else {
            self.pieces.clone()
        }
    }
}

/// What the trim makes of one ring: the wall's and the cap's contacts in
/// pieces, the meridian, and which way each contact's parameter turns.
struct Made {
    wall: Cut,
    wall_turn: f64,
    cap: Cut,
    cap_turn: f64,
    meridian: Shape,
}

/// `a` with the rim wedges of `b` cut from it, or fused on to it where
/// `fill` is set, by trimming the faces they sit on, keeping as they stand
/// the faces the general fuse sets aside when it applies `whole` (`b`
/// itself, or a tool `b`'s lumps are part of); `None` where a condition of
/// the module's fails.
pub(super) fn ring_cut(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fill: bool,
    whole: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(mut solid) = read_solid(model, a)? else {
        return Ok(None);
    };
    let Some(wedges) = lumps_apart(model, b, tol)? else {
        return Ok(None);
    };
    let mut placed = Vec::with_capacity(wedges.len());
    for wedge in &wedges {
        let Some(ring) = read_ring(model, wedge, tol)? else {
            return Ok(None);
        };
        let Some(on) = place(model, &mut solid, ring, fill, tol)? else {
            return Ok(None);
        };
        placed.push(on);
    }
    // No face hosts two seats on one loop, and no two rings share a rim.
    for (i, p) in placed.iter().enumerate() {
        for q in &placed[i + 1..] {
            if (p.cap_host, p.cap_loop) == (q.cap_host, q.cap_loop)
                || p.rims.iter().any(|r| q.rims.contains(r))
            {
                return Ok(None);
            }
        }
    }
    // Each wall's seam split where the wall's contact crosses it.
    let mut splits = Vec::with_capacity(placed.len());
    for p in &placed {
        let Some((curve, range, _)) = curve_of(model, &p.seam) else {
            return Ok(None);
        };
        let (s, e) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
        let (a0, a1) = (p.ring.axis.along(s), p.ring.axis.along(e));
        let at = s + (e - s) * (a0 / (a0 - a1));
        let Some(vertex) = super::split_at(model, &mut solid, at, p.ring.tolerance, tol)? else {
            return Ok(None);
        };
        splits.push(vertex);
    }
    // The faces the general fuse would set aside applying the whole tool,
    // which it passes through as they stand: none is one a seat lies on,
    // and for a fill, which keeps the tool's own, none is a face of `b`.
    let mut aside = vec![false; solid.faces.len()];
    if let Ok(Some(sides)) = crate::aside::set_aside(model, a, whole, tol) {
        if fill {
            let own = ogeom_topo::explore_unique(model, b, ShapeType::Face)?;
            if sides[1]
                .faces
                .iter()
                .any(|(face, _)| own.iter().any(|f| f.is_same(face)))
            {
                return Ok(None);
            }
        }
        for (face, _) in &sides[0].faces {
            if let Some(i) = solid.faces.iter().position(|f| f.is_same(face)) {
                aside[i] = true;
            }
        }
    }
    if placed
        .iter()
        .any(|p| aside[p.cap_host] || aside[p.wall_host])
    {
        return Ok(None);
    }
    build(model, &solid, &placed, &splits, &aside, fill, a, b, tol)
}

/// The loop of the wall face the ring's run lies in, where the run starts
/// in it and how long it is, once the seam is split: the seam's piece up to
/// the rim, the rim, and the piece back.
fn wall_run(solid: &Solid, p: &Placed, split: TShapeId) -> Option<(usize, usize, usize)> {
    for (l, walked) in solid.loops[p.wall_host].iter().enumerate() {
        let n = walked.len();
        let on = |k: usize| p.rims.contains(&walked[k % n].edge.node());
        let Some(start) = (0..n).find(|&k| on(k) && !on(k + n - 1)) else {
            continue;
        };
        let len = p.rims.len();
        let (before, after) = (&walked[(start + n - 1) % n], &walked[(start + len) % n]);
        if before.edge.node() != after.edge.node() || before.from != split || after.to != split {
            return None;
        }
        return Some((l, (start + n - 1) % n, len + 2));
    }
    None
}

/// Make the faces the rings call for and the solid they bound, with the
/// faces' history; `None` where a face made does not close or its chart
/// breaks.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn build(
    model: &mut Model,
    solid: &Solid,
    placed: &[Placed],
    splits: &[Shape],
    aside: &[bool],
    fill: bool,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut copies = Copies {
        kept: FastSet::default(),
        vertices: FastMap::default(),
        edges: FastMap::default(),
    };
    for (i, rings) in solid.loops.iter().enumerate() {
        if !aside[i] {
            continue;
        }
        for walked in rings.iter().flatten() {
            copies.kept.insert(walked.edge.node());
            for end in model.node(&walked.edge).map_or(&[][..], |n| n.children()) {
                copies.kept.insert(end.node());
            }
        }
    }
    // The runs each ring takes from its two faces.
    let mut runs: Vec<(usize, usize, usize)> = Vec::with_capacity(placed.len());
    for (p, split) in placed.iter().zip(splits) {
        let Some(run) = wall_run(solid, p, split.node()) else {
            return Ok(None);
        };
        runs.push(run);
    }
    let dropped: FastSet<TShapeId> = placed
        .iter()
        .zip(&runs)
        .flat_map(|(p, &(l, start, len))| {
            let walked = &solid.loops[p.wall_host][l];
            (0..len)
                .map(move |k| walked[(start + k) % walked.len()].edge.node())
                .chain(p.rims.iter().copied())
        })
        .collect();
    let Some(halves) = halved(model, solid, aside, &dropped, &mut copies, tol)? else {
        return Ok(None);
    };
    let mut made: Vec<Made> = Vec::with_capacity(placed.len());
    for ((p, split), run) in placed.iter().zip(splits).zip(&runs) {
        let Some(m) = cut_contacts(model, p, split, &mut copies, tol)? else {
            return Ok(None);
        };
        if !image_contacts(model, solid, p, &m, *run, tol)? {
            return Ok(None);
        }
        made.push(m);
    }
    // An edge of the solid as the result holds it, walked as `walked` walks
    // it: its halves where it is closed.
    let held =
        |model: &mut Model, copies: &mut Copies, walked: &Walked| -> OgeomResult<Vec<Shape>> {
            let id = walked.edge.node();
            if let Some(parts) = halves.get(&id) {
                return Ok(if walked.edge.orientation() == Orientation::Reversed {
                    parts.iter().rev().map(Shape::reversed).collect()
                } else {
                    parts.clone()
                });
            }
            let edge = copies.edge(model, &walked.edge.oriented(Orientation::Forward))?;
            Ok(vec![edge.oriented(walked.edge.orientation())])
        };
    let mut history = History::new();
    let mut faces: Vec<Shape> = Vec::with_capacity(solid.faces.len() + placed.len());
    let mut checked: Vec<Shape> = Vec::new();
    for (i, face) in solid.faces.iter().enumerate() {
        if aside[i] {
            history.copy(face, face.clone());
            faces.push(face.clone());
            continue;
        }
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            return Ok(None);
        };
        let mut data = (**data).clone();
        data.triangulation = None;
        let mut changed = solid.split.contains(&i);
        let mut wires = Vec::with_capacity(solid.loops[i].len());
        for (l, walked) in solid.loops[i].iter().enumerate() {
            let mut edges: Vec<Shape> = Vec::with_capacity(walked.len() + 2);
            let cap_of = placed
                .iter()
                .zip(&made)
                .find(|(p, _)| p.cap_host == i && p.cap_loop == l);
            if let Some((p, m)) = cap_of {
                let Some(sense) = turning(model, &p.ring.axis, &walked[0].edge, tol)? else {
                    return Ok(None);
                };
                edges.extend(m.cap.walked(sense, m.cap_turn));
                changed = true;
            } else {
                // Each position a ring's run takes, and the ring whose run
                // starts there; walked from a position none takes.
                let n = walked.len();
                let mut taken = vec![false; n];
                let mut starts: FastMap<usize, usize> = FastMap::default();
                for (r, (p, &(at_loop, start, len))) in placed.iter().zip(&runs).enumerate() {
                    if p.wall_host != i || at_loop != l {
                        continue;
                    }
                    for k in 0..len {
                        let q = (start + k) % n;
                        if taken[q] {
                            return Ok(None);
                        }
                        taken[q] = true;
                    }
                    starts.insert(start, r);
                }
                let Some(first) = (0..n).find(|&q| !taken[q]) else {
                    return Ok(None);
                };
                for k in 0..n {
                    let q = (first + k) % n;
                    if let Some(&r) = starts.get(&q) {
                        let rim = &walked[(q + 1) % n].edge;
                        let Some(sense) = turning(model, &placed[r].ring.axis, rim, tol)? else {
                            return Ok(None);
                        };
                        edges.extend(made[r].wall.walked(sense, made[r].wall_turn));
                        changed = true;
                    }
                    if !taken[q] {
                        let w = &walked[q];
                        changed |= halves.contains_key(&w.edge.node());
                        edges.extend(held(model, &mut copies, w)?);
                    }
                }
            }
            if face.orientation() == Orientation::Reversed {
                edges = edges.iter().rev().map(Shape::reversed).collect();
            }
            wires.push(make_wire(model, &edges, tol)?.shape);
        }
        if changed {
            data.natural_restriction = false;
        }
        let image = model.add_face(data, &wires)?.oriented(face.orientation());
        if changed {
            history.modify(face, image.clone());
            checked.push(image.clone());
        } else {
            history.copy(face, image.clone());
        }
        faces.push(image);
    }
    // Each blend, its rings in pieces, facing out of the material: turned
    // for a cut.
    for (p, m) in placed.iter().zip(&made) {
        let Some(face) = blend_face(model, p, m, fill, tol)? else {
            return Ok(None);
        };
        history.modify(&p.ring.blend, face.clone());
        history.delete(&p.ring.cap);
        history.delete(&p.ring.wall);
        checked.push(face.clone());
        faces.push(face);
    }
    if !walked_twice(model, &faces)? {
        return Ok(None);
    }
    for face in &checked {
        if !ogeom_algo::chart_breaks(model, face, &|_| false, tol)?.is_empty() {
            return Ok(None);
        }
    }
    let shell = model.add_shell(&faces)?;
    let result = model.add_solid(&[shell])?;
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Some(Built::new(result, history)))
}

/// Every closed edge of the solid the faces made again keep, cut at its
/// middle as the general fuse cuts one holding a single vertex; `None` where
/// a face set aside keeps it whole.
fn halved(
    model: &mut Model,
    solid: &Solid,
    aside: &[bool],
    dropped: &FastSet<TShapeId>,
    copies: &mut Copies,
    tol: Tolerances,
) -> OgeomResult<Option<FastMap<TShapeId, Vec<Shape>>>> {
    let mut halves: FastMap<TShapeId, Vec<Shape>> = FastMap::default();
    for (i, rings) in solid.loops.iter().enumerate() {
        if aside[i] {
            continue;
        }
        for walked in rings.iter().flatten() {
            let id = walked.edge.node();
            if walked.from != walked.to || dropped.contains(&id) || halves.contains_key(&id) {
                continue;
            }
            if copies.kept.contains(&id) {
                return Ok(None);
            }
            let edge = walked.edge.oriented(Orientation::Forward);
            let Some((curve, range, own)) = curve_of(model, &edge) else {
                return Ok(None);
            };
            let Some((start, _)) = ogeom_algo::edge_vertices(model, &edge)? else {
                return Ok(None);
            };
            let start = copies.vertex(model, &start)?;
            let middle = f64::midpoint(range.0, range.1);
            let mid = vertex_on(model, &curve, middle, own, tol)?;
            let Some(cut) = cut_ring(
                model,
                &edge,
                &[range.0, middle, range.1],
                &[start, mid],
                &|r| {
                    matches!(
                        r,
                        EdgeRepr::Curve3d { .. } | EdgeRepr::PCurve { .. } | EdgeRepr::Seam { .. }
                    )
                },
            )?
            else {
                return Ok(None);
            };
            halves.insert(id, cut);
        }
    }
    Ok(Some(halves))
}

/// A ring's contacts in pieces and its meridian, on the vertices the result
/// holds: the wall's contact cut where the seam and the meridian meet it,
/// and at its middle where they are one, the cap's at its middle.
fn cut_contacts(
    model: &mut Model,
    p: &Placed,
    split: &Shape,
    copies: &mut Copies,
    tol: Tolerances,
) -> OgeomResult<Option<Made>> {
    let ring = &p.ring;
    let Some(NodeData::Face(blend)) = model.node(&ring.blend).map(|n| n.data()) else {
        return Ok(None);
    };
    let blend_surface = blend.surface;
    let keep = move |r: &EdgeRepr| on_surface(r, blend_surface);
    let (Some((wall_curve, wall_range, wall_tol)), Some((cap_curve, cap_range, cap_tol))) = (
        curve_of(model, &ring.wall_ring),
        curve_of(model, &ring.cap_ring),
    ) else {
        return Ok(None);
    };
    let (Some((wall_vertex, _)), Some((cap_vertex, _))) = (
        ogeom_algo::edge_vertices(model, &ring.wall_ring)?,
        ogeom_algo::edge_vertices(model, &ring.cap_ring)?,
    ) else {
        return Ok(None);
    };
    let (Some((meets, meets_tol)), Some((at_split, split_tol))) = (
        super::vertex_at(model, &wall_vertex),
        super::vertex_at(model, split),
    ) else {
        return Ok(None);
    };
    let start = wall_curve.point_at(wall_range.0, tol)?;
    if start.distance(meets) > meets_tol + wall_tol {
        return Ok(None);
    }
    let split_copy = copies.vertex(model, split)?;
    // The wall's contact, cut where the seam and the meridian meet it,
    // and at its middle where they are one.
    let apart = at_split.distance(meets);
    let (wall_cuts, wall_at, meridian_start, rotated) = if apart <= split_tol + meets_tol {
        let middle = wall_range.0 + core::f64::consts::PI;
        let mid = vertex_on(model, &wall_curve, middle, wall_tol, tol)?;
        (
            vec![wall_range.0, middle, wall_range.1],
            vec![split_copy.clone(), mid],
            split_copy.clone(),
            false,
        )
    } else if apart > split_tol + meets_tol + crate::sub_shapes::reach(tol) {
        let Some(t) = turn_at(&wall_curve, wall_range.0, at_split, tol)? else {
            return Ok(None);
        };
        if t <= wall_range.0 || t >= wall_range.1 {
            return Ok(None);
        }
        let Some(data) = model.node(&wall_vertex).and_then(|n| n.data().as_vertex()) else {
            return Ok(None);
        };
        let own = model.add_vertex(data.clone());
        (
            vec![wall_range.0, t, wall_range.1],
            vec![own.clone(), split_copy.clone()],
            own,
            true,
        )
    } else {
        return Ok(None);
    };
    let Some(wall_pieces) = cut_ring(model, &ring.wall_ring, &wall_cuts, &wall_at, &keep)? else {
        return Ok(None);
    };
    let Some(data) = model.node(&cap_vertex).and_then(|n| n.data().as_vertex()) else {
        return Ok(None);
    };
    let cap_own = model.add_vertex(data.clone());
    let middle = cap_range.0 + core::f64::consts::PI;
    let cap_mid = vertex_on(model, &cap_curve, middle, cap_tol, tol)?;
    let Some(cap_pieces) = cut_ring(
        model,
        &ring.cap_ring,
        &[cap_range.0, middle, cap_range.1],
        &[cap_own.clone(), cap_mid],
        &keep,
    )?
    else {
        return Ok(None);
    };
    // The meridian, between the contacts' vertices where it meets them.
    let Some(data) = model
        .node(&ring.meridian)
        .and_then(|n| n.data().as_edge())
        .cloned()
    else {
        return Ok(None);
    };
    let Some(EdgeRepr::Curve3d { range, .. }) = data.curve3d().cloned() else {
        return Ok(None);
    };
    let Some(data) = piece(&data, range, range, &keep) else {
        return Ok(None);
    };
    let Some((first, _)) = ogeom_algo::edge_vertices(model, &ring.meridian)? else {
        return Ok(None);
    };
    let ends = if first.node() == wall_vertex.node() {
        [meridian_start, cap_own]
    } else {
        [cap_own, meridian_start]
    };
    let meridian = model.add_edge(data, &ends)?;
    let (Some(wall_turn), Some(cap_turn)) = (
        turning(model, &ring.axis, &ring.wall_ring, tol)?,
        turning(model, &ring.axis, &ring.cap_ring, tol)?,
    ) else {
        return Ok(None);
    };
    let mut around = wall_pieces.clone();
    if rotated {
        around.rotate_left(1);
    }
    Ok(Some(Made {
        wall: Cut {
            pieces: wall_pieces,
            around,
        },
        wall_turn,
        cap: Cut {
            around: cap_pieces.clone(),
            pieces: cap_pieces,
        },
        cap_turn,
        meridian,
    }))
}

/// The contacts' images on the faces of the solid they cross: the cap's on
/// its plane, the wall's on its cylinder in the turn the rim's own images
/// run over; `false` where one has no image in closed form.
fn image_contacts(
    model: &mut Model,
    solid: &Solid,
    p: &Placed,
    m: &Made,
    (l, start, len): (usize, usize, usize),
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some(NodeData::Face(cap)) = model.node(&solid.faces[p.cap_host]).map(|n| n.data()) else {
        return Ok(false);
    };
    let cap_surface = cap.surface;
    let Some(NodeData::Face(wall)) = model.node(&solid.faces[p.wall_host]).map(|n| n.data()) else {
        return Ok(false);
    };
    let wall_surface = wall.surface;
    let (Some(cap_geometry), Some(wall_geometry)) = (
        model.geometry().surface(cap_surface).cloned(),
        model.geometry().surface(wall_surface).cloned(),
    ) else {
        return Ok(false);
    };
    let walked = &solid.loops[p.wall_host][l];
    let mut span = (f64::INFINITY, f64::NEG_INFINITY);
    for k in 1..len - 1 {
        let edge = &walked[(start + k) % walked.len()].edge;
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            return Ok(false);
        };
        let Some(EdgeRepr::PCurve { curve, range, .. }) = data.pcurve_on(wall_surface) else {
            return Ok(false);
        };
        let Some(image) = model.geometry().pcurve(*curve) else {
            return Ok(false);
        };
        for t in [range.0, f64::midpoint(range.0, range.1), range.1] {
            let u = image.point_at(t, tol)?.x;
            span = (span.0.min(u), span.1.max(u));
        }
    }
    use ogeom_geom::Surface as _;
    let period = wall_geometry
        .is_periodic_u()
        .then(|| wall_geometry.domain().0)
        .map(|(lo, hi)| hi - lo);
    let Some(period) = period else {
        return Ok(false);
    };
    for (edges, surface, geometry, turn_into) in [
        (&m.cap.pieces, cap_surface, &cap_geometry, None),
        (&m.wall.pieces, wall_surface, &wall_geometry, Some(span)),
    ] {
        for edge in edges {
            let Some((curve, range, _)) = curve_of(model, edge) else {
                return Ok(false);
            };
            let Some(mut image) = ogeom_intersect::exact_pcurve_over(&curve, range, geometry, tol)
            else {
                return Ok(false);
            };
            if let Some((lo, hi)) = turn_into {
                let u = image.point_at(f64::midpoint(range.0, range.1), tol)?.x;
                let k = ((lo - u) / period).ceil();
                let shifted = u + k * period;
                if shifted > hi {
                    return Ok(false);
                }
                if k != 0.0 {
                    image = image.transformed(
                        &ogeom_math::Transform2::translation(ogeom_math::Vector2::new(
                            k * period,
                            0.0,
                        )),
                        tol,
                    )?;
                }
            }
            ogeom_algo::attach_pcurve(model, edge, image, surface, Location::identity(), range)?;
        }
    }
    Ok(true)
}

/// A ring's blend with its rings in pieces, facing out of the material:
/// turned for a cut, as the wedge faces for a fill.
fn blend_face(
    model: &mut Model,
    p: &Placed,
    m: &Made,
    fill: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let ring = &p.ring;
    let Some(NodeData::Face(data)) = model.node(&ring.blend).map(|n| n.data()) else {
        return Ok(None);
    };
    let mut data = (**data).clone();
    data.triangulation = None;
    data.natural_restriction = false;
    let Some(mut rings) = walk_face(model, &ring.blend)? else {
        return Ok(None);
    };
    let [walked] = rings.as_mut_slice() else {
        return Ok(None);
    };
    let mut edges: Vec<Shape> = Vec::with_capacity(walked.len() + 2);
    for w in walked.iter() {
        let id = w.edge.node();
        if id == ring.wall_ring.node() {
            edges.extend(m.wall.instead_of(&w.edge));
        } else if id == ring.cap_ring.node() {
            edges.extend(m.cap.instead_of(&w.edge));
        } else if id == ring.meridian.node() {
            edges.push(m.meridian.oriented(w.edge.orientation()));
        } else {
            return Ok(None);
        }
    }
    // The node holds the loop as the wedge's face walks it where that
    // face is used forward; for a cut the blend is then used the other way
    // round, facing the other way.
    if ring.blend.orientation() == Orientation::Reversed {
        edges = edges.iter().rev().map(Shape::reversed).collect();
    }
    let used = match (ring.blend.orientation(), fill) {
        (own, true) => own,
        (Orientation::Forward, false) => Orientation::Reversed,
        (_, false) => Orientation::Forward,
    };
    let wire = make_wire(model, &edges, tol)?.shape;
    let face = model.add_face(data, &[wire])?.oriented(used);
    Ok(Some(face))
}

/// Whether the faces walk every edge they hold twice, once each way.
fn walked_twice(model: &Model, faces: &[Shape]) -> OgeomResult<bool> {
    let mut uses: FastMap<TShapeId, (usize, usize)> = FastMap::default();
    for face in faces {
        let Some(rings) = walk_face(model, face)? else {
            return Ok(false);
        };
        for walked in rings.iter().flatten() {
            let slot = uses.entry(walked.edge.node()).or_default();
            if walked.edge.orientation() == Orientation::Reversed {
                slot.1 += 1;
            } else {
                slot.0 += 1;
            }
        }
    }
    Ok(uses.values().all(|&(f, r)| f == 1 && r == 1))
}

#[cfg(test)]
mod tests {
    use super::hull_distance;

    #[test]
    fn a_box_s_shadow_stands_off_the_axis_by_its_nearest_side() {
        let square = [(2.0, -1.0), (4.0, -1.0), (4.0, 1.0), (2.0, 1.0)];
        assert!((hull_distance(&square) - 2.0).abs() < 1e-12);
        let around = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
        assert_eq!(hull_distance(&around), 0.0);
        let corner = [(3.0, 4.0), (6.0, 8.0), (3.0, 8.0)];
        assert!((hull_distance(&corner) - 5.0).abs() < 1e-12);
    }
}
