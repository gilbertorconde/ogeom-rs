//! A corner one ball rounds, its creases each built as it is alone (ADR
//! 0001, amendment "a corner one ball rounds, from its bands alone"): three
//! convex creases of one radius meeting at a vertex of three faces, each a
//! plane or a wall holding the material inside it, where one ball a radius
//! in from all three faces fits.
//!
//! Locally the solid is convex there, and the round is the rolling ball's:
//! what no ball within the solid reaches. A point of a crease's spandrel
//! (between its faces and the ball swept along it) is reached by no ball
//! within those two faces, so by none within the solid: each crease's
//! round alone, carried past its ends as it is alone, is part of the
//! corner's round. What is left of it is the corner's cell: the points
//! whose nearest centre of a ball within the solid is the corner ball's
//! (the cone of the three faces' normals at its centre) farther than a
//! radius from it. Each crease's alone round is a tool, and the cell (that
//! cone, held within the three faces' tangent planes moved a radius out,
//! less the ball) is one more; they are applied as the one-boolean
//! fallback's.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::Curved;
use super::defeatured::Rounder;
use super::fallback::{Sense, Tool};

/// `once` rounded on `solid` where it is three convex creases of one
/// radius meeting at a vertex of exactly three edges and three faces, each
/// face a plane or a cylinder or sphere holding the material inside it,
/// where one ball a radius in from the three faces touches each within its
/// trim and no face of the solid stands nearer than the radius; `None`
/// otherwise, or where a crease's round alone, a tool or the boolean will
/// not build, or the result is not valid.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] where the work is cancelled, and reading the
/// solid.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let [(_, r), (_, r1), (_, r2)] = once else {
        return Ok(None);
    };
    let r = *r;
    if *r1 != r || *r2 != r || !(r.is_finite() && r > tol.confusion()) {
        return Ok(None);
    }
    let Some(curved) = super::curved::read(model, solid, tol)? else {
        return Ok(None);
    };
    let ids: Vec<_> = once.iter().map(|(e, _)| super::id(e)).collect();
    let mut ends = Vec::with_capacity(3);
    for e in &ids {
        let Some(edge) = curved.edges.get(e) else {
            return Ok(None);
        };
        ends.push(edge.ends);
    }
    // The vertex all three share, with no other edge.
    let shared: Vec<_> = [ends[0].0, ends[0].1]
        .into_iter()
        .filter(|v| ends.iter().all(|(a, b)| a == v || b == v))
        .collect();
    let [v] = shared[..] else {
        return Ok(None);
    };
    let at_v = curved
        .edges
        .values()
        .filter(|edge| edge.ends.0 == v || edge.ends.1 == v)
        .count();
    if at_v != 3 {
        return Ok(None);
    }
    let mut faces: Vec<usize> = Vec::new();
    for e in &ids {
        for f in super::trimmed::faces_of(&curved, *e) {
            if !faces.contains(&f) {
                faces.push(f);
            }
        }
    }
    let [fa, fb, fc] = faces[..] else {
        return Ok(None);
    };
    let faces = [fa, fb, fc];
    if faces.iter().any(|f| !holds_inside(&curved, *f)) {
        return Ok(None);
    }
    let corner = curved.points[&v];
    let Some(centre) = ball_centre(&curved, faces, corner, r) else {
        return Ok(None);
    };
    let mut normals = [Vector::ZERO; 3];
    for (k, f) in faces.iter().enumerate() {
        let Some((_, n)) = past(&curved, *f, centre) else {
            return Ok(None);
        };
        normals[k] = n;
    }
    let mut scratch = model.clone();
    if !fits(&mut scratch, solid, &curved, faces, centre, normals, r, tol)? {
        return Ok(None);
    }
    // Each crease's round alone: what it takes off, its band's faces named
    // for its edge.
    let mut tools = Vec::with_capacity(4);
    for (edge, radius) in once {
        let alone = match round(&mut scratch, solid, &[(edge.clone(), *radius)], tol) {
            Ok(Some(alone)) => alone,
            Err(e @ OgeomError::Cancelled) => return Err(e),
            _ => return Ok(None),
        };
        let (Ok(taken), Ok(added)) = (
            ogeom_bool::cut(&mut scratch, solid, &alone.shape, tol),
            ogeom_bool::cut(&mut scratch, &alone.shape, solid, tol),
        ) else {
            return Ok(None);
        };
        if !explore_unique(&scratch, &added.shape, ShapeType::Solid)?.is_empty()
            || explore_unique(&scratch, &taken.shape, ShapeType::Solid)?.is_empty()
        {
            return Ok(None);
        }
        let mut named = Vec::new();
        for band in alone.history.generated(edge) {
            if scratch.kind_of(band)? != ShapeType::Face {
                continue;
            }
            for face in taken.history.trace(band) {
                named.push((face.clone(), edge.clone()));
            }
        }
        if named.is_empty() {
            return Ok(None);
        }
        tools.push(Tool {
            solid: taken.shape,
            sense: Sense::Cut,
            faces: named,
        });
    }
    // The corner's cell: the cone of the faces' normals at the centre, held
    // within their tangent planes at the ball's touches moved a radius out
    // (the solid stands within them there), less the ball.
    let mut planes: Vec<(Vector, f64)> = Vec::with_capacity(6);
    for k in 0..3 {
        let (ni, nj, nk) = (normals[k], normals[(k + 1) % 3], normals[(k + 2) % 3]);
        let m = ni.cross(nj);
        let out = if m.dot(nk) > 0.0 { -m } else { m };
        let len = out.magnitude();
        if len <= tol.angular() {
            return Ok(None);
        }
        let out = out / len;
        planes.push((out, out.dot(centre - Point::ORIGIN)));
        planes.push((ni, ni.dot(centre - Point::ORIGIN) + 2.0 * r));
    }
    let Some(cell) = super::fallback::polytope(&mut scratch, &planes, tol)? else {
        return Ok(None);
    };
    let frame = Frame::new(
        centre,
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        tol,
    )?;
    let ball = ogeom_algo::make_sphere(&mut scratch, frame, r, tol)?.shape;
    let Ok(cut) = ogeom_bool::cut(&mut scratch, &cell, &ball, tol) else {
        return Ok(None);
    };
    let vertex = curved.shapes[&v].clone();
    let mut named = Vec::new();
    for face in explore_unique(&scratch, &cut.shape, ShapeType::Face)? {
        let on_ball = scratch
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| scratch.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Sphere(_)));
        if on_ball {
            named.push((face, vertex.clone()));
        }
    }
    if named.is_empty() {
        return Ok(None);
    }
    tools.push(Tool {
        solid: cut.shape,
        sense: Sense::Cut,
        faces: named,
    });
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(&scratch, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    *model = scratch;
    Ok(Some(built))
}

/// Whether face `f` is a plane, or a cylinder or sphere holding the
/// material inside it: the solid stands within its tangent planes near it.
fn holds_inside(curved: &Curved, f: usize) -> bool {
    let face = &curved.faces[f];
    match face.surface {
        SurfaceGeometry::Plane(_) => true,
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Sphere(_) => !face.reversed,
        _ => false,
    }
}

/// How far `p` stands past face `f`'s surface along its outward normal,
/// and that normal, for a plane, a cylinder or a sphere; `None` otherwise
/// or on an axis or centre.
fn past(curved: &Curved, f: usize, p: Point) -> Option<(f64, Vector)> {
    let face = &curved.faces[f];
    let sign = if face.reversed { -1.0 } else { 1.0 };
    match &face.surface {
        SurfaceGeometry::Plane(plane) => {
            let plane = plane.plane();
            let n = plane.normal().vector() * sign;
            Some(((p - plane.origin()).dot(n), n))
        }
        SurfaceGeometry::Cylinder(c) => {
            let cylinder = c.cylinder();
            let axis = cylinder.axis();
            let (o, a) = (axis.location, axis.direction.vector());
            let radial = (p - o) - a * (p - o).dot(a);
            let m = radial.magnitude();
            (m > 0.0).then(|| (sign * (m - cylinder.radius()), radial / m * sign))
        }
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let q = p - sphere.centre();
            let m = q.magnitude();
            (m > 0.0).then(|| (sign * (m - sphere.radius()), q / m * sign))
        }
        _ => None,
    }
}

/// The centre a radius in from all three faces, by Newton's method from
/// the corner moved in along its faces' normals: the distance past each
/// face is the radius in, and its gradient the face's outward normal at the
/// centre's foot. `None` where it does not converge.
fn ball_centre(curved: &Curved, faces: [usize; 3], corner: Point, r: f64) -> Option<Point> {
    let mut inward = Vector::ZERO;
    for f in faces {
        inward -= past(curved, f, corner)?.1;
    }
    let len = inward.magnitude();
    if len <= f64::EPSILON {
        return None;
    }
    let mut c = corner + inward * (r * 3.0_f64.sqrt() / len);
    for _ in 0..64 {
        let mut rows = [Vector::ZERO; 3];
        let mut rhs = [0.0; 3];
        for (k, f) in faces.iter().enumerate() {
            let (d, n) = past(curved, *f, c)?;
            rows[k] = n;
            rhs[k] = -(d + r);
        }
        let det = rows[0].dot(rows[1].cross(rows[2]));
        if det.abs() <= 1e-12 {
            return None;
        }
        // Cramer's rule on rows · step = rhs.
        let step = (rows[1].cross(rows[2]) * rhs[0]
            + rows[2].cross(rows[0]) * rhs[1]
            + rows[0].cross(rows[1]) * rhs[2])
            / det;
        c += step;
        if step.magnitude() <= 1e-14 * (1.0 + r) {
            let residual = faces
                .iter()
                .map(|f| past(curved, *f, c).map_or(f64::INFINITY, |(d, _)| (d + r).abs()))
                .fold(0.0, f64::max);
            return (residual <= 1e-12 * (1.0 + r)).then_some(c);
        }
    }
    None
}

/// Whether the ball about `centre` touches each of the three faces within
/// its trim (its foot there on the face) and no face of the solid stands
/// nearer than its radius.
#[allow(clippy::too_many_arguments, reason = "one corner's reading")]
fn fits(
    scratch: &mut Model,
    solid: &Shape,
    curved: &Curved,
    faces: [usize; 3],
    centre: Point,
    normals: [Vector; 3],
    r: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let options = ogeom_intersect::ExtremaOptions::default();
    for (k, f) in faces.iter().enumerate() {
        let foot = ogeom_algo::make_vertex(scratch, centre + normals[k] * r).shape;
        let off = ogeom_algo::distance_between_shapes(
            scratch,
            &foot,
            &curved.faces[*f].shape,
            options,
            tol,
        )?;
        if off.distance > 10.0 * tol.confusion() {
            return Ok(false);
        }
    }
    let at = ogeom_algo::make_vertex(scratch, centre).shape;
    for face in explore_unique(scratch, solid, ShapeType::Face)? {
        let off = ogeom_algo::distance_between_shapes(scratch, &at, &face, options, tol)?;
        if off.distance < r - 10.0 * tol.confusion() {
            return Ok(false);
        }
    }
    Ok(true)
}
