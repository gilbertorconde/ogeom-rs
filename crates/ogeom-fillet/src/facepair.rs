//! The blend between two faces that share no edge.
//!
//! A rolling ball does not care whether the solid has an edge where the two
//! supports would meet. It cares where they *would* meet — for two planes,
//! their own line of intersection — and rolls in the corner that line
//! defines. So a face-face blend is the edge blend seated on a line the
//! solid does not have: found from the planes, cut back to the stretch both
//! faces actually reach, and handed to the same wedge construction.
//!
//! A step is the shape that names the case: a tall block beside a low one,
//! the tall one's wall and the low one's lid facing each other across a
//! corner that belongs to neither.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

use crate::fillet::seated_fillet;
use crate::support::Seat;

/// Blend two faces of one solid with a rolling ball of the given radius.
///
/// The two faces need not touch. What they must do is face each other
/// across a corner: their planes must meet, both must reach the stretch of
/// that meeting line the blend will sit on, and the material must fill the
/// dihedral between them — which is asked of the solid rather than assumed
/// from the normals, because normals cannot tell a step from a slot.
///
/// Planar supports only. A curved face-face blend needs the marching seat
/// — the spine that is the two offset surfaces' own intersection — which is
/// recorded as owed in `docs/PLAN.md` rather than guessed at here.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// either face is not planar, the planes are parallel, the faces do not
/// both reach the meeting line, or the radius is not a usable length.
pub fn blend_faces(
    model: &mut Model,
    solid: &Shape,
    a: &Shape,
    b: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a blend of radius {radius} rounds nothing");
    }
    let (plane_a, normal_a) = planar_face_of(model, a, tol)?;
    let (plane_b, normal_b) = planar_face_of(model, b, tol)?;

    // The line the two planes meet on: direction from the normals' cross,
    // a point from the two plane equations plus the direction as a third.
    let along = normal_a.cross(normal_b);
    let magnitude = along.magnitude();
    if magnitude <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the two faces are parallel; there is no corner between them"
        );
    }
    let along = along / magnitude;
    let seed = meet(plane_a, normal_a, plane_b, normal_b, along, tol)?;

    // Cut the line back to the stretch both faces reach: each face's own
    // vertices, projected onto it, give the span it can seat a blend over.
    let span = |face: &Shape| -> OgeomResult<(f64, f64)> {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            let Some(point) = model
                .node(&vertex)
                .and_then(|n| n.data().as_vertex().map(|v| v.point))
            else {
                continue;
            };
            let placed = vertex.transform(model.datums())?.apply(point);
            let t = (placed - seed).dot(along);
            lo = lo.min(t);
            hi = hi.max(t);
        }
        if !lo.is_finite() {
            ogeom_bail!(Construction, "a face with no vertices seats nothing");
        }
        Ok((lo, hi))
    };
    let (a0, a1) = span(a)?;
    let (b0, b1) = span(b)?;
    let (lo, hi) = (a0.max(b0), a1.min(b1));
    if hi - lo <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the two faces do not both reach the line their planes meet on, \
             so there is no stretch to seat a blend over"
        );
    }

    // Which way does the corner turn? The solid answers, and the question
    // has to be asked in the right place: the quadrant opposite both
    // normals is material either way — that is what makes both faces
    // outward-facing — so it tells a convex corner from a concave one not
    // at all. The *side* quadrants do. Around a convex edge the material is
    // the opposite quadrant alone; around a concave one, a step's inner
    // corner, it is three of the four, and a side probe lands in it.
    let unit = |v: Vector| -> OgeomResult<Vector> {
        let m = v.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "the faces meet too sharply to seat a blend");
        }
        Ok(v / m)
    };
    let middle = seed + along * f64::midpoint(lo, hi);
    let step = (hi - lo).min(radius) * 1e-3 + tol.confusion();
    let mut convex = true;
    for side in [unit(normal_a - normal_b)?, unit(normal_b - normal_a)?] {
        if matches!(
            ogeom_algo::classify_in_solid_exact(model, solid, middle + side * step, tol)?,
            ogeom_algo::Containment::In
        ) {
            convex = false;
        }
    }

    let seat = Seat {
        start: seed + along * lo,
        end: seed + along * hi,
        along,
        normals: [normal_a, normal_b],
        faces: [a.clone(), b.clone()],
        convex,
    };
    seated_fillet(model, solid, &seat, radius, None, tol)
}

/// A face's plane origin and its outward normal, refusing anything curved.
fn planar_face_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<(Point, Vector)> {
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
        ogeom_bail!(
            Construction,
            "a face-face blend between curved supports needs the marching \
             seat; this is the planar form"
        );
    };
    let placement = face.transform(model.datums())?;
    let origin = placement.apply(plane.plane().frame().origin());
    let mut normal = placement.apply_vector(plane.plane().normal().vector());
    if face.orientation() == ogeom_topo::Orientation::Reversed {
        normal = -normal;
    }
    let magnitude = normal.magnitude();
    if magnitude <= tol.angular() {
        ogeom_bail!(Construction, "a face with no normal faces nothing");
    }
    Ok((origin, normal / magnitude))
}

/// A point on both planes: the one nearest the two origins' midpoint, found
/// by solving the two plane equations with the meeting direction as the
/// third.
fn meet(
    origin_a: Point,
    normal_a: Vector,
    origin_b: Point,
    normal_b: Vector,
    along: Vector,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let rows = [normal_a, normal_b, along];
    let rhs = [
        normal_a.dot(origin_a.to_vector()),
        normal_b.dot(origin_b.to_vector()),
        along.dot(Point::midpoint(origin_a, origin_b).to_vector()),
    ];
    let det = rows[0].dot(rows[1].cross(rows[2]));
    if det.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "the two planes do not meet in a line");
    }
    // The inverse of a three-by-three whose rows are these: its columns are
    // the cross products of the other two rows, over the determinant.
    let c0 = rows[1].cross(rows[2]);
    let c1 = rows[2].cross(rows[0]);
    let c2 = rows[0].cross(rows[1]);
    let v = (c0 * rhs[0] + c1 * rhs[1] + c2 * rhs[2]) / det;
    Ok(Point::ORIGIN + v)
}
