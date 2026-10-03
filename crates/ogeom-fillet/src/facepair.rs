//! The blend between two faces of one solid, which need not share an edge.
//!
//! A rolling ball does not care whether the solid has an edge where the two
//! supports would meet. It cares where they *would* meet (for two planes,
//! their own line of intersection) and rolls in the corner that line
//! defines. So a face-face blend between planes is the edge blend seated on
//! a line the solid does not have: found from the planes, cut back to the
//! stretch both faces actually reach, and handed to the same wedge
//! construction.
//!
//! A step is the shape that names the case: a tall block beside a low one,
//! the tall one's wall and the low one's lid facing each other across a
//! corner that belongs to neither.
//!
//! Curved faces that meet along edges of the solid are blended along those
//! edges by the edge blend, exact or marched. Curved faces that share no
//! edge are blended where the ball's section is the same all along the
//! corner (surfaces sharing a direction or an axis) by sweeping the corner
//! between the round and the crease from that section, and otherwise by
//! marching the ball round where the two surfaces cross and building the
//! corner from the marched band, as the edge blend does along an edge.
//! Either corner is cut off or fused on.

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
/// across a corner, and the material must fill that corner (the ball rolls
/// inside it and the blend takes the corner off) or leave it open (the
/// ball rolls in the open and the blend fills it), which is asked of the
/// solid rather than assumed from the normals, because normals cannot tell
/// a step from a slot.
///
/// - Two planes: their meeting line must cross the stretch both faces
///   reach, and the blend is a cylinder of `radius` along it.
/// - A curved face and another face sharing edges of the solid: the corner
///   is those edges, and the blend is the edge blend of all of them at
///   once ([`fillet_edges`](crate::fillet_edges)): a cylinder or torus
///   where the seat has a closed form, a B-spline band fitted through the
///   marched ball otherwise (B-spline, cone, sphere and torus hosts
///   included).
/// - A curved face and another face sharing no edge, their surfaces
///   sharing a direction (planes and cylinders along it) or an axis (planes
///   square to it, cylinders, cones, spheres and tori about it): the ball
///   touches each face along a line or a circle, and the blend is a
///   cylinder or a torus of `radius` over the run both faces reach, or all
///   the way round, ended in the ball's section through the end of the
///   solid's edges along the crease where those end first. The corner it takes off or fills is the section between
///   the ball's arc and where the two surfaces cross, swept over the run;
///   the faces are trimmed to their lines of contact by the boolean that
///   applies it.
/// - A curved face and another face sharing no edge whose surfaces share
///   no direction or axis (a B-spline face, a cylinder at a slant to a
///   plane): the ball is marched along where the two surfaces cross, the
///   surfaces carried on past the faces for it (a plane's or a cylinder's
///   window widened, a B-spline patch continued, each written back as the
///   face's own surface). The blend is a B-spline band fitted through the
///   ball's arcs, and the corner between the band and the crease, bounded
///   by the two surfaces, is cut off or fused on as the edge blend's is:
///   all the way round where the ball touches both faces round a closed
///   seat, and otherwise over the stretch where it touches both, capped at
///   each end by a plane face in the ball's section where its line of
///   contact leaves a face (a crease running off the solid, or a face that
///   holds part of the seat), or, where the solid's edges along the crease
///   end first (a crease running off a block's side while the ball still
///   rolls over its top), in the section through their end. The band rides
///   both faces within a tenth of a degree.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, if the radius is not a usable length, and:
///
/// - for two planes, if they are parallel or the faces do not both reach
///   the meeting line;
/// - for a curved face, if the two are one face or either is not a face
///   of the solid;
/// - for a curved face sharing edges with the other, as
///   [`fillet_edges`](crate::fillet_edges);
/// - for a curved face sharing no edge with the other, if their surfaces
///   do not cross; and, for each side of the faces the ball could roll on,
///   if no ball touches both faces there, its round would be a torus
///   crossing its own axis, or the solid beside the middle of the round is
///   not what that side needs (material behind the faces, open in front of
///   them). It is also refused where a ball seats on both sides;
/// - for a marched seat, on a side, if no ball seats along where the
///   surfaces cross, more than one crossing holds one touching both faces,
///   the ball leaves the two faces at different places round its seat and
///   never both at once, its line of contact enters a face more than once,
///   or the march stops while the ball still touches both faces; and if a
///   face is a trimmed or offset surface.
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
    if !(is_planar(model, a)? && is_planar(model, b)?) {
        return curved_blend(model, solid, a, b, radius, tol);
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
    // normals is material either way (that is what makes both faces
    // outward-facing), so it tells a convex corner from a concave one not
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
pub(crate) fn planar_face_of(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Point, Vector)> {
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Construction, "expected a planar face");
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
pub(crate) fn meet(
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

/// Whether a face lies on a plane.
fn is_planar(model: &Model, face: &Shape) -> OgeomResult<bool> {
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "expected a face");
    };
    Ok(matches!(
        model.geometry().surface(data.surface),
        Some(SurfaceGeometry::Plane(_))
    ))
}

/// The blend between two faces at least one of which is curved.
///
/// Faces that meet along edges of the solid have their corner there: the
/// ball rolls along those edges, and the blend is the edge blend of all of
/// them at once, exact where the seat has a closed form and marched where
/// it does not. Faces that share no edge are blended from one section
/// swept along the corner where their surfaces share a direction or an
/// axis, and from the ball marched round where the surfaces cross
/// otherwise.
fn curved_blend(
    model: &mut Model,
    solid: &Shape,
    a: &Shape,
    b: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if a.is_same(b) {
        ogeom_bail!(Construction, "a face is not blended against itself");
    }
    let own = explore_unique(model, solid, ShapeType::Face)?;
    for (name, face) in [("first", a), ("second", b)] {
        if !own.iter().any(|f| f.is_same(face)) {
            ogeom_bail!(
                Construction,
                "the {name} face is not a face of the solid being blended"
            );
        }
    }
    let theirs = explore_unique(model, b, ShapeType::Edge)?;
    let shared: Vec<Shape> = explore_unique(model, a, ShapeType::Edge)?
        .into_iter()
        .filter(|e| {
            theirs
                .iter()
                .any(|t| crate::support::same_occurrence(model, t, e, tol))
        })
        .collect();
    if !shared.is_empty() {
        return crate::fillet::fillet_edges(model, solid, &shared, radius, tol);
    }

    // Which corner the ball rounds is read from both sides: behind both
    // faces it rolls in the material and the round takes a convex corner
    // off; in front of both it rolls in the open and fills a concave one.
    // A side counts only where the solid agrees, the corner beside the
    // round being material for the first and open for the second.
    // The corner is swept from one section where the surfaces share a
    // direction or an axis, and marched along where they cross otherwise.
    let swept = crate::sheet_curved::closed_form_between(model, [a, b], tol)?;
    let hosts = if swept {
        None
    } else {
        Some(crate::pair_marched::hosts_of(model, [a, b], radius, tol)?)
    };
    let mut seated = Vec::new();
    let mut misses: Vec<(&str, String)> = Vec::new();
    for behind in [true, false] {
        let side = if behind { "behind" } else { "in front of" };
        let built = match &hosts {
            None => crate::sheet_curved::face_seat(model, [a, b], radius, behind, tol)
                .and_then(|seat| crate::sheet_curved::corner_wedge(model, &seat, solid, tol))
                .map(|(wedge, probe)| (Corner::Swept(wedge), probe)),
            Some(hosts) => crate::pair_marched::pair_seat(
                model,
                [a, b],
                hosts,
                radius,
                behind,
                Some(solid),
                tol,
            )
            .and_then(|seat| match seat.probe() {
                Some(probe) => Ok((Corner::Marched(Box::new(seat)), probe)),
                None => {
                    ogeom_bail!(Invariant, "a marched seat on the crease has no probe")
                }
            }),
        };
        let (wedge, probe) = match built {
            Ok(built) => built,
            Err(
                ogeom_core::OgeomError::Construction(why) | ogeom_core::OgeomError::NotDone(why),
            ) => {
                misses.push((side, why.to_string()));
                continue;
            }
            Err(e) => return Err(e),
        };
        let material = matches!(
            ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?,
            ogeom_algo::Containment::In
        );
        if material == behind {
            seated.push((wedge, behind));
        } else {
            let found = if material { "material" } else { "open" };
            misses.push((side, format!("the solid beside the round is {found}")));
        }
    }
    let (wedge, behind) = match seated.len() {
        1 => seated.remove(0),
        0 => {
            // A reason that holds on both sides is said once.
            let said = if misses.len() == 2 && misses[0].1 == misses[1].1 {
                misses[0].1.clone()
            } else {
                misses
                    .iter()
                    .map(|(side, why)| format!("{side} them, {why}"))
                    .collect::<Vec<_>>()
                    .join("; ")
            };
            ogeom_bail!(
                Construction,
                "no ball of radius {radius} rounds a corner between the two faces: {said}"
            )
        }
        _ => ogeom_bail!(
            Construction,
            "a ball of radius {radius} rounds a corner both behind and in front of the two \
             faces; which one to blend is ambiguous"
        ),
    };
    match wedge {
        Corner::Swept(wedge) if behind => ogeom_bool::cut(model, solid, &wedge, tol),
        Corner::Swept(wedge) => ogeom_bool::fuse(model, solid, &wedge, tol),
        Corner::Marched(seat) => {
            crate::pair_marched::apply_to_solid(model, solid, *seat, radius, behind, tol)
        }
    }
}

/// The corner a curved face blend takes off or fills: the swept solid
/// itself, or the marched seat its wedge is built from.
enum Corner {
    Swept(Shape),
    Marched(Box<crate::pair_marched::PairSeat>),
}
