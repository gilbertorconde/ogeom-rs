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

/// The stations along a two-plane run where each contact is asked to lie
/// on its face, the middles of as many equal stretches as the swept
/// check's first reading. A contact found off its face proves the round
/// would stand on no face there; a face notched between two stations is
/// not seen.
const CONTACT_STATIONS: u32 = 16;

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
///   included). Where the shared edges stop short of closing round a loop
///   because another face stands over the crease between them (a wedge
///   fused into the corner), while the ball touches both faces all the way
///   round a closed seat, the blend is that closed round, built as between
///   faces sharing no edge: the other face is taken in whole where every
///   ball clears it.
/// - A curved face and another face sharing no edge, their surfaces
///   sharing a direction (planes and cylinders along it) or an axis (planes
///   square to it, cylinders, cones, spheres and tori about it): the ball
///   touches each face along a line or a circle, and the blend is a
///   cylinder or a torus of `radius` over the run both faces reach, or all
///   the way round, ended in the ball's section through the end of the
///   solid's edges along the crease where those end first. The corner it
///   takes off or fills is the section between the ball's arc and where
///   the two surfaces cross, swept over the run; the faces are trimmed to
///   their lines of contact by the boolean that applies it. Where the
///   surfaces cross in more than one corner the ball rolls in (a drum's
///   flat and its wall), every such corner is blended. A strip between the
///   faces (an older, smaller round) that lies wholly in the corner goes
///   with it.
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
///   the meeting line, or the ball touching both touches a face's plane
///   off the face at a station of the run (a face narrower than the ball's
///   setback), or another face of the solid reaches into its lune, the
///   part of it between its centre and its arc, at a station of the run (a
///   strip between them, such as a chamfer, that the ball crosses rather
///   than takes off whole). A face reaching into the back of the ball only,
///   behind its centre (the far side of a part thinner than the ball), is
///   clear of the round;
/// - for a curved face, if the two are one face or either is not a face
///   of the solid;
/// - for a curved face sharing edges with the other, as
///   [`fillet_edges`](crate::fillet_edges), and, where the ball touches
///   both faces all the way round a closed seat the shared edges do not
///   close, if it reaches into another face of the solid at a station of
///   the run;
/// - for a curved face sharing no edge with the other, if their surfaces
///   do not cross; and, for each side of the faces the ball could roll on,
///   if no ball touches both faces there, its round would be a torus
///   crossing its own axis, or the solid beside the middle of the round is
///   not what that side needs (material behind the faces, open in front of
///   them), or another face of the solid reaches into the ball's lune at a
///   station of the run (a chamfer or a wedge between them it crosses),
///   other than a face the round runs off at the run's ends or along the
///   crease, or the corners the ball rolls in on a side stand close enough
///   that their rounds could meet. It is also refused where a ball seats on
///   both sides;
/// - for a marched seat, on a side, if no ball seats along where the
///   surfaces cross, more than one crossing holds one touching both faces,
///   the ball leaves the two faces at different places round its seat and
///   never both at once, its line of contact enters a face more than once,
///   or the march stops while the ball still touches both faces; and if a
///   face is a trimmed or offset surface.
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where whether
/// another face comes inside the ball between two stations cannot be told
/// within the halvings allowed.
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

    // The ball touching both planes, behind both on a convex corner and in
    // front of both on a concave one: its centre `o` off the line has
    // `n_a . o = n_b . o = -r` (or `r`), so `o = -r (n_a + n_b) / (1 + n_a . n_b)`.
    // It is the only ball touching both planes on that side, so where another
    // face of the solid reaches inside it (a strip between the faces, such as
    // a chamfer, crossing its arc) no ball rolls in the corner touching the
    // two faces alone, and the wedge would leave that face standing between
    // two pieces of round. The planes are not parallel, so `1 + n_a . n_b`
    // is positive.
    let k = normal_a.dot(normal_b);
    let side = if convex { -1.0 } else { 1.0 };
    let centre = seed + (normal_a + normal_b) * (side * radius / (1.0 + k));
    let contacts = [
        centre - normal_a * (side * radius),
        centre - normal_b * (side * radius),
    ];
    if let Some(at) = crate::sheet_curved::face_inside_rolling_ball(
        model,
        solid,
        [a, b],
        [centre, seed],
        contacts,
        along,
        radius,
        (lo, hi - lo),
        tol,
    )? {
        ogeom_bail!(
            Construction,
            "no ball of radius {radius} rolls between the two faces: {}",
            stands_inside(at)
        );
    }

    // Each contact must lie on its own face along the run. A contact off
    // the face (the face narrower than the ball's setback) touches only
    // the face's plane, and the round's rail there would stand on no face:
    // a fill hanging past the face's far edge, or a cut running into the
    // face beyond it.
    let deflection = ogeom_mesh::Deflection {
        chord: (radius * 1e-3).max(tol.confusion() * 1e2),
        ..ogeom_mesh::Deflection::default()
    };
    for (face, contact, name) in [(a, contacts[0], "first"), (b, contacts[1], "second")] {
        for k in 0..CONTACT_STATIONS {
            let t = (hi - lo).mul_add((f64::from(k) + 0.5) / f64::from(CONTACT_STATIONS), lo);
            let at = contact + along * t;
            if ogeom_algo::classify_on_face(model, face, at, deflection, tol)?
                == ogeom_algo::Containment::Out
            {
                ogeom_bail!(
                    Construction,
                    "no ball of radius {radius} rolls between the two faces: the ball touching \
                     both touches the {name} face's plane at ({:.6}, {:.6}, {:.6}), off the face",
                    at.x,
                    at.y,
                    at.z
                );
            }
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
    seated_fillet(model, solid, &seat, radius, None, &[], tol)
}

/// Why a ball touching both faces rounds nothing: it reaches into another
/// face of the solid at `at`.
fn stands_inside(at: Point) -> String {
    format!(
        "the ball touching both cuts into another face of the solid at ({:.6}, {:.6}, {:.6}), \
         a face standing between them (a strip such as a chamfer) that a round of the two \
         faces alone would leave standing between two pieces of round",
        at.x, at.y, at.z
    )
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
    // Shared edges that close round a loop are the whole corner. Open ones
    // end where the crease leaves a face, or where another face of the
    // solid stands over it (a wedge fused into the corner) while the ball
    // still rolls on the two faces past it. Where the ball rolls touching
    // both all the way round one closed seat, its round is the blend
    // whatever stands over the crease: built where every ball clears the
    // other faces, refused where one reaches into a face. Otherwise the
    // shared edges are the corner.
    if !shared.is_empty() {
        if closes_round(model, &shared)? {
            return crate::fillet::fillet_edges(model, solid, &shared, radius, tol);
        }
        let held = crate::sheet_curved::surfaces_of(model, [a, b])?;
        let reads = match read_sides(model, solid, a, b, radius, tol) {
            Ok(reads) => reads,
            Err(ogeom_core::OgeomError::Construction(_) | ogeom_core::OgeomError::NotDone(_)) => {
                Vec::new()
            }
            Err(e) => return Err(e),
        };
        let mut round = reads
            .into_iter()
            .filter(|read| read.closed && read.corners == 1)
            .filter_map(|read| match read.made {
                Made::Rounds(corner) => Some(Ok((corner, read.behind))),
                Made::Crosses(at) => Some(Err(at)),
                Made::Misses(_) => None,
            });
        if let (Some(found), None) = (round.next(), round.next()) {
            let (corner, behind) = match found {
                Ok(seated) => seated,
                Err(at) => ogeom_bail!(
                    Construction,
                    "no ball of radius {radius} rounds a corner between the two faces: the \
                     ball rolls on both all the way round, and {}",
                    stands_inside(at)
                ),
            };
            return apply(model, solid, corner, behind, radius, tol);
        }
        crate::sheet_curved::put_back(model, held);
        return crate::fillet::fillet_edges(model, solid, &shared, radius, tol);
    }

    // The faces' windows, widened to read the ball, are put back where no
    // blend is made.
    let held = crate::sheet_curved::surfaces_of(model, [a, b])?;
    let blended = blend_apart(model, solid, a, b, radius, tol);
    if blended.is_err() {
        crate::sheet_curved::put_back(model, held);
    }
    blended
}

/// The blend between two faces that share no edge.
fn blend_apart(
    model: &mut Model,
    solid: &Shape,
    a: &Shape,
    b: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let reads = read_sides(model, solid, a, b, radius, tol)?;
    let mut seated = Vec::new();
    let mut misses: Vec<(&str, String)> = Vec::new();
    for read in reads {
        let side = if read.behind { "behind" } else { "in front of" };
        match read.made {
            Made::Rounds(corner) => seated.push((corner, read.behind)),
            Made::Crosses(at) => misses.push((side, stands_inside(at))),
            Made::Misses(why) => misses.push((side, why)),
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
    apply(model, solid, wedge, behind, radius, tol)
}

/// What a ball on one side of two faces makes of the corner between them.
enum Made {
    /// It rolls on the two faces alone, beside the corner that side needs.
    Rounds(Corner),
    /// It rolls beside the corner that side needs, but reaches into another
    /// face of the solid at this point.
    Crosses(Point),
    /// It rounds no corner on that side, for this reason.
    Misses(String),
}

/// One side's reading: behind both faces or in front of them, whether the
/// ball's seat closes round a loop touching both faces all the way, how
/// many corners the ball rolls in there, and what it makes of them.
struct SideRead {
    behind: bool,
    closed: bool,
    corners: usize,
    made: Made,
}

/// The ball read on each side of two faces that share no edge, or whose
/// shared edges leave part of the crease: behind both, where it rolls in
/// the material, then in front of both. The corner is swept from one
/// section where the surfaces share a direction or an axis, and marched
/// along where they cross otherwise. A side counts only where the solid
/// agrees, the corner beside the round being material for the first and
/// open for the second.
fn read_sides(
    model: &mut Model,
    solid: &Shape,
    a: &Shape,
    b: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<SideRead>> {
    let swept = crate::sheet_curved::closed_form_between(model, [a, b], tol)?;
    let hosts = if swept {
        None
    } else {
        Some(crate::pair_marched::hosts_of(model, [a, b], radius, tol)?)
    };
    let mut reads = Vec::with_capacity(2);
    for behind in [true, false] {
        let built = match &hosts {
            None => crate::sheet_curved::face_seats(model, [a, b], radius, behind, tol).and_then(
                |seats| {
                    let closed = seats.iter().all(crate::sheet_curved::FaceSeat::closed);
                    let mut corners = Vec::with_capacity(seats.len());
                    for seat in &seats {
                        let (wedge, [corner, lune], run) =
                            crate::sheet_curved::corner_wedge(model, seat, solid, tol)?;
                        let inside = crate::sheet_curved::face_inside_seat(
                            model,
                            solid,
                            [a, b],
                            seat,
                            run,
                            tol,
                        )?;
                        // A lune no face stands in is all material or all
                        // open, as the ball's side is, whatever a strip
                        // between the faces (an older, smaller round) leaves
                        // of the corner.
                        let probe = if inside.is_none() { lune } else { corner };
                        corners.push((wedge, probe, inside));
                    }
                    Ok((Corners::Swept(corners), closed))
                },
            ),
            Some(hosts) => crate::pair_marched::pair_seat(
                model,
                [a, b],
                hosts,
                radius,
                behind,
                Some(solid),
                tol,
            )
            .and_then(|seat| {
                let Some(probe) = seat.probe() else {
                    ogeom_bail!(Invariant, "a marched seat on the crease has no probe");
                };
                let closed = seat.closed();
                let inside =
                    crate::pair_marched::face_inside(model, solid, [a, b], &seat, radius, tol)?;
                Ok((Corners::Marched(Box::new(seat), probe, inside), closed))
            }),
        };
        let (corners, closed) = match built {
            Ok(built) => built,
            Err(
                ogeom_core::OgeomError::Construction(why) | ogeom_core::OgeomError::NotDone(why),
            ) => {
                reads.push(SideRead {
                    behind,
                    closed: false,
                    corners: 0,
                    made: Made::Misses(why.to_string()),
                });
                continue;
            }
            Err(e) => return Err(e),
        };
        // Every corner on the side must round for the side to.
        let read = |model: &mut Model, probe: Point, inside: Option<Point>| -> OgeomResult<_> {
            let material = matches!(
                ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?,
                ogeom_algo::Containment::In
            );
            Ok(match (material == behind, inside) {
                (true, Some(at)) => Some(Made::Crosses(at)),
                (true, None) => None,
                (false, _) => {
                    let found = if material { "material" } else { "open" };
                    Some(Made::Misses(format!(
                        "the solid beside the round is {found}"
                    )))
                }
            })
        };
        let count = match &corners {
            Corners::Swept(corners) => corners.len(),
            Corners::Marched(..) => 1,
        };
        let made = match corners {
            Corners::Swept(corners) => {
                let mut wedges = Vec::with_capacity(corners.len());
                let mut fails = None;
                for (wedge, probe, inside) in corners {
                    if let Some(made) = read(model, probe, inside)? {
                        fails = Some(made);
                        break;
                    }
                    wedges.push(wedge);
                }
                fails.unwrap_or(Made::Rounds(Corner::Swept(wedges)))
            }
            Corners::Marched(seat, probe, inside) => {
                read(model, probe, inside)?.unwrap_or(Made::Rounds(Corner::Marched(seat)))
            }
        };
        reads.push(SideRead {
            behind,
            closed,
            corners: count,
            made,
        });
    }
    Ok(reads)
}

/// The corners read on one side before the solid is asked: each swept
/// corner with its probe and any point of another face inside its ball's
/// lune, or the marched seat with its own.
enum Corners {
    Swept(Vec<(Shape, Point, Option<Point>)>),
    Marched(Box<crate::pair_marched::PairSeat>, Point, Option<Point>),
}

/// Whether edges close round loops: every vertex ending one of them ends
/// two.
fn closes_round(model: &Model, edges: &[Shape]) -> OgeomResult<bool> {
    let mut ends: Vec<(ogeom_topo::TShapeId, usize)> = Vec::new();
    for edge in edges {
        let Some((first, last)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(false);
        };
        for vertex in [first, last] {
            match ends.iter_mut().find(|(node, _)| *node == vertex.node()) {
                Some((_, count)) => *count += 1,
                None => ends.push((vertex.node(), 1)),
            }
        }
    }
    Ok(ends.iter().all(|&(_, count)| count == 2))
}

/// The corner taken off (`behind`) or filled.
fn apply(
    model: &mut Model,
    solid: &Shape,
    wedge: Corner,
    behind: bool,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    match wedge {
        Corner::Swept(wedges) => {
            // Corners that stand apart, taken off or filled one after the
            // other.
            let mut built: Option<Built> = None;
            for wedge in &wedges {
                let from = built.as_ref().map_or(solid, |b| &b.shape).clone();
                let next = if behind {
                    ogeom_bool::cut(model, &from, wedge, tol)?
                } else {
                    ogeom_bool::fuse(model, &from, wedge, tol)?
                };
                built = Some(match built {
                    None => next,
                    Some(earlier) => {
                        let history = earlier.history.then(&next.history);
                        Built::new(next.shape, history)
                    }
                });
            }
            let Some(built) = built else {
                ogeom_bail!(Invariant, "a swept blend has no corner to take off or fill");
            };
            Ok(built)
        }
        Corner::Marched(seat) => {
            crate::pair_marched::apply_to_solid(model, solid, *seat, radius, behind, tol)
        }
    }
}

/// The corners a curved face blend takes off or fills: the swept solids
/// themselves, one for each corner the faces share, or the marched seat its
/// wedge is built from.
enum Corner {
    Swept(Vec<Shape>),
    Marched(Box<crate::pair_marched::PairSeat>),
}
