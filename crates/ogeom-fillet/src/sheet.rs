//! Rounds on sheets: where two faces of an open shell meet at an edge, and
//! in the corner between two faces of separate shapes.
//!
//! A solid's fillet is a wedge cut away or fused on. A sheet has no volume
//! for a boolean to take a wedge from, so here the faces are rebuilt: each
//! face the ball touches loses the strip between the corner and its line of
//! contact, its boundary is re-routed along that line, and the round shares
//! the edge there with it. The shared edges are built once and used by both
//! faces, so the result's topology is exact by construction rather than
//! sewn within a tolerance.
//!
//! Between two planes the rolling ball's envelope is a cylinder whose axis
//! is where the two planes' offsets meet, and the lines of contact are
//! straight: every curve and every trim here is closed-form. Curved faces
//! are rounded in `sheet_curved`: exactly where the supports share a
//! direction or an axis, marched otherwise.

use core::f64::consts::PI;

use ogeom_algo::{
    Built, History, attach_pcurve, edge_vertices, make_edge_between, make_half_space,
    make_natural_face, make_vertex, make_wire, shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, CylinderSurface, PlaneSurface, SurfaceGeometry};
use ogeom_geom::{Curve3d as _, Surface as _};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{
    FaceData, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId, explore_unique,
};

use crate::facepair::{meet, planar_face_of};
use crate::support::{edge_curve, face_from_edges, segment_between};

/// Round edges of a sheet where two of its faces meet, with a rolling ball
/// of `radius`.
///
/// The ball rolls on the concave side of each corner: the side where the
/// two faces make an angle under a half turn. Each of the two faces is
/// rebuilt on its own surface without the strip between the edge and the
/// line the ball touches it along, and the round, tangent to both, fills
/// the gap, sharing an edge with each. Its normal follows the sheet's:
/// where the faces' normals point into the corner, the round's points to
/// the ball's centre.
///
/// The round's surface:
///
/// - between two planes, or planes and cylinders all along the straight
///   edge's direction, a cylinder of `radius` along it;
/// - between planes square to an axis and cylinders, cones, spheres and
///   tori about it, meeting along a circle or arc about the axis, a torus
///   about the axis with `radius` for its tube;
/// - otherwise (B-spline, swept or other faces, or an edge with no closed
///   form) a B-spline surface fitted through the ball's arcs as it is
///   marched along an open edge, its borders the lines of contact, within
///   two tenths of a micron at unit scale; the edges it shares are widened
///   to say so.
///
/// The round ends in the ball's section at each end of an open edge (the
/// plane square to a straight edge, the half-plane through the axis of an
/// arc, the march's section elsewhere), and those ends become part of the
/// sheet's free boundary. A closed edge is rounded all the way round.
///
/// Edges are rounded one after another, each on the sheet the previous one
/// left, so two rounds may trim the same face from different sides.
///
/// History: the sheet is modified into the result; every rebuilt face and
/// every shortened boundary edge is modified into its new self; each
/// rounded edge is deleted and generates its round.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, where:
///
/// - `radius` is not a positive length, or `edges` is empty;
/// - `sheet` is not a shell, or is placed;
/// - an edge is not an edge of the sheet with a face on each side (a free
///   edge, an edge of three faces, an edge an earlier round in the same
///   call consumed);
/// - a face beside an edge is placed, or runs along it twice;
/// - the two faces continue each other, fold onto each other, or are
///   oriented inconsistently across the edge;
/// - no ball of `radius` seats in the corner, a ball does not fit inside a
///   curved face it rolls on, it touches the faces at the ends of a
///   diameter, or its round would be a torus crossing its own axis;
/// - the edge is closed and its round would have to be marched;
/// - the ball does not seat between the faces at an end of a marched edge;
/// - at an end of the edge, a face's boundary does not leave it along an
///   edge of the sheet's free boundary lying in the round's end section (for
///   two planes, square to the edge and straight), or other faces meet there
///   too;
/// - the ball sets back past the far end of a face's boundary edge, or some
///   other part of a face's boundary comes within reach of the strip the
///   round takes (checked against a sampled distance held low by the
///   sampling step, so an edge only near the strip may be refused).
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a march
/// gives too few stations or its band misses its fit target.
pub fn fillet_sheet_edges(
    model: &mut Model,
    sheet: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    usable_radius(radius, tol)?;
    if edges.is_empty() {
        ogeom_bail!(Construction, "no edges were named to round");
    }
    let kind = model.kind_of(sheet)?;
    if kind != ShapeType::Shell {
        ogeom_bail!(
            Construction,
            "a sheet's edges are rounded on a shell; got a {kind:?}"
        );
    }
    if !sheet.location().is_identity() {
        ogeom_bail!(
            Construction,
            "a placed sheet is not rounded; bake its placement into its geometry first"
        );
    }
    let mut current = sheet.clone();
    let mut steps = Vec::with_capacity(edges.len());
    for edge in edges {
        if model.kind_of(edge)? != ShapeType::Edge {
            ogeom_bail!(Construction, "a sheet is rounded along its edges");
        }
        let step = round_sheet_edge(model, &current, edge, radius, tol)?;
        current = step.shape.clone();
        steps.push(step.history);
    }
    Ok(Built::new(current, History::chain(&steps)))
}

/// Round the corner between two faces of separate shapes with a rolling
/// ball of `radius`.
///
/// The ball rolls on the side each face's normal points to: it touches
/// the front of both. Reverse a face to roll it on that face's other side.
///
/// The faces' surfaces must share a direction or an axis, which gives the
/// round a closed form:
///
/// - two planes meet on a line, and planes and cylinders all along one
///   direction share it: the ball's centre runs along that direction, it
///   touches each face along a line, and the round is a cylinder of
///   `radius` running over the stretch both faces reach (where each face's
///   line of contact crosses it), ending square to the direction;
/// - planes square to an axis and cylinders, cones, spheres and tori about
///   it share the axis: the ball's centre turns about it, it touches each
///   face along a circle, and the round is a torus about the axis with
///   `radius` for its tube, running over the turn both faces reach and
///   ending in half-planes through the axis, or all the way round where
///   both lines of contact close on their faces.
///
/// With `trim`, each face is cut back to its line of contact (everything of
/// it on the round's side of that line goes) and the result is one shell
/// of the two trimmed faces and the round, sharing an edge along each line
/// of contact; where a face's line of contact runs past the round's ends
/// (or closes round on itself), the rest of it stays free boundary.
/// Without it the faces are left as they are and the result is
/// the round alone, its normal following the faces' as in the shell.
///
/// History: each face generates the round; with `trim`, each face is
/// modified into its trimmed self.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, where:
///
/// - `radius` is not a positive length;
/// - an argument is not a face, or the two are one face;
/// - a face is placed;
/// - the faces' surfaces share no direction or axis (two B-spline faces,
///   a cylinder and a plane at a slant to it, and the like);
/// - two planes are parallel, or no ball of `radius` touches the front of
///   both surfaces, or it does not fit inside a curved face, or more than
///   one does and touches both faces;
/// - the ball touching both surfaces misses a face, or touches the two at
///   the ends of a diameter, or its round would be a torus crossing its
///   own axis;
/// - for two planes, a face does not reach its line of contact (it lies
///   wholly beyond it, or wholly between it and the corner);
/// - a line of contact crosses its face more than once;
/// - the stretches the two faces reach do not overlap;
/// - one face's line of contact closes on itself across the face's own
///   boundary (round a seam) while the round spans only part of the turn;
/// - with `trim`, cutting a face back leaves it in more than one piece, or
///   leaves a placed face.
pub fn fillet_faces(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    radius: f64,
    trim: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    usable_radius(radius, tol)?;
    let mut planar = true;
    for (name, face) in [("first", a), ("second", b)] {
        let kind = model.kind_of(face)?;
        if kind != ShapeType::Face {
            ogeom_bail!(
                Construction,
                "fillet_faces rounds between two faces; the {name} argument is a {kind:?}"
            );
        }
        planar &= matches!(surface_of(model, face)?, SurfaceGeometry::Plane(_));
    }
    if a.is_same(b) {
        ogeom_bail!(Construction, "a face is not rounded against itself");
    }
    if !planar {
        return crate::sheet_curved::fillet_faces(model, a, b, radius, trim, tol);
    }
    let (origin_a, normal_a) = planar_face_of(model, a, tol)?;
    let (origin_b, normal_b) = planar_face_of(model, b, tol)?;
    let cross = normal_a.cross(normal_b);
    let magnitude = cross.magnitude();
    if magnitude <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the two faces are parallel; their planes meet nowhere for a ball to roll along"
        );
    }
    let along = cross / magnitude;
    let foot = meet(origin_a, normal_a, origin_b, normal_b, along, tol)?;

    // The ball in front of both: a centre `radius` in front of each plane.
    // Its contact with one plane then lies on the side the other plane's
    // normal points to, and that direction, within the plane and square to
    // the meeting line, is the in-plane part of the other's normal.
    let cosine = normal_a.dot(normal_b);
    let inward = [
        (normal_b - normal_a * cosine).normalized(tol)?,
        (normal_a - normal_b * cosine).normalized(tol)?,
    ];
    let setback = Round::between([foot, foot + along], inward, radius, tol)?.setback;
    let reaches = [
        reach(model, a, foot, inward[0] * setback, along, tol)?,
        reach(model, b, foot, inward[1] * setback, along, tol)?,
    ];
    let lo = reaches[0].span.0.max(reaches[1].span.0);
    let hi = reaches[0].span.1.min(reaches[1].span.1);
    if hi - lo <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the two faces reach the round over stretches of the corner that do not overlap"
        );
    }
    let round = Round::between([foot + along * lo, foot + along * hi], inward, radius, tol)?;
    let (s, e) = (round.start, 1 - round.start);

    let mut history = History::new();
    if !trim {
        let vertices = round
            .contacts
            .map(|row| row.map(|p| make_vertex(model, p).shape));
        let rails = [
            segment_between(
                model,
                (&vertices[0][s], round.contacts[0][s]),
                (&vertices[0][e], round.contacts[0][e]),
                tol,
            )?,
            segment_between(
                model,
                (&vertices[1][s], round.contacts[1][s]),
                (&vertices[1][e], round.contacts[1][e]),
                tol,
            )?,
        ];
        let face = round_face(model, &round, &vertices, &rails, normal_a, tol)?;
        history.generate(a, face.clone());
        history.generate(b, face.clone());
        return Ok(Built::new(face, history));
    }

    // Each trimmed face meets the round along the stretch of its contact
    // edge the round spans. Where the round ends short of an end of that
    // edge, the edge is split there and the rest stays free boundary.
    let mut flats = Vec::with_capacity(2);
    let mut rows = Vec::with_capacity(2);
    let mut runs = Vec::with_capacity(2);
    for (f, (reached, normal)) in reaches.iter().zip([normal_a, normal_b]).enumerate() {
        let flat = Flat::read(model, &reached.face)?;
        if flat.normal.dot(normal) <= 0.0 {
            ogeom_bail!(
                Invariant,
                "trimming a face turned it over; the round would face the wrong way"
            );
        }
        let (w, i) = flat.locate(&reached.edge)?;
        let Some((p, q)) = edge_vertices(model, &flat.loops[w][i])? else {
            ogeom_bail!(Construction, "the face's line of contact has no vertices");
        };
        let (pp, qp) = (vertex_point(model, &p)?, vertex_point(model, &q)?);
        // The end of the round each end of the contact edge is nearer.
        let p_first = (pp - foot).dot(along) <= (qp - foot).dot(along);
        let (k_p, k_q) = if p_first { (0, 1) } else { (1, 0) };
        let mut at = |vertex: &Shape, point: Point, k: usize| {
            if point.distance(round.contacts[f][k]) <= tol.confusion() {
                vertex.clone()
            } else {
                make_vertex(model, round.contacts[f][k]).shape
            }
        };
        let (at_p, at_q) = (at(&p, pp, k_p), at(&q, qp, k_q));
        rows.push(if p_first { [at_p, at_q] } else { [at_q, at_p] });
        runs.push(((w, i), (p, pp), (q, qp), (k_p, k_q)));
        flats.push(flat);
    }
    let vertices = [rows[0].clone(), rows[1].clone()];
    let rails = [
        flats[0].segment(
            model,
            (&vertices[0][s], round.contacts[0][s]),
            (&vertices[0][e], round.contacts[0][e]),
            tol,
        )?,
        flats[1].segment(
            model,
            (&vertices[1][s], round.contacts[1][s]),
            (&vertices[1][e], round.contacts[1][e]),
            tol,
        )?,
    ];
    let blend = round_face(model, &round, &vertices, &rails, normal_a, tol)?;

    let mut shell_faces = Vec::with_capacity(3);
    for (f, (flat, ((w, i), (p, pp), (q, qp), (k_p, k_q)))) in flats.iter().zip(runs).enumerate() {
        let mut pieces = Vec::with_capacity(3);
        if !vertices[f][k_p].is_same(&p) {
            pieces.push(flat.segment(
                model,
                (&p, pp),
                (&vertices[f][k_p], round.contacts[f][k_p]),
                tol,
            )?);
        }
        pieces.push(if k_p == s {
            rails[f].clone()
        } else {
            rails[f].reversed()
        });
        if !vertices[f][k_q].is_same(&q) {
            pieces.push(flat.segment(
                model,
                (&vertices[f][k_q], round.contacts[f][k_q]),
                (&q, qp),
                tol,
            )?);
        }
        let face = flat.rebuilt(model, w, i, 1, pieces, tol)?;
        history.modify(if f == 0 { a } else { b }, face.clone());
        shell_faces.push(face);
    }
    shell_faces.insert(1, blend.clone());
    let shell = model.add_shell(&shell_faces)?;
    history.generate(a, blend.clone());
    history.generate(b, blend);
    Ok(Built::new(shell, history))
}

/// Refuse a radius that rounds nothing.
fn usable_radius(radius: f64, tol: Tolerances) -> OgeomResult<()> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a round of radius {radius} rounds nothing");
    }
    Ok(())
}

/// The rolling ball's envelope between two planes over one stretch of the
/// line they meet on.
struct Round {
    /// The stretch's ends on the meeting line, in the caller's order.
    ends: [Point; 2],
    /// Where the ball touches each face at each end: `contacts[face][end]`.
    contacts: [[Point; 2]; 2],
    /// The ball's centre at each end.
    centres: [Point; 2],
    /// The direction the cylinder's axis runs: the arc from the first
    /// face's contact to the second's turns positively about it.
    axis: Vector,
    /// The end the cylinder's height starts at.
    start: usize,
    /// The ball's radius.
    radius: f64,
    /// The angle the arc sweeps: the corner's supplement.
    sweep: f64,
    /// How far each line of contact stands from the corner, along its face.
    setback: f64,
}

impl Round {
    /// The ball seated in the corner the two faces make along the stretch
    /// from `ends[0]` to `ends[1]`, each face leaving the corner along its
    /// unit `sides` direction, square to the stretch.
    ///
    /// The corner's opening `θ` is the angle between the two directions.
    /// The ball's centre lies on their bisector at `r / sin(θ/2)` from the
    /// corner, and touches each face `r / tan(θ/2)` along it.
    fn between(
        ends: [Point; 2],
        sides: [Vector; 2],
        radius: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let [first, second] = sides;
        let opening = first.cross(second).magnitude().atan2(first.dot(second));
        if opening <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the two faces fold onto each other; no ball fits between them"
            );
        }
        if opening >= PI - tol.angular() {
            ogeom_bail!(
                Construction,
                "the two faces continue each other across the corner; there is no corner to round"
            );
        }
        let half = opening / 2.0;
        let setback = radius / half.tan();
        let depth = radius / half.sin();
        let bisector = (first + second).normalized(tol)?;
        let travel = ends[1] - ends[0];
        let length = travel.magnitude();
        if length <= tol.confusion() {
            ogeom_bail!(Construction, "the stretch to round has no length");
        }
        let along = travel / length;
        let centres = ends.map(|p| p + bisector * depth);
        let contacts = [
            ends.map(|p| p + first * setback),
            ends.map(|p| p + second * setback),
        ];
        let turn = (contacts[0][0] - centres[0])
            .cross(contacts[1][0] - centres[0])
            .dot(along);
        let (axis, start) = if turn > 0.0 { (along, 0) } else { (-along, 1) };
        Ok(Self {
            ends,
            contacts,
            centres,
            axis,
            start,
            radius,
            sweep: PI - opening,
            setback,
        })
    }

    /// The stretch's length.
    fn length(&self) -> f64 {
        self.ends[0].distance(self.ends[1])
    }
}

/// The round's face, sharing its lines of contact with the faces it meets:
/// `vertices[face][end]` are the contacts' vertices and `rails[face]` the
/// edge along each line of contact, running from the end the cylinder's
/// height starts at. Oriented so its normal agrees with `normal`, the first
/// face's, along their shared edge.
fn round_face(
    model: &mut Model,
    round: &Round,
    vertices: &[[Shape; 2]; 2],
    rails: &[Shape; 2],
    normal: Vector,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (s, e) = (round.start, 1 - round.start);
    let base = round.centres[s];
    let radial = (round.contacts[0][s] - base) / round.radius;
    let axis_dir = Direction::new(round.axis, tol)?;
    let radial_dir = Direction::new(radial, tol)?;
    let length = round.length();
    // An arc of the ball's circle at height `h` up the axis, from the first
    // face's contact to the second's: its frame is the cylinder's moved up
    // the axis, so its trim on the cylinder is a straight line.
    let arc = |model: &mut Model, h: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let frame = Frame::new(base + round.axis * h, axis_dir, radial_dir, tol)?;
        let circle = Circle::new(frame, round.radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, round.sweep), from, to, tol)?.shape)
    };
    let frame = Frame::new(base, axis_dir, radial_dir, tol)?;
    let surface = CylinderSurface::new(Cylinder::new(frame, round.radius, tol)?, (0.0, length))?;
    let edges = vec![
        arc(model, 0.0, &vertices[0][s], &vertices[1][s])?,
        rails[1].clone(),
        arc(model, length, &vertices[0][e], &vertices[1][e])?.reversed(),
        rails[0].reversed(),
    ];
    let face = face_from_edges(model, surface.into(), &edges, tol)?;
    // The cylinder's own normal points away from its axis.
    Ok(if radial.dot(normal) > 0.0 {
        face
    } else {
        face.reversed()
    })
}

/// A planar face as a rebuild reads it: its surface, its normal as the
/// face is oriented, and its boundary as stored.
struct Flat {
    face: Shape,
    surface_id: SurfaceId,
    surface: SurfaceGeometry,
    frame: Frame,
    normal: Vector,
    tolerance: ogeom_core::Tolerance,
    wires: Vec<Shape>,
    loops: Vec<Vec<Shape>>,
}

impl Flat {
    fn read(model: &Model, face: &Shape) -> OgeomResult<Self> {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        if !face.location().is_identity() || !data.location.is_identity() {
            ogeom_bail!(
                Construction,
                "a placed face is not rounded; bake its placement into its geometry first"
            );
        }
        let (surface_id, tolerance) = (data.surface, data.tolerance);
        let surface = surface_of(model, face)?;
        let SurfaceGeometry::Plane(plane) = &surface else {
            ogeom_bail!(
                Construction,
                "only planar faces are rounded here; a face on a {:?} surface is not",
                surface.kind()
            );
        };
        let frame = plane.plane().frame();
        let mut normal = frame.z().vector();
        if face.orientation() == Orientation::Reversed {
            normal = -normal;
        }
        let wires = model.children_of(&face.oriented(Orientation::Forward))?;
        let mut loops = Vec::with_capacity(wires.len());
        for wire in &wires {
            loops.push(model.children_of(&wire.oriented(Orientation::Forward))?);
        }
        Ok(Self {
            face: face.clone(),
            surface_id,
            surface,
            frame,
            normal,
            tolerance,
            wires,
            loops,
        })
    }

    /// Where `edge` stands in the boundary: its wire and its place in it.
    fn locate(&self, edge: &Shape) -> OgeomResult<(usize, usize)> {
        let mut hits = Vec::new();
        for (w, ring) in self.loops.iter().enumerate() {
            for (i, e) in ring.iter().enumerate() {
                if e.is_same(edge) {
                    hits.push((w, i));
                }
            }
        }
        match hits.as_slice() {
            [one] => Ok(*one),
            [] => ogeom_bail!(Construction, "the edge is not on the face's boundary"),
            _ => ogeom_bail!(
                Construction,
                "the face's boundary runs along the edge twice (a seam); it has no one side \
                 to round"
            ),
        }
    }

    /// A straight edge between two vertices in the face's plane, its exact
    /// trim on that plane attached.
    fn segment(
        &self,
        model: &mut Model,
        from: (&Shape, Point),
        to: (&Shape, Point),
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let edge = segment_between(model, from, to, tol)?;
        let (curve, range) = edge_curve(model, &edge, tol)?;
        let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &self.surface, tol) else {
            ogeom_bail!(
                Invariant,
                "a straight edge in a face's plane has no trim on that plane"
            );
        };
        attach_pcurve(
            model,
            &edge,
            pcurve,
            self.surface_id,
            Location::identity(),
            range,
        )?;
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) {
            data.assert_same_parameter(true);
        }
        Ok(edge)
    }

    /// The face again, on its own surface, with `count` edges of wire `w`
    /// from `first` on (cyclically) replaced by `pieces`, which must walk
    /// from where the replaced run began to where it ended.
    fn rebuilt(
        &self,
        model: &mut Model,
        w: usize,
        first: usize,
        count: usize,
        pieces: Vec<Shape>,
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let ring = &self.loops[w];
        let n = ring.len();
        let mut edges: Vec<Shape> = (count..n).map(|j| ring[(first + j) % n].clone()).collect();
        edges.extend(pieces);
        let wire = make_wire(model, &edges, tol)?
            .shape
            .oriented(self.wires[w].orientation());
        let mut wires = self.wires.clone();
        wires[w] = wire;
        let mut data = FaceData::new(self.surface_id, Location::identity());
        data.tolerance = self.tolerance;
        Ok(model
            .add_face(data, &wires)?
            .oriented(self.face.orientation()))
    }

    /// The unit direction, in the plane and square to the edge at `(w, i)`,
    /// in which the face lies beside it.
    ///
    /// Read from the winding of the edge's loop in the plane's chart: the
    /// face lies left of an outer loop walked counter-clockwise and right
    /// of a hole's. The outer loop is the one enclosing the most area.
    fn inward(&self, model: &Model, w: usize, i: usize, tol: Tolerances) -> OgeomResult<Vector> {
        let chart = |p: Point| {
            let local = self.frame.to_local(p);
            (local.x, local.y)
        };
        let mut areas = Vec::with_capacity(self.loops.len());
        for ring in &self.loops {
            let mut polygon = Vec::new();
            for edge in ring {
                let mut samples = edge_samples(model, edge, 8, tol)?;
                if edge.orientation() == Orientation::Reversed {
                    samples.reverse();
                }
                samples.pop();
                polygon.extend(samples.into_iter().map(chart));
            }
            let mut area = 0.0;
            for (k, p) in polygon.iter().enumerate() {
                let q = polygon[(k + 1) % polygon.len()];
                area += p.0 * q.1 - p.1 * q.0;
            }
            areas.push(area);
        }
        let outer = areas
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.abs().total_cmp(&y.1.abs()))
            .map_or(0, |(k, _)| k);
        let Some((from, to)) = edge_vertices(model, &self.loops[w][i])? else {
            ogeom_bail!(Construction, "the edge has no vertices");
        };
        let tangent = vertex_point(model, &to)? - vertex_point(model, &from)?;
        let left = self.frame.z().vector().cross(tangent);
        let on_left = (areas[w] > 0.0) == (w == outer);
        let side = if on_left { left } else { -left };
        side.normalized(tol)
    }
}

/// Points along an edge in its stored direction, placed in the world.
fn edge_samples(
    model: &Model,
    edge: &Shape,
    count: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    let (curve, (lo, hi)) = edge_curve(model, edge, tol)?;
    (0..=count)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let t = lo + (hi - lo) * (k as f64) / (count as f64);
            curve.point_at(t, tol)
        })
        .collect()
}

/// A vertex's point, placed in the world.
fn vertex_point(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(point) = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex().map(|v| v.point))
    else {
        ogeom_bail!(Dangling, "vertex is not in this model");
    };
    Ok(vertex.transform(model.datums())?.apply(point))
}

/// A face's surface, cloned out of the model.
fn surface_of(model: &Model, face: &Shape) -> OgeomResult<SurfaceGeometry> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    Ok(surface.clone())
}

/// What of a face lies beyond its line of contact, the edge it gains along
/// that line, and the stretch of the line it spans, measured from `foot`
/// along `along`.
struct Reach {
    face: Shape,
    edge: Shape,
    span: (f64, f64),
}

/// Cut `face` back to the line through `foot + offset` along `along`,
/// keeping the side `offset` points to.
fn reach(
    model: &mut Model,
    face: &Shape,
    foot: Point,
    offset: Vector,
    along: Vector,
    tol: Tolerances,
) -> OgeomResult<Reach> {
    let contact = foot + offset;
    let inward = offset.normalized(tol)?;
    let mut size = 1.0_f64;
    for corner in shape_bounds(model, face, tol)?.corners() {
        size = size.max(corner.distance(contact));
    }
    let size = size * 4.0;
    let plane = Plane::through(contact, Direction::new(inward, tol)?);
    let wall = make_natural_face(
        model,
        PlaneSurface::over(plane, (-size, size), (-size, size))?.into(),
    )?
    .shape;
    let keep = make_half_space(model, &wall, contact + inward * size * 0.5, tol)?.shape;
    let kept = ogeom_bool::common(model, face, &keep, tol)?.shape;
    let pieces = explore_unique(model, &kept, ShapeType::Face)?;
    let reach = tol.confusion() * 10.0;
    let off_line = |p: Point| {
        let d = p - contact;
        (d - along * d.dot(along)).magnitude()
    };
    let mut found = Vec::new();
    for piece in &pieces {
        for edge in explore_unique(model, piece, ShapeType::Edge)? {
            let (curve, (lo, hi)) = edge_curve(model, &edge, tol)?;
            let Curve::Line(line) = &curve else {
                continue;
            };
            if line.axis().direction.vector().cross(along).magnitude() > tol.angular() {
                continue;
            }
            let (p, q) = (curve.point_at(lo, tol)?, curve.point_at(hi, tol)?);
            if off_line(p) > reach || off_line(q) > reach {
                continue;
            }
            let (tp, tq) = ((p - foot).dot(along), (q - foot).dot(along));
            found.push((piece.clone(), edge, (tp.min(tq), tp.max(tq))));
        }
    }
    match (pieces.len(), found.len()) {
        (0, _) => ogeom_bail!(
            Construction,
            "a face lies wholly between the corner and the line the ball touches it along; \
             it does not reach the round"
        ),
        (_, 0) => ogeom_bail!(
            Construction,
            "a face lies wholly beyond the line the ball would touch it along; the ball \
             does not reach it"
        ),
        (1, 1) => {
            let (face, edge, span) = found.remove(0);
            Ok(Reach { face, edge, span })
        }
        _ => ogeom_bail!(
            Construction,
            "the line the ball touches a face along crosses it more than once; the stretch \
             the round spans is ambiguous"
        ),
    }
}

/// One face beside a sheet edge being rounded: where the edge sits in its
/// boundary and the boundary edges that leave the edge's ends.
struct Side {
    flat: Flat,
    wire: usize,
    at: usize,
    /// The direction, in the face and square to the edge, the face lies in.
    inward: Vector,
    /// The end of the edge (0 its start, 1 its end) the loop arrives at.
    arrives: usize,
    /// The boundary edge arriving at the rounded edge, as the loop walks it.
    before: Shape,
    /// The one leaving it.
    after: Shape,
    /// Where `before` starts.
    from: (Shape, Point),
    /// Where `after` ends.
    to: (Shape, Point),
}

/// Round one edge of a shell.
#[allow(clippy::too_many_lines, reason = "one rebuild, checked then assembled")]
fn round_sheet_edge(
    model: &mut Model,
    shell: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let faces = explore_unique(model, shell, ShapeType::Face)?;
    let mut face_edges = Vec::with_capacity(faces.len());
    for face in &faces {
        face_edges.push(explore_unique(model, face, ShapeType::Edge)?);
    }
    let users_of = |e: &Shape| -> Vec<usize> {
        face_edges
            .iter()
            .enumerate()
            .filter(|(_, edges)| edges.iter().any(|x| x.is_same(e)))
            .map(|(k, _)| k)
            .collect()
    };
    let users = users_of(edge);
    match users.len() {
        0 => ogeom_bail!(
            Construction,
            "the edge is not an edge of the sheet (an earlier round in the same call may \
             have consumed it)"
        ),
        1 => ogeom_bail!(
            Construction,
            "the edge bounds one face of the sheet only; a round needs a face on each side"
        ),
        2 => {}
        n => ogeom_bail!(
            Construction,
            "the edge is shared by {n} faces of the sheet; a round has two sides"
        ),
    }
    let (curve, _) = edge_curve(model, edge, tol)?;
    let planar = [&faces[users[0]], &faces[users[1]]]
        .into_iter()
        .map(|face| surface_of(model, face))
        .collect::<OgeomResult<Vec<_>>>()?
        .iter()
        .all(|s| matches!(s, SurfaceGeometry::Plane(_)));
    let Curve::Line(line) = &curve else {
        return crate::sheet_curved::round_sheet_edge(
            model,
            shell,
            edge,
            [&faces[users[0]], &faces[users[1]]],
            &|e| users_of(e).len(),
            radius,
            tol,
        );
    };
    if !planar {
        return crate::sheet_curved::round_sheet_edge(
            model,
            shell,
            edge,
            [&faces[users[0]], &faces[users[1]]],
            &|e| users_of(e).len(),
            radius,
            tol,
        );
    }
    let along = line.axis().direction.vector();
    let Some((v0, v1)) = edge_vertices(model, &edge.oriented(Orientation::Forward))? else {
        ogeom_bail!(Construction, "the edge has no vertices");
    };
    if v0.is_same(&v1) {
        ogeom_bail!(
            Construction,
            "the edge closes on itself; a straight edge between planes has two ends"
        );
    }
    let ends = [vertex_point(model, &v0)?, vertex_point(model, &v1)?];
    let end_vertices = [v0.clone(), v1.clone()];

    let mut sides = Vec::with_capacity(2);
    for &k in &users {
        let flat = Flat::read(model, &faces[k])?;
        let (wire, at) = flat.locate(edge)?;
        let ring = &flat.loops[wire];
        let n = ring.len();
        if n < 3 {
            ogeom_bail!(
                Construction,
                "a face bounded by {n} edges has no boundary edges to end a round on"
            );
        }
        let before = ring[(at + n - 1) % n].clone();
        let after = ring[(at + 1) % n].clone();
        let Some((arrive, _)) = edge_vertices(model, &ring[at])? else {
            ogeom_bail!(Construction, "the edge has no vertices");
        };
        let arrives = usize::from(!arrive.is_same(&end_vertices[0]));
        let (Some((from, _)), Some((_, to))) = (
            edge_vertices(model, &before)?,
            edge_vertices(model, &after)?,
        ) else {
            ogeom_bail!(Construction, "a boundary edge has no vertices");
        };
        let inward = flat.inward(model, wire, at, tol)?;
        let from_point = vertex_point(model, &from)?;
        let to_point = vertex_point(model, &to)?;
        sides.push(Side {
            flat,
            wire,
            at,
            inward,
            arrives,
            before,
            after,
            from: (from, from_point),
            to: (to, to_point),
        });
    }

    // Which side of the sheet the corner is concave on, and whether the
    // two faces agree about it: each face's normal points toward the other
    // face, or both point away.
    let toward = [
        sides[0].flat.normal.dot(sides[1].inward),
        sides[1].flat.normal.dot(sides[0].inward),
    ];
    if (toward[0] > 0.0) != (toward[1] > 0.0) {
        ogeom_bail!(
            Construction,
            "the two faces are oriented inconsistently across the edge; the round has no \
             side to face"
        );
    }
    let round = Round::between(ends, [sides[0].inward, sides[1].inward], radius, tol)?;

    // The ends: square to the edge, straight, free, and long enough.
    for side in &sides {
        for (boundary, near, far) in [
            (&side.before, ends[side.arrives], side.from.1),
            (&side.after, ends[1 - side.arrives], side.to.1),
        ] {
            let (curve, _) = edge_curve(model, boundary, tol)?;
            let Curve::Line(line) = &curve else {
                ogeom_bail!(
                    Construction,
                    "a face's boundary leaves an end of the edge along a curve; the round \
                     ends on straight boundary edges only"
                );
            };
            if line.axis().direction.vector().dot(along).abs() > tol.angular() {
                ogeom_bail!(
                    Construction,
                    "a face's boundary leaves an end of the edge at an angle to it; the round \
                     ends square to the edge, so only square ends are rounded"
                );
            }
            if users_of(boundary).len() != 1 {
                ogeom_bail!(
                    Construction,
                    "the round would run into another face of the sheet at an end of the \
                     edge; only an edge whose ends lie on the sheet's free boundary is rounded"
                );
            }
            let span = (far - near).dot(side.inward);
            if span <= round.setback + tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the round sets back {} from the edge, past the far end of a face's \
                     boundary edge {span} long",
                    round.setback
                );
            }
        }
    }
    // Nothing else meets the edge's ends.
    for vertex in &end_vertices {
        let mut meeting = 0;
        for edge_of in explore_unique(model, shell, ShapeType::Edge)? {
            if let Some((p, q)) = edge_vertices(model, &edge_of)?
                && (p.is_same(vertex) || q.is_same(vertex))
            {
                meeting += 1;
            }
        }
        if meeting != 3 {
            ogeom_bail!(
                Construction,
                "{meeting} edges of the sheet meet at an end of the edge; only an edge whose \
                 ends lie on the sheet's free boundary, with one boundary edge of each face \
                 leaving it, is rounded"
            );
        }
    }
    // Nothing else of either face comes into the strip the round takes.
    let length = round.length();
    for side in &sides {
        for (w, ring) in side.flat.loops.iter().enumerate() {
            for (i, other) in ring.iter().enumerate() {
                let n = ring.len();
                if w == side.wire
                    && (i == side.at || i == (side.at + 1) % n || i == (side.at + n - 1) % n)
                {
                    continue;
                }
                if !clear_of_strip(
                    model,
                    other,
                    ends[0],
                    (ends[1] - ends[0]) / length,
                    side.inward,
                    (length, round.setback),
                    tol,
                )? {
                    ogeom_bail!(
                        Construction,
                        "another part of a face's boundary comes within the strip the round \
                         takes from it; the round would cross it"
                    );
                }
            }
        }
    }

    // Assembly: the contact vertices and rails, the round, the two faces
    // rebuilt around them, and the shell with the faces exchanged.
    let (s, e) = (round.start, 1 - round.start);
    let vertices = round
        .contacts
        .map(|row| row.map(|p| make_vertex(model, p).shape));
    let rails = [
        sides[0].flat.segment(
            model,
            (&vertices[0][s], round.contacts[0][s]),
            (&vertices[0][e], round.contacts[0][e]),
            tol,
        )?,
        sides[1].flat.segment(
            model,
            (&vertices[1][s], round.contacts[1][s]),
            (&vertices[1][e], round.contacts[1][e]),
            tol,
        )?,
    ];
    let blend = round_face(model, &round, &vertices, &rails, sides[0].flat.normal, tol)?;

    let mut history = History::new();
    let mut replaced = Vec::with_capacity(2);
    for (f, side) in sides.iter().enumerate() {
        let (k_in, k_out) = (side.arrives, 1 - side.arrives);
        let lead = side.flat.segment(
            model,
            (&side.from.0, side.from.1),
            (&vertices[f][k_in], round.contacts[f][k_in]),
            tol,
        )?;
        let rail = if k_in == s {
            rails[f].clone()
        } else {
            rails[f].reversed()
        };
        let tail = side.flat.segment(
            model,
            (&vertices[f][k_out], round.contacts[f][k_out]),
            (&side.to.0, side.to.1),
            tol,
        )?;
        let n = side.flat.loops[side.wire].len();
        let face = side.flat.rebuilt(
            model,
            side.wire,
            (side.at + n - 1) % n,
            3,
            vec![lead.clone(), rail, tail.clone()],
            tol,
        )?;
        history.modify(&side.flat.face, face.clone());
        history.modify(&side.before, lead);
        history.modify(&side.after, tail);
        replaced.push((side.flat.face.clone(), face));
    }
    let mut shell_faces: Vec<Shape> = faces
        .iter()
        .map(|face| {
            replaced
                .iter()
                .find(|(old, _)| old.is_same(face))
                .map_or_else(|| face.clone(), |(_, new)| new.clone())
        })
        .collect();
    shell_faces.push(blend.clone());
    let rounded = model.add_shell(&shell_faces)?;
    history.modify(shell, rounded.clone());
    history.delete(edge);
    history.delete(&v0);
    history.delete(&v1);
    history.generate(edge, blend);
    Ok(Built::new(rounded, history))
}

/// Whether an edge keeps clear of the strip a round takes from a face: the
/// rectangle from the corner at `origin`, `extent.0` along `along` and
/// `extent.1` along `inward`, widened by the confusion distance.
///
/// Exact for a straight edge. A curved one is held to its bounding box,
/// which can refuse an edge that only comes near the strip but never
/// passes one that enters it.
fn clear_of_strip(
    model: &Model,
    edge: &Shape,
    origin: Point,
    along: Vector,
    inward: Vector,
    extent: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<bool> {
    let margin = tol.confusion();
    let lo = (-margin, -margin);
    let hi = (extent.0 + margin, extent.1 + margin);
    let flat = |p: Point| {
        let d = p - origin;
        (d.dot(along), d.dot(inward))
    };
    let (curve, (t0, t1)) = edge_curve(model, edge, tol)?;
    if let Curve::Line(_) = curve {
        let a = flat(curve.point_at(t0, tol)?);
        let b = flat(curve.point_at(t1, tol)?);
        return Ok(!segment_meets_box(a, b, lo, hi));
    }
    let (mut min, mut max) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for corner in shape_bounds(model, edge, tol)?.corners() {
        let p = flat(corner);
        min = (min.0.min(p.0), min.1.min(p.1));
        max = (max.0.max(p.0), max.1.max(p.1));
    }
    Ok(max.0 < lo.0 || min.0 > hi.0 || max.1 < lo.1 || min.1 > hi.1)
}

/// Whether the segment from `a` to `b` meets the box from `lo` to `hi`:
/// the segment clipped against each slab in turn.
fn segment_meets_box(a: (f64, f64), b: (f64, f64), lo: (f64, f64), hi: (f64, f64)) -> bool {
    let d = (b.0 - a.0, b.1 - a.1);
    let (mut enter, mut leave) = (0.0_f64, 1.0_f64);
    for (p, q) in [
        (-d.0, a.0 - lo.0),
        (d.0, hi.0 - a.0),
        (-d.1, a.1 - lo.1),
        (d.1, hi.1 - a.1),
    ] {
        if p.abs() <= f64::MIN_POSITIVE {
            if q < 0.0 {
                return false;
            }
            continue;
        }
        let r = q / p;
        if p < 0.0 {
            enter = enter.max(r);
        } else {
            leave = leave.min(r);
        }
        if enter > leave {
            return false;
        }
    }
    true
}
