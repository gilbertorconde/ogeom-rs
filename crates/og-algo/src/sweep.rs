//! Sweeping: dragging a shape through space to make one of higher dimension.
//!
//! A vertex sweeps into an edge, an edge into a face, a wire into a shell, a
//! face into a solid. One rule, applied at every level, which is why a prism
//! over a face falls out of the prism over its edges rather than being built
//! separately.
//!
//! # The top is the bottom, moved
//!
//! The far end of a prism is not a copy of the near end. It is the *same*
//! topology node at a different [`Location`] — the shape triple's whole reason
//! for existing (`docs/DATA_MODEL.md` §2). A copy would double the geometry, and
//! then a later edit would have to find and fix both. Sharing means the two ends
//! of a prism cannot drift apart, because there is only one of them.
//!
//! It also means an assembly of a thousand identical extrusions holds one
//! profile and a thousand placements, which is the case the location chain was
//! designed for.
//!
//! # History
//!
//! A swept edge is *both* consumed and generative: it survives as the bottom of
//! the prism and it generates the lateral face. Recording only one of those is
//! the classic way to break downstream naming — a reference to "that edge"
//! resolves to nothing, or a reference to "the face from that edge" does.

use std::collections::HashMap;

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve2d, ExtrusionSurface, Line2d, PlanarCurve, Surface};
use og_math::{Point2, Transform, Vector};
use og_topo::{EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, TShapeId};

use crate::build::{make_face_on, make_shell, make_solid, make_wire};
use crate::history::{Built, History};

/// Roles a sweep assigns.
pub mod roles {
    use og_core::Role;

    /// The face the sweep started from.
    pub const SWEEP_BOTTOM: Role = Role::op_defined(20);
    /// The face the sweep ended at.
    pub const SWEEP_TOP: Role = Role::op_defined(21);
    /// A face swept out by one edge of the profile.
    pub const SWEEP_SIDE: Role = Role::op_defined(22);
    /// An edge swept out by one vertex of the profile.
    pub const SWEEP_RAIL: Role = Role::op_defined(23);
}

/// Extrude a shape along `vector`.
///
/// A face becomes a solid, a wire becomes a shell, an edge becomes a face. The
/// result's history reports each input as generating what it swept out, and the
/// profile itself as surviving into the near end.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `vector` has no
/// length, if the shape is of a kind that cannot be swept, or if an edge of the
/// profile has no 3D curve to sweep.
pub fn make_prism(
    model: &mut Model,
    profile: &Shape,
    vector: Vector,
    tol: Tolerances,
) -> OgResult<Built> {
    if !vector.is_finite() || vector.magnitude() <= tol.confusion() {
        og_bail!(
            Construction,
            "a prism needs a direction to travel; {vector:?} has no length"
        );
    }
    model.begin_operation();

    // One datum for the whole sweep, so every entity at the far end shares a
    // single placement rather than each carrying its own copy of the same
    // transform. Comparing two far-end shapes is then a comparison of one
    // identifier, which is what makes instance detection cheap.
    let datum = model.add_datum(Transform::translation(vector));
    let displacement = Location::of(datum);

    let rails = &mut Rails::new();
    match model.kind_of(profile)? {
        ShapeType::Face => prism_over_face(model, rails, profile, &displacement, vector, tol),
        ShapeType::Wire => {
            let (faces, history) =
                prism_over_wire(model, rails, profile, &displacement, vector, tol)?;
            let shell = make_shell(model, &faces)?.shape;
            Ok(Built::new(shell, history))
        }
        ShapeType::Edge => {
            let (face, history) =
                prism_over_edge(model, rails, profile, &displacement, vector, tol)?;
            Ok(Built::new(face, history))
        }
        other => og_bail!(
            Construction,
            "a {other:?} cannot be swept into anything; sweep an edge, a wire or \
             a face"
        ),
    }
}

/// A face swept into a solid.
fn prism_over_face(
    model: &mut Model,
    rails: &mut Rails,
    face: &Shape,
    displacement: &Location,
    vector: Vector,
    tol: Tolerances,
) -> OgResult<Built> {
    // Which side of the profile the material lands on is decided by the sweep,
    // not by which way the profile was handed over. A profile facing against
    // the sweep does not describe a different solid — it describes the same one
    // from the other side — so it is turned round here and everything below
    // proceeds as if it had faced along all along.
    //
    // Left unturned, both end caps present the wrong side: the mesh still
    // closes and the shell is still closed, so nothing topological notices, and
    // the volume comes back short by twice the caps' contribution.
    let normal = face_normal(model, face, tol)?;
    let travel = vector.magnitude();
    let along = normal.dot(vector) / travel;
    if along.abs() <= tol.angular() {
        og_bail!(
            Construction,
            "the sweep runs along the profile's own surface, so it encloses no \
             volume; a face swept within its own plane is not a solid"
        );
    }
    let profile = if along < 0.0 {
        face.reversed()
    } else {
        face.clone()
    };

    let mut history = History::new();
    let mut faces = Vec::new();

    for wire in model.children_of(&profile)? {
        let (sides, wire_history) =
            prism_over_wire(model, rails, &wire, displacement, vector, tol)?;
        history = history.then(&wire_history);
        faces.extend(sides);
    }

    // The near end faces backwards, because the solid is on the far side of it.
    // Getting this wrong makes a solid that is inside out along one face, and
    // the volume comes out short by exactly that face's contribution rather
    // than obviously wrong.
    let bottom = profile.reversed();
    let top = profile.moved(displacement);
    model.set_derived(&bottom, std::slice::from_ref(face), roles::SWEEP_BOTTOM)?;
    model.set_derived(&top, std::slice::from_ref(face), roles::SWEEP_TOP)?;
    history.generate(face, top.clone());
    faces.push(bottom);
    faces.push(top);

    let shell = make_shell(model, &faces)?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    history.generate(face, solid.clone());
    Ok(Built::new(solid, history))
}

/// Every face a wire sweeps out.
fn prism_over_wire(
    model: &mut Model,
    rails: &mut Rails,
    wire: &Shape,
    displacement: &Location,
    vector: Vector,
    tol: Tolerances,
) -> OgResult<(Vec<Shape>, History)> {
    let mut faces = Vec::new();
    let mut history = History::new();
    for edge in model.ordered_children_of(wire)? {
        let (face, edge_history) = prism_over_edge(model, rails, &edge, displacement, vector, tol)?;
        history = history.then(&edge_history);
        faces.push(face);
    }
    if faces.is_empty() {
        og_bail!(Construction, "a wire with no edges sweeps out nothing");
    }
    Ok((faces, history))
}

/// One edge swept into one face.
///
/// The face's surface is the extrusion of the edge's own 3D curve, so the
/// lateral surface is exact for whatever the edge was — a line gives a plane, an
/// arc gives a cylinder, a spline gives an extruded spline — rather than
/// everything becoming a plane through an approximation.
fn prism_over_edge(
    model: &mut Model,
    rails: &mut Rails,
    edge: &Shape,
    displacement: &Location,
    vector: Vector,
    tol: Tolerances,
) -> OgResult<(Shape, History)> {
    let Some(node) = model.node(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        og_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        og_bail!(
            Construction,
            "an edge with no curve in space has no shape to sweep; a degenerate \
             edge sweeps out nothing and has to be handled by its face, not here"
        );
    };
    let Some(geometry) = model.geometry().curve(*curve).cloned() else {
        og_bail!(Dangling, "curve is not in this model");
    };

    let placement = edge.transform(model.datums())?;
    let (lo, hi) = *range;
    let travel = vector.magnitude();
    let direction = og_math::Direction::new(placement.inverse()?.apply_vector(vector), tol)?;
    let surface = model
        .geometry_mut()
        .add_surface(ExtrusionSurface::new(geometry, direction, travel)?.into());

    // The extrusion's `u` is the *curve's* own parameter, and the curve does
    // not care which way the wire walks it. So a reversed occurrence is
    // traversed from `hi` to `lo`, and the rail its walk starts at stands at
    // `u = hi`, not at `u = lo`.
    //
    // Pinning the rails to `lo` and `hi` regardless — which is what this did —
    // puts each rail's pcurve on the wrong side of the parameter rectangle, and
    // the boundary comes out as a bow tie enclosing nothing. The face then
    // fails to triangulate outright, while the topology looks perfect: the wire
    // closes, the shell closes, and every edge is used twice.
    let reversed = edge.orientation() == Orientation::Reversed;
    let (u_start, u_end) = if reversed { (hi, lo) } else { (lo, hi) };

    // The four sides of the extrusion's parameter rectangle: the edge along the
    // bottom, the same edge displaced along the top, and the two vertical rails
    // its endpoints sweep out.
    let bottom = edge.clone();
    let top = edge.moved(displacement);
    let start_rail = rail(model, rails, edge, displacement, vector, false, tol)?;
    let end_rail = rail(model, rails, edge, displacement, vector, true, tol)?;

    pcurve(model, &bottom, surface, (lo, 0.0), (hi, 0.0), tol)?;
    pcurve(model, &top, surface, (lo, travel), (hi, travel), tol)?;
    if start_rail.is_same(&end_rail) {
        // A closed profile edge — a full circle — starts and ends at one
        // vertex, so its two rails are one edge appearing at both `u = lo` and
        // `u = hi`. That is a seam, and it needs both pcurves: giving it one
        // would leave the face's boundary running up the same side twice and
        // enclosing nothing. Which pcurve is which is decided by the ring
        // below: the rail is walked forward at `u_end` and backward at
        // `u_start`.
        seam_pcurves(
            model,
            &start_rail,
            surface,
            ((u_end, 0.0), (u_end, travel)),
            ((u_start, 0.0), (u_start, travel)),
            tol,
        )?;
    } else {
        pcurve(
            model,
            &start_rail,
            surface,
            (u_start, 0.0),
            (u_start, travel),
            tol,
        )?;
        pcurve(
            model,
            &end_rail,
            surface,
            (u_end, 0.0),
            (u_end, travel),
            tol,
        )?;
    }

    // Round the rectangle: along the bottom, up the far rail, back along the
    // top, down the near rail.
    let ring = [
        bottom.clone(),
        end_rail.clone(),
        top.reversed(),
        start_rail.reversed(),
    ];
    let boundary = make_wire(model, &ring, tol)?.shape;
    let built = make_face_on(model, surface, std::slice::from_ref(&boundary), tol)?.shape;

    // The extrusion's normal is the curve's tangent crossed with the sweep, so
    // it follows the *curve* and not the wire's walk of it. An edge the wire
    // walks backwards therefore makes a face whose default side points into the
    // solid, and the occurrence has to be reversed to present the other one.
    // Every profile with a mixed wire — four of a box's six faces — has some of
    // each, so this cannot be decided once for the profile.
    let face = if reversed { built.reversed() } else { built };
    model.set_derived(&face, std::slice::from_ref(edge), roles::SWEEP_SIDE)?;

    let mut history = History::new();
    // Both, not either. The edge survives as the bottom of the prism *and*
    // makes the lateral face; recording only one is how a reference to "that
    // edge" or to "the face from that edge" ends up resolving to nothing.
    history.generate(edge, face.clone());
    history.generate(edge, top);
    Ok((face, history))
}

/// The direction a face presents, in space.
///
/// Sampled at the mean of its boundary in parameter space, which for a planar
/// profile is exact everywhere and for a curved one is representative: a
/// profile whose normal turns past perpendicular to the sweep somewhere across
/// its own extent sweeps into a solid that passes through itself, and one
/// sample is enough to decide which side the material lands on in every case
/// this can build. A face with no boundary at all covers its whole surface, so
/// the middle of the domain is the point to ask about.
fn face_normal(model: &Model, face: &Shape, tol: Tolerances) -> OgResult<Vector> {
    let Some(node) = model.node(face) else {
        og_bail!(Dangling, "face is not in this model");
    };
    let Some(data) = node.data().as_face() else {
        og_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        og_bail!(Dangling, "face refers to a surface not in this model");
    };

    let mut sum = (0.0, 0.0);
    let mut count = 0_u32;
    // The outer wire is the first, and it alone bounds the region; a hole would
    // only pull the sample towards a point the face does not cover.
    for edge in match model.children_of(face)?.first() {
        Some(outer) => model.children_of(outer)?,
        None => Vec::new(),
    } {
        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let (id, range) = match edge_data.pcurve_for(data.surface, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
            Some(EdgeRepr::Seam { forward, range, .. }) => (*forward, *range),
            _ => continue,
        };
        let Some(pcurve) = model.geometry().pcurve(id) else {
            og_bail!(Dangling, "pcurve is not in this model");
        };
        for at in [range.0, f64::midpoint(range.0, range.1), range.1] {
            let p = pcurve.point_at(at, tol)?;
            sum = (sum.0 + p.x, sum.1 + p.y);
            count += 1;
        }
    }

    let ((ua, ub), (va, vb)) = surface.domain();
    let (u, v) = if count == 0 {
        (f64::midpoint(ua, ub), f64::midpoint(va, vb))
    } else {
        let n = f64::from(count);
        (sum.0 / n, sum.1 / n)
    };
    let normal = surface.normal_at(u, v, tol)?;

    let placed = face
        .transform(model.datums())?
        .apply_vector(normal.vector());
    Ok(if face.orientation() == Orientation::Reversed {
        -placed
    } else {
        placed
    })
}

/// The edge one endpoint of the profile sweeps out.
///
/// Shared between the two faces that meet along it — the previous edge's sweep
/// and this one's — which is what makes the shell close. Building a rail per
/// face instead leaves every one used once and the prism open along every
/// corner.
fn rail(
    model: &mut Model,
    rails: &mut Rails,
    edge: &Shape,
    displacement: &Location,
    vector: Vector,
    at_end: bool,
    tol: Tolerances,
) -> OgResult<Shape> {
    let Some((start, end)) = crate::build::edge_vertices(model, edge)? else {
        og_bail!(
            Construction,
            "an unbounded edge has no endpoints to sweep into rails"
        );
    };
    let base = if at_end { end } else { start };
    let raised = base.moved(displacement);

    // A rail between the same two vertices already exists if a neighbouring
    // edge swept it. Reusing it is not an optimization: two rails between one
    // pair of vertices would leave each used once, and the shell open.
    if let Some(existing) = rails.get(&base.node()) {
        return Ok(existing.clone());
    }

    let Some(node) = model.node(&base) else {
        og_bail!(Dangling, "vertex is not in this model");
    };
    let Some(data) = node.data().as_vertex() else {
        og_bail!(Construction, "vertex node holds no point");
    };
    let from = base.transform(model.datums())?.apply(data.point);

    let line = og_geom::LineCurve::segment(from, from + vector, tol)?;
    let built = crate::build::make_edge_between(
        model,
        line.into(),
        (0.0, vector.magnitude()),
        &base,
        &raised,
        tol,
    )?;
    model.set_derived(&built.shape, std::slice::from_ref(&base), roles::SWEEP_RAIL)?;
    rails.insert(base.node(), built.shape.clone());
    Ok(built.shape)
}

/// The rails built so far in one sweep, keyed by the vertex each rose from.
///
/// Threaded through rather than looked up in the model, because "is there
/// already an edge between these two vertices" is a question the model cannot
/// answer without a search, and the answer is only ever about *this* sweep.
type Rails = HashMap<TShapeId, Shape>;

/// Attach a seam edge's two pcurves, one for each side of the rectangle it
/// bounds twice.
fn seam_pcurves(
    model: &mut Model,
    edge: &Shape,
    surface: og_topo::SurfaceId,
    forward: ((f64, f64), (f64, f64)),
    reversed: ((f64, f64), (f64, f64)),
    tol: Tolerances,
) -> OgResult<()> {
    let flat = |p: (f64, f64)| Point2::new(p.0, p.1);
    let length = flat(forward.0).distance(flat(forward.1));
    let first = model
        .geometry_mut()
        .add_pcurve(Line2d::segment(flat(forward.0), flat(forward.1), tol)?.into());
    let second = model
        .geometry_mut()
        .add_pcurve(Line2d::segment(flat(reversed.0), flat(reversed.1), tol)?.into());

    let Some(node) = model.node_mut(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data_mut() else {
        og_bail!(Construction, "edge node holds no edge data");
    };
    data.add(EdgeRepr::Seam {
        forward: first,
        reversed: second,
        surface,
        location: Location::identity(),
        range: (0.0, length),
    });
    Ok(())
}

/// Attach a straight pcurve between two points of a surface's parameter space.
fn pcurve(
    model: &mut Model,
    edge: &Shape,
    surface: og_topo::SurfaceId,
    from: (f64, f64),
    to: (f64, f64),
    tol: Tolerances,
) -> OgResult<()> {
    let (a, b) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
    let curve: PlanarCurve = Line2d::segment(a, b, tol)?.into();
    // Keyed by the occurrence's own placement. The bottom and top of a prism
    // are one edge node at two locations, running along two different lines of
    // the same parameter space; attached without the placement they would be
    // indistinguishable and the face would collapse onto one of them.
    crate::build::attach_pcurve(
        model,
        edge,
        curve,
        surface,
        edge.location().clone(),
        (0.0, a.distance(b)),
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::build::is_shell_closed;
    use crate::mass::volume_properties;
    use approx::assert_relative_eq;
    use og_math::{Frame, Point};
    use og_mesh::{Deflection, triangulate};
    use og_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn deflection(chord: f64) -> Deflection {
        Deflection {
            chord,
            ..Deflection::default()
        }
    }

    /// One face of a box of `side`, named by its role.
    fn box_face(model: &mut Model, side: f64, role: og_core::Role) -> Shape {
        let built = crate::make_box(model, Frame::WORLD, (side, side, side), T).unwrap();
        explore_unique(model, &built.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| model.provenance_of(f).and_then(og_core::Provenance::role) == Some(role))
            .expect("the box has a face with that role")
    }

    /// A square face in the xy plane, one unit on a side from the origin.
    fn square(model: &mut Model, side: f64) -> Shape {
        box_face(model, side, crate::primitive::roles::FACE_MAX_Z)
    }

    #[test]
    fn a_profile_facing_away_from_the_sweep_gives_the_same_solid_as_one_facing_along_it() {
        // The defect this pins: the `-Z` face of a box has all four of its
        // edges reversed within its wire, and the `+Z` face has none. Sweeping
        // either along `+Z` describes the same solid, so the two had better
        // agree about it — in face count, in mesh closure and in volume.
        for (role, centre) in [
            // The `+Z` face sits at z = 2 and sweeps to z = 5; the `-Z` face
            // sits at z = 0 and sweeps to z = 3.
            (
                crate::primitive::roles::FACE_MAX_Z,
                Point::new(1.0, 1.0, 3.5),
            ),
            (
                crate::primitive::roles::FACE_MIN_Z,
                Point::new(1.0, 1.0, 1.5),
            ),
        ] {
            let mut model = Model::new();
            let face = box_face(&mut model, 2.0, role);
            let built = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 3.0), T).unwrap();

            let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
            assert_eq!(counts(ShapeType::Face), 6, "{role:?}");
            assert_eq!(counts(ShapeType::Edge), 12, "{role:?}");

            // Every face triangulates: a lateral face whose boundary ring runs
            // up the same side twice encloses nothing and fails outright.
            for face in explode(&model, &built.shape) {
                og_mesh::triangulate_face(&model, &face, deflection(0.01), T)
                    .unwrap_or_else(|e| panic!("{role:?}: a face would not triangulate: {e}"));
            }

            let mesh = triangulate(&model, &built.shape, deflection(0.01), T).unwrap();
            assert!(mesh.is_closed(), "{role:?}: the mesh has a slit in it");
            // Positive, and 2 * 2 * 3. A cap left facing inward keeps the mesh
            // closed and takes its own contribution out of the volume twice,
            // which is wrong by an amount nothing else reports.
            assert_relative_eq!(mesh.volume(), 12.0, epsilon = 1e-9);

            let props = volume_properties(&model, &built.shape, deflection(0.01), T).unwrap();
            assert_relative_eq!(props.mass, 12.0, epsilon = 1e-9);
            assert!(
                props.centre.distance(centre) < 1e-9,
                "{role:?}: got {:?}",
                props.centre
            );

            assert!(
                crate::check_tessellation(&model, &built.shape, deflection(0.01), T)
                    .unwrap()
                    .is_valid(),
                "{role:?}: the mesh disagrees with the topology"
            );
        }
    }

    #[test]
    fn every_face_of_a_box_sweeps_into_a_solid_of_the_right_volume() {
        // Four of the six have their wire's edges mixed — some forward, some
        // reversed — which is the case a per-face flip would not have caught.
        use crate::primitive::roles;
        let roles = [
            (roles::FACE_MIN_X, Vector::new(-3.0, 0.0, 0.0)),
            (roles::FACE_MAX_X, Vector::new(3.0, 0.0, 0.0)),
            (roles::FACE_MIN_Y, Vector::new(0.0, -3.0, 0.0)),
            (roles::FACE_MAX_Y, Vector::new(0.0, 3.0, 0.0)),
            (roles::FACE_MIN_Z, Vector::new(0.0, 0.0, -3.0)),
            (roles::FACE_MAX_Z, Vector::new(0.0, 0.0, 3.0)),
        ];
        for (role, vector) in roles {
            let mut model = Model::new();
            let face = box_face(&mut model, 2.0, role);
            let built = make_prism(&mut model, &face, vector, T).unwrap();
            let mesh = triangulate(&model, &built.shape, deflection(0.01), T).unwrap();
            assert!(mesh.is_closed(), "{role:?}: the mesh has a slit in it");
            assert_relative_eq!(mesh.volume(), 12.0, epsilon = 1e-9);
        }
    }

    /// Every face below a shape.
    fn explode(model: &Model, shape: &Shape) -> Vec<Shape> {
        og_topo::explore(model, shape, og_topo::Filter::OfType(ShapeType::Face)).unwrap()
    }

    #[test]
    fn a_square_swept_upward_is_a_box() {
        let mut model = Model::new();
        let face = square(&mut model, 2.0);
        let built = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 3.0), T).unwrap();

        let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
        assert_eq!(counts(ShapeType::Face), 6);
        assert_eq!(counts(ShapeType::Edge), 12);
        assert_eq!(counts(ShapeType::Vertex), 8);

        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(is_shell_closed(&model, &shell).unwrap());

        let props = volume_properties(&model, &built.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(props.mass, 12.0, epsilon = 1e-9);
    }

    #[test]
    fn the_far_end_is_the_same_topology_at_a_different_place() {
        // The point of using a location rather than a copy: one profile, two
        // placements. A copy would double the geometry and let the two ends
        // drift apart under a later edit.
        let mut model = Model::new();
        let face = square(&mut model, 1.0);
        let before = model.node_count();
        let built = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 1.0), T).unwrap();

        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();
        let ends: Vec<&Shape> = faces.iter().filter(|f| f.is_partner(&face)).collect();
        assert_eq!(ends.len(), 2, "both ends share the profile's node");
        assert!(
            !ends[0].is_same(ends[1]),
            "and are still distinct, because their placements differ"
        );

        // Four side faces, four rails, one wire, a shell and a solid — but no
        // second copy of the profile's four edges or four vertices.
        assert!(
            model.node_count() - before < 20,
            "sweeping copied more than it should have: {} new nodes",
            model.node_count() - before
        );
    }

    #[test]
    fn a_swept_edge_is_reported_as_both_surviving_and_generating() {
        // A swept edge is consumed into the bottom of the prism *and* makes the
        // lateral face. Recording one and not the other is how a reference to
        // "that edge", or to "the face it made", resolves to nothing.
        let mut model = Model::new();
        let face = square(&mut model, 1.0);
        let edge = model
            .children_of(&model.children_of(&face).unwrap()[0])
            .unwrap()[0]
            .clone();

        let built = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 1.0), T).unwrap();
        let generated = built.history.generated(&edge);
        assert_eq!(
            generated.len(),
            2,
            "the lateral face and the displaced edge, got {generated:?}"
        );
        assert!(!built.history.is_deleted(&edge), "the edge survives");
    }

    #[test]
    fn an_arc_sweeps_into_a_cylindrical_face_not_a_flat_one() {
        // The lateral surface is the extrusion of the edge's own curve, so it is
        // exact for whatever the edge was. Approximating every side as a plane
        // would make a swept arc visibly faceted and its area wrong.
        let mut model = Model::new();
        let (radius, height) = (2.0_f64, 5.0);
        let cylinder = crate::make_cylinder(&mut model, Frame::WORLD, radius, 1.0, T).unwrap();
        let rim = explore_unique(&model, &cylinder.shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .and_then(og_topo::EdgeData::curve3d)
                    .is_some_and(|r| matches!(r, EdgeRepr::Curve3d { range, .. } if range.1 > 6.0))
            })
            .expect("the cylinder has a full circular rim");

        let built = make_prism(&mut model, &rim, Vector::new(0.0, 0.0, height), T).unwrap();
        assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::Face);

        let mesh = triangulate(&model, &built.shape, deflection(0.005), T).unwrap();
        let area = mesh.area();
        let exact = std::f64::consts::TAU * radius * height;
        assert!(
            area < exact,
            "an inscribed area cannot exceed the surface's"
        );
        assert!(area > exact * 0.999, "{area} against {exact}");
    }

    #[test]
    fn a_wire_sweeps_into_an_open_shell() {
        let mut model = Model::new();
        let face = square(&mut model, 2.0);
        let wire = model.children_of(&face).unwrap()[0].clone();

        let built = make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 3.0), T).unwrap();
        assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::Shell);
        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Face)
                .unwrap()
                .len(),
            4,
            "one side per edge, and no ends"
        );
    }

    #[test]
    fn a_sweep_that_goes_nowhere_is_refused() {
        let mut model = Model::new();
        let face = square(&mut model, 1.0);
        for vector in [
            Vector::ZERO,
            Vector::new(f64::NAN, 0.0, 0.0),
            Vector::new(0.0, 0.0, f64::INFINITY),
        ] {
            assert!(make_prism(&mut model, &face, vector, T).is_err());
        }
    }

    #[test]
    fn a_face_swept_within_its_own_plane_is_refused() {
        // It encloses no volume, and the two ends would land on top of each
        // other. Building it anyway gives a solid whose faces all have area and
        // which measures zero, which is the shape of answer that gets trusted.
        let mut model = Model::new();
        let face = square(&mut model, 1.0);
        let err = make_prism(&mut model, &face, Vector::new(1.0, 1.0, 0.0), T).unwrap_err();
        assert!(
            err.to_string().contains("encloses no volume"),
            "unexpected message: {err}"
        );
        // A wire has no side for the sweep to lie in, so the same vector is
        // fine there — it makes a perfectly good open shell.
        let wire = model.children_of(&face).unwrap()[0].clone();
        assert!(make_prism(&mut model, &wire, Vector::new(1.0, 1.0, 0.0), T).is_ok());
    }

    #[test]
    fn a_vertex_is_not_something_this_sweeps() {
        // A vertex sweeps into an edge, which is a real operation — but it is
        // not one this returns, and claiming otherwise by returning something
        // of the wrong kind would be worse than saying so.
        let mut model = Model::new();
        let vertex = model.add_point(Point::ORIGIN);
        assert!(make_prism(&mut model, &vertex, Vector::Z, T).is_err());
    }
}
