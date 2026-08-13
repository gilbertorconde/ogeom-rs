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

use core::f64::consts::TAU;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d, ExtrusionSurface, Line2d, PlanarCurve, Transformable};
use ogeom_math::{Axis, Circle, Direction, Frame, Point, Point2, Transform, Vector};
use ogeom_topo::{EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, TShapeId};

use crate::build::{make_face_on, make_shell, make_solid, make_wire};
use crate::history::{Built, History};

/// Roles a sweep assigns.
pub mod roles {
    use ogeom_core::Role;

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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `vector` has no
/// length, if the shape is of a kind that cannot be swept, or if an edge of the
/// profile has no 3D curve to sweep.
pub fn make_prism(
    model: &mut Model,
    profile: &Shape,
    vector: Vector,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !vector.is_finite() || vector.magnitude() <= tol.confusion() {
        ogeom_bail!(
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
        other => ogeom_bail!(
            Construction,
            "a {other:?} cannot be swept into anything; sweep an edge, a wire or \
             a face"
        ),
    }
}

/// Sweep a planar face into a tapered prism: every wall leans by `taper`.
///
/// The draft-prism semantics: each section is the profile's own offset at
/// the rate the taper names — the outer loop outward, holes inward — so a
/// positive taper widens the far end and narrows every hole, and each wall
/// makes exactly `taper` with the travel. The far ring is genuinely new
/// topology: a taper breaks the plain prism's the-top-is-the-bottom-moved
/// invariant, so corners re-solve where the offset lines meet, straight
/// walls come out as exact tilted planes, and a full circular hole's wall
/// is an exact cone.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// profile is not a planar face, the travel is not square to it, a profile
/// edge is neither straight nor a full circle (the curved-wall taper needs a
/// fitted ruling — docs/PARITY.md, offset.sweeps), or the taper collapses a
/// loop over the height.
pub fn make_prism_tapered(
    model: &mut Model,
    profile: &Shape,
    vector: Vector,
    taper: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve;

    if !taper.is_finite() || taper.abs() <= tol.angular() {
        return make_prism(model, profile, vector, tol);
    }
    if taper.abs() >= core::f64::consts::FRAC_PI_2 - tol.angular() {
        ogeom_bail!(Construction, "a taper of {taper} flattens the prism");
    }
    let travel = vector.magnitude();
    if !travel.is_finite() || travel <= tol.confusion() {
        ogeom_bail!(Construction, "a sweep along {vector:?} goes nowhere");
    }
    if model.kind_of(profile)? != ShapeType::Face {
        ogeom_bail!(
            Construction,
            "a tapered prism encloses volume; sweep a planar face"
        );
    }
    let Some(plane) = crate::build::find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a tapered prism sweeps a planar face");
    };
    let mut normal = plane.normal().vector();
    if normal.cross(vector / travel).magnitude() > tol.angular().max(1e-9) {
        ogeom_bail!(
            Construction,
            "an oblique tapered sweep is ambiguous about its own sections; \
             the travel must run square to the profile"
        );
    }
    // Everything below leans on the travel direction, not the face's own
    // normal sign.
    if normal.dot(vector) < 0.0 {
        normal = -normal;
    }
    let up = Direction::new(normal, tol)?;
    let spread = travel * taper.tan();

    let mut history = History::new();
    let mut faces: Vec<Shape> = Vec::new();
    let mut near_wires: Vec<Shape> = Vec::new();
    let mut far_wires: Vec<Shape> = Vec::new();
    for wire in model.children_of(profile)? {
        let edges = model.ordered_children_of(&wire)?;
        // A loop of one closed circle tapers as a cone; a loop of straight
        // edges tapers as tilted planes with re-mitred corners.
        let lone_circle = edges.len() == 1 && {
            let (curve, _) = edge_geometry(model, &edges[0])?;
            matches!(curve, Curve::Circle(_))
        };
        if lone_circle {
            let (curve, range) = edge_geometry(model, &edges[0])?;
            let Curve::Circle(c) = curve else {
                unreachable!("just matched")
            };
            let circle = c.circle();
            // Which way is away from the material, by measurement: probe a
            // little outside the ring and ask the face. Windings are not to
            // be trusted — a hole's wire may come wound either way.
            let start = c.point_at(range.0, tol)?;
            let radial = (start - circle.centre()) / circle.radius();
            let sigma = away_sign(model, profile, start, radial, circle.radius(), tol)?;
            let far_radius = sigma.mul_add(spread, circle.radius());
            if far_radius <= tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the taper collapses a circular loop of radius {} over \
                     this height",
                    circle.radius()
                );
            }
            let near_frame = Frame::new(circle.centre(), up, circle.frame().x(), tol)?;
            let far_frame = Frame::new(circle.centre() + vector, up, circle.frame().x(), tol)?;
            let near_curve: Curve =
                ogeom_geom::CircleCurve::new(Circle::new(near_frame, circle.radius(), tol)?).into();
            let near_domain = near_curve.domain();
            let near = crate::build::make_edge(model, near_curve, near_domain, tol)?.shape;
            let far_curve: Curve =
                ogeom_geom::CircleCurve::new(Circle::new(far_frame, far_radius, tol)?).into();
            let far_domain = far_curve.domain();
            let far = crate::build::make_edge(model, far_curve, far_domain, tol)?.shape;
            // The wall cone: reference radius at the near plane, the radius
            // running to the far one; the surface's own normal points away
            // from the axis, which is outward exactly when the material is
            // inside the ring.
            let slope = (far_radius - circle.radius()) / travel;
            let cone = ogeom_math::Cone::new(near_frame, circle.radius(), slope.atan(), tol)?;
            let pad = travel * 0.1;
            let surface: ogeom_geom::SurfaceGeometry =
                ogeom_geom::ConeSurface::new(cone, (-pad, travel + pad))?.into();
            let band = crate::build::make_revolution_band(model, &surface, &near, &far, tol)?;
            let wall = if sigma > 0.0 { band } else { band.reversed() };
            model.set_derived(&wall, std::slice::from_ref(&edges[0]), roles::SWEEP_SIDE)?;
            history.generate(&edges[0], wall.clone());
            faces.push(wall);
            near_wires.push(make_wire(model, std::slice::from_ref(&near), tol)?.shape);
            far_wires.push(make_wire(model, std::slice::from_ref(&far), tol)?.shape);
            continue;
        }

        // Straight loops: offset every line away from the material and
        // re-mitre the corners in the far plane.
        let mut corners_near: Vec<Point> = Vec::new();
        let mut aways: Vec<Vector> = Vec::new();
        let mut dirs: Vec<Vector> = Vec::new();
        for edge in &edges {
            let (curve, range) = edge_geometry(model, edge)?;
            let Curve::Line(_) = curve else {
                ogeom_bail!(
                    Construction,
                    "a tapered wall over an edge that is neither straight nor \
                     a full circle needs a fitted ruling — docs/PARITY.md, \
                     offset.sweeps"
                );
            };
            let reversed = edge.orientation() == Orientation::Reversed;
            let (t0, t1) = if reversed {
                (range.1, range.0)
            } else {
                (range.0, range.1)
            };
            let from = curve.point_at(t0, tol)?;
            let to = curve.point_at(t1, tol)?;
            let dir = (to - from) / from.distance(to);
            corners_near.push(from);
            dirs.push(dir);
            aways.push(dir.cross(up.vector()));
        }
        let count = corners_near.len();
        if count < 3 {
            ogeom_bail!(Construction, "a straight loop needs at least three edges");
        }
        // The loop's material side, measured once and applied to every edge:
        // the wire's own winding is not to be trusted for holes.
        let flip = {
            let mid = corners_near[0] + dirs[0] * (corners_near[0].distance(corners_near[1]) / 2.0);
            let scale = corners_near[0].distance(corners_near[1]);
            away_sign(model, profile, mid, aways[0], scale, tol)?
        };
        if flip < 0.0 {
            for a in &mut aways {
                *a = -*a;
            }
        }
        // Far corners: each is where the two neighbouring offset lines meet,
        // in the far plane.
        let mut corners_far: Vec<Point> = Vec::with_capacity(count);
        for i in 0..count {
            let prev = (i + count - 1) % count;
            let (d0, d1) = (dirs[prev], dirs[i]);
            let (a0, a1) = (aways[prev], aways[i]);
            let p0 = corners_near[i] + a0 * spread + vector;
            let p1 = corners_near[i] + a1 * spread + vector;
            let cross = d0.cross(d1);
            let m = cross.magnitude();
            let far = if m <= tol.angular() {
                // Collinear neighbours offset to the same line.
                p1
            } else {
                // Solve p0 + s d0 = p1 + t d1 in the loop's own plane.
                let w = p1 - p0;
                let s = w.cross(d1).dot(cross) / (m * m);
                p0 + d0 * s
            };
            corners_far.push(far);
        }
        for (i, far) in corners_far.iter().enumerate() {
            let next = corners_far[(i + 1) % count];
            let d = next - *far;
            if d.magnitude() <= tol.confusion() || d.dot(dirs[i]) <= 0.0 {
                ogeom_bail!(
                    Construction,
                    "the taper collapses the profile's loop over this height"
                );
            }
        }

        let near_vertices: Vec<Shape> = corners_near
            .iter()
            .map(|p| crate::build::make_vertex(model, *p).shape)
            .collect();
        let far_vertices: Vec<Shape> = corners_far
            .iter()
            .map(|p| crate::build::make_vertex(model, *p).shape)
            .collect();
        let segment =
            |model: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| -> OgeomResult<Shape> {
                let line = ogeom_geom::LineCurve::segment(from.1, to.1, tol)?;
                let curve: Curve = line.into();
                let domain = curve.domain();
                Ok(crate::build::make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape)
            };
        let mut near_edges = Vec::with_capacity(count);
        let mut far_edges = Vec::with_capacity(count);
        let mut rails = Vec::with_capacity(count);
        for i in 0..count {
            let next = (i + 1) % count;
            near_edges.push(segment(
                model,
                (&near_vertices[i], corners_near[i]),
                (&near_vertices[next], corners_near[next]),
            )?);
            far_edges.push(segment(
                model,
                (&far_vertices[i], corners_far[i]),
                (&far_vertices[next], corners_far[next]),
            )?);
            rails.push(segment(
                model,
                (&near_vertices[i], corners_near[i]),
                (&far_vertices[i], corners_far[i]),
            )?);
        }
        for (i, edge) in edges.iter().enumerate() {
            let next = (i + 1) % count;
            // The wall's own plane: through the near edge, leaning with the
            // rails. Its normal is the in-plane away tilted by the taper,
            // which faces off the material by construction.
            let outward = {
                let lean = aways[i] * taper.cos() - up.vector() * taper.sin();
                Direction::new(lean, tol)?
            };
            let wall_plane = ogeom_math::Plane::through(corners_near[i], outward);
            let mut reach = travel + 1.0_f64;
            for p in [
                corners_near[i],
                corners_near[next],
                corners_far[i],
                corners_far[next],
            ] {
                reach = reach.max(p.distance(corners_near[i]) * 2.0);
            }
            let surface: ogeom_geom::SurfaceGeometry =
                ogeom_geom::PlaneSurface::over(wall_plane, (-reach, reach), (-reach, reach))?
                    .into();
            let wall = crate::build::make_face_with_pcurves(
                model,
                surface,
                &[vec![
                    near_edges[i].clone(),
                    rails[next].clone(),
                    far_edges[i].reversed(),
                    rails[i].reversed(),
                ]],
                tol,
            )?
            .shape;
            model.set_derived(&wall, std::slice::from_ref(edge), roles::SWEEP_SIDE)?;
            history.generate(edge, wall.clone());
            faces.push(wall);
        }
        near_wires.push(make_wire(model, &near_edges, tol)?.shape);
        far_wires.push(make_wire(model, &far_edges, tol)?.shape);
    }

    // Caps: the near one facing backwards, both rebuilt over the fresh rings
    // so every wall shares its vertices.
    let near_surface: ogeom_geom::SurfaceGeometry = {
        let base = ogeom_math::Plane::through(plane.origin(), up);
        ogeom_geom::PlaneSurface::over(base, (-1e6, 1e6), (-1e6, 1e6))?.into()
    };
    let far_surface: ogeom_geom::SurfaceGeometry = {
        let lifted = ogeom_math::Plane::through(plane.origin() + vector, up);
        ogeom_geom::PlaneSurface::over(lifted, (-1e6, 1e6), (-1e6, 1e6))?.into()
    };
    let bottom = crate::build::make_face(model, near_surface, &near_wires, tol)?
        .shape
        .reversed();
    let top = crate::build::make_face(model, far_surface, &far_wires, tol)?.shape;
    attach_cap_pcurves(model, &bottom, tol)?;
    attach_cap_pcurves(model, &top, tol)?;
    model.set_derived(&bottom, std::slice::from_ref(profile), roles::SWEEP_BOTTOM)?;
    model.set_derived(&top, std::slice::from_ref(profile), roles::SWEEP_TOP)?;
    history.generate(profile, top.clone());
    faces.push(bottom);
    faces.push(top);

    let sewn = crate::sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !crate::build::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the tapered prism did not close");
    }
    let solid = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?.shape;
    history.generate(profile, solid.clone());
    Ok(Built::new(solid, history))
}

/// The sign that points `candidate` away from the face's material at `at`,
/// probed against the face's own trim at a few scales of `scale`.
fn away_sign(
    model: &Model,
    face: &Shape,
    at: Point,
    candidate: Vector,
    scale: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    for factor in [1e-3, 1e-2, 5e-2] {
        let eps = scale * factor;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        for sign in [1.0_f64, -1.0] {
            let probe = at + candidate * (sign * eps);
            if crate::classify_on_face(model, face, probe, deflection, tol)?
                == crate::Containment::In
            {
                // Material on this side: away is the other one.
                return Ok(-sign);
            }
        }
    }
    ogeom_bail!(
        Construction,
        "cannot read which side of the profile the material is on"
    )
}

/// An edge's 3D curve and range, cloned out of the model.
fn edge_geometry(model: &Model, edge: &Shape) -> OgeomResult<(ogeom_geom::Curve, (f64, f64))> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Construction, "an edge holds no data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "an edge has no curve");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    Ok((geometry.clone(), *range))
}

/// Attach exact pcurves to every edge of a planar cap that lacks one.
fn attach_cap_pcurves(model: &mut Model, cap: &Shape, tol: Tolerances) -> OgeomResult<()> {
    let cap_id = {
        let Some(node) = model.node(cap) else {
            ogeom_bail!(Dangling, "the cap just built is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "the cap holds no face data");
        };
        data.surface
    };
    let Some(surface) = model.geometry().surface(cap_id).cloned() else {
        ogeom_bail!(Dangling, "the cap's surface is not in this model");
    };
    for edge in ogeom_topo::explore(model, cap, ogeom_topo::Filter::OfType(ShapeType::Edge))? {
        let (curve, range) = edge_geometry(model, &edge)?;
        let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) else {
            ogeom_bail!(Construction, "a cap edge has no closed-form pcurve");
        };
        crate::build::attach_pcurve(model, &edge, pcurve, cap_id, Location::identity(), range)?;
    }
    Ok(())
}

/// A face swept into a solid.
fn prism_over_face(
    model: &mut Model,
    rails: &mut Rails,
    face: &Shape,
    displacement: &Location,
    vector: Vector,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // Which side of the profile the material lands on is decided by the sweep,
    // not by which way the profile was handed over. A profile facing against
    // the sweep does not describe a different solid — it describes the same one
    // from the other side — so it is turned round here and everything below
    // proceeds as if it had faced along all along.
    //
    // Left unturned, both end caps present the wrong side: the mesh still
    // closes and the shell is still closed, so nothing topological notices, and
    // the volume comes back short by twice the caps' contribution.
    let (_, normal) = crate::measure::face_normal(model, face, tol)?;
    let travel = vector.magnitude();
    let along = normal.dot(vector) / travel;
    if along.abs() <= tol.angular() {
        ogeom_bail!(
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
) -> OgeomResult<(Vec<Shape>, History)> {
    let mut faces = Vec::new();
    let mut history = History::new();
    for edge in model.ordered_children_of(wire)? {
        let (face, edge_history) = prism_over_edge(model, rails, &edge, displacement, vector, tol)?;
        history = history.then(&edge_history);
        faces.push(face);
    }
    if faces.is_empty() {
        ogeom_bail!(Construction, "a wire with no edges sweeps out nothing");
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
) -> OgeomResult<(Shape, History)> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(
            Construction,
            "an edge with no curve in space has no shape to sweep; a degenerate \
             edge sweeps out nothing and has to be handled by its face, not here"
        );
    };
    let Some(geometry) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };

    // The surface is built from the edge's curve *where the edge actually is*.
    // Built from the stored curve instead, every side wall of a profile that
    // has been placed lands back at the placement's origin, while its two ends
    // — which are the profile itself, at its own location — land correctly, so
    // the mesh comes apart along every lateral face at once.
    //
    // A placement may carry a uniform scale, and a scale rescales a curve's
    // parameter with it. The edge's range is in the stored curve's parameter
    // and the surface's `u` is in the placed one's, so the range is carried
    // across by where it sits in the domain rather than by its value.
    let placement = edge.transform(model.datums())?;
    let stored = geometry.domain();
    let geometry = geometry.transformed(&placement, tol)?;
    let placed = geometry.domain();
    let (lo, hi) = (
        rescale(range.0, stored, placed),
        rescale(range.1, stored, placed),
    );
    let travel = vector.magnitude();
    let direction = ogeom_math::Direction::new(vector, tol)?;
    // A straight profile edge sweeps a plane, and the plane is built so its
    // chart *is* the extrusion's — origin on the line, x along it, y along
    // the travel — so every pcurve below serves either surface unchanged.
    // Naming the plane it actually made is what lets the boolean's
    // same-domain resolution meet a prism wall as the plane it is.
    let canonical: Option<ogeom_geom::SurfaceGeometry> = if let ogeom_geom::Curve::Line(line) =
        &geometry
        && line.axis().direction.vector().dot(vector).abs() <= tol.angular() * travel
    {
        // Chart identity needs the travel perpendicular to the line; an
        // oblique sweep's plane exists too, but with a sheared chart the
        // pcurves below would no longer describe.
        let axis = line.axis();
        let normal = ogeom_math::Direction::new(axis.direction.vector().cross(vector), tol).ok();
        normal
            .map(|n| -> OgeomResult<ogeom_geom::SurfaceGeometry> {
                let frame = ogeom_math::Frame::new(axis.location, n, axis.direction, tol)?;
                let plane = ogeom_math::Plane::new(frame);
                let margin = (hi - lo).abs().max(travel) * 0.1 + 1.0;
                Ok(ogeom_geom::PlaneSurface::over(
                    plane,
                    (lo.min(hi) - margin, lo.max(hi) + margin),
                    (-margin, travel + margin),
                )?
                .into())
            })
            .transpose()?
    } else if let ogeom_geom::Curve::Circle(c) = &geometry
        && !c.is_reversed()
        && c.circle().frame().z().vector().dot(direction.vector()) >= 1.0 - tol.angular()
    {
        // A circular profile edge swept along its own axis is a cylinder,
        // and on the circle's own frame the chart *is* the extrusion's —
        // u the circle's angle, v the travel — so the pcurves below serve
        // either surface unchanged, and the boolean's same-domain
        // resolution meets a prism wall as the cylinder it is.
        let circle = c.circle();
        let margin = travel * 0.1 + 1.0;
        Some(
            ogeom_geom::CylinderSurface::new(
                ogeom_math::Cylinder::new(circle.frame(), circle.radius(), tol)?,
                (-margin, travel + margin),
            )?
            .into(),
        )
    } else {
        None
    };
    let surface = model.geometry_mut().add_surface(match canonical {
        Some(exact) => exact,
        None => ExtrusionSurface::new(geometry, direction, travel)?.into(),
    });

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

/// The turn a revolution makes, carried down to every entity it builds.
struct Turn {
    /// What it turns about.
    axis: Axis,
    /// How far, in radians.
    angle: f64,
    /// Whether the turn closes on itself, so the two ends are one seam rather
    /// than two faces.
    full: bool,
    /// Where the far end sits. The identity for a full turn, because there is
    /// no far end — it is the near end again.
    displacement: Location,
}

/// Revolve a shape about `axis`, through `angle` radians.
///
/// A face becomes a solid, a wire becomes a shell, an edge becomes a face —
/// the same rule as [`make_prism`], turning instead of travelling. A full turn
/// closes on itself: its two ends are the *same* profile, meeting along a seam,
/// which is the topology [`make_cylinder`](crate::make_cylinder) produces for
/// the identical solid. A partial turn has two distinct ends and caps them.
///
/// The seam is not an implementation detail to be avoided. Building a full turn
/// as two half-turns would sidestep it and give a different face count for the
/// same solid, and a boolean later has to split along seams like any other
/// edge.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the angle is
/// not in `(0, 2pi]`, if the shape is of a kind that cannot be swept, if an
/// edge of the profile has no 3D curve, if the profile meets the axis anywhere
/// but at its ends — it would sweep through itself — or if the profile lies in
/// a surface the turn runs along, which encloses no volume.
pub fn make_revolution(
    model: &mut Model,
    profile: &Shape,
    axis: Axis,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !angle.is_finite() || angle <= tol.angular() || angle > TAU + tol.angular() {
        ogeom_bail!(
            Construction,
            "a revolution turns through (0, 2pi]; {angle} does not"
        );
    }
    let angle = angle.min(TAU);
    let full = TAU - angle <= tol.angular();
    model.begin_operation();

    let turn = Turn {
        axis,
        angle,
        full,
        // One datum for the whole turn, so every entity at the far end shares a
        // single placement — and none at all for a full turn, whose far end is
        // its near end.
        displacement: if full {
            Location::identity()
        } else {
            Location::of(model.add_datum(Transform::rotation(axis, angle)))
        },
    };

    let rails = &mut Rails::new();
    match model.kind_of(profile)? {
        ShapeType::Face => revolution_over_face(model, rails, profile, &turn, tol),
        ShapeType::Wire => {
            let (faces, history) = revolution_over_wire(model, rails, profile, &turn, tol)?;
            let shell = make_shell(model, &faces)?.shape;
            Ok(Built::new(shell, history))
        }
        ShapeType::Edge => {
            let (face, history) = revolution_over_edge(model, rails, profile, &turn, tol)?;
            let Some(face) = face else {
                ogeom_bail!(
                    Construction,
                    "an edge lying along the axis turns onto itself and sweeps \
                     out no face"
                );
            };
            Ok(Built::new(face, history))
        }
        other => ogeom_bail!(
            Construction,
            "a {other:?} cannot be revolved into anything; revolve an edge, a \
             wire or a face"
        ),
    }
}

/// A face revolved into a solid.
fn revolution_over_face(
    model: &mut Model,
    rails: &mut Rails,
    face: &Shape,
    turn: &Turn,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // The same question the prism asks, with the sweep direction read off the
    // turn: at the profile, revolving moves it along the tangent to its circle
    // about the axis, so that tangent is what its normal is compared against.
    let (point, normal) = crate::measure::face_normal(model, face, tol)?;
    let tangent = turn
        .axis
        .direction
        .cross_with(point - turn.axis.project(point));
    let reach = tangent.magnitude();
    if reach <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the profile sits on the axis, so revolving it sweeps out nothing"
        );
    }
    let along = normal.dot(tangent) / reach;
    if along.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the turn runs along the profile's own surface, so it encloses no \
             volume; a face revolved within its own plane is not a solid"
        );
    }
    // The sweep's material side follows the profile wire's own walk, so the
    // walk is normalized to one hand — measured from the traversal itself,
    // as the loop's area vector against the sweep tangent. The face's
    // stated normal cannot answer this: it speaks the carrier's chart,
    // and one loop reads as either hand depending on which way the chart
    // was laid down.
    let hand = {
        let mut area = ogeom_math::Vector::ZERO;
        let wires = model.ordered_children_of(face)?;
        let Some(outer) = wires.first() else {
            ogeom_bail!(Construction, "the profile has no boundary to revolve");
        };
        let mut walk: Vec<Point> = Vec::new();
        for edge in model.ordered_children_of(outer)? {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "a profile edge holds no data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "a profile edge has no curve");
            };
            let Some(stored) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            let placed = stored
                .clone()
                .transformed(&edge.transform(model.datums())?, tol)?;
            let flipped = edge.orientation() == Orientation::Reversed;
            for k in 0..8 {
                let f = f64::from(k) / 8.0;
                let f = if flipped { 1.0 - f } else { f };
                let t = (range.1 - range.0).mul_add(f, range.0);
                walk.push(placed.point_at(t, tol)?);
            }
        }
        for k in 0..walk.len() {
            let (a, b) = (walk[k], walk[(k + 1) % walk.len()]);
            area += (a - point).cross(b - point);
        }
        area.dot(tangent)
    };
    let profile = if hand < 0.0 {
        face.reversed()
    } else {
        face.clone()
    };

    let mut history = History::new();
    let mut faces = Vec::new();
    for wire in model.children_of(&profile)? {
        let (sides, wire_history) = revolution_over_wire(model, rails, &wire, turn, tol)?;
        history = history.then(&wire_history);
        faces.extend(sides);
    }

    if turn.full {
        // A full turn has no ends. The profile is an interior cross-section of
        // the result, not a face of it — so it is *deleted*, while still
        // generating everything its edges and vertices swept out. Reporting it
        // as surviving would leave a reference resolving to a face that is not
        // on the solid.
        history.delete(face);
    } else {
        let bottom = profile.reversed();
        let top = profile.moved(&turn.displacement);
        model.set_derived(&bottom, std::slice::from_ref(face), roles::SWEEP_BOTTOM)?;
        model.set_derived(&top, std::slice::from_ref(face), roles::SWEEP_TOP)?;
        history.generate(face, top.clone());
        faces.push(bottom);
        faces.push(top);
    }

    let shell = make_shell(model, &faces)?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    history.generate(face, solid.clone());
    Ok(Built::new(solid, history))
}

/// Every face a wire revolves out.
fn revolution_over_wire(
    model: &mut Model,
    rails: &mut Rails,
    wire: &Shape,
    turn: &Turn,
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, History)> {
    let mut faces = Vec::new();
    let mut history = History::new();
    for edge in model.ordered_children_of(wire)? {
        let (face, edge_history) = revolution_over_edge(model, rails, &edge, turn, tol)?;
        history = history.then(&edge_history);
        // An edge lying along the axis turns onto itself. It contributes no
        // face, which is what makes a rectangle with one side on the axis
        // revolve into a cylinder — three faces — rather than into a cylinder
        // with a fourth face of no area down its middle.
        faces.extend(face);
    }
    if faces.is_empty() {
        ogeom_bail!(
            Construction,
            "the whole wire lies along the axis, so it revolves out nothing"
        );
    }
    Ok((faces, history))
}

/// One edge revolved into one face.
///
/// The rectangle is transposed from the prism's: a revolution's `u` is the
/// angle turned and its `v` is the generating curve's own parameter, so the
/// profile edge runs up the *sides* of the parameter rectangle and the circles
/// its endpoints sweep run across the top and bottom.
fn revolution_over_edge(
    model: &mut Model,
    rails: &mut Rails,
    edge: &Shape,
    turn: &Turn,
    tol: Tolerances,
) -> OgeomResult<(Option<Shape>, History)> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(
            Construction,
            "an edge with no curve in space has no shape to revolve; a \
             degenerate edge sweeps out nothing and has to be handled by its \
             face, not here"
        );
    };
    let Some(geometry) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };

    // As for the prism: the surface is the curve where the edge actually is,
    // and the range is carried across the placement's effect on the parameter.
    let placement = edge.transform(model.datums())?;
    let stored = geometry.domain();
    let geometry = geometry.transformed(&placement, tol)?;
    let placed = geometry.domain();
    let (lo, hi) = (
        rescale(range.0, stored, placed),
        rescale(range.1, stored, placed),
    );
    if axis_relation(&geometry, (lo, hi), turn.axis, tol)? == AxisRelation::On {
        return Ok((None, History::new()));
    }

    // A profile line perpendicular to the axis stays at one height as it
    // turns: the face it sweeps is a region of a *plane*, and naming that
    // plane — rather than dressing it as a revolution with a polar chart, a
    // seam and a degenerate centre edge — is what lets a revolved rectangle
    // match `make_cylinder` face for face and edge for edge.
    if let ogeom_geom::Curve::Line(line) = &geometry
        && line
            .axis()
            .direction
            .vector()
            .dot(turn.axis.direction.vector())
            .abs()
            <= tol.angular()
    {
        return flat_revolution(model, rails, edge, &geometry, (lo, hi), turn, tol);
    }

    let surface = model
        .geometry_mut()
        .add_surface(ogeom_geom::RevolutionSurface::new(geometry, turn.axis, turn.angle)?.into());

    let reversed = edge.orientation() == Orientation::Reversed;
    let (v_start, v_end) = if reversed { (hi, lo) } else { (lo, hi) };
    let (near, far) = (0.0, turn.angle);

    // The circles the two ends sweep. A full turn brings each back to where it
    // started, so it is one closed edge; a partial turn leaves an arc between
    // the endpoint and its rotated copy.
    let start_rail = revolved_rail(model, rails, edge, turn, false, tol)?;
    let end_rail = revolved_rail(model, rails, edge, turn, true, tol)?;

    if start_rail.is_same(&end_rail) {
        // A closed profile edge — revolving a circle makes a torus — returns to
        // one vertex, so its two rails are one edge bounding the face across
        // both the bottom and the top of the parameter rectangle. That is a
        // seam in `v`, and it needs both pcurves for the same reason a seam in
        // `u` does.
        seam_pcurves(
            model,
            &start_rail,
            surface,
            ((near, v_start), (far, v_start)),
            ((near, v_end), (far, v_end)),
            tol,
        )?;
    } else {
        pcurve(
            model,
            &start_rail,
            surface,
            (near, v_start),
            (far, v_start),
            tol,
        )?;
        pcurve(model, &end_rail, surface, (near, v_end), (far, v_end), tol)?;
    }

    // The profile edge itself runs up the sides. Both pcurves follow the
    // curve's own parameterization — `lo` to `hi` — because a pcurve describes
    // the edge and not the wire's walk of it.
    let displaced = edge.moved(&turn.displacement);
    if turn.full {
        // The two sides are one edge appearing twice, at `u = 0` and at
        // `u = 2pi`. Which pcurve applies is decided by the occurrence's
        // orientation, and the ring below puts the walk that goes *up* the far
        // side on whichever occurrence carries the edge's own direction.
        let (forward, reversed_side) = if reversed { (near, far) } else { (far, near) };
        seam_pcurves(
            model,
            edge,
            surface,
            ((forward, lo), (forward, hi)),
            ((reversed_side, lo), (reversed_side, hi)),
            tol,
        )?;
    } else {
        pcurve(model, edge, surface, (near, lo), (near, hi), tol)?;
        pcurve(model, &displaced, surface, (far, lo), (far, hi), tol)?;
    }

    // Round the rectangle: across the bottom in the direction of the turn, up
    // the far side of the profile, back across the top, down the near side.
    let ring = [
        start_rail.clone(),
        displaced.clone(),
        end_rail.reversed(),
        edge.reversed(),
    ];
    let boundary = make_wire(model, &ring, tol)?.shape;
    let built = make_face_on(model, surface, std::slice::from_ref(&boundary), tol)?.shape;

    // A revolution's normal is the *turn's* tangent crossed with the curve's,
    // because the angle is `u` and the curve is `v`. The prism's is the other
    // way round, so the two disagree by a sign for the same walk — and an
    // occurrence the wire walks forwards is the one that has to be reversed
    // here. It is not the surface that decides which side is material; it is
    // which way the profile's wire goes round.
    let face = if reversed { built } else { built.reversed() };
    model.set_derived(&face, std::slice::from_ref(edge), roles::SWEEP_SIDE)?;

    let mut history = History::new();
    // Both, not either — as for the prism. The edge makes the lateral face
    // *and* survives: on a partial turn as the far side, and on a full turn as
    // the seam, which is the same edge occurring twice on one face rather than
    // an edge that ceased to exist.
    history.generate(edge, face.clone());
    if !turn.full {
        history.generate(edge, displaced);
    }
    Ok((Some(face), history))
}

/// The planar face a radial profile line sweeps: a disc, an annulus or a
/// pie, on the plane it actually turns in.
fn flat_revolution(
    model: &mut Model,
    rails: &mut Rails,
    edge: &Shape,
    geometry: &ogeom_geom::Curve,
    range: (f64, f64),
    turn: &Turn,
    tol: Tolerances,
) -> OgeomResult<(Option<Shape>, History)> {
    use ogeom_geom::Curve3d as _;
    let (lo, hi) = range;
    let axis_dir = turn.axis.direction;
    let at = geometry.point_at(lo, tol)?;
    let height = (at - turn.axis.location).dot(axis_dir.vector());
    let foot = turn.axis.location + axis_dir.vector() * height;
    let radius_of = |p: ogeom_math::Point| (p - foot).magnitude();
    let far = geometry.point_at(hi, tol)?;
    let reach = radius_of(at).max(radius_of(far)) * 2.0 + 1.0;
    let plane = ogeom_math::Plane::new(ogeom_math::Frame::about(foot, axis_dir));
    let plane_surface: ogeom_geom::SurfaceGeometry =
        ogeom_geom::PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
    let surface = model.geometry_mut().add_surface(plane_surface.clone());

    let start_rail = revolved_rail(model, rails, edge, turn, false, tol)?;
    let end_rail = revolved_rail(model, rails, edge, turn, true, tol)?;
    let is_degenerate = |model: &Model, e: &Shape| -> bool {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|d| d.curve3d().is_none())
    };

    // Exact Cartesian pcurves for whichever edges bound the face; the
    // degenerate centre of the old polar chart simply has no place here.
    let attach = |model: &mut Model, occurrence: &Shape| -> OgeomResult<()> {
        let Some(data) = model.node(occurrence).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "a flat revolution edge holds no data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "a flat revolution edge has no curve");
        };
        let (curve, range) = (*curve, *range);
        let Some(stored) = model.geometry().curve(curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let placed = stored
            .clone()
            .transformed(&occurrence.transform(model.datums())?, tol)?;
        let Some(pc) = ogeom_intersect::exact_pcurve_of(&placed, &plane_surface, tol) else {
            ogeom_bail!(
                Construction,
                "a flat revolution edge has no closed-form pcurve on its plane"
            );
        };
        crate::build::attach_pcurve(
            model,
            occurrence,
            pc,
            surface,
            occurrence.location().clone(),
            range,
        )
    };

    let reversed = edge.orientation() == Orientation::Reversed;
    let built = if turn.full {
        let mut wires = Vec::new();
        for rail in [&start_rail, &end_rail] {
            if is_degenerate(model, rail) {
                continue;
            }
            attach(model, rail)?;
            wires.push(make_wire(model, std::slice::from_ref(rail), tol)?.shape);
        }
        if wires.is_empty() {
            ogeom_bail!(Construction, "a flat revolution swept out no boundary");
        }
        make_face_on(model, surface, &wires, tol)?.shape
    } else {
        let displaced = edge.moved(&turn.displacement);
        let mut ring: Vec<Shape> = Vec::new();
        if !is_degenerate(model, &start_rail) {
            attach(model, &start_rail)?;
            ring.push(start_rail.clone());
        }
        attach(model, &displaced)?;
        ring.push(displaced.clone());
        if !is_degenerate(model, &end_rail) {
            attach(model, &end_rail)?;
            ring.push(end_rail.reversed());
        }
        attach(model, edge)?;
        ring.push(edge.reversed());
        let boundary = make_wire(model, &ring, tol)?.shape;
        make_face_on(model, surface, std::slice::from_ref(&boundary), tol)?.shape
    };

    // The material side is the profile wire's business, as for every sweep;
    // the plane's own normal relates to the revolution's by the sign of the
    // line's outward sense.
    let outward_sense = {
        let radial = if radius_of(far) >= radius_of(at) {
            far - foot
        } else {
            at - foot
        };
        let d = (far - at).dot(radial);
        d < 0.0
    };
    let face = if reversed != outward_sense {
        built.reversed()
    } else {
        built
    };
    model.set_derived(&face, std::slice::from_ref(edge), roles::SWEEP_SIDE)?;

    let mut history = History::new();
    history.generate(edge, face.clone());
    if !turn.full {
        history.generate(edge, edge.moved(&turn.displacement));
    }
    Ok((Some(face), history))
}

/// The circle or arc one endpoint of the profile sweeps out.
///
/// Shared between the two faces that meet along it, exactly as the prism's
/// rails are — building one per face would leave every rail used once and the
/// solid open along every corner.
///
/// An endpoint *on* the axis sweeps out nothing, and gets a degenerate edge: it
/// still bounds the face across its side of the parameter rectangle, and
/// leaving it out would leave the boundary open there with nothing to trim to.
fn revolved_rail(
    model: &mut Model,
    rails: &mut Rails,
    edge: &Shape,
    turn: &Turn,
    at_end: bool,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let Some((start, end)) = crate::build::edge_vertices(model, edge)? else {
        ogeom_bail!(
            Construction,
            "an unbounded edge has no endpoints to sweep into rails"
        );
    };
    let base = if at_end { end } else { start };
    if let Some(existing) = rails.get(&base.node()) {
        return Ok(existing.clone());
    }

    let Some(node) = model.node(&base) else {
        ogeom_bail!(Dangling, "vertex is not in this model");
    };
    let Some(data) = node.data().as_vertex() else {
        ogeom_bail!(Construction, "vertex node holds no point");
    };
    let from = base.transform(model.datums())?.apply(data.point);

    // A full turn brings the endpoint back to itself, so the rail is one closed
    // edge named twice by the same vertex — which is what keeps "walk to the
    // end" meaningful all the way round.
    let raised = if turn.full {
        base.clone()
    } else {
        base.moved(&turn.displacement)
    };

    let radius = from - turn.axis.project(from);
    let built = if radius.magnitude() <= tol.confusion() {
        let mut data = ogeom_topo::EdgeData::new();
        data.degenerate = true;
        model.add_edge(data, &[base.clone(), raised])?
    } else {
        // The circle's `x` points from the axis out to the endpoint, so its
        // angle parameter *is* the revolution's `u`: at zero it lands on the
        // endpoint exactly, rather than merely nearby.
        let frame = Frame::new(
            turn.axis.project(from),
            turn.axis.direction,
            Direction::new(radius, tol)?,
            tol,
        )?;
        let circle = Circle::new(frame, radius.magnitude(), tol)?;
        crate::build::make_edge_between(
            model,
            ogeom_geom::CircleCurve::new(circle).into(),
            (0.0, turn.angle),
            &base,
            &raised,
            tol,
        )?
        .shape
    };

    model.set_derived(&built, std::slice::from_ref(&base), roles::SWEEP_RAIL)?;
    rails.insert(base.node(), built.clone());
    Ok(built)
}

/// How many places along an edge are checked against the axis.
const AXIS_SAMPLES: usize = 32;

/// Where an edge stands in relation to the axis it is to be turned about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AxisRelation {
    /// Clear of it, except possibly at its ends — which sweep out poles.
    Clear,
    /// Lying along it, so it sweeps out nothing at all.
    On,
}

/// Where an edge stands in relation to the axis, refusing the case that cannot
/// be built.
///
/// An edge *crossing* the axis sweeps out a surface that passes through itself,
/// and the solid built on it is wrong in a way nothing downstream can detect:
/// its volume is finite and plausible and counts part of space twice. That is
/// refused. An edge lying *along* the axis sweeps out nothing — the inner side
/// of a rectangle that revolves into a cylinder — and gets no face. An edge
/// merely touching the axis at an end sweeps out a pole, which is ordinary.
///
/// Solved, not sampled. An earlier version sampled the radial vector at
/// thirty-two places and watched for its direction reversing, with the note
/// that deciding exactly where a curve meets a line is the intersector's work
/// and the intersector did not exist yet. It does now, and the decision is
/// exact in two layers: the intersector names every point where the curve
/// meets the axis within tolerance, and the extrema machinery names every
/// stationary nearest approach — which is what catches the case the sampling
/// never could, a profile grazing the axis *tangentially* between samples,
/// where the radial direction never reverses and no sample lands on the
/// touch.
fn axis_relation(
    curve: &ogeom_geom::Curve,
    range: (f64, f64),
    axis: Axis,
    tol: Tolerances,
) -> OgeomResult<AxisRelation> {
    // The cheap layer stays for what it answers exactly: an edge every one of
    // whose samples sits on the axis is a straight edge lying along it.
    let mut points = Vec::with_capacity(AXIS_SAMPLES + 1);
    for i in 0..=AXIS_SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = i as f64 / AXIS_SAMPLES as f64;
        let at = range.0 + (range.1 - range.0) * t;
        points.push(curve.point_at(at, tol)?);
    }
    let on_axis = |p: &ogeom_math::Point| (*p - axis.project(*p)).magnitude() <= tol.confusion();
    if points.iter().all(on_axis) {
        return Ok(AxisRelation::On);
    }

    // The axis as a bounded line covering the edge's whole extent along it,
    // with room to spare either side.
    let along = |p: &ogeom_math::Point| (*p - axis.location).dot(axis.direction.vector());
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for p in &points {
        lo = lo.min(along(p));
        hi = hi.max(along(p));
    }
    let margin = (hi - lo).max(1.0);
    let line: ogeom_geom::Curve =
        ogeom_geom::LineCurve::over(axis, lo - margin, hi + margin)?.into();
    let trimmed: ogeom_geom::Curve = if (range.0, range.1) == curve.domain() {
        curve.clone()
    } else {
        ogeom_geom::TrimmedCurve::new(curve.clone(), range.0, range.1, tol)?.into()
    };

    // A contact at an end is a pole and is ordinary; anywhere else the
    // revolution would pinch or pass through itself.
    let (start, end) = (points[0], points[AXIS_SAMPLES]);
    let interior = |p: ogeom_math::Point| {
        p.distance(start) > tol.confusion() && p.distance(end) > tol.confusion()
    };

    let hits = ogeom_intersect::intersect_curves(
        &trimmed,
        &line,
        ogeom_intersect::CurveCurveOptions::default(),
        tol,
    )?;
    if !hits.overlaps.is_empty() {
        // Lying along the axis in part but not whole: the part that veers off
        // still meets the axis away from the ends. Whole-edge overlap was the
        // all-samples case above.
        ogeom_bail!(
            Construction,
            "the profile touches the axis away from its ends; revolving it \
             would sweep a surface through itself. Split the profile where it \
             meets the axis"
        );
    }
    if hits.crossings.iter().any(|c| interior(c.point)) {
        ogeom_bail!(
            Construction,
            "the profile passes through the axis away from its ends; \
             revolving it would sweep a surface through itself. Split the \
             profile where it meets the axis"
        );
    }

    // The tangential graze: no crossing, but a stationary nearest approach
    // reaching the axis at an interior parameter.
    let near = ogeom_intersect::extrema_curve_curve(
        &trimmed,
        &line,
        ogeom_intersect::ExtremaOptions::default(),
        tol,
    )?;
    if near
        .approaches
        .iter()
        .any(|a| a.distance <= tol.confusion() && interior(a.point_a))
    {
        ogeom_bail!(
            Construction,
            "the profile touches the axis away from its ends; revolving it \
             would sweep a surface through itself. Split the profile where it \
             meets the axis"
        );
    }
    Ok(AxisRelation::Clear)
}

/// Carry a parameter from one domain to the corresponding place in another.
///
/// A rigid motion leaves a curve's parameterization alone; a uniform scale
/// stretches it, because a line's parameter is a length. Rather than knowing
/// which curve types do which, the parameter is placed by where it sits between
/// the domain's ends — which is the same affine map in both cases, and the
/// identity when the two domains agree.
fn rescale(u: f64, from: (f64, f64), to: (f64, f64)) -> f64 {
    let span = from.1 - from.0;
    if span.abs() <= f64::MIN_POSITIVE {
        return to.0;
    }
    to.0 + (to.1 - to.0) * (u - from.0) / span
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
) -> OgeomResult<Shape> {
    let Some((start, end)) = crate::build::edge_vertices(model, edge)? else {
        ogeom_bail!(
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
        ogeom_bail!(Dangling, "vertex is not in this model");
    };
    let Some(data) = node.data().as_vertex() else {
        ogeom_bail!(Construction, "vertex node holds no point");
    };
    let from = base.transform(model.datums())?.apply(data.point);

    let line = ogeom_geom::LineCurve::segment(from, from + vector, tol)?;
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
    surface: ogeom_topo::SurfaceId,
    forward: ((f64, f64), (f64, f64)),
    reversed: ((f64, f64), (f64, f64)),
    tol: Tolerances,
) -> OgeomResult<()> {
    let flat = |p: (f64, f64)| Point2::new(p.0, p.1);
    let length = flat(forward.0).distance(flat(forward.1));
    let first = model
        .geometry_mut()
        .add_pcurve(Line2d::segment(flat(forward.0), flat(forward.1), tol)?.into());
    let second = model
        .geometry_mut()
        .add_pcurve(Line2d::segment(flat(reversed.0), flat(reversed.1), tol)?.into());

    let Some(node) = model.node_mut(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data_mut() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
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
    surface: ogeom_topo::SurfaceId,
    from: (f64, f64),
    to: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<()> {
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
    use ogeom_math::{Frame, Point};
    use ogeom_mesh::{Deflection, triangulate};
    use ogeom_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn deflection(chord: f64) -> Deflection {
        Deflection {
            chord,
            ..Deflection::default()
        }
    }

    /// One face of a box of `side`, named by its role.
    fn box_face(model: &mut Model, side: f64, role: ogeom_core::Role) -> Shape {
        let built = crate::make_box(model, Frame::WORLD, (side, side, side), T).unwrap();
        explore_unique(model, &built.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                model
                    .provenance_of(f)
                    .and_then(ogeom_core::Provenance::role)
                    == Some(role)
            })
            .expect("the box has a face with that role")
    }

    /// A square face in the xy plane, one unit on a side from the origin.
    fn square(model: &mut Model, side: f64) -> Shape {
        box_face(model, side, crate::primitive::roles::FACE_MAX_Z)
    }

    #[test]
    fn a_tapered_square_prism_is_the_frustum_the_closed_form_names() {
        let mut model = Model::new();
        let profile = square(&mut model, 10.0);
        let taper = 5.0_f64.to_radians();
        let built =
            crate::make_prism_tapered(&mut model, &profile, Vector::new(0.0, 0.0, 10.0), taper, T)
                .unwrap();
        let diagnosis = crate::check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

        // The far face measures 10 plus twice the height times the tangent
        // per side: its corner vertex says so directly.
        let d = 10.0 * taper.tan();
        let has_far_corner = explore_unique(&model, &built.shape, ShapeType::Vertex)
            .unwrap()
            .into_iter()
            .any(|v| {
                model
                    .node(&v)
                    .and_then(|n| n.data().as_vertex().map(|data| data.point))
                    .is_some_and(|p| {
                        (p.z - 20.0).abs() < 1e-9
                            && (p.x + d).abs() < 1e-9
                            && (p.y + d).abs() < 1e-9
                    })
            });
        assert!(has_far_corner, "the far ring widened by the taper");

        // All-planar, so the mesh integrates the frustum exactly.
        let (a0, a1) = (100.0, (10.0 + 2.0 * d) * (10.0 + 2.0 * d));
        let expected = 10.0 / 3.0 * (a1.mul_add(1.0, a0) + (a0 * a1).sqrt());
        let measured = volume_properties(&model, &built.shape, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!(
            (measured - expected).abs() < 1e-6,
            "tapered prism volume {measured} against {expected}"
        );
        assert!(!built.history.generated(&profile).is_empty());
    }

    #[test]
    fn a_hole_tapers_with_its_profile_into_a_cone() {
        let mut model = Model::new();
        // A 10 mm square with an off-centre round hole.
        let plane = ogeom_math::Plane::new(Frame::WORLD);
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ];
        let outer = crate::build::make_polygon(&mut model, &corners, true, T)
            .unwrap()
            .shape;
        let hole_centre = Point::new(3.5, 6.0, 0.0);
        let hole_r = 2.0;
        let circle = Circle::new(
            Frame::new(
                hole_centre,
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            hole_r,
            T,
        )
        .unwrap();
        let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
        let domain = curve.domain();
        let ring = crate::build::make_edge(&mut model, curve, domain, T)
            .unwrap()
            .shape;
        let hole = make_wire(&mut model, std::slice::from_ref(&ring), T)
            .unwrap()
            .shape;
        let surface: ogeom_geom::SurfaceGeometry =
            ogeom_geom::PlaneSurface::over(plane, (-20.0, 20.0), (-20.0, 20.0))
                .unwrap()
                .into();
        let outer_edges = model.ordered_children_of(&outer).unwrap();
        let profile = crate::build::make_face_with_pcurves(
            &mut model,
            surface,
            &[outer_edges, vec![ring.clone()]],
            T,
        )
        .unwrap()
        .shape;
        let _ = hole;

        let taper = 5.0_f64.to_radians();
        let built =
            crate::make_prism_tapered(&mut model, &profile, Vector::new(0.0, 0.0, 10.0), taper, T)
                .unwrap();
        let diagnosis = crate::check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

        // The hole's wall is a genuine cone, narrowing with height wherever
        // the hole sits in the profile.
        let cones = explore_unique(&model, &built.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                model
                    .node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface))
                    .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Cone(_)))
            })
            .count();
        assert_eq!(cones, 1, "the hole wall is a cone");

        let pi = core::f64::consts::PI;
        let d = 10.0 * taper.tan();
        let (a0, a1) = (100.0, (10.0 + 2.0 * d) * (10.0 + 2.0 * d));
        let outer_frustum = 10.0 / 3.0 * (a1.mul_add(1.0, a0) + (a0 * a1).sqrt());
        let r1 = hole_r - d;
        let hole_frustum = pi * 10.0 / 3.0 * (hole_r.mul_add(hole_r, hole_r * r1) + r1 * r1);
        let expected = outer_frustum - hole_frustum;
        let measured = volume_properties(&model, &built.shape, deflection(1e-3), T)
            .unwrap()
            .mass;
        assert!(
            (measured - expected).abs() < 5e-2,
            "holed tapered prism volume {measured} against {expected}"
        );
    }

    #[test]
    fn a_taper_that_collapses_a_hole_is_refused_by_name() {
        let mut model = Model::new();
        let plane = ogeom_math::Plane::new(Frame::WORLD);
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ];
        let outer = crate::build::make_polygon(&mut model, &corners, true, T)
            .unwrap()
            .shape;
        let circle = Circle::new(
            Frame::new(
                Point::new(5.0, 5.0, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            1.0,
            T,
        )
        .unwrap();
        let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
        let domain = curve.domain();
        let ring = crate::build::make_edge(&mut model, curve, domain, T)
            .unwrap()
            .shape;
        let hole = make_wire(&mut model, std::slice::from_ref(&ring), T)
            .unwrap()
            .shape;
        let surface: ogeom_geom::SurfaceGeometry =
            ogeom_geom::PlaneSurface::over(plane, (-20.0, 20.0), (-20.0, 20.0))
                .unwrap()
                .into();
        let outer_edges = model.ordered_children_of(&outer).unwrap();
        let profile = crate::build::make_face_with_pcurves(
            &mut model,
            surface,
            &[outer_edges, vec![ring.clone()]],
            T,
        )
        .unwrap()
        .shape;
        let _ = hole;
        let err = crate::make_prism_tapered(
            &mut model,
            &profile,
            Vector::new(0.0, 0.0, 10.0),
            8.0_f64.to_radians(),
            T,
        )
        .unwrap_err();
        assert!(err.to_string().contains("collapses"), "{err}");
    }

    #[test]
    fn a_curved_profile_edge_is_refused_by_name() {
        let mut model = Model::new();
        let plane = ogeom_math::Plane::new(Frame::WORLD);
        let ellipse = ogeom_math::Ellipse::new(Frame::WORLD, 4.0, 2.0, T).unwrap();
        let curve: ogeom_geom::Curve = ogeom_geom::EllipseCurve::new(ellipse).into();
        let domain = curve.domain();
        let ring = crate::build::make_edge(&mut model, curve, domain, T)
            .unwrap()
            .shape;
        let wire = make_wire(&mut model, std::slice::from_ref(&ring), T)
            .unwrap()
            .shape;
        let surface: ogeom_geom::SurfaceGeometry =
            ogeom_geom::PlaneSurface::over(plane, (-10.0, 10.0), (-10.0, 10.0))
                .unwrap()
                .into();
        let profile =
            crate::build::make_face_with_pcurves(&mut model, surface, &[vec![ring.clone()]], T)
                .unwrap()
                .shape;
        let _ = wire;
        let err = crate::make_prism_tapered(
            &mut model,
            &profile,
            Vector::new(0.0, 0.0, 5.0),
            5.0_f64.to_radians(),
            T,
        )
        .unwrap_err();
        assert!(err.to_string().contains("fitted ruling"), "{err}");
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
                ogeom_mesh::triangulate_face(&model, &face, deflection(0.01), T)
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
        ogeom_topo::explore(model, shape, ogeom_topo::Filter::OfType(ShapeType::Face)).unwrap()
    }

    /// A square profile in the xz plane, `offset` out from the z axis, `side`
    /// on a side, built as a face so it can be revolved.
    ///
    /// The corners run counter-clockwise about `-y`, so that is the plane's
    /// normal: a face whose wire winds against its own normal is inside out,
    /// and would sweep into a solid that measures negative — which is a
    /// property of the profile, not of the sweep.
    fn ring_profile(model: &mut Model, offset: f64, side: f64) -> Shape {
        let frame = Frame::new(
            Point::new(offset, 0.0, 0.0),
            -ogeom_math::Direction::Y,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let corners = [
            Point::new(offset, 0.0, 0.0),
            Point::new(offset + side, 0.0, 0.0),
            Point::new(offset + side, 0.0, side),
            Point::new(offset, 0.0, side),
        ];
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| model.add_vertex(ogeom_topo::VertexData::new(*p)))
            .collect();
        let edges: Vec<Shape> = (0..4)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                crate::build::make_edge_between(
                    model,
                    ogeom_geom::LineCurve::segment(a, b, T).unwrap().into(),
                    (0.0, a.distance(b)),
                    &vertices[i],
                    &vertices[(i + 1) % 4],
                    T,
                )
                .unwrap()
                .shape
            })
            .collect();
        let wire = crate::make_wire(model, &edges, T).unwrap().shape;
        let surface = model
            .geometry_mut()
            .add_surface(ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)).into());
        for (i, edge) in edges.iter().enumerate() {
            let (a, b) = (corners[i], corners[(i + 1) % 4]);
            let flat = |p: ogeom_math::Point| {
                let l = frame.to_local(p);
                Point2::new(l.x, l.y)
            };
            crate::attach_pcurve(
                model,
                edge,
                Line2d::segment(flat(a), flat(b), T).unwrap().into(),
                surface,
                ogeom_topo::Location::identity(),
                (0.0, a.distance(b)),
            )
            .unwrap();
        }
        crate::make_face_on(model, surface, std::slice::from_ref(&wire), T)
            .unwrap()
            .shape
    }

    #[test]
    fn a_square_revolved_a_full_turn_is_a_ring_that_agrees_with_itself() {
        // The case the reverted draft got wrong: correct topology, a closed
        // shell, per-face triangulations matching Pappus — and twelve unshared
        // triangle edges at the seam, because the two sides of a face closed in
        // `u` did not weld together.
        let (offset, side) = (3.0_f64, 2.0_f64);
        let mut model = Model::new();
        let profile = ring_profile(&mut model, offset, side);
        let built = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap();

        let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
        assert_eq!(counts(ShapeType::Face), 4, "one per profile edge, no caps");
        assert_eq!(
            counts(ShapeType::Edge),
            6,
            "a rail per profile vertex, and a seam only on the cylindrical \
             walls — the flat annuli are plane faces bounded by their rails \
             alone"
        );
        assert_eq!(counts(ShapeType::Vertex), 4, "a full turn adds none");

        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(is_shell_closed(&model, &shell).unwrap());
        assert!(
            crate::check(&model, &built.shape, T).unwrap().is_valid(),
            "{}",
            crate::check(&model, &built.shape, T).unwrap()
        );

        for face in explode(&model, &built.shape) {
            ogeom_mesh::triangulate_face(&model, &face, deflection(0.01), T)
                .unwrap_or_else(|e| panic!("a face would not triangulate: {e}"));
        }

        // Pappus: the volume is the profile's area times the distance its
        // centroid travels.
        let exact = side * side * TAU * (offset + side / 2.0);
        let found = crate::check_tessellation(&model, &built.shape, deflection(0.005), T).unwrap();
        assert!(found.is_valid(), "the mesh came apart: {found}");

        let mesh = triangulate(&model, &built.shape, deflection(0.005), T).unwrap();
        assert!(mesh.is_closed(), "the mesh has a slit in it");
        assert!(mesh.volume() > 0.0, "the solid is inside out");
        // Not bounded below by the exact value the way a convex solid's mesh
        // is: chords across the *inner* wall cut into the hole rather than into
        // the material, so they add volume where the outer wall's take it away.
        assert_relative_eq!(mesh.volume(), exact, max_relative = 1e-3);
    }

    #[test]
    fn a_square_revolved_part_way_has_two_ends_and_the_volume_of_that_wedge() {
        let (offset, side) = (3.0_f64, 2.0_f64);
        let angle = std::f64::consts::FRAC_PI_2;
        let mut model = Model::new();
        let profile = ring_profile(&mut model, offset, side);
        let built = make_revolution(&mut model, &profile, Axis::Z, angle, T).unwrap();

        let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
        assert_eq!(counts(ShapeType::Face), 6, "four sides and two ends");

        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(is_shell_closed(&model, &shell).unwrap());

        let exact = side * side * angle * (offset + side / 2.0);
        let mesh = triangulate(&model, &built.shape, deflection(0.005), T).unwrap();
        assert!(mesh.is_closed(), "the mesh has a slit in it");
        assert!(mesh.volume() > 0.0, "the solid is inside out");
        assert_relative_eq!(mesh.volume(), exact, max_relative = 1e-3);
        assert!(
            crate::check_tessellation(&model, &built.shape, deflection(0.005), T)
                .unwrap()
                .is_valid()
        );
    }

    /// A quadrilateral profile in the xz plane, wound counter-clockwise about
    /// `-y` so its wire agrees with its own normal.
    fn profile_from(model: &mut Model, corners: &[Point]) -> Shape {
        let frame = Frame::new(
            corners[0],
            -ogeom_math::Direction::Y,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let n = corners.len();
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| model.add_vertex(ogeom_topo::VertexData::new(*p)))
            .collect();
        let surface = model
            .geometry_mut()
            .add_surface(ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)).into());
        let flat = |p: ogeom_math::Point| {
            let l = frame.to_local(p);
            Point2::new(l.x, l.y)
        };

        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            let (a, b) = (corners[i], corners[(i + 1) % n]);
            let edge = crate::build::make_edge_between(
                model,
                ogeom_geom::LineCurve::segment(a, b, T).unwrap().into(),
                (0.0, a.distance(b)),
                &vertices[i],
                &vertices[(i + 1) % n],
                T,
            )
            .unwrap()
            .shape;
            crate::attach_pcurve(
                model,
                &edge,
                Line2d::segment(flat(a), flat(b), T).unwrap().into(),
                surface,
                ogeom_topo::Location::identity(),
                (0.0, a.distance(b)),
            )
            .unwrap();
            edges.push(edge);
        }
        let wire = crate::make_wire(model, &edges, T).unwrap().shape;
        crate::make_face_on(model, surface, std::slice::from_ref(&wire), T)
            .unwrap()
            .shape
    }

    #[test]
    fn a_rectangle_with_a_side_on_the_axis_revolves_into_a_cylinder_face_for_face() {
        // The claim the seam decision rests on: the same solid gets the same
        // counts whichever way it was built. Each lateral face is one face
        // closed on itself at a seam rather than two halves; a side lying
        // along the axis turns onto itself and contributes no face; and a
        // radial side sweeps a *plane* — the sweep names it as one, so the
        // caps are plane faces bounded by their rim circles alone, with no
        // seam and no degenerate centre, exactly as `make_cylinder` builds
        // them. Faces, edges and vertices all agree.
        let (radius, height) = (2.0_f64, 5.0_f64);
        let mut model = Model::new();
        let profile = profile_from(
            &mut model,
            &[
                Point::new(0.0, 0.0, 0.0),
                Point::new(radius, 0.0, 0.0),
                Point::new(radius, 0.0, height),
                Point::new(0.0, 0.0, height),
            ],
        );
        let revolved = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap();
        let primitive = crate::make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();

        let counts = |shape: &Shape, kind| explore_unique(&model, shape, kind).unwrap().len();
        assert_eq!(
            counts(&revolved.shape, ShapeType::Face),
            counts(&primitive.shape, ShapeType::Face),
            "a side and two caps, the same as make_cylinder"
        );
        assert_eq!(counts(&revolved.shape, ShapeType::Face), 3);
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            assert_eq!(
                counts(&revolved.shape, kind),
                counts(&primitive.shape, kind),
                "canonical caps carry a rim circle and nothing else: {kind:?}"
            );
        }

        let shell = explore_unique(&model, &revolved.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(is_shell_closed(&model, &shell).unwrap());
        assert!(
            crate::check(&model, &revolved.shape, T).unwrap().is_valid(),
            "{}",
            crate::check(&model, &revolved.shape, T).unwrap()
        );
        assert!(
            crate::check_tessellation(&model, &revolved.shape, deflection(0.005), T)
                .unwrap()
                .is_valid()
        );

        let exact = std::f64::consts::PI * radius * radius * height;
        let mesh = triangulate(&model, &revolved.shape, deflection(0.005), T).unwrap();
        assert!(mesh.is_closed());
        assert!(mesh.volume() > 0.0, "the solid is inside out");
        assert!(
            mesh.volume() < exact,
            "an inscribed volume cannot exceed it"
        );
        // Against the primitive at the same deflection, not against a bound
        // pulled out of the air: both inscribe the same cylinder with the same
        // chord, so they should agree to far better than either agrees with the
        // exact value.
        let reference = triangulate(&model, &primitive.shape, deflection(0.005), T).unwrap();
        assert_relative_eq!(mesh.volume(), reference.volume(), max_relative = 1e-6);
        assert!(
            mesh.volume() > exact * 0.995,
            "{} against {exact}",
            mesh.volume()
        );
    }

    #[test]
    fn a_triangle_touching_the_axis_revolves_into_a_cone() {
        // The endpoint on the axis sweeps out nothing, so its rail is a
        // degenerate edge — an apex. Leaving it out would leave the flank's
        // boundary open along one side of its parameter rectangle with nothing
        // for the triangulator to trim to.
        let (radius, height) = (3.0_f64, 4.0_f64);
        let mut model = Model::new();
        let profile = profile_from(
            &mut model,
            &[
                Point::new(0.0, 0.0, 0.0),
                Point::new(radius, 0.0, 0.0),
                Point::new(0.0, 0.0, height),
            ],
        );
        let built = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap();

        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Face)
                .unwrap()
                .len(),
            2,
            "a flank and one cap; the side on the axis sweeps out nothing"
        );
        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(is_shell_closed(&model, &shell).unwrap());
        assert!(
            crate::check_tessellation(&model, &built.shape, deflection(0.005), T)
                .unwrap()
                .is_valid()
        );

        let exact = std::f64::consts::PI * radius * radius * height / 3.0;
        let mesh = triangulate(&model, &built.shape, deflection(0.005), T).unwrap();
        assert!(mesh.volume() > 0.0, "the solid is inside out");
        assert!(mesh.volume() < exact);
        assert!(
            mesh.volume() > exact * 0.99,
            "{} against {exact}",
            mesh.volume()
        );
    }

    #[test]
    fn a_disc_revolved_a_full_turn_is_a_torus_seamed_both_ways() {
        // The case that decides whether seam handling is general: the profile
        // edge is *closed*, so the circle its one vertex sweeps bounds the face
        // across both the top and the bottom of the parameter rectangle. That
        // is a seam in `v` on a face already seamed in `u`, and the result has
        // to come out with the same counts `make_torus` gives for the same
        // solid.
        let (major, minor) = (5.0_f64, 2.0_f64);
        let mut model = Model::new();

        let frame = Frame::new(
            Point::new(major, 0.0, 0.0),
            -ogeom_math::Direction::Y,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let circle = ogeom_math::Circle::new(frame, minor, T).unwrap();
        let start = model.add_vertex(ogeom_topo::VertexData::new(Point::new(
            major + minor,
            0.0,
            0.0,
        )));
        let edge = crate::build::make_edge_between(
            &mut model,
            ogeom_geom::CircleCurve::new(circle).into(),
            (0.0, TAU),
            &start,
            &start,
            T,
        )
        .unwrap()
        .shape;
        let surface = model
            .geometry_mut()
            .add_surface(ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)).into());
        crate::attach_pcurve(
            &mut model,
            &edge,
            ogeom_geom::Circle2d::new(
                ogeom_math::Circle2::new(
                    ogeom_math::Frame2::new(Point2::ORIGIN, ogeom_math::Direction2::X),
                    minor,
                    T,
                )
                .unwrap(),
            )
            .into(),
            surface,
            ogeom_topo::Location::identity(),
            (0.0, TAU),
        )
        .unwrap();
        let wire = crate::make_wire(&mut model, std::slice::from_ref(&edge), T)
            .unwrap()
            .shape;
        let profile = crate::make_face_on(&mut model, surface, std::slice::from_ref(&wire), T)
            .unwrap()
            .shape;

        let built = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap();
        let primitive = crate::make_torus(&mut model, Frame::WORLD, major, minor, T).unwrap();

        let counts = |shape: &Shape, kind| explore_unique(&model, shape, kind).unwrap().len();
        for kind in [ShapeType::Face, ShapeType::Edge, ShapeType::Vertex] {
            assert_eq!(
                counts(&built.shape, kind),
                counts(&primitive.shape, kind),
                "{kind:?} count differs from make_torus's"
            );
        }
        assert_eq!(
            counts(&built.shape, ShapeType::Edge),
            2,
            "one seam each way"
        );

        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(is_shell_closed(&model, &shell).unwrap());
        assert!(
            crate::check_tessellation(&model, &built.shape, deflection(0.02), T)
                .unwrap()
                .is_valid()
        );

        let exact = 2.0 * std::f64::consts::PI * std::f64::consts::PI * major * minor * minor;
        let mesh = triangulate(&model, &built.shape, deflection(0.02), T).unwrap();
        assert!(mesh.is_closed(), "the mesh has a slit in it");
        assert!(mesh.volume() > 0.0, "the solid is inside out");
        assert!(
            mesh.volume() > exact * 0.99 && mesh.volume() < exact,
            "{} against {exact}",
            mesh.volume()
        );
    }

    #[test]
    fn a_full_turn_consumes_the_profile_face_but_not_its_edges() {
        // The profile of a full turn is an interior cross-section of the
        // result: no face of the solid is it, so it is deleted. Its edges are a
        // different matter — each survives as the seam of the face it made, and
        // reporting them deleted would break a reference to an edge that is
        // still right there.
        let mut model = Model::new();
        let profile = ring_profile(&mut model, 3.0, 2.0);
        let edge = model
            .children_of(&model.children_of(&profile).unwrap()[0])
            .unwrap()[0]
            .clone();

        let built = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap();
        assert!(
            built.history.is_deleted(&profile),
            "the profile is interior"
        );
        assert!(!built.history.is_deleted(&edge), "its edges are not");
        assert_eq!(
            built.history.generated(&edge).len(),
            1,
            "the lateral face it made"
        );

        // A partial turn keeps the profile as its near cap, so nothing is
        // deleted at all.
        let mut model = Model::new();
        let profile = ring_profile(&mut model, 3.0, 2.0);
        let partial = make_revolution(&mut model, &profile, Axis::Z, 1.0, T).unwrap();
        assert!(!partial.history.is_deleted(&profile));
    }

    #[test]
    fn a_profile_crossing_the_axis_is_refused_rather_than_swept_through_itself() {
        // The solid would have a finite, plausible volume that counts part of
        // space twice, and nothing downstream could tell.
        let mut model = Model::new();
        let profile = profile_from(
            &mut model,
            &[
                Point::new(-1.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 1.0),
                Point::new(-1.0, 0.0, 1.0),
            ],
        );
        let err = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap_err();
        assert!(
            err.to_string().contains("through itself"),
            "unexpected message: {err}"
        );

        // And the crossing is a third of the way along the bottom edge, which
        // no evenly spaced sample lands on. Testing the distance to the axis
        // would have missed it; testing which side of the axis the profile is
        // on does not.
        let mut model = Model::new();
        let grazing = profile_from(
            &mut model,
            &[
                Point::new(-1.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 3.0),
                Point::new(-1.0, 0.0, 3.0),
            ],
        );
        assert!(make_revolution(&mut model, &grazing, Axis::Z, TAU, T).is_err());
    }

    #[test]
    fn a_profile_grazing_the_axis_between_samples_is_refused_exactly() {
        // The case the sampled check could never see, and the reason the
        // exact one replaced it. The bottom of this profile is the quadratic
        // Bezier x(t) = (1 - 3t)^2: it dips to touch the axis tangentially at
        // t = 1/3 — not on any evenly spaced sample grid — and comes back
        // without ever changing side, so the radial direction never reverses
        // either. Sampling saw a profile clear of the axis; the extrema layer
        // sees the stationary approach that reaches it, and the revolution
        // would pinch to a point mid-face there.
        let mut model = Model::new();
        let frame = Frame::new(
            Point::new(1.0, 0.0, 0.0),
            -ogeom_math::Direction::Y,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let surface = model
            .geometry_mut()
            .add_surface(ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)).into());
        let flat = |p: ogeom_math::Point| {
            let l = frame.to_local(p);
            Point2::new(l.x, l.y)
        };

        let controls = [
            Point::new(1.0, 0.0, 0.0),
            Point::new(-2.0, 0.0, 0.5),
            Point::new(4.0, 0.0, 1.0),
        ];
        let corners = [
            controls[0],
            controls[2],
            Point::new(5.0, 0.0, 1.0),
            Point::new(5.0, 0.0, 0.0),
        ];
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| model.add_vertex(ogeom_topo::VertexData::new(*p)))
            .collect();

        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap();
        let dip: ogeom_geom::Curve =
            ogeom_geom::BSplineCurve::new(knots.clone(), controls.to_vec(), T)
                .unwrap()
                .into();
        let mut edges = vec![
            crate::build::make_edge_between(
                &mut model,
                dip,
                (0.0, 1.0),
                &vertices[0],
                &vertices[1],
                T,
            )
            .unwrap()
            .shape,
        ];
        crate::attach_pcurve(
            &mut model,
            &edges[0],
            ogeom_geom::BSpline2d::new(knots, controls.iter().map(|p| flat(*p)).collect(), T)
                .unwrap()
                .into(),
            surface,
            ogeom_topo::Location::identity(),
            (0.0, 1.0),
        )
        .unwrap();
        for i in 1..corners.len() {
            let (a, b) = (corners[i], corners[(i + 1) % corners.len()]);
            let edge = crate::build::make_edge_between(
                &mut model,
                ogeom_geom::LineCurve::segment(a, b, T).unwrap().into(),
                (0.0, a.distance(b)),
                &vertices[i],
                &vertices[(i + 1) % corners.len()],
                T,
            )
            .unwrap()
            .shape;
            crate::attach_pcurve(
                &mut model,
                &edge,
                Line2d::segment(flat(a), flat(b), T).unwrap().into(),
                surface,
                ogeom_topo::Location::identity(),
                (0.0, a.distance(b)),
            )
            .unwrap();
            edges.push(edge);
        }
        let wire = crate::build::make_wire(&mut model, &edges, T)
            .unwrap()
            .shape;
        let profile = crate::build::make_face_on(&mut model, surface, &[wire], T)
            .unwrap()
            .shape;

        let err = make_revolution(&mut model, &profile, Axis::Z, TAU, T).unwrap_err();
        assert!(
            err.to_string().contains("touches the axis"),
            "unexpected message: {err}"
        );
    }

    #[test]
    fn a_turn_that_goes_nowhere_or_too_far_is_refused() {
        let mut model = Model::new();
        let profile = ring_profile(&mut model, 3.0, 2.0);
        for angle in [0.0, -1.0, TAU * 1.5, f64::NAN, f64::INFINITY] {
            assert!(
                make_revolution(&mut model, &profile, Axis::Z, angle, T).is_err(),
                "accepted {angle}"
            );
        }
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
                    .and_then(ogeom_topo::EdgeData::curve3d)
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
    fn a_profile_that_has_been_placed_sweeps_where_it_actually_sits() {
        // A placed profile's edges arrive at a location, and the lateral
        // surface is built from the edge's *stored* curve. Building it without
        // the placement puts every side wall back at the origin.
        let mut model = Model::new();
        let face = square(&mut model, 2.0);
        let moved = crate::transformed(
            &mut model,
            &face,
            Transform::translation(Vector::new(10.0, 0.0, 0.0)),
        )
        .unwrap()
        .shape;

        let built = make_prism(&mut model, &moved, Vector::new(0.0, 0.0, 3.0), T).unwrap();
        let mesh = triangulate(&model, &built.shape, deflection(0.01), T).unwrap();
        assert!(mesh.is_closed(), "the mesh has a slit in it");
        assert_relative_eq!(mesh.volume(), 12.0, epsilon = 1e-9);

        let props = volume_properties(&model, &built.shape, deflection(0.01), T).unwrap();
        assert!(
            props.centre.distance(Point::new(11.0, 1.0, 3.5)) < 1e-9,
            "got {:?}",
            props.centre
        );
    }

    #[test]
    fn a_profile_placed_with_a_scale_sweeps_at_the_size_it_is_now() {
        // A placement may carry a uniform scale, and a scale stretches a line's
        // parameter with it, because that parameter is a length. The edge's
        // range is in the stored curve's parameter and the lateral surface's
        // `u` is in the placed one's, so a range copied across unchanged would
        // trim the surface at the wrong place — here, at half of it.
        let mut model = Model::new();
        let face = square(&mut model, 2.0);
        let scaled = crate::transformed(
            &mut model,
            &face,
            Transform::scaling(Point::ORIGIN, 2.0, T).unwrap(),
        )
        .unwrap()
        .shape;

        let built = make_prism(&mut model, &scaled, Vector::new(0.0, 0.0, 3.0), T).unwrap();
        let mesh = triangulate(&model, &built.shape, deflection(0.01), T).unwrap();
        assert!(mesh.is_closed(), "the mesh has a slit in it");
        // A four-by-four square, three tall.
        assert_relative_eq!(mesh.volume(), 48.0, epsilon = 1e-9);
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
