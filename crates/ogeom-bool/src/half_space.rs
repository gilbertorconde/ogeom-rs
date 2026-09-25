//! Half spaces bounded by curved faces, resolved into finite solids.
//!
//! A half space is one face and the side of it `make_half_space` was told
//! the material is on. The boolean acts on solids, so a half space is
//! resolved against the other operand into the finite solid that agrees
//! with it wherever the other operand reaches: the part of the material
//! side inside a box past the other operand's bound.
//!
//! A closed surface divides all of space. A drum, a cone, a ball or a
//! ring is resolved into its own primitive where the material is inside
//! it, and into a box with the primitive cut out where the material is
//! outside; the primitive's surface is the face's own, so the result is
//! exact on it. A drum and a cone are read as unbounded along their axes,
//! whatever height their surfaces declare.
//!
//! An open surface (a spline patch, an extrusion, a revolution, a trimmed
//! or offset surface of those) divides space only across its extent. It is
//! resolved into the prism the face sweeps travelling into the material,
//! which is the material side wherever each line along that direction
//! crosses the face once and the other operand lies within the face's
//! reach. Where either fails the half space is refused by name.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface as _, SurfaceGeometry, Transformable as _};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// The face bounding a half space, where `shape` is one: a solid of one
/// shell of one face, whose shell is open (the face does not close on
/// itself) or closes inside out (the material outside a ball or a ring).
/// A ball or ring the right way out is an ordinary solid.
pub(crate) fn half_space_face(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    if model.kind_of(shape)? != ShapeType::Solid {
        return Ok(None);
    }
    let shells = ogeom_topo::explore_unique(model, shape, ShapeType::Shell)?;
    let faces = ogeom_topo::explore_unique(model, shape, ShapeType::Face)?;
    if shells.len() != 1 || faces.len() != 1 {
        return Ok(None);
    }
    if !ogeom_algo::is_shell_closed(model, &shells[0])? {
        return Ok(Some(faces[0].clone()));
    }
    // Closed on itself: a half space only if it encloses the outside, which
    // its drawn volume says by its sign.
    let coarse = ogeom_mesh::Deflection::default();
    let mesh = ogeom_mesh::triangulate(model, shape, coarse, tol)?;
    Ok((mesh.volume() < 0.0).then(|| faces[0].clone()))
}

/// The surface a face lies on, in world space, with trims and offsets of
/// the closed kinds seen through to the closed surface they are.
fn world_surface(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        ogeom_bail!(Construction, "the half space's face holds no face data");
    };
    let Some(stored) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;
    let placed = stored.transformed(&placement, tol)?;
    Ok(seen_through(&placed, tol).unwrap_or(placed))
}

/// A trimmed surface's basis, and an offset of a drum, ball or ring as the
/// drum, ball or ring it is; `None` where there is nothing to see through.
fn seen_through(surface: &SurfaceGeometry, tol: Tolerances) -> Option<SurfaceGeometry> {
    match surface {
        SurfaceGeometry::Trimmed(t) => {
            Some(seen_through(t.basis(), tol).unwrap_or_else(|| t.basis().clone()))
        }
        SurfaceGeometry::Offset(o) => {
            let basis = seen_through(o.basis(), tol).unwrap_or_else(|| o.basis().clone());
            let d = o.distance();
            match basis {
                SurfaceGeometry::Cylinder(c) => {
                    let cyl = c.cylinder();
                    let grown =
                        ogeom_math::Cylinder::new(cyl.frame(), cyl.radius() + d, tol).ok()?;
                    Some(
                        ogeom_geom::CylinderSurface::new(grown, c.domain().1)
                            .ok()?
                            .into(),
                    )
                }
                SurfaceGeometry::Sphere(s) => {
                    let ball = s.sphere();
                    let grown =
                        ogeom_math::Sphere::new(ball.frame(), ball.radius() + d, tol).ok()?;
                    Some(ogeom_geom::SphereSurface::new(grown).into())
                }
                SurfaceGeometry::Torus(t) => {
                    let ring = t.torus();
                    let grown = ogeom_math::Torus::new(
                        ring.frame(),
                        ring.major_radius(),
                        ring.minor_radius() + d,
                        tol,
                    )
                    .ok()?;
                    Some(ogeom_geom::TorusSurface::new(grown).into())
                }
                _ => None,
            }
        }
        _ => None,
    }
}

/// The other operand's bound in `frame`'s coordinates: its eight corners
/// carried over, low and high.
fn local_extent(frame: &Frame, bound: &ogeom_math::Aabb) -> Option<(Point, Point)> {
    let (lo, hi) = (bound.low()?, bound.high()?);
    let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for k in 0..8 {
        let corner = Point::new(
            if k & 1 == 0 { lo.x } else { hi.x },
            if k & 2 == 0 { lo.y } else { hi.y },
            if k & 4 == 0 { lo.z } else { hi.z },
        );
        let local = frame.to_local(corner);
        low = Point::new(low.x.min(local.x), low.y.min(local.y), low.z.min(local.z));
        high = Point::new(
            high.x.max(local.x),
            high.y.max(local.y),
            high.z.max(local.z),
        );
    }
    Some((low, high))
}

/// A box in `frame` from `low` to `high` in its coordinates.
fn framed_box(
    model: &mut Model,
    frame: &Frame,
    low: Point,
    high: Point,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let corner = Frame::new(frame.to_world(low), frame.z(), frame.x(), tol)?;
    let size = (high.x - low.x, high.y - low.y, high.z - low.z);
    Ok(ogeom_algo::make_box(model, corner, size, tol)?.shape)
}

/// A vector in `frame`'s coordinates, in the world's.
fn world_vector(frame: &Frame, v: Vector) -> Vector {
    frame.x().vector() * v.x + frame.y().vector() * v.y + frame.z().vector() * v.z
}

/// `frame` moved along its own `z` to `height`.
fn raised(frame: &Frame, height: f64, tol: Tolerances) -> OgeomResult<Frame> {
    Frame::new(
        frame.to_world(Point::new(0.0, 0.0, height)),
        frame.z(),
        frame.x(),
        tol,
    )
}

/// A half space bounded by a curved face, resolved against `other` into
/// the finite solid that agrees with it wherever `other` reaches.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// the surface is open and `other` reaches past it, or it folds back over
/// itself along the direction into the material; or where `other` lies
/// wholly beyond a cone's apex, on the side the cone does not reach.
pub(crate) fn resolved(
    model: &mut Model,
    face: &Shape,
    other: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let surface = world_surface(model, face, tol)?;
    let (at, outward) = ogeom_algo::face_normal(model, face, tol)?;
    let bound = ogeom_algo::shape_bounds(model, other, tol)?;
    if bound.low().is_none() {
        ogeom_bail!(
            Construction,
            "the other argument has no bound to fill against"
        );
    }
    let margin = bound.diagonal().max(tol.confusion() * 1e3) * 0.5 + tol.confusion() * 1e3;
    // Whether the material is inside the closed surface: its outward
    // normal leads away from the axis or centre.
    let inside_of = |away: Vector| away.dot(outward) > 0.0;
    match surface {
        SurfaceGeometry::Cylinder(c) => {
            let cylinder = c.cylinder();
            let frame = cylinder.frame();
            let axis = cylinder.axis();
            let rel = at - axis.location;
            let radial = rel - axis.direction.vector() * rel.dot(axis.direction.vector());
            let Some((low, high)) = local_extent(&frame, &bound) else {
                ogeom_bail!(
                    Construction,
                    "the other argument has no bound to fill against"
                );
            };
            let r = cylinder.radius();
            if inside_of(radial) {
                let start = raised(&frame, low.z - margin, tol)?;
                let length = high.z - low.z + 2.0 * margin;
                Ok(ogeom_algo::make_cylinder(model, start, r, length, tol)?.shape)
            } else {
                let pad = Vector::new(margin, margin, margin);
                let low = Point::new(low.x.min(-r), low.y.min(-r), low.z) - pad;
                let high = Point::new(high.x.max(r), high.y.max(r), high.z) + pad;
                let block = framed_box(model, &frame, low, high, tol)?;
                let start = raised(&frame, low.z - margin, tol)?;
                let drum =
                    ogeom_algo::make_cylinder(model, start, r, high.z - low.z + 2.0 * margin, tol)?
                        .shape;
                Ok(crate::cut(model, &block, &drum, tol)?.shape)
            }
        }
        SurfaceGeometry::Sphere(s) => {
            let ball = s.sphere();
            let solid = ogeom_algo::make_sphere(model, ball.frame(), ball.radius(), tol)?.shape;
            if inside_of(at - ball.centre()) {
                return Ok(solid);
            }
            let frame = ball.frame();
            let Some((low, high)) = local_extent(&frame, &bound) else {
                ogeom_bail!(
                    Construction,
                    "the other argument has no bound to fill against"
                );
            };
            let r = ball.radius();
            let pad = Vector::new(margin, margin, margin);
            let low = Point::new(low.x.min(-r), low.y.min(-r), low.z.min(-r)) - pad;
            let high = Point::new(high.x.max(r), high.y.max(r), high.z.max(r)) + pad;
            let block = framed_box(model, &frame, low, high, tol)?;
            Ok(crate::cut(model, &block, &solid, tol)?.shape)
        }
        SurfaceGeometry::Torus(t) => {
            let ring = t.torus();
            let frame = ring.frame();
            let local = frame.to_local(at);
            let flat = Vector::new(local.x, local.y, 0.0);
            let tube = if flat.magnitude() > tol.confusion() {
                flat * (ring.major_radius() / flat.magnitude())
            } else {
                Vector::new(ring.major_radius(), 0.0, 0.0)
            };
            let away = world_vector(&frame, Vector::new(local.x, local.y, local.z) - tube);
            let solid = ogeom_algo::make_torus(
                model,
                frame,
                ring.major_radius(),
                ring.minor_radius(),
                tol,
            )?
            .shape;
            if inside_of(away) {
                return Ok(solid);
            }
            let Some((low, high)) = local_extent(&frame, &bound) else {
                ogeom_bail!(
                    Construction,
                    "the other argument has no bound to fill against"
                );
            };
            let reach = ring.major_radius() + ring.minor_radius();
            let pad = Vector::new(margin, margin, margin);
            let low = Point::new(low.x.min(-reach), low.y.min(-reach), low.z.min(-reach)) - pad;
            let high = Point::new(high.x.max(reach), high.y.max(reach), high.z.max(reach)) + pad;
            let block = framed_box(model, &frame, low, high, tol)?;
            Ok(crate::cut(model, &block, &solid, tol)?.shape)
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let frame = cone.frame();
            let slope = cone.half_angle().tan();
            let apex = -cone.reference_radius() / slope;
            let radius_at = |h: f64| (cone.reference_radius() + h * slope).max(0.0);
            let local = frame.to_local(at);
            let radial = world_vector(&frame, Vector::new(local.x, local.y, 0.0));
            let Some((low, high)) = local_extent(&frame, &bound) else {
                ogeom_bail!(
                    Construction,
                    "the other argument has no bound to fill against"
                );
            };
            let (h0, h1) = ((low.z - 2.0 * margin).max(apex), high.z + 2.0 * margin);
            if h1 <= apex + tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the other argument lies wholly beyond the cone's apex, where \
                     its surface does not reach"
                );
            }
            let start = raised(&frame, h0, tol)?;
            let frustum =
                ogeom_algo::make_cone(model, start, radius_at(h0), radius_at(h1), h1 - h0, tol)?
                    .shape;
            if inside_of(radial) {
                return Ok(frustum);
            }
            let reach = radius_at(h1);
            let pad = Vector::new(margin, margin, margin);
            let low = Point::new(low.x.min(-reach), low.y.min(-reach), low.z) - pad;
            let high = Point::new(high.x.max(reach), high.y.max(reach), high.z) + pad;
            let block = framed_box(model, &frame, low, high, tol)?;
            Ok(crate::cut(model, &block, &frustum, tol)?.shape)
        }
        _ => swept_into_material(model, face, outward, at, &bound, margin, tol),
    }
}

/// An open surface's half space over the other operand: the prism the face
/// sweeps into the material, far enough to pass the other operand.
fn swept_into_material(
    model: &mut Model,
    face: &Shape,
    outward: Vector,
    at: Point,
    bound: &ogeom_math::Aabb,
    margin: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let into = Direction::new(-outward, tol)?.vector();
    let face = &bounded(model, face, tol)?;
    let fine = ogeom_mesh::Deflection::default();
    let mesh = ogeom_mesh::triangulate_face(model, face, fine, tol)?;
    // Every line along `into` must cross the face once: no drawn triangle
    // may turn edge-on to it or face back along it.
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let n = (b - a).cross(c - a);
        let m = n.magnitude();
        if m <= tol.confusion() * tol.confusion() {
            continue;
        }
        if n.dot(outward) / m <= 1e-3 {
            ogeom_bail!(
                Construction,
                "the half space's surface folds back along the direction into its \
                 material, so a sweep of it is not its side"
            );
        }
    }
    // And the other operand must lie across the face's reach: every corner
    // of its bound on a line along `into` that meets the face.
    let (Some(lo), Some(hi)) = (bound.low(), bound.high()) else {
        ogeom_bail!(
            Construction,
            "the other argument has no bound to fill against"
        );
    };
    let meets = |p: Point| -> bool {
        mesh.triangles.iter().any(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            crosses(p, into, a, b, c)
        })
    };
    let mut deepest = 0.0_f64;
    for k in 0..8 {
        let corner = Point::new(
            if k & 1 == 0 { lo.x } else { hi.x },
            if k & 2 == 0 { lo.y } else { hi.y },
            if k & 4 == 0 { lo.z } else { hi.z },
        );
        if !meets(corner) {
            ogeom_bail!(
                Construction,
                "the other argument reaches past the half space's surface, which \
                 divides space only across its own extent"
            );
        }
        deepest = deepest.max((corner - at).dot(into));
    }
    let reach = mesh
        .positions
        .iter()
        .map(|p| (at - *p).dot(into))
        .fold(0.0_f64, f64::max);
    let length = deepest + reach + margin;
    Ok(ogeom_algo::make_prism(model, face, into * length, tol)?.shape)
}

/// A face with edges round it: the face itself where it has them, and a
/// face bounded by its surface's whole domain where it is the surface's
/// natural face, whose sweep would have no sides. The four edges run along
/// the domain's sides, each fitted through the surface and carrying the
/// straight side as its image.
fn bounded(model: &mut Model, face: &Shape, tol: Tolerances) -> OgeomResult<Shape> {
    if !model.children_of(face)?.is_empty() {
        return Ok(face.clone());
    }
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        ogeom_bail!(Construction, "the half space's face holds no face data");
    };
    let id = data.surface;
    let Some(surface) = model.geometry().surface(id).cloned() else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    if [u0, u1, v0, v1]
        .iter()
        .any(|x| !x.is_finite() || x.abs() > 1e6)
    {
        ogeom_bail!(
            Construction,
            "the half space's surface reaches without bound; a face trimmed to              the part of it that matters divides space across that part"
        );
    }
    let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
    let mut vertices = Vec::with_capacity(4);
    for &(u, v) in &corners {
        let p = surface.point_at(u, v, tol)?;
        vertices.push(ogeom_algo::make_vertex(model, p).shape);
    }
    const SAMPLES: u32 = 32;
    let mut edges = Vec::with_capacity(4);
    for k in 0..4 {
        let (from, to) = (corners[k], corners[(k + 1) % 4]);
        let params: Vec<f64> = (0..=SAMPLES)
            .map(|i| f64::from(i) / f64::from(SAMPLES))
            .collect();
        let mut points = Vec::with_capacity(params.len());
        for &t in &params {
            let (u, v) = (from.0 + (to.0 - from.0) * t, from.1 + (to.1 - from.1) * t);
            points.push(surface.point_at(u, v, tol)?);
        }
        let fitted = ogeom_geom::fit::fit_points_at(&params, &points, 3, tol.confusion(), tol)?;
        let edge = ogeom_algo::make_edge_between(
            model,
            fitted.curve.into(),
            (0.0, 1.0),
            &vertices[k],
            &vertices[(k + 1) % 4],
            tol,
        )?
        .shape;
        let side = ogeom_geom::Line2d::segment(
            ogeom_math::Point2::new(from.0, from.1),
            ogeom_math::Point2::new(to.0, to.1),
            tol,
        )?;
        let length = (to.0 - from.0).hypot(to.1 - from.1);
        ogeom_algo::attach_pcurve(
            model,
            &edge,
            side.into(),
            id,
            ogeom_topo::Location::identity(),
            (0.0, length),
        )?;
        edges.push(edge);
    }
    let wire = ogeom_algo::make_wire(model, &edges, tol)?.shape;
    let built = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape;
    let placed = built.moved(face.location());
    Ok(if face.orientation() == ogeom_topo::Orientation::Reversed {
        placed.reversed()
    } else {
        placed
    })
}

/// Whether the line through `p` along `d` crosses the triangle `a b c`.
fn crosses(p: Point, d: Vector, a: Point, b: Point, c: Point) -> bool {
    let (e1, e2) = (b - a, c - a);
    let h = d.cross(e2);
    let det = e1.dot(h);
    if det.abs() <= f64::EPSILON {
        return false;
    }
    let s = p - a;
    let u = s.dot(h) / det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let q = s.cross(e1);
    let v = d.dot(q) / det;
    v >= 0.0 && u + v <= 1.0
}
