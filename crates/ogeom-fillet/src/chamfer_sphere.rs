//! The chamfer on a rim where a plane meets a sphere.
//!
//! The seat is revolved, as a cylinder's or a cone's rim is: a planar cap
//! square to an axis, and a sphere whose centre stands on that axis. In
//! the rim's meridian half-plane, `(rho, h)` with `h` along the cap's
//! outward normal from the rim, the cap is the line `h = 0` and the wall
//! the sphere's great circle about `(0, c)`, `c` the centre's height, of
//! radius `S`. The rim stands at `(R, 0)`, at the latitude
//! `psi0 = atan2(-c, R)`, and a point of the wall at latitude `psi` at
//! `(S cos psi, c + S sin psi)`.
//!
//! The cap's leg runs `on_cap` along the cap, square to the rim; the
//! wall's runs `on_wall` down the great circle, the geodesic square to the
//! rim, to the latitude `psi0 - tau on_wall / S`. The bevel is the cone
//! (a cylinder where the radii agree) through the two contact rings, and
//! the wedge between it, the cap's annulus and the sphere's band is the
//! section turned about the axis: the whole turn for a closed rim, or the
//! arcs a boolean split one into, the arc's own sweep for an arc that
//! stops short of its circle.

use crate::support::{
    edge_curve, face_from_edges, revolved_annulus, same_occurrence, segment_between,
};
use ogeom_algo::{make_edge_between, make_revolution_band, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, ConeSurface, Curve, Curve3d as _, CylinderSurface, PlaneSurface, SphereSurface,
    SurfaceGeometry, Transformable as _,
};
use ogeom_math::{Circle, Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};
use ogeom_topo::{Filter, Model, NodeData, Orientation, Shape, ShapeType, explore};

/// A rim between a planar cap and a sphere centred on the rim's axis.
pub(crate) struct SphereSeat {
    /// The rim's centre.
    pub centre: Point,
    /// The rim's radius.
    pub radius: f64,
    /// The cap's outward unit normal.
    pub up: Vector,
    /// The rim frame's `x`, the angle every ring and chart is counted from.
    pub x_ref: Direction,
    /// The sphere centre's height along `up` from the rim's plane.
    pub height: f64,
    /// The sphere's radius.
    pub sphere_radius: f64,
    /// `+1` with the material inside the sphere, `-1` outside it.
    pub sigma: f64,
    /// `+1` when the wall runs away from the cap's outward side, the convex
    /// seats, and `-1` alongside it, the concave ones.
    pub tau: f64,
    /// The planar cap.
    pub cap_face: Shape,
    /// The sphere.
    pub wall_face: Shape,
    /// The sphere's own surface, placed: the chart the solid's face is on.
    pub wall_surface: SurfaceGeometry,
}

impl SphereSeat {
    /// The rim's latitude on the sphere, about `up`.
    fn rim_latitude(&self) -> f64 {
        (-self.height).atan2(self.radius)
    }

    /// The latitude a leg of `on_wall` down the wall from the rim reaches.
    pub fn wall_latitude(&self, on_wall: f64) -> f64 {
        (-self.tau).mul_add(on_wall / self.sphere_radius, self.rim_latitude())
    }

    /// The wall's point at `latitude`, in the meridian: `(rho, h)`.
    fn wall_at(&self, latitude: f64) -> (f64, f64) {
        let s = self.sphere_radius;
        (s * latitude.cos(), s.mul_add(latitude.sin(), self.height))
    }

    /// How far a leg may run down the wall before its ring closes on the
    /// sphere's pole.
    pub fn to_pole(&self) -> f64 {
        self.sphere_radius
            * self
                .tau
                .mul_add(self.rim_latitude(), core::f64::consts::FRAC_PI_2)
    }

    /// The cap contact's distance off the axis for a leg of `on_cap`.
    pub fn cap_rho(&self, on_cap: f64) -> f64 {
        (self.sigma * self.tau).mul_add(-on_cap, self.radius)
    }

    /// The other leg of a bevel leaving the leg `distance` along `face` at
    /// `angle` from that face's tangent in the rim's section, into the
    /// wedge's side: where the line meets the other face. Returned as
    /// `(on_wall, on_cap)`.
    pub fn angle_legs(
        &self,
        on_cap_face: bool,
        distance: f64,
        angle: f64,
        tol: Tolerances,
    ) -> OgeomResult<(f64, f64)> {
        let st = self.sigma * self.tau;
        let never = || {
            ogeom_core::ogeom_err!(
                Construction,
                "the bevel at that angle never meets the edge's other face"
            )
        };
        if on_cap_face {
            // From the cap's contact, back toward the rim and down toward
            // the wall's side; where that first meets the great circle.
            let c = (self.cap_rho(distance), 0.0);
            let dir = (st * angle.cos(), -self.tau * angle.sin());
            let rel = (c.0, c.1 - self.height);
            let b = rel.0.mul_add(dir.0, rel.1 * dir.1);
            let k = rel.0.mul_add(rel.0, rel.1 * rel.1) - self.sphere_radius.powi(2);
            let disc = b.mul_add(b, -k);
            if disc < 0.0 {
                return Err(never());
            }
            // Inside the sphere one root lies ahead; outside it the first
            // ahead is where the line reaches the wall.
            let roots = [-b - disc.sqrt(), -b + disc.sqrt()];
            let Some(t) = roots.into_iter().find(|t| *t > tol.confusion()) else {
                return Err(never());
            };
            let p = (t.mul_add(dir.0, c.0), t.mul_add(dir.1, c.1));
            if p.0 <= tol.confusion() {
                return Err(never());
            }
            let latitude = (p.1 - self.height).atan2(p.0);
            let on_wall = -self.tau * (latitude - self.rim_latitude()) * self.sphere_radius;
            if on_wall <= tol.confusion() {
                return Err(never());
            }
            Ok((on_wall, distance))
        } else {
            // From the wall's contact, along the great circle's tangent back
            // toward the rim, turned toward the wedge's side; where that
            // meets the cap's line.
            let latitude = self.wall_latitude(distance);
            let w = self.wall_at(latitude);
            let back = (-self.tau * latitude.sin(), self.tau * latitude.cos());
            let into = (-st * latitude.cos(), -st * latitude.sin());
            let dir = (
                back.0.mul_add(angle.cos(), into.0 * angle.sin()),
                back.1.mul_add(angle.cos(), into.1 * angle.sin()),
            );
            if dir.1.abs() <= tol.angular() {
                return Err(never());
            }
            let t = -w.1 / dir.1;
            if t <= tol.confusion() {
                return Err(never());
            }
            let rho = t.mul_add(dir.0, w.0);
            let on_cap = st * (self.radius - rho);
            if on_cap <= tol.confusion() || rho <= tol.confusion() {
                return Err(never());
            }
            Ok((distance, on_cap))
        }
    }
}

/// Read a rim's seat between a planar cap and a sphere: `None` where the
/// rim's faces are not one plane and one sphere.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// plane is not square to the rim's axis through its centre, the sphere is
/// not centred on the axis through the rim, or the side the wall runs from
/// the cap cannot be read at the rim.
pub(crate) fn sphere_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    tol: Tolerances,
) -> OgeomResult<Option<SphereSeat>> {
    let rim_circle = rim.circle();
    let rim_centre = rim_circle.centre();
    let rim_radius = rim_circle.radius();
    let rim_axis = rim_circle.frame().z().vector();
    let mut cap: Option<(Shape, Vector)> = None;
    let mut wall: Option<(Shape, Point, f64, f64, SurfaceGeometry)> = None;
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
                        "the rim's cap is not the perpendicular plane through it"
                    );
                }
                cap = Some((face.clone(), normal));
            }
            Some(SurfaceGeometry::Sphere(s)) => {
                let sphere = s.sphere();
                let centre = placement.apply(sphere.centre());
                let to_centre = centre - rim_centre;
                let along = to_centre.dot(rim_axis);
                let off_axis = (to_centre - rim_axis * along).magnitude();
                if off_axis > tol.confusion()
                    || (rim_radius.hypot(along) - sphere.radius()).abs() > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's sphere is not centred on the rim's axis \
                         through it"
                    );
                }
                let sigma = if reversed { -1.0 } else { 1.0 };
                let placed = SurfaceGeometry::Sphere(*s).transformed(&placement, tol)?;
                wall = Some((face.clone(), centre, sphere.radius(), sigma, placed));
            }
            _ => return Ok(None),
        }
    }
    let (
        Some((cap_face, up_raw)),
        Some((wall_face, sphere_centre, sphere_radius, sigma, wall_surface)),
    ) = (cap, wall)
    else {
        return Ok(None);
    };
    let up = up_raw / up_raw.magnitude();
    let mut seat = SphereSeat {
        centre: rim_centre,
        radius: rim_radius,
        up,
        x_ref: rim_circle.frame().x(),
        height: (sphere_centre - rim_centre).dot(up),
        sphere_radius,
        sigma,
        tau: 0.0,
        cap_face,
        wall_face,
        wall_surface,
    };
    // Which side of the cap's plane the wall runs, read at the edge itself:
    // just off the rim along the great circle each way, the wall stands on
    // one side only.
    let (curve, range) = edge_curve(model, edge, tol)?;
    let at = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    let radial = {
        let off = at - rim_centre;
        let off = off - up * off.dot(up);
        off / off.magnitude()
    };
    let eps = rim_radius * 1e-2;
    let deflection = ogeom_mesh::Deflection {
        chord: eps * 0.1,
        ..ogeom_mesh::Deflection::default()
    };
    let probe = |latitude: f64| {
        let (rho, h) = seat.wall_at(latitude);
        rim_centre + up * h + radial * rho
    };
    let stands = |p: Point| -> OgeomResult<bool> {
        Ok(
            ogeom_algo::classify_on_face(model, &seat.wall_face, p, deflection, tol)?
                == ogeom_algo::Containment::In,
        )
    };
    let step = eps / sphere_radius;
    let rim_latitude = seat.rim_latitude();
    seat.tau = match (
        stands(probe(rim_latitude + step))?,
        stands(probe(rim_latitude - step))?,
    ) {
        (true, false) => -1.0,
        (false, true) => 1.0,
        _ => ogeom_bail!(
            Construction,
            "the sphere stands on both sides of the rim's cap, or on neither"
        ),
    };
    Ok(Some(seat))
}

/// The bevel's wedge on a sphere's rim: legs `on_wall` down the sphere's
/// great circle and `on_cap` along the cap, the bevel the cone through
/// the two contact rings, over the whole turn or over the arc `span`
/// (its start about the seat's `x_ref` and its sweep).
pub(crate) fn sphere_bevel(
    model: &mut Model,
    seat: &SphereSeat,
    on_wall: f64,
    on_cap: f64,
    span: Option<(f64, f64)>,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let conf = tol.confusion();
    for distance in [on_wall, on_cap] {
        if !distance.is_finite() || distance <= conf {
            ogeom_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    let cap_rho = seat.cap_rho(on_cap);
    if cap_rho <= conf {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_cap} along the cap swallows the axis of a rim \
             of radius {}",
            seat.radius
        );
    }
    if on_wall >= seat.to_pole() {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_wall} down the sphere runs past its pole"
        );
    }
    let latitude = seat.wall_latitude(on_wall);
    let (wall_rho, wall_h) = seat.wall_at(latitude);
    if wall_rho <= conf {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_wall} down the sphere runs past its pole"
        );
    }
    // The bevel's chord must stay on the wedge's side of the sphere: from
    // the wall's contact `W` toward the cap's `C`, the line meets the great
    // circle again at `t* = -2 (W - O).(C - W) / |C - W|²` along it, which
    // must not fall between the two.
    let chord = (cap_rho - wall_rho, -wall_h);
    let spoke = (wall_rho, wall_h - seat.height);
    let length2 = chord.0.mul_add(chord.0, chord.1 * chord.1);
    let again = -2.0 * spoke.0.mul_add(chord.0, spoke.1 * chord.1) / length2;
    if again * length2.sqrt() > conf && again < 1.0 {
        ogeom_bail!(
            Construction,
            "the bevel between the chamfer's contacts runs through the sphere"
        );
    }
    // The wedge's outward side on the bevel: square to the chord, away
    // from the rim. A cone's natural normal points off the axis.
    let bevel_turned = {
        let normal = (chord.1, -chord.0);
        let to_rim = (seat.radius - cap_rho, 0.0);
        let normal = if normal.0.mul_add(to_rim.0, normal.1 * to_rim.1) > 0.0 {
            (-normal.0, -normal.1)
        } else {
            normal
        };
        normal.0 < 0.0
    };
    let geometry = Section {
        cap_rho,
        wall_rho,
        wall_h,
        rim_latitude: seat.rim_latitude(),
        wall_latitude: latitude,
        bevel_turned,
    };
    match span {
        None => whole_turn(model, seat, &geometry, tol),
        Some((start, sweep)) => sector(model, seat, &geometry, start, sweep, tol),
    }
}

/// The section of the wedge in the meridian.
struct Section {
    cap_rho: f64,
    wall_rho: f64,
    wall_h: f64,
    rim_latitude: f64,
    wall_latitude: f64,
    /// Whether the bevel's face turns against its cone's natural normal.
    bevel_turned: bool,
}

/// The bevel's surface: the cone through the cap's contact ring and the
/// wall's, a cylinder where their radii agree, framed by `frame` at the
/// cap's level, its chart kept clear of its apex.
fn bevel_surface(frame: Frame, section: &Section, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let slope = (section.wall_rho - section.cap_rho) / section.wall_h;
    let rows = (section.wall_h.min(0.0), section.wall_h.max(0.0));
    let depth = section.wall_h.abs();
    if slope.atan().abs() <= tol.angular() {
        let pad = 0.1 * depth;
        return Ok(CylinderSurface::new(
            Cylinder::new(frame, section.cap_rho, tol)?,
            (rows.0 - pad, rows.1 + pad),
        )?
        .into());
    }
    let apex = -section.cap_rho / slope;
    let clear = (apex - rows.0).abs().min((apex - rows.1).abs());
    let pad = (0.1 * depth).min(0.5 * clear);
    let cone = Cone::new(frame, section.cap_rho, slope.atan(), tol)?;
    Ok(ConeSurface::new(cone, (rows.0 - pad, rows.1 + pad))?.into())
}

/// The sphere framed about the seat's axis, its angles counted from `x`.
fn wall_surface(seat: &SphereSeat, x: Direction, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let frame = Frame::new(
        seat.centre + seat.up * seat.height,
        Direction::new(seat.up, tol)?,
        x,
        tol,
    )?;
    Ok(SphereSurface::new(Sphere::new(frame, seat.sphere_radius, tol)?).into())
}

/// The wedge over the whole turn: the sphere's band, the cap's annulus and
/// the bevel's band, between three rings.
fn whole_turn(
    model: &mut Model,
    seat: &SphereSeat,
    section: &Section,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let ring = |model: &mut Model, origin: Point, r: f64| -> OgeomResult<Shape> {
        let frame = Frame::new(origin, Direction::new(seat.up, tol)?, seat.x_ref, tol)?;
        let curve = Curve::Circle(CircleCurve::new(Circle::new(frame, r, tol)?));
        let domain = curve.domain();
        Ok(ogeom_algo::make_edge(model, curve, domain, tol)?.shape)
    };
    let apex_ring = ring(model, seat.centre, seat.radius)?;
    let wall_ring = ring(
        model,
        seat.centre + seat.up * section.wall_h,
        section.wall_rho,
    )?;
    let cap_ring = ring(model, seat.centre, section.cap_rho)?;

    let wall_band = {
        let surface = wall_surface(seat, seat.x_ref, tol)?;
        let band = make_revolution_band(model, &surface, &wall_ring, &apex_ring, tol)?;
        if seat.sigma * seat.tau < 0.0 {
            band.reversed()
        } else {
            band
        }
    };
    let annulus = revolved_annulus(
        model,
        seat.centre,
        seat.up,
        seat.tau,
        (seat.radius, section.cap_rho),
        (seat.radius + section.wall_h.abs() + section.cap_rho) * 2.0,
        (&apex_ring, &cap_ring),
        tol,
    )?;
    let bevel_band = {
        let frame = Frame::new(seat.centre, Direction::new(seat.up, tol)?, seat.x_ref, tol)?;
        let surface = bevel_surface(frame, section, tol)?;
        let band = make_revolution_band(model, &surface, &wall_ring, &cap_ring, tol)?;
        if section.bevel_turned {
            band.reversed()
        } else {
            band
        }
    };
    Ok(vec![wall_band, annulus, bevel_band])
}

/// The wedge over an arc of the rim, from `start` about the seat's `x_ref`
/// through `sweep`: the three revolved patches and the two section faces
/// closing its ends in the meridian half-planes.
fn sector(
    model: &mut Model,
    seat: &SphereSeat,
    section: &Section,
    start: f64,
    sweep: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let up = Direction::new(seat.up, tol)?;
    let seat_frame = Frame::new(seat.centre, up, seat.x_ref, tol)?;
    let dir_at = |theta: f64| -> Vector {
        seat_frame.x().vector() * theta.cos() + seat_frame.y().vector() * theta.sin()
    };
    let (theta0, theta1) = (start, start + sweep);
    let x0 = Direction::new(dir_at(theta0), tol)?;
    let wedge_frame_at = |origin: Point| Frame::new(origin, up, x0, tol);
    let wall_level = seat.centre + seat.up * section.wall_h;
    let apex_at = |theta: f64| seat.centre + dir_at(theta) * seat.radius;
    let wall_at = |theta: f64| wall_level + dir_at(theta) * section.wall_rho;
    let cap_at = |theta: f64| seat.centre + dir_at(theta) * section.cap_rho;

    let va0 = make_vertex(model, apex_at(theta0)).shape;
    let va1 = make_vertex(model, apex_at(theta1)).shape;
    let vw0 = make_vertex(model, wall_at(theta0)).shape;
    let vw1 = make_vertex(model, wall_at(theta1)).shape;
    let vc0 = make_vertex(model, cap_at(theta0)).shape;
    let vc1 = make_vertex(model, cap_at(theta1)).shape;

    let ring_arc = |model: &mut Model,
                    origin: Point,
                    r: f64,
                    from: &Shape,
                    to: &Shape|
     -> OgeomResult<Shape> {
        let circle = Circle::new(wedge_frame_at(origin)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };
    // The great circle's arc in the meridian at `theta`, from the rim to
    // the wall's contact.
    let sphere_centre = seat.centre + seat.up * seat.height;
    let meridian =
        |model: &mut Model, theta: f64, rim: &Shape, contact: &Shape| -> OgeomResult<Shape> {
            let dir = dir_at(theta);
            let (lo, hi, from, to) = if section.wall_latitude < section.rim_latitude {
                (section.wall_latitude, section.rim_latitude, contact, rim)
            } else {
                (section.rim_latitude, section.wall_latitude, rim, contact)
            };
            let x = Direction::new(dir * lo.cos() + seat.up * lo.sin(), tol)?;
            let z = Direction::new(dir.cross(seat.up), tol)?;
            let circle = Circle::new(
                Frame::new(sphere_centre, z, x, tol)?,
                seat.sphere_radius,
                tol,
            )?;
            let curve = Curve::Circle(CircleCurve::new(circle));
            let edge = make_edge_between(model, curve, (0.0, hi - lo), from, to, tol)?.shape;
            Ok(if from.is_same(rim) {
                edge
            } else {
                edge.reversed()
            })
        };

    let arc_apex = ring_arc(model, seat.centre, seat.radius, &va0, &va1)?;
    let arc_wall = ring_arc(model, wall_level, section.wall_rho, &vw0, &vw1)?;
    let arc_cap = ring_arc(model, seat.centre, section.cap_rho, &vc0, &vc1)?;
    let meridian0 = meridian(model, theta0, &va0, &vw0)?;
    let meridian1 = meridian(model, theta1, &va1, &vw1)?.reversed();
    let radial0 = segment_between(model, (&vc0, cap_at(theta0)), (&va0, apex_at(theta0)), tol)?;
    let radial1 = segment_between(model, (&va1, apex_at(theta1)), (&vc1, cap_at(theta1)), tol)?;
    let bevel0 = segment_between(model, (&vw0, wall_at(theta0)), (&vc0, cap_at(theta0)), tol)?;
    let bevel1 = segment_between(model, (&vw1, wall_at(theta1)), (&vc1, cap_at(theta1)), tol)?;
    let walked = |edges: Vec<Shape>, back: bool| -> Vec<Shape> {
        if back {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        }
    };

    // The sphere's patch, on the solid's own chart of the sphere, so the
    // melt reads the two as one. In a chart whose angle runs about `up`
    // and whose latitude runs up it, the ring of contact ring forward and
    // rim back is anticlockwise where the contact stands below the rim;
    // that is its turn about the outward normal, which the solid's chart
    // keeps where its own normal points out and turns where it points in.
    let wall_patch = {
        let surface = seat.wall_surface.clone();
        let inward = {
            use ogeom_geom::Surface as _;
            let ((u, _), (v0, v1)) = surface.domain();
            let v = f64::midpoint(v0, v1);
            let normal = surface.normal_at(u, v, tol)?.vector();
            normal.dot(surface.point_at(u, v, tol)? - sphere_centre) < 0.0
        };
        let edges = walked(
            vec![
                arc_wall.clone(),
                meridian1.clone(),
                arc_apex.reversed(),
                meridian0.clone(),
            ],
            (seat.tau < 0.0) != inward,
        );
        let face = face_from_edges(model, surface, &edges, tol)?;
        if (seat.sigma * seat.tau < 0.0) != inward {
            face.reversed()
        } else {
            face
        }
    };

    // The cap's sector between the rim and the contact ring, as the rim
    // fillet's.
    let reach = (seat.radius + section.wall_h.abs() + section.cap_rho + 1.0) * 2.0;
    let cap_sector = {
        let plane = Plane::through(seat.centre, Direction::new(seat.up * seat.tau, tol)?);
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let rim_wider = seat.radius > section.cap_rho;
        let edges = walked(
            vec![
                arc_apex.clone(),
                radial1.clone(),
                arc_cap.reversed(),
                radial0.clone(),
            ],
            rim_wider != (seat.tau > 0.0),
        );
        face_from_edges(model, surface, &edges, tol)?
    };

    // The bevel's patch: its chart's height runs up the axis, so the ring
    // of wall contact forward and cap contact back is anticlockwise where
    // the wall's contact stands below the cap.
    let bevel_patch = {
        let surface = bevel_surface(wedge_frame_at(seat.centre)?, section, tol)?;
        let edges = walked(
            vec![
                arc_wall.clone(),
                bevel1.clone(),
                arc_cap.reversed(),
                bevel0.reversed(),
            ],
            section.wall_h > 0.0,
        );
        let face = face_from_edges(model, surface, &edges, tol)?;
        if section.bevel_turned {
            face.reversed()
        } else {
            face
        }
    };

    // The ends: the section in the meridian half-planes, outward along
    // the rim's travel, each ring walked back where its corners turn
    // against `outward`.
    let tangent_at = |theta: f64| seat.up.cross(dir_at(theta));
    let end_cap = |model: &mut Model,
                   theta: f64,
                   outward: Vector,
                   edges: Vec<Shape>,
                   [p, q, r]: [Point; 3]|
     -> OgeomResult<Shape> {
        let plane = Plane::through(apex_at(theta), Direction::new(outward, tol)?);
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let edges = walked(edges, (q - p).cross(r - p).dot(outward) < 0.0);
        face_from_edges(model, surface, &edges, tol)
    };
    let cap_start = end_cap(
        model,
        theta0,
        -tangent_at(theta0),
        vec![meridian0.clone(), bevel0.clone(), radial0.clone()],
        [apex_at(theta0), wall_at(theta0), cap_at(theta0)],
    )?;
    let cap_end = end_cap(
        model,
        theta1,
        tangent_at(theta1),
        vec![meridian1.clone(), radial1.clone(), bevel1.reversed()],
        [wall_at(theta1), apex_at(theta1), cap_at(theta1)],
    )?;
    Ok(vec![
        wall_patch,
        cap_sector,
        bevel_patch,
        cap_start,
        cap_end,
    ])
}

/// Whether the rims on the seat's circle run its whole turn between them:
/// the edges of `solid` on `circle` that one face on the seat's sphere and
/// one on its cap's plane share, as a boolean leaves a rim it split into
/// arcs.
pub(crate) fn rim_closes(
    model: &Model,
    solid: &Shape,
    seat: &SphereSeat,
    circle: &Circle,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let reach = tol.confusion() * 100.0;
    let sphere_centre = seat.centre + seat.up * seat.height;
    let mut on_sphere: Vec<Shape> = Vec::new();
    let mut on_cap: Vec<Shape> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(surface) = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
        else {
            continue;
        };
        let placement = face.transform(model.datums())?;
        let kept = match surface {
            SurfaceGeometry::Sphere(s) => {
                let sphere = s.sphere();
                placement.apply(sphere.centre()).distance(sphere_centre) <= reach
                    && (sphere.radius() - seat.sphere_radius).abs() <= reach
            }
            SurfaceGeometry::Plane(p) => {
                let normal = placement.apply_vector(p.plane().normal().vector());
                let origin = placement.apply(p.plane().frame().origin());
                normal.cross(seat.up).magnitude() <= tol.angular()
                    && (origin - seat.centre).dot(seat.up).abs() <= reach
            }
            _ => false,
        };
        if !kept {
            continue;
        }
        let edges = ogeom_topo::explore_unique(model, &face, ShapeType::Edge)?;
        if matches!(surface, SurfaceGeometry::Sphere(_)) {
            on_sphere.extend(edges);
        } else {
            on_cap.extend(edges);
        }
    }
    let shared: Vec<Shape> = on_sphere
        .into_iter()
        .filter(|e| on_cap.iter().any(|c| c.is_same(e)))
        .collect();
    let frame = Frame::new(seat.centre, Direction::new(seat.up, tol)?, seat.x_ref, tol)?;
    crate::support::edges_cover(model, shared, circle, &frame, tol)
}
