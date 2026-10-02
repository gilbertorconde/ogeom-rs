//! Rounds on sheets where a face beside the corner is curved.
//!
//! Two families of seat have closed forms. Where both supports are swept
//! along one direction (planes parallel to it, cylinders whose axes run
//! along it) every section square to that direction is the same, and the
//! round is a cylinder along it touching each support along a line. Where
//! both supports turn about one axis (planes square to it, cylinders and
//! cones about it, spheres and tori centred on it) every half-plane through
//! the axis holds the same section, and the round is a torus about it
//! touching each support along a circle. In both the ball's section is a
//! circle of the radius tangent to two lines or circles in a plane, solved
//! exactly.
//!
//! Every other seat is marched: the ball's centre and contact points are
//! solved station by station along the rounded edge, and the round is the
//! surface fitted through the ball's arcs, its borders the lines of
//! contact.
//!
//! Each face the ball touches is rebuilt on its own surface with its
//! boundary re-routed along its line of contact, and the round shares that
//! edge with it.

use core::f64::consts::{PI, TAU};

use ogeom_algo::{
    Built, History, attach_pcurve, attach_seam, edge_vertices, make_edge, make_edge_between,
    make_face_on, make_vertex, make_wire, project_on_curve, project_on_surface, shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, CylinderSurface, LineCurve, PlanarCurve, Surface as _,
    SurfaceGeometry, TorusSurface, TrimmedCurve,
};
use ogeom_intersect::Marching;
use ogeom_math::{
    Axis, Axis2, Circle, Cylinder, Direction, Direction2, Frame, Point, Point2, Torus, Transform2,
    Vector, Vector2,
};
use ogeom_topo::{
    EdgeRepr, FaceData, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId,
    explore_unique,
};

use crate::march::{MarchedBlend, Sides, march_blend_sided, seat_section};
use crate::marched::{band_fit_target, fit_open_band, touching_ball};
use crate::support::{edge_curve, face_from_edges, segment_between};

/// A pair of coordinates in a section plane.
type P2 = (f64, f64);

fn sub2(a: P2, b: P2) -> P2 {
    (a.0 - b.0, a.1 - b.1)
}

fn add2(a: P2, b: P2) -> P2 {
    (a.0 + b.0, a.1 + b.1)
}

fn scale2(a: P2, k: f64) -> P2 {
    (a.0 * k, a.1 * k)
}

fn dot2(a: P2, b: P2) -> f64 {
    a.0.mul_add(b.0, a.1 * b.1)
}

fn cross2(a: P2, b: P2) -> f64 {
    a.0.mul_add(b.1, -(a.1 * b.0))
}

fn norm2(a: P2) -> f64 {
    a.0.hypot(a.1)
}

/// A quarter turn counter-clockwise.
const fn perp2(a: P2) -> P2 {
    (-a.1, a.0)
}

/// How the round's sections are laid out in space.
#[derive(Debug, Clone, Copy)]
enum Layout {
    /// Each section is the plane square to `along` at a distance `s` along
    /// it from `origin`, with coordinates along `x` and `y`.
    Extruded {
        origin: Point,
        along: Vector,
        x: Vector,
        y: Vector,
    },
    /// Each section is the half-plane through the axis at the angle `s`
    /// from `x` toward `y`, with coordinates the distance from the axis and
    /// the height along it. `x`, `y` and `axis` are right-handed.
    Revolved {
        origin: Point,
        axis: Vector,
        x: Vector,
        y: Vector,
    },
}

impl Layout {
    /// Sections square to `along`, with `origin` at station zero.
    fn extruded(origin: Point, along: Vector, tol: Tolerances) -> OgeomResult<Self> {
        let along = along.normalized(tol)?;
        let frame = Frame::new(
            origin,
            Direction::new(along, tol)?,
            any_square_to(along, tol)?,
            tol,
        )?;
        Ok(Self::Extruded {
            origin,
            along,
            x: frame.x().vector(),
            y: frame.y().vector(),
        })
    }

    /// Sections through the axis at `origin` along `axis`, angles measured
    /// from `x` (which need not be square to the axis).
    fn revolved(origin: Point, axis: Vector, x: Vector, tol: Tolerances) -> OgeomResult<Self> {
        let frame = Frame::new(
            origin,
            Direction::new(axis, tol)?,
            Direction::new(x, tol)?,
            tol,
        )?;
        Ok(Self::Revolved {
            origin,
            axis: frame.z().vector(),
            x: frame.x().vector(),
            y: frame.y().vector(),
        })
    }

    /// Where along the round a point stands: its distance along the
    /// direction, or its angle about the axis.
    fn station(&self, p: Point) -> f64 {
        match *self {
            Self::Extruded { origin, along, .. } => (p - origin).dot(along),
            Self::Revolved { origin, x, y, .. } => {
                let d = p - origin;
                d.dot(y).atan2(d.dot(x))
            }
        }
    }

    /// The radial direction of the half-plane at angle `s`.
    fn radial(x: Vector, y: Vector, s: f64) -> Vector {
        let (sin, cos) = s.sin_cos();
        x * cos + y * sin
    }

    /// A point's coordinates in its own section.
    fn flat(&self, p: Point) -> P2 {
        match *self {
            Self::Extruded { origin, x, y, .. } => {
                let d = p - origin;
                (d.dot(x), d.dot(y))
            }
            Self::Revolved { origin, axis, .. } => {
                let d = p - origin;
                let h = d.dot(axis);
                ((d - axis * h).magnitude(), h)
            }
        }
    }

    /// A vector at station `s` in section coordinates (its part in the
    /// section).
    fn flat_vector(&self, v: Vector, s: f64) -> P2 {
        match *self {
            Self::Extruded { x, y, .. } => (v.dot(x), v.dot(y)),
            Self::Revolved { axis, x, y, .. } => (v.dot(Self::radial(x, y, s)), v.dot(axis)),
        }
    }

    /// The point at section coordinates `q` in the section at station `s`.
    fn lift(&self, q: P2, s: f64) -> Point {
        match *self {
            Self::Extruded {
                origin,
                along,
                x,
                y,
            } => origin + x * q.0 + y * q.1 + along * s,
            Self::Revolved { origin, axis, x, y } => {
                origin + Self::radial(x, y, s) * q.0 + axis * q.1
            }
        }
    }

    /// A section vector lifted into space at station `s`.
    fn lift_vector(&self, v: P2, s: f64) -> Vector {
        match *self {
            Self::Extruded { x, y, .. } => x * v.0 + y * v.1,
            Self::Revolved { axis, x, y, .. } => Self::radial(x, y, s) * v.0 + axis * v.1,
        }
    }

    /// Whether stations go round once and come back.
    const fn turns(&self) -> bool {
        matches!(self, Self::Revolved { .. })
    }
}

/// A unit direction square to `v`.
fn any_square_to(v: Vector, tol: Tolerances) -> OgeomResult<Direction> {
    let trial = if v.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    Direction::new(trial - v * trial.dot(v), tol)
}

/// A support's trace in every section: a line or a circle.
#[derive(Debug, Clone, Copy)]
enum Profile {
    Line { at: P2, normal: P2 },
    Circle { centre: P2, radius: f64 },
}

impl Profile {
    /// The profile's own normal at `q` on it: a line's stated normal, a
    /// circle's outward one.
    fn normal_at(&self, q: P2) -> P2 {
        match *self {
            Self::Line { normal, .. } => normal,
            Self::Circle { centre, .. } => {
                let d = sub2(q, centre);
                scale2(d, 1.0 / norm2(d))
            }
        }
    }

    /// The point of the profile nearest the ball's centre `c`, where the
    /// ball of radius `r` on the `side` of it touches.
    fn foot(&self, c: P2, side: f64, r: f64) -> P2 {
        match *self {
            Self::Line { normal, .. } => sub2(c, scale2(normal, side * r)),
            Self::Circle { centre, radius } => {
                let d = sub2(c, centre);
                add2(centre, scale2(d, radius / norm2(d)))
            }
        }
    }

    /// Where a ball of radius `r` on the `side` of the profile has its
    /// centre: the profile offset by `r`, or nothing where a circle's
    /// offset closes up.
    fn offset(&self, side: f64, r: f64, tol: Tolerances) -> Option<Self> {
        match *self {
            Self::Line { at, normal } => Some(Self::Line {
                at: add2(at, scale2(normal, side * r)),
                normal,
            }),
            Self::Circle { centre, radius } => {
                let reach = side.mul_add(r, radius);
                (reach > tol.confusion()).then_some(Self::Circle {
                    centre,
                    radius: reach,
                })
            }
        }
    }

    /// Length along the profile from `from` to `to`, both on it, the short
    /// way round.
    fn length_between(&self, from: P2, to: P2) -> f64 {
        match *self {
            Self::Line { .. } => norm2(sub2(to, from)),
            Self::Circle { centre, radius } => {
                let (a, b) = (sub2(from, centre), sub2(to, centre));
                radius * cross2(a, b).abs().atan2(dot2(a, b))
            }
        }
    }
}

/// Where two profiles cross.
fn crossings(a: Profile, b: Profile, tol: Tolerances) -> Vec<P2> {
    match (a, b) {
        (Profile::Line { at: p, normal: n }, Profile::Line { at: q, normal: m }) => {
            let (d, e) = (perp2(n), perp2(m));
            let det = cross2(d, e);
            if det.abs() <= tol.angular() {
                return Vec::new();
            }
            let t = cross2(sub2(q, p), e) / det;
            vec![add2(p, scale2(d, t))]
        }
        (Profile::Line { at, normal }, Profile::Circle { centre, radius })
        | (Profile::Circle { centre, radius }, Profile::Line { at, normal }) => {
            let d = perp2(normal);
            let foot = add2(at, scale2(d, dot2(sub2(centre, at), d)));
            let h = norm2(sub2(centre, foot));
            if h > radius + tol.confusion() {
                return Vec::new();
            }
            let w = (radius * radius - h * h).max(0.0).sqrt();
            if w <= tol.confusion() {
                return vec![foot];
            }
            vec![add2(foot, scale2(d, w)), sub2(foot, scale2(d, w))]
        }
        (
            Profile::Circle {
                centre: c0,
                radius: r0,
            },
            Profile::Circle {
                centre: c1,
                radius: r1,
            },
        ) => {
            let d = sub2(c1, c0);
            let gap = norm2(d);
            if gap <= tol.confusion() || gap > r0 + r1 + tol.confusion() || gap < (r0 - r1).abs() {
                return Vec::new();
            }
            let along = r0.mul_add(r0, -(r1 * r1)) / (2.0 * gap) + gap / 2.0;
            let u = scale2(d, 1.0 / gap);
            let base = add2(c0, scale2(u, along));
            let w = (r0 * r0 - along * along).max(0.0).sqrt();
            if w <= tol.confusion() {
                return vec![base];
            }
            let v = perp2(u);
            vec![add2(base, scale2(v, w)), sub2(base, scale2(v, w))]
        }
    }
}

/// A support's profile in the layout's sections, where its trace is the
/// same line or circle in every one.
fn profile_of(surface: &SurfaceGeometry, layout: &Layout, tol: Tolerances) -> Option<Profile> {
    let near_axis = |p: Point, origin: Point, axis: Vector| {
        let d = p - origin;
        (d - axis * d.dot(axis)).magnitude() <= tol.confusion() * 10.0
    };
    match (*layout, surface) {
        (Layout::Extruded { along, .. }, SurfaceGeometry::Plane(plane)) => {
            let frame = plane.plane().frame();
            let n = frame.z().vector();
            if n.dot(along).abs() > tol.angular() {
                return None;
            }
            let normal = layout.flat_vector(n, 0.0);
            Some(Profile::Line {
                at: layout.flat(frame.origin()),
                normal: scale2(normal, 1.0 / norm2(normal)),
            })
        }
        (Layout::Extruded { along, .. }, SurfaceGeometry::Cylinder(c)) => {
            let cylinder = c.cylinder();
            let axis = cylinder.axis();
            if axis.direction.vector().cross(along).magnitude() > tol.angular() {
                return None;
            }
            Some(Profile::Circle {
                centre: layout.flat(axis.location),
                radius: cylinder.radius(),
            })
        }
        (Layout::Revolved { origin, axis, .. }, SurfaceGeometry::Plane(plane)) => {
            let frame = plane.plane().frame();
            let n = frame.z().vector();
            if n.cross(axis).magnitude() > tol.angular() {
                return None;
            }
            Some(Profile::Line {
                at: (0.0, (frame.origin() - origin).dot(axis)),
                normal: (0.0, n.dot(axis).signum()),
            })
        }
        (Layout::Revolved { origin, axis, .. }, SurfaceGeometry::Cylinder(c)) => {
            let cylinder = c.cylinder();
            let own = cylinder.axis();
            if own.direction.vector().cross(axis).magnitude() > tol.angular()
                || !near_axis(own.location, origin, axis)
            {
                return None;
            }
            Some(Profile::Line {
                at: (cylinder.radius(), 0.0),
                normal: (1.0, 0.0),
            })
        }
        (Layout::Revolved { origin, axis, .. }, SurfaceGeometry::Cone(c)) => {
            let cone = c.cone();
            let frame = cone.frame();
            let z = frame.z().vector();
            if z.cross(axis).magnitude() > tol.angular() || !near_axis(frame.origin(), origin, axis)
            {
                return None;
            }
            let rise = z.dot(axis).signum();
            let slope = cone.half_angle().tan();
            let length = slope.hypot(1.0);
            Some(Profile::Line {
                at: (cone.reference_radius(), (frame.origin() - origin).dot(axis)),
                normal: (rise / length, -slope / length),
            })
        }
        (Layout::Revolved { origin, axis, .. }, SurfaceGeometry::Sphere(s)) => {
            let sphere = s.sphere();
            if !near_axis(sphere.centre(), origin, axis) {
                return None;
            }
            Some(Profile::Circle {
                centre: (0.0, (sphere.centre() - origin).dot(axis)),
                radius: sphere.radius(),
            })
        }
        (Layout::Revolved { origin, axis, .. }, SurfaceGeometry::Torus(t)) => {
            let torus = t.torus();
            let frame = torus.frame();
            if frame.z().vector().cross(axis).magnitude() > tol.angular()
                || !near_axis(frame.origin(), origin, axis)
            {
                return None;
            }
            Some(Profile::Circle {
                centre: (torus.major_radius(), (frame.origin() - origin).dot(axis)),
                radius: torus.minor_radius(),
            })
        }
        _ => None,
    }
}

/// A surface's unit normal at the point of it nearest `p`.
fn surface_normal_near(
    surface: &SurfaceGeometry,
    p: Point,
    tol: Tolerances,
) -> OgeomResult<Vector> {
    let found = project_on_surface(surface, p, 24, tol)?;
    let (u, v) = found.parameters;
    let (du, dv) = surface.d1_at(u, v, tol)?;
    du.cross(dv).normalized(tol)
}

/// The sign that turns a profile's own normal into the face's, read at a
/// point `p` of the face where the face's normal is `normal`.
fn profile_sign(profile: &Profile, layout: &Layout, p: Point, normal: Vector) -> OgeomResult<f64> {
    let q = layout.flat(p);
    let along = dot2(
        layout.flat_vector(normal, layout.station(p)),
        profile.normal_at(q),
    );
    if along.abs() < 0.9 {
        ogeom_bail!(
            Invariant,
            "a support's normal does not lie in the round's section ({along})"
        );
    }
    Ok(along.signum())
}

/// The ball's section: its centre and where it touches each support.
#[derive(Debug, Clone, Copy)]
struct Section {
    centre: P2,
    contacts: [P2; 2],
}

impl Section {
    /// Which way the short arc from the first contact to the second turns
    /// about the centre, and how far.
    fn turn(&self) -> (f64, f64) {
        let a = sub2(self.contacts[0], self.centre);
        let b = sub2(self.contacts[1], self.centre);
        let c = cross2(a, b);
        (c.signum(), c.abs().atan2(dot2(a, b)))
    }

    /// The direction along each support away from the round: against the
    /// arc's tangent where it leaves that support.
    fn away(&self) -> [P2; 2] {
        let (sign, _) = self.turn();
        let a = sub2(self.contacts[0], self.centre);
        let b = sub2(self.contacts[1], self.centre);
        let into_first = scale2(perp2(a), sign);
        let into_second = scale2(perp2(b), -sign);
        [
            scale2(into_first, -1.0 / norm2(into_first)),
            scale2(into_second, -1.0 / norm2(into_second)),
        ]
    }
}

/// Every ball of `radius` on the given sides of two profiles.
fn seatings(
    profiles: [Profile; 2],
    sides: [f64; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Section>> {
    let (Some(a), Some(b)) = (
        profiles[0].offset(sides[0], radius, tol),
        profiles[1].offset(sides[1], radius, tol),
    ) else {
        ogeom_bail!(
            Construction,
            "a ball of radius {radius} does not fit inside a curved face it would roll on"
        );
    };
    let mut found = Vec::new();
    for centre in crossings(a, b, tol) {
        let contacts = [
            profiles[0].foot(centre, sides[0], radius),
            profiles[1].foot(centre, sides[1], radius),
        ];
        if norm2(sub2(contacts[0], contacts[1])) <= tol.confusion() {
            continue;
        }
        found.push(Section { centre, contacts });
    }
    Ok(found)
}

/// A seat with a closed form: the layout and the section every station
/// shares.
#[derive(Debug, Clone, Copy)]
struct Exact {
    layout: Layout,
    section: Section,
    radius: f64,
}

impl Exact {
    /// Check the section makes a round: the arc neither collapses nor
    /// spans half the ball, and a torus does not cross its own axis.
    fn checked(self, tol: Tolerances) -> OgeomResult<Self> {
        let (_, sweep) = self.section.turn();
        if sweep <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the ball touches the two faces at one point; there is no corner to round"
            );
        }
        if sweep >= PI - tol.angular() {
            ogeom_bail!(
                Construction,
                "the ball touches the two faces at the ends of one diameter; which half of it \
                 rounds the corner is undecided"
            );
        }
        if self.layout.turns() {
            if self.section.centre.0 - self.radius <= tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the ball's centre runs within its radius of the axis; the round would be a \
                     torus crossing its own axis"
                );
            }
            for contact in self.section.contacts {
                if contact.0 <= tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "the ball touches a face on the axis; its line of contact has no length"
                    );
                }
            }
        }
        Ok(self)
    }

    fn contact(&self, f: usize, s: f64) -> Point {
        self.layout.lift(self.section.contacts[f], s)
    }

    fn centre(&self, s: f64) -> Point {
        self.layout.lift(self.section.centre, s)
    }

    /// The line of contact on the `f`th support, parameterized by station.
    fn contact_curve(&self, f: usize, tol: Tolerances) -> OgeomResult<Curve> {
        self.contact_curve_from(f, 0.0, tol)
    }

    /// The line of contact on the `f`th support, parameterized by station
    /// less `from`: a circle's parameter starts at the station `from`.
    fn contact_curve_from(&self, f: usize, from: f64, tol: Tolerances) -> OgeomResult<Curve> {
        let q = self.section.contacts[f];
        match self.layout {
            Layout::Extruded { along, .. } => Ok(Curve::Line(LineCurve::new(Axis {
                location: self.layout.lift(q, from),
                direction: Direction::new(along, tol)?,
            }))),
            Layout::Revolved {
                origin, axis, x, y, ..
            } => {
                let frame = Frame::new(
                    origin + axis * q.1,
                    Direction::new(axis, tol)?,
                    Direction::new(Layout::radial(x, y, from), tol)?,
                    tol,
                )?;
                Ok(Curve::Circle(ogeom_geom::CircleCurve::new(Circle::new(
                    frame, q.0, tol,
                )?)))
            }
        }
    }

    /// The line of contact on the `f`th support over the stations `run`,
    /// as a curve and its range.
    fn contact_run(
        &self,
        f: usize,
        run: (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<(Curve, (f64, f64))> {
        Ok((
            self.contact_curve_from(f, run.0, tol)?,
            (0.0, run.1 - run.0),
        ))
    }

    /// The ball's arc at station `s`, from the first contact to the second,
    /// over its range.
    fn arc(&self, s: f64, tol: Tolerances) -> OgeomResult<(Curve, (f64, f64))> {
        let centre = self.centre(s);
        let a = self.contact(0, s) - centre;
        let b = self.contact(1, s) - centre;
        let (_, sweep) = self.section.turn();
        let frame = Frame::new(
            centre,
            Direction::new(a.cross(b), tol)?,
            Direction::new(a, tol)?,
            tol,
        )?;
        let circle = Circle::new(frame, self.radius, tol)?;
        Ok((
            Curve::Circle(ogeom_geom::CircleCurve::new(circle)),
            (0.0, sweep),
        ))
    }

    /// The round's surface, its window over the stations `span` where it
    /// has one.
    fn surface(&self, span: (f64, f64), tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
        let c = self.section.centre;
        match self.layout {
            Layout::Extruded { along, x, .. } => {
                let frame = Frame::new(
                    self.layout.lift(c, 0.0),
                    Direction::new(along, tol)?,
                    Direction::new(x, tol)?,
                    tol,
                )?;
                Ok(CylinderSurface::new(Cylinder::new(frame, self.radius, tol)?, span)?.into())
            }
            Layout::Revolved {
                origin, axis, x, ..
            } => {
                let frame = Frame::new(
                    origin + axis * c.1,
                    Direction::new(axis, tol)?,
                    Direction::new(x, tol)?,
                    tol,
                )?;
                Ok(TorusSurface::new(Torus::new(frame, c.0, self.radius, tol)?).into())
            }
        }
    }

    /// How far each face's strip reaches along it from the corner `crease`
    /// (in section coordinates) to its line of contact.
    fn strip_width(&self, profile: &Profile, f: usize, crease: P2) -> f64 {
        profile.length_between(crease, self.section.contacts[f])
    }
}

/// A seat that has to be marched, with the band fitted through it.
struct Marched {
    blend: MarchedBlend,
    band: ogeom_geom::BSplineSurface,
    fit_target: f64,
}

impl Marched {
    /// The band's corner on the `f`th border at the `k`th end: control
    /// points the clamped borders interpolate.
    fn corner(&self, f: usize, k: usize) -> Point {
        let grid = self.band.grid();
        let (nu, nv) = (grid.u_count(), grid.v_count());
        let i = if f == 0 { 0 } else { nu - 1 };
        let j = if k == 0 { 0 } else { nv - 1 };
        grid.points()[i * nv + j].point()
    }

    /// The band's border along the `f`th line of contact, or its row at the
    /// `k`th end, as curves off the control net.
    fn border(&self, f: usize, tol: Tolerances) -> OgeomResult<Curve> {
        let grid = self.band.grid();
        let (nu, nv) = (grid.u_count(), grid.v_count());
        let i = if f == 0 { 0 } else { nu - 1 };
        let control: Vec<Point> = (0..nv).map(|j| grid.points()[i * nv + j].point()).collect();
        Ok(Curve::BSpline(ogeom_geom::BSplineCurve::new(
            self.band.v_knots().clone(),
            control,
            tol,
        )?))
    }

    fn row(&self, k: usize, tol: Tolerances) -> OgeomResult<Curve> {
        let grid = self.band.grid();
        let (nu, nv) = (grid.u_count(), grid.v_count());
        let j = if k == 0 { 0 } else { nv - 1 };
        let control: Vec<Point> = (0..nu).map(|i| grid.points()[i * nv + j].point()).collect();
        Ok(Curve::BSpline(ogeom_geom::BSplineCurve::new(
            self.band.u_knots().clone(),
            control,
            tol,
        )?))
    }
}

/// The seat of one round.
enum Seat {
    Exact(Exact, [Profile; 2]),
    Marched(Marched),
}

impl Seat {
    /// Where the ball touches the `f`th face at the `k`th end of the run
    /// from station `span.0` to `span.1`.
    fn end_contact(&self, f: usize, k: usize, span: (f64, f64)) -> Point {
        match self {
            Self::Exact(exact, _) => exact.contact(f, if k == 0 { span.0 } else { span.1 }),
            Self::Marched(marched) => marched.corner(f, k),
        }
    }

    /// How far a point may sit from the boundary it is said to lie on.
    fn reach(&self, tol: Tolerances) -> f64 {
        match self {
            Self::Exact(..) => tol.confusion() * 10.0,
            Self::Marched(marched) => marched.fit_target,
        }
    }
}

/// One face beside the rounded edge, as the rebuild reads it.
struct Beside {
    face: Shape,
    surface_id: SurfaceId,
    surface: SurfaceGeometry,
    /// The face's normal at the probe point on the edge.
    normal: Vector,
    /// The direction, square to the edge in the face, the face lies in.
    inward: Vector,
    /// The face's wire holding the edge, and where in it.
    wire: usize,
    at: usize,
    loops: Vec<Vec<Shape>>,
    tolerance: Tolerance,
}

/// A face's surface id, surface and tolerance; refused where it is placed.
fn unplaced_face(
    model: &Model,
    face: &Shape,
) -> OgeomResult<(SurfaceId, SurfaceGeometry, Tolerance)> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    if !face.location().is_identity() || !data.location.is_identity() {
        ogeom_bail!(
            Construction,
            "a placed face is not rounded; bake its placement into its geometry first"
        );
    }
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    Ok((data.surface, surface.clone(), data.tolerance))
}

/// The face's oriented normal at the point of its surface nearest `p`.
fn face_normal_at(
    surface: &SurfaceGeometry,
    face: &Shape,
    p: Point,
    tol: Tolerances,
) -> OgeomResult<Vector> {
    let n = surface_normal_near(surface, p, tol)?;
    Ok(if face.orientation() == Orientation::Reversed {
        -n
    } else {
        n
    })
}

/// Which way, square to the edge with tangent `tangent` at `p`, the face
/// lies: probed on the face itself at growing distances.
#[allow(clippy::too_many_arguments, reason = "one probe, all its data")]
fn probe_inward(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    p: Point,
    tangent: Vector,
    normal: Vector,
    span: f64,
    tol: Tolerances,
) -> OgeomResult<Vector> {
    let raw = normal.cross(tangent).normalized(tol)?;
    for scale in [1e-3, 1e-2, 5e-2] {
        let eps = span * scale;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        for dir in [raw, -raw] {
            let probe = project_on_surface(surface, p + dir * eps, 32, tol)?.point;
            if ogeom_algo::classify_on_face(model, face, probe, deflection, tol)?
                == ogeom_algo::Containment::In
            {
                return Ok(dir);
            }
        }
    }
    ogeom_bail!(
        Construction,
        "cannot read which way a face extends from the edge; the face is thinner than the \
         probe can resolve"
    );
}

/// The pcurve an edge occurrence carries on `surface`, as stored: one
/// image, or a seam's two.
enum Image {
    One(PlanarCurve),
    Seam(PlanarCurve, PlanarCurve),
}

fn image_on(model: &Model, edge: &Shape, surface: SurfaceId) -> OgeomResult<Option<Image>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let get = |id| model.geometry().pcurve(id).cloned();
    Ok(match data.pcurve_for(surface, edge.location()) {
        Some(EdgeRepr::PCurve { curve, .. }) => get(*curve).map(Image::One),
        Some(EdgeRepr::Seam {
            forward, reversed, ..
        }) => match (get(*forward), get(*reversed)) {
            (Some(f), Some(r)) => Some(Image::Seam(f, r)),
            _ => None,
        },
        _ => None,
    })
}

/// The image an occurrence uses: a seam's forward image for a forward
/// occurrence, its reversed one otherwise.
fn occurrence_image(model: &Model, edge: &Shape, surface: SurfaceId) -> OgeomResult<PlanarCurve> {
    match image_on(model, edge, surface)? {
        Some(Image::One(c)) => Ok(c),
        Some(Image::Seam(f, r)) => Ok(if edge.orientation() == Orientation::Reversed {
            r
        } else {
            f
        }),
        None => ogeom_bail!(
            Construction,
            "a boundary edge of a face has no image in the face's chart"
        ),
    }
}

/// Mark an edge's descriptions as sharing one parameter.
fn same_parameter(model: &mut Model, edge: &Shape) {
    if let Some(NodeData::Edge(data)) = model.node_mut(edge).map(|n| n.data_mut()) {
        data.assert_same_parameter(true);
    }
}

/// An edge on the curve of `old` (an occurrence, placed as it stands) over
/// `range` of its parameter, between `from` and `to`, carrying the
/// occurrence's own image on `surface` over the same range.
fn piece_of(
    model: &mut Model,
    old: &Shape,
    surface: SurfaceId,
    range: (f64, f64),
    from: &Shape,
    to: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let forward = old.oriented(Orientation::Forward);
    let (curve, _) = edge_curve(model, &forward, tol)?;
    let image = image_on(model, &forward, surface)?;
    let tolerance = model
        .node(old)
        .and_then(|n| n.data().as_edge())
        .map_or(Tolerance::MIN, |d| d.tolerance);
    let edge = make_edge_between(model, curve, range, from, to, tol)?.shape;
    match image {
        Some(Image::One(c)) => {
            attach_pcurve(model, &edge, c, surface, Location::identity(), range)?;
        }
        Some(Image::Seam(f, r)) => {
            attach_seam(model, &edge, f, r, surface, Location::identity(), range)?;
        }
        None => ogeom_bail!(
            Construction,
            "a boundary edge of a face has no image in the face's chart"
        ),
    }
    same_parameter(model, &edge);
    model.widen(&edge, tolerance)?;
    Ok(edge)
}

/// The image of `curve` over `range` on `surface`: exact where it has a
/// closed form, fitted otherwise with the edge and its vertices widened to
/// hold the fit; put on the branch of a periodic chart nearest `near` at
/// the image's start.
fn attach_image(
    model: &mut Model,
    edge: &Shape,
    surface_id: SurfaceId,
    surface: &SurfaceGeometry,
    near: Option<Point2>,
    tol: Tolerances,
) -> OgeomResult<()> {
    let (curve, range) = edge_curve(model, edge, tol)?;
    let mut image =
        if let Some(exact) = ogeom_intersect::exact_pcurve_over(&curve, range, surface, tol) {
            exact
        } else {
            let (fitted, _, _, worst, _) =
                ogeom_algo::pcurve_fit::fit_projected_pcurve(&curve, range, surface, tol)?;
            if worst > tol.confusion() {
                let widened = Tolerance::new(worst + tol.confusion())?;
                model.widen(edge, widened)?;
                if let Some((a, b)) = edge_vertices(model, edge)? {
                    model.widen(&a, widened)?;
                    model.widen(&b, widened)?;
                }
            }
            fitted
        };
    if let Some(target) = near {
        let ((ua, ub), (va, vb)) = surface.domain();
        let start = image.point_at(range.0, tol)?;
        let shift = |periodic: bool, span: f64, d: f64| {
            if periodic && span > 0.0 {
                (d / span).round() * span
            } else {
                0.0
            }
        };
        let offset = Vector2::new(
            shift(surface.is_periodic_u(), ub - ua, target.x - start.x),
            shift(surface.is_periodic_v(), vb - va, target.y - start.y),
        );
        if offset.x != 0.0 || offset.y != 0.0 {
            image = image.transformed(&Transform2::translation(offset), tol)?;
        }
    }
    attach_pcurve(model, edge, image, surface_id, Location::identity(), range)?;
    same_parameter(model, edge);
    Ok(())
}

/// A face again on its own surface, every edge of `swap` replaced by its
/// pieces (forward edges, in the replaced edge's stored order), each
/// occurrence keeping its orientation.
fn rebuilt_face(
    model: &mut Model,
    face: &Shape,
    surface_id: SurfaceId,
    tolerance: Tolerance,
    swap: &[(Shape, Vec<Shape>)],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let wires = model.children_of(&face.oriented(Orientation::Forward))?;
    let mut rebuilt = Vec::with_capacity(wires.len());
    for wire in &wires {
        let ring = model.children_of(&wire.oriented(Orientation::Forward))?;
        let mut edges = Vec::with_capacity(ring.len() + 2);
        let mut touched = false;
        for edge in ring {
            match swap.iter().find(|(old, _)| old.is_same(&edge)) {
                Some((_, pieces)) => {
                    touched = true;
                    if edge.orientation() == Orientation::Reversed {
                        edges.extend(pieces.iter().rev().map(Shape::reversed));
                    } else {
                        edges.extend(pieces.iter().cloned());
                    }
                }
                None => edges.push(edge),
            }
        }
        if touched {
            rebuilt.push(
                make_wire(model, &edges, tol)?
                    .shape
                    .oriented(wire.orientation()),
            );
        } else {
            rebuilt.push(wire.clone());
        }
    }
    let mut data = FaceData::new(surface_id, Location::identity());
    data.tolerance = tolerance;
    Ok(model.add_face(data, &rebuilt)?.oriented(face.orientation()))
}

/// Orient a round so its normal agrees with the faces it meets: away from
/// the ball where `away` holds, toward it otherwise. `ball` finds the
/// centre of the ball touching the round at a point.
fn oriented_round(
    model: &Model,
    round: Shape,
    away: bool,
    ball: impl Fn(Point) -> Point,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (p, n) = ogeom_algo::face_normal(model, &round, tol)?;
    let outward = n.dot(p - ball(p)) > 0.0;
    Ok(if outward == away {
        round
    } else {
        round.reversed()
    })
}

/// Sampled distance from every point of an edge to a curve, held low by
/// half the longest step between samples: a lower bound on the true least
/// distance, the distance changing no faster than the walk along the edge.
fn least_distance(model: &Model, edge: &Shape, to: &Curve, tol: Tolerances) -> OgeomResult<f64> {
    const SAMPLES: usize = 64;
    let (curve, (lo, hi)) = edge_curve(model, edge, tol)?;
    let mut least = f64::INFINITY;
    let mut step = 0.0_f64;
    let mut previous: Option<Point> = None;
    for k in 0..=SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * (k as f64) / (SAMPLES as f64);
        let p = curve.point_at(t, tol)?;
        least = least.min(project_on_curve(to, p, 32, tol)?.distance);
        if let Some(q) = previous {
            step = step.max(p.distance(q));
        }
        previous = Some(p);
    }
    // A chord is shorter than its arc; a tenth more covers the bend
    // between samples of any edge smooth enough to round against.
    Ok(least - step * 0.55)
}

/// Round one edge of a shell where a face beside it is curved, or the
/// edge is.
#[allow(clippy::too_many_lines, reason = "one rebuild, checked then assembled")]
pub(crate) fn round_sheet_edge(
    model: &mut Model,
    shell: &Shape,
    edge: &Shape,
    faces: [&Shape; 2],
    users_of: &dyn Fn(&Shape) -> usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (curve, range) = edge_curve(model, edge, tol)?;
    let Some((v0, v1)) = edge_vertices(model, &edge.oriented(Orientation::Forward))? else {
        ogeom_bail!(Construction, "the edge has no vertices");
    };
    let closed = v0.is_same(&v1);
    let end_vertices = [v0.clone(), v1.clone()];
    let mid = f64::midpoint(range.0, range.1);
    let p = curve.point_at(mid, tol)?;
    let tangent = curve.d1_at(mid, tol)?.normalized(tol)?;
    let span = {
        let a = curve.point_at(range.0, tol)?;
        let b = curve.point_at(range.1, tol)?;
        a.distance(p).max(b.distance(p)).max(radius)
    };

    let mut sides = Vec::with_capacity(2);
    for face in faces {
        let (surface_id, surface, tolerance) = unplaced_face(model, face)?;
        let normal = face_normal_at(&surface, face, p, tol)?;
        let inward = probe_inward(model, face, &surface, p, tangent, normal, span, tol)?;
        let wires = model.children_of(&face.oriented(Orientation::Forward))?;
        let mut loops = Vec::with_capacity(wires.len());
        for wire in &wires {
            loops.push(model.children_of(&wire.oriented(Orientation::Forward))?);
        }
        let mut hits = Vec::new();
        for (w, ring) in loops.iter().enumerate() {
            for (i, e) in ring.iter().enumerate() {
                if e.is_same(edge) {
                    hits.push((w, i));
                }
            }
        }
        let [(wire, at)] = hits.as_slice() else {
            ogeom_bail!(
                Construction,
                "a face's boundary runs along the edge twice (a seam); it has no one side to \
                 round"
            );
        };
        sides.push(Beside {
            face: face.clone(),
            surface_id,
            surface,
            normal,
            inward,
            wire: *wire,
            at: *at,
            loops,
            tolerance,
        });
    }

    // The corner is concave on the side each face's normal leans toward
    // the other face, or on the other side for both.
    let toward = [
        sides[0].normal.dot(sides[1].inward),
        sides[1].normal.dot(sides[0].inward),
    ];
    if toward[0].abs() <= tol.angular() && toward[1].abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the two faces continue each other across the edge; there is no corner to round"
        );
    }
    if (toward[0] > 0.0) != (toward[1] > 0.0) {
        ogeom_bail!(
            Construction,
            "the two faces are oriented inconsistently across the edge; the round has no \
             side to face"
        );
    }
    let concave = toward[0].signum();

    // The seat: a closed form where the layout fits both faces, marched
    // otherwise.
    let layout = match &curve {
        Curve::Line(line) => Some(Layout::extruded(
            line.axis().location,
            line.axis().direction.vector(),
            tol,
        )?),
        Curve::Circle(c) => {
            let frame = c.circle().frame();
            let (x, y) = (frame.x().vector(), frame.y().vector());
            Some(Layout::revolved(c.circle().centre(), x.cross(y), x, tol)?)
        }
        _ => None,
    };
    let exact = match layout {
        Some(layout) => match (
            profile_of(&sides[0].surface, &layout, tol),
            profile_of(&sides[1].surface, &layout, tol),
        ) {
            (Some(a), Some(b)) => Some((layout, [a, b])),
            _ => None,
        },
        None => None,
    };
    let seat = if let Some((layout, profiles)) = exact {
        let crease = layout.flat(p);
        let station = layout.station(p);
        let mut ball_sides = [0.0; 2];
        let mut inward = [(0.0, 0.0); 2];
        for f in 0..2 {
            ball_sides[f] = concave * profile_sign(&profiles[f], &layout, p, sides[f].normal)?;
            inward[f] = layout.flat_vector(sides[f].inward, station);
        }
        let chosen = seatings(profiles, ball_sides, radius, tol)?
            .into_iter()
            .filter(|s| {
                (0..2).all(|f| dot2(sub2(s.contacts[f], crease), inward[f]) > tol.confusion())
            })
            .min_by(|a, b| norm2(sub2(a.centre, crease)).total_cmp(&norm2(sub2(b.centre, crease))));
        let Some(section) = chosen else {
            ogeom_bail!(
                Construction,
                "no ball of radius {radius} seats in the corner between the two faces"
            );
        };
        let exact = Exact {
            layout,
            section,
            radius,
        }
        .checked(tol)?;
        for (f, side) in sides.iter().enumerate() {
            let q = exact.contact(f, station);
            if project_on_surface(&side.surface, q, 24, tol)?.distance > tol.confusion() * 10.0 {
                ogeom_bail!(Invariant, "a closed-form contact misses its face's surface");
            }
        }
        Seat::Exact(exact, profiles)
    } else {
        if closed {
            ogeom_bail!(
                Construction,
                "a closed edge is rounded where its round has a closed form only (faces turning \
                 about one axis); this one would have to be marched"
            );
        }
        Seat::Marched(marched_seat(&sides, &curve, range, concave, radius, tol)?)
    };
    // The stations of the edge's two ends, the second reached from the
    // first the way the edge runs: a whole turn on for a closed edge.
    let stations = match &seat {
        Seat::Exact(exact, _) => {
            let layout = exact.layout;
            let first = curve.point_at(range.0, tol)?;
            let s0 = layout.station(first);
            let s1 = layout.station(curve.point_at(range.1, tol)?);
            match layout {
                Layout::Revolved { axis, x, y, .. } => {
                    let around = axis.cross(Layout::radial(x, y, s0));
                    let rising = curve.d1_at(range.0, tol)?.dot(around) > 0.0;
                    let swept = match (closed, rising) {
                        (true, true) => TAU,
                        (true, false) => -TAU,
                        (false, true) => (s1 - s0).rem_euclid(TAU),
                        (false, false) => -(s0 - s1).rem_euclid(TAU),
                    };
                    (s0, s0 + swept)
                }
                Layout::Extruded { .. } => (s0, s1),
            }
        }
        Seat::Marched(_) => range,
    };

    // The faces' boundaries at the ends: the edge leading in and the edge
    // leading out of the rounded one in each face's loop, cut where the
    // ball touches them.
    struct Cut {
        old: Shape,
        surface: SurfaceId,
        start: Option<(usize, usize)>,
        end: Option<(usize, usize)>,
    }
    let mut cuts: Vec<Cut> = Vec::new();
    let mut lead_outs: Vec<Shape> = Vec::new();
    for (f, side) in sides.iter().enumerate() {
        let ring = &side.loops[side.wire];
        let n = ring.len();
        if n == 1 {
            continue;
        }
        let occurrence = &ring[side.at];
        let k_in = usize::from(occurrence.orientation() == Orientation::Reversed);
        let k_out = if closed { k_in } else { 1 - k_in };
        for (lead, k) in [
            (&ring[(side.at + n - 1) % n], k_in),
            (&ring[(side.at + 1) % n], k_out),
        ] {
            if users_of(lead) != 1 {
                ogeom_bail!(
                    Construction,
                    "the round would run into another face of the sheet at an end of the edge; \
                     only an edge whose ends lie on the sheet's free boundary is rounded"
                );
            }
            let Some((s, e)) = edge_vertices(model, &lead.oriented(Orientation::Forward))? else {
                ogeom_bail!(Construction, "a boundary edge has no vertices");
            };
            let (at_start, at_end) = (s.is_same(&end_vertices[k]), e.is_same(&end_vertices[k]));
            if at_start == at_end {
                ogeom_bail!(
                    Construction,
                    "a face's boundary edge leaving an end of the edge closes on itself there"
                );
            }
            let forward = lead.oriented(Orientation::Forward);
            if !lead_outs.iter().any(|l| l.is_same(&forward)) {
                lead_outs.push(forward.clone());
            }
            let slot = match cuts.iter().position(|c| c.old.is_same(&forward)) {
                Some(slot) => slot,
                None => {
                    cuts.push(Cut {
                        old: forward.clone(),
                        surface: side.surface_id,
                        start: None,
                        end: None,
                    });
                    cuts.len() - 1
                }
            };
            if at_start {
                cuts[slot].start = Some((f, k));
            } else {
                cuts[slot].end = Some((f, k));
            }
        }
    }
    // Nothing else meets the edge's ends.
    for (k, vertex) in end_vertices.iter().enumerate() {
        if closed && k == 1 {
            break;
        }
        let mut meeting: Vec<Shape> = Vec::new();
        for edge_of in explore_unique(model, shell, ShapeType::Edge)? {
            if let Some((a, b)) = edge_vertices(model, &edge_of)?
                && (a.is_same(vertex) || b.is_same(vertex))
                && !meeting.iter().any(|m| m.is_same(&edge_of))
            {
                meeting.push(edge_of);
            }
        }
        let expected = 1 + lead_outs
            .iter()
            .filter(|l| {
                edge_vertices(model, l)
                    .ok()
                    .flatten()
                    .is_some_and(|(a, b)| a.is_same(vertex) || b.is_same(vertex))
            })
            .count();
        if meeting.len() != expected {
            ogeom_bail!(
                Construction,
                "{} edges of the sheet meet at an end of the edge; only an edge whose ends lie \
                 on the sheet's free boundary, with one boundary edge of each face leaving it, \
                 is rounded",
                meeting.len()
            );
        }
    }
    // Each cut lands on its edge, inside it.
    let reach = seat.reach(tol);
    let mut cut_at: Vec<(Option<f64>, Option<f64>)> = Vec::with_capacity(cuts.len());
    for cut in &cuts {
        let (lead_curve, lead_range) = edge_curve(model, &cut.old, tol)?;
        let bounded = Curve::Trimmed(Box::new(TrimmedCurve::new(
            lead_curve,
            lead_range.0,
            lead_range.1,
            tol,
        )?));
        let mut found = (None, None);
        for (slot, which) in [(cut.start, 0), (cut.end, 1)] {
            let Some((f, k)) = slot else { continue };
            let point = seat.end_contact(f, k, stations);
            let hit = project_on_curve(&bounded, point, 64, tol)?;
            let far = if which == 0 {
                lead_range.1
            } else {
                lead_range.0
            };
            let near = if which == 0 {
                lead_range.0
            } else {
                lead_range.1
            };
            let length = (lead_range.1 - lead_range.0).abs();
            let reached = (hit.parameter - near).abs();
            if reached >= length - tol.parametric().max(length * 1e-9)
                || (hit.parameter - far).abs() <= tol.parametric()
            {
                ogeom_bail!(
                    Construction,
                    "the round sets back past the far end of a face's boundary edge leaving the \
                     edge"
                );
            }
            if hit.distance > reach {
                ogeom_bail!(
                    Construction,
                    "a face's boundary leaves an end of the edge outside the round's end \
                     section ({} from where the ball touches); the round ends in its section, \
                     so only a boundary running along it is rounded",
                    hit.distance
                );
            }
            if which == 0 {
                found.0 = Some(hit.parameter);
            } else {
                found.1 = Some(hit.parameter);
            }
        }
        if let (Some(a), Some(b)) = found
            && a >= b
        {
            ogeom_bail!(
                Construction,
                "the rounds at the two ends of the edge cross on a face's boundary edge"
            );
        }
        cut_at.push(found);
    }
    // Nothing else of either face comes within the strip the round takes.
    let bounded_edge = Curve::Trimmed(Box::new(TrimmedCurve::new(
        curve.clone(),
        range.0,
        range.1,
        tol,
    )?));
    for (f, side) in sides.iter().enumerate() {
        let width = match &seat {
            Seat::Exact(exact, profiles) => {
                exact.strip_width(&profiles[f], f, exact.layout.flat(p))
            }
            Seat::Marched(marched) => marched_strip_width(marched, side, f, &curve, tol)?,
        };
        for ring in &side.loops {
            for other in ring {
                let forward = other.oriented(Orientation::Forward);
                if other.is_same(edge) || lead_outs.iter().any(|l| l.is_same(&forward)) {
                    continue;
                }
                if least_distance(model, other, &bounded_edge, tol)?
                    <= width + reach + tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "another part of a face's boundary comes within the strip the round \
                         takes from it; the round could cross it"
                    );
                }
            }
        }
    }

    // Assembly: the contact vertices, the lines of contact with their
    // images on each face, the round, the cut boundary edges, the faces
    // rebuilt around them.
    let mut corners: Vec<[Shape; 2]> = Vec::with_capacity(2);
    for f in 0..2 {
        let a = make_vertex(model, seat.end_contact(f, 0, stations)).shape;
        let b = if closed {
            a.clone()
        } else {
            make_vertex(model, seat.end_contact(f, 1, stations)).shape
        };
        if let Seat::Marched(marched) = &seat {
            let widened = Tolerance::new(marched.fit_target)?;
            model.widen(&a, widened)?;
            model.widen(&b, widened)?;
        }
        corners.push([a, b]);
    }
    let corners = [corners[0].clone(), corners[1].clone()];
    // Each rail is stored running up the stations; `along_edge` is it as
    // the rounded edge runs, which is how it takes the edge's place.
    let rising = stations.0 <= stations.1;
    let run = if rising {
        stations
    } else {
        (stations.1, stations.0)
    };
    let mut rails = Vec::with_capacity(2);
    let mut along_edge = Vec::with_capacity(2);
    for (f, side) in sides.iter().enumerate() {
        let (from, to) = if rising { (0, 1) } else { (1, 0) };
        let (rail_curve, rail_range) = match &seat {
            Seat::Exact(exact, _) => exact.contact_run(f, run, tol)?,
            Seat::Marched(marched) => {
                let border = marched.border(f, tol)?;
                let domain = border.domain();
                (border, domain)
            }
        };
        let rail = make_edge_between(
            model,
            rail_curve,
            rail_range,
            &corners[f][from],
            &corners[f][to],
            tol,
        )?
        .shape;
        if let Seat::Marched(marched) = &seat {
            model.widen(&rail, Tolerance::new(marched.fit_target)?)?;
        }
        let edge_image = occurrence_image(model, &side.loops[side.wire][side.at], side.surface_id)?;
        let near = edge_image.point_at(if rising { range.0 } else { range.1 }, tol)?;
        attach_image(
            model,
            &rail,
            side.surface_id,
            &side.surface,
            Some(near),
            tol,
        )?;
        along_edge.push(if rising {
            rail.clone()
        } else {
            rail.reversed()
        });
        rails.push(rail);
    }
    let rails = [rails[0].clone(), rails[1].clone()];
    let ordered = if rising {
        corners.clone()
    } else {
        [
            [corners[0][1].clone(), corners[0][0].clone()],
            [corners[1][1].clone(), corners[1][0].clone()],
        ]
    };
    let away = concave < 0.0;
    let round = match &seat {
        Seat::Exact(exact, _) => {
            if closed {
                let surface = exact.surface(run, tol)?;
                let band =
                    ogeom_algo::make_revolution_band(model, &surface, &rails[0], &rails[1], tol)?;
                oriented_round(
                    model,
                    band,
                    away,
                    |q| exact.centre(exact.layout.station(q)),
                    tol,
                )?
            } else {
                exact_open_round(model, exact, run, &rails, &ordered, away, tol)?
            }
        }
        Seat::Marched(marched) => marched_round(model, marched, &rails, &corners, away, tol)?,
    };

    let mut history = History::new();
    let mut swaps: Vec<(Shape, Vec<Shape>)> = Vec::with_capacity(cuts.len() + 1);
    for (cut, found) in cuts.iter().zip(&cut_at) {
        let (_, lead_range) = edge_curve(model, &cut.old, tol)?;
        let Some((s, e)) = edge_vertices(model, &cut.old)? else {
            ogeom_bail!(Construction, "a boundary edge has no vertices");
        };
        let (from, lo) = match (cut.start, found.0) {
            (Some((f, k)), Some(t)) => (corners[f][k].clone(), t),
            _ => (s, lead_range.0),
        };
        let (to, hi) = match (cut.end, found.1) {
            (Some((f, k)), Some(t)) => (corners[f][k].clone(), t),
            _ => (e, lead_range.1),
        };
        let piece = piece_of(model, &cut.old, cut.surface, (lo, hi), &from, &to, tol)?;
        history.modify(&cut.old, piece.clone());
        swaps.push((cut.old.clone(), vec![piece]));
    }
    let mut replaced = Vec::with_capacity(2);
    for (f, side) in sides.iter().enumerate() {
        let mut swap = swaps.clone();
        swap.push((
            edge.oriented(Orientation::Forward),
            vec![along_edge[f].clone()],
        ));
        let face = rebuilt_face(
            model,
            &side.face,
            side.surface_id,
            side.tolerance,
            &swap,
            tol,
        )?;
        if side.surface.is_periodic_u() || side.surface.is_periodic_v() {
            let wires = model.children_of(&face.oriented(Orientation::Forward))?;
            ogeom_algo::chain_wire_branches(model, side.surface_id, &wires, tol)?;
        }
        history.modify(&side.face, face.clone());
        replaced.push((side.face.clone(), face));
    }
    let faces_now = explore_unique(model, shell, ShapeType::Face)?;
    let mut shell_faces: Vec<Shape> = faces_now
        .iter()
        .map(|face| {
            replaced
                .iter()
                .find(|(old, _)| old.is_same(face))
                .map_or_else(|| face.clone(), |(_, new)| new.clone())
        })
        .collect();
    shell_faces.push(round.clone());
    let rounded = model.add_shell(&shell_faces)?;
    history.modify(shell, rounded.clone());
    history.delete(edge);
    history.delete(&v0);
    if !closed {
        history.delete(&v1);
    }
    history.generate(edge, round);
    Ok(Built::new(rounded, history))
}

/// The round over an open run of an exact seat from station `span.0` to
/// `span.1`: its two lines of contact (`rails`, each running from the
/// first station to the last) and the ball's arcs at the ends.
fn exact_open_round(
    model: &mut Model,
    exact: &Exact,
    span: (f64, f64),
    rails: &[Shape; 2],
    corners: &[[Shape; 2]; 2],
    away: bool,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut arcs = Vec::with_capacity(2);
    for (k, s) in [span.0, span.1].into_iter().enumerate() {
        let (arc, arc_range) = exact.arc(s, tol)?;
        arcs.push(
            make_edge_between(model, arc, arc_range, &corners[0][k], &corners[1][k], tol)?.shape,
        );
    }
    let surface = exact.surface(span, tol)?;
    let edges = vec![
        arcs[0].clone(),
        rails[1].clone(),
        arcs[1].reversed(),
        rails[0].reversed(),
    ];
    let face = face_from_edges(model, surface, &edges, tol)?;
    oriented_round(
        model,
        face,
        away,
        |q| exact.centre(exact.layout.station(q)),
        tol,
    )
}

/// The marched seat between the two faces along the edge's curve over
/// `range`, with exact sections at both ends.
fn marched_seat(
    sides: &[Beside],
    curve: &Curve,
    range: (f64, f64),
    concave: f64,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Marched> {
    let guide = Curve::Trimmed(Box::new(TrimmedCurve::new(
        curve.clone(),
        range.0,
        range.1,
        tol,
    )?));
    let sign = |side: &Beside| -> i8 {
        let flipped = side.face.orientation() == Orientation::Reversed;
        if (concave > 0.0) != flipped { 1 } else { -1 }
    };
    let ball_sides = Sides {
        first: sign(&sides[0]),
        second: sign(&sides[1]),
    };
    let (first, second) = (&sides[0].surface, &sides[1].surface);
    let walked = march_blend_sided(
        first,
        second,
        radius,
        &guide,
        ball_sides,
        Marching::default(),
        tol,
    )?;
    if walked.len() < 2 {
        ogeom_bail!(
            NotDone,
            "the rolling ball's march along the edge gave too few stations to fit a round"
        );
    }
    let mut order: Vec<usize> = (0..walked.len()).collect();
    order.sort_by(|&a, &b| walked.along[a].total_cmp(&walked.along[b]));
    let margin = (range.1 - range.0) * 1e-6;
    let inner: Vec<usize> = order
        .into_iter()
        .filter(|&i| walked.along[i] > range.0 + margin && walked.along[i] < range.1 - margin)
        .collect();
    let (Some(&head), Some(&tail)) = (inner.first(), inner.last()) else {
        ogeom_bail!(
            NotDone,
            "the rolling ball's march stayed at the edge's ends; there is nothing to fit"
        );
    };
    let near = |i: usize| -> [f64; 4] {
        [
            walked.on_first[i].0,
            walked.on_first[i].1,
            walked.on_second[i].0,
            walked.on_second[i].1,
        ]
    };
    let mut blend = MarchedBlend {
        spine: Vec::with_capacity(inner.len() + 2),
        on_first: Vec::with_capacity(inner.len() + 2),
        on_second: Vec::with_capacity(inner.len() + 2),
        touch_first: Vec::with_capacity(inner.len() + 2),
        touch_second: Vec::with_capacity(inner.len() + 2),
        along: Vec::with_capacity(inner.len() + 2),
        sides: ball_sides,
        stopped: walked.stopped,
    };
    let push_exact = |blend: &mut MarchedBlend, at: f64, seed: [f64; 4]| -> OgeomResult<()> {
        let x = seat_section(first, second, radius, &guide, ball_sides, at, seed, tol).map_err(
            |_| {
                ogeom_core::ogeom_err!(
                    Construction,
                    "the ball does not seat between the faces at an end of the edge"
                )
            },
        )?;
        let p1 = first.point_at(x[0], x[1], tol)?;
        let p2 = second.point_at(x[2], x[3], tol)?;
        let (du, dv) = first.d1_at(x[0], x[1], tol)?;
        let n = du.cross(dv).normalized(tol)?;
        blend
            .spine
            .push(p1 + n * (f64::from(ball_sides.first) * radius));
        blend.on_first.push((x[0], x[1]));
        blend.on_second.push((x[2], x[3]));
        blend.touch_first.push(p1);
        blend.touch_second.push(p2);
        blend.along.push(x[4]);
        Ok(())
    };
    push_exact(&mut blend, range.0, near(head))?;
    for &i in &inner {
        blend.spine.push(walked.spine[i]);
        blend.on_first.push(walked.on_first[i]);
        blend.on_second.push(walked.on_second[i]);
        blend.touch_first.push(walked.touch_first[i]);
        blend.touch_second.push(walked.touch_second[i]);
        blend.along.push(walked.along[i]);
    }
    push_exact(&mut blend, range.1, near(tail))?;
    let band = fit_open_band(&blend, radius, [false, false], tol)?;
    Ok(Marched {
        blend,
        band,
        fit_target: band_fit_target(tol),
    })
}

/// How far a marched seat's strip reaches along a face from the edge: the
/// longest walk on the face's surface from the edge to the line of contact
/// at any station, along the straight line between them in the chart.
fn marched_strip_width(
    marched: &Marched,
    side: &Beside,
    f: usize,
    edge: &Curve,
    tol: Tolerances,
) -> OgeomResult<f64> {
    const STEPS: usize = 16;
    let surface = &side.surface;
    let ((ua, ub), (va, vb)) = surface.domain();
    let blend = &marched.blend;
    let mut widest = 0.0_f64;
    for i in 0..blend.len() {
        let touch = if f == 0 {
            blend.on_first[i]
        } else {
            blend.on_second[i]
        };
        let crease =
            project_on_surface(surface, edge.point_at(blend.along[i], tol)?, 24, tol)?.parameters;
        let unwrap = |periodic: bool, span: f64, from: f64, to: f64| {
            if periodic && span > 0.0 {
                to - ((to - from) / span).round() * span
            } else {
                to
            }
        };
        let target = (
            unwrap(surface.is_periodic_u(), ub - ua, crease.0, touch.0),
            unwrap(surface.is_periodic_v(), vb - va, crease.1, touch.1),
        );
        let mut length = 0.0;
        let mut previous = surface.point_at(crease.0, crease.1, tol)?;
        for k in 1..=STEPS {
            #[allow(clippy::cast_precision_loss)]
            let w = (k as f64) / (STEPS as f64);
            let q = surface.point_at(
                crease.0 + (target.0 - crease.0) * w,
                crease.1 + (target.1 - crease.1) * w,
                tol,
            )?;
            length += q.distance(previous);
            previous = q;
        }
        widest = widest.max(length);
    }
    Ok(widest * 1.05)
}

/// The marched round's face: the band bounded by its borders along the
/// lines of contact and its rows at the ends.
fn marched_round(
    model: &mut Model,
    marched: &Marched,
    rails: &[Shape; 2],
    corners: &[[Shape; 2]; 2],
    away: bool,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let band: SurfaceGeometry = marched.band.clone().into();
    let (u_dom, v_dom) = band.domain();
    let band_id = model.geometry_mut().add_surface(band.clone());
    let widened = Tolerance::new(marched.fit_target)?;
    let mut arcs = Vec::with_capacity(2);
    for (k, (start, end)) in corners[0].iter().zip(&corners[1]).enumerate() {
        let row = marched.row(k, tol)?;
        let arc = make_edge_between(model, row, u_dom, start, end, tol)?.shape;
        model.widen(&arc, widened)?;
        arcs.push(arc);
    }
    let iso = |at: Point2, along: Direction2, lo: f64, hi: f64| -> OgeomResult<PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(Axis2::new(at, along), lo - 1.0, hi + 1.0)?.into())
    };
    for (f, u) in [(0, u_dom.0), (1, u_dom.1)] {
        attach_pcurve(
            model,
            &rails[f],
            iso(Point2::new(u, 0.0), Direction2::Y, v_dom.0, v_dom.1)?,
            band_id,
            Location::identity(),
            v_dom,
        )?;
    }
    for (k, v) in [(0, v_dom.0), (1, v_dom.1)] {
        attach_pcurve(
            model,
            &arcs[k],
            iso(Point2::new(0.0, v), Direction2::X, u_dom.0, u_dom.1)?,
            band_id,
            Location::identity(),
            u_dom,
        )?;
    }
    let wire = make_wire(
        model,
        &[
            arcs[0].clone(),
            rails[1].clone(),
            arcs[1].reversed(),
            rails[0].reversed(),
        ],
        tol,
    )?
    .shape;
    let face = make_face_on(model, band_id, std::slice::from_ref(&wire), tol)?.shape;
    let spine = marched.blend.spine.clone();
    oriented_round(model, face, away, |q| touching_ball(&spine, q), tol)
}

/// A closed-form seat between two faces: the ball, where it touches each
/// face, and the run of stations both faces reach.
pub(crate) struct FaceSeat {
    exact: Exact,
    stretches: Vec<Stretch>,
    /// The stations the round spans, or `None` for the whole turn.
    span: Option<(f64, f64)>,
    profiles: [Profile; 2],
    surfaces: Vec<SurfaceGeometry>,
    reads: Vec<(SurfaceId, Tolerance, Point, Vector)>,
}

/// The ball of `radius` touching both faces, on the side each face's
/// normal points to, or behind both with `behind`, where the two faces'
/// surfaces share a direction or an axis.
pub(crate) fn face_seat(
    model: &mut Model,
    faces: [&Shape; 2],
    radius: f64,
    behind: bool,
    tol: Tolerances,
) -> OgeomResult<FaceSeat> {
    let mut surfaces = Vec::with_capacity(2);
    let mut reads = Vec::with_capacity(2);
    for face in faces {
        let (id, surface, tolerance) = unplaced_face(model, face)?;
        let (p, n) = ogeom_algo::face_normal(model, face, tol)?;
        surfaces.push(surface);
        reads.push((id, tolerance, p, n));
    }
    let Some((layout, profiles)) = shared_layout(&surfaces, tol)? else {
        ogeom_bail!(
            Construction,
            "a round between two faces that share no edge is built where it has a closed \
             form: planes and cylinders along one direction, or planes, cylinders, cones, \
             spheres and tori about one axis; these two faces share neither"
        );
    };
    let side = if behind { "back" } else { "front" };
    let mut ball_sides = [0.0; 2];
    for f in 0..2 {
        let (_, _, p, n) = reads[f];
        let flip = if behind { -1.0 } else { 1.0 };
        ball_sides[f] = profile_sign(&profiles[f], &layout, p, n)? * flip;
    }
    let candidates = seatings(profiles, ball_sides, radius, tol)?;
    if candidates.is_empty() {
        ogeom_bail!(
            Construction,
            "no ball of radius {radius} touches the {side} of both faces' surfaces"
        );
    }

    // The stations each face spans, a little widened: its line of contact
    // is cut where it leaves the face, and stays inside the window of the
    // face's surface, which reaches a little past the face.
    let mut extents = [(0.0, 0.0); 2];
    for (extent, face) in extents.iter_mut().zip(faces) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for corner in shape_bounds(model, face, tol)?.corners() {
            let s = layout.station(corner);
            lo = lo.min(s);
            hi = hi.max(s);
        }
        let pad = ((hi - lo).abs() * 1e-3).max(tol.confusion() * 1e3);
        *extent = (lo - pad, hi + pad);
    }

    let mut seated = Vec::new();
    let (mut unmade, mut made) = (None, false);
    for section in candidates {
        let exact = match (Exact {
            layout,
            section,
            radius,
        })
        .checked(tol)
        {
            Ok(exact) => exact,
            Err(why) => {
                unmade = Some(why);
                continue;
            }
        };
        made = true;
        let mut found = Vec::with_capacity(2);
        for (f, face) in faces.iter().enumerate() {
            match contact_stretch(model, face, &exact, f, extents[f], tol)? {
                Some(stretch) => found.push(stretch),
                None => break,
            }
        }
        if found.len() == 2 {
            seated.push((exact, found));
        }
    }
    // Every ball the surfaces seat fails to make a round: that is why.
    if let (Some(why), false) = (unmade, made) {
        return Err(why);
    }
    let (exact, stretches) = match seated.len() {
        0 => ogeom_bail!(
            Construction,
            "a ball of radius {radius} touching the {side} of both faces' surfaces does not \
             touch both faces; it misses at least one of them"
        ),
        1 => seated.remove(0),
        _ => ogeom_bail!(
            Construction,
            "more than one ball of radius {radius} touches the {side} of both faces; which \
             corner to round is ambiguous"
        ),
    };

    // The run the round spans: where both lines of contact lie on their
    // faces.
    let closed = [stretches[0].closed, stretches[1].closed];
    let span = match closed {
        [true, true] => None,
        [false, false] => {
            let (a0, a1) = stretches[0].span;
            let (mut b0, mut b1) = stretches[1].span;
            if layout.turns() {
                let best = [-TAU, 0.0, TAU]
                    .into_iter()
                    .max_by(|x, y| {
                        let o = |k: f64| a1.min(b1 + k) - a0.max(b0 + k);
                        o(*x).total_cmp(&o(*y))
                    })
                    .unwrap_or(0.0);
                b0 += best;
                b1 += best;
            }
            let (lo, hi) = (a0.max(b0), a1.min(b1));
            if hi - lo <= tol.parametric() {
                ogeom_bail!(
                    Construction,
                    "the two faces reach the round over stretches of the corner that do not \
                     overlap"
                );
            }
            Some((lo, hi))
        }
        // One line closes on its face and the other runs across its own:
        // the round spans the open one's run.
        [true, false] => Some(stretches[1].span),
        [false, true] => Some(stretches[0].span),
    };
    Ok(FaceSeat {
        exact,
        stretches,
        span,
        profiles,
        surfaces,
        reads,
    })
}

/// The solid between a seat's round and the corner it rounds: the section
/// bounded by the two faces' surfaces, from where they cross to where the
/// ball touches each, and by the ball's arc, swept over the run. With it, a
/// point inside it beside the middle of the round, where the solid being
/// blended says whether the corner is material.
pub(crate) fn corner_wedge(
    model: &mut Model,
    seat: &FaceSeat,
    tol: Tolerances,
) -> OgeomResult<(Shape, Point)> {
    let exact = &seat.exact;
    let layout = exact.layout;
    let section = exact.section;
    let (s0, turn) = match seat.span {
        Some((lo, hi)) => (lo, hi - lo),
        None => (0.0, TAU),
    };
    let Some(crease) = crossings(seat.profiles[0], seat.profiles[1], tol)
        .into_iter()
        .min_by(|p, q| norm2(sub2(*p, section.centre)).total_cmp(&norm2(sub2(*q, section.centre))))
    else {
        ogeom_bail!(
            Construction,
            "the two faces' surfaces do not meet, so there is no corner between them for \
             the round to take off or fill"
        );
    };
    for contact in section.contacts {
        if norm2(sub2(contact, crease)) <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "the ball touches a face where the two surfaces cross; there is no corner \
                 between the round and the crease"
            );
        }
    }
    let apex = layout.lift(crease, s0);
    let ends = [exact.contact(0, s0), exact.contact(1, s0)];
    let apex_v = make_vertex(model, apex).shape;
    let end_v = [
        make_vertex(model, ends[0]).shape,
        make_vertex(model, ends[1]).shape,
    ];
    let normal = layout
        .lift_vector((1.0, 0.0), s0)
        .cross(layout.lift_vector((0.0, 1.0), s0));
    let normal = Direction::new(normal, tol)?;
    // Each leg runs along its face's profile between the crease and the
    // ball's touch: a segment, or the short arc of a circle.
    let mut legs = Vec::with_capacity(2);
    for f in 0..2 {
        let leg = match seat.profiles[f] {
            Profile::Line { .. } => {
                segment_between(model, (&apex_v, apex), (&end_v[f], ends[f]), tol)?
            }
            Profile::Circle { centre, radius } => {
                let centre = layout.lift(centre, s0);
                let (from, to) = (apex - centre, ends[f] - centre);
                let axis = from.cross(to);
                let sweep = axis.magnitude().atan2(from.dot(to));
                let frame = Frame::new(
                    centre,
                    Direction::new(axis, tol)?,
                    Direction::new(from, tol)?,
                    tol,
                )?;
                let circle = Curve::Circle(ogeom_geom::CircleCurve::new(Circle::new(
                    frame, radius, tol,
                )?));
                make_edge_between(model, circle, (0.0, sweep), &apex_v, &end_v[f], tol)?.shape
            }
        };
        legs.push(leg);
    }
    let (arc, range) = exact.arc(s0, tol)?;
    let arc = make_edge_between(model, arc, range, &end_v[0], &end_v[1], tol)?.shape;
    let reach = (apex.distance(ends[0]).max(apex.distance(ends[1])) + exact.radius) * 2.0;
    let plane = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(apex, normal),
        (-reach, reach),
        (-reach, reach),
    )?;
    let face = face_from_edges(
        model,
        plane.into(),
        &[legs[0].clone(), arc, legs[1].reversed()],
        tol,
    )?;
    let wedge = match layout {
        Layout::Extruded { along, .. } => {
            ogeom_algo::make_prism(model, &face, along * turn, tol)?.shape
        }
        Layout::Revolved { origin, axis, .. } => {
            let axis = Axis {
                location: origin,
                direction: Direction::new(axis, tol)?,
            };
            ogeom_algo::make_revolution(model, &face, axis, turn, tol)?.shape
        }
    };

    // Just off the middle of the ball's arc, toward the crease.
    let bisector = add2(
        sub2(section.contacts[0], section.centre),
        sub2(section.contacts[1], section.centre),
    );
    let on_arc = add2(
        section.centre,
        scale2(bisector, exact.radius / norm2(bisector)),
    );
    let gap = norm2(sub2(crease, on_arc));
    let toward = scale2(sub2(crease, on_arc), 1.0 / gap);
    let probe = add2(on_arc, scale2(toward, gap.min(exact.radius) * 0.1));
    Ok((wedge, layout.lift(probe, s0 + turn / 2.0)))
}

/// Round the corner between two faces of separate shapes where one is
/// curved and the seat has a closed form.
#[allow(clippy::too_many_lines, reason = "one rebuild, checked then assembled")]
pub(crate) fn fillet_faces(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    radius: f64,
    trim: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let FaceSeat {
        exact,
        stretches,
        span,
        profiles: _,
        surfaces,
        reads,
    } = face_seat(model, [a, b], radius, false, tol)?;
    let layout = exact.layout;

    let mut history = History::new();
    let away = false;
    if !trim {
        let round = match span {
            Some((lo, hi)) => {
                let mut corners = Vec::with_capacity(2);
                let mut rails = Vec::with_capacity(2);
                for f in 0..2 {
                    let v0 = make_vertex(model, exact.contact(f, lo)).shape;
                    let v1 = make_vertex(model, exact.contact(f, hi)).shape;
                    let (curve, range) = exact.contact_run(f, (lo, hi), tol)?;
                    let rail = make_edge_between(model, curve, range, &v0, &v1, tol)?.shape;
                    corners.push([v0, v1]);
                    rails.push(rail);
                }
                exact_open_round(
                    model,
                    &exact,
                    (lo, hi),
                    &[rails[0].clone(), rails[1].clone()],
                    &[corners[0].clone(), corners[1].clone()],
                    away,
                    tol,
                )?
            }
            None => {
                let start = layout.station(stretches[0].start);
                let mut rings = Vec::with_capacity(2);
                for f in 0..2 {
                    let curve = exact.contact_curve(f, tol)?;
                    let v = make_vertex(model, exact.contact(f, start)).shape;
                    rings.push(
                        make_edge_between(model, curve, (start, start + TAU), &v, &v, tol)?.shape,
                    );
                }
                let surface = exact.surface((0.0, TAU), tol)?;
                let band =
                    ogeom_algo::make_revolution_band(model, &surface, &rings[0], &rings[1], tol)?;
                oriented_round(model, band, away, |q| exact.centre(layout.station(q)), tol)?
            }
        };
        history.generate(a, round.clone());
        history.generate(b, round.clone());
        return Ok(Built::new(round, history));
    }

    // Each kept piece with its line of contact cut to the run, and the
    // round between them.
    let mut kept = Vec::with_capacity(2);
    let mut rails = Vec::with_capacity(2);
    let mut corners = Vec::with_capacity(2);
    for f in 0..2 {
        let stretch = &stretches[f];
        let (surface_id, tolerance, _, _) = reads[f];
        match span {
            None => {
                kept.push(stretch.piece.clone());
                rails.push(stretch.edge.clone());
                let Some((v, _)) = edge_vertices(model, &stretch.edge)? else {
                    ogeom_bail!(Construction, "a line of contact has no vertices");
                };
                corners.push([v.clone(), v]);
            }
            Some((lo, hi)) if stretch.closed => {
                let (piece, rail, ends) = closed_cut_to_run(
                    model,
                    stretch,
                    &exact,
                    f,
                    (lo, hi),
                    (surface_id, &surfaces[f], tolerance),
                    tol,
                )?;
                kept.push(piece);
                rails.push(rail);
                corners.push(ends);
            }
            Some((lo, hi)) => {
                let edge = stretch.edge.oriented(Orientation::Forward);
                let (edge_curve_of, edge_range) = edge_curve(model, &edge, tol)?;
                let Some((s, e)) = edge_vertices(model, &edge)? else {
                    ogeom_bail!(Construction, "a line of contact has no vertices");
                };
                let bounded = Curve::Trimmed(Box::new(TrimmedCurve::new(
                    edge_curve_of,
                    edge_range.0,
                    edge_range.1,
                    tol,
                )?));
                // The run's ends on this edge, in its own parameter.
                let mut ends = [0.0; 2];
                for (k, station) in [lo, hi].into_iter().enumerate() {
                    ends[k] =
                        project_on_curve(&bounded, exact.contact(f, station), 64, tol)?.parameter;
                }
                let rising = ends[1] > ends[0];
                let (t_lo, t_hi) = if rising {
                    (ends[0], ends[1])
                } else {
                    (ends[1], ends[0])
                };
                let mut cuts = vec![(edge_range.0, s.clone())];
                let at_lo = if (t_lo - edge_range.0).abs() <= tol.parametric() {
                    s.clone()
                } else {
                    let v =
                        make_vertex(model, exact.contact(f, if rising { lo } else { hi })).shape;
                    cuts.push((t_lo, v.clone()));
                    v
                };
                let at_hi = if (edge_range.1 - t_hi).abs() <= tol.parametric() {
                    e.clone()
                } else {
                    let v =
                        make_vertex(model, exact.contact(f, if rising { hi } else { lo })).shape;
                    cuts.push((t_hi, v.clone()));
                    v
                };
                cuts.push((edge_range.1, e.clone()));
                let (piece, rail) = if cuts.len() == 2 {
                    (stretch.piece.clone(), edge.clone())
                } else {
                    let mut pieces = Vec::with_capacity(3);
                    let mut rail = None;
                    for w in cuts.windows(2) {
                        let made = piece_of(
                            model,
                            &edge,
                            surface_id,
                            (w[0].0, w[1].0),
                            &w[0].1,
                            &w[1].1,
                            tol,
                        )?;
                        if (w[0].0 - t_lo).abs() <= tol.parametric() {
                            rail = Some(made.clone());
                        }
                        pieces.push(made);
                    }
                    let Some(rail) = rail else {
                        ogeom_bail!(Invariant, "the run's stretch of a line of contact was lost");
                    };
                    let face = rebuilt_face(
                        model,
                        &stretch.piece,
                        surface_id,
                        tolerance,
                        &[(edge.clone(), pieces)],
                        tol,
                    )?;
                    (face, rail)
                };
                kept.push(piece);
                // Each rail runs from the run's start to its end.
                if rising {
                    rails.push(rail);
                    corners.push([at_lo, at_hi]);
                } else {
                    rails.push(rail.reversed());
                    corners.push([at_hi, at_lo]);
                }
            }
        }
    }
    let rails = [rails[0].clone(), rails[1].clone()];
    let corners = [corners[0].clone(), corners[1].clone()];
    let round = match span {
        Some(run) => exact_open_round(model, &exact, run, &rails, &corners, away, tol)?,
        None => {
            let surface = exact.surface((0.0, TAU), tol)?;
            let band =
                ogeom_algo::make_revolution_band(model, &surface, &rails[0], &rails[1], tol)?;
            oriented_round(model, band, away, |q| exact.centre(layout.station(q)), tol)?
        }
    };
    let shell = model.add_shell(&[kept[0].clone(), round.clone(), kept[1].clone()])?;
    history.modify(a, kept[0].clone());
    history.modify(b, kept[1].clone());
    history.generate(a, round.clone());
    history.generate(b, round);
    Ok(Built::new(shell, history))
}

/// A kept piece whose line of contact closes on itself as a loop of its
/// own, cut to the run `lo..hi` the round spans: the loop restated as the
/// run (the rail, running from `lo` to `hi`) and the rest of the turn,
/// which stays free boundary. The rail and the run's end vertices come
/// back with the rebuilt piece.
fn closed_cut_to_run(
    model: &mut Model,
    stretch: &Stretch,
    exact: &Exact,
    f: usize,
    run: (f64, f64),
    host: (SurfaceId, &SurfaceGeometry, Tolerance),
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape, [Shape; 2])> {
    let (surface_id, surface, tolerance) = host;
    let edge = stretch.edge.oriented(Orientation::Forward);
    let alone = model
        .children_of(&stretch.piece.oriented(Orientation::Forward))?
        .into_iter()
        .any(|wire| {
            model
                .children_of(&wire.oriented(Orientation::Forward))
                .is_ok_and(|ring| ring.len() == 1 && ring[0].is_same(&edge))
        });
    if !alone {
        ogeom_bail!(
            Construction,
            "a face's line of contact closes on itself across the face's own boundary; the \
             round spanning part of its turn has nowhere to end on it"
        );
    }
    let (lo, hi) = run;
    let v_lo = make_vertex(model, exact.contact(f, lo)).shape;
    let v_hi = make_vertex(model, exact.contact(f, hi)).shape;
    let curve = exact.contact_curve_from(f, lo, tol)?;
    let rail = make_edge_between(model, curve.clone(), (0.0, hi - lo), &v_lo, &v_hi, tol)?.shape;
    let rest = make_edge_between(model, curve.clone(), (hi - lo, TAU), &v_hi, &v_lo, tol)?.shape;
    for made in [&rail, &rest] {
        attach_image(model, made, surface_id, surface, None, tol)?;
    }
    // The pieces run as the contact curve does; the old loop may not.
    let (old_curve, old_range) = edge_curve(model, &edge, tol)?;
    let at = old_curve.point_at(old_range.0, tol)?;
    let alike = old_curve
        .d1_at(old_range.0, tol)?
        .dot(curve.d1_at(exact.layout.station(at) - lo, tol)?)
        > 0.0;
    let pieces = if alike {
        vec![rail.clone(), rest]
    } else {
        vec![rest.reversed(), rail.reversed()]
    };
    let face = rebuilt_face(
        model,
        &stretch.piece,
        surface_id,
        tolerance,
        &[(edge, pieces)],
        tol,
    )?;
    if surface.is_periodic_u() || surface.is_periodic_v() {
        let wires = model.children_of(&face.oriented(Orientation::Forward))?;
        ogeom_algo::chain_wire_branches(model, surface_id, &wires, tol)?;
    }
    Ok((face, rail, [v_lo, v_hi]))
}

/// The layout both supports' profiles fit, if any: along a cylinder's
/// axis, or about an axis a revolved support (or a sphere's centre and a
/// plane's normal, or two spheres' centres) states.
fn shared_layout(
    surfaces: &[SurfaceGeometry],
    tol: Tolerances,
) -> OgeomResult<Option<(Layout, [Profile; 2])>> {
    let mut layouts = Vec::new();
    for surface in surfaces {
        if let SurfaceGeometry::Cylinder(c) = surface {
            let axis = c.cylinder().axis();
            layouts.push(Layout::extruded(
                axis.location,
                axis.direction.vector(),
                tol,
            )?);
        }
    }
    for surface in surfaces {
        let frame = match surface {
            SurfaceGeometry::Cylinder(c) => Some(c.cylinder().frame()),
            SurfaceGeometry::Cone(c) => Some(c.cone().frame()),
            SurfaceGeometry::Torus(t) => Some(t.torus().frame()),
            _ => None,
        };
        if let Some(frame) = frame {
            layouts.push(Layout::revolved(
                frame.origin(),
                frame.z().vector(),
                frame.x().vector(),
                tol,
            )?);
        }
    }
    let centres: Vec<(Point, Vector)> = surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Sphere(s) => {
                Some((s.sphere().centre(), s.sphere().frame().x().vector()))
            }
            _ => None,
        })
        .collect();
    for (centre, x) in &centres {
        for surface in surfaces {
            match surface {
                SurfaceGeometry::Plane(p) => {
                    let n = p.plane().frame().z().vector();
                    layouts.push(Layout::revolved(
                        *centre,
                        n,
                        any_square_to(n, tol)?.vector(),
                        tol,
                    )?);
                }
                SurfaceGeometry::Sphere(s)
                    if (s.sphere().centre() - *centre).magnitude() > tol.confusion() =>
                {
                    let axis = s.sphere().centre() - *centre;
                    let reference = if axis.cross(*x).magnitude() > tol.angular() {
                        *x
                    } else {
                        any_square_to(axis, tol)?.vector()
                    };
                    layouts.push(Layout::revolved(*centre, axis, reference, tol)?);
                }
                _ => {}
            }
        }
    }
    for layout in layouts {
        if let (Some(a), Some(b)) = (
            profile_of(&surfaces[0], &layout, tol),
            profile_of(&surfaces[1], &layout, tol),
        ) {
            return Ok(Some((layout, [a, b])));
        }
    }
    Ok(None)
}

/// Where a face's line of contact runs on it.
struct Stretch {
    /// The piece of the face on the far side of the line from the round.
    piece: Shape,
    /// The edge the line of contact cut it along.
    edge: Shape,
    /// Whether it closes on itself.
    closed: bool,
    /// The stations it spans, rising.
    span: (f64, f64),
    /// Where it starts.
    start: Point,
}

/// The face itself where none of its edges is placed; otherwise the face
/// again with each placed edge restated where it stands, on its placed
/// curve between the same (placed) vertices, with its image on the face.
fn with_edges_unplaced(model: &mut Model, face: &Shape, tol: Tolerances) -> OgeomResult<Shape> {
    let placed: Vec<Shape> = explore_unique(model, face, ShapeType::Edge)?
        .into_iter()
        .filter(|e| !e.location().is_identity())
        .map(|e| e.oriented(Orientation::Forward))
        .collect();
    if placed.is_empty() {
        return Ok(face.clone());
    }
    let (surface_id, _, tolerance) = unplaced_face(model, face)?;
    let mut swap = Vec::with_capacity(placed.len());
    for edge in placed {
        let (_, range) = edge_curve(model, &edge, tol)?;
        let Some((from, to)) = edge_vertices(model, &edge)? else {
            ogeom_bail!(Construction, "a boundary edge has no vertices");
        };
        let restated = piece_of(model, &edge, surface_id, range, &from, &to, tol)?;
        swap.push((edge, vec![restated]));
    }
    rebuilt_face(model, face, surface_id, tolerance, &swap, tol)
}

/// Cut a face along its line of contact over `extent` (or the whole turn)
/// and keep the piece away from the round; `None` where the line misses
/// the face.
fn contact_stretch(
    model: &mut Model,
    face: &Shape,
    exact: &Exact,
    f: usize,
    extent: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Stretch>> {
    let layout = exact.layout;
    let curve = exact.contact_curve(f, tol)?;
    let range = if layout.turns() { (0.0, TAU) } else { extent };
    let line = make_edge(model, curve, range, tol)?.shape;
    let face = &with_edges_unplaced(model, face, tol)?;
    let Ok(split) = ogeom_heal::split_face(
        model,
        face,
        face,
        std::slice::from_ref(&line),
        ogeom_heal::Projection::OnFace,
        tol,
    ) else {
        return Ok(None);
    };
    let cut: Vec<Shape> = split
        .history
        .generated(&line)
        .iter()
        .filter(|s| model.kind_of(s).is_ok_and(|k| k == ShapeType::Edge))
        .cloned()
        .collect();
    let edge = match cut.as_slice() {
        [] => return Ok(None),
        [one] => one.clone(),
        _ => ogeom_bail!(
            Construction,
            "the line the ball touches a face along crosses it more than once; the stretch the \
             round spans is ambiguous"
        ),
    };
    let (stretch_curve, stretch_range) = edge_curve(model, &edge, tol)?;
    let start = stretch_curve.point_at(stretch_range.0, tol)?;
    let end = stretch_curve.point_at(stretch_range.1, tol)?;
    let middle = stretch_curve.point_at(f64::midpoint(stretch_range.0, stretch_range.1), tol)?;
    let Some((first_vertex, last_vertex)) = edge_vertices(model, &edge)? else {
        ogeom_bail!(Construction, "a line of contact has no vertices");
    };
    // A full turn cut across a seam can come back as an open edge whose
    // ends meet: closed all the same.
    let closed = first_vertex.is_same(&last_vertex)
        || (layout.turns() && start.distance(end) <= tol.confusion() * 10.0);
    let (s0, sm, s1) = (
        layout.station(start),
        layout.station(middle),
        layout.station(end),
    );
    let span = if closed {
        (s0, s0 + TAU)
    } else if layout.turns() {
        let up = |from: f64, to: f64| from + (to - from).rem_euclid(TAU);
        let (m, e) = (up(s0, sm), up(s0, s1));
        if m < e { (s0, e) } else { (s1, up(s1, s0)) }
    } else {
        (s0.min(s1), s0.max(s1))
    };
    // Probed half way along, clear of any seam the stretch starts on.
    let probe_at = sm;

    // The piece away from the round: probed just off the line on the far
    // side from the round.
    let away = layout.lift_vector(exact.section.away()[f], probe_at);
    let on = exact.contact(f, probe_at);
    let (surface_id, surface, tolerance) = unplaced_face(model, face)?;
    let pieces = explore_unique(model, &split.shape, ShapeType::Face)?;
    let reach = if closed {
        exact.radius
    } else {
        (span.1 - span.0).abs().min(exact.radius)
    };
    for scale in [1e-3, 1e-2, 5e-2] {
        let eps = reach * scale;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        let probe = project_on_surface(&surface, on + away * eps, 32, tol)?.point;
        let mut holding = Vec::new();
        for piece in &pieces {
            if ogeom_algo::classify_on_face(model, piece, probe, deflection, tol)?
                == ogeom_algo::Containment::In
            {
                holding.push(piece.clone());
            }
        }
        let [piece] = holding.as_slice() else {
            continue;
        };
        if !explore_unique(model, piece, ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(&edge))
        {
            return Ok(None);
        }
        if closed && !first_vertex.is_same(&last_vertex) {
            let (piece, edge) = closed_up(
                model,
                piece,
                &edge,
                [&first_vertex, &last_vertex],
                exact,
                f,
                (surface_id, &surface, tolerance),
                tol,
            )?;
            return Ok(Some(Stretch {
                piece,
                edge,
                closed,
                span,
                start,
            }));
        }
        return Ok(Some(Stretch {
            piece: piece.clone(),
            edge,
            closed,
            span,
            start,
        }));
    }
    Ok(None)
}

/// A piece whose line of contact goes a full turn but ends on a second
/// vertex where it started, with whatever sliver of boundary the cut left
/// between the two: the line restated as one closed edge on one new vertex
/// where the ball touches, the sliver dropped, and every other edge at
/// either old vertex moved onto the new one.
#[allow(clippy::too_many_arguments, reason = "one repair, all its data")]
fn closed_up(
    model: &mut Model,
    piece: &Shape,
    edge: &Shape,
    ends: [&Shape; 2],
    exact: &Exact,
    f: usize,
    host: (SurfaceId, &SurfaceGeometry, Tolerance),
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let (surface_id, surface, tolerance) = host;
    let at = vertex_point(model, ends[0])?;
    if at.distance(vertex_point(model, ends[1])?) > tol.confusion() * 10.0 {
        ogeom_bail!(Invariant, "a closed line of contact's ends do not meet");
    }
    // Started on the layout's own column where the cut stands on it within
    // the confusion distance, so a band between this ring and another
    // started there finds one column for its seam.
    let s0 = if exact.contact(f, 0.0).distance(at) <= tol.confusion() {
        0.0
    } else {
        exact.layout.station(at)
    };
    let fresh = make_vertex(model, exact.contact(f, s0)).shape;
    let ring = make_edge_between(
        model,
        exact.contact_curve(f, tol)?,
        (s0, s0 + TAU),
        &fresh,
        &fresh,
        tol,
    )?
    .shape;
    let old_image = occurrence_image(model, &edge.oriented(Orientation::Forward), surface_id)?;
    let (old_curve, old_range) = edge_curve(model, edge, tol)?;
    let near = old_image.point_at(old_range.0, tol)?;
    attach_image(model, &ring, surface_id, surface, Some(near), tol)?;
    // The ring runs as the contact curve does; the old edge may not.
    let runs_alike = old_curve
        .d1_at(old_range.0, tol)?
        .dot(exact.contact_curve(f, tol)?.d1_at(s0, tol)?)
        > 0.0;
    let mut swap = vec![(
        edge.oriented(Orientation::Forward),
        vec![if runs_alike {
            ring.clone()
        } else {
            ring.reversed()
        }],
    )];
    let old_end = |v: &Shape| ends.iter().any(|e| e.is_same(v));
    for other in explore_unique(model, piece, ShapeType::Edge)? {
        let forward = other.oriented(Orientation::Forward);
        if forward.is_same(edge) {
            continue;
        }
        let Some((a, b)) = edge_vertices(model, &forward)? else {
            continue;
        };
        let (a_moves, b_moves) = (old_end(&a), old_end(&b));
        if a_moves && b_moves {
            // The sliver between the two ends.
            swap.push((forward, Vec::new()));
            continue;
        }
        if !a_moves && !b_moves {
            continue;
        }
        let (_, range) = edge_curve(model, &forward, tol)?;
        let from = if a_moves { fresh.clone() } else { a };
        let to = if b_moves { fresh.clone() } else { b };
        let moved = piece_of(model, &forward, surface_id, range, &from, &to, tol)?;
        swap.push((forward, vec![moved]));
    }
    let face = rebuilt_face(model, piece, surface_id, tolerance, &swap, tol)?;
    Ok((face, ring))
}

/// A vertex's point, placed.
fn vertex_point(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(point) = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex().map(|v| v.point))
    else {
        ogeom_bail!(Dangling, "vertex is not in this model");
    };
    Ok(vertex.transform(model.datums())?.apply(point))
}
