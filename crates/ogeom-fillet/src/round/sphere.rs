//! P3's exact vertex blend for rounds: the sphere where one ball touches
//! all three faces of a trihedral vertex (ADR 0001, §4.4).
//!
//! Equal rounds on the three edges of a convex (or a concave) trihedral
//! vertex roll balls whose spines meet at one point, a radius from each
//! plane: the centre of the ball that touches all three. Each strip's end
//! section, in the plane square to its edge through that centre, is a
//! great circle of the ball, and its two ends are the ball's feet on the
//! strip's two faces. The vertex blend is the spherical triangle between
//! the three feet, bounded by those arcs.
//!
//! Its face is charted exactly where one of the planes is square to the
//! other two (a box's corner, a prism's corner on a square end): with the
//! sphere's pole at that plane's foot, the two arcs through it are
//! meridians, the third is the equator, and every pcurve is a straight line
//! of the chart, the pole a degenerate edge along its row. A trihedral
//! vertex with no such plane takes the general patch.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, SphereSurface, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Direction2, Frame, Point, Point2, Sphere, Vector, Vector2};
use ogeom_topo::{EdgeData, Location, Model, Shape};

/// The ball touching the three planes of a trihedral vertex, and how its
/// triangle is charted.
#[derive(Debug, Clone)]
pub(crate) struct SphereCorner {
    /// The ball's centre.
    pub centre: Point,
    /// Its radius.
    pub radius: f64,
    /// Its foot on each plane, in the order the planes were given.
    pub feet: [Point; 3],
    /// Whether the vertex is convex: the blend faces away from the centre.
    pub convex: bool,
    /// The plane whose foot is the chart's pole: square to the other two.
    pole: usize,
    /// The other two, the first at longitude zero and the second east of
    /// it.
    equator: [usize; 2],
}

/// One arc of the triangle: the curve, its range, and the feet it runs
/// between (by plane), from the first to the second.
pub(crate) type Arc = (Curve, (f64, f64), [usize; 2]);

impl SphereCorner {
    /// The ball of `radius` touching the three planes, each an outward unit
    /// normal and a point on it, inside a convex vertex's material or in a
    /// concave vertex's open corner; `None` where the planes meet in no
    /// single point, or no plane stands square to the other two, so the
    /// triangle has no exact chart.
    pub(crate) fn new(
        planes: [(Vector, Point); 3],
        radius: f64,
        convex: bool,
        tol: Tolerances,
    ) -> OgeomResult<Option<Self>> {
        if !radius.is_finite() || radius <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "a corner ball of radius {radius} rounds nothing"
            );
        }
        let n = planes.map(|(n, _)| n);
        // Inside a convex vertex the centre stands a radius behind each
        // plane, in a concave one a radius before it.
        let off = if convex { -radius } else { radius };
        let rhs = planes.map(|(n, p)| n.dot(p.to_vector()) + off);
        let det = n[0].dot(n[1].cross(n[2]));
        if det.abs() <= tol.angular() {
            return Ok(None);
        }
        let c = (n[1].cross(n[2]) * rhs[0] + n[2].cross(n[0]) * rhs[1] + n[0].cross(n[1]) * rhs[2])
            / det;
        let centre = Point::from_vector(c);
        let feet = [0, 1, 2].map(|i| centre - n[i] * off);
        let square = |i: usize, j: usize| n[i].dot(n[j]).abs() <= tol.angular();
        let Some(pole) = (0..3).find(|&i| (0..3).filter(|&j| j != i).all(|j| square(i, j))) else {
            return Ok(None);
        };
        let others: Vec<usize> = (0..3).filter(|&j| j != pole).collect();
        let z = feet[pole] - centre;
        let east = |a: usize, b: usize| z.cross(feet[a] - centre).dot(feet[b] - centre) > 0.0;
        let equator = if east(others[0], others[1]) {
            [others[0], others[1]]
        } else {
            [others[1], others[0]]
        };
        Ok(Some(Self {
            centre,
            radius,
            feet,
            convex,
            pole,
            equator,
        }))
    }

    /// The direction from the centre to plane `i`'s foot.
    fn toward(&self, i: usize, tol: Tolerances) -> OgeomResult<Direction> {
        Direction::new(self.feet[i] - self.centre, tol)
    }

    /// The sphere, charted with its pole at the pole plane's foot and its
    /// longitude zero at the first equator foot.
    pub(crate) fn surface(&self, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
        let frame = Frame::new(
            self.centre,
            self.toward(self.pole, tol)?,
            self.toward(self.equator[0], tol)?,
            tol,
        )?;
        Ok(SphereSurface::new(Sphere::new(frame, self.radius, tol)?).into())
    }

    /// The longitude of the second equator foot.
    fn east(&self) -> f64 {
        let (a, b) = (
            self.feet[self.equator[0]] - self.centre,
            self.feet[self.equator[1]] - self.centre,
        );
        a.cross(b).magnitude().atan2(a.dot(b))
    }

    /// The triangle's three arcs: along the equator between the equator
    /// feet, and down each meridian from the pole foot's plane to an
    /// equator foot (each run from the equator foot up to the pole).
    pub(crate) fn arcs(&self, tol: Tolerances) -> OgeomResult<[Arc; 3]> {
        let z = self.toward(self.pole, tol)?;
        let [a, b] = self.equator;
        let equator = Circle::new(
            Frame::new(self.centre, z, self.toward(a, tol)?, tol)?,
            self.radius,
            tol,
        )?;
        let meridian = |i: usize| -> OgeomResult<Curve> {
            let x = self.toward(i, tol)?;
            let w = Direction::new(x.vector().cross(z.vector()), tol)?;
            Ok(Curve::Circle(CircleCurve::new(Circle::new(
                Frame::new(self.centre, w, x, tol)?,
                self.radius,
                tol,
            )?)))
        };
        let quarter = core::f64::consts::FRAC_PI_2;
        Ok([
            (
                Curve::Circle(CircleCurve::new(equator)),
                (0.0, self.east()),
                [a, b],
            ),
            (meridian(a)?, (0.0, quarter), [a, self.pole]),
            (meridian(b)?, (0.0, quarter), [b, self.pole]),
        ])
    }

    /// The triangle's face on its sphere, bounded by `edges`, the arcs of
    /// [`SphereCorner::arcs`] in that order as edges between vertices at
    /// the feet (each running as its arc does), each given its line in the
    /// sphere's chart, and a degenerate edge along the pole's row; facing
    /// away from the centre at a convex vertex and toward it at a concave
    /// one.
    pub(crate) fn face(
        &self,
        model: &mut Model,
        edges: [&Shape; 3],
        pole: &Shape,
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let surface = self.surface(tol)?;
        let id = model.geometry_mut().add_surface(surface);
        let east = self.east();
        let quarter = core::f64::consts::FRAC_PI_2;
        let line = |from: (f64, f64),
                    along: (f64, f64),
                    to: f64|
         -> OgeomResult<ogeom_geom::PlanarCurve> {
            Ok(ogeom_geom::Line2d::over(
                ogeom_math::Axis2::new(
                    Point2::new(from.0, from.1),
                    Direction2::new(Vector2::new(along.0, along.1), tol)?,
                ),
                0.0,
                to,
            )?
            .into())
        };
        let [equator, west, east_meridian] = edges;
        for (edge, pcurve, range) in [
            (equator, line((0.0, 0.0), (1.0, 0.0), east)?, (0.0, east)),
            (west, line((0.0, 0.0), (0.0, 1.0), quarter)?, (0.0, quarter)),
            (
                east_meridian,
                line((east, 0.0), (0.0, 1.0), quarter)?,
                (0.0, quarter),
            ),
        ] {
            ogeom_algo::attach_pcurve(model, edge, pcurve, id, Location::identity(), range)?;
        }
        // The pole: a degenerate edge along the chart's top row, back from
        // the east meridian to the west.
        let collapsed = {
            let mut data = EdgeData::new();
            data.degenerate = true;
            model.add_edge(data, &[pole.clone(), pole.clone()])?
        };
        ogeom_algo::attach_pcurve(
            model,
            &collapsed,
            line((east, quarter), (-1.0, 0.0), east)?,
            id,
            Location::identity(),
            (0.0, east),
        )?;
        let wire = ogeom_algo::make_wire(
            model,
            &[
                equator.clone(),
                east_meridian.clone(),
                collapsed,
                west.reversed(),
            ],
            tol,
        )?
        .shape;
        let face = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape;
        Ok(if self.convex { face } else { face.reversed() })
    }
}

/// The ball of `radius` touching the three planes, each an outward unit
/// normal and a point on it, inside a convex vertex's material or in a
/// concave vertex's open corner: its centre and its foot on each plane, in
/// the order the planes were given; `None` where the planes meet in no
/// single point.
pub(crate) fn corner_ball(
    planes: [(Vector, Point); 3],
    radius: f64,
    convex: bool,
    tol: Tolerances,
) -> Option<(Point, [Point; 3])> {
    let n = planes.map(|(n, _)| n);
    let off = if convex { -radius } else { radius };
    let rhs = planes.map(|(n, p)| n.dot(p.to_vector()) + off);
    let det = n[0].dot(n[1].cross(n[2]));
    if det.abs() <= tol.angular() {
        return None;
    }
    let c =
        (n[1].cross(n[2]) * rhs[0] + n[2].cross(n[0]) * rhs[1] + n[0].cross(n[1]) * rhs[2]) / det;
    let centre = Point::from_vector(c);
    Some((centre, [0, 1, 2].map(|i| centre - n[i] * off)))
}

/// The ball about `centre` of `radius`, charted for the spherical triangle
/// between `feet`: its longitude half a turn from the triangle's middle and
/// its equator through it, so the triangle, less than a hemisphere, stands
/// clear of both poles and of the chart's seam, whatever the corner.
pub(crate) fn triangle_sphere(
    centre: Point,
    radius: f64,
    feet: [Point; 3],
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let middle = feet
        .iter()
        .fold(Vector::ZERO, |sum, f| sum + (*f - centre) / radius);
    let away = Direction::new(-middle, tol)?;
    let a = away.vector();
    let helper = if a.x.abs() <= a.y.abs() && a.x.abs() <= a.z.abs() {
        Vector::new(1.0, 0.0, 0.0)
    } else if a.y.abs() <= a.z.abs() {
        Vector::new(0.0, 1.0, 0.0)
    } else {
        Vector::new(0.0, 0.0, 1.0)
    };
    let up = Direction::new(a.cross(helper), tol)?;
    let frame = Frame::new(centre, up, away, tol)?;
    Ok(SphereSurface::new(Sphere::new(frame, radius, tol)?).into())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::SphereCorner;
    use ogeom_core::Tolerances;
    use ogeom_geom::Curve3d as _;
    use ogeom_math::{Point, Vector};

    const T: Tolerances = Tolerances::millimetres();

    /// The box corner at the origin's top: the top (+z at 10), the front
    /// (-y at 0) and the left (-x at 0).
    fn corner() -> [(Vector, Point); 3] {
        [
            (Vector::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, 10.0)),
            (Vector::new(0.0, -1.0, 0.0), Point::ORIGIN),
            (Vector::new(-1.0, 0.0, 0.0), Point::ORIGIN),
        ]
    }

    /// The ball sits a radius in from each face, each foot is the centre's
    /// projection on its face, each arc runs between its feet a radius off
    /// the centre, and the triangle's area is an eighth of the sphere's.
    #[test]
    fn the_corner_ball_touches_three_faces_and_its_triangle_is_an_octant() {
        let r = 2.0;
        let ball = SphereCorner::new(corner(), r, true, T).unwrap().unwrap();
        assert!(
            ball.centre.distance(Point::new(2.0, 2.0, 8.0)) < 1e-12,
            "{:?}",
            ball.centre
        );
        assert!(ball.feet[0].distance(Point::new(2.0, 2.0, 10.0)) < 1e-12);
        assert!(ball.feet[1].distance(Point::new(2.0, 0.0, 8.0)) < 1e-12);
        assert!(ball.feet[2].distance(Point::new(0.0, 2.0, 8.0)) < 1e-12);
        let arcs = ball.arcs(T).unwrap();
        for (curve, (a, b), [from, to]) in &arcs {
            assert!(curve.point_at(*a, T).unwrap().distance(ball.feet[*from]) < 1e-12);
            assert!(curve.point_at(*b, T).unwrap().distance(ball.feet[*to]) < 1e-12);
            let mid = curve.point_at(0.5 * (a + b), T).unwrap();
            assert!((mid.distance(ball.centre) - r).abs() < 1e-12);
        }
        let mut m = ogeom_topo::Model::new();
        let v = ball.feet.map(|p| ogeom_algo::make_vertex(&mut m, p).shape);
        let edge = |m: &mut ogeom_topo::Model, k: usize| {
            let (curve, range, [from, to]) = arcs[k].clone();
            ogeom_algo::make_edge_between(m, curve, range, &v[from], &v[to], T)
                .unwrap()
                .shape
        };
        let (e0, e1, e2) = (edge(&mut m, 0), edge(&mut m, 1), edge(&mut m, 2));
        let face = ball.face(&mut m, [&e0, &e1, &e2], &v[0], T).unwrap();
        let fine = ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        };
        let area = ogeom_algo::surface_properties(&m, &face, fine, T)
            .unwrap()
            .mass;
        let octant = core::f64::consts::PI * r * r / 2.0;
        assert!((area - octant).abs() < 1e-6, "{area} against {octant}");
    }

    /// In a concave corner the ball stands a radius out from each face,
    /// and where no plane is square to the other two there is no exact
    /// chart.
    #[test]
    fn a_concave_corner_s_ball_stands_out_and_an_oblique_corner_has_no_chart() {
        let inward = corner().map(|(n, p)| (-n, p));
        let ball = SphereCorner::new(inward, 1.5, false, T).unwrap().unwrap();
        assert!(
            ball.centre.distance(Point::new(1.5, 1.5, 8.5)) < 1e-12,
            "{:?}",
            ball.centre
        );
        let mut oblique = corner();
        let lean = Vector::new(-1.0, -0.3, 0.3);
        oblique[2].0 = lean / lean.magnitude();
        assert!(SphereCorner::new(oblique, 1.0, true, T).unwrap().is_none());
    }

    /// A prism's corner on its square end, its sides leaning apart: the
    /// end is square to both, its foot the pole, and the triangle between
    /// the pole and the two side feet covers `r^2 phi`, `phi` the angle
    /// between the sides' normals.
    #[test]
    fn a_prism_corner_on_its_square_end_is_charted_and_covers_r2_phi() {
        let mut leaning = corner();
        let lean = Vector::new(-1.0, -0.4, 0.0);
        leaning[2].0 = lean / lean.magnitude();
        let r = 1.5;
        let ball = SphereCorner::new(leaning, r, true, T).unwrap().unwrap();
        let arcs = ball.arcs(T).unwrap();
        let mut m = ogeom_topo::Model::new();
        let v = ball.feet.map(|p| ogeom_algo::make_vertex(&mut m, p).shape);
        let mut edges = Vec::new();
        for (curve, range, [from, to]) in arcs {
            edges.push(
                ogeom_algo::make_edge_between(&mut m, curve, range, &v[from], &v[to], T)
                    .unwrap()
                    .shape,
            );
        }
        let face = ball
            .face(&mut m, [&edges[0], &edges[1], &edges[2]], &v[0], T)
            .unwrap();
        let fine = ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        };
        let area = ogeom_algo::surface_properties(&m, &face, fine, T)
            .unwrap()
            .mass;
        let phi = leaning[1].0.dot(leaning[2].0).clamp(-1.0, 1.0).acos();
        assert!(
            (area - r * r * phi).abs() < 1e-6,
            "{area} against {}",
            r * r * phi
        );
    }
}
