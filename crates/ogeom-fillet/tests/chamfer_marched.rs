//! Chamfers on edges only the marching machinery seats: a bore's rim on a
//! drum.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// A drum of radius 10 and height 30 with a bore of radius 3 across it
/// along y at height 15.
fn holed_drum(m: &mut Model) -> Shape {
    let drum = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 0.0), 10.0, 30.0, T)
        .unwrap()
        .shape;
    let f = Frame::new(Point::new(0.0, -20.0, 15.0), Direction::Y, Direction::X, T).unwrap();
    let bore = ogeom_algo::make_cylinder(m, f, 3.0, 40.0, T).unwrap().shape;
    ogeom_bool::cut(m, &drum, &bore, T).unwrap().shape
}

/// The edge of `s` passing nearest `p`.
fn edge_near(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// How far `p` stands off the drum's wall (radius 10 about z) and off the
/// bore's (radius 3 about the line x = 0, z = 15).
fn off_drum(p: Point) -> f64 {
    (p.x.hypot(p.y) - 10.0).abs()
}

fn off_bore(p: Point) -> f64 {
    (p.x.hypot(p.z - 15.0) - 3.0).abs()
}

/// The nearest point of the crease to `p`, refined on the crease's own
/// parameter: where `p` stands in the plane square to the crease; and the
/// crease's unit tangent there.
fn crease_near(p: Point) -> (Point, Vector) {
    let at = |t: f64| {
        let x = 3.0 * t.cos();
        Point::new(x, (100.0 - x * x).sqrt(), 3.0f64.mul_add(t.sin(), 15.0))
    };
    let n = 4096;
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..n {
        let t = core::f64::consts::TAU * f64::from(i) / f64::from(n);
        let d = at(t).distance(p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let step = core::f64::consts::TAU / f64::from(n);
    let (mut lo, mut hi) = (best.1 - step, best.1 + step);
    for _ in 0..80 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if at(a).distance(p) < at(b).distance(p) {
            hi = b;
        } else {
            lo = a;
        }
    }
    let t = f64::midpoint(lo, hi);
    let h = 1e-6;
    let tangent = at(t + h) - at(t - h);
    (at(t), tangent / tangent.magnitude())
}

fn nearest_on_crease(p: Point) -> Point {
    crease_near(p).0
}

/// The holed drum's bore rim chamfered by `spec` (handed the drum's face
/// and the bore's), checked: back come its bevel's two contacts, each as
/// 32 points along each of the bevel faces' edges on that face, on the
/// drum and on the bore; and a bevel's surface. The boolean leaves the rim
/// in two pieces of one loop: with `whole`, both are chamfered, running on
/// into each other; else the one nearest a point on it, its bevel ending
/// where the rim runs on sharp.
fn bore_contacts(
    spec: impl Fn(&Shape, &Shape) -> Chamfer,
    whole: bool,
) -> (Vec<Point>, Vec<Point>, ogeom_geom::SurfaceGeometry) {
    let mut m = Model::new();
    let s = holed_drum(&mut m);
    let rim = edge_near(&m, &s, [0.0, 9.5, 18.0]);
    let beside: Vec<Shape> = explore_unique(&m, &s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(&m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(&rim))
        })
        .collect();
    let surface = |m: &Model, f: &Shape| {
        m.node(f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| m.geometry().surface(d.surface))
            .cloned()
            .unwrap()
    };
    let is_drum = |f: &Shape| match surface(&m, f) {
        ogeom_geom::SurfaceGeometry::Cylinder(c) => (c.cylinder().radius() - 10.0).abs() < 1e-9,
        _ => false,
    };
    let drum = beside
        .iter()
        .find(|f| is_drum(f))
        .expect("the drum")
        .clone();
    let bore = beside
        .iter()
        .find(|f| !is_drum(f))
        .expect("the bore")
        .clone();
    let rims: Vec<Shape> = if whole {
        let of = |f: &Shape| explore_unique(&m, f, ShapeType::Edge).unwrap();
        let on_bore = of(&bore);
        of(&drum)
            .into_iter()
            .filter(|e| {
                on_bore.iter().any(|x| x.is_same(e))
                    && ogeom_algo::shape_bounds(&m, e, T).unwrap().low().unwrap().y > 0.0
            })
            .collect()
    } else {
        vec![rim.clone()]
    };
    let specs: Vec<(Shape, Chamfer)> = rims
        .iter()
        .map(|e| (e.clone(), spec(&drum, &bore)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for e in &rims {
        assert!(built.history.is_deleted(e));
    }
    let bevels: Vec<Shape> = explore_unique(&m, &built.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| matches!(surface(&m, f), ogeom_geom::SurfaceGeometry::BSpline(_)))
        .collect();
    assert!(!bevels.is_empty(), "the bevel");
    if whole {
        assert_eq!(bevels.len(), rims.len(), "a bevel on each piece");
    }
    let bevel = bevels[0].clone();
    let (mut on_drum, mut on_bore) = (Vec::new(), Vec::new());
    for e in bevels
        .iter()
        .flat_map(|b| explore_unique(&m, b, ShapeType::Edge).unwrap())
    {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
            continue;
        };
        let points: Vec<Point> = (0..32)
            .map(|i| {
                c.point_at((b - a).mul_add(f64::from(i) / 32.0, a), T)
                    .unwrap()
            })
            .collect();
        let middle = c.point_at(0.5 * (a + b), T).unwrap();
        if off_drum(middle) <= FIT {
            on_drum.extend(points);
        } else if off_bore(middle) <= FIT {
            on_bore.extend(points);
        }
    }
    assert!(
        !on_drum.is_empty() && !on_bore.is_empty(),
        "a contact on each face"
    );
    (on_drum, on_bore, surface(&m, &bevel))
}

/// The band's fit target: how far a rail may stand off its face.
const FIT: f64 = 2e-4;

/// A leg of `length` measured along a section of curvature at most `kappa`
/// spans a chord between `length (1 - (kappa length)² / 24)` and `length`.
fn leg_fits(chord: f64, length: f64, kappa: f64) {
    let shortest = length * (1.0 - (kappa * length).powi(2) / 24.0);
    assert!(
        chord <= length + FIT && chord >= shortest - FIT,
        "a chord of {chord} for a leg of {length}"
    );
}

/// How sharply a cylinder of radius `r`'s section by the plane square to
/// the crease bends within `length` of the crease. The crease's tangent
/// lies in the face, so at the crease that plane holds the face's normal: a
/// normal section, of curvature at most `1 / r`. Along the leg the normal
/// turns by at most `length / r`, which leans the plane off it by as much,
/// and by Meusnier the section's curvature is the normal one over the
/// lean's cosine.
fn kappa(r: f64, length: f64) -> f64 {
    1.0 / (r * (length / r).cos())
}

#[test]
fn a_bore_rim_on_a_drum_chamfers_symmetric() {
    let (on_drum, on_bore, _) = bore_contacts(|_, _| Chamfer::Symmetric(0.5), false);
    for p in on_drum {
        leg_fits(p.distance(nearest_on_crease(p)), 0.5, kappa(10.0, 0.5));
    }
    for p in on_bore {
        leg_fits(p.distance(nearest_on_crease(p)), 0.5, kappa(3.0, 0.5));
    }
}

#[test]
fn a_bore_rim_on_a_drum_chamfers_by_two_distances() {
    let (on_drum, on_bore, _) = bore_contacts(
        |drum, _| Chamfer::Distances {
            face: drum.clone(),
            on_face: 0.8,
            on_other: 0.4,
        },
        false,
    );
    for p in on_drum {
        leg_fits(p.distance(nearest_on_crease(p)), 0.8, kappa(10.0, 0.8));
    }
    for p in on_bore {
        leg_fits(p.distance(nearest_on_crease(p)), 0.4, kappa(3.0, 0.4));
    }
}

#[test]
fn a_bore_rim_on_a_drum_chamfers_by_distance_and_angle() {
    let angle = 40.0_f64.to_radians();
    let (on_drum, _, bevel) = bore_contacts(
        |drum, _| Chamfer::Angle {
            face: drum.clone(),
            distance: 0.6,
            angle,
        },
        false,
    );
    for on_drum in on_drum {
        leg_fits(
            on_drum.distance(nearest_on_crease(on_drum)),
            0.6,
            kappa(10.0, 0.6),
        );
        // The bevel's line through the drum's contact in the plane square
        // to the crease there leaves the drum's section at the angle: its
        // sine against the drum's normal within that plane.
        let (_, along) = crease_near(on_drum);
        let foot = ogeom_geom::project::project_on_surface(&bevel, on_drum, 16, T).unwrap();
        assert!(foot.distance <= FIT, "the contact off the bevel");
        let n = bevel
            .normal_at(foot.parameters.0, foot.parameters.1, T)
            .unwrap()
            .vector();
        let ruling = n.cross(along);
        let normal = Vector::new(on_drum.x, on_drum.y, 0.0) / 10.0;
        let normal = normal - along * normal.dot(along);
        let normal = normal / normal.magnitude();
        let off = (ruling.dot(normal) / ruling.magnitude()).abs().asin();
        assert!((off - angle).abs() <= 1e-3, "{} degrees", off.to_degrees());
    }
}

/// Both pieces of the bore's rim chamfered 0.5: they run on into each other
/// round the whole loop, a bevel on each piece, meeting where the fitted
/// seam turns at its join. Each contact on its face, each leg its length
/// along its face's section.
#[test]
fn a_bore_rim_on_a_drum_chamfers_all_round() {
    let (on_drum, on_bore, _) = bore_contacts(|_, _| Chamfer::Symmetric(0.5), true);
    for p in on_drum {
        leg_fits(p.distance(nearest_on_crease(p)), 0.5, kappa(10.0, 0.5));
    }
    for p in on_bore {
        leg_fits(p.distance(nearest_on_crease(p)), 0.5, kappa(3.0, 0.5));
    }
}
