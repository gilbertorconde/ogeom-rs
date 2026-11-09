//! A pad of a closed periodic spline profile: its top edge chamfered and
//! rounded, against the swept section.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with, fillet_edges};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const TAU: f64 = core::f64::consts::TAU;

/// The profile: a closed periodic cubic through eight points round the
/// origin in `z = 0`.
fn profile() -> ogeom_geom::Curve {
    let mut points: Vec<Point> = [
        (6.0, 0.0),
        (4.0, 4.0),
        (0.0, 5.0),
        (-4.0, 3.0),
        (-6.0, 0.0),
        (-3.0, -4.0),
        (0.0, -5.0),
        (4.0, -3.0),
    ]
    .iter()
    .map(|&(x, y)| Point::new(x, y, 0.0))
    .collect();
    points.push(points[0]);
    ogeom_geom::fit::fit_points_closed(&points, 3, 1e-6, T)
        .unwrap()
        .curve
        .into()
}

/// The pad 2 high on the profile: the solid, its top edge, its top face
/// and its wall.
fn pad(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let curve = profile();
    let domain = curve.domain();
    let edge = ogeom_algo::make_edge(m, curve, domain, T).unwrap().shape;
    let wire = ogeom_algo::make_wire(m, &[edge], T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let solid = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, 2.0), T)
        .unwrap()
        .shape;
    let rim = explore_unique(m, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom_algo::shape_bounds(m, e, T).unwrap();
            (b.low().unwrap().z - 2.0).abs() < 1e-6
                && (b.high().unwrap().z - 2.0).abs() < 1e-6
                && b.high().unwrap().x - b.low().unwrap().x > 1.0
        })
        .expect("the top edge");
    let beside: Vec<Shape> = explore_unique(m, &solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(&rim))
        })
        .collect();
    let is_plane = |f: &Shape| {
        matches!(
            m.node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| m.geometry().surface(d.surface)),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        )
    };
    let top = beside
        .iter()
        .find(|f| is_plane(f))
        .expect("the top")
        .clone();
    let wall = beside
        .iter()
        .find(|f| !is_plane(f))
        .expect("the wall")
        .clone();
    (solid, rim, top, wall)
}

/// The profile's length, by Gauss–Legendre over 512 pieces.
fn length() -> f64 {
    let curve = profile();
    let (lo, hi) = curve.domain();
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let mut length = 0.0;
    for k in 0..512 {
        let a = (hi - lo).mul_add(f64::from(k) / 512.0, lo);
        let b = (hi - lo).mul_add(f64::from(k + 1) / 512.0, lo);
        for (x, w) in nodes {
            let t = (b - a).mul_add(x, a + b) * 0.5;
            length += w * (b - a) * 0.5 * curve.d1_at(t, T).unwrap().magnitude();
        }
    }
    length
}

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// The volume a section of `area`, its centroid `inward` in from the
/// wall, takes off swept round the top edge: the edge closes once round,
/// its tangent turning through a whole turn, so each strip of the section
/// at depth `h` in sweeps the length less a turn's worth of `h`.
fn swept(area: f64, inward: f64) -> f64 {
    area * TAU.mul_add(-inward, length())
}

fn close(v: f64, expected: f64) {
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

fn after(build: impl Fn(&mut Model, &Shape, &Shape, &Shape, &Shape) -> Shape) -> (f64, f64) {
    let mut m = Model::new();
    let (solid, rim, top, wall) = pad(&mut m);
    let before = vol(&m, &solid);
    let built = build(&mut m, &solid, &rim, &top, &wall);
    let diagnosis = ogeom_algo::check(&m, &built, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    (before, vol(&m, &built))
}

#[test]
fn a_periodic_pad_chamfers_symmetric() {
    let (before, v) = after(|m, s, rim, _, _| {
        chamfer_edges_with(m, s, &[(rim.clone(), Chamfer::Symmetric(0.5))], T)
            .unwrap()
            .shape
    });
    close(v, before - swept(0.5 * 0.5 * 0.5, 0.5 / 3.0));
}

#[test]
fn a_periodic_pad_chamfers_by_two_distances() {
    let (before, v) = after(|m, s, rim, top, _| {
        let spec = Chamfer::Distances {
            face: top.clone(),
            on_face: 0.8,
            on_other: 0.4,
        };
        chamfer_edges_with(m, s, &[(rim.clone(), spec)], T)
            .unwrap()
            .shape
    });
    close(v, before - swept(0.5 * 0.8 * 0.4, 0.8 / 3.0));
}

#[test]
fn a_periodic_pad_rounds_its_top_edge() {
    // The corner between the legs and the quarter circle: its area
    // r² (1 - π/4), its centroid r (5/6 - π/4) / (1 - π/4) in from the
    // wall.
    let r = 0.5;
    let quarter = 1.0 - core::f64::consts::FRAC_PI_4;
    let (before, v) = after(|m, s, rim, _, _| {
        fillet_edges(m, s, std::slice::from_ref(rim), r, T)
            .unwrap()
            .shape
    });
    let inward = r * (5.0 / 6.0 - core::f64::consts::FRAC_PI_4) / quarter;
    close(v, before - swept(r * r * quarter, inward));
}
