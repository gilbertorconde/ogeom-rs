//! Fitted surfaces measured between their samples:
//! each result sampled densely and compared with the geometry it stands
//! for, independently of the points it was fitted to.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Circle, Frame, Point};
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The spline surfaces of every face of `shape`.
fn spline_surfaces(model: &Model, shape: &Shape) -> Vec<ogeom_geom::BSplineSurface> {
    let mut out = Vec::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
        let NodeData::Face(data) = model.node(&face).unwrap().data() else {
            continue;
        };
        if let Some(SurfaceGeometry::BSpline(patch)) = model.geometry().surface(data.surface) {
            out.push(patch.clone());
        }
    }
    out
}

/// The parameter `k / n` of the way across `range`.
fn across(range: (f64, f64), k: usize, n: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let f = k as f64 / n as f64;
    range.0 + (range.1 - range.0) * f
}

/// The distance from `p` to the circle of `radius` about the origin in the
/// XY plane.
fn to_flat_circle(p: Point, radius: f64) -> f64 {
    (p.x.hypot(p.y) - radius).hypot(p.z)
}

#[test]
fn a_filling_holds_its_curved_border_between_samples() {
    // A half disc's arc below the x axis, closed by three lines through a
    // height of one.
    let mut model = Model::new();
    let arc = ogeom_geom::CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap());
    let corners = [
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 1.0),
        Point::new(-1.0, 0.0, 1.0),
        Point::new(-1.0, 0.0, 0.0),
    ];
    let bottom = ogeom_algo::make_edge(&mut model, arc.into(), (PI, 2.0 * PI), T)
        .unwrap()
        .shape;
    let line = |model: &mut Model, a: Point, b: Point| {
        ogeom_algo::make_edge(
            model,
            ogeom_geom::LineCurve::segment(a, b, T).unwrap().into(),
            (0.0, a.distance(b)),
            T,
        )
        .unwrap()
        .shape
    };
    let right = line(&mut model, corners[0], corners[1]);
    let top = line(&mut model, corners[1], corners[2]);
    let left = line(&mut model, corners[2], corners[3]);
    let tolerance = 1e-4;
    let filled =
        ogeom_offset::make_filling(&mut model, &[bottom, right, top, left], 4, tolerance, T)
            .unwrap()
            .shape;
    let patches = spline_surfaces(&model, &filled);
    let patch = &patches[0];
    let (ud, vd) = patch.domain();
    let mut worst = 0.0_f64;
    for i in 0..=2000 {
        let p = patch.point_at(across(ud, i, 2000), vd.0, T).unwrap();
        worst = worst.max(to_flat_circle(p, 1.0));
    }
    eprintln!("filling border off its arc by {worst}");
    assert!(worst <= tolerance, "the border strays {worst} from its arc");
}
