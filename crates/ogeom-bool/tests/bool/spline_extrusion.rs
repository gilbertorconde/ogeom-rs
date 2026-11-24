//! Coincident contacts of a spline cut leave a complete floor.
#![allow(clippy::unwrap_used, reason = "test constructions")]
use ogeom_core::Tolerances;
use ogeom_geom::{BSplineCurve, Curve, PlaneSurface};
use ogeom_math::{Frame, KnotVector, Plane, Point, Vector};
use ogeom_topo::Model;
const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_coplanar_spline_cut_keeps_its_floor() {
    let mut model = Model::new();
    let poles = [
        (-179.129269, -19.793100000000056),
        (-89.41785586436728, -35.29058632234713),
        (-27.470414982562755, -28.96334532492927),
        (-7.1318, -11.6753),
        (-2.931494952547587, 1.6093846456950238),
        (12.822659840391081, 1.8499786355516328),
        (17.652333784763506, -13.961291221885196),
        (18.16767826291577, -29.091017034056588),
        (22.847037645776954, -69.9317858803962),
        (-169.29083535933546, -73.82738406508365),
        (-217.71690988535858, -43.78406244359132),
        (-212.97543058987898, -25.29362538754298),
        (-179.129269, -19.793100000000056),
    ];
    let control: Vec<_> = poles
        .iter()
        .rev()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let curve = BSplineCurve::new(
        KnotVector::clamped_uniform(3, control.len()).unwrap(),
        control,
        T,
    )
    .unwrap();
    let edge = ogeom_algo::make_edge(&mut model, Curve::BSpline(curve), (0.0, 1.0), T)
        .unwrap()
        .shape;
    let wire = ogeom_algo::make_wire(&mut model, &[edge], T).unwrap().shape;
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let tool = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 50.0), T)
        .unwrap()
        .shape;
    let corners = [
        (-59.53236991800965, -4.4894737953398645),
        (-59.274113, -32.003345),
        (5.430181999999988, -35.096714),
        (5.430181999999988, 0.9563679999999977),
    ];
    let points: Vec<_> = corners
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let wire = ogeom_algo::make_polygon(&mut model, &points, true, T)
        .unwrap()
        .shape;
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let target = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 24.5), T)
        .unwrap()
        .shape;
    let result = ogeom_bool::cut(&mut model, &target, &tool, T).unwrap();
    let check = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    let mass =
        ogeom_algo::volume_properties(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    // Independent profile sampling and convex clipping at 16k/64k/256k
    // stations extrapolates to this volume. The 1e-5 mm³ bound exceeds its
    // last extrapolation residual while resolving a missing floor.
    assert!((mass - 28840.851225099).abs() < 1e-5, "{mass}");
    assert!(
        ogeom_mesh::triangulate(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .is_closed()
    );
}
