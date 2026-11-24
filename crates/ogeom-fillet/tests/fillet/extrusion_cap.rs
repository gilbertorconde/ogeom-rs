//! A spline cut preserves the exact cap of a rounded prism.
#![allow(clippy::unwrap_used, reason = "test constructions")]
use ogeom_core::Tolerances;
use ogeom_geom::{BSplineCurve, Curve, PlaneSurface, SurfaceGeometry, Transformable as _};
use ogeom_math::{Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom_topo::{Model, ShapeType, explore_unique};
const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_spline_cut_preserves_the_rounded_prisms_exact_cap() {
    let mut model = Model::new();
    let outline = [
        Point::new(5.43018, -35.461197, 0.0),
        Point::new(5.43018, 0.9563679999999977, 0.0),
        Point::new(-57.6881, -0.309752, 0.0),
        Point::new(-63.996323, -35.461197, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(&mut model, &outline, true, T)
        .unwrap()
        .shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    let face = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let target = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 24.5), T)
        .unwrap()
        .shape;
    let corner = Point::new(-57.6881, -0.309752, 24.5);
    let edges: Vec<_> = explore_unique(&model, &target, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|edge| {
            ogeom_algo::edge_vertices(&model, edge)
                .unwrap()
                .is_some_and(|(a, b)| {
                    [a, b].iter().any(|v| {
                        v.transform(model.datums())
                            .unwrap()
                            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
                            .distance(corner)
                            < T.confusion()
                    })
                })
        })
        .collect();
    assert_eq!(edges.len(), 3);
    let rounded = ogeom_fillet::fillet_edges(&mut model, &target, &edges, 5.0, T)
        .unwrap()
        .shape;
    let poles = [
        (-179.129269, -19.793100000000056),
        (-95.0316, -36.0137),
        (-27.3093783690672, -29.017896380949523),
        (-6.99477577057121, -11.647169254618543),
        (-2.9549141810512656, 1.6062356066273102),
        (12.823546881498215, 1.8499141341258323),
        (17.652645346179227, -13.96122200785365),
        (18.167615150405485, -29.091026102155464),
        (22.84704112905679, -69.93178584283447),
        (-169.29083468245173, -73.82738389906461),
        (-217.71691005887297, -43.78406246931849),
        (-212.97543059260772, -25.293625390522354),
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
    let result = ogeom_bool::cut(&mut model, &rounded, &tool, T).unwrap();
    assert!(
        ogeom_algo::check(&model, &result.shape, T)
            .unwrap()
            .is_valid()
    );
    let mut caps = 0;
    let mut worst = 0.0_f64;
    for face in explore_unique(&model, &result.shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        let surface = model
            .geometry()
            .surface(data.surface)
            .unwrap()
            .transformed(&face.transform(model.datums()).unwrap(), T)
            .unwrap();
        let SurfaceGeometry::Plane(plane) = surface else {
            continue;
        };
        let plane = plane.plane();
        if plane.normal().dot(Direction::Z).abs() < 1.0 - T.angular()
            || (plane.origin().z - 24.5).abs() > T.confusion()
        {
            continue;
        }
        caps += 1;
        for v in explore_unique(&model, &face, ShapeType::Vertex).unwrap() {
            let p = v
                .transform(model.datums())
                .unwrap()
                .apply(model.node(&v).unwrap().data().as_vertex().unwrap().point);
            worst = worst.max((p.z - 24.5).abs());
        }
    }
    assert!(caps > 0);
    // An intersection may already have an off-plane geometric endpoint
    // while its topological vertex still lies on the exact cap. A later
    // Boolean uses that curve to split edges, so check the geometry too.
    for edge in explore_unique(&model, &result.shape, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(ogeom_topo::EdgeRepr::Curve3d {
            curve,
            range,
            location,
        }) = data.curve3d()
        else {
            continue;
        };
        use ogeom_geom::Curve3d as _;
        let curve = model
            .geometry()
            .curve(*curve)
            .unwrap()
            .transformed(&location.composed(model.datums()).unwrap(), T)
            .unwrap()
            .transformed(&edge.transform(model.datums()).unwrap(), T)
            .unwrap();
        let ends = [
            curve.point_at(range.0, T).unwrap(),
            curve.point_at(range.1, T).unwrap(),
        ];
        if let Some((a, b)) = ogeom_algo::edge_vertices(&model, &edge).unwrap() {
            for v in [a, b] {
                let p = v
                    .transform(model.datums())
                    .unwrap()
                    .apply(model.node(&v).unwrap().data().as_vertex().unwrap().point);
                if (p.z - 24.5).abs() < T.confusion() {
                    let q = ends
                        .iter()
                        .min_by(|a, b| a.square_distance(p).total_cmp(&b.square_distance(p)))
                        .unwrap();
                    worst = worst.max((q.z - 24.5).abs());
                }
            }
        }
    }
    // The original cap is an exact plane. A large remote portion of the
    // cutting profile does not relax its model-confusion accuracy budget.
    assert!(worst < T.confusion(), "cap endpoint off z=24.5 by {worst}");
}
