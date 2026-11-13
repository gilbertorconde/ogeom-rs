//! Exact contacts, carried sections and spherical trim in edge rounds.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn cylinder_contacts_leave_the_middle_wall_between_both_bands() {
    for (wall, pitch, height, radius) in [(3.58_f64, 5.5_f64, 14.0, 8.0), (4.85, 8.5, 12.0, 4.0)] {
        let mut model = Model::new();
        let mut lobes = Vec::new();
        for i in 0..3 {
            let frame = Frame::new(
                Point::new(f64::from(i) * pitch, 0.0, 0.0),
                Direction::Z,
                Direction::X,
                T,
            )
            .unwrap();
            lobes.push(
                ogeom_algo::make_cylinder(&mut model, frame, wall, height, T)
                    .unwrap()
                    .shape,
            );
        }
        let two = ogeom_bool::fuse(&mut model, &lobes[0], &lobes[1], T)
            .unwrap()
            .shape;
        let solid = ogeom_bool::fuse(&mut model, &two, &lobes[2], T)
            .unwrap()
            .shape;
        let y = -(wall * wall - pitch * pitch / 4.0).sqrt();
        let mut edges = Vec::new();
        for x in [pitch / 2.0, 3.0 * pitch / 2.0] {
            let edge = explore_unique(&model, &solid, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .find(|edge| {
                    ogeom_algo::edge_vertices(&model, edge)
                        .unwrap()
                        .is_some_and(|(a, b)| {
                            let p = |v: &ogeom_topo::Shape| {
                                v.transform(model.datums())
                                    .unwrap()
                                    .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
                            };
                            let (a, b) = (p(&a), p(&b));
                            (a.x - x).abs() < 1e-7
                                && (b.x - x).abs() < 1e-7
                                && (a.y - y).abs() < 1e-7
                                && (b.y - y).abs() < 1e-7
                                && (a.z - b.z).abs() > height / 2.0
                        })
                })
                .expect("the middle lobe has both concave vertical seams");
            edges.push(edge);
        }
        // A tangent ball centre in this cross-section is at distance wall+radius
        // from both cylinder axes. Its two contacts on the middle circle stay
        // on opposite sides of the bottom: they leave a positive angular gap.
        let ball_y = -((wall + radius).powi(2) - pitch * pitch / 4.0).sqrt();
        let gap = 2.0 * (pitch / (2.0 * (wall + radius))).asin();
        assert!(gap > 0.0);
        let built = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, radius, T).unwrap();
        let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
        assert!(check.is_valid(), "{:?}", check.problems);
        let found = ogeom_algo::volume_properties(
            &model,
            &built.shape,
            ogeom_mesh::Deflection::default(),
            T,
        )
        .unwrap()
        .mass;
        // Integrate the cross-section independently: three discs less their
        // two lenses, plus two additions bounded by tangent circle arcs.
        let primitive =
            |r: f64, x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).asin());
        let half = pitch / 2.0;
        let contact = wall * half / (wall + radius);
        let lens =
            2.0 * wall * wall * (half / wall).acos() - pitch * (wall * wall - half * half).sqrt();
        let addition = 2.0
            * (-primitive(wall, half) + primitive(wall, contact)
                - ball_y * (half - contact)
                - primitive(radius, 0.0)
                + primitive(radius, contact - half));
        let volume =
            (3.0 * core::f64::consts::PI * wall * wall - 2.0 * lens + 2.0 * addition) * height;
        // All supports and strips are analytic cylinders, so mass is exact
        // to numerical quadrature, here within 1e-6 mm³.
        assert!(
            (found - volume).abs() < 1e-6,
            "{found} vs analytic {volume}"
        );
        // The new cylinders' axes are the exact external tangency centres,
        // not offsets along the tangent planes at the original seams.
        for x in [pitch / 2.0, 3.0 * pitch / 2.0] {
            let expected = Point::new(x, ball_y, height / 2.0);
            let generated: Vec<_> =
                ogeom_topo::explore_unique(&model, &built.shape, ShapeType::Face)
                    .unwrap()
                    .into_iter()
                    .filter(|f| {
                        let data = model.node(f).unwrap().data().as_face().unwrap();
                        let Some(ogeom_geom::SurfaceGeometry::Cylinder(c)) =
                            model.geometry().surface(data.surface)
                        else {
                            return false;
                        };
                        let axis = c.cylinder().axis();
                        let transform = f.transform(model.datums()).unwrap();
                        let origin = transform.apply(axis.location);
                        let direction = transform.apply_vector(axis.direction.vector());
                        (c.cylinder().radius() - radius).abs() < T.confusion()
                            && (expected - origin).cross(direction).magnitude() < T.confusion()
                    })
                    .collect();
            assert_eq!(
                generated.len(),
                1,
                "each seam has its own exact ball cylinder"
            );
        }
    }
}

/// A cylinder Ø20×26 and a 12.6×10×21 box. The shelf's exposed width at
/// y=2.5 is 12.6-sqrt(10²-2.5²)=2.91754, under the requested setback 3.
#[test]
#[ignore = "heavy"]
fn narrow_cylinder_shelf() {
    let mut model = Model::new();
    let cylinder = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 26.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(0.0, -5.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom_algo::make_box(&mut model, frame, (12.6, 10.0, 21.0), T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &cylinder, &block, T)
        .unwrap()
        .shape;
    let edges: Vec<_> = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        v.transform(model.datums())
                            .unwrap()
                            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
                    };
                    let (a, b) = (p(&a), p(&b));
                    ((a.x - 12.6).abs() < 1e-7
                        && (b.x - 12.6).abs() < 1e-7
                        && (a.z - 21.0).abs() < 1e-7
                        && (b.z - 21.0).abs() < 1e-7)
                        || ((a.z - 26.0).abs() < 1e-7 && (b.z - 26.0).abs() < 1e-7)
                })
        })
        .collect();
    assert!(edges.len() >= 2);
    let result = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 3.0, T).unwrap();
    let check = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    // Independent samples of the exact rolling envelope, including both
    // sides of each transition; no curve or surface from the builder.
    for y in [
        -5.0_f64, -2.80001, -2.8, -2.79999, -1.3, 0.0, 1.3, 2.79999, 2.8, 2.80001, 5.0,
    ] {
        let rho = 9.6_f64.hypot(y);
        let over = (10.0 - rho).max(0.0);
        let height = (9.0 - over * over).sqrt();
        let centre = ogeom_math::Point::new(9.6, y, 21.0 - height);
        let cap = if over > 0.0 {
            Point::new(96.0 / rho, 10.0 * y / rho, 21.0)
        } else {
            Point::new(9.6, y, 21.0)
        };
        let b = (cap - centre) / 3.0;
        let turn = b.x.acos();
        let across = (b - ogeom_math::Vector::X * b.x) / turn.sin();
        for k in 0..=16 {
            let t = turn * f64::from(k) / 16.0;
            let p = centre + (ogeom_math::Vector::X * t.cos() + across * t.sin()) * 3.0;
            // Ask each *trimmed* face directly. The generic shape-distance
            // scan can miss a spline's foot within its parameter epsilon
            // of a run's end, precisely where these transition probes sit.
            let mut d = f64::INFINITY;
            for f in explore_unique(&model, &result.shape, ShapeType::Face).unwrap() {
                use ogeom_geom::Transformable as _;
                let data = model.node(&f).unwrap().data().as_face().unwrap();
                let surface = model
                    .geometry()
                    .surface(data.surface)
                    .unwrap()
                    .clone()
                    .transformed(&f.transform(model.datums()).unwrap(), T)
                    .unwrap();
                let q = ogeom_algo::project_on_surface(&surface, p, 32, T).unwrap();
                if q.distance >= d {
                    continue;
                }
                let trim = ogeom_algo::classify::FaceTrim::of(
                    &model,
                    &f,
                    ogeom_mesh::Deflection {
                        chord: T.confusion(),
                        ..Default::default()
                    },
                    T,
                )
                .unwrap();
                if trim
                    .holds(&model, q.parameters, T.confusion(), T)
                    .unwrap()
                    .is_inside_or_on()
                {
                    d = d.min(q.distance);
                }
            }
            assert!(d < 2e-6, "exact envelope at {p:?} is {d} off the result");
        }
    }
    // A cavity between the supports intersects the new envelope but leaves
    // both contact rails intact. The swept-ball guard must still refuse it.
    let frame = Frame::new(Point::new(11.7, 0.2, 20.0), Direction::Z, Direction::X, T).unwrap();
    let hole = ogeom_algo::make_sphere(&mut model, frame, 0.3, T)
        .unwrap()
        .shape;
    let drilled = ogeom_bool::cut(&mut model, &solid, &hole, T).unwrap();
    let selected: Vec<_> = edges
        .iter()
        .map(|e| drilled.history.trace(e)[0].clone())
        .collect();
    let refusal =
        ogeom_fillet::fillet_edges(&mut model, &drilled.shape, &selected, 3.0, T).unwrap_err();
    assert!(
        refusal
            .to_string()
            .contains("enters a carried shelf's ball"),
        "{refusal}"
    );
}

/// No source spline: a quadrilateral prism, three incident edges, R=5.
/// Surface samples for an oracle must lie on the trimmed boundary as well
/// as on the supporting sphere. Every emitted node is measured, and the
/// resulting mesh must still be closed.
#[test]
fn oblique_corner_mesh() {
    use ogeom_math::{Plane, Vector};
    for degrees in [0.0_f64, 180.0, 270.0, 360.0, 450.0, 720.0, -90.0] {
        let (sin, cos) = degrees.to_radians().sin_cos();
        let rotate = |x: f64, y: f64, z: f64| Point::new(x * cos - y * sin, x * sin + y * cos, z);
        let mut model = Model::new();
        let outline = [
            (5.43018, -35.461197),
            (5.43018, 0.9563679999999977),
            (-57.6881, -0.309752),
            (-63.996323, -35.461197),
        ];
        let points: Vec<_> = outline.iter().map(|&(x, y)| rotate(x, y, 0.0)).collect();
        let wire = ogeom_algo::make_polygon(&mut model, &points, true, T)
            .unwrap()
            .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::Z),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let solid = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 24.5), T)
            .unwrap()
            .shape;
        let corner = rotate(-57.6881, -0.309752, 24.5);
        let edges: Vec<_> = explore_unique(&model, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                ogeom_algo::edge_vertices(&model, e)
                    .unwrap()
                    .is_some_and(|(a, b)| {
                        [&a, &b].into_iter().any(|v| {
                            v.transform(model.datums())
                                .unwrap()
                                .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
                                .distance(corner)
                                < 1e-7
                        })
                    })
            })
            .collect();
        assert_eq!(edges.len(), 3);
        let result = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 5.0, T).unwrap();
        let mesh = ogeom_mesh::triangulate(
            &model,
            &result.shape,
            ogeom_mesh::Deflection::with_chord(0.1).unwrap(),
            T,
        )
        .unwrap();
        assert!(
            mesh.is_closed(),
            "a corrected trim must retain a closed mesh"
        );
        let mut worst = (0.0, Point::ORIGIN);
        for p in &mesh.positions {
            let vertex = ogeom_algo::make_vertex(&mut model, *p).shape;
            let found = ogeom_algo::distance_between_shapes(
                &model,
                &vertex,
                &result.shape,
                ogeom_intersect::ExtremaOptions::default(),
                T,
            )
            .unwrap();
            if found.distance > worst.0 {
                worst = (found.distance, *p);
            }
        }
        assert!(
            worst.0 < T.confusion(),
            "at rotation {degrees}, an emitted node lies outside its own trim: {worst:?}"
        );
    }
}
