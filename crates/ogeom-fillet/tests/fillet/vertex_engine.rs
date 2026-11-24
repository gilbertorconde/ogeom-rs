//! Vertex-only geometry, ancestry and engine ownership.
#![allow(
    clippy::unwrap_used,
    reason = "test constructions and source inspection"
)]

use ogeom_core::Tolerances;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A vertex-only round sheds its radius-sized corner block less one ball
/// octant. Rounding the three incident edges would shed long strips too.
#[test]
fn an_isolated_round_owns_its_sphere_and_three_planar_caps() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let vertex = explore_unique(&model, &solid, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| model.node(v).unwrap().data().as_vertex().unwrap().point == Point::ORIGIN)
        .unwrap();
    let radius: f64 = 2.0;
    let built = ogeom_fillet::round_vertex(&mut model, &solid, &vertex, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let generated = built.history.generated(&vertex);
    assert_eq!(generated.len(), 4, "one sphere and the three flush caps");
    let mut caps = [false; 3];
    let mut spheres = 0;
    for face in generated {
        let data = model.node(face).unwrap().data().as_face().unwrap();
        match model.geometry().surface(data.surface).unwrap() {
            SurfaceGeometry::Sphere(sphere) => {
                spheres += 1;
                assert!(
                    sphere
                        .sphere()
                        .frame()
                        .origin()
                        .distance(Point::new(radius, radius, radius))
                        < T.confusion()
                );
                assert!((sphere.sphere().radius() - radius).abs() < T.confusion());
            }
            SurfaceGeometry::Plane(_) => {
                let (at, normal) = ogeom_algo::face_normal(&model, face, T).unwrap();
                for (k, component) in [normal.x, normal.y, normal.z].iter().enumerate() {
                    if *component < -1.0 + T.angular() {
                        assert!(([at.x, at.y, at.z][k] - radius).abs() < T.confusion());
                        caps[k] = true;
                    }
                }
            }
            other => panic!("a corner generated an unexpected surface: {other:?}"),
        }
    }
    assert_eq!(spheres, 1);
    assert_eq!(caps, [true; 3]);
    let volume = ogeom_algo::volume_properties(
        &model,
        &built.shape,
        ogeom_mesh::Deflection::with_chord(1e-4).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    let expected = 8000.0 - radius.powi(3) * (1.0 - core::f64::consts::PI / 6.0);
    // The exact sphere and planes admit geometric integration; 1e-6 mm³
    // allows floating-point quadrature error at this scale.
    assert!(
        (volume - expected).abs() < 1e-6,
        "{volume} against {expected}"
    );
}

/// The cavity intersects the corner's cell, but no cavity face is wholly
/// inside it: every face extends past at least one setback plane.
#[test]
fn a_foreign_shell_crossing_the_setback_planes_is_refused_by_name() {
    for sharp in [false, true] {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 20.0), T)
            .unwrap()
            .shape;
        let cavity = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(1.0, 1.0, 1.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (6.0, 6.0, 6.0),
            T,
        )
        .unwrap()
        .shape;
        let solid = ogeom_bool::cut(&mut model, &solid, &cavity, T)
            .unwrap()
            .shape;
        let vertex = explore_unique(&model, &solid, ShapeType::Vertex)
            .unwrap()
            .into_iter()
            .find(|v| model.node(v).unwrap().data().as_vertex().unwrap().point == Point::ORIGIN)
            .unwrap();
        let setbacks = ogeom_fillet::Setbacks {
            per_edge: Vec::new(),
            default: 4.0,
            continuity: ogeom_geom::Continuity::G1,
        };
        let got = if sharp {
            let edge = explore_unique(&model, &solid, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .find(|e| {
                    let Some((a, b)) = ogeom_algo::edge_vertices(&model, e).unwrap() else {
                        return false;
                    };
                    let p = model.node(&a).unwrap().data().as_vertex().unwrap().point;
                    let q = model.node(&b).unwrap().data().as_vertex().unwrap().point;
                    (a.is_same(&vertex) || b.is_same(&vertex)) && (p.z - q.z).abs() > 1.0
                })
                .unwrap();
            ogeom_fillet::round_vertex_setback_sharp(
                &mut model, &solid, &vertex, &edge, 0.5, &setbacks, T,
            )
        } else {
            ogeom_fillet::round_vertex_setback(&mut model, &solid, &vertex, 0.5, &setbacks, T)
        };
        assert!(
            matches!(got, Err(ogeom_core::OgeomError::Construction(ref why)) if why.to_string().contains("another feature of the solid lies inside the setback corner's cell")),
            "sharp {sharp}: {got:?}"
        );
    }
}

fn function<'a>(source: &'a str, name: &str) -> &'a str {
    let start = source.find(&format!("fn {name}(")).unwrap();
    let body = &source[start..];
    let open = body.find('{').unwrap();
    let mut depth = 0;
    for (offset, c) in body[open..].char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    return &body[..=open + offset];
                }
            }
            _ => {}
        }
    }
    panic!("unclosed function {name}")
}

/// The unlogged vertex entries must not acquire an old constructor behind
/// an engine attempt. Geometry tests alone cannot distinguish that fallback.
#[test]
fn vertex_entries_do_not_reach_legacy_constructors() {
    assert!(include_str!("../../src/lib.rs").contains("pub use engine::vertex::round_vertex;"));
    for (source, name) in [
        (include_str!("../../src/setback.rs"), "setback_own"),
        (include_str!("../../src/setback.rs"), "cornered_own"),
        (include_str!("../../src/setback.rs"), "cornered_by_engine"),
        (include_str!("../../src/sharp.rs"), "sharp_own"),
    ] {
        let body = function(source, name);
        for legacy in [
            "fillet::fillet_edges(",
            "close_corner(",
            "close_corner_ends(",
            "close_sharp_corner(",
            "round_vertex_with(",
        ] {
            assert!(!body.contains(legacy), "{name} reaches {legacy}");
        }
    }
    let engine = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/engine");
    for file in std::fs::read_dir(engine).unwrap() {
        let file = file.unwrap();
        let name = file.file_name();
        let name = name.to_str().unwrap();
        if !name.starts_with("vertex") && name != "point_corner.rs" {
            continue;
        }
        let source = std::fs::read_to_string(file.path()).unwrap();
        for legacy in [
            "crate::corner::",
            "crate::corner_curved::",
            "crate::fillet::",
            "crate::sharp::",
            "close_corner_ends(",
            "close_sharp_corner(",
        ] {
            assert!(!source.contains(legacy), "{name} reaches {legacy}");
        }
    }
}

fn volume(model: &Model, shape: &ogeom_topo::Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn vertex_at(model: &Model, shape: &ogeom_topo::Shape, at: Point) -> ogeom_topo::Shape {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
            p.distance(at) < 1e-9
        })
        .unwrap()
}

/// Each corner of a box rounded alone generates its own sphere, centred a
/// radius in from every face at that corner, and its three planar caps:
/// the vertex's ancestry names its patch, not the whole solid. Composed
/// with all twelve bands, each vertex still generates exactly one face of
/// the result: composing histories keeps each corner's own patch.
#[test]
fn generated_corner_faces_belong_to_their_own_vertices() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let vertices = explore_unique(&model, &solid, ShapeType::Vertex).unwrap();
    for vertex in &vertices {
        let single = ogeom_fillet::round_vertex(&mut model, &solid, vertex, 1.0, T).unwrap();
        let point = model
            .node(vertex)
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point;
        let inward = |x: f64| if x < 5.0 { 1.0 } else { 9.0 };
        let centre = Point::new(inward(point.x), inward(point.y), inward(point.z));
        assert_eq!(single.history.generated(vertex).len(), 4);
        let mut spheres = 0;
        for patch in single.history.generated(vertex) {
            assert_eq!(model.kind_of(patch).unwrap(), ShapeType::Face);
            let surface = model.node(patch).unwrap().data().as_face().unwrap().surface;
            match model.geometry().surface(surface) {
                Some(SurfaceGeometry::Sphere(sphere)) => {
                    spheres += 1;
                    assert!(sphere.sphere().centre().distance(centre) < T.confusion());
                    assert!((sphere.sphere().radius() - 1.0).abs() < T.confusion());
                }
                Some(SurfaceGeometry::Plane(_)) => {}
                _ => panic!("a box corner has only its sphere and caps"),
            }
        }
        assert_eq!(spheres, 1);
    }
    let edges = explore_unique(&model, &solid, ShapeType::Edge).unwrap();
    let fillet = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 1.0, T).unwrap();
    let faces = explore_unique(&model, &fillet.shape, ShapeType::Face).unwrap();
    for vertex in &vertices {
        let patches = fillet
            .history
            .generated(vertex)
            .iter()
            .filter(|image| faces.iter().any(|f| f.is_same(image)))
            .count();
        assert_eq!(patches, 1, "composed history lost a corner patch");
    }
}

/// A box's eight corners rounded one after another, each corner's block
/// and ball clear of the others, take off exactly eight corner blocks less
/// their ball octants, `r^3 (1 - pi / 6)` each, at every radius up to balls
/// a fifth of a millimetre apart.
#[test]
fn a_box_s_corners_rounded_in_turn_take_off_their_corner_blocks_exactly() {
    for radius in [1.0_f64, 2.0, 4.9] {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let vertices = explore_unique(&model, &solid, ShapeType::Vertex).unwrap();
        let mut rounded = solid.clone();
        for vertex in &vertices {
            rounded = ogeom_fillet::round_vertex(&mut model, &rounded, vertex, radius, T)
                .unwrap()
                .shape;
        }
        let diagnosis = ogeom_algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "r {radius}: {:?}", diagnosis.problems);
        let want = 1000.0 - 8.0 * radius.powi(3) * (1.0 - core::f64::consts::PI / 6.0);
        let got = volume(&model, &rounded);
        assert!(
            (got - want).abs() <= want * 1e-9,
            "r {radius}: {got} against {want}"
        );
    }
}

/// The rounded corners of every labelling agree to within a part in a
/// billion of the solid.
fn labellings_agree(rounded: &[f64]) {
    for (index, v) in rounded.iter().enumerate() {
        assert!(
            (v - rounded[0]).abs() < 1e-9 * rounded[0],
            "labelling {index} builds {v}, labelling 0 {}",
            rounded[0]
        );
    }
}

/// An oblique corner of a sheared block, the block built from its three
/// sides in every order: the same corner, labelled six ways, rounds to the
/// same solid.
#[test]
fn an_oblique_corner_rounds_the_same_under_every_labelling() {
    use ogeom_math::Vector;
    let sides = [
        Vector::new(20.0, 0.0, 0.0),
        Vector::new(6.0, 20.0, 0.0),
        Vector::new(3.6, 6.0, 20.0),
    ];
    let orders = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut rounded = Vec::new();
    for order in orders {
        let mut model = Model::new();
        let block =
            ogeom_algo::make_parallelepiped(&mut model, Point::ORIGIN, order.map(|k| sides[k]), T)
                .unwrap()
                .shape;
        let vertex = vertex_at(&model, &block, Point::ORIGIN);
        let built = ogeom_fillet::round_vertex(&mut model, &block, &vertex, 2.0, T)
            .unwrap_or_else(|e| panic!("order {order:?}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &built.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "order {order:?}: {:?}",
            diagnosis.problems
        );
        rounded.push(volume(&model, &built.shape));
    }
    labellings_agree(&rounded);
}

/// A square pyramid's apex, its faces listed starting from each slope and in
/// both orders, each ring starting at another corner: the same apex,
/// labelled eight ways, rounds to the same solid.
#[test]
fn a_square_pyramid_apex_rounds_the_same_under_every_labelling() {
    let apex = Point::new(0.0, 0.0, 15.0);
    let points = [
        Point::new(-10.0, -10.0, 0.0),
        Point::new(10.0, -10.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(-10.0, 10.0, 0.0),
        apex,
    ];
    let mut rounded = Vec::new();
    for start in 0..4 {
        for backwards in [false, true] {
            let mut rings: Vec<Vec<usize>> = (0..4)
                .map(|step| {
                    let k = (start + step) % 4;
                    let ring = [k, (k + 1) % 4, 4];
                    // Each ring starting at another of its corners.
                    (0..3).map(|i| ring[(i + start) % 3]).collect()
                })
                .collect();
            rings.push((0..4).map(|i| [0, 3, 2, 1][(i + start) % 4]).collect());
            if backwards {
                rings.reverse();
            }
            let mut model = Model::new();
            let pyramid = ogeom_algo::make_polyhedron(&mut model, &points, &rings, T)
                .unwrap()
                .shape;
            let vertex = vertex_at(&model, &pyramid, apex);
            let built = ogeom_fillet::round_vertex(&mut model, &pyramid, &vertex, 2.0, T)
                .unwrap_or_else(|e| panic!("start {start}, backwards {backwards}: {e}"));
            let diagnosis = ogeom_algo::check(&model, &built.shape, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "start {start}, backwards {backwards}: {:?}",
                diagnosis.problems
            );
            rounded.push(volume(&model, &built.shape));
        }
    }
    labellings_agree(&rounded);
}
