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
    assert!(include_str!("../src/lib.rs").contains("pub use engine::vertex::round_vertex;"));
    for (source, name) in [
        (include_str!("../src/setback.rs"), "setback_own"),
        (include_str!("../src/setback.rs"), "cornered_own"),
        (include_str!("../src/setback.rs"), "cornered_by_engine"),
        (include_str!("../src/sharp.rs"), "sharp_own"),
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
