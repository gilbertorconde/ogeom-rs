//! Solids fused with their own mirror image, through STEP and IGES: the
//! image's circles and ellipses stand in frames the mirror left-handed, its
//! faces under a placement that mirrors, and its seams meet the original's
//! at poles and along shared rims. Each reads back as the solid it was,
//! split or unified.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn corner(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The part of `solid` inside the box at `keep`, fused with its mirror
/// image in the plane through the origin square to `normal`.
fn fused_with_its_image(
    model: &mut Model,
    solid: &Shape,
    keep: (Frame, (f64, f64, f64)),
    normal: Direction,
) -> Shape {
    let keep = ogeom::algo::make_box(model, keep.0, keep.1, T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(model, solid, &keep, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(model, &half).unwrap().shape;
    let mirror = Transform::plane_mirror(Point::ORIGIN, normal);
    let image = ogeom::algo::transformed(model, &copy, mirror)
        .unwrap()
        .shape;
    ogeom::boolean::fuse(model, &half, &image, T).unwrap().shape
}

/// A prism `h` high on the ellipse of radii `radii` about the z axis, its
/// major axis along y: the ellipse's frame has its `x` along y and its `y`
/// along -x, so the mirror in x = 0 runs its angle backwards.
fn elliptic_prism(model: &mut Model, radii: (f64, f64), h: f64) -> Shape {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::Y, T).unwrap();
    let ellipse = ogeom::math::Ellipse::new(frame, radii.0, radii.1, T).unwrap();
    let rim = ogeom::algo::make_edge(
        model,
        ogeom::geom::Curve::Ellipse(ogeom::geom::EllipseCurve::new(ellipse)),
        (0.0, 2.0 * core::f64::consts::PI),
        T,
    )
    .unwrap()
    .shape;
    let wire = ogeom::algo::make_wire(model, &[rim], T).unwrap().shape;
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD));
    let base = ogeom::algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom::algo::make_prism(model, &base, Vector::new(0.0, 0.0, h), T)
        .unwrap()
        .shape
}

/// The exchange formats a solid goes through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    Step,
    Iges,
}

/// `solid` written as `format` and read back: valid, drawing closed, and
/// measuring `closed_form`. The volume is an exact integral over the faces
/// read back, so it meets the closed form to rounding; 1e-9 relative leaves
/// room for the sums.
fn round_trips(format: Format, name: &str, model: &Model, solid: &Shape, closed_form: f64) {
    let mut document = ogeom::doc::Document::over(model.clone());
    document.add_part("part", solid.clone());
    let (document, solids) = match format {
        Format::Step => {
            let text = ogeom::io::write_step(&document, T).unwrap();
            let import = ogeom::io::read_step(&text, T).unwrap();
            (import.document, import.solids)
        }
        Format::Iges => {
            let text = ogeom::io::write_iges(&document, T).unwrap();
            let import = ogeom::io::read_iges(&text, T).unwrap();
            (import.document, import.solids)
        }
    };
    let name = format!("{name} through {format:?}");
    assert_eq!(solids.len(), 1, "{name}");
    let (back, solid) = (document.model(), &solids[0]);
    let diagnosis = ogeom::algo::check(back, solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let mesh = ogeom::mesh::triangulate(back, solid, Deflection::default(), T).unwrap();
    assert!(mesh.is_closed(), "{name}: the solid read back draws closed");
    let volume = ogeom::algo::volume_properties(back, solid, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (volume - closed_form).abs() < closed_form * 1e-9,
        "{name}: {volume} against {closed_form}"
    );
}

/// A ball about the x axis halved at its equator, a torus halved through
/// its middle plane, and a prism on an ellipse halved across its long
/// axis, each fused with its mirror image, then unified: every one, split
/// and unified, reads back from STEP and from IGES as the solid it was. The
/// unified prism's wall closes round the ellipse with its seam on the
/// mirror plane, and IGES carries that wall as a spline patch closed where
/// the ellipse starts, not periodic: no face on that patch can run across
/// its closing line, so it goes through STEP alone.
#[test]
fn solids_fused_with_their_mirror_image_round_trip() {
    use core::f64::consts::PI;
    let (ball, (major, minor), (a, b, h)) = (5.0, (10.0, 3.0), (4.0, 2.0, 5.0));
    let cases: [(&str, f64); 3] = [
        ("ball", 4.0 / 3.0 * PI * ball * ball * ball),
        ("torus", 2.0 * PI * PI * major * minor * minor),
        ("elliptic prism", PI * a * b * h),
    ];
    for (name, closed_form) in cases {
        let mut model = Model::new();
        let fused = match name {
            "ball" => {
                let polar = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
                let solid = ogeom::algo::make_sphere(&mut model, polar, ball, T)
                    .unwrap()
                    .shape;
                let keep = (corner(0.0, -6.0, -6.0), (6.0, 12.0, 12.0));
                fused_with_its_image(&mut model, &solid, keep, Direction::X)
            }
            "torus" => {
                let solid = ogeom::algo::make_torus(&mut model, Frame::WORLD, major, minor, T)
                    .unwrap()
                    .shape;
                let keep = (corner(-15.0, -15.0, 0.0), (30.0, 30.0, 5.0));
                fused_with_its_image(&mut model, &solid, keep, Direction::Z)
            }
            _ => {
                let solid = elliptic_prism(&mut model, (a, b), h);
                let keep = (corner(0.0, -6.0, -1.0), (6.0, 12.0, h + 2.0));
                fused_with_its_image(&mut model, &solid, keep, Direction::X)
            }
        };
        for format in [Format::Step, Format::Iges] {
            round_trips(
                format,
                &format!("{name}, split"),
                &model,
                &fused,
                closed_form,
            );
        }
        let unified = ogeom::heal::unify_same_domain(&mut model, &fused, T)
            .unwrap()
            .0
            .shape;
        let formats: &[Format] = if name == "elliptic prism" {
            &[Format::Step]
        } else {
            &[Format::Step, Format::Iges]
        };
        for &format in formats {
            round_trips(
                format,
                &format!("{name}, unified"),
                &model,
                &unified,
                closed_form,
            );
        }
    }
}

/// A ball of radius 5 bounded by one face whose seam is two edges running
/// away from a vertex on the equator, one to each pole: the file drops the
/// pole edges, and each seam's second use leaves a pole. It reads back
/// whole.
#[test]
fn a_seam_of_two_edges_running_from_one_vertex_round_trips() {
    let path = format!(
        "{}/../../tests/corpus/sphere_seam_two_ways.ogeom",
        env!("CARGO_MANIFEST_DIR")
    );
    let (model, roots) = ogeom::io::native::read(&std::fs::read_to_string(path).unwrap()).unwrap();
    let ball = 4.0 / 3.0 * core::f64::consts::PI * 125.0;
    for format in [Format::Step, Format::Iges] {
        round_trips(format, "ball", &model, &roots[0], ball);
    }
}

/// Balls about each axis, as built and mirrored: a ball's poles stand on
/// its axis, and its seam and pole edges with them, wherever the axis
/// points.
#[test]
fn balls_about_any_axis_round_trip() {
    let ball = 4.0 / 3.0 * core::f64::consts::PI * 125.0;
    let mirror = Transform::plane_mirror(Point::ORIGIN, Direction::X);
    for (axis, x) in [
        (Direction::Z, Direction::X),
        (Direction::X, Direction::Y),
        (Direction::Y, Direction::Z),
    ] {
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN, axis, x, T).unwrap();
        let solid = ogeom::algo::make_sphere(&mut model, frame, 5.0, T)
            .unwrap()
            .shape;
        let copy = ogeom::algo::copied(&mut model, &solid).unwrap().shape;
        let image = ogeom::algo::transformed(&mut model, &copy, mirror)
            .unwrap()
            .shape;
        for format in [Format::Step, Format::Iges] {
            round_trips(
                format,
                &format!("ball about {axis:?}"),
                &model,
                &solid,
                ball,
            );
            round_trips(
                format,
                &format!("mirrored ball about {axis:?}"),
                &model,
                &image,
                ball,
            );
        }
    }
}

/// An elliptic prism and its half on x > 0, each placed by a mirror: the
/// mirror turns every face's chart over, its planes written about their
/// normals and its wall as a spline patch, chart and all. Each reads back
/// with every face facing out of the material.
#[test]
fn a_mirrored_prism_of_planes_and_a_swept_wall_round_trips() {
    let (a, b, h) = (4.0, 2.0, 5.0);
    let prism = core::f64::consts::PI * a * b * h;
    let mut model = Model::new();
    let whole = elliptic_prism(&mut model, (a, b), h);
    let keep = ogeom::algo::make_box(&mut model, corner(0.0, -6.0, -1.0), (6.0, 12.0, h + 2.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &keep, T)
        .unwrap()
        .shape;
    let mirror = Transform::plane_mirror(Point::ORIGIN, Direction::X);
    for (name, solid, closed_form) in [("prism", whole, prism), ("half", half, prism / 2.0)] {
        let copy = ogeom::algo::copied(&mut model, &solid).unwrap().shape;
        let image = ogeom::algo::transformed(&mut model, &copy, mirror)
            .unwrap()
            .shape;
        for format in [Format::Step, Format::Iges] {
            round_trips(
                format,
                &format!("mirrored {name}"),
                &model,
                &image,
                closed_form,
            );
        }
    }
}

/// Each face of `solid` on a surface of its own, restated by `restate`
/// over the surface it had.
fn restated_surfaces(
    model: &mut Model,
    solid: &Shape,
    restate: impl Fn(usize, ogeom::geom::SurfaceGeometry) -> ogeom::geom::SurfaceGeometry,
) {
    use ogeom::topo::{NodeData, ShapeType};
    let faces = ogeom::topo::explore_unique(model, solid, ShapeType::Face).unwrap();
    for (k, face) in faces.iter().enumerate() {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
            continue;
        };
        let held = model.geometry().surface(data.surface).unwrap().clone();
        *model.geometry_mut().surface_mut(data.surface).unwrap() = restate(k, held);
    }
}

/// A box with every other face on a trimmed plane, and a drum whose wall
/// is stated as an offset of a narrower drum, each placed by a mirror.
/// STEP states a trimmed or offset surface over its basis, in the
/// right-handed frame on the basis's axis, so a mirrored basis's normal is
/// turned back there as a bare one's is: the faces keep their sense, and
/// the offset runs the other way along the normal it is read against.
#[test]
fn mirrored_faces_on_trimmed_and_offset_surfaces_round_trip_through_step() {
    use ogeom::geom::{OffsetSurface, Surface as _, SurfaceGeometry, TrimmedSurface};
    let mirror = Transform::plane_mirror(Point::ORIGIN, Direction::X);
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, corner(1.0, -2.0, 0.0), (3.0, 4.0, 5.0), T)
        .unwrap()
        .shape;
    restated_surfaces(&mut model, &block, |k, surface| {
        if k % 2 == 1 {
            return surface;
        }
        let (u, v) = surface.domain();
        SurfaceGeometry::Trimmed(Box::new(TrimmedSurface::new(surface, u, v, T).unwrap()))
    });
    let (r, h, d) = (3.0, 5.0, 0.5);
    let drum = ogeom::algo::make_cylinder(&mut model, corner(1.0, 0.0, 0.0), r, h, T)
        .unwrap()
        .shape;
    restated_surfaces(&mut model, &drum, |_, surface| match surface {
        SurfaceGeometry::Cylinder(c) => {
            let narrower = ogeom::math::Cylinder::new(c.cylinder().frame(), r - d, T).unwrap();
            let basis = ogeom::geom::CylinderSurface::new(narrower, c.domain().1).unwrap();
            SurfaceGeometry::Offset(Box::new(OffsetSurface::new(basis.into(), d).unwrap()))
        }
        other => other,
    });
    for (name, solid, closed_form) in [
        ("box", block, 3.0 * 4.0 * 5.0),
        ("drum", drum, core::f64::consts::PI * r * r * h),
    ] {
        let copy = ogeom::algo::copied(&mut model, &solid).unwrap().shape;
        let image = ogeom::algo::transformed(&mut model, &copy, mirror)
            .unwrap()
            .shape;
        round_trips(Format::Step, name, &model, &solid, closed_form);
        round_trips(
            Format::Step,
            &format!("mirrored {name}"),
            &model,
            &image,
            closed_form,
        );
    }
}

/// Every planar face of a box and of its mirror image, as written to STEP,
/// keeps its interior on the left of its outer bound about the face's
/// normal (its plane's axis, turned where the face's sense is false), as
/// ISO 10303-42 has a face's bounds: the bound's loop read through its
/// oriented edges' vertices, run the other way where the bound says so, and
/// its area vector (Newell's) taken against that normal.
#[test]
fn a_mirrored_solid_bounds_its_faces_interior_on_the_left_in_step() {
    use ogeom::core::FastMap;
    let mirror = Transform::plane_mirror(Point::ORIGIN, Direction::X);
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, corner(1.0, -2.0, 0.0), (3.0, 4.0, 5.0), T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(&mut model, &block).unwrap().shape;
    let image = ogeom::algo::transformed(&mut model, &copy, mirror)
        .unwrap()
        .shape;
    for (name, solid) in [("box", block), ("mirrored box", image)] {
        let mut document = ogeom::doc::Document::over(model.clone());
        document.add_part("part", solid);
        let text = ogeom::io::write_step(&document, T).unwrap();
        // Each entity's keyword and its arguments' references and flags.
        let mut entities: FastMap<u64, String> = FastMap::default();
        for statement in text.split(";\n") {
            let statement = statement.trim();
            if let Some(rest) = statement.strip_prefix('#')
                && let Some((id, body)) = rest.split_once('=')
            {
                entities.insert(id.trim().parse().unwrap(), body.to_string());
            }
        }
        let refs = |body: &str| -> Vec<u64> {
            body.split('#')
                .skip(1)
                .map(|s| {
                    s.chars()
                        .take_while(char::is_ascii_digit)
                        .collect::<String>()
                })
                .map(|s| s.parse().unwrap())
                .collect()
        };
        let reals = |body: &str| -> Vec<f64> {
            let inner = &body[body.rfind('(').unwrap() + 1..body.find(')').unwrap()];
            inner
                .split(',')
                .map(|x| x.trim().parse().unwrap())
                .collect()
        };
        let point = |id: u64| reals(&entities[&id]);
        let mut planar = 0;
        for body in entities.values().filter(|b| b.starts_with("ADVANCED_FACE")) {
            let face_refs = refs(body);
            let surface = &entities[face_refs.last().unwrap()];
            if !surface.starts_with("PLANE") {
                continue;
            }
            planar += 1;
            let axis = point(refs(&entities[&refs(surface)[0]])[1]);
            let sense = if body.trim_end().ends_with(".T.)") {
                1.0
            } else {
                -1.0
            };
            let outer = face_refs[..face_refs.len() - 1]
                .iter()
                .map(|b| &entities[b])
                .find(|b| b.starts_with("FACE_OUTER_BOUND"))
                .unwrap();
            let along = outer.trim_end().ends_with(".T.)");
            let mut corners: Vec<Vec<f64>> = Vec::new();
            for oriented in refs(&entities[&refs(outer)[0]]) {
                let oriented_body = &entities[&oriented];
                let forward = oriented_body.trim_end().ends_with(".T.)");
                let edge = refs(&entities[&refs(oriented_body)[0]]);
                let vertex = if forward { edge[0] } else { edge[1] };
                corners.push(point(refs(&entities[&vertex])[0]));
            }
            if !along {
                corners.reverse();
            }
            let mut area = [0.0f64; 3];
            for (i, p) in corners.iter().enumerate() {
                let q = &corners[(i + 1) % corners.len()];
                area[0] += (p[1] - q[1]) * (p[2] + q[2]);
                area[1] += (p[2] - q[2]) * (p[0] + q[0]);
                area[2] += (p[0] - q[0]) * (p[1] + q[1]);
            }
            let facing = sense * (area[0] * axis[0] + area[1] * axis[1] + area[2] * axis[2]);
            assert!(facing > 0.0, "{name}: {body} winds against its normal");
        }
        assert_eq!(planar, 6, "{name}");
    }
}

/// `text`, an IGES deck, with every loop (508) that fits on one parameter
/// line and holds edges alone, each with no parameter-space curve, started
/// `by` edges later: the same loop, walked from another edge. The entries
/// keep their text, so the line keeps its length and every pointer holds.
fn loops_started_later(text: &str, by: usize) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        let body = line.trim_end_matches(['\n', '\r']);
        let rotated = (body.len() == 80 && &body[72..73] == "P" && body.starts_with("508,"))
            .then(|| {
                let params = body[..64].trim_end();
                let inner = params.strip_prefix("508,")?.strip_suffix(';')?;
                let fields: Vec<&str> = inner.split(',').collect();
                let n: usize = fields.first()?.parse().ok()?;
                let entries: Vec<&[&str]> = fields[1..].chunks(5).collect();
                if entries.len() != n
                    || entries
                        .iter()
                        .any(|e| e.len() != 5 || e[0] != "0" || e[4] != "0")
                {
                    return None;
                }
                let mut turned = entries.clone();
                turned.rotate_left(by % n);
                let joined: Vec<String> = turned.iter().map(|e| e.join(",")).collect();
                Some(format!(
                    "{:<64}{}",
                    format!("508,{n},{};", joined.join(",")),
                    &body[64..]
                ))
            })
            .flatten();
        match rotated {
            Some(rotated) => {
                out.push_str(&rotated);
                out.push_str(&line[body.len()..]);
            }
            None => out.push_str(line),
        }
    }
    out
}

/// A prism on an ellipse, cut to its half on x > 0, goes out of IGES with
/// its wall on a spline patch closed where the ellipse starts, and one edge
/// of the wall lies along that closing line. Read with every loop started
/// at each of its edges in turn, as another system may write it, the half
/// reads back the same: a loop that opens on the closing line meets its
/// last edge across it.
#[test]
fn a_loop_opening_on_a_closing_line_reads_back_from_iges() {
    let (a, b, h) = (4.0, 2.0, 5.0);
    let mut model = Model::new();
    let whole = elliptic_prism(&mut model, (a, b), h);
    let keep = ogeom::algo::make_box(&mut model, corner(0.0, -6.0, -1.0), (6.0, 12.0, h + 2.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &keep, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("part", half);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    let closed_form = core::f64::consts::PI * a * b * h / 2.0;
    for by in 0..4 {
        let deck = loops_started_later(&text, by);
        assert!(
            by == 0 || deck != text,
            "the deck holds loops to start later"
        );
        let import = ogeom::io::read_iges(&deck, T).unwrap();
        let (back, solid) = (import.document.model(), &import.solids[0]);
        let diagnosis = ogeom::algo::check(back, solid, T).unwrap();
        assert!(diagnosis.is_valid(), "started {by} later: {diagnosis}");
        let volume = ogeom::algo::volume_properties(back, solid, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!(
            (volume - closed_form).abs() < closed_form * 1e-9,
            "started {by} later: {volume} against {closed_form}"
        );
    }
}

/// `text`, a STEP file, with every edge's stated pcurves dropped (an
/// `EDGE_CURVE` naming its `SURFACE_CURVE` or `SEAM_CURVE` names the curve
/// in space under it instead), so a reader derives each one, and with
/// every `EDGE_LOOP` started `by` edges later.
fn bare_loops_started_later(text: &str, by: usize) -> String {
    use ogeom::core::FastMap;
    let statements: Vec<&str> = text.split(";\n").collect();
    let mut under: FastMap<String, String> = FastMap::default();
    for statement in &statements {
        let statement = statement.trim();
        for keyword in ["SURFACE_CURVE('',", "SEAM_CURVE('',"] {
            if let Some((id, body)) = statement.split_once('=')
                && let Some(rest) = body.strip_prefix(keyword)
            {
                let curve = rest.split(',').next().unwrap().to_string();
                under.insert(id.trim().to_string(), curve);
            }
        }
    }
    let mut out = Vec::with_capacity(statements.len());
    for statement in statements {
        let mut statement = statement.to_string();
        if let Some((id, body)) = statement.clone().split_once('=') {
            if let Some(rest) = body.strip_prefix("EDGE_CURVE('',") {
                let mut parts: Vec<String> = rest.split(',').map(str::to_string).collect();
                if let Some(curve) = under.get(&parts[2]) {
                    parts[2].clone_from(curve);
                }
                statement = format!("{id}=EDGE_CURVE('',{}", parts.join(","));
            } else if let Some(rest) = body.strip_prefix("EDGE_LOOP('',(") {
                let list = rest.strip_suffix("))").unwrap();
                let mut uses: Vec<&str> = list.split(',').collect();
                let n = uses.len();
                uses.rotate_left(by % n);
                statement = format!("{id}=EDGE_LOOP('',({}))", uses.join(","));
            }
        }
        out.push(statement);
    }
    out.join(";\n")
}

/// The same half prism through STEP, its wall the closed spline patch it
/// reads back from IGES as: written with no pcurve, so the reader derives
/// every one, and every loop started at each of its edges in turn. It reads
/// back the same each time.
#[test]
fn a_loop_opening_on_a_closing_line_reads_back_from_step() {
    let (a, b, h) = (4.0, 2.0, 5.0);
    let mut model = Model::new();
    let whole = elliptic_prism(&mut model, (a, b), h);
    let keep = ogeom::algo::make_box(&mut model, corner(0.0, -6.0, -1.0), (6.0, 12.0, h + 2.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &keep, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("part", half);
    let iges = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
    let mut patched = iges.document;
    let solid = iges.solids[0].clone();
    patched.add_part("patched", solid);
    let text = ogeom::io::write_step(&patched, T).unwrap();
    let closed_form = core::f64::consts::PI * a * b * h / 2.0;
    for by in 0..4 {
        let file = bare_loops_started_later(&text, by);
        let import = ogeom::io::read_step(&file, T).unwrap();
        let (back, solid) = (import.document.model(), import.solids.last().unwrap());
        let diagnosis = ogeom::algo::check(back, solid, T).unwrap();
        assert!(diagnosis.is_valid(), "started {by} later: {diagnosis}");
        let volume = ogeom::algo::volume_properties(back, solid, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!(
            (volume - closed_form).abs() < closed_form * 1e-9,
            "started {by} later: {volume} against {closed_form}"
        );
    }
}
