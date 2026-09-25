//! STEP's other curves and surfaces: the swept and offset surfaces and the
//! open conics written exactly and read back as themselves, and the kinds
//! a writer may spell differently (trims, polylines, composites, the knot
//! forms, trimmed, offset and replicated surfaces) read as the geometry
//! they name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::collections::BTreeMap;

use ogeom::core::Tolerances;
use ogeom::geom::{Curve, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> Deflection {
    Deflection::with_chord(1e-3).unwrap()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

fn round_trip(model: Model, solid: &Shape) -> (f64, ogeom::io::StepImport) {
    let before = volume(&model, solid);
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("part", solid.clone());
    let text = ogeom::io::write_step(&document, T).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    (before, import)
}

fn surfaces(model: &Model, shape: &Shape) -> Vec<SurfaceGeometry> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .map(|f| {
            let data = model.node(&f).unwrap().data().as_face().unwrap().clone();
            model.geometry().surface(data.surface).unwrap().clone()
        })
        .collect()
}

fn curves(model: &Model, shape: &Shape) -> Vec<Curve> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter_map(|e| {
            let data = model.node(&e)?.data().as_edge()?.clone();
            let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return None;
            };
            model.geometry().curve(*curve).cloned()
        })
        .collect()
}

/// A profile with a spline side turned half round the z axis: the spline
/// sweeps a surface of revolution, written as one and read back as one.
#[test]
fn a_surface_of_revolution_round_trips_exactly() {
    let mut model = Model::new();
    let pts = [
        Point::new(2.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 5.0),
        Point::new(2.0, 0.0, 5.0),
    ];
    let v: Vec<Shape> = pts
        .iter()
        .map(|p| ogeom::algo::make_vertex(&mut model, *p).shape)
        .collect();
    let line = |model: &mut Model, a: usize, b: usize| {
        let c: Curve = ogeom::geom::LineCurve::segment(pts[a], pts[b], T)
            .unwrap()
            .into();
        let d = ogeom::geom::Curve3d::domain(&c);
        ogeom::algo::make_edge_between(model, c, d, &v[a], &v[b], T)
            .unwrap()
            .shape
    };
    let bottom = line(&mut model, 0, 1);
    let spline = {
        let fitted =
            ogeom::geom::fit::fit_points(&[pts[1], Point::new(4.5, 0.0, 2.5), pts[2]], 2, 1e-9, T)
                .unwrap();
        let c = Curve::BSpline(fitted.curve);
        let d = ogeom::geom::Curve3d::domain(&c);
        ogeom::algo::make_edge_between(&mut model, c, d, &v[1], &v[2], T)
            .unwrap()
            .shape
    };
    let top = line(&mut model, 2, 3);
    let inner = line(&mut model, 3, 0);
    let plane =
        ogeom::geom::PlaneSurface::new(ogeom::math::Plane::through(Point::ORIGIN, Direction::Y));
    let profile = ogeom::algo::make_face_with_pcurves(
        &mut model,
        plane.into(),
        &[vec![bottom, spline, top, inner]],
        T,
    )
    .unwrap()
    .shape;
    let axis = ogeom::math::Axis::new(Point::ORIGIN, Direction::Z);
    let solid = ogeom::algo::make_revolution(&mut model, &profile, axis, core::f64::consts::PI, T)
        .unwrap()
        .shape;
    assert!(
        surfaces(&model, &solid)
            .iter()
            .any(|s| matches!(s, SurfaceGeometry::Revolution(_)))
    );
    let (before, import) = round_trip(model, &solid);
    let back = import.document.model();
    let read = &import.solids[0];
    assert!(
        surfaces(back, read)
            .iter()
            .any(|s| matches!(s, SurfaceGeometry::Revolution(_))),
        "the spline's sweep reads back as a surface of revolution"
    );
    let after = volume(back, read);
    assert!(
        (after - before).abs() < 1e-6 * before,
        "{before} -> {after}"
    );
}

/// A prism over a parabolic segment: its curved wall is the parabola
/// swept straight, written as a linear extrusion of a parabola.
#[test]
fn a_parabola_and_its_extrusion_round_trip_exactly() {
    let mut model = Model::new();
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::Y, T).unwrap();
    let parabola = ogeom::math::Parabola::new(frame, 1.0, T).unwrap();
    // Opening along +y: the arc from x = -2 to x = 2, closed by y = 1.
    let arc: Curve = ogeom::geom::ParabolaCurve::new(parabola, 10.0)
        .unwrap()
        .into();
    let (a, b) = (Point::new(2.0, 1.0, 0.0), Point::new(-2.0, 1.0, 0.0));
    let va = ogeom::algo::make_vertex(&mut model, a).shape;
    let vb = ogeom::algo::make_vertex(&mut model, b).shape;
    let t_of = |p: Point| frame.to_local(p).y;
    let (ta, tb) = (t_of(a), t_of(b));
    let (lo, hi, from, to) = if ta < tb {
        (ta, tb, &va, &vb)
    } else {
        (tb, ta, &vb, &va)
    };
    let curved = ogeom::algo::make_edge_between(&mut model, arc, (lo, hi), from, to, T)
        .unwrap()
        .shape;
    let chord: Curve = ogeom::geom::LineCurve::segment(b, a, T).unwrap().into();
    let d = ogeom::geom::Curve3d::domain(&chord);
    let straight = ogeom::algo::make_edge_between(&mut model, chord, d, &vb, &va, T)
        .unwrap()
        .shape;
    let wire = ogeom::algo::make_wire_unordered(&mut model, &[curved, straight], T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &wire, ShapeType::Edge).unwrap();
    let plane =
        ogeom::geom::PlaneSurface::new(ogeom::math::Plane::through(Point::ORIGIN, Direction::Z));
    let profile = ogeom::algo::make_face_with_pcurves(&mut model, plane.into(), &[edges], T)
        .unwrap()
        .shape;
    let solid = ogeom::algo::make_prism(&mut model, &profile, Vector::Z * 3.0, T)
        .unwrap()
        .shape;
    let (before, import) = round_trip(model, &solid);
    let back = import.document.model();
    let read = &import.solids[0];
    assert!(
        curves(back, read)
            .iter()
            .any(|c| matches!(c, Curve::Parabola(_))),
        "the parabola reads back as one"
    );
    // The segment's area is two thirds of its bounding rectangle, 8/3,
    // over a height of three.
    assert!((before - 8.0).abs() < 1e-6, "{before}");
    let after = volume(back, read);
    assert!(
        (after - before).abs() < 1e-6 * before,
        "{before} -> {after}"
    );
}

/// A written cube's STEP text, entity by entity.
struct Deck {
    lines: BTreeMap<u64, String>,
    head: String,
    tail: String,
}

impl Deck {
    fn cube() -> Self {
        let mut model = Model::new();
        let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let mut document = ogeom::doc::Document::over(model);
        document.add_part("cube", cube);
        let text = ogeom::io::write_step(&document, T).unwrap();
        let (head, rest) = text.split_once("DATA;\n").unwrap();
        let (body, tail) = rest.split_once("ENDSEC;\nEND").unwrap();
        let mut lines = BTreeMap::new();
        for line in body.lines() {
            let (id, _) = line.split_once('=').unwrap();
            lines.insert(id[1..].parse::<u64>().unwrap(), line.to_owned());
        }
        Self {
            lines,
            head: format!("{head}DATA;\n"),
            tail: format!("ENDSEC;\nEND{tail}"),
        }
    }

    fn next_id(&self) -> u64 {
        self.lines.keys().last().copied().unwrap_or(0) + 1
    }

    fn add(&mut self, body: &str) -> u64 {
        let id = self.next_id();
        self.lines.insert(id, format!("#{id}={body};"));
        id
    }

    fn body(&self, id: u64) -> &str {
        let line = &self.lines[&id];
        &line[line.find('=').unwrap() + 1..]
    }

    fn refs(&self, id: u64) -> Vec<u64> {
        let body = self.body(id);
        let mut out = Vec::new();
        let mut rest = body;
        while let Some(at) = rest.find('#') {
            let digits: String = rest[at + 1..]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect();
            out.push(digits.parse().unwrap());
            rest = &rest[at + 1 + digits.len()..];
        }
        out
    }

    fn find(&self, keyword: &str) -> Vec<u64> {
        self.lines
            .iter()
            .filter(|(_, l)| l.contains(&format!("={keyword}(")))
            .map(|(id, _)| *id)
            .collect()
    }

    fn replace_in(&mut self, id: u64, from: &str, to: &str) {
        let line = self.lines[&id].replacen(from, to, 1);
        self.lines.insert(id, line);
    }

    fn text(&self) -> String {
        let mut out = self.head.clone();
        for line in self.lines.values() {
            out.push_str(line);
            out.push('\n');
        }
        out.push_str(&self.tail);
        out
    }

    fn read_volume(&self) -> f64 {
        let import = ogeom::io::read_step(&self.text(), T).unwrap();
        assert_eq!(import.solids.len(), 1, "{:?}", import.report.warnings);
        volume(import.document.model(), &import.solids[0])
    }

    /// The first straight edge: its EDGE_CURVE, its LINE, and its two
    /// vertices' CARTESIAN_POINTs.
    fn first_line_edge(&self) -> (u64, u64, u64, u64) {
        for edge in self.find("EDGE_CURVE") {
            let r = self.refs(edge);
            let (va, vb, curve) = (r[0], r[1], r[2]);
            if self.body(curve).starts_with("LINE(") {
                return (edge, curve, self.refs(va)[0], self.refs(vb)[0]);
            }
        }
        panic!("the cube has a straight edge");
    }

    /// The numbers in an entity's last parenthesised list.
    fn coords(&self, id: u64) -> Vec<f64> {
        let body = self.body(id);
        let inner = &body[body.rfind('(').unwrap() + 1..body.find(')').unwrap()];
        inner
            .split(',')
            .map(|x| x.trim().parse().unwrap())
            .collect()
    }

    /// The top face's ADVANCED_FACE and its PLANE: the plane through z = 10.
    fn top_face(&self) -> (u64, u64) {
        for face in self.find("ADVANCED_FACE") {
            let plane = *self.refs(face).last().unwrap();
            if !self.body(plane).starts_with("PLANE(") {
                continue;
            }
            let placement = self.refs(plane)[0];
            let point = self.refs(placement)[0];
            let z = self.coords(point)[2];
            if (z - 10.0).abs() < 1e-9 {
                return (face, plane);
            }
        }
        panic!("the cube has a top face at z = 10");
    }
}

#[test]
fn edges_spelled_as_trims_polylines_composites_and_knot_forms_read() {
    for spelling in ["trimmed", "polyline", "composite", "bezier"] {
        let mut deck = Deck::cube();
        let (edge, line, pa, pb) = deck.first_line_edge();
        let replacement = match spelling {
            "trimmed" => deck.add(&format!(
                "TRIMMED_CURVE('',#{line},(#{pa}),(#{pb}),.T.,.CARTESIAN.)"
            )),
            "polyline" => deck.add(&format!("POLYLINE('',(#{pa},#{pb}))")),
            "composite" => {
                let (ca, cb) = (deck.coords(pa), deck.coords(pb));
                let mid = deck.add(&format!(
                    "CARTESIAN_POINT('',({},{},{}))",
                    (ca[0] + cb[0]) / 2.0,
                    (ca[1] + cb[1]) / 2.0,
                    (ca[2] + cb[2]) / 2.0
                ));
                let first = deck.add(&format!(
                    "TRIMMED_CURVE('',#{line},(#{pa}),(#{mid}),.T.,.CARTESIAN.)"
                ));
                let second = deck.add(&format!(
                    "TRIMMED_CURVE('',#{line},(#{mid}),(#{pb}),.T.,.CARTESIAN.)"
                ));
                let s1 = deck.add(&format!(
                    "COMPOSITE_CURVE_SEGMENT(.CONTINUOUS.,.T.,#{first})"
                ));
                let s2 = deck.add(&format!(
                    "COMPOSITE_CURVE_SEGMENT(.CONTINUOUS.,.T.,#{second})"
                ));
                deck.add(&format!("COMPOSITE_CURVE('',(#{s1},#{s2}),.F.)"))
            }
            _ => deck.add(&format!(
                "BEZIER_CURVE('',1,(#{pa},#{pb}),.UNSPECIFIED.,.F.,.F.)"
            )),
        };
        deck.replace_in(edge, &format!("#{line},"), &format!("#{replacement},"));
        let got = deck.read_volume();
        assert!((got - 1000.0).abs() < 1e-6, "{spelling}: {got}");
    }
}

#[test]
fn faces_on_trimmed_offset_and_replicated_planes_read() {
    for spelling in ["trimmed", "offset", "replica", "composite"] {
        let mut deck = Deck::cube();
        let (face, plane) = deck.top_face();
        let replacement = match spelling {
            "trimmed" => deck.add(&format!(
                "RECTANGULAR_TRIMMED_SURFACE('',#{plane},-5.,15.,-5.,15.,.T.,.T.)"
            )),
            "offset" => {
                // The plane half a unit below, offset half a unit up.
                let p = deck.add("CARTESIAN_POINT('',(0.,0.,9.5))");
                let z = deck.add("DIRECTION('',(0.,0.,1.))");
                let x = deck.add("DIRECTION('',(1.,0.,0.))");
                let frame = deck.add(&format!("AXIS2_PLACEMENT_3D('',#{p},#{z},#{x})"));
                let lower = deck.add(&format!("PLANE('',#{frame})"));
                deck.add(&format!("OFFSET_SURFACE('',#{lower},0.5,.F.)"))
            }
            "composite" => {
                // Two patches across the top, meeting at x = 5: a bilinear
                // spline patch, and the plane itself trimmed to its window.
                let corners: Vec<u64> = [(-1.0, -1.0), (-1.0, 11.0), (5.0, -1.0), (5.0, 11.0)]
                    .iter()
                    .map(|(x, y)| deck.add(&format!("CARTESIAN_POINT('',({x:?},{y:?},10.))")))
                    .collect();
                let spline = deck.add(&format!(
                    "B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#{},#{}),(#{},#{})),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)",
                    corners[0], corners[1], corners[2], corners[3]
                ));
                let window = deck.add(&format!(
                    "RECTANGULAR_TRIMMED_SURFACE('',#{plane},5.,11.,-1.,11.,.T.,.T.)"
                ));
                let a = deck.add(&format!(
                    "SURFACE_PATCH('',#{spline},.DISCONTINUOUS.,.DISCONTINUOUS.,.T.,.T.)"
                ));
                let b = deck.add(&format!(
                    "SURFACE_PATCH('',#{window},.DISCONTINUOUS.,.DISCONTINUOUS.,.T.,.T.)"
                ));
                deck.add(&format!(
                    "RECTANGULAR_COMPOSITE_SURFACE('',((#{a}),(#{b})))"
                ))
            }
            _ => {
                // The floor's plane, moved up ten.
                let p = deck.add("CARTESIAN_POINT('',(0.,0.,0.))");
                let z = deck.add("DIRECTION('',(0.,0.,1.))");
                let x = deck.add("DIRECTION('',(1.,0.,0.))");
                let frame = deck.add(&format!("AXIS2_PLACEMENT_3D('',#{p},#{z},#{x})"));
                let floor = deck.add(&format!("PLANE('',#{frame})"));
                let up = deck.add("CARTESIAN_POINT('',(0.,0.,10.))");
                let operator = deck.add(&format!(
                    "CARTESIAN_TRANSFORMATION_OPERATOR_3D('',$,$,#{up},$,$)"
                ));
                deck.add(&format!("SURFACE_REPLICA('',#{floor},#{operator})"))
            }
        };
        deck.replace_in(face, &format!("#{plane},"), &format!("#{replacement},"));
        let got = deck.read_volume();
        assert!((got - 1000.0).abs() < 1e-6, "{spelling}: {got}");
    }
}
