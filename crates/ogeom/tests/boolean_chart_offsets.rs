//! Booleans against parts whose charts carry the small slips an imported
//! or edited part does: a spline fillet whose small ring closes its chart
//! only to the projection's accuracy, a torus whose trim sits a period off
//! its window both ways, an edge whose pcurve lags its curve, and a part
//! whose edges carry a few microns of slop. Each reads valid, and each
//! boolean against it is valid and measures what the common part says.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, make_torus, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A pin of radius `major + minor` standing `height` above its rounded
/// foot, as a STEP file: a bottom disc of radius `major`, a quarter of a
/// tube of radius `minor` round a circle of radius `major` carrying it out
/// to the wall, the wall and a top disc. The quarter tube is a rational
/// spline closed round the axis over `[-pi, pi]`, its tube angle over
/// `[-pi/2, 0]`, with one seam, and the file states no pcurves: the reader
/// projects every edge into it. Round the small bottom ring a radian of
/// the chart is `major` of the surface.
fn pin_step(major: f64, minor: f64, height: f64) -> String {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let outer = major + minor;
    // The tube's quarter profile in (radius, z), from the bottom ring out to
    // the wall, and the unit square a full circle's nine controls stand on,
    // starting at angle pi.
    let profile = [(major, -minor), (outer, -minor), (outer, 0.0)];
    let square = [
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
    ];
    let mut out = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\n\
         FILE_NAME('','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n\
         ENDSEC;\nDATA;\n",
    );
    let mut rows = Vec::new();
    let mut weights = Vec::new();
    for (i, (rho, z)) in profile.iter().enumerate() {
        let mut row = Vec::new();
        let mut wrow = Vec::new();
        for (j, (cx, cy)) in square.iter().enumerate() {
            let id = 1 + 9 * i + j;
            out.push_str(&format!(
                "#{id}=CARTESIAN_POINT('',({:?},{:?},{:?}));\n",
                rho * cx,
                rho * cy,
                z
            ));
            row.push(format!("#{id}"));
            let w = if i == 1 { s } else { 1.0 } * if j % 2 == 1 { s } else { 1.0 };
            wrow.push(format!("{w:?}"));
        }
        rows.push(format!("({})", row.join(",")));
        weights.push(format!("({})", wrow.join(",")));
    }
    out.push_str(&format!(
        "#30=(BOUNDED_SURFACE()B_SPLINE_SURFACE(2,2,({}),.UNSPECIFIED.,.F.,.F.,.F.)\
         B_SPLINE_SURFACE_WITH_KNOTS((3,3),(3,2,2,2,3),(-1.57079632679549,0.),\
         (-3.14159265358979,-1.5707963267949,0.,1.5707963267949,3.14159265358979),\
         .UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(({}))\
         REPRESENTATION_ITEM('')SURFACE());\n",
        rows.join(","),
        weights.join(",")
    ));
    out.push_str(&format!(
        "#40=DIRECTION('',(0.,0.,1.));\n#41=DIRECTION('',(1.,0.,0.));\n\
         #42=DIRECTION('',(0.,1.,0.));\n#43=DIRECTION('',(-1.,0.,0.));\n\
         #50=CARTESIAN_POINT('',({m:?},0.,{d:?}));\n#51=VERTEX_POINT('',#50);\n\
         #52=CARTESIAN_POINT('',({o:?},0.,0.));\n#53=VERTEX_POINT('',#52);\n\
         #54=CARTESIAN_POINT('',({o:?},0.,{height:?}));\n#55=VERTEX_POINT('',#54);\n\
         #60=CARTESIAN_POINT('',(0.,0.,{d:?}));\n#61=AXIS2_PLACEMENT_3D('',#60,#40,#41);\n\
         #62=CIRCLE('',#61,{major:?});\n#63=EDGE_CURVE('',#51,#51,#62,.T.);\n\
         #64=CARTESIAN_POINT('',({m:?},0.,0.));\n#65=AXIS2_PLACEMENT_3D('',#64,#42,#43);\n\
         #66=CIRCLE('',#65,{minor:?});\n#67=EDGE_CURVE('',#51,#53,#66,.T.);\n\
         #68=CARTESIAN_POINT('',(0.,0.,0.));\n#69=AXIS2_PLACEMENT_3D('',#68,#40,#41);\n\
         #70=CIRCLE('',#69,{outer:?});\n#71=EDGE_CURVE('',#53,#53,#70,.T.);\n\
         #72=CARTESIAN_POINT('',(0.,0.,{height:?}));\n#73=AXIS2_PLACEMENT_3D('',#72,#40,#41);\n\
         #74=CIRCLE('',#73,{outer:?});\n#75=EDGE_CURVE('',#55,#55,#74,.T.);\n\
         #76=VECTOR('',#40,1.);\n#77=LINE('',#52,#76);\n#78=EDGE_CURVE('',#53,#55,#77,.T.);\n\
         #80=ORIENTED_EDGE('',*,*,#63,.T.);\n#81=ORIENTED_EDGE('',*,*,#67,.T.);\n\
         #82=ORIENTED_EDGE('',*,*,#71,.F.);\n#83=ORIENTED_EDGE('',*,*,#67,.F.);\n\
         #84=EDGE_LOOP('',(#80,#81,#82,#83));\n#85=FACE_OUTER_BOUND('',#84,.T.);\n\
         #86=ADVANCED_FACE('',(#85),#30,.F.);\n\
         #90=PLANE('',#61);\n#92=ORIENTED_EDGE('',*,*,#63,.F.);\n#93=EDGE_LOOP('',(#92));\n\
         #94=FACE_OUTER_BOUND('',#93,.T.);\n#95=ADVANCED_FACE('',(#94),#90,.F.);\n\
         #101=CYLINDRICAL_SURFACE('',#69,{outer:?});\n\
         #102=ORIENTED_EDGE('',*,*,#71,.T.);\n#103=ORIENTED_EDGE('',*,*,#78,.T.);\n\
         #104=ORIENTED_EDGE('',*,*,#75,.F.);\n#105=ORIENTED_EDGE('',*,*,#78,.F.);\n\
         #106=EDGE_LOOP('',(#102,#103,#104,#105));\n#107=FACE_OUTER_BOUND('',#106,.T.);\n\
         #108=ADVANCED_FACE('',(#107),#101,.T.);\n\
         #110=PLANE('',#73);\n#111=ORIENTED_EDGE('',*,*,#75,.T.);\n#112=EDGE_LOOP('',(#111));\n\
         #113=FACE_OUTER_BOUND('',#112,.T.);\n#114=ADVANCED_FACE('',(#113),#110,.T.);\n\
         #120=CLOSED_SHELL('',(#86,#95,#108,#114));\n#121=MANIFOLD_SOLID_BREP('',#120);\n\
         ENDSEC;\nEND-ISO-10303-21;\n",
        m = -major,
        d = -minor,
        o = -outer,
    ));
    out
}

/// The pin's volume: the foot, a solid of revolution whose radius runs
/// `major + sqrt(minor^2 - z^2)` over `z` in `[-minor, 0]`, and the drum
/// above it.
fn pin_volume(major: f64, minor: f64, height: f64) -> f64 {
    let foot = PI
        * (major * major * minor
            + major * PI * minor * minor / 2.0
            + 2.0 * minor * minor * minor / 3.0);
    foot + PI * (major + minor).powi(2) * height
}

/// The three booleans of a part and a tool, each valid, measured against
/// the common one: what a cut leaves and a fuse holds.
fn all_three(model: &mut Model, part: &Shape, tool: &Shape, part_volume: f64) -> f64 {
    let tool_volume = volume(model, tool);
    let common = ogeom::boolean::common(model, part, tool, T).unwrap().shape;
    let cut = ogeom::boolean::cut(model, part, tool, T).unwrap().shape;
    let fused = ogeom::boolean::fuse(model, part, tool, T).unwrap().shape;
    for (name, shape) in [("common", &common), ("cut", &cut), ("fuse", &fused)] {
        let found = check(model, shape, T).unwrap();
        assert!(found.problems.is_empty(), "{name}: {:?}", found.problems);
    }
    let shared = volume(model, &common);
    let left = volume(model, &cut);
    let held = volume(model, &fused);
    let slack = 1e-6 * (part_volume + tool_volume);
    assert!(
        (left - (part_volume - shared)).abs() <= slack,
        "cut {left} against {part_volume} less {shared}"
    );
    assert!(
        (held - (part_volume + tool_volume - shared)).abs() <= slack,
        "fuse {held} against {part_volume} and {tool_volume} less {shared}"
    );
    shared
}

/// The pin's foot and height: a ring a sixth of a millimetre in radius
/// under a fillet of three quarters of one.
const PIN: (f64, f64, f64) = (0.16, 0.75, 2.0);

/// The pin as read, valid and of its own volume.
fn read_pin() -> (ogeom::doc::Document, Shape, f64) {
    let (major, minor, height) = PIN;
    let import = ogeom::io::step::read_step(&pin_step(major, minor, height), T).unwrap();
    assert!(import.report.untrimmed_faces.is_empty());
    let pin = import.solids[0].clone();
    let document = import.document;
    let model = document.model();
    assert!(check(model, &pin, T).unwrap().problems.is_empty());
    let whole = volume(model, &pin);
    let expected = pin_volume(major, minor, height);
    assert!(
        (whole - expected).abs() <= 1e-9 * expected,
        "the pin as read: {whole} against {expected}"
    );
    (document, pin, whole)
}

/// A spline fillet whose small ring closes its chart only to the
/// projection's accuracy takes part in every boolean.
///
/// The ring's image ends where the projection of its vertex stopped. A
/// projection that stops once its residual (the gap times the chart's
/// speed) is under the confusion leaves that end several confusions short
/// along a ring a sixth of a millimetre in radius, microradians off the
/// seam column, and the fillet's boundary then reads as not closing round the
/// seam: every boolean that touches the pin refuses with a dangling strand,
/// wherever the tool stands, though the pin reads valid.
#[test]
fn a_spline_fillet_with_a_small_ring_takes_part_in_booleans() {
    let (major, minor, _) = PIN;
    let (mut document, pin, whole) = read_pin();
    let model = document.model_mut();

    // Across the top, far from the fillet: a slab of the drum, in closed form.
    let lid = make_box(
        model,
        Frame::WORLD.with_origin(Point::new(-2.0, -2.0, 1.5)),
        (4.0, 4.0, 1.0),
        T,
    )
    .unwrap()
    .shape;
    let shared = all_three(model, &pin, &lid, whole);
    let slab = PI * (major + minor).powi(2) * 0.5;
    assert!(
        (shared - slab).abs() <= 1e-6 * slab,
        "{shared} against {slab}"
    );

    // A slot through the foot across the fillet's seam, and a drum through
    // the fillet beside the seam.
    let slot = make_box(
        model,
        Frame::WORLD.with_origin(Point::new(-0.6, -2.0, -1.0)),
        (0.4, 4.0, 0.6),
        T,
    )
    .unwrap()
    .shape;
    assert!(all_three(model, &pin, &slot, whole) > 0.0);
    let drill = make_cylinder(
        model,
        Frame::WORLD.with_origin(Point::new(-0.8, 0.0, -1.5)),
        0.2,
        4.0,
        T,
    )
    .unwrap()
    .shape;
    assert!(all_three(model, &pin, &drill, whole) > 0.0);
}

/// The same fillet under tools at more places: from the side across the
/// seam, from below off it, a shallow pocket in the foot, and drums down
/// the axis and through the wall.
#[test]
#[ignore = "heavy"]
fn a_spline_fillet_with_a_small_ring_takes_tools_anywhere() {
    let (mut document, pin, whole) = read_pin();
    let model = document.model_mut();
    let boxes = [
        ((-2.0, -0.25, -0.5), (1.6, 0.5, 1.0)),
        ((0.3, -0.25, -1.0), (1.0, 0.5, 1.0)),
        ((-0.5, -0.5, -0.4), (1.0, 1.0, 0.3)),
    ];
    for ((x, y, z), size) in boxes {
        let tool = make_box(
            model,
            Frame::WORLD.with_origin(Point::new(x, y, z)),
            size,
            T,
        )
        .unwrap()
        .shape;
        let shared = all_three(model, &pin, &tool, whole);
        assert!(
            shared > 0.0,
            "the box at ({x}, {y}, {z}) takes some of the pin"
        );
    }
    for (x, y, radius) in [(0.0, 0.0, 0.1), (0.5, 0.5, 0.3)] {
        let tool = make_cylinder(
            model,
            Frame::WORLD.with_origin(Point::new(x, y, -1.5)),
            radius,
            4.0,
            T,
        )
        .unwrap()
        .shape;
        let shared = all_three(model, &pin, &tool, whole);
        assert!(shared > 0.0, "the drum at ({x}, {y}) takes some of the pin");
    }
}

/// A torus whose trim sits a period off its window in both directions is
/// cut as the same torus at home.
///
/// A file's torus face can run its rims' circles past a turn and its tube
/// angle below zero, and the reader keeps the chart where the file put it.
/// A section's image, folded into the window, then reaches the trim only
/// by a shift of a period along both directions at once. Tried one
/// direction at a time it never lands inside, the section is left off the
/// face, and the cut does not close.
#[test]
fn a_torus_chart_a_period_off_both_ways_is_cut_as_at_home() {
    let cut_through = |shifted: bool| -> (f64, f64) {
        let mut model = Model::new();
        let ring = make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
            .unwrap()
            .shape;
        if shifted {
            let away = ogeom::math::Transform2::translation(ogeom::math::Vector2::new(
                2.0 * PI,
                -2.0 * PI,
            ));
            let mut ids = Vec::new();
            for edge in explore_unique(&model, &ring, ShapeType::Edge).unwrap() {
                let data = model.node(&edge).unwrap().data().as_edge().unwrap();
                for repr in &data.representations {
                    match repr {
                        EdgeRepr::PCurve { curve, .. } => ids.push(*curve),
                        EdgeRepr::Seam {
                            forward, reversed, ..
                        } => ids.extend([*forward, *reversed]),
                        _ => {}
                    }
                }
            }
            ids.sort();
            ids.dedup();
            for id in ids {
                let held = model.geometry_mut().pcurve_mut(id).unwrap();
                *held = held.transformed(&away, T).unwrap();
            }
            assert!(check(&model, &ring, T).unwrap().problems.is_empty());
        }
        let whole = volume(&model, &ring);
        let mut shared = 0.0;
        for (x, y, radius) in [(10.0, 0.0, 2.0), (-7.0, 7.5, 1.5)] {
            let tool = make_cylinder(
                &mut model,
                Frame::WORLD.with_origin(Point::new(x, y, -5.0)),
                radius,
                10.0,
                T,
            )
            .unwrap()
            .shape;
            shared += all_three(&mut model, &ring, &tool, whole);
        }
        (whole, shared)
    };
    let (home, home_shared) = cut_through(false);
    let (away, away_shared) = cut_through(true);
    assert!((home - away).abs() <= 1e-9 * home, "{away} against {home}");
    assert!(home_shared > 0.0);
    assert!(
        (away_shared - home_shared).abs() <= 1e-6 * home_shared,
        "{away_shared} against {home_shared}"
    );
}

/// An edge whose pcurve lies on its curve but lags it, cut short, states
/// the lag at its new ends.
///
/// The checker holds a pcurve to the curve's points, not to its pace, so a
/// pcurve running microns behind its curve's parameter reads valid on the
/// whole edge. A piece of it ends at one parameter on both, and the lag
/// then stands between the piece's vertex and its image in the chart: the
/// piece must state it, or the cut it comes from reads broken.
#[test]
fn a_piece_of_an_edge_whose_pcurve_lags_states_the_lag() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    // The top face's edge along `x` at `y = 0`: its image there is redrawn
    // along the same chart line, half a micron behind at its middle: five
    // times the edge's tolerance, and within the arrangement's weld.
    let top = explore_unique(&model, &block, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            ogeom::algo::shape_bounds(&model, f, T)
                .unwrap()
                .low()
                .is_some_and(|low| low.z > 9.0)
        })
        .unwrap();
    let Some(NodeData::Face(face)) = model.node(&top).map(|n| n.data()) else {
        panic!("the top is a face");
    };
    let plane = face.surface;
    let edge = explore_unique(&model, &top, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom::algo::shape_bounds(&model, e, T).unwrap();
            b.high().is_some_and(|h| h.x > 9.0 && h.y.abs() < 1e-3)
        })
        .unwrap();
    let data = model.node(&edge).unwrap().data().as_edge().unwrap();
    let (id, range) = data
        .representations
        .iter()
        .find_map(|r| match r {
            EdgeRepr::PCurve {
                curve,
                surface,
                range,
                ..
            } if *surface == plane => Some((*curve, *range)),
            _ => None,
        })
        .unwrap();
    use ogeom::geom::Curve2d as _;
    let image = model.geometry().pcurve(id).unwrap().clone();
    let (a, b) = (
        image.point_at(range.0, T).unwrap(),
        image.point_at(range.1, T).unwrap(),
    );
    let behind = 1e-6 / a.distance(b);
    let middle = ogeom::math::Point2::new(
        (b.x - a.x).mul_add(-behind, f64::midpoint(a.x, b.x)),
        (b.y - a.y).mul_add(-behind, f64::midpoint(a.y, b.y)),
    );
    let knots = ogeom::math::KnotVector::new(
        vec![range.0, range.0, range.0, range.1, range.1, range.1],
        2,
    )
    .unwrap();
    *model.geometry_mut().pcurve_mut(id).unwrap() =
        ogeom::geom::BSpline2d::new(knots, vec![a, middle, b], T)
            .unwrap()
            .into();
    assert!(
        check(&model, &block, T).unwrap().problems.is_empty(),
        "a pcurve on its curve reads valid whatever its pace"
    );
    // A notch across the middle of that edge.
    let notch = make_box(
        &mut model,
        Frame::WORLD.with_origin(Point::new(4.0, -1.0, 9.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap()
    .shape;
    let cut = ogeom::boolean::cut(&mut model, &block, &notch, T)
        .unwrap()
        .shape;
    let found = check(&model, &cut, T).unwrap();
    assert!(found.problems.is_empty(), "{:?}", found.problems);
    let left = volume(&model, &cut);
    assert!((left - 998.0).abs() <= 1e-6, "{left}");
}

/// A tool standing in a part's hole, its surfaces crossing the part's,
/// leaves the part as it was though the part's edges carry microns of slop.
///
/// Every piece of both reads outside the other, and the boolean asks the
/// two boundaries, sampled, whether they share volume. The samples are read
/// at the part's own looseness; drawn at ten thousand times it, the trims
/// they are read against take tens of millimetres a chord, every ray lands
/// beside some trim, and the question goes unanswered.
#[test]
fn a_tool_in_a_loose_part_s_hole_leaves_it_whole() {
    let mut model = Model::new();
    let slab = make_box(&mut model, Frame::WORLD, (30.0, 30.0, 10.0), T)
        .unwrap()
        .shape;
    let window = make_box(
        &mut model,
        Frame::WORLD.with_origin(Point::new(10.0, 10.0, -1.0)),
        (10.0, 10.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let frame = ogeom::boolean::cut(&mut model, &slab, &window, T)
        .unwrap()
        .shape;
    let slop = ogeom::core::Tolerance::new(3e-3).unwrap();
    for kind in [ShapeType::Edge, ShapeType::Vertex] {
        for sub in explore_unique(&model, &frame, kind).unwrap() {
            match model.node_mut(&sub).unwrap().data_mut() {
                NodeData::Edge(data) => data.widen(slop),
                NodeData::Vertex(data) => data.widen(slop),
                _ => {}
            }
        }
    }
    assert!(check(&model, &frame, T).unwrap().problems.is_empty());
    let whole = volume(&model, &frame);
    assert!((whole - 8000.0).abs() <= 1e-6, "{whole}");
    let post = make_box(
        &mut model,
        Frame::WORLD.with_origin(Point::new(13.0, 13.0, -2.0)),
        (4.0, 4.0, 14.0),
        T,
    )
    .unwrap()
    .shape;
    assert!(all_three(&mut model, &frame, &post, whole).abs() <= 1e-9);
}
