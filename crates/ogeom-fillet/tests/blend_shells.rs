//! What a round or a bevel built through the public API takes off or adds,
//! measured on the solids it gives: each edge walked once each way, every
//! pcurve within its edge's stated tolerance, the meshes closed and true to
//! the exact measures, and the volumes against closed forms.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::{FastMap, Tolerances};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// How many edges two faces of `shape` walk the same way.
fn walked_one_way(model: &Model, shape: &Shape) -> usize {
    let mut walks: FastMap<_, Vec<Orientation>> = FastMap::default();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
        for wire in model.children_of(&face).unwrap() {
            for edge in model.children_of(&wire).unwrap() {
                if !model
                    .node(&edge)
                    .unwrap()
                    .data()
                    .as_edge()
                    .unwrap()
                    .degenerate
                {
                    walks
                        .entry(edge.node())
                        .or_default()
                        .push(edge.orientation());
                }
            }
        }
    }
    walks
        .values()
        .filter(|w| w.len() == 2 && w[0] == w[1])
        .count()
}

/// The edge of `shape` whose middle stands nearest `at`.
fn edge_at(model: &Model, shape: &Shape, at: Point) -> Shape {
    let off = |e: &Shape| {
        let b = ogeom_algo::shape_bounds(model, e, T).unwrap();
        b.low().unwrap().midpoint(b.high().unwrap()).distance(at)
    };
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| off(a).total_cmp(&off(b)))
        .unwrap()
}

/// The edge of `shape` whose curve is the first `kind` accepts.
fn edge_on(model: &Model, shape: &Shape, kind: impl Fn(&ogeom_geom::Curve) -> bool) -> Shape {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::placed_edge_curve(model, e, T)
                .unwrap()
                .is_some_and(|(curve, _)| kind(&curve))
        })
        .unwrap()
}

/// The volume of a solid, through the mesh where a face is fitted.
fn volume(model: &Model, solid: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, solid, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// Every edge of `shape` stands within its stated tolerance of each pcurve
/// it carries, measured along its whole length.
fn assert_pcurves_within(model: &Model, shape: &Shape, name: &str) {
    for e in explore_unique(model, shape, ShapeType::Edge).unwrap() {
        let stated = model
            .node(&e)
            .unwrap()
            .data()
            .as_edge()
            .unwrap()
            .tolerance
            .get();
        if let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &e, T).unwrap() {
            assert!(
                gap <= stated,
                "{name}: a pcurve stands {gap:e} off an edge stating {stated:e}"
            );
        }
    }
}

/// What `blended` takes off `solid` and what it adds to it, each as the
/// solids of a boolean.
fn changes(model: &mut Model, solid: &Shape, blended: &Shape) -> (Shape, Shape) {
    let taken = ogeom_bool::cut(model, solid, blended, T).unwrap().shape;
    let added = ogeom_bool::cut(model, blended, solid, T).unwrap().shape;
    (taken, added)
}

/// Whether `shape` holds no solid.
fn empty(model: &Model, shape: &Shape) -> bool {
    explore_unique(model, shape, ShapeType::Solid)
        .unwrap()
        .is_empty()
}

/// A straight edge's round and bevel and a drum rim's round walk each edge
/// of the solid they give once each way, and what each takes off walks
/// each edge so too and holds the volume the blend removes: a quarter
/// round's section, or a half square's, along the edge or round the rim at
/// the section's centroid.
#[test]
fn a_straight_round_and_bevel_and_a_drum_rim_round_walk_each_edge_once_each_way() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
        .unwrap()
        .shape;
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 4.0, T)
        .unwrap()
        .shape;
    let ridge = edge_at(&model, &block, Point::new(5.0, 0.0, 4.0));
    let rim = edge_at(&model, &drum, Point::new(0.0, 0.0, 4.0));
    let pi = core::f64::consts::PI;
    let section = 1.0 - pi / 4.0;
    let centroid = (5.0 / 6.0 - pi / 4.0) / section;
    let cases = [
        ("round", &block, &ridge, false, section * 10.0),
        ("bevel", &block, &ridge, true, 0.5 * 10.0),
        (
            "rim",
            &drum,
            &rim,
            false,
            2.0 * pi * (3.0 - centroid) * section,
        ),
    ];
    for (name, solid, edge, bevel, removed) in cases {
        let blended = if bevel {
            ogeom_fillet::chamfer_edge(&mut model, solid, edge, 1.0, T)
        } else {
            ogeom_fillet::fillet_edge(&mut model, solid, edge, 1.0, T)
        }
        .unwrap()
        .shape;
        assert_eq!(walked_one_way(&model, &blended), 0, "{name}");
        let (taken, added) = changes(&mut model, solid, &blended);
        assert!(empty(&model, &added), "{name}: a convex blend adds nothing");
        assert_eq!(walked_one_way(&model, &taken), 0, "{name}: what it takes");
        for got in [
            volume(&model, &taken),
            volume(&model, solid) - volume(&model, &blended),
        ] {
            assert!(
                (got - removed).abs() < removed * 1e-6,
                "{name}: {got} against {removed}"
            );
        }
    }
}

/// A straight round, on a convex edge and on a concave one, gives a solid
/// whose every pcurve stands within the tolerance its edge states.
#[test]
fn straight_rounds_hold_their_pcurves_within_their_edges() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
        .unwrap()
        .shape;
    let step = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 3.0, 8.0), T)
        .unwrap()
        .shape;
    let ell = ogeom_bool::fuse(&mut model, &block, &step, T)
        .unwrap()
        .shape;
    let ridge = edge_at(&model, &block, Point::new(5.0, 0.0, 4.0));
    let crease = edge_at(&model, &ell, Point::new(5.0, 3.0, 4.0));
    for (name, solid, edge) in [("convex", &block, &ridge), ("concave", &ell, &crease)] {
        let blended = ogeom_fillet::fillet_edge(&mut model, solid, edge, 1.0, T)
            .unwrap()
            .shape;
        assert!(ogeom_algo::check(&model, &blended, T).unwrap().is_valid());
        assert_pcurves_within(&model, &blended, name);
    }
}

/// The volume, validity and face, edge and vertex counts of a solid.
fn summary(model: &Model, shape: &Shape) -> (f64, bool, [usize; 3]) {
    let valid = ogeom_algo::check(model, shape, T).unwrap().is_valid();
    let count = |kind| explore_unique(model, shape, kind).unwrap().len();
    (
        volume(model, shape),
        valid,
        [
            count(ShapeType::Face),
            count(ShapeType::Edge),
            count(ShapeType::Vertex),
        ],
    )
}

/// The spandrel a round of radius `r` takes off the edge of a block along
/// x at y = 0, z = `top`, from x = 0 to `length`: the square corner of side
/// `r` less the quarter disc about (y, z) = (r, top - r), swept along x.
fn spandrel(model: &mut Model, length: f64, top: f64, r: f64) -> Shape {
    use ogeom_math::{Circle, Plane};
    let p = |y: f64, z: f64| Point::new(0.0, y, z);
    let [corner, on_top, on_front] = [p(0.0, top), p(r, top), p(0.0, top - r)]
        .map(|q| (ogeom_algo::make_vertex(model, q).shape, q));
    let segment = |model: &mut Model, a: &(Shape, Point), b: &(Shape, Point)| {
        let line = ogeom_geom::LineCurve::segment(a.1, b.1, T).unwrap();
        let curve = ogeom_geom::Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        ogeom_algo::make_edge_between(model, curve, domain, &a.0, &b.0, T)
            .unwrap()
            .shape
    };
    let along_top = segment(model, &corner, &on_top);
    let down_front = segment(model, &on_front, &corner);
    // The arc about (r, top - r) from the front contact up to the top's.
    let frame = Frame::new(p(r, top - r), -Direction::X, -Direction::Y, T).unwrap();
    let arc = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(
        Circle::new(frame, r, T).unwrap(),
    ));
    let arc = ogeom_algo::make_edge_between(
        model,
        arc,
        (0.0, core::f64::consts::FRAC_PI_2),
        &on_front.0,
        &on_top.0,
        T,
    )
    .unwrap()
    .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(p(0.0, 0.0), Direction::X),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(
        model,
        plane.into(),
        &[vec![along_top, arc.reversed(), down_front]],
        T,
    )
    .unwrap()
    .shape;
    let prism = ogeom_algo::make_prism(model, &face, Vector::new(length, 0.0, 0.0), T)
        .unwrap()
        .shape;
    // The prism's far end is its profile placed; the seated cut reads a
    // tool's own nodes, so the placement is baked in.
    ogeom_algo::baked_where_placed(model, &prism, T)
        .unwrap()
        .shape
}

/// A round's spandrel cuts the block it was built for by trimming the
/// faces it sits on. The same spandrel on the block with a blind hole
/// reaching up into it from below, whose wall its blend would cross, is
/// left to the general fuse. Either way the solid is the general fuse's.
#[test]
fn a_seated_spandrel_is_cut_by_trimming_only_where_nothing_reaches_into_it() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
        .unwrap()
        .shape;
    let tool = spandrel(&mut model, 10.0, 4.0, 1.0);
    let section = 1.0 - core::f64::consts::FRAC_PI_4;
    assert!((volume(&model, &tool) - 10.0 * section).abs() < 10.0 * section * 1e-9);
    // Up to a tenth below the ridge, a tenth in from the front: inside the
    // corner the round takes away.
    let at = Frame::new(Point::new(5.0, 0.1, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom_algo::make_cylinder(&mut model, at, 0.05, 4.9, T)
        .unwrap()
        .shape;
    let drilled = ogeom_bool::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    for (solid, seated) in [(&block, 1), (&drilled, 0)] {
        let mut general = model.clone();
        let want = ogeom_bool::cut(&mut general, solid, &tool, T).unwrap();
        let (got, work) = ogeom_bool::work::counted(|| {
            ogeom_bool::cut_seated(&mut model, solid, &tool, T).unwrap()
        });
        assert_eq!(work.seated_cuts, seated);
        let (want, got) = (summary(&general, &want.shape), summary(&model, &got.shape));
        assert!(got.1, "the cut is valid");
        assert_eq!((got.1, got.2), (want.1, want.2));
        assert!(
            (got.0 - want.0).abs() <= want.0 * 1e-9,
            "{} against {}",
            got.0,
            want.0
        );
    }
}

/// A stadium prism from `z0` up `height`: two runs of `length` along x
/// joined by semicircles of radius `r`, centred on (0, 0) and (`length`,
/// 0).
fn stadium(model: &mut Model, length: f64, r: f64, z0: f64, height: f64) -> Shape {
    use ogeom_math::{Circle, Plane};
    let corner = |x: f64, y: f64| Point::new(x, y, z0);
    let [tl, tr, br, bl] = [(0.0, r), (length, r), (length, -r), (0.0, -r)]
        .map(|(x, y)| ogeom_algo::make_vertex(model, corner(x, y)).shape);
    let arc = |model: &mut Model, cx: f64, x: Direction, from: &Shape, to: &Shape| {
        let frame = Frame::new(corner(cx, 0.0), Direction::Z, x, T).unwrap();
        let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(
            Circle::new(frame, r, T).unwrap(),
        ));
        ogeom_algo::make_edge_between(model, curve, (0.0, core::f64::consts::PI), from, to, T)
            .unwrap()
            .shape
    };
    let segment = |model: &mut Model, a: (&Shape, Point), b: (&Shape, Point)| {
        let line = ogeom_geom::LineCurve::segment(a.1, b.1, T).unwrap();
        let curve = ogeom_geom::Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        ogeom_algo::make_edge_between(model, curve, domain, a.0, b.0, T)
            .unwrap()
            .shape
    };
    let top = segment(model, (&tl, corner(0.0, r)), (&tr, corner(length, r)));
    let right = arc(model, length, -Direction::Y, &br, &tr);
    let bottom = segment(model, (&br, corner(length, -r)), (&bl, corner(0.0, -r)));
    let left = arc(model, 0.0, Direction::Y, &tl, &bl);
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(corner(0.0, 0.0), Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(
        model,
        plane.into(),
        &[vec![left, bottom.reversed(), right, top.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// An arc of a rim rounded, whichever side of the rim the cap and the
/// material stand on, gives a solid walking each edge once each way, and
/// what the round takes off or adds walks each edge so too and holds a
/// quarter round's section swept round the half turn at its centroid's
/// radius: inside the rim where the wall faces out, outside it where the
/// wall faces the axis or the rim is concave.
#[test]
fn rim_arc_rounds_walk_each_edge_once_each_way() {
    let mut model = Model::new();
    let (length, r) = (10.0, 5.0);
    let post = stadium(&mut model, length, r, 0.0, 4.0);
    let under =
        |z: f64| Frame::new(Point::new(-10.0, -10.0, z), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom_algo::make_box(&mut model, under(0.0), (30.0, 20.0, 4.0), T)
        .unwrap()
        .shape;
    let bore = stadium(&mut model, length, r, -1.0, 6.0);
    let slot = ogeom_bool::cut(&mut model, &block, &bore, T).unwrap().shape;
    let plate = ogeom_algo::make_box(&mut model, under(-4.0), (30.0, 20.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_bool::fuse(&mut model, &plate, &post, T)
        .unwrap()
        .shape;
    let pi = core::f64::consts::PI;
    let section = 1.0 - pi / 4.0;
    let centroid = (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi);
    let inside = pi * (r - centroid) * section;
    let outside = pi * (r + centroid) * section;
    let at = |z: f64| Point::new(length + r / 2.0, 0.0, z);
    let cases = [
        ("post top", &post, at(4.0), inside),
        ("post foot", &post, at(0.0), inside),
        ("slot top", &slot, at(4.0), outside),
        ("slot foot", &slot, at(0.0), outside),
        ("boss foot", &boss, at(0.0), outside),
    ];
    for (name, solid, near, expected) in cases {
        let rim = edge_at(&model, solid, near);
        let blended = ogeom_fillet::fillet_edge(&mut model, solid, &rim, 1.0, T)
            .unwrap()
            .shape;
        assert_eq!(walked_one_way(&model, &blended), 0, "{name}");
        let (taken, added) = changes(&mut model, solid, &blended);
        // A rim either takes or adds, never both.
        let change = match (empty(&model, &taken), empty(&model, &added)) {
            (false, true) => taken,
            (true, false) => added,
            _ => panic!("{name}: the round both takes and adds, or neither"),
        };
        assert_eq!(walked_one_way(&model, &change), 0, "{name}: its change");
        let moved = (volume(&model, &blended) - volume(&model, solid)).abs();
        for got in [volume(&model, &change), moved] {
            assert!(
                (got - expected).abs() < expected * 1e-6,
                "{name}: {got} against {expected}"
            );
        }
    }
}

/// A bore down a ball beside its pole, rounded where it leaves the far
/// side: the ball rolling round the rim touches the sphere along a rail
/// that passes over the pole, so the sphere's leg, the part of it the
/// round takes off, runs from the rail to the pole with the rim cut from
/// it. Its seam to the pole runs clear of the rim however the bore is
/// turned about the ball's axis: what the round takes off meshes closed to
/// what it measures, the leg's mesh covers its area, and the blend is
/// valid and measures the same at each turn.
#[test]
fn a_pole_leg_meshes_to_its_area_as_the_bore_turns() {
    let reference = pole_leg_blend(0.0);
    for turn in [0.9, 3.2, 5.6] {
        let whole = pole_leg_blend(turn);
        assert!(
            (whole - reference).abs() < reference * 1e-6,
            "turn {turn}: {whole} against {reference}"
        );
    }
}

/// The bore of [`a_pole_leg_meshes_to_its_area_as_the_bore_turns`] turned
/// so the rim's middle is its point nearest the pole: the side probe there
/// stands a tenth of a millimetre either side of the rim, half a
/// millimetre from the pole, and the ball's touch on the sphere lies
/// across the pole from it. The blend seats at another point of the rim
/// and measures what it does at the other turns.
#[test]
fn a_pole_leg_blends_with_the_rim_s_middle_nearest_the_pole() {
    let reference = pole_leg_blend(0.0);
    let whole = pole_leg_blend(4.7);
    assert!(
        (whole - reference).abs() < reference * 1e-6,
        "{whole} against {reference}"
    );
}

/// The far rim of a bore of radius 1.5 down a ball of radius 10, its axis
/// 2 from the ball's and `turn` round it, rounded at radius 1: checks the
/// mesh of what the round takes off, the sphere leg's mesh and the blend's
/// validity, and gives the blend's volume.
fn pole_leg_blend(turn: f64) -> f64 {
    let chord = 1e-2;
    let deflection = ogeom_mesh::Deflection::with_chord(chord).unwrap();
    let (c, s) = (turn.cos(), turn.sin());
    let mut model = Model::new();
    let ball = ogeom_algo::make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(
        Point::new(2.0 * c, 2.0 * s, -20.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let bore = ogeom_algo::make_cylinder(&mut model, frame, 1.5, 40.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom_bool::cut(&mut model, &ball, &bore, T).unwrap().shape;
    // The rim on the ball's top, on the bore's side away from the axis.
    let top = (100.0_f64 - 3.5 * 3.5).sqrt();
    let rim = edge_at(&model, &drilled, Point::new(3.5 * c, 3.5 * s, top));
    let blended = ogeom_fillet::fillet_edge(&mut model, &drilled, &rim, 1.0, T)
        .unwrap()
        .shape;
    assert!(
        ogeom_algo::check(&model, &blended, T).unwrap().is_valid(),
        "turn {turn}"
    );
    let (taken, added) = changes(&mut model, &drilled, &blended);
    assert!(
        empty(&model, &added),
        "turn {turn}: a convex rim adds nothing"
    );
    let exact = volume(&model, &taken);
    let mesh = ogeom_mesh::triangulate(&model, &taken, deflection, T).unwrap();
    assert!(mesh.is_closed(), "turn {turn}");
    assert!(
        (mesh.volume() - exact).abs() < exact * 0.05,
        "turn {turn}: the taken mesh {} against {exact}",
        mesh.volume()
    );
    let leg = explore_unique(&model, &taken, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|face| {
            model
                .node(face)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|g| matches!(g, ogeom_geom::SurfaceGeometry::Sphere(_)))
        })
        .unwrap();
    let area = ogeom_algo::surface_properties(&model, &leg, deflection, T)
        .unwrap()
        .mass;
    let leg_mesh = ogeom_mesh::triangulate_face(&model, &leg, deflection, T).unwrap();
    let meshed: f64 = leg_mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| leg_mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude() * 0.5
        })
        .sum();
    // A flat facet whose corners lie on a sphere of radius 10, sagging
    // under it by at most the chord, falls short of the sphere's area over
    // it by under chord / 20 of it; along the boundary each chord of an
    // edge stands within the chord of its arc, the sliver between them
    // under two thirds of the chord times the arc's length.
    let perimeter = ogeom_algo::linear_properties(&model, &leg, deflection, T)
        .unwrap()
        .mass;
    let within = area * chord / 20.0 + 2.0 / 3.0 * chord * perimeter;
    assert!(
        (meshed - area).abs() <= within,
        "turn {turn}: leg mesh {meshed} against {area}, within {within}"
    );
    volume(&model, &blended)
}

/// Where two bores cross inside a block, what the round of each edge where
/// they meet takes off meshes to what it measures. Round a closed loop
/// each leg is the patch of one bore's wall between two closed rings, both
/// on one branch of the wall's chart, so the mesh covers the patch and not
/// the wall round between them.
#[test]
fn a_band_round_crossing_bores_meshes_to_its_volume() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let bore = |model: &mut Model, at: Point, axis: Direction, x: Direction, r: f64| {
        let frame = Frame::new(at, axis, x, T).unwrap();
        ogeom_algo::make_cylinder(model, frame, r, 22.0, T)
            .unwrap()
            .shape
    };
    let a = bore(
        &mut model,
        Point::new(-1.0, 10.0, 10.0),
        Direction::X,
        Direction::Y,
        4.0,
    );
    let b = bore(
        &mut model,
        Point::new(10.0, -1.0, 12.0),
        Direction::Y,
        Direction::Z,
        3.0,
    );
    let once = ogeom_bool::cut(&mut model, &block, &a, T).unwrap().shape;
    let bored = ogeom_bool::cut(&mut model, &once, &b, T).unwrap().shape;
    let meeting: Vec<Shape> = explore_unique(&model, &bored, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            ogeom_algo::placed_edge_curve(&model, e, T)
                .unwrap()
                .is_some_and(|(curve, _)| matches!(curve, ogeom_geom::Curve::BSpline(_)))
        })
        .collect();
    assert_eq!(meeting.len(), 3);
    let chord = 1e-2;
    let deflection = ogeom_mesh::Deflection::with_chord(chord).unwrap();
    for edge in &meeting {
        let blended = ogeom_fillet::fillet_edge(&mut model, &bored, edge, 1.0, T)
            .unwrap()
            .shape;
        let (taken, added) = changes(&mut model, &bored, &blended);
        assert!(empty(&model, &added), "a convex crease adds nothing");
        let exact = volume(&model, &taken);
        let mesh = ogeom_mesh::triangulate(&model, &taken, deflection, T).unwrap();
        assert!(mesh.is_closed());
        // Every point of the boundary stands within the chord of the mesh,
        // so the two differ by at most the chord times the boundary's area:
        // a mesh covering a wall round between two rings would not.
        let area = ogeom_algo::surface_properties(&model, &taken, deflection, T)
            .unwrap()
            .mass;
        assert!(
            (mesh.volume() - exact).abs() <= chord * area,
            "mesh {} against {exact}, within {}",
            mesh.volume(),
            chord * area
        );
    }
}

/// A cylinder's axis: a point on it and its unit direction.
#[derive(Clone, Copy)]
struct Axis {
    at: Point,
    along: Vector,
}

impl Axis {
    /// The distance of `p` from the axis.
    fn rho(self, p: Point) -> f64 {
        let d = p - self.at;
        (d - self.along * d.dot(self.along)).magnitude()
    }

    /// The point of the cylinder of radius `radius` about the axis nearest
    /// `p`, off the axis: `p` moved along its radial line.
    fn foot(self, p: Point, radius: f64) -> Point {
        let d = p - self.at;
        let radial = d - self.along * d.dot(self.along);
        p + radial * (radius / radial.magnitude() - 1.0)
    }
}

/// One side of the region where a rolling ball's centre may stand: outside
/// a cylinder of the given radius, or on the positive side of a plane
/// (through `at`, unit normal `normal`).
#[derive(Clone, Copy)]
enum Side {
    Outside(Axis, f64),
    Above(Point, Vector),
}

impl Side {
    /// How far `c` stands inside this side (negative outside it).
    fn value(self, c: Point) -> f64 {
        match self {
            Side::Outside(axis, radius) => axis.rho(c) - radius,
            Side::Above(at, normal) => (c - at).dot(normal),
        }
    }

    /// The point of this side's boundary nearest `p`.
    fn foot(self, p: Point) -> Point {
        match self {
            Side::Outside(axis, radius) => axis.foot(p, radius),
            Side::Above(at, normal) => p - normal * (p - at).dot(normal),
        }
    }
}

/// A curve along which two sides' boundaries cross, by parameter, with a
/// first guess of the parameter of its point nearest a given point.
struct Spine {
    at: Box<dyn Fn(f64) -> Point>,
    guess: Box<dyn Fn(Point) -> f64>,
}

/// The region a rolling ball's centre may stand in, where it meets the
/// solid at most on its surface: the sides all hold, and the curves where
/// two sides' boundaries cross bound it.
struct Centres {
    sides: Vec<Side>,
    spines: Vec<Spine>,
}

impl Centres {
    /// The distance from `p` to the region. A closed set's nearest point
    /// to a point outside it lies on its boundary: on one side's boundary,
    /// where the nearest point of that whole boundary (the plane's foot,
    /// or the radial one for a cylinder, its only minimum) lies within the
    /// others, or else on a spine. Each spine's nearest point is searched
    /// within half a radian of the guess, past which every spine of these
    /// cases stands further from the band than three radii.
    fn distance(&self, p: Point) -> f64 {
        if self.sides.iter().all(|s| s.value(p) >= 0.0) {
            return 0.0;
        }
        let mut best = f64::INFINITY;
        for (i, side) in self.sides.iter().enumerate() {
            let foot = side.foot(p);
            let within = self
                .sides
                .iter()
                .enumerate()
                .all(|(j, s)| j == i || s.value(foot) >= -1e-12);
            if within {
                best = best.min(p.distance(foot));
            }
        }
        for spine in &self.spines {
            let g = (spine.guess)(p);
            let off = |t: f64| (spine.at)(t).distance(p);
            // Golden section to a parameter step of 1e-13: the distance is
            // then within rounding of the spine's nearest.
            let (mut a, mut b) = (g - 0.5, g + 0.5);
            let k = 0.5 * (5.0_f64.sqrt() - 1.0);
            let (mut c, mut d) = (b - k * (b - a), a + k * (b - a));
            let (mut fc, mut fd) = (off(c), off(d));
            while b - a > 1e-13 {
                if fc < fd {
                    (b, d, fd) = (d, c, fc);
                    c = b - k * (b - a);
                    fc = off(c);
                } else {
                    (a, c, fc) = (c, d, fd);
                    d = a + k * (b - a);
                    fd = off(d);
                }
            }
            best = best.min(fc.min(fd));
        }
        best
    }
}

/// The lower and upper bounds of `∫∫ ρ dρ dz` over the part of the box
/// `ρ0..ρ1`, `z0..z1` of a half-plane where `inside` is positive. `inside`
/// is 1-Lipschitz in the half-plane's own coordinates, which are lengths:
/// a cell whose centre's value exceeds its half-diagonal in size lies
/// wholly on that side. Cells are halved until decided or `least` wide;
/// undecided cells of that width count toward the upper bound alone.
fn section(
    inside: &dyn Fn(f64, f64) -> f64,
    (r0, r1): (f64, f64),
    (z0, z1): (f64, f64),
    least: f64,
) -> (f64, f64) {
    let (mut lo, mut hi) = (0.0, 0.0);
    let mut cells = vec![(r0, r1, z0, z1)];
    while let Some((a, b, c, d)) = cells.pop() {
        let half = 0.5 * (b - a).hypot(d - c);
        let value = inside(0.5 * (a + b), 0.5 * (c + d));
        let weight = 0.5 * (b * b - a * a) * (d - c);
        if value > half + 1e-9 {
            lo += weight;
            hi += weight;
        } else if value < -half - 1e-9 {
        } else if b - a <= least {
            hi += weight;
        } else {
            let (m, n) = (0.5 * (a + b), 0.5 * (c + d));
            cells.extend([(a, m, c, n), (m, b, c, n), (a, m, n, d), (m, b, n, d)]);
        }
    }
    (lo, hi)
}

/// The volume of a band by its sections in the half-planes about `axis`
/// (`e1`, `e2` square to it and to each other) at the angles `angles`, each
/// weighted as the quadrature over them weighs it: each section's box,
/// in distance from the axis and along it, and `inside` in space. Gives
/// the lower and upper bounds of the quadrature's sum.
fn band_volume(
    axis: Axis,
    (e1, e2): (Vector, Vector),
    angles: &[(f64, f64)],
    sections: &dyn Fn(f64) -> ((f64, f64), (f64, f64)),
    inside: &dyn Fn(Point) -> f64,
    least: f64,
) -> (f64, f64) {
    let (mut lo, mut hi) = (0.0, 0.0);
    for &(angle, weight) in angles {
        let out = e1 * angle.cos() + e2 * angle.sin();
        let at = |rho: f64, along: f64| inside(axis.at + axis.along * along + out * rho);
        let (rho, along) = sections(angle);
        let (l, h) = section(&at, rho, along, least);
        lo += weight * l;
        hi += weight * h;
    }
    (lo, hi)
}

/// An oracle's bounds on a band's volume: at the angles summed, and at
/// half as many.
type Oracle<'a> = dyn Fn() -> [(f64, f64); 2] + 'a;

/// The midpoints of `n` equal parts of `from..to`, each with its part's
/// width.
fn midpoints(from: f64, to: f64, n: u32) -> Vec<(f64, f64)> {
    let w = (to - from) / f64::from(n);
    (0..n)
        .map(|k| (from + w * (f64::from(k) + 0.5), w))
        .collect()
}

/// The rolling ball's band against the round built: the volume of what the
/// round adds or takes, measured on the boolean's solid through a mesh of
/// chord `chord`, within the sum of the bounds it can stray by of the
/// oracle's, sampled from the ball's envelope. The oracle's sections are
/// bounded rigorously; the quadrature over the angle is held to the
/// difference between `angles` and the half as many `coarse`; the round's
/// band strays by at most its fit target, and the mesh by its chord, over
/// the change's surface.
fn assert_band_volume(
    name: &str,
    model: &Model,
    change: &Shape,
    band: [(f64, f64); 2],
    chord: f64,
) {
    let deflection = ogeom_mesh::Deflection::with_chord(chord).unwrap();
    let got = ogeom_algo::volume_properties(model, change, deflection, T)
        .unwrap()
        .mass;
    let area = ogeom_algo::surface_properties(model, change, deflection, T)
        .unwrap()
        .mass;
    let [(lo, hi), (coarse_lo, coarse_hi)] = band;
    let mid = 0.5 * (lo + hi);
    let quadrature = (mid - 0.5 * (coarse_lo + coarse_hi)).abs();
    // The band's fit target (2e-4 at millimetre tolerances).
    let fit = 2e-4;
    let within = 0.5 * (hi - lo) + quadrature + (fit + chord) * area;
    eprintln!(
        "{name}: built {got}, oracle {mid} in [{lo}, {hi}], quadrature {quadrature}, within {within}"
    );
    assert!(
        (got - mid).abs() <= within,
        "{name}: the round moves {got}, the rolling ball {mid}, within {within}"
    );
}

/// A marched band rounded gives a solid walking each edge once each way,
/// every pcurve within its edge's tolerance, that passes the checker; the
/// round only adds, at a concave crease, or only takes, at a convex one;
/// and what it moves is the rolling ball's band: a closed band round a
/// branch's foot, a band running out through a block's walls at both ends,
/// and a band pinched where two equal drums turn tangent.
///
/// The oracle is the ball's envelope in the creases' own geometry, taken
/// whole (the drums and the floor unbounded): at a concave crease the
/// points outside the solid that no ball of the radius clear of it
/// reaches, at a convex one the points inside that no ball within it
/// reaches, kept to the solid's walls and to the crease asked. A point's
/// distance to the balls' centres is exact (the faces' feet, and the
/// nearest point of the curves where their offsets cross), so every bound
/// is 1-Lipschitz, and each section about the band's axis is bounded by
/// cells it decides; the sections are summed over the angle.
#[test]
#[ignore = "heavy"]
fn marched_rounds_walk_each_edge_once_each_way() {
    use core::f64::consts::{FRAC_PI_2, PI};
    use ogeom_geom::Curve;
    let mut model = Model::new();
    let across = |x: f64, y: f64, z: f64| {
        Frame::new(Point::new(x, y, z), Direction::X, Direction::Y, T).unwrap()
    };
    let main = ogeom_algo::make_cylinder(&mut model, across(-20.0, 0.0, 0.0), 10.0, 40.0, T)
        .unwrap()
        .shape;
    let branch = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let branched = ogeom_bool::fuse(&mut model, &main, &branch, T)
        .unwrap()
        .shape;
    let foot = edge_on(&model, &branched, |c| {
        use ogeom_geom::Curve3d as _;
        let (lo, hi) = c.domain();
        matches!(c, Curve::BSpline(_))
            && c.point_at(lo, T)
                .and_then(|p| c.point_at(hi, T).map(|q| p.distance(q)))
                .is_ok_and(|d| d <= 1e-6)
    });

    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let tilt = 0.35_f64;
    let (tc, ts) = (tilt.cos(), tilt.sin());
    let axis = Direction::new(Vector::new(0.0, tc, -ts), T).unwrap();
    let drum_at = Point::new(10.0, -5.0, 11.5);
    let drum = ogeom_algo::make_cylinder(
        &mut model,
        Frame::new(drum_at, axis, Direction::X, T).unwrap(),
        4.0,
        30.0,
        T,
    )
    .unwrap()
    .shape;
    let grooved = ogeom_bool::cut(&mut model, &block, &drum, T).unwrap().shape;
    let groove = edge_at(&model, &grooved, Point::new(10.0, 14.84, 0.0));

    let upright = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let lying = ogeom_algo::make_cylinder(&mut model, across(-20.0, 0.0, 10.0), 5.0, 40.0, T)
        .unwrap()
        .shape;
    let crossed = ogeom_bool::fuse(&mut model, &upright, &lying, T)
        .unwrap()
        .shape;
    let pinch = edge_on(&model, &crossed, |c| matches!(c, Curve::Ellipse(_)));

    // Sections no wider than a ten-thousandth decide the bounds.
    let least = 1e-4;
    let z_axis = Axis {
        at: Point::ORIGIN,
        along: Vector::Z,
    };

    // The closed band: outside the main drum (radius 10 about x) and the
    // branch (5 about z), the ball of 1.5 rolling between them.
    let closed = || {
        let r = 1.5;
        let main = Axis {
            at: Point::ORIGIN,
            along: Vector::X,
        };
        let centres = Centres {
            sides: vec![
                Side::Outside(main, 10.0 + r),
                Side::Outside(z_axis, 5.0 + r),
            ],
            spines: vec![Spine {
                at: Box::new(move |t: f64| {
                    let (x, y) = (6.5 * t.cos(), 6.5 * t.sin());
                    Point::new(x, y, (11.5_f64 * 11.5 - y * y).sqrt())
                }),
                guess: Box::new(|p: Point| p.y.atan2(p.x)),
            }],
        };
        // Above the axis: the branch's other crossing, below, is not the
        // solid's.
        let inside = move |p: Point| {
            (main.rho(p) - 10.0)
                .min(z_axis.rho(p) - 5.0)
                .min(centres.distance(p) - r)
                .min(p.z)
        };
        let sections = move |t: f64| {
            let z = (100.0 - 25.0 * t.sin().powi(2)).sqrt();
            ((5.0, 5.0 + 3.0 * r), (z - 3.0 * r, z + 3.0 * r))
        };
        // Symmetric across both planes through the branch's axis: a
        // quarter turn, four times. The sections repeat with the turn, so
        // the midpoint sum converges faster than any power.
        let quarter = |n| -> Vec<(f64, f64)> {
            midpoints(0.0, FRAC_PI_2, n)
                .into_iter()
                .map(|(t, w)| (t, 4.0 * w))
                .collect()
        };
        let e = (Vector::X, Vector::Y);
        [
            band_volume(z_axis, e, &quarter(16), &sections, &inside, least),
            band_volume(z_axis, e, &quarter(8), &sections, &inside, least),
        ]
    };

    // The run-out band: above the floor z = 0 and outside the tilted drum
    // (radius 4), the ball of 1 rolling between them, kept within the
    // block.
    let run_out = || {
        let r = 1.0;
        let drum = Axis {
            at: drum_at,
            along: Vector::new(0.0, tc, -ts),
        };
        // Square to the drum's axis: up the axis's plane, and along x.
        let (e1, e2) = (Vector::new(0.0, ts, tc), Vector::X);
        let spine = move |t: f64| {
            let out = e1 * t.cos() + e2 * t.sin();
            let along = (drum_at.z + (4.0 + r) * out.z - r) / ts;
            drum_at + drum.along * along + out * (4.0 + r)
        };
        let centres = Centres {
            sides: vec![
                Side::Above(Point::new(0.0, 0.0, r), Vector::Z),
                Side::Outside(drum, 4.0 + r),
            ],
            spines: vec![Spine {
                at: Box::new(spine),
                guess: Box::new(move |p: Point| {
                    let d = p - drum_at;
                    d.dot(e2).atan2(d.dot(e1))
                }),
            }],
        };
        let inside = move |p: Point| {
            p.z.min(drum.rho(p) - 4.0)
                .min(centres.distance(p) - r)
                .min(20.0 - p.y)
                .min(p.x)
                .min(20.0 - p.x)
                .min(10.0 - p.z)
        };
        // Where the drum's line at angle `t` meets the floor.
        let sections = move |t: f64| {
            let along = (drum_at.z + 4.0 * e1.z * t.cos()) / ts;
            ((4.0, 4.0 + 4.0 * r), (along - 12.0 * r, along + 4.0 * r))
        };
        // Symmetric across x = 10: half the run, twice. Past 1.6 from the
        // lowest line the crease stands beyond the back wall by more than
        // three radii.
        let half = |n| -> Vec<(f64, f64)> {
            midpoints(PI, PI + 1.6, n)
                .into_iter()
                .map(|(t, w)| (t, 2.0 * w))
                .collect()
        };
        [
            band_volume(drum, (e1, e2), &half(64), &sections, &inside, least),
            band_volume(drum, (e1, e2), &half(32), &sections, &inside, least),
        ]
    };

    // The pinched band: outside two drums of radius 5, upright about z and
    // lying about x at z = 10, the ball of 1 rolling between them. Their
    // crossing is two ellipses in the planes z - 10 = ±x, each split at
    // the pinches (0, ±5, 10) into arcs, and each arc at its middle, y = 0,
    // where its tangent is along y: the edge asked runs from a pinch to the
    // middle, on the sides of x, y and z = 10 its middle stands, and its
    // band ends in the plane y = 0, its section there.
    let pinched = |sx: f64, sy: f64, sz: f64| {
        let r = 1.0;
        let lying = Axis {
            at: Point::new(0.0, 0.0, 10.0),
            along: Vector::X,
        };
        let ellipse = |s: f64| Spine {
            at: Box::new(move |t: f64| {
                Point::new(6.0 * t.cos(), 6.0 * t.sin(), 10.0 + s * 6.0 * t.cos())
            }),
            guess: Box::new(|p: Point| p.y.atan2(p.x)),
        };
        let centres = Centres {
            sides: vec![
                Side::Outside(z_axis, 5.0 + r),
                Side::Outside(lying, 5.0 + r),
            ],
            spines: vec![ellipse(1.0), ellipse(-1.0)],
        };
        let inside = move |p: Point| {
            (z_axis.rho(p) - 5.0)
                .min(lying.rho(p) - 5.0)
                .min(centres.distance(p) - r)
                .min(sx * p.x)
                .min(sy * p.y)
                .min(sz * (p.z - 10.0))
        };
        let sections = move |t: f64| {
            let z = 10.0 + sz * 5.0 * t.cos();
            ((5.0, 5.0 + 3.0 * r), (z - 3.0 * r, z + 3.0 * r))
        };
        // From the arc's middle to a pinch.
        let half = |n| midpoints(0.0, FRAC_PI_2, n);
        let e = (Vector::X * sx, Vector::Y * sy);
        [
            band_volume(z_axis, e, &half(32), &sections, &inside, least),
            band_volume(z_axis, e, &half(16), &sections, &inside, least),
        ]
    };

    let middle = |model: &Model, edge: &Shape| {
        use ogeom_geom::Curve3d as _;
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(model, edge, T)
            .unwrap()
            .unwrap();
        c.point_at(0.5 * (a + b), T).unwrap()
    };
    let arc = middle(&model, &pinch);
    let ends = ogeom_algo::edge_vertices(&model, &pinch).unwrap().unwrap();
    let mut ends =
        [ends.0, ends.1].map(|v| model.node(&v).unwrap().data().as_vertex().unwrap().point);
    ends.sort_by(|a, b| a.x.abs().total_cmp(&b.x.abs()));
    let [at_pinch, at_middle] = ends;
    assert!(
        at_pinch.x.abs() < 1e-6 && (at_pinch.y.abs() - 5.0).abs() < 1e-6,
        "{at_pinch:?} is a pinch"
    );
    assert!(
        at_middle.y.abs() < 1e-6,
        "{at_middle:?} is the arc's middle"
    );
    let (sx, sy, sz) = (arc.x.signum(), arc.y.signum(), (arc.z - 10.0).signum());

    let cases: [(&str, &Shape, &Shape, f64, &Oracle<'_>); 3] = [
        ("closed band", &branched, &foot, 1.5, &closed),
        ("run-out band", &grooved, &groove, 1.0, &run_out),
        ("pinched band", &crossed, &pinch, 1.0, &|| {
            pinched(sx, sy, sz)
        }),
    ];
    for (name, solid, edge, radius, oracle) in cases {
        let blended = ogeom_fillet::fillet_edge(&mut model, solid, edge, radius, T)
            .unwrap()
            .shape;
        assert_eq!(walked_one_way(&model, &blended), 0, "{name}");
        assert_pcurves_within(&model, &blended, name);
        let found = ogeom_algo::check(&model, &blended, T).unwrap();
        assert!(found.is_usable(), "{name}: the blended solid: {found}");
        let (taken, added) = changes(&mut model, solid, &blended);
        let change = match (empty(&model, &taken), empty(&model, &added)) {
            (false, true) => taken,
            (true, false) => added,
            _ => panic!("{name}: the round both takes and adds, or neither"),
        };
        assert_band_volume(name, &model, &change, oracle(), 1e-4);
    }
}

/// A trapezoidal block's top front edge rounded, its band ending square on
/// one end face and on the other carried on to the oblique end face's
/// plane, which cuts it: the solid is valid, and the round takes off its
/// quarter section's area times the length of the section's centroid line
/// between the two end planes (a prism cut by planes holds its section
/// times its centroid line's length between them).
#[test]
fn a_band_cut_by_an_oblique_end_face_is_valid_and_holds_its_section() {
    use ogeom_math::Plane;
    let mut model = Model::new();
    // The trapezoid (0, 0), (10, 0), (7, 4), (0, 4), 2 thick: the oblique
    // end's plane is 0.8 x + 0.6 y = 8.
    let corners =
        [(0.0, 0.0), (10.0, 0.0), (7.0, 4.0), (0.0, 4.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let outline = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(&mut model, plane.into(), &[outline], T)
        .unwrap()
        .shape;
    let block = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 2.0), T)
        .unwrap()
        .shape;
    let r = 0.5;
    let edge = edge_at(&model, &block, Point::new(5.0, 0.0, 2.0));
    let blended = ogeom_fillet::fillet_edge(&mut model, &block, &edge, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom_algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let pi = core::f64::consts::PI;
    let section = (1.0 - pi / 4.0) * r * r;
    // The section's centroid stands this far in from the front and down
    // from the top; the oblique plane meets its line at x = 10 - 0.75 c.
    let c = r * (5.0 / 6.0 - pi / 4.0) / (1.0 - pi / 4.0);
    let removed = section * (10.0 - 0.75 * c);
    let got = volume(&model, &block).abs() - volume(&model, &blended).abs();
    assert!(
        (got - removed).abs() < removed * 1e-6,
        "{got} against {removed}"
    );
}
