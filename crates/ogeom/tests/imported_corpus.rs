//! M3's closing argument: the corpus imported, healed, measured, and
//! operated on — real files in, real modelling out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::topo::{Filter, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// Healing sweeps the whole corpus: every part reads, every shell closes,
/// every part measures, and each volume pins to its known figure — the
/// values ftc_07 and ftc_11 arbitrated against the kernel's exact ray
/// classifier and ctc_01 against an orientation-free even-odd grid. A loose
/// relative band absorbs future mesh refinements; an orientation or unit
/// mistake moves a volume by whole factors and cannot hide inside it.
#[test]
fn the_corpus_heals_and_every_part_measures() {
    let files: [(&str, f64); 11] = [
        ("nist_ctc_01_asme1_rd.stp", 14_643_073.4),
        ("nist_ctc_02_asme1_rc.stp", 47_101_909.1),
        ("nist_ctc_03_asme1_rc.stp", 331_884.8),
        ("nist_ctc_04_asme1_rd.stp", 17_519_472.2),
        ("nist_ctc_05_asme1_rd.stp", 12_694_721.1),
        ("nist_ftc_06_asme1_rd.stp", 3_289_989.2),
        ("nist_ftc_07_asme1_rd.stp", 1_726_286.9),
        ("nist_ftc_08_asme1_rc.stp", 503_596.4),
        ("nist_ftc_09_asme1_rd.stp", 136_453.8),
        ("nist_ftc_10_asme1_rb.stp", 188_388.8),
        ("nist_ftc_11_asme1_rb.stp", 5_122.3),
    ];
    let fine = ogeom::mesh::Deflection {
        chord: 1e-2,
        ..ogeom::mesh::Deflection::default()
    };
    for (name, expected) in files {
        let text = corpus(name);
        let mut import = ogeom::io::read_step(&text, T).unwrap();
        let solid = import.solids[0].clone();
        // Healing is idempotent where nothing is broken and surgery where
        // something is; either way the shell must close.
        let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
            .map_or(solid, |h| h.0.shape);
        let shell = explore_unique(import.document.model(), &healed, ShapeType::Shell)
            .unwrap()
            .remove(0);
        assert!(
            ogeom::algo::is_shell_closed(import.document.model(), &shell).unwrap(),
            "{name}: the healed shell closes"
        );
        let props = ogeom::algo::volume_properties(import.document.model(), &healed, fine, T)
            .unwrap_or_else(|e| panic!("{name}: does not measure: {e}"));
        eprintln!("REPORT {name}: volume {:.3} mm^3", props.mass);
        assert!(
            (props.mass - expected).abs() <= expected * 1e-2,
            "{name}: volume {:.3} strays from its pinned {expected:.1}",
            props.mass
        );
    }
}

/// A boolean over imported geometry, history checked: the milestone's
/// operations run on the world's parts, not only on this kernel's own.
#[test]
fn an_imported_part_takes_a_boolean_cut() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let mut import = ogeom::io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
        .unwrap()
        .0
        .shape;
    let fine = ogeom::mesh::Deflection {
        chord: 1e-2,
        ..ogeom::mesh::Deflection::default()
    };
    let before = ogeom::algo::volume_properties(import.document.model(), &healed, fine, T)
        .unwrap()
        .mass;

    // A square post cut down through the plate's solid ring — the part has a
    // large central pocket, and a post through fresh air cuts nothing, as an
    // earlier version of this test discovered only once the result's mesh
    // first became measurable.
    let frame = ogeom::math::Frame::new(
        ogeom::math::Point::new(20.0, -4.0, -3.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let post =
        ogeom::algo::make_box(import.document.model_mut(), frame, (8.0, 8.0, 6.0), T).unwrap();
    let result = ogeom::boolean::cut(import.document.model_mut(), &healed, &post.shape, T).unwrap();

    // The cut runs on the imported part: the result is a solid whose shell
    // closes, built from pieces of the world's geometry and this kernel's,
    // and it *measures* — the mesh welds across the file's slop because the
    // slop is recorded on the edges and vertices and the weld honours it.
    let shell = explore_unique(import.document.model(), &result.shape, ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(ogeom::algo::is_shell_closed(import.document.model(), &shell).unwrap());
    let props =
        ogeom::algo::volume_properties(import.document.model(), &result.shape, fine, T).unwrap();
    eprintln!("REPORT cut imported: {before:.3} -> {:.3} mm^3", props.mass);
    assert!(props.mass < before, "the cut removed material");
    assert!(props.mass > 0.0);
    // The post removes at most its own volume, and the difference cannot
    // exceed it: the bound the arithmetic itself provides.
    assert!(
        before - props.mass <= 8.0 * 8.0 * 6.0,
        "the cut removed more than the post could"
    );

    // History carried through: the imported solid is recorded as modified
    // into the result, and at least one of its faces was split or consumed.
    assert_eq!(
        result.history.modified(&healed),
        std::slice::from_ref(&result.shape)
    );
    let touched = explore(
        import.document.model(),
        &healed,
        Filter::OfType(ShapeType::Face),
    )
    .unwrap()
    .iter()
    .filter(|f| result.history.is_affected(f))
    .count();
    assert!(touched > 0, "faces of the imported part appear in history");
}

/// A body is bounded by what it is trimmed to, not by the carriers under it.
///
/// `shape_bounds` promised a guarantee and delivered the carriers: an
/// imported plane reported its own window, which spans kilometres, and a
/// cylinder its height domain. A consumer asking how big a screw is got a
/// box a billion millimetres across and had to fall back to the hull of the
/// topological vertices — which does not contain the body either, since a
/// button head's apex is a bulge between its rims, three millimetres past
/// every vertex the head has.
///
/// Two things were wrong. An edge reported the bound of its whole curve
/// rather than of the range it uses, so a segment on a line that runs to
/// the ends of the world reported the ends of the world. And a curved face
/// reported its whole surface, when its boundary is bounded by the edges
/// below it and all the face itself has to add is where its surface bulges
/// past that boundary.
#[test]
fn a_body_is_bounded_by_what_it_is_trimmed_to() {
    let deflection = ogeom::mesh::Deflection::default();
    // The same screw in two encodings, five bodies in one and one in the
    // other, and the second is where the figures below come from: an M5×16
    // button head, 18.75 long over a head 9.5 across.
    for (file, bodies) in [("m5x16_bhcs.step", 5), ("m5x16_bhcs_loops.step", 1)] {
        let text = corpus(file);
        let import = ogeom::io::read_step(&text, T).unwrap();
        let model = import.document.model();
        assert_eq!(import.solids.len(), bodies, "{file}: its bodies");
        for (which, solid) in import.solids.iter().enumerate() {
            let bounds = ogeom::algo::shape_bounds(model, solid, T).unwrap();
            let mesh = ogeom::mesh::triangulate(model, solid, deflection, T).unwrap();
            let mut hull = ogeom::math::Aabb::default();
            for at in &mesh.positions {
                hull = hull.with_point(*at);
            }
            assert!(
                bounds.contains_box(&hull),
                "{file} body {which}: the bound holds the mesh, {:?}..{:?} against {:?}..{:?}",
                bounds.low(),
                bounds.high(),
                hull.low(),
                hull.high()
            );
            // And holds little else. Where the bound stands proud of the mesh it
            // is the mesh that is short — a chord across an arc falls inside it
            // — so the margin is the chord's, either way round.
            let (low, high) = (bounds.low().unwrap(), bounds.high().unwrap());
            let (near, far) = (hull.low().unwrap(), hull.high().unwrap());
            for (bound, meshed) in [
                (low.x, near.x),
                (low.y, near.y),
                (low.z, near.z),
                (high.x, far.x),
                (high.y, far.y),
                (high.z, far.z),
            ] {
                assert!(
                    (bound - meshed).abs() <= deflection.chord,
                    "{file} body {which}: the bound stands {} from the mesh",
                    (bound - meshed).abs()
                );
            }
            // The screw itself, and not the carriers under it.
            let size = high - low;
            assert!(
                size.x < 20.0 && (size.y - 9.5).abs() < 0.1 && (size.z - 9.5).abs() < 0.1,
                "{file} body {which}: an M5×16 button head, {size:?}"
            );
        }
    }
}

/// A sliver's boundary is drawn finely enough to be a boundary.
///
/// The face is a quarter-arc forty-five millimetres long and eighteen
/// microns wide, between two nearly concentric circles, lifted out of a
/// community printer assembly. At the default tenth of a millimetre the
/// sagitta of each bounding arc is twenty-nine microns — wider than the
/// region itself — so the inner polyline crosses the outer one and what
/// reaches the triangulator is not a region. It answered with sixteen
/// triangles in fifteen disconnected pieces, and the holes between them
/// were what kept the body it belongs to from meshing closed.
///
/// `V - E + F` is the test because it is the question: one for a disc,
/// `1 - holes` for a face with inner loops, and anything else means the
/// pieces are not joined. Fifteen is fifteen fragments.
#[test]
fn a_sliver_face_is_drawn_fine_enough_to_triangulate_whole() {
    use std::collections::HashMap;
    let text = corpus("sliver_face_falls_apart.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");

    let wires = model.ordered_children_of(&faces[0]).unwrap().len();
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let euler = mesh.positions.len() as i64 - uses.len() as i64 + mesh.triangles.len() as i64;
    assert_eq!(
        euler,
        2 - wires as i64,
        "one connected piece, not {} of them: {} verts, {} tris",
        euler,
        mesh.positions.len(),
        mesh.triangles.len()
    );

    // And drawn finer than asked, which is the point: the caller's chord
    // would have crossed the boundary with itself.
    assert!(
        mesh.triangles.len() >= 32,
        "refined past the caller's chord: {} triangles",
        mesh.triangles.len()
    );
}

/// A run along a cone's apex row is boundary, whatever angle it spans.
///
/// The face is a cone sector bounded by two rulings into the apex and one
/// arc, the rulings standing exactly a quarter turn apart. In the chart the
/// two apex ends are distinct points a quarter period apart; in space they
/// are one vertex. The ring has to keep both — the run between them along
/// the degenerate row is the face's own boundary — and a test that read
/// "apart" as "more than a quarter period" dropped the second, cut the
/// corner through the face, and lost the triangle at the apex. The face
/// was still a disc; it was a third smaller than it should be, and the
/// body it belongs to had a hole exactly one triangle wide.
///
/// Area is the assertion because area is what went missing: 0.842 mm² with
/// the corner cut, 1.104 with it kept.
#[test]
fn a_cone_s_apex_run_is_kept_at_a_quarter_turn() {
    let text = corpus("cone_apex_quarter_turn.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let area: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude() * 0.5
        })
        .sum();
    assert!(
        (area - 1.10391).abs() < 0.01,
        "the apex triangle is drawn: area {area:.5}, not 0.84195"
    );
    // Both apex ends survive into the mesh: two vertices at the apex row.
    let (_, (va, _)) = ogeom::geom::Surface::domain(
        model
            .geometry()
            .surface(
                model
                    .node(&faces[0])
                    .unwrap()
                    .data()
                    .as_face()
                    .unwrap()
                    .surface,
            )
            .unwrap(),
    );
    let at_apex = mesh
        .parameters
        .iter()
        .filter(|(_, v)| (v - va).abs() < 1e-9 || (v + 1.0).abs() < 1e-9)
        .count();
    assert!(
        at_apex >= 2,
        "both rulings reach the apex row: {at_apex} vertices there"
    );
}

/// A ring folds across a closed chart's join, periodic or not.
///
/// The face lies on a B-spline tube that closes on itself in `u` — the
/// same points at `u = 0` and `u = 1` — without being periodic, and its
/// trim crosses that join twice. Walking the ring, the fold onto the branch
/// that continues it engaged on periodicity alone, so on this surface it
/// never fired: consecutive edges' images stood a whole chart apart, the
/// ring jumped the width of the chart twice, and the triangulator drew six
/// pieces. Closure, not periodicity, is the test — the same distinction the
/// projected-fit unwrap learned — and folded, the face is one piece.
#[test]
fn a_ring_folds_across_a_closed_chart_s_join() {
    use std::collections::HashMap;
    let text = corpus("closed_tube_face_crosses_its_join.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let euler = mesh.positions.len() as i64 - uses.len() as i64 + mesh.triangles.len() as i64;
    assert_eq!(
        euler,
        1,
        "one piece, not six: {} verts, {} tris",
        mesh.positions.len(),
        mesh.triangles.len()
    );
    let area: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude() * 0.5
        })
        .sum();
    assert!(
        (area - 106.5).abs() < 1.0,
        "the face's own area, not the pieces' overlap: {area:.2}"
    );
}

/// An inner loop thinner than a micron is a slit, not a hole.
///
/// The face is a plane with twenty inner loops. Fifteen are holes. Five run
/// out along two arcs and back along two splines fitted to the same arcs —
/// three millimetres long, a fifth of a micron wide, enclosing nothing —
/// and read as holes they are a tangle the triangulator cannot classify:
/// it drew the face with thirty-two holes. Measured in space, they are
/// slits, and dropped; the fifteen real holes stay.
///
/// `V - E + F` is `1 - holes` for a face with inner loops: fifteen holes
/// give −14. Thirty-two gave −31.
#[test]
fn a_slit_loop_is_not_a_hole() {
    use std::collections::HashMap;
    let text = corpus("slit_loops_are_not_holes.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    assert_eq!(
        model.ordered_children_of(&faces[0]).unwrap().len(),
        21,
        "twenty inner loops as read"
    );
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let euler = mesh.positions.len() as i64 - uses.len() as i64 + mesh.triangles.len() as i64;
    assert_eq!(
        euler,
        1 - 15,
        "fifteen holes, not thirty-two: {} verts, {} tris",
        mesh.positions.len(),
        mesh.triangles.len()
    );
}

/// A long bore is drawn round between its cross holes, not square.
///
/// The face is a cylinder 2.1 mm in radius and four hundred long, crossed
/// by holes wider than itself. It never sags along its axis, so sag gave
/// the grid one interior row, and the Delaunay triangulation bridged two
/// hundred millimetres from each rim to that row with triangles a quarter
/// turn wide — each sagging less than the three chords the repair pass
/// fires at, so they stayed. Grid cells are held to a bounded aspect now:
/// rows close enough that no triangle can reach across more than a few
/// columns.
///
/// Every ten-millimetre band of the bore holds vertices at fifteen
/// distinct whole degrees; the square bore held three.
#[test]
fn a_long_bore_is_round_between_its_holes() {
    use std::collections::{BTreeMap, BTreeSet};
    let text = corpus("long_bore_between_cross_holes.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    // The bore's axis is `z`; a band's angles are the whole degrees its
    // vertices sit at about that axis.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "whole degrees and centimetres"
    )]
    let whole = |x: f64| x.round() as i64;
    let mut bands: BTreeMap<i64, BTreeSet<i64>> = BTreeMap::new();
    for p in &mesh.positions {
        let degrees = whole(p.y.atan2(p.x).to_degrees());
        bands
            .entry(whole((p.z / 10.0).floor()))
            .or_default()
            .insert(degrees);
    }
    let (band, angles) = bands
        .iter()
        .map(|(z, a)| (*z, a.len()))
        .min_by_key(|b| b.1)
        .unwrap();
    assert!(
        angles >= 12,
        "the band at z {}..{} mm holds {angles} distinct angles: a square, not a bore",
        band * 10,
        band * 10 + 10
    );
}

/// A bore's inside is held to the angular deflection its rims are.
///
/// A viewer scaling its chord to a body's size hands a long extrusion a
/// chord of a third of a millimetre, and holds the tangent's turn to half a
/// radian. Its edges are drawn to both, so the bore's rims come out
/// thirteen-sided; the interior grid, held to the chord alone, came out
/// seven-sided, and the bore changed shape a chord in from each rim. The
/// grid's cells are now held to the normal's turn as well.
///
/// The same bore as above, at that viewer's deflection: every ten
/// millimetres holds vertices at twelve or more whole degrees.
#[test]
fn a_bore_s_inside_is_as_round_as_its_rims() {
    use std::collections::{BTreeMap, BTreeSet};
    let text = corpus("long_bore_between_cross_holes.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    let deflection = ogeom::mesh::Deflection {
        chord: 0.31,
        angular: 0.5,
        ..ogeom::mesh::Deflection::default()
    };
    let mesh = ogeom::mesh::triangulate_face(model, &faces[0], deflection, T).unwrap();
    #[allow(
        clippy::cast_possible_truncation,
        reason = "whole degrees and centimetres"
    )]
    let whole = |x: f64| x.round() as i64;
    let mut bands: BTreeMap<i64, BTreeSet<i64>> = BTreeMap::new();
    for p in &mesh.positions {
        bands
            .entry(whole((p.z / 10.0).floor()))
            .or_default()
            .insert(whole(p.y.atan2(p.x).to_degrees()));
    }
    let fewest = bands.values().map(BTreeSet::len).min().unwrap();
    assert!(
        fewest >= 12,
        "a band holds only {fewest} distinct angles: the inside is coarser than the rims"
    );
}

/// A chart's size is how far it reaches, not where it sits.
///
/// The face is a cylinder twenty millimetres tall whose axis point the
/// file placed half a metre away, so its chart spans `v` from −500 000 to
/// −499 980. The scale a degenerate triangle was measured against was the
/// difference between the smallest and largest coordinate over both axes
/// — half a million — and at that scale a quarter of a chart unit was a
/// hair: every cell of the grid was dropped and the face drew as two
/// triangles. The scale is the region's span now.
#[test]
fn a_chart_far_from_its_origin_is_not_degenerate() {
    let text = corpus("chart_far_from_its_origin.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let area: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude() * 0.5
        })
        .sum();
    // 2π · 2.25 · 20, less the chord's inscribed deficit.
    assert!(
        (area - 282.7).abs() < 2.0,
        "the whole face, not two triangles of it: area {area:.2} over {} triangles",
        mesh.triangles.len()
    );
}

/// Half a radian of angular deflection: what a viewer that keeps a circle
/// at thirteen segments asks for, and coarser than the default's eleven
/// degrees.
fn half_a_radian() -> ogeom::mesh::Deflection {
    ogeom::mesh::Deflection {
        angular: 0.5,
        ..ogeom::mesh::Deflection::default()
    }
}

/// `V - E + F` of a face's mesh, counting each undirected edge once.
fn euler_of(mesh: &ogeom::topo::Triangulation) -> i64 {
    use std::collections::HashMap;
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    mesh.positions.len() as i64 - uses.len() as i64 + mesh.triangles.len() as i64
}

/// An annulus narrower than its rims' sag is drawn finer, not refused.
///
/// Forty microns wide between rims of 2.845 and 2.805 mm: at half a radian
/// each rim is a sixteen-gon sagging fifty-five microns, the two polygons
/// cross, and the first pass encloses nothing. A first pass that came back
/// in *fragments* was already drawn again with finer edges; one that came
/// back *empty* was refused before it could be. Empty is short too.
#[test]
fn an_annulus_narrower_than_its_rims_sag_is_drawn_finer() {
    let text = corpus("annulus_narrower_than_its_rims_sag.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh = ogeom::mesh::triangulate_face(model, &faces[0], half_a_radian(), T)
        .expect("drawn finer, not refused");
    assert_eq!(euler_of(&mesh), 0, "one hole: an annulus");
    let area: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude() * 0.5
        })
        .sum();
    // π · (2.8448² − 2.8054²), less what the inscribed polygons leave out.
    assert!(
        (area - 0.699).abs() < 0.03,
        "the annulus's own area: {area:.4} over {} triangles",
        mesh.triangles.len()
    );
}

/// A repair point that lands on a vertex already there is not inserted.
///
/// A turned part with a three-edged B-spline patch whose bottom row
/// collapses to a point. Meshed whole at half a radian, three grid points
/// on that patch sit on a diagonal, the middle one a rounding off the
/// line, and the sliver they make has its centre at that middle point to
/// the last bits. The sag repair inserted the centre, round after round,
/// each a hair on the last; the degenerate filter dropped the hairs and
/// left a hole, and the solid was open by six edges.
#[test]
fn a_sliver_on_a_diagonal_of_the_grid_leaves_the_solid_closed() {
    let text = corpus("sliver_on_a_diagonal_of_the_grid.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let mesh = ogeom::mesh::triangulate(model, &import.solids[0], half_a_radian(), T).unwrap();
    assert!(
        mesh.is_closed(),
        "a hole where the hairs were: {} triangles",
        mesh.triangles.len()
    );
}

/// A grid point on a boundary segment is not inserted.
///
/// A turned part with a B-spline patch and the torus across one of its
/// edges, an edge that runs diagonally across the patch's chart. Meshed
/// whole at half a radian, a grid point falls exactly on that segment —
/// the midpoint of two grid corners the ring joins — and even-odd counting
/// calls it inside; inserted, it split the constraint on the patch alone,
/// and the torus was drawn to the unsplit edge: a T-junction, and the
/// solid open by six edges.
#[test]
fn a_grid_point_on_a_diagonal_boundary_leaves_the_solid_closed() {
    let text = corpus("grid_point_on_a_diagonal_boundary.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let mesh = ogeom::mesh::triangulate(model, &import.solids[0], half_a_radian(), T).unwrap();
    assert!(
        mesh.is_closed(),
        "a T-junction on the patch's edge: {} triangles",
        mesh.triangles.len()
    );
}
