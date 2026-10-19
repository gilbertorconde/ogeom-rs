//! The corpus imported, healed, measured, and operated on: real files in,
//! real modelling out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom::core::Tolerances;
use ogeom::topo::{Filter, ShapeType, explore, explore_unique};
use walks::edges_walked_one_way;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// Healing sweeps the whole corpus: every part reads, every shell closes,
/// every part measures, and each volume pins to its known figure: the
/// values ftc_07 and ftc_11 arbitrated against the kernel's exact ray
/// classifier and ctc_01 against an orientation-free even-odd grid, and
/// ftc_06 with its two balls standing on the plate, as their own loops and
/// the concave fillets at their feet say, none of its faces inside out. A loose
/// relative band absorbs mesh refinements. An orientation or unit
/// mistake moves a volume by whole factors and cannot hide inside it.
#[test]
#[ignore = "heavy"]
fn the_corpus_heals_and_every_part_measures() {
    let files: [(&str, f64); 11] = [
        ("nist_ctc_01_asme1_rd.stp", 14_643_073.4),
        ("nist_ctc_02_asme1_rc.stp", 47_101_909.1),
        ("nist_ctc_03_asme1_rc.stp", 331_884.8),
        ("nist_ctc_04_asme1_rd.stp", 17_519_472.2),
        ("nist_ctc_05_asme1_rd.stp", 12_694_721.1),
        ("nist_ftc_06_asme1_rd.stp", 3_323_546.1),
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
        // something is. Either way the shell must close.
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

/// A boolean over imported geometry, history checked: the operations run
/// on the world's parts, not only on this kernel's own.
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

    // A square post cut down through the plate's solid ring. The part has a
    // large central pocket, and a post through fresh air would cut nothing.
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
    // and it *measures*: the mesh welds across the file's slop because the
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
/// An imported plane's own window spans kilometres and a cylinder's height
/// domain is unbounded, so a bound read from the carriers is a box a
/// billion millimetres across. The hull of the topological vertices does
/// not contain the body either: a button head's apex is a bulge between
/// its rims, three millimetres past every vertex the head has.
///
/// An edge bounds the range it uses, not its whole curve, so a segment on a
/// line that runs to the ends of the world does not report the ends of the
/// world. A curved face is bounded by the edges below it, plus where its
/// surface bulges past that boundary, not by its whole surface.
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
            // is the mesh that is short (a chord across an arc falls inside it),
            // so the margin is the chord's, either way round.
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
/// microns wide, between two nearly concentric circles. At the default
/// tenth of a millimetre the sagitta of each bounding arc is twenty-nine
/// microns (wider than the region itself), so at the caller's chord the
/// inner polyline crosses the outer one and what reaches the triangulator
/// is not a region: it comes back as a scatter of disconnected pieces, and
/// the holes between them keep the body from meshing closed.
///
/// `V - E + F` is the test because it is the question: one for a disc,
/// `1 - holes` for a face with inner loops, and anything else means the
/// pieces are not joined.
#[test]
fn a_sliver_face_is_drawn_fine_enough_to_triangulate_whole() {
    use ogeom::core::FastMap;
    let text = corpus("sliver_face_falls_apart.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");

    let wires = model.ordered_children_of(&faces[0]).unwrap().len();
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
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
/// two apex ends are distinct points a quarter period apart. In space they
/// are one vertex. The ring has to keep both (the run between them along
/// the degenerate row is the face's own boundary): reading "apart" as
/// "more than a quarter period" drops the second, cuts the corner through
/// the face, and loses the triangle at the apex. The face is still a disc,
/// a third smaller than it should be, and the body it belongs to has a
/// hole exactly one triangle wide.
///
/// Area is the assertion because area is what goes missing: 0.842 mm² with
/// the corner cut, 1.104 with it kept.
#[test]
fn a_cone_s_apex_run_is_kept_at_a_quarter_turn() {
    let text = corpus("cone_apex_quarter_turn.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    // The angular deflection the areas below are measured at: the
    // inscribed area moves with the deflection, so the default is not
    // relied on.
    let deflection = ogeom::mesh::Deflection {
        angular: 0.2,
        ..ogeom::mesh::Deflection::default()
    };
    let mesh = ogeom::mesh::triangulate_face(model, &faces[0], deflection, T).unwrap();
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
/// The face lies on a B-spline tube that closes on itself in `u` (the
/// same points at `u = 0` and `u = 1`) without being periodic, and its
/// trim crosses that join twice. Walking the ring, the fold onto the branch
/// that continues it engages on closure, not on periodicity. On periodicity
/// alone it never fires here: consecutive edges' images stand a whole chart
/// apart, the ring jumps the width of the chart twice, and the triangulator
/// draws six pieces. Folded, the face is one piece.
#[test]
fn a_ring_folds_across_a_closed_chart_s_join() {
    use ogeom::core::FastMap;
    let text = corpus("closed_tube_face_crosses_its_join.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
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

/// An inner loop walked out and back through the same vertices, thinner
/// than a micron, is a slit, not a hole.
///
/// The face is a plane with twenty inner loops. Fifteen are holes. Five run
/// out along two arcs and back along two splines fitted to the same arcs (
/// three millimetres long, about six tenths of a micron wide on average,
/// enclosing nothing).
/// Read as holes they are a tangle the triangulator cannot classify, and
/// the face draws with thirty-two holes. Each edge out is answered by one
/// back between the same two vertices, and measured in space they are
/// under a micron: slits, and dropped. The fifteen real holes stay.
///
/// `V - E + F` is `1 - holes` for a face with inner loops: fifteen holes
/// give −14.
#[test]
fn a_slit_loop_is_not_a_hole() {
    use ogeom::core::FastMap;
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
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
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
/// by holes wider than itself. It never sags along its axis, so sag alone
/// gives the grid one interior row, and a Delaunay triangulation bridges
/// two hundred millimetres from each rim to that row with triangles a
/// quarter turn wide, each sagging less than the three chords the repair
/// pass fires at. Grid cells are held to a bounded aspect: rows close
/// enough that no triangle can reach across more than a few columns.
///
/// Every ten-millimetre band of the bore holds vertices at fifteen
/// distinct whole degrees, where a square bore holds three.
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
/// thirteen-sided. An interior grid held to the chord alone comes out
/// seven-sided, and the bore changes shape a chord in from each rim, so
/// the grid's cells are held to the normal's turn as well.
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
/// −499 980. Measured against the difference between the smallest and
/// largest coordinate over both axes (half a million), a quarter of a
/// chart unit is a hair: every cell of the grid is dropped and the face
/// draws as two triangles. The scale a degenerate triangle is measured
/// against is the region's span.
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

/// Half a radian of angular deflection, spelled out: the fixtures below
/// are measured at this value, whatever the default is.
fn half_a_radian() -> ogeom::mesh::Deflection {
    ogeom::mesh::Deflection {
        angular: 0.5,
        ..ogeom::mesh::Deflection::default()
    }
}

/// `V - E + F` of a face's mesh, counting each undirected edge once.
fn euler_of(mesh: &ogeom::topo::Triangulation) -> i64 {
    use ogeom::core::FastMap;
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
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
/// cross, and the first pass encloses nothing. A first pass that comes
/// back in *fragments* is drawn again with finer edges, and so is one that
/// comes back *empty*: empty is short too.
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
/// the last bits. Inserting the centre puts a hair on the last point,
/// round after round. The degenerate filter drops the hairs and leaves a
/// hole, and the solid is open by six edges.
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

/// Two halves of a cone on one surface share their seam ruling, whose
/// chart image stands at u = 0; the half that runs from 2pi back to pi
/// needs it a period over. Walked into the apex, the ruling's two columns
/// are half a turn from the last point either way, and left where it was
/// it wound the ring onto the other half: both halves meshed over one
/// side, and the part tessellated open.
#[test]
fn halves_of_a_cone_sharing_their_seam_ruling_mesh_closed() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let mesh = ogeom::mesh::triangulate(
        model,
        &import.solids[0],
        ogeom::mesh::Deflection::default(),
        T,
    )
    .unwrap();
    assert!(mesh.is_closed());
}

/// A grid point on a boundary segment is not inserted.
///
/// A turned part with a B-spline patch and the torus across one of its
/// edges, an edge that runs diagonally across the patch's chart. Meshed
/// whole at half a radian, a grid point falls exactly on that segment (
/// the midpoint of two grid corners the ring joins), and even-odd counting
/// calls it inside. Inserted, it splits the constraint on the patch alone
/// while the torus is drawn to the unsplit edge: a T-junction, and the
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

/// A crossing the boundary alone can tell of.
///
/// A cylinder wall with six slanted slots, each slot's two sides different
/// curves between the same two points. At half a radian the two sides'
/// polylines cross, and a face drawn from them has two holes more than it
/// has: `V - E + F` of −7 for a face with six holes, which is −5. A count
/// of triangles against boundary points sees nothing: the crossing costs a
/// handful of triangles and the face's hundreds of interior points bury
/// the difference. The boundary is triangulated on its own first, where
/// the count is exact, and the slots' edges are drawn finer.
#[test]
fn a_crossing_is_told_by_the_boundary_alone() {
    let text = corpus("slots_cross_at_half_a_radian.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    assert_eq!(
        model.ordered_children_of(&faces[0]).unwrap().len(),
        7,
        "six slots"
    );
    let mesh = ogeom::mesh::triangulate_face(model, &faces[0], half_a_radian(), T).unwrap();
    assert_eq!(
        euler_of(&mesh),
        1 - 6,
        "six holes, not eight: {} verts, {} tris",
        mesh.positions.len(),
        mesh.triangles.len()
    );
}

/// A spike on the boundary bounds nothing and is not drawn.
///
/// A plane whose wire runs out along an edge to a point and straight back
/// over the same curve. Kept, the spike triangulates to two hairs (one
/// vertex and one triangle more than a boundary that encloses a region
/// has), and the face is two pieces by count. Stripped before the
/// triangulation, the face is one.
#[test]
fn a_spike_on_the_boundary_is_stripped() {
    let text = corpus("a_spike_on_the_boundary.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let mesh =
        ogeom::mesh::triangulate_face(model, &faces[0], ogeom::mesh::Deflection::default(), T)
            .unwrap();
    assert_eq!(
        euler_of(&mesh),
        1,
        "one piece: {} verts, {} tris",
        mesh.positions.len(),
        mesh.triangles.len()
    );
    assert!(
        mesh.positions.len() <= 12,
        "the spike's tip is not a vertex: {} verts",
        mesh.positions.len()
    );
}

/// A face is drawn on its surface, not as a quilt of fins off it.
///
/// A fillet strip a centimetre long and a couple of millimetres wide,
/// turning a quarter turn across its width, whose chart runs `u` over a
/// fiftieth of a unit and `v` over one. Delaunay in the raw chart joins
/// points along the strip across columns rather than to the row beside
/// them, and the triangles it makes are slivers in space that lift
/// folded: flat across a bend the surface takes in between, their normals
/// pointing where none of their vertices' do, so half the face shades as
/// creases. Triangulated in the chart scaled to the surface's own metric,
/// every triangle lies with its vertices' normals.
#[test]
fn a_narrow_chart_s_triangles_lie_on_the_surface() {
    let text = corpus("fillet_strip_with_a_narrow_chart.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let deflection = ogeom::mesh::Deflection {
        chord: 0.31,
        angular: 0.5,
        ..ogeom::mesh::Deflection::default()
    };
    let mesh = ogeom::mesh::triangulate_face(model, &faces[0], deflection, T).unwrap();
    // Each triangle against the mean of its vertices' normals. A fin
    // stands near ninety degrees off, a facet on a fillet at most a few
    // tens.
    let total: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude()
        })
        .sum();
    let mut worst = 0.0_f64;
    let mut fins = 0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let own = (b - a).cross(c - a);
        if own.magnitude() < total * 1e-6 {
            continue;
        }
        let theirs =
            mesh.normals[t[0] as usize] + mesh.normals[t[1] as usize] + mesh.normals[t[2] as usize];
        let angle = (own.dot(theirs) / (own.magnitude() * theirs.magnitude()))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees();
        worst = worst.max(angle);
        if angle > 60.0 {
            fins += 1;
        }
    }
    assert_eq!(
        fins,
        0,
        "{fins} of {} triangles stand off the surface; the worst by {worst:.1} degrees",
        mesh.triangles.len()
    );
    assert!(worst < 45.0, "the worst triangle leans {worst:.1} degrees");
}

/// A face narrower than the chord is drawn to a fraction of its width.
///
/// A thread flank a tenth of a millimetre wide, drawn at a chord three
/// times that: the boundary sags by more than the face is wide between
/// its points, and every triangle across the width (a boundary chord
/// two millimetres long against a point a hundredth of a millimetre in)
/// stands off the surface by that sag. The face's edges are drawn to a
/// quarter of its width instead, and the triangles lie on it.
#[test]
fn a_face_narrower_than_the_chord_draws_its_edges_finer() {
    let text = corpus("thread_flank_narrower_than_a_chord.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 1, "the fixture is the one face");
    let deflection = ogeom::mesh::Deflection {
        chord: 0.31,
        angular: 0.5,
        ..ogeom::mesh::Deflection::default()
    };
    let mesh = ogeom::mesh::triangulate_face(model, &faces[0], deflection, T).unwrap();
    let (mut total, mut fins) = (0.0_f64, 0.0_f64);
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let own = (b - a).cross(c - a);
        let area = own.magnitude() * 0.5;
        total += area;
        let theirs =
            mesh.normals[t[0] as usize] + mesh.normals[t[1] as usize] + mesh.normals[t[2] as usize];
        if own.magnitude() > 0.0 && theirs.magnitude() > 0.0 {
            let angle = (own.dot(theirs) / (own.magnitude() * theirs.magnitude()))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            if angle > 45.0 {
                fins += area;
            }
        }
    }
    assert!(
        fins < total * 0.01,
        "{fins:.2} of {total:.2} mm² stands off the surface over {} triangles",
        mesh.triangles.len()
    );
}

/// A plane cutting a drum almost along its axis leaves an ellipse six and
/// a half metres by 1.8 millimetres, and the eight faces round one corner
/// of the fixture are bounded by short arcs of such ellipses. Their
/// vertices sit a couple of microns off the curve, where the closed-form
/// inversion can read them nine millimetres along it, run the edges the
/// long way round the ellipse, and draw three faces out to 10.8 m from a
/// 1 × 3 mm corner. Every face's mesh stays within its vertices' span,
/// grown by a millimetre for an arc's own bulge.
#[test]
fn faces_bounded_by_arcs_of_an_eccentric_ellipse_mesh_on_themselves() {
    let text = corpus("eccentric_ellipse_edges.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let faces = explore_unique(model, &import.solids[0], ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 8);
    for (index, face) in faces.iter().enumerate() {
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for vertex in explore_unique(model, face, ShapeType::Vertex).unwrap() {
            let data = model.node(&vertex).unwrap().data().as_vertex().unwrap();
            let p = vertex.transform(model.datums()).unwrap().apply(data.point);
            for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[k] = lo[k].min(c);
                hi[k] = hi[k].max(c);
            }
        }
        let mesh =
            ogeom::mesh::triangulate_face(model, face, ogeom::mesh::Deflection::default(), T)
                .unwrap();
        assert!(!mesh.triangles.is_empty(), "face {index} draws");
        for p in &mesh.positions {
            for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                assert!(
                    c >= lo[k] - 1.0 && c <= hi[k] + 1.0,
                    "face {index} draws a point off itself: {p:?}"
                );
            }
        }
    }
}

/// Re-anchoring a solid's periodic rings rebuilds every face that uses a
/// moved ring, and a rebuilt face keeps the side its flag named: a plate's
/// underside, holed by rings the heal moves, still faces down afterwards,
/// so a boolean or a volume read from it counts the plate as material.
#[test]
fn a_face_rebuilt_by_ring_reanchoring_keeps_its_side() {
    use ogeom::algo::{Containment, classify_on_face};
    use ogeom::math::Point;
    let text = corpus("nist_ctc_01_asme1_rd.stp");
    let mut import = ogeom::io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    // A point on the plate's underside, material above it.
    let at = Point::new(73.3, -91.16, -50.0);
    let facing = |model: &ogeom::topo::Model, shape: &ogeom::topo::Shape| -> f64 {
        let face = explore(model, shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .into_iter()
            .find(|f| {
                classify_on_face(model, f, at, ogeom::mesh::Deflection::default(), T).unwrap()
                    == Containment::In
            })
            .expect("the underside holds the point");
        let (_, normal) = ogeom::algo::face_normal(model, &face, T).unwrap();
        normal.z
    };
    assert!(
        facing(import.document.model(), &solid) < 0.0,
        "imported facing down"
    );
    let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
        .unwrap()
        .0
        .shape;
    assert!(
        facing(import.document.model(), &healed) < 0.0,
        "the rebuilt underside faces down still"
    );
}

/// Faces turned against their surfaces read with their loops as the file
/// walks them about the face's own normal: every edge of every part read
/// is walked once each way.
#[test]
fn every_edge_reads_walked_once_each_way() {
    for name in [
        "nist_ctc_01_asme1_rd.stp",
        "nist_ftc_09_asme1_rd.stp",
        "m5x16_bhcs.step",
        "sliver_on_a_diagonal_of_the_grid.step",
    ] {
        let import = ogeom::io::read_step(&corpus(name), T).unwrap();
        assert!(
            !import
                .report
                .warnings
                .iter()
                .any(|w| w.contains("walked back")),
            "{name}"
        );
        for solid in &import.solids {
            assert_eq!(
                edges_walked_one_way(import.document.model(), solid),
                0,
                "{name}"
            );
        }
    }
}

/// Every reversed face's bounds in a file turned the other way: the loops
/// of each face turned against its surface walked backward.
fn bounds_walked_back(text: &str) -> String {
    let reference = |s: &str| s.trim().trim_start_matches('#').parse::<u64>().ok();
    let mut bounds = ogeom::core::FastSet::default();
    for line in text.lines() {
        if line.contains("=ADVANCED_FACE(") && line.ends_with(",.F.);") {
            let list = line.split_once("',(").unwrap().1.split_once(')').unwrap().0;
            bounds.extend(list.split(',').filter_map(reference));
        }
    }
    text.lines()
        .map(|line| {
            let id = line.split_once('=').and_then(|(id, _)| reference(id));
            if id.is_some_and(|id| bounds.contains(&id)) {
                if line.ends_with(",.T.);") {
                    return line.replace(",.T.);", ",.F.);");
                }
                return line.replace(",.F.);", ",.T.);");
            }
            line.to_string()
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// A file whose faces turned against their surfaces walk their loops
/// backward reads as if it walked them right: each such face's loops are
/// walked back, every edge is walked once each way, and the part has the
/// faces and volume it has read from the file written right.
#[test]
fn faces_walking_their_loops_backward_are_walked_back() {
    let coarse = ogeom::mesh::Deflection::default();
    let right = corpus("m5x16_bhcs.step");
    let wrong = bounds_walked_back(&right);
    assert_ne!(right, wrong);
    let read = |text: &str| {
        let import = ogeom::io::read_step(text, T).unwrap();
        let model = import.document.model();
        let measured: Vec<(usize, f64)> = import
            .solids
            .iter()
            .map(|solid| {
                assert_eq!(edges_walked_one_way(model, solid), 0);
                (
                    explore_unique(model, solid, ShapeType::Face).unwrap().len(),
                    ogeom::algo::volume_properties(model, solid, coarse, T)
                        .unwrap()
                        .mass,
                )
            })
            .collect();
        let walked_back = import
            .report
            .warnings
            .iter()
            .any(|w| w.contains("walked back"));
        (measured, walked_back)
    };
    let (want, _) = read(&right);
    let (got, walked_back) = read(&wrong);
    assert!(walked_back, "the report says the loops were walked back");
    assert_eq!(want.len(), got.len());
    for ((faces, volume), (got_faces, got_volume)) in want.iter().zip(&got) {
        assert_eq!(faces, got_faces);
        assert!(
            (volume - got_volume).abs() <= volume * 1e-9,
            "{volume} read right, {got_volume} walked back"
        );
    }
}

/// A part as read, its faces listed from three different first faces: the
/// moments are taken from a point of the first, and a boundary left open
/// would move the volume by the open flux times the step. Its fitted
/// edges' pcurves meet their curves at every inner sample and stray up to
/// 1.5e-4 off them towards their ends; the strips there close it, and the
/// volume agrees to a part in a hundred billion.
#[test]
#[ignore = "heavy"]
fn a_read_part_s_volume_does_not_move_with_its_first_face() {
    use ogeom::algo::{make_shell, make_solid, shape_bounds, volume_properties};
    let mut import = ogeom::io::read_step(&corpus("nist_ctc_02_asme1_rc.stp"), T).unwrap();
    let part = import.solids[0].clone();
    let model = import.document.model_mut();
    let diagonal = shape_bounds(model, &part, T).unwrap().diagonal();
    let fine = ogeom::mesh::Deflection::with_chord(diagonal * 1e-5).unwrap();
    let faces = explore_unique(model, &part, ShapeType::Face).unwrap();
    let mut volumes = Vec::new();
    for first in [0, faces.len() / 3, 2 * faces.len() / 3] {
        let mut listed = faces.clone();
        listed.rotate_left(first);
        let shell = make_shell(model, &listed).unwrap().shape;
        let solid = make_solid(model, &[shell]).unwrap().shape;
        let measured = volume_properties(model, &solid, fine, T).unwrap();
        assert_eq!(measured.deflection, 0.0, "integrated, not meshed");
        volumes.push(measured.mass);
    }
    for volume in &volumes[1..] {
        assert!(
            (volume - volumes[0]).abs() <= volumes[0] * 1e-11,
            "the volume moves with the first face: {volumes:?}"
        );
    }
}

/// Faces meeting at a vertex place it at the vertex's own point, bit for
/// bit, so meshes joined on exact coordinates close. One of this part's
/// faces has a corner whose chart parameter comes out subnormal-tiny; the
/// triangulation holds it at zero, and the corner is still the vertex's.
#[test]
fn faces_meeting_at_a_vertex_place_it_at_one_point() {
    use ogeom::core::FastMap;
    let import = ogeom::io::read_step(&corpus("stepconf_sled_base.step"), T).unwrap();
    let model = import.document.model();
    let solid = &import.solids[0];
    let bits = |p: ogeom::math::Point| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    for chord in [0.0125, 0.05, 0.1, 0.5] {
        let deflection = ogeom::mesh::Deflection {
            chord,
            ..ogeom::mesh::Deflection::default()
        };
        let (faces, _) = ogeom::mesh::face_meshes_for(model, solid, deflection, None, T).unwrap();
        let mut sides: FastMap<([u64; 3], [u64; 3]), usize> = FastMap::default();
        for (_, mesh) in &faces {
            let mesh = mesh.as_ref().unwrap();
            for t in &mesh.triangles {
                for k in 0..3 {
                    let p = bits(mesh.positions[t[k] as usize]);
                    let q = bits(mesh.positions[t[(k + 1) % 3] as usize]);
                    if p != q {
                        *sides.entry((p.min(q), p.max(q))).or_default() += 1;
                    }
                }
            }
        }
        let open = sides.values().filter(|&&n| n == 1).count();
        assert_eq!(open, 0, "chord {chord}: {open} sides with one triangle");
    }
}

/// A bowl 380 mm across, its rim a surface of revolution: at a chord of
/// 0.0033 mm the grid needs more than 512 steps round it. Each direction
/// may take as many as the deflection's `max_segments`, refined evenly
/// along its length. The rim stands within the chord of its surface, and
/// the bowl's doubly curved face, whose grid cells may sag both ways at
/// once, within half as much again.
#[test]
fn a_wide_revolved_rim_meets_a_fine_chord() {
    use ogeom::geom::{Surface as _, SurfaceKind, Transformable as _};
    let import = ogeom::io::read_step(&corpus("stepconf_bowl_01.step"), T).unwrap();
    let model = import.document.model();
    let part = import.solids.first().or(import.shells.first()).unwrap();
    let chord = 0.0033;
    let deflection = ogeom::mesh::Deflection {
        chord,
        ..ogeom::mesh::Deflection::default()
    };
    let (faces, _) = ogeom::mesh::face_meshes_for(model, part, deflection, None, T).unwrap();
    let mut rims = 0;
    for (face, mesh) in &faces {
        let data = model.node(face).unwrap().data().as_face().unwrap();
        let surface = model.geometry().surface(data.surface).unwrap().clone();
        if surface.kind() != SurfaceKind::Revolution {
            continue;
        }
        rims += 1;
        let surface = surface
            .transformed(&face.transform(model.datums()).unwrap(), T)
            .unwrap();
        let mesh = mesh.as_ref().unwrap();
        let mut worst = 0.0_f64;
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
            let middle = ogeom::math::Point::from_vector(
                (a.to_vector() + b.to_vector() + c.to_vector()) / 3.0,
            );
            let off = ogeom::algo::project_on_surface(&surface, middle, 8, T)
                .unwrap()
                .distance;
            worst = worst.max(off);
        }
        assert!(worst <= 1.5 * chord, "a middle {worst} off a revolved face");
    }
    assert!(rims > 0, "the bowl has a revolved rim");
}

/// The edge of `solid` passing nearest `at` while running along `along`
/// there, read off samples of each edge's curve.
fn edge_along(
    model: &ogeom::topo::Model,
    solid: &ogeom::topo::Shape,
    at: ogeom::math::Point,
    along: ogeom::math::Vector,
) -> ogeom::topo::Shape {
    use ogeom::geom::Curve3d as _;
    let mut best: Option<(f64, ogeom::topo::Shape)> = None;
    for edge in explore_unique(model, solid, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let placed = edge.transform(model.datums()).unwrap();
        let points: Vec<_> = (0..=256)
            .map(|i| {
                let t = (range.1 - range.0).mul_add(f64::from(i) / 256.0, range.0);
                placed.apply(curve.point_at(t, T).unwrap())
            })
            .collect();
        for pair in points.windows(2) {
            let step = pair[1] - pair[0];
            let length = step.magnitude();
            if length <= 0.0 || (step.dot(along) / length).abs() < 0.9 {
                continue;
            }
            let off = pair[0].distance(at).min(pair[1].distance(at));
            if best.as_ref().is_none_or(|(held, _)| off < *held) {
                best = Some((off, edge.clone()));
            }
        }
    }
    let (off, edge) = best.unwrap();
    assert!(off < 0.1, "no edge runs along {along:?} through {at:?}");
    edge
}

/// A solid's volume, its bounds, and whether it checks.
/// Each box of a half carriage fused with its mirror image across the
/// plane its swept walls, rounds and bore are cut in half by, then
/// unified: every cylinder and extrusion face the plane split is one face
/// again (the bore a whole tube closed by a seam), the image's halves
/// carried into the part's charts, and the solid is unchanged.
#[test]
fn curved_walls_split_by_the_mirror_plane_unify() {
    use ogeom::geom::{Surface as _, SurfaceKind};
    let (mut model, roots) =
        ogeom::io::native::read(&corpus("curved_walls_across_the_mirror_plane.ogeom")).unwrap();
    // How many faces lie on each kind of surface, and how many edges on
    // the mirror plane two faces on one kind of curved surface share.
    let census = |model: &ogeom::topo::Model, solid: &ogeom::topo::Shape| {
        let mut kinds = (0, 0);
        let mut users: ogeom::core::FastMap<_, Vec<SurfaceKind>> = ogeom::core::FastMap::default();
        for face in explore_unique(model, solid, ShapeType::Face).unwrap() {
            let Some(ogeom::topo::NodeData::Face(data)) = model.node(&face).map(|n| n.data())
            else {
                continue;
            };
            let kind = model.geometry().surface(data.surface).unwrap().kind();
            match kind {
                SurfaceKind::Cylinder => kinds.0 += 1,
                SurfaceKind::Extrusion => kinds.1 += 1,
                _ => {}
            }
            for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
                let b = ogeom::algo::tight_bounds(model, &edge, T).unwrap();
                if b.low().unwrap().x.abs() <= 1e-6 && b.high().unwrap().x.abs() <= 1e-6 {
                    users.entry(edge.node()).or_default().push(kind);
                }
            }
        }
        let split = users
            .values()
            .filter(|k| matches!(k.as_slice(), [a, b] if a == b && *a != SurfaceKind::Plane))
            .count();
        (kinds, split)
    };
    // Cylinders and extrusions fused, then unified: the sweeps' box holds
    // two rounds and two swept walls on the plane and a pair of walls
    // split in the part itself; the bore's box holds the bore.
    let expected = [((6, 10), (4, 7)), ((2, 0), (1, 0))];
    for (part, (fused_kinds, unified_kinds)) in roots.iter().zip(expected) {
        let (half, _) = measured(&model, part);
        let copy = ogeom::algo::copied(&mut model, part).unwrap().shape;
        let mirror = ogeom::math::Transform::plane_mirror(
            ogeom::math::Point::ORIGIN,
            ogeom::math::Direction::X,
        );
        let image = ogeom::algo::transformed(&mut model, &copy, mirror)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, part, &image, T)
            .unwrap()
            .shape;
        let (kinds, split) = census(&model, &fused);
        assert_eq!(kinds, fused_kinds);
        assert!(split > 0, "the fuse leaves faces split on the plane");

        let (unified, _) = ogeom::heal::unify_same_domain(&mut model, &fused, T).unwrap();
        let (whole, _) = measured(&model, &unified.shape);
        assert_eq!(census(&model, &unified.shape), (unified_kinds, 0));
        assert!(
            (whole - 2.0 * half).abs() <= half * 1e-8,
            "the whole measures {whole} against twice the half's {half}"
        );
    }
}

fn measured(model: &ogeom::topo::Model, solid: &ogeom::topo::Shape) -> (f64, ogeom::math::Aabb) {
    assert!(ogeom::algo::check(model, solid, T).unwrap().is_valid());
    let volume =
        ogeom::algo::volume_properties(model, solid, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    (volume, ogeom::algo::tight_bounds(model, solid, T).unwrap())
}

/// What square convex edges of these lengths lose to a round of `radius`:
/// the corner outside the ball's quarter circle, along each.
fn square_round_loss(model: &ogeom::topo::Model, edges: &[ogeom::topo::Shape], radius: f64) -> f64 {
    let length: f64 = edges
        .iter()
        .map(|e| {
            ogeom::algo::linear_properties(model, e, ogeom::mesh::Deflection::default(), T)
                .unwrap()
                .mass
        })
        .sum();
    (1.0 - std::f64::consts::FRAC_PI_4) * radius * radius * length
}

/// What two rounds of `radius` meeting across a square inside corner take
/// past their caps: each runs on to the plane bisecting the corner, and the
/// section a distance w in from its wall reaches w along, so each adds the
/// corner outside its quarter circle weighted by w, r^3 (5/6 - pi/4).
fn square_corner_cell(radius: f64) -> f64 {
    2.0 * (5.0 / 6.0 - std::f64::consts::FRAC_PI_4) * radius.powi(3)
}

/// A rib's top edges rounded where one of its walls is swept from a
/// fitted profile: a straight wall, the swept one, and a round end, the
/// three meeting at shallow corners. The swept wall ends where its
/// profile does, and the crease along its top is a fitted cubic that does
/// not continue the way the wall does; the ball's run-out past either end
/// stands on both. Each round loses the corner outside its quarter circle
/// along its edge and nothing else, whichever order the edges are named in,
/// and the rib keeps its bounds.
#[test]
fn a_rib_rounds_along_a_wall_swept_from_a_fitted_profile() {
    use ogeom::math::{Point, Vector};
    let (model, roots) =
        ogeom::io::native::read(&corpus("rib_round_on_a_swept_spline_wall.ogeom")).unwrap();
    let rib = &roots[0];
    let straight = edge_along(
        &model,
        rib,
        Point::new(-13.13, 9.51, -14.4),
        Vector::new(0.0, 0.0, 1.0),
    );
    let swept = edge_along(
        &model,
        rib,
        Point::new(-16.0, 9.51, -20.46),
        Vector::new(-0.8, 0.0, -0.6),
    );
    let end = edge_along(
        &model,
        rib,
        Point::new(-24.36, 9.51, -26.69),
        Vector::new(0.26, 0.0, -0.97),
    );
    let radius = 0.85;
    let (before, bounds) = measured(&model, rib);
    let loss = square_round_loss(
        &model,
        &[straight.clone(), swept.clone(), end.clone()],
        radius,
    );
    for named in [
        [straight.clone(), swept.clone(), end.clone()],
        [end.clone(), straight.clone(), swept.clone()],
    ] {
        let mut model = model.clone();
        let rounded = ogeom::fillet::fillet_edges(&mut model, rib, &named, radius, T).unwrap();
        let (after, held) = measured(&model, &rounded.shape);
        assert!(
            (before - after - loss).abs() <= loss * 0.02,
            "the rounds took {} where their sections over their edges are {loss}",
            before - after
        );
        assert!(
            held.diagonal() - bounds.diagonal() <= 1e-6,
            "the rib keeps its bounds"
        );
    }
    // Alone, the swept wall's round runs out past both its ends over its
    // neighbours' sharp corners: it takes its own section and some more.
    let mut alone = model.clone();
    let rounded =
        ogeom::fillet::fillet_edges(&mut alone, rib, std::slice::from_ref(&swept), radius, T)
            .unwrap();
    let (after, _) = measured(&alone, &rounded.shape);
    let own = square_round_loss(&model, std::slice::from_ref(&swept), radius);
    assert!(
        before - after >= own && before - after <= own * 1.5,
        "the swept wall's round took {} against its own section's {own}",
        before - after
    );
}

/// Two top edges rounded where they meet at an inside corner of the top
/// face, on a part whose second edge's curve runs a fraction of a micron
/// past the vertex it ends at (the vertex's tolerance says so). Each round
/// is capped at the vertex itself, where the other's round is tangent to
/// the cap's plane; capped where the curve ends, the cap cut a sliver off
/// the other round's end. Either order rounds both, each losing the corner
/// outside its quarter circle along its edge, and the two meeting across
/// the square inside corner in the plane that bisects it: each runs on past
/// its cap to that plane, which takes the cell past both caps that both
/// rounds hold, two times r^3 (5/6 - pi/4).
#[test]
fn rounds_meeting_at_a_vertex_their_edge_overruns_are_capped_there() {
    use ogeom::math::{Point, Vector};
    let (model, roots) =
        ogeom::io::native::read(&corpus("inner_corner_edge_past_its_vertex.ogeom")).unwrap();
    let part = &roots[0];
    let down = edge_along(
        &model,
        part,
        Point::new(21.6, -4.477, 0.0),
        Vector::new(0.0, 0.0, 1.0),
    );
    let across = edge_along(
        &model,
        part,
        Point::new(21.6, -5.7, -1.5),
        Vector::new(0.0, 1.0, 0.0),
    );
    let (before, bounds) = measured(&model, part);
    let loss =
        square_round_loss(&model, &[down.clone(), across.clone()], 1.0) + square_corner_cell(1.0);
    for named in [
        [down.clone(), across.clone()],
        [across.clone(), down.clone()],
    ] {
        let mut model = model.clone();
        let rounded = ogeom::fillet::fillet_edges(&mut model, part, &named, 1.0, T).unwrap();
        let (after, held) = measured(&model, &rounded.shape);
        assert!(
            (before - after - loss).abs() <= loss * 1e-3,
            "the rounds took {} where their sections over their edges are {loss}",
            before - after
        );
        assert!((held.diagonal() - bounds.diagonal()).abs() <= 1e-6);
    }
}

/// The same two top edges on a part where the inside corner's vertex
/// stands where the second edge's line starts, 0.24 µm off the wall the
/// first edge runs down, which the second edge's line runs into that far:
/// measured from the vertex, neither line runs past it. Each round is
/// capped where its line meets the other edge's, flush with the other
/// round; capped at the vertex, the second round cut a sliver off the
/// first's end, or off the wall's corner left after it, which the
/// boolean's arrangement could not place. Either order rounds both, each
/// losing the corner outside its quarter circle along its edge, and the
/// two meeting in the plane bisecting the square inside corner, as above.
#[test]
fn rounds_meeting_at_a_vertex_off_one_edges_wall_are_capped_at_the_wall() {
    use ogeom::math::{Point, Vector};
    let (model, roots) =
        ogeom::io::native::read(&corpus("inner_corner_vertex_off_the_wall.ogeom")).unwrap();
    let part = &roots[0];
    let down = edge_along(
        &model,
        part,
        Point::new(21.6, -4.477, 4.95),
        Vector::new(0.0, 0.0, 1.0),
    );
    let across = edge_along(
        &model,
        part,
        Point::new(21.6, -5.7, -1.5),
        Vector::new(0.0, 1.0, 0.0),
    );
    let (before, bounds) = measured(&model, part);
    let loss =
        square_round_loss(&model, &[down.clone(), across.clone()], 1.0) + square_corner_cell(1.0);
    for named in [
        [down.clone(), across.clone()],
        [across.clone(), down.clone()],
    ] {
        let mut model = model.clone();
        let rounded = ogeom::fillet::fillet_edges(&mut model, part, &named, 1.0, T).unwrap();
        let (after, held) = measured(&model, &rounded.shape);
        assert!(
            (before - after - loss).abs() <= loss * 1e-3,
            "the rounds took {} where their sections over their edges are {loss}",
            before - after
        );
        assert!((held.diagonal() - bounds.diagonal()).abs() <= 1e-6);
    }
}

/// A floor edge rounded where, at its end, the floor gives way to a ramp
/// rising at a little over 45 degrees and the wall the round rolls on
/// carries on above the ramp. The round's cylinder runs on past the end
/// along the wall until the ramp has cut its whole section away, a radius
/// over the ramp's slope, and stops on the ramp: no cap stands at the
/// edge's end. Past the end, the section above the ramp standing h over
/// the floor is the corner outside the quarter circle above h, and summed
/// along the ramp it is r^3 (5/6 - pi/4) over the slope, on top of the
/// corner outside the quarter circle along the edge.
#[test]
fn a_round_runs_out_along_its_wall_onto_a_ramp_past_its_end() {
    use ogeom::math::{Point, Vector};
    let (mut model, roots) =
        ogeom::io::native::read(&corpus("floor_edge_ending_at_a_ramp.ogeom")).unwrap();
    let part = &roots[0];
    let edge = edge_along(
        &model,
        part,
        Point::new(-4.2, 2.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    );
    let radius = 5.0;
    // The ramp's outward normal is (0, 0.717..., -0.697...).
    let slope = 0.717_058_867_531_502 / 0.697_012_611_431_414_5;
    let end = ogeom::algo::shape_bounds(&model, &edge, T)
        .unwrap()
        .high()
        .unwrap()
        .y;
    let (before, bounds) = measured(&model, part);
    let loss = square_round_loss(&model, std::slice::from_ref(&edge), radius)
        + radius.powi(3) * (5.0 / 6.0 - std::f64::consts::FRAC_PI_4) / slope;
    let rounded =
        ogeom::fillet::fillet_edges(&mut model, part, std::slice::from_ref(&edge), radius, T)
            .unwrap();
    let (after, held) = measured(&model, &rounded.shape);
    assert!(
        (before - after - loss).abs() <= loss * 1e-3,
        "the round took {} where its section along the edge and past it is {loss}",
        before - after
    );
    assert!((held.diagonal() - bounds.diagonal()).abs() <= 1e-6);
    let [band] = rounded.history.generated(&edge) else {
        panic!("the round is one face");
    };
    let reach = ogeom::algo::shape_bounds(&model, band, T)
        .unwrap()
        .high()
        .unwrap()
        .y;
    assert!(
        (reach - (end + radius / slope)).abs() <= 1e-3,
        "the round reaches y = {reach}"
    );
    for face in explore_unique(&model, &rounded.shape, ShapeType::Face).unwrap() {
        let b = ogeom::algo::shape_bounds(&model, &face, T).unwrap();
        let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
        assert!(
            !((lo.y - end).abs() <= 1e-3 && (hi.y - end).abs() <= 1e-3 && hi.x > -9.2),
            "a cap stands at the edge's end"
        );
    }
}

/// A crease of the rib's top plane `y = 9.51`: its length, how far it
/// turns towards the material along it, and its two ends with the unit
/// direction leaving each along the crease. Which side the material lies
/// is read off a point beside the crease's middle, under the top.
fn rib_crease(
    model: &ogeom::topo::Model,
    rib: &ogeom::topo::Shape,
    edge: &ogeom::topo::Shape,
) -> (f64, f64, [(ogeom::math::Point, ogeom::math::Vector); 2]) {
    use ogeom::geom::Curve3d as _;
    let data = model.node(edge).unwrap().data().as_edge().unwrap();
    let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        unreachable!("a rib crease has its curve");
    };
    let curve = model.geometry().curve(*curve).unwrap();
    let (lo, hi) = *range;
    let samples = 4096;
    let at = |i: u32| (hi - lo).mul_add(f64::from(i) / f64::from(samples), lo);
    let (mut length, mut turn) = (0.0, 0.0);
    let heading = |t: f64| {
        let d = curve.d1_at(t, T).unwrap();
        d.z.atan2(d.x)
    };
    for i in 0..samples {
        length += curve
            .point_at(at(i + 1), T)
            .unwrap()
            .distance(curve.point_at(at(i), T).unwrap());
        let mut step = heading(at(i + 1)) - heading(at(i));
        if step > std::f64::consts::PI {
            step -= std::f64::consts::TAU;
        } else if step < -std::f64::consts::PI {
            step += std::f64::consts::TAU;
        }
        turn += step;
    }
    let middle = f64::midpoint(lo, hi);
    let d = curve.d1_at(middle, T).unwrap();
    let left = ogeom::math::Vector::new(-d.z, 0.0, d.x) / d.magnitude();
    let beside =
        curve.point_at(middle, T).unwrap() + left * 0.05 - ogeom::math::Vector::new(0.0, 0.3, 0.0);
    let on_left = ogeom::algo::classify_in_solid_exact(model, rib, beside, T).unwrap()
        == ogeom::algo::Containment::In;
    let end = |t: f64, sign: f64| {
        let d = curve.d1_at(t, T).unwrap();
        (curve.point_at(t, T).unwrap(), d / d.magnitude() * sign)
    };
    (
        length,
        if on_left { turn } else { -turn },
        [end(lo, 1.0), end(hi, -1.0)],
    )
}

/// The turn between two creases at the end they share: the angle between
/// the one carried on past it and the other leaving it.
fn turn_between(
    a: &[(ogeom::math::Point, ogeom::math::Vector); 2],
    b: &[(ogeom::math::Point, ogeom::math::Vector); 2],
) -> f64 {
    for (p, leaving_a) in a {
        for (q, leaving_b) in b {
            if p.distance(*q) < 1e-6 {
                return (-*leaving_a).dot(*leaving_b).clamp(-1.0, 1.0).acos();
            }
        }
    }
    unreachable!("the creases share an end")
}

/// The rib's three top creases rounded together, meeting at two shallow
/// turns: the straight wall turns back behind the swept one between them,
/// and the round end's rim turns out from the swept wall. Each round takes
/// its section along its crease less its moment about the wall times its
/// crease's turn towards the material. Across the inside turn the two
/// rounds run on past their caps to the plane bisecting it, and take the
/// cell past both caps, twice the moment times tan(t/2); across the
/// outside turn they meet in the plane bisecting it, where their sections
/// agree to the fit, and leave the same twice the moment times tan(t/2) of
/// their capped wedges' overlap. The rib keeps its bounds.
#[test]
fn a_rib_s_creases_meeting_at_shallow_turns_take_their_sections_and_their_corners() {
    use ogeom::math::{Point, Vector};
    let (model, roots) =
        ogeom::io::native::read(&corpus("rib_round_on_a_swept_spline_wall.ogeom")).unwrap();
    let rib = &roots[0];
    let straight = edge_along(
        &model,
        rib,
        Point::new(-13.13, 9.51, -14.4),
        Vector::new(0.0, 0.0, 1.0),
    );
    let swept = edge_along(
        &model,
        rib,
        Point::new(-16.0, 9.51, -20.46),
        Vector::new(-0.8, 0.0, -0.6),
    );
    let end = edge_along(
        &model,
        rib,
        Point::new(-24.36, 9.51, -26.69),
        Vector::new(0.26, 0.0, -0.97),
    );
    let radius: f64 = 0.85;
    let section = (1.0 - std::f64::consts::FRAC_PI_4) * radius * radius;
    let moment = (5.0 / 6.0 - std::f64::consts::FRAC_PI_4) * radius.powi(3);
    let creases = [&straight, &swept, &end].map(|e| rib_crease(&model, rib, e));
    let bands: f64 = creases
        .iter()
        .map(|(length, turn, _)| section.mul_add(*length, -moment * turn))
        .sum();
    let inside = turn_between(&creases[0].2, &creases[1].2);
    let outside = turn_between(&creases[1].2, &creases[2].2);
    let want = 2.0f64.mul_add(
        moment * ((inside / 2.0).tan() - (outside / 2.0).tan()),
        bands,
    );
    let (before, bounds) = measured(&model, rib);
    for named in [
        [straight.clone(), swept.clone(), end.clone()],
        [end.clone(), straight.clone(), swept.clone()],
    ] {
        let mut model = model.clone();
        let rounded = ogeom::fillet::fillet_edges(&mut model, rib, &named, radius, T).unwrap();
        let (after, held) = measured(&model, &rounded.shape);
        assert!(
            (before - after - want).abs() <= want * 1e-4,
            "the rounds took {} where their sections and corners are {want}",
            before - after
        );
        assert!(held.diagonal() - bounds.diagonal() <= 1e-6);
    }
}

/// A rib's top edges rounded as one chain carried on along their tangent
/// continuations: two straight runs that continue up under a slanted
/// roof, and a round end whose rim is split at the wall's seam into a
/// long arc and a short one. The runs under the roof stop on it, cut off
/// by its plane: at a point of the section a distance b in from the rib's
/// top the roof stands h - k b up the run, k its slope, so each takes the
/// section (1 - pi/4) r^2 times h less k times the section's centroid's
/// distance from the top, r (10 - 3 pi) / (3 (4 - pi)). The rim stops
/// where the rib's walls leave it, so the short arc's round runs its own
/// arc and takes the section turned about the axis at the centroid's
/// radius over that arc, alone or with the rest. Each of the three adds
/// that to what the rest of the chain takes, whichever order the edges
/// are named in.
#[test]
fn a_rib_s_chain_rounds_on_along_its_tangent_continuations() {
    use ogeom::math::{Point, Vector};
    let (model, roots) =
        ogeom::io::native::read(&corpus("rib_chain_under_a_slanted_roof.ogeom")).unwrap();
    let rib = &roots[0];
    let along = |at: (f64, f64), dir: (f64, f64)| {
        edge_along(
            &model,
            rib,
            Point::new(at.0, 9.51, at.1),
            Vector::new(dir.0, 0.0, dir.1),
        )
    };
    let picked = [
        along((-9.21, -5.88), (0.0, 1.0)),
        along((-13.13, -7.84), (0.0, 1.0)),
        along((-4.605, -11.76), (1.0, 0.0)),
        along((-4.605, -15.68), (1.0, 0.0)),
        along((-16.004, -20.458), (-0.806, -0.592)),
        along((-24.472, -26.004), (-0.068, 0.998)),
    ];
    let short_arc = along((-20.191, -29.293), (-0.982, -0.191));
    let under_roof = [
        along((-13.13, 1.043), (0.0, 1.0)),
        along((-9.21, 1.043), (0.0, 1.0)),
    ];
    let radius = 0.85;
    let pi = std::f64::consts::PI;
    let section = (1.0 - pi / 4.0) * radius * radius;
    let centroid = radius * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
    // The roof's outward normal is (0, 0.717..., -0.697...).
    let slope: f64 = 0.717_058_867_531_502 / 0.697_012_611_431_414_5;
    let length = |edge: &ogeom::topo::Shape| {
        ogeom::algo::linear_properties(&model, edge, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let roofed: f64 = under_roof
        .iter()
        .map(|e| section * slope.mul_add(-centroid, length(e)))
        .sum();
    // The rim's radius is 3.6, and a circle's parameter is its angle.
    let turn = match model
        .node(&short_arc)
        .unwrap()
        .data()
        .as_edge()
        .unwrap()
        .curve3d()
    {
        Some(ogeom::topo::EdgeRepr::Curve3d { range, .. }) => range.1 - range.0,
        _ => unreachable!("the short arc has its curve"),
    };
    let arc = section * (3.6 - centroid) * turn;
    let (before, bounds) = measured(&model, rib);
    let took = |named: &[ogeom::topo::Shape]| {
        let mut model = model.clone();
        let rounded = ogeom::fillet::fillet_edges(&mut model, rib, named, radius, T).unwrap();
        let (after, held) = measured(&model, &rounded.shape);
        assert!((held.diagonal() - bounds.diagonal()).abs() <= 1e-6);
        before - after
    };
    let alone = took(std::slice::from_ref(&short_arc));
    assert!(
        (alone - arc).abs() <= arc * 1e-3,
        "the short arc's round took {alone} where its section over its arc is {arc}"
    );
    let rest = took(&picked);
    let mut chain: Vec<_> = picked.to_vec();
    chain.push(short_arc);
    chain.extend(under_roof);
    let added = roofed + arc;
    for order in [false, true] {
        if order {
            chain.reverse();
        }
        let whole = took(&chain);
        assert!(
            (whole - rest - added).abs() <= added * 1e-3,
            "the continuations took {} where their sections are {added}",
            whole - rest
        );
    }
}

/// The rib's lower top crease, along a wall swept from a fitted profile,
/// rounded alone and with the rest of the rib's chain. It leaves the round
/// end's rim with the wall turned a degree and a half back behind the
/// rim's, and ends where the outline turns back into the material along a
/// straight top edge. Alone, its round stops at both ends of the crease,
/// since past either the crease carries on into the material: the corner's
/// section along it, less its moment about the wall times the crease's turn
/// towards the material. With the chain it meets the rim's round and the
/// straight edge's across those two re-entrant turns, each adding the cell
/// past both caps, twice the moment times tan(t/2), t the turn.
#[test]
fn a_rib_s_lower_spline_crease_rounds_alone_and_with_its_chain() {
    use ogeom::math::{Point, Vector};
    let (model, roots) =
        ogeom::io::native::read(&corpus("rib_chain_under_a_slanted_roof.ogeom")).unwrap();
    let rib = &roots[0];
    let along = |at: (f64, f64), dir: (f64, f64)| {
        edge_along(
            &model,
            rib,
            Point::new(at.0, 9.51, at.1),
            Vector::new(dir.0, 0.0, dir.1),
        )
    };
    let straight_end = along((-4.605, -15.68), (1.0, 0.0));
    let short_arc = along((-20.191, -29.293), (-0.982, -0.191));
    let chain = [
        along((-9.21, -5.88), (0.0, 1.0)),
        along((-13.13, -7.84), (0.0, 1.0)),
        along((-4.605, -11.76), (1.0, 0.0)),
        straight_end.clone(),
        along((-16.004, -20.458), (-0.806, -0.592)),
        along((-24.472, -26.004), (-0.068, 0.998)),
        short_arc.clone(),
        along((-13.13, 1.043), (0.0, 1.0)),
        along((-9.21, 1.043), (0.0, 1.0)),
    ];
    let spline = along((-19.06, -28.92), (0.935, 0.355));
    let radius: f64 = 0.85;
    let section = (1.0 - std::f64::consts::FRAC_PI_4) * radius * radius;
    let moment = (5.0 / 6.0 - std::f64::consts::FRAC_PI_4) * radius.powi(3);
    let (length, turn, ends) = rib_crease(&model, rib, &spline);
    let band = section.mul_add(length, -moment * turn);
    let cells = [&short_arc, &straight_end]
        .map(|e| turn_between(&ends, &rib_crease(&model, rib, e).2))
        .iter()
        .map(|t| 2.0 * moment * (t / 2.0).tan())
        .sum::<f64>();
    let (before, bounds) = measured(&model, rib);
    let took = |named: &[ogeom::topo::Shape]| {
        let mut model = model.clone();
        let rounded = ogeom::fillet::fillet_edges(&mut model, rib, named, radius, T).unwrap();
        let (after, held) = measured(&model, &rounded.shape);
        assert!(held.diagonal() - bounds.diagonal() <= 1e-6);
        before - after
    };
    let alone = took(std::slice::from_ref(&spline));
    assert!(
        (alone - band).abs() <= band * 1e-4,
        "the spline's round took {alone} where its section along it is {band}"
    );
    let rest = took(&chain);
    let mut whole = chain.to_vec();
    whole.push(spline);
    let added = took(&whole) - rest;
    assert!(
        (added - band - cells).abs() <= (band + cells) * 1e-4,
        "the spline's round added {added} to the chain where its band and the corners \
         take {}",
        band + cells
    );
}

/// A pad fused onto a part where the pad's curved side touches a wall
/// swept from a quadratic profile, at the point where both their profiles
/// end on their apexes and the pad's foot stands on the wall's end. The
/// section between the two curved faces runs into that touch ever more
/// shallowly; it ends on the touch, where the faces' boundary edges meet
/// it, and the fuse closes. What the part gains is the pad less what of
/// it the part already held.
#[test]
fn a_pad_touching_a_swept_wall_where_both_profiles_end_fuses() {
    let (mut model, roots) =
        ogeom::io::native::read(&corpus("pad_arc_touching_a_swept_wall_at_its_end.ogeom")).unwrap();
    let (part, pad) = (&roots[0], &roots[1]);
    let (before, _) = measured(&model, part);
    let (pad_volume, _) = measured(&model, pad);
    let mut held = model.clone();
    let fused = ogeom::boolean::fuse(&mut model, part, pad, T).unwrap();
    let (after, _) = measured(&model, &fused.shape);
    let common = ogeom::boolean::common(&mut held, part, pad, T).unwrap();
    let (shared, _) = measured(&held, &common.shape);
    // The pad's section is the corner its quadratic side cuts off a 2 mm
    // square: a third of the square's half, 2/3 mm², swept 3.92 mm. Its
    // faces are planes and one swept quadratic, which the mesh measures
    // to well under a ten-thousandth.
    assert!((pad_volume - 2.0 / 3.0 * 3.92).abs() <= 1e-4);
    // Inclusion and exclusion, to a hundredth of what the pad adds: the
    // meshes measure the part's spline faces to a few millionths.
    let gained = pad_volume - shared;
    assert!(
        (after - before - gained).abs() <= gained * 1e-2,
        "the part gained {} where the pad stands {gained} outside it",
        after - before
    );
}

/// A part standing on its mirror plane fused with its own mirror image. Its
/// half-round channel is centred on that plane, so the two halves meet in
/// the channel's two halves and in the faces they share on the plane, and
/// the image's edges, on left-handed frames, split the part's channel as
/// contacts. The slanted bore's walls carry fitted edges ending where the
/// bore's chart closes, which the image's bake refits. The halves share
/// only faces on the plane, so their volumes add.
#[test]
fn a_part_on_its_mirror_plane_fuses_with_its_mirror_image() {
    let (mut model, roots) =
        ogeom::io::native::read(&corpus("half_channel_on_the_mirror_plane.ogeom")).unwrap();
    let part = &roots[0];
    let (half, _) = measured(&model, part);
    let copy = ogeom::algo::copied(&mut model, part).unwrap().shape;
    let mirror =
        ogeom::math::Transform::plane_mirror(ogeom::math::Point::ORIGIN, ogeom::math::Direction::X);
    let image = ogeom::algo::transformed(&mut model, &copy, mirror)
        .unwrap()
        .shape;
    let fused = ogeom::boolean::fuse(&mut model, part, &image, T).unwrap();
    let (whole, bounds) = measured(&model, &fused.shape);
    assert_eq!(edges_walked_one_way(&model, &fused.shape), 0);
    // Both volumes are the same exact integrals over the same faces and
    // their images: rounding, and the fitted edges' micron.
    assert!(
        (whole - 2.0 * half).abs() <= half * 1e-8,
        "the whole measures {whole} against twice the half's {half}"
    );
    // The part reaches the plane from one side; the whole spans it.
    assert!(bounds.low().unwrap().x < -3.0 && bounds.high().unwrap().x > 3.0);
}
