//! IGES decks the way surface exporters write them: untrimmed patches
//! standing alone, trims given only in a surface of revolution's
//! parameters over a circular-arc generatrix, copious data (106) drawn as
//! paths, wireframe curves, and subfigures placed singly and in arrays.
//! Each reads to the length, area or volume its closed form gives.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::surface_properties;
use ogeom::core::Tolerances;
use ogeom::mesh::Deflection;

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// One entity of a deck: its type, form, whether another entity owns it
/// (the subordinate switch), and its parameters after the type number.
struct Item {
    kind: i64,
    form: i64,
    owned: bool,
    /// The directory entry of the transformation matrix (124) placing it,
    /// or 0.
    placed: usize,
    /// The entity use flag: 0 geometry, 5 a parameter-space curve.
    used: u8,
    params: String,
}

fn item(kind: i64, form: i64, owned: bool, params: &str) -> Item {
    Item {
        kind,
        form,
        owned,
        placed: 0,
        used: 0,
        params: params.to_owned(),
    }
}

impl Item {
    fn placed_by(mut self, de: usize) -> Self {
        self.placed = de;
        self
    }

    fn used_as(mut self, flag: u8) -> Self {
        self.used = flag;
        self
    }
}

/// A record padded to 72 columns with its section letter and number.
fn record(body: &str, section: char, seq: usize) -> String {
    format!("{body:<72}{section}{seq:>7}\n")
}

/// A millimetre deck of `items`, the `i`th at directory entry `2 i + 1`.
fn deck(items: &[Item]) -> String {
    let mut out = String::new();
    out += &record("surface deck", 'S', 1);
    out += &record(
        "1H,,1H;,4Htest,8Htest.igs,6Hogeom,6Hogeom,32,38,6,308,15,4Htest,",
        'G',
        1,
    );
    out += &record(
        "1.0,2,2HMM,1,0.01,15H20261021.120000,1D-06,10.0,2Hme,",
        'G',
        2,
    );
    out += &record("4Hhere,11,0,15H20261021.120000;", 'G', 3);
    let mut directory = String::new();
    let mut parameters = String::new();
    let mut p_seq = 1;
    for (i, it) in items.iter().enumerate() {
        let de = 2 * i + 1;
        let text = format!("{},{};", it.kind, it.params);
        // Split at commas into pieces of at most 64 columns.
        let mut lines: Vec<String> = Vec::new();
        let mut current = String::new();
        for piece in text.split_inclusive(',') {
            if current.len() + piece.len() > 64 {
                lines.push(std::mem::take(&mut current));
            }
            current += piece;
        }
        lines.push(current);
        let status = format!("00{:02}{:02}00", u8::from(it.owned), it.used);
        let field = |s: &str| format!("{s:>8}");
        let first: String = [
            it.kind.to_string(),
            p_seq.to_string(),
            "0".into(),
            "0".into(),
            "0".into(),
            "0".into(),
            it.placed.to_string(),
            "0".into(),
            status,
        ]
        .iter()
        .map(|s| field(s))
        .collect();
        let second: String = [
            it.kind.to_string(),
            "0".into(),
            "0".into(),
            lines.len().to_string(),
            it.form.to_string(),
            String::new(),
            String::new(),
            String::new(),
            "0".into(),
        ]
        .iter()
        .map(|s| field(s))
        .collect();
        directory += &record(&first, 'D', de);
        directory += &record(&second, 'D', de + 1);
        for line in lines {
            parameters += &format!("{line:<64}{de:>8}P{p_seq:>7}\n");
            p_seq += 1;
        }
    }
    out += &directory;
    out += &parameters;
    out += &record(
        &format!("S      1G      3D{:>7}P{:>7}", 2 * items.len(), p_seq - 1),
        'T',
        1,
    );
    out
}

/// The chord deflection a curved face's area is measured to.
const CHORD: f64 = 1e-5;

fn area(model: &ogeom::topo::Model, shape: &ogeom::topo::Shape) -> f64 {
    let fine = Deflection {
        chord: CHORD,
        ..Deflection::default()
    };
    surface_properties(model, shape, fine, T).unwrap().mass
}

/// How far a curved face's measured area may fall short of its own: a
/// chord sagging `CHORD` under an arc of radius `rho` is shorter than the
/// arc by under `2 CHORD / rho` of it, so each curved direction of the
/// mesh loses under that share of the area; twice over is a bound.
fn mesh_bound(area: f64, rho: f64) -> f64 {
    area * 4.0 * CHORD / rho
}

/// The sheets' total area, and the same read back through this crate's own
/// writer.
fn areas_and_round_trip(import: &ogeom::io::IgesImport) -> (f64, f64) {
    let model = import.document.model();
    let here: f64 = import.sheets.iter().map(|s| area(model, s)).sum();
    let text = ogeom::io::write_iges(&import.document, T).unwrap();
    let again = ogeom::io::read_iges(&text, T).unwrap();
    let there: f64 = again
        .sheets
        .iter()
        .map(|s| area(again.document.model(), s))
        .sum();
    (here, there)
}

/// The bilinear B-spline patch (128) over the square `[0, side]²` at height
/// `z`.
fn flat_patch(side: f64, z: f64, owned: bool) -> Item {
    item(
        128,
        0,
        owned,
        &format!(
            "1,1,1,1,0,0,1,0,0,0.,0.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,\
             0.,0.,{z},{side},0.,{z},0.,{side},{z},{side},{side},{z},0.,1.,0.,1."
        ),
    )
}

/// An untrimmed B-spline patch and an untrimmed surface of revolution,
/// both independent, read as faces bounded by their own extents; a patch a
/// trimmed surface already takes is not read a second time.
#[test]
fn independent_untrimmed_surfaces_read_as_their_whole_faces() {
    let items = [
        // D1: a square patch standing alone.
        flat_patch(10.0, 0.0, false),
        // D3: the drum of radius 2 and height 6 about z, a whole turn.
        item(120, 0, false, &format!("5,7,0.,{}", 2.0 * PI)),
        // D5, D7: its axis and generatrix.
        item(110, 0, true, "0.,0.,0.,0.,0.,1."),
        item(110, 0, true, "2.,0.,0.,2.,0.,6."),
        // D9: a trimmed surface with no boundary over D11, a patch whose
        // own status says it stands alone.
        item(144, 0, false, "11,0,0,0"),
        flat_patch(10.0, 20.0, false),
    ];
    let import = ogeom::io::read_iges(&deck(&items), T).unwrap();
    let (total, again) = areas_and_round_trip(&import);
    let want = 100.0 + 2.0 * PI * 2.0 * 6.0 + 100.0;
    assert_eq!(import.sheets.len(), 3, "{:?}", import.report);
    let bound = mesh_bound(want, 2.0);
    assert!(
        (total - want).abs() < bound,
        "{total} against {want}, bound {bound}"
    );
    assert!(
        (again - want).abs() < bound,
        "written and read back: {again}"
    );
    assert!(
        import.report.skipped.is_empty(),
        "{:?}",
        import.report.skipped
    );
}

/// A torus patch spelt as a surface of revolution (120) about the x axis
/// of a circular arc (100) in the xy plane about (0, R) of radius r, its
/// trim given only in the surface's parameters, as a closed copious-data
/// path (106 form 63): the arc's angle `θ` across, the angle turned `φ`
/// along. A point of the patch stands `R + r sin θ` from the axis, so the
/// patch's area is `r Δφ (R Δθ - r (cos θ₂ - cos θ₁))`. The arc's own
/// angle starts where it does: a quarter turn, and three quarters (where
/// the start's angle, read either way round, must be the format's).
#[test]
fn a_trim_on_a_surface_of_revolution_of_an_arc_lifts_through_the_arc_s_angle() {
    let (big, small) = (10.0, 2.0);
    for start in [PI / 2.0, 1.5 * PI] {
        let end = start + PI / 2.0;
        let at = |t: f64| (small * t.cos(), big + small * t.sin());
        let (s, e) = (at(start), at(end));
        let (t1, t2) = (start + 0.2, end - 0.3);
        let (f1, f2) = (0.1, 1.2);
        let items = [
            // D1: the trimmed surface over D3, bounded by D5.
            item(144, 0, false, "3,1,0,5"),
            // D3: the arc D9 turned about the x axis D7 a quarter turn.
            item(120, 0, true, &format!("7,9,0.,{}", PI / 2.0)),
            // D5: the boundary, in the surface's parameters only.
            item(142, 0, true, "0,3,11,0,1"),
            item(110, 0, true, "0.,0.,0.,1.,0.,0."),
            item(
                100,
                0,
                true,
                &format!("0.,0.,{big},{},{},{},{}", s.0, s.1, e.0, e.1),
            ),
            // D11: the rectangle (θ, φ) as a closed path in the plane.
            item(
                106,
                63,
                true,
                &format!("1,4,0.,{t1},{f1},{t2},{f1},{t2},{f2},{t1},{f2}"),
            ),
        ];
        let import =
            ogeom::io::read_iges(&deck(&items), T).unwrap_or_else(|e| panic!("start {start}: {e}"));
        assert_eq!(import.sheets.len(), 1, "start {start}: {:?}", import.report);
        let (got, again) = areas_and_round_trip(&import);
        let want = small * (f2 - f1) * (big * (t2 - t1) - small * (t2.cos() - t1.cos()));
        let bound = mesh_bound(want, small);
        assert!(
            (got - want).abs() < bound,
            "start {start}: {got} against {want}, bound {bound}"
        );
        assert!(
            (again - want).abs() < bound,
            "start {start}: written and read back, {again}"
        );
    }
}

/// A trim given only in model space as copious data (106 form 12, points
/// as triples) closing a triangle on a square patch: the triangle's area.
#[test]
fn a_copious_data_path_bounds_a_trimmed_surface() {
    let items = [
        item(144, 0, false, "3,1,0,5"),
        flat_patch(10.0, 0.0, true),
        item(142, 0, true, "0,3,0,7,2"),
        item(106, 12, true, "2,4,0.,0.,0.,8.,0.,0.,0.,6.,0.,0.,0.,0."),
    ];
    let import = ogeom::io::read_iges(&deck(&items), T).unwrap();
    assert_eq!(import.sheets.len(), 1, "{:?}", import.report);
    assert!(import.report.warnings.is_empty(), "{:?}", import.report);
    let (got, again) = areas_and_round_trip(&import);
    assert!((got - 24.0).abs() < 24.0 * 1e-12, "{got}");
    assert!(
        (again - 24.0).abs() < 24.0 * 1e-12,
        "written and read back: {again}"
    );
}

fn length(model: &ogeom::topo::Model, shape: &ogeom::topo::Shape) -> (f64, ogeom::math::Point) {
    let fine = Deflection {
        chord: CHORD,
        ..Deflection::default()
    };
    let props = ogeom::algo::linear_properties(model, shape, fine, T).unwrap();
    (props.mass, props.centre)
}

fn count(
    model: &ogeom::topo::Model,
    shape: &ogeom::topo::Shape,
    kind: ogeom::topo::ShapeType,
) -> usize {
    ogeom::topo::explore(model, shape, ogeom::topo::Filter::OfType(kind))
        .unwrap()
        .into_iter()
        .map(|s| s.node())
        .collect::<std::collections::BTreeSet<_>>()
        .len()
}

/// Independent curves read as wires: a line, a whole circle closing on
/// one vertex, a polyline, and a composite curve whose matrix lifts its
/// pieces 5 off the plane (the pieces' own matrices first, then the
/// composite's) with a point marking its start. A parameter-space curve
/// and a curve another entity owns are not wireframe. Each wire is its
/// curve's length, and comes back the same through the writer.
#[test]
fn independent_curves_read_as_wireframe() {
    let items = [
        // D1: a line of length 10.
        item(110, 0, false, "0.,0.,0.,10.,0.,0."),
        // D3: a whole circle of radius 3 about (0, 20).
        item(100, 0, false, "0.,0.,20.,3.,20.,3.,20."),
        // D5: the composite of D9, D11 and D13, placed by D7.
        item(102, 0, false, "3,9,11,13").placed_by(7),
        // D7: a translation by 5 along z.
        item(124, 0, true, "1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,5."),
        // D9: the point the composite starts at.
        item(116, 0, true, "0.,0.,0.,0"),
        item(110, 0, true, "0.,0.,0.,4.,0.,0."),
        item(110, 0, true, "4.,0.,0.,4.,3.,0."),
        // D15: a path of two segments, 3 and 4 long.
        item(106, 12, false, "2,3,0.,-10.,0.,3.,-10.,0.,3.,-6.,0."),
        // D17: a curve in some surface's parameters, standing alone.
        item(110, 0, false, "0.,0.,0.,1.,1.,0.").used_as(5),
        // D19: a curve another entity owns.
        item(110, 0, true, "0.,0.,0.,0.,0.,1."),
    ];
    let import = ogeom::io::read_iges(&deck(&items), T).unwrap();
    assert!(import.report.warnings.is_empty(), "{:?}", import.report);
    assert_eq!(
        import.report.skipped,
        [("type 110".to_owned(), 2)].into_iter().collect()
    );
    assert!(import.solids.is_empty() && import.sheets.is_empty());
    let model = import.document.model();
    let mut measured: Vec<(f64, ogeom::math::Point, usize, usize)> = import
        .wires
        .iter()
        .map(|w| {
            let (l, c) = length(model, w);
            let edges = count(model, w, ogeom::topo::ShapeType::Edge);
            let vertices = count(model, w, ogeom::topo::ShapeType::Vertex);
            (l, c, edges, vertices)
        })
        .collect();
    measured.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.z.total_cmp(&b.1.z)));
    assert_eq!(measured.len(), 4, "{measured:?}");
    let circle = 6.0 * PI;
    let wants = [(7.0, 2, 3), (7.0, 2, 3), (10.0, 1, 2), (circle, 1, 1)];
    for (got, want) in measured.iter().zip(wants) {
        let bound = if want.0 == circle {
            mesh_bound(circle, 3.0)
        } else {
            1e-12
        };
        assert!((got.0 - want.0).abs() <= bound, "{got:?} against {want:?}");
        assert_eq!((got.2, got.3), (want.1, want.2), "{got:?}");
    }
    // The polyline lies in z = 0; the composite stands 5 above it, its
    // centre the length-weighted middle of its pieces.
    assert!(
        measured[0]
            .1
            .distance(ogeom::math::Point::new(16.5 / 7.0, -62.0 / 7.0, 0.0))
            < 1e-12
    );
    assert!(
        measured[1]
            .1
            .distance(ogeom::math::Point::new(20.0 / 7.0, 4.5 / 7.0, 5.0))
            < 1e-12
    );

    // Written as parts and read back, each wire is the same length.
    let mut lengths: Vec<f64> = measured.iter().map(|m| m.0).collect();
    let mut document = import.document;
    for (i, wire) in import.wires.iter().enumerate() {
        document.add_part(format!("wire-{i}"), wire.clone());
    }
    let again = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
    assert!(again.report.skipped.is_empty(), "{:?}", again.report);
    let mut back: Vec<f64> = again
        .wires
        .iter()
        .map(|w| length(again.document.model(), w).0)
        .collect();
    lengths.sort_by(f64::total_cmp);
    back.sort_by(f64::total_cmp);
    assert_eq!(back.len(), 4);
    for (a, b) in lengths.iter().zip(&back) {
        assert!(
            (a - b).abs() <= mesh_bound(circle, 3.0),
            "{lengths:?} against {back:?}"
        );
    }
}

/// Subfigures placed singly and in arrays. A rectangular array (412) of a
/// definition holding an untrimmed unit patch, doubled, three columns 10
/// apart by two rows 20 apart, the grid turned a quarter about its corner
/// at (100, 0) and position 2 left out: five patches of area 4, counted up
/// each column. A circular array (414) of a definition holding a 2-cube
/// (150), four positions round a circle of radius 10 about (0, 50) of
/// which it keeps 1 and 3. A singular instance (408) of the same cube
/// under a quarter turn about z with its translation (0, -100, 0), which
/// the turn does not move.
#[test]
fn subfigure_arrays_place_their_base_at_each_kept_position() {
    let items = [
        // D1: the definition holding the patch D3.
        item(308, 0, true, "0,4HTILE,1,3"),
        flat_patch(1.0, 0.0, true),
        // D5: the definition holding the block D7.
        item(308, 0, true, "0,4HCUBE,1,7"),
        item(150, 0, true, "2.,2.,2.,0.,0.,0.,1.,0.,0.,0.,0.,1."),
        // D9: the rectangular array, leaving out position 2.
        item(
            412,
            0,
            false,
            &format!("1,2.,100.,0.,0.,3,2,10.,20.,{},1,1,2", PI / 2.0),
        ),
        // D11: the circular array, keeping positions 1 and 3.
        item(
            414,
            0,
            false,
            &format!("5,4,0.,50.,0.,10.,0.,{},2,0,1,3", PI / 2.0),
        ),
        // D13: the singular instance, turned by D15.
        item(408, 0, false, "5,0.,-100.,0.,1.").placed_by(15),
        item(124, 0, true, "0.,-1.,0.,0.,1.,0.,0.,0.,0.,0.,1.,0."),
    ];
    let import = ogeom::io::read_iges(&deck(&items), T).unwrap();
    assert!(import.report.warnings.is_empty(), "{:?}", import.report);
    assert!(import.report.skipped.is_empty(), "{:?}", import.report);
    let model = import.document.model();
    let fine = Deflection::default();
    let close = |a: ogeom::math::Point, b: [f64; 3]| {
        a.distance(ogeom::math::Point::new(b[0], b[1], b[2])) < 1e-9
    };

    let mut patches: Vec<(f64, ogeom::math::Point)> = import
        .sheets
        .iter()
        .map(|s| {
            let props = surface_properties(model, s, fine, T).unwrap();
            (props.mass, props.centre)
        })
        .collect();
    patches.sort_by(|a, b| a.1.x.total_cmp(&b.1.x).then(a.1.y.total_cmp(&b.1.y)));
    let want = [
        [81.0, 11.0, 0.0],
        [81.0, 21.0, 0.0],
        [101.0, 1.0, 0.0],
        [101.0, 11.0, 0.0],
        [101.0, 21.0, 0.0],
    ];
    assert_eq!(patches.len(), want.len(), "{patches:?}");
    for ((area, centre), want) in patches.iter().zip(want) {
        assert!(
            (area - 4.0).abs() < 1e-9 && close(*centre, want),
            "{patches:?}"
        );
    }

    let mut cubes: Vec<(f64, ogeom::math::Point)> = import
        .solids
        .iter()
        .map(|s| {
            let props = ogeom::algo::volume_properties(model, s, fine, T).unwrap();
            (props.mass, props.centre)
        })
        .collect();
    cubes.sort_by(|a, b| a.1.y.total_cmp(&b.1.y).then(a.1.x.total_cmp(&b.1.x)));
    let want = [[-1.0, -99.0, 1.0], [-9.0, 51.0, 1.0], [11.0, 51.0, 1.0]];
    assert_eq!(cubes.len(), want.len(), "{cubes:?}");
    for ((volume, centre), want) in cubes.iter().zip(want) {
        assert!(
            (volume - 8.0).abs() < 1e-9 && close(*centre, want),
            "{cubes:?}"
        );
    }
}

/// Network subfigures and a group as an array's base. A network
/// definition (320) holds an untrimmed unit patch; its physical instance
/// (420) places it at (100, 0, 0) scaled by 3 on every axis (y's and z's
/// scale defaulted to x's): one patch of area 9 centred at (101.5, 1.5).
/// A logical instance (a schematic's) and one scaling its axes unlike
/// place nothing and say why. A rectangular array (412) of a group (402,
/// form 7) holding a unit patch and a 2-cube (150), two columns 10 apart
/// from (0, -50, 0), places both members at each position.
#[test]
fn network_subfigures_and_grouped_array_bases_place_their_members() {
    let items = [
        // D1: the network definition holding the patch D3, physical, with
        // no connect points.
        item(320, 0, true, "0,4HPIPE,1,3,2,3HP01,0,0"),
        flat_patch(1.0, 0.0, true),
        // D5: its instance, scaled by 3.
        item(420, 0, false, "1,100.,0.,0.,3.,,,0,3HP01,0,0"),
        // D7: a logical instance.
        item(420, 0, false, "1,0.,0.,0.,1.,1.,1.,1,3HP02,0,0"),
        // D9: an instance scaling y twice x.
        item(420, 0, false, "1,0.,200.,0.,1.,2.,1.,0,3HP03,0,0"),
        // D11: the group of the patch D13 and the block D15.
        item(402, 7, true, "2,13,15"),
        flat_patch(1.0, 0.0, true),
        item(150, 0, true, "2.,2.,2.,0.,0.,0.,1.,0.,0.,0.,0.,1."),
        // D17: the array of the group, every position kept.
        item(412, 0, false, "11,1.,0.,-50.,0.,2,1,10.,0.,0.,0"),
    ];
    let import = ogeom::io::read_iges(&deck(&items), T).unwrap();
    let warnings = &import.report.warnings;
    assert_eq!(warnings.len(), 2, "{:?}", import.report);
    assert!(warnings[0].starts_with("D7:") && warnings[0].contains("logical"));
    assert!(warnings[1].starts_with("D9:") && warnings[1].contains("unlike"));
    assert!(import.report.skipped.is_empty(), "{:?}", import.report);
    let model = import.document.model();
    let fine = Deflection::default();
    let close = |a: ogeom::math::Point, b: [f64; 3]| {
        a.distance(ogeom::math::Point::new(b[0], b[1], b[2])) < 1e-9
    };

    let mut patches: Vec<(f64, ogeom::math::Point)> = import
        .sheets
        .iter()
        .map(|s| {
            let props = surface_properties(model, s, fine, T).unwrap();
            (props.mass, props.centre)
        })
        .collect();
    patches.sort_by(|a, b| a.1.x.total_cmp(&b.1.x));
    let want = [
        (1.0, [0.5, -49.5, 0.0]),
        (1.0, [10.5, -49.5, 0.0]),
        (9.0, [101.5, 1.5, 0.0]),
    ];
    assert_eq!(patches.len(), want.len(), "{patches:?}");
    for ((area, centre), (want_area, want_centre)) in patches.iter().zip(want) {
        assert!(
            (area - want_area).abs() < 1e-9 && close(*centre, want_centre),
            "{patches:?}"
        );
    }

    let mut cubes: Vec<(f64, ogeom::math::Point)> = import
        .solids
        .iter()
        .map(|s| {
            let props = ogeom::algo::volume_properties(model, s, fine, T).unwrap();
            (props.mass, props.centre)
        })
        .collect();
    cubes.sort_by(|a, b| a.1.x.total_cmp(&b.1.x));
    let want = [[1.0, -49.0, 1.0], [11.0, -49.0, 1.0]];
    assert_eq!(cubes.len(), want.len(), "{cubes:?}");
    for ((volume, centre), want) in cubes.iter().zip(want) {
        assert!(
            (volume - 8.0).abs() < 1e-9 && close(*centre, want),
            "{cubes:?}"
        );
    }
}

/// Where each of `points` stands, sorted.
fn positions(model: &ogeom::topo::Model, points: &[ogeom::topo::Shape]) -> Vec<[f64; 3]> {
    let mut at: Vec<[f64; 3]> = points
        .iter()
        .map(|p| {
            let own = model.node(p).unwrap().data().as_vertex().unwrap().point;
            let q = p.transform(model.datums()).unwrap().apply(own);
            [q.x, q.y, q.z]
        })
        .collect();
    at.sort_by(|a, b| a.partial_cmp(b).unwrap());
    at
}

/// Stand-alone points read as vertices: one where the file puts it, one
/// carried by its matrix, and one a subfigure definition holds, placed by
/// its instance. A point that is annotation, or that another entity owns,
/// is not one. Written as parts and read back, each stands where it did.
#[test]
fn stand_alone_points_read_and_write_as_vertices() {
    let items = [
        // D1: a point at (1, 2, 3).
        item(116, 0, false, "1.,2.,3.,0"),
        // D3: a point at the origin, moved by D5 to (10, 0, -4).
        item(116, 0, false, "0.,0.,0.,0").placed_by(5),
        item(124, 0, true, "1.,0.,0.,10.,0.,1.,0.,0.,0.,0.,1.,-4."),
        // D7: a point drawn as annotation.
        item(116, 0, false, "5.,5.,5.,0").used_as(1),
        // D9: a point another entity owns.
        item(116, 0, true, "6.,6.,6.,0"),
        // D11: a definition holding the point D13, instanced by D15 at
        // (0, 0, 100) scaled by 2.
        item(308, 0, true, "0,4HMARK,1,13"),
        item(116, 0, true, "1.,1.,1.,0"),
        item(408, 0, false, "11,0.,0.,100.,2."),
    ];
    let import = ogeom::io::read_iges(&deck(&items), T).unwrap();
    assert!(import.report.warnings.is_empty(), "{:?}", import.report);
    assert_eq!(
        import.report.skipped,
        [("type 116".to_owned(), 2)].into_iter().collect()
    );
    assert!(import.solids.is_empty() && import.sheets.is_empty() && import.wires.is_empty());
    let want = vec![[1.0, 2.0, 3.0], [2.0, 2.0, 102.0], [10.0, 0.0, -4.0]];
    assert_eq!(positions(import.document.model(), &import.points), want);

    let mut document = import.document;
    for (i, point) in import.points.iter().enumerate() {
        document.add_part(format!("point-{i}"), point.clone());
    }
    let again = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
    assert!(again.report.skipped.is_empty(), "{:?}", again.report);
    assert_eq!(positions(again.document.model(), &again.points), want);
}
