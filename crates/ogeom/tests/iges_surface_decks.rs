//! IGES decks the way surface exporters write them: untrimmed patches
//! standing alone, trims given only in a surface of revolution's
//! parameters over a circular-arc generatrix, and copious data (106) drawn
//! as paths. Each reads to the area its closed form gives.
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
    params: String,
}

fn item(kind: i64, form: i64, owned: bool, params: &str) -> Item {
    Item {
        kind,
        form,
        owned,
        params: params.to_owned(),
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
        let status = if it.owned { "00010000" } else { "00000000" };
        let field = |s: &str| format!("{s:>8}");
        let first: String = [
            it.kind.to_string(),
            p_seq.to_string(),
            "0".into(),
            "0".into(),
            "0".into(),
            "0".into(),
            "0".into(),
            "0".into(),
            status.into(),
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
