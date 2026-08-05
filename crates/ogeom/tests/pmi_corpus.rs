//! AP242 semantic PMI: the values, not the leader lines — read from NIST's
//! own annotated part and carried through the writer intact.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::doc::{MeasureKind, Pmi};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// One dimension, flattened for comparison.
type DimensionRow = (String, Vec<String>, Option<String>, usize);
/// One tolerance, flattened for comparison.
type ToleranceRow = (String, String, String, Vec<String>, usize);
/// One datum, flattened for comparison.
type DatumRow = (String, usize);

/// The comparable core of a PMI set: everything except node identities,
/// which are model-scoped and change across a round trip.
fn summary(pmi: &Pmi) -> (Vec<DimensionRow>, Vec<ToleranceRow>, Vec<DatumRow>) {
    let f = |v: f64| format!("{v:.9}");
    let mut dims: Vec<_> = pmi
        .dimensions
        .iter()
        .map(|d| {
            (
                format!("{}|{:?}|{}", d.name, d.kind, d.location),
                d.values.iter().map(|&v| f(v)).collect::<Vec<_>>(),
                d.plus
                    .map(f)
                    .zip(d.minus.map(f))
                    .map(|(p, m)| format!("{p}/{m}")),
                d.items.len(),
            )
        })
        .collect();
    dims.sort();
    let mut tols: Vec<_> = pmi
        .tolerances
        .iter()
        .map(|t| {
            (
                t.kind.clone(),
                t.name.clone(),
                f(t.magnitude),
                t.datums.clone(),
                t.items.len(),
            )
        })
        .collect();
    tols.sort();
    let mut datums: Vec<_> = pmi
        .datums
        .iter()
        .map(|d| (d.label.clone(), d.items.len()))
        .collect();
    datums.sort();
    (dims, tols, datums)
}

#[test]
fn the_nist_part_reads_its_annotations_semantically() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let pmi = import.document.pmi();

    // Eight dimensional characteristics, values in millimetres.
    assert_eq!(pmi.dimensions.len(), 8);
    let diameters: Vec<f64> = pmi
        .dimensions
        .iter()
        .filter(|d| d.name == "diameter")
        .map(|d| d.values[0])
        .collect();
    assert_eq!(diameters, [35.0, 35.0, 20.0, 20.0, 35.0, 35.0, 25.0]);
    // The one angular dimension: sixty degrees, plus/minus half a degree,
    // radians throughout.
    let angle = pmi
        .dimensions
        .iter()
        .find(|d| d.kind == MeasureKind::Angle)
        .expect("the part dimensions an angle");
    assert!((angle.values[0] - 60_f64.to_radians()).abs() < 1e-12);
    assert!((angle.plus.unwrap() - 0.5_f64.to_radians()).abs() < 1e-12);

    // Six geometric tolerances, magnitudes and datum references as NIST
    // published them.
    let mut kinds: Vec<(&str, f64, &[String])> = pmi
        .tolerances
        .iter()
        .map(|t| (t.kind.as_str(), t.magnitude, t.datums.as_slice()))
        .collect();
    kinds.sort_by(|a, b| a.0.cmp(b.0).then(a.1.total_cmp(&b.1)));
    let a = ["A".to_string()];
    let none: [String; 0] = [];
    assert_eq!(
        kinds,
        [
            ("flatness", 0.2, &none[..]),
            ("perpendicularity", 1.5, &a[..]),
            ("position", 0.75, &a[..]),
            ("position", 0.75, &a[..]),
            ("surface_profile", 0.5, &a[..]),
            ("surface_profile", 1.25, &a[..]),
        ]
    );

    // Three datums: A, B, C.
    let labels: Vec<&str> = pmi.datums.iter().map(|d| d.label.as_str()).collect();
    assert_eq!(labels, ["A", "B", "C"]);
}

#[test]
fn pmi_round_trips_through_the_writer() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let first = ogeom::io::read_step(&text, T).unwrap();
    let written = ogeom::io::write_step(&first.document, T).unwrap();
    let second = ogeom::io::read_step(&written, T).unwrap();
    assert_eq!(
        summary(first.document.pmi()),
        summary(second.document.pmi())
    );
}
