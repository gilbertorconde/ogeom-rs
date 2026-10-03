//! Oblique drills through a machined part whose slot corners are partial
//! drums, their trims stated two turns below their charts, and whose thin
//! bores lie all but tangent inside the drill's wall. Each cut and common
//! is valid, the two add up to the part, and the common matches the
//! drill's share of the part's mesh as read. So does a drill lying in a
//! thin plate's underside between the rounds at its ends.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, Triangulation};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap()
        .mass
}

/// How much of a closed triangle mesh lies inside an unbounded cylinder:
/// the length the mesh encloses along each ray parallel to the axis,
/// integrated over the cylinder's disc by the midpoint rule in polar
/// coordinates. Each crossing adds its height along the axis signed by
/// which way the triangle faces, so the sum is the length inside, in any
/// order.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "grid indices"
)]
fn mesh_inside_cylinder(
    mesh: &Triangulation,
    origin: Point,
    (x, y, z): (Vector, Vector, Vector),
    radius: f64,
    (rings, spokes): (usize, usize),
) -> f64 {
    // Every vertex in the cylinder's own frame: (across, across, along).
    let local: Vec<[f64; 3]> = mesh
        .positions
        .iter()
        .map(|p| {
            let d = *p - origin;
            [d.dot(x), d.dot(y), d.dot(z)]
        })
        .collect();
    // Triangles bucketed by their extent across the axis, so a ray meets
    // few.
    let cell = 2.0;
    let side = (2.0 * radius / cell).ceil() as usize + 1;
    let bucket = |v: f64| ((v + radius) / cell).floor().clamp(0.0, (side - 1) as f64) as usize;
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); side * side];
    for (t, tri) in mesh.triangles.iter().enumerate() {
        let p = tri.map(|i| local[i as usize]);
        let (al, ah) = (
            p[0][0].min(p[1][0]).min(p[2][0]),
            p[0][0].max(p[1][0]).max(p[2][0]),
        );
        let (bl, bh) = (
            p[0][1].min(p[1][1]).min(p[2][1]),
            p[0][1].max(p[1][1]).max(p[2][1]),
        );
        if ah < -radius || bh < -radius || al > radius || bl > radius {
            continue;
        }
        for i in bucket(al)..=bucket(ah) {
            for k in bucket(bl)..=bucket(bh) {
                buckets[i * side + k].push(t);
            }
        }
    }
    let mut volume = 0.0;
    let (dr, dt) = (
        radius / rings as f64,
        core::f64::consts::TAU / spokes as f64,
    );
    for i in 0..rings {
        let r = (i as f64 + 0.5) * dr;
        for j in 0..spokes {
            // An odd offset keeps the rays off the lines the part's edges
            // run along.
            let t = (j as f64 + 0.25 * (i % 4) as f64 + 0.123) * dt;
            let (a, b) = (r * t.cos(), r * t.sin());
            let mut length = 0.0;
            for &tri in &buckets[bucket(a) * side + bucket(b)] {
                let [p, q, s] = mesh.triangles[tri].map(|i| local[i as usize]);
                // Barycentric coordinates of (a, b) in the projected
                // triangle.
                let det = (q[0] - p[0]) * (s[1] - p[1]) - (s[0] - p[0]) * (q[1] - p[1]);
                if det.abs() < 1e-300 {
                    continue;
                }
                let u = ((a - p[0]) * (s[1] - p[1]) - (s[0] - p[0]) * (b - p[1])) / det;
                let v = ((q[0] - p[0]) * (b - p[1]) - (a - p[0]) * (q[1] - p[1])) / det;
                if u < 0.0 || v < 0.0 || u + v > 1.0 {
                    continue;
                }
                let h = p[2] + u * (q[2] - p[2]) + v * (s[2] - p[2]);
                // The sign of the normal's component along the axis:
                // positive where the ray leaves the solid.
                length += h * det.signum();
            }
            volume += length * r * dr * dt;
        }
    }
    volume
}

/// The part healed as a caller imports it, its volume, and its meshes as
/// read at two chords, coarse then fine.
fn part() -> (ogeom::io::StepImport, Shape, f64, [Triangulation; 2]) {
    let mut import = ogeom::io::read_step(&corpus("nist_ftc_06_asme1_rd.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let meshes = [COARSE, FINE].map(|chord| {
        ogeom::mesh::triangulate(model, &solid, Deflection::with_chord(chord).unwrap(), T).unwrap()
    });
    let part = ogeom::heal::reanchor_periodic_rings(model, &solid, T)
        .unwrap()
        .0
        .shape;
    let whole = volume(model, &part);
    (import, part, whole, meshes)
}

const COARSE: f64 = 0.01;
const FINE: f64 = 0.002;

/// Cut and common of the part by a drill of `radius` from `start` along
/// `axis` for `length`, its frame's x square to the axis and to world x:
/// both valid, adding up to the part, the common within `share` of the
/// drill's share of the part as its meshes measure it.
fn drill(start: Point, axis: Vector, radius: f64, length: f64, share: f64) {
    let (mut import, part, whole, meshes) = part();
    let model = import.document.model_mut();
    let z = axis.normalized(T).unwrap();
    let x = z.cross(Vector::new(1.0, 0.0, 0.0)).normalized(T).unwrap();
    let y = z.cross(x);
    let frame = Frame::new(
        start,
        Direction::new(z, T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(model, frame, radius, length, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(model, &part, &drill, T).unwrap().shape;
    let common = ogeom::boolean::common(model, &part, &drill, T)
        .unwrap()
        .shape;
    for result in [&cut, &common] {
        assert!(ogeom::algo::check(model, result, T).unwrap().is_valid());
    }
    let (a, b) = (volume(model, &cut), volume(model, &common));
    assert!(
        (a + b - whole).abs() < whole * 1e-4,
        "{a} + {b} against {whole}"
    );
    // A mesh stands off the part's curved faces by up to its chord, so
    // what it measures inside the drill misses by an amount proportional
    // to the chord: extrapolated from two chords to none. Sampled half as
    // finely each way, either measure moves by under 6e-5 of itself; the
    // coarse one alone misses the common by up to 3.5e-4, the extrapolated
    // one by under 5e-6.
    let [coarse, fine] =
        meshes.map(|mesh| mesh_inside_cylinder(&mesh, start, (x, y, z), radius, (300, 1200)));
    let inside = fine + (fine - coarse) * FINE / (COARSE - FINE);
    assert!(
        (b - inside).abs() < b * share,
        "common {b} against {inside} measured on the meshes as read"
    );
}

/// The drill's wall passes a few hundredths of a millimetre outside a thin
/// bore where the bore runs all but tangent inside it, and crosses the
/// slot's corner drums.
#[test]
fn an_oblique_drill_along_a_thin_bore_cuts() {
    drill(
        Point::new(
            -14.478_610_818_124_423,
            3.999_371_635_181_901_8,
            -277.138_401_510_994_7,
        ),
        Vector::new(
            -0.565_016_635_381_368_8,
            0.528_780_602_945_684_7,
            0.633_361_883_673_714_3,
        ),
        57.088_766_757_544_09,
        171.306_147_621_762_6,
        2e-5,
    );
}

/// A wider drill from below crossing the slot's corner drums.
#[test]
fn an_oblique_drill_through_the_slot_corners_cuts() {
    drill(
        Point::new(
            -106.231_280_769_777_24,
            -121.575_356_939_727_43,
            -467.682_084_637_834_9,
        ),
        Vector::new(
            0.451_720_459_210_370_6,
            0.616_916_458_860_542_9,
            0.644_486_409_206_619_3,
        ),
        59.733_482_940_789_8,
        594.045_655_070_655_9,
        2e-5,
    );
}

/// A drill along a thin plate, its lowest line lying in the plate's
/// underside, which runs into a round at each end of the plate. Each round
/// is tangent to the underside within a few microns of the drill, so the
/// drill's section with it is a figure eight crossing the round's tangent
/// edge twice at its double point, or two loops crossing that edge a
/// hundredth of a millimetre apart where the edge cannot be told from the
/// drill between them, which is one junction. Cut and common are valid, cut
/// and common add up to the part and the drill less the part and the common
/// add up to the drill. The common matches the drill's share of the part's
/// mesh as read: at a tangency the mesh's chord moves that share by the
/// root of the chord, so no extrapolation holds. Over chords from 2e-3 to
/// 1.25e-4 and grids up to sixteen times finer it lands from 3e-5 below
/// the common to 7e-5 above it, and 4e-5 above on the chord and grid used
/// here.
#[test]
fn a_drill_lying_in_a_plate_s_underside_between_its_rounds_cuts() {
    let mut import = ogeom::io::read_step(&corpus("nist_ctc_03_asme1_rc.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let mesh =
        ogeom::mesh::triangulate(model, &solid, Deflection::with_chord(FINE).unwrap(), T).unwrap();
    let part = ogeom::heal::reanchor_periodic_rings(model, &solid, T)
        .unwrap()
        .0
        .shape;
    let whole = volume(model, &part);
    let (radius, length) = (23.616_669_978_788_245, 586.769_672_450_339_8);
    let start = Point::new(
        -473.877_236_225_169_9,
        210.972_399_999_999_96,
        99.816_669_978_788_24,
    );
    let (x, z) = (Vector::new(0.0, 0.0, 1.0), Vector::new(1.0, 0.0, 0.0));
    let frame = Frame::new(
        start,
        Direction::new(z, T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(model, frame, radius, length, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(model, &part, &drill, T).unwrap().shape;
    let common = ogeom::boolean::common(model, &part, &drill, T)
        .unwrap()
        .shape;
    let rest = ogeom::boolean::cut(model, &drill, &part, T).unwrap().shape;
    for result in [&cut, &common, &rest] {
        assert!(ogeom::algo::check(model, result, T).unwrap().is_valid());
    }
    let (a, b, c) = (
        volume(model, &cut),
        volume(model, &common),
        volume(model, &rest),
    );
    assert!(
        (a + b - whole).abs() < whole * 1e-6,
        "{a} + {b} against {whole}"
    );
    let bore = core::f64::consts::PI * radius * radius * length;
    assert!(
        (c + b - bore).abs() < bore * 1e-6,
        "{c} + {b} against the drill's {bore}"
    );
    let inside = mesh_inside_cylinder(&mesh, start, (x, z.cross(x), z), radius, (300, 1200));
    assert!(
        (b - inside).abs() < b * 5e-5,
        "common {b} against {inside} measured on the mesh as read"
    );
}
