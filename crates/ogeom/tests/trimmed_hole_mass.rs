//! A block with a through bore whose face stands on a trimmed cylinder
//! rather than the cylinder itself, as a saved model can read back. The
//! closed forms leave such a face to its own mesh while the block's faces
//! stay exact, and the sum must still give the block less the bore: its
//! volume, its centre and its inertia, against the closed form.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{MassProperties, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{SurfaceGeometry, TrimmedSurface};
use ogeom::math::{Direction, Frame, Matrix3, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;
const SIZE: (f64, f64, f64) = (40.0, 30.0, 10.0);
const R: f64 = 3.0;

/// Every cylindrical face of `shape` moved onto its surface trimmed to its
/// full turn and `v`; how many surfaces were moved.
fn trim_cylinders(model: &mut Model, shape: &Shape, v: (f64, f64)) -> usize {
    let mut ids = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let NodeData::Face(data) = model.node(&face).unwrap().data() else {
            continue;
        };
        if matches!(
            model.geometry().surface(data.surface),
            Some(SurfaceGeometry::Cylinder(_))
        ) {
            ids.push(data.surface);
        }
    }
    ids.sort_unstable();
    ids.dedup();
    for &id in &ids {
        let surface = model.geometry_mut().surface_mut(id).unwrap();
        let trimmed = TrimmedSurface::new(surface.clone(), (0.0, 2.0 * PI), v, T).unwrap();
        *surface = trimmed.into();
    }
    ids.len()
}

/// The 40 x 30 x 10 block with a radius-3 bore along `axis` (0, 1 or 2)
/// through `at` (the bore's centre line passes through it), the bore's
/// face on a trimmed cylinder.
fn holed_block(axis: usize, at: Point) -> (Model, Shape) {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, SIZE, T)
        .unwrap()
        .shape;
    let extent = [SIZE.0, SIZE.1, SIZE.2][axis];
    let (direction, across) = match axis {
        0 => (Direction::X, Direction::Y),
        1 => (Direction::Y, Direction::Z),
        _ => (Direction::Z, Direction::X),
    };
    let mut start = [at.x, at.y, at.z];
    start[axis] = -1.0;
    let origin = Point::new(start[0], start[1], start[2]);
    let frame = Frame::new(origin, direction, across, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, frame, R, extent + 2.0, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    assert_eq!(trim_cylinders(&mut model, &part, (0.0, extent + 2.0)), 1);
    (model, part)
}

/// The inertia about the origin of a body of mass `m` centred at `c` with
/// inertia `own` about its centre.
fn about_origin(m: f64, c: Vector, own: [f64; 3]) -> [[f64; 3]; 3] {
    let d = [c.x, c.y, c.z];
    let dd = c.dot(c);
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            out[i][j] = delta * own[i] + m * (delta * dd - d[i] * d[j]);
        }
    }
    out
}

/// The block less the bore in closed form: mass, centre and the inertia
/// about the centre.
fn expected(axis: usize, at: Point) -> (f64, Point, [[f64; 3]; 3]) {
    let (a, b, c) = SIZE;
    let box_mass = a * b * c;
    let box_centre = Vector::new(a / 2.0, b / 2.0, c / 2.0);
    let box_own = [
        box_mass * (b * b + c * c) / 12.0,
        box_mass * (a * a + c * c) / 12.0,
        box_mass * (a * a + b * b) / 12.0,
    ];
    let length = [a, b, c][axis];
    let bore_mass = PI * R * R * length;
    let mut centre = [at.x, at.y, at.z];
    centre[axis] = length / 2.0;
    let bore_centre = Vector::new(centre[0], centre[1], centre[2]);
    let transverse = bore_mass * (3.0 * R * R + length * length) / 12.0;
    let mut bore_own = [transverse; 3];
    bore_own[axis] = bore_mass * R * R / 2.0;

    let mass = box_mass - bore_mass;
    let first = box_centre * box_mass - bore_centre * bore_mass;
    let g = first / mass;
    let outer = about_origin(box_mass, box_centre, box_own);
    let inner = about_origin(bore_mass, bore_centre, bore_own);
    // Back from the origin to the centre: the parallel axis theorem
    // reversed.
    let shift = about_origin(mass, g, [0.0; 3]);
    let mut inertia = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            inertia[i][j] = outer[i][j] - inner[i][j] - shift[i][j];
        }
    }
    (mass, Point::ORIGIN + g, inertia)
}

fn check_against_closed_form(axis: usize, at: Point, found: &MassProperties) {
    let (mass, centre, inertia) = expected(axis, at);
    let what = format!("bore along axis {axis} through {at:?}: {found:?}");
    // The bore's face is meshed, to a tenth of the chord asked for.
    assert!(
        (found.mass - mass).abs() < 1e-4 * mass,
        "volume {} against {mass}; {what}",
        found.mass
    );
    assert!(
        found.centre.distance(centre) < 1e-3,
        "centre {:?} against {centre:?}; {what}",
        found.centre
    );
    let scale = inertia[0][0].max(inertia[1][1]).max(inertia[2][2]);
    let Matrix3 { rows } = found.inertia;
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (rows[i][j] - inertia[i][j]).abs() < 1e-4 * scale,
                "inertia [{i}][{j}] {} against {}; {what}",
                rows[i][j],
                inertia[i][j]
            );
        }
    }
}

fn chord() -> Deflection {
    Deflection::with_chord(1e-2).unwrap()
}

/// The block centred on its bore along y: symmetric about y = 15, so the
/// centre must sit there.
#[test]
fn a_block_with_a_trimmed_cylinder_bore_has_its_centre_on_its_symmetry_plane() {
    let at = Point::new(20.0, 0.0, 5.0);
    let (model, part) = holed_block(1, at);
    let found = volume_properties(&model, &part, chord(), T).unwrap();
    assert!(
        (found.centre.y - 15.0).abs() < 1e-3,
        "centre {:?}",
        found.centre
    );
    check_against_closed_form(1, at, &found);
}

/// Bores off centre and along every axis: the centre and the inertia come
/// out as the block's less the bore's.
#[test]
fn a_trimmed_cylinder_bore_anywhere_gives_the_closed_form_moments() {
    for (axis, at) in [
        (1, Point::new(12.0, 0.0, 4.0)),
        (1, Point::new(31.0, 0.0, 6.5)),
        (0, Point::new(0.0, 9.0, 5.0)),
        (0, Point::new(0.0, 22.0, 3.5)),
        (2, Point::new(10.0, 8.0, 0.0)),
        (2, Point::new(27.0, 19.0, 0.0)),
    ] {
        let (model, part) = holed_block(axis, at);
        let found = volume_properties(&model, &part, chord(), T).unwrap();
        check_against_closed_form(axis, at, &found);
    }
}
