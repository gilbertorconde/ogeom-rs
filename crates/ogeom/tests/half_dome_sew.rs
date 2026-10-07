//! Two half domes sewn along their meridian seams leave their two half rims
//! free, whether the quarter circle they turn is centred on the axis or a
//! little inside the confusion distance off it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{copied, is_shell_closed, make_edge, make_revolution, sew, transformed};
use ogeom::core::Tolerances;
use ogeom::geom::CircleCurve;
use ogeom::math::{Axis, Circle, Direction, Frame, Point, Transform, Vector};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The quarter circle from (10, 0, 0) up to (0, 0, 10) in the XZ plane,
/// its centre `off` along -X and -Z from the origin and its radius grown
/// to still reach both ends within the confusion distance.
fn quarter(model: &mut Model, off: f64) -> Shape {
    let centre = Point::new(-off, 0.0, -off);
    let radius = 10.0 + off;
    let frame = Frame::new(
        centre,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let circle = CircleCurve::new(Circle::new(frame, radius, T).unwrap());
    make_edge(model, circle.into(), (0.0, core::f64::consts::FRAC_PI_2), T)
        .unwrap()
        .shape
}

/// The quarter circle turned by pi about the Z axis pointing `up` or down.
fn half_dome(model: &mut Model, off: f64, up: bool) -> Shape {
    let edge = quarter(model, off);
    let direction = if up { Direction::Z } else { -Direction::Z };
    let axis = Axis {
        location: Point::ORIGIN,
        direction,
    };
    make_revolution(model, &edge, axis, core::f64::consts::PI, T)
        .unwrap()
        .shape
}

fn faces(model: &Model, shapes: &[Shape]) -> Vec<Shape> {
    shapes
        .iter()
        .flat_map(|s| explore_unique(model, s, ShapeType::Face).unwrap())
        .collect()
}

/// Sew the faces and require the open dome: one shell, the two meridian
/// seams joined, the two half rims free.
fn sews_open(model: &mut Model, faces: &[Shape], what: &str) {
    let sewn = sew(model, faces, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "{what}");
    assert_eq!(sewn.joined, 2, "{what}");
    assert_eq!(sewn.free_edges.len(), 2, "{what}");
    assert!(!is_shell_closed(model, &sewn.shells[0]).unwrap(), "{what}");
}

#[test]
fn half_domes_about_opposite_axes_sew_into_an_open_dome() {
    for off in [0.0, 4.051_126_4e-8] {
        let mut model = Model::new();
        let a = half_dome(&mut model, off, true);
        let b = half_dome(&mut model, off, false);
        let faces = faces(&model, &[a, b]);
        sews_open(&mut model, &faces, &format!("centre {off:e} off the axis"));
    }
}

#[test]
fn a_half_dome_and_its_mirror_image_sew_into_an_open_dome() {
    for off in [0.0, 4.051_126_4e-8] {
        let mut model = Model::new();
        let a = half_dome(&mut model, off, true);
        let copy = copied(&mut model, &a).unwrap().shape;
        let b = transformed(
            &mut model,
            &copy,
            Transform::plane_mirror(Point::ORIGIN, Direction::Y),
        )
        .unwrap()
        .shape;
        let faces = faces(&model, &[a, b]);
        sews_open(&mut model, &faces, &format!("centre {off:e} off the axis"));
    }
}

/// The two half domes as the native format holds them, the quarter
/// circle's centre about 4e-8 off the axis.
const HALVES: &str = "\
ogeom 2
units 1.0
datum 0:0 -1.0 -1.2246467991473532e-16 0.0 1.2246467991473532e-16 -1.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0
datum 1:0 -1.0 1.2246467991473532e-16 0.0 -1.2246467991473532e-16 -1.0 -0.0 -0.0 0.0 1.0 1.0 0.0 0.0 0.0
curve 0:0 circle -4.0511263676989984e-8 0.0 -4.051126789651678e-8 1.0 0.0 4.051126773240051e-9 -4.051126773240051e-9 0.0 1.0 0.0 -1.0 0.0 10.000000040511264 0
curve 1:0 circle 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 10.0 0
curve 2:0 circle -4.0511263676989984e-8 0.0 -4.051126789651678e-8 1.0 0.0 4.051126773240051e-9 -4.051126773240051e-9 0.0 1.0 0.0 -1.0 0.0 10.000000040511264 0
curve 3:0 circle 0.0 0.0 0.0 1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 -1.0 10.0 0
pcurve 0:0 line2 0.0 0.0 1.0 0.0 0.0 3.141592653589793
pcurve 1:0 line2 0.0 1.5707963186926437 1.0 0.0 0.0 3.141592653589793
pcurve 2:0 line2 0.0 0.0 0.0 1.0 0.0 1.5707963186926437
pcurve 3:0 line2 3.141592653589793 0.0 0.0 1.0 0.0 1.5707963186926437
pcurve 4:0 line2 0.0 0.0 1.0 0.0 0.0 3.141592653589793
pcurve 5:0 line2 0.0 1.5707963186926437 1.0 0.0 0.0 3.141592653589793
pcurve 6:0 line2 0.0 0.0 0.0 1.0 0.0 1.5707963186926437
pcurve 7:0 line2 3.141592653589793 0.0 0.0 1.0 0.0 1.5707963186926437
surface 0:0 revolution 0.0 0.0 0.0 0.0 0.0 1.0 0.0 3.141592653589793 circle -4.0511263676989984e-8 0.0 -4.051126789651678e-8 1.0 0.0 4.051126773240051e-9 -4.051126773240051e-9 0.0 1.0 0.0 -1.0 0.0 10.000000040511264 0
surface 1:0 revolution 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 3.141592653589793 circle -4.0511263676989984e-8 0.0 -4.051126789651678e-8 1.0 0.0 4.051126773240051e-9 -4.051126773240051e-9 0.0 1.0 0.0 -1.0 0.0 10.000000040511264 0
entity 1 primitive 0 0
entity 2 primitive 1 0
entity 3 derived 1 1047 0
entity 4 primitive 1 0
entity 5 derived 1 1047 0
entity 6 primitive 1 0
entity 7 derived 1 1046 1 1
entity 8 primitive 1 0
entity 9 primitive 2 0
entity 10 derived 2 1047 0
entity 11 primitive 2 0
entity 12 derived 2 1047 0
entity 13 primitive 2 0
entity 14 derived 2 1046 1 8
node 0:0 vertex 1e-7 10.0 0.0 0.0 0
node 1:0 vertex 1e-7 0.0 0.0 10.0 0
node 2:0 edge 1e-7 0 0 3 2 0:0/F 1:0/F
r curve3d 0:0 - 0.0 1.5707963186926435
r pcurve 2:0 0:0 - 0.0 1.5707963186926437
r pcurve 3:0 0:0 0:0^1 0.0 1.5707963186926437
node 3:0 edge 1e-7 0 0 2 2 0:0/F 0:0/F/0:0^1
r curve3d 1:0 - 0.0 3.141592653589793
r pcurve 0:0 0:0 - 0.0 3.141592653589793
node 4:0 edge 1e-7 0 1 1 2 1:0/F 1:0/F/0:0^1
r pcurve 1:0 0:0 - 0.0 3.141592653589793
node 5:0 wire 4 3:0/F 2:0/F/0:0^1 4:0/R 2:0/R
node 6:0 face 1e-7 0:0 - 0 - 1 5:0/F
node 7:0 shell 1 6:0/R
node 8:0 vertex 1e-7 10.0 0.0 0.0 0
node 9:0 vertex 1e-7 0.0 0.0 10.0 0
node 10:0 edge 1e-7 0 0 3 2 8:0/F 9:0/F
r curve3d 2:0 - 0.0 1.5707963186926435
r pcurve 6:0 1:0 - 0.0 1.5707963186926437
r pcurve 7:0 1:0 1:0^1 0.0 1.5707963186926437
node 11:0 edge 1e-7 0 0 2 2 8:0/F 8:0/F/1:0^1
r curve3d 3:0 - 0.0 3.141592653589793
r pcurve 4:0 1:0 - 0.0 3.141592653589793
node 12:0 edge 1e-7 0 1 1 2 9:0/F 9:0/F/1:0^1
r pcurve 5:0 1:0 - 0.0 3.141592653589793
node 13:0 wire 4 11:0/F 10:0/F/1:0^1 12:0/R 10:0/R
node 14:0 face 1e-7 1:0 - 0 - 1 13:0/F
node 15:0 shell 1 14:0/R
node 16:0 compound 2 7:0/F 15:0/F
identity 2:0 1
identity 3:0 3
identity 4:0 5
identity 6:0 7
identity 10:0 8
identity 11:0 10
identity 12:0 12
identity 14:0 14
operation 2
root 16:0/F
";

#[test]
fn half_domes_read_from_the_native_format_sew_into_an_open_dome() {
    let mut model = Model::new();
    let root = ogeom::io::native::read_into(&mut model, HALVES)
        .unwrap()
        .shapes[0]
        .clone();
    let faces = explore_unique(&model, &root, ShapeType::Face).unwrap();
    sews_open(&mut model, &faces, "read");
}
