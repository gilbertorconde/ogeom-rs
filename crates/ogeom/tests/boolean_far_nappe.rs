//! A cone face on the far nappe of its surface, past the apex where the
//! chart's radius `R + v·tan α` runs negative, cut by planes, cylinders, a
//! ball and boxes. Its chart's window starts half a turn round from the
//! surface's own zero, and the face's trim closes at the apex's pole.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// The conical surface sits at the origin with radius 1 and half-angle π/4,
/// so its apex is at `z = -1`; the face runs from there down to its rim of
/// radius 1 at `z = -2`, where the chart's radius is -1. The solid is the
/// cone of radius 1 and height 1, of volume π/3.
const FAR_NAPPE_CONE: &str = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('far nappe','2026-10-28',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CONICAL_SURFACE('',#4,1.,0.7853981633974483);
#6=CARTESIAN_POINT('',(0.,0.,-2.));
#7=AXIS2_PLACEMENT_3D('',#6,#2,#3);
#8=PLANE('',#7);
#9=CARTESIAN_POINT('',(0.,0.,-1.));
#10=CARTESIAN_POINT('',(1.,0.,-2.));
#11=VERTEX_POINT('',#9);
#12=VERTEX_POINT('',#10);
#13=DIRECTION('',(0.7071067811865476,0.,-0.7071067811865476));
#14=VECTOR('',#13,1.);
#15=LINE('',#9,#14);
#16=CIRCLE('',#7,1.);
#17=EDGE_CURVE('',#11,#12,#15,.T.);
#18=EDGE_CURVE('',#12,#12,#16,.T.);
#19=ORIENTED_EDGE('',*,*,#17,.T.);
#20=ORIENTED_EDGE('',*,*,#18,.T.);
#21=ORIENTED_EDGE('',*,*,#17,.F.);
#22=EDGE_LOOP('',(#19,#20,#21));
#23=FACE_OUTER_BOUND('',#22,.T.);
#24=ADVANCED_FACE('',(#23),#5,.T.);
#25=ORIENTED_EDGE('',*,*,#18,.F.);
#26=EDGE_LOOP('',(#25));
#27=FACE_OUTER_BOUND('',#26,.T.);
#28=ADVANCED_FACE('',(#27),#8,.F.);
#29=CLOSED_SHELL('',(#24,#28));
#30=MANIFOLD_SOLID_BREP('',#29);
ENDSEC;
END-ISO-10303-21;
"#;

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn frame(at: [f64; 3], z: [f64; 3], x: [f64; 3]) -> Frame {
    let direction = |v: [f64; 3]| Direction::new(Vector::new(v[0], v[1], v[2]), T).unwrap();
    Frame::new(
        Point::new(at[0], at[1], at[2]),
        direction(z),
        direction(x),
        T,
    )
    .unwrap()
}

/// The far-nappe cone, read into a model of its own.
fn far_nappe_cone() -> (Model, Shape) {
    let import = ogeom::io::read_step(FAR_NAPPE_CONE, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    (import.document.model().clone(), import.solids[0].clone())
}

/// The cut, the common and the fuse of `a` with `b`, and the cut of `b`
/// by `a`, each checked valid, by volume.
fn four_volumes(model: &mut Model, a: &Shape, b: &Shape) -> [f64; 4] {
    let shapes = [
        ogeom::boolean::cut(model, a, b, T).unwrap().shape,
        ogeom::boolean::common(model, a, b, T).unwrap().shape,
        ogeom::boolean::fuse(model, a, b, T).unwrap().shape,
        ogeom::boolean::cut(model, b, a, T).unwrap().shape,
    ];
    shapes.map(|shape| {
        let diagnosis = ogeom::algo::check(model, &shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{diagnosis}");
        volume(model, &shape)
    })
}

/// Slabs across the axis, a coaxial drum and a half-space through the axis
/// against the far-nappe cone, each meeting it where the closed forms say.
///
/// A slab's plane meets the face in one circle that never reaches the apex.
/// Each cut keeps a frustum `πh(R² + Rr + r²)/3` or the cone above it, and
/// the drum of radius 1/2 keeps the drum up to `z = -1.5` and the cone of
/// radius 1/2 above. The face's trim spans the chart from half a turn, so
/// a section's half that starts in the surface's first turn is folded a
/// turn on to lie inside it.
#[test]
fn a_cone_face_past_its_apex_cuts_to_its_closed_forms() {
    let (mut model, cone) = far_nappe_cone();
    let frustum = |h: f64, big: f64, small: f64| {
        PI * h * small.mul_add(small, big.mul_add(big, big * small)) / 3.0
    };
    let whole = PI / 3.0;
    let cases = [
        (
            "slab above z = -1.5",
            ogeom::algo::make_box(
                &mut model,
                frame([-2.0, -2.0, -1.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                (4.0, 4.0, 2.0),
                T,
            ),
            32.0,
            frustum(0.5, 0.5, 0.0),
        ),
        (
            "slab below z = -1.5",
            ogeom::algo::make_box(
                &mut model,
                frame([-2.0, -2.0, -3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                (4.0, 4.0, 1.5),
                T,
            ),
            24.0,
            frustum(0.5, 1.0, 0.5),
        ),
        (
            "slab z in [-1.8, -1.2]",
            ogeom::algo::make_box(
                &mut model,
                frame([-2.0, -2.0, -1.8], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                (4.0, 4.0, 0.6),
                T,
            ),
            9.6,
            frustum(0.6, 0.8, 0.2),
        ),
        (
            "coaxial drum of radius 1/2",
            ogeom::algo::make_cylinder(
                &mut model,
                frame([0.0, 0.0, -3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                0.5,
                3.0,
                T,
            ),
            0.75 * PI,
            0.125f64.mul_add(PI, frustum(0.5, 0.5, 0.0)),
        ),
        (
            "half-space x >= 0",
            ogeom::algo::make_box(
                &mut model,
                frame([0.0, -2.0, -3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                (2.0, 4.0, 3.0),
                T,
            ),
            24.0,
            whole / 2.0,
        ),
    ];
    for (name, tool, tool_volume, common) in cases {
        let tool = tool.unwrap().shape;
        let got = four_volumes(&mut model, &cone, &tool);
        let want = [
            whole - common,
            common,
            whole + tool_volume - common,
            tool_volume - common,
        ];
        for (got, want) in got.iter().zip(want) {
            assert!(
                (got - want).abs() <= 1e-9 * want.max(1.0),
                "{name}: {got:?} against {want}"
            );
        }
    }
}

/// Tools meeting the far-nappe cone where no closed form is at hand: a
/// half-space off the axis (a hyperbola), a tilted slab (an ellipse), a drum
/// off the axis and one across it, a ball on the rim and a box standing on
/// an edge. The near-nappe cone `make_cone` builds is the same solid on the
/// other nappe of a surface of its own, and each result matches its twin's.
#[test]
fn a_cone_face_past_its_apex_cuts_as_its_near_nappe_twin() {
    let (mut model, far) = far_nappe_cone();
    let near = ogeom::algo::make_cone(
        &mut model,
        frame([0.0, 0.0, -2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        1.0,
        0.0,
        1.0,
        T,
    )
    .unwrap()
    .shape;
    let half = 0.5f64.sqrt();
    let tools = [
        ogeom::algo::make_box(
            &mut model,
            frame([0.3, -2.0, -3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
            (2.0, 4.0, 3.0),
            T,
        ),
        ogeom::algo::make_box(
            &mut model,
            frame([3.0, -3.0, -2.2], [0.2, 0.0, 1.0], [0.0, 1.0, 0.0]),
            (6.0, 6.0, 3.0),
            T,
        ),
        ogeom::algo::make_cylinder(
            &mut model,
            frame([0.6, 0.0, -3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
            0.3,
            3.0,
            T,
        ),
        ogeom::algo::make_cylinder(
            &mut model,
            frame([-2.0, 0.1, -1.6], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            0.2,
            4.0,
            T,
        ),
        ogeom::algo::make_sphere(
            &mut model,
            frame([0.5, 0.0, -2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
            0.4,
            T,
        ),
        ogeom::algo::make_box(
            &mut model,
            frame([0.0, 0.0, -1.5], [half, 0.0, half], [0.0, 1.0, 0.0]),
            (1.0, 1.0, 1.0),
            T,
        ),
    ];
    for (k, tool) in tools.into_iter().enumerate() {
        let tool = tool.unwrap().shape;
        let got = four_volumes(&mut model, &far, &tool);
        let want = four_volumes(&mut model, &near, &tool);
        for (got, want) in got.iter().zip(want) {
            assert!(
                (got - want).abs() <= 1e-6 * want.max(1.0),
                "tool {k}: {got} against the twin's {want}"
            );
        }
    }
}
