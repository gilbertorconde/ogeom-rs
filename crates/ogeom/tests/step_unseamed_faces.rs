#![allow(
    clippy::unwrap_used,
    reason = "test code; a failed unwrap is a failed test"
)]

//! Faces a STEP file bounds by rings alone, leaving the rest of their
//! boundary (a pole, a seam) for the reader to put in.

use ogeom::core::Tolerances;
use ogeom::math::Point;
use ogeom::mesh::{Deflection, face_meshes_for};

const T: Tolerances = Tolerances::millimetres();

/// A STEP file of one solid: the context entities, then `body`, whose
/// entities are numbered from 100 and whose last names the closed shell.
fn file(body: &[String]) -> String {
    let mut out = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test'),'2;1');\n\
         FILE_NAME('t.step','2026-01-01T00:00:00',(''),(''),'','','');\n\
         FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n\
         #1=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n\
         #2=( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );\n\
         #3=( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() );\n\
         #4=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-7),#1,'distance_accuracy_value','');\n\
         #5=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#4)) \
         GLOBAL_UNIT_ASSIGNED_CONTEXT((#1,#2,#3)) REPRESENTATION_CONTEXT('','3D') );\n\
         #6=APPLICATION_CONTEXT('automotive design');\n\
         #8=PRODUCT_CONTEXT('',#6,'mechanical');\n\
         #9=PRODUCT_DEFINITION_CONTEXT('part definition',#6,'design');\n\
         #10=PRODUCT('p','p','',(#8));\n#11=PRODUCT_DEFINITION_FORMATION('','',#10);\n\
         #12=PRODUCT_DEFINITION('design','',#11,#9);\n#13=PRODUCT_DEFINITION_SHAPE('','',#12);\n\
         #14=CARTESIAN_POINT('',(0.,0.,0.));\n#15=DIRECTION('',(0.,0.,1.));\n\
         #16=DIRECTION('',(1.,0.,0.));\n#17=AXIS2_PLACEMENT_3D('',#14,#15,#16);\n",
    );
    for (i, entity) in body.iter().enumerate() {
        out.push_str(&format!("#{}={entity};\n", 100 + i));
    }
    let shell = 99 + body.len();
    out.push_str(&format!(
        "#90=MANIFOLD_SOLID_BREP('',#{shell});\n\
         #91=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#17,#90),#5);\n\
         #92=SHAPE_DEFINITION_REPRESENTATION(#13,#91);\nENDSEC;\nEND-ISO-10303-21;\n"
    ));
    out
}

/// Entities numbered as [`file`] numbers them, and named by what they are.
#[derive(Default)]
struct Body(Vec<String>);

impl Body {
    fn add(&mut self, entity: impl Into<String>) -> usize {
        self.0.push(entity.into());
        99 + self.0.len()
    }

    fn point(&mut self, p: (f64, f64, f64)) -> usize {
        self.add(format!(
            "CARTESIAN_POINT('',({:?},{:?},{:?}))",
            p.0, p.1, p.2
        ))
    }

    /// An axis placement at `origin` along z, x along x.
    fn frame(&mut self, origin: (f64, f64, f64)) -> usize {
        let o = self.point(origin);
        let z = self.add("DIRECTION('',(0.,0.,1.))");
        let x = self.add("DIRECTION('',(1.,0.,0.))");
        self.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{z},#{x})"))
    }

    /// A closed circle edge of `radius` about z at height `z`, its vertex
    /// on +x.
    fn ring(&mut self, radius: f64, z: f64) -> usize {
        let frame = self.frame((0.0, 0.0, z));
        let circle = self.add(format!("CIRCLE('',#{frame},{radius:?})"));
        let at = self.point((radius, 0.0, z));
        let vertex = self.add(format!("VERTEX_POINT('',#{at})"));
        self.add(format!("EDGE_CURVE('',#{vertex},#{vertex},#{circle},.T.)"))
    }

    /// A face of `surface` bounded by `rings`, each walked as `forward`
    /// says, the face along its surface as `same_sense` says.
    fn face(&mut self, surface: usize, rings: &[(usize, bool)], same_sense: bool) -> usize {
        let bounds: Vec<String> = rings
            .iter()
            .enumerate()
            .map(|(i, &(edge, forward))| {
                let sense = if forward { ".T." } else { ".F." };
                let used = self.add(format!("ORIENTED_EDGE('',*,*,#{edge},{sense})"));
                let lp = self.add(format!("EDGE_LOOP('',(#{used}))"));
                let kind = if i == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                };
                format!("#{}", self.add(format!("{kind}('',#{lp},.T.)")))
            })
            .collect();
        let sense = if same_sense { ".T." } else { ".F." };
        self.add(format!(
            "ADVANCED_FACE('',({}),#{surface},{sense})",
            bounds.join(",")
        ))
    }

    /// A horizontal plane at height `z`.
    fn plane(&mut self, z: f64) -> usize {
        let frame = self.frame((0.0, 0.0, z));
        self.add(format!("PLANE('',#{frame})"))
    }
}

/// Every face meshed at `chord`: none refused or empty, every vertex on the
/// true surface `off` measures against, and the volume the triangles
/// enclose, short of the true one by the sag of a few thousand square
/// millimetres of curved face at most.
fn meshed(text: &str, chord: f64, off: impl Fn(Point) -> f64) -> f64 {
    let import = ogeom::io::read_step(text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let model = import.document.model();
    let deflection = Deflection {
        chord,
        ..Deflection::default()
    };
    let (faces, _) = face_meshes_for(model, &import.solids[0], deflection, None, T).unwrap();
    let mut volume = 0.0;
    for (_, mesh) in faces {
        let mesh = mesh.unwrap();
        assert!(!mesh.triangles.is_empty(), "a face meshed to nothing");
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
            volume += a.to_vector().dot(b.to_vector().cross(c.to_vector())) / 6.0;
            for p in [a, b, c] {
                assert!(off(p) <= 1e-6, "chord {chord}: {p:?} stands {} off", off(p));
            }
        }
    }
    volume
}

/// A ball face bounded by its equator alone, walked so its material lies
/// above: the upper hemisphere, the pole the file never wrote put in on
/// that side, not the lower one.
#[test]
fn a_ball_bounded_by_its_equator_is_the_half_its_ring_walks_round() {
    let mut body = Body::default();
    let equator = body.ring(20.0, 0.0);
    let floor = body.plane(0.0);
    let base = body.face(floor, &[(equator, false)], false);
    let frame = body.frame((0.0, 0.0, 0.0));
    let ball = body.add(format!("SPHERICAL_SURFACE('',#{frame},20.)"));
    let dome = body.face(ball, &[(equator, true)], true);
    body.add(format!("CLOSED_SHELL('',(#{base},#{dome}))"));
    let text = file(&body.0);
    let want = 2.0 / 3.0 * core::f64::consts::PI * 20f64.powi(3);
    for chord in [0.0125, 0.1, 0.5] {
        let off = |p: Point| {
            if p.z.abs() < 1e-9 {
                (p.x.hypot(p.y) - 20.0).max(0.0)
            } else if p.z < 0.0 {
                f64::INFINITY
            } else {
                (p.to_vector().magnitude() - 20.0).abs()
            }
        };
        let v = meshed(&text, chord, off);
        assert!(
            v > want - 3000.0 * chord && v < want,
            "chord {chord}: {v} against {want}"
        );
    }
}
