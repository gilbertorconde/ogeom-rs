//! Rolling corners (`FilletCorner::Rolling`): every vertex where rounds
//! meet closed by the ball rolling through it, on spheres and tori, never
//! a patch; and refused by name where the ball's sweep does not close it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{FilletCorner, fillet_edges, fillet_edges_corner};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn block(m: &mut Model, at: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    ogeom_algo::make_box(
        m,
        Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap(),
        size,
        T,
    )
    .unwrap()
    .shape
}

/// Blocks fused, their coplanar faces unified and their collinear edges
/// merged, as printCAD hands a stepped part over.
fn fused(m: &mut Model, blocks: &[Shape]) -> Shape {
    let mut whole = blocks[0].clone();
    for b in &blocks[1..] {
        whole = ogeom_bool::fuse(m, &whole, b, T).unwrap().shape;
    }
    let unified = ogeom_heal::unify_same_domain(m, &whole, T).unwrap().0.shape;
    let merged = ogeom_heal::merge_edges(m, &unified, T).unwrap().0.shape;
    explore_unique(m, &merged, ShapeType::Solid).unwrap()[0].clone()
}

/// The printCAD stairs: three blocks 14 deep in `y`, stepping up in `x`.
/// Its profile has six convex corners and two concave ones.
fn stairs(m: &mut Model) -> Shape {
    let a = block(m, (-15.0, -7.0, 0.0), (30.0, 14.0, 7.0));
    let b = block(m, (0.0, -7.0, 7.0), (15.0, 14.0, 7.0));
    let c = block(m, (7.0, -7.0, 14.0), (8.0, 14.0, 7.0));
    fused(m, &[a, b, c])
}

/// An L-bracket's profile (one step) 10 deep in `y`.
fn step(m: &mut Model) -> Shape {
    let a = block(m, (-10.0, -5.0, 0.0), (20.0, 10.0, 5.0));
    let b = block(m, (0.0, -5.0, 5.0), (10.0, 10.0, 5.0));
    fused(m, &[a, b])
}

fn edges_of(m: &Model, s: &Shape) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge).unwrap()
}

fn surface_of(m: &Model, face: &Shape) -> SurfaceGeometry {
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().clone()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The faces' surfaces by kind: planes, cylinders, spheres, tori, other.
fn kinds(m: &Model, s: &Shape) -> [usize; 5] {
    let mut n = [0; 5];
    for face in explore_unique(m, s, ShapeType::Face).unwrap() {
        let k = match surface_of(m, &face) {
            SurfaceGeometry::Plane(_) => 0,
            SurfaceGeometry::Cylinder(_) => 1,
            SurfaceGeometry::Sphere(_) => 2,
            SurfaceGeometry::Torus(_) => 3,
            _ => 4,
        };
        n[k] += 1;
    }
    n
}

/// The largest angle between the two faces' normals, and the largest gap
/// off either surface, over points read along every edge's own curve and
/// projected onto each face's surface: no pcurve is read.
fn worst_join(m: &Model, s: &Shape) -> (f64, f64) {
    worst_join_below(m, s, f64::INFINITY)
}

/// [`worst_join`] over the edges of `s` touching a curved face that do not
/// lie wholly above `z = ceiling`: the joins a request made, not the edges
/// it left sharp or a band's flush end in a face above.
fn worst_join_below(m: &Model, s: &Shape, ceiling: f64) -> (f64, f64) {
    let (mut angle, mut gap) = (0.0_f64, 0.0_f64);
    for e in edges_of(m, s) {
        let mut faces: Vec<Shape> = Vec::new();
        for f in ancestors_of(m, s, &e, ShapeType::Face).unwrap() {
            if !faces.iter().any(|g| g.is_same(&f)) {
                faces.push(f);
            }
        }
        // A seam has one face either side.
        let [f, g] = faces.as_slice() else {
            continue;
        };
        let flat = |face: &Shape| matches!(surface_of(m, face), SurfaceGeometry::Plane(_));
        if ceiling.is_finite() && flat(f) && flat(g) {
            continue;
        }
        let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        let (p, q) = (curve.point_at(a, T).unwrap(), curve.point_at(b, T).unwrap());
        if p.z.min(q.z) > ceiling {
            continue;
        }
        for k in 1..16 {
            let p = curve
                .point_at(a + (b - a) * f64::from(k) / 16.0, T)
                .unwrap();
            let mut normal = |face: &Shape| -> Vector {
                let surface = surface_of(m, face);
                let foot = ogeom_geom::project_on_surface(&surface, p, 16, T).unwrap();
                gap = gap.max(foot.distance);
                let (u, v) = foot.parameters;
                let (du, dv) = surface.d1_at(u, v, T).unwrap();
                let n = du.cross(dv);
                n / n.magnitude()
            };
            let n = normal(f);
            let o = normal(g);
            angle = angle.max(n.cross(o).magnitude().asin());
        }
    }
    (angle, gap)
}

/// What says two results are one solid: every vertex's point and the
/// faces' surface kinds, to the bit.
fn signature(m: &Model, s: &Shape) -> (Vec<[u64; 3]>, [usize; 5]) {
    let mut points: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    points.sort_unstable();
    (points, kinds(m, s))
}

/// The stairs rounded on all 24 edges at radius `r = 1`, each vertex
/// rolled.
///
/// The exact volume, by parts:
/// - the blocks, `30·14·7 + 15·14·7 + 8·14·7 = 5194`;
/// - each straight band's spandrel `(1 − π/4)·r²·L`, taken off along the
///   convex bands and put on along the concave ones. Every band ends a
///   radius short of a convex corner and two short of a concave one's
///   cap, where the torus takes over; the eight bands along `y` run 12,
///   and on each cap the eight profile bands run 28, 19, 6, 5, 5, 5, 13
///   and 5. Convex: `6·12 + 2·86 = 244`; concave: `2·12 = 24`;
/// - each of the 12 convex corners, the cube `r³` less the ball's
///   octant: `r³·(1 − π/6)` off;
/// - each of the 4 elbows (the two concave profile corners on each cap):
///   in the cell `2r × r × 2r` the elbow's straight bands leave between
///   them, the slice at depth `t` under the cap's band centres keeps the
///   square less a quarter disc of radius `2r − √(r² − t²)` about the
///   concave band's axis, so the cell keeps
///   `∫₀ʳ (4r² − (π/4)(2r − √(r² − t²))²) dt = r³(4 − 7π/6 + π²/4)`
///   of the `3r³` it held: `r³(1 − 7π/6 + π²/4)` on.
///
/// The total, `4966 + 157π/3 + π²`, is the same as the ball's rolling
/// body read by slices (the cap's offset region `E` of area `281 − 2π` and
/// perimeter `86 + 2π`, `12·(A + P + π) + 2∫₀¹ (A + P·ρ + π·ρ²) dt`).
fn stairs_volume(r: f64) -> f64 {
    let blocks = 5194.0;
    let spandrel = (1.0 - PI / 4.0) * r * r;
    let straight = -spandrel * 244.0 + spandrel * 24.0;
    let spheres = -12.0 * r.powi(3) * (1.0 - PI / 6.0);
    let elbows = 4.0 * r.powi(3) * (1.0 - 7.0 * PI / 6.0 + PI * PI / 4.0);
    blocks + straight + spheres + elbows
}

#[test]
fn the_stairs_roll_every_corner_on_spheres_and_tori() {
    let mut m = Model::new();
    let solid = stairs(&mut m);
    assert_eq!(
        explore_unique(&m, &solid, ShapeType::Face).unwrap().len(),
        10
    );
    let edges = edges_of(&m, &solid);
    assert_eq!(edges.len(), 24);
    let built = fillet_edges_corner(&mut m, &solid, &edges, 1.0, FilletCorner::Rolling, T).unwrap();
    let rolled = &built.shape;
    assert!(ogeom_algo::check(&m, rolled, T).unwrap().is_valid());
    // Every face a plane, a cylinder, a sphere or a torus: the ten faces,
    // a band per edge, a sphere per convex corner and a torus per elbow.
    assert_eq!(kinds(&m, rolled), [10, 24, 12, 4, 0]);
    // No two faces left on one surface across an edge.
    let (_, merged) = ogeom_heal::unify_same_domain(&mut m, rolled, T).unwrap();
    assert_eq!(merged, 0, "faces left split on one surface");
    let (angle, gap) = worst_join(&m, rolled);
    assert!(angle <= 1e-6, "a join turns {angle:.2e} rad");
    assert!(gap <= 1e-6, "an edge stands {gap:.2e} off its faces");
    let v = volume(&m, rolled);
    let exact = stairs_volume(1.0);
    assert!((exact - (4966.0 + 157.0 * PI / 3.0 + PI * PI)).abs() < 1e-9);
    assert!(
        (v - exact).abs() <= 1e-6 * exact,
        "volume {v}, exact {exact}"
    );
}

#[test]
fn the_stairs_roll_alike_in_any_edge_order() {
    let build = |order: &dyn Fn(Vec<Shape>) -> Vec<Shape>| {
        let mut m = Model::new();
        let solid = stairs(&mut m);
        let edges = order(edges_of(&m, &solid));
        let built =
            fillet_edges_corner(&mut m, &solid, &edges, 1.0, FilletCorner::Rolling, T).unwrap();
        signature(&m, &built.shape)
    };
    let forward = build(&|e| e);
    let backward = build(&|mut e| {
        e.reverse();
        e
    });
    let shuffled = build(&|e| {
        let (odd, even): (Vec<_>, Vec<_>) =
            e.into_iter().enumerate().partition(|(i, _)| i % 2 == 1);
        odd.into_iter().chain(even).map(|(_, s)| s).collect()
    });
    assert_eq!(forward, backward);
    assert_eq!(forward, shuffled);
}

#[test]
fn a_setback_corner_is_fillet_edges_itself() {
    let run = |corner: Option<FilletCorner>| {
        let mut m = Model::new();
        let solid = stairs(&mut m);
        let edges = edges_of(&m, &solid);
        let built = match corner {
            Some(c) => fillet_edges_corner(&mut m, &solid, &edges, 1.0, c, T),
            None => fillet_edges(&mut m, &solid, &edges, 1.0, T),
        };
        built
            .map(|b| signature(&m, &b.shape))
            .map_err(|e| e.to_string())
    };
    assert_eq!(run(Some(FilletCorner::Setback)), run(None));
    assert_eq!(FilletCorner::default(), FilletCorner::Setback);
}

fn refusal(m: &mut Model, solid: &Shape, edges: &[Shape]) -> String {
    match fillet_edges_corner(m, solid, edges, 1.0, FilletCorner::Rolling, T) {
        Ok(_) => panic!("a rolling corner was built where none closes"),
        Err(e) => e.to_string(),
    }
}

/// The edges of `s` meeting at the vertex at `p`.
fn edges_at(m: &Model, s: &Shape, p: Point) -> Vec<Shape> {
    edges_of(m, s)
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            [a, b].iter().any(|v| {
                m.node(v)
                    .unwrap()
                    .data()
                    .as_vertex()
                    .unwrap()
                    .point
                    .distance(p)
                    <= 1e-9
            })
        })
        .collect()
}

#[test]
fn a_vertex_of_four_edges_is_refused_by_name() {
    let mut m = Model::new();
    let points = [
        Point::new(-5.0, -5.0, 0.0),
        Point::new(5.0, -5.0, 0.0),
        Point::new(5.0, 5.0, 0.0),
        Point::new(-5.0, 5.0, 0.0),
        Point::new(0.0, 0.0, 6.0),
    ];
    let rings = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 4],
        vec![1, 2, 4],
        vec![2, 3, 4],
        vec![3, 0, 4],
    ];
    let (built, _) =
        ogeom_algo::make_faceted_solid(&mut m, &points, &rings, &[None; 5], T).unwrap();
    let solid = explore_unique(&m, &built.shape, ShapeType::Solid).unwrap()[0].clone();
    let edges = edges_of(&m, &solid);
    let why = refusal(&mut m, &solid, &edges);
    assert!(
        why.contains("(0.000000, 0.000000, 6.000000)") && why.contains("joins 4 edges"),
        "{why}"
    );
}

#[test]
fn a_vertex_on_a_curved_face_is_refused_by_name() {
    let mut m = Model::new();
    let body = block(&mut m, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    let drill = ogeom_algo::make_cylinder(
        &mut m,
        Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
        4.0,
        12.0,
        T,
    )
    .unwrap()
    .shape;
    let notched = ogeom_bool::cut(&mut m, &body, &drill, T).unwrap().shape;
    let solid = explore_unique(&m, &notched, ShapeType::Solid).unwrap()[0].clone();
    let at = Point::new(6.0, 10.0, 10.0);
    let edges = edges_at(&m, &solid, at);
    assert_eq!(edges.len(), 3);
    let why = refusal(&mut m, &solid, &edges);
    assert!(
        why.contains("(6.000000, 10.000000, 10.000000)") && why.contains("curved"),
        "{why}"
    );
}

#[test]
fn a_concave_edge_leaning_to_the_face_across_it_is_refused_by_name() {
    let mut m = Model::new();
    let solid = step(&mut m);
    // The front cut off by a plane leaning 15° back from upright: the
    // concave edge along `y` then meets it aslant.
    let a = 15.0_f64.to_radians();
    let normal = Vector::new(0.0, -a.cos(), a.sin());
    let x = Vector::new(1.0, 0.0, 0.0);
    let y = normal.cross(x);
    let on = Point::new(0.0, -4.0, 5.0);
    let origin = on - x * 50.0 - y * 50.0;
    let knife = ogeom_algo::make_box(
        &mut m,
        Frame::new(
            origin,
            Direction::new(normal, T).unwrap(),
            Direction::new(x, T).unwrap(),
            T,
        )
        .unwrap(),
        (100.0, 100.0, 50.0),
        T,
    )
    .unwrap()
    .shape;
    let cut = ogeom_bool::cut(&mut m, &solid, &knife, T).unwrap().shape;
    let solid = explore_unique(&m, &cut, ShapeType::Solid).unwrap()[0].clone();
    // The concave edge's front end: on the leaning face, at `x = 0`,
    // `z = 5`.
    let front = Point::new(0.0, -4.0, 5.0);
    let edges = edges_at(&m, &solid, front);
    assert_eq!(edges.len(), 3);
    let why = refusal(&mut m, &solid, &edges);
    assert!(
        why.contains("(0.000000, -4.000000, 5.000000)") && why.contains("square"),
        "{why}"
    );
}

#[test]
fn a_convex_and_a_concave_round_beside_a_sharp_edge_are_refused_by_name() {
    let mut m = Model::new();
    let solid = step(&mut m);
    let at = Point::new(0.0, -5.0, 5.0);
    // The concave edge along `y` and the convex front edge of the upper
    // block's wall, the front edge of the lower block's top left sharp.
    let edges: Vec<Shape> = edges_at(&m, &solid, at)
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            let d = c.point_at(b, T).unwrap() - c.point_at(a, T).unwrap();
            d.x.abs() < 1e-9
        })
        .collect();
    assert_eq!(edges.len(), 2);
    let why = refusal(&mut m, &solid, &edges);
    assert!(
        why.contains("(0.000000, -5.000000, 5.000000)") && why.contains("left sharp"),
        "{why}"
    );
}

/// A boss 6 square and 6 high standing on a plate 20 square and 4 thick.
fn boss(m: &mut Model) -> Shape {
    let plate = block(m, (-10.0, -10.0, 0.0), (20.0, 20.0, 4.0));
    let boss = block(m, (-3.0, -3.0, 4.0), (6.0, 6.0, 6.0));
    fused(m, &[plate, boss])
}

/// The boss's four upright edges and the four edges of its foot.
fn boss_edges(m: &Model, s: &Shape) -> Vec<Shape> {
    edges_of(m, s)
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let at = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
            let (p, q) = (at(&a), at(&b));
            let on_boss = |p: Point| p.x.abs() <= 3.0 + 1e-9 && p.y.abs() <= 3.0 + 1e-9;
            let upright = (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9;
            let foot = (p.z - 4.0).abs() < 1e-9 && (q.z - 4.0).abs() < 1e-9;
            on_boss(p) && on_boss(q) && (upright || foot)
        })
        .collect()
}

/// The boss rounded at `r` on its upright edges and its foot.
///
/// The exact volume, by slices in `z`:
/// - the plate, `1600`, as it is;
/// - above `z = 4 + r` the boss's square with its corners rounded,
///   `(36 − 4(1 − π/4) r²)` over a height of `6 − r`;
/// - below that, at `s = 4 + r − z`, the foot's fill stands
///   `w = r − √(r² − s²)` off every wall, and round each corner the torus's
///   inner face stands `r + w` off the upright band's axis: the slice is
///   the rounded square `Q` grown by `w`, `A_Q + P_Q w + π w²`, with
///   `A_Q = 32 + π` and `P_Q = 16 + 2π` at `r = 1`. Over the slab,
///   `∫ w = 1 − π/4` and `∫ w² = 5/3 − π/2`.
///
/// The total, `1808 + 17π/3 − π²`, is also the blocks less the upright
/// spandrels (`4 × 5` long) plus the foot's (`4 × 4` long, between the
/// tori), each `(1 − π/4)·r²·L`, and four elbows each taking
/// `r³(1 − 7π/6 + π²/4)`: the stairs' elbow with its sign turned, as
/// rounding the complement turns it.
fn boss_volume() -> f64 {
    let slab = (32.0 + PI) + (16.0 + 2.0 * PI) * (1.0 - PI / 4.0) + PI * (5.0 / 3.0 - PI / 2.0);
    1600.0 + 5.0 * (32.0 + PI) + slab
}

#[test]
fn a_boss_s_foot_rolls_round_its_upright_edges_on_tori() {
    let mut m = Model::new();
    let solid = boss(&mut m);
    let edges = boss_edges(&m, &solid);
    assert_eq!(edges.len(), 8);
    let built = fillet_edges_corner(&mut m, &solid, &edges, 1.0, FilletCorner::Rolling, T).unwrap();
    let rolled = &built.shape;
    assert!(ogeom_algo::check(&m, rolled, T).unwrap().is_valid());
    // The eleven faces, a band per edge and a torus per corner of the foot;
    // the upright bands end flush in the boss's top.
    assert_eq!(kinds(&m, rolled), [11, 8, 0, 4, 0]);
    let (_, merged) = ogeom_heal::unify_same_domain(&mut m, rolled, T).unwrap();
    assert_eq!(merged, 0, "faces left split on one surface");
    // Every edge but the plate's, the boss's top and the upright bands'
    // flush ends there is a tangent join.
    let (angle, gap) = worst_join_below(&m, rolled, 9.5);
    assert!(angle <= 1e-6, "a join turns {angle:.2e} rad");
    assert!(gap <= 1e-6, "an edge stands {gap:.2e} off its faces");
    let v = volume(&m, rolled);
    let exact = boss_volume();
    assert!((exact - (1808.0 + 17.0 * PI / 3.0 - PI * PI)).abs() < 1e-9);
    assert!(
        (v - exact).abs() <= 1e-6 * exact,
        "volume {v}, exact {exact}"
    );
}

#[test]
fn a_boss_s_foot_rolls_alike_in_any_edge_order() {
    let build = |order: &dyn Fn(Vec<Shape>) -> Vec<Shape>| {
        let mut m = Model::new();
        let solid = boss(&mut m);
        let edges = order(boss_edges(&m, &solid));
        let built =
            fillet_edges_corner(&mut m, &solid, &edges, 1.0, FilletCorner::Rolling, T).unwrap();
        signature(&m, &built.shape)
    };
    let forward = build(&|e| e);
    let backward = build(&|mut e| {
        e.reverse();
        e
    });
    let shuffled = build(&|e| {
        let (odd, even): (Vec<_>, Vec<_>) =
            e.into_iter().enumerate().partition(|(i, _)| i % 2 == 1);
        odd.into_iter().chain(even).map(|(_, s)| s).collect()
    });
    assert_eq!(forward, backward);
    assert_eq!(forward, shuffled);
}
