//! Rounds where a flat wall runs tangent into a curved one, the face
//! between them narrowing to a cusp at the tangent point.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The edge of `s` passing nearest `q`, 64 samples along each.
fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(m, s, deflection, T)
        .unwrap()
        .mass
}

/// Every vertex's point, sorted, as bits.
fn points(m: &Model, s: &Shape) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    out
}

/// A 40 x 26 x 5 plate, x from -20 to 20 and y from -6 to 20, fused with a
/// boss of radius 6 about the z axis, 10 high: the plate's side at y = -6
/// runs tangent into the boss's wall at (0, -6).
fn tangent_boss(m: &mut Model) -> Shape {
    let plate = ogeom_algo::make_box(m, at(-20.0, -6.0, 0.0), (40.0, 26.0, 5.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 0.0), 6.0, 10.0, T)
        .unwrap()
        .shape;
    ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape
}

/// How far past the plate's side `y = -6` the result stands above the
/// plate's top, outside the boss: read on every face's mesh, each point
/// classified against the side's plane. Nought where nothing hangs over.
fn overhang(m: &Model, s: &Shape) -> f64 {
    let mut worst = 0.0_f64;
    for face in explore_unique(m, s, ShapeType::Face).unwrap() {
        let mesh = ogeom_mesh::triangulate_face(
            m,
            &face,
            ogeom_mesh::Deflection {
                chord: 1e-3,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap();
        for p in &mesh.positions {
            if p.z >= 5.0 - 1e-9 && p.x.hypot(p.y) > 6.0 + 1e-6 {
                worst = worst.max(-6.0 - p.y);
            }
        }
    }
    worst
}

/// Every edge of the boss's foot, the circle of radius 6 at `z = 5`.
fn foot_arcs(m: &Model, s: &Shape) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return false;
            };
            [a, 0.5 * (a + b), b].iter().all(|t| {
                let p = c.point_at(*t, T).unwrap();
                (p.z - 5.0).abs() < 1e-9 && (p.x.hypot(p.y) - 6.0).abs() < 1e-9
            })
        })
        .collect()
}

/// The volume of the foot's concave band of radius 1 about the boss's axis
/// standing past the plane `y = -6`, over the whole turn: the band's
/// section (the corner square less the quarter disc, of height
/// `1 - sqrt(1 - (7 - rho)^2)` at `rho`) swept over the angles where
/// `rho sin(phi) < -6`, `pi - 2 asin(6 / rho)` of them, at radius `rho`.
fn overhang_volume() -> f64 {
    let n = 2_000_000_u32;
    let mut sum = 0.0;
    for k in 0..n {
        let rho = 6.0 + (f64::from(k) + 0.5) / f64::from(n);
        let height = 1.0 - (1.0 - (7.0 - rho) * (7.0 - rho)).sqrt();
        sum += rho * height * (std::f64::consts::PI - 2.0 * (6.0 / rho).asin());
    }
    sum / f64::from(n)
}

/// The boss's foot alone, rounded at 1. Within about 3.5 of (0, -6, 5) the
/// plate's top between the foot and its side narrows below the band's
/// width, so the band runs past the side; it is cut off at the side's
/// plane, its fill ending flush with the side. In either order, the half
/// foot beside the cusp and the whole foot: valid, no point past the side,
/// and the volume the plate and boss plus the band (Pappus) less its part
/// past the side.
#[test]
fn a_foot_round_on_a_face_narrowing_to_a_cusp_is_cut_off_at_the_face_s_side() {
    use std::f64::consts::PI;
    let left = Point::new(-4.242_640_687, 4.242_640_687, 5.0);
    let area = 1.0 - PI / 4.0;
    let offset = (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let base = 40.0 * 26.0 * 5.0 + PI * 36.0 * 5.0;
    let past = overhang_volume();
    for (whole, turn, share) in [(false, PI, 0.5), (true, 2.0 * PI, 1.0)] {
        let want = base + area * (6.0 + offset) * turn - past * share;
        let mut first: Option<f64> = None;
        for flip in [false, true] {
            let mut m = Model::new();
            let s = tangent_boss(&mut m);
            let mut edges = if whole {
                foot_arcs(&m, &s)
            } else {
                vec![edge_near(&m, &s, left)]
            };
            if flip {
                edges.reverse();
            }
            let rounded = ogeom_fillet::fillet_edges(&mut m, &s, &edges, 1.0, T)
                .unwrap_or_else(|e| panic!("whole {whole}, flip {flip}: {e}"))
                .shape;
            let found = ogeom_algo::check(&m, &rounded, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            let out = overhang(&m, &rounded);
            let got = volume(&m, &rounded);
            eprintln!(
                "whole {whole}, flip {flip}: {out} past the side, volume {got} against {want}"
            );
            assert!(out <= 1e-6, "{out} past the plate's side");
            assert!((got - want).abs() < 2e-4, "{got} against {want}");
            // The foot's arcs are rounded in the order asked, as before the
            // cut; the two orders agree on the solid's volume.
            match first {
                None => first = Some(got),
                Some(held) => assert!(
                    (held - got).abs() < 1e-9 * got,
                    "flip {flip}: {held} against {got}"
                ),
            }
        }
    }
}

/// The foot on the side x < 0 (concave) and the plate's top edge on the
/// side x > 0 (convex), a tangent chain running on through the cusp at
/// (0, -6, 5). The two bands lie either side of the vertex; the foot's is
/// cut off at the plate's side as it is alone. Both orders give the same
/// solid, valid, every point of it inside the plate and boss or in the
/// foot's fill over the plate's top, and the volume the plate and boss
/// plus the foot's band less its part past the side, less the edge's band
/// over its 20.
#[test]
fn a_tangent_chain_through_the_cusp_rounds_either_side_of_it_in_either_order() {
    use std::f64::consts::PI;
    let r = 1.0_f64;
    let area = r * r * (1.0 - PI / 4.0);
    let offset = r * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let want = 40.0 * 26.0 * 5.0 + PI * 36.0 * 5.0 + area * PI * (6.0 + offset)
        - 0.5 * overhang_volume()
        - area * 20.0;
    let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
    for flip in [false, true] {
        let mut m = Model::new();
        let s = tangent_boss(&mut m);
        let foot = edge_near(&m, &s, Point::new(-4.242_640_687, 4.242_640_687, 5.0));
        let side = edge_near(&m, &s, Point::new(10.0, -6.0, 5.0));
        let asked = if flip {
            vec![side, foot]
        } else {
            vec![foot, side]
        };
        let rounded = ogeom_fillet::fillet_edges(&mut m, &s, &asked, r, T)
            .unwrap_or_else(|e| panic!("flip {flip}: {e}"))
            .shape;
        let found = ogeom_algo::check(&m, &rounded, T).unwrap();
        assert!(found.is_valid(), "flip {flip}: {:?}", found.problems);
        let got = volume(&m, &rounded);
        eprintln!("flip {flip}: volume {got} against {want}");
        assert!(
            (got - want).abs() < 2e-4,
            "flip {flip}: {got} against {want}"
        );
        let inside = |p: Point| -> bool {
            let e = 1e-6;
            let plate = p.x.abs() <= 20.0 + e
                && p.y >= -6.0 - e
                && p.y <= 20.0 + e
                && p.z >= -e
                && p.z <= 5.0 + e;
            let rho = p.x.hypot(p.y);
            let boss = rho <= 6.0 + e && p.z >= -e && p.z <= 10.0 + e;
            let fill = p.x <= e
                && p.y >= -6.0 - e
                && p.z >= 5.0 - e
                && p.z <= 5.0 + r + e
                && rho >= 6.0 - e
                && rho <= 6.0 + r + e
                && (rho - 6.0 - r).hypot(p.z - 5.0 - r) >= r - 1e-4;
            plate || boss || fill
        };
        for face in explore_unique(&m, &rounded, ShapeType::Face).unwrap() {
            let mesh = ogeom_mesh::triangulate_face(
                &m,
                &face,
                ogeom_mesh::Deflection {
                    chord: 1e-3,
                    ..ogeom_mesh::Deflection::default()
                },
                T,
            )
            .unwrap();
            for p in &mesh.positions {
                assert!(
                    inside(*p),
                    "flip {flip}: outside the solid and the fill: {p:?}"
                );
            }
        }
        let key = (got.to_bits(), points(&m, &rounded));
        match &first {
            None => first = Some(key),
            Some(held) => assert_eq!(held, &key, "flip {flip}: another solid"),
        }
    }
}

/// printCAD's corpus case: a 40 x 40 x 5 plate under a boss of radius 6 at
/// (14, 0), 8 high, its wall tangent to the plate's side x = 20 at
/// (20, 0, 5).
fn printcad_boss(m: &mut Model) -> Shape {
    let plate = ogeom_algo::make_box(m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 5.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(m, at(14.0, 0.0, 5.0), 6.0, 8.0, T)
        .unwrap()
        .shape;
    ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape
}

/// The boss's foot and the plate's top edge at x = 20: two arcs of the
/// foot, two halves of the edge, one tangent chain.
fn printcad_chain(m: &Model, s: &Shape) -> Vec<Shape> {
    let mut chain = Vec::new();
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
        let mid = c.point_at(0.5 * (a + b), T).unwrap();
        let on_top = [p, q, mid].iter().all(|x| (x.z - 5.0).abs() < 1e-6);
        let foot = ((mid.x - 14.0).hypot(mid.y) - 6.0).abs() < 1e-6;
        let side = [p, q, mid].iter().all(|x| (x.x - 20.0).abs() < 1e-6);
        if on_top && (foot || side) {
            chain.push(e);
        }
    }
    assert_eq!(
        chain.len(),
        4,
        "two arcs of the foot, two halves of the edge"
    );
    chain
}

/// How far `(x, z)` lies from the arc about `(cx, cz)` of radius `r` over
/// the angles `a0..a1` (from +x towards +z).
fn arc_distance(x: f64, z: f64, (cx, cz): (f64, f64), r: f64, (a0, a1): (f64, f64)) -> f64 {
    let mut a = (z - cz).atan2(x - cx);
    if a < a0 - 1e-12 {
        a += std::f64::consts::TAU;
    }
    if a <= a1 + 1e-12 {
        return ((x - cx).hypot(z - cz) - r).abs();
    }
    let end = |t: f64| (x - cx - r * t.cos()).hypot(z - cz - r * t.sin());
    end(a0).min(end(a1))
}

/// How far along y from the tangent point the top between the plate's
/// side and the boss's foot is narrower than `2 r`: there `6 - sqrt(36 -
/// y^2) = 2 r`, at `y = 2 sqrt(r (6 - r))`.
fn printcad_s_reach(r: f64) -> f64 {
    2.0 * (r * (6.0 - r)).sqrt()
}

/// How far `p` lies from the round printCAD's chain asks for. Nearer the
/// tangent point than `printcad_s_reach`, in the section at `y`, the step
/// between the side and the wall (a top `w = 6 - sqrt(36 - y^2)` wide) is
/// rounded by an S: a convex arc of radius `r` off the side turning
/// `theta`, `cos theta = 1 - w/2r`, and a concave one of radius `r` onto
/// the wall, point symmetric about `(20 - w/2, 5)`. Past it, the two bands
/// on their own: the edge's quarter round, the top, the foot's torus.
fn printcad_deviation(p: Point, r: f64) -> f64 {
    use std::f64::consts::PI;
    if p.y.abs() < printcad_s_reach(r) {
        let w = 6.0 - (36.0 - p.y * p.y).sqrt();
        let theta = (1.0 - w / (2.0 * r)).acos();
        let rise = r * theta.sin();
        let convex = arc_distance(p.x, p.z, (20.0 - r, 5.0 - rise), r, (0.0, theta));
        let concave = arc_distance(p.x, p.z, (20.0 - w + r, 5.0 + rise), r, (PI, PI + theta));
        return convex.min(concave);
    }
    let rho = (p.x - 14.0).hypot(p.y);
    let edge = arc_distance(p.x, p.z, (20.0 - r, 5.0 - r), r, (0.0, 0.5 * PI));
    let foot = arc_distance(rho, p.z, (6.0 + r, 5.0 + r), r, (PI, 1.5 * PI));
    let top = if p.x <= 20.0 - r && rho >= 6.0 + r {
        (p.z - 5.0).abs()
    } else {
        f64::INFINITY
    };
    edge.min(foot).min(top)
}

/// The volume of the round `printcad_deviation` measures from, by sections
/// across y. Past the S's reach: the plate and boss, less the edge's
/// quarter round over its length there, plus the foot's torus band (its
/// section, of height `r - sqrt(r^2 - (6 + r - rho)^2)` at `rho`, swept
/// over the angles at `rho` outside the slab). Within it the S takes as
/// much from the plate's corner as it adds at the wall's foot, being point
/// symmetric with the step it rounds: the plate and boss alone.
fn printcad_volume(r: f64) -> f64 {
    use std::f64::consts::PI;
    let reach = printcad_s_reach(r);
    let n = 2_000_000_u32;
    let mut band = 0.0;
    for k in 0..n {
        let rho = 6.0 + r * (f64::from(k) + 0.5) / f64::from(n);
        let height = r - (r * r - (6.0 + r - rho).powi(2)).sqrt();
        band += height * rho * (2.0 * PI - 2.0 * (reach / rho).asin());
    }
    band *= r / f64::from(n);
    40.0 * 40.0 * 5.0 + PI * 36.0 * 8.0 - (1.0 - PI / 4.0) * r * r * (40.0 - 2.0 * reach) + band
}

/// A solid's volume bits and sorted vertex bits, to compare orders by.
type SolidKey = (u64, Vec<[u64; 3]>);

/// What a cusp round is checked against.
struct CuspCase<'a> {
    /// How many fitted faces close it: on each side of the tangent point
    /// rounded, the S's strip and a patch either side of it.
    patches: usize,
    /// Carries a point into printCAD's frame, where the side is x = 20,
    /// the wall's axis (14, 0) and the top z = 5, for `printcad_deviation`.
    frame: &'a dyn Fn(Point) -> Point,
    /// Whether a point lies inside the solid asked on or the foot's fill.
    inside: &'a dyn Fn(Point) -> bool,
    /// Whether a point stands over the side near the tangent point, where
    /// no face may look down.
    over_side: &'a dyn Fn(Point) -> bool,
    /// The volume of the round `printcad_deviation` measures from.
    want: f64,
}

/// Asserts every check of a cusp round at `r` on `rounded`: valid; its
/// patches tangent to every neighbour within 2e-4 rad; every point of
/// every face inside the solid or the foot's fill; no face looking down
/// over the side; and the volume within the patches' own distance from the
/// S (each mesh triangle's area times its corners' mean distance, summed)
/// of the case's own. Gives the solid's volume bits and vertices.
fn assert_cusp_round(
    m: &Model,
    rounded: &Shape,
    r: f64,
    case: &CuspCase<'_>,
    label: &str,
) -> SolidKey {
    let found = ogeom_algo::check(m, rounded, T).unwrap();
    if !found.is_valid() {
        for pr in &found.problems {
            let f = &pr.at;
            if let Some(d) = m.node(f).and_then(|n| n.data().as_face()) {
                let kind = format!("{:?}", m.geometry().surface(d.surface).unwrap());
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for v in explore_unique(m, f, ShapeType::Vertex).unwrap() {
                    let q = m.node(&v).unwrap().data().as_vertex().unwrap().point;
                    for (i, c) in [q.x, q.y, q.z].iter().enumerate() {
                        lo[i] = lo[i].min(*c);
                        hi[i] = hi[i].max(*c);
                    }
                }
                eprintln!("BADFACE {:?} {lo:.3?}..{hi:.3?} {}", f.orientation(), &kind[..kind.len().min(40)]);
            }
        }
    }
    assert!(found.is_valid(), "{label}: {:?}", found.problems);
    let patches: Vec<Shape> = explore_unique(m, rounded, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let d = m.node(f).unwrap().data().as_face().unwrap();
            matches!(
                m.geometry().surface(d.surface),
                Some(ogeom_geom::SurfaceGeometry::BSpline(_))
            )
        })
        .collect();
    if std::env::var_os("DOWNWHERE").is_some() {
        for f in &patches {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for v in explore_unique(m, f, ShapeType::Vertex).unwrap() {
                let q = m.node(&v).unwrap().data().as_vertex().unwrap().point;
                for (i, c) in [q.x, q.y, q.z].iter().enumerate() {
                    lo[i] = lo[i].min(*c);
                    hi[i] = hi[i].max(*c);
                }
            }
            let d = m.node(f).unwrap().data().as_face().unwrap();
            let kind = format!("{:?}", m.geometry().surface(d.surface).unwrap());
            let edges = explore_unique(m, f, ShapeType::Edge).unwrap();
            eprintln!("FITTED face box {lo:.3?} .. {hi:.3?} edges {} surface {:?}", edges.len(), d.surface);
            let _ = &kind;
            if edges.len() == 3 && hi[1] - lo[1] < 0.1 {
                for e in &edges {
                    match ogeom_algo::placed_edge_curve(m, e, T).unwrap() {
                        Some((c, (a, b))) => eprintln!("  SLIVER edge {:?} -> {:?} mid {:?}", c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap(), c.point_at(0.5*(a+b), T).unwrap()),
                        None => eprintln!("  SLIVER edge without curve"),
                    }
                }
            }
        }
    }
    let mut tangent = 0.0_f64;
    for patch in &patches {
        for c in ogeom_fillet::analyse_blend(m, rounded, patch, 17, T).unwrap() {
            tangent = tangent.max(c.tangency_error);
        }
    }
    let (mut farthest, mut departure, mut down) = (0.0_f64, 0.0_f64, 1.0_f64);
    for face in explore_unique(m, rounded, ShapeType::Face).unwrap() {
        let mesh = ogeom_mesh::triangulate_face(
            m,
            &face,
            ogeom_mesh::Deflection {
                chord: 1e-3,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap();
        for p in &mesh.positions {
            assert!(
                (case.inside)(*p),
                "{label}: outside the solid and the fill: {p:?}"
            );
        }
        for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
            if (case.over_side)(*p) {
                if n.z < -2e-4 && std::env::var_os("DOWNWHERE").is_some() {
                    let d = m.node(&face).unwrap().data().as_face().unwrap();
                    eprintln!("DOWN {:.2e} at {:.4} {:.4} {:.4} on {:?}", n.z, p.x, p.y, p.z, std::mem::discriminant(m.geometry().surface(d.surface).unwrap()));
                }
                down = down.min(n.z);
            }
        }
        if !patches.iter().any(|q| q.node() == face.node()) {
            continue;
        }
        let off: Vec<f64> = mesh
            .positions
            .iter()
            .map(|p| printcad_deviation((case.frame)(*p), r))
            .collect();
        farthest = off.iter().copied().fold(farthest, f64::max);
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let mean = t.iter().map(|&i| off[i as usize]).sum::<f64>() / 3.0;
            departure += 0.5 * (b - a).cross(c - a).magnitude() * mean;
        }
    }
    let got = volume(m, rounded);
    eprintln!(
        "{label}: {} patches tangent within {tangent:e}, lowest normal z {down:e}, at most \
         {farthest:e} from the S, {departure:e} by area; volume {got} against {} ({:e})",
        patches.len(),
        case.want,
        got - case.want
    );
    if std::env::var_os("DOWNWHERE").is_none() {
        assert_eq!(patches.len(), case.patches, "{label}: fitted faces");
    }
    assert!(tangent < 2e-4, "{label}: tangent within {tangent}");
    assert!(down > -5e-3, "{label}: a face looks down, {down}");
    assert!(
        (got - case.want).abs() <= departure + 1e-4,
        "{label}: {got} against {}, past the patches' {departure} from the S",
        case.want
    );
    (got.to_bits(), points(m, rounded))
}

/// printCAD's chain rounded at `r`, its edges asked in reverse if `flip`,
/// every check of `assert_cusp_round` asserted where it builds, or why not.
fn printcad_round(r: f64, flip: bool) -> Result<SolidKey, ogeom_core::OgeomError> {
    let mut m = Model::new();
    let s = printcad_boss(&mut m);
    let mut chain = printcad_chain(&m, &s);
    if flip {
        chain.reverse();
    }
    let rounded = ogeom_fillet::fillet_edges(&mut m, &s, &chain, r, T)?.shape;
    let inside = |p: Point| -> bool {
        let e = 1e-6;
        let plate = p.x.abs() <= 20.0 + e && p.y.abs() <= 20.0 + e && p.z >= -e && p.z <= 5.0 + e;
        let rho = (p.x - 14.0).hypot(p.y);
        let boss = rho <= 6.0 + e && p.z >= 5.0 - e && p.z <= 13.0 + e;
        let fill = p.x <= 20.0 + e
            && p.z >= 5.0 - e
            && p.z <= 5.0 + r + e
            && rho >= 6.0 - e
            && rho <= 6.0 + r + e
            && (rho - 6.0 - r).hypot(p.z - 5.0 - r) >= r - 1e-4;
        plate || boss || fill
    };
    let over_side = |p: Point| p.x > 20.0 - 2.0 * r && p.z > 5.0 - r && p.y.abs() < 6.0;
    let case = CuspCase {
        patches: 6,
        frame: &|p| p,
        inside: &inside,
        over_side: &over_side,
        want: printcad_volume(r),
    };
    Ok(assert_cusp_round(
        &m,
        &rounded,
        r,
        &case,
        &format!("printCAD r {r}, flip {flip}"),
    ))
}

/// printCAD's chain at its own radius 1 rounded: the edge's round turns
/// into the foot's through an S either side of (20, 0, 5), passing every
/// check of `assert_cusp_round`, one solid to the bit in either order.
/// Before 0.18.0 the foot's band ran on past the plate's side (material at
/// (20.95, 0, 5.001), outside the solid); 0.18.0 and 0.18.1 refuse it.
#[test]
fn a_foot_and_an_edge_through_the_cusp_where_a_boss_runs_tangent_into_a_plate_turn_through_an_s() {
    let first = printcad_round(1.0, false).unwrap_or_else(|e| panic!("in order: {e}"));
    let second = printcad_round(1.0, true).unwrap_or_else(|e| panic!("reversed: {e}"));
    assert!(first == second, "the two orders give two solids");
}

/// printCAD's chain at other radii: each order either builds a solid that
/// passes every check of `assert_cusp_round`, or is refused by name; the
/// two orders agree. A quarter of the boss's radius and up is refused for
/// leaving the S no room.
#[test]
fn printcad_s_chain_at_other_radii_closes_soundly_or_is_refused_by_name() {
    for r in [0.25, 0.5, 0.75, 1.25, 1.5] {
        let got: Vec<Result<SolidKey, String>> = [false, true]
            .iter()
            .map(|&flip| match printcad_round(r, flip) {
                Ok(key) => Ok(key),
                Err(
                    e @ (ogeom_core::OgeomError::NotDone(_)
                    | ogeom_core::OgeomError::Construction(_)),
                ) => Err(e.to_string()),
                Err(e) => panic!("r {r}, flip {flip}: refused without a name: {e}"),
            })
            .collect();
        eprintln!(
            "r {r}: {}",
            if got[0].is_ok() { "built" } else { "refused" }
        );
        match (&got[0], &got[1]) {
            (Ok(a), Ok(b)) => assert!(a == b, "r {r}: the two orders give two solids"),
            (Err(_), Err(_)) => {}
            _ => panic!("r {r}: one order builds, the other is refused: {got:?}"),
        }
    }
}

/// The volume of the round the arm and hub ask at `r`, by sections across
/// x: the arm and the hub less their overlap (the half disc under the
/// arm's top), less the side edge's quarter round past the S's reach, plus
/// the foot's torus band over its quarter turn from (6, 0) to (0, 6), save
/// where it stands nearer the hub's side of the S's reach (`x <
/// printcad_s_reach`). Within that reach the S adds nothing, as in
/// `printcad_volume`.
fn hub_volume(r: f64) -> f64 {
    use std::f64::consts::PI;
    let reach = printcad_s_reach(r);
    let n = 2_000_000_u32;
    let mut band = 0.0;
    for k in 0..n {
        let rho = 6.0 + r * (f64::from(k) + 0.5) / f64::from(n);
        let height = r - (r * r - (6.0 + r - rho).powi(2)).sqrt();
        band += height * rho * (reach / rho).acos();
    }
    band *= r / f64::from(n);
    40.0 * 12.0 * 5.0 + 270.0 * PI - (1.0 - PI / 4.0) * r * r * (40.0 - reach) + band
}

/// The issue's arm and hub rounded at `r`, its edges asked in reverse if
/// `flip`, every check of `assert_cusp_round` asserted where it builds.
fn hub_round(r: f64, flip: bool) -> Result<SolidKey, ogeom_core::OgeomError> {
    let mut m = Model::new();
    let arm = ogeom_algo::make_box(&mut m, at(0.0, -6.0, 0.0), (40.0, 12.0, 5.0), T)
        .unwrap()
        .shape;
    let hub = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 0.0), 6.0, 10.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &arm, &hub, T).unwrap().shape;
    let mut edges = vec![
        edge_near(&m, &s, Point::new(20.0, 6.0, 5.0)),
        edge_near(&m, &s, Point::new(4.242_640_687, 4.242_640_687, 5.0)),
    ];
    if flip {
        edges.reverse();
    }
    let rounded = ogeom_fillet::fillet_edges(&mut m, &s, &edges, r, T)?.shape;
    let inside = |p: Point| -> bool {
        let e = 1e-6;
        let arm =
            p.x >= -e && p.x <= 40.0 + e && p.y.abs() <= 6.0 + e && p.z >= -e && p.z <= 5.0 + e;
        let rho = p.x.hypot(p.y);
        let hub = rho <= 6.0 + e && p.z >= -e && p.z <= 10.0 + e;
        let fill = p.x >= -e
            && p.y >= -e
            && p.y <= 6.0 + e
            && p.z >= 5.0 - e
            && p.z <= 5.0 + r + e
            && rho >= 6.0 - e
            && rho <= 6.0 + r + e
            && (rho - 6.0 - r).hypot(p.z - 5.0 - r) >= r - 1e-4;
        arm || hub || fill
    };
    let over_side = |p: Point| p.y > 6.0 - 2.0 * r && p.z > 5.0 - r && p.x < 6.0;
    let case = CuspCase {
        patches: 3,
        frame: &|p: Point| Point::new(14.0 + p.y, p.x, p.z),
        inside: &inside,
        over_side: &over_side,
        want: hub_volume(r),
    };
    Ok(assert_cusp_round(
        &m,
        &rounded,
        r,
        &case,
        &format!("hub r {r}, flip {flip}"),
    ))
}

/// The issue's arm and hub: an arm x 0..40, y -6..6, z 0..5 fused with a
/// hub of radius 6 about z, 10 high, the arm's sides tangent to the hub's
/// wall. Its top side edge (convex) and the arc where its top meets the hub
/// (concave) leave (0, 6, 5) the same way, the arm's top between them a
/// cusp: at radius 1 the side edge's round turns into the foot's through
/// an S on that one side of the tangent point, passing every check of
/// `assert_cusp_round`, one solid to the bit in either order. 0.17.1 to
/// 0.18.1 refuse it.
#[test]
fn the_arm_s_side_edge_and_its_arc_into_the_hub_turn_through_an_s_in_either_order() {
    let first = hub_round(1.0, false).unwrap_or_else(|e| panic!("in order: {e}"));
    let second = hub_round(1.0, true).unwrap_or_else(|e| panic!("reversed: {e}"));
    assert!(first == second, "the two orders give two solids");
}

/// printCAD's plate with one side of the tangent point asked only: the
/// foot's arc and the plate's top edge on the side y > 0, the edges past
/// the tangent point left sharp. A one-sided S beside sharp edges is not
/// closed yet: refused by name in either order, not built.
#[test]
fn printcad_s_chain_asked_on_one_side_of_the_tangent_point_is_refused_by_name() {
    for flip in [false, true] {
        let mut m = Model::new();
        let s = printcad_boss(&mut m);
        let mut chain: Vec<Shape> = printcad_chain(&m, &s)
            .into_iter()
            .filter(|e| {
                let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
                c.point_at(0.5 * (a + b), T).unwrap().y > 0.0
            })
            .collect();
        assert_eq!(chain.len(), 2, "an arc of the foot and half the edge");
        if flip {
            chain.reverse();
        }
        match ogeom_fillet::fillet_edges(&mut m, &s, &chain, 1.0, T) {
            Err(e @ ogeom_core::OgeomError::Construction(_)) => {
                let why = e.to_string();
                assert!(why.contains("did not close"), "flip {flip}: {why}");
            }
            Err(e) => panic!("flip {flip}: refused otherwise: {e}"),
            Ok(_) => panic!("flip {flip}: built"),
        }
    }
}

#[test]
fn zz_probe() {
    let Ok(rs) = std::env::var("CUSPR") else { return };
    for r in rs.split(',').map(|x| x.parse::<f64>().unwrap()) {
        for (name, f) in [("pc", printcad_round as fn(f64, bool) -> _), ("hub", hub_round)] {
            match f(r, false) {
                Ok(_) => eprintln!("PROBE {name} r {r}: built"),
                Err(e) => eprintln!("PROBE {name} r {r}: {e}"),
            }
        }
    }
}
