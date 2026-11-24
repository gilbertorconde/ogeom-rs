//! Overflow handling for chamfers (`chamfer_edges_overflow`): a bevel
//! reaching geometry other than its own two faces. Blend into round runs a
//! boss's foot bevel on over the straight round it reaches; trim through is
//! the construction the chamfer built before the option set; with both off
//! the request is refused by name. Hold edge does not apply to a bevel.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};

use ogeom_core::Tolerances;
use ogeom_fillet::{
    Chamfer, ChamferCorner, Overflow, chamfer_edges, chamfer_edges_overflow, fillet_edges,
};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The block's size, its edge round's radius, and the boss's radius and
/// height.
const BLOCK: (f64, f64, f64) = (40.0, 30.0, 15.0);
const EDGE_ROUND: f64 = 5.0;
const BOSS: f64 = 6.0;
const BOSS_HIGH: f64 = 15.0;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [0.0, 0.3, 0.7, 1.0]
                .iter()
                .all(|t| on(c.point_at(a + (b - a) * t, T).unwrap()))
        })
        .collect()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn surface_of(m: &Model, face: &Shape) -> SurfaceGeometry {
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().clone()
}

fn kind(g: &SurfaceGeometry) -> &'static str {
    match g {
        SurfaceGeometry::Plane(_) => "plane",
        SurfaceGeometry::Cylinder(_) => "cylinder",
        SurfaceGeometry::Cone(_) => "cone",
        SurfaceGeometry::BSpline(_) => "fitted",
        _ => "other",
    }
}

/// Each face's surface kind and area, sorted.
fn faces(m: &Model, s: &Shape) -> Vec<(&'static str, f64)> {
    let mut all: Vec<(&'static str, f64)> = explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let area = ogeom_algo::surface_properties(m, f, ogeom_mesh::Deflection::default(), T)
                .unwrap()
                .mass;
            (kind(&surface_of(m, f)), area)
        })
        .collect();
    all.sort_by(|a, b| a.0.cmp(b.0).then(a.1.total_cmp(&b.1)));
    all
}

/// The block 40 by 30 by 15, its top edge along `y = 0` rounded at 5, and
/// a boss of radius 6, 15 tall, standing on the top with its axis at
/// `(20, boss_y)`; with the boss's foot rim and the top face.
fn bossed(m: &mut Model, boss_y: f64) -> (Shape, Vec<Shape>, Shape) {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), BLOCK, T)
        .unwrap()
        .shape;
    let top = edges_where(m, &block, |p| {
        p.y.abs() < 1e-9 && (p.z - BLOCK.2).abs() < 1e-9
    });
    let rounded = fillet_edges(m, &block, &top, EDGE_ROUND, T).unwrap().shape;
    let boss = ogeom_algo::make_cylinder(m, at(20.0, boss_y, BLOCK.2), BOSS, BOSS_HIGH, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(m, &rounded, &boss, T).unwrap().shape;
    let foot = edges_where(m, &s, |p| {
        (p.z - BLOCK.2).abs() < 1e-9 && ((p.x - 20.0).hypot(p.y - boss_y) - BOSS).abs() < 1e-9
    });
    assert!(!foot.is_empty(), "the foot rim");
    let floor = ancestors_of(m, &s, &foot[0], ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface_of(m, f), SurfaceGeometry::Plane(_)))
        .unwrap();
    (s, foot, floor)
}

/// The independent reference: the fill a foot bevel `out` along the floor
/// and `rise` up the wall adds, by its definition, integrated without the
/// construction. At each angle about the boss's axis the bevel is the
/// straight line in the plane through the axis from the wall at `rise`
/// over the top to the floor straight below the circle `out` from the
/// wall, on the top or on the edge round past its tangent line `y = 5`. A
/// column's height is the line's over the floor there; the columns are
/// integrated by Gauss-Legendre quadrature over the ring, split where the
/// leg on the floor leaves the top and where a column crosses the tangent
/// line.
struct Reference {
    boss_y: f64,
    out: f64,
    rise: f64,
}

impl Reference {
    /// The floor's height at `y`: the top, or the edge round past its
    /// tangent line.
    fn floor(y: f64) -> f64 {
        if y >= EDGE_ROUND {
            BLOCK.2
        } else {
            let d = EDGE_ROUND - y;
            (BLOCK.2 - EDGE_ROUND) + (EDGE_ROUND * EDGE_ROUND - d * d).max(0.0).sqrt()
        }
    }

    /// The column at `rho` from the axis, at angle `a` from `-y`.
    fn column(&self, rho: f64, a: f64) -> f64 {
        let y_at = |r: f64| self.boss_y - r * a.cos();
        let leg = Self::floor(y_at(BOSS + self.out));
        let top = BLOCK.2 + self.rise;
        let line = top + (leg - top) * (rho - BOSS) / self.out;
        (line - Self::floor(y_at(rho))).max(0.0)
    }

    fn volume(&self, n: usize) -> f64 {
        let g = gauss(n);
        let reach = BOSS + self.out;
        let gap = self.boss_y - EDGE_ROUND;
        let half = if gap < reach {
            (gap / reach).acos()
        } else {
            0.0
        };
        let mut sum = 0.0;
        for (a0, a1) in [(-half, half), (half, TAU - half)] {
            if a1 - a0 <= 0.0 {
                continue;
            }
            let ha = 0.5 * (a1 - a0);
            for (ga, wa) in &g {
                let a = a0 + ha * (ga + 1.0);
                // Where the column crosses the tangent line, if it does.
                let cross = gap / a.cos();
                let mut cuts = vec![BOSS];
                if a.cos() > 0.0 && cross > BOSS && cross < reach {
                    cuts.push(cross);
                }
                cuts.push(reach);
                for w in cuts.windows(2) {
                    let hr = 0.5 * (w[1] - w[0]);
                    for (gr, wr) in &g {
                        let rho = w[0] + hr * (gr + 1.0);
                        sum += wa * wr * self.column(rho, a) * rho * hr * ha;
                    }
                }
            }
        }
        sum
    }
}

/// Gauss-Legendre nodes and weights on [-1, 1], by Newton on the
/// polynomial's roots.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
    let nf = n as f64;
    (0..n)
        .map(|i| {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
            let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
            let mut d = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
                    let kf = k as f64;
                    let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                    p0 = p1;
                    p1 = p2;
                }
                d = nf * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / d;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            (x, 2.0 / ((1.0 - x * x) * d * d))
        })
        .collect()
}

/// The largest angle between the normals either side of each edge where
/// the fitted piece meets the cone, sampled along the edge's own curve and
/// projected onto each face's surface; with the number of such edges.
fn worst_transition(m: &Model, s: &Shape) -> (f64, usize) {
    let mut worst = 0.0_f64;
    let mut count = 0;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(m, a), surface_of(m, b));
        let mut kinds = [kind(&sa), kind(&sb)];
        kinds.sort_unstable();
        if kinds != ["cone", "fitted"] {
            continue;
        }
        count += 1;
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 1..20 {
            let p = c.point_at(lo + (hi - lo) * f64::from(k) / 20.0, T).unwrap();
            let normal = |g: &SurfaceGeometry| -> Vector {
                let foot = ogeom_geom::project::project_on_surface(g, p, 16, T).unwrap();
                let (u, v) = foot.parameters;
                g.normal_at(u, v, T).unwrap().vector()
            };
            let angle = normal(&sa).dot(normal(&sb)).abs().clamp(0.0, 1.0).acos();
            worst = worst.max(angle);
        }
    }
    (worst, count)
}

/// Whether every edge between the fitted piece and the edge round lies on
/// the edge round (there is at least one), and no planar face of the
/// result reaches past the tangent line at the top's height: the bevel
/// meets the round itself, with no ledge left over it.
fn meets_the_round_without_a_step(m: &Model, s: &Shape) -> bool {
    let mut met = 0;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let kinds: Vec<&str> = faces.iter().map(|f| kind(&surface_of(m, f))).collect();
        if !(kinds.contains(&"fitted") && kinds.contains(&"cylinder")) {
            continue;
        }
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        let on_round = (0..=10).all(|k| {
            let p = c.point_at(lo + (hi - lo) * f64::from(k) / 10.0, T).unwrap();
            ((p.y - EDGE_ROUND).hypot(p.z - (BLOCK.2 - EDGE_ROUND)) - EDGE_ROUND).abs() < 1e-6
        });
        if on_round && c.point_at(0.5 * (lo + hi), T).unwrap().y < EDGE_ROUND {
            met += 1;
        }
    }
    let ledge = explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .iter()
        .any(|f| {
            let bb = ogeom_algo::shape_bounds(m, f, T).unwrap();
            let (l, h) = (bb.low().unwrap(), bb.high().unwrap());
            matches!(surface_of(m, f), SurfaceGeometry::Plane(_))
                && (l.z - BLOCK.2).abs() < 1e-6
                && (h.z - BLOCK.2).abs() < 1e-6
                && l.y < EDGE_ROUND - 1e-6
        });
    met > 0 && !ledge
}

/// The boss at `y = 12`, its foot bevelled 3 along the top and 3 up the
/// wall: the leg on the top reaches 2 past the edge round's tangent line
/// `y = 5`. Blend into round runs the bevel on over the edge round: the
/// result is valid and closed, the bevel meets the round with no ledge,
/// the fitted piece meets the cone G1 across both transition edges, and the
/// fill matches the reference within its bound. The reference is checked
/// against Pappus first on a boss clear of the edge round (the cone's
/// fill, `2 pi (6 + 1) 4.5`); its own error bound is the change from 48 to
/// 96 points per range. The same holds for the legs 3 along the top and 2
/// up the wall, asked as distances on the top.
#[test]
fn a_boss_s_foot_bevel_blends_into_the_edge_round_it_reaches() {
    let pappus = TAU * (BOSS + 1.0) * 4.5;
    let clear = Reference {
        boss_y: 20.0,
        out: 3.0,
        rise: 3.0,
    }
    .volume(32);
    assert!(
        (clear - pappus).abs() < 1e-9,
        "the reference {clear} against Pappus {pappus}"
    );

    for (out, rise) in [(3.0, 3.0), (3.0, 2.0)] {
        let mut m = Model::new();
        let (s, foot, floor) = bossed(&mut m, 12.0);
        let before = volume(&m, &s);
        let specs: Vec<(Shape, Chamfer)> = foot
            .iter()
            .map(|e| {
                let spec = if out == rise {
                    Chamfer::Symmetric(out)
                } else {
                    Chamfer::Distances {
                        face: floor.clone(),
                        on_face: out,
                        on_other: rise,
                    }
                };
                (e.clone(), spec)
            })
            .collect();
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{out}/{rise}: {:?}", found.problems);
        for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
            assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
        }
        for e in &foot {
            assert!(!built.history.generated(e).is_empty(), "the rim's bevel");
        }
        assert!(
            meets_the_round_without_a_step(&m, &built.shape),
            "{out}/{rise}: the bevel meets the round with no ledge"
        );
        let (crease, transitions) = worst_transition(&m, &built.shape);
        assert_eq!(transitions, 2, "{out}/{rise}: the two transition edges");
        assert!(
            crease < 1e-5,
            "{out}/{rise}: the transition creases by {crease}"
        );

        let reference = Reference {
            boss_y: 12.0,
            out,
            rise,
        };
        let (coarse, fine) = (reference.volume(48), reference.volume(96));
        let bound = (fine - coarse).abs().max(1e-9);
        let got = volume(&m, &built.shape) - before;
        assert!(
            (got - fine).abs() < 1e-6 + bound,
            "{out}/{rise}: the fill {got} against the reference {fine} (within {bound})"
        );
        // Not the cone's fill: the bevel runs down over the edge round.
        let cone = TAU * (BOSS + out / 3.0) * (0.5 * out * rise);
        assert!((got - cone).abs() > 1e-2, "{got} against the cone's {cone}");
    }
}

/// Blend into round off: the request goes on to trim through, and the
/// result is the solid the chamfer built before overflow handling: the
/// whole cone's fill, `2 pi 7 4.5`, its floor ring carried on flat as a
/// ledge over the edge round, the faces those of that solid (recorded from
/// it, each to 1e-6 in area). Both off: refused by name. Hold edge off
/// changes nothing; trim through off alone leaves blend into round to win.
#[test]
fn blend_into_round_off_trims_through_as_before_and_both_off_refuse() {
    let mut m = Model::new();
    let (s, foot, _) = bossed(&mut m, 12.0);
    let before = volume(&m, &s);
    let specs: Vec<(Shape, Chamfer)> = foot
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(3.0)))
        .collect();
    let corner = ChamferCorner::Facet;
    let mut blend_off = Overflow::default();
    blend_off.blend_into_round = false;
    let trimmed = chamfer_edges_overflow(&mut m, &s, &specs, corner, blend_off, T).unwrap();
    assert!(ogeom_algo::check(&m, &trimmed.shape, T).unwrap().is_valid());
    let want = TAU * 7.0 * 4.5;
    let got = volume(&m, &trimmed.shape) - before;
    assert!(
        (got - want).abs() < 1e-6,
        "the cone's fill {got} against {want}"
    );
    let before_overflow = [
        ("cone", 199.929_732_217),
        ("cylinder", 314.159_265_359),
        ("cylinder", 452.389_342_117),
        ("plane", 15.455_599_585),
        ("plane", 113.097_335_529),
        ("plane", 400.0),
        ("plane", 444.634_954_085),
        ("plane", 444.634_954_085),
        ("plane", 600.0),
        ("plane", 760.986_594_644),
        ("plane", 1200.0),
    ];
    let now = faces(&m, &trimmed.shape);
    assert_eq!(now.len(), before_overflow.len(), "{now:?}");
    for ((k, a), (wk, wa)) in now.iter().zip(before_overflow) {
        assert!(*k == wk && (a - wa).abs() < 1e-6, "{now:?}");
    }
    assert!(
        !meets_the_round_without_a_step(&m, &trimmed.shape),
        "the ledge"
    );

    let err = chamfer_edges_overflow(&mut m, &s, &specs, corner, Overflow::none(), T).unwrap_err();
    assert!(err.to_string().contains("blend into round"), "{err}");

    let blended = chamfer_edges(&mut m, &s, &foot, 3.0, T).unwrap();
    let mut hold_off = Overflow::default();
    hold_off.hold_edge = false;
    hold_off.hold_edge_at_end = false;
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    for options in [hold_off, trim_off] {
        let other = chamfer_edges_overflow(&mut m, &s, &specs, corner, options, T).unwrap();
        assert_eq!(
            faces(&m, &other.shape).len(),
            faces(&m, &blended.shape).len()
        );
        assert!((volume(&m, &other.shape) - volume(&m, &blended.shape)).abs() < 1e-9);
    }
}

/// A plate 10 × 10 × 2, its top front edge chamfered 1 on the top and 3
/// down the 2 mm front: the bevel runs past the front and on into the
/// floor as a notch, trimming through it (ADR 0001, amendment "M5b"). Blend
/// into round does not take it, so every option set with trim through on
/// builds the notch, `200 - 10 (4/3)`, face for face; with trim through
/// off it is refused by name.
#[test]
fn a_notch_trims_through_and_trim_through_off_refuses_it() {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (10.0, 10.0, 2.0), T)
        .unwrap()
        .shape;
    let top = explore_unique(&m, &s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| ogeom_algo::face_normal(&m, f, T).unwrap().1.z > 0.5)
        .unwrap();
    let edge = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && (p.z - 2.0).abs() < 1e-9);
    let specs = vec![(
        edge[0].clone(),
        Chamfer::Distances {
            face: top,
            on_face: 1.0,
            on_other: 3.0,
        },
    )];
    let corner = ChamferCorner::Facet;
    let want = 200.0 - 10.0 * 4.0 / 3.0;
    let default =
        chamfer_edges_overflow(&mut m, &s, &specs, corner, Overflow::default(), T).unwrap();
    let mut blend_off = Overflow::default();
    blend_off.blend_into_round = false;
    let trimmed = chamfer_edges_overflow(&mut m, &s, &specs, corner, blend_off, T).unwrap();
    for built in [&default, &trimmed] {
        assert!((volume(&m, &built.shape) - want).abs() <= 1e-9 * want);
    }
    assert_eq!(
        faces(&m, &default.shape).len(),
        faces(&m, &trimmed.shape).len()
    );
    for (a, b) in faces(&m, &default.shape)
        .iter()
        .zip(faces(&m, &trimmed.shape))
    {
        assert!(a.0 == b.0 && (a.1 - b.1).abs() < 1e-9);
    }

    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let err = chamfer_edges_overflow(&mut m, &s, &specs, corner, trim_off, T).unwrap_err();
    assert!(
        err.to_string().contains("trim through is off"),
        "refused by name: {err}"
    );
}

/// The boss at `y = 12` with its foot bevelled `d` on both legs, built with
/// every overflow behaviour on: valid and closed, its fill against the
/// reference within the reference's bound. With `past`, its leg on the
/// floor reaches past the tangent line `y = 5` by more than the booleans'
/// widened edges: the bevel meets the edge round with no ledge, the fitted
/// piece meeting the cone G1 across both transition edges.
fn foot_bevel_builds(d: f64, past: bool) {
    let mut m = Model::new();
    let (s, foot, _) = bossed(&mut m, 12.0);
    let before = volume(&m, &s);
    let specs: Vec<(Shape, Chamfer)> = foot
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(d)))
        .collect();
    let built = chamfer_edges_overflow(
        &mut m,
        &s,
        &specs,
        ChamferCorner::Facet,
        Overflow::default(),
        T,
    )
    .unwrap_or_else(|e| panic!("{d}: {e}"));
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{d}: {:?}", found.problems);
    for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(
            ogeom_algo::is_shell_closed(&m, &shell).unwrap(),
            "{d}: closed"
        );
    }
    let reference = Reference {
        boss_y: 12.0,
        out: d,
        rise: d,
    };
    let (coarse, fine) = (reference.volume(48), reference.volume(96));
    let bound = (fine - coarse).abs().max(1e-9);
    let got = volume(&m, &built.shape) - before;
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "{d}: the fill {got} against the reference {fine} (within {bound})"
    );
    if past {
        assert!(
            meets_the_round_without_a_step(&m, &built.shape),
            "{d}: the bevel meets the round with no ledge"
        );
        let (crease, transitions) = worst_transition(&m, &built.shape);
        assert_eq!(transitions, 2, "{d}: the two transition edges");
        assert!(crease < 1e-5, "{d}: the transition creases by {crease}");
    }
}

/// The foot bevel short of the tangent line (0.9, the cone alone) and past
/// it (1.1 to 2.7, blended into the edge round), where a reading of the
/// largest bevel steps from its first refusal just above 1 straight to
/// 2.71: every size between builds.
#[test]
fn a_foot_bevel_builds_either_side_of_the_tangent_line() {
    for d in [0.9, 1.1, 1.9, 2.7, 2.71] {
        foot_bevel_builds(d, d > 1.0);
    }
}

/// [`a_foot_bevel_builds_either_side_of_the_tangent_line`] at every tenth
/// from 1.1 to 2.6.
#[test]
#[ignore = "heavy: sixteen foot bevels and their references"]
fn a_foot_bevel_builds_at_every_tenth_past_the_tangent_line() {
    for k in 11..=26 {
        foot_bevel_builds(f64::from(k) / 10.0, true);
    }
}

/// The foot bevel 1 on both legs: its leg on the floor, the circle of
/// radius 7 about the boss's axis, touches the tangent line `y = 5` 7 from
/// the axis. The floor and the edge round would be left wedges closing to
/// nothing there, and neither blend into round nor trim through builds
/// them: refused by a name that says so, through `chamfer_edges_overflow`
/// with blend into round on or off and through `chamfer_edges`. Just past
/// the line, where the leg crosses it at under a hundredth of a radian
/// (`d - 1` below about `3.5e-4`), the bevel either builds, valid and
/// matching the reference, or is refused by the same name; a build there
/// may fold the round's sliver under the bevel into an edge widened to it.
#[test]
fn a_foot_bevel_grazing_the_tangent_line_is_refused_by_name() {
    let named = |e: &ogeom_core::OgeomError| {
        let text = e.to_string();
        text.contains("grazes the tangent line of an existing round") && text.contains("chamfer")
    };
    let mut m = Model::new();
    let (s, foot, _) = bossed(&mut m, 12.0);
    let specs: Vec<(Shape, Chamfer)> = foot
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(1.0)))
        .collect();
    let mut blend_off = Overflow::default();
    blend_off.blend_into_round = false;
    for options in [Overflow::default(), blend_off] {
        let err = chamfer_edges_overflow(&mut m, &s, &specs, ChamferCorner::Facet, options, T)
            .unwrap_err();
        assert!(named(&err), "{err}");
    }
    let err = chamfer_edges(&mut m, &s, &foot, 1.0, T).unwrap_err();
    assert!(named(&err), "{err}");

    for past in [1e-6, 1e-5, 1e-4] {
        let d = 1.0 + past;
        let mut m = Model::new();
        let (s, foot, _) = bossed(&mut m, 12.0);
        let specs: Vec<(Shape, Chamfer)> = foot
            .iter()
            .map(|e| (e.clone(), Chamfer::Symmetric(d)))
            .collect();
        match chamfer_edges_overflow(
            &mut m,
            &s,
            &specs,
            ChamferCorner::Facet,
            Overflow::default(),
            T,
        ) {
            Ok(_) => foot_bevel_builds(d, false),
            Err(e) => assert!(named(&e), "{d}: {e}"),
        }
    }
}

/// The notch of [`a_notch_trims_through_and_trim_through_off_refuses_it`]
/// with the three edges at the corner `(10, 10, 0)` bevelled 0.5 in the
/// same request. The corner's bottom edges are seated on the floor the
/// notch trims, but the notch's bevel is not: trim through off still
/// refuses the request by name, under both corners. With trim
/// through on it builds, the notch and the corner taking the plate's
/// material apart.
#[test]
fn trim_through_off_holds_per_blend_when_another_edge_owns_the_face() {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (10.0, 10.0, 2.0), T)
        .unwrap()
        .shape;
    let top = explore_unique(&m, &s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| ogeom_algo::face_normal(&m, f, T).unwrap().1.z > 0.5)
        .unwrap();
    let notch = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && (p.z - 2.0).abs() < 1e-9);
    let mut specs = vec![(
        notch[0].clone(),
        Chamfer::Distances {
            face: top,
            on_face: 1.0,
            on_other: 3.0,
        },
    )];
    let corner = edges_where(&m, &s, |p| {
        let at_x = (p.x - 10.0).abs() < 1e-9;
        let at_y = (p.y - 10.0).abs() < 1e-9;
        let at_z = p.z.abs() < 1e-9;
        (at_x && at_y) || (at_x && at_z) || (at_y && at_z)
    });
    assert_eq!(corner.len(), 3, "the corner's three edges");
    specs.extend(corner.iter().map(|e| (e.clone(), Chamfer::Symmetric(0.5))));
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    for corner in [ChamferCorner::Point, ChamferCorner::Facet] {
        let built = chamfer_edges_overflow(&mut m, &s, &specs, corner, Overflow::default(), T)
            .unwrap_or_else(|e| panic!("{corner:?}: trim through builds it: {e}"));
        assert!(volume(&m, &built.shape) < 200.0 - 10.0 * 4.0 / 3.0 + 1e-9);
        let err = chamfer_edges_overflow(&mut m, &s, &specs, corner, trim_off, T)
            .map(|_| ())
            .unwrap_err();
        assert!(
            err.to_string().contains("trim through is off"),
            "{corner:?}: refused by name: {err}"
        );
    }
}

/// The block 40 by 22 by 15, both its top long edges rounded at 5, and the
/// boss of radius 6 standing on the top at `(20, 11)`; with the boss's
/// foot rim.
fn bossed_between_rounds(m: &mut Model) -> (Shape, Vec<Shape>) {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (40.0, 22.0, 15.0), T)
        .unwrap()
        .shape;
    let top = edges_where(m, &block, |p| {
        (p.y.abs() < 1e-9 || (p.y - 22.0).abs() < 1e-9) && (p.z - 15.0).abs() < 1e-9
    });
    assert_eq!(top.len(), 2);
    let rounded = fillet_edges(m, &block, &top, EDGE_ROUND, T).unwrap().shape;
    let boss = ogeom_algo::make_cylinder(m, at(20.0, 11.0, 15.0), BOSS, BOSS_HIGH, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(m, &rounded, &boss, T).unwrap().shape;
    let foot = edges_where(m, &s, |p| {
        (p.z - 15.0).abs() < 1e-9 && ((p.x - 20.0).hypot(p.y - 11.0) - BOSS).abs() < 1e-9
    });
    assert!(!foot.is_empty(), "the foot rim");
    (s, foot)
}

/// A foot bevel whose floor ring trim through carries on flat as a ledge
/// over an edge round changes no face of the solid, but lays the floor on
/// past its side over nothing: only trim through builds it, so with trim
/// through off it is refused by name. The foot between two rounds
/// (which blend into round does not take) under every option set with
/// trim through off and both corners, and one arc of the one-round foot
/// alone with every option off.
#[test]
fn trim_through_off_refuses_a_foot_ledge_over_an_edge_round() {
    let named = |e: &ogeom_core::OgeomError| e.to_string().contains("trim through is off");
    let mut m = Model::new();
    let (s, foot) = bossed_between_rounds(&mut m);
    let specs: Vec<(Shape, Chamfer)> = foot
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(3.0)))
        .collect();
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    for corner in [ChamferCorner::Point, ChamferCorner::Facet] {
        for options in [trim_off, Overflow::none()] {
            let err = chamfer_edges_overflow(&mut m, &s, &specs, corner, options, T)
                .map(|_| ())
                .unwrap_err();
            assert!(named(&err), "{corner:?}: {err}");
        }
    }

    let mut m = Model::new();
    let (s, foot, _) = bossed(&mut m, 12.0);
    assert!(foot.len() > 1, "the rim's arcs");
    for arc in &foot {
        let specs = vec![(arc.clone(), Chamfer::Symmetric(3.0))];
        let err = chamfer_edges_overflow(
            &mut m,
            &s,
            &specs,
            ChamferCorner::Facet,
            Overflow::none(),
            T,
        )
        .map(|_| ())
        .unwrap_err();
        assert!(
            named(&err) || err.to_string().contains("blend into round"),
            "{err}"
        );
    }
}
