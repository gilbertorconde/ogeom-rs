//! Overflow handling (`Overflow`): a round reaching geometry other than its
//! own two faces. Blend into round bends a boss's foot round into the
//! straight round it reaches; each behaviour switched off passes the
//! request to the next, or it is refused by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_core::Tolerances;
use ogeom_fillet::{
    Blend, ChamferCorner, FilletCorner, Overflow, RadiusLaw, RoundProfile, blend_edges,
    blend_edges_overflow, fillet_edge_distances, fillet_edge_distances_overflow, fillet_edges,
    fillet_edges_corner, fillet_edges_corner_overflow, fillet_edges_overflow,
    fillet_edges_variable, fillet_edges_variable_overflow,
};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The block's size, the edge round's radius, the boss's radius and
/// height, and the foot round's radius.
const BLOCK: (f64, f64, f64) = (40.0, 30.0, 15.0);
const EDGE_ROUND: f64 = 5.0;
const BOSS: f64 = 6.0;
const BOSS_HIGH: f64 = 15.0;
const FOOT: f64 = 3.0;

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

/// The block 40 by 30 by 15, its top edge along `y = 0` rounded at 5, and
/// a boss of radius 6, 15 tall, standing on the top with its axis at
/// `(20, boss_y)`; with the boss's foot rim.
fn bossed(m: &mut Model, boss_y: f64) -> (Shape, Vec<Shape>) {
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
    (s, foot)
}

/// The independent reference: the fill a ball of the foot round's radius
/// leaves in the boss's foot, as the morphological closing defines it. A
/// point above the floor and outside the boss is in the fill where no ball
/// clear of the solid reaches it: where its distance from the set `Free`
/// of centres clear of the solid (at least the foot round's radius from the
/// boss's wall and from the floor and the edge round) is more than that
/// radius. `Free` is bounded near the foot by the cylinder about the boss's
/// axis and the floor's offset (the plane 3 above the top, and the
/// cylinder of radius 8 about the edge round's axis past its tangent line),
/// which meet along the spine. Its distance is the least of the
/// projections onto either boundary that land on the other's side, and the
/// distance from the spine, by sampling and golden-section search.
struct Closing {
    boss_y: f64,
}

impl Closing {
    /// The floor's height at `y`: the top, or the edge round past its
    /// tangent line `y = 5`.
    fn floor(y: f64) -> f64 {
        if y >= EDGE_ROUND {
            BLOCK.2
        } else {
            let d = EDGE_ROUND - y;
            (BLOCK.2 - EDGE_ROUND) + (EDGE_ROUND * EDGE_ROUND - d * d).max(0.0).sqrt()
        }
    }

    /// The signed clearance of a centre from the floor less the radius:
    /// non-negative where a ball there clears the floor and the edge round.
    fn above(c: Point) -> f64 {
        if c.y >= EDGE_ROUND {
            c.z - BLOCK.2 - FOOT
        } else {
            let (dy, dz) = (c.y - EDGE_ROUND, c.z - (BLOCK.2 - EDGE_ROUND));
            dy.hypot(dz) - EDGE_ROUND - FOOT
        }
    }

    fn rho(&self, p: Point) -> f64 {
        (p.x - 20.0).hypot(p.y - self.boss_y)
    }

    /// The spine at angle `t` about the boss's axis.
    fn spine(&self, t: f64) -> Point {
        let reach = BOSS + FOOT;
        let (x, y) = (20.0 + reach * t.cos(), self.boss_y + reach * t.sin());
        let z = if y >= EDGE_ROUND {
            BLOCK.2 + FOOT
        } else {
            let dy = EDGE_ROUND - y;
            let r = EDGE_ROUND + FOOT;
            (BLOCK.2 - EDGE_ROUND) + (r * r - dy * dy).sqrt()
        };
        Point::new(x, y, z)
    }

    fn distance_to_free(&self, p: Point) -> f64 {
        let reach = BOSS + FOOT;
        let free = |c: Point| self.rho(c) >= reach - 1e-12 && Self::above(c) >= -1e-12;
        if free(p) {
            return 0.0;
        }
        let mut best = f64::INFINITY;
        // Onto the cylinder about the boss's axis.
        let rho = self.rho(p);
        if rho > 0.0 {
            let (ux, uy) = ((p.x - 20.0) / rho, (p.y - self.boss_y) / rho);
            let q = Point::new(20.0 + reach * ux, self.boss_y + reach * uy, p.z);
            if Self::above(q) >= 0.0 {
                best = best.min(p.distance(q));
            }
        }
        // Onto the floor's offset: the plane over the top, the cylinder
        // past the tangent line.
        let flat = Point::new(p.x, p.y, BLOCK.2 + FOOT);
        if flat.y >= EDGE_ROUND && self.rho(flat) >= reach {
            best = best.min(p.distance(flat));
        }
        let (dy, dz) = (p.y - EDGE_ROUND, p.z - (BLOCK.2 - EDGE_ROUND));
        let m = dy.hypot(dz);
        if m > 0.0 {
            let k = (EDGE_ROUND + FOOT) / m;
            let q = Point::new(p.x, EDGE_ROUND + dy * k, (BLOCK.2 - EDGE_ROUND) + dz * k);
            if q.y <= EDGE_ROUND && self.rho(q) >= reach {
                best = best.min(p.distance(q));
            }
        }
        // The spine, near the point's own angle.
        let t0 = (p.y - self.boss_y).atan2(p.x - 20.0);
        let d = |t: f64| p.distance(self.spine(t));
        let (mut lo, mut hi) = (t0 - 1.2, t0 + 1.2);
        let mut at = t0;
        let mut least = f64::INFINITY;
        for k in 0..=96 {
            let t = lo + (hi - lo) * f64::from(k) / 96.0;
            if d(t) < least {
                least = d(t);
                at = t;
            }
        }
        lo = at - 2.4 / 96.0;
        hi = at + 2.4 / 96.0;
        let g = 0.5 * (5.0_f64.sqrt() - 1.0);
        for _ in 0..80 {
            let (a, b) = (hi - g * (hi - lo), lo + g * (hi - lo));
            if d(a) < d(b) {
                hi = b;
            } else {
                lo = a;
            }
        }
        best.min(d(0.5 * (lo + hi)))
    }

    /// The fill's height over the floor at `(x, y)`, by bisection on the
    /// height where the distance from `Free` comes down to the radius.
    fn column(&self, x: f64, y: f64) -> f64 {
        let z0 = Self::floor(y);
        if self.distance_to_free(Point::new(x, y, z0)) <= FOOT {
            return 0.0;
        }
        let (mut lo, mut hi) = (0.0, FOOT + 1.0);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.distance_to_free(Point::new(x, y, z0 + mid)) > FOOT {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The fill's volume: the columns integrated over the ring between the
    /// boss and the ball's reach, `rho = 9 - 3 cos s` (smooth where the
    /// round meets the wall), each angular range Gauss-Legendre of `n`
    /// points, split where the ball's seat changes.
    fn volume(&self, n: usize) -> f64 {
        let g = gauss(n);
        let reach = BOSS + FOOT;
        let gap = self.boss_y - EDGE_ROUND;
        let half = if gap < reach {
            (gap / reach).acos()
        } else {
            0.0
        };
        // Angles about the boss's axis from -y: the edge round's side.
        let ranges = [(-half, half), (half, TAU - half)];
        let mut sum = 0.0;
        for (a0, a1) in ranges {
            if a1 - a0 <= 0.0 {
                continue;
            }
            let ha = 0.5 * (a1 - a0);
            for (ga, wa) in &g {
                let a = a0 + ha * (ga + 1.0);
                let (ux, uy) = (a.sin(), -a.cos());
                for (gs, ws) in &g {
                    let s = FRAC_PI_2 * 0.25 * (gs + 1.0) * 2.0;
                    let rho = FOOT.mul_add(-s.cos(), reach);
                    let drho = FOOT * s.sin() * FRAC_PI_2 * 0.5;
                    let h = self.column(20.0 + rho * ux, self.boss_y + rho * uy);
                    sum += wa * ws * h * rho * drho * ha;
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

/// The largest angle between the normals of the two faces either side of
/// each edge bounding a face on a fitted surface, sampled along the edge's
/// own curve and projected onto each face's surface.
fn worst_crease(m: &Model, s: &Shape) -> f64 {
    let mut worst = 0.0_f64;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(m, a), surface_of(m, b));
        let fitted = |g: &SurfaceGeometry| matches!(g, SurfaceGeometry::BSpline(_));
        if !fitted(&sa) && !fitted(&sb) {
            continue;
        }
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 1..20 {
            let p = c.point_at(lo + (hi - lo) * f64::from(k) / 20.0, T).unwrap();
            let normal = |g: &SurfaceGeometry| -> Vector {
                let foot = ogeom_geom::project::project_on_surface(g, p, 16, T).unwrap();
                let (u, v) = foot.parameters;
                g.normal_at(u, v, T).unwrap().vector()
            };
            let (na, nb) = (normal(&sa), normal(&sb));
            let angle = na.dot(nb).abs().clamp(0.0, 1.0).acos();
            worst = worst.max(angle);
        }
    }
    worst
}

/// The issue's case: the boss at `y = 12`, its foot round reaching 3 past
/// the edge round's tangent line `y = 5`. The round rolls on over the edge
/// round: the result is valid and closed, every edge of the fitted piece
/// joins its neighbours G1, and its volume is the block's with the edge
/// round and the boss, plus the fill the closing reference gives. That
/// reference is checked against Pappus first on a boss clear of the edge
/// round (the torus's fill, `2 pi (6 + 3 (10 - 3 pi) / (12 - 3 pi)) 9 (1 -
/// pi / 4)`), and its own error bound is the change from 48 to 96 points
/// per range.
#[test]
fn a_boss_s_foot_round_blends_into_the_edge_round_it_reaches() {
    let clear = Closing { boss_y: 20.0 };
    let pappus = TAU
        * (BOSS + FOOT * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI))
        * FOOT
        * FOOT
        * (1.0 - PI / 4.0);
    let self_check = clear.volume(32);
    assert!(
        (self_check - pappus).abs() < 1e-6,
        "the reference {self_check} against Pappus {pappus}"
    );

    let mut m = Model::new();
    let (s, foot) = bossed(&mut m, 12.0);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &foot, FOOT, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
    }
    for e in &foot {
        assert!(!built.history.generated(e).is_empty(), "the rim's round");
    }

    let reference = Closing { boss_y: 12.0 };
    let (coarse, fine) = (reference.volume(48), reference.volume(96));
    let bound = (fine - coarse).abs().max(1e-9);
    let got = volume(&m, &built.shape) - before;
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "the fill {got} against the reference {fine} (within {bound})"
    );
    // Not the torus's fill: the round bends over the edge round.
    assert!(
        (got - pappus).abs() > 1e-2,
        "{got} against the torus's {pappus}"
    );

    let crease = worst_crease(&m, &built.shape);
    assert!(crease < 1e-5, "the round's joins crease by {crease}");
}

/// Blend into round off: the request goes on to trim through, the engine's
/// own construction, which refuses this foot by name, as before; every
/// behaviour off refuses it naming blend into round.
#[test]
fn blend_into_round_off_falls_to_trim_through_or_is_refused() {
    let mut m = Model::new();
    let (s, foot) = bossed(&mut m, 12.0);
    let mut off = Overflow::default();
    off.blend_into_round = false;
    let err = fillet_edges_overflow(&mut m, &s, &foot, FOOT, off, T).unwrap_err();
    assert!(
        err.to_string().contains("does not reach so far"),
        "trim through's refusal: {err}"
    );
    let err = fillet_edges_overflow(&mut m, &s, &foot, FOOT, Overflow::none(), T).unwrap_err();
    assert!(err.to_string().contains("blend into round"), "{err}");
    // Trim through off alone: blend into round still wins.
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let a = fillet_edges_overflow(&mut m, &s, &foot, FOOT, trim_off, T).unwrap();
    let b = fillet_edges(&mut m, &s, &foot, FOOT, T).unwrap();
    assert!((volume(&m, &a.shape) - volume(&m, &b.shape)).abs() < 1e-9);
}

/// A plate with a boss, its top rim rounded at 1 beside a pocket sunk into
/// the boss: the round reaches the pocket's walls and trims through them
/// (as `non_local_rims` measures). With trim through off it is refused by
/// name, hold edge not taking a rim; a local round with every behaviour off
/// builds as with all on.
#[test]
fn trim_through_off_refuses_a_round_reaching_another_face() {
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(&mut m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 4.0), 5.0, 5.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let pocket = ogeom_algo::make_box(&mut m, at(4.3, -0.25, 6.0), (0.5, 0.5, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &s, &pocket, T).unwrap().shape;
    let rim = edges_where(&m, &s, |p| {
        (p.z - 9.0).abs() < 1e-6 && (p.x.hypot(p.y) - 5.0).abs() < 1e-6
    });
    assert!(
        fillet_edges(&mut m, &s, &rim, 1.0, T).is_ok(),
        "trimmed through"
    );
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let err = fillet_edges_overflow(&mut m, &s, &rim, 1.0, trim_off, T).unwrap_err();
    assert!(
        err.to_string().contains("trim through is off"),
        "refused by name: {err}"
    );

    let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (10.0, 8.0, 6.0), T)
        .unwrap()
        .shape;
    let top = edges_where(&m, &block, |p| p.y.abs() < 1e-9 && (p.z - 6.0).abs() < 1e-9);
    let local = fillet_edges_overflow(&mut m, &block, &top, 2.0, Overflow::none(), T).unwrap();
    let want = 480.0 - 10.0 * 4.0 * (1.0 - PI / 4.0);
    assert!((volume(&m, &local.shape) - want).abs() < 1e-9);
}

/// One round call under test, taking the option set.
type RoundCall<'a> =
    dyn FnMut(&mut Model, Overflow) -> ogeom_core::OgeomResult<ogeom_algo::Built> + 'a;

/// The volume of what `built` builds on `m`.
fn measured(
    m: &mut Model,
    built: impl FnOnce(&mut Model) -> ogeom_core::OgeomResult<ogeom_algo::Built>,
) -> f64 {
    let b = built(m).unwrap();
    volume(m, &b.shape)
}

/// Blend into round off, every behaviour off, and trim through off alone,
/// for one round request: each option off gives the next option's result
/// or a refusal naming it, as [`fillet_edges_overflow`] does on the foot.
fn options_off_as_setback(m: &mut Model, round: &mut RoundCall<'_>, want: f64) {
    let mut off = Overflow::default();
    off.blend_into_round = false;
    let err = round(m, off).unwrap_err();
    assert!(
        err.to_string().contains("does not reach so far"),
        "trim through's refusal: {err}"
    );
    let err = round(m, Overflow::none()).unwrap_err();
    assert!(err.to_string().contains("blend into round"), "{err}");
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let built = round(m, trim_off).unwrap();
    assert!((volume(m, &built.shape) - want).abs() < 1e-9);
}

/// The foot with a rolling corner: no vertex where rounds meet, so the
/// round is the setback call's, with its options. The default is
/// `fillet_edges_corner`'s solid; blend into round off is refused as the
/// setback call refuses it.
#[test]
fn a_rolling_corner_takes_the_overflow_options() {
    let mut m = Model::new();
    let (s, foot) = bossed(&mut m, 12.0);
    let setback = fillet_edges(&mut m, &s, &foot, FOOT, T).unwrap();
    let want = volume(&m, &setback.shape);
    let rolled = fillet_edges_corner(&mut m, &s, &foot, FOOT, FilletCorner::Rolling, T).unwrap();
    let by_default = fillet_edges_corner_overflow(
        &mut m,
        &s,
        &foot,
        FOOT,
        FilletCorner::Rolling,
        Overflow::default(),
        T,
    )
    .unwrap();
    assert!((volume(&m, &rolled.shape) - want).abs() < 1e-9);
    assert!((volume(&m, &by_default.shape) - want).abs() < 1e-9);
    options_off_as_setback(
        &mut m,
        &mut |m, o| fillet_edges_corner_overflow(m, &s, &foot, FOOT, FilletCorner::Rolling, o, T),
        want,
    );
}

/// The foot rounded by a law holding 3 all along (linear and by points),
/// by two equal distances of the elliptical profile and by a plain round
/// in `blend_edges`: before the option set each was refused, the variable
/// and the blend calls rounding the foot without blend into round. With
/// the volume of the foot's solid before, and of the constant round's.
fn level_feet() -> (Model, Vec<ogeom_algo::Built>, f64, f64) {
    let mut m = Model::new();
    let (s, foot) = bossed(&mut m, 12.0);
    let before = volume(&m, &s);
    let constant = fillet_edges(&mut m, &s, &foot, FOOT, T).unwrap();
    let want = volume(&m, &constant.shape);
    let floor = ancestors_of(&m, &s, &foot[0], ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface_of(&m, f), SurfaceGeometry::Plane(_)))
        .unwrap();
    let distances: Vec<(Shape, Blend)> = foot
        .iter()
        .map(|e| {
            (
                e.clone(),
                Blend::RoundDistances {
                    face: floor.clone(),
                    on_face: FOOT,
                    on_other: FOOT,
                    profile: RoundProfile::Elliptical,
                },
            )
        })
        .collect();
    let rounds: Vec<(Shape, Blend)> = foot
        .iter()
        .map(|e| (e.clone(), Blend::Round(FOOT)))
        .collect();
    let level = RadiusLaw::Points(vec![(0.0, FOOT), (0.4, FOOT), (1.0, FOOT)]);
    let built = vec![
        fillet_edges_variable(
            &mut m,
            &s,
            &foot,
            &RadiusLaw::Linear {
                start: FOOT,
                end: FOOT,
            },
            T,
        )
        .unwrap(),
        fillet_edges_variable(&mut m, &s, &foot, &level, T).unwrap(),
        blend_edges(&mut m, &s, &distances, ChamferCorner::Facet, T).unwrap(),
        blend_edges(&mut m, &s, &rounds, ChamferCorner::Facet, T).unwrap(),
    ];
    (m, built, before, want)
}

/// Each of [`level_feet`] is now the constant round blended into the edge
/// round: the same volume within 1e-9, valid and closed, and G1 at every
/// join of its fitted piece.
#[test]
fn a_level_law_and_equal_distances_blend_into_the_edge_round() {
    let (m, built, _, want) = level_feet();
    for b in &built {
        let found = ogeom_algo::check(&m, &b.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for shell in explore_unique(&m, &b.shape, ShapeType::Shell).unwrap() {
            assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
        }
        let v = volume(&m, &b.shape);
        assert!((v - want).abs() < 1e-9, "{v} against the round's {want}");
        let crease = worst_crease(&m, &b.shape);
        assert!(crease < 1e-5, "the round's joins crease by {crease}");
    }
}

/// Each of [`level_feet`]'s fills against the closing reference, within
/// its bound (the change from 48 to 96 points per range) plus 1e-6.
#[test]
#[ignore = "heavy: the sampled reference takes most of a minute"]
fn a_level_law_and_equal_distances_fill_as_the_closing() {
    let (m, built, before, _) = level_feet();
    let reference = Closing { boss_y: 12.0 };
    let (coarse, fine) = (reference.volume(48), reference.volume(96));
    let bound = (fine - coarse).abs().max(1e-9);
    for b in &built {
        let got = volume(&m, &b.shape) - before;
        assert!(
            (got - fine).abs() < 1e-6 + bound,
            "the fill {got} against the reference {fine} (within {bound})"
        );
    }
}

/// The level law's and the equal distances' options off, as the setback
/// call's; a law that varies along the foot is not blended (blend into
/// round is built for one radius), so it keeps trim through's refusal
/// with every option on, and with every option off; two different
/// distances likewise keep theirs.
#[test]
fn variable_and_two_distance_feet_take_the_options_off_as_the_round() {
    let mut m = Model::new();
    let (s, foot) = bossed(&mut m, 12.0);
    let want = measured(&mut m, |m| fillet_edges(m, &s, &foot, FOOT, T));
    let level = RadiusLaw::Linear {
        start: FOOT,
        end: FOOT,
    };
    options_off_as_setback(
        &mut m,
        &mut |m, o| fillet_edges_variable_overflow(m, &s, &foot, &level, o, T),
        want,
    );
    let floor = ancestors_of(&m, &s, &foot[0], ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface_of(&m, f), SurfaceGeometry::Plane(_)))
        .unwrap();
    let asked = |on_other: f64| -> Vec<(Shape, Blend)> {
        foot.iter()
            .map(|e| {
                (
                    e.clone(),
                    Blend::RoundDistances {
                        face: floor.clone(),
                        on_face: FOOT,
                        on_other,
                        profile: RoundProfile::Elliptical,
                    },
                )
            })
            .collect()
    };
    let equal = asked(FOOT);
    options_off_as_setback(
        &mut m,
        &mut |m, o| blend_edges_overflow(m, &s, &equal, ChamferCorner::Facet, o, T),
        want,
    );

    let varying = RadiusLaw::Points(vec![(0.0, FOOT), (0.5, 2.5), (1.0, FOOT)]);
    let unequal = asked(2.0);
    for o in [Overflow::default(), Overflow::none()] {
        let err = fillet_edges_variable_overflow(&mut m, &s, &foot, &varying, o, T).unwrap_err();
        assert!(err.to_string().contains("does not reach so far"), "{err}");
        let err =
            blend_edges_overflow(&mut m, &s, &unequal, ChamferCorner::Facet, o, T).unwrap_err();
        assert!(err.to_string().contains("sets back past the side"), "{err}");
    }
}

/// The plate's rim beside a pocket: a varying law, a level one, a round
/// of two distances and the same in `blend_edges` all trim through the
/// pocket's walls by default, each as its call built it before the option
/// set, and are refused by name with trim through off.
#[test]
fn trim_through_off_refuses_variable_and_two_distance_rounds_reaching_a_face() {
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(&mut m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 4.0), 5.0, 5.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let pocket = ogeom_algo::make_box(&mut m, at(4.3, -0.25, 6.0), (0.5, 0.5, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &s, &pocket, T).unwrap().shape;
    let rim = edges_where(&m, &s, |p| {
        (p.z - 9.0).abs() < 1e-6 && (p.x.hypot(p.y) - 5.0).abs() < 1e-6
    });
    let top = ancestors_of(&m, &s, &rim[0], ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface_of(&m, f), SurfaceGeometry::Plane(_)))
        .unwrap();
    let varying = RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 0.8), (1.0, 1.0)]);
    let level = RadiusLaw::Linear {
        start: 1.0,
        end: 1.0,
    };
    let blends: Vec<(Shape, Blend)> = rim
        .iter()
        .map(|e| {
            (
                e.clone(),
                Blend::RoundDistances {
                    face: top.clone(),
                    on_face: 1.0,
                    on_other: 0.7,
                    profile: RoundProfile::Elliptical,
                },
            )
        })
        .collect();
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let calls: [&mut RoundCall<'_>; 4] = [
        &mut |m, o| fillet_edges_variable_overflow(m, &s, &rim, &varying, o, T),
        &mut |m, o| fillet_edges_variable_overflow(m, &s, &rim, &level, o, T),
        &mut |m, o| {
            fillet_edge_distances_overflow(
                m,
                &s,
                &rim[0],
                &top,
                1.0,
                0.7,
                RoundProfile::Elliptical,
                o,
                T,
            )
        },
        &mut |m, o| blend_edges_overflow(m, &s, &blends, ChamferCorner::Facet, o, T),
    ];
    let before = [
        measured(&mut m, |m| fillet_edges_variable(m, &s, &rim, &varying, T)),
        measured(&mut m, |m| fillet_edges_variable(m, &s, &rim, &level, T)),
        measured(&mut m, |m| {
            fillet_edge_distances(m, &s, &rim[0], &top, 1.0, 0.7, RoundProfile::Elliptical, T)
        }),
        measured(&mut m, |m| {
            blend_edges(m, &s, &blends, ChamferCorner::Facet, T)
        }),
    ];
    for (call, want) in calls.into_iter().zip(before) {
        let got = call(&mut m, Overflow::default()).unwrap();
        assert!((volume(&m, &got.shape) - want).abs() < 1e-9);
        let err = call(&mut m, trim_off).unwrap_err();
        assert!(err.to_string().contains("trim through is off"), "{err}");
    }
}

/// A local round of each new call builds the same solid with every
/// behaviour off as with every one on: a block's top edge, rounded by a
/// level and a varying law, by two distances, and in `blend_edges`.
#[test]
fn local_rounds_build_alike_with_every_option_off() {
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (10.0, 8.0, 6.0), T)
        .unwrap()
        .shape;
    let top = edges_where(&m, &block, |p| p.y.abs() < 1e-9 && (p.z - 6.0).abs() < 1e-9);
    let face = ancestors_of(&m, &block, &top[0], ShapeType::Face).unwrap()[0].clone();
    let law = RadiusLaw::Linear {
        start: 1.0,
        end: 2.0,
    };
    let blends = vec![(top[0].clone(), Blend::Round(2.0))];
    let calls: [&mut RoundCall<'_>; 5] = [
        &mut |m, o| fillet_edges_variable_overflow(m, &block, &top, &law, o, T),
        &mut |m, o| {
            fillet_edges_variable_overflow(
                m,
                &block,
                &top,
                &RadiusLaw::Linear {
                    start: 2.0,
                    end: 2.0,
                },
                o,
                T,
            )
        },
        &mut |m, o| {
            fillet_edge_distances_overflow(
                m,
                &block,
                &top[0],
                &face,
                2.0,
                1.0,
                RoundProfile::Elliptical,
                o,
                T,
            )
        },
        &mut |m, o| blend_edges_overflow(m, &block, &blends, ChamferCorner::Facet, o, T),
        &mut |m, o| fillet_edges_corner_overflow(m, &block, &top, 2.0, FilletCorner::Rolling, o, T),
    ];
    for call in calls {
        let on = call(&mut m, Overflow::default()).unwrap();
        let off = call(&mut m, Overflow::none()).unwrap();
        assert!((volume(&m, &on.shape) - volume(&m, &off.shape)).abs() < 1e-12);
    }
}
