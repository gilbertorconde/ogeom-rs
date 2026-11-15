//! Swept-section volumes for the rolling rounds' exact oracles: a ball's
//! section clipped to half-spaces, and the adaptive rule that integrates
//! the volume it sweeps.

use ogeom_math::Point;

/// The area of a ball's section kept within the end faces' planes `caps`
/// (each `n . p <= d`), and the volume rate it sweeps: the section, in the
/// plane square to the centre's path (`c`, `c'`, `c''`), is the triangle
/// from the seam's crossing `x` to the contacts `touch` and `foot` less the
/// circular segment the ball's arc cuts off it, and a plane turning with
/// the path sweeps `(|c'| - n'.(p - c)) dA` at each point `p` of it. The
/// arc is read as chords, each chord's sliver of the ball taken off where
/// both of its ends are kept (a sliver the clip cuts is a part in ten
/// million of a section).
pub fn clipped_section(
    caps: &[(ogeom_math::Vector, f64)],
    (c, d, dd): (Point, ogeom_math::Vector, ogeom_math::Vector),
    [x, touch, foot]: [Point; 3],
    radius: f64,
) -> (f64, f64) {
    use ogeom_math::Vector;
    const CHORDS: usize = 1024;
    let kept = |p: Point| caps.iter().all(|(n, d)| n.dot(p - Point::ORIGIN) <= *d);
    let (p1, q1) = ((touch - c) / radius, (foot - c) / radius);
    let axis = p1.cross(q1) / p1.cross(q1).magnitude();
    let sweep = p1.dot(q1).clamp(-1.0, 1.0).acos();
    #[allow(clippy::cast_precision_loss)]
    let angle = |i: usize| sweep * (i as f64) / (CHORDS as f64);
    let on_arc = |t: f64| c + (p1 * t.cos() + axis.cross(p1) * t.sin()) * radius;
    let mut ring: Vec<Point> = vec![x];
    ring.extend((0..=CHORDS).map(|i| on_arc(angle(i))));
    for (n, off) in caps {
        let side = |p: Point| n.dot(p - Point::ORIGIN) - off;
        let mut clipped = Vec::with_capacity(ring.len() + 2);
        for (i, &p) in ring.iter().enumerate() {
            let q = ring[(i + 1) % ring.len()];
            let (sp, sq) = (side(p), side(q));
            if sp <= 0.0 {
                clipped.push(p);
            }
            if (sp < 0.0 && sq > 0.0) || (sp > 0.0 && sq < 0.0) {
                clipped.push(p + (q - p) * (sp / (sp - sq)));
            }
        }
        ring = clipped;
    }
    let normal = d / d.magnitude();
    let (mut area, mut moment) = (0.0, Vector::ZERO);
    for i in 1..ring.len().saturating_sub(1) {
        let (o, p, q) = (ring[0], ring[i], ring[i + 1]);
        let piece = (p - o).cross(q - o).dot(normal) / 2.0;
        area += piece;
        moment += ((o - Point::ORIGIN) + (p - Point::ORIGIN) + (q - Point::ORIGIN)) * (piece / 3.0);
    }
    if area < 0.0 {
        area = -area;
        moment = -moment;
    }
    // The chords cut through the ball: each kept sliver comes off.
    #[allow(clippy::cast_precision_loss)]
    let step = sweep / (CHORDS as f64);
    let sliver = radius * radius * (step - step.sin()) / 2.0;
    let lever = 4.0 * radius * (step / 2.0).sin().powi(3) / (3.0 * (step - step.sin()));
    for i in 0..CHORDS {
        let (a0, a1) = (angle(i), angle(i + 1));
        if kept(on_arc(a0)) && kept(on_arc(a1)) {
            let mid = on_arc(0.5 * (a0 + a1));
            area -= sliver;
            moment -= ((c - Point::ORIGIN) + (mid - c) * (lever / radius)) * sliver;
        }
    }
    let speed = d.magnitude();
    let turning = (dd - normal * dd.dot(normal)) * (1.0 / speed);
    (
        area,
        speed * area - turning.dot(moment - (c - Point::ORIGIN) * area),
    )
}

/// The integral of `f` over `[lo, hi]` by adaptive Simpson on sixty-four
/// pieces, each halved until its two halves agree with it to 1e-13 (or
/// forty times): a clip puts kinks in a swept rate, which a fixed rule
/// would not see.
pub fn adaptive(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    #[allow(clippy::too_many_arguments, reason = "one step of adaptive Simpson")]
    fn simpson(
        f: &dyn Fn(f64) -> f64,
        a: f64,
        b: f64,
        fa: f64,
        fm: f64,
        fb: f64,
        whole: f64,
        depth: u32,
    ) -> f64 {
        let m = 0.5 * (a + b);
        let (l, r) = (0.5 * (a + m), 0.5 * (m + b));
        let (fl, fr) = (f(l), f(r));
        let left = (m - a) * (fa + 4.0 * fl + fm) / 6.0;
        let right = (b - m) * (fm + 4.0 * fr + fb) / 6.0;
        if depth == 0 || (left + right - whole).abs() <= 1e-13 {
            return left + right + (left + right - whole) / 15.0;
        }
        simpson(f, a, m, fa, fl, fm, left, depth - 1)
            + simpson(f, m, b, fm, fr, fb, right, depth - 1)
    }
    // Sixty-four pieces first: a clip narrow in the parameter (a post's
    // width) would read nought at three samples and stop there.
    (0..64)
        .map(|k| {
            let (a, b) = (
                lo + (hi - lo) * f64::from(k) / 64.0,
                lo + (hi - lo) * f64::from(k + 1) / 64.0,
            );
            let (fa, fm, fb) = (f(a), f(0.5 * (a + b)), f(b));
            let whole = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
            simpson(f, a, b, fa, fm, fb, whole, 40)
        })
        .sum()
}
