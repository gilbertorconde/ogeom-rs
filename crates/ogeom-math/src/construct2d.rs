//! The classical 2D constructions: circles tangent to three entities,
//! tangent lines, and bisector curves — the straightedge-and-compass
//! repertoire, solved algebraically.
//!
//! One linearization carries the whole tangency family. A circle with
//! centre `c` and radius `r` is tangent to a target once a *side* is
//! chosen, and with that side fixed every constraint is linear in
//! `(c, r, Q)` where `Q = |c|² − r²`:
//!
//! - a target circle `(cᵢ, rᵢ)` on side `sᵢ`: `Q − 2cᵢ·c − 2sᵢrᵢr = rᵢ² − |cᵢ|²`;
//! - a point is a zero-radius circle;
//! - a line with unit normal `n` and offset `d` on side `t`: `n·c − t·r = d`,
//!   with no `Q` at all.
//!
//! Three constraints give three linear equations in at most four unknowns;
//! the solution family is a line, and re-imposing `Q = |c|² − r²` is a
//! quadratic along it. Enumerating the sides, solving, and *verifying every
//! candidate against the literal tangency distances* — the linearization can
//! manufacture roots the geometry rejects — yields exactly the classical
//! solution sets, Apollonius's eight included.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::conic::{Circle2, Ellipse2, Hyperbola2, Parabola2};
use crate::direction::Direction2;
use crate::frame::{Axis2, Frame2};
use crate::point::Point2;
use crate::vector::Vector2;

/// An entity a construction can be tangent to, or equidistant from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Target2 {
    /// A point — a zero-radius circle for tangency, itself for distance.
    Point(Point2),
    /// An unbounded line.
    Line(Axis2),
    /// A circle.
    Circle(Circle2),
}

impl Target2 {
    /// The distance from `p` to this target's own locus — for a circle, the
    /// distance to its *boundary*.
    #[must_use]
    pub fn distance_to(&self, p: Point2) -> f64 {
        match self {
            Self::Point(q) => p.distance(*q),
            Self::Line(axis) => axis.distance_to(p),
            Self::Circle(c) => (p.distance(c.centre()) - c.radius()).abs(),
        }
    }
}

/// How a tangent circle stands to one of its targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placement {
    /// Passes through a point target.
    Through,
    /// Touches a line target.
    Tangent,
    /// Touches a circle target from outside — the circles exclude each
    /// other.
    Outside,
    /// The solution contains the target circle.
    Enclosing,
    /// The target circle contains the solution.
    Enclosed,
}

/// One tangent circle, with its standing toward each target in order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentCircle {
    /// The solution.
    pub circle: Circle2,
    /// How it stands to each target, in the order they were given.
    pub placements: [Placement; 3],
}

/// A row of the linearized system: coefficients of `(cx, cy, r, Q)` and the
/// right-hand side.
type Row = ([f64; 4], f64);

fn rows_for(target: &Target2, side: f64) -> Row {
    match target {
        Target2::Point(p) => {
            // Q − 2p·c = −|p|²
            ([-2.0 * p.x, -2.0 * p.y, 0.0, 1.0], -(p.x * p.x + p.y * p.y))
        }
        Target2::Circle(c) => {
            let centre = c.centre();
            let r = c.radius();
            (
                [-2.0 * centre.x, -2.0 * centre.y, -2.0 * side * r, 1.0],
                r * r - (centre.x * centre.x + centre.y * centre.y),
            )
        }
        Target2::Line(axis) => {
            let n = normal_of(axis);
            let d = n.dot(axis.location.to_vector());
            ([n.x, n.y, -side, 0.0], d)
        }
    }
}

/// The unit left normal of a line.
fn normal_of(axis: &Axis2) -> Vector2 {
    let d = axis.direction.vector();
    Vector2::new(-d.y, d.x)
}

/// The sides to enumerate for one target: points have no side.
fn sides_of(target: &Target2) -> &'static [f64] {
    match target {
        Target2::Point(_) => &[1.0],
        _ => &[1.0, -1.0],
    }
}

/// Circles tangent to all three targets — the Apollonius family and its
/// degenerate relatives, every candidate verified against the literal
/// tangency distances before it is returned.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if two
/// targets coincide, which asks a different, underdetermined question.
pub fn circles_tangent_to_three(
    targets: &[Target2; 3],
    tol: Tolerances,
) -> OgeomResult<Vec<TangentCircle>> {
    for i in 0..3 {
        for j in i + 1..3 {
            if targets_coincide(&targets[i], &targets[j], tol) {
                ogeom_bail!(
                    Construction,
                    "targets {i} and {j} coincide; the tangency family is underdetermined"
                );
            }
        }
    }
    let mut out: Vec<TangentCircle> = Vec::new();
    for &s0 in sides_of(&targets[0]) {
        for &s1 in sides_of(&targets[1]) {
            for &s2 in sides_of(&targets[2]) {
                let rows = [
                    rows_for(&targets[0], s0),
                    rows_for(&targets[1], s1),
                    rows_for(&targets[2], s2),
                ];
                for candidate in solve_rows(&rows, tol) {
                    admit(&mut out, candidate, targets, tol);
                }
            }
        }
    }
    Ok(out)
}

/// Circles of a fixed radius tangent to two targets — the same machinery
/// with the radius row supplied.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// radius is not finite and positive, or the targets coincide.
pub fn circles_of_radius_tangent_to_two(
    radius: f64,
    targets: &[Target2; 2],
    tol: Tolerances,
) -> OgeomResult<Vec<TangentCircle>> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a tangent circle of radius {radius} is not a circle"
        );
    }
    if targets_coincide(&targets[0], &targets[1], tol) {
        ogeom_bail!(Construction, "the two targets coincide");
    }
    let mut out: Vec<TangentCircle> = Vec::new();
    let radius_row: Row = ([0.0, 0.0, 1.0, 0.0], radius);
    for &s0 in sides_of(&targets[0]) {
        for &s1 in sides_of(&targets[1]) {
            let rows = [
                rows_for(&targets[0], s0),
                rows_for(&targets[1], s1),
                radius_row,
            ];
            for candidate in solve_rows(&rows, tol) {
                let three = [targets[0], targets[1], targets[1]];
                let mut kept = out.clone();
                admit(&mut kept, candidate, &three, tol);
                // Re-derive placements for the two real targets only.
                if kept.len() > out.len() {
                    let solution = kept[kept.len() - 1].circle;
                    let placements = [
                        placement_of(&solution, &targets[0], tol),
                        placement_of(&solution, &targets[1], tol),
                        placement_of(&solution, &targets[1], tol),
                    ];
                    out.push(TangentCircle {
                        circle: solution,
                        placements,
                    });
                }
            }
        }
    }
    Ok(out)
}

/// The up-to-four lines tangent to two circles: the external pair where the
/// circles lie on the same side, the internal pair where they straddle.
#[must_use]
pub fn lines_tangent_to_two_circles(a: &Circle2, b: &Circle2, tol: Tolerances) -> Vec<Axis2> {
    let e = b.centre() - a.centre();
    let distance = e.magnitude();
    if distance <= tol.confusion() {
        return Vec::new();
    }
    let along = e / distance;
    let across = Vector2::new(-along.y, along.x);
    let mut out = Vec::new();
    // Unit normal n with n·(ca − cb) = s_a·ra − s_b·rb, d = n·ca − s_a·ra.
    for (sa, sb) in [(1.0, 1.0), (1.0, -1.0)] {
        let k = (sa * a.radius() - sb * b.radius()) / distance;
        if k.abs() > 1.0 - tol.angular() {
            continue;
        }
        let across_part = (1.0 - k * k).sqrt();
        for flip in [1.0, -1.0] {
            let n = along * -k + across * (across_part * flip);
            let d = n.dot(a.centre().to_vector()) - sa * a.radius();
            // The line's own frame: direction perpendicular to n, located at
            // the foot nearest the midpoint of the centres.
            let mid = a.centre() + e * 0.5;
            let foot = mid - n * (n.dot(mid.to_vector()) - d);
            if let Ok(direction) = Direction2::new(Vector2::new(n.y, -n.x), tol) {
                out.push(Axis2::new(foot, direction));
            }
        }
    }
    out
}

/// Solve three rows for `(c, r)` candidates: direct where `Q` is absent,
/// the line-family-plus-quadratic where it is present.
fn solve_rows(rows: &[Row; 3], tol: Tolerances) -> Vec<(Point2, f64)> {
    let uses_q = rows.iter().any(|(coeffs, _)| coeffs[3] != 0.0);
    if uses_q {
        solve_with_q(rows, tol)
    } else {
        solve_linear(rows, tol)
    }
}

/// All-lines: three equations in `(cx, cy, r)`.
fn solve_linear(rows: &[Row; 3], _tol: Tolerances) -> Vec<(Point2, f64)> {
    let m = nalgebra::Matrix3::new(
        rows[0].0[0],
        rows[0].0[1],
        rows[0].0[2],
        rows[1].0[0],
        rows[1].0[1],
        rows[1].0[2],
        rows[2].0[0],
        rows[2].0[1],
        rows[2].0[2],
    );
    let b = nalgebra::Vector3::new(rows[0].1, rows[1].1, rows[2].1);
    let Some(solution) = m.lu().solve(&b) else {
        return Vec::new();
    };
    vec![(Point2::new(solution[0], solution[1]), solution[2])]
}

/// With `Q` present: the 3×4 system's solution line, cut by the quadratic
/// `|c|² − r² − Q = 0`. The null vector comes from the four 3×3 minors —
/// the generalized cross product — and the particular solution from the
/// best-conditioned 3×3 subsystem with the remaining unknown pinned to
/// zero.
fn solve_with_q(rows: &[Row; 3], tol: Tolerances) -> Vec<(Point2, f64)> {
    let m = [rows[0].0, rows[1].0, rows[2].0];
    let b = [rows[0].1, rows[1].1, rows[2].1];

    // Minor j: the determinant with column j removed, alternating sign.
    let minor = |skip: usize| -> f64 {
        let cols: Vec<usize> = (0..4).filter(|c| *c != skip).collect();

        nalgebra::Matrix3::new(
            m[0][cols[0]],
            m[0][cols[1]],
            m[0][cols[2]],
            m[1][cols[0]],
            m[1][cols[1]],
            m[1][cols[2]],
            m[2][cols[0]],
            m[2][cols[1]],
            m[2][cols[2]],
        )
        .determinant()
    };
    let null: [f64; 4] = [minor(0), -minor(1), minor(2), -minor(3)];
    let biggest = null.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    if biggest <= 1e-12 {
        // Rank below three: a degenerate side pattern.
        return Vec::new();
    }

    // Particular solution: pin the unknown whose removal leaves the
    // best-conditioned square system.
    let pin = (0..4)
        .max_by(|a, b| {
            minor(*a)
                .abs()
                .partial_cmp(&minor(*b).abs())
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .unwrap_or(3);
    let cols: Vec<usize> = (0..4).filter(|c| *c != pin).collect();
    let square = nalgebra::Matrix3::new(
        m[0][cols[0]],
        m[0][cols[1]],
        m[0][cols[2]],
        m[1][cols[0]],
        m[1][cols[1]],
        m[1][cols[2]],
        m[2][cols[0]],
        m[2][cols[1]],
        m[2][cols[2]],
    );
    let rhs = nalgebra::Vector3::new(b[0], b[1], b[2]);
    let Some(solved) = square.lu().solve(&rhs) else {
        return Vec::new();
    };
    let mut particular = [0.0f64; 4];
    for (slot, col) in cols.iter().enumerate() {
        particular[*col] = solved[slot];
    }

    // g(λ) = |c|² − r² − Q along x = particular + λ·null.
    let (px, py, pr, pq) = (particular[0], particular[1], particular[2], particular[3]);
    let (nx, ny, nr, nq) = (null[0], null[1], null[2], null[3]);
    let a2 = nx * nx + ny * ny - nr * nr;
    let a1 = 2.0 * (px * nx + py * ny - pr * nr) - nq;
    let a0 = px * px + py * py - pr * pr - pq;

    let mut lambdas = Vec::new();
    if a2.abs() <= 1e-14 * (a1.abs().max(a0.abs()).max(1.0)) {
        if a1.abs() > 1e-14 {
            lambdas.push(-a0 / a1);
        }
    } else {
        // The quadratic's vertex is always a candidate: a tangent
        // configuration's double root sits exactly there, and rounding
        // renders its discriminant a hair negative. The literal tangency
        // verification downstream rejects the vertex whenever it is not a
        // real solution, so offering it costs nothing and loses nothing.
        lambdas.push(-a1 / (2.0 * a2));
        let disc = a1.mul_add(a1, -4.0 * a2 * a0);
        if disc > 0.0 {
            let root = disc.sqrt();
            lambdas.push((-a1 + root) / (2.0 * a2));
            lambdas.push((-a1 - root) / (2.0 * a2));
        }
    }
    lambdas
        .into_iter()
        .map(|l| (Point2::new(px + l * nx, py + l * ny), pr + l * nr))
        .filter(|(_, r)| r.is_finite() && *r > tol.confusion())
        .collect()
}

/// Verify a candidate against the literal tangency distances and admit it
/// once.
fn admit(
    out: &mut Vec<TangentCircle>,
    (centre, radius): (Point2, f64),
    targets: &[Target2; 3],
    tol: Tolerances,
) {
    let slack = tol.confusion() * 1e3 * radius.max(1.0);
    for target in targets {
        let touch = match target {
            Target2::Point(p) => (centre.distance(*p) - radius).abs(),
            Target2::Line(axis) => (axis.distance_to(centre) - radius).abs(),
            Target2::Circle(c) => {
                let d = centre.distance(c.centre());
                (d - (radius + c.radius()))
                    .abs()
                    .min((d - (radius - c.radius()).abs()).abs())
            }
        };
        if touch > slack {
            return;
        }
    }
    if out.iter().any(|held| {
        held.circle.centre().distance(centre) <= slack
            && (held.circle.radius() - radius).abs() <= slack
    }) {
        return;
    }
    let Ok(circle) = Circle2::new(Frame2::new(centre, Direction2::X), radius, tol) else {
        return;
    };
    let placements = [
        placement_of(&circle, &targets[0], tol),
        placement_of(&circle, &targets[1], tol),
        placement_of(&circle, &targets[2], tol),
    ];
    out.push(TangentCircle { circle, placements });
}

fn placement_of(circle: &Circle2, target: &Target2, tol: Tolerances) -> Placement {
    match target {
        Target2::Point(_) => Placement::Through,
        Target2::Line(_) => Placement::Tangent,
        Target2::Circle(c) => {
            let d = circle.centre().distance(c.centre());
            let slack = tol.confusion() * 1e3 * circle.radius().max(1.0);
            if (d - (circle.radius() + c.radius())).abs() <= slack {
                Placement::Outside
            } else if circle.radius() >= c.radius()
                && (d - (circle.radius() - c.radius())).abs() <= slack
            {
                Placement::Enclosing
            } else {
                Placement::Enclosed
            }
        }
    }
}

fn targets_coincide(a: &Target2, b: &Target2, tol: Tolerances) -> bool {
    match (a, b) {
        (Target2::Point(p), Target2::Point(q)) => p.is_equal(*q, tol),
        (Target2::Circle(c), Target2::Circle(d)) => {
            c.centre().is_equal(d.centre(), tol)
                && (c.radius() - d.radius()).abs() <= tol.confusion()
        }
        (Target2::Line(a), Target2::Line(b)) => {
            let na = normal_of(a);
            let nb = normal_of(b);
            na.cross(nb).abs() <= tol.angular() && a.distance_to(b.location) <= tol.confusion()
        }
        _ => false,
    }
}

// --- bisectors ---------------------------------------------------------------

/// The locus of points equidistant from two targets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Bisector2 {
    /// A single line — two points, or two parallel lines.
    Line(Axis2),
    /// The two angle bisectors of intersecting lines.
    Pair([Axis2; 2]),
    /// Point against line, or line against circle: a parabola.
    Parabola(Parabola2),
    /// A point inside a circle, or nested circles: an ellipse.
    Ellipse(Ellipse2),
    /// A point outside a circle, or circles of unequal radius: a hyperbola.
    /// The equidistant locus is the branch on the frame's `+x` side — toward
    /// the point, or toward the smaller circle; the mirror branch comes with
    /// the conic but bisects nothing.
    Hyperbola(Hyperbola2),
}

/// The bisector of two targets: the equidistant locus, as the conic it is.
///
/// For circle targets the distance is to the *boundary*, which is what makes
/// the answer a conic with the centres as foci.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// targets coincide, or a point lies on a line or circle target — the locus
/// then degenerates to something that is not a curve.
pub fn bisector(a: &Target2, b: &Target2, tol: Tolerances) -> OgeomResult<Bisector2> {
    if targets_coincide(a, b, tol) {
        ogeom_bail!(Construction, "coincident targets bisect everywhere");
    }
    // Normalize the order so each pair is handled once.
    match (a, b) {
        (Target2::Point(p), Target2::Point(q)) => {
            let mid = *p + (*q - *p) * 0.5;
            let direction = Direction2::new(perp(*q - *p), tol)?;
            Ok(Bisector2::Line(Axis2::new(mid, direction)))
        }

        (Target2::Line(l), Target2::Line(m)) => {
            let nl = normal_of(l);
            let nm = normal_of(m);
            let dl = nl.dot(l.location.to_vector());
            let dm = nm.dot(m.location.to_vector());
            if nl.cross(nm).abs() <= tol.angular() {
                // Parallel: the midline. Align the normals first.
                let (nm, dm) = if nl.dot(nm) < 0.0 {
                    (-nm, -dm)
                } else {
                    (nm, dm)
                };
                let _ = nm;
                let offset = f64::midpoint(dl, dm);
                let foot = Point2::new(nl.x * offset, nl.y * offset);
                return Ok(Bisector2::Line(Axis2::new(foot, l.direction)));
            }
            // Intersecting: n_l·x − d_l = ±(n_m·x − d_m).
            let apex = intersect_lines(nl, dl, nm, dm)?;
            let d1 = Direction2::new(l.direction.vector() + m.direction.vector(), tol)
                .or_else(|_| Direction2::new(perp(l.direction.vector()), tol))?;
            let d2 = Direction2::new(perp(d1.vector()), tol)?;
            Ok(Bisector2::Pair([
                Axis2::new(apex, d1),
                Axis2::new(apex, d2),
            ]))
        }

        (Target2::Point(p), Target2::Line(l)) | (Target2::Line(l), Target2::Point(p)) => {
            let n = normal_of(l);
            let signed = n.dot(*p - l.location);
            if signed.abs() <= tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the point lies on the line; the locus degenerates"
                );
            }
            // Focus at the point, directrix the line: apex midway, opening
            // from the line toward the point.
            let foot = *p - n * signed;
            let apex = foot + (*p - foot) * 0.5;
            let x = Direction2::new(*p - foot, tol)?;
            let frame = Frame2::new(apex, x);
            Ok(Bisector2::Parabola(Parabola2::new(
                frame,
                signed.abs() / 2.0,
                tol,
            )?))
        }

        (Target2::Point(p), Target2::Circle(c)) | (Target2::Circle(c), Target2::Point(p)) => {
            let spread = p.distance(c.centre());
            let r = c.radius();
            if (spread - r).abs() <= tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the point lies on the circle; the locus degenerates"
                );
            }
            foci_conic(c.centre(), *p, r, spread, tol)
        }

        (Target2::Line(l), Target2::Circle(c)) | (Target2::Circle(c), Target2::Line(l)) => {
            let n = normal_of(l);
            let signed = n.dot(c.centre() - l.location);
            if signed.abs() <= c.radius() + tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the line meets the circle; the equidistant locus is not one conic"
                );
            }
            // |x − centre| − r = distance to line, on the circle's side:
            // a parabola with the centre as focus and the line shifted r
            // toward the circle... away from it, as the directrix.
            let toward = if signed > 0.0 { n } else { -n };
            let directrix_foot = l.location + perp_foot_shift(l, c.centre()) - toward * c.radius();
            let focus = c.centre();
            let foot_to_focus = focus - directrix_foot;
            let apex = directrix_foot + foot_to_focus * 0.5;
            let x = Direction2::new(foot_to_focus, tol)?;
            Ok(Bisector2::Parabola(Parabola2::new(
                Frame2::new(apex, x),
                foot_to_focus.magnitude() / 2.0,
                tol,
            )?))
        }

        (Target2::Circle(c1), Target2::Circle(c2)) => {
            let spread = c1.centre().distance(c2.centre());
            if spread <= tol.confusion() {
                // Concentric: the midway circle.
                let radius = f64::midpoint(c1.radius(), c2.radius());
                let circle = Circle2::new(Frame2::new(c1.centre(), Direction2::X), radius, tol)?;
                let _ = circle;
                ogeom_bail!(
                    Construction,
                    "concentric circles bisect on a circle; ask for it as one"
                );
            }
            if (c1.radius() - c2.radius()).abs() <= tol.confusion() {
                // Equal radii: the perpendicular bisector of the centres.
                let mid = c1.centre() + (c2.centre() - c1.centre()) * 0.5;
                let direction = Direction2::new(perp(c2.centre() - c1.centre()), tol)?;
                return Ok(Bisector2::Line(Axis2::new(mid, direction)));
            }
            // | |x−c1| − |x−c2| | = |r1 − r2|: a hyperbola with the centres
            // as foci.
            let difference = (c1.radius() - c2.radius()).abs();
            if difference >= spread - tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "one circle encloses the other too deeply; the locus degenerates"
                );
            }
            let centre = c1.centre() + (c2.centre() - c1.centre()) * 0.5;
            let x = Direction2::new(c2.centre() - c1.centre(), tol)?;
            let a_half = difference / 2.0;
            let c_half = spread / 2.0;
            let b_half = (c_half * c_half - a_half * a_half).sqrt();
            Ok(Bisector2::Hyperbola(Hyperbola2::new(
                Frame2::new(centre, x),
                a_half,
                b_half,
                tol,
            )?))
        }
    }
}

/// The conic with foci at `f1`, `f2` where the boundary-distance equality
/// gives `|x−f1| ± |x−f2| = r`: an ellipse when the point sits inside the
/// circle, a hyperbola outside.
fn foci_conic(
    circle_centre: Point2,
    point: Point2,
    r: f64,
    spread: f64,
    tol: Tolerances,
) -> OgeomResult<Bisector2> {
    let centre = circle_centre + (point - circle_centre) * 0.5;
    let x = Direction2::new(point - circle_centre, tol)?;
    let a_half = r / 2.0;
    let c_half = spread / 2.0;
    if spread < r {
        // Inside: |x−centre| + |x−p| = r, an ellipse.
        let b_half = (a_half * a_half - c_half * c_half).sqrt();
        Ok(Bisector2::Ellipse(Ellipse2::new(
            Frame2::new(centre, x),
            a_half,
            b_half,
            tol,
        )?))
    } else {
        // Outside: |x−centre| − |x−p| = ±r, a hyperbola.
        let b_half = (c_half * c_half - a_half * a_half).sqrt();
        Ok(Bisector2::Hyperbola(Hyperbola2::new(
            Frame2::new(centre, x),
            a_half,
            b_half,
            tol,
        )?))
    }
}

fn perp(v: Vector2) -> Vector2 {
    Vector2::new(-v.y, v.x)
}

/// The component of `to − axis.location` along the line — the foot offset.
fn perp_foot_shift(axis: &Axis2, to: Point2) -> Vector2 {
    let along = axis.direction.vector();
    along * along.dot(to - axis.location)
}

fn intersect_lines(n1: Vector2, d1: f64, n2: Vector2, d2: f64) -> OgeomResult<Point2> {
    let det = n1.x * n2.y - n1.y * n2.x;
    if det.abs() <= f64::MIN_POSITIVE {
        ogeom_bail!(Construction, "parallel lines do not meet");
    }
    Ok(Point2::new(
        (d1 * n2.y - d2 * n1.y) / det,
        (n1.x * d2 - n2.x * d1) / det,
    ))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    fn circle(x: f64, y: f64, r: f64) -> Circle2 {
        Circle2::new(Frame2::new(Point2::new(x, y), Direction2::X), r, T).unwrap()
    }

    /// Every returned circle touches every target, by measurement.
    fn assert_tangent(solutions: &[TangentCircle], targets: &[Target2; 3]) {
        assert!(!solutions.is_empty(), "the construction found nothing");
        for s in solutions {
            for target in targets {
                let gap = match target {
                    Target2::Point(p) => (s.circle.centre().distance(*p) - s.circle.radius()).abs(),
                    Target2::Line(l) => {
                        (l.distance_to(s.circle.centre()) - s.circle.radius()).abs()
                    }
                    Target2::Circle(c) => {
                        let d = s.circle.centre().distance(c.centre());
                        (d - (s.circle.radius() + c.radius()))
                            .abs()
                            .min((d - (s.circle.radius() - c.radius()).abs()).abs())
                    }
                };
                assert!(gap < 1e-9, "tangency gap {gap} on {target:?} for {s:?}");
            }
        }
    }

    #[test]
    fn three_points_give_the_circumcircle() {
        let targets = [
            Target2::Point(Point2::new(0.0, 0.0)),
            Target2::Point(Point2::new(4.0, 0.0)),
            Target2::Point(Point2::new(0.0, 3.0)),
        ];
        let found = circles_tangent_to_three(&targets, T).unwrap();
        assert_eq!(found.len(), 1);
        // The 3-4-5 right triangle's circumradius is the hypotenuse over two.
        assert!((found[0].circle.radius() - 2.5).abs() < 1e-9);
        assert_tangent(&found, &targets);
    }

    #[test]
    fn three_lines_give_the_incircle_and_excircles() {
        // The 3-4-5 right triangle: incircle radius 1, three excircles.
        let targets = [
            Target2::Line(Axis2::new(Point2::new(0.0, 0.0), Direction2::X)),
            Target2::Line(Axis2::new(Point2::new(0.0, 0.0), Direction2::Y)),
            Target2::Line(Axis2::new(
                Point2::new(4.0, 0.0),
                Direction2::new(Vector2::new(-4.0, 3.0), T).unwrap(),
            )),
        ];
        let found = circles_tangent_to_three(&targets, T).unwrap();
        assert_eq!(found.len(), 4, "incircle and three excircles: {found:?}");
        assert!(
            found.iter().any(|s| (s.circle.radius() - 1.0).abs() < 1e-9),
            "the incircle of 3-4-5 has radius 1"
        );
        assert_tangent(&found, &targets);
    }

    #[test]
    fn apollonius_three_circles_yields_eight() {
        // The classical configuration: three mutually external circles in
        // general position give all eight Apollonius circles.
        let targets = [
            Target2::Circle(circle(0.0, 0.0, 1.0)),
            Target2::Circle(circle(6.0, 0.0, 1.5)),
            Target2::Circle(circle(2.5, 5.0, 2.0)),
        ];
        let found = circles_tangent_to_three(&targets, T).unwrap();
        assert_eq!(found.len(), 8, "Apollonius promises eight: {}", found.len());
        assert_tangent(&found, &targets);
        // Among them, one touches all three from outside and one encloses
        // all three.
        assert!(
            found
                .iter()
                .any(|s| s.placements == [Placement::Outside; 3])
        );
        assert!(
            found
                .iter()
                .any(|s| s.placements == [Placement::Enclosing; 3])
        );
    }

    #[test]
    fn mixed_targets_and_fixed_radius_answer() {
        let targets = [
            Target2::Point(Point2::new(1.0, 2.0)),
            Target2::Line(Axis2::new(Point2::new(0.0, -1.0), Direction2::X)),
            Target2::Circle(circle(5.0, 3.0, 1.0)),
        ];
        let found = circles_tangent_to_three(&targets, T).unwrap();
        assert_tangent(&found, &targets);

        let two = [
            Target2::Line(Axis2::new(Point2::new(0.0, 0.0), Direction2::X)),
            Target2::Circle(circle(0.0, 5.0, 1.0)),
        ];
        let sized = circles_of_radius_tangent_to_two(2.0, &two, T).unwrap();
        assert!(!sized.is_empty());
        for s in &sized {
            assert!((s.circle.radius() - 2.0).abs() < 1e-9);
            let d0 = Target2::distance_to(&two[0], s.circle.centre());
            let d1 = Target2::distance_to(&two[1], s.circle.centre());
            assert!((d0 - 2.0).abs() < 1e-9 && (d1 - 2.0).abs() < 1e-9, "{s:?}");
        }
    }

    #[test]
    fn bitangent_lines_touch_both_circles() {
        let a = circle(0.0, 0.0, 2.0);
        let b = circle(8.0, 0.0, 1.0);
        let lines = lines_tangent_to_two_circles(&a, &b, T);
        assert_eq!(lines.len(), 4, "external pair and internal pair");
        for line in &lines {
            assert!((line.distance_to(a.centre()) - 2.0).abs() < 1e-9);
            assert!((line.distance_to(b.centre()) - 1.0).abs() < 1e-9);
        }
    }

    /// Sample a bisector and assert the defining property: equidistance
    /// from both targets, measured literally.
    fn assert_equidistant(bisector: &Bisector2, a: &Target2, b: &Target2) {
        let probes: Vec<Point2> = match bisector {
            Bisector2::Line(axis) => (-5..=5)
                .map(|i| axis.location + axis.direction.vector() * f64::from(i))
                .collect(),
            Bisector2::Pair(axes) => axes
                .iter()
                .flat_map(|axis| {
                    (-3..=3).map(move |i| axis.location + axis.direction.vector() * f64::from(i))
                })
                .collect(),
            Bisector2::Parabola(p) => (-5..=5)
                .map(|i| {
                    let t = f64::from(i);
                    let frame = p.frame();
                    frame.origin()
                        + frame.x().vector() * (t * t / (4.0 * p.focal()))
                        + frame.y().vector() * t
                })
                .collect(),
            Bisector2::Ellipse(e) => (0..12)
                .map(|i| {
                    let t = core::f64::consts::TAU * f64::from(i) / 12.0;
                    let frame = e.frame();
                    frame.origin()
                        + frame.x().vector() * (e.major_radius() * t.cos())
                        + frame.y().vector() * (e.minor_radius() * t.sin())
                })
                .collect(),
            // Only the +x branch is the boundary bisector.
            Bisector2::Hyperbola(h) => (-3..=3)
                .map(|i| {
                    let t = 0.6 * f64::from(i);
                    let frame = h.frame();
                    frame.origin()
                        + frame.x().vector() * (h.major_radius() * t.cosh())
                        + frame.y().vector() * (h.minor_radius() * t.sinh())
                })
                .collect(),
        };
        for p in probes {
            let (da, db) = (a.distance_to(p), b.distance_to(p));
            // A hyperbola carries both branches; each point serves one side.
            assert!(
                (da - db).abs() < 1e-9,
                "not equidistant at {p:?}: {da} vs {db} for {bisector:?}"
            );
        }
    }

    #[test]
    fn bisectors_are_equidistant_loci() {
        let point = Target2::Point(Point2::new(1.0, 1.0));
        let other = Target2::Point(Point2::new(-1.0, 2.0));
        let line = Target2::Line(Axis2::new(Point2::new(0.0, -2.0), Direction2::X));
        let small = Target2::Circle(circle(0.0, 0.0, 5.0));
        let far = Target2::Circle(circle(12.0, 0.0, 2.0));

        assert_equidistant(&bisector(&point, &other, T).unwrap(), &point, &other);
        assert_equidistant(&bisector(&point, &line, T).unwrap(), &point, &line);
        // Point inside the circle: an ellipse.
        let inside = bisector(&point, &small, T).unwrap();
        assert!(matches!(inside, Bisector2::Ellipse(_)), "{inside:?}");
        assert_equidistant(&inside, &point, &small);
        // Unequal circles: a hyperbola.
        let between = bisector(&small, &far, T).unwrap();
        assert!(matches!(between, Bisector2::Hyperbola(_)), "{between:?}");
        assert_equidistant(&between, &small, &far);
        // Intersecting lines: the two angle bisectors.
        let slanted = Target2::Line(Axis2::new(
            Point2::new(0.0, -2.0),
            Direction2::new(Vector2::new(1.0, 1.0), T).unwrap(),
        ));
        let pair = bisector(&line, &slanted, T).unwrap();
        assert!(matches!(pair, Bisector2::Pair(_)), "{pair:?}");
        assert_equidistant(&pair, &line, &slanted);
    }
}
