//! The general surface/surface intersector: seed, then walk.
//!
//! Where two surfaces meet has no closed form in general: the curve is
//! transcendental, and `docs/DATA_MODEL.md` §9 is blunt about the consequence:
//! there is no exact answer to be exact about, which is why the topology carries
//! tolerances. What there is instead is a curve that can be *followed*, one
//! corrected step at a time, to a stated accuracy.
//!
//! # Two problems, kept apart on purpose
//!
//! **Finding a branch** and **following one** fail in completely different ways,
//! and lumping them together is how an intersector comes to look better than it
//! is. A tracer that follows one branch beautifully while never noticing the
//! second reports a smooth, accurate, *wrong* answer, and the obvious accuracy
//! measure, "is every point on both surfaces", scores it perfectly.
//!
//! So [`seeds`] and [`trace`] are separate, separately testable, and separately
//! measured. Seeding is polyhedral: both surfaces are sampled into triangles and
//! the triangle pairs that cross give starting points. A branch thinner than
//! the grid crosses no pair of triangles; where it meets a spline's border it
//! is found as a piercing, and a closed loop inside both surfaces is found
//! from the approach of the surfaces inside it. A branch running along a
//! border of either chart, where one surface ends in the other, is found by
//! sampling the border itself. Where the surfaces lie close
//! along a whole curve only part of that search is made, and there the grid
//! is the knob.
//!
//! # Following the curve
//!
//! At a point on both surfaces the intersection runs along the cross product of
//! the two normals: the one direction that stays in both tangent planes. Step
//! along it and you leave both surfaces slightly; a Newton correction brings you
//! back.
//!
//! The correction has four unknowns (two parameters on each surface) and three
//! equations, `A(u1,v1) = B(u2,v2)`. That is deliberately one short, because the
//! solution set *is* the curve and pinning it to a point needs one more
//! condition. The fourth is a plane across the direction of travel: it says how
//! far along to land, and it is what turns "somewhere on the curve" into "the
//! next point".
//!
//! # What it reports about itself
//!
//! Whether the curve closed, and whether it ran out of steps. A polyline that
//! stopped because it hit a limit is not the same answer as one that stopped
//! because the curve ended, and a caller that cannot tell them apart will treat
//! a truncated branch as a complete one.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface, SurfaceGeometry};
use ogeom_math::{Point, Vector, solve};

/// How hard to look, and how closely to follow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Marching {
    /// How far the polyline may sit from the true curve, in space.
    pub chord: f64,
    /// How finely each surface is sampled when looking for branches.
    ///
    /// A branch narrower than one cell crosses no pair of cells, and is found
    /// only from a border it pierces or from the approach of the surfaces
    /// near it. Raising this costs time quadratically and makes a thin branch
    /// cross cells.
    pub grid: usize,
    /// A ceiling on the points in one branch, so a curve that will not close
    /// cannot run forever.
    pub max_points: usize,
}

impl Default for Marching {
    fn default() -> Self {
        Self {
            chord: 1e-4,
            grid: 24,
            max_points: 20_000,
        }
    }
}

impl Marching {
    /// Check the settings are usable.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the chord is
    /// not positive, the grid is too coarse to hold a triangle, or no points are
    /// allowed.
    pub fn validate(&self) -> OgeomResult<()> {
        if !self.chord.is_finite() || self.chord <= 0.0 {
            ogeom_bail!(Construction, "a chord of {} is not a distance", self.chord);
        }
        if self.grid < 2 {
            ogeom_bail!(Construction, "a sampling grid needs at least two steps");
        }
        if self.max_points < 2 {
            ogeom_bail!(Construction, "a branch needs at least two points");
        }
        Ok(())
    }
}

/// A point that lies on both surfaces, with where it is on each.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// Parameters on the first surface.
    pub on_a: (f64, f64),
    /// Parameters on the second.
    pub on_b: (f64, f64),
    /// Where that is in space.
    pub point: Point,
}

/// Why a traced branch stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stopped {
    /// It came back to where it started.
    Closed,
    /// It reached the edge of one surface's domain.
    LeftTheDomain,
    /// The correction stopped converging: a tangency or a singular point.
    ///
    /// Reported rather than pushed through. Marching past a point where the two
    /// normals are parallel is how a tracer jumps onto the wrong branch, and a
    /// wrong branch is a plausible answer to a different question.
    Stalled,
    /// It hit [`Marching::max_points`].
    ///
    /// Distinct from every other reason, because this one means the answer is
    /// *incomplete* rather than finished.
    RanOut,
}

/// One traced branch.
#[derive(Debug, Clone, PartialEq)]
pub struct Traced {
    /// The points along it, in order.
    pub points: Vec<Point>,
    /// Where each point is on the first surface.
    pub on_a: Vec<(f64, f64)>,
    /// Where each point is on the second.
    pub on_b: Vec<(f64, f64)>,
    /// Why it stopped.
    pub stopped: Stopped,
}

impl Traced {
    /// Whether the branch is finished rather than truncated.
    #[must_use]
    pub const fn complete(&self) -> bool {
        !matches!(self.stopped, Stopped::RanOut)
    }

    /// Whether it closed on itself.
    #[must_use]
    pub const fn closed(&self) -> bool {
        matches!(self.stopped, Stopped::Closed)
    }
}

/// Starting points on the intersection, one per branch found.
///
/// Polyhedral: both surfaces are sampled into triangles, the pairs that cross
/// give approximate points, and each is corrected onto both surfaces exactly.
/// Points that land on the same spot are merged, so a branch crossing many
/// cells yields one seed rather than dozens. A spline's borders are pierced
/// by the other surface, and the pairs that come close without crossing are
/// searched for a crossing or for an approach a small loop runs round.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable.
pub fn seeds(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Contact>> {
    Ok(seeded(a, b, options, tol)?.0)
}

/// [`seeds`], and the stretches of either surface's borders lying in the
/// other, found on the same sampling: see [`border_stretches`].
fn seeded(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<(Vec<Contact>, Vec<Traced>)> {
    options.validate()?;
    let (mesh_a, mesh_b) = (sample(a, options.grid, tol), sample(b, options.grid, tol));
    let apart = span(a).min(span(b)) / f64::from(u32::try_from(options.grid).unwrap_or(1));
    // The bins reach as far as the surfaces can bow from their cells, so
    // the pairs that come close without crossing are found as well.
    let bow = |mesh: &[Cell]| mesh.iter().map(|c| c.sag).fold(0.0, f64::max);
    let near = CellBins::over(
        &mesh_b,
        BOW.mul_add(bow(&mesh_a) + bow(&mesh_b), options.chord),
    );
    let (facing_a, facing_b): (Vec<Facing>, Vec<Facing>) = (
        mesh_a.iter().map(Facing::of).collect(),
        mesh_b.iter().map(Facing::of).collect(),
    );

    let mut found: Vec<Seed> = Vec::new();
    // The pairs of cells that come close without crossing, and the grid
    // squares of the first surface some pair crosses in.
    let mut close: Vec<(usize, usize)> = Vec::new();
    let mut crossed = vec![false; mesh_a.len().div_ceil(2)];
    let mut candidates: Vec<usize> = Vec::new();
    for (i, cell_a) in mesh_a.iter().enumerate() {
        // Cheap rejection first: most pairs are nowhere near each other,
        // and the segment test below is far from free. The bins hand over
        // only the cells whose boxes could meet this one's, in their own
        // order, so what follows sees what the full scan would have.
        near.candidates(cell_a, &mut candidates);
        for &j in &candidates {
            let cell_b = &mesh_b[j];
            let crossing = overlap(cell_a, cell_b, options.chord)
                .then(|| triangles_cross(cell_a, cell_b))
                .flatten();
            let Some(guess) = crossing else {
                if may_touch(
                    (cell_a, &facing_a[i]),
                    (cell_b, &facing_b[j]),
                    options.chord,
                ) {
                    close.push((i, j));
                }
                continue;
            };
            crossed[i / 2] = true;
            let start = [cell_a.at.0, cell_a.at.1, cell_b.at.0, cell_b.at.1];
            let Some(contact) = correct(a, b, start, guess, None, options.chord, tol) else {
                continue;
            };
            // One seed per branch, not one per cell it passes through. The
            // reach is the *finer* surface's grid: the coarser surface's
            // cells say nothing about how far apart branches can be, and a
            // plane's clamped domain spans a million units.
            if merges(&found, contact, apart) {
                continue;
            }
            found.push(Seed::at(a, b, contact, tol));
        }
    }
    // A branch that runs in from a border at a grazing angle can be
    // thinner than the sampling's sag, and a short one cut off in a corner
    // of both domains can lie inside a single pair of cells without
    // crossing them: no pair of cells crosses on it. Where it meets the
    // border it is a curve piercing a surface, which is found exactly:
    // each border of each spline, and each rim of a drum or a cone, is
    // intersected with the other surface, and every piercing seeds.
    for (from_a, border_of, other) in [(true, a, b), (false, b, a)] {
        let borders = spline_borders(border_of, tol)
            .into_iter()
            .chain(rims(border_of, tol));
        for (border, at) in borders {
            let Ok(met) = crate::intersect_curve_surface(
                &border,
                other,
                crate::CurveSurfaceOptions::default(),
                tol,
            ) else {
                continue;
            };
            for piercing in met.crossings {
                let on_border = at(piercing.on_curve);
                let start = if from_a {
                    [
                        on_border.0,
                        on_border.1,
                        piercing.on_surface.0,
                        piercing.on_surface.1,
                    ]
                } else {
                    [
                        piercing.on_surface.0,
                        piercing.on_surface.1,
                        on_border.0,
                        on_border.1,
                    ]
                };
                let Some(contact) = correct(a, b, start, piercing.point, None, options.chord, tol)
                else {
                    continue;
                };
                if merges(&found, contact, apart) {
                    continue;
                }
                found.push(Seed::at(a, b, contact, tol));
            }
        }
    }
    // A closed loop inside both surfaces, smaller than a cell and touching
    // no border, crosses no pair of cells and pierces no border. Inside it
    // the separation of the surfaces has an extremum, where the normals are
    // collinear and the gap is small; the cells that come close without
    // crossing seed a search for it (Sederberg and Meyers 1988), and the
    // loop is found by stepping out from it until the separation changes
    // sign.
    let extra = extremum_seeds(a, b, (&mesh_a, &mesh_b), &close, &crossed, options, tol);
    for contact in extra {
        if !merges(&found, contact, apart) {
            found.push(Seed::at(a, b, contact, tol));
        }
    }
    let lying = border_stretches(a, b, (&mesh_a, &mesh_b), options, tol)?;
    Ok((found.into_iter().map(|s| s.contact).collect(), lying))
}

/// How many times a cell's sag the surface is allowed to bow from the
/// cell when deciding whether two cells could touch: the sag is measured
/// at the cell's middle only.
const BOW: f64 = 2.0;

/// The sine of the angle between two cells' normals that always passes:
/// a floor under the turning estimated from the sags.
const FACING: f64 = 0.02;

/// Whether two cells that do not cross could still hold a loop between
/// them: their boxes, grown by how far each surface bows from its cell,
/// meet; each cell's corners straddle or come within that bow of the
/// other's plane; and their normals are near enough to collinear for a
/// stationary approach, where the surfaces' normals agree, to lie between
/// them.
fn may_touch(a: (&Cell, &Facing), b: (&Cell, &Facing), chord: f64) -> bool {
    let ((a, facing_a), (b, facing_b)) = (a, b);
    let slack = BOW.mul_add(a.sag + b.sag, chord);
    if !overlap(a, b, slack) {
        return false;
    }
    let (Some(na), Some(nb)) = (facing_a.normal, facing_b.normal) else {
        return false;
    };
    if na.cross(nb).magnitude() > FACING + facing_a.turn + facing_b.turn {
        return false;
    }
    let within = |corners: &[Point; 3], origin: Point, n: Vector| {
        let side = corners.map(|p| (p - origin).dot(n));
        let low = side.iter().copied().fold(f64::INFINITY, f64::min);
        let high = side.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        low <= slack && high >= -slack
    };
    within(&a.corners, b.corners[0], nb) && within(&b.corners, a.corners[0], na)
}

/// Which way a cell faces, and how far the surface's normal strays from
/// the cell's within it.
struct Facing {
    normal: Option<Vector>,
    turn: f64,
}

impl Facing {
    fn of(cell: &Cell) -> Self {
        // An arc of chord `l` and sag `s` turns through about `8 s / l`.
        Self {
            normal: cell_normal(cell),
            turn: 8.0 * BOW * cell.sag / cell_size(cell).max(f64::MIN_POSITIVE),
        }
    }
}

/// A cell's unit normal, `None` for a cell with no area.
fn cell_normal(cell: &Cell) -> Option<Vector> {
    let [p, q, r] = cell.corners;
    let n = (q - p).cross(r - p);
    let length = n.magnitude();
    (length > f64::MIN_POSITIVE).then(|| n * (1.0 / length))
}

/// How near two cells come, for ranking: the nearer of each one's middle
/// to the other's triangle.
fn closeness(a: &Cell, b: &Cell) -> f64 {
    point_to_triangle(centroid(a), b.corners).min(point_to_triangle(centroid(b), a.corners))
}

/// The distance from a point to a triangle.
fn point_to_triangle(p: Point, t: [Point; 3]) -> f64 {
    let edges = (0..3).map(|k| distance_to_segment(p, t[k], t[(k + 1) % 3]));
    let rim = edges.fold(f64::INFINITY, f64::min);
    let n = (t[1] - t[0]).cross(t[2] - t[0]);
    let area = n.magnitude();
    if area <= f64::MIN_POSITIVE {
        return rim;
    }
    let n = n * (1.0 / area);
    let foot = p - n * (p - t[0]).dot(n);
    let inside = (0..3).all(|k| {
        let (from, to) = (t[k], t[(k + 1) % 3]);
        (to - from).cross(foot - from).dot(n) >= 0.0
    });
    if inside { p.distance(foot) } else { rim }
}

/// The middle of a cell.
fn centroid(cell: &Cell) -> Point {
    let [p, q, r] = cell.corners;
    Point::new(
        (p.x + q.x + r.x) / 3.0,
        (p.y + q.y + r.y) / 3.0,
        (p.z + q.z + r.z) / 3.0,
    )
}

/// Whether a parameter lies in the grid square a cell was cut from (the
/// box of its corners' parameters, which takes in the square's diagonal),
/// or within `around` squares of it.
fn in_square(cell: &Cell, at: (f64, f64), around: f64) -> bool {
    let (u, v) = cell.params.iter().fold(
        (
            (f64::INFINITY, f64::NEG_INFINITY),
            (f64::INFINITY, f64::NEG_INFINITY),
        ),
        |((u0, u1), (v0, v1)), &(u, v)| ((u0.min(u), u1.max(u)), (v0.min(v), v1.max(v))),
    );
    let (du, dv) = ((u.1 - u.0) * around, (v.1 - v.0) * around);
    (u.0 - du..=u.1 + du).contains(&at.0) && (v.0 - dv..=v.1 + dv).contains(&at.1)
}

/// A cell's longest edge.
fn cell_size(cell: &Cell) -> f64 {
    let [p, q, r] = cell.corners;
    p.distance(q).max(q.distance(r)).max(r.distance(p))
}

/// Contacts found near the cell pairs `close`, where the surfaces come
/// close without their cells crossing: a point of a crossing the grid
/// stepped over, or a point of a loop around a stationary approach where
/// the surfaces pass through each other. `crossed` marks the first
/// surface's grid squares some pair of cells crosses in.
fn extremum_seeds(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    (mesh_a, mesh_b): (&[Cell], &[Cell]),
    close: &[(usize, usize)],
    crossed: &[bool],
    options: Marching,
    tol: Tolerances,
) -> Vec<Contact> {
    let middle = |c: &Cell| {
        let [p, q, r] = c.params;
        ((p.0 + q.0 + r.0) / 3.0, (p.1 + q.1 + r.1) / 3.0)
    };
    // One search per grid square of the first surface, from the cell of
    // the second nearest it: the two triangles of a square share it, and
    // a square holds one approach unless the second surface folds back
    // within a cell.
    let mut nearest: std::collections::BTreeMap<usize, (f64, usize, usize)> =
        std::collections::BTreeMap::new();
    for &(i, j) in close {
        let d = closeness(&mesh_a[i], &mesh_b[j]);
        let best = nearest.entry(i / 2).or_insert((d, i, j));
        if d < best.0 {
            *best = (d, i, j);
        }
    }
    // The squares no found branch crosses in go first, then the nearest.
    let mut order: Vec<(f64, usize, usize)> = nearest.into_values().collect();
    order.sort_by(|x, y| {
        crossed[x.1 / 2]
            .cmp(&crossed[y.1 / 2])
            .then(x.0.total_cmp(&y.0))
            .then(x.1.cmp(&y.1))
    });
    let mut out: Vec<Contact> = Vec::new();
    let mut stationary: Vec<Contact> = Vec::new();
    // Where a search found nothing to follow. Where the surfaces lie close
    // along a whole curve (a blend against the face it rolls on) every
    // square is close, and a search in each would cost more than the grid
    // did; the squares around such a search are passed over.
    let mut settled: Vec<(f64, f64)> = Vec::new();
    for (_, i, j) in order {
        let (cell_a, cell_b) = (&mesh_a[i], &mesh_b[j]);
        if settled.iter().any(|&at| in_square(cell_a, at, 1.0)) {
            continue;
        }
        // A pair of squares an earlier search already landed in has its
        // approach.
        if stationary
            .iter()
            .any(|known| in_square(cell_a, known.on_a, 0.0) && in_square(cell_b, known.on_b, 0.0))
        {
            continue;
        }
        let (from_a, from_b) = (middle(cell_a), middle(cell_b));
        // The surfaces meeting near the squares answer first, and cheaply:
        // a point of a crossing the grid stepped over seeds, and a touch
        // has nothing to follow. Only where they do not meet is the
        // approach between them looked for.
        let guess = centroid(cell_a).midpoint(centroid(cell_b));
        if let Some(contact) = correct(
            a,
            b,
            [from_a.0, from_a.1, from_b.0, from_b.1],
            guess,
            None,
            options.chord,
            tol,
        ) {
            if passes_through(a, b, contact, tol) {
                out.push(contact);
            } else {
                settled.push(contact.on_a);
            }
            continue;
        }
        let Some((ua, va, ub, vb)) =
            crate::extrema::stationary_surface_surface(a, b, from_a, from_b, STATIONARY_STEPS, tol)
        else {
            settled.push(from_a);
            continue;
        };
        let (Ok(pa), Ok(pb)) = (a.point_at(ua, va, tol), b.point_at(ub, vb, tol)) else {
            continue;
        };
        let reach = tol.confusion() * 100.0;
        if stationary.iter().any(|known| {
            known.point.distance(pa) <= reach
                && b.point_at(known.on_b.0, known.on_b.1, tol)
                    .is_ok_and(|q| q.distance(pb) <= reach)
        }) {
            continue;
        }
        stationary.push(Contact {
            on_a: (ua, va),
            on_b: (ub, vb),
            point: pa,
        });
        // A loop thinner than the cells passes through the other surface
        // by no more than the cells bow.
        let gap = pa.distance(pb);
        if gap > BOW.mul_add(cell_a.sag + cell_b.sag, options.chord) {
            continue;
        }
        let at = Contact {
            on_a: (ua, va),
            on_b: (ub, vb),
            point: pa.midpoint(pb),
        };
        if gap <= tol.confusion() {
            // On both surfaces already: the least distance is nought on the
            // intersection itself. A point of a crossing the grid stepped
            // over seeds; a touch has nothing to follow.
            if !passes_through(a, b, at, tol) {
                settled.push(at.on_a);
            } else if let Some(contact) =
                correct(a, b, [ua, va, ub, vb], at.point, None, options.chord, tol)
            {
                out.push(contact);
            }
            continue;
        }
        // Where the normals agree and the gap is within the chord, it is a
        // touch at the accuracy the trace is held to: any loop about it is
        // shallower than a traced polyline's own error.
        if gap <= options.chord {
            settled.push(at.on_a);
            continue;
        }
        let size = cell_size(cell_a).max(cell_size(cell_b));
        if let Some(contact) = loop_around(a, b, at, size, options.chord, tol) {
            out.push(contact);
        }
    }
    out
}

/// How many Newton steps the search for a stationary approach takes from
/// the middle of a square. From that near, a regular approach converges in
/// a handful; one that has not by then lies where the surfaces are all but
/// tangent along a curve, and has no isolated approach to find.
const STATIONARY_STEPS: usize = 12;

/// The signed separation of the first surface from the second at `on_a`:
/// the gap to its nearest point on the second, along the second's normal
/// there. `near` is where on the second to start looking, and is moved to
/// what was found.
fn separation(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    near: &mut (f64, f64),
    tol: Tolerances,
) -> Option<f64> {
    let p = a.point_at(on_a.0, on_a.1, tol).ok()?;
    let (on_b, q) = nearest_on(b, *near, p, tol)?;
    let n = normal_at(b, on_b, tol)?;
    *near = on_b;
    Some((p - q).dot(n))
}

/// The parameter step that moves a unit along the tangent `t`, for a
/// surface whose derivatives are `du` and `dv`: the first fundamental form
/// inverted. `None` where the derivatives are parallel.
fn params_along(du: Vector, dv: Vector, t: Vector) -> Option<(f64, f64)> {
    let (g11, g12, g22) = (du.dot(du), du.dot(dv), dv.dot(dv));
    let det = g11.mul_add(g22, -(g12 * g12));
    if det <= f64::MIN_POSITIVE {
        return None;
    }
    let (x, y) = (du.dot(t), dv.dot(t));
    Some((
        g22.mul_add(x, -(g12 * y)) / det,
        g11.mul_add(y, -(g12 * x)) / det,
    ))
}

/// Whether the surfaces pass through each other at a point on both: the
/// separation has opposite signs, each beyond the confusion distance, a
/// little to either side of the intersection on the first surface. At a
/// crossing the separation grows with the sine of the angle between the
/// normals; where the surfaces only touch, rounding leaves the normals a
/// whisker apart, and the separation keeps one sign on both sides.
fn passes_through(a: &SurfaceGeometry, b: &SurfaceGeometry, at: Contact, tol: Tolerances) -> bool {
    let Some(along) = tangent_at(a, b, at, tol) else {
        return false;
    };
    let (Some(na), Ok((du, dv))) = (
        normal_at(a, at.on_a, tol),
        a.d1_at(at.on_a.0, at.on_a.1, tol),
    ) else {
        return false;
    };
    let sine = crossing_sine(a, b, at.on_a, at.on_b, tol);
    let Some(d) = (sine > f64::MIN_POSITIVE)
        .then(|| params_along(du, dv, na.cross(along)))
        .flatten()
    else {
        return false;
    };
    let h = 4.0 * tol.confusion() / sine;
    let side = |s: f64| {
        let (u, v) = clamp(a, d.0.mul_add(s, at.on_a.0), d.1.mul_add(s, at.on_a.1));
        let mut near = at.on_b;
        separation(a, b, (u, v), &mut near, tol)
    };
    match (side(h), side(-h)) {
        (Some(p), Some(q)) => {
            p.abs() > tol.confusion() && q.abs() > tol.confusion() && p.signum() != q.signum()
        }
        _ => false,
    }
}

/// A point of the loop around a stationary approach `at` where the
/// surfaces pass through each other, if there is one within twice `size`
/// of it.
///
/// Steps out from the approach on the first surface along eight directions
/// of its tangent plane, at distances doubling from a thousandth of `size`,
/// until the separation changes sign; the sign change is halved down and
/// corrected onto both surfaces.
fn loop_around(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Contact,
    size: f64,
    within: f64,
    tol: Tolerances,
) -> Option<Contact> {
    const DIRECTIONS: u32 = 8;
    const DOUBLINGS: i32 = 11;
    const HALVINGS: usize = 40;
    let mut near = at.on_b;
    let here = separation(a, b, at.on_a, &mut near, tol)?;
    if here == 0.0 {
        return None;
    }
    let (du, dv) = a.d1_at(at.on_a.0, at.on_a.1, tol).ok()?;
    let normal = unit_normal(du, dv, tol)?;
    let e1 = du * (1.0 / du.magnitude());
    let e2 = normal.cross(e1);
    let step = |d: (f64, f64), r: f64| {
        let (u, v) = (d.0.mul_add(r, at.on_a.0), d.1.mul_add(r, at.on_a.1));
        let held = clamp(a, u, v);
        let wrapped = |x: f64, y: f64, periodic: bool| periodic || x == y;
        (wrapped(u, held.0, a.is_periodic_u()) && wrapped(v, held.1, a.is_periodic_v()))
            .then_some(held)
    };
    for k in 0..DIRECTIONS {
        let angle = core::f64::consts::TAU * f64::from(k) / f64::from(DIRECTIONS);
        let d = params_along(du, dv, e1 * angle.cos() + e2 * angle.sin())?;
        let mut near = at.on_b;
        let (mut inside, mut r) = (0.0, size * 2f64.powi(-DOUBLINGS + 1));
        let mut across = None;
        while r <= 2.0 * size {
            let Some(uv) = step(d, r) else { break };
            let Some(there) = separation(a, b, uv, &mut near, tol) else {
                break;
            };
            if there.signum() != here.signum() {
                across = Some(r);
                break;
            }
            inside = r;
            r *= 2.0;
        }
        let Some(mut outside) = across else {
            continue;
        };
        for _ in 0..HALVINGS {
            let mid = f64::midpoint(inside, outside);
            let Some(uv) = step(d, mid) else { break };
            let Some(there) = separation(a, b, uv, &mut near, tol) else {
                break;
            };
            if there.signum() == here.signum() {
                inside = mid;
            } else {
                outside = mid;
            }
        }
        let Some(uv) = step(d, outside) else {
            continue;
        };
        let Ok(guess) = a.point_at(uv.0, uv.1, tol) else {
            continue;
        };
        if let Some(contact) = correct(a, b, [uv.0, uv.1, near.0, near.1], guess, None, within, tol)
        {
            return Some(contact);
        }
    }
    None
}

/// A seed kept, with the direction the intersection runs through it.
struct Seed {
    contact: Contact,
    along: Option<Vector>,
}

impl Seed {
    fn at(a: &SurfaceGeometry, b: &SurfaceGeometry, contact: Contact, tol: Tolerances) -> Self {
        Self {
            contact,
            along: tangent_at(a, b, contact, tol),
        }
    }
}

/// Whether a contact is a seed already kept, met again further along its
/// branch: within `apart` of it, and off the line the intersection runs
/// through it by less than a quarter of the way there. Two branches side by
/// side, or a loop beside another the same drum cuts, lie across that line
/// rather than along it, and each keeps a seed of its own however close.
/// Where the surfaces are all but tangent the intersection has no direction
/// to ask, and nearness alone decides.
fn merges(found: &[Seed], contact: Contact, apart: f64) -> bool {
    found.iter().any(|seed| {
        let offset = contact.point - seed.contact.point;
        let gap = offset.magnitude();
        gap <= apart
            && seed
                .along
                .is_none_or(|along| offset.cross(along).magnitude() <= 0.25 * gap)
    })
}

/// Past this a chart's bound is a stand-in for no end at all.
const FAR: f64 = 1e6;

/// The sine of the angle between two surfaces' normals at or below which
/// they are taken as tangent where they meet: a branch along which they
/// stay that close to tangent is a contact, not a crossing.
pub(crate) const TANGENTIAL_SINE: f64 = 3e-2;

/// Whether a surface meeting another along a curve the two share leaves
/// it there at an angle, rather than touching it: `sine` is the sine
/// between their normals on the curve, `probe` the distance in from the
/// curve to the nearest sample that says which side of the other the
/// surface keeps, `standing` how far that sample stands off the other, and
/// `width` the distance from the other within which a sample says nothing.
///
/// Leaving at a sine `s`, the surface rises `s d` a distance `d` in; it
/// also bows, by what is left of `standing`. It leaves at an angle where
/// its rise is resolved at the probe, `s probe > width`, and is no less
/// than its bow there, `s probe ≥ standing - s probe`: the side it keeps
/// is then the one its angle takes it to. One touching the other has a
/// sine of naught on the curve, or the little a fit leaves, and stands off
/// it by its bow, so it fails one or the other. Past [`TANGENTIAL_SINE`],
/// the sine below which two surfaces are taken as tangent wherever they
/// meet, it always leaves.
pub(crate) fn leaves_at_an_angle(sine: f64, probe: f64, standing: f64, width: f64) -> bool {
    sine > TANGENTIAL_SINE || (sine * probe > width && 2.0 * sine * probe >= standing)
}

/// A border of a surface as a curve, and the map from the curve's
/// parameter to the surface's.
type Border = (ogeom_geom::Curve, Box<dyn Fn(f64) -> (f64, f64)>);

/// The open borders of a spline surface, or of a window of one.
fn spline_borders(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Border> {
    let spline = match surface {
        SurfaceGeometry::BSpline(spline) => spline,
        SurfaceGeometry::Trimmed(trimmed) => match trimmed.basis() {
            SurfaceGeometry::BSpline(spline) => spline,
            _ => return Vec::new(),
        },
        _ => return Vec::new(),
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let whole = ogeom_geom::Surface::domain(spline);
    // An iso curve runs the spline's whole other direction; on a window it
    // is cut to the window's.
    let cut = |c: ogeom_geom::BSplineCurve, (lo, hi): (f64, f64), (a, b): (f64, f64)| {
        let c = ogeom_geom::Curve::BSpline(c);
        if lo <= a && hi >= b {
            return Some(c);
        }
        ogeom_geom::TrimmedCurve::new(c, lo, hi, tol)
            .ok()
            .map(Into::into)
    };
    let mut out: Vec<Border> = Vec::new();
    if !surface.is_closed_u(tol) {
        for u in [u0, u1] {
            if let Some(c) = spline
                .iso_u_curve(u, tol)
                .ok()
                .and_then(|c| cut(c, (v0, v1), whole.1))
            {
                out.push((c, Box::new(move |t| (u, t))));
            }
        }
    }
    if !surface.is_closed_v(tol) {
        for v in [v0, v1] {
            if let Some(c) = spline
                .iso_v_curve(v, tol)
                .ok()
                .and_then(|c| cut(c, (u0, u1), whole.0))
            {
                out.push((c, Box::new(move |t| (t, v))));
            }
        }
    }
    out
}

/// The rims of a drum or a cone: the circles at the ends of its height,
/// where those ends are a part's height and not a chart's clamp. On a
/// window of one less than a turn wide, the arcs of those circles the
/// window spans.
fn rims(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Border> {
    let is_drum = |s: &SurfaceGeometry| {
        matches!(s, SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)) && s.is_closed_u(tol)
    };
    let (drum, windowed) = match surface {
        SurfaceGeometry::Trimmed(trimmed) if is_drum(trimmed.basis()) => (trimmed.basis(), true),
        _ if is_drum(surface) => (surface, false),
        _ => return Vec::new(),
    };
    let turn = drum.domain().0.1 - drum.domain().0.0;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let mut out: Vec<Border> = Vec::new();
    for v in [v0, v1] {
        if !v.is_finite() || v.abs() > FAR {
            continue;
        }
        // Three points a third of a turn apart: the circle through them
        // starts at the first and runs the way the angle does, so its
        // parameter is the angle from the window's start. A cone's apex is
        // no circle, and is passed over.
        let at = |k: f64| drum.point_at(turn.mul_add(k / 3.0, u0), v, tol);
        let (Ok(p0), Ok(p1), Ok(p2)) = (at(0.0), at(1.0), at(2.0)) else {
            continue;
        };
        let Ok(circle) = ogeom_math::Circle::through(p0, p1, p2, tol) else {
            continue;
        };
        let circle: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
        let rim = if windowed {
            let Ok(arc) = ogeom_geom::TrimmedCurve::new(circle, 0.0, u1 - u0, tol) else {
                continue;
            };
            arc.into()
        } else {
            circle
        };
        out.push((rim, Box::new(move |t| (u0 + t, v))));
    }
    out
}

/// A border of a surface's chart: the parameter it holds (`u` where
/// `holds_u`), the value it holds it at, the range of the other parameter
/// along it, and whether it closes on itself.
#[derive(Debug, Clone, Copy)]
struct ChartBorder {
    holds_u: bool,
    at: f64,
    range: (f64, f64),
    closed: bool,
}

impl ChartBorder {
    /// The surface's parameters at `t` along the border.
    const fn params(&self, t: f64) -> (f64, f64) {
        if self.holds_u {
            (self.at, t)
        } else {
            (t, self.at)
        }
    }

    /// The surface's point at `t` along the border, and the derivative
    /// along it there.
    fn point_d1(
        &self,
        surface: &SurfaceGeometry,
        t: f64,
        tol: Tolerances,
    ) -> Option<(Point, Vector)> {
        let (u, v) = self.params(t);
        let (p, du, dv) = surface.point_d1_at(u, v, tol).ok()?;
        Some((p, if self.holds_u { dv } else { du }))
    }
}

/// The borders of a surface's chart in each direction that neither repeats
/// nor closes, where the bound and the range along it are a part's and not
/// a chart's stand-in for no end.
fn chart_borders(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<ChartBorder> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let real = |x: f64| x.is_finite() && x.abs() <= FAR;
    let open_u = !surface.is_periodic_u() && !surface.is_closed_u(tol);
    let open_v = !surface.is_periodic_v() && !surface.is_closed_v(tol);
    let mut out = Vec::new();
    for (holds_u, open, ends, range, closed) in [
        (true, open_u, (u0, u1), (v0, v1), !open_v),
        (false, open_v, (v0, v1), (u0, u1), !open_u),
    ] {
        if !open || !real(range.0) || !real(range.1) || range.1 <= range.0 {
            continue;
        }
        for at in [ends.0, ends.1] {
            if real(at) {
                out.push(ChartBorder {
                    holds_u,
                    at,
                    range,
                    closed,
                });
            }
        }
    }
    out
}

/// The stretches of either surface's chart borders that lie in the other
/// surface and leave it at an angle, each as a traced branch.
///
/// A border lying in the other surface is a curve of the section, and one
/// neither seeding finds: no pair of cells crosses on it (the cells beside
/// it lie on one side of the other surface, or meet it only along an
/// edge), the border pierces the other surface nowhere, and a walk along
/// it stops at once, as every step leaves the chart. So each border is
/// sampled at the grid's spacing. A sample lies in the other surface where
/// a cell of the other comes within the distance the surface bows from its
/// cells, and the foot projected from that cell is within the confusion
/// distance, the distance a corrected seed is held to. Two neighbouring
/// samples in it make a stretch, carried by bisection out to where the
/// border leaves the other surface, then stepped along at the chord as a
/// walk is, every point projected again. A stretch stands only where every
/// point lies in the other surface and the surface leaves the other there
/// at an angle, as [`leaves_at_an_angle`] has it, against the grid's next
/// row in and the confusion distance: a border the surface leaves
/// tangentially, as a blend's edge leaves its host, is a contact, and is
/// left to the marcher and the contact handling. A stretch both surfaces'
/// borders hold is stated once.
fn border_stretches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    (mesh_a, mesh_b): (&[Cell], &[Cell]),
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Traced>> {
    let reach = options.chord.max(tol.confusion()) * 8.0;
    let mut out: Vec<Traced> = Vec::new();
    for (from_a, border_of, other, cells) in [(true, a, b, mesh_b), (false, b, a, mesh_a)] {
        let borders = chart_borders(border_of, tol);
        if borders.is_empty() || cells.is_empty() {
            continue;
        }
        for border in borders {
            for piece in stretches(border_of, other, cells, border, options, tol)? {
                let piece = if from_a {
                    piece
                } else {
                    Traced {
                        on_a: piece.on_b,
                        on_b: piece.on_a,
                        ..piece
                    }
                };
                let n = piece.points.len();
                let middle = if n > 2 {
                    piece.points[n / 2]
                } else {
                    piece.points[0].midpoint(piece.points[n - 1])
                };
                if out.iter().any(|other| passes_near(other, middle, reach)) {
                    continue;
                }
                out.push(piece);
            }
        }
    }
    Ok(out)
}

/// The stretches of `border` of `surface` lying in `other`, sampled into
/// `cells`, as traced branches whose first parameters are on `surface`:
/// see [`border_stretches`].
fn stretches(
    surface: &SurfaceGeometry,
    other: &SurfaceGeometry,
    cells: &[Cell],
    border: ChartBorder,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Traced>> {
    let grid = options.grid;
    let (lo, hi) = border.range;
    #[allow(clippy::cast_precision_loss)]
    let at = |k: usize| (hi - lo).mul_add(k as f64 / grid as f64, lo);
    let samples: Vec<Option<Point>> = (0..=grid)
        .map(|k| {
            let (u, v) = border.params(at(k));
            surface.point_at(u, v, tol).ok()
        })
        .collect();
    let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in samples.iter().flatten() {
        low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
        high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
    }
    // A border that is a single point (a pole, an apex) has no length to
    // lie along.
    if samples.iter().flatten().count() < 2 || (high - low).magnitude() <= tol.confusion() {
        return Ok(Vec::new());
    }
    // The cells the samples could lie on: those whose boxes come within
    // the distance the surface bows from them of the samples' box.
    let near: Vec<usize> = (0..cells.len())
        .filter(|&j| {
            let (cell, reach) = (&cells[j], BOW.mul_add(cells[j].sag, options.chord));
            cell.low.x <= high.x + reach
                && low.x <= cell.high.x + reach
                && cell.low.y <= high.y + reach
                && low.y <= cell.high.y + reach
                && cell.low.z <= high.z + reach
                && low.z <= cell.high.z + reach
        })
        .collect();
    if near.is_empty() {
        return Ok(Vec::new());
    }
    let lying = Lying {
        surface,
        other,
        cells,
        near,
        border,
        options,
        tol,
    };
    let feet: Vec<Option<(f64, f64)>> = samples
        .iter()
        .map(|p| {
            let p = (*p)?;
            lying.foot(p, lying.nearest_cell(p)?)
        })
        .collect();
    lying.stretches(&feet)
}

/// One border of a surface asked against another surface, and the cells of
/// the other near it: see [`border_stretches`].
struct Lying<'s> {
    surface: &'s SurfaceGeometry,
    other: &'s SurfaceGeometry,
    cells: &'s [Cell],
    near: Vec<usize>,
    border: ChartBorder,
    options: Marching,
    tol: Tolerances,
}

impl Lying<'_> {
    /// The stretches lying in the other surface, from the feet of the
    /// border's samples there.
    fn stretches(&self, feet: &[Option<(f64, f64)>]) -> OgeomResult<Vec<Traced>> {
        let border = self.border;
        let grid = feet.len() - 1;
        let (lo, hi) = border.range;
        #[allow(clippy::cast_precision_loss)]
        let at = |k: usize| (hi - lo).mul_add(k as f64 / grid as f64, lo);
        let mut out = Vec::new();
        if border.closed && feet.iter().all(Option::is_some) {
            if let Some(mut whole) = self.walked((lo, feet[0]), (hi, feet[grid]))? {
                whole.points.pop();
                whole.on_a.pop();
                whole.on_b.pop();
                whole.points.push(whole.points[0]);
                whole.on_a.push(whole.on_a[0]);
                whole.on_b.push(whole.on_b[0]);
                whole.stopped = Stopped::Closed;
                out.push(whole);
            }
            return Ok(out);
        }
        let mut k = 0;
        while k < grid {
            if feet[k].is_none() || feet[k + 1].is_none() {
                k += 1;
                continue;
            }
            let mut end = k + 1;
            while end < grid && feet[end + 1].is_some() {
                end += 1;
            }
            let start = if k == 0 {
                (lo, feet[0])
            } else {
                self.leaves((at(k), feet[k]), at(k - 1))
            };
            let finish = if end == grid {
                (hi, feet[grid])
            } else {
                self.leaves((at(end), feet[end]), at(end + 1))
            };
            if let Some(piece) = self.walked(start, finish)? {
                out.push(piece);
            }
            k = end + 1;
        }
        Ok(out)
    }

    /// The parameters of the cell corner of the other surface nearest `p`,
    /// among the cells it could lie on: those within the distance the
    /// surface bows from them.
    fn nearest_cell(&self, p: Point) -> Option<(f64, f64)> {
        let (_, cell) = self
            .near
            .iter()
            .map(|&j| &self.cells[j])
            .filter_map(|cell| {
                let reach = BOW.mul_add(cell.sag, self.options.chord);
                let boxed = (cell.low.x - reach..=cell.high.x + reach).contains(&p.x)
                    && (cell.low.y - reach..=cell.high.y + reach).contains(&p.y)
                    && (cell.low.z - reach..=cell.high.z + reach).contains(&p.z);
                if !boxed {
                    return None;
                }
                let d = point_to_triangle(p, cell.corners);
                (d <= reach).then_some((d, cell))
            })
            .min_by(|x, y| x.0.total_cmp(&y.0))?;
        let corner = (0..3).min_by(|&i, &j| {
            p.square_distance(cell.corners[i])
                .total_cmp(&p.square_distance(cell.corners[j]))
        })?;
        Some(cell.params[corner])
    }

    /// The foot of `p` on the other surface, from `guess`, where it lies
    /// within the confusion distance.
    fn foot(&self, p: Point, guess: (f64, f64)) -> Option<(f64, f64)> {
        let found = ogeom_geom::project_on_surface_from(self.other, p, guess, self.tol).ok()?;
        (found.distance <= self.tol.confusion()).then_some(found.parameters)
    }

    /// The foot at `t` along the border, from `guess` and failing that
    /// from the nearest cell.
    fn foot_at(&self, t: f64, guess: Option<(f64, f64)>) -> Option<(Point, (f64, f64))> {
        let (u, v) = self.border.params(t);
        let p = self.surface.point_at(u, v, self.tol).ok()?;
        let found = guess
            .and_then(|g| self.foot(p, g))
            .or_else(|| self.foot(p, self.nearest_cell(p)?))?;
        Some((p, found))
    }

    /// Where the border leaves the other surface between `inside`, whose
    /// foot is known, and `outside`: bisected down to the parametric
    /// confusion.
    fn leaves(&self, inside: (f64, Option<(f64, f64)>), outside: f64) -> (f64, Option<(f64, f64)>) {
        let (mut on, mut foot) = inside;
        let mut off = outside;
        for _ in 0..64 {
            if (off - on).abs() <= self.tol.parametric() {
                break;
            }
            let mid = f64::midpoint(on, off);
            match self.foot_at(mid, foot) {
                Some((_, found)) => {
                    on = mid;
                    foot = Some(found);
                }
                None => off = mid,
            }
        }
        (on, foot)
    }

    /// The surface's sampling grid's next row in from the border at `t`,
    /// where the border stands at `p` with its foot `foot` on the other:
    /// how far that row stands from `p`, and from the other surface. It is
    /// the nearest sample saying which side of the other the surface keeps.
    fn probe(&self, t: f64, p: Point, foot: (f64, f64)) -> Option<(f64, f64)> {
        let ((u0, u1), (v0, v1)) = self.surface.domain();
        let (lo, hi) = if self.border.holds_u {
            (u0, u1)
        } else {
            (v0, v1)
        };
        #[allow(clippy::cast_precision_loss)]
        let step = (hi - lo) / self.options.grid as f64;
        let at = self.border.at;
        let inside = if (at - lo).abs() <= (hi - at).abs() {
            at + step
        } else {
            at - step
        };
        let (u, v) = self.border.params(t);
        let (u, v) = if self.border.holds_u {
            (inside, v)
        } else {
            (u, inside)
        };
        let q = self.surface.point_at(u, v, self.tol).ok()?;
        let off = ogeom_geom::project_on_surface_from(self.other, q, foot, self.tol).ok()?;
        Some((q.distance(p), off.distance))
    }

    /// The border from `start` to `finish` stepped along at the chord as a
    /// walk is, `None` unless every point lies in the other surface and
    /// the two leave each other at an angle there.
    fn walked(
        &self,
        start: (f64, Option<(f64, f64)>),
        finish: (f64, Option<(f64, f64)>),
    ) -> OgeomResult<Option<Traced>> {
        let (border, chord, tol) = (self.border, self.options.chord, self.tol);
        let ((t0, foot0), (t1, _)) = (start, finish);
        if t1 <= t0 {
            return Ok(None);
        }
        let extent = span(self.surface).max(span(self.other));
        let ceiling = extent / 8.0;
        let mut step = (chord * extent).sqrt().clamp(tol.confusion(), ceiling);
        let mut traced = Traced {
            points: Vec::new(),
            on_a: Vec::new(),
            on_b: Vec::new(),
            stopped: Stopped::LeftTheDomain,
        };
        let keep = |traced: &mut Traced, t: f64, guess: Option<(f64, f64)>| -> Option<(f64, f64)> {
            let (p, foot) = self.foot_at(t, guess)?;
            let ns = normal_at(self.surface, border.params(t), tol)?;
            let no = normal_at(self.other, foot, tol)?;
            let sine = ns.cross(no).magnitude();
            if sine <= TANGENTIAL_SINE {
                let (probe, standing) = self.probe(t, p, foot)?;
                if !leaves_at_an_angle(sine, probe, standing, tol.confusion()) {
                    return None;
                }
            }
            traced.points.push(p);
            traced.on_a.push(border.params(t));
            traced.on_b.push(foot);
            Some(foot)
        };
        let (Some(mut foot), Some((mut here, mut heading))) = (
            keep(&mut traced, t0, foot0),
            border.point_d1(self.surface, t0, tol),
        ) else {
            return Ok(None);
        };
        let mut t = t0;
        while t < t1 {
            ogeom_core::progress::checkpoint()?;
            let speed = heading.magnitude();
            if speed <= tol.confusion() {
                return Ok(None);
            }
            let mut taken = None;
            for _ in 0..40 {
                let next = (t + step / speed).min(t1);
                let Some((there, ahead)) = border.point_d1(self.surface, next, tol) else {
                    return Ok(None);
                };
                let turn = heading
                    .normalized(tol)
                    .ok()
                    .zip(ahead.normalized(tol).ok())
                    .map_or(0.0, |(d, e)| d.dot(e).clamp(-1.0, 1.0).acos());
                let sag = here.distance(there) * turn / 8.0;
                if sag <= chord || step <= tol.confusion() * 8.0 {
                    let scale = if sag > 0.0 {
                        (chord / sag).sqrt().clamp(0.5, 2.0)
                    } else {
                        2.0
                    };
                    taken = Some((
                        next,
                        there,
                        ahead,
                        (step * scale).clamp(tol.confusion(), ceiling),
                    ));
                    break;
                }
                step *= (chord / sag).sqrt().clamp(0.25, 0.9);
            }
            let Some((next, there, ahead, following)) = taken else {
                return Ok(None);
            };
            let Some(landed) = keep(&mut traced, next, Some(foot)) else {
                return Ok(None);
            };
            if traced.points.len() >= self.options.max_points {
                return Ok(None);
            }
            (t, here, heading, step, foot) = (next, there, ahead, following, landed);
        }
        Ok(Some(traced))
    }
}

/// Every branch of the intersection: seed, trace each, and keep the distinct
/// ones.
///
/// A branch crossing many sampling cells produces many seeds, and tracing from
/// any of them gives the same curve. So a seed already lying on something
/// traced is dropped rather than followed again, which is what makes the
/// *number* of branches returned meaningful, and it is the number a boolean
/// will act on.
///
/// Where the surfaces stand within the confusion distance `ε` of each other
/// across a band rather than along a curve, the section is that band to
/// tolerance. Two curves of it in one band, with nowhere between them that
/// the surfaces part by more than `ε`, may come back as one: each sample
/// lies on both surfaces within `ε`, and each point of the curve left out
/// joins one returned by a straight way that does too. A branch is placed
/// across the section only to within `ε` over the sine of the angle the
/// surfaces cross at, so where they cross shallowly its samples may zigzag
/// within that band and its length exceed the curve's.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable.
pub fn branches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Traced>> {
    let (found, lying) = seeded(a, b, options, tol)?;
    // A stretch of a border lying in the other surface is a branch already
    // traced, and a seed on it is not followed again.
    let mut out: Vec<Traced> = lying;
    let mut again: Vec<Traced> = Vec::new();
    let reach = options.chord.max(tol.confusion()) * 8.0;
    let follow = |seed: Contact, out: &mut Vec<Traced>, again: &mut Vec<Traced>| {
        // A seed that will not trace (a tangency) is reported by being
        // absent rather than by an error, since the other branches are still
        // real answers.
        let Ok(branch) = trace(a, b, seed, options, tol) else {
            return;
        };
        if branch.points.len() < 2 || is_fragment(&branch, options) {
            return;
        }
        for piece in cut_at_edge_touches(branch, a, b, tol) {
            // A branch whose middle lies on one already traced is that
            // branch again, reached from a seed its trace stopped short of,
            // or one that turned a corner where the section crosses itself
            // and ran on along another: held to find that crossing from.
            // A branch traced before that lies along this piece all its
            // length is, the other way round, this piece again: its walk
            // stopped at a touch this one ran over and was cut at, or at the
            // closure of a chart this one crossed. So a piece whose middle
            // lies only on such branches, and which runs on past them, is
            // the longer trace of their stretch and takes their place.
            let along: Vec<bool> = out
                .iter()
                .map(|other| other.points.iter().all(|&p| passes_near(&piece, p, reach)))
                .collect();
            let n = piece.points.len();
            let middle = if n > 2 {
                piece.points[n / 2]
            } else {
                piece.points[0].midpoint(piece.points[n - 1])
            };
            let on_middle: Vec<bool> = out
                .iter()
                .zip(&along)
                .filter(|(other, _)| passes_near(other, middle, reach))
                .map(|(_, &along)| along)
                .collect();
            let runs_past = || {
                piece.points.iter().any(|&p| {
                    !out.iter()
                        .zip(&along)
                        .any(|(other, &along)| along && passes_near(other, p, reach))
                })
            };
            if !on_middle.is_empty() && (on_middle.contains(&false) || !runs_past()) {
                again.push(piece);
                continue;
            }
            for (other, along) in core::mem::take(out).into_iter().zip(along) {
                if along {
                    again.push(other);
                } else {
                    out.push(other);
                }
            }
            out.push(piece);
        }
    };
    for seed in found {
        // Already on something followed.
        if out
            .iter()
            .any(|branch| passes_near(branch, seed.point, reach))
        {
            continue;
        }
        follow(seed, &mut out, &mut again);
    }
    // A touch of an edge is where two arcs of the section end, one either
    // side. A seed there lies on the arc already followed to it and is
    // skipped, and the arc on the other side is found only if some seed lies
    // on it: so each touch ending one arc alone is followed from, each way,
    // and the arc it runs back along is that arc again.
    let mut followed: Vec<Point> = Vec::new();
    loop {
        let alone: Vec<Contact> = lone_touches(&out, a, b, tol)
            .into_iter()
            .filter(|touch| {
                followed
                    .iter()
                    .all(|p| p.distance(touch.point) > tol.confusion())
            })
            .collect();
        if alone.is_empty() {
            break;
        }
        for touch in alone {
            followed.push(touch.point);
            follow(touch, &mut out, &mut again);
        }
    }
    let stitched = stitch_stalled(out, a, b, options, tol);
    Ok(split_at_touches(stitched, again, a, b, options, tol))
}

/// The touches of an edge that end one branch alone, as contacts to follow
/// the section from.
fn lone_touches(
    branches: &[Traced],
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Vec<Contact> {
    let pair = SurfacePair { a, b };
    let ends: Vec<(Point, [f64; 4])> = branches
        .iter()
        .filter(|branch| !branch.closed() && !branch.points.is_empty())
        .flat_map(|branch| {
            let last = branch.points.len() - 1;
            [0, last].map(|i| {
                let (on_a, on_b) = (branch.on_a[i], branch.on_b[i]);
                (branch.points[i], [on_a.0, on_a.1, on_b.0, on_b.1])
            })
        })
        .collect();
    let mut alone: Vec<Contact> = Vec::new();
    for &(point, state) in &ends {
        let touch = [(a, 0_usize), (b, 2)]
            .into_iter()
            .find_map(|(surface, first)| {
                let ((ua, ub), (va, vb)) = surface.domain();
                [
                    (first, (ua, ub), surface.is_periodic_u()),
                    (first + 1, (va, vb), surface.is_periodic_v()),
                ]
                .into_iter()
                .filter(|&(_, _, periodic)| !periodic)
                .flat_map(|(k, (lo, hi), _)| [(k, lo), (k, hi)])
                .filter(|&(k, bound)| (state[k] - bound).abs() <= tol.parametric())
                .find_map(|(k, bound)| touch_on_edge(&pair, &state, k, bound, tol))
            });
        let Some(at) = touch else {
            continue;
        };
        if a.point_at(at[0], at[1], tol)
            .is_ok_and(|p| p.distance(point) > tol.confusion())
        {
            continue;
        }
        let meeting = ends
            .iter()
            .filter(|(other, _)| other.distance(point) <= tol.confusion())
            .count();
        if meeting == 1 {
            alone.push(Contact {
                on_a: (at[0], at[1]),
                on_b: (at[2], at[3]),
                point,
            });
        }
    }
    alone
}

/// A branch cut where it runs over a touch of either surface's edge: the
/// pieces between the touches, each ending on the touch.
///
/// A walk ends at a touch of an edge where its step lands on the edge
/// ([`end_at_touch`]). A step can as well land beyond the touch, back inside
/// the domain, and the walk runs on over it. A section is then split at a
/// touch or not by where the steps fell, and the same stretch comes back
/// from two seeds as two different curves: one ending at the touch, one
/// running over it. So a branch is cut at every touch it runs over.
///
/// A touch run over is beside a sample where one parameter turns back
/// toward the edge it bounds. It is solved there by [`touch_on_edge`], and
/// kept where it lies within the steps either side of the sample.
fn cut_at_edge_touches(
    branch: Traced,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Vec<Traced> {
    let pair = SurfacePair { a, b };
    let closed = branch.closed();
    let n = branch.points.len();
    if n < 3 {
        return vec![branch];
    }
    let state = |i: usize| {
        [
            branch.on_a[i].0,
            branch.on_a[i].1,
            branch.on_b[i].0,
            branch.on_b[i].1,
        ]
    };
    // A loop's last sample is its first.
    let count = if closed { n - 1 } else { n };
    let period = |surface: &SurfaceGeometry, along_u: bool| {
        let ((ua, ub), (va, vb)) = surface.domain();
        if along_u {
            surface.is_periodic_u().then_some(ub - ua)
        } else {
            surface.is_periodic_v().then_some(vb - va)
        }
    };
    let periods = [
        period(a, true),
        period(a, false),
        period(b, true),
        period(b, false),
    ];
    let mut touches: Vec<(usize, [f64; 4], Point)> = Vec::new();
    for i in 0..count {
        let (before, after) = if closed {
            ((i + count - 1) % count, (i + 1) % count)
        } else if i == 0 || i + 1 == n {
            continue;
        } else {
            (i - 1, i + 1)
        };
        let here = state(i);
        let (back, ahead) = (
            branch.points[i].distance(branch.points[before]),
            branch.points[i].distance(branch.points[after]),
        );
        if back <= 0.0 || ahead <= 0.0 {
            continue;
        }
        let stride = back.max(ahead);
        // A periodic parameter is read on the turn of the middle sample's.
        let unwrap = |x: [f64; 4]| -> [f64; 4] {
            core::array::from_fn(|j| match periods[j] {
                Some(period) => x[j] - ((x[j] - here[j]) / period).round() * period,
                None => x[j],
            })
        };
        let around = [unwrap(state(before)), here, unwrap(state(after))];
        'bounds: for (surface, first) in [(a, 0_usize), (b, 2)] {
            let ((ua, ub), (va, vb)) = surface.domain();
            for (k, (lo, hi), periodic) in [
                (first, (ua, ub), surface.is_periodic_u()),
                (first + 1, (va, vb), surface.is_periodic_v()),
            ] {
                if periodic || !(hi - lo).is_finite() {
                    continue;
                }
                for bound in [lo, hi] {
                    let off = |x: &[f64; 4]| (x[k] - bound).abs();
                    if off(&here) >= off(&around[0]) || off(&here) > off(&around[2]) {
                        continue;
                    }
                    // A touch within the steps either side lies on an edge
                    // passing within them.
                    let mut on_edge = here;
                    on_edge[k] = bound;
                    let edge = surface.point_at(on_edge[first], on_edge[first + 1], tol);
                    if edge.is_ok_and(|p| p.distance(branch.points[i]) > stride) {
                        continue;
                    }
                    let start = turning(around, (back, ahead), k);
                    let Some(at) = touch_on_edge(&pair, &start, k, bound, tol) else {
                        continue;
                    };
                    let Ok(point) = a.point_at(at[0], at[1], tol) else {
                        continue;
                    };
                    if point.distance(branch.points[i]) <= stride {
                        touches.push((i, at, point));
                        break 'bounds;
                    }
                }
            }
        }
    }
    if touches.is_empty() {
        return vec![branch];
    }
    // The samples with each touch in its place, and where the touches are.
    let mut samples: Vec<(Point, [f64; 4])> = Vec::with_capacity(n + touches.len());
    let mut cuts: Vec<usize> = Vec::new();
    let mut next = touches.iter().peekable();
    for i in 0..count {
        let Some(&(_, at, point)) = next.next_if(|t| t.0 == i) else {
            samples.push((branch.points[i], state(i)));
            continue;
        };
        let (before, after) = (
            branch.points[if i == 0 { count - 1 } else { i - 1 }],
            branch.points[(i + 1) % n],
        );
        if point.distance(branch.points[i]) <= tol.confusion() {
            cuts.push(samples.len());
            samples.push((point, at));
        } else if (point - branch.points[i]).dot(after - before) > 0.0 {
            samples.push((branch.points[i], state(i)));
            cuts.push(samples.len());
            samples.push((point, at));
        } else {
            cuts.push(samples.len());
            samples.push((point, at));
            samples.push((branch.points[i], state(i)));
        }
    }
    let piece = |run: &[(Point, [f64; 4])], stopped: Stopped| Traced {
        points: run.iter().map(|s| s.0).collect(),
        on_a: run.iter().map(|s| (s.1[0], s.1[1])).collect(),
        on_b: run.iter().map(|s| (s.1[2], s.1[3])).collect(),
        stopped,
    };
    let mut pieces = Vec::new();
    if closed {
        // Round the loop from its first touch back to it.
        let first = cuts[0];
        let mut round: Vec<(Point, [f64; 4])> = samples[first..].to_vec();
        round.extend_from_slice(&samples[..=first]);
        let mut from = 0;
        for cut in cuts
            .iter()
            .skip(1)
            .map(|c| c - first)
            .chain([round.len() - 1])
        {
            if cut > from {
                pieces.push(piece(&round[from..=cut], Stopped::LeftTheDomain));
            }
            from = cut;
        }
    } else {
        let mut from = 0;
        for cut in cuts.iter().copied().chain([samples.len() - 1]) {
            if cut > from {
                pieces.push(piece(&samples[from..=cut], Stopped::LeftTheDomain));
            }
            from = cut;
        }
        // The ends the trace stopped at are its own: a branch truncated or
        // stalled at either end is reported so in its first and last piece.
        if let Some(first) = pieces.first_mut() {
            first.stopped = branch.stopped;
        }
        if let Some(last) = pieces.last_mut() {
            last.stopped = branch.stopped;
        }
    }
    pieces
}

/// Where unknown `k` turns along three samples of a branch, to second
/// order: the states of the quadratics through them, as functions of the
/// length along their chords, `back` behind the middle one and `ahead`
/// past it, at the quadratic's vertex in `k`, held between the outer two.
fn turning(around: [[f64; 4]; 3], (back, ahead): (f64, f64), k: usize) -> [f64; 4] {
    let [prior, here, next] = around;
    let (rise, fall) = ((here[k] - prior[k]) / back, (next[k] - here[k]) / ahead);
    let t = if fall == rise {
        0.0
    } else {
        ((rise * (ahead + back) / (rise - fall) - back) / 2.0).clamp(-back, ahead)
    };
    let weights = [
        t * (t - ahead) / (back * (back + ahead)),
        (t + back) * (ahead - t) / (back * ahead),
        (t + back) * t / (ahead * (ahead + back)),
    ];
    core::array::from_fn(|j| prior[j] * weights[0] + here[j] * weights[1] + next[j] * weights[2])
}

/// Below this sine the surfaces count as tangent at a point: the
/// branch-point certificate a stall end must carry to participate in
/// stitching.
pub(crate) const BRANCH_POINT_SINE: f64 = 0.05;

/// The sine of the normal angle at a contact: the transversality measure.
fn crossing_sine(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> f64 {
    let Ok(na) = a.normal_at(on_a.0, on_a.1, tol) else {
        return 0.0;
    };
    let Ok(nb) = b.normal_at(on_b.0, on_b.1, tol) else {
        return 0.0;
    };
    na.vector().cross(nb.vector()).magnitude()
}

/// Whether a stalled trace is a fragment rather than a curve.
///
/// Coincident or near-coincident surfaces defeat the tangency check at a seed:
/// rounding in the corrected parameters leaves the two normals a whisker apart,
/// the walk takes a couple of steps, and then stalls where the arithmetic gives
/// out. What comes back lies on both surfaces perfectly and describes nothing:
/// identical spheres yield such fragments, each a few points long.
///
/// A stalled branch shorter than a handful of chords carries no information the
/// seed did not, so it is noise from a degenerate configuration and dropped. A
/// *real* stalled branch (one that ran into a genuine tangency) has length
/// behind it and is kept, because a truncated real answer is still an answer.
///
/// The marcher is deliberately not a coincidence detector: for the pairs with
/// closed forms, [`surface_surface`](crate::surface_surface) answers
/// [`Same`](crate::Meeting::Same), and that check belongs before this one.
fn is_fragment(branch: &Traced, options: Marching) -> bool {
    if branch.stopped != Stopped::Stalled {
        return false;
    }
    let length: f64 = branch
        .points
        .windows(2)
        .map(|pair| pair[0].distance(pair[1]))
        .sum();
    length < options.chord * 10.0
}

/// Whether a traced branch passes within a distance of a point.
///
/// Measured against the polyline's *segments*, not its vertices. The vertices
/// are a marching step apart (far more than the chord tolerance), so a seed
/// sitting neatly between two of them looks distant from both, and comparing to
/// vertices alone would report one circle many times.
fn passes_near(branch: &Traced, p: Point, reach: f64) -> bool {
    branch
        .points
        .windows(2)
        .any(|pair| distance_to_segment(p, pair[0], pair[1]) <= reach)
}

/// Distance from a point to a segment.
pub(crate) fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.square_magnitude();
    if length <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

/// Follow the intersection from a starting point, in both directions.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the two surfaces
/// are tangent at the seed, where there is no single direction to follow.
pub fn trace(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    options.validate()?;
    if tangent_at(a, b, from, tol).is_none() {
        ogeom_bail!(
            NotDone,
            "the surfaces are tangent here, so the intersection has no single \
             direction to follow; that is a branch point and needs the seed \
             moved off it"
        );
    }

    // Forwards first. If it closes, that is the whole branch and there is
    // nothing behind. A backward half that closes is the whole branch too,
    // whatever stopped the forward one, turned to run the way that set off.
    let ahead = walk(a, b, from, 1.0, options, tol)?;
    if ahead.stopped == Stopped::Closed {
        return Ok(ahead);
    }
    let behind = walk(a, b, from, -1.0, options, tol)?;
    if behind.stopped == Stopped::Closed {
        return Ok(reversed(behind));
    }
    let last_step = |walked: &[Point]| -> f64 {
        walked
            .windows(2)
            .last()
            .map_or(0.0, |w| w[0].distance(w[1]))
    };
    let steps = last_step(&ahead.points).max(last_step(&behind.points));

    // Join them, the backward half reversed and its shared first point dropped.
    let mut points = behind.points;
    let mut on_a = behind.on_a;
    let mut on_b = behind.on_b;
    points.reverse();
    on_a.reverse();
    on_b.reverse();
    points.pop();
    on_a.pop();
    on_b.pop();
    points.extend(ahead.points);
    on_a.extend(ahead.on_a);
    on_b.extend(ahead.on_b);

    // The worse of the two reasons: a branch truncated at either end is
    // truncated.
    let mut stopped = if ahead.stopped == Stopped::RanOut || behind.stopped == Stopped::RanOut {
        Stopped::RanOut
    } else if ahead.stopped == Stopped::Stalled || behind.stopped == Stopped::Stalled {
        Stopped::Stalled
    } else {
        Stopped::LeftTheDomain
    };
    // A loop cut at a seam. A closed patch is clamped, not periodic: a
    // section that runs round it (a rim's circle on a converted drum)
    // is walked from the seed to the seam one way and to the seam the
    // other, each walk stopping a fraction of a step short of it, and the
    // two ends meet where the surface closes on itself. That is the whole
    // loop, and it is closed: left open, the arrangement downstream would
    // hold a circle with two ends at one point and find no face piece to keep.
    // The ends are within a couple of the walks' own last steps of each
    // other, and the loop is closed exactly on its first point.
    if stopped == Stopped::LeftTheDomain && points.len() > 3 {
        let gap = points[0].distance(points[points.len() - 1]);
        if gap <= (steps * 2.0).max(tol.confusion() * 10.0) {
            // Both walks may have landed on the seam itself, one point.
            if gap <= tol.confusion() {
                points.pop();
                on_a.pop();
                on_b.pop();
            }
            points.push(points[0]);
            on_a.push(on_a[0]);
            on_b.push(on_b[0]);
            stopped = Stopped::Closed;
        }
    }
    Ok(Traced {
        points,
        on_a,
        on_b,
        stopped,
    })
}

/// A traced curve run the other way.
pub(crate) fn reversed(mut traced: Traced) -> Traced {
    traced.points.reverse();
    traced.on_a.reverse();
    traced.on_b.reverse();
    traced
}

/// Two surfaces, as a condition for the walker: four unknowns and three
/// equations saying the two points coincide.
///
/// The intersector's own walk goes through [`crate::walk`] like everything
/// else, and what is *not* generic lives here: the domain clamps a surface
/// pair needs, and the tangent, which the intersector computes from the two
/// normals rather than from the null space so that it can refuse a crossing
/// too shallow to be more than the correction's own noise.
struct SurfacePair<'s> {
    a: &'s SurfaceGeometry,
    b: &'s SurfaceGeometry,
}

impl crate::walk::Condition for SurfacePair<'_> {
    fn unknowns(&self) -> usize {
        4
    }

    fn position(&self, x: &[f64], tol: Tolerances) -> Option<Point> {
        self.a.point_at(x[0], x[1], tol).ok()
    }

    fn position_gradient(&self, x: &[f64], tol: Tolerances) -> Option<Vec<Vector>> {
        let (au, av) = self.a.d1_at(x[0], x[1], tol).ok()?;
        // The point is taken from the first surface, so it does not move with
        // the second's parameters at all.
        Some(vec![au, av, Vector::ZERO, Vector::ZERO])
    }

    fn system(&self, x: &[f64], tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
        Some(self.system_at(x, tol)?.0)
    }

    fn system_at(
        &self,
        x: &[f64],
        tol: Tolerances,
    ) -> Option<((Vec<f64>, Vec<Vec<f64>>), Point, Vec<Vector>)> {
        let (pa, au, av) = self.a.point_d1_at(x[0], x[1], tol).ok()?;
        let (pb, bu, bv) = self.b.point_d1_at(x[2], x[3], tol).ok()?;
        let gap = pa - pb;
        Some((
            (
                vec![gap.x, gap.y, gap.z],
                vec![
                    vec![au.x, av.x, -bu.x, -bv.x],
                    vec![au.y, av.y, -bu.y, -bv.y],
                    vec![au.z, av.z, -bu.z, -bv.z],
                ],
            ),
            pa,
            vec![au, av, Vector::ZERO, Vector::ZERO],
        ))
    }

    fn clamp(&self, x: &mut [f64]) {
        let (ua, va) = clamp(self.a, x[0], x[1]);
        let (ub, vb) = clamp(self.b, x[2], x[3]);
        x[0] = ua;
        x[1] = va;
        x[2] = ub;
        x[3] = vb;
    }

    fn outside(&self, x: &[f64], tol: Tolerances) -> bool {
        outside(self.a, (x[0], x[1]), tol) || outside(self.b, (x[2], x[3]), tol)
    }

    fn near_edge(&self, x: &[f64]) -> bool {
        near_edge(self.a, (x[0], x[1])) || near_edge(self.b, (x[2], x[3]))
    }

    fn extent(&self) -> f64 {
        span(self.a).max(span(self.b))
    }

    fn tangent_is_oriented(&self) -> bool {
        // The cross product of the two normals, whose sign is the surfaces'
        // own and whose flip at a tangency is what stops the march.
        true
    }

    fn tangent(&self, x: &[f64], tol: Tolerances) -> Option<Vector> {
        tangent_at(
            self.a,
            self.b,
            Contact {
                on_a: (x[0], x[1]),
                on_b: (x[2], x[3]),
                point: Point::ORIGIN,
            },
            tol,
        )
    }

    fn tangent_from(
        &self,
        x: &[f64],
        jacobian: &[Vec<f64>],
        gradient: &[Vector],
        tol: Tolerances,
    ) -> Option<Vector> {
        // The system's columns are the two surfaces' first derivatives, the
        // second's negated: the normals the tangent is made of.
        if jacobian.len() != 3 || jacobian.iter().any(|row| row.len() != 4) {
            return self.tangent(x, tol);
        }
        let _ = gradient;
        let column = |c: usize| Vector::new(jacobian[0][c], jacobian[1][c], jacobian[2][c]);
        let na = unit_normal(column(0), column(1), tol)?;
        let nb = unit_normal(-column(2), -column(3), tol)?;
        tangent_of_normals(na, nb, tol)
    }
}

/// Walk one way from a seed.
fn walk(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    sense: f64,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    let pair = SurfacePair { a, b };
    let start = [from.on_a.0, from.on_a.1, from.on_b.0, from.on_b.1];
    let (mut walked, past) = crate::walk::walk_one_way_past(&pair, &start, sense, options, tol)?;
    if walked.stopped == Stopped::LeftTheDomain
        && !past.is_some_and(|past| end_at_touch(&pair, &mut walked, &past, tol))
    {
        land_on_edge(&pair, &mut walked, tol);
    }
    Ok(Traced {
        on_a: walked.states.iter().map(|x| (x[0], x[1])).collect(),
        on_b: walked.states.iter().map(|x| (x[2], x[3])).collect(),
        points: walked.points,
        stopped: walked.stopped,
    })
}

/// End a walk that left a surface's domain on the edge it left through.
///
/// The walk stops a fraction of a step short of the edge, and a section
/// ending there misses the edge by that much: across a patch boundary the
/// section on the next patch starts on the edge, and the two ends stand
/// apart by more than either curve's tolerance. The last point is carried
/// onto the nearest non-periodic bound of either surface, the pair still
/// meeting there, and kept if it lies ahead within a couple of steps.
fn land_on_edge(pair: &SurfacePair<'_>, walked: &mut crate::walk::Walked, tol: Tolerances) {
    use crate::walk::Condition as _;
    let n = walked.states.len();
    if n < 2 {
        return;
    }
    let (last, before) = (&walked.states[n - 1], &walked.states[n - 2]);
    let step = walked.points[n - 1].distance(walked.points[n - 2]);
    let heading = walked.points[n - 1] - walked.points[n - 2];
    let mut bounds: Vec<(f64, usize, f64)> = Vec::new();
    for (surface, first) in [(pair.a, 0_usize), (pair.b, 2)] {
        let ((ua, ub), (va, vb)) = surface.domain();
        for (k, (lo, hi), periodic) in [
            (first, (ua, ub), surface.is_periodic_u()),
            (first + 1, (va, vb), surface.is_periodic_v()),
        ] {
            let span = hi - lo;
            if periodic || !span.is_finite() || span <= 0.0 {
                continue;
            }
            // Only a bound the walk was heading for, and one that is an
            // edge: a sphere's pole bounds its chart at a point.
            let moving = last[k] - before[k];
            let across = if k == first {
                surface.domain().1
            } else {
                surface.domain().0
            };
            let collapses = |bound: f64| {
                let at = |t: f64| {
                    if k == first {
                        surface.point_at(bound, t, tol)
                    } else {
                        surface.point_at(t, bound, tol)
                    }
                };
                match (at(across.0), at(f64::midpoint(across.0, across.1))) {
                    (Ok(p), Ok(q)) => p.distance(q) <= tol.confusion(),
                    _ => true,
                }
            };
            for (bound, toward) in [(lo, moving < 0.0), (hi, moving > 0.0)] {
                if (toward || (last[k] - bound).abs() <= span * 1e-4) && !collapses(bound) {
                    bounds.push(((last[k] - bound).abs() / span, k, bound));
                }
            }
        }
    }
    bounds.sort_by(|x, y| x.0.total_cmp(&y.0));
    for (_, k, bound) in bounds {
        let system = |x: &[f64; 4]| {
            let mut residual = [f64::INFINITY; 4];
            let mut jacobian = [[0.0; 4]; 4];
            let Some(((rows, matrix), _, _)) = pair.system_at(x, tol) else {
                return (residual, jacobian);
            };
            residual[..3].copy_from_slice(&rows);
            for (to, row) in jacobian.iter_mut().zip(&matrix) {
                to.copy_from_slice(row);
            }
            residual[3] = x[k] - bound;
            jacobian[3][k] = 1.0;
            (residual, jacobian)
        };
        let criteria = solve::Criteria {
            residual: tol.confusion() * 0.01,
            step: tol.parametric(),
            max_iterations: 40,
        };
        let Ok(start) = <[f64; 4]>::try_from(last.as_slice()) else {
            return;
        };
        let Ok((at, norm, _, _)) = solve::newton_system_fixed(system, start, criteria) else {
            continue;
        };
        if norm > tol.confusion() {
            continue;
        }
        let within = [(pair.a, 0_usize), (pair.b, 2)]
            .iter()
            .all(|(surface, first)| {
                let ((ua, ub), (va, vb)) = surface.domain();
                [
                    (*first, ua, ub, surface.is_periodic_u()),
                    (first + 1, va, vb, surface.is_periodic_v()),
                ]
                .iter()
                .all(|&(i, lo, hi, periodic)| {
                    periodic || (at[i] >= lo - tol.parametric() && at[i] <= hi + tol.parametric())
                })
            });
        if !within {
            continue;
        }
        let Ok(point) = pair.a.point_at(at[0], at[1], tol) else {
            continue;
        };
        let ahead = point - walked.points[n - 1];
        if ahead.magnitude() > (step * 2.0).max(tol.confusion())
            || ahead.dot(heading) < -tol.confusion() * step
        {
            continue;
        }
        if ahead.magnitude() <= tol.confusion() {
            walked.states[n - 1] = at.to_vec();
            walked.points[n - 1] = point;
        } else {
            walked.states.push(at.to_vec());
            walked.points.push(point);
        }
        return;
    }
}

/// End a walk whose last step landed on a surface's edge it only touches,
/// at the touch: whether it did.
///
/// A section touching an edge, a drum's rim lying in a plane tangent to the
/// other surface, runs along the edge to second order: within the
/// confusion distance of the other surface along a stretch of the edge
/// about `√(2·ε·ρ)` long, `ρ` the radius their parting curves at. A step
/// landing there is held on the edge with its residual within the
/// correction's, and reads as having left. Solving the section onto the
/// edge from there is singular, its root a double one: the solve stops
/// where the residual, the square of the distance left, is small, short of
/// the touch by the square root of the residual. The touch is solved as one
/// instead, from the landing, by [`touch_on_edge`], and kept where it lies
/// ahead of the walk's last point.
fn end_at_touch(
    pair: &SurfacePair<'_>,
    walked: &mut crate::walk::Walked,
    (landed, beyond): &crate::walk::Past,
    tol: Tolerances,
) -> bool {
    let n = walked.points.len();
    if n == 0 || landed.len() != 4 {
        return false;
    }
    let heading = *beyond - walked.points[n - 1];
    for (surface, first) in [(pair.a, 0_usize), (pair.b, 2)] {
        let ((ua, ub), (va, vb)) = surface.domain();
        for (k, (lo, hi), periodic) in [
            (first, (ua, ub), surface.is_periodic_u()),
            (first + 1, (va, vb), surface.is_periodic_v()),
        ] {
            if periodic {
                continue;
            }
            for bound in [lo, hi] {
                if (landed[k] - bound).abs() > tol.parametric() {
                    continue;
                }
                let Some(at) = touch_on_edge(pair, landed, k, bound, tol) else {
                    continue;
                };
                let Ok(point) = pair.a.point_at(at[0], at[1], tol) else {
                    continue;
                };
                let ahead = point - walked.points[n - 1];
                if ahead.dot(heading) <= 0.0 {
                    continue;
                }
                if ahead.magnitude() <= tol.confusion() {
                    walked.states[n - 1] = at.to_vec();
                    walked.points[n - 1] = point;
                } else {
                    walked.states.push(at.to_vec());
                    walked.points.push(point);
                }
                return true;
            }
        }
    }
    false
}

/// Where the section touches a surface's edge without crossing it, solved
/// from a state beside it: the state there, or `None` if the edge has no
/// such touch there.
///
/// `k` is the unknown the edge holds at `bound`. Along the edge the
/// bordered surface's point runs as a curve of its free parameter, and the
/// section touches the edge where that curve touches the other surface, as
/// [`crate::curve_surface::touch_solved`] solves it. The touch is kept only
/// within the confusion distance, and only where the section leaves the
/// edge into the domain: a section outside it, touching the edge from
/// there, meets the domain at the touch alone, and a walk along the edge
/// beside it is the correction held on the edge, not the section.
fn touch_on_edge(
    pair: &SurfacePair<'_>,
    from: &[f64],
    k: usize,
    bound: f64,
    tol: Tolerances,
) -> Option<[f64; 4]> {
    let (edged, other, first, rest) = if k < 2 {
        (pair.a, pair.b, 0_usize, 2_usize)
    } else {
        (pair.b, pair.a, 2, 0)
    };
    let along_u = k == first + 1;
    let free = if along_u { first } else { first + 1 };
    let place = |s: f64| if along_u { (s, bound) } else { (bound, s) };
    let along = |s: f64| {
        let (u, v) = place(s);
        let (p, su, sv) = edged.point_d1_at(u, v, tol).ok()?;
        let (suu, _, svv) = edged.d2_at(u, v, tol).ok()?;
        Some(if along_u { (p, su, suu) } else { (p, sv, svv) })
    };
    let start = [from[free], from[rest], from[rest + 1]];
    let found = crate::curve_surface::touch_solved(along, other, start, tol)?;
    let (u, v) = place(found[0]);
    let within = |surface: &SurfaceGeometry, (u, v): (f64, f64)| {
        let ((ua, ub), (va, vb)) = surface.domain();
        [
            (u, ua, ub, surface.is_periodic_u()),
            (v, va, vb, surface.is_periodic_v()),
        ]
        .iter()
        .all(|&(x, lo, hi, periodic)| {
            periodic || (x >= lo - tol.parametric() && x <= hi + tol.parametric())
        })
    };
    if !within(edged, (u, v)) || !within(other, (found[1], found[2])) {
        return None;
    }
    let (p, su, sv) = edged.point_d1_at(u, v, tol).ok()?;
    let (suu, _, svv) = edged.d2_at(u, v, tol).ok()?;
    let (q, ou, ov) = other.point_d1_at(found[1], found[2], tol).ok()?;
    if p.distance(q) > tol.confusion() {
        return None;
    }
    // Along the edge the gap to the other surface grows as `g·s²/2`, `g`
    // the edge's normal curvature less the other surface's along it, and
    // off the edge as `a·t`, `a` the edge's cross derivative along the
    // normal. So the section stands `−g·s² / (2·a)` across the edge: inside
    // where that has the sign of the way into the domain. Where `g` or `a`
    // vanishes the contact is of higher order, and its side is not read.
    let (c, c2, across) = if along_u {
        (su, suu, sv)
    } else {
        (sv, svv, su)
    };
    let (ouu, ouv, ovv) = other.d2_at(found[1], found[2], tol).ok()?;
    let normal = ou.cross(ov);
    let metric = [[ou.dot(ou), ou.dot(ov)], [ou.dot(ov), ov.dot(ov)]];
    let [cu, cv] = solve::linear_fixed(metric, [c.dot(ou), c.dot(ov)])?;
    let bend = (c2 - (ouu * (cu * cu) + ouv * (2.0 * cu * cv) + ovv * (cv * cv))).dot(normal);
    let ((ua, ub), (va, vb)) = edged.domain();
    let (lo, hi) = if along_u { (va, vb) } else { (ua, ub) };
    let inward = if (bound - hi).abs() <= (bound - lo).abs() {
        -1.0
    } else {
        1.0
    };
    let side = -bend * across.dot(normal) * inward;
    if side <= 0.0 || side.is_nan() {
        return None;
    }
    let mut at = [0.0; 4];
    at[first] = u;
    at[first + 1] = v;
    at[rest] = found[1];
    at[rest + 1] = found[2];
    Some(at)
}

/// The sine of the shallowest crossing angle the marcher will follow.
///
/// One microradian, and the number is set by the *correction*, not by taste.
/// `correct` accepts a residual up to the confusion tolerance, so the two
/// parameter points of a contact can disagree by that much in space, and on
/// coincident or near-coincident surfaces, that disagreement shows up as a
/// spurious angle between the two computed normals of about the residual over
/// the local feature size. A gate below that floor reads the correction's own
/// noise as a direction and marches along it: identical spheres come back as
/// confident little curves that exist nowhere but in rounding.
///
/// So below this angle the marcher cannot tell an ultra-shallow crossing from
/// coincidence, and refuses both rather than guessing. A genuine crossing
/// shallower than a microradian is also one the Newton correction cannot
/// reliably follow (its travel constraint becomes numerically dependent on
/// the surface-gap rows at exactly the same rate), so the gate refuses what
/// could not have been followed anyway.
const SHALLOWEST: f64 = 1e-6;

/// The direction the intersection runs at a contact.
///
/// The cross product of the two normals: the one direction lying in both
/// tangent planes. `None` where the normals are parallel to within
/// [`SHALLOWEST`]: the surfaces are tangent or coincident there, and the
/// intersection has no direction the marcher can trust.
fn tangent_at(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Contact,
    tol: Tolerances,
) -> Option<Vector> {
    let na = normal_at(a, at.on_a, tol)?;
    let nb = normal_at(b, at.on_b, tol)?;
    tangent_of_normals(na, nb, tol)
}

/// The intersection tangent from the two unit normals, as [`tangent_at`]
/// decides it.
fn tangent_of_normals(na: Vector, nb: Vector, tol: Tolerances) -> Option<Vector> {
    let cross = na.cross(nb);
    let length = cross.magnitude();
    // The decision is made through intervals rather than a bare compare:
    // each normal component carries the correction's stated residual as an
    // uncertainty, the squared cross magnitude is computed as an enclosure,
    // and only a crossing *certainly* above the floor is followed. A sine
    // inside the enclosure's undecided band is exactly the case the floor
    // exists for (the correction's own noise masquerading as an angle),
    // and it is refused with a certificate instead of a guess.
    let floor = tol.angular().max(SHALLOWEST);
    let widen = |value: f64| ogeom_math::Interval::about(value, tol.confusion());
    let (ax, ay, az) = (widen(na.x), widen(na.y), widen(na.z));
    let (bx, by, bz) = (widen(nb.x), widen(nb.y), widen(nb.z));
    let cx = ay.mul(&bz).sub(&az.mul(&by));
    let cy = az.mul(&bx).sub(&ax.mul(&bz));
    let cz = ax.mul(&by).sub(&ay.mul(&bx));
    let magnitude2 = cx.square().add(&cy.square()).add(&cz.square());
    let above = magnitude2.sub(&ogeom_math::Interval::point(floor * floor));
    if above.certain_sign() != Some(ogeom_core::Sign::Positive) || length <= f64::MIN_POSITIVE {
        return None;
    }
    Some(cross * (1.0 / length))
}

/// A surface's unit normal at a parameter.
fn normal_at(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> Option<Vector> {
    let (du, dv) = surface.d1_at(at.0, at.1, tol).ok()?;
    unit_normal(du, dv, tol)
}

/// The unit normal of two tangents, `None` where they are parallel.
fn unit_normal(du: Vector, dv: Vector, tol: Tolerances) -> Option<Vector> {
    let cross = du.cross(dv);
    let length = cross.magnitude();
    if length <= tol.confusion() {
        return None;
    }
    Some(cross * (1.0 / length))
}

/// A state [`correct`] measured, the system's Jacobian and residual there,
/// and the first surface's two derivatives.
type Measured = ([f64; 4], [[f64; 4]; 4], [f64; 4], [Vector; 2]);

/// Bring a parameter guess onto both surfaces.
///
/// Three equations say the two surface points coincide; the fourth says how far
/// along the direction of travel to land. Without that fourth the system is
/// underdetermined (its solution set *is* the curve), and Newton would wander
/// along it instead of converging to a point.
///
/// `constraint` is `(anchor, direction, distance)`. Without one, the guess
/// itself is used as the anchor and the direction is the intersection tangent,
/// which is what seeding wants: land anywhere on the curve near here. The
/// point it lands on stands within `within` of the section as far as the
/// solve can tell, as [`crate::walk::settling_residual`] sets out: where the
/// surfaces cross at a shallow angle the residual is taken further down.
fn correct(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    start: [f64; 4],
    guess: Point,
    constraint: Option<(Point, Vector, f64)>,
    within: f64,
    tol: Tolerances,
) -> Option<Contact> {
    let (anchor, along, reach) = match constraint {
        Some(given) => given,
        None => {
            // No direction to travel: hold the guess still along whichever way
            // the curve runs, so the solve slides onto the curve rather than
            // along it.
            let at = Contact {
                on_a: (start[0], start[1]),
                on_b: (start[2], start[3]),
                point: guess,
            };
            (guess, tangent_at(a, b, at, tol).unwrap_or(Vector::X), 0.0)
        }
    };

    // The state and Jacobian last measured, and the first surface's
    // derivatives there: what the point moves by with each unknown.
    let last: core::cell::Cell<Option<Measured>> = core::cell::Cell::new(None);
    let system = |x: &[f64; 4]| {
        last.set(None);
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        // Where either surface cannot be evaluated the residual is
        // infinite, so the damped step backs off rather than reading a
        // made-up point as a root.
        let (Ok((pa, au, av)), Ok((pb, bu, bv))) =
            (a.point_d1_at(ua, va, tol), b.point_d1_at(ub, vb, tol))
        else {
            return ([f64::INFINITY; 4], [[0.0; 4]; 4]);
        };

        let gap = pa - pb;
        let residual = [gap.x, gap.y, gap.z, (pa - anchor).dot(along) - reach];
        let jacobian = [
            [au.x, av.x, -bu.x, -bv.x],
            [au.y, av.y, -bu.y, -bv.y],
            [au.z, av.z, -bu.z, -bv.z],
            [au.dot(along), av.dot(along), 0.0, 0.0],
        ];
        last.set(Some((*x, jacobian, residual, [au, av])));
        (residual, jacobian)
    };

    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let mut found = solve::newton_system_fixed(system, start, criteria).ok()?;
    if found.1 > tol.confusion() {
        return None;
    }
    if last.get().is_none_or(|(x, _, _, _)| x != found.0) {
        system(&found.0);
    }
    if let Some((_, jacobian, residual, [au, av])) = last.get() {
        let gradient = [au, av, Vector::ZERO, Vector::ZERO];
        if let Some(settled) = crate::walk::settle(
            system,
            (found.0, &residual),
            (&jacobian, &gradient),
            within,
            criteria,
        ) {
            found.0 = settled;
        }
    }
    let (ua, va) = clamp(a, found.0[0], found.0[1]);
    let (ub, vb) = clamp(b, found.0[2], found.0[3]);
    Some(Contact {
        on_a: (ua, va),
        on_b: (ub, vb),
        point: a.point_at(ua, va, tol).ok()?,
    })
}

/// Hold a parameter inside a surface's domain.
///
/// A periodic direction wraps instead, so a curve crossing a cylinder's seam
/// keeps going rather than stopping at a boundary that is not one.
pub(crate) fn clamp(surface: &SurfaceGeometry, u: f64, v: f64) -> (f64, f64) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let fold = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if !periodic {
            return x.clamp(lo, hi);
        }
        let span = hi - lo;
        if span <= 0.0 {
            return x;
        }
        lo + (x - lo).rem_euclid(span)
    };
    (
        fold(u, ua, ub, surface.is_periodic_u()),
        fold(v, va, vb, surface.is_periodic_v()),
    )
}

/// Whether a parameter sits close enough to a non-periodic edge that a stalled
/// walk there means the edge rather than a singularity.
///
/// The band is a fraction of the domain's own span: a walk stalls within a
/// step of the boundary, and the step is far larger than the strict band
/// [`outside`] uses to decide a point has actually crossed.
fn near_edge(surface: &SurfaceGeometry, at: (f64, f64)) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let close = |x: f64, lo: f64, hi: f64, periodic: bool| {
        !periodic && {
            let band = (hi - lo).abs() * 1e-4;
            x <= lo + band || x >= hi - band
        }
    };
    close(at.0, ua, ub, surface.is_periodic_u()) || close(at.1, va, vb, surface.is_periodic_v())
}

/// Whether a parameter has left a surface's domain, in a direction that has one.
fn outside(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let past = |x: f64, lo: f64, hi: f64, periodic: bool| {
        !periodic && (x <= lo + tol.parametric() || x >= hi - tol.parametric())
    };
    past(at.0, ua, ub, surface.is_periodic_u()) || past(at.1, va, vb, surface.is_periodic_v())
}

/// A surface's rough size, for spacing seeds.
pub(crate) fn span(surface: &SurfaceGeometry) -> f64 {
    let ((ua, ub), (va, vb)) = surface.domain();
    let tol = Tolerances::millimetres();
    let corners = [(ua, va), (ub, va), (ua, vb), (ub, vb)];
    let mut low = Point::new(f64::MAX, f64::MAX, f64::MAX);
    let mut high = Point::new(f64::MIN, f64::MIN, f64::MIN);
    for (u, v) in corners {
        if let Ok(p) = surface.point_at(u, v, tol) {
            low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
            high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
        }
    }
    let size = (high - low).magnitude();
    if size.is_finite() && size > 0.0 {
        size
    } else {
        1.0
    }
}

/// One sampled triangle of a surface, with the parameters it came from.
#[derive(Debug, Clone)]
pub(crate) struct Cell {
    pub(crate) corners: [Point; 3],
    pub(crate) at: (f64, f64),
    pub(crate) low: Point,
    pub(crate) high: Point,
    /// How far the surface bows from the flat cell: its middle's distance
    /// from the middle of the diagonal the two triangles share.
    pub(crate) sag: f64,
    /// The parameters of the three corners, in `corners` order.
    pub(crate) params: [(f64, f64); 3],
}

/// Sample a surface into triangles.
pub(crate) fn sample(surface: &SurfaceGeometry, grid: usize, tol: Tolerances) -> Vec<Cell> {
    sample_by(surface, (grid, grid), tol)
}

/// Sample a surface into triangles, `counts` cells along `u` and along `v`.
pub(crate) fn sample_by(
    surface: &SurfaceGeometry,
    counts: (usize, usize),
    tol: Tolerances,
) -> Vec<Cell> {
    let ((ua, ub), (va, vb)) = surface.domain();
    // An unbounded domain would put the samples a billion units apart and find
    // nothing. Clamped to something a real model lives inside.
    let limit = 1.0e6;
    let (ua, ub) = (ua.max(-limit), ub.min(limit));
    let (va, vb) = (va.max(-limit), vb.min(limit));

    let mut out = Vec::new();
    #[allow(clippy::cast_precision_loss)]
    let (nu, nv) = (counts.0 as f64, counts.1 as f64);
    let at = |s: f64, t: f64| {
        let (u, v) = (ua + (ub - ua) * s, va + (vb - va) * t);
        surface.point_at(u, v, tol).ok().map(|p| ((u, v), p))
    };
    // Each grid corner is evaluated once and shared by the cells around it,
    // a row of corners at a time.
    #[allow(clippy::cast_precision_loss)]
    let row = |i: usize| -> Vec<_> {
        (0..=counts.1)
            .map(|j| at(i as f64 / nu, j as f64 / nv))
            .collect()
    };
    let mut below = row(0);
    for i in 0..counts.0 {
        let above = row(i + 1);
        for j in 0..counts.1 {
            #[allow(clippy::cast_precision_loss)]
            let (s0, s1) = (i as f64 / nu, (i + 1) as f64 / nu);
            #[allow(clippy::cast_precision_loss)]
            let (t0, t1) = (j as f64 / nv, (j + 1) as f64 / nv);
            let (Some((p00, a00)), Some((p10, a10)), Some((p01, a01)), Some((p11, a11))) =
                (below[j], above[j], below[j + 1], above[j + 1])
            else {
                continue;
            };
            let sag = at(f64::midpoint(s0, s1), f64::midpoint(t0, t1))
                .map_or(0.0, |(_, middle)| middle.distance(a00.midpoint(a11)));
            for (corners, params) in [
                ([a00, a10, a11], [p00, p10, p11]),
                ([a00, a11, a01], [p00, p11, p01]),
            ] {
                let low = Point::new(
                    corners.iter().map(|p| p.x).fold(f64::MAX, f64::min),
                    corners.iter().map(|p| p.y).fold(f64::MAX, f64::min),
                    corners.iter().map(|p| p.z).fold(f64::MAX, f64::min),
                );
                let high = Point::new(
                    corners.iter().map(|p| p.x).fold(f64::MIN, f64::max),
                    corners.iter().map(|p| p.y).fold(f64::MIN, f64::max),
                    corners.iter().map(|p| p.z).fold(f64::MIN, f64::max),
                );
                out.push(Cell {
                    corners,
                    at: p00,
                    low,
                    high,
                    sag,
                    params,
                });
            }
        }
        below = above;
    }
    out
}

/// Whether two cells' boxes come within a margin of each other.
/// Sampled cells binned by their boxes on a coarse grid, for finding the
/// cells whose boxes could meet a given one without testing them all.
struct CellBins {
    low: Point,
    size: f64,
    counts: [usize; 3],
    bins: Vec<Vec<usize>>,
    margin: f64,
}

impl CellBins {
    /// At most this many bins a side: a grid on a clamped plane's million
    /// units would otherwise be all empty bins.
    const MOST: usize = 48;

    fn over(cells: &[Cell], margin: f64) -> Self {
        let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for c in cells {
            low = Point::new(low.x.min(c.low.x), low.y.min(c.low.y), low.z.min(c.low.z));
            high = Point::new(
                high.x.max(c.high.x),
                high.y.max(c.high.y),
                high.z.max(c.high.z),
            );
        }
        let extent = (high - low).magnitude();
        #[allow(clippy::cast_precision_loss)]
        let size = if extent.is_finite() && extent > 0.0 {
            (extent / Self::MOST as f64).max(margin)
        } else {
            1.0
        };
        let count = |lo: f64, hi: f64| -> usize {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let n = ((hi - lo) / size).floor() as usize + 1;
            n.clamp(1, Self::MOST + 1)
        };
        let counts = if extent.is_finite() {
            [
                count(low.x, high.x),
                count(low.y, high.y),
                count(low.z, high.z),
            ]
        } else {
            [1, 1, 1]
        };
        let mut bins = vec![Vec::new(); counts[0] * counts[1] * counts[2]];
        let mut this = Self {
            low,
            size,
            counts,
            bins: Vec::new(),
            margin,
        };
        for (i, c) in cells.iter().enumerate() {
            this.each_bin(c.low, c.high, |k| bins[k].push(i));
        }
        this.bins = bins;
        this
    }

    /// Every bin a box grown by the margin touches.
    fn each_bin(&self, low: Point, high: Point, mut visit: impl FnMut(usize)) {
        let index = |x: f64, lo: f64, n: usize| -> usize {
            if !x.is_finite() {
                return if x > 0.0 { n - 1 } else { 0 };
            }
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let k = ((x - lo) / self.size).floor().max(0.0) as usize;
            k.min(n - 1)
        };
        let m = self.margin;
        let [nx, ny, nz] = self.counts;
        let (x0, x1) = (
            index(low.x - m, self.low.x, nx),
            index(high.x + m, self.low.x, nx),
        );
        let (y0, y1) = (
            index(low.y - m, self.low.y, ny),
            index(high.y + m, self.low.y, ny),
        );
        let (z0, z1) = (
            index(low.z - m, self.low.z, nz),
            index(high.z + m, self.low.z, nz),
        );
        for x in x0..=x1 {
            for y in y0..=y1 {
                for z in z0..=z1 {
                    visit((x * ny + y) * nz + z);
                }
            }
        }
    }

    /// The cells whose boxes could meet `cell`'s, lowest index first.
    fn candidates(&self, cell: &Cell, out: &mut Vec<usize>) {
        out.clear();
        self.each_bin(cell.low, cell.high, |k| {
            out.extend_from_slice(&self.bins[k])
        });
        out.sort_unstable();
        out.dedup();
    }
}

fn overlap(a: &Cell, b: &Cell, margin: f64) -> bool {
    a.low.x <= b.high.x + margin
        && b.low.x <= a.high.x + margin
        && a.low.y <= b.high.y + margin
        && b.low.y <= a.high.y + margin
        && a.low.z <= b.high.z + margin
        && b.low.z <= a.high.z + margin
}

/// A point where two triangles cross, if they do.
///
/// Each triangle's edges are tested against the other's plane and then against
/// the triangle itself. Only an approximate answer is needed: it is a seed, and
/// the Newton correction that follows is what makes it a point on the curve.
fn triangles_cross(a: &Cell, b: &Cell) -> Option<Point> {
    for (edges, target) in [(a, b), (b, a)] {
        for k in 0..3 {
            let (from, to) = (edges.corners[k], edges.corners[(k + 1) % 3]);
            if let Some(hit) = segment_meets_triangle(from, to, target.corners) {
                return Some(hit);
            }
        }
    }
    None
}

/// Where a segment crosses a triangle.
pub(crate) fn segment_meets_triangle(from: Point, to: Point, t: [Point; 3]) -> Option<Point> {
    let direction = to - from;
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let h = direction.cross(e2);
    let determinant = e1.dot(h);
    if determinant.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let inverse = 1.0 / determinant;
    let s = from - t[0];
    let u = inverse * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = inverse * direction.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let along = inverse * e2.dot(q);
    if !(0.0..=1.0).contains(&along) {
        return None;
    }
    Some(from + direction * along)
}

/// One arc between branch points, cut from a stalled fragment.
struct Arc {
    points: Vec<Point>,
    on_a: Vec<(f64, f64)>,
    on_b: Vec<(f64, f64)>,
    /// The branch-point cluster each end attaches to, if any.
    head_bp: Option<usize>,
    tail_bp: Option<usize>,
}

impl Arc {
    fn length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|pair| pair[0].distance(pair[1]))
            .sum()
    }

    /// The direction the arc leaves its end, read across a window deep
    /// enough to stand clear of any residual wander.
    fn outgoing(&self, tail: bool) -> Option<Vector> {
        let n = self.points.len();
        if n < 2 {
            return None;
        }
        let window = (n - 1).min(24);
        let (at, back) = if tail {
            (n - 1, n - 1 - window)
        } else {
            (0, window)
        };
        let out = self.points[at] - self.points[back];
        let m = out.magnitude();
        (m > f64::MIN_POSITIVE).then(|| out / m)
    }
}

/// Whether a stalled branch is transversal *somewhere* in its interior:
/// the certificate that it is a curve passing branch points rather than
/// tangential-contact debris. A plane resting on a torus produces
/// fragments tangent along their whole length; a real curve through a
/// branch point is tangent only in passing.
fn interior_is_transversal(
    branch: &Traced,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    let n = branch.points.len();
    if n < 5 {
        return false;
    }
    [n / 4, n / 2, 3 * n / 4]
        .into_iter()
        .any(|i| crossing_sine(a, b, branch.on_a[i], branch.on_b[i], tol) > BRANCH_POINT_SINE)
}

/// Join stalled fragments that meet at branch points into the curves they
/// belong to.
///
/// Where the two normals become parallel the intersection has no single
/// direction: the walk stalls there, wanders in place while the correction
/// gives out, and a curve that passes *through* the singularity comes back
/// as fragments with parked ends. The reassembly is geometric, not
/// bookkeeping: cluster the near-tangent stall ends into branch points;
/// cut every fragment at its branch-point visits, which trims the wander
/// and separates the arcs a walk-through glued together; drop debris and
/// duplicate coverage; then, at each branch point, pair arc ends whose
/// tangents continue one another (a smooth curve crosses the singularity
/// collinearly, and the crossing curve turns through the crossing angle)
/// and chain the pairs into whole curves, closing the loops that close.
///
/// Only fragments transversal somewhere in their interior participate:
/// tangential contact along a whole curve is a different phenomenon and
/// keeps its honest fragments.
fn stitch_stalled(
    found: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<Traced> {
    let reach = options.chord.max(tol.confusion()) * 60.0;
    const CONTINUES: f64 = 0.5;

    let (candidates, mut out): (Vec<Traced>, Vec<Traced>) = found.into_iter().partition(|branch| {
        branch.stopped == Stopped::Stalled && interior_is_transversal(branch, a, b, tol)
    });
    if candidates.is_empty() {
        return out;
    }

    // Branch points: the near-tangent stall ends, clustered.
    let mut bps: Vec<Point> = Vec::new();
    for branch in &candidates {
        let n = branch.points.len();
        for at in [0, n - 1] {
            if crossing_sine(a, b, branch.on_a[at], branch.on_b[at], tol) < BRANCH_POINT_SINE {
                let p = branch.points[at];
                if !bps.iter().any(|held| held.distance(p) <= reach) {
                    bps.push(p);
                }
            }
        }
    }
    if bps.is_empty() {
        out.extend(candidates);
        return out;
    }
    let bp_of =
        |p: Point| -> Option<usize> { bps.iter().position(|held| held.distance(p) <= reach) };

    // Cut each fragment at its branch-point visits: maximal runs of points
    // clear of every branch point become arcs, attached to the branch
    // points beside them.
    let mut arcs: Vec<Arc> = Vec::new();
    for branch in &candidates {
        let n = branch.points.len();
        let mut run_start: Option<usize> = None;
        for i in 0..=n {
            let near = i < n && bp_of(branch.points[i]).is_some();
            match (run_start, near, i == n) {
                (None, false, false) => run_start = Some(i),
                (Some(s), true, _) | (Some(s), _, true) => {
                    let e = i;
                    if e > s + 1 {
                        let head_bp = if s > 0 {
                            bp_of(branch.points[s - 1])
                        } else {
                            None
                        };
                        let tail_bp = if e < n { bp_of(branch.points[e]) } else { None };
                        arcs.push(Arc {
                            points: branch.points[s..e].to_vec(),
                            on_a: branch.on_a[s..e].to_vec(),
                            on_b: branch.on_b[s..e].to_vec(),
                            head_bp,
                            tail_bp,
                        });
                    }
                    run_start = None;
                }
                _ => {}
            }
        }
    }

    // Debris and duplicate coverage out; longest first so the fuller
    // tracing of a doubly-walked arc is the one kept.
    arcs.retain(|arc| arc.length() > options.chord * 10.0 && arc.points.len() >= 4);
    arcs.sort_by(|x, y| {
        y.length()
            .partial_cmp(&x.length())
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let mut kept: Vec<Arc> = Vec::new();
    'candidate: for arc in arcs {
        let n = arc.points.len();
        for probe in [n / 4, n / 2, 3 * n / 4] {
            let p = arc.points[probe];
            if kept.iter().any(|held| {
                held.points
                    .windows(2)
                    .any(|pair| distance_to_segment(p, pair[0], pair[1]) <= reach)
            }) {
                continue 'candidate;
            }
        }
        kept.push(arc);
    }

    // Pair arc ends at each branch point by tangent continuation: the
    // smooth curve runs straight through, the crossing one turns.
    let ends: Vec<(usize, bool, usize, Vector)> = kept
        .iter()
        .enumerate()
        .flat_map(|(i, arc)| {
            [(false, arc.head_bp), (true, arc.tail_bp)]
                .into_iter()
                .filter_map(move |(tail, bp)| Some((i, tail, bp?, arc.outgoing(tail)?)))
        })
        .collect();
    let mut partner: Vec<Option<usize>> = vec![None; ends.len()];
    for bp in 0..bps.len() {
        loop {
            let mut best: Option<(usize, usize, f64)> = None;
            for x in 0..ends.len() {
                if partner[x].is_some() || ends[x].2 != bp {
                    continue;
                }
                for y in (x + 1)..ends.len() {
                    if partner[y].is_some() || ends[y].2 != bp {
                        continue;
                    }
                    let score = -ends[x].3.dot(ends[y].3);
                    if score > CONTINUES && best.is_none_or(|(_, _, held)| score > held) {
                        best = Some((x, y, score));
                    }
                }
            }
            let Some((x, y, _)) = best else { break };
            partner[x] = Some(y);
            partner[y] = Some(x);
        }
    }

    // Chain the arcs through the pairings into whole curves.
    let end_index = |arc: usize, tail: bool| -> Option<usize> {
        ends.iter().position(|e| e.0 == arc && e.1 == tail)
    };
    let mut used = vec![false; kept.len()];
    for start in 0..kept.len() {
        if used[start] {
            continue;
        }
        // Walk backwards first to a free entry, unless the chain loops.
        let mut first = start;
        let mut first_reversed = false;
        let mut seen_back = vec![false; kept.len()];
        loop {
            seen_back[first] = true;
            // Forward traversal enters an arc at its head, reversed at its
            // tail, so the entry end is named by the orientation flag.
            let Some(entry) = end_index(first, first_reversed) else {
                break;
            };
            let Some(p) = partner[entry] else { break };
            let (prev, prev_tail, _, _) = ends[p];
            if seen_back[prev] {
                break; // the chain is a loop; any start serves
            }
            first = prev;
            // The previous arc *leaves* through the paired end: leaving at
            // its tail means it runs forward.
            first_reversed = !prev_tail;
        }

        // Now walk forwards from `first`, consuming arcs.
        let mut points: Vec<Point> = Vec::new();
        let mut on_a: Vec<(f64, f64)> = Vec::new();
        let mut on_b: Vec<(f64, f64)> = Vec::new();
        let mut current = first;
        let mut reversed = first_reversed;
        let mut closed = false;
        loop {
            used[current] = true;
            let arc = &kept[current];
            type Run = (Vec<Point>, Vec<(f64, f64)>, Vec<(f64, f64)>);
            let (pts, pa, pb): Run = if reversed {
                (
                    arc.points.iter().rev().copied().collect(),
                    arc.on_a.iter().rev().copied().collect(),
                    arc.on_b.iter().rev().copied().collect(),
                )
            } else {
                (arc.points.clone(), arc.on_a.clone(), arc.on_b.clone())
            };
            // Insert the branch point itself at the junction.
            if !points.is_empty() {
                let joint_bp = if reversed { arc.tail_bp } else { arc.head_bp };
                if let Some(bp) = joint_bp {
                    points.push(bps[bp]);
                    on_a.push(pa[0]);
                    on_b.push(pb[0]);
                }
            }
            points.extend(pts);
            on_a.extend(pa);
            on_b.extend(pb);

            let leaving = end_index(current, !reversed);
            let Some(l) = leaving else { break };
            let Some(p) = partner[l] else { break };
            let (next, next_tail, _, _) = ends[p];
            if used[next] {
                closed = next == first;
                break;
            }
            current = next;
            reversed = next_tail;
        }
        if closed && points.len() > 3 {
            let bridge = points[0];
            let ba = on_a[0];
            let bb = on_b[0];
            points.push(bridge);
            on_a.push(ba);
            on_b.push(bb);
        }
        out.push(Traced {
            points,
            on_a,
            on_b,
            stopped: if closed {
                Stopped::Closed
            } else {
                Stopped::Stalled
            },
        });
    }
    out
}

/// A point where the two surfaces touch, their normals parallel, and the
/// section crosses itself.
#[derive(Debug, Clone, Copy)]
struct Touch {
    point: Point,
    on_a: (f64, f64),
    on_b: (f64, f64),
    /// How far the second surface stands off the first along its normal
    /// there, within the confusion distance.
    apart: f64,
    /// How far the solve for the touch may land from it: where a gap the
    /// size of its residual opens, `√(2r / |λ|)` for the least eigenvalue
    /// `λ` of `Q`.
    blur: f64,
    /// The two lines the section crosses itself along, where they stand
    /// far enough apart to be followed one by one: see [`touch_arms`].
    arms: Option<[Vector; 2]>,
    /// How far round the touch a sample counts as on it: see
    /// [`settle_touches`].
    radius: f64,
    /// How far round the touch a stalled end or a step between its arms is
    /// taken as heading for it: see [`settle_touches`].
    around: f64,
    /// How fast the surfaces part across an arm, `√(−det Q)` of
    /// [`touch_arms`]: a sample `ρ` out along an arm with the surfaces
    /// within `ε` of each other stands up to `ε / (ρ · spread)` off it.
    spread: f64,
    /// How far from the touch a straight chord to it stays on both
    /// surfaces within the confusion distance `ε`: `√(2ε / k)`, `k` the
    /// largest bend of either surface or of their difference there.
    snap: f64,
}

/// The two lines along which the section leaves a touch, as unit vectors
/// in the shared tangent plane.
///
/// Near the touch each surface stands off the tangent plane by half its
/// second fundamental form, so the section is where the two forms agree:
/// the null lines of their difference `Q`. Where `Q` is indefinite there
/// are two, and the section crosses itself along them; where it is
/// definite the touch is an isolated point, and where it is singular the
/// surfaces touch along a curve. `None` unless there are two; with them,
/// how fast the surfaces part across either, `√(−det Q)`, the largest bend
/// of either surface or of `Q`, the spectral radius of the three forms,
/// and the eigenvalues of `Q`, least first.
fn touch_arms(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> Option<([Vector; 2], f64, f64, [f64; 2])> {
    let n = normal_at(a, on_a, tol)?;
    let (au, _) = a.d1_at(on_a.0, on_a.1, tol).ok()?;
    let e1 = au - n * au.dot(n);
    let e1 = e1 / e1.magnitude();
    let e2 = n.cross(e1);
    // The second fundamental form of one surface along `n`, in the basis
    // (e1, e2): each basis vector written in the chart through the first
    // fundamental form, then the second derivatives read off along `n`.
    let form = |surface: &SurfaceGeometry, (u, v): (f64, f64)| -> Option<[f64; 3]> {
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let (suu, suv, svv) = surface.d2_at(u, v, tol).ok()?;
        let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = e.mul_add(g, -(f * f));
        if det <= f64::MIN_POSITIVE {
            return None;
        }
        let chart = |w: Vector| {
            let (p, q) = (su.dot(w), sv.dot(w));
            (g.mul_add(p, -(f * q)) / det, e.mul_add(q, -(f * p)) / det)
        };
        let (x, y) = (chart(e1), chart(e2));
        let second = |s: (f64, f64), t: (f64, f64)| {
            n.dot(suu) * s.0 * t.0
                + n.dot(suv) * s.0.mul_add(t.1, s.1 * t.0)
                + n.dot(svv) * s.1 * t.1
        };
        Some([second(x, x), second(x, y), second(y, y)])
    };
    let (fa, fb) = (form(a, on_a)?, form(b, on_b)?);
    let q = [fa[0] - fb[0], fa[1] - fb[1], fa[2] - fb[2]];
    // Along the line at angle θ to e1 the difference is A + B cos 2θ +
    // C sin 2θ = A + R cos(2θ − φ), zero at 2θ = φ ± acos(−A / R): two
    // lines where |A| < R, the form indefinite.
    let mean = f64::midpoint(q[0], q[2]);
    let (half_gap, skew) = ((q[0] - q[2]) * 0.5, q[1]);
    let swing = half_gap.hypot(skew);
    if mean.abs().partial_cmp(&swing) != Some(core::cmp::Ordering::Less) {
        return None;
    }
    let (phi, open) = (skew.atan2(half_gap), (-mean / swing).acos());
    let line = |theta: f64| e1 * theta.cos() + e2 * theta.sin();
    let arms = [line((phi + open) * 0.5), line((phi - open) * 0.5)];
    let bend = |m: [f64; 3]| f64::midpoint(m[0], m[2]).abs() + ((m[0] - m[2]) * 0.5).hypot(m[1]);
    Some((
        arms,
        swing.mul_add(swing, -(mean * mean)).sqrt(),
        bend(fa).max(bend(fb)).max(bend(q)),
        [mean - swing, mean + swing],
    ))
}

/// Whether the section crosses itself at a touch, not just passes it: the
/// surfaces stand `apart` along the normal there, so the section is where
/// half of `Q` makes that up, a hyperbola whose branches pass
/// `√(2 |apart| / |λ|)` from the touch, `λ` the eigenvalue of `Q` of the
/// sign of `apart`. Two surfaces within the confusion distance of each
/// other at a near tangency leave such a point where nothing crosses; the
/// section crosses at the touch where the branches pass within the chord
/// `near` of it, or where `apart` is within the residual the solve for the
/// touch settles to, so no different from none.
fn crosses(
    touch: &Touch,
    (_, _, _, eigen): ([Vector; 2], f64, f64, [f64; 2]),
    near: f64,
    tol: Tolerances,
) -> bool {
    let lambda = if touch.apart >= 0.0 {
        eigen[1]
    } else {
        eigen[0]
    };
    touch.apart.abs() <= TOUCH_RESIDUAL * tol.confusion()
        || (lambda.abs() > 0.0 && 2.0 * touch.apart.abs() / lambda.abs() <= near * near)
}

/// How far a solve for a touch may land from it, from the eigenvalues of
/// `Q` there: see [`Touch::blur`].
fn blur_of(eigen: [f64; 2], tol: Tolerances) -> f64 {
    let least = eigen[0].abs().min(eigen[1].abs());
    if least > 0.0 {
        (2.0 * TOUCH_RESIDUAL * tol.confusion() / least).sqrt()
    } else {
        0.0
    }
}

/// Whether two solved touches are one: they stand within the chord `near`,
/// or within where either solve may land.
fn same_touch(t: &Touch, u: &Touch, near: f64) -> bool {
    t.point.distance(u.point) <= near.max(t.blur).max(u.blur)
}

/// Whether two arms stand far enough apart that a sample `radius` out
/// along one lies more than twice `near`, the march's chord, off the
/// other, so each can be told from the other and followed on its own.
fn apart(arms: [Vector; 2], radius: f64, near: f64) -> bool {
    radius * arms[0].cross(arms[1]).magnitude() > 2.0 * near
}

/// The residual a solve for a touch settles to, as a fraction of the
/// confusion distance.
const TOUCH_RESIDUAL: f64 = 1e-3;

/// Where the two surfaces touch near a contact: on both, normals parallel.
///
/// Posed as a pair of points standing apart by `t` along the first
/// surface's normal, with the second surface's tangent plane square to that
/// normal: five equations in the four parameters and `t`, regular where the
/// surfaces bend apart differently across the touch, as they do where a
/// section crosses itself. A touch is a solution with `t` inside the
/// confusion distance; a pair standing further apart there is a near miss,
/// which the march follows as it is. `None` where the solve does not settle.
fn touching_point(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> Option<Touch> {
    let eval = |x: &[f64]| -> Option<[f64; 5]> {
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        let pa = a.point_at(ua, va, tol).ok()?;
        let pb = b.point_at(ub, vb, tol).ok()?;
        let na = normal_at(a, (ua, va), tol)?;
        let (bu, bv) = b.d1_at(ub, vb, tol).ok()?;
        let (lu, lv) = (bu.magnitude(), bv.magnitude());
        if lu <= tol.confusion() || lv <= tol.confusion() {
            return None;
        }
        let gap = pa + na * x[4] - pb;
        Some([gap.x, gap.y, gap.z, na.dot(bu) / lu, na.dot(bv) / lv])
    };
    // Forward differences: the Jacobian's own error slows the iteration and
    // moves no root. They are taken at the point the residual measures,
    // held to the charts, and step back into a chart where a step forward
    // would leave it: a touch on a bound that does not repeat (a profile
    // ending on its apex) draws the iteration onto the bound and past it,
    // and a step clamped back onto the bound would read that parameter as
    // moving nothing, the Jacobian as singular.
    const STEP: f64 = 1e-7;
    let step_from = |k: usize, at: f64| -> f64 {
        let surface = if k < 2 { a } else { b };
        let ((_, ub), (_, vb)) = surface.domain();
        let (periodic, high) = if matches!(k, 0 | 2) {
            (surface.is_periodic_u(), ub)
        } else {
            (surface.is_periodic_v(), vb)
        };
        if periodic || at + STEP <= high {
            STEP
        } else {
            -STEP
        }
    };
    let system = |x: &[f64; 5]| {
        let failed = || ([f64::INFINITY; 5], [[0.0; 5]; 5]);
        let mut at = *x;
        (at[0], at[1]) = clamp(a, x[0], x[1]);
        (at[2], at[3]) = clamp(b, x[2], x[3]);
        let Some(here) = eval(&at) else {
            return failed();
        };
        let mut jacobian = [[0.0; 5]; 5];
        for k in 0..5 {
            let step = if k < 4 { step_from(k, at[k]) } else { STEP };
            let mut moved = at;
            moved[k] += step;
            let Some(there) = eval(&moved) else {
                return failed();
            };
            for (row, (t, h)) in jacobian.iter_mut().zip(there.iter().zip(here.iter())) {
                row[k] = (t - h) / step;
            }
        }
        (here, jacobian)
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * TOUCH_RESIDUAL,
        step: tol.parametric() * 1e-3,
        max_iterations: 50,
    };
    let start = [on_a.0, on_a.1, on_b.0, on_b.1, 0.0];
    let (x, norm, _, _) = solve::newton_system_fixed(system, start, criteria).ok()?;
    if norm > tol.confusion() || x[4].abs() > tol.confusion() {
        return None;
    }
    let (on_a, on_b) = (clamp(a, x[0], x[1]), clamp(b, x[2], x[3]));
    let point = a.point_at(on_a.0, on_a.1, tol).ok()?;
    if b.point_at(on_b.0, on_b.1, tol).ok()?.distance(point) > tol.confusion() {
        return None;
    }
    Some(Touch {
        point,
        on_a,
        on_b,
        apart: x[4],
        blur: 0.0,
        arms: None,
        radius: 0.0,
        around: 0.0,
        spread: 0.0,
        snap: 0.0,
    })
}

/// A rectangle of a surface's chart, sampled on a 3 × 3 grid: where the
/// surface lies, within how far it bows from the samples, and which way it
/// faces.
#[derive(Debug, Clone, Copy)]
struct Patch {
    u: (f64, f64),
    v: (f64, f64),
    /// The box of the samples, grown by the bow.
    low: Point,
    high: Point,
    /// The normal at the middle, and the angle every normal of the patch
    /// stays within of it.
    axis: Vector,
    spread: f64,
    /// The diagonal of the box of the samples.
    size: f64,
    /// The normals at the samples.
    normals: [Vector; 9],
    /// Whether the patch is a piece of a plane: its samples in one plane
    /// and its normals one direction, so no quarter of it says more.
    flat: bool,
    /// Across `u` and across `v`: how far the samples run, and how far the
    /// normals turn.
    across: [(f64, f64); 2],
}

impl Patch {
    /// The patch over `u` × `v`, or `None` where the surface cannot be
    /// evaluated or has no normal there.
    fn over(
        surface: &SurfaceGeometry,
        u: (f64, f64),
        v: (f64, f64),
        tol: Tolerances,
    ) -> Option<Self> {
        let mut points = [[Point::ORIGIN; 3]; 3];
        let mut normals = [[Vector::new(0.0, 0.0, 0.0); 3]; 3];
        for (i, row) in points.iter_mut().enumerate() {
            for (j, point) in row.iter_mut().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let (s, t) = (i as f64 * 0.5, j as f64 * 0.5);
                let (pu, pv) = (u.0 + (u.1 - u.0) * s, v.0 + (v.1 - v.0) * t);
                *point = surface.point_at(pu, pv, tol).ok()?;
                normals[i][j] = normal_at(surface, (pu, pv), tol)?;
            }
        }
        // How far the surface bows from the samples: the middle samples'
        // distance from the straight blend of the ones either side, the
        // sag of a quadratic through them, taken twice over as the
        // seeding's cells take it.
        let mid = |p: Point, q: Point| p + (q - p) * 0.5;
        let mut bow: f64 = 0.0;
        for (k, row) in points.iter().enumerate() {
            bow = bow
                .max(points[1][k].distance(mid(points[0][k], points[2][k])))
                .max(row[1].distance(mid(row[0], row[2])));
        }
        let grow = BOW * bow;
        let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in points.iter().flatten() {
            low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
            high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
        }
        let pad = Vector::new(grow, grow, grow);
        let axis = normals[1][1];
        // The normals turn between samples by as much again as they turn
        // across them, so the spread is taken twice over too.
        let mut turn: f64 = 0.0;
        for n in normals.iter().flatten() {
            turn = turn.max(axis.dot(*n).clamp(-1.0, 1.0).acos());
        }
        let mut flat = [axis; 9];
        for (k, n) in normals.iter().flatten().enumerate() {
            flat[k] = *n;
        }
        let angle = |m: Vector, n: Vector| m.dot(n).clamp(-1.0, 1.0).acos();
        // Taken through the middle sample, as a patch over a whole closed
        // chart ends where it starts.
        let mut across = [(0.0_f64, 0.0_f64); 2];
        for k in 0..3 {
            let (p, n) = (&points, &normals);
            across[0].0 = across[0]
                .0
                .max(p[0][k].distance(p[1][k]) + p[1][k].distance(p[2][k]));
            across[0].1 = across[0]
                .1
                .max(angle(n[0][k], n[1][k]) + angle(n[1][k], n[2][k]));
            across[1].0 = across[1]
                .0
                .max(p[k][0].distance(p[k][1]) + p[k][1].distance(p[k][2]));
            across[1].1 = across[1]
                .1
                .max(angle(n[k][0], n[k][1]) + angle(n[k][1], n[k][2]));
        }
        Some(Self {
            u,
            v,
            low: low - pad,
            high: high + pad,
            axis,
            spread: BOW * turn,
            size: (high - low).magnitude(),
            normals: flat,
            flat: grow == 0.0 && turn == 0.0,
            across,
        })
    }

    /// The way the patch is halved next, if it still runs further than
    /// `near` any way, and how much that says: whether the normals turn
    /// that way, and by how much, or else how far it runs. Where two
    /// patches could touch, the normals that turn across a patch are what
    /// must part for the pair to be set aside, so the patch is halved
    /// across the way they turn most, a thin strip along a line of
    /// near-tangency parting from its neighbours without being cut along
    /// its length; the ways the normals do not turn are halved once those
    /// that do run within `near`.
    fn cut(&self, near: f64) -> Option<(usize, (bool, f64))> {
        if self.flat {
            return None;
        }
        let open = |k: usize| self.across[k].0 > near;
        // Halved across the way that says most, or quartered where the
        // other says at least half as much.
        let pick = |says: &dyn Fn(usize) -> f64| {
            let k = (0..2)
                .filter(|&k| open(k) && says(k) > 0.0)
                .max_by(|&x, &y| says(x).total_cmp(&says(y)))?;
            let other = 1 - k;
            let way = if open(other) && says(other) * 2.0 >= says(k) {
                2
            } else {
                k
            };
            Some((way, says(k)))
        };
        if let Some((way, turn)) = pick(&|k| self.across[k].1) {
            return Some((way, (true, turn)));
        }
        pick(&|k| self.across[k].0).map(|(way, run)| (way, (false, run)))
    }

    /// The two halves of the patch across `u` (`way` 0) or `v` (1), or its
    /// four quarters (2).
    fn halves(&self, way: usize, surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Self> {
        let (um, vm) = self.middle();
        let (u0, u1) = ((self.u.0, um), (um, self.u.1));
        let (v0, v1) = ((self.v.0, vm), (vm, self.v.1));
        let parts = match way {
            0 => vec![(u0, self.v), (u1, self.v)],
            1 => vec![(self.u, v0), (self.u, v1)],
            _ => vec![(u0, v0), (u1, v0), (u0, v1), (u1, v1)],
        };
        parts
            .into_iter()
            .filter_map(|(u, v)| Self::over(surface, u, v, tol))
            .collect()
    }

    /// The middle of the patch in its chart.
    fn middle(&self) -> (f64, f64) {
        (
            f64::midpoint(self.u.0, self.u.1),
            f64::midpoint(self.v.0, self.v.1),
        )
    }
}

/// Whether no normal of one patch can be parallel to a normal of the
/// other although their cones meet: some component of the cross product
/// of the normals keeps one sign over every pair of samples, its value
/// at the middles further from zero than it strays across the samples,
/// taken twice over as the bow is.
fn normals_part(pa: &Patch, pb: &Patch) -> bool {
    let middle = pa.normals[4].cross(pb.normals[4]);
    let mut stray = [0.0_f64; 3];
    for na in &pa.normals {
        for nb in &pb.normals {
            let c = na.cross(*nb);
            stray[0] = stray[0].max((c.x - middle.x).abs());
            stray[1] = stray[1].max((c.y - middle.y).abs());
            stray[2] = stray[2].max((c.z - middle.z).abs());
        }
    }
    [middle.x, middle.y, middle.z]
        .iter()
        .zip(stray)
        .any(|(m, d)| m.abs() > BOW * d)
}

/// A surface's patches, and for each the indices of its halves once cut.
type PatchTree = (Vec<Patch>, Vec<Option<Vec<usize>>>);

/// The crossings where the two surfaces touch, found whatever the march's
/// stride: both charts are cut into patches, and each pair is set aside
/// where their boxes do not meet, or where no normal of one can be
/// parallel to a normal of the other, their normal cones apart; every
/// other pair is halved, across the way the normals turn most (see
/// [`Patch::cut`]), until both stand within the chord, and the touch is
/// solved from the middles of what is left. Only touches the section crosses itself at, with two arms, are
/// kept. `None` where the charts are not bounded or the pairs to weigh
/// would outnumber the points a branch may hold, so the search is not
/// complete.
fn subdivided_touches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Option<Vec<Touch>> {
    search_touches(a, b, options, tol, true).or_else(|| search_touches(a, b, options, tol, false))
}

/// The search of [`subdivided_touches`], each patch halved across the way
/// its normals turn most where `across` holds (see [`Patch::cut`]), and
/// otherwise quartered, the larger first, until it stands within the
/// chord. Halving parts a strip along a line of near-tangency from its
/// neighbours without cutting it along its length; where the normals turn
/// alike every way, as about many crossings close together, quartering
/// weighs fewer pairs, and the search is run so where halving does not
/// complete.
fn search_touches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
    across: bool,
) -> Option<Vec<Touch>> {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    let whole = |surface: &SurfaceGeometry| {
        let (u, v) = surface.domain();
        [u.0, u.1, v.0, v.1]
            .iter()
            .all(|x| x.is_finite())
            .then(|| Patch::over(surface, u, v, tol))
            .flatten()
    };
    let cut_of = |p: &Patch| {
        if across {
            p.cut(near)
        } else {
            (!p.flat && p.size > near).then_some((2, (false, p.size)))
        }
    };
    // Each surface's patches as a tree, each halved once however many
    // patches of the other it is weighed against.
    let mut trees: [PatchTree; 2] = [(vec![whole(a)?], vec![None]), (vec![whole(b)?], vec![None])];
    let halve = |tree: &mut PatchTree, surface: &SurfaceGeometry, k: usize| -> Vec<usize> {
        if let Some(children) = &tree.1[k] {
            return children.clone();
        }
        let start = tree.0.len();
        let Some((way, _)) = cut_of(&tree.0[k]) else {
            return Vec::new();
        };
        let halves = tree.0[k].halves(way, surface, tol);
        let count = halves.len();
        tree.0.extend(halves);
        tree.1.extend(core::iter::repeat_with(|| None).take(count));
        let children: Vec<usize> = (start..start + count).collect();
        tree.1[k] = Some(children.clone());
        children
    };
    let mut pairs = vec![(0_usize, 0_usize)];
    let mut touches: Vec<Touch> = Vec::new();
    let mut weighed = 0_usize;
    while let Some((ia, ib)) = pairs.pop() {
        weighed += 1;
        if weighed > options.max_points {
            return None;
        }
        let (pa, pb) = (trees[0].0[ia], trees[1].0[ib]);
        let meet = pa.low.x <= pb.high.x
            && pb.low.x <= pa.high.x
            && pa.low.y <= pb.high.y
            && pb.low.y <= pa.high.y
            && pa.low.z <= pb.high.z
            && pb.low.z <= pa.high.z;
        if !meet {
            continue;
        }
        let facing = pa.axis.dot(pb.axis).abs().clamp(-1.0, 1.0).acos();
        if facing > pa.spread + pb.spread || normals_part(&pa, &pb) {
            continue;
        }
        let (cut_a, cut_b) = (cut_of(&pa), cut_of(&pb));
        if cut_a.is_none() && cut_b.is_none() {
            // A touch already found within the pair is the one there:
            // touches stand further apart than the chord.
            let held = touches.iter().any(|t| {
                t.point.x >= pa.low.x - near
                    && t.point.x <= pa.high.x + near
                    && t.point.y >= pa.low.y - near
                    && t.point.y <= pa.high.y + near
                    && t.point.z >= pa.low.z - near
                    && t.point.z <= pa.high.z + near
            });
            if held {
                continue;
            }
            // A flat patch is started from beside the other's middle.
            let (start_a, start_b) = match (pa.flat, pb.flat) {
                (true, false) => {
                    let at = b.point_at(pb.middle().0, pb.middle().1, tol).ok();
                    let on = at.and_then(|p| nearest_on(a, pa.middle(), p, tol));
                    (on.map_or(pa.middle(), |(uv, _)| uv), pb.middle())
                }
                (false, true) => {
                    let at = a.point_at(pa.middle().0, pa.middle().1, tol).ok();
                    let on = at.and_then(|p| nearest_on(b, pb.middle(), p, tol));
                    (pa.middle(), on.map_or(pb.middle(), |(uv, _)| uv))
                }
                _ => (pa.middle(), pb.middle()),
            };
            let Some(mut touch) = touching_point(a, b, start_a, start_b, tol) else {
                continue;
            };
            let Some(form) = touch_arms(a, b, touch.on_a, touch.on_b, tol) else {
                continue;
            };
            if !crosses(&touch, form, near, tol) || !apart(form.0, reach, near) {
                continue;
            }
            touch.arms = Some(form.0);
            touch.blur = blur_of(form.3, tol);
            if !touches.iter().any(|t| same_touch(t, &touch, near)) {
                touches.push(touch);
            }
            continue;
        }
        // The patch whose cut says more is halved: one across which the
        // normals turn most, or else the longest.
        let says = |c: Option<(usize, (bool, f64))>| c.map(|(_, w)| w);
        if says(cut_a) >= says(cut_b) {
            for q in halve(&mut trees[0], a, ia) {
                pairs.push((q, ib));
            }
        } else {
            for q in halve(&mut trees[1], b, ib) {
                pairs.push((ia, q));
            }
        }
    }
    Some(touches)
}

/// A parameter pair moved by whole periods to stand nearest `near`, so a
/// point joined onto a run of samples continues it in the chart.
fn beside(surface: &SurfaceGeometry, at: (f64, f64), near: (f64, f64)) -> (f64, f64) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let shift = |x: f64, to: f64, span: f64, periodic: bool| {
        if !periodic || !span.is_finite() || span <= 0.0 {
            return x;
        }
        x + ((to - x) / span).round() * span
    };
    (
        shift(at.0, near.0, ub - ua, surface.is_periodic_u()),
        shift(at.1, near.1, vb - va, surface.is_periodic_v()),
    )
}

/// Branches cut where the surfaces touch and the section crosses itself.
///
/// Where two surfaces touch at a point and bend apart differently across
/// it (a drill lying inside a ring's outer wall and touching it), the
/// section is a figure eight: two loops meeting at the touch, four arms
/// leaving it. The walk has no single direction there. It runs straight
/// through on one arm, so one branch goes round both loops, or it turns
/// onto the next arm and closes one loop, or it stalls a short way off
/// the touch. As one branch the figure eight is no ring a face can be
/// split by.
///
/// So each touch is found exactly from the low-angle stretches of the
/// branches, those traced `again` (a walk that turned a corner and ran on
/// over a branch already held) among them, and a branch crossing itself
/// there, turning a corner on it, ending on it or stopping short of it is
/// cut there: the samples within reach of it go and the arcs either side
/// end on it, a stalled end heading for it carried on to it. A branch
/// running once straight through each touch it passes is a smooth curve
/// another crosses there, and stays whole; where the touch's two arms are
/// known, straight through means in along one arm and on along it. An arc
/// starting and ending on one
/// touch is a whole loop, cut again so each piece is open with two
/// distinct ends: off its middle sample, which on a loop symmetric about
/// the touch lands on whatever else lies on the plane of symmetry, a seam
/// crossing another seam among them. An arm of a touch that nothing runs
/// along is then followed out from the touch, and a branch or arc lying
/// wholly on a longer one is dropped, or cut back where its ends do.
/// Branches nowhere near a touch, and those running along a tangency, are
/// returned as they came.
///
/// An arc shorter than ten chords, or of fewer than four samples, is
/// dropped as debris, an arm cut that short among them: a branch running
/// off the edge of a surface within a few chords of a touch loses the
/// stub between the two.
fn split_at_touches(
    found: Vec<Traced>,
    again: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<Traced> {
    let near = options.chord.max(tol.confusion());
    let mut touches: Vec<Touch> = Vec::new();
    let solved_ends: Vec<[Option<usize>; 2]> = found
        .iter()
        .map(|branch| find_touches(branch, &mut touches, (a, b), options, tol))
        .collect();
    for branch in &again {
        find_touches(branch, &mut touches, (a, b), options, tol);
    }
    if touches.is_empty() {
        return found;
    }
    settle_touches(&mut touches, (a, b), options, tol);
    // The march strides up to its longest step, and two touches within one
    // stride of each other may stand in one step, the section's direction
    // turning about twice in it and reading the same at both ends: a touch
    // between them goes unfound, and cutting at those found loses the
    // stretch the step covered.
    let stride = found
        .iter()
        .chain(&again)
        .flat_map(|branch| branch.points.windows(2).map(|w| w[0].distance(w[1])))
        .fold(0.0, f64::max);
    let crowded = touches.iter().enumerate().any(|(i, t)| {
        t.arms.is_some()
            && touches
                .iter()
                .enumerate()
                .any(|(j, u)| j != i && t.point.distance(u.point) <= stride)
    });
    // Every crossing of the pair is then found by cutting both charts into
    // patches, whatever the march's stride: a crossing found alone may
    // have another within a stride of it just the same. Where crossings
    // crowd and that search does not complete, the branches are cut as
    // they were traced.
    let crossing = touches.iter().any(|t| t.arms.is_some());
    let complete = crossing
        && subdivided_touches(a, b, options, tol).is_some_and(|crossings| {
            for touch in crossings {
                if !touches.iter().any(|t| same_touch(t, &touch, near)) {
                    touches.push(touch);
                }
            }
            true
        });
    if crowded && !complete {
        return split_at_touches_as_traced(found, a, b, options, tol);
    }
    if complete {
        settle_touches(&mut touches, (a, b), options, tol);
    }

    // Branches no touch cuts stay as they are, the others are cut into
    // arcs; each notes the arms of the touches it runs along.
    let mut kept: Vec<Traced> = Vec::with_capacity(found.len());
    let mut arcs: Vec<Traced> = Vec::new();
    let mut covered: Vec<(usize, usize)> = Vec::new();
    for (branch, solved) in found.into_iter().zip(solved_ends) {
        let cut = cut_at_touches(branch, solved, &touches, (a, b), options, tol);
        covered.extend(cut.covered);
        if cut.cut {
            arcs.extend(cut.pieces);
        } else {
            kept.extend(cut.pieces);
        }
    }
    let (joins, walked_out) = join_arms(&touches, &mut covered, &[], (a, b), options, tol);
    arcs.extend(joins);
    if walked_out.is_empty() && arcs.is_empty() {
        return kept;
    }
    let solved_out: Vec<[Option<usize>; 2]> = walked_out
        .iter()
        .map(|arm| find_touches(arm, &mut touches, (a, b), options, tol))
        .collect();
    settle_touches(&mut touches, (a, b), options, tol);
    let walked_any = !walked_out.is_empty();
    for (arm, solved) in walked_out.into_iter().zip(solved_out) {
        arcs.extend(cut_at_touches(arm, solved, &touches, (a, b), options, tol).pieces);
    }

    // An arc lying wholly on another, every sample clear of the touches
    // within the chord of it, is that stretch twice: longest first, each
    // is kept only if some sample of it lies off everything kept so far,
    // and cut back where its ends run along what is kept. A branch no
    // touch cut is weighed too where arms were followed out, as one may
    // run along it. Two arcs between the same two touches, leaving each
    // along the same arm, are one stretch: near the touches neither may
    // have a sample clear of them to weigh.
    let mut weighed = arcs;
    if walked_any {
        weighed.append(&mut kept);
    }
    weighed.sort_by(|x, y| arc_length(y).total_cmp(&arc_length(x)));
    let weigh = |arc: Traced, kept: &mut Vec<Traced>, touches: &[Touch]| {
        let across = across_section(&arc, touches, near, (a, b), tol);
        if kept.iter().any(|other| {
            lies_on(&arc, other, near, &across, touches) || same_stretch(&arc, other, near, touches)
        }) {
            return;
        }
        if let Some(arc) = trimmed(arc, kept, near, &across, touches) {
            kept.push(arc);
        }
    };
    for arc in weighed {
        weigh(arc, &mut kept, &touches);
    }
    // With every touch found, an arm the arcs still do not run along is
    // section left out: joined or walked out once more from what is kept.
    if complete {
        let mut covered = arms_covered(&kept, &touches);
        let ends = loose_ends(&kept, &touches, near);
        let (joins, walks) = join_arms(&touches, &mut covered, &ends, (a, b), options, tol);
        let mut more = joins;
        for arm in walks {
            more.extend(cut_at_touches(arm, [None, None], &touches, (a, b), options, tol).pieces);
        }
        more.sort_by(|x, y| arc_length(y).total_cmp(&arc_length(x)));
        for arc in more {
            weigh(arc, &mut kept, &touches);
        }
    }
    // An arc is cut back only against those kept before it, and an end on
    // a touch is never cut back: where a loose end of one runs along an arc
    // kept after it that ends on a touch, the stretch between is traced
    // twice. Each is cut back once more against all the others.
    let mut k = 0;
    while k < kept.len() {
        let arc = kept.remove(k);
        let across = across_section(&arc, &touches, near, (a, b), tol);
        if let Some(arc) = trimmed(arc, &kept, near, &across, &touches) {
            kept.insert(k, arc);
            k += 1;
        }
    }
    kept
}

/// The arms of the touches `covered` does not list, joined to a facing
/// touch along the section, or else walked out from the touch: the joins,
/// each noted as covering both arms it runs along, and the walks.
fn join_arms(
    touches: &[Touch],
    covered: &mut Vec<(usize, usize)>,
    ends: &[(Contact, Vector)],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> (Vec<Traced>, Vec<Traced>) {
    let near = options.chord.max(tol.confusion());
    let carry = near * 60.0 * 40.0;
    let mut arcs: Vec<Traced> = Vec::new();
    // Two touches facing each other along arms nothing runs along are
    // joined along the section between them: a stretch between two
    // crossings may be too shallow for a walk to set out on, and is
    // followed by samples placed where the chord calls for them, each
    // settling on the section within half its chord.
    let touch_contact = |t: &Touch| Contact {
        on_a: t.on_a,
        on_b: t.on_b,
        point: t.point,
    };
    for t in 0..touches.len() {
        if touches[t].arms.is_none() {
            continue;
        }
        for h in 0..4 {
            if covered.contains(&(t, h)) {
                continue;
            }
            let facing = (0..touches.len())
                .filter(|&u| u != t && touches[u].arms.is_some())
                .filter(|&u| half_arm(&touches[t], touches[u].point) == Some(h))
                .filter_map(|u| {
                    let back = half_arm(&touches[u], touches[t].point)?;
                    (!covered.contains(&(u, back))).then_some((u, back))
                })
                .min_by(|&(x, _), &(y, _)| {
                    let d = |u: usize| touches[u].point.distance(touches[t].point);
                    d(x).total_cmp(&d(y))
                });
            let from = touches[t].point;
            let ahead = |q: Point| {
                touches.iter().enumerate().any(|(v, w)| {
                    v != t
                        && half_arm(&touches[t], w.point) == Some(h)
                        && w.point.distance(from) < q.distance(from)
                })
            };
            // An arc that stopped short of the touch, heading at it with no
            // other touch on the way, is joined to it the same way.
            let loose = ends
                .iter()
                .filter(|(e, heading)| {
                    half_arm(&touches[t], e.point) == Some(h)
                        && heading.dot(from - e.point) > 0.0
                        && !ahead(e.point)
                })
                .min_by(|x, y| {
                    x.0.point
                        .distance(from)
                        .total_cmp(&y.0.point.distance(from))
                });
            let target = match (facing, loose) {
                (Some((u, _)), Some((e, _)))
                    if e.point.distance(from) < touches[u].point.distance(from) =>
                {
                    (*e, None)
                }
                (Some((u, back)), _) => (touch_contact(&touches[u]), Some((u, back))),
                (None, Some((e, _))) => (*e, None),
                (None, None) => continue,
            };
            let (p, q) = (touch_contact(&touches[t]), target.0);
            let q = Contact {
                on_a: beside(a, q.on_a, p.on_a),
                on_b: beside(b, q.on_b, p.on_b),
                point: q.point,
            };
            let mut budget = options.max_points;
            let Some(between) = fill((a, b), p, q, near, &mut budget, tol) else {
                continue;
            };
            let leaves = between
                .first()
                .is_none_or(|c| half_arm(&touches[t], c.point) == Some(h))
                && target.1.is_none_or(|(u, back)| {
                    between
                        .last()
                        .is_none_or(|c| half_arm(&touches[u], c.point) == Some(back))
                });
            if !leaves {
                continue;
            }
            let mut stretch = Traced {
                points: vec![p.point],
                on_a: vec![p.on_a],
                on_b: vec![p.on_b],
                stopped: Stopped::Stalled,
            };
            for c in between.iter().chain([&q]) {
                stretch.points.push(c.point);
                stretch.on_a.push(c.on_a);
                stretch.on_b.push(c.on_b);
            }
            covered.push((t, h));
            covered.extend(target.1);
            arcs.push(stretch);
        }
    }
    // An arm of a touch nothing runs along is followed out from the touch:
    // a walk turning a corner there, or stalling, leaves the arm it did not
    // take untraced, and every seed on that arm may have walked onto
    // another one first.
    let mut walked_out: Vec<Traced> = Vec::new();
    for (t, touch) in touches.iter().enumerate() {
        let Some(arms) = touch.arms else { continue };
        for (h, half) in halves(arms).into_iter().enumerate() {
            if covered.contains(&(t, h)) {
                continue;
            }
            // The walk starts where the section has a direction it can
            // follow: out from the touch's radius, twice as far each time,
            // as far as a stall is carried onto a touch.
            let Some((start, tangent)) = core::iter::successors(Some(touch.radius), |r| {
                Some(r * 2.0).filter(|&r| r <= carry)
            })
            .find_map(|out| {
                let start = along_arm(a, b, touch, half, (out, arms), options.chord, tol)?;
                Some((start, tangent_at(a, b, start, tol)?))
            }) else {
                continue;
            };
            let sense = if tangent.dot(half) >= 0.0 { 1.0 } else { -1.0 };
            let Ok(walked) = walk(a, b, start, sense, options, tol) else {
                continue;
            };
            let mut arm = Traced {
                points: vec![touch.point],
                on_a: vec![beside(a, touch.on_a, start.on_a)],
                on_b: vec![beside(b, touch.on_b, start.on_b)],
                stopped: walked.stopped,
            };
            arm.points.extend(walked.points);
            arm.on_a.extend(walked.on_a);
            arm.on_b.extend(walked.on_b);
            walked_out.push(arm);
        }
    }
    (arcs, walked_out)
}

/// The arms of the touches the arcs run along: each run of samples within
/// a touch's radius, or each step passing within it, and the samples
/// either side clear of it.
fn arms_covered(arcs: &[Traced], touches: &[Touch]) -> Vec<(usize, usize)> {
    let mut covered = Vec::new();
    for arc in arcs {
        let m = arc.points.len();
        for (t, touch) in touches.iter().enumerate() {
            if touch.arms.is_none() {
                continue;
            }
            let inside = |k: usize| arc.points[k].distance(touch.point) <= touch.radius;
            for k in 0..m {
                let passes = inside(k)
                    || (k + 1 < m
                        && distance_to_segment(touch.point, arc.points[k], arc.points[k + 1])
                            <= touch.radius);
                if !passes {
                    continue;
                }
                let before = (0..=k).rev().find(|&j| !inside(j));
                let after = (k + 1..m).find(|&j| !inside(j));
                for j in [before, after].into_iter().flatten() {
                    if let Some(h) = half_arm(touch, arc.points[j]) {
                        covered.push((t, h));
                    }
                }
            }
        }
    }
    covered
}

/// The ends of arcs that stop on no touch and meet no other end, each with
/// the way the arc runs out there.
fn loose_ends(arcs: &[Traced], touches: &[Touch], near: f64) -> Vec<(Contact, Vector)> {
    let mut ends = Vec::new();
    let all: Vec<(usize, Point)> = arcs
        .iter()
        .enumerate()
        .filter(|(_, arc)| !arc.points.is_empty())
        .flat_map(|(i, arc)| {
            [
                (2 * i, arc.points[0]),
                (2 * i + 1, arc.points[arc.points.len() - 1]),
            ]
        })
        .collect();
    for (i, arc) in arcs.iter().enumerate() {
        let m = arc.points.len();
        if m < 2 {
            continue;
        }
        for (end, (k, before)) in [(0, 1), (m - 1, m - 2)].into_iter().enumerate() {
            let p = arc.points[k];
            if touches.iter().any(|t| p.distance(t.point) <= t.radius)
                || all
                    .iter()
                    .any(|&(j, q)| j != 2 * i + end && q.distance(p) <= near)
            {
                continue;
            }
            let contact = Contact {
                on_a: arc.on_a[k],
                on_b: arc.on_b[k],
                point: p,
            };
            ends.push((contact, p - arc.points[before]));
        }
    }
    ends
}

/// Branches cut at the touches as they were traced, where the touches stand
/// within a stride of the march of each other (see [`split_at_touches`]):
/// each touch found from the low-angle samples of transversal branches,
/// within `reach` of a sample, a branch cut there unless it runs straight
/// through, judged by the directions in and out, and the arcs either side
/// carried onto it along the straight way.
fn split_at_touches_as_traced(
    found: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<Traced> {
    let reach = options.chord.max(tol.confusion()) * 60.0;
    // How far a touch solved from a low-angle sample may stand from it and
    // still be the one it lies beside, and how far short of a touch any
    // stalled end is carried onto it: the walk gives out where the arms
    // close in on each other. A touch solved from a stalled end itself is
    // that end's at any distance: a stall cut back to the last sample
    // clear of it stands a whole stride short, and `carry_onto` follows
    // the section out to the touch or leaves the end where it is. Both ends
    // of a branch marked stalled count as stalled, though one may have
    // left the domain: `carry_onto` leaves such an end where it is unless
    // the section runs on from it, ahead, to the touch.
    let carry = reach * 40.0;

    let mut touches: Vec<Touch> = Vec::new();
    // The touch solved from each stalled branch's ends, head then tail.
    let mut solved_ends: Vec<[Option<usize>; 2]> = vec![[None; 2]; found.len()];
    for (branch, solved) in found.iter().zip(&mut solved_ends) {
        let n = branch.points.len();
        if n < 5 || !interior_is_transversal(branch, a, b, tol) {
            continue;
        }
        let sines: Vec<f64> = (0..n)
            .map(|i| crossing_sine(a, b, branch.on_a[i], branch.on_b[i], tol))
            .collect();
        for i in 0..n {
            let low = sines[i] < BRANCH_POINT_SINE
                && (i == 0 || sines[i] <= sines[i - 1])
                && (i + 1 == n || sines[i] <= sines[i + 1]);
            if !low
                || touches
                    .iter()
                    .any(|t| t.point.distance(branch.points[i]) <= reach)
            {
                continue;
            }
            let Some(touch) = touching_point(a, b, branch.on_a[i], branch.on_b[i], tol) else {
                continue;
            };
            let stall = branch.stopped == Stopped::Stalled && (i == 0 || i + 1 == n);
            if !stall && touch.point.distance(branch.points[i]) > carry {
                continue;
            }
            let at = touches
                .iter()
                .position(|t| t.point.distance(touch.point) <= reach)
                .unwrap_or_else(|| {
                    touches.push(touch);
                    touches.len() - 1
                });
            if stall {
                solved[usize::from(i != 0)] = Some(at);
            }
        }
    }
    if touches.is_empty() {
        return found;
    }
    let touch_near = |p: Point, within: f64| -> Option<usize> {
        touches.iter().position(|t| t.point.distance(p) <= within)
    };

    let mut out = Vec::with_capacity(found.len());
    for (branch, solved) in found.into_iter().zip(solved_ends) {
        let n = branch.points.len();
        // A sample on a touch, or either end of a step passing within
        // reach of one: the walk strides straight arms in long steps.
        let mut visits: Vec<Option<usize>> = branch
            .points
            .iter()
            .map(|p| touch_near(*p, reach))
            .collect();
        for i in 1..n {
            let (p, q) = (branch.points[i - 1], branch.points[i]);
            if let Some(t) = touches
                .iter()
                .position(|t| distance_to_segment(t.point, p, q) <= reach)
            {
                visits[i - 1].get_or_insert(t);
                visits[i].get_or_insert(t);
            }
        }
        // The touch a stalled end heads for: one within `carry` of it, or
        // failing that the one solved from the end itself. Where the end
        // will not carry onto the first, the one solved from it is tried
        // next: two crossings a little apart on one branch both lie within
        // `carry` of its end, and the one listed first may be the other.
        let solved_at = |at: usize| if at == 0 { solved[0] } else { solved[1] };
        let end_touch = |at: usize| -> Option<usize> {
            if branch.stopped != Stopped::Stalled {
                return None;
            }
            touch_near(branch.points[at], carry).or_else(|| solved_at(at))
        };
        let end_fallback = |at: usize| -> Option<usize> {
            let first = end_touch(at)?;
            solved_at(at).filter(|&t| t != first)
        };
        let stall_end = |at: usize| end_touch(at).is_some();
        let touched =
            visits.iter().any(Option::is_some) || (n > 0 && (stall_end(0) || stall_end(n - 1)));
        if !touched || !interior_is_transversal(&branch, a, b, tol) {
            out.push(branch);
            continue;
        }
        let closed = branch.closed();
        // A loop is read from the start of a pass over a touch, so each arc
        // of it is one run of samples clear of every touch.
        let order: Vec<usize> = if closed {
            let last = if branch.points[0].distance(branch.points[n - 1]) <= tol.confusion() {
                n - 1
            } else {
                n
            };
            let Some(first) =
                (0..last).find(|&i| visits[i].is_some() && visits[(i + last - 1) % last].is_none())
            else {
                out.push(branch);
                continue;
            };
            (0..last).map(|k| (first + k) % last).collect()
        } else {
            (0..n).collect()
        };
        let mut runs: Vec<(Vec<usize>, Option<usize>, Option<usize>)> = Vec::new();
        let mut run: Vec<usize> = Vec::new();
        let mut head: Option<usize> = None;
        for (k, &i) in order.iter().enumerate() {
            if let Some(t) = visits[i] {
                if !run.is_empty() {
                    runs.push((core::mem::take(&mut run), head, Some(t)));
                }
                head = Some(t);
                continue;
            }
            if run.is_empty() && k == 0 && !closed && stall_end(i) {
                head = end_touch(i);
            }
            run.push(i);
        }
        if !run.is_empty() {
            let last = *run.last().unwrap_or(&0);
            let tail = if closed {
                visits[order[0]]
            } else if last == n - 1 {
                end_touch(last)
            } else {
                None
            };
            runs.push((run, head, tail));
        }
        // A branch passing each touch once, running straight through,
        // stays whole: the crossing is the caller's to split like any other.
        let unit = |v: Vector| {
            let m = v.magnitude();
            (m > tol.confusion()).then(|| v / m)
        };
        let pass = |before: &[usize], after: &[usize], t: usize| -> (usize, bool) {
            let (Some(&i), Some(&o)) = (before.last(), after.first()) else {
                return (t, false);
            };
            let into = unit(touches[t].point - branch.points[i]);
            let onward = unit(branch.points[o] - touches[t].point);
            (
                t,
                matches!((into, onward), (Some(x), Some(y)) if x.dot(y) > 0.5),
            )
        };
        let mut passes: Vec<(usize, bool)> = runs
            .windows(2)
            .filter_map(|w| Some(pass(&w[0].0, &w[1].0, w[0].2?)))
            .collect();
        if closed
            && let (Some(last), Some(first)) = (runs.last(), runs.first())
            && let Some(t) = last.2
        {
            passes.push(pass(&last.0, &first.0, t));
        }
        let mut met: Vec<usize> = passes.iter().map(|(t, _)| *t).collect();
        met.sort_unstable();
        let repeated = met.windows(2).any(|w| w[0] == w[1]);
        let straight = passes.iter().all(|(_, s)| *s);
        let stopped_short = !closed && (stall_end(0) || stall_end(n - 1));
        if !repeated && straight && !stopped_short {
            out.push(branch);
            continue;
        }

        let mut arcs: Vec<Traced> = Vec::new();
        for (run, head, tail) in runs {
            let mut arc = Traced {
                points: run.iter().map(|&i| branch.points[i]).collect(),
                on_a: run.iter().map(|&i| branch.on_a[i]).collect(),
                on_b: run.iter().map(|&i| branch.on_b[i]).collect(),
                stopped: Stopped::Stalled,
            };
            if arc.points.len() < 2 {
                continue;
            }
            let fallback = |at_end: bool, index: usize| {
                (at_end && !closed && visits[index].is_none())
                    .then(|| end_fallback(index))
                    .flatten()
            };
            let ends = [
                (tail, false, fallback(run.last() == Some(&(n - 1)), n - 1)),
                (head, true, fallback(run.first() == Some(&0), 0)),
            ];
            for (end, at_head, next) in ends {
                if let Some(t) = end
                    && !carry_onto_straight(&mut arc, &touches[t], at_head, (a, b), options, tol)
                    && let Some(t) = next
                {
                    carry_onto_straight(&mut arc, &touches[t], at_head, (a, b), options, tol);
                }
            }
            let length: f64 = arc.points.windows(2).map(|w| w[0].distance(w[1])).sum();
            if length <= options.chord * 10.0 || arc.points.len() < 4 {
                continue;
            }
            let m = arc.points.len();
            if head.is_some() && head == tail && arc.points[0].distance(arc.points[m - 1]) <= reach
            {
                // The golden section of the samples, clear of the middle
                // and quarters where a symmetric loop's features fall.
                let middle = m * 382 / 1000;
                let part = |r: core::ops::RangeInclusive<usize>| Traced {
                    points: arc.points[r.clone()].to_vec(),
                    on_a: arc.on_a[r.clone()].to_vec(),
                    on_b: arc.on_b[r].to_vec(),
                    stopped: Stopped::Stalled,
                };
                arcs.push(part(0..=middle));
                arcs.push(part(middle..=m - 1));
            } else {
                arcs.push(arc);
            }
        }
        out.extend(arcs);
    }
    out
}

/// An arc's end carried onto a touch along the straight way there, as
/// [`split_at_touches_as_traced`] carries it: samples `reach` apart down to
/// the touch, then halving to an eighth of it, each corrected onto both
/// surfaces. Nothing is added where the end does not head for the touch,
/// a sample further than `reach` from it will not settle onto the section,
/// or the samples would take the arc past the points a branch may hold.
/// Whether the end was carried onto the touch.
fn carry_onto_straight(
    arc: &mut Traced,
    touch: &Touch,
    at_head: bool,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> bool {
    let reach = options.chord.max(tol.confusion()) * 60.0;
    let n = arc.points.len();
    let (end, inner) = if at_head { (0, 1) } else { (n - 1, n - 2) };
    let from = arc.points[end];
    let gap = touch.point - from;
    let distance = gap.magnitude();
    let mut added: Vec<Contact> = Vec::new();
    if distance > tol.confusion() {
        let along = gap * (1.0 / distance);
        let heading = from - arc.points[inner];
        if distance > reach && heading.dot(along) < 0.5 * heading.magnitude() {
            return false;
        }
        // How far short of the touch each sample stands: steps of at most
        // `reach` down to it, then halving to an eighth of it, where the
        // two arms still stand far enough apart to correct onto the right
        // one and the last straight stretch bows from the arm well inside
        // the chord.
        let pieces = (distance / reach).ceil().max(1.0);
        #[allow(clippy::cast_precision_loss)]
        if pieces + 3.0 > options.max_points.saturating_sub(n) as f64 {
            return false;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = pieces as usize;
        #[allow(clippy::cast_precision_loss)]
        let mut short: Vec<f64> = (1..count)
            .map(|k| distance - distance * k as f64 / pieces)
            .collect();
        let mut close = short.last().copied().unwrap_or(distance).min(reach);
        for _ in 0..3 {
            close *= 0.5;
            short.push(close);
        }
        let (mut on_a, mut on_b) = (arc.on_a[end], arc.on_b[end]);
        for remaining in short {
            let s = distance - remaining;
            let guess = from + along * s;
            let start = [on_a.0, on_a.1, on_b.0, on_b.1];
            let settled = correct(
                a,
                b,
                start,
                guess,
                Some((from, along, s)),
                options.chord,
                tol,
            )
            .filter(|c| c.point.distance(guess) <= remaining * 0.5);
            let Some(c) = settled else {
                // Beside the touch a sample that will not settle is left
                // out; further off, the arm is not there to follow.
                if remaining > reach {
                    return false;
                }
                continue;
            };
            on_a = beside(a, c.on_a, on_a);
            on_b = beside(b, c.on_b, on_b);
            added.push(Contact {
                on_a,
                on_b,
                point: c.point,
            });
        }
        added.push(Contact {
            on_a: beside(a, touch.on_a, on_a),
            on_b: beside(b, touch.on_b, on_b),
            point: touch.point,
        });
    } else {
        let (on_a, on_b) = (arc.on_a[end], arc.on_b[end]);
        arc.points[end] = touch.point;
        arc.on_a[end] = beside(a, touch.on_a, on_a);
        arc.on_b[end] = beside(b, touch.on_b, on_b);
        return true;
    }
    if at_head {
        added.reverse();
        arc.points.splice(0..0, added.iter().map(|c| c.point));
        arc.on_a.splice(0..0, added.iter().map(|c| c.on_a));
        arc.on_b.splice(0..0, added.iter().map(|c| c.on_b));
    } else {
        arc.points.extend(added.iter().map(|c| c.point));
        arc.on_a.extend(added.iter().map(|c| c.on_a));
        arc.on_b.extend(added.iter().map(|c| c.on_b));
    }
    true
}

/// Whether every sample of `arc` clear of the touches lies within `near` of
/// `other`, the middles of its chords within twice that: the samples of
/// both lie on the section, and the chords of a march stand within its
/// chord of it, so a stretch traced twice holds each of its samples within
/// the chord of the other tracing. Across the section that is widened to
/// where the confusion distance leaves it: see [`across_section`]. Round a
/// touch the samples are cut and carried, not marched, and are not
/// weighed.
fn lies_on(arc: &Traced, other: &Traced, near: f64, across: &[Across], touches: &[Touch]) -> bool {
    let mut last = 0;
    let mut weighed = false;
    // The samples, and the middle of each chord: an arc joining two touches
    // may have no sample clear of them.
    let middles = arc.points.windows(2).zip(across.windows(2)).map(|(w, r)| {
        let step = w[1] - w[0];
        let length = step.magnitude();
        let way = (length > f64::MIN_POSITIVE).then(|| step * (1.0 / length));
        (
            w[0] + step * 0.5,
            2.0 * near,
            (2.0 * r[0].0.min(r[1].0), way),
        )
    });
    let samples = arc.points.iter().zip(across).map(|(p, r)| (*p, near, *r));
    for (p, within, wide) in samples.chain(middles) {
        if touches.iter().any(|t| t.point.distance(p) <= t.radius) {
            continue;
        }
        if !holds(other, p, within, wide, &mut last) {
            return false;
        }
        weighed = true;
    }
    weighed
}

/// How far across the section two tracings of it may stand apart at each
/// sample of `arc`, and the way the section runs there: a sample is held
/// to both surfaces within the confusion distance `ε`, and a step `δ`
/// across the section opens a gap of `δ · sin θ`, `θ` the angle the
/// surfaces cross at, so each tracing lies within `ε / sin θ` of the
/// section across it, and two within twice that of each other. Round a
/// touch, where the surfaces cross at no angle, and where the angle is not
/// known, no further than `near`.
fn across_section(
    arc: &Traced,
    touches: &[Touch],
    near: f64,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    tol: Tolerances,
) -> Vec<Across> {
    (0..arc.points.len())
        .map(|k| {
            let p = arc.points[k];
            let none = (near, None);
            if touches.iter().any(|t| t.point.distance(p) <= t.radius) {
                return none;
            }
            let (Some(na), Some(nb)) = (
                normal_at(a, arc.on_a[k], tol),
                normal_at(b, arc.on_b[k], tol),
            ) else {
                return none;
            };
            let cross = na.cross(nb);
            let sine = cross.magnitude();
            if sine <= f64::MIN_POSITIVE {
                return none;
            }
            (
                near.max(2.0 * tol.confusion() / sine),
                Some(cross * (1.0 / sine)),
            )
        })
        .collect()
}

/// How far across the section a point may stand off another tracing of
/// it, and the way the section runs there: see [`across_section`].
type Across = (f64, Option<Vector>);

/// Whether two arcs are one stretch between the same two touches: each
/// runs from one touch to the other, leaving each along the same half arm.
fn same_stretch(arc: &Traced, other: &Traced, near: f64, touches: &[Touch]) -> bool {
    let ends = |t: &Traced| -> Option<[(usize, usize); 2]> {
        let m = t.points.len();
        if m < 2 {
            return None;
        }
        let at = |end: usize, beside: usize| {
            let k = touches
                .iter()
                .position(|u| u.point.distance(t.points[end]) <= near.max(u.blur))?;
            Some((k, half_arm(&touches[k], t.points[beside])?))
        };
        Some([at(0, 1)?, at(m - 1, m - 2)?])
    };
    match (ends(arc), ends(other)) {
        (Some(x), Some(y)) => x == y || x == [y[1], y[0]],
        _ => false,
    }
}

/// Whether a point lies within `near` of a polyline, or beside a segment
/// of it within `across` (see [`beside_segment`]), searched outward from
/// the segment `last` a point beside it was found on, which it is left at:
/// points taken in order along one curve are found at once.
fn holds(other: &Traced, p: Point, near: f64, across: Across, last: &mut usize) -> bool {
    let segments = other.points.len().saturating_sub(1);
    let on = |k: usize| beside_segment(p, other.points[k], other.points[k + 1], near, across);
    let from = (*last).min(segments.saturating_sub(1));
    let found = (0..segments).find_map(|step| {
        let ahead = from + step;
        if ahead < segments && on(ahead) {
            return Some(ahead);
        }
        let behind = from.checked_sub(step)?;
        on(behind).then_some(behind)
    });
    found.is_some_and(|k| {
        *last = k;
        true
    })
}

/// Whether a point lies within `near` of the segment `from`–`to`, or, with
/// the way `along` the section runs at the point, stands off a point of the
/// segment within `near` along the section and within `across` across it.
fn beside_segment(p: Point, from: Point, to: Point, near: f64, (across, along): Across) -> bool {
    if distance_to_segment(p, from, to) <= near {
        return true;
    }
    let Some(along) = along else {
        return false;
    };
    let step = to - from;
    let length = step.square_magnitude();
    if across <= near || length <= f64::MIN_POSITIVE {
        return false;
    }
    let t = (p - from).dot(step) / length;
    if !(0.0..=1.0).contains(&t) {
        return false;
    }
    let off = p - (from + step * t);
    let ahead = off.dot(along);
    ahead.abs() <= near && (off - along * ahead).magnitude() <= across
}

/// An arc with the samples at either end that lie on arcs already kept cut
/// back, so a stretch traced twice where two tracings overlap at their
/// ends is kept once. Each cut end keeps the last sample on the other
/// tracing, so the two still meet; an end on a touch is kept.
fn trimmed(
    arc: Traced,
    kept: &[Traced],
    near: f64,
    across: &[Across],
    touches: &[Touch],
) -> Option<Traced> {
    let n = arc.points.len();
    // How many samples from one end run along arcs already kept, and the
    // arc the last of them lies on: each sample is looked for first on the
    // arc the one before was found on.
    let run = |order: &mut dyn Iterator<Item = usize>| -> (usize, Option<usize>) {
        let mut holder: Option<(usize, usize)> = None;
        let mut count = 0;
        for k in order {
            let p = arc.points[k];
            if touches.iter().any(|t| t.point.distance(p) <= t.radius) {
                break;
            }
            let on_holder = holder.and_then(|(o, mut last)| {
                holds(&kept[o], p, near, across[k], &mut last).then_some((o, last))
            });
            let found = on_holder.or_else(|| {
                kept.iter().enumerate().find_map(|(o, other)| {
                    let mut last = 0;
                    holds(other, p, near, across[k], &mut last).then_some((o, last))
                })
            });
            if found.is_none() {
                break;
            }
            holder = found;
            count += 1;
        }
        (count, holder.map(|(o, _)| o))
    };
    let (head, head_on) = run(&mut (0..n));
    if head == n {
        return None;
    }
    let (back, tail_on) = run(&mut (0..n).rev());
    let tail = n - 1 - back;
    if tail < head {
        return None;
    }
    let mut points = arc.points[head..=tail].to_vec();
    let mut on_a = arc.on_a[head..=tail].to_vec();
    let mut on_b = arc.on_b[head..=tail].to_vec();
    // Each cut end meets the arc it ran along: at that arc's own end where
    // it stops between the last sample cut and the first kept, so the two
    // do not run side by side, and otherwise at the last sample cut.
    let meet = |cut: usize, kept_at: usize, on: Option<usize>| -> (Point, (f64, f64), (f64, f64)) {
        let (u, h) = (arc.points[kept_at], arc.points[cut]);
        on.map(|o| &kept[o])
            .and_then(|other| {
                let m = other.points.len();
                [0, m - 1]
                    .into_iter()
                    .filter(|&e| {
                        let wide = (across[cut].0.min(across[kept_at].0), across[cut].1);
                        beside_segment(other.points[e], u, h, near, wide)
                    })
                    .min_by(|&x, &y| {
                        other.points[x]
                            .distance(u)
                            .total_cmp(&other.points[y].distance(u))
                    })
                    .map(|e| (other.points[e], other.on_a[e], other.on_b[e]))
            })
            .unwrap_or((h, arc.on_a[cut], arc.on_b[cut]))
    };
    if head > 0 {
        let (p, a, b) = meet(head - 1, head, head_on);
        points.insert(0, p);
        on_a.insert(0, a);
        on_b.insert(0, b);
    }
    if tail + 1 < n {
        let (p, a, b) = meet(tail + 1, tail, tail_on);
        points.push(p);
        on_a.push(a);
        on_b.push(b);
    }
    (points.len() >= 2).then_some(Traced {
        points,
        on_a,
        on_b,
        stopped: arc.stopped,
    })
}

/// Whether the chord of the step into sample `i` stands off the section by
/// more than `near`, the march's chord: the section solved on the plane
/// square to the step through its middle lands further than that from it,
/// so the step is not one along a single curve.
fn off_the_section(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    branch: &Traced,
    i: usize,
    near: f64,
    tol: Tolerances,
) -> bool {
    let (p, q) = (branch.points[i - 1], branch.points[i]);
    let step = q - p;
    let length = step.magnitude();
    if length <= near {
        return false;
    }
    let along = step / length;
    let middle = p + step * 0.5;
    let start = [
        branch.on_a[i - 1].0,
        branch.on_a[i - 1].1,
        branch.on_b[i - 1].0,
        branch.on_b[i - 1].1,
    ];
    correct(
        a,
        b,
        start,
        middle,
        Some((p, along, length * 0.5)),
        near,
        tol,
    )
    .is_some_and(|landed| landed.point.distance(middle) > near)
}

/// The four halves of a touch's two arms, each pointing away from it.
fn halves(arms: [Vector; 2]) -> [Vector; 4] {
    [arms[0], -arms[0], arms[1], -arms[1]]
}

/// The half arm of a touch a point lies nearest in direction.
fn half_arm(touch: &Touch, p: Point) -> Option<usize> {
    let arms = touch.arms?;
    let off = p - touch.point;
    let m = off.magnitude();
    if m <= f64::MIN_POSITIVE {
        return None;
    }
    halves(arms)
        .iter()
        .map(|h| h.dot(off) / m)
        .enumerate()
        .max_by(|x, y| x.1.total_cmp(&y.1))
        .map(|(k, _)| k)
}

/// Whether a sample of the section lies on one half arm of a touch: the
/// section solved on the plane square to that half through the sample
/// lands within `near` of it, or, with `blurred`, within the band the
/// confusion distance leaves round the arm that far out.
fn on_arm(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    (half, blurred): (usize, bool),
    p: Point,
    near: f64,
    tol: Tolerances,
) -> bool {
    let Some(arms) = touch.arms else {
        return false;
    };
    let direction = halves(arms)[half];
    let along = (p - touch.point).dot(direction);
    let band = if blurred && touch.spread > 0.0 {
        tol.confusion() / (along * touch.spread)
    } else {
        0.0
    };
    along > 0.0
        && along_arm(a, b, touch, direction, (along, arms), near, tol)
            .is_some_and(|landed| landed.point.distance(p) <= near.max(band))
}

/// Each touch's radius, arms and snap, once all are known: a sample counts
/// as on the touch within the reach of [`split_at_touches`], and a stalled
/// end or a step between its arms as heading for it within the `carry` of
/// [`find_touches`], each held to half the way to the nearest other touch
/// where that is less, so nothing counts as on or heading for two; the
/// arms are kept where they stand apart at the radius.
fn settle_touches(
    touches: &mut [Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    let points: Vec<Point> = touches.iter().map(|t| t.point).collect();
    for (i, touch) in touches.iter_mut().enumerate() {
        let nearest = points
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, p)| p.distance(touch.point))
            .fold(f64::INFINITY, f64::min);
        let form = touch_arms(a, b, touch.on_a, touch.on_b, tol);
        let spread = form.map_or(0.0, |(_, spread, _, _)| spread);
        touch.snap = form.map_or(0.0, |(_, _, bend, _)| {
            if bend > 0.0 {
                (2.0 * tol.confusion() / bend).sqrt()
            } else {
                0.0
            }
        });
        touch.radius = reach.min(nearest * 0.5);
        touch.around = (reach * 40.0).min(nearest * 0.5);
        touch.spread = spread;
        touch.arms = form
            .map(|(arms, _, _, _)| arms)
            .filter(|&arms| apart(arms, touch.radius, near));
    }
}

/// The length of a traced polyline.
fn arc_length(arc: &Traced) -> f64 {
    arc.points.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// The section a distance `reach` out from a touch along one half of an
/// arm: corrected onto both surfaces on the plane square to the arm there,
/// from parameters moved out along it on each chart. `None` where the
/// correction does not settle nearer that arm than the other.
fn along_arm(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    touch: &Touch,
    half: Vector,
    (reach, arms): (f64, [Vector; 2]),
    within: f64,
    tol: Tolerances,
) -> Option<Contact> {
    let moved = |surface: &SurfaceGeometry, (u, v): (f64, f64)| -> Option<(f64, f64)> {
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = e.mul_add(g, -(f * f));
        if det <= f64::MIN_POSITIVE {
            return None;
        }
        let (p, q) = (su.dot(half) * reach, sv.dot(half) * reach);
        Some((
            u + g.mul_add(p, -(f * q)) / det,
            v + e.mul_add(q, -(f * p)) / det,
        ))
    };
    let (ua, va) = moved(a, touch.on_a)?;
    let (ub, vb) = moved(b, touch.on_b)?;
    let guess = touch.point + half * reach;
    let landed = correct(
        a,
        b,
        [ua, va, ub, vb],
        guess,
        Some((touch.point, half, reach)),
        within,
        tol,
    )?;
    // The other arm crosses that plane `reach · tan α` off this one.
    let apart = arms[0].cross(arms[1]).magnitude();
    (landed.point.distance(guess) * 2.0 < reach * apart).then_some(landed)
}

/// The touches a branch's low-angle samples are solved onto, added to
/// `touches` unless one is already there, and the touch solved from each of
/// its ends where it stalled, head then tail.
fn find_touches(
    branch: &Traced,
    touches: &mut Vec<Touch>,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> [Option<usize>; 2] {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    // How far a touch solved from a low-angle sample may stand from it and
    // still be the one it lies beside, and how far short of a touch any
    // stalled end is carried onto it: the walk gives out where the arms
    // close in on each other. A touch solved from a stalled end itself is
    // that end's at any distance: a stall cut back to the last sample
    // clear of it stands a whole stride short, and `carry_onto` follows
    // the section out to the touch or leaves the end where it is. Both ends
    // of a branch marked stalled count as stalled, though one may have
    // left the domain: `carry_onto` leaves such an end where it is unless
    // the section runs on from it, ahead, to the touch.
    let carry = reach * 40.0;
    let mut solved = [None; 2];
    let n = branch.points.len();
    if n < 5 {
        return solved;
    }
    // A branch with no transversal stretch to certify it may be tangential
    // contact rather than a curve, and only crossings are taken from it.
    let transversal = interior_is_transversal(branch, a, b, tol);
    let sines: Vec<f64> = (0..n)
        .map(|i| crossing_sine(a, b, branch.on_a[i], branch.on_b[i], tol))
        .collect();
    for i in 0..n {
        let low = sines[i] < BRANCH_POINT_SINE
            && (i == 0 || sines[i] <= sines[i - 1])
            && (i + 1 == n || sines[i] <= sines[i + 1]);
        if !low
            || touches
                .iter()
                .any(|t| t.point.distance(branch.points[i]) <= near)
        {
            continue;
        }
        let Some(mut touch) = touching_point(a, b, branch.on_a[i], branch.on_b[i], tol) else {
            continue;
        };
        // A touch whose section passes it without crossing there is none.
        let form = touch_arms(a, b, touch.on_a, touch.on_b, tol);
        if form.is_some_and(|form| !crosses(&touch, form, near, tol)) {
            continue;
        }
        touch.blur = form.map_or(0.0, |form| blur_of(form.3, tol));
        touch.arms = form
            .map(|(arms, _, _, _)| arms)
            .filter(|&arms| apart(arms, reach, near));
        if !transversal && touch.arms.is_none() {
            continue;
        }
        let stall = branch.stopped == Stopped::Stalled && (i == 0 || i + 1 == n);
        if !stall && touch.point.distance(branch.points[i]) > carry {
            continue;
        }
        // Two solves land on one touch to within the chord the march
        // resolves, or where either may land; touches further apart are
        // distinct, however close.
        let at = touches
            .iter()
            .position(|t| same_touch(t, &touch, near))
            .unwrap_or_else(|| {
                touches.push(touch);
                touches.len() - 1
            });
        if stall {
            solved[usize::from(i != 0)] = Some(at);
        }
    }
    solved
}

/// A branch as [`cut_at_touches`] leaves it.
struct Cut {
    /// The branch whole, or the arcs cut from it.
    pieces: Vec<Traced>,
    /// Whether it was cut.
    cut: bool,
    /// The touches it runs into or through and the half arm it runs along
    /// at each, by index into the touches and into [`halves`].
    covered: Vec<(usize, usize)>,
}

/// A branch with each touch a step of it passes over, both ends of the
/// step clear of the touch's radius, put in as a sample where the step
/// passes it; with each sample, the touch it was put in for.
///
/// A step passes over a touch it passes within the touch's radius, or
/// within its `snap`, held to its `around`: the straight way from the
/// step to a touch that near lies on both surfaces within the confusion
/// distance, so the walk, which places a sample across the section only
/// to within that distance over the sine of the crossing angle, may stride
/// past the touch that far off it.
fn with_touches_on_steps(
    branch: Traced,
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    counts: &dyn Fn(&Touch) -> bool,
) -> (Traced, Vec<Option<usize>>) {
    let n = branch.points.len();
    let mut out = Traced {
        points: Vec::with_capacity(n),
        on_a: Vec::with_capacity(n),
        on_b: Vec::with_capacity(n),
        stopped: branch.stopped,
    };
    let mut on_touch = Vec::with_capacity(n);
    for i in 0..n {
        if i > 0 {
            let (p, q) = (branch.points[i - 1], branch.points[i]);
            let step = q - p;
            let length2 = step.square_magnitude();
            let mut passed: Vec<(f64, usize)> = touches
                .iter()
                .enumerate()
                .filter(|(_, t)| {
                    let over = (t.point - p).dot(step);
                    counts(t)
                        && t.point.distance(p) > t.radius
                        && t.point.distance(q) > t.radius
                        && over > 0.0
                        && over < length2
                        && distance_to_segment(t.point, p, q) <= t.radius.max(t.snap.min(t.around))
                })
                .map(|(k, t)| ((t.point - p).dot(step) / length2.max(f64::MIN_POSITIVE), k))
                .collect();
            passed.sort_by(|x, y| x.0.total_cmp(&y.0));
            for (_, k) in passed {
                let t = &touches[k];
                out.points.push(t.point);
                out.on_a.push(beside(a, t.on_a, branch.on_a[i - 1]));
                out.on_b.push(beside(b, t.on_b, branch.on_b[i - 1]));
                on_touch.push(Some(k));
            }
        }
        out.points.push(branch.points[i]);
        out.on_a.push(branch.on_a[i]);
        out.on_b.push(branch.on_b[i]);
        on_touch.push(None);
    }
    (out, on_touch)
}

/// One branch cut at the touches it visits, as [`split_at_touches`] cuts
/// it. `solved` holds the touches solved from its ends, head then tail.
fn cut_at_touches(
    branch: Traced,
    solved: [Option<usize>; 2],
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> Cut {
    let near = options.chord.max(tol.confusion());
    // A branch with no transversal stretch to certify it is cut only at
    // crossings, as only crossings are found from it.
    let transversal = interior_is_transversal(&branch, a, b, tol);
    let counts = |t: &Touch| transversal || t.arms.is_some();
    let touch_near = |p: Point| -> Option<usize> {
        (0..touches.len())
            .filter(|&t| counts(&touches[t]) && touches[t].point.distance(p) <= touches[t].around)
            .min_by(|&x, &y| {
                touches[x]
                    .point
                    .distance(p)
                    .total_cmp(&touches[y].point.distance(p))
            })
    };
    // A step passing over touches, its ends clear of them, gets each as a
    // sample where it passes: a straight arm walked in long strides passes
    // several touches a step apart, and each is cut at on its own.
    let (branch, on_touch) = with_touches_on_steps(branch, touches, (a, b), &counts);
    let n = branch.points.len();
    // A sample on a touch, or either end of a step passing within its
    // radius: the walk strides straight arms in long steps.
    let mut visits: Vec<Option<usize>> = branch
        .points
        .iter()
        .map(|p| {
            touches
                .iter()
                .position(|t| counts(t) && t.point.distance(*p) <= t.radius)
        })
        .collect();
    for i in 1..n {
        let (p, q) = (branch.points[i - 1], branch.points[i]);
        if let Some(t) = touches
            .iter()
            .position(|t| counts(t) && distance_to_segment(t.point, p, q) <= t.radius)
        {
            visits[i - 1].get_or_insert(t);
            visits[i].get_or_insert(t);
        }
    }
    // A step from one arm of a touch onto another passes over the touch
    // however far off it the step runs: a walk turning a corner there cuts
    // it by as much as its step. One stepping from an arm of one touch onto
    // an arm of another passes over both where its chord leaves the
    // section, as a step along one arm through the other touch does not.
    // Both ends must lie on the arms, so a stretch of the section near a
    // touch but on neither is not taken for one.
    let arm_of = |t: usize, p: Point| -> Option<usize> {
        let touch = &touches[t];
        if !counts(touch) || touch.point.distance(p) > touch.around {
            return None;
        }
        let h = half_arm(touch, p)?;
        on_arm((a, b), touch, (h, true), p, near, tol).then_some(h)
    };
    for i in 1..n {
        if visits[i - 1].is_some() || visits[i].is_some() {
            continue;
        }
        let (p, q) = (branch.points[i - 1], branch.points[i]);
        let mut ends: [Option<(usize, usize)>; 2] = [None; 2];
        let mut over = None;
        for t in 0..touches.len() {
            let (hp, hq) = (arm_of(t, p), arm_of(t, q));
            if let (Some(hp), Some(hq)) = (hp, hq) {
                if hp != hq {
                    over = Some((t, t));
                }
                ends = [Some((t, hp)), Some((t, hq))];
                break;
            }
            ends[0] = ends[0].or(hp.map(|h| (t, h)));
            ends[1] = ends[1].or(hq.map(|h| (t, h)));
        }
        if over.is_none()
            && let [Some((tp, _)), Some((tq, _))] = ends
            && tp != tq
            && off_the_section((a, b), &branch, i, near, tol)
        {
            over = Some((tp, tq));
        }
        if let Some((tp, tq)) = over {
            visits[i - 1] = Some(tp);
            visits[i] = Some(tq);
        }
    }
    // The touch a stalled end heads for: one within `carry` of it, or
    // failing that the one solved from the end itself. Where the end
    // will not carry onto the first, the one solved from it is tried
    // next: two crossings a little apart on one branch both lie within
    // `carry` of its end, and the one listed first may be the other.
    let solved_at = |at: usize| if at == 0 { solved[0] } else { solved[1] };
    let end_touch = |at: usize| -> Option<usize> {
        if branch.stopped != Stopped::Stalled {
            return None;
        }
        touch_near(branch.points[at]).or_else(|| solved_at(at))
    };
    let end_fallback = |at: usize| -> Option<usize> {
        let first = end_touch(at)?;
        solved_at(at).filter(|&t| t != first)
    };
    let stall_end = |at: usize| end_touch(at).is_some();
    let touched =
        visits.iter().any(Option::is_some) || (n > 0 && (stall_end(0) || stall_end(n - 1)));
    if !touched {
        return Cut {
            pieces: vec![branch],
            cut: false,
            covered: Vec::new(),
        };
    }
    let closed = branch.closed();
    // The half arm a branch leaves each touch along: the sample beside each
    // run of samples on it, and an end carried onto it.
    let mut covered: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        let Some(t) = visits[i] else { continue };
        let beside = [
            i.checked_sub(1).or_else(|| closed.then_some(n - 1)),
            (i + 1 < n).then_some(i + 1).or_else(|| closed.then_some(0)),
        ];
        for j in beside.into_iter().flatten() {
            // The sample beside, clear of the touch or on another one a
            // step put beside it.
            if visits[j].is_none_or(|u| u != t && on_touch[j].is_some())
                && let Some(h) = half_arm(&touches[t], branch.points[j])
            {
                covered.push((t, h));
            }
        }
    }
    for at in [0, n - 1] {
        if !closed
            && visits[at].is_none()
            && let Some(t) = end_touch(at)
            && let Some(h) = half_arm(&touches[t], branch.points[at])
        {
            covered.push((t, h));
        }
    }
    // A loop is read from the start of a pass over a touch, so each arc
    // of it is one run of samples clear of every touch.
    let order: Vec<usize> = if closed {
        let last = if branch.points[0].distance(branch.points[n - 1]) <= tol.confusion() {
            n - 1
        } else {
            n
        };
        let Some(first) =
            (0..last).find(|&i| visits[i].is_some() && visits[(i + last - 1) % last].is_none())
        else {
            return Cut {
                pieces: vec![branch],
                cut: false,
                covered,
            };
        };
        (0..last).map(|k| (first + k) % last).collect()
    } else {
        (0..n).collect()
    };
    let mut runs: Vec<(Vec<usize>, Option<usize>, Option<usize>)> = Vec::new();
    let mut run: Vec<usize> = Vec::new();
    let mut head: Option<usize> = None;
    let mut before: Option<usize> = None;
    for (k, &i) in order.iter().enumerate() {
        if let Some(t) = visits[i] {
            if !run.is_empty() {
                runs.push((core::mem::take(&mut run), head, Some(t)));
            } else if let Some(j) = before
                && let (Some(u), Some(v)) = (on_touch[j], on_touch[i])
                && u != v
            {
                // Two touches a step put side by side: the stretch between
                // them is that much of the step.
                runs.push((vec![j, i], Some(u), Some(v)));
            }
            head = Some(t);
            before = Some(i);
            continue;
        }
        before = Some(i);
        if run.is_empty() && k == 0 && !closed && stall_end(i) {
            head = end_touch(i);
        }
        run.push(i);
    }
    if !run.is_empty() {
        let last = *run.last().unwrap_or(&0);
        let tail = if closed {
            visits[order[0]]
        } else if last == n - 1 {
            end_touch(last)
        } else {
            None
        };
        runs.push((run, head, tail));
    }
    // A branch passing each touch once, running straight through,
    // stays whole: the crossing is the caller's to split like any other.
    let unit = |v: Vector| {
        let m = v.magnitude();
        (m > tol.confusion()).then(|| v / m)
    };
    // Where the touch's two arms are known, the branch runs straight
    // through when it comes in along one arm and goes on along it, the
    // same way; it turns a corner when it goes on along the other. Which
    // arm a sample is on is read from the section's own direction there,
    // the cross of the normals, not the chord to the touch, which a long
    // stride tilts.
    let pass = |before: &[usize], after: &[usize], t: usize| -> (usize, bool) {
        // The samples either side clear of the touch: a stretch between two
        // touches starts on one.
        let clear = |&&k: &&usize| branch.points[k].distance(touches[t].point) > touches[t].radius;
        let (Some(&i), Some(&o)) = (before.iter().rev().find(clear), after.iter().find(clear))
        else {
            return (t, false);
        };
        let into = unit(touches[t].point - branch.points[i]);
        let onward = unit(branch.points[o] - touches[t].point);
        let (Some(into), Some(onward)) = (into, onward) else {
            return (t, false);
        };
        let straight = match touches[t].arms {
            Some(arms) => {
                let heading = |k: usize, chord: Vector| {
                    let at = Contact {
                        on_a: branch.on_a[k],
                        on_b: branch.on_b[k],
                        point: branch.points[k],
                    };
                    tangent_at(a, b, at, tol).unwrap_or(chord)
                };
                let line = |w: Vector| w.dot(arms[1]).abs() > w.dot(arms[0]).abs();
                line(heading(i, into)) == line(heading(o, onward)) && into.dot(onward) > 0.0
            }
            None => into.dot(onward) > 0.5,
        };
        (t, straight)
    };
    // Runs either side of a step over two touches meet neither straight.
    let mut passes: Vec<(usize, bool)> = runs
        .windows(2)
        .filter_map(|w| {
            let t = w[0].2?;
            Some(if w[1].1 == Some(t) {
                pass(&w[0].0, &w[1].0, t)
            } else {
                (t, false)
            })
        })
        .collect();
    if closed
        && let (Some(last), Some(first)) = (runs.last(), runs.first())
        && let Some(t) = last.2
    {
        passes.push(pass(&last.0, &first.0, t));
    }
    let mut met: Vec<usize> = passes.iter().map(|(t, _)| *t).collect();
    met.sort_unstable();
    let repeated = met.windows(2).any(|w| w[0] == w[1]);
    let straight = passes.iter().all(|(_, s)| *s);
    // An open branch ending on a touch or short of one is cut there too:
    // the samples it wanders round the touch on go, and its end is carried
    // onto the touch.
    let ends_at_touch = !closed
        && n > 0
        && (stall_end(0) || stall_end(n - 1) || visits[0].is_some() || visits[n - 1].is_some());
    if !repeated && straight && !ends_at_touch {
        return Cut {
            pieces: vec![branch],
            cut: false,
            covered,
        };
    }

    let mut arcs: Vec<Traced> = Vec::new();
    for (run, head, tail) in runs {
        let mut arc = Traced {
            points: run.iter().map(|&i| branch.points[i]).collect(),
            on_a: run.iter().map(|&i| branch.on_a[i]).collect(),
            on_b: run.iter().map(|&i| branch.on_b[i]).collect(),
            stopped: Stopped::Stalled,
        };
        // A run between two touches is the stretch between two crossings,
        // whatever few samples a long stride left on it: one sample is
        // joined to both touches along the section.
        let between = match (head, tail) {
            (Some(h), Some(t)) if h != t => Some((h, t)),
            _ => None,
        };
        if arc.points.len() == 1
            && let Some((h, t)) = between
        {
            let touch_contact = |touch: &Touch, near_to: Contact| Contact {
                on_a: beside(a, touch.on_a, near_to.on_a),
                on_b: beside(b, touch.on_b, near_to.on_b),
                point: touch.point,
            };
            let sample = Contact {
                on_a: arc.on_a[0],
                on_b: arc.on_b[0],
                point: arc.points[0],
            };
            let (from, to) = (
                touch_contact(&touches[h], sample),
                touch_contact(&touches[t], sample),
            );
            let mut budget = options.max_points;
            let joined = fill((a, b), from, sample, near, &mut budget, tol).zip(fill(
                (a, b),
                sample,
                to,
                near,
                &mut budget,
                tol,
            ));
            if let Some((before, after)) = joined {
                let all: Vec<Contact> = core::iter::once(from)
                    .chain(before)
                    .chain([sample])
                    .chain(after)
                    .chain([to])
                    .collect();
                arcs.push(Traced {
                    points: all.iter().map(|c| c.point).collect(),
                    on_a: all.iter().map(|c| c.on_a).collect(),
                    on_b: all.iter().map(|c| c.on_b).collect(),
                    stopped: Stopped::Stalled,
                });
            }
            continue;
        }
        if arc.points.len() < 2 {
            continue;
        }
        let fallback = |at_end: bool, index: usize| {
            (at_end && !closed && visits[index].is_none())
                .then(|| end_fallback(index))
                .flatten()
        };
        let ends = [
            (tail, false, fallback(run.last() == Some(&(n - 1)), n - 1)),
            (head, true, fallback(run.first() == Some(&0), 0)),
        ];
        for (end, at_head, next) in ends {
            let Some(t) = end else { continue };
            // Samples by the touch that stand off its arms by more than the
            // chord are the walk's wandering round it, where the surfaces
            // stand within the confusion distance of each other across the
            // corner, and go before the end is carried on.
            let mut trial = arc.clone();
            if touches[t].arms.is_some() {
                while trial.points.len() > 2 {
                    let at = if at_head { 0 } else { trial.points.len() - 1 };
                    let p = trial.points[at];
                    let wanders = p.distance(touches[t].point) <= touches[t].around
                        && half_arm(&touches[t], p).is_some_and(|h| {
                            !on_arm((a, b), &touches[t], (h, false), p, near, tol)
                        });
                    if !wanders {
                        break;
                    }
                    trial.points.remove(at);
                    trial.on_a.remove(at);
                    trial.on_b.remove(at);
                }
            }
            let carried = carry_onto(&mut trial, &touches[t], at_head, (a, b), options, tol)
                || next.is_some_and(|t| {
                    carry_onto(&mut trial, &touches[t], at_head, (a, b), options, tol)
                });
            if carried {
                arc = trial;
                continue;
            }
            // An end that will not carry onto the touch keeps the samples
            // the walk took past it, as the walk left them.
            let at = if at_head { run.first() } else { run.last() };
            let Some(&at) = at else { continue };
            let Some(k) = order.iter().position(|&i| i == at) else {
                continue;
            };
            let m = order.len();
            let step = |k: usize| -> Option<usize> {
                if at_head {
                    k.checked_sub(1).or_else(|| closed.then_some(m - 1))
                } else {
                    (k + 1 < m).then_some(k + 1).or_else(|| closed.then_some(0))
                }
            };
            let mut j = step(k);
            let mut seen = 0;
            while let Some(jj) = j
                && visits[order[jj]].is_some()
                && seen < m
            {
                let i = order[jj];
                if at_head {
                    arc.points.insert(0, branch.points[i]);
                    arc.on_a.insert(0, branch.on_a[i]);
                    arc.on_b.insert(0, branch.on_b[i]);
                } else {
                    arc.points.push(branch.points[i]);
                    arc.on_a.push(branch.on_a[i]);
                    arc.on_b.push(branch.on_b[i]);
                }
                j = step(jj);
                seen += 1;
            }
        }
        let length: f64 = arc.points.windows(2).map(|w| w[0].distance(w[1])).sum();
        let joins = between.is_some_and(|(h, t)| {
            let m = arc.points.len();
            arc.points[0].distance(touches[h].point) <= tol.confusion()
                && arc.points[m - 1].distance(touches[t].point) <= tol.confusion()
        });
        if !joins
            && (length <= options.chord * 10.0 || (arc.points.len() < 4 && length <= near * 60.0))
        {
            continue;
        }
        let m = arc.points.len();
        if let (Some(h), Some(t)) = (head, tail)
            && h == t
            && arc.points[0].distance(arc.points[m - 1]) <= touches[h].radius
        {
            // The golden section of the samples, clear of the middle
            // and quarters where a symmetric loop's features fall.
            let middle = m * 382 / 1000;
            let part = |r: core::ops::RangeInclusive<usize>| Traced {
                points: arc.points[r.clone()].to_vec(),
                on_a: arc.on_a[r.clone()].to_vec(),
                on_b: arc.on_b[r].to_vec(),
                stopped: Stopped::Stalled,
            };
            arcs.push(part(0..=middle));
            arcs.push(part(middle..=m - 1));
        } else {
            arcs.push(arc);
        }
    }
    // The arcs cut from it run along the arms they leave the touches by,
    // where they reach the touches.
    let mut covered = Vec::new();
    for arc in &arcs {
        let m = arc.points.len();
        for (end, beside) in [(0, 1), (m - 1, m.saturating_sub(2))] {
            if m < 2 {
                continue;
            }
            if let Some(t) = touches
                .iter()
                .position(|t| t.point.distance(arc.points[end]) <= tol.confusion())
                && let Some(h) = half_arm(&touches[t], arc.points[beside])
            {
                covered.push((t, h));
            }
        }
    }
    Cut {
        pieces: arcs,
        cut: true,
        covered,
    }
}

/// Samples of the section between two of its points, each chord between
/// neighbours standing within `near` of it: the section solved on the plane
/// square to a chord through its middle, the chord halved where that lands
/// further off. `None` where a solve does not settle within half the chord
/// of its middle, or more samples than `budget` would be needed.
fn fill(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    p: Contact,
    q: Contact,
    near: f64,
    budget: &mut usize,
    tol: Tolerances,
) -> Option<Vec<Contact>> {
    let chord = q.point - p.point;
    let length = chord.magnitude();
    if length <= near {
        return Some(Vec::new());
    }
    let middle = p.point + chord * 0.5;
    // The solve starts from the points of both surfaces nearest the middle
    // of the chord, so it lands on the stretch of section beside it: from
    // an end standing on a touch, where the surfaces are tangent and the
    // solve is singular, it may land on another stretch crossing there.
    let feet = nearest_on(a, p.on_a, middle, tol).zip(nearest_on(b, p.on_b, middle, tol));
    let start = feet.map_or(
        [p.on_a.0, p.on_a.1, p.on_b.0, p.on_b.1],
        |((ua, _), (ub, _))| [ua.0, ua.1, ub.0, ub.1],
    );
    // Where the solve gives out, as it does where the surfaces part too
    // slowly to tell apart, the middle of the chord lying within the
    // confusion distance of both surfaces is on the section as far as
    // anything can be.
    let on_both = || -> Option<Contact> {
        let ((on_a, pa), (on_b, pb)) = feet?;
        (pa.distance(middle) <= tol.confusion() && pb.distance(middle) <= tol.confusion())
            .then_some(Contact {
                on_a,
                on_b,
                point: middle,
            })
    };
    let landed = correct(
        a,
        b,
        start,
        middle,
        Some((p.point, chord / length, length * 0.5)),
        near,
        tol,
    )
    .filter(|c| c.point.distance(middle) <= length * 0.5)
    .or_else(on_both)?;
    let off = landed.point.distance(middle);
    if off > length * 0.5 {
        return None;
    }
    // A chord whose middle lies within the confusion distance of both
    // surfaces is on the section to tolerance, though the section solved
    // through it lands further off: an end settled onto a touch stands off
    // the exact crossing by up to the touch's own tolerance.
    if off <= near || on_both().is_some() {
        return Some(Vec::new());
    }
    *budget = budget.checked_sub(1)?;
    let m = Contact {
        on_a: beside(a, landed.on_a, p.on_a),
        on_b: beside(b, landed.on_b, p.on_b),
        point: landed.point,
    };
    let mut out = fill((a, b), p, m, near, budget, tol)?;
    out.push(m);
    out.extend(fill((a, b), m, q, near, budget, tol)?);
    Some(out)
}

/// Whether a sample beside a touch runs along the arm it lies out on: the
/// section's direction there, the cross of the normals, nearer that arm's
/// line than the other's. A sample of a walk wandering across the corner
/// does not.
fn along_its_arm(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    at: Contact,
    tol: Tolerances,
) -> bool {
    let (Some(arms), Some(h)) = (touch.arms, half_arm(touch, at.point)) else {
        return false;
    };
    let Some(tangent) = tangent_at(a, b, at, tol) else {
        return false;
    };
    let (own, other) = (arms[h / 2], arms[1 - h / 2]);
    tangent.dot(own).abs() >= tangent.dot(other).abs()
}

/// An arc's end carried onto a touch: samples along the straight way there
/// corrected onto both surfaces, then the touch itself. The arc stays as it
/// was where the end does not head for the touch, a sample further than
/// `reach` from it will not settle onto the section, or the samples would
/// take the arc past the points a branch may hold. Whether the end was
/// carried onto the touch.
fn carry_onto(
    arc: &mut Traced,
    touch: &Touch,
    at_head: bool,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> bool {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    let carry = reach * 40.0;
    let n = arc.points.len();
    let (end, inner) = if at_head { (0, 1) } else { (n - 1, n - 2) };
    let from = arc.points[end];
    let gap = touch.point - from;
    let distance = gap.magnitude();
    let mut added: Vec<Contact> = Vec::new();
    let end_contact = Contact {
        on_a: arc.on_a[end],
        on_b: arc.on_b[end],
        point: from,
    };
    if distance > tol.confusion()
        && distance <= touch.snap
        && along_its_arm((a, b), touch, end_contact, tol)
    {
        // The straight way to the touch lies on both surfaces to within
        // the confusion distance.
        let (on_a, on_b) = (arc.on_a[end], arc.on_b[end]);
        added.push(Contact {
            on_a: beside(a, touch.on_a, on_a),
            on_b: beside(b, touch.on_b, on_b),
            point: touch.point,
        });
    } else if distance > tol.confusion() {
        let along = gap * (1.0 / distance);
        let heading = from - arc.points[inner];
        if distance > reach && heading.dot(along) < 0.5 * heading.magnitude() {
            return false;
        }
        // How far short of the touch each sample stands: steps of at most
        // `reach` down to it, then halving to an eighth of it, where the
        // two arms still stand far enough apart to correct onto the right
        // one and the last straight stretch bows from the arm well inside
        // the chord. An end further off than `carry` is followed in to
        // `carry` short of the touch by samples placed where the chord
        // calls for them.
        let span = distance.min(carry);
        let pieces = (span / reach).ceil().max(1.0);
        #[allow(clippy::cast_precision_loss)]
        if pieces + 3.0 > options.max_points.saturating_sub(n) as f64 {
            return false;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = pieces as usize;
        #[allow(clippy::cast_precision_loss)]
        let mut short: Vec<f64> = (1..count)
            .map(|k| span - span * k as f64 / pieces)
            .collect();
        if distance > carry {
            short.insert(0, carry);
        }
        let mut close = short.last().copied().unwrap_or(distance).min(reach);
        for _ in 0..3 {
            close *= 0.5;
            short.push(close);
        }
        let (mut on_a, mut on_b) = (arc.on_a[end], arc.on_b[end]);
        for remaining in short {
            let s = distance - remaining;
            let guess = from + along * s;
            let start = [on_a.0, on_a.1, on_b.0, on_b.1];
            let settled = correct(
                a,
                b,
                start,
                guess,
                Some((from, along, s)),
                options.chord,
                tol,
            )
            .filter(|c| c.point.distance(guess) <= remaining * 0.5);
            let Some(c) = settled else {
                // Beside the touch a sample that will not settle is left
                // out; further off, the arm is not there to follow.
                if remaining > reach {
                    return false;
                }
                continue;
            };
            on_a = beside(a, c.on_a, on_a);
            on_b = beside(b, c.on_b, on_b);
            let c = Contact {
                on_a,
                on_b,
                point: c.point,
            };
            if added.is_empty() && remaining >= carry && distance > carry {
                let mut budget = options.max_points.saturating_sub(n);
                let Some(between) = fill((a, b), end_contact, c, near, &mut budget, tol) else {
                    return false;
                };
                added.extend(between);
            }
            added.push(c);
        }
        added.push(Contact {
            on_a: beside(a, touch.on_a, on_a),
            on_b: beside(b, touch.on_b, on_b),
            point: touch.point,
        });
    } else {
        let (on_a, on_b) = (arc.on_a[end], arc.on_b[end]);
        arc.points[end] = touch.point;
        arc.on_a[end] = beside(a, touch.on_a, on_a);
        arc.on_b[end] = beside(b, touch.on_b, on_b);
        return true;
    }
    if at_head {
        added.reverse();
        arc.points.splice(0..0, added.iter().map(|c| c.point));
        arc.on_a.splice(0..0, added.iter().map(|c| c.on_a));
        arc.on_b.splice(0..0, added.iter().map(|c| c.on_b));
    } else {
        arc.points.extend(added.iter().map(|c| c.point));
        arc.on_a.extend(added.iter().map(|c| c.on_a));
        arc.on_b.extend(added.iter().map(|c| c.on_b));
    }
    true
}

/// Newton projection of a point onto a surface, warm-started.
pub(crate) fn nearest_on(
    surface: &SurfaceGeometry,
    seed: (f64, f64),
    target: Point,
    tol: Tolerances,
) -> Option<((f64, f64), Point)> {
    let (mut u, mut v) = seed;
    for _ in 0..16 {
        let (u_ok, v_ok) = surface.normalize_parameters(u, v, tol).ok()?;
        u = u_ok;
        v = v_ok;
        let p = surface.point_at(u, v, tol).ok()?;
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let r = p - target;
        let (a11, a12, a22) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = a11.mul_add(a22, -(a12 * a12));
        if det.abs() <= f64::MIN_POSITIVE {
            break;
        }
        let (b1, b2) = (-su.dot(r), -sv.dot(r));
        let du = b1.mul_add(a22, -(b2 * a12)) / det;
        let dv = a11.mul_add(b2, -(a12 * b1)) / det;
        u += du;
        v += dv;
        if du.hypot(dv) < 1e-14 {
            break;
        }
    }
    let (u, v) = surface.normalize_parameters(u, v, tol).ok()?;
    Some(((u, v), surface.point_at(u, v, tol).ok()?))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::print_stdout)]
mod tests {
    use super::*;
    use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere};

    const T: Tolerances = Tolerances::millimetres();

    fn cylinder(origin: Point, axis: Vector, radius: f64, height: (f64, f64)) -> SurfaceGeometry {
        let frame = Frame::new(
            origin,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), height)
            .unwrap()
            .into()
    }

    /// A unit drum spelt as its exact clamped patch, which closes round
    /// without repeating, cut by a plane touching one of its end rims off
    /// the patch's seam: `x cos a + y sin a + c d = 1`, `a` the touch's
    /// angle and `d` the depth into the drum from that rim. The section is
    /// the ellipse of semi-axes 1 and `b = √(1 + 1/c²)`, touching the rim.
    /// A walk from a seed on one side of the touch stops on it, and one
    /// from the other side may run over it and round the seam: the arc the
    /// first traced is then a stretch of the second, and is returned once,
    /// and the stretch from the touch to the seam is returned too. So the
    /// arcs' polylines together run as long as the ellipse's perimeter `P`,
    /// by Simpson's rule, less what chords sagging at most the march's
    /// chord `h` cut off: a chord of a stretch of curvature `κ` sagging `s`
    /// is shorter than it by `κ·s·L/3` to leading order, `L` its length, at
    /// most `b·h·L/3` here, and at most `b·h·P/3` all round.
    #[test]
    fn a_plane_touching_a_spline_drums_rim_is_traced_once() {
        let drum: SurfaceGeometry = SurfaceGeometry::from(
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (0.0, 2.0)).unwrap(),
        )
        .to_bspline(T)
        .unwrap()
        .into();
        for c in [1.5_f64, 2.5, 4.0] {
            let b = (1.0 + 1.0 / (c * c)).sqrt();
            let quarter = core::f64::consts::FRAC_PI_2 / 512.0;
            let perimeter = 4.0 * quarter / 3.0
                * (0..=512)
                    .map(|k| {
                        let t = quarter * f64::from(k);
                        let weight = if k == 0 || k == 512 {
                            1.0
                        } else if k % 2 == 1 {
                            4.0
                        } else {
                            2.0
                        };
                        weight * t.sin().hypot(b * t.cos())
                    })
                    .sum::<f64>();
            for k in 0..12 {
                let a = f64::from(k) / 12.0 * core::f64::consts::TAU + 0.1;
                for (height, up) in [(0.0, 1.0), (2.0, -1.0)] {
                    let plane: SurfaceGeometry = PlaneSurface::over(
                        Plane::through(
                            Point::new(a.cos(), a.sin(), height),
                            Direction::new(Vector::new(a.cos(), a.sin(), up * c), T).unwrap(),
                        ),
                        (-8.0, 8.0),
                        (-8.0, 8.0),
                    )
                    .unwrap()
                    .into();
                    let found = branches(&drum, &plane, Marching::default(), T).unwrap();
                    let length: f64 = found
                        .iter()
                        .flat_map(|branch| branch.points.windows(2))
                        .map(|pair| pair[0].distance(pair[1]))
                        .sum();
                    // A polyline through samples of a convex loop runs
                    // shorter than the loop. A stretch returned twice adds
                    // its length, and one missed takes its length away.
                    let cut_off = b * Marching::default().chord * perimeter / 3.0;
                    assert!(
                        length < perimeter && perimeter - length <= cut_off,
                        "c {c}, touch at {a} on the rim at {height}: {} arcs run {length}, \
                         the ellipse {perimeter}, chords cut off at most {cut_off}",
                        found.len()
                    );
                }
            }
        }
    }

    /// A drum of radius 0.5 along x, its axis at y = `c` and z = 1.5,
    /// ending in the plane x = 2.5 tangent to a bore of radius 2.5 along z.
    /// The drum's wall meets the bore in a loop whose x is greatest, 2.5,
    /// where it crosses y = 0: at z = 1.5 ± √(0.25 − c²), two points of the
    /// drum's end rim, which the loop touches and turns back from. So the
    /// loop is two arcs between the touches, each ending on both, with
    /// either surface first and the rim at either end of the drum's height.
    /// At c = 0.45 and 0.475 the arc below y = 0 stands within 0.0005 and
    /// 0.00013 of the rim. A walk ending short of a touch or running over
    /// it leaves an arc's end microns to millimetres off the rim's point,
    /// or the loop whole.
    #[test]
    fn a_section_touching_a_drum_rim_ends_on_the_touch() {
        let bore: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.5, T).unwrap(), (0.0, 3.0))
                .unwrap()
                .into();
        for c in [0.2, 0.45, 0.475] {
            let lift = (0.25_f64 - c * c).sqrt();
            let touches = [
                Point::new(2.5, 0.0, 1.5 - lift),
                Point::new(2.5, 0.0, 1.5 + lift),
            ];
            for (origin, axis) in [
                (Point::new(1.0, c, 1.5), Direction::X),
                (
                    Point::new(2.5, c, 1.5),
                    Direction::new(-Vector::X, T).unwrap(),
                ),
            ] {
                let frame = Frame::new(origin, axis, Direction::Z, T).unwrap();
                let drum: SurfaceGeometry =
                    CylinderSurface::new(Cylinder::new(frame, 0.5, T).unwrap(), (0.0, 1.5))
                        .unwrap()
                        .into();
                for (a, b) in [(&bore, &drum), (&drum, &bore)] {
                    let found = branches(a, b, Marching::default(), T).unwrap();
                    assert_eq!(found.len(), 2, "{c}, {axis:?}: the loop's two arcs");
                    for touch in touches {
                        let ends = found
                            .iter()
                            .flat_map(|arc| [arc.points[0], arc.points[arc.points.len() - 1]])
                            .filter(|end| end.distance(touch) <= T.confusion())
                            .count();
                        assert_eq!(ends, 2, "{c}, {axis:?}: the arcs ending on {touch:?}");
                    }
                }
            }
        }
    }

    /// Two drums across each other, the second's foot a little above the
    /// first's axis and the first's end a little past the second's axis.
    /// Their crossing is one loop, and the domains cut it into two arcs:
    /// a long one, and one 0.09 long between the first drum's rim and the
    /// second's, inside one pair of cells and crossing neither. It is
    /// seeded where the rims pierce the other drum.
    #[test]
    fn a_short_arc_cut_off_between_two_rims_is_found() {
        let drum = |origin: Point, axis: Direction, x: Direction, radius: f64, height: f64| {
            let frame = Frame::new(origin, axis, x, T).unwrap();
            SurfaceGeometry::from(
                CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (0.0, height))
                    .unwrap(),
            )
        };
        let a = drum(
            Point::new(-1.419_880_124_370_016_3, -1.0, 3.0),
            Direction::X,
            Direction::Y,
            1.773_938_116_486_465_8,
            2.5,
        );
        let b = drum(
            Point::new(0.510_316_713_993_857_4, 1.5, 2.230_205_316_817_85),
            Direction::Z,
            Direction::X,
            1.040_566_819_682_834_4,
            4.0,
        );
        let found = branches(&a, &b, Marching::default(), T).unwrap();
        // The short arc runs between the two rims, all of it past
        // x = 1.02; the long one reaches back to x = -0.24.
        assert!(
            found
                .iter()
                .any(|branch| branch.points.iter().all(|p| p.x > 1.02)),
            "{} branches, none the short arc",
            found.len()
        );
    }

    /// The patch over `0 ≤ x ≤ 1`, `0 ≤ y ≤ 2` standing `z = s y + q y²`
    /// off the plane `z = 0`: linear along `x`, its quadratic in `y` with
    /// the middle control point on the start's heading half way. Its
    /// border `y = 0` lies in the plane, and it leaves the plane there at
    /// the sine `s / √(1 + s²)`, bowing on at `q`.
    fn leaving(s: f64, q: f64) -> SurfaceGeometry {
        let points: Vec<Point> = [0.0, 1.0]
            .into_iter()
            .flat_map(|x| {
                [(0.0, 0.0), (1.0, s), (2.0, 2.0 * s + 4.0 * q)].map(|(y, z)| Point::new(x, y, z))
            })
            .collect();
        ogeom_geom::BSplineSurface::new(
            ogeom_math::KnotVector::clamped_uniform(1, 2).unwrap(),
            ogeom_math::KnotVector::clamped_uniform(2, 3).unwrap(),
            &ogeom_math::ControlGrid::new(points, 2, 3).unwrap(),
            T,
        )
        .unwrap()
        .into()
    }

    /// A border lying in the other surface is a stretch of the section
    /// however shallowly the surface leaves it, so long as the angle takes
    /// it off: a sine `s` puts it `s d` off a distance `d` in, past the
    /// confusion distance by the grid's first row in, here `d = 1/12`. One
    /// leaving tangentially, a sine of naught there and rising only as it
    /// bows, is a contact and gives no stretch, and so is one leaving at
    /// `2e-4` (the little a fit of a tangent blend leaves) while bowing at a
    /// curvature of one, which stands it off two hundred times further there.
    #[test]
    fn a_border_left_shallowly_is_a_stretch_and_one_left_tangentially_is_not() {
        let ground = plane(Point::ORIGIN, Vector::Z);
        let options = Marching::default();
        for (s, q, stretches) in [
            (0.002, 0.005, 1),
            (0.01, 0.0, 1),
            (0.02, 0.02, 1),
            (0.0, 0.02, 0),
            (2e-4, 0.5, 0),
        ] {
            let patch = leaving(s, q);
            for (a, b) in [(&ground, &patch), (&patch, &ground)] {
                let mesh = (sample(a, options.grid, T), sample(b, options.grid, T));
                let found = border_stretches(a, b, (&mesh.0, &mesh.1), options, T).unwrap();
                assert_eq!(found.len(), stretches, "sine {s}, bow {q}: {found:?}");
                for piece in &found {
                    let ends = [piece.points[0], piece.points[piece.points.len() - 1]];
                    for p in &piece.points {
                        assert!(p.y.abs().max(p.z.abs()) <= T.confusion(), "sine {s}: {p:?}");
                    }
                    let mut xs = ends.map(|p| p.x);
                    xs.sort_by(f64::total_cmp);
                    assert!(
                        xs[0].abs() <= T.confusion() && (xs[1] - 1.0).abs() <= T.confusion(),
                        "sine {s}: the stretch runs {xs:?}"
                    );
                }
            }
        }
    }

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-8.0, 8.0),
            (-8.0, 8.0),
        )
        .unwrap()
        .into()
    }

    /// How far a point is from a quadric, in closed form.
    fn off(surface: &SurfaceGeometry, p: Point) -> f64 {
        match surface {
            SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
            SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
            SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
            _ => 0.0,
        }
    }

    /// The worst distance from a traced branch to either surface.
    fn deviation(a: &SurfaceGeometry, b: &SurfaceGeometry, traced: &Traced) -> f64 {
        traced
            .points
            .iter()
            .map(|p| off(a, *p).abs().max(off(b, *p).abs()))
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn a_plane_through_a_bent_strip_seeds_both_branches() {
        // A cubic strip bent into an arch crosses a level plane twice. The
        // plane's domain is the unbounded one a face carries, and its
        // sampling cells span a million units; merging seeds at *its*
        // spacing would call the two crossings one branch and trace only
        // one of them. Seeds merge at the finer surface's spacing instead.
        use ogeom_geom::BSplineSurface;
        use ogeom_math::{ControlGrid, KnotVector};
        let mut points = Vec::new();
        for i in 0..7 {
            let a = core::f64::consts::PI * f64::from(i) / 6.0;
            for j in 0..2 {
                points.push(Point::new(2.0 * a.cos(), f64::from(j), 2.0 * a.sin()));
            }
        }
        let grid = ControlGrid::new(points, 7, 2).unwrap();
        let strip: SurfaceGeometry = BSplineSurface::new(
            KnotVector::clamped_uniform(3, 7).unwrap(),
            KnotVector::clamped_uniform(1, 2).unwrap(),
            &grid,
            T,
        )
        .unwrap()
        .into();
        let level: SurfaceGeometry = PlaneSurface::over(
            Plane::through(Point::new(0.0, 0.0, 1.0), Direction::Z),
            (-1.0e9, 1.0e9),
            (-1.0e9, 1.0e9),
        )
        .unwrap()
        .into();
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };
        let found = branches(&strip, &level, options, T).unwrap();
        assert_eq!(
            found.len(),
            2,
            "the arch crosses the level twice: {}",
            found.len()
        );
        for branch in &found {
            assert!(!branch.closed());
            for p in &branch.points {
                assert!((p.z - 1.0).abs() < 1e-4, "on the level: {p:?}");
            }
        }
    }

    #[test]
    fn two_crossed_cylinders_are_traced_onto_both_of_them() {
        // The case with no closed form, and the one the analytic module
        // explicitly refuses: two cylinders on perpendicular axes meet in a
        // quartic space curve. Every point of the trace must be on both.
        let a = cylinder(Point::ORIGIN, Vector::Z, 1.0, (-4.0, 4.0));
        let b = cylinder(Point::ORIGIN, Vector::X, 1.0, (-4.0, 4.0));
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };

        // Two equal cylinders crossing at right angles meet in two ellipses,
        // the Steinmetz solid's seams, crossing where the cylinders touch at
        // (0, -1, 0) and (0, 1, 0): four arcs, each from one touch to the
        // other.
        let found = branches(&a, &b, options, T).unwrap();
        assert_eq!(found.len(), 4, "the two ellipses' halves");

        let touches = [Point::new(0.0, -1.0, 0.0), Point::new(0.0, 1.0, 0.0)];
        let mut worst = 0.0_f64;
        for branch in &found {
            let ends = [branch.points[0], branch.points[branch.points.len() - 1]];
            for touch in touches {
                assert!(
                    ends.iter().any(|end| end.distance(touch) < 1e-9),
                    "an arc from {:?} to {:?} misses the touch {touch:?}",
                    ends[0],
                    ends[1]
                );
            }
            assert!(
                branch.points.len() > 100,
                "a branch of only {} points",
                branch.points.len()
            );
            worst = worst.max(deviation(&a, &b, branch));
        }
        println!(
            "crossed cylinders: {} branches, worst deviation {worst:e}",
            found.len()
        );
        assert!(worst < 1e-7, "traced off the surfaces by {worst:e}");
    }

    /// A cylinder lying inside a torus's outer equator and touching it at
    /// (13, 0, 0), where the torus's chart has its corner, meets it in a
    /// figure eight. Traced either way round, the section comes back as the
    /// two loops' four halves, each ending on the touch, and the loops'
    /// length is the section's.
    #[test]
    fn a_figure_eight_is_cut_at_its_double_point() {
        let radius = 4.353_623_591_855_474;
        let torus: SurfaceGeometry = ogeom_geom::TorusSurface::new(
            ogeom_math::Torus::new(Frame::WORLD, 10.0, 3.0, T).unwrap(),
        )
        .into();
        let frame = Frame::new(
            Point::new(13.0 - radius, 0.0, -10.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let drill: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (0.0, 20.0))
                .unwrap()
                .into();
        let options = Marching {
            chord: 6e-6,
            ..Marching::default()
        };
        // The section's length: the drill's wall unrolled, z the root of 9
        // less the square of the distance from the tube's centre circle, at
        // each angle round the drill where that is positive, four times
        // over for the two loops' upper and lower halves.
        let steps = 200_000;
        let mut expected = 0.0;
        let at = |t: f64| {
            let (x, y) = (13.0 - radius + radius * t.cos(), radius * t.sin());
            let off = x.hypot(y) - 10.0;
            ((9.0 - off * off).max(0.0).sqrt(), x, y)
        };
        for k in 0..steps {
            let (t0, t1) = (
                core::f64::consts::TAU * f64::from(k) / f64::from(steps),
                core::f64::consts::TAU * f64::from(k + 1) / f64::from(steps),
            );
            let ((z0, x0, y0), (z1, x1, y1)) = (at(t0), at(t1));
            if z0 > 0.0 || z1 > 0.0 {
                expected += 2.0 * Point::new(x0, y0, z0).distance(Point::new(x1, y1, z1));
            }
        }
        let touch = Point::new(13.0, 0.0, 0.0);
        for (a, b) in [(&torus, &drill), (&drill, &torus)] {
            let found = branches(a, b, options, T).unwrap();
            assert_eq!(found.len(), 4, "the figure eight's four halves");
            let mut length = 0.0;
            for branch in &found {
                let ends = [branch.points[0], branch.points[branch.points.len() - 1]];
                assert!(
                    ends.iter().filter(|end| end.distance(touch) < 1e-9).count() == 1,
                    "a half from {:?} to {:?} does not end once on the touch",
                    ends[0],
                    ends[1]
                );
                length += branch
                    .points
                    .windows(2)
                    .map(|w| w[0].distance(w[1]))
                    .sum::<f64>();
                assert!(deviation(a, b, branch) < 1e-6, "traced off the surfaces");
            }
            assert!(
                (length - expected).abs() < 1e-3,
                "the halves run {length} against the section's {expected}"
            );
        }
    }

    #[test]
    fn unequal_crossed_cylinders_meet_in_two_curves_as_well() {
        // The non-degenerate cousin. Equal radii put the two curves through
        // each other at the tangency points, so getting the count right there
        // says less than getting it right here, where there is no singularity
        // for a tracer to be lucky about.
        let a = cylinder(Point::ORIGIN, Vector::Z, 1.0, (-4.0, 4.0));
        let b = cylinder(Point::ORIGIN, Vector::X, 1.6, (-4.0, 4.0));
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };

        let found = branches(&a, &b, options, T).unwrap();
        assert_eq!(found.len(), 2);
        for branch in &found {
            assert!(branch.closed());
            assert!(deviation(&a, &b, branch) < 1e-7);
        }
    }

    #[test]
    fn a_traced_circle_agrees_with_the_circle_it_should_be() {
        // A sphere cut by a plane through its centre is a circle of known
        // radius, and the marcher does not know that. Tracing it and checking
        // against the closed form is the strongest single check there is: it
        // tests the tracer against an answer derived independently of it.
        let s = sphere(Point::ORIGIN, 3.0);
        let cut = plane(Point::ORIGIN, Vector::Z);
        let options = Marching {
            chord: 1e-6,
            ..Marching::default()
        };

        let found = seeds(&s, &cut, options, T).unwrap();
        assert!(!found.is_empty());
        let branch = trace(&s, &cut, found[0], options, T).unwrap();

        assert!(
            branch.closed(),
            "a plane through a sphere gives a closed loop"
        );
        for p in &branch.points {
            let radius = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (radius - 3.0).abs() < 1e-7,
                "a point at radius {radius} on a circle of 3"
            );
            assert!(p.z.abs() < 1e-7, "off the cutting plane by {}", p.z);
        }
    }

    /// A branch the backward half closes is the loop, whatever stopped the
    /// forward half. An oblique plane cuts a drum in an ellipse, which takes
    /// a different number of steps round each way from most seeds. From one
    /// that goes round in fewer backwards, with points enough only for that,
    /// the forward half runs out and the backward one closes: joined, the
    /// two ran past the seed a second time and were reported truncated.
    #[test]
    fn a_branch_closed_going_backwards_is_the_loop() {
        let drum = cylinder(Point::new(0.0, 0.0, -4.0), Vector::Z, 1.0, (0.0, 8.0));
        let cut = plane(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0));
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };
        let found = seeds(&drum, &cut, options, T).unwrap();
        let round = trace(&drum, &cut, found[0], options, T).unwrap();
        let Some((from, ahead, behind)) = (0..round.points.len())
            .step_by(37)
            .map(|k| {
                let from = Contact {
                    on_a: round.on_a[k],
                    on_b: round.on_b[k],
                    point: round.points[k],
                };
                let ahead = walk(&drum, &cut, from, 1.0, options, T).unwrap();
                let behind = walk(&drum, &cut, from, -1.0, options, T).unwrap();
                (from, ahead, behind)
            })
            .find(|(_, ahead, behind)| {
                ahead.stopped == Stopped::Closed
                    && behind.stopped == Stopped::Closed
                    && behind.points.len() < ahead.points.len()
            })
        else {
            panic!("a seed the backward walk closes from in fewer points");
        };

        let options = Marching {
            max_points: behind.points.len(),
            ..options
        };
        let branch = trace(&drum, &cut, from, options, T).unwrap();
        assert_eq!(branch.stopped, Stopped::Closed, "the backward half closed");
        assert_eq!(branch.points.len(), behind.points.len(), "once round");
        assert_eq!(branch.points[0], branch.points[branch.points.len() - 1]);
        assert!(
            (branch.points[1] - from.point).dot(ahead.points[1] - from.point) > 0.0,
            "it runs the way the forward half set off"
        );
    }

    #[test]
    fn a_branch_that_leaves_the_surface_says_so_rather_than_stopping_quietly() {
        // A truncated branch and a finished one are different answers, and a
        // caller that cannot tell them apart treats the first as the second.
        let s = sphere(Point::ORIGIN, 3.0);
        let cut = plane(Point::new(0.0, 0.0, 0.0), Vector::Z);
        let options = Marching {
            chord: 1e-4,
            max_points: 8,
            ..Marching::default()
        };
        let found = seeds(&s, &cut, options, T).unwrap();
        let branch = trace(&s, &cut, found[0], options, T).unwrap();
        assert_eq!(branch.stopped, Stopped::RanOut);
        assert!(!branch.complete(), "a truncated branch is not complete");
    }

    #[test]
    fn tangent_surfaces_are_refused_rather_than_followed_onto_a_guess() {
        // Where the normals are parallel the intersection has no single
        // direction, and marching through such a point is how a tracer changes
        // branch without noticing. A sphere resting on a plane is the case.
        let s = sphere(Point::new(0.0, 0.0, 3.0), 3.0);
        let ground = plane(Point::ORIGIN, Vector::Z);
        let touch = Contact {
            on_a: (0.0, -core::f64::consts::FRAC_PI_2),
            on_b: (0.0, 0.0),
            point: Point::ORIGIN,
        };
        let err = trace(&s, &ground, touch, Marching::default(), T).unwrap_err();
        assert!(err.to_string().contains("tangent"), "unexpected: {err}");
    }

    #[test]
    fn the_number_of_branches_is_the_number_there_are() {
        // The failure the accuracy measure cannot see. Every point of one
        // circle is on both surfaces, so returning one of two scores perfectly;
        // the count is the only thing that catches it.
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };

        // A sphere cut by a plane off its centre: one circle.
        let one = branches(
            &sphere(Point::ORIGIN, 3.0),
            &plane(Point::new(0.0, 0.0, 1.0), Vector::Z),
            options,
            T,
        )
        .unwrap();
        assert_eq!(one.len(), 1, "one plane through a sphere cuts one circle");
        assert!(one[0].closed());

        // A sphere and a coaxial cylinder narrower than it: two circles, one
        // above and one below.
        let two = branches(
            &sphere(Point::ORIGIN, 3.0),
            &cylinder(Point::ORIGIN, Vector::Z, 1.5, (-4.0, 4.0)),
            options,
            T,
        )
        .unwrap();
        assert_eq!(two.len(), 2, "a coaxial cylinder cuts a sphere twice");
        for branch in &two {
            assert!(branch.closed(), "each is a closed circle");
        }
        // And they are on opposite sides, rather than the same one twice.
        let heights: Vec<f64> = two.iter().map(|b| b.points[0].z).collect();
        assert!(
            heights[0] * heights[1] < 0.0,
            "both branches came back on the same side: {heights:?}"
        );
    }

    /// A thin drum at right angles to a long drill, its axis passing 46
    /// from the drill's (radius 56), runs in through the drill's wall and
    /// out again: two loops, about the two points where the drum's axis
    /// pierces the wall, 64 apart. That is closer than a drill cell's long
    /// side, and nearer than either surface's length over the grid.
    #[test]
    fn a_thin_drum_through_a_long_drill_s_wall_cuts_it_twice() {
        let (radius, offset) = (56.0, 46.0);
        let drill = cylinder(
            Point::new(-900.0, 0.0, 0.0),
            Vector::X,
            radius,
            (0.0, 1800.0),
        );
        let drum = cylinder(
            Point::new(5.0, 0.0, offset),
            Vector::Y,
            5.0,
            (-1000.0, 1000.0),
        );
        let found = branches(&drum, &drill, Marching::default(), T).unwrap();
        assert_eq!(found.len(), 2, "the drum passes through the wall twice");
        let pierce = (radius * radius - offset * offset).sqrt();
        let mut sides: Vec<f64> = found
            .iter()
            .map(|branch| {
                assert!(branch.closed(), "each is a loop round the drum");
                let mean = branch.points.iter().map(|p| p.y).sum::<f64>();
                #[allow(clippy::cast_precision_loss)]
                let mean = mean / branch.points.len() as f64;
                assert!(
                    (mean.abs() - pierce).abs() < 5.0,
                    "a loop about y = {mean}, the axis pierces at {pierce}"
                );
                mean
            })
            .collect();
        sides.sort_by(f64::total_cmp);
        assert!(sides[0] < 0.0 && sides[1] > 0.0, "{sides:?}");
    }

    #[test]
    fn a_branch_thinner_than_the_sampling_is_found_from_the_surfaces_extremum() {
        // Two spheres barely overlapping meet in a circle far smaller than a
        // coarse grid's cell: no pair of cells crosses on it, and a sphere
        // has no border to pierce. The approach of the two surfaces inside
        // the circle seeds it at any grid, and a finer grid still finds it
        // by crossing cells.
        let a = sphere(Point::ORIGIN, 3.0);
        let b = sphere(Point::new(5.98, 0.0, 0.0), 3.0);
        let radius = (9.0_f64 - 2.99 * 2.99).sqrt();
        for grid in [6, 24, 120] {
            let options = Marching {
                grid,
                ..Marching::default()
            };
            assert!(
                !seeds(&a, &b, options, T).unwrap().is_empty(),
                "grid {grid}"
            );
            let found = branches(&a, &b, options, T).unwrap();
            assert_eq!(found.len(), 1, "one circle at grid {grid}");
            assert!(found[0].closed());
            assert!(deviation(&a, &b, &found[0]) < 1e-7);
            for p in &found[0].points {
                assert!((p.x - 2.99).abs() < 1e-6, "in the mid plane: {p:?}");
                assert!(
                    (p.y.hypot(p.z) - radius).abs() < 1e-6,
                    "on the circle: {p:?}"
                );
            }
        }
    }

    /// A bicubic patch over the square `[0, 10]`, level at z = 0 but for a
    /// bump raised by one control point near the middle, its top at
    /// `peak`.
    fn bumped(peak: f64) -> SurfaceGeometry {
        use ogeom_geom::BSplineSurface;
        use ogeom_math::{ControlGrid, KnotVector};
        const N: usize = 41;
        // A raised control point lifts its patch by four ninths of its own
        // height: the cubic basis peaks at two thirds each way.
        let lift = peak * 9.0 / 4.0;
        let mut points = Vec::new();
        for i in 0..N {
            for j in 0..N {
                #[allow(clippy::cast_precision_loss)]
                let (x, y) = (10.0 * i as f64 / 40.0, 10.0 * j as f64 / 40.0);
                let z = if (i, j) == (19, 22) { lift } else { 0.0 };
                points.push(Point::new(x, y, z));
            }
        }
        BSplineSurface::new(
            KnotVector::clamped_uniform(3, N).unwrap(),
            KnotVector::clamped_uniform(3, N).unwrap(),
            &ControlGrid::new(points, N, N).unwrap(),
            T,
        )
        .unwrap()
        .into()
    }

    /// A bilinear patch level at `height`, wider than the bumped patch.
    fn level_patch(height: f64) -> SurfaceGeometry {
        use ogeom_geom::BSplineSurface;
        use ogeom_math::{ControlGrid, KnotVector};
        let points = vec![
            Point::new(-1.0, -1.0, height),
            Point::new(-1.0, 11.0, height),
            Point::new(11.0, -1.0, height),
            Point::new(11.0, 11.0, height),
        ];
        BSplineSurface::new(
            KnotVector::clamped_uniform(1, 2).unwrap(),
            KnotVector::clamped_uniform(1, 2).unwrap(),
            &ControlGrid::new(points, 2, 2).unwrap(),
            T,
        )
        .unwrap()
        .into()
    }

    #[test]
    fn a_small_bump_through_a_spline_patch_cuts_a_loop_inside_both() {
        // The bump's top pokes a hundredth of its height through a level
        // patch: a closed loop far smaller than a grid cell, touching no
        // border of either patch.
        let peak = 0.5;
        let bump = bumped(peak);
        let level = level_patch(peak * 0.99);
        let found = branches(&bump, &level, Marching::default(), T).unwrap();
        assert_eq!(found.len(), 1, "one loop round the bump's top");
        let loop_ = &found[0];
        assert!(loop_.closed());
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for (p, &(u, v)) in loop_.points.iter().zip(&loop_.on_a) {
            assert!((p.z - peak * 0.99).abs() < 1e-7, "on the level: {p:?}");
            assert!(bump.point_at(u, v, T).unwrap().distance(*p) < 1e-7);
            low = low.min(p.x);
            high = high.max(p.x);
        }
        assert!(
            high - low < 10.0 / 24.0,
            "thinner than a cell: {}",
            high - low
        );

        // Just touching: the top meets the level within tolerance, and a
        // touch has no curve to follow.
        let touching = level_patch(peak);
        assert!(
            branches(&bump, &touching, Marching::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_drill_grazing_a_ball_cuts_one_small_loop() {
        // A drill along x passes a thousandth deeper than tangent to a ball:
        // a loop a few hundredths across, inside one cell of either surface.
        let ball = sphere(Point::ORIGIN, 3.0);
        let drill = cylinder(Point::new(0.0, 0.0, 3.999), Vector::X, 1.0, (-4.0, 4.0));
        let found = branches(&drill, &ball, Marching::default(), T).unwrap();
        assert_eq!(found.len(), 1, "one loop where the drill dips in");
        assert!(found[0].closed());
        assert!(deviation(&drill, &ball, &found[0]) < 1e-7);
        for p in &found[0].points {
            assert!(p.x.abs() < 0.2 && p.y.abs() < 0.1, "near the graze: {p:?}");
        }
    }

    #[test]
    fn settings_that_could_not_work_are_refused() {
        let a = sphere(Point::ORIGIN, 1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        for options in [
            Marching {
                chord: 0.0,
                ..Marching::default()
            },
            Marching {
                grid: 1,
                ..Marching::default()
            },
            Marching {
                max_points: 1,
                ..Marching::default()
            },
        ] {
            assert!(seeds(&a, &b, options, T).is_err());
        }
    }
}
