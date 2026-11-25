//! Seeds: starting points on the intersection, from both surfaces sampled into cells.

use super::*;

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
    Ok(seeded(a, b, options, tol)?.found)
}

/// What [`seeded`] finds on one sampling of a pair.
pub(super) struct Seeded {
    /// One seed per branch, as [`seeds`] returns them.
    pub(super) found: Vec<Contact>,
    /// The contacts within reach of a seed where the surfaces are tangent,
    /// and of no other: see [`merges`].
    pub(super) beside_tangency: Vec<Contact>,
    /// The stretches of either surface's borders lying in the other: see
    /// [`border_stretches`].
    pub(super) lying: Vec<Traced>,
}

/// [`seeds`], the contacts held back beside a tangency, and the stretches
/// of either surface's borders lying in the other, found on the same
/// sampling.
pub(super) fn seeded(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Seeded> {
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
    let mut beside_tangency: Vec<Contact> = Vec::new();
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
            // A contact held on a chart's bound slides along it by its
            // residual over the sine the section crosses it at, and is
            // solved again on it before it is weighed against the seeds.
            let contact = on_its_bound(a, b, contact, options.chord, tol);
            // One seed per branch, not one per cell it passes through. The
            // reach is the *finer* surface's grid: the coarser surface's
            // cells say nothing about how far apart branches can be, and a
            // plane's clamped domain spans a million units.
            match merges(&found, contact, apart) {
                Merged::Along => {}
                Merged::BesideTangency => beside_tangency.push(contact),
                Merged::Not => found.push(Seed::at(a, b, contact, tol)),
            }
        }
    }
    // A branch that runs in from a border at a grazing angle can be
    // thinner than the sampling's sag, and a short one cut off in a corner
    // of both domains can lie inside a single pair of cells without
    // crossing them: no pair of cells crosses on it. Where it meets the
    // border it is a curve piercing a surface, which is found exactly:
    // each border of each spline and of each extrusion, and each rim of a
    // drum or a cone, is intersected with the other surface, and every
    // piercing seeds.
    for (from_a, border_of, other) in [(true, a, b), (false, b, a)] {
        let borders = spline_borders(border_of, tol)
            .into_iter()
            .chain(rims(border_of, tol))
            .chain(extrusion_borders(border_of, tol));
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
                let contact = on_its_bound(a, b, contact, options.chord, tol);
                match merges(&found, contact, apart) {
                    Merged::Along => {}
                    Merged::BesideTangency => beside_tangency.push(contact),
                    Merged::Not => found.push(Seed::at(a, b, contact, tol)),
                }
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
        let contact = on_its_bound(a, b, contact, options.chord, tol);
        match merges(&found, contact, apart) {
            Merged::Along => {}
            Merged::BesideTangency => beside_tangency.push(contact),
            Merged::Not => found.push(Seed::at(a, b, contact, tol)),
        }
    }
    let lying = border_stretches(a, b, (&mesh_a, &mesh_b), options, tol)?;
    Ok(Seeded {
        found: found.into_iter().map(|s| s.contact).collect(),
        beside_tangency,
        lying,
    })
}

/// How many times a cell's sag the surface is allowed to bow from the
/// cell when deciding whether two cells could touch: the sag is measured
/// at the cell's middle only.
pub(super) const BOW: f64 = 2.0;

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
pub(super) fn point_to_triangle(p: Point, t: [Point; 3]) -> f64 {
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
pub(crate) const STATIONARY_STEPS: usize = 12;

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
    let sine = crossing_sine(a, b, at.on_a, at.on_b, tol);
    sine > f64::MIN_POSITIVE && parts_either_side(a, b, at, 4.0 * tol.confusion() / sine, tol)
}

/// Whether the first surface stands off the second by more than the
/// confusion distance, on opposite sides, a distance `h` either way across
/// the intersection from `at`.
fn parts_either_side(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Contact,
    h: f64,
    tol: Tolerances,
) -> bool {
    let Some(along) = tangent_at(a, b, at, tol) else {
        return false;
    };
    let (Some(na), Ok((du, dv))) = (
        normal_at(a, at.on_a, tol),
        a.d1_at(at.on_a.0, at.on_a.1, tol),
    ) else {
        return false;
    };
    let Some(d) = params_along(du, dv, na.cross(along)) else {
        return false;
    };
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

/// Whether two surfaces that are their own geometry, not a fit of another,
/// cross at `at` by more than the confusion distance `w` every marched
/// point stands within of both, rather than touch there.
///
/// A distance `x` across the intersection on the first surface, the first
/// stands off the second by `s x + k x² / 2`: `s` the sine between their
/// normals, `k` the difference of their normal curvatures across it.
/// - Touching along a curve, a marched point within `w` of both lies `x`
///   off the touch with `k x² / 2 ≤ w`, where the sine reads `k x`: so
///   `s² ≤ 2 k w`.
/// - Crossing, the two part on one side and bend back towards each other
///   on the other, `s² / 2k` apart at the deepest, `s / k` across. That is
///   past `w` exactly where `s² > 2 k w`: the one bound tells a touch from
///   a crossing that can be measured, both ways.
///
/// Past the bound the separation is measured either side: `s / k` across,
/// or `4 w / s` where that is nearer, there `4 w − 8 k w² / s²`, at least
/// `2 w`. A touch keeps one side either way, so no term the model leaves
/// out reads one as a crossing.
pub(crate) fn crosses_measurably(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Contact,
    tol: Tolerances,
) -> bool {
    let w = tol.confusion();
    let sine = crossing_sine(a, b, at.on_a, at.on_b, tol);
    let (Some(along), Some(na)) = (tangent_at(a, b, at, tol), normal_at(a, at.on_a, tol)) else {
        return false;
    };
    let across = na.cross(along);
    let bending = |s: &SurfaceGeometry, (u, v): (f64, f64)| -> Option<f64> {
        let curvature = s.curvature_at(u, v, tol).ok()?;
        let k = curvature.normal_curvature(across)?;
        Some(if curvature.normal.vector().dot(na) < 0.0 {
            -k
        } else {
            k
        })
    };
    let (Some(ka), Some(kb)) = (bending(a, at.on_a), bending(b, at.on_b)) else {
        return false;
    };
    let k = (ka - kb).abs();
    if sine * sine <= 2.0 * k * w {
        return false;
    }
    let near = 4.0 * w / sine;
    let h = if k * near > sine { sine / k } else { near };
    parts_either_side(a, b, at, h, tol)
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

/// Whether a contact is a seed already kept: see [`merges`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Merged {
    /// Met again further along a kept seed's branch.
    Along,
    /// Within reach of a kept seed where the surfaces are tangent, and of
    /// no other.
    BesideTangency,
    /// A seed of its own.
    Not,
}

/// Whether a contact is a seed already kept, met again further along its
/// branch: within `apart` of it, and off the line the intersection runs
/// through it by less than a quarter of the way there. Two branches side by
/// side, or a loop beside another the same drum cuts, lie across that line
/// rather than along it, and each keeps a seed of its own however close.
/// Where the surfaces are all but tangent the intersection has no direction
/// to ask, and nearness alone decides; but [`trace`] follows nothing from
/// such a seed, and the branches leaving a node where the section crosses
/// itself pass within `apart` of it, so a contact held by it alone is held
/// back rather than dropped: see [`branches_from_touches`].
fn merges(found: &[Seed], contact: Contact, apart: f64) -> Merged {
    let mut beside_tangency = false;
    for seed in found {
        let offset = contact.point - seed.contact.point;
        let gap = offset.magnitude();
        if gap > apart {
            continue;
        }
        match seed.along {
            Some(along) if offset.cross(along).magnitude() <= 0.25 * gap => return Merged::Along,
            Some(_) => {}
            None => beside_tangency = true,
        }
    }
    if beside_tangency {
        Merged::BesideTangency
    } else {
        Merged::Not
    }
}

/// A contact held on a bound of either chart that does not repeat, solved
/// on that bound to within `within` of where the section crosses it, as
/// [`settled_from`] solves; as it was where that solve does not settle.
///
/// The solve that lands a contact stops once the surfaces stand within a
/// residual of each other, and a point held on a bound slides along it by
/// that residual over the sine the section crosses the bound at: beside a
/// tangency, where the sine falls towards zero, further than the chord
/// from where the section leaves the chart.
fn on_its_bound(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    contact: Contact,
    within: f64,
    tol: Tolerances,
) -> Contact {
    let x = [
        contact.on_a.0,
        contact.on_a.1,
        contact.on_b.0,
        contact.on_b.1,
    ];
    let ((ua, ub), (va, vb)) = a.domain();
    let ((uc, ud), (vc, vd)) = b.domain();
    let bounds = [
        (a.is_periodic_u(), ua, ub),
        (a.is_periodic_v(), va, vb),
        (b.is_periodic_u(), uc, ud),
        (b.is_periodic_v(), vc, vd),
    ];
    let Some(held) = (0..4).find(|&k| {
        let (periodic, lo, hi) = bounds[k];
        !periodic && (x[k] == lo || x[k] == hi)
    }) else {
        return contact;
    };
    let mut row = [0.0; 4];
    row[held] = 1.0;
    let on = |y: &[f64; 4], _: Point, _: [Vector; 2]| (y[held] - x[held], row);
    match settled_from((a, b), x, &on, within, tol) {
        Landing::At(landed) => landed,
        _ => contact,
    }
}

/// Past this a chart's bound is a stand-in for no end at all.
pub(super) const FAR: f64 = 1e6;

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
