//! Branches found on a border: spline and extrusion borders, rims, and chart borders.

use super::*;

/// A border of a surface as a curve, and the map from the curve's
/// parameter to the surface's.
type Border = (ogeom_geom::Curve, Box<dyn Fn(f64) -> (f64, f64)>);

/// The open borders of a spline surface, or of a window of one.
pub(super) fn spline_borders(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Border> {
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

/// The open borders of a linear extrusion, or of a window of one: the
/// rulings at the ends of its `u` and the profile carried to the ends of
/// its `v`, each at the extrusion's own parameter along it.
pub(super) fn extrusion_borders(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Border> {
    use ogeom_geom::Curve3d as _;
    let extrusion = match surface {
        SurfaceGeometry::Extrusion(e) => e,
        SurfaceGeometry::Trimmed(trimmed) => match trimmed.basis() {
            SurfaceGeometry::Extrusion(e) => e,
            _ => return Vec::new(),
        },
        _ => return Vec::new(),
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    if ![u0, u1, v0, v1]
        .iter()
        .all(|x| x.is_finite() && x.abs() <= FAR)
    {
        return Vec::new();
    }
    let profile = extrusion.curve();
    let d = extrusion.direction();
    let mut out: Vec<Border> = Vec::new();
    if !surface.is_closed_u(tol) {
        for u in [u0, u1] {
            let Ok(at) = profile.point_at(u, tol) else {
                continue;
            };
            let axis = ogeom_math::Axis {
                location: at,
                direction: d,
            };
            if let Ok(line) = ogeom_geom::LineCurve::over(axis, v0, v1) {
                out.push((line.into(), Box::new(move |t| (u, t))));
            }
        }
    }
    let (p0, p1) = profile.domain();
    for v in [v0, v1] {
        let moved = ogeom_math::Transform::translation(d.vector() * v);
        let Ok(carried) = ogeom_geom::Transformable::transformed(profile, &moved, tol) else {
            continue;
        };
        let iso = if u0 <= p0 + tol.parametric() && u1 >= p1 - tol.parametric() {
            Some(carried)
        } else {
            ogeom_geom::TrimmedCurve::new(carried, u0, u1, tol)
                .ok()
                .map(Into::into)
        };
        if let Some(iso) = iso {
            out.push((iso, Box::new(move |t| (t, v))));
        }
    }
    out
}

/// The rims of a drum or a cone: the circles at the ends of its height,
/// where those ends are a part's height and not a chart's clamp. On a
/// window of one less than a turn wide, the arcs of those circles the
/// window spans.
pub(super) fn rims(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Border> {
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
pub(super) fn border_stretches(
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
