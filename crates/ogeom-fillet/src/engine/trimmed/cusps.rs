//! Cusps and lobes at vertex patches, and the conic sections between them.

use super::*;

/// A vertex patch's ring: its vertex, each side's sheet edge, whether the
/// ring walks it forward and the strip it ends, and the region the patch is
/// proved within.
pub(super) type PatchRing = (Id, Vec<(usize, bool, usize)>, ogeom_algo::HalfSpaces);

/// A vertex where a face narrows to nothing between a convex straight
/// round and a concave rim's round running tangent to it (a boss tangent
/// to a plate's side, its foot and the plate's top edge there rounded):
/// the top's loop passes the vertex twice, the rounds' contacts crossing on
/// it either side. Each straight strip ends where the boss's meridian plane
/// through its crossing cuts it, so it runs on tangent into the rim's
/// round's end in that plane: one S side from the side face over the
/// crossing to the wall. The hole between the two S sides, the wall along
/// the rim's contact and the side along the straight rounds' closes on one
/// patch.
pub(super) struct Cusp {
    pub(super) v: Id,
    /// The top, the side and the wall.
    pub(super) faces: [usize; 3],
    /// Each straight strip's side of the vertex: two, or one where the top
    /// ends at the vertex.
    pub(super) lobes: Vec<Lobe>,
    /// One lobe's S's feet: on the wall from the lobe's end there down to
    /// the vertex, and on the side from the vertex to the lobe's end there;
    /// its hole is the lobe and its two feet. Empty for two lobes, whose
    /// hole is the two and the strips' contacts between them.
    pub(super) feet: Vec<(Curve, Curve)>,
    /// Two lobes: where the wall's own seam leaves the rim's contact
    /// between them, where it has one there.
    pub(super) foot: Option<Point>,
}

/// A side of a cusp's hole: its curve over `range`, whether the walk runs
/// with it, and its edges, each the sheet's edge over a range of the curve
/// with the surface across it, in the curve's order.
pub(super) struct HoleSpec {
    pub(super) curve: Curve,
    pub(super) range: (f64, f64),
    pub(super) forward: bool,
    pub(super) pieces: Vec<(usize, (f64, f64), SurfaceGeometry)>,
}

impl HoleSpec {
    /// The side walked the other way.
    pub(super) fn back(self) -> Self {
        Self {
            forward: !self.forward,
            ..self
        }
    }

    /// Where the walk starts.
    pub(super) fn start(&self, tol: Tolerances) -> OgeomResult<Point> {
        self.curve.point_at(
            if self.forward {
                self.range.0
            } else {
                self.range.1
            },
            tol,
        )
    }

    /// The patch's side, its edges the sheet's `edges`.
    pub(super) fn side(&self, edges: &[Shape]) -> crate::round::transfinite::HoleSide {
        let (a, b) = self.range;
        let along = |t: f64| {
            let x = (t - a) / (b - a);
            if self.forward { x } else { 1.0 - x }
        };
        let mut walked: Vec<&(usize, (f64, f64), SurfaceGeometry)> = self.pieces.iter().collect();
        if !self.forward {
            walked.reverse();
        }
        let first = walked[0];
        crate::round::transfinite::HoleSide {
            edge: edges[first.0].clone(),
            curve: self.curve.clone(),
            range: self.range,
            forward: self.forward,
            support: first.2.clone(),
            pieces: if walked.len() == 1 {
                Vec::new()
            } else {
                walked
                    .iter()
                    .map(|(k, (t0, t1), support)| {
                        let (x, y) = (along(*t0), along(*t1));
                        crate::round::transfinite::SidePiece {
                            edge: edges[*k].clone(),
                            range: (*t0, *t1),
                            forward: self.forward,
                            along: (x.min(y), x.max(y)),
                            support: support.clone(),
                        }
                    })
                    .collect()
            },
        }
    }
}

/// The sheet's edge on `curve` over `range` from point `ends.0` to
/// `ends.1`, kept as the strip's end between them.
pub(super) fn piece_of(
    edges: &mut Vec<SheetCurve>,
    made: &mut FastMap<Key, usize>,
    curve: &Curve,
    ends: (usize, usize),
    range: (f64, f64),
) -> usize {
    edges.push(SheetCurve {
        ends,
        curve: curve.clone(),
        range,
    });
    made.insert(
        Key::Section(ends.0.min(ends.1), ends.0.max(ends.1)),
        edges.len() - 1,
    );
    edges.len() - 1
}

/// One straight strip's side of a [`Cusp`].
pub(super) struct Lobe {
    /// The strip.
    pub(super) line: usize,
    /// The rim's strip it crosses on the top: the rim's one closed edge,
    /// or its arc on this side.
    pub(super) rim: usize,
    /// Its crossing with the rim's contact on the top.
    pub(super) top: Point,
    /// Where the meridian plane through the crossing cuts its contact on
    /// the side, and the rim's on the wall.
    pub(super) side: Point,
    pub(super) wall: Point,
    /// That plane.
    pub(super) plane: SurfaceGeometry,
}

/// The [`Cusp`]s among `curved`'s vertices.
///
/// # Errors
///
/// As making the planes.
pub(super) fn cusps_of(
    curved: &Curved,
    bevels: &[Strip],
    corners: &FastMap<(usize, Id), Corner>,
    tol: Tolerances,
) -> OgeomResult<Vec<Cusp>> {
    let mut found = Vec::new();
    let mut keyed: Vec<&Id> = curved.points.keys().collect();
    keyed.sort_unstable();
    for v in keyed {
        let pinched_top =
            (0..curved.faces.len()).find(|fi| corners.contains_key(&(*fi, pinched(*v))));
        let ends_at = |bi: &usize| {
            let e = &curved.edges[&bevels[*bi].edge];
            e.ends.0 == *v || e.ends.1 == *v
        };
        let at: Vec<usize> = (0..bevels.len()).filter(ends_at).collect();
        let rims: Vec<usize> = at
            .iter()
            .copied()
            .filter(|bi| matches!(bevels[*bi].surface, SurfaceGeometry::Torus(_)))
            .collect();
        let lines: Vec<usize> = at
            .iter()
            .copied()
            .filter(|bi| {
                matches!(curved.edges[&bevels[*bi].edge].curve, Curve::Line(_))
                    && matches!(bevels[*bi].surface, SurfaceGeometry::Cylinder(_))
            })
            .collect();
        let closed = |bi: &usize| {
            let e = &curved.edges[&bevels[*bi].edge];
            e.ends.0 == e.ends.1
        };
        // Two lobes where the top passes the vertex twice; one where it
        // ends there, a straight strip and an arc of the rim leaving it on
        // the one side.
        let (top, strips): (usize, Vec<usize>) = match (pinched_top, lines.as_slice()) {
            (Some(top), [l0, l1]) => {
                let one = matches!(rims.as_slice(), [r] if closed(r));
                let two = matches!(rims.as_slice(), [a, b] if !closed(a) && !closed(b));
                if !(one || two) || at.len() != 2 + rims.len() {
                    continue;
                }
                (top, vec![*l0, *l1])
            }
            (pinched_top, [l0]) => {
                let ([r], 2) = (rims.as_slice(), at.len()) else {
                    continue;
                };
                if closed(r) {
                    continue;
                }
                let (fl, fr) = (bevels[*l0].faces, bevels[*r].faces);
                let Some(top) = fl.iter().copied().find(|f| fr.contains(f)) else {
                    continue;
                };
                // The top's loop passes the vertex once, between these two;
                // or twice, a pass there with its edges left sharp, the top
                // running on past the vertex.
                let at_v = curved
                    .edges
                    .iter()
                    .filter(|(e, c)| {
                        (c.ends.0 == *v || c.ends.1 == *v) && faces_of(curved, **e).contains(&top)
                    })
                    .count();
                let passes = if pinched_top.is_some() { 2 } else { 1 };
                if pinched_top.is_some_and(|p| p != top)
                    || at_v != 2 * passes
                    || !corners.contains_key(&(top, *v))
                {
                    continue;
                }
                // A cusp: the top narrows to nothing, the two edges leaving
                // the vertex the same way. Leaving it opposite ways they
                // only continue each other (a stadium's side into its end).
                let leaving = |bi: usize| -> OgeomResult<Option<Vector>> {
                    let e = &curved.edges[&bevels[bi].edge];
                    let (t, sign) = if e.ends.0 == *v {
                        (e.range.0, 1.0)
                    } else {
                        (e.range.1, -1.0)
                    };
                    let d = e.curve.d1_at(t, tol)? * sign;
                    let m = d.magnitude();
                    Ok((m > tol.confusion()).then(|| d / m))
                };
                let (Some(dl), Some(dr)) = (leaving(*l0)?, leaving(*r)?) else {
                    continue;
                };
                if dl.dot(dr) < 1.0 - tol.angular() {
                    continue;
                }
                (top, vec![*l0])
            }
            _ => continue,
        };
        let other = |bi: usize| {
            let f = bevels[bi].faces;
            f.contains(&top)
                .then(|| if f[0] == top { f[1] } else { f[0] })
        };
        let (Some(wall), Some(side)) = (other(rims[0]), other(strips[0])) else {
            continue;
        };
        if strips.iter().any(|l| other(*l) != Some(side))
            || rims.iter().any(|r| other(*r) != Some(wall))
        {
            continue;
        }
        let crossings = [corners.get(&(top, *v)), corners.get(&(top, pinched(*v)))];
        let mut lobes = Vec::with_capacity(2);
        for bi in strips.iter().copied() {
            let e = bevels[bi].edge;
            let Some(crossing) = crossings.iter().flatten().find(|c| c.on(e).is_some()) else {
                break;
            };
            let t = crossing.point;
            let Some(rim) = rims
                .iter()
                .copied()
                .find(|r| crossing.on(bevels[*r].edge).is_some())
            else {
                break;
            };
            let r = &bevels[rim];
            let (SurfaceGeometry::Torus(torus), Curve::Circle(ring)) =
                (&r.surface, &r.contacts[usize::from(r.faces[1] == wall)])
            else {
                break;
            };
            let torus = torus.torus();
            let (o, z) = (torus.frame().origin(), torus.frame().z().vector());
            let (ring_at, ring_radius) = (ring.circle().frame().origin(), ring.circle().radius());
            let radial = {
                let d = t - o;
                (d - z * d.dot(z)).normalized(tol)?
            };
            let normal = z.cross(radial);
            let Curve::Line(line) = &bevels[bi].contacts[usize::from(bevels[bi].faces[1] == side)]
            else {
                break;
            };
            let axis = line.axis();
            let (at, along) = (axis.location, axis.direction.vector());
            if normal.dot(along).abs() <= tol.angular() {
                break;
            }
            let frame = Frame::new(
                t,
                Direction::new(normal, tol)?,
                Direction::new(z, tol)?,
                tol,
            )?;
            lobes.push(Lobe {
                line: bi,
                rim,
                top: t,
                side: at + along * (normal.dot(t - at) / normal.dot(along)),
                wall: ring_at + radial * ring_radius,
                plane: PlaneSurface::new(Plane::new(frame)).into(),
            });
        }
        if lobes.len() != strips.len() {
            continue;
        }
        // Each lobe's S, section by section square to the straight edge: at
        // `s` from the vertex the top stands open `w` from the side's line
        // out to the wall's foot along the radius, and the S is a convex arc
        // off the side turning `t` and a concave one onto the wall, point
        // symmetric, `cos t = 1 - u`, `u` that opening there over the one at
        // the lobe's end on that face: its foot on the side `sin t` of the lobe's depth there
        // below the top, on the wall `sin t` of its height there above. At
        // the vertex both feet stand on the top's edge; at the lobe each
        // runs into its strip's contact tangentially.
        let at = curved.points[v];
        let (Some(n_side), Some(n_wall), Some(n_top)) = (
            outward(curved, side, at),
            outward(curved, wall, at),
            outward(curved, top, at),
        ) else {
            continue;
        };
        if n_side.dot(n_wall) < 1.0 - tol.angular() {
            continue;
        }
        // Only a one-sided cusp's lobe takes its S's feet.
        let one_sided: &[Lobe] = if lobes.len() == 1 { &lobes } else { &[] };
        let mut feet = Vec::with_capacity(one_sided.len());
        for lobe in one_sided {
            let r = &bevels[lobe.rim];
            let Curve::Circle(rim_circle) = &curved.edges[&r.edge].curve else {
                break;
            };
            let (c, big) = (
                rim_circle.circle().frame().origin(),
                rim_circle.circle().radius(),
            );
            let along = {
                let d = lobe.side - at;
                (d - n_top * d.dot(n_top) - n_side * d.dot(n_side)).normalized(tol)?
            };
            // The top's width at `s`, from the side in along the top to the
            // wall's foot.
            let width = |s: f64| -> Option<f64> {
                let e = at + along * s;
                let d = e - c;
                let (b, k) = (-d.dot(n_side), d.dot(d) - big * big);
                let disc = b * b - k;
                (disc >= -tol.confusion()).then(|| {
                    let root = disc.max(0.0).sqrt();
                    let t = -b - root;
                    if t >= -tol.confusion() {
                        t.max(0.0)
                    } else {
                        -b + root
                    }
                })
            };
            let s_side = (lobe.side - at).dot(along);
            let on_top = lobe.wall - n_top * (lobe.wall - at).dot(n_top);
            let s_wall = (on_top - at).dot(along);
            // How far the top stands open at `s`: from the side's line
            // there out to the wall's foot, along the radius. Nought at the
            // vertex, growing as the square of `s` there, and read at every
            // station, past the foot's reach along the top too.
            let open = |s: f64| ((at + along * s) - c).magnitude() - big;
            let (w_side, w_wall) = (open(s_side), open(s_wall));
            if width(s_wall).is_none()
                || s_side <= tol.confusion()
                || s_wall <= tol.confusion()
                || w_side <= tol.confusion()
                || w_wall <= tol.confusion()
            {
                break;
            }
            let side_off = lobe.side - (at + along * s_side);
            let wall_up = n_top * (lobe.wall - on_top).dot(n_top);
            let lift = |s: f64, full: f64| -> f64 {
                let u = (open(s) / full).clamp(0.0, 1.0);
                (1.0 - u).acos().sin()
            };
            let count = 128_u32;
            let mut on_wall = Vec::with_capacity(count as usize + 1);
            let mut on_side = Vec::with_capacity(count as usize + 1);
            for k in 0..=count {
                let f = f64::from(k) / f64::from(count);
                // From the lobe to the vertex on the wall, from the vertex to
                // the lobe on the side.
                let sw = s_wall * (1.0 - f);
                let foot = at + along * sw - n_side * width(sw).unwrap_or(0.0);
                on_wall.push(foot + wall_up * lift(sw, w_wall));
                let ss = s_side * f;
                on_side.push(at + along * ss + side_off * lift(ss, w_side));
            }
            let n = on_wall.len() - 1;
            on_wall[0] = lobe.wall;
            on_wall[n] = at;
            on_side[0] = at;
            on_side[n] = lobe.side;
            let fit = crate::marched::band_fit_target(tol);
            let wall_foot = ogeom_geom::fit::fit_points(&on_wall, 3, 0.1 * fit, tol)?;
            let side_foot = ogeom_geom::fit::fit_points(&on_side, 3, 0.1 * fit, tol)?;
            feet.push((
                Curve::BSpline(wall_foot.curve),
                Curve::BSpline(side_foot.curve),
            ));
        }
        if feet.len() != one_sided.len() {
            continue;
        }
        let foot = (lobes.len() == 2)
            .then(|| {
                curved
                    .edges
                    .iter()
                    .find(|(e, c)| {
                        (c.ends.0 == *v || c.ends.1 == *v) && faces_of(curved, **e) == [wall]
                    })
                    .and_then(|(e, _)| corners.get(&(wall, *v)).filter(|c| c.on(*e).is_some()))
                    .map(|c| c.point)
            })
            .flatten();
        found.push(Cusp {
            v: *v,
            faces: [top, side, wall],
            lobes,
            feet,
            foot,
        });
    }
    Ok(found)
}

/// A rational quadratic arc's three control points, each with its weight.
type Conic = [(Point, f64); 3];

/// The arc of the circle about `centre` from `a` to `b` the shorter way,
/// as a rational quadratic: the two ends and where their tangents meet,
/// weighted by the cosine of half its turn.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where the arc turns half way round or more, or its ends meet.
pub(super) fn arc_control(
    centre: Point,
    a: Point,
    b: Point,
    tol: Tolerances,
) -> OgeomResult<Conic> {
    let (x, w) = ((a - centre).normalized(tol)?, (b - centre).normalized(tol)?);
    let radius = (a - centre).magnitude();
    let half = ((1.0 + x.dot(w).clamp(-1.0, 1.0)) / 2.0).sqrt();
    if half <= tol.angular() {
        ogeom_core::ogeom_bail!(Construction, "an arc of a half turn");
    }
    let apex = centre + (x + w).normalized(tol)? * (radius / half);
    Ok([(a, 1.0), (apex, half), (b, 1.0)])
}

/// A straight strip's end from `a` to `b`, its corners at stations
/// `along` apart on its axis (`centre` the axis's point at `a`'s): the
/// strip cut by the plane through both of normal `normal`, or where none is
/// given the one holding square to the axis the chord of the section at
/// `a`. That is the arc there carried along the axis onto the plane, an
/// affine map keeping a rational quadratic arc's weights, so exact.
///
/// # Errors
///
/// As [`arc_control`], and where the plane runs along the axis.
pub(super) fn skew_control(
    centre: Point,
    along: Vector,
    a: Point,
    b: Point,
    normal: Option<Vector>,
    tol: Tolerances,
) -> OgeomResult<Conic> {
    let back = b - along;
    let [_, (apex, half), _] = arc_control(centre, a, back, tol)?;
    let normal = normal.unwrap_or_else(|| (b - a).cross(along.cross(back - a)));
    if normal.dot(along).abs() <= tol.angular() * normal.magnitude() * along.magnitude() {
        ogeom_core::ogeom_bail!(Construction, "a strip cut along its own axis has no end");
    }
    let onto = |p: Point| p + along * (normal.dot(a - p) / normal.dot(along));
    Ok([(a, 1.0), (onto(apex), half), (b, 1.0)])
}

/// The rational quadratics `pieces`, each starting where the last ends, as
/// one curve, a knot doubled at each join, and where each piece ends on it:
/// the first over `(0, 1)`, each next over as long as keeps the curve's
/// speed through the join (the pieces meeting tangent there, the curve is
/// then smooth to first order in its parameter, as a patch fitted along it
/// needs).
///
/// # Errors
///
/// As making the curve.
pub(super) fn rational_quadratics(
    pieces: &[Conic],
    tol: Tolerances,
) -> OgeomResult<(Curve, Vec<f64>)> {
    let mut knots = vec![0.0; 3];
    let mut ends = Vec::with_capacity(pieces.len());
    let mut control = Vec::with_capacity(2 * pieces.len() + 1);
    let (mut at, mut speed) = (0.0, 0.0);
    for (k, piece) in pieces.iter().enumerate() {
        let [(p0, _), (p1, w1), (p2, _)] = *piece;
        let span = if k == 0 {
            1.0
        } else {
            2.0 * w1 * (p1 - p0).magnitude() / speed
        };
        if !(span.is_finite() && span > 0.0) {
            ogeom_core::ogeom_bail!(Construction, "a piece of a side without length");
        }
        for &(p, w) in &piece[usize::from(k > 0)..] {
            control.push(ogeom_math::Weighted::new(p, w, tol)?);
        }
        at += span;
        speed = 2.0 * w1 * (p2 - p1).magnitude() / span;
        knots.extend_from_slice(&[at, at]);
        if k + 1 == pieces.len() {
            knots.push(at);
        }
        ends.push(at);
    }
    let knots = ogeom_math::KnotVector::new(knots, 2)?;
    Ok((
        Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?),
        ends,
    ))
}

/// The points a torus strip's skew end is interpolated through.
const SKEW_SAMPLES: u32 = 128;

/// A torus strip's end from `a` to `b`, standing at two stations round its
/// axis: the curve running evenly in both of the torus's angles from one
/// to the other, the shorter way round in each, interpolated through
/// [`SKEW_SAMPLES`] of its points (each on the torus) at even parameters;
/// and its range.
///
/// # Errors
///
/// As the interpolation.
pub(super) fn torus_skew(
    torus: &ogeom_math::Torus,
    a: Point,
    b: Point,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64))> {
    use core::f64::consts::{PI, TAU};
    let frame = torus.frame();
    let (o, x, y, z) = (
        frame.origin(),
        frame.x().vector(),
        frame.y().vector(),
        frame.z().vector(),
    );
    let (big, small) = (torus.major_radius(), torus.minor_radius());
    let angles = |p: Point| {
        let d = p - o;
        let u = d.dot(y).atan2(d.dot(x));
        let rho = d.dot(x).hypot(d.dot(y));
        (u, d.dot(z).atan2(rho - big))
    };
    let shorter = |from: f64, to: f64| {
        let d = (to - from).rem_euclid(TAU);
        if d > PI { d - TAU } else { d }
    };
    let ((u0, v0), (u1, v1)) = (angles(a), angles(b));
    let (du, dv) = (shorter(u0, u1), shorter(v0, v1));
    let mut points = Vec::with_capacity(SKEW_SAMPLES as usize + 1);
    for k in 0..=SKEW_SAMPLES {
        let t = f64::from(k) / f64::from(SKEW_SAMPLES);
        let (u, v) = (u0 + du * t, v0 + dv * t);
        let rho = small.mul_add(v.cos(), big);
        points.push(o + (x * u.cos() + y * u.sin()) * rho + z * (small * v.sin()));
    }
    points[0] = a;
    points[SKEW_SAMPLES as usize] = b;
    let curve = ogeom_algo::interpolate(&points, 3, ogeom_algo::Spacing::Uniform, tol)?;
    let range = curve.domain();
    Ok((Curve::BSpline(curve), range))
}

/// A piece of a curve: the curve, its range, and whether the first point
/// asked for stands at the range's low end.
pub(super) type Met = (Curve, (f64, f64), bool);

/// The piece of where `surface` meets the plane `cut` between `p` and `q`,
/// the shorter way along it: its curve, its range, and whether `p` stands
/// at the range's low end. `None` where the two have no closed-form meeting
/// through both points.
pub(super) fn met_between(
    surface: &SurfaceGeometry,
    cut: &SurfaceGeometry,
    p: Point,
    q: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Met>> {
    let Ok(ogeom_intersect::surface::Meeting::Along(curves)) =
        ogeom_intersect::surface::surface_surface(surface, cut, tol)
    else {
        return Ok(None);
    };
    for curve in curves {
        let domain = curve.domain();
        if !(domain.0.is_finite() && domain.1.is_finite()) {
            continue;
        }
        let (tp, tq) = (
            param_on(&curve, domain, p, tol)?,
            param_on(&curve, domain, q, tol)?,
        );
        if curve.point_at(tp, tol)?.distance(p) > tol.confusion()
            || curve.point_at(tq, tol)?.distance(q) > tol.confusion()
        {
            continue;
        }
        let period = curve.is_periodic().then_some(domain.1 - domain.0);
        let (lo, hi, p_low) = if tp <= tq {
            (tp, tq, true)
        } else {
            (tq, tp, false)
        };
        return Ok(Some(match period {
            Some(period) if hi - lo > period / 2.0 => (curve, (hi, lo + period), !p_low),
            _ => (curve, (lo, hi), p_low),
        }));
    }
    Ok(None)
}

/// The point every one of `surfaces` passes through near `start`: Newton's
/// steps on the first three's signed distances (each surface's normal at
/// the foot of the point), the rest checked to pass there too. `None`
/// where the steps do not settle within the confusion distance, three
/// normals there span no corner, or another surface misses the point.
pub(super) fn shared_point(
    surfaces: &[&SurfaceGeometry],
    start: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    const STEPS: usize = 50;
    if surfaces.len() < 3 {
        return Ok(None);
    }
    let off = |surface: &SurfaceGeometry, x: Point| -> OgeomResult<(f64, Vector)> {
        let foot = ogeom_geom::project::project_on_surface(surface, x, 16, tol)?;
        let n = ogeom_geom::Surface::normal_at(surface, foot.parameters.0, foot.parameters.1, tol)?
            .vector();
        Ok(((x - foot.point).dot(n), n))
    };
    let mut x = start;
    for _ in 0..STEPS {
        let (fa, a) = off(surfaces[0], x)?;
        let (fb, b) = off(surfaces[1], x)?;
        let (fc, c) = off(surfaces[2], x)?;
        if fa.abs().max(fb.abs()).max(fc.abs()) <= 1e-3 * tol.confusion() {
            break;
        }
        let det = a.dot(b.cross(c));
        if det.abs() <= 1e-12 {
            return Ok(None);
        }
        x -= (b.cross(c) * fa + c.cross(a) * fb + a.cross(b) * fc) / det;
    }
    for surface in surfaces {
        if off(surface, x)?.0.abs() > tol.confusion() {
            return Ok(None);
        }
    }
    Ok(Some(x))
}

/// Whether the segment from `p` to `q` lies on every one of `surfaces`,
/// to the confusion distance at each sample.
pub(super) fn segment_on(
    surfaces: &[&SurfaceGeometry],
    p: Point,
    q: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    const SAMPLES: u32 = 4;
    for i in 1..SAMPLES {
        let at = p + (q - p) * (f64::from(i) / f64::from(SAMPLES));
        for surface in surfaces {
            let foot = ogeom_geom::project::project_on_surface(surface, at, 16, tol)?;
            if foot.point.distance(at) > tol.confusion() {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The most a sheet's points merge over: a fitted contact wider than this
/// keeps its corners apart.
pub(super) const MERGE_CAP: f64 = 1e-5;

/// How many points a closed ridge is sampled at to make it periodic.
pub(super) const RIDGE_SAMPLES: usize = 256;

/// A closed curve as a periodic cubic B-spline through `n` points evenly
/// spaced in its parameter, and how far it strays from the curve between
/// them. The uniform periodic cubic's point at its `j`-th knot is a sixth of
/// its controls `j`, `j + 2` and four of `j + 1`: solved by Jacobi's steps,
/// which the four on the diagonal against the two off it make converge.
pub(super) fn periodic_through(
    curve: &Curve,
    n: usize,
    tol: Tolerances,
) -> OgeomResult<(ogeom_geom::BSplineCurve, f64)> {
    let (lo, hi) = curve.domain();
    #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
    let at = |j: usize| lo + (hi - lo) * (j as f64) / (n as f64);
    let samples: Vec<Point> = (0..n)
        .map(|j| curve.point_at(at(j), tol))
        .collect::<OgeomResult<_>>()?;
    // Control `j + 1` sits under sample `j`.
    let mut control: Vec<Vector> = (0..n)
        .map(|j| samples[(j + n - 1) % n].to_vector())
        .collect();
    for _ in 0..200 {
        let mut worst: f64 = 0.0;
        let next: Vec<Vector> = (0..n)
            .map(|i| {
                let j = (i + n - 1) % n;
                let v = (samples[j].to_vector() * 6.0 - control[j] - control[(i + 1) % n]) / 4.0;
                worst = worst.max((v - control[i]).magnitude());
                v
            })
            .collect();
        control = next;
        if worst <= 1e-3 * tol.confusion() {
            break;
        }
    }
    let ring: Vec<Point> = control.into_iter().map(Point::from_vector).collect();
    let periodic = ogeom_geom::BSplineCurve::periodic(&ring, 3, tol)?;
    let (p0, _) = ogeom_geom::Curve3d::domain(&periodic);
    let mut off: f64 = 0.0;
    for j in 0..n {
        #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
        let u = p0 + j as f64 + 0.5;
        let p = ogeom_geom::Curve3d::point_at(&periodic, u, tol)?;
        let near = param_on(curve, (at(j), at(j + 1)), p, tol)?;
        off = off.max(curve.point_at(near, tol)?.distance(p));
    }
    Ok((periodic, off))
}
