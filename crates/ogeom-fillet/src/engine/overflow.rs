//! Blend into round (ADR 0001, amendment "overflow handling"): a round
//! whose ball, rolling on its floor, reaches an existing round of the
//! body rolls on over that round instead of trimming it; and a bevel whose
//! leg on the floor reaches such a round runs on over it.
//!
//! Built here for a boss's foot: the closed rim between a planar floor and
//! a cylindrical wall standing square on it, where the foot round's contact
//! circle crosses a straight edge of the floor along which a convex
//! cylindrical round meets the floor tangentially. Where the contact lies
//! on the floor the round is the torus its ball sweeps; past the tangent
//! line the ball rests on the wall and the existing round, and its centre
//! follows the intersection of their offsets (the cylinder about the wall's
//! axis and the one about the round's), so the two pieces meet G1 along the
//! section where the ball touches the tangent line. The fill between the
//! round, the wall, the floor and the existing round is one closed tool,
//! fused into the solid.
//!
//! A bevel of the same foot is the cone between its two legs over the
//! floor. Past the tangent line its leg on the floor is where the wall's
//! offset meets the existing round, and its leg on the wall where the
//! floor's offset meets the wall, so each section is the straight ruling
//! between them in the plane through the wall's axis; the two pieces share
//! that ruling where the leg crosses the tangent line, and meet G1 there
//! since the existing round leaves the floor tangentially.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{
    CircleCurve, ConeSurface, Curve, Curve3d as _, LineCurve, Surface as _, SurfaceGeometry,
    TorusSurface,
};
use ogeom_math::{Aabb, Circle, Cone, Direction, Frame, Point, Torus, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, SheetCurve, SheetRings, read};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// Stations along the piece over the existing round.
const STATIONS: usize = 128;
/// Points across each of its sections.
const ACROSS: usize = 65;

/// What the blend sweeps along the foot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Size {
    /// A round of this radius: the ball's arc between its contacts.
    Round(f64),
    /// A bevel, `out` along the floor and `rise` up the wall: the straight
    /// ruling between its legs.
    Bevel {
        /// From the wall along the floor.
        out: f64,
        /// Up the wall from the floor.
        rise: f64,
    },
}

impl Size {
    /// How far the blend reaches along the floor from the wall.
    const fn out(self) -> f64 {
        match self {
            Self::Round(r) => r,
            Self::Bevel { out, .. } => out,
        }
    }

    /// How far it reaches up the wall over the floor.
    const fn rise(self) -> f64 {
        match self {
            Self::Round(r) => r,
            Self::Bevel { rise, .. } => rise,
        }
    }
}

/// A foot blend reaching a straight round tangent to its floor.
#[derive(Debug, Clone)]
pub(crate) struct Reach {
    curved: Curved,
    edges: Vec<Id>,
    wall: usize,
    floor: usize,
    round: usize,
    line: Id,
    /// The foot circle's centre, on the floor.
    base: Point,
    /// The floor's outward normal.
    up: Vector,
    /// Square to the tangent line, from the boss toward it, in the floor.
    across: Vector,
    /// Along the tangent line: `up x across`.
    along: Vector,
    /// The wall's radius.
    wall_radius: f64,
    /// From the foot circle's centre to the tangent line.
    gap: f64,
    /// The existing round's radius.
    round_radius: f64,
    /// The foot blend's size.
    size: Size,
    /// How far it reaches along the floor: the round's radius, or the
    /// bevel's leg there.
    radius: f64,
}

impl Reach {
    fn radial(&self, theta: f64) -> Vector {
        self.across * theta.cos() + self.along * theta.sin()
    }

    /// The angle either side of `across` within which the ball rests on
    /// the existing round.
    fn half_angle(&self) -> f64 {
        (self.gap / (self.wall_radius + self.radius)).acos()
    }

    /// How far past the tangent line the ball's centre stands at `theta`.
    fn past(&self, theta: f64) -> f64 {
        (self.wall_radius + self.radius).mul_add(theta.cos(), -self.gap)
    }

    /// The ball centre's height above the floor at `theta`, or the height
    /// of the bevel's leg on the wall.
    fn height(&self, theta: f64) -> f64 {
        if let Size::Bevel { rise, .. } = self.size {
            return rise;
        }
        let s = self.past(theta);
        if s <= 0.0 {
            return self.radius;
        }
        let reach = self.round_radius + self.radius;
        (reach * reach - s * s).sqrt() - self.round_radius
    }

    fn centre(&self, theta: f64) -> Point {
        self.base
            + self.radial(theta) * (self.wall_radius + self.radius)
            + self.up * self.height(theta)
    }

    fn wall_touch(&self, theta: f64) -> Point {
        self.base + self.radial(theta) * self.wall_radius + self.up * self.height(theta)
    }

    /// Where the ball touches the floor or the existing round at `theta`,
    /// or where the bevel's leg there lies: on the wall's offset, straight
    /// below `on_line` on the floor or on the existing round.
    fn floor_touch(&self, theta: f64) -> Point {
        if let Size::Bevel { .. } = self.size {
            let s = self.past(theta);
            if s <= 0.0 {
                return self.on_line(theta);
            }
            let r = self.round_radius;
            return self.on_line(theta) + self.up * ((r * r - s * s).sqrt() - r);
        }
        let c = self.centre(theta);
        if self.past(theta) <= 0.0 {
            return c - self.up * self.radius;
        }
        let line = self.base + self.across * self.gap;
        let axis = line - self.up * self.round_radius + self.along * (c - line).dot(self.along);
        let w = c - axis;
        axis + w * (self.round_radius / w.magnitude())
    }

    fn foot(&self, theta: f64) -> Point {
        self.base + self.radial(theta) * self.wall_radius
    }

    /// The point where the ball's contact meets the tangent line at the
    /// seat change.
    fn on_line(&self, theta: f64) -> Point {
        self.base + self.radial(theta) * (self.wall_radius + self.radius)
    }

    /// The section at `theta`, from the wall's touch to the floor's: the
    /// ball's arc or the bevel's ruling, at `n` points.
    fn section(&self, theta: f64, n: usize, tol: Tolerances) -> Option<Vec<Point>> {
        let (w, f) = (self.wall_touch(theta), self.floor_touch(theta));
        #[expect(clippy::cast_precision_loss, reason = "a few points")]
        let at = |k: usize| k as f64 / (n - 1) as f64;
        if let Size::Bevel { .. } = self.size {
            return Some((0..n).map(|k| w + (f - w) * at(k)).collect());
        }
        let c = self.centre(theta);
        let a = (w - c) / self.radius;
        let b = (f - c) / self.radius;
        let cross = a.cross(b);
        let m = cross.magnitude();
        if m <= tol.angular() {
            return None;
        }
        let axis = cross / m;
        let sweep = a.dot(b).clamp(-1.0, 1.0).acos();
        Some(
            (0..n)
                .map(|k| {
                    let t = sweep * at(k);
                    c + (a * t.cos() + axis.cross(a) * t.sin()) * self.radius
                })
                .collect(),
        )
    }

    /// A point inside the blend's face at `theta`, and the fill's outward
    /// direction there: away from the corner the fill stands in.
    fn blend_mid(&self, theta: f64) -> (Point, Vector) {
        let (w, f) = (self.wall_touch(theta), self.floor_touch(theta));
        if let Size::Bevel { .. } = self.size {
            let d = f - w;
            let radial = self.radial(theta);
            let n = radial * -d.dot(self.up) + self.up * d.dot(radial);
            let n = n / n.magnitude();
            let mid = w + d * 0.5;
            let away = mid - self.foot(theta);
            return (mid, if n.dot(away) < 0.0 { -n } else { n });
        }
        let c = self.centre(theta);
        let m = (w - c) + (f - c);
        let m = m / m.magnitude();
        (c + m * self.radius, -m)
    }
}

/// The plane face and the cylindrical face `faces` holds, in that order.
fn plane_and_cylinder(curved: &Curved, faces: &[usize]) -> Option<(usize, usize)> {
    let [a, b] = faces[..] else {
        return None;
    };
    match (&curved.faces[a].surface, &curved.faces[b].surface) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Cylinder(_)) => Some((a, b)),
        (SurfaceGeometry::Cylinder(_), SurfaceGeometry::Plane(_)) => Some((b, a)),
        _ => None,
    }
}

/// `edges` read as a boss's foot rim blended at `size` whose blend
/// reaches one straight round tangent to the floor; `None` for any other
/// request, or where anything else of the solid stands within reach of the
/// fill.
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<Option<Reach>> {
    let eps = tol.confusion();
    let ang = tol.angular();
    let (radius, rise) = (size.out(), size.rise());
    if edges.is_empty() || [radius, rise].iter().any(|d| !d.is_finite() || *d <= eps) {
        return Ok(None);
    }
    // A plane and a cylinder either side of the first edge, before the
    // solid is read.
    let kinds: Vec<Option<SurfaceGeometry>> =
        ogeom_topo::ancestors_of(model, solid, &edges[0], ogeom_topo::ShapeType::Face)?
            .iter()
            .map(|f| {
                model
                    .node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface))
                    .cloned()
            })
            .collect();
    let plane = kinds
        .iter()
        .filter(|k| matches!(k, Some(SurfaceGeometry::Plane(_))))
        .count();
    let cylinder = kinds
        .iter()
        .filter(|k| matches!(k, Some(SurfaceGeometry::Cylinder(_))))
        .count();
    if kinds.len() != 2 || plane != 1 || cylinder != 1 {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut asked: Vec<Id> = edges.iter().map(id).collect();
    asked.sort_unstable();
    asked.dedup();
    if asked.iter().any(|e| !curved.edges.contains_key(e)) {
        return Ok(None);
    }
    let Some((floor, wall)) = plane_and_cylinder(&curved, &faces_of(&curved, asked[0])) else {
        return Ok(None);
    };
    // Every edge the floor and the wall share is asked, and nothing else.
    let mut shared: Vec<Id> = curved
        .edges
        .keys()
        .copied()
        .filter(|e| {
            let mut f = faces_of(&curved, *e);
            f.sort_unstable();
            let mut pair = vec![floor, wall];
            pair.sort_unstable();
            f == pair
        })
        .collect();
    shared.sort_unstable();
    if shared != asked {
        return Ok(None);
    }
    let SurfaceGeometry::Plane(plane) = &curved.faces[floor].surface else {
        return Ok(None);
    };
    let SurfaceGeometry::Cylinder(cylinder) = &curved.faces[wall].surface else {
        return Ok(None);
    };
    let plane = plane.plane();
    let up = plane.normal().vector()
        * if curved.faces[floor].reversed {
            -1.0
        } else {
            1.0
        };
    let axis = cylinder.cylinder().axis();
    if axis.direction.vector().cross(up).magnitude() > ang {
        return Ok(None);
    }
    let wall_radius = cylinder.cylinder().radius();
    let base = axis.location - up * (axis.location - plane.origin()).dot(up);
    // The foot circle closes round the axis on the floor.
    for e in &asked {
        let edge = &curved.edges[e];
        for k in 0..=4 {
            let t = edge.range.0 + (edge.range.1 - edge.range.0) * f64::from(k) / 4.0;
            let p = edge.curve.point_at(t, tol)?;
            let q = p - base;
            if q.dot(up).abs() > eps * 10.0 || (q.magnitude() - wall_radius).abs() > eps * 10.0 {
                return Ok(None);
            }
        }
    }
    // A boss: the wall's material outside its axis, standing above the
    // floor.
    let some_foot = curved.edges[&asked[0]]
        .curve
        .point_at(curved.edges[&asked[0]].range.0, tol)?;
    let out = outward(&curved, wall, some_foot).unwrap_or(Vector::ZERO);
    if out.dot(some_foot - base) <= 0.0 {
        return Ok(None);
    }
    let above = curved.faces[wall].loops.iter().flatten().any(|(e, _)| {
        let c = &curved.edges[e];
        [c.ends.0, c.ends.1]
            .iter()
            .any(|v| (curved.points[v] - base).dot(up) > rise)
    });
    if !above {
        return Ok(None);
    }
    // The floor's straight edges along which a convex round meets it
    // tangentially, those the contact circle crosses.
    let mut crossed: Vec<(Id, usize, Vector, Vector, f64, f64)> = Vec::new();
    for (e, _) in curved.faces[floor].loops.iter().flatten() {
        if asked.contains(e) {
            continue;
        }
        let edge = &curved.edges[e];
        let Curve::Line(line) = &edge.curve else {
            continue;
        };
        let others: Vec<usize> = faces_of(&curved, *e)
            .into_iter()
            .filter(|f| *f != floor)
            .collect();
        let [round] = others[..] else {
            continue;
        };
        let SurfaceGeometry::Cylinder(c) = &curved.faces[round].surface else {
            continue;
        };
        let along = line.axis().direction.vector();
        let round_axis = c.cylinder().axis();
        if round_axis.direction.vector().cross(along).magnitude() > ang {
            continue;
        }
        let round_radius = c.cylinder().radius();
        let p = line.axis().location;
        let off = p - round_axis.location;
        let off = off - along * off.dot(along);
        if (off - up * round_radius).magnitude() > eps * 10.0 {
            continue;
        }
        if outward(&curved, round, p).is_none_or(|n| n.dot(up) < 1.0 - ang) {
            continue;
        }
        let to = p - base;
        let to = to - up * to.dot(up);
        let to = to - along * to.dot(along);
        let gap = to.magnitude();
        if gap <= wall_radius + eps * 10.0 {
            continue;
        }
        // A bevel's leg stays on the existing round, short of where it
        // turns down from the floor.
        let on_round = match size {
            Size::Round(_) => true,
            Size::Bevel { .. } => wall_radius + radius - gap < round_radius - eps * 10.0,
        };
        if gap < wall_radius + radius - eps * 10.0 && on_round {
            let across = to / gap;
            crossed.push((*e, round, across, up.cross(across), gap, round_radius));
        }
    }
    let [(line, round, across, along, gap, round_radius)] = crossed[..] else {
        return Ok(None);
    };
    let reach = Reach {
        curved,
        edges: asked,
        wall,
        floor,
        round,
        line,
        base,
        up,
        across,
        along,
        wall_radius,
        gap,
        round_radius,
        size,
        radius,
    };
    if !clear(model, &reach, tol)? {
        return Ok(None);
    }
    Ok(Some(reach))
}

/// Whether nothing of the solid but the wall, the floor and the existing
/// round stands within the box the fill can take up: within the cylinder
/// of the contacts' greatest distance from the axis, between the lowest
/// contact and the ball's top.
fn clear(model: &Model, reach: &Reach, tol: Tolerances) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let half = reach.half_angle();
    let mut far = reach.wall_radius + reach.radius;
    let mut low = 0.0_f64;
    for k in 0..=64 {
        let theta = half * (2.0 * f64::from(k) / 64.0 - 1.0);
        let t = reach.floor_touch(theta) - reach.base;
        low = low.min(t.dot(reach.up));
        far = far.max((t - reach.up * t.dot(reach.up)).magnitude());
    }
    let margin = reach.radius * 1e-3 + eps * 10.0;
    let far = far + margin;
    let (lo, hi) = (low - margin, reach.size.rise().max(reach.radius) + margin);
    let mut corners = Vec::with_capacity(16);
    for h in [lo, hi] {
        for k in 0..8 {
            let a = core::f64::consts::FRAC_PI_4 * f64::from(k);
            // The octagon round the circle of radius `far`.
            let r = far / core::f64::consts::FRAC_PI_8.cos();
            corners.push(reach.base + reach.radial(a) * r + reach.up * h);
        }
    }
    let room = Aabb::of_points(&corners);
    let c = &reach.curved;
    for (fi, face) in c.faces.iter().enumerate() {
        if fi == reach.wall || fi == reach.floor || fi == reach.round {
            continue;
        }
        if ogeom_algo::shape_bounds(model, &face.shape, tol)?.intersects(&room) {
            return Ok(false);
        }
    }
    for fi in [reach.wall, reach.floor, reach.round] {
        for (e, _) in c.faces[fi].loops.iter().flatten() {
            if reach.edges.contains(e) || *e == reach.line {
                continue;
            }
            // A seam of the wall is the wall's own.
            if faces_of(c, *e) == [reach.wall] {
                continue;
            }
            if ogeom_algo::shape_bounds(model, &c.edges[e].shape, tol)?.intersects(&room) {
                return Ok(false);
            }
        }
    }
    let line = &c.edges[&reach.line];
    Ok(!room.contains(c.points[&line.ends.0]) && !room.contains(c.points[&line.ends.1]))
}

/// The arc of `radius` about `centre` in the plane square to `z`, from the
/// direction `x` through `sweep`.
fn arc(
    centre: Point,
    z: Vector,
    x: Vector,
    radius: f64,
    sweep: f64,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64))> {
    let frame = Frame::new(
        centre,
        Direction::new(z, tol)?,
        Direction::new(x, tol)?,
        tol,
    )?;
    Ok((
        CircleCurve::new(Circle::new(frame, radius, tol)?).into(),
        (0.0, sweep),
    ))
}

fn segment(a: Point, b: Point, tol: Tolerances) -> OgeomResult<(Curve, (f64, f64))> {
    Ok((
        Curve::Line(LineCurve::segment(a, b, tol)?),
        (0.0, a.distance(b)),
    ))
}

/// The piece of the blend over the existing round: its section at each
/// station from its wall touch to its touch on the round, fitted.
fn strip(reach: &Reach, tol: Tolerances) -> OgeomResult<Option<ogeom_geom::BSplineSurface>> {
    let half = reach.half_angle();
    let mut stations: Vec<Vec<Point>> = Vec::with_capacity(STATIONS + 1);
    for i in 0..=STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
        let theta = half * (2.0 * i as f64 / STATIONS as f64 - 1.0);
        let Some(row) = reach.section(theta, ACROSS, tol) else {
            return Ok(None);
        };
        stations.push(row);
    }
    let target = tol.confusion() * 0.5;
    let Ok(fitted) = ogeom_geom::fit::fit_surface_grid_chordal(&stations, 3, target, tol) else {
        return Ok(None);
    };
    Ok((fitted.error <= tol.confusion()).then_some(fitted.curve))
}

/// The boundary iso-curve of `surface` between `points[ends.0]` and
/// `points[ends.1]`, as a sheet curve over those ends.
fn boundary_from(
    surface: &ogeom_geom::BSplineSurface,
    points: &[Point],
    ends: (usize, usize),
    tol: Tolerances,
) -> OgeomResult<Option<SheetCurve>> {
    let (from, to) = (points[ends.0], points[ends.1]);
    let ((u0, u1), (v0, v1)) = surface.domain();
    let mut candidates = Vec::with_capacity(4);
    for u in [u0, u1] {
        candidates.push(surface.iso_u_curve(u, tol)?);
    }
    for v in [v0, v1] {
        candidates.push(surface.iso_v_curve(v, tol)?);
    }
    let reach = tol.confusion() * 10.0;
    for c in candidates {
        let (a, b) = c.domain();
        let (pa, pb) = (c.point_at(a, tol)?, c.point_at(b, tol)?);
        if pa.distance(from) <= reach && pb.distance(to) <= reach {
            return Ok(Some(SheetCurve {
                ends,
                curve: c.into(),
                range: (a, b),
            }));
        }
        if pb.distance(from) <= reach && pa.distance(to) <= reach {
            return Ok(Some(SheetCurve {
                ends: (ends.1, ends.0),
                curve: c.into(),
                range: (a, b),
            }));
        }
    }
    Ok(None)
}

/// The parameters' periods of `surface`, where they close.
fn periods(surface: &SurfaceGeometry) -> (Option<f64>, Option<f64>) {
    let turn = core::f64::consts::TAU;
    match surface {
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => (Some(turn), None),
        SurfaceGeometry::Torus(_) => (Some(turn), Some(turn)),
        SurfaceGeometry::Sphere(_) => (Some(turn), None),
        _ => (None, None),
    }
}

/// One face of a tool: its surface, its loop by curve, and a point of it
/// with the tool's outward direction there.
struct ToolFace {
    surface: SurfaceGeometry,
    curves: Vec<usize>,
    at: Point,
    out: Vector,
}

/// The closed solid `faces` bound over `points` and `curves`, each loop
/// walked so its face lies on its left about the surface's normal (by the
/// loop's signed area in the surface's parameters), each face turned so
/// it faces `out`. `None` where a loop will not close, the shell is open,
/// or the solid holds no volume.
fn closed(
    model: &mut Model,
    points: &[Point],
    curves: &[SheetCurve],
    faces: &[ToolFace],
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let mut sheet: Vec<SheetRings> = Vec::with_capacity(faces.len());
    for face in faces {
        let first = face.curves[0];
        let mut walk = vec![(first, true)];
        let mut at = curves[first].ends.1;
        let mut left: Vec<usize> = face.curves[1..].to_vec();
        while !left.is_empty() {
            let Some(k) = left
                .iter()
                .position(|c| curves[*c].ends.0 == at || curves[*c].ends.1 == at)
            else {
                return Ok(None);
            };
            let c = left.remove(k);
            let forward = curves[c].ends.0 == at;
            at = if forward {
                curves[c].ends.1
            } else {
                curves[c].ends.0
            };
            walk.push((c, forward));
        }
        if at != curves[first].ends.0 {
            return Ok(None);
        }
        // The loop's signed area in the surface's parameters, each period
        // unwrapped along the walk.
        let (pu, pv) = periods(&face.surface);
        let mut uv: Vec<(f64, f64)> = Vec::new();
        for (c, forward) in &walk {
            let (lo, hi) = curves[*c].range;
            for k in 0..24 {
                let s = f64::from(k) / 24.0;
                let t = if *forward {
                    lo + (hi - lo) * s
                } else {
                    hi - (hi - lo) * s
                };
                let p = curves[*c].curve.point_at(t, tol)?;
                let foot = ogeom_geom::project::project_on_surface(&face.surface, p, 16, tol)?;
                let (mut u, mut v) = foot.parameters;
                if let Some(&(lu, lv)) = uv.last() {
                    if let Some(period) = pu {
                        u += period * ((lu - u) / period).round();
                    }
                    if let Some(period) = pv {
                        v += period * ((lv - v) / period).round();
                    }
                }
                uv.push((u, v));
            }
        }
        let n = uv.len();
        let area: f64 = (0..n)
            .map(|i| {
                let (a, b) = (uv[i], uv[(i + 1) % n]);
                a.0 * b.1 - b.0 * a.1
            })
            .sum();
        if area.abs() <= f64::EPSILON {
            return Ok(None);
        }
        if area < 0.0 {
            walk = walk.iter().rev().map(|(c, f)| (*c, !f)).collect();
        }
        let foot = ogeom_geom::project::project_on_surface(&face.surface, face.at, 16, tol)?;
        let (u, v) = foot.parameters;
        let (du, dv) = face.surface.d1_at(u, v, tol)?;
        sheet.push(SheetRings {
            surface: face.surface.clone(),
            rings: vec![walk],
            reversed: du.cross(dv).dot(face.out) < 0.0,
        });
    }
    let Ok(built) = super::curved::build(model, points, curves, &sheet, tol) else {
        return Ok(None);
    };
    let shell = ogeom_algo::make_shell(model, &built.faces)?.shape;
    if !ogeom_algo::is_shell_closed(model, &shell)? {
        return Ok(None);
    }
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    let volume =
        ogeom_algo::volume_properties(model, &solid, ogeom_mesh::Deflection::default(), tol)?.mass;
    if volume <= 0.0 {
        return Ok(None);
    }
    Ok(Some((solid, built.faces)))
}

/// The blend's face over the floor: the torus the ball sweeps, or the cone
/// between the bevel's legs.
fn over_floor(reach: &Reach, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let (r, rr, up) = (reach.radius, reach.wall_radius, reach.up);
    let Size::Bevel { rise, .. } = reach.size else {
        return Ok(TorusSurface::new(Torus::new(
            Frame::new(
                reach.base + up * r,
                Direction::new(up, tol)?,
                Direction::new(reach.across, tol)?,
                tol,
            )?,
            rr + r,
            r,
            tol,
        )?)
        .into());
    };
    // The cone's radius falls from `rr + r` on the floor to `rr` at the
    // leg on the wall; its apex stands past that leg.
    let frame = Frame::new(
        reach.base,
        Direction::new(up, tol)?,
        Direction::new(reach.across, tol)?,
        tol,
    )?;
    let apex = rise * (rr + r) / r;
    let pad = (0.1 * rise).min(0.5 * (apex - rise));
    let cone = Cone::new(frame, rr + r, (-r / rise).atan(), tol)?;
    Ok(ConeSurface::new(cone, (-pad, rise + pad))?.into())
}

/// The fill of `reach` as one closed tool, with its two blend faces.
fn tool(
    model: &mut Model,
    reach: &Reach,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, [Shape; 2])>> {
    let r = reach.radius;
    let rise = reach.size.rise();
    let rr = reach.wall_radius;
    let half = reach.half_angle();
    let (t0, t1) = (-half, half);
    let up = reach.up;
    let Some(band) = strip(reach, tol)? else {
        return Ok(None);
    };
    // The points at each seat change: the foot, the wall's touch, the
    // touch on the tangent line.
    let points = vec![
        reach.foot(t0),
        reach.wall_touch(t0),
        reach.on_line(t0),
        reach.foot(t1),
        reach.wall_touch(t1),
        reach.on_line(t1),
    ];
    let (f0, w0, p0, f1, w1, p1) = (0, 1, 2, 3, 4, 5);
    let around = core::f64::consts::TAU - 2.0 * half;
    let mut curves: Vec<SheetCurve> = Vec::new();
    let mut push = |ends: (usize, usize), (curve, range): (Curve, (f64, f64))| {
        curves.push(SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    let x1 = reach.radial(t1);
    // Round the floor: the foot, the wall's touch and the floor's touch.
    let foot_a = push((f1, f0), arc(reach.base, up, x1, rr, around, tol)?);
    let top_a = push(
        (w1, w0),
        arc(reach.base + up * rise, up, x1, rr, around, tol)?,
    );
    let out_a = push((p1, p0), arc(reach.base, up, x1, rr + r, around, tol)?);
    // Over the existing round.
    let x0 = reach.radial(t0);
    let foot_b = push((f0, f1), arc(reach.base, up, x0, rr, 2.0 * half, tol)?);
    let Some(wall_side) = boundary_from(&band, &points, (w0, w1), tol)? else {
        return Ok(None);
    };
    let wall_b = push(wall_side.ends, (wall_side.curve, wall_side.range));
    let Some(round_side) = boundary_from(&band, &points, (p0, p1), tol)? else {
        return Ok(None);
    };
    let round_b = push(round_side.ends, (round_side.curve, round_side.range));
    let line_b = push((p0, p1), segment(points[p0], points[p1], tol)?);
    // At each seat change: up the wall, out along the floor, and the
    // section between them.
    let mut ups = [0; 2];
    let mut rays = [0; 2];
    let mut quarters = [0; 2];
    for (k, (theta, f, w, p)) in [(t0, f0, w0, p0), (t1, f1, w1, p1)].into_iter().enumerate() {
        let radial = reach.radial(theta);
        ups[k] = push((f, w), segment(points[f], points[w], tol)?);
        rays[k] = push((f, p), segment(points[f], points[p], tol)?);
        let section = match reach.size {
            Size::Round(_) => arc(
                reach.centre(theta),
                radial.cross(up),
                -radial,
                r,
                core::f64::consts::FRAC_PI_2,
                tol,
            )?,
            Size::Bevel { .. } => segment(points[w], points[p], tol)?,
        };
        quarters[k] = push((w, p), section);
    }
    let curved = &reach.curved;
    let wall = curved.faces[reach.wall].surface.clone();
    let floor = curved.faces[reach.floor].surface.clone();
    let round = curved.faces[reach.round].surface.clone();
    let opposite = core::f64::consts::PI;
    let blend = over_floor(reach, tol)?;
    // A point of each face and the tool's outward direction there.
    let (blend_at, blend_out) = reach.blend_mid(opposite);
    let (band_at, band_out) = reach.blend_mid(0.0);
    let line_at = reach.base + reach.across * reach.gap;
    let round_axis = line_at - up * reach.round_radius;
    let on_round = {
        let t = reach.floor_touch(0.0) - round_axis;
        let l = line_at - round_axis;
        let m = t / t.magnitude() + l / l.magnitude();
        m / m.magnitude()
    };
    let faces = vec![
        ToolFace {
            surface: wall.clone(),
            curves: vec![foot_a, ups[0], top_a, ups[1]],
            at: reach.foot(opposite) + up * (rise * 0.5),
            out: -reach.radial(opposite),
        },
        ToolFace {
            surface: floor.clone(),
            curves: vec![foot_a, rays[0], out_a, rays[1]],
            at: reach.base + reach.radial(opposite) * (rr + r * 0.5),
            out: -up,
        },
        ToolFace {
            surface: blend,
            curves: vec![top_a, quarters[0], out_a, quarters[1]],
            at: blend_at,
            out: blend_out,
        },
        ToolFace {
            surface: wall,
            curves: vec![foot_b, ups[1], wall_b, ups[0]],
            at: reach.foot(0.0) + up * (reach.height(0.0) * 0.5),
            out: -reach.across,
        },
        ToolFace {
            surface: floor,
            curves: vec![foot_b, rays[1], line_b, rays[0]],
            at: reach.base + reach.across * (0.5 * (rr + reach.gap)),
            out: -up,
        },
        ToolFace {
            surface: round,
            curves: vec![round_b, line_b],
            at: round_axis + on_round * reach.round_radius,
            out: -on_round,
        },
        ToolFace {
            surface: band.into(),
            curves: vec![wall_b, quarters[1], round_b, quarters[0]],
            at: band_at,
            out: band_out,
        },
    ];
    let Some((solid, made)) = closed(model, &points, &curves, &faces, tol)? else {
        return Ok(None);
    };
    Ok(Some((solid, [made[2].clone(), made[6].clone()])))
}

/// The foot blend of `edges` at `size` blended into the round it
/// reaches, fused into `solid` as one tool; `None` where `reach` does not
/// hold the request or the tool will not build, apply, or prove.
///
/// # Errors
///
/// As reading the solid and the booleans, other than their refusals.
pub(crate) fn blended(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(reach) = reach(model, solid, edges, size, tol)? else {
        return Ok(None);
    };
    let Some((fill, [over, band])) = tool(model, &reach, tol)? else {
        return Ok(None);
    };
    let deflection = ogeom_mesh::Deflection::default();
    let before = ogeom_algo::volume_properties(model, solid, deflection, tol)?.mass;
    let added = ogeom_algo::volume_properties(model, &fill, deflection, tol)?.mass;
    let mut named = Vec::with_capacity(2 * edges.len());
    for e in edges {
        named.push((over.clone(), e.clone()));
        named.push((band.clone(), e.clone()));
    }
    let tools = [super::fallback::Tool {
        solid: fill,
        sense: super::fallback::Sense::Fuse,
        faces: named,
    }];
    let Some(built) = super::fallback::applied(model, solid, edges, &tools, tol)? else {
        return Ok(None);
    };
    let built = super::merged(model, solid, built, tol)?;
    // The fill stands outside the solid, touching it only on the faces it
    // was closed on: the fuse adds its whole volume.
    let after = ogeom_algo::volume_properties(model, &built.shape, deflection, tol)?.mass;
    if (after - before - added).abs() > 1e-6 * (1.0 + before.abs()) {
        return Ok(None);
    }
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    Ok(Some(built))
}

/// Whether `edges` at `size` are a foot blend that blend into round
/// takes ([`reach`]).
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn reach_found(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<bool> {
    Ok(reach(model, solid, edges, size, tol)?.is_some())
}

/// The edges of `specs` and the one bevel they ask, where every spec asks
/// the same legs on a rim between a cylinder's wall and a plane square to
/// it; `None` for any other request.
///
/// # Errors
///
/// As reading the seats, other than a spec's own refusal, which the
/// engine's construction reports.
pub(crate) fn bevel_size(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, crate::Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<Shape>, Size)>> {
    let eps = tol.confusion();
    let mut legs: Option<(f64, f64)> = None;
    for (edge, spec) in specs {
        let Ok(crate::chamfer::Legs::Rim(seat, on_wall, on_cap)) =
            crate::chamfer::legs(model, solid, edge, spec, tol)
        else {
            return Ok(None);
        };
        if seat.slope.abs() > tol.angular() {
            return Ok(None);
        }
        match legs {
            None => legs = Some((on_wall, on_cap)),
            Some((w, c)) if (w - on_wall).abs() <= eps && (c - on_cap).abs() <= eps => {}
            Some(_) => return Ok(None),
        }
    }
    let Some((rise, out)) = legs else {
        return Ok(None);
    };
    let edges = specs.iter().map(|(e, _)| e.clone()).collect();
    Ok(Some((edges, Size::Bevel { out, rise })))
}

/// Whether `built`, the rounds of `edges` on `solid`, changed a face of
/// `solid` other than the edges' own and those at their ends: a face the
/// rounds reached and trimmed through. A face is unchanged where the result
/// holds a face of its area and its bounds on whose surface a point of it
/// lies.
///
/// # Errors
///
/// As exploring the shapes and measuring the faces.
pub(crate) fn changes_others(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    built: &Built,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_topo::{ShapeType, ancestors_of, explore_unique};
    let eps = tol.confusion() * 10.0;
    let deflection = ogeom_mesh::Deflection::default();
    let mut own: Vec<Shape> = Vec::new();
    for edge in edges {
        own.extend(ancestors_of(model, solid, edge, ShapeType::Face)?);
        for vertex in explore_unique(model, edge, ShapeType::Vertex)? {
            own.extend(ancestors_of(model, solid, &vertex, ShapeType::Face)?);
        }
    }
    let read = |face: &Shape| -> OgeomResult<(f64, Aabb)> {
        Ok((
            ogeom_algo::surface_properties(model, face, deflection, tol)?.mass,
            ogeom_algo::shape_bounds(model, face, tol)?,
        ))
    };
    let same_box = |a: &Aabb, b: &Aabb| match (a.low(), a.high(), b.low(), b.high()) {
        (Some(al), Some(ah), Some(bl), Some(bh)) => al.distance(bl).max(ah.distance(bh)) <= eps,
        _ => false,
    };
    let made: Vec<(Shape, f64, Aabb)> = explore_unique(model, &built.shape, ShapeType::Face)?
        .into_iter()
        .map(|f| read(&f).map(|(area, bounds)| (f, area, bounds)))
        .collect::<OgeomResult<_>>()?;
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if own.iter().any(|o| o.is_same(&face)) {
            continue;
        }
        let (area, bounds) = read(&face)?;
        let (p, _) = ogeom_algo::face_normal(model, &face, tol)?;
        let mut kept = false;
        for (g, now, around) in &made {
            if (area - now).abs() > 1e-9 * (1.0 + area) || !same_box(&bounds, around) {
                continue;
            }
            let Some(surface) = model
                .node(g)
                .and_then(|x| x.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
            else {
                continue;
            };
            let on = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?.distance;
            if on <= eps {
                kept = true;
                break;
            }
        }
        if !kept {
            return Ok(true);
        }
    }
    Ok(false)
}
