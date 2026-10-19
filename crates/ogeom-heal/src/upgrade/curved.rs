//! Unification on curved carriers: cylinders, cones, spheres and tori on
//! one axis with one size, and extrusions of one basis curve (a line, a
//! conic or a B-spline) along one line.
//!
//! Two charts of one such carrier differ by a map that is affine in each
//! parameter alone: the angle about the axis turned and possibly run the
//! other way, the height (or latitude, or tube angle, or sweep) shifted and
//! possibly run the other way. The map is read off the two frames exactly
//! and then measured over the face it carries; a pair the frames do not
//! relate so is not one carrier here, however close the surfaces pass.

use core::ops::Range;

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction2, Frame, Point, Point2, Transform, Transform2, Vector, Vector2};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, SurfaceId, TShapeId};

/// How many points each pcurve is read at, ends included.
const SAMPLES: u32 = 17;

/// Where a curved face lies: its surface as stored, what places that in
/// the world, the side the face takes, the box its boundary spans in the
/// stored chart, and whether its boundary walks a seam.
#[derive(Debug, Clone)]
pub(super) struct Curved {
    pub(super) id: SurfaceId,
    surface: SurfaceGeometry,
    placement: Transform,
    reversed: bool,
    window: (Point2, Point2),
    seamed: bool,
}

/// The carrier of a face on a cylinder, cone, sphere, torus or extrusion,
/// or `None` for any other face, a face placed by a scaling, or one whose
/// boundary is not all pcurves or seams on its surface.
pub(super) fn carrier_of(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Curved>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        return Ok(None);
    };
    if !matches!(
        surface,
        SurfaceGeometry::Cylinder(_)
            | SurfaceGeometry::Cone(_)
            | SurfaceGeometry::Sphere(_)
            | SurfaceGeometry::Torus(_)
            | SurfaceGeometry::Extrusion(_)
    ) {
        return Ok(None);
    }
    let placement = face.transform(model.datums())?;
    if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
        return Ok(None);
    }
    let reversed = face.orientation() == Orientation::Reversed;
    let mut lo = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut seamed = false;
    for wire in model.children_of(face)? {
        for edge in model.ordered_children_of(&wire)? {
            let Some((pcurve, range, seam)) = use_pcurve(model, &edge, data.surface, reversed)
            else {
                return Ok(None);
            };
            seamed |= seam;
            for p in samples(&pcurve, range, false, tol)? {
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
    }
    if !(lo.x.is_finite() && lo.y.is_finite() && hi.x.is_finite() && hi.y.is_finite()) {
        return Ok(None);
    }
    Ok(Some(Curved {
        id: data.surface,
        surface,
        placement,
        reversed,
        window: (lo, hi),
        seamed,
    }))
}

/// Whether two curved carriers are one surface facing one way.
pub(super) fn same_carrier(a: &Curved, b: &Curved, tol: Tolerances) -> bool {
    chart_map(a, b, tol).is_some_and(|map| faces_one_way(a, b, &map))
}

/// [`same_carrier`], for two faces neither of which walks a seam, related
/// by a plain map (see [`chart_map_of`]): the pairs a cluster falls back to
/// where the whole cluster does not merge.
pub(super) fn same_plain_carrier(a: &Curved, b: &Curved, tol: Tolerances) -> bool {
    !a.seamed
        && !b.seamed
        && chart_map_of(a, b, tol).is_some_and(|(map, plain)| plain && faces_one_way(a, b, &map))
}

/// Whether two faces related by `map` take the same side of their
/// carrier. Each face's outward normal is its placement applied to its
/// stored chart's `Su × Sv`, turned back where the face is used reversed;
/// in placed tangents that is `Su × Sv` turned again by a mirroring
/// placement, and the placed tangents of one chart are those of the other
/// times the map's Jacobian, whose determinant is the sign of its flips.
fn faces_one_way(a: &Curved, b: &Curved, map: &ChartMap) -> bool {
    let turned = |c: &Curved| c.reversed != !c.placement.preserves_handedness();
    (turned(a) != (map.flip.0 != map.flip.1)) == turned(b)
}

/// A map between two charts of one carrier: each parameter run one way or
/// the other, then shifted.
#[derive(Debug, Clone, Copy)]
struct ChartMap {
    flip: (bool, bool),
    shift: Vector2,
}

impl ChartMap {
    const IDENTITY: Self = Self {
        flip: (false, false),
        shift: Vector2::ZERO,
    };

    fn apply(&self, p: Point2) -> Point2 {
        Point2::new(
            if self.flip.0 { -p.x } else { p.x } + self.shift.x,
            if self.flip.1 { -p.y } else { p.y } + self.shift.y,
        )
    }

    /// The same map as a planar transform, built from exact reflections in
    /// the chart's axes so a pcurve carried by it lands where `apply` puts
    /// its points.
    fn transform(&self) -> Transform2 {
        let mut map = Transform2::translation(self.shift);
        if self.flip.0 {
            map = map * Transform2::line_mirror(Point2::ORIGIN, Direction2::X);
        }
        if self.flip.1 {
            map = map * Transform2::line_mirror(Point2::ORIGIN, Direction2::Y);
        }
        map
    }
}

/// A stored frame's origin and axes carried into the world, as vectors: a
/// mirroring placement leaves them left-handed, which the chart formulas
/// read the same way.
struct Axes {
    origin: Point,
    x: Vector,
    y: Vector,
    z: Vector,
}

fn axes(frame: &Frame, placement: &Transform) -> Axes {
    Axes {
        origin: placement.apply(frame.origin()),
        x: placement.apply_vector(frame.x().vector()),
        y: placement.apply_vector(frame.y().vector()),
        z: placement.apply_vector(frame.z().vector()),
    }
}

/// The angle map between two frames sharing an axis line, for charts
/// whose angle `u` names the direction `cos u·x + sin u·y`: `from`'s angle
/// `u` is `to`'s `±u + turn`, run backwards where the two frames see the
/// turn about the axis opposite ways. Returns that flip, the turn, and
/// whether the axes run opposite ways. `None` where the axes are not
/// parallel.
fn about_one_axis(from: &Axes, to: &Axes, tol: Tolerances) -> Option<(bool, f64, bool)> {
    let along = from.z.dot(to.z);
    if (along.abs() - 1.0).abs() > tol.angular() {
        return None;
    }
    let turn = from.x.dot(to.y).atan2(from.x.dot(to.x));
    let backwards = from.x.cross(from.y).dot(to.x.cross(to.y)) < 0.0;
    Some((backwards, turn, along < 0.0))
}

/// How far `offset` runs along `axis`, where it runs along it alone.
fn along_axis(offset: Vector, axis: Vector, tol: Tolerances) -> Option<f64> {
    let run = offset.dot(axis);
    ((offset - axis * run).magnitude() <= tol.confusion()).then_some(run)
}

/// The map taking `from`'s chart onto `to`'s, where the two are one
/// carrier, and `None` otherwise.
fn chart_map(from: &Curved, to: &Curved, tol: Tolerances) -> Option<ChartMap> {
    chart_map_of(from, to, tol).map(|(map, _)| map)
}

/// [`chart_map`], and whether the map is plain: two extrusions whose basis
/// curves are one moved along the sweep with their parameter kept, or any
/// other pair of carriers. A basis run backwards, shifted along itself or
/// of another kind than a line, circle or B-spline is not plain.
fn chart_map_of(from: &Curved, to: &Curved, tol: Tolerances) -> Option<(ChartMap, bool)> {
    use SurfaceGeometry as S;
    let close = |a: f64, b: f64| (a - b).abs() <= tol.confusion();
    let mut plain = true;
    let map = match (&from.surface, &to.surface) {
        (S::Cylinder(a), S::Cylinder(b)) => {
            let (a, b) = (a.cylinder(), b.cylinder());
            if !close(a.radius(), b.radius()) {
                return None;
            }
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            let rise = along_axis(fa.origin - fb.origin, fb.z, tol)?;
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, rise),
            }
        }
        (S::Cone(a), S::Cone(b)) => {
            let (a, b) = (a.cone(), b.cone());
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            let rise = along_axis(fa.origin - fb.origin, fb.z, tol)?;
            // The radius `r + h·tan(a)` at `from`'s height `h` is the
            // radius at `to`'s height `±h + rise`.
            let (slope_a, slope_b) = (a.half_angle().tan(), b.half_angle().tan());
            let slope_b_seen = if opposed { -slope_b } else { slope_b };
            if (slope_a - slope_b_seen).abs() > tol.angular()
                || !close(
                    a.reference_radius(),
                    rise.mul_add(slope_b, b.reference_radius()),
                )
            {
                return None;
            }
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, rise),
            }
        }
        (S::Sphere(a), S::Sphere(b)) => {
            let (a, b) = (a.sphere(), b.sphere());
            if !close(a.radius(), b.radius()) {
                return None;
            }
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            if fa.origin.distance(fb.origin) > tol.confusion() {
                return None;
            }
            // Latitude is the angle off the equator, so it runs the other
            // way where the axes do.
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, 0.0),
            }
        }
        (S::Torus(a), S::Torus(b)) => {
            let (a, b) = (a.torus(), b.torus());
            if !close(a.major_radius(), b.major_radius())
                || !close(a.minor_radius(), b.minor_radius())
            {
                return None;
            }
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            if fa.origin.distance(fb.origin) > tol.confusion() {
                return None;
            }
            // The tube angle rises along the axis, so it too runs the other
            // way where the axes do.
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, 0.0),
            }
        }
        (S::Extrusion(a), S::Extrusion(b)) => {
            use ogeom_geom::Transformable as _;
            let basis_a = a.curve().transformed(&from.placement, tol).ok()?;
            let basis_b = b.curve().transformed(&to.placement, tol).ok()?;
            let da = from.placement.apply_vector(a.direction().vector());
            let db = to.placement.apply_vector(b.direction().vector());
            let along = da.dot(db);
            if (along.abs() - 1.0).abs() > tol.angular() {
                return None;
            }
            // `b`'s basis at `t` is `a`'s at `s` moved `lift` along `da`,
            // `t` being `±s + shift`, so `a`'s `(s, v)` is `b`'s
            // `(±s + shift, ±(v - lift))`.
            let moved = match translate_of(&basis_a, &basis_b, da, tol) {
                Some(lift) => BasisMap {
                    flip: false,
                    shift: 0.0,
                    lift,
                },
                None => {
                    plain = false;
                    basis_map(&basis_a, &basis_b, da, tol)?
                }
            };
            ChartMap {
                flip: (moved.flip, along < 0.0),
                shift: Vector2::new(
                    moved.shift,
                    if along < 0.0 { moved.lift } else { -moved.lift },
                ),
            }
        }
        _ => return None,
    };
    // The map is read off the frames; it must also carry every point of
    // the face where the other chart has it.
    let (lo, hi) = from.window;
    const GRID: u32 = 4;
    for i in 0..=GRID {
        for j in 0..=GRID {
            let (fu, fv) = (
                f64::from(i) / f64::from(GRID),
                f64::from(j) / f64::from(GRID),
            );
            let at = Point2::new(
                (hi.x - lo.x).mul_add(fu, lo.x),
                (hi.y - lo.y).mul_add(fv, lo.y),
            );
            let (Some(p), Some(q)) = (
                world_point(from, at, tol),
                world_point(to, map.apply(at), tol),
            ) else {
                return None;
            };
            if p.distance(q) > tol.confusion() {
                return None;
            }
        }
    }
    Some((map, plain))
}

/// The point at `(u, v)` of a carrier's chart, placed in the world, read
/// off the unbounded surface so a chart point past the stored window still
/// names its point.
fn world_point(c: &Curved, at: Point2, tol: Tolerances) -> Option<Point> {
    use ogeom_math::elementary as e;
    let local = match &c.surface {
        SurfaceGeometry::Cylinder(s) => e::cylinder_at(&s.cylinder(), at.x, at.y).point,
        SurfaceGeometry::Cone(s) => e::cone_at(&s.cone(), at.x, at.y).point,
        SurfaceGeometry::Sphere(s) => e::sphere_at(&s.sphere(), at.x, at.y).point,
        SurfaceGeometry::Torus(s) => e::torus_at(&s.torus(), at.x, at.y).point,
        SurfaceGeometry::Extrusion(s) => {
            s.curve().point_at(at.x, tol).ok()? + s.direction().vector() * at.y
        }
        _ => return None,
    };
    Some(c.placement.apply(local))
}

/// How far along `d` the curve `b` is `a` moved, where it is `a` moved
/// along `d` with its parameter kept: a line, circle or B-spline whose
/// defining points all moved by one multiple of `d` and whose shape and
/// parameter range are the same. Moving every defining point by at most
/// the confusion distance moves every point of the curve by no more.
fn translate_of(a: &Curve, b: &Curve, d: Vector, tol: Tolerances) -> Option<f64> {
    let close = |x: f64, y: f64| (x - y).abs() <= tol.confusion();
    let same_range = {
        let ((a0, a1), (b0, b1)) = (a.domain(), b.domain());
        (a0 - b0).abs() <= tol.parametric() && (a1 - b1).abs() <= tol.parametric()
    };
    if !same_range {
        return None;
    }
    match (a, b) {
        (Curve::Line(la), Curve::Line(lb)) => {
            let (aa, ab) = (la.axis(), lb.axis());
            let (a0, a1) = a.domain();
            let reach = a0.abs().max(a1.abs());
            let turned = (aa.direction.vector() - ab.direction.vector()).magnitude() * reach;
            if turned > tol.confusion() {
                return None;
            }
            along_axis(ab.location - aa.location, d, tol)
        }
        (Curve::Circle(ca), Curve::Circle(cb)) => {
            let (ra, rb) = (ca.circle(), cb.circle());
            if ca.is_reversed() != cb.is_reversed() || !close(ra.radius(), rb.radius()) {
                return None;
            }
            let (fa, fb) = (ra.frame(), rb.frame());
            let r = ra.radius();
            if (fa.x().vector() - fb.x().vector()).magnitude() * r > tol.confusion()
                || (fa.y().vector() - fb.y().vector()).magnitude() * r > tol.confusion()
            {
                return None;
            }
            along_axis(fb.origin() - fa.origin(), d, tol)
        }
        (Curve::BSpline(sa), Curve::BSpline(sb)) => {
            let (ka, kb) = (sa.knots(), sb.knots());
            let (pa, pb) = (sa.control_points(), sb.control_points());
            if ka.degree() != kb.degree()
                || ka.knots().len() != kb.knots().len()
                || ka
                    .knots()
                    .iter()
                    .zip(kb.knots())
                    .any(|(x, y)| (x - y).abs() > tol.parametric())
                || a.is_periodic() != b.is_periodic()
                || pa.len() != pb.len()
                || pa.is_empty()
            {
                return None;
            }
            let mut lift: Option<f64> = None;
            for (p, q) in pa.iter().zip(pb) {
                if (p.weight - q.weight).abs() > 1e-12 * p.weight.abs() || p.weight <= 0.0 {
                    return None;
                }
                let moved = along_axis(q.point() - p.point(), d, tol)?;
                match lift {
                    None => lift = Some(moved),
                    Some(first) if close(first, moved) => {}
                    Some(_) => return None,
                }
            }
            lift
        }
        _ => None,
    }
}

/// How one extrusion's basis curve is another's moved along the sweep: the
/// point of `b` at `t` is the point of `a` at `s` moved `lift` along the
/// sweep, where `t` is `s` (or `-s` where `flip`) plus `shift`.
#[derive(Debug, Clone, Copy)]
struct BasisMap {
    flip: bool,
    shift: f64,
    lift: f64,
}

/// How `b` is `a` moved along `d`, its parameter run either way and
/// shifted: the general form of [`translate_of`].
///
/// Each kind is compared through the vectors that state it, so the map is
/// exact and not a fit:
/// - a line `L + r·D`, by its direction (one way or the other) and its
///   origin, moved along itself and along `d`;
/// - a circle or ellipse `O + Q·cos r + R·sin r`, by its axis vectors turned
///   by the angle `γ` and its sine run either way;
/// - a hyperbola `O + Q·cosh r + R·sinh r` or a parabola `O + Q·r² + R·r`,
///   by the even part `Q` alike and the odd part `R` run either way;
/// - a B-spline, by its knots shifted, or mirrored, and its weighted
///   control points moved alike, in the same order or the reverse.
///
/// Each vector's error is weighed by the largest value its function takes
/// over `b`'s parameter range, so no point of `b` strays further than the
/// confusion distance per term. `None` where `b` is not `a` so moved, or
/// for any other kind of curve.
fn basis_map(a: &Curve, b: &Curve, d: Vector, tol: Tolerances) -> Option<BasisMap> {
    let conf = tol.confusion();
    match (a, b) {
        (Curve::Line(la), Curve::Line(lb)) => {
            let (aa, ab) = (la.axis(), lb.axis());
            let (da, db) = (aa.direction.vector(), ab.direction.vector());
            // A line swept along itself sweeps no surface.
            let across = da.cross(d).magnitude();
            if across <= tol.angular() {
                return None;
            }
            let e = if da.dot(db) < 0.0 { -1.0 } else { 1.0 };
            let (b0, b1) = b.domain();
            if (db - da * e).magnitude() * b0.abs().max(b1.abs()) > conf {
                return None;
            }
            // `Lb - La = c·Da + lift·d`, so `b(r) = a(e·r + c) + lift·d`.
            let offset = ab.location - aa.location;
            let k = da.dot(d);
            let (p, q) = (offset.dot(da), offset.dot(d));
            let det = across * across;
            let c = k.mul_add(-q, p) / det;
            let lift = k.mul_add(-p, q) / det;
            if (offset - da * c - d * lift).magnitude() > conf {
                return None;
            }
            Some(BasisMap {
                flip: e < 0.0,
                shift: -e * c,
                lift,
            })
        }
        (Curve::Circle(_) | Curve::Ellipse(_), Curve::Circle(_) | Curve::Ellipse(_)) => {
            let (oa, qa, ra, sa) = trig_parts(a)?;
            let (ob, qb, rb, sb) = trig_parts(b)?;
            // `a`'s angle `σa·s` is `ε·σb·t + γ`: `Qb = Qa·cos γ + Ra·sin γ`
            // and `Rb = ε·(Ra·cos γ - Qa·sin γ)`.
            let gamma = (qb.dot(ra) / ra.dot(ra)).atan2(qb.dot(qa) / qa.dot(qa));
            let (sin, cos) = gamma.sin_cos();
            let turned = ra * cos - qa * sin;
            let e = if rb.dot(turned) < 0.0 { -1.0 } else { 1.0 };
            if (qb - (qa * cos + ra * sin)).magnitude() > conf
                || (rb - turned * e).magnitude() > conf
            {
                return None;
            }
            let lift = along_axis(ob - oa, d, tol)?;
            // `s = E·t + σa·γ` with `E = σa·ε·σb`, so `t = E·(s - σa·γ)`.
            let big_e = sa * e * sb;
            Some(BasisMap {
                flip: big_e < 0.0,
                shift: -big_e * sa * gamma,
                lift,
            })
        }
        (Curve::Hyperbola(_), Curve::Hyperbola(_)) | (Curve::Parabola(_), Curve::Parabola(_)) => {
            let (oa, qa, ra, sa, ka) = even_odd_parts(a)?;
            let (ob, qb, rb, sb, kb) = even_odd_parts(b)?;
            // `b`'s raw parameter at its ends, and the most its even and
            // odd functions reach there.
            let (b0, b1) = b.domain();
            let reach = sb.mul_add(b0, kb).abs().max(sb.mul_add(b1, kb).abs());
            let (even, odd) = match b {
                Curve::Hyperbola(_) => (reach.cosh(), reach.sinh()),
                _ => (reach * reach, reach),
            };
            // `a`'s raw parameter is `ε` times `b`'s.
            let e = if rb.dot(ra) < 0.0 { -1.0 } else { 1.0 };
            if (qb - qa).magnitude() * even > conf || (rb - ra * e).magnitude() * odd > conf {
                return None;
            }
            let lift = along_axis(ob - oa, d, tol)?;
            // `σa·s + κa = ε·(σb·t + κb)`, so `s = E·t + σa·(ε·κb - κa)`
            // with `E = σa·ε·σb`.
            let big_e = sa * e * sb;
            let offset = sa * e.mul_add(kb, -ka);
            Some(BasisMap {
                flip: big_e < 0.0,
                shift: -big_e * offset,
                lift,
            })
        }
        (Curve::BSpline(sa), Curve::BSpline(sb)) => {
            let (ka, kb) = (sa.knots(), sb.knots());
            let (pa, pb) = (sa.control_points(), sb.control_points());
            if ka.degree() != kb.degree()
                || ka.knots().len() != kb.knots().len()
                || a.is_periodic() != b.is_periodic()
                || pa.len() != pb.len()
                || pa.is_empty()
            {
                return None;
            }
            let (ka, kb) = (ka.knots(), kb.knots());
            // `b`'s knots are `a`'s shifted by `δ`, or mirrored in `δ`; the
            // control points follow in the same order or the reverse.
            for flip in [false, true] {
                let at = |i: usize| if flip { ka[ka.len() - 1 - i] } else { ka[i] };
                let delta = if flip { kb[0] + at(0) } else { kb[0] - at(0) };
                let knots_agree = (0..kb.len()).all(|i| {
                    let mapped = if flip { delta - at(i) } else { at(i) + delta };
                    (kb[i] - mapped).abs() <= tol.parametric()
                });
                if !knots_agree {
                    continue;
                }
                let mut lift: Option<f64> = None;
                let mut agree = true;
                for (i, q) in pb.iter().enumerate() {
                    let p = &pa[if flip { pa.len() - 1 - i } else { i }];
                    if (p.weight - q.weight).abs() > 1e-12 * p.weight.abs() || p.weight <= 0.0 {
                        agree = false;
                        break;
                    }
                    let Some(moved) = along_axis(q.point() - p.point(), d, tol) else {
                        agree = false;
                        break;
                    };
                    match lift {
                        None => lift = Some(moved),
                        Some(first) if (first - moved).abs() <= conf => {}
                        Some(_) => {
                            agree = false;
                            break;
                        }
                    }
                }
                if let (true, Some(lift)) = (agree, lift) {
                    return Some(BasisMap {
                        flip,
                        shift: delta,
                        lift,
                    });
                }
            }
            None
        }
        _ => None,
    }
}

/// A circle or ellipse as `O + Q·cos(σ·u) + R·sin(σ·u)`: its centre, its
/// axis vectors scaled by its radii, and `σ`, which is `-1` where it runs
/// backwards.
fn trig_parts(c: &Curve) -> Option<(Point, Vector, Vector, f64)> {
    let sign = |reversed: bool| if reversed { -1.0 } else { 1.0 };
    match c {
        Curve::Circle(c) => {
            let (circle, f) = (c.circle(), c.circle().frame());
            let r = circle.radius();
            Some((
                circle.centre(),
                f.x().vector() * r,
                f.y().vector() * r,
                sign(c.is_reversed()),
            ))
        }
        Curve::Ellipse(c) => {
            let (ellipse, f) = (c.ellipse(), c.ellipse().frame());
            Some((
                ellipse.centre(),
                f.x().vector() * ellipse.major_radius(),
                f.y().vector() * ellipse.minor_radius(),
                sign(c.is_reversed()),
            ))
        }
        _ => None,
    }
}

/// A hyperbola as `O + Q·cosh r + R·sinh r`, or a parabola as
/// `O + Q·r² + R·r`, with `r = σ·u + κ`: its centre or apex, its even and
/// odd vectors, and `σ` and `κ`, which run a reversed curve back over its
/// range.
fn even_odd_parts(c: &Curve) -> Option<(Point, Vector, Vector, f64, f64)> {
    let (lo, hi) = c.domain();
    let run = |reversed: bool| {
        if reversed {
            (-1.0, lo + hi)
        } else {
            (1.0, 0.0)
        }
    };
    match c {
        Curve::Hyperbola(h) => {
            let (curve, f) = (h.hyperbola(), h.hyperbola().frame());
            let (s, k) = run(h.is_reversed());
            Some((
                curve.centre(),
                f.x().vector() * curve.major_radius(),
                f.y().vector() * curve.minor_radius(),
                s,
                k,
            ))
        }
        Curve::Parabola(p) => {
            let (curve, f) = (p.parabola(), p.parabola().frame());
            let (s, k) = run(p.is_reversed());
            Some((
                curve.apex(),
                f.x().vector() / (4.0 * curve.focal()),
                f.y().vector(),
                s,
                k,
            ))
        }
        _ => None,
    }
}

/// The plain pcurve `edge` carries on `surface`, where it carries one.
fn pcurve_on(model: &Model, edge: &Shape, surface: SurfaceId) -> Option<(PlanarCurve, (f64, f64))> {
    let data = model.node(edge)?.data().as_edge()?;
    let EdgeRepr::PCurve { curve, range, .. } = data.pcurve_for(surface, edge.location())? else {
        return None;
    };
    Some((model.geometry().pcurve(*curve)?.clone(), *range))
}

/// The pcurve a face's use of `edge` reads on `surface`: the edge's plain
/// pcurve, or the side of its seam the use walks, which is the forward side
/// where the face's own node walks the edge forward. `face_reversed` is
/// whether the face is used reversed, which `edge`'s orientation carries.
/// The flag says whether the edge is a seam there.
fn use_pcurve(
    model: &Model,
    edge: &Shape,
    surface: SurfaceId,
    face_reversed: bool,
) -> Option<(PlanarCurve, (f64, f64), bool)> {
    let data = model.node(edge)?.data().as_edge()?;
    let (id, range, seam) = match data.pcurve_for(surface, edge.location())? {
        EdgeRepr::PCurve { curve, range, .. } => (*curve, *range, false),
        EdgeRepr::Seam {
            forward,
            reversed,
            range,
            ..
        } => {
            let node_forward = (edge.orientation() == Orientation::Forward) != face_reversed;
            (
                if node_forward { *forward } else { *reversed },
                *range,
                true,
            )
        }
        _ => return None,
    };
    Some((model.geometry().pcurve(id)?.clone(), range, seam))
}

/// Whether `edge` carries any parameterization on `surface` already.
fn has_slot(model: &Model, edge: &Shape, surface: SurfaceId) -> bool {
    model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .is_some_and(|d| d.pcurve_for(surface, edge.location()).is_some())
}

/// `pcurve` read along `range` at evenly spaced parameters, backwards where
/// `backwards`.
fn samples(
    pcurve: &PlanarCurve,
    range: (f64, f64),
    backwards: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
    let mut out = Vec::with_capacity(SAMPLES as usize);
    for k in 0..SAMPLES {
        let f = f64::from(k) / f64::from(SAMPLES - 1);
        let f = if backwards { 1.0 - f } else { f };
        out.push(pcurve.point_at((range.1 - range.0).mul_add(f, range.0), tol)?);
    }
    Ok(out)
}

/// The chart's period in each direction, where it repeats or closes.
fn periods(surface: &SurfaceGeometry, tol: Tolerances) -> (Option<f64>, Option<f64>) {
    let ((ua, ub), (va, vb)) = surface.domain();
    (
        (surface.is_periodic_u() || surface.is_closed_u(tol)).then_some(ub - ua),
        (surface.is_periodic_v() || surface.is_closed_v(tol)).then_some(vb - va),
    )
}

/// Whether two chart points that name one point of the carrier are one
/// chart point too, and not a period apart.
fn same_sheet(a: Point2, b: Point2, periods: (Option<f64>, Option<f64>)) -> bool {
    periods.0.is_none_or(|p| (a.x - b.x).abs() < p / 4.0)
        && periods.1.is_none_or(|p| (a.y - b.y).abs() < p / 4.0)
}

/// Twice the area a chain of chart points sweeps about `anchor`: summed
/// over closed rings, what they enclose, positive counter-clockwise.
fn swept(points: &[Point2], anchor: Point2) -> f64 {
    points
        .windows(2)
        .map(|w| (w[0] - anchor).cross(w[1] - anchor))
        .sum()
}

/// The cluster's members, as indices into the faces, each with its
/// carrier; the first is the keeper, whose chart the merged face takes.
/// `wires` holds each wire the members walk, as a run of the walks.
pub(super) struct Members<'a> {
    pub(super) members: &'a [usize],
    pub(super) carriers: &'a [&'a Curved],
    pub(super) wires: &'a [Range<usize>],
}

/// Carry a cluster of faces on one curved carrier into its first face's
/// chart, and chain the merged boundary into rings, the outer one first,
/// each a list of indices into `walks`. `walks` holds each member's edges
/// as it walks them, seen from the side they all face, tagged with the
/// member's index into the faces; `kept` marks those no two members share.
///
/// Each use of an edge reads its own pcurve, a seam's use the side it
/// walks. A member's wire is read as one line through its chart, each
/// pcurve moved by whole periods to start where the one before it ended:
/// an edge two faces on one surface share carries one pcurve, which only
/// one of them reads on its own sheet. Each member's chart is mapped onto
/// the keeper's, the period it lands on settled by an edge it shares with
/// a member already placed; the members' pcurves are carried across by
/// that map, which is exact, and measured against the members' own lifts.
/// An edge two members share dissolves where its two sides land on one
/// sheet of the chart, and stays as a seam where they land a period apart
/// (two halves of a bore make a whole one). A member's own seam stays a
/// seam, its two sides carried across. The rings are chained through the
/// chart, so a seam's two sides are told apart, and seams of two members
/// that meet end to end run on as one seam line.
///
/// A member on the keeper's own surface, placed where the keeper is, keeps
/// the sheet its stored pcurves give it.
///
/// The cluster is left as it is (`None`) where a member walks an edge twice
/// other than down a seam's two sides, a member's seam stores its sides
/// other than a period apart or its wire reads them elsewhere, a shared
/// edge's sides land neither together nor a period apart, a ring does not
/// close in the chart (a band round the axis with no seam), a seam would
/// run along its period or change direction at a vertex, the rings do not
/// wind as one outer ring round its holes, or the keeper's window
/// cannot be grown to hold the others'. Nothing is changed in the model:
/// the edits come back with the rings, for the caller to make once the
/// merged face is built.
pub(super) fn merge(
    model: &mut Model,
    cluster: &Members<'_>,
    walks: &[(usize, Shape)],
    kept: &[bool],
    reversed: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Merged>> {
    let Members {
        members,
        carriers,
        wires,
    } = *cluster;
    let keeper = carriers[0];
    let periods = periods(&keeper.surface, tol);
    let slot_of: FastMap<usize, usize> = members.iter().enumerate().map(|(k, &i)| (i, k)).collect();

    // Each walk's member, as a slot in `members`, and the pcurve its use
    // reads there.
    let mut by_node: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    let mut own: Vec<(usize, PlanarCurve, (f64, f64))> = Vec::with_capacity(walks.len());
    for (w, (i, edge)) in walks.iter().enumerate() {
        let slot = slot_of[i];
        let Some((pcurve, range, seam)) =
            use_pcurve(model, edge, carriers[slot].id, carriers[slot].reversed)
        else {
            return Ok(None);
        };
        let uses = by_node.entry(edge.node()).or_default();
        if uses
            .iter()
            .any(|&o| walks[o].0 == *i && !(seam && walks[o].1.orientation() != edge.orientation()))
        {
            return Ok(None);
        }
        uses.push(w);
        own.push((slot, pcurve, range));
    }
    let (shared, own_seams) = {
        let mut shared: Vec<(TShapeId, usize, usize)> = Vec::new();
        let mut own_seams: Vec<(usize, usize)> = Vec::new();
        for (node, uses) in &by_node {
            if let &[a, b] = uses.as_slice() {
                if walks[a].0 == walks[b].0 {
                    own_seams.push((a.min(b), a.max(b)));
                } else {
                    shared.push((*node, a, b));
                }
            }
        }
        shared.sort_by_key(|&(_, a, _)| a);
        own_seams.sort_unstable();
        (shared, own_seams)
    };

    // Each walk read in its member's chart, in the direction it is walked,
    // and the whole periods that put it where the walk before it ended.
    let mut raw: Vec<Vec<Point2>> = Vec::with_capacity(walks.len());
    for (w, (_, edge)) in walks.iter().enumerate() {
        let (_, pcurve, range) = &own[w];
        let backwards = edge.orientation() == Orientation::Reversed;
        raw.push(samples(pcurve, *range, backwards, tol)?);
    }
    let mut wrap = vec![Vector2::ZERO; walks.len()];
    for run in wires {
        let Some(first) = run.clone().next() else {
            continue;
        };
        let member_periods = self::periods(&carriers[own[first].0].surface, tol);
        let mut end: Option<Point2> = None;
        for w in run.clone() {
            let (Some(&start), Some(&last)) = (raw[w].first(), raw[w].last()) else {
                return Ok(None);
            };
            if let Some(end) = end {
                wrap[w] = whole_periods(end - start, member_periods);
            }
            end = Some(last + wrap[w]);
        }
        // The wire stays on the sheet most of its walks were stored on.
        let mut counts: Vec<(Vector2, usize)> = Vec::new();
        for w in run.clone() {
            match counts.iter_mut().find(|(v, _)| *v == wrap[w]) {
                Some((_, n)) => *n += 1,
                None => counts.push((wrap[w], 1)),
            }
        }
        let mut most = counts[0];
        for &(v, n) in &counts[1..] {
            if n > most.1 {
                most = (v, n);
            }
        }
        for w in run.clone() {
            wrap[w] -= most.0;
        }
    }

    // Each member's chart placed on the keeper's, reached across the edges
    // the members share.
    let mut maps: Vec<Option<ChartMap>> = vec![None; members.len()];
    maps[0] = Some(ChartMap::IDENTITY);
    let mut placed = 1;
    while placed < members.len() {
        let before = placed;
        for &(_, a, b) in &shared {
            let (sa, sb) = (own[a].0, own[b].0);
            let (from, onto, unplaced) = match (maps[sa], maps[sb]) {
                (Some(m), None) => (b, (a, m), sb),
                (None, Some(m)) => (a, (b, m), sa),
                _ => continue,
            };
            let Some(base) = chart_map(carriers[unplaced], keeper, tol) else {
                return Ok(None);
            };
            // A member stored on the keeper's own surface, placed where the
            // keeper is, is in the keeper's chart already: the pcurves its
            // edges keep there say which sheet it lies on.
            if carriers[unplaced].id == keeper.id
                && base.flip == (false, false)
                && base.shift == Vector2::ZERO
            {
                maps[unplaced] = Some(ChartMap::IDENTITY);
                placed += 1;
                continue;
            }
            // Both sides read at the middle of the edge.
            let mid = |w: usize| -> OgeomResult<Point2> {
                let (_, pcurve, range) = &own[w];
                Ok(pcurve.point_at(0.5 * (range.0 + range.1), tol)? + wrap[w])
            };
            let there = onto.1.apply(mid(onto.0)?);
            let here = base.apply(mid(from)?);
            let shift = whole_periods(there - here, periods);
            maps[unplaced] = Some(ChartMap {
                flip: base.flip,
                shift: base.shift + shift,
            });
            placed += 1;
        }
        if placed == before {
            return Ok(None);
        }
    }
    let maps: Vec<ChartMap> = maps.into_iter().flatten().collect();

    // Every walk read in the keeper's chart, each point measured against
    // the member's own lift.
    let mut lines: Vec<Vec<Point2>> = Vec::with_capacity(walks.len());
    for (w, points) in raw.iter().enumerate() {
        let slot = own[w].0;
        let mut line = Vec::with_capacity(points.len());
        for &p in points {
            let q = maps[slot].apply(p + wrap[w]);
            if slot != 0 {
                let (Some(theirs), Some(ours)) = (
                    world_point(carriers[slot], p, tol),
                    world_point(keeper, q, tol),
                ) else {
                    return Ok(None);
                };
                if theirs.distance(ours) > tol.confusion() {
                    return Ok(None);
                }
            }
            line.push(q);
        }
        lines.push(line);
    }

    // Each shared edge dissolves, or is a seam a whole period across.
    let mut kept = kept.to_vec();
    let mut seams: Vec<(TShapeId, usize, usize, Vector2)> = Vec::new();
    for &(node, a, b) in &shared {
        let pairs = || lines[a].iter().zip(lines[b].iter().rev());
        if pairs().all(|(p, q)| same_sheet(*p, *q, periods)) {
            continue;
        }
        let (Some(&p), Some(&q)) = (lines[a].first(), lines[b].last()) else {
            return Ok(None);
        };
        let across = whole_periods(p - q, periods);
        if across == Vector2::ZERO || !pairs().all(|(p, q)| same_sheet(*p, *q + across, periods)) {
            return Ok(None);
        }
        // A seam runs straight across the period it closes. An edge that
        // runs along it, such as an arc between two bands' seams that start
        // at different angles, would make the seam a staircase that meshes,
        // booleans and offsets do not read as a seam.
        let spread = |at: fn(&Point2) -> f64| {
            let (lo, hi) = lines[a]
                .iter()
                .map(at)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                    (lo.min(x), hi.max(x))
                });
            hi - lo
        };
        if (across.x != 0.0 && spread(|p| p.x) > tol.parametric())
            || (across.y != 0.0 && spread(|p| p.y) > tol.parametric())
        {
            return Ok(None);
        }
        kept[a] = true;
        kept[b] = true;
        seams.push((node, a, b, across));
    }
    // A member's own seam stores its two sides a whole period apart, and
    // its wire reads both where they are stored: two sides on one sheet
    // close no band, and a wire that reads them otherwise is not the face
    // the seam bounds.
    let mut seam_line: FastMap<TShapeId, (bool, bool)> = seams
        .iter()
        .map(|&(node, _, _, across)| (node, (across.x != 0.0, across.y != 0.0)))
        .collect();
    for &(a, b) in &own_seams {
        let (Some(&p), Some(&q)) = (raw[a].first(), raw[b].last()) else {
            return Ok(None);
        };
        let across = whole_periods(p - q, periods);
        if across == Vector2::ZERO
            || wrap[a] != wrap[b]
            || !raw[a]
                .iter()
                .zip(raw[b].iter().rev())
                .all(|(p, q)| same_sheet(*p, *q + across, periods))
        {
            return Ok(None);
        }
        seam_line.insert(walks[a].1.node(), (across.x != 0.0, across.y != 0.0));
    }

    // Chain the kept walks into rings, head to tail through the chart: of
    // the walks leaving a vertex, the one starting where the last ended.
    let boundary: Vec<usize> = (0..walks.len()).filter(|&w| kept[w]).collect();
    if boundary.is_empty() {
        return Ok(None);
    }
    let mut ends: Vec<(TShapeId, TShapeId)> = Vec::with_capacity(boundary.len());
    let mut leaving: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    for (k, &w) in boundary.iter().enumerate() {
        let Some((start, end)) = ogeom_algo::edge_vertices(model, &walks[w].1)? else {
            return Ok(None);
        };
        leaving.entry(start.node()).or_default().push(k);
        ends.push((start.node(), end.node()));
    }
    let first_of = |k: usize| lines[boundary[k]].first().copied();
    let last_of = |k: usize| lines[boundary[k]].last().copied();
    let mut chained = vec![false; boundary.len()];
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for first in 0..boundary.len() {
        if chained[first] {
            continue;
        }
        let mut ring = Vec::new();
        let mut k = first;
        loop {
            chained[k] = true;
            ring.push(boundary[k]);
            let Some(end) = last_of(k) else {
                return Ok(None);
            };
            let mut next = leaving
                .get(&ends[k].1)
                .map_or(&[][..], Vec::as_slice)
                .iter()
                .copied()
                .filter(|&n| first_of(n).is_some_and(|start| same_sheet(start, end, periods)));
            let (Some(next), None) = (next.next(), next.next()) else {
                return Ok(None);
            };
            if next == first {
                break;
            }
            if chained[next] {
                return Ok(None);
            }
            k = next;
        }
        rings.push(ring);
    }

    // A seam line runs one way: where a ring walks from one seam edge
    // straight on into another closing the same period, both are walked
    // alike, head to tail. Two seam edges running opposite ways from the
    // vertex between them make a seam line its readers do not follow.
    for ring in &rings {
        for (k, &w) in ring.iter().enumerate() {
            let next = ring[(k + 1) % ring.len()];
            let (here, there) = (walks[w].1.node(), walks[next].1.node());
            if here == there {
                continue;
            }
            if let (Some(a), Some(b)) = (seam_line.get(&here), seam_line.get(&there))
                && a == b
                && walks[w].1.orientation() != walks[next].1.orientation()
            {
                return Ok(None);
            }
        }
    }

    // The side the members' boundaries wind about, all alike.
    let anchor = lines[0][0];
    let mut winding = 0.0f64;
    for slot in 0..members.len() {
        let area: f64 = lines
            .iter()
            .zip(&own)
            .filter(|(_, o)| o.0 == slot)
            .map(|(line, _)| swept(line, anchor))
            .sum();
        if area == 0.0 || (winding != 0.0 && area.signum() != winding) {
            return Ok(None);
        }
        winding = area.signum();
    }

    // One merged ring winds as the members do, the rest the other way.
    let mut outer = None;
    for (r, ring) in rings.iter().enumerate() {
        let mut points: Vec<Point2> = Vec::new();
        for &w in ring {
            points.extend_from_slice(&lines[w]);
        }
        points.push(points[0]);
        let area = swept(&points, anchor);
        if area == 0.0 {
            return Ok(None);
        }
        if area.signum() == winding && outer.replace(r).is_some() {
            return Ok(None);
        }
    }
    let Some(outer) = outer else {
        return Ok(None);
    };

    // The keeper's window grown to hold every member's, carried across.
    let Some(grown) = grown_window(keeper, carriers, &maps, tol) else {
        return Ok(None);
    };

    // The boundary's pcurves on the keeper's surface: its own where it has
    // them, the members' carried across where the edge has none there yet,
    // and where it has one already, that one only if it lies where the
    // carried one does.
    let carried = |w: usize| -> OgeomResult<(PlanarCurve, (f64, f64))> {
        let (slot, pcurve, range) = &own[w];
        let map = maps[*slot].transform() * Transform2::translation(wrap[w]);
        Ok((pcurve.transformed(&map, tol)?, *range))
    };
    let lies_along = |there: &(PlanarCurve, (f64, f64)),
                      ours: &(PlanarCurve, (f64, f64))|
     -> OgeomResult<bool> {
        let theirs = samples(&there.0, there.1, false, tol)?;
        let ours = samples(&ours.0, ours.1, false, tol)?;
        Ok(theirs
            .iter()
            .zip(&ours)
            .all(|(p, q)| same_sheet(*p, *q, periods)))
    };
    // The walk of a seam's two that reads its forward side, as the merged
    // face's own node walks it.
    let forward_walk = |a: usize, b: usize| {
        if (walks[a].1.orientation() == Orientation::Forward) != reversed {
            (a, b)
        } else {
            (b, a)
        }
    };
    let seam_nodes: Vec<TShapeId> = seams.iter().map(|s| s.0).collect();
    let mut attach: Vec<(Shape, PlanarCurve, (f64, f64))> = Vec::new();
    for (w, (_, edge)) in walks.iter().enumerate() {
        if !kept[w] || own[w].0 == 0 || seam_nodes.contains(&edge.node()) {
            continue;
        }
        if own_seams.iter().any(|&(a, b)| a == w || b == w) {
            continue;
        }
        let ours = carried(w)?;
        if has_slot(model, edge, keeper.id) {
            let Some(there) = pcurve_on(model, edge, keeper.id) else {
                return Ok(None);
            };
            if !lies_along(&there, &ours)? {
                return Ok(None);
            }
            continue;
        }
        attach.push((edge.clone(), ours.0, ours.1));
    }
    // A seam's Forward use, as the merged face's own node walks it, reads
    // the side that use walks. A seam two members make reads one side off
    // the other a whole period across, so both sides are one curve; a
    // member's own seam carries both of its sides across.
    let mut seamed: Vec<(Shape, PlanarCurve, PlanarCurve, (f64, f64))> = Vec::new();
    for &(_, a, b, across) in &seams {
        let (pa, range) = carried(a)?;
        let pb = pa.transformed(&Transform2::translation(-across), tol)?;
        let (forward, backward) = if forward_walk(a, b).0 == a {
            (pa, pb)
        } else {
            (pb, pa)
        };
        seamed.push((walks[a].1.clone(), forward, backward, range));
    }
    for &(a, b) in &own_seams {
        if own[a].0 == 0 {
            continue;
        }
        let (f, r) = forward_walk(a, b);
        let ((forward, range), (backward, _)) = (carried(f)?, carried(r)?);
        let edge = &walks[a].1;
        if has_slot(model, edge, keeper.id) {
            let (Some((there_f, ..)), Some((there_r, ..))) = (
                use_pcurve(model, &walks[f].1, keeper.id, reversed),
                use_pcurve(model, &walks[r].1, keeper.id, reversed),
            ) else {
                return Ok(None);
            };
            if !lies_along(&(there_f, range), &(forward, range))?
                || !lies_along(&(there_r, range), &(backward, range))?
            {
                return Ok(None);
            }
            continue;
        }
        seamed.push((edge.clone(), forward, backward, range));
    }

    rings.swap(0, outer);
    Ok(Some(Merged {
        rings,
        edits: Edits {
            surface: keeper.id,
            grown,
            attach,
            seamed,
        },
    }))
}

/// What [`merge`] found: the merged boundary's rings, the outer one first,
/// as indices into the walks, and the edits that put the keeper's surface
/// and the boundary's pcurves in its chart.
pub(super) struct Merged {
    pub(super) rings: Vec<Vec<usize>>,
    pub(super) edits: Edits,
}

/// The model edits a merge makes, held back until the merged face is built.
pub(super) struct Edits {
    surface: SurfaceId,
    grown: Option<SurfaceGeometry>,
    attach: Vec<(Shape, PlanarCurve, (f64, f64))>,
    seamed: Vec<(Shape, PlanarCurve, PlanarCurve, (f64, f64))>,
}

impl Edits {
    /// Grow the keeper's window, give the boundary its pcurves on the
    /// keeper's surface, and make the seams seams there.
    pub(super) fn apply(self, model: &mut Model) -> OgeomResult<()> {
        let Self {
            surface,
            grown,
            attach,
            seamed,
        } = self;
        if let Some(grown) = grown
            && let Some(held) = model.geometry_mut().surface_mut(surface)
        {
            *held = grown;
        }
        for (edge, pcurve, range) in attach {
            if !has_slot(model, &edge, surface) {
                ogeom_algo::attach_pcurve(
                    model,
                    &edge,
                    pcurve,
                    surface,
                    edge.location().clone(),
                    range,
                )?;
            }
        }
        for (edge, forward, backward, range) in seamed {
            if let Some(node) = model.node_mut(&edge)
                && let NodeData::Edge(data) = node.data_mut()
            {
                data.representations.retain(|r| {
                    r.surface() != Some(surface) || r.location() != Some(edge.location())
                });
            }
            ogeom_algo::attach_seam(
                model,
                &edge,
                forward,
                backward,
                surface,
                edge.location().clone(),
                range,
            )?;
        }
        Ok(())
    }
}

/// The whole periods nearest `gap` in each direction the chart repeats.
fn whole_periods(gap: Vector2, periods: (Option<f64>, Option<f64>)) -> Vector2 {
    let wrap = |g: f64, period: Option<f64>| period.map_or(0.0, |p| (g / p).round() * p);
    Vector2::new(wrap(gap.x, periods.0), wrap(gap.y, periods.1))
}

/// The keeper's surface with its bounded window grown to hold every
/// member's carried into its chart: `Some(None)` where it holds them
/// already, `None` where a member reaches past what the keeper's chart
/// has (an extrusion's basis curve ends there).
fn grown_window(
    keeper: &Curved,
    carriers: &[&Curved],
    maps: &[ChartMap],
    tol: Tolerances,
) -> Option<Option<SurfaceGeometry>> {
    let ((ua, ub), (va, vb)) = keeper.surface.domain();
    let (pu, _) = periods(&keeper.surface, tol);
    let (mut lo, mut hi) = (va, vb);
    for (c, map) in carriers.iter().zip(maps).skip(1) {
        let ((a, b), (c0, c1)) = c.surface.domain();
        let corners = [map.apply(Point2::new(a, c0)), map.apply(Point2::new(b, c1))];
        if pu.is_none()
            && corners
                .iter()
                .any(|p| p.x < ua - tol.parametric() || p.x > ub + tol.parametric())
        {
            return None;
        }
        for p in corners {
            (lo, hi) = (lo.min(p.y), hi.max(p.y));
        }
    }
    if (lo, hi) == (va, vb) {
        return Some(None);
    }
    match &keeper.surface {
        SurfaceGeometry::Cylinder(s) => ogeom_geom::CylinderSurface::new(s.cylinder(), (lo, hi))
            .ok()
            .map(|s| Some(s.into())),
        SurfaceGeometry::Cone(s) => ogeom_geom::ConeSurface::new(s.cone(), (lo, hi))
            .ok()
            .map(|s| Some(s.into())),
        SurfaceGeometry::Extrusion(s) => {
            ogeom_geom::ExtrusionSurface::over(s.curve().clone(), s.direction(), (lo, hi))
                .ok()
                .map(|s| Some(s.into()))
        }
        // A sphere's latitude and a torus's tube angle span their whole
        // range already.
        _ => Some(None),
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_geom::{Reversible as _, Transformable as _};
    use ogeom_math::{Axis, Direction, KnotVector};

    const T: Tolerances = Tolerances::millimetres();

    /// Whether `map` takes `a` onto `b` moved along `d`: at every sampled
    /// parameter `t` of `b` whose `s` (the parameter the map names) `a`
    /// holds too, `b(t)` is `a(s) + lift·d`, read off the curves' own
    /// evaluation. At least half the samples must be compared; the copies
    /// overlap that much.
    fn carries(a: &Curve, b: &Curve, d: Vector, map: BasisMap) -> bool {
        let (lo, hi) = b.domain();
        let (a0, a1) = a.domain();
        let mut compared = 0;
        for k in 0..=16 {
            let t = (hi - lo).mul_add(f64::from(k) / 16.0, lo);
            let s = if map.flip {
                map.shift - t
            } else {
                t - map.shift
            };
            if !a.is_periodic() && (s < a0 || s > a1) {
                continue;
            }
            compared += 1;
            let there = a.point_at(s, T).unwrap() + d * map.lift;
            if b.point_at(t, T).unwrap().distance(there) > 1e-9 {
                return false;
            }
        }
        compared >= 8
    }

    /// Each kind of basis curve against copies of itself run backwards,
    /// mirrored in a plane along the sweep, turned about the sweep or
    /// shifted along itself, each lifted along the sweep: the map found
    /// carries one onto the other point for point. A copy moved across
    /// the sweep, or an ellipse turned by other than a half turn, is no
    /// such copy.
    #[test]
    fn a_basis_moved_along_the_sweep_is_found_whichever_way_it_runs() {
        let d = Vector::new(0.0, 0.0, 1.0);
        let lift = Transform::translation(d * 2.5);
        let mirror = Transform::plane_mirror(Point::ORIGIN, Direction::X);
        // The conics' frame sits on the mirror plane, and the turns are
        // about its own axis.
        let centre = Point::new(0.0, -0.25, 0.0);
        let turn = |angle: f64| Transform::rotation(Axis::new(centre, Direction::Z), angle);
        let frame = Frame::new(centre, Direction::Z, Direction::Y, T).unwrap();
        let line = Curve::Line(
            ogeom_geom::LineCurve::over(
                Axis::new(Point::new(0.5, 2.0, 0.0), Direction::X),
                -2.0,
                3.0,
            )
            .unwrap(),
        );
        let circle = Curve::Circle(ogeom_geom::CircleCurve::new(
            ogeom_math::Circle::new(frame, 2.0, T).unwrap(),
        ));
        let ellipse = Curve::Ellipse(ogeom_geom::EllipseCurve::new(
            ogeom_math::Ellipse::new(frame, 4.0, 2.0, T).unwrap(),
        ));
        let hyperbola = Curve::Hyperbola(
            ogeom_geom::HyperbolaCurve::over(
                ogeom_math::Hyperbola::new(frame, 2.0, 1.5, T).unwrap(),
                -0.5,
                1.2,
            )
            .unwrap(),
        );
        let parabola = Curve::Parabola(
            ogeom_geom::ParabolaCurve::over(
                ogeom_math::Parabola::new(frame, 0.75, T).unwrap(),
                -1.0,
                3.0,
            )
            .unwrap(),
        );
        let spline = Curve::BSpline(
            ogeom_geom::BSplineCurve::new(
                KnotVector::new(vec![0.0, 0.0, 0.0, 0.3, 1.0, 1.0, 1.0], 2).unwrap(),
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(1.0, 2.0, 0.0),
                    Point::new(3.0, 1.0, 0.0),
                    Point::new(4.0, 3.0, 0.0),
                ],
                T,
            )
            .unwrap(),
        );
        let moved = |c: &Curve, t: &Transform| c.transformed(t, T).unwrap();
        let pairs = [
            (line.clone(), moved(&line.reversed(), &lift)),
            (line.clone(), moved(&line, &(lift * mirror))),
            (circle.clone(), moved(&circle, &(lift * turn(0.7)))),
            (circle.clone(), moved(&circle.reversed(), &(lift * mirror))),
            (ellipse.clone(), moved(&ellipse, &(lift * mirror))),
            (ellipse.clone(), moved(&ellipse.reversed(), &lift)),
            (
                ellipse.clone(),
                moved(&ellipse, &(lift * turn(core::f64::consts::PI))),
            ),
            (hyperbola.clone(), moved(&hyperbola.reversed(), &lift)),
            (hyperbola.clone(), moved(&hyperbola, &(lift * mirror))),
            (parabola.clone(), moved(&parabola.reversed(), &lift)),
            (parabola.clone(), moved(&parabola, &(lift * mirror))),
            (spline.clone(), moved(&spline.reversed(), &lift)),
        ];
        for (a, b) in &pairs {
            let map = basis_map(a, b, d, T);
            assert!(
                map.is_some_and(|map| carries(a, b, d, map)),
                "{a:?} against {b:?}: {map:?}"
            );
        }
        let across = Transform::translation(Vector::new(0.0, 0.5, 2.5));
        for (a, b) in [
            (ellipse.clone(), moved(&ellipse, &(lift * turn(0.7)))),
            (hyperbola.clone(), moved(&hyperbola, &across)),
            (spline.clone(), moved(&spline.reversed(), &across)),
        ] {
            assert!(basis_map(&a, &b, d, T).is_none(), "{a:?} against {b:?}");
        }
    }
}
