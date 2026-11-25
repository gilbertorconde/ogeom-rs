//! Faces from closed planar wires: which wire bounds which face, and which
//! way each one runs.
//!
//! Every edge is read as its exact image in the plane (a line, a conic, a
//! B-spline's projected control net), so nesting, crossing and orientation
//! are decided on the curves themselves and never on a sampled polygon:
//!
//! - Two wires meet where their curves' closed-form or Newton-polished
//!   intersections say they do. Between those points each stretch of the
//!   smaller wire lies wholly on one side of the larger, and one point of it
//!   decides which, by the signed crossings of a ray with the larger
//!   wire's curves. A ray through a vertex of that wire, or grazing one of
//!   its curves, is replaced by another.
//! - Orientation is the sign of the enclosed area, Green's theorem over each
//!   edge: in closed form for lines, circles and ellipses, and by a
//!   Gauss rule per knot span that is exact for polynomial B-splines up to
//!   degree ten. Other curves integrate adaptively to rounding.
//! - Only wires whose bounding boxes overlap are compared, found by one
//!   sweep along the plane's `x`: disjoint wires cost nothing.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, Line2d, PlanarCurve, PlaneSurface, SurfaceGeometry,
    Trimmed2d,
};
use ogeom_intersect::{CurveCurveOptions, exact_pcurve_over, intersect_curves_2d};
use ogeom_math::{Axis2, Direction2, Frame, Handedness, Plane, Point, Point2, Vector2};
use ogeom_topo::{Model, Orientation, Shape, ShapeType};

use crate::build::{
    is_wire_closed, make_face_on, placed_edge_curve, trimmed_where_bare, walked_back,
};
use crate::measure::curve_bounds_over;

/// The faces [`faces_from_planar_wires`] built, and the wires bounding each.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlanarFaces {
    /// One face per region, in the order of their outer wires' indices.
    pub faces: Vec<Shape>,
    /// Per face, its outer wire's index into the given wires, then its
    /// holes' in increasing order.
    pub wire_of: Vec<Vec<usize>>,
}

/// The faces closed planar `wires` bound under the even-odd rule.
///
/// A wire inside an odd number of the others is a hole in the innermost
/// one around it; a wire inside an even number, none included, is the outer
/// boundary of a new face. So an island drawn inside a hole is a face of
/// its own, and so is anything inside that island's holes. Each face lies
/// on `plane`, faces its normal, and walks its outer wire counterclockwise
/// about that normal and its holes clockwise, whichever way the wires were
/// given and in whatever order. A wire walked the wrong way is replaced in
/// the face by a fresh wire walking the same edges back; `wire_of` still
/// names the given one.
///
/// Wires of two different faces may touch at isolated points, a wire
/// tangent to another included: two discs meeting at a point, or an island
/// tangent to the hole it sits in. Each wire is still wholly inside or
/// outside the other, and nesting follows from that; the faces are built as
/// they stand, without merging or splitting at the touching point. Two
/// wires of one face may not touch (a hole tangent to its outer wire, two
/// holes meeting): that face's boundary would meet itself.
///
/// Containment, crossing and orientation are decided on the exact curves,
/// never on samples of them: see the module documentation.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction),
/// naming the wire (or the two wires) at fault, if a given shape is not a
/// wire, a wire is open, encloses no area, has an edge off `plane` or with
/// no exact image in it, crosses itself, or if two wires cross, run along
/// each other for a stretch, or touch while bounding one face.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if no ray from
/// a point of one wire crossed another cleanly.
pub fn faces_from_planar_wires(
    model: &mut Model,
    plane: Plane,
    wires: &[Shape],
    tol: Tolerances,
) -> OgeomResult<PlanarFaces> {
    Ok(arranged(model, plane, wires, tol)?.0)
}

/// [`faces_from_planar_wires`], with how many points it classified.
fn arranged(
    model: &mut Model,
    plane: Plane,
    wires: &[Shape],
    tol: Tolerances,
) -> OgeomResult<(PlanarFaces, usize)> {
    let plane = right_handed(plane, tol)?;
    let surface = SurfaceGeometry::Plane(PlaneSurface::new(plane));
    let rings = wires
        .iter()
        .enumerate()
        .map(|(i, wire)| ring_of(model, i, wire, plane, &surface, tol))
        .collect::<OgeomResult<Vec<_>>>()?;
    let mut probe = Probe {
        margin: tol.confusion(),
        tests: 0,
        tol,
    };
    for (i, ring) in rings.iter().enumerate() {
        probe.simple(i, ring)?;
        // A wire closing on itself within the confusion distance encloses
        // no more than this.
        if ring.area.abs() <= tol.confusion() * tol.confusion().max(ring.bound.diagonal()) {
            ogeom_bail!(Construction, "wire {i} encloses no area");
        }
    }

    // Which wires each lies inside.
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); rings.len()];
    let mut touching = Vec::new();
    for (a, b) in overlapping_pairs(&rings, probe.margin) {
        let (inner, outer) = if rings[a].area.abs() <= rings[b].area.abs() {
            (a, b)
        } else {
            (b, a)
        };
        let met = probe.meet((inner, &rings[inner]), (outer, &rings[outer]))?;
        if met.inside {
            around[inner].push(outer);
        }
        if met.touching {
            touching.push((a, b));
        }
    }

    // Non-crossing wires nest as a forest: the containers of a wire are a
    // chain, and the innermost of them is the one inside all the others.
    let depth: Vec<usize> = around.iter().map(Vec::len).collect();
    // Each wire's face: its own where it is an outer wire, its parent's
    // where it is a hole.
    let mut face_of: Vec<usize> = (0..rings.len()).collect();
    let mut holes: Vec<Vec<usize>> = vec![Vec::new(); rings.len()];
    for (i, containers) in around.iter().enumerate() {
        if depth[i].is_multiple_of(2) {
            continue;
        }
        let Some(&parent) = containers.iter().max_by_key(|&&c| depth[c]) else {
            continue;
        };
        if depth[parent] + 1 != depth[i] {
            ogeom_bail!(
                Construction,
                "wire {i} sits inside wires {containers:?}, which do not nest; \
                 two of them cross"
            );
        }
        holes[parent].push(i);
        face_of[i] = parent;
    }

    // A face whose own wires touch has a boundary meeting itself, which no
    // valid face has; wires of two faces may touch.
    for (a, b) in touching {
        if face_of[a] == face_of[b] {
            ogeom_bail!(
                Construction,
                "wires {a} and {b} touch, and both would bound one face; a face's \
                 boundary may not meet itself"
            );
        }
    }

    let mut faces = Vec::new();
    let mut wire_of = Vec::new();
    for outer in (0..rings.len()).filter(|&i| depth[i].is_multiple_of(2)) {
        let mut bounding = vec![outer];
        bounding.extend(holes[outer].iter().copied());
        let mut wound = Vec::with_capacity(bounding.len());
        for (k, &i) in bounding.iter().enumerate() {
            let counterclockwise = rings[i].area > 0.0;
            wound.push(if counterclockwise == (k == 0) {
                wires[i].clone()
            } else {
                walked_back(model, &wires[i], tol)?
            });
        }
        let id = model.geometry_mut().add_surface(surface.clone());
        let face = make_face_on(model, id, &wound, tol)?.shape;
        // Every edge has an exact image in the plane, read above.
        trimmed_where_bare(model, &face, tol)?;
        faces.push(face);
        wire_of.push(bounding);
    }
    Ok((PlanarFaces { faces, wire_of }, probe.tests))
}

/// `plane` with a right-handed frame and the same normal, so that the
/// plane's chart runs counterclockwise about the normal and a face on it
/// faces that normal.
fn right_handed(plane: Plane, tol: Tolerances) -> OgeomResult<Plane> {
    let frame = plane.frame();
    if frame.handedness() == Handedness::Right {
        return Ok(plane);
    }
    Ok(Plane::new(Frame::from_axes(
        frame.origin(),
        frame.x(),
        -frame.y(),
        frame.z(),
        tol,
    )?))
}

/// An axis-aligned box in the plane's chart.
#[derive(Debug, Clone, Copy)]
struct Bound {
    lo: Point2,
    hi: Point2,
}

impl Bound {
    fn of(points: impl IntoIterator<Item = Point2>) -> Self {
        let mut bound = Self {
            lo: Point2::new(f64::INFINITY, f64::INFINITY),
            hi: Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        };
        for p in points {
            bound.lo = Point2::new(bound.lo.x.min(p.x), bound.lo.y.min(p.y));
            bound.hi = Point2::new(bound.hi.x.max(p.x), bound.hi.y.max(p.y));
        }
        bound
    }

    fn joined(self, other: Self) -> Self {
        Self::of([self.lo, self.hi, other.lo, other.hi])
    }

    fn overlaps(&self, other: &Self, margin: f64) -> bool {
        self.lo.x <= other.hi.x + margin
            && other.lo.x <= self.hi.x + margin
            && self.lo.y <= other.hi.y + margin
            && other.lo.y <= self.hi.y + margin
    }

    fn diagonal(&self) -> f64 {
        (self.hi - self.lo).magnitude()
    }
}

/// One edge of a wire, in the plane's chart.
#[derive(Debug, Clone)]
struct Piece {
    /// The edge's exact image, restricted to the stretch the edge uses.
    curve: PlanarCurve,
    /// That stretch, in the image's parameter (the edge's own).
    range: (f64, f64),
    /// Whether the wire walks the edge from `range.1` to `range.0`.
    reversed: bool,
    /// Where the stretch starts and ends, in the image's parameter order.
    ends: [Point2; 2],
    bound: Bound,
}

/// A closed wire in the plane's chart.
#[derive(Debug, Clone)]
struct Ring {
    pieces: Vec<Piece>,
    /// The signed area it encloses: positive walked counterclockwise.
    area: f64,
    bound: Bound,
}

/// Wire `i` read into `plane`'s chart.
fn ring_of(
    model: &Model,
    i: usize,
    wire: &Shape,
    plane: Plane,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Ring> {
    if model.kind_of(wire)? != ShapeType::Wire {
        ogeom_bail!(Construction, "shape {i} is not a wire");
    }
    if !is_wire_closed(model, wire, tol)? {
        ogeom_bail!(
            Construction,
            "wire {i} is open; an open wire bounds no region"
        );
    }
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };
    let mut pieces = Vec::new();
    for (k, edge) in model.ordered_children_of(wire)?.iter().enumerate() {
        let Some((curve, range)) = placed_edge_curve(model, edge, tol)? else {
            ogeom_bail!(Construction, "wire {i}: edge {k} has no curve in space");
        };
        if !lies_in(&curve, range, plane, tol)? {
            ogeom_bail!(
                Construction,
                "wire {i} does not lie in the plane: its edge {k} stands off it"
            );
        }
        let Some(image) = exact_pcurve_over(&curve, range, surface, tol) else {
            ogeom_bail!(
                Construction,
                "wire {i}: edge {k}'s curve has no exact image in the plane"
            );
        };
        let curve2 = windowed(image, range, tol)?;
        let ends = [
            curve2.point_at(range.0, tol)?,
            curve2.point_at(range.1, tol)?,
        ];
        let bound = {
            let space = curve_bounds_over(&curve, range, tol)?;
            let (Some(lo), Some(hi)) = (space.low(), space.high()) else {
                ogeom_bail!(Construction, "wire {i}: edge {k} has no extent");
            };
            Bound::of(
                [
                    Point::new(lo.x, lo.y, lo.z),
                    Point::new(hi.x, lo.y, lo.z),
                    Point::new(lo.x, hi.y, lo.z),
                    Point::new(hi.x, hi.y, lo.z),
                    Point::new(lo.x, lo.y, hi.z),
                    Point::new(hi.x, lo.y, hi.z),
                    Point::new(lo.x, hi.y, hi.z),
                    Point::new(hi.x, hi.y, hi.z),
                ]
                .map(flat),
            )
            .joined(Bound::of(ends))
        };
        pieces.push(Piece {
            curve: curve2,
            range,
            reversed: edge.orientation() == Orientation::Reversed,
            ends,
            bound,
        });
    }
    let Some(first) = pieces.first() else {
        ogeom_bail!(Construction, "wire {i} has no edges");
    };
    // About a point of the wire itself: the sum is the same about any point,
    // and about one far off it is a difference of large products.
    let anchor = first.ends[0];
    let bound = pieces
        .iter()
        .map(|p| p.bound)
        .reduce(Bound::joined)
        .unwrap_or(first.bound);
    let mut area = 0.0;
    for piece in &pieces {
        let swept = swept_area(piece, anchor, bound.diagonal(), tol)?;
        area += if piece.reversed { -swept } else { swept };
    }
    Ok(Ring {
        pieces,
        area,
        bound,
    })
}

/// Whether `curve` over `range` lies in `plane`: decided from what fixes
/// the curve (a segment's ends, a conic's frame, a B-spline's control net,
/// which holds the curve in its hull). A curve drawn on a surface is
/// left to its exact image, which exists only where that surface is this
/// plane.
fn lies_in(curve: &Curve, range: (f64, f64), plane: Plane, tol: Tolerances) -> OgeomResult<bool> {
    let normal = plane.frame().z();
    let near = |p: Point| plane.distance_to(p) <= tol.confusion();
    let conic = |frame: Frame| frame.z().is_parallel(normal, tol) && near(frame.origin());
    Ok(match curve {
        Curve::Line(_) => {
            near(curve.point_at(range.0, tol)?) && near(curve.point_at(range.1, tol)?)
        }
        Curve::Circle(c) => conic(c.circle().frame()),
        Curve::Ellipse(e) => conic(e.ellipse().frame()),
        Curve::Hyperbola(h) => conic(h.hyperbola().frame()),
        Curve::Parabola(p) => conic(p.parabola().frame()),
        Curve::BSpline(b) => b.control_points().iter().all(|w| near(w.point())),
        Curve::Trimmed(t) => {
            let window = t.domain();
            let range = if t.is_reversed() {
                let sum = window.0 + window.1;
                (sum - range.1, sum - range.0)
            } else {
                range
            };
            return lies_in(t.basis(), range, plane, tol);
        }
        Curve::Offset(o) => {
            o.reference().is_parallel(normal, tol) && lies_in(o.basis(), range, plane, tol)?
        }
        Curve::OnSurface(_) => true,
        Curve::Helix(_) => false,
    })
}

/// `image` restricted to `range`, through any trim it came with: a trim
/// shares its basis's parameter, so one window serves both.
fn windowed(image: PlanarCurve, range: (f64, f64), tol: Tolerances) -> OgeomResult<PlanarCurve> {
    let basis = match image {
        PlanarCurve::Trimmed(t) if !t.is_reversed() => t.basis().clone(),
        other => other,
    };
    let (lo, hi) = basis.domain();
    if range.0 <= lo + tol.parametric() && range.1 >= hi - tol.parametric() {
        return Ok(basis);
    }
    Ok(Trimmed2d::new(basis, range.0, range.1, tol)?.into())
}

/// Half the integral of `(P - anchor) x P'` along `piece` from `range.0` to
/// `range.1`: the area the stretch sweeps about `anchor`, by Green's
/// theorem.
fn swept_area(piece: &Piece, anchor: Point2, scale: f64, tol: Tolerances) -> OgeomResult<f64> {
    let curve = &piece.curve;
    let (t0, t1) = piece.range;
    let basis = match curve {
        PlanarCurve::Trimmed(t) if !t.is_reversed() => t.basis(),
        other => other,
    };
    // A circle or an ellipse about its centre sweeps at a constant rate:
    // `(P - c) x P'` is the product of its radii, signed by the way it runs.
    let conic = |centre: Point2| -> OgeomResult<f64> {
        let [p0, p1] = piece.ends;
        let rate = (p0 - centre).cross(curve.d1_at(t0, tol)?);
        Ok(0.5 * ((centre - anchor).cross(p1 - p0) + rate * (t1 - t0)))
    };
    match basis {
        PlanarCurve::Line(_) => {
            let [p0, p1] = piece.ends;
            return Ok(0.5 * (p0 - anchor).cross(p1 - anchor));
        }
        PlanarCurve::Circle(c) => return conic(c.circle().centre()),
        PlanarCurve::Ellipse(e) => return conic(e.ellipse().centre()),
        _ => {}
    }
    let integrand = |t: f64| -> OgeomResult<f64> {
        Ok((curve.point_at(t, tol)? - anchor).cross(curve.d1_at(t, tol)?))
    };
    // Split at the knots, where a B-spline's pieces join.
    let mut breaks = vec![t0];
    let mut polynomial = false;
    if let PlanarCurve::BSpline(b) = basis {
        breaks.extend(
            b.knots()
                .distinct()
                .into_iter()
                .map(|(k, _)| k)
                .filter(|&k| k > t0 && k < t1),
        );
        // `(P - a) x P'` is a polynomial of degree `2p - 1` on a span,
        // within the ten-point rule's nineteen.
        polynomial = !b.is_rational() && b.knots().degree() <= 10;
    }
    breaks.push(t1);
    let mut total = 0.0;
    for span in breaks.windows(2) {
        let (a, b) = (span[0], span[1]);
        if polynomial {
            for (t, w) in ogeom_math::integrate::gauss_legendre_rule(a, b) {
                total += w * integrand(t)?;
            }
        } else {
            let mut failed = None;
            let value = ogeom_math::integrate::integrate(
                |t| {
                    integrand(t).unwrap_or_else(|e| {
                        failed.get_or_insert(e);
                        0.0
                    })
                },
                a,
                b,
                // Rounding on an area of the wire's own size.
                (1e-14 * scale * scale).max(f64::MIN_POSITIVE),
            )?;
            if let Some(e) = failed {
                return Err(e);
            }
            total += value;
        }
    }
    Ok(0.5 * total)
}

/// Every pair of rings whose bounds overlap, by one sweep along `x`.
fn overlapping_pairs(rings: &[Ring], margin: f64) -> Vec<(usize, usize)> {
    let mut order: Vec<usize> = (0..rings.len()).collect();
    order.sort_by(|&a, &b| rings[a].bound.lo.x.total_cmp(&rings[b].bound.lo.x));
    let mut active: Vec<usize> = Vec::new();
    let mut pairs = Vec::new();
    for &i in &order {
        let here = rings[i].bound;
        active.retain(|&j| rings[j].bound.hi.x + margin >= here.lo.x);
        pairs.extend(
            active
                .iter()
                .filter(|&&j| rings[j].bound.overlaps(&here, margin))
                .map(|&j| (j.min(i), j.max(i))),
        );
        active.push(i);
    }
    pairs
}

/// Where a point stands against a ring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    In,
    Out,
    On,
}

/// How two rings that do not cross stand.
#[derive(Debug, Clone, Copy)]
struct Met {
    /// The smaller lies inside the larger.
    inside: bool,
    /// They meet at isolated points.
    touching: bool,
}

/// How far round each ray turns from the last: the golden angle, so no
/// two of the few tried line up with each other or with a drawing's axes.
const RAY_TURN: f64 = 2.399_963_229_728_653;

/// How many rays to try before giving up.
const RAYS: u32 = 12;

/// The sine of the angle under which a ray must cross a curve to be
/// counted, rather than replaced by another ray.
const GRAZING: f64 = 1e-6;

/// The predicates, and a count of the points classified.
struct Probe {
    margin: f64,
    tests: usize,
    tol: Tolerances,
}

impl Probe {
    fn options(&self) -> CurveCurveOptions {
        CurveCurveOptions::for_tolerances(self.tol)
    }

    /// Refuses a ring whose edges cross or run over each other anywhere
    /// but where two of them share an end.
    fn simple(&self, i: usize, ring: &Ring) -> OgeomResult<()> {
        let tol = self.tol;
        for (k, a) in ring.pieces.iter().enumerate() {
            for b in &ring.pieces[k + 1..] {
                if !a.bound.overlaps(&b.bound, self.margin) {
                    continue;
                }
                let at_shared_end = |p: Point2| {
                    a.ends.iter().any(|e| (p - *e).magnitude() <= self.margin)
                        && b.ends.iter().any(|e| (p - *e).magnitude() <= self.margin)
                };
                let met = intersect_curves_2d(&a.curve, &b.curve, self.options(), tol)?;
                for overlap in &met.overlaps {
                    if !self.is_point(&a.curve, overlap.on_a)? {
                        ogeom_bail!(Construction, "wire {i} runs back over itself");
                    }
                }
                if met.crossings.iter().any(|c| !at_shared_end(c.point)) {
                    ogeom_bail!(Construction, "wire {i} crosses itself");
                }
            }
        }
        Ok(())
    }

    /// Whether the stretch `on` of `curve` is no longer than the margin.
    fn is_point(&self, curve: &PlanarCurve, on: (f64, f64)) -> OgeomResult<bool> {
        let tol = self.tol;
        let start = curve.point_at(on.0, tol)?;
        let middle = curve.point_at(f64::midpoint(on.0, on.1), tol)?;
        let end = curve.point_at(on.1, tol)?;
        Ok((end - start).magnitude() <= self.margin && (middle - start).magnitude() <= self.margin)
    }

    /// Whether ring `inner` lies inside ring `outer`, and whether the two
    /// touch, refusing the pair where they cross or run along each other.
    ///
    /// The inner ring is cut where the two meet; between the cuts each
    /// stretch lies on one side of the outer ring, and the middle of each
    /// says which.
    fn meet(&mut self, inner: (usize, &Ring), outer: (usize, &Ring)) -> OgeomResult<Met> {
        let tol = self.tol;
        let ((i, small), (j, large)) = (inner, outer);
        let mut stretches: Vec<(&Piece, (f64, f64))> = Vec::new();
        let mut met_any = false;
        for piece in &small.pieces {
            let mut cuts = vec![piece.range.0, piece.range.1];
            if piece.bound.overlaps(&large.bound, self.margin) {
                for other in &large.pieces {
                    if !piece.bound.overlaps(&other.bound, self.margin) {
                        continue;
                    }
                    let met = intersect_curves_2d(&piece.curve, &other.curve, self.options(), tol)?;
                    for overlap in &met.overlaps {
                        if !self.is_point(&piece.curve, overlap.on_a)? {
                            ogeom_bail!(
                                Construction,
                                "wires {} and {} run along each other",
                                i.min(j),
                                i.max(j)
                            );
                        }
                        cuts.extend([overlap.on_a.0, overlap.on_a.1]);
                    }
                    cuts.extend(met.crossings.iter().map(|c| c.on_a));
                    met_any |= !met.is_empty();
                }
            }
            cuts.sort_by(f64::total_cmp);
            stretches.extend(
                cuts.windows(2)
                    .filter(|w| w[1] - w[0] > tol.parametric())
                    .map(|w| (piece, (w[0], w[1]))),
            );
        }
        if !met_any {
            // Clear of each other: any one point of the inner ring decides.
            stretches.truncate(1);
        }
        let mut found: Option<Side> = None;
        for (piece, (a, b)) in stretches {
            let p = piece.curve.point_at(f64::midpoint(a, b), tol)?;
            match self.side_of(p, large)? {
                Side::On => {}
                side if found.is_some_and(|f| f != side) => {
                    ogeom_bail!(Construction, "wires {} and {} cross", i.min(j), i.max(j))
                }
                side => found = Some(side),
            }
        }
        match found {
            Some(side) => Ok(Met {
                inside: side == Side::In,
                touching: met_any,
            }),
            None => ogeom_bail!(
                Construction,
                "wires {} and {} run along each other",
                i.min(j),
                i.max(j)
            ),
        }
    }

    /// Where `p` stands against `ring`, by the signed crossings of a ray
    /// from it with the ring's curves.
    fn side_of(&mut self, p: Point2, ring: &Ring) -> OgeomResult<Side> {
        let tol = self.tol;
        self.tests += 1;
        if !ring.bound.overlaps(&Bound::of([p]), self.margin) {
            return Ok(Side::Out);
        }
        let reach = (p - ring.bound.lo)
            .magnitude()
            .max((p - ring.bound.hi).magnitude())
            + ring.bound.diagonal()
            + 1.0;
        // A ray this near a vertex may meet both edges there or neither.
        let clearance = (100.0 * self.margin).max(1e-9 * ring.bound.diagonal());
        'ray: for k in 0..RAYS {
            let angle = 0.3 + f64::from(k) * RAY_TURN;
            let d = Vector2::new(angle.cos(), angle.sin());
            for piece in &ring.pieces {
                for v in piece.ends {
                    let along = (v - p).dot(d).clamp(0.0, reach);
                    let off = (v - (p + d * along)).magnitude();
                    if along == 0.0 && off <= self.margin {
                        return Ok(Side::On);
                    }
                    if off <= clearance {
                        continue 'ray;
                    }
                }
            }
            let ray: PlanarCurve =
                Line2d::over(Axis2::new(p, Direction2::new(d, tol)?), 0.0, reach)?.into();
            let ray_bound = Bound::of([p, p + d * reach]);
            let mut winding = 0_i64;
            for piece in &ring.pieces {
                if !piece.bound.overlaps(&ray_bound, self.margin) {
                    continue;
                }
                let met = intersect_curves_2d(&ray, &piece.curve, self.options(), tol)?;
                if let Some(overlap) = met.overlaps.first() {
                    if overlap.on_a.0.min(overlap.on_a.1) <= self.margin {
                        return Ok(Side::On);
                    }
                    continue 'ray;
                }
                for crossing in &met.crossings {
                    if crossing.on_a <= self.margin {
                        return Ok(Side::On);
                    }
                    let mut tangent = piece.curve.d1_at(crossing.on_b, tol)?;
                    if piece.reversed {
                        tangent = -tangent;
                    }
                    let speed = tangent.magnitude();
                    let sine = if speed > 0.0 {
                        d.cross(tangent) / speed
                    } else {
                        0.0
                    };
                    if sine.abs() < GRAZING {
                        continue 'ray;
                    }
                    winding += if sine > 0.0 { 1 } else { -1 };
                }
            }
            return Ok(if winding == 0 { Side::Out } else { Side::In });
        }
        ogeom_bail!(
            NotDone,
            "every ray from a point of a wire met another wire at a vertex or a tangent"
        )
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests;
