//! Minimum distance between shapes.
//!
//! *Elsewhere* this is `BRepExtrema_DistShapeShape`. The geometry-level
//! extrema in `ogeom-intersect` answer where two curves or surfaces come
//! nearest; this module assembles those answers for topology, where a shape
//! is vertices, edges and faces and the nearest approach may land on any of
//! them.
//!
//! # The assembly argument
//!
//! The nearest distance between two shapes is attained either at an interior
//! stationary approach of a pair of elements, or on some element's boundary,
//! and an element's boundary is itself an element: a face's boundary is its
//! edges, an edge's boundary is its vertices. So walking every pair of
//! elements (vertex against vertex, edge and face; edge against edge and
//! face; face against face) with stationary approaches for the interiors and
//! projections for the points covers every candidate, and the geometry level
//! is allowed to answer "no interior approach" honestly because the pair that
//! owns the boundary case is in the same sweep.
//!
//! A face's interior approach is accepted only where its foot lands inside
//! the face's trimming; a foot outside or too near the boundary is dropped,
//! because the true nearest point of that configuration is on an edge and the
//! edge pairs find it exactly.
//!
//! # What the distance is between
//!
//! Boundaries. A shape strictly inside another reports the gap between their
//! boundaries, not zero: whether a point is *inside* a solid is
//! [`classify_in_solid_exact`](crate::classify_in_solid_exact)'s question,
//! and conflating the two would make this answer wrong for the shells it is
//! right for.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, SurfaceGeometry, Transformable, TrimmedCurve};
use ogeom_math::{Point, Point2, Transform, Vector};
use ogeom_mesh::{Deflection, face_boundary, inside_boundary};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

use crate::classify::parametric_band;
use crate::measure::{project_on_curve, project_on_surface};
use ogeom_intersect::ExtremaOptions;

/// One pair of nearest points, with the elements they lie on.
#[derive(Debug, Clone)]
pub struct ClosestPair {
    /// The nearest point on the first shape.
    pub point_a: Point,
    /// The nearest point on the second.
    pub point_b: Point,
    /// The vertex, edge or face of the first shape the point lies on.
    pub support_a: Shape,
    /// The same for the second shape.
    pub support_b: Shape,
}

/// The minimum distance between two shapes, with everywhere it is attained.
#[derive(Debug, Clone)]
pub struct ShapeDistance {
    /// The distance.
    pub distance: f64,
    /// Every pair of nearest points found within tolerance of the minimum.
    /// Parallel walls meet at a representative pair, not at every point of
    /// the overlap.
    pub pairs: Vec<ClosestPair>,
}

/// One element of a shape, with its geometry carried into world space.
enum Element {
    Vertex(Shape, Point),
    Edge(Shape, Box<Curve>),
    Face(Shape, Box<Prepared>),
}

/// A face's surface twice over: in world space for the extrema, and local
/// with its placement and rings for the trim test. Parameters on a
/// transformed surface need not match the rings, which live in the stored
/// surface's parameter space, so the trim question is always asked locally.
struct Prepared {
    world: SurfaceGeometry,
    local: SurfaceGeometry,
    to_local: Transform,
    rings: Vec<Vec<Point2>>,
    /// The face's placement, carrying `local` into world space.
    placement: Transform,
    /// The parameter rectangle `world` is restricted to, when it is.
    span: Option<((f64, f64), (f64, f64))>,
}

/// The minimum distance between two shapes' boundaries.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either shape
/// has no vertices, edges or faces to measure to;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve.
pub fn distance_between_shapes(
    model: &Model,
    a: &Shape,
    b: &Shape,
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<ShapeDistance> {
    let ea = elements(model, a, tol)?;
    let eb = elements(model, b, tol)?;
    distance_between_elements(&ea, &eb, options, tol)
}

/// The elements of a shape (its vertices, edges and faces, placed) for
/// measuring it against others many times: see [`distance_between_elements`].
pub(crate) struct Elements(Vec<Element>);

impl Elements {
    pub(crate) fn of(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self(elements(model, shape, tol)?))
    }

    /// Whether every point of these elements is proven farther than `gap`
    /// from every point of `other`.
    ///
    /// Each element is held by the convex hull of a few points: a vertex by
    /// itself, a B-spline curve or patch by its control net, a line by its
    /// ends, a plane by its corners, any other curve by the corners of its
    /// bound. Two hulls are apart where their boxes are, or where their
    /// shadows on one axis are (a patch's normal, or the line between the
    /// pieces with a curve's chord taken out). Pairs not yet apart are
    /// refined by halving the larger piece, until every pair is apart or a
    /// pair cannot be refined within a fixed budget. `false` proves
    /// nothing: the pair is then measured.
    ///
    /// A face counts as its whole surface over the rectangle its rings
    /// enclose, the same restricted surface the distance is measured on,
    /// so the hulls hold every candidate the measurement could find.
    pub(crate) fn apart_by_more_than(&self, other: &Self, gap: f64, tol: Tolerances) -> bool {
        // Without a vertex on each side a measurement may find no candidate
        // at all and refuse; that refusal is the measurement's to make.
        let has_vertex = |e: &Self| e.0.iter().any(|e| matches!(e, Element::Vertex(..)));
        if !has_vertex(self) || !has_vertex(other) {
            return false;
        }
        let mut arena = Arena::default();
        let (Some(mine), Some(theirs)) = (arena.roots(&self.0, tol), arena.roots(&other.0, tol))
        else {
            return false;
        };
        let mut stack: Vec<(usize, usize)> = mine
            .iter()
            .flat_map(|&i| theirs.iter().map(move |&j| (i, j)))
            .collect();
        let mut visits = 0_usize;
        while let Some((i, j)) = stack.pop() {
            visits += 1;
            if visits > COVER_VISITS {
                return false;
            }
            if arena.nodes[i].apart_from(&arena.nodes[j], gap) {
                continue;
            }
            let (big, small) = if arena.nodes[i].size >= arena.nodes[j].size {
                (i, j)
            } else {
                (j, i)
            };
            let halves = match arena.halves(big, gap, tol) {
                Some(h) => Some((h, small)),
                None => arena.halves(small, gap, tol).map(|h| (h, big)),
            };
            let Some(([a, b], with)) = halves else {
                return false;
            };
            stack.push((a, with));
            stack.push((b, with));
        }
        true
    }
}

/// How many pairs of pieces [`Elements::apart_by_more_than`] looks at before
/// it leaves the pair to be measured.
const COVER_VISITS: usize = 20_000;

/// A piece of one element's geometry that can be held and halved.
enum Piece<'a> {
    /// Points whose hull holds the piece, with nothing finer to offer.
    Fixed(Vec<Point>),
    /// A non-periodic B-spline curve, held by its control polygon.
    Spline(Box<ogeom_geom::BSplineCurve>),
    /// A line or conic over a range.
    Arc(&'a Curve, (f64, f64)),
    /// A B-spline patch in its face's frame, held by its control net.
    Patch(Box<ogeom_geom::BSplineSurface>, &'a Transform),
    /// A plane over a parameter rectangle, held by its corners.
    Flat(
        &'a ogeom_geom::PlaneSurface,
        (f64, f64),
        (f64, f64),
        &'a Transform,
    ),
}

struct Node<'a> {
    piece: Piece<'a>,
    /// World points whose convex hull holds the piece.
    hull: Vec<Point>,
    bound: ogeom_math::Aabb,
    /// The bound's diagonal.
    size: f64,
    /// Across a patch or a plane: the normal of its corners.
    normal: Option<Vector>,
    /// Along a curve: from its start to its end.
    chord: Option<Vector>,
    /// The two halves once asked for; `Some(None)` where there are none.
    halves: Option<Option<[usize; 2]>>,
}

impl Node<'_> {
    /// Whether the two hulls are proven farther than `gap` apart.
    fn apart_from(&self, other: &Self, gap: f64) -> bool {
        if self.bound.distance_to_box(&other.bound) > gap {
            return true;
        }
        let (Some(here), Some(there)) = (self.bound.centre(), other.bound.centre()) else {
            return false;
        };
        let between = there - here;
        let square = |v: Vector| {
            let length = v.magnitude();
            (length > 0.0 && length.is_finite()).then(|| v / length)
        };
        // The line between the pieces with a curve's own direction taken
        // out, so two runs side by side are told apart across their gap.
        let across = |chord: Option<Vector>| {
            let chord = square(chord?)?;
            square(between - chord * between.dot(chord))
        };
        [
            self.normal.and_then(square),
            other.normal.and_then(square),
            square(between),
            across(self.chord),
            across(other.chord),
        ]
        .into_iter()
        .flatten()
        .any(|axis| {
            let shadow = |hull: &[Point]| {
                hull.iter()
                    .map(|p| (*p - Point::ORIGIN).dot(axis))
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
                        (lo.min(d), hi.max(d))
                    })
            };
            let (a, b) = (shadow(&self.hull), shadow(&other.hull));
            b.0 - a.1 > gap || a.0 - b.1 > gap
        })
    }
}

#[derive(Default)]
struct Arena<'a> {
    nodes: Vec<Node<'a>>,
}

impl<'a> Arena<'a> {
    /// One node per element, or `None` if an element has no hull to start
    /// from.
    fn roots(&mut self, elements: &'a [Element], tol: Tolerances) -> Option<Vec<usize>> {
        let mut out = Vec::with_capacity(elements.len());
        for element in elements {
            let piece = match element {
                Element::Vertex(_, p) => Piece::Fixed(vec![*p]),
                Element::Edge(_, curve) => curve_piece(curve, tol)?,
                Element::Face(_, face) => face_piece(face, tol)?,
            };
            out.push(self.add(piece, tol)?);
        }
        Some(out)
    }

    fn add(&mut self, piece: Piece<'a>, tol: Tolerances) -> Option<usize> {
        let (hull, normal, chord) = hull_of(&piece, tol)?;
        // An empty or unmeasurable hull would prove a distance it does not
        // hold.
        if hull.is_empty()
            || !hull
                .iter()
                .all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite())
        {
            return None;
        }
        let bound = ogeom_math::Aabb::of_points(&hull);
        self.nodes.push(Node {
            piece,
            size: bound.diagonal(),
            hull,
            bound,
            normal,
            chord,
            halves: None,
        });
        Some(self.nodes.len() - 1)
    }

    /// The node's two halves, made once; `None` where it cannot be halved
    /// or is no larger than `gap` already.
    fn halves(&mut self, i: usize, gap: f64, tol: Tolerances) -> Option<[usize; 2]> {
        if let Some(known) = self.nodes[i].halves {
            return known;
        }
        let found = if self.nodes[i].size <= gap {
            None
        } else {
            split(&self.nodes[i].piece, tol).and_then(|(a, b)| {
                let a = self.add(a, tol)?;
                let b = self.add(b, tol)?;
                Some([a, b])
            })
        };
        self.nodes[i].halves = Some(found);
        found
    }
}

/// An edge's curve as a piece; `None` where it has no bound.
fn curve_piece(curve: &Curve, tol: Tolerances) -> Option<Piece<'_>> {
    use ogeom_geom::Curve3d as _;
    let range = curve.domain();
    let whole = || {
        crate::measure::curve_bounds(curve, tol)
            .ok()
            .map(|bound| Piece::Fixed(bound.corners()))
    };
    match curve {
        Curve::BSpline(spline) if !spline.is_periodic() => {
            Some(Piece::Spline(Box::new(spline.clone())))
        }
        Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) => Some(Piece::Arc(curve, range)),
        Curve::Trimmed(trimmed) if !trimmed.is_reversed() => match trimmed.basis() {
            Curve::BSpline(spline) => spline
                .segment(range, tol)
                .map_or_else(|_| whole(), |piece| Some(Piece::Spline(Box::new(piece)))),
            basis @ (Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_)) => {
                Some(Piece::Arc(basis, range))
            }
            _ => whole(),
        },
        _ => whole(),
    }
}

/// A face's restricted surface as a piece; `None` for a surface this does
/// not hold.
fn face_piece(face: &Prepared, tol: Tolerances) -> Option<Piece<'_>> {
    match &face.local {
        SurfaceGeometry::BSpline(spline) => {
            let patch = match face.span {
                Some((u, v)) => spline.segment(u, v, tol).ok()?,
                None => spline.clone(),
            };
            Some(Piece::Patch(Box::new(patch), &face.placement))
        }
        SurfaceGeometry::Plane(plane) => {
            let (u, v) = face.span?;
            Some(Piece::Flat(plane, u, v, &face.placement))
        }
        _ => None,
    }
}

/// The points whose hull holds a piece, its normal where it is a surface,
/// and its chord where it is a curve.
type Hull = (Vec<Point>, Option<Vector>, Option<Vector>);

fn hull_of(piece: &Piece<'_>, tol: Tolerances) -> Option<Hull> {
    use ogeom_geom::{Curve3d as _, Surface as _};
    // The normal of four corners, as the cross of the diagonals.
    let normal = |[a, b, c, d]: [Point; 4]| (d - a).cross(c - b);
    Some(match piece {
        Piece::Fixed(points) => (points.clone(), None, None),
        Piece::Spline(spline) => {
            let hull: Vec<Point> = spline.control_points().iter().map(|w| w.point()).collect();
            let chord = *hull.last()? - *hull.first()?;
            (hull, None, Some(chord))
        }
        Piece::Arc(curve, range) => {
            let (start, end) = (
                curve.point_at(range.0, tol).ok()?,
                curve.point_at(range.1, tol).ok()?,
            );
            let hull = match curve {
                Curve::Line(_) => vec![start, end],
                _ => crate::measure::curve_bounds_over(curve, *range, tol)
                    .ok()?
                    .corners(),
            };
            (hull, None, Some(end - start))
        }
        Piece::Patch(patch, placement) => {
            let grid = patch.grid();
            let (nu, nv) = (grid.u_count(), grid.v_count());
            let at = |i, j| grid.get(i, j).map(|w| placement.apply(w.point()));
            let corners = [
                at(0, 0)?,
                at(nu - 1, 0)?,
                at(0, nv - 1)?,
                at(nu - 1, nv - 1)?,
            ];
            let hull = grid
                .points()
                .iter()
                .map(|w| placement.apply(w.point()))
                .collect();
            (hull, Some(normal(corners)), None)
        }
        Piece::Flat(plane, u, v, placement) => {
            let mut corners = [Point::ORIGIN; 4];
            for (corner, (a, b)) in
                corners
                    .iter_mut()
                    .zip([(u.0, v.0), (u.1, v.0), (u.0, v.1), (u.1, v.1)])
            {
                *corner = placement.apply(plane.point_at(a, b, tol).ok()?);
            }
            (corners.to_vec(), Some(normal(corners)), None)
        }
    })
}

/// A piece cut into two halves of its parameter range.
fn split<'a>(piece: &Piece<'a>, tol: Tolerances) -> Option<(Piece<'a>, Piece<'a>)> {
    use ogeom_geom::{Curve3d as _, Surface as _};
    let mid = |(a, b): (f64, f64)| 0.5 * (a + b);
    match piece {
        Piece::Fixed(_) => None,
        Piece::Spline(spline) => {
            let (a, b) = spline.domain();
            let m = mid((a, b));
            Some((
                Piece::Spline(Box::new(spline.segment((a, m), tol).ok()?)),
                Piece::Spline(Box::new(spline.segment((m, b), tol).ok()?)),
            ))
        }
        Piece::Arc(curve, (a, b)) => {
            let m = mid((*a, *b));
            Some((Piece::Arc(curve, (*a, m)), Piece::Arc(curve, (m, *b))))
        }
        Piece::Patch(patch, placement) => {
            let (u, v) = patch.domain();
            // Across the direction whose control polygons run longer.
            let grid = patch.grid();
            let length = |along_u: bool| {
                let (across, along) = if along_u {
                    (grid.v_count(), grid.u_count())
                } else {
                    (grid.u_count(), grid.v_count())
                };
                let at = |k: usize, o: usize| {
                    if along_u {
                        grid.get(k, o)
                    } else {
                        grid.get(o, k)
                    }
                };
                let mut longest = 0.0_f64;
                for o in 0..across {
                    let mut run = 0.0;
                    for k in 1..along {
                        if let (Some(p), Some(q)) = (at(k - 1, o), at(k, o)) {
                            run += p.point().distance(q.point());
                        }
                    }
                    longest = longest.max(run);
                }
                longest
            };
            let pieces = if length(true) >= length(false) {
                let m = mid(u);
                (
                    patch.segment((u.0, m), v, tol),
                    patch.segment((m, u.1), v, tol),
                )
            } else {
                let m = mid(v);
                (
                    patch.segment(u, (v.0, m), tol),
                    patch.segment(u, (m, v.1), tol),
                )
            };
            Some((
                Piece::Patch(Box::new(pieces.0.ok()?), placement),
                Piece::Patch(Box::new(pieces.1.ok()?), placement),
            ))
        }
        Piece::Flat(plane, u, v, placement) => {
            if u.1 - u.0 >= v.1 - v.0 {
                let m = mid(*u);
                Some((
                    Piece::Flat(plane, (u.0, m), *v, placement),
                    Piece::Flat(plane, (m, u.1), *v, placement),
                ))
            } else {
                let m = mid(*v);
                Some((
                    Piece::Flat(plane, *u, (v.0, m), placement),
                    Piece::Flat(plane, *u, (m, v.1), placement),
                ))
            }
        }
    }
}

/// [`distance_between_shapes`] on elements already gathered.
pub(crate) fn distance_between_prepared(
    a: &Elements,
    b: &Elements,
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<ShapeDistance> {
    distance_between_elements(&a.0, &b.0, options, tol)
}

fn distance_between_elements(
    ea: &[Element],
    eb: &[Element],
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<ShapeDistance> {
    if ea.is_empty() || eb.is_empty() {
        ogeom_bail!(Construction, "a shape with no elements has no distance");
    }

    let mut candidates: Vec<(f64, ClosestPair)> = Vec::new();
    for element_a in ea {
        for element_b in eb {
            approach(element_a, element_b, options, tol, &mut candidates)?;
        }
    }
    // A candidate that is not a number measures nothing; it neither wins
    // nor, compared equal, lets every other through.
    candidates.retain(|(d, _)| d.is_finite());
    let Some(least) = candidates.iter().map(|(d, _)| *d).min_by(f64::total_cmp) else {
        ogeom_bail!(
            NotDone,
            "no candidate approach was found between these shapes"
        );
    };

    let mut pairs: Vec<ClosestPair> = Vec::new();
    for (d, pair) in candidates {
        if d - least > tol.confusion() {
            continue;
        }
        // The same nearest pair arrives from several element pairs: a corner
        // is on a vertex, three edges and three faces at once. Keep the first
        // at each location.
        if pairs.iter().any(|known| {
            known.point_a.distance(pair.point_a) <= tol.confusion() * 1e2
                && known.point_b.distance(pair.point_b) <= tol.confusion() * 1e2
        }) {
            continue;
        }
        pairs.push(pair);
    }
    Ok(ShapeDistance {
        distance: least,
        pairs,
    })
}

/// Candidate approaches between one pair of elements.
fn approach(
    a: &Element,
    b: &Element,
    options: ExtremaOptions,
    tol: Tolerances,
    out: &mut Vec<(f64, ClosestPair)>,
) -> OgeomResult<()> {
    let mut push = |distance: f64, pa: Point, pb: Point, sa: &Shape, sb: &Shape| {
        out.push((
            distance,
            ClosestPair {
                point_a: pa,
                point_b: pb,
                support_a: sa.clone(),
                support_b: sb.clone(),
            },
        ));
    };
    match (a, b) {
        (Element::Vertex(sa, pa), Element::Vertex(sb, pb)) => {
            push(pa.distance(*pb), *pa, *pb, sa, sb);
        }
        (Element::Vertex(sa, pa), Element::Edge(sb, curve)) => {
            let foot = project_on_curve(curve, *pa, 64, tol)?;
            push(foot.distance, *pa, foot.point, sa, sb);
        }
        (Element::Edge(sa, curve), Element::Vertex(sb, pb)) => {
            let foot = project_on_curve(curve, *pb, 64, tol)?;
            push(foot.distance, foot.point, *pb, sa, sb);
        }
        (Element::Vertex(sa, pa), Element::Face(sb, face)) => {
            let foot = project_on_surface(&face.world, *pa, 32, tol)?;
            if inside_trim(face, foot.point, tol)? {
                push(foot.distance, *pa, foot.point, sa, sb);
            }
        }
        (Element::Face(sa, face), Element::Vertex(sb, pb)) => {
            let foot = project_on_surface(&face.world, *pb, 32, tol)?;
            if inside_trim(face, foot.point, tol)? {
                push(foot.distance, foot.point, *pb, sa, sb);
            }
        }
        (Element::Edge(sa, ca), Element::Edge(sb, cb)) => {
            let found = ogeom_intersect::extrema_curve_curve(ca, cb, options, tol)?;
            for near in &found.approaches {
                push(near.distance, near.point_a, near.point_b, sa, sb);
            }
        }
        (Element::Edge(sa, curve), Element::Face(sb, face)) => {
            let found = ogeom_intersect::extrema_curve_surface(curve, &face.world, options, tol)?;
            for near in &found.approaches {
                if inside_trim(face, near.point_b, tol)? {
                    push(near.distance, near.point_a, near.point_b, sa, sb);
                }
            }
        }
        (Element::Face(sa, face), Element::Edge(sb, curve)) => {
            let found = ogeom_intersect::extrema_curve_surface(curve, &face.world, options, tol)?;
            for near in &found.approaches {
                if inside_trim(face, near.point_b, tol)? {
                    push(near.distance, near.point_b, near.point_a, sa, sb);
                }
            }
        }
        (Element::Face(sa, fa), Element::Face(sb, fb)) => {
            let found =
                ogeom_intersect::extrema_surface_surface(&fa.world, &fb.world, options, tol)?;
            for near in &found.approaches {
                if inside_trim(fa, near.point_a, tol)? && inside_trim(fb, near.point_b, tol)? {
                    push(near.distance, near.point_a, near.point_b, sa, sb);
                }
            }
        }
    }
    Ok(())
}

/// Whether a world-space point on a face's surface lands inside its trimming.
///
/// Asked in the stored surface's own parameter space: the point is carried
/// into the face's frame and projected there, exactly as `classify_on_face`
/// does it, because rings and world-surface parameters need not agree under a
/// placement that scales. Too near the boundary counts as outside: the edge
/// pairs own that candidate and answer it exactly.
fn inside_trim(face: &Prepared, world_point: Point, tol: Tolerances) -> OgeomResult<bool> {
    let local = face.to_local.apply(world_point);
    let projection = project_on_surface(&face.local, local, 32, tol)?;
    let (u, v) = projection.parameters;
    let at = Point2::new(u, v);
    let band = parametric_band(&face.local, (u, v), tol.confusion() + RING_CHORD, tol);
    if band.meets_rings(&face.rings, at) {
        return Ok(false);
    }
    Ok(inside_boundary(&face.rings, at))
}

/// The rings' polylining error, spatially: the trim test's uncertainty band.
const RING_CHORD: f64 = 1e-3;

/// Every element of a shape, with world-space geometry.
fn elements(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<Element>> {
    let mut out = Vec::new();
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let Some(node) = model.node(&vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        let placed = vertex.transform(model.datums())?.apply(data.point);
        out.push(Element::Vertex(vertex, placed));
    }
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(node) = model.node(&edge) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            // A degenerate edge has no extent of its own; its vertex and its
            // face carry its geometry.
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let placement = edge.transform(model.datums())?;
        let trimmed: Curve = if (range.0, range.1) == {
            use ogeom_geom::Curve3d as _;
            geometry.domain()
        } {
            geometry.clone()
        } else {
            TrimmedCurve::new(geometry.clone(), range.0, range.1, tol)?.into()
        };
        out.push(Element::Edge(
            edge,
            Box::new(trimmed.transformed(&placement, tol)?),
        ));
    }
    let ring_deflection = Deflection {
        chord: RING_CHORD,
        angular: 0.05,
        ..Deflection::default()
    };
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        let rings = face_boundary(model, &face, ring_deflection, tol)?;
        // The stored surface may declare an enormous domain (a plane spans
        // ±1e9), and the extrema layer rightly refuses to sample that. The
        // face only uses what its rings enclose, so the surface handed over
        // is trimmed to their parameter bound, with a margin for the rings'
        // own polylining, before being carried into world space.
        let span = ring_span(surface, &rings, tol);
        let restricted = match span {
            Some((u, v)) => ogeom_geom::TrimmedSurface::new(surface.clone(), u, v, tol)?.into(),
            None => surface.clone(),
        };
        out.push(Element::Face(
            face,
            Box::new(Prepared {
                world: restricted.transformed(&placement, tol)?,
                local: surface.clone(),
                to_local: placement.inverse()?,
                rings,
                placement,
                span,
            }),
        ));
    }
    Ok(out)
}

/// The parameter rectangle a face's rings enclose, the surface restricted
/// to it being all of the face there is to measure; `None` for a face with
/// no rings, which uses its surface's whole domain.
///
/// The margin is proportional to the used span: the exact boundary lies
/// within the rings' polylining of it, and `inside_trim` already treats the
/// near-boundary band as the edges' territory, so the margin only has to
/// keep the whole face inside the restriction; it does not have to be
/// tight.
fn ring_span(
    surface: &SurfaceGeometry,
    rings: &[Vec<Point2>],
    tol: Tolerances,
) -> Option<((f64, f64), (f64, f64))> {
    use ogeom_geom::Surface as _;
    let ((ua, ub), (va, vb)) = surface.domain();
    let mut u = (f64::INFINITY, f64::NEG_INFINITY);
    let mut v = (f64::INFINITY, f64::NEG_INFINITY);
    for ring in rings {
        for p in ring {
            u = (u.0.min(p.x), u.1.max(p.x));
            v = (v.0.min(p.y), v.1.max(p.y));
        }
    }
    if u.0 > u.1 || v.0 > v.1 {
        return None;
    }
    let margin_u = (u.1 - u.0).mul_add(0.05, tol.parametric());
    let margin_v = (v.1 - v.0).mul_add(0.05, tol.parametric());
    let lo_u = (u.0 - margin_u).max(ua);
    let hi_u = (u.1 + margin_u).min(ub);
    let lo_v = (v.0 - margin_v).max(va);
    let hi_v = (v.1 + margin_v).min(vb);
    (lo_u < hi_u && lo_v < hi_v).then_some(((lo_u, hi_u), (lo_v, hi_v)))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{make_box, make_cylinder, make_sphere};
    use ogeom_math::{Direction, Frame, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn frame_at(origin: Point) -> Frame {
        Frame::new(origin, Direction::Z, Direction::X, T).unwrap()
    }

    #[test]
    fn parallel_box_walls_meet_at_the_gap_between_them() {
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(5.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        assert!((found.distance - 3.0).abs() < 1e-9, "{}", found.distance);
        assert!(!found.pairs.is_empty());
        for pair in &found.pairs {
            assert!((pair.point_a.distance(pair.point_b) - found.distance).abs() < 1e-9);
        }
    }

    #[test]
    fn diagonal_boxes_meet_corner_to_corner() {
        // Offset along all three axes: the nearest points are two vertices,
        // exactly the candidates the geometry level declines to invent.
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(3.0, 3.0, 3.0)),
            (1.0, 1.0, 1.0),
            T,
        )
        .unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        let exact = (3.0_f64 * 4.0).sqrt(); // corner (1,1,1) to corner (3,3,3)
        assert!((found.distance - exact).abs() < 1e-9);
        let pair = &found.pairs[0];
        assert!(pair.point_a.is_equal(Point::new(1.0, 1.0, 1.0), T));
        assert!(pair.point_b.is_equal(Point::new(3.0, 3.0, 3.0), T));
        assert_eq!(model.kind_of(&pair.support_a).unwrap(), ShapeType::Vertex);
        assert_eq!(model.kind_of(&pair.support_b).unwrap(), ShapeType::Vertex);
    }

    #[test]
    fn a_sphere_over_a_box_measures_to_the_top_face() {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let ball = make_sphere(&mut model, frame_at(Point::new(2.0, 2.0, 4.0)), 1.0, T).unwrap();
        let found = distance_between_shapes(
            &model,
            &block.shape,
            &ball.shape,
            ExtremaOptions::default(),
            T,
        )
        .unwrap();
        assert!((found.distance - 2.0).abs() < 1e-7, "{}", found.distance);
        let pair = &found.pairs[0];
        assert!(pair.point_a.is_equal(Point::new(2.0, 2.0, 1.0), T));
        assert!(pair.point_b.is_equal(Point::new(2.0, 2.0, 3.0), T));
    }

    #[test]
    fn parallel_cylinders_meet_wall_to_wall() {
        // The nearest locus is a pair of facing rulings: a family at the
        // geometry level, a representative pair here, with the distance exact.
        let mut model = Model::new();
        let a = make_cylinder(&mut model, Frame::WORLD, 1.0, 4.0, T).unwrap();
        let b =
            make_cylinder(&mut model, frame_at(Point::new(5.0, 0.0, 0.0)), 1.0, 4.0, T).unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        assert!((found.distance - 3.0).abs() < 1e-7, "{}", found.distance);
    }

    #[test]
    fn touching_boxes_report_zero() {
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(2.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        assert!(found.distance < 1e-9, "{}", found.distance);
    }

    #[test]
    fn a_box_inside_a_box_measures_boundary_to_boundary() {
        // Containment is the classifier's question. Distance is between
        // boundaries, and the gap between nested walls is what comes back.
        let mut model = Model::new();
        let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
        let inner = make_box(
            &mut model,
            frame_at(Point::new(2.0, 2.0, 2.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let found = distance_between_shapes(
            &model,
            &outer.shape,
            &inner.shape,
            ExtremaOptions::default(),
            T,
        )
        .unwrap();
        assert!((found.distance - 2.0).abs() < 1e-9, "{}", found.distance);
    }

    #[test]
    fn a_rotated_box_measures_edge_to_edge() {
        // Roll one box forty-five degrees about x and lift it: what faces the
        // top of the lower box is a single edge, and the nearest pair is that
        // edge against the top face.
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let tilted = Frame::new(
            Point::new(2.0, 2.0, 3.0),
            Direction::new(Vector::new(0.0, 1.0, 1.0), T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        let b = make_box(&mut model, tilted, (1.0, 1.0, 1.0), T).unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        // The tilted box's lowest feature is the edge its roll brings down to
        // z = 3 - 1/sqrt(2), facing the top face at z = 1.
        let exact = 2.0 - core::f64::consts::FRAC_1_SQRT_2;
        assert!((found.distance - exact).abs() < 1e-7, "{}", found.distance);
    }
}
