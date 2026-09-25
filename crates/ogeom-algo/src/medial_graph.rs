//! The medial axis of any planar face: holes, reflex corners and curved
//! boundaries included, each branch its exact bisector.
//!
//! The boundary is taken apart into *sites*: every edge's open interior,
//! and every reflex corner as a point (a convex corner is where its two
//! edges' bisector reaches the boundary, not a site of its own). The axis
//! is the set of points with two or more nearest sites, and where two
//! sites meet it is a curve with a closed form: two segments' bisector is a
//! line, a segment and a point's a parabola, and a circular arc speaks as
//! a circle, so an arc against a segment is a parabola too and an arc
//! against a point or an arc is a hyperbola or an ellipse. A site on any
//! other curve has no closed form, and its branches are fitted through
//! points solved onto the bisector exactly.
//!
//! Which sites meet where (the axis's topology) is read off a sampling:
//! the Delaunay triangulation of points along every site has, inside the
//! region, circumcentres whose three corners name their sites. A centre
//! naming two sites is a point on their branch, one naming three a branch
//! point. Each branch point is then solved exactly (equidistant from its
//! sites, by Newton), and each branch is its exact bisector between the
//! branch points the sampling says it runs between. The result is checked:
//! along every branch no third site may come nearer than the branch's own
//! two, and every sampled centre must lie on some branch. A result that
//! fails is sampled again, finer.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point, Point2, Vector2};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape};
use spade::{DelaunayTriangulation, Point2 as SpadePoint, Triangulation as _};
use std::collections::HashMap;

/// A boundary element a medial branch keeps its distance from.
#[derive(Debug, Clone)]
pub enum MedialSite {
    /// The open interior of one curve along the boundary: an edge, or
    /// neighbouring edges on one line or one circle, which the boundary
    /// splits and the axis does not see.
    Edges(Vec<Shape>),
    /// A reflex corner, as a point.
    Vertex(Shape),
}

/// A branch point of the medial axis, or a branch's end on the boundary.
#[derive(Debug, Clone)]
pub struct MedialVertex {
    /// Where it stands.
    pub point: Point,
    /// The radius of the largest disc centred here inside the face: the
    /// distance to the nearest boundary. Zero where a branch reaches a
    /// convex corner.
    pub clearance: f64,
}

/// One branch: the curve of points equidistant from two sites.
#[derive(Debug, Clone)]
pub struct MedialBranch {
    /// The branch's curve, exact wherever the two sites have a closed-form
    /// bisector (a line, a parabola, a hyperbola or an ellipse), fitted
    /// otherwise.
    pub curve: ogeom_geom::Curve,
    /// The portion of `curve` the branch covers.
    pub range: (f64, f64),
    /// The vertices at `range.0` and `range.1`, by index.
    pub ends: [usize; 2],
    /// The two sites the branch is equidistant from.
    pub sites: [MedialSite; 2],
}

/// The medial axis of a planar face as a graph of exact branches.
#[derive(Debug, Clone)]
pub struct MedialGraph {
    /// Branch points and boundary ends.
    pub vertices: Vec<MedialVertex>,
    /// The branches between them.
    pub branches: Vec<MedialBranch>,
    /// The largest amount, measured along the branches, by which a site
    /// other than a branch's own two comes nearer than they do: zero for a
    /// branch that is all medial, and held under the tolerance asked for.
    pub deviation: f64,
    frame: Frame,
    sites: Vec<Site>,
    branch_sites: Vec<[usize; 2]>,
}

impl MedialGraph {
    /// The clearance at parameter `t` of branch `branch`: its distance from
    /// either of its two sites.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// the branch does not exist or the curve cannot be evaluated there.
    pub fn clearance_at(&self, branch: usize, t: f64, tol: Tolerances) -> OgeomResult<f64> {
        let Some(b) = self.branches.get(branch) else {
            ogeom_bail!(Construction, "no branch {branch}");
        };
        let p = self.frame.to_local(b.curve.point_at(t, tol)?);
        let site = &self.sites[self.branch_sites[branch][0]];
        site.distance(Point2::new(p.x, p.y), tol)
    }
}

/// How many times the sampling may be refined before the construction
/// gives up.
const REFINEMENTS: usize = 4;

/// The medial axis of a planar face: every point with two or more nearest
/// boundary elements, as exact branches between branch points.
///
/// Holes, reflex corners and circular arcs are handled exactly; an edge on
/// any other curve bisects along a fitted branch held to `tolerance`. The
/// result's own check is reported in [`MedialGraph::deviation`].
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// face is not planar or its boundary cannot be read;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the branches
/// still fail their check after every refinement of the sampling.
pub fn medial_graph(
    model: &Model,
    face: &Shape,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<MedialGraph> {
    if !tolerance.is_finite() || tolerance <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a medial axis to {tolerance} is not a distance"
        );
    }
    let boundary = Boundary::read(model, face, tol)?;
    let mut spacing = boundary.extent / 256.0;
    let mut last = None;
    for _ in 0..REFINEMENTS {
        match build(&boundary, spacing, tolerance, tol) {
            Ok(graph) => return Ok(graph),
            Err(e) => last = Some(e),
        }
        spacing *= 0.5;
    }
    ogeom_bail!(
        NotDone,
        "the medial axis did not settle under refinement: {}",
        last.map_or_else(String::new, |e| e.to_string())
    )
}

// --- the boundary, in the face's own plane ------------------------------

#[derive(Debug, Clone)]
enum SiteKind {
    Segment {
        a: Point2,
        b: Point2,
    },
    /// A circular arc: centre, radius, start angle and signed sweep.
    Arc {
        centre: Point2,
        radius: f64,
        start: f64,
        sweep: f64,
    },
    Point {
        at: Point2,
    },
    /// Any other curve, placed, over its range in the direction of travel,
    /// with a dense polyline for seeding.
    Curve {
        curve: Box<ogeom_geom::Curve>,
        range: (f64, f64),
        reversed: bool,
        polyline: Vec<(f64, Point2)>,
        frame: Frame,
    },
}

#[derive(Debug, Clone)]
struct Site {
    kind: SiteKind,
    origin: MedialSite,
}

impl Site {
    fn foot(&self, p: Point2, tol: Tolerances) -> OgeomResult<Point2> {
        Ok(match &self.kind {
            SiteKind::Segment { a, b } => {
                let d = *b - *a;
                let s = ((p - *a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
                *a + d * s
            }
            SiteKind::Arc {
                centre,
                radius,
                start,
                sweep,
            } => {
                let v = p - *centre;
                let theta = v.y.atan2(v.x);
                let within = if *sweep >= 0.0 {
                    (theta - start).rem_euclid(core::f64::consts::TAU) <= *sweep
                } else {
                    (start - theta).rem_euclid(core::f64::consts::TAU) <= -sweep
                };
                let on = |angle: f64| *centre + Vector2::new(angle.cos(), angle.sin()) * *radius;
                if within && v.magnitude() > 0.0 {
                    on(theta)
                } else {
                    let (s, e) = (on(*start), on(start + sweep));
                    if s.distance(p) <= e.distance(p) { s } else { e }
                }
            }
            SiteKind::Point { at } => *at,
            SiteKind::Curve {
                curve,
                range,
                polyline,
                frame,
                ..
            } => {
                let seed = polyline
                    .iter()
                    .min_by(|x, y| x.1.distance(p).total_cmp(&y.1.distance(p)))
                    .map_or(range.0, |s| s.0);
                let lifted = lift(frame, p);
                let mut t = seed;
                for _ in 0..30 {
                    let d = curve.d1_at(t, tol)?;
                    let q = curve.point_at(t, tol)?;
                    let step = (lifted - q).dot(d) / d.dot(d).max(1e-300);
                    let next = (t + step).clamp(range.0.min(range.1), range.0.max(range.1));
                    if (next - t).abs() <= 1e-14 * (1.0 + t.abs()) {
                        t = next;
                        break;
                    }
                    t = next;
                }
                flat(frame, curve.point_at(t, tol)?)
            }
        })
    }

    fn distance(&self, p: Point2, tol: Tolerances) -> OgeomResult<f64> {
        Ok(self.foot(p, tol)?.distance(p))
    }

    /// Points along the site at about `spacing` apart, its ends left out:
    /// a site is its open interior.
    fn samples(&self, spacing: f64, tol: Tolerances) -> OgeomResult<Vec<Point2>> {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = |length: f64| (length / spacing).ceil().max(4.0) as usize;
        Ok(match &self.kind {
            SiteKind::Segment { a, b } => {
                let n = count(a.distance(*b));
                (0..n)
                    .map(|i| {
                        #[allow(clippy::cast_precision_loss)]
                        let f = (i as f64 + 0.5) / n as f64;
                        a.lerp(*b, f)
                    })
                    .collect()
            }
            SiteKind::Arc {
                centre,
                radius,
                start,
                sweep,
            } => {
                let n = count(radius * sweep.abs());
                (0..n)
                    .map(|i| {
                        #[allow(clippy::cast_precision_loss)]
                        let angle = start + sweep * (i as f64 + 0.5) / n as f64;
                        *centre + Vector2::new(angle.cos(), angle.sin()) * *radius
                    })
                    .collect()
            }
            SiteKind::Point { at } => vec![*at],
            SiteKind::Curve {
                curve,
                range,
                polyline,
                frame,
                ..
            } => {
                let length: f64 = polyline.windows(2).map(|w| w[0].1.distance(w[1].1)).sum();
                let n = count(length);
                let mut out = Vec::with_capacity(n);
                for i in 0..n {
                    #[allow(clippy::cast_precision_loss)]
                    let f = (i as f64 + 0.5) / n as f64;
                    let t = range.0 + (range.1 - range.0) * f;
                    out.push(flat(frame, curve.point_at(t, tol)?));
                }
                out
            }
        })
    }

    /// The site as the circle or line its exact bisectors speak: a segment
    /// by its line and inward normal, an arc by its circle and the side the
    /// region lies on, a point as a circle of no radius.
    fn primitive(&self) -> Option<Primitive> {
        match &self.kind {
            SiteKind::Segment { a, b } => {
                let d = *b - *a;
                let m = d.magnitude();
                // Inward is to the left: the boundary is wound so.
                Some(Primitive::Line {
                    normal: Vector2::new(-d.y / m, d.x / m),
                    through: *a,
                })
            }
            SiteKind::Arc {
                centre,
                radius,
                sweep,
                ..
            } => Some(Primitive::Circle {
                centre: *centre,
                radius: *radius,
                // Turning left the arc holds the region inside its circle.
                outside: *sweep < 0.0,
            }),
            SiteKind::Point { at } => Some(Primitive::Circle {
                centre: *at,
                radius: 0.0,
                outside: true,
            }),
            SiteKind::Curve { .. } => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Primitive {
    Line {
        normal: Vector2,
        through: Point2,
    },
    Circle {
        centre: Point2,
        radius: f64,
        outside: bool,
    },
}

/// The face's boundary as sites, with its corners and a fine outline.
struct Boundary {
    frame: Frame,
    sites: Vec<Site>,
    /// Convex corners: the point and the two sites meeting there.
    convex: Vec<(Point2, usize, usize)>,
    /// Reflex corners: the point's own site and its two edges' sites.
    reflex: Vec<(usize, usize, usize)>,
    /// Every loop as a fine polyline, for the inside test.
    outline: Vec<Vec<Point2>>,
    extent: f64,
}

fn lift(frame: &Frame, p: Point2) -> Point {
    frame.origin() + frame.x().vector() * p.x + frame.y().vector() * p.y
}

fn flat(frame: &Frame, p: Point) -> Point2 {
    let l = frame.to_local(p);
    Point2::new(l.x, l.y)
}

/// One edge as travelled round its loop: its site and its tangents in and
/// out, with points along it for the outline.
struct Travel {
    site: SiteKind,
    edges: Vec<Shape>,
    /// The vertices the travel leaves and arrives at.
    from: Shape,
    to: Shape,
    enter: Vector2,
    leave: Vector2,
    start: Point2,
    points: Vec<Point2>,
}

impl Boundary {
    fn read(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Self> {
        let Some(data) = model.node(face).and_then(|n| match n.data() {
            NodeData::Face(d) => Some(d.clone()),
            _ => None,
        }) else {
            ogeom_bail!(Construction, "the shape is not a face");
        };
        let Some(ogeom_geom::SurfaceGeometry::Plane(plane)) =
            model.geometry().surface(data.surface)
        else {
            ogeom_bail!(
                Construction,
                "the medial axis is computed for planar faces; this face's \
                 surface is not a plane"
            );
        };
        let placement = face.transform(model.datums())?;
        let frame = {
            let f = plane.plane().frame();
            Frame::new(
                placement.apply(f.origin()),
                Direction::new(placement.apply_vector(f.z().vector()), tol)?,
                Direction::new(placement.apply_vector(f.x().vector()), tol)?,
                tol,
            )?
        };

        let mut loops: Vec<Vec<Travel>> = Vec::new();
        for wire in model.ordered_children_of(face)? {
            let mut travel = Vec::new();
            for edge in model.ordered_children_of(&wire)? {
                travel.push(read_edge(model, &edge, &frame, tol)?);
            }
            if !travel.is_empty() {
                loops.push(travel);
            }
        }
        if loops.is_empty() {
            ogeom_bail!(Construction, "the face has no boundary to read");
        }
        // Wind the loops so the region lies to the left of every one: the
        // outer loop counter-clockwise and each hole clockwise, whatever
        // sense each was stored in. The outer is the loop of largest area.
        let areas: Vec<f64> = loops
            .iter()
            .map(|l| signed_area(&l.iter().flat_map(|t| t.points.clone()).collect::<Vec<_>>()))
            .collect();
        let outer = (0..loops.len())
            .max_by(|&a, &b| areas[a].abs().total_cmp(&areas[b].abs()))
            .unwrap_or(0);
        for (i, l) in loops.iter_mut().enumerate() {
            let wound_wrong = if i == outer {
                areas[i] < 0.0
            } else {
                areas[i] > 0.0
            };
            if wound_wrong {
                l.reverse();
                for t in l.iter_mut() {
                    reverse_travel(t);
                }
            }
        }

        // Neighbours on one line or one circle are one site: the boundary
        // split them, and a join between them is no corner.
        for l in &mut loops {
            merge_continuations(l, tol);
        }

        let mut sites: Vec<Site> = Vec::new();
        let mut convex = Vec::new();
        let mut reflex = Vec::new();
        let mut outline = Vec::new();
        let mut extent = 0.0_f64;
        let mut all: Vec<Point2> = Vec::new();
        for l in &loops {
            let first = sites.len();
            for t in l {
                sites.push(Site {
                    kind: t.site.clone(),
                    origin: MedialSite::Edges(t.edges.clone()),
                });
            }
            let n = l.len();
            for i in 0..n {
                let j = (i + 1) % n;
                let (into, out) = (l[i].leave, l[j].enter);
                let turn = into.cross(out);
                let at = l[j].start;
                if turn > tol.angular() {
                    convex.push((at, first + i, first + j));
                } else if turn < -tol.angular() || into.dot(out) < 0.0 {
                    let vertex = l[j].from.clone();
                    let point_site = sites.len();
                    sites.push(Site {
                        kind: SiteKind::Point { at },
                        origin: MedialSite::Vertex(vertex),
                    });
                    reflex.push((point_site, first + i, first + j));
                }
            }
            let ring: Vec<Point2> = l.iter().flat_map(|t| t.points.clone()).collect();
            all.extend(ring.iter().copied());
            outline.push(ring);
        }
        if let (Some(lo), Some(hi)) = (
            all.iter()
                .copied()
                .reduce(|a, b| Point2::new(a.x.min(b.x), a.y.min(b.y))),
            all.iter()
                .copied()
                .reduce(|a, b| Point2::new(a.x.max(b.x), a.y.max(b.y))),
        ) {
            extent = lo.distance(hi);
        }
        if extent <= tol.confusion() {
            ogeom_bail!(Construction, "the face has no extent");
        }
        Ok(Self {
            frame,
            sites,
            convex,
            reflex,
            outline,
            extent,
        })
    }

    fn inside(&self, p: Point2) -> bool {
        let mut inside = false;
        for ring in &self.outline {
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                if (a.y > p.y) != (b.y > p.y) {
                    let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                    if x > p.x {
                        inside = !inside;
                    }
                }
            }
        }
        inside
    }
}

/// The traversal start vertex of an edge (its first vertex as used), and
/// the vertex it arrives at.
fn edge_ends(model: &Model, edge: &Shape) -> OgeomResult<Option<(Shape, Shape)>> {
    let Some((a, b)) = crate::edge_vertices(model, edge)? else {
        return Ok(None);
    };
    Ok(Some(if edge.orientation() == Orientation::Reversed {
        (b, a)
    } else {
        (a, b)
    }))
}

fn reverse_travel(t: &mut Travel) {
    t.points.reverse();
    let (enter, leave) = (-t.leave, -t.enter);
    t.enter = enter;
    t.leave = leave;
    t.start = t.points[0];
    t.site = match &t.site {
        SiteKind::Segment { a, b } => SiteKind::Segment { a: *b, b: *a },
        SiteKind::Arc {
            centre,
            radius,
            start,
            sweep,
        } => SiteKind::Arc {
            centre: *centre,
            radius: *radius,
            start: start + sweep,
            sweep: -sweep,
        },
        SiteKind::Point { at } => SiteKind::Point { at: *at },
        SiteKind::Curve {
            curve,
            range,
            reversed,
            polyline,
            frame,
        } => SiteKind::Curve {
            curve: curve.clone(),
            range: (range.1, range.0),
            reversed: !reversed,
            polyline: polyline.iter().rev().copied().collect(),
            frame: *frame,
        },
    };
    core::mem::swap(&mut t.from, &mut t.to);
}

fn read_edge(model: &Model, edge: &Shape, frame: &Frame, tol: Tolerances) -> OgeomResult<Travel> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Construction, "a boundary edge holds no data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "a boundary edge carries no curve");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Construction, "a boundary curve is not in this model");
    };
    use ogeom_geom::Transformable as _;
    let placed = geometry
        .clone()
        .transformed(&edge.transform(model.datums())?, tol)?;
    let reversed = edge.orientation() == Orientation::Reversed;
    let (t0, t1) = if reversed {
        (range.1, range.0)
    } else {
        (range.0, range.1)
    };
    let at = |t: f64| -> OgeomResult<Point2> { Ok(flat(frame, placed.point_at(t, tol)?)) };
    let heading = |t: f64| -> OgeomResult<Vector2> {
        let d = placed.d1_at(t, tol)?;
        let l = frame.to_local(frame.origin() + d);
        let v = Vector2::new(l.x, l.y);
        let v = v.normalized(tol)?;
        Ok(if reversed { -v } else { v })
    };
    let mut points = Vec::with_capacity(65);
    for i in 0..64 {
        let t = t0 + (t1 - t0) * f64::from(i) / 64.0;
        points.push(at(t)?);
    }
    let (start, end) = (at(t0)?, at(t1)?);
    let site = match &placed {
        ogeom_geom::Curve::Line(_) => SiteKind::Segment { a: start, b: end },
        ogeom_geom::Curve::Circle(c) => {
            let centre = flat(frame, c.circle().frame().origin());
            let radius = c.circle().radius();
            let angle = |p: Point2| (p.y - centre.y).atan2(p.x - centre.x);
            let a0 = angle(start);
            let mid = at(f64::midpoint(t0, t1))?;
            // The sweep from start through the midpoint to the end: turning
            // left is positive.
            let turning = (mid - start).cross(end - mid);
            let travel = (angle(end) - a0).rem_euclid(core::f64::consts::TAU);
            let closed = start.distance(end) <= tol.confusion();
            let left = if closed {
                heading(t0)?.cross(start - centre) < 0.0
                    || (start - centre).cross(heading(t0)?) > 0.0
            } else {
                turning > 0.0
                    || (turning.abs() <= tol.confusion() && (mid - centre).cross(end - start) > 0.0)
            };
            let sweep = if closed {
                if left {
                    core::f64::consts::TAU
                } else {
                    -core::f64::consts::TAU
                }
            } else if left {
                if travel <= 0.0 {
                    core::f64::consts::TAU
                } else {
                    travel
                }
            } else {
                -(core::f64::consts::TAU - travel)
            };
            SiteKind::Arc {
                centre,
                radius,
                start: a0,
                sweep,
            }
        }
        other => {
            let mut polyline = Vec::with_capacity(129);
            for i in 0..=128 {
                let t = t0 + (t1 - t0) * f64::from(i) / 128.0;
                polyline.push((t, at(t)?));
            }
            SiteKind::Curve {
                curve: Box::new(other.clone()),
                range: (t0, t1),
                reversed,
                polyline,
                frame: *frame,
            }
        }
    };
    let Some((from, to)) = edge_ends(model, edge)? else {
        ogeom_bail!(Construction, "a boundary edge has no vertices");
    };
    Ok(Travel {
        site,
        edges: vec![edge.clone()],
        from,
        to,
        enter: heading(t0)?,
        leave: heading(t1)?,
        start,
        points,
    })
}

/// Fuse neighbouring travels that continue one another: two segments on
/// one line, or two arcs on one circle turning the same way. A loop that
/// is one circle all round becomes a single full-turn arc.
fn merge_continuations(l: &mut Vec<Travel>, tol: Tolerances) {
    let same = |a: &SiteKind, b: &SiteKind| -> Option<SiteKind> {
        match (a, b) {
            (SiteKind::Segment { a: p, b: q }, SiteKind::Segment { a: r, b: t }) => {
                let d1 = (*q - *p).normalized(tol).ok()?;
                let d2 = (*t - *r).normalized(tol).ok()?;
                (d1.cross(d2).abs() <= tol.angular()
                    && d1.dot(d2) > 0.0
                    && q.distance(*r) <= tol.confusion())
                .then_some(SiteKind::Segment { a: *p, b: *t })
            }
            (
                SiteKind::Arc {
                    centre: c1,
                    radius: r1,
                    start,
                    sweep: s1,
                },
                SiteKind::Arc {
                    centre: c2,
                    radius: r2,
                    sweep: s2,
                    ..
                },
            ) => (c1.distance(*c2) <= tol.confusion()
                && (r1 - r2).abs() <= tol.confusion()
                && s1.signum() == s2.signum())
            .then(|| SiteKind::Arc {
                centre: *c1,
                radius: *r1,
                start: *start,
                sweep: (s1 + s2).clamp(-core::f64::consts::TAU, core::f64::consts::TAU),
            }),
            _ => None,
        }
    };
    let mut i = 0;
    while l.len() > 1 && i < l.len() {
        let j = (i + 1) % l.len();
        if let Some(kind) = same(&l[i].site, &l[j].site) {
            let next = l.remove(j);
            let i = if j < i { i - 1 } else { i };
            let t = &mut l[i];
            t.site = kind;
            t.edges.extend(next.edges);
            t.leave = next.leave;
            t.to = next.to;
            t.points.extend(next.points);
            continue;
        }
        i += 1;
    }
}

fn signed_area(ring: &[Point2]) -> f64 {
    let mut sum = 0.0;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        sum += a.x * b.y - b.x * a.y;
    }
    sum / 2.0
}

// --- exact bisectors ------------------------------------------------------

/// Where two sites' distances agree, as a curve with a parameter.
#[derive(Debug, Clone, Copy)]
enum Bisector {
    Line {
        origin: Point2,
        dir: Vector2,
    },
    Parabola {
        apex: Point2,
        x: Vector2,
        focal: f64,
    },
    Hyperbola {
        centre: Point2,
        x: Vector2,
        a: f64,
        b: f64,
    },
    Ellipse {
        centre: Point2,
        x: Vector2,
        a: f64,
        b: f64,
    },
}

fn perp(v: Vector2) -> Vector2 {
    Vector2::new(-v.y, v.x)
}

impl Bisector {
    fn between(p: Primitive, q: Primitive) -> Option<Self> {
        const EPS: f64 = 1e-12;
        match (p, q) {
            (
                Primitive::Line {
                    normal: n1,
                    through: a1,
                },
                Primitive::Line {
                    normal: n2,
                    through: a2,
                },
            ) => {
                let m = n1 - n2;
                let mm = m.dot(m);
                if mm <= EPS {
                    return None;
                }
                let k = n1.dot(a1.to_vector()) - n2.dot(a2.to_vector());
                Some(Self::Line {
                    origin: Point2::ORIGIN + m * (k / mm),
                    dir: perp(m) * (1.0 / mm.sqrt()),
                })
            }
            (
                Primitive::Line { normal, through },
                Primitive::Circle {
                    centre,
                    radius,
                    outside,
                },
            )
            | (
                Primitive::Circle {
                    centre,
                    radius,
                    outside,
                },
                Primitive::Line { normal, through },
            ) => {
                // n·(p − a) = σ(|p − c| − r)  ⇔  |p − c| = σ n·(p − a*),
                // a* = a − σ r n: focus c, directrix through a*.
                let sigma = if outside { 1.0 } else { -1.0 };
                let x = normal * sigma;
                let on_directrix = through - normal * (sigma * radius);
                let focal = x.dot(centre - on_directrix) / 2.0;
                if focal <= EPS {
                    return Some(Self::Line {
                        origin: centre,
                        dir: x,
                    });
                }
                Some(Self::Parabola {
                    apex: centre - x * focal,
                    x,
                    focal,
                })
            }
            (
                Primitive::Circle {
                    centre: c1,
                    radius: r1,
                    outside: o1,
                },
                Primitive::Circle {
                    centre: c2,
                    radius: r2,
                    outside: o2,
                },
            ) => {
                let mid = c1.lerp(c2, 0.5);
                let span = c2 - c1;
                let e = span.magnitude() / 2.0;
                if o1 == o2 {
                    // |p − c1| − |p − c2| = r1 − r2.
                    let delta = r1 - r2;
                    if delta.abs() <= EPS {
                        if e <= EPS {
                            return None;
                        }
                        return Some(Self::Line {
                            origin: mid,
                            dir: perp(span) * (1.0 / span.magnitude()),
                        });
                    }
                    let a = delta.abs() / 2.0;
                    if e <= a + EPS {
                        return None;
                    }
                    let toward = if delta > 0.0 { span } else { -span };
                    Some(Self::Hyperbola {
                        centre: mid,
                        x: toward * (1.0 / toward.magnitude()),
                        a,
                        b: (e * e - a * a).sqrt(),
                    })
                } else {
                    // |p − c1| + |p − c2| = r1 + r2.
                    let a = (r1 + r2) / 2.0;
                    if a <= e + EPS {
                        return None;
                    }
                    let x = if e > EPS {
                        span * (1.0 / span.magnitude())
                    } else {
                        Vector2::new(1.0, 0.0)
                    };
                    Some(Self::Ellipse {
                        centre: mid,
                        x,
                        a,
                        b: (a * a - e * e).sqrt(),
                    })
                }
            }
        }
    }

    fn param(&self, p: Point2) -> f64 {
        match *self {
            Self::Line { origin, dir } => (p - origin).dot(dir),
            Self::Parabola { apex, x, .. } => (p - apex).dot(perp(x)),
            Self::Hyperbola { centre, x, b, .. } => ((p - centre).dot(perp(x)) / b).asinh(),
            Self::Ellipse { centre, x, a, b } => {
                let d = p - centre;
                (d.dot(perp(x)) / b).atan2(d.dot(x) / a)
            }
        }
    }

    fn at(&self, t: f64) -> Point2 {
        match *self {
            Self::Line { origin, dir } => origin + dir * t,
            Self::Parabola { apex, x, focal } => apex + x * (t * t / (4.0 * focal)) + perp(x) * t,
            Self::Hyperbola { centre, x, a, b } => {
                centre + x * (a * t.cosh()) + perp(x) * (b * t.sinh())
            }
            Self::Ellipse { centre, x, a, b } => {
                centre + x * (a * t.cos()) + perp(x) * (b * t.sin())
            }
        }
    }

    fn periodic(&self) -> bool {
        matches!(self, Self::Ellipse { .. })
    }

    fn curve(
        &self,
        frame: &Frame,
        range: (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<ogeom_geom::Curve> {
        let plane_frame = |origin: Point2, x: Vector2| -> OgeomResult<Frame> {
            Frame::new(
                lift(frame, origin),
                frame.z(),
                Direction::new(frame.x().vector() * x.x + frame.y().vector() * x.y, tol)?,
                tol,
            )
        };
        Ok(match *self {
            Self::Line { .. } => ogeom_geom::LineCurve::segment(
                lift(frame, self.at(range.0)),
                lift(frame, self.at(range.1)),
                tol,
            )?
            .into(),
            Self::Parabola { apex, x, focal } => ogeom_geom::ParabolaCurve::over(
                ogeom_math::Parabola::new(plane_frame(apex, x)?, focal, tol)?,
                range.0,
                range.1,
            )?
            .into(),
            Self::Hyperbola { centre, x, a, b } => ogeom_geom::HyperbolaCurve::over(
                ogeom_math::Hyperbola::new(plane_frame(centre, x)?, a, b, tol)?,
                range.0,
                range.1,
            )?
            .into(),
            Self::Ellipse { centre, x, a, b } => ogeom_geom::EllipseCurve::new(
                ogeom_math::Ellipse::new(plane_frame(centre, x)?, a, b, tol)?,
            )
            .into(),
        })
    }
}

// --- the construction -----------------------------------------------------

/// A branch point being found: the sites it is equidistant from and where.
struct Node {
    at: Point2,
    sites: Vec<usize>,
    boundary: bool,
}

#[allow(
    clippy::too_many_lines,
    reason = "one construction, read top to bottom"
)]
fn build(
    boundary: &Boundary,
    spacing: f64,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<MedialGraph> {
    let sites = &boundary.sites;
    // Sample every site, and triangulate the samples.
    let mut owner: Vec<usize> = Vec::new();
    let mut dt: DelaunayTriangulation<SpadePoint<f64>> = DelaunayTriangulation::new();
    let mut index_of: HashMap<spade::handles::FixedVertexHandle, usize> = HashMap::new();
    for (s, site) in sites.iter().enumerate() {
        for p in site.samples(spacing, tol)? {
            let handle = dt.insert(SpadePoint::new(p.x, p.y)).map_err(|e| {
                ogeom_core::ogeom_err!(Construction, "sampling the boundary: {e:?}")
            })?;
            index_of.insert(handle, owner.len());
            owner.push(s);
        }
    }

    // Circumcentres inside the region, labelled by the sites they touch.
    let mut pair_samples: HashMap<(usize, usize), Vec<Point2>> = HashMap::new();
    let mut triples: Vec<(Point2, Vec<usize>)> = Vec::new();
    for face in dt.inner_faces() {
        let corners = face.vertices();
        let pts: Vec<Point2> = corners
            .iter()
            .map(|v| Point2::new(v.position().x, v.position().y))
            .collect();
        let Some(centre) = circumcentre(pts[0], pts[1], pts[2]) else {
            continue;
        };
        if !boundary.inside(centre) {
            continue;
        }
        let mut labels: Vec<usize> = corners
            .iter()
            .filter_map(|v| index_of.get(&v.fix()).map(|&i| owner[i]))
            .collect();
        labels.sort_unstable();
        labels.dedup();
        // Two sites meeting at one foot (a smooth join) are one boundary
        // there, not two: drop the pair.
        let mut distinct: Vec<usize> = Vec::new();
        for &s in &labels {
            let fs = sites[s].foot(centre, tol)?;
            let mut keep = true;
            for &t in &distinct {
                if sites[t].foot(centre, tol)?.distance(fs) <= spacing * 2.0 {
                    keep = false;
                }
            }
            if keep {
                distinct.push(s);
            }
        }
        // A branch point keeps every site it touches, feet shared or not:
        // where a reflex corner's point and its own edge share a foot, the
        // axis hands over from the one's branch to the other's there.
        if labels.len() >= 3 && distinct.len() >= 2 {
            triples.push((centre, labels));
        } else if distinct.len() == 2 {
            pair_samples
                .entry((distinct[0], distinct[1]))
                .or_default()
                .push(centre);
        }
    }

    // Branch points: nearby triples merged, then solved exactly.
    let mut nodes: Vec<Node> = Vec::new();
    let merge = spacing * 4.0;
    for (at, labels) in triples {
        if let Some(node) = nodes
            .iter_mut()
            .find(|n| !n.boundary && n.at.distance(at) <= merge)
        {
            for s in labels {
                if !node.sites.contains(&s) {
                    node.sites.push(s);
                }
            }
            continue;
        }
        nodes.push(Node {
            at,
            sites: labels,
            boundary: false,
        });
    }
    for node in &mut nodes {
        // A branch point touching an arc the region lies inside may be the
        // arc's own centre, where the arc's distance has no gradient for
        // Newton to follow: tried first, and kept where every site stands
        // at one distance from it.
        let centre_of = node.sites.iter().find_map(|&s| match sites[s].kind {
            SiteKind::Arc {
                centre,
                sweep,
                radius,
                ..
            } if sweep > 0.0 && centre.distance(node.at) <= merge => Some((centre, radius)),
            _ => None,
        });
        if let Some((centre, radius)) = centre_of {
            let mut agrees = true;
            for &s in &node.sites {
                agrees &= (sites[s].distance(centre, tol)? - radius).abs() <= tol.confusion();
            }
            if agrees {
                node.at = centre;
                continue;
            }
        }
        node.at = solve_node(sites, &node.sites, node.at, spacing, tol)?;
    }
    // A branch point is only where its sites truly stand at one distance:
    // a centre beside a branch can name a third site a sampling's width
    // farther, and no point near it is equidistant from all three. The
    // sites farther than the nearest are dropped, and a point left with
    // fewer than three is no branch point.
    let mut solved: Vec<Node> = Vec::with_capacity(nodes.len());
    for mut node in nodes {
        let mut distances = Vec::with_capacity(node.sites.len());
        for &s in &node.sites {
            distances.push((s, sites[s].distance(node.at, tol)?));
        }
        let nearest = distances.iter().map(|d| d.1).fold(f64::INFINITY, f64::min);
        node.sites = distances
            .iter()
            .filter(|(_, d)| *d <= nearest + tol.confusion() * 10.0)
            .map(|(s, _)| *s)
            .collect();
        let whole_circle = node.sites.len() == 1
            && matches!(sites[node.sites[0]].kind, SiteKind::Arc { sweep, .. } if sweep >= core::f64::consts::TAU - 1e-9);
        if node.sites.len() >= 3 || whole_circle {
            solved.push(node);
        }
    }
    let mut nodes = solved;
    // A full circle holds its own centre, equidistant from all of it.
    for (s, site) in sites.iter().enumerate() {
        if let SiteKind::Arc {
            centre,
            radius,
            sweep,
            ..
        } = site.kind
            && sweep >= core::f64::consts::TAU - 1e-9
            && boundary.inside(centre)
            && !nodes.iter().any(|n| n.at.distance(centre) <= merge)
        {
            let nearer = sites.iter().enumerate().any(|(t, other)| {
                t != s
                    && other
                        .distance(centre, tol)
                        .is_ok_and(|d| d < radius - tolerance)
            });
            if !nearer {
                nodes.push(Node {
                    at: centre,
                    sites: vec![s],
                    boundary: false,
                });
            }
        }
    }
    // Boundary ends: convex corners and reflex corners.
    for &(at, s1, s2) in &boundary.convex {
        nodes.push(Node {
            at,
            sites: vec![s1, s2],
            boundary: true,
        });
    }
    for &(point_site, s1, s2) in &boundary.reflex {
        let SiteKind::Point { at } = sites[point_site].kind else {
            continue;
        };
        nodes.push(Node {
            at,
            sites: vec![point_site, s1, s2],
            boundary: true,
        });
    }

    let mut vertices: Vec<MedialVertex> = Vec::with_capacity(nodes.len());
    for node in &nodes {
        let clearance = if node.boundary {
            0.0
        } else {
            sites[node.sites[0]].distance(node.at, tol)?
        };
        vertices.push(MedialVertex {
            point: lift(&boundary.frame, node.at),
            clearance,
        });
    }

    // Branches: each pair's bisector between consecutive ends along it that
    // the sampling puts material between.
    let mut branches: Vec<MedialBranch> = Vec::new();
    let mut branch_sites: Vec<[usize; 2]> = Vec::new();
    let mut pairs: Vec<(usize, usize)> = pair_samples.keys().copied().collect();
    pairs.sort_unstable();
    for pair in pairs {
        let samples = &pair_samples[&pair];
        let ends: Vec<usize> = (0..nodes.len())
            .filter(|&i| nodes[i].sites.contains(&pair.0) && nodes[i].sites.contains(&pair.1))
            .collect();
        let exact = match (sites[pair.0].primitive(), sites[pair.1].primitive()) {
            (Some(p), Some(q)) => Bisector::between(p, q),
            _ => None,
        };
        let Some(bisector) = exact else {
            // No closed form: a fitted branch between exactly two ends.
            if ends.len() != 2 {
                if samples
                    .iter()
                    .all(|s| nodes.iter().any(|n| n.at.distance(*s) <= merge))
                {
                    continue;
                }
                ogeom_bail!(
                    NotDone,
                    "a fitted branch found {} ends where it needs two",
                    ends.len()
                );
            }
            let curve = fitted_branch(
                sites,
                pair,
                (nodes[ends[0]].at, nodes[ends[1]].at),
                samples,
                &boundary.frame,
                tolerance,
                tol,
            )?;
            let range = curve.domain();
            branches.push(MedialBranch {
                curve,
                range,
                ends: [ends[0], ends[1]],
                sites: [sites[pair.0].origin.clone(), sites[pair.1].origin.clone()],
            });
            branch_sites.push([pair.0, pair.1]);
            continue;
        };
        let unwrap = |t: f64, about: f64| -> f64 {
            if bisector.periodic() {
                let tau = core::f64::consts::TAU;
                about + (t - about + core::f64::consts::PI).rem_euclid(tau) - core::f64::consts::PI
            } else {
                t
            }
        };
        let about = bisector.param(samples[0]);
        let mut at_ends: Vec<(f64, usize)> = ends
            .iter()
            .map(|&i| (unwrap(bisector.param(nodes[i].at), about), i))
            .collect();
        at_ends.sort_by(|a, b| a.0.total_cmp(&b.0));
        let params: Vec<f64> = samples
            .iter()
            .map(|s| unwrap(bisector.param(*s), about))
            .collect();
        for w in at_ends.windows(2) {
            let ((t0, i0), (t1, i1)) = (w[0], w[1]);
            if (t1 - t0).abs() <= 1e-12 {
                continue;
            }
            // Samples strictly between the two ends, clear of both.
            let between = params.iter().zip(samples).any(|(t, s)| {
                *t > t0
                    && *t < t1
                    && s.distance(nodes[i0].at) > spacing
                    && s.distance(nodes[i1].at) > spacing
            });
            if !between {
                continue;
            }
            let curve = bisector.curve(&boundary.frame, (t0, t1), tol)?;
            let range = if matches!(bisector, Bisector::Line { .. }) {
                curve.domain()
            } else {
                (t0, t1)
            };
            branches.push(MedialBranch {
                curve,
                range,
                ends: [i0, i1],
                sites: [sites[pair.0].origin.clone(), sites[pair.1].origin.clone()],
            });
            branch_sites.push([pair.0, pair.1]);
        }
    }

    // The check: along every branch no third site nearer than its own two,
    // and every sampled centre on some branch.
    let mut deviation = 0.0_f64;
    for (b, branch) in branches.iter().enumerate() {
        let [s1, s2] = branch_sites[b];
        for k in 1..16 {
            let t = branch.range.0 + (branch.range.1 - branch.range.0) * f64::from(k) / 16.0;
            let p = flat(&boundary.frame, branch.curve.point_at(t, tol)?);
            let own = sites[s1].distance(p, tol)?;
            let other = sites[s2].distance(p, tol)?;
            deviation = deviation.max((own - other).abs());
            for (s, site) in sites.iter().enumerate() {
                if s == s1 || s == s2 {
                    continue;
                }
                let d = site.distance(p, tol)?;
                deviation = deviation.max(own - d);
            }
        }
    }
    if deviation > tolerance {
        ogeom_bail!(
            NotDone,
            "a branch strays {deviation} nearer another site than its own"
        );
    }
    for (pair, samples) in &pair_samples {
        for s in samples {
            if nodes.iter().any(|n| n.at.distance(*s) <= merge) {
                continue;
            }
            let on_some = branches.iter().enumerate().any(|(b, _)| {
                let [a, c] = branch_sites[b];
                (a, c) == *pair
            });
            if !on_some {
                ogeom_bail!(NotDone, "sampled axis points lie on no branch");
            }
        }
    }
    Ok(MedialGraph {
        vertices,
        branches,
        deviation,
        frame: boundary.frame,
        sites: sites.clone(),
        branch_sites,
    })
}

fn circumcentre(a: Point2, b: Point2, c: Point2) -> Option<Point2> {
    let d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if d.abs() <= 1e-300 {
        return None;
    }
    let (a2, b2, c2) = (
        a.x * a.x + a.y * a.y,
        b.x * b.x + b.y * b.y,
        c.x * c.x + c.y * c.y,
    );
    Some(Point2::new(
        (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d,
        (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d,
    ))
}

/// The point equidistant from a branch point's sites, by Newton from the
/// sampled estimate: the first site's distance against each other's.
fn solve_node(
    sites: &[Site],
    labels: &[usize],
    guess: Point2,
    spacing: f64,
    tol: Tolerances,
) -> OgeomResult<Point2> {
    if labels.len() < 3 {
        return Ok(guess);
    }
    // A handover: two of the sites share a foot (a reflex corner and its
    // own edge, a line and the arc it runs on into), so their distances
    // agree all along the normal there, and the axis passes from one to
    // the other where that normal meets a third site's bisector.
    let mut feet = Vec::with_capacity(labels.len());
    for &s in labels {
        feet.push(sites[s].foot(guess, tol)?);
    }
    for i in 0..labels.len() {
        for j in i + 1..labels.len() {
            if feet[i].distance(feet[j]) > spacing * 2.0 {
                continue;
            }
            let Some(k) = (0..labels.len())
                .find(|&k| k != i && k != j && feet[k].distance(feet[i]) > spacing * 2.0)
            else {
                continue;
            };
            let joint = shared_joint(&sites[labels[i]], &sites[labels[j]]);
            let normal = shared_normal(&sites[labels[i]], &sites[labels[j]], joint, guess, tol);
            let third = &sites[labels[k]];
            let mut t = (guess - joint).dot(normal);
            for _ in 0..60 {
                let p = joint + normal * t;
                let f = third.distance(p, tol)? - t;
                let h = (t.abs() + 1.0) * 1e-7;
                let g = (third.distance(joint + normal * (t + h), tol)? - (t + h) - f) / h;
                if g.abs() <= 1e-300 {
                    break;
                }
                let step = f / g;
                t -= step;
                if step.abs() <= tol.confusion() * 1e-4 {
                    break;
                }
            }
            return Ok(joint + normal * t);
        }
    }
    let gradient = |s: usize, p: Point2| -> OgeomResult<(f64, Vector2)> {
        let foot = sites[s].foot(p, tol)?;
        let d = p - foot;
        let m = d.magnitude();
        Ok((
            m,
            if m > 0.0 {
                d * (1.0 / m)
            } else {
                Vector2::new(0.0, 0.0)
            },
        ))
    };
    let mut p = guess;
    for _ in 0..60 {
        let (d0, g0) = gradient(labels[0], p)?;
        // Least squares over every other site: a branch point of four
        // sites (a square's centre) is still one point.
        let (mut jtj, mut jtf) = ([[0.0_f64; 2]; 2], [0.0_f64; 2]);
        let mut worst = 0.0_f64;
        for &s in &labels[1..] {
            let (d, g) = gradient(s, p)?;
            let f = d0 - d;
            let row = g0 - g;
            worst = worst.max(f.abs());
            jtj[0][0] += row.x * row.x;
            jtj[0][1] += row.x * row.y;
            jtj[1][0] += row.y * row.x;
            jtj[1][1] += row.y * row.y;
            jtf[0] += row.x * f;
            jtf[1] += row.y * f;
        }
        if worst <= tol.confusion() * 1e-3 {
            return Ok(p);
        }
        let det = jtj[0][0] * jtj[1][1] - jtj[0][1] * jtj[1][0];
        if det.abs() <= 1e-300 {
            break;
        }
        let dx = (jtj[1][1] * jtf[0] - jtj[0][1] * jtf[1]) / det;
        let dy = (jtj[0][0] * jtf[1] - jtj[1][0] * jtf[0]) / det;
        p = Point2::new(p.x - dx, p.y - dy);
    }
    Ok(p)
}

/// Where two sites meet: a corner's own point, or the ends of the two
/// that coincide.
fn shared_joint(a: &Site, b: &Site) -> Point2 {
    let ends = |site: &Site| -> Vec<Point2> {
        match &site.kind {
            SiteKind::Segment { a, b } => vec![*a, *b],
            SiteKind::Arc {
                centre,
                radius,
                start,
                sweep,
            } => {
                let at = |angle: f64| *centre + Vector2::new(angle.cos(), angle.sin()) * *radius;
                vec![at(*start), at(start + sweep)]
            }
            SiteKind::Point { at } => vec![*at],
            SiteKind::Curve { polyline, .. } => {
                vec![polyline[0].1, polyline[polyline.len() - 1].1]
            }
        }
    };
    let (ea, eb) = (ends(a), ends(b));
    let mut best = (f64::INFINITY, ea[0]);
    for p in &ea {
        for q in &eb {
            if p.distance(*q) < best.0 {
                best = (p.distance(*q), if eb.len() == 1 { *q } else { *p });
            }
        }
    }
    best.1
}

/// The inward normal two sites share at their common foot: a segment's
/// own, an arc's radius, and otherwise the direction from the foot to the
/// estimate.
fn shared_normal(a: &Site, b: &Site, joint: Point2, guess: Point2, tol: Tolerances) -> Vector2 {
    let own = |site: &Site| -> Option<Vector2> {
        match site.kind {
            SiteKind::Segment { a, b } => {
                let d = b - a;
                Vector2::new(-d.y, d.x).normalized(tol).ok()
            }
            SiteKind::Arc { centre, sweep, .. } => {
                let radial = (centre - joint).normalized(tol).ok()?;
                Some(if sweep > 0.0 { radial } else { -radial })
            }
            _ => None,
        }
    };
    own(a)
        .or_else(|| own(b))
        .or_else(|| (guess - joint).normalized(tol).ok())
        .unwrap_or(Vector2::new(0.0, 1.0))
}

/// A branch with no closed form: points solved onto the bisector exactly,
/// ordered from one end to the other, and fitted.
fn fitted_branch(
    sites: &[Site],
    pair: (usize, usize),
    ends: (Point2, Point2),
    samples: &[Point2],
    frame: &Frame,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::Curve> {
    let f = |p: Point2| -> OgeomResult<(f64, Vector2)> {
        let (fa, fb) = (sites[pair.0].foot(p, tol)?, sites[pair.1].foot(p, tol)?);
        let (da, db) = (p - fa, p - fb);
        let (ma, mb) = (da.magnitude(), db.magnitude());
        let g = da * (1.0 / ma.max(1e-300)) - db * (1.0 / mb.max(1e-300));
        Ok((ma - mb, g))
    };
    let chord = ends.1 - ends.0;
    let mut placed: Vec<(f64, Point2)> = Vec::with_capacity(samples.len());
    for s in samples {
        let mut p = *s;
        for _ in 0..30 {
            let (value, g) = f(p)?;
            let gg = g.dot(g);
            if value.abs() <= tol.confusion() * 1e-3 || gg <= 1e-300 {
                break;
            }
            p -= g * (value / gg);
        }
        placed.push(((p - ends.0).dot(chord), p));
    }
    placed.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut points: Vec<Point> = vec![lift(frame, ends.0)];
    for (_, p) in &placed {
        if p.distance(ends.0) > tol.confusion() * 10.0
            && p.distance(ends.1) > tol.confusion() * 10.0
        {
            points.push(lift(frame, *p));
        }
    }
    points.push(lift(frame, ends.1));
    let fitted = ogeom_geom::fit::fit_points(&points, 3, tolerance * 0.5, tol)?;
    Ok(ogeom_geom::Curve::BSpline(fitted.curve))
}
