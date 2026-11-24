//! Boolean operations.
//!
//! Fuse, common, cut, section and split are one algorithm (general fuse)
//! plus a selection predicate: every face of each solid is split by the
//! other solid's faces, the pieces are classified against the other solid,
//! and a filter keeps the pieces the operation wants.
//!
//! # One pipeline, any analytic face
//!
//! The pipeline runs on solids whose faces are any surface the intersectors
//! answer for (planes, cylinders, spheres, cones, tori), with the planar
//! case as nothing more than the case where every curve is a line. Sections
//! come from [`intersect_surfaces`](ogeom_intersect::intersect_surfaces), exact
//! with same-parameter pcurves where the projection has a closed form and
//! *marched and fitted to a stated tolerance* where it does not. Paves (the
//! points where sections meet boundary edges and each other) come from the
//! exact curve/curve intersection. Every face is then split in its own
//! parameter space by an arrangement of *strands*: polyline scaffolding, each
//! naming the exact sub-curve it stands for, so the combinatorics run on
//! polylines and the rebuilt result is exact curves, pcurves attached both
//! sides, sewn back into shared topology.
//!
//! Pieces are classified against the other solid by the exact ray classifier
//! and the filters select: fuse keeps what is outside, common what is inside,
//! cut flips the tool's contribution. Same-domain contact (a piece lying
//! *on* the other boundary) is resolved by whether the two solids' outward
//! normals agree there, where the contact's edges can be imaged in the
//! shared chart, and refused by name where they cannot. A tangential touch
//! splits nothing and is carried as the contact curve it is.
//!
//! An edge of one operand lying in a face of the other is an on-on contact:
//! it splits that face as a section would, and the pieces either side are
//! classified on their own. A face's outline for those tests holds every
//! side of its chart, a side that closes to one point (a pole) included.
//!
//! # Sheets
//!
//! A sheet (a face, a shell or a compound of those) runs the same pipeline
//! with one side left out: its faces are split by the other argument's, and
//! the other argument's faces only cut. [`common`] and [`cut`] classify the
//! sheet's pieces against a solid or a half space and keep those inside or
//! outside it. [`split_sheet`] keeps every piece, unclassified. The kept
//! pieces are sewn into shells with no closure asked of them.

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

mod along_edges;
mod arrange;
mod aside;
mod ball_chart;
mod bins;
mod box_tree;
mod defeature;
mod fill;
mod fuse_all;
mod gather;
mod half_space;
mod mesh_region;
mod seam_join;
mod seated;
mod section_face;
mod sections;
mod self_union;
mod sub_shapes;
#[cfg(feature = "test-hooks")]
#[doc(hidden)]
pub mod work;

pub use along_edges::note_meeting_along_edge;

pub use defeature::remove_faces;
pub use fuse_all::fuse_all;
pub use mesh_region::cut_by_mesh_region;
pub use ogeom_algo::POLE_TOUCHING;
// The pole re-chart moved to `ogeom_algo`; re-exported here for one release.
pub use ogeom_algo::{POLE_MARGIN, Recharted, rechart_clear_of_poles};
pub use section_face::section_face;
pub use self_union::self_union;

use ogeom_algo::{
    Built, Containment, History, is_shell_closed, make_edge_between, make_face_on, make_vertex,
    make_wire, sew_around, sew_around_within, shape_bounds,
};
use ogeom_core::diag::{self, Channel};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve2d as _;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry, Transformable};
use ogeom_math::{Point, Point2, Vector, Weighted};
use ogeom_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, explore, explore_unique,
};

use arrange::{Strand, Traversal, assemble as arrange_pieces, inside_many, inside_rings};

use fill::*;
use gather::*;
use sections::*;

/// Parameter-space chord for the polyline scaffolding.
const SCAFFOLD_CHORD: f64 = 1e-3;

/// How far a curved strand's chords may stand off its curve on a plane, and
/// a face outline's chords off its pcurves on any surface. A piece is
/// classified by a point inside its polyline outline, and a piece thinner
/// than its outline's chords bow (the cap annulus between two drums a
/// thousandth apart in radius) has that point outside the true piece. A
/// section's piece is kept on a face by its middle read against the face's
/// outline, and a piece standing nearer the rim than the rim's chords bow
/// (an arc between two crossings of a bowed seam) reads on the wrong side.
const STRAND_SAG: f64 = 1e-4;

/// Samples enough that a curve's chords stand within [`STRAND_SAG`] of it in
/// the chart: from `count`, raised as the measured sag at the chords'
/// middles asks, up to 512. With the count, the sag measured at it, or
/// infinity where the count was raised past the last measure.
fn sag_count(
    at: &dyn Fn(f64) -> OgeomResult<Point2>,
    lo: f64,
    hi: f64,
    mut count: usize,
) -> OgeomResult<(usize, f64)> {
    const MOST: usize = 512;
    let mut measured = f64::INFINITY;
    for _ in 0..3 {
        if count >= MOST {
            break;
        }
        let mut sag = 0.0_f64;
        #[allow(clippy::cast_precision_loss)]
        let step = (hi - lo) / count as f64;
        let mut previous = at(lo)?;
        for i in 1..=count {
            #[allow(clippy::cast_precision_loss)]
            let t = (i as f64).mul_add(step, lo);
            let here = at(t)?;
            let middle = at(t - step / 2.0)?;
            sag = sag.max(middle.distance(Point2::midpoint(previous, here)));
            previous = here;
        }
        if sag <= STRAND_SAG {
            measured = sag;
            break;
        }
        // Sag falls as the square of the step.
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let raised = ((count as f64) * (sag / STRAND_SAG).sqrt()).ceil() as usize;
        count = raised.clamp(count + 1, MOST);
    }
    Ok((count, measured))
}

/// [`sag_count`] for a circle, whose chords' sag is known: `r(1 - cos(θ/2))`
/// over a step of `θ`. `None` for any other curve.
fn circle_sag_count(pcurve: &PlanarCurve, lo: f64, hi: f64, count: usize) -> Option<usize> {
    const MOST: usize = 512;
    let PlanarCurve::Circle(circle) = pcurve else {
        return None;
    };
    let radius = circle.circle().radius();
    if radius <= STRAND_SAG {
        return Some(count);
    }
    let step = 2.0 * (1.0 - STRAND_SAG / radius).acos();
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let needed = ((hi - lo).abs() / step).ceil().min(MOST as f64) as usize;
    Some(count.max(needed))
}

/// How far the chords of a circle sampled in `count` equal steps over
/// `lo..hi` stand off it: `r(1 - cos(θ/2))` over a step of `θ`. `None` for
/// any other curve.
fn circle_sag(pcurve: &PlanarCurve, lo: f64, hi: f64, count: usize) -> Option<f64> {
    let PlanarCurve::Circle(circle) = pcurve else {
        return None;
    };
    #[allow(clippy::cast_precision_loss)]
    let step = (hi - lo).abs() / count.max(1) as f64;
    Some(circle.circle().radius() * (1.0 - (step / 2.0).cos()))
}

/// Parameter-space node snap for the arrangement.
const PARAM_SNAP: f64 = 1e-6;

/// The widest chart gap a seam column's end may stand from the face's
/// other edges and still meet them within their radii in space.
const SEAM_HAIR: f64 = 1e-3;

/// How many chords a curved strand of a planar piece is resampled into.
const FINE: usize = 1024;

/// How many of a resampled arc's chords are passed over together where a
/// ray cannot reach them.
const FINE_BLOCK: usize = 32;

/// A planar piece's rings resampled finely off their exact curves. The
/// piece's own polylines bow inside a curved strand by more than the width
/// of a sliver beside it, so a probe picked inside them can stand outside
/// the piece. These rings hold it to a few microns.
///
/// A circular strand's samples are drawn as a ray asks for them: a ray
/// meets a circle in at most two places, and the samples far from its
/// height cannot be crossed by it.
struct FineRings<'a> {
    rings: Vec<Vec<FinePart<'a>>>,
    tol: Tolerances,
}

/// One strand of a [`FineRings`] ring.
enum FinePart<'a> {
    /// Its samples, in the ring's direction.
    Points(Vec<Point2>),
    /// A circle's samples, drawn when first asked for.
    Arc(FineArc<'a>),
}

/// A circular strand's [`FINE`] chords, each sample drawn once, when a ray
/// first needs it.
struct FineArc<'a> {
    curve: &'a PlanarCurve,
    circle: ogeom_math::Circle2,
    from: f64,
    to: f64,
    wrap: Option<(f64, f64)>,
    reversed: bool,
    /// The samples drawn so far; empty until the first is.
    drawn: core::cell::RefCell<Vec<Option<Point2>>>,
}

impl FineArc<'_> {
    /// The parameter sample `i` is drawn at.
    fn parameter(&self, i: usize) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let t = self.from + (self.to - self.from) * (i as f64 / FINE as f64);
        self.wrap.map_or(t, |domain| fold(t, domain))
    }

    /// Sample `i`. A circle is drawn at any parameter.
    fn point(&self, i: usize, tol: Tolerances) -> Point2 {
        if let Some(Some(p)) = self.drawn.borrow().get(i) {
            return *p;
        }
        let p = self
            .curve
            .point_at(self.parameter(i), tol)
            .unwrap_or(Point2::new(f64::NAN, f64::NAN));
        let mut drawn = self.drawn.borrow_mut();
        if drawn.is_empty() {
            drawn.resize(FINE + 1, None);
        }
        drawn[i] = Some(p);
        p
    }

    /// How far the circle's coordinates stand from its centre's, over and
    /// above its rounding.
    fn slack(&self) -> f64 {
        let centre = self.circle.centre();
        1e-9 * (1.0 + self.circle.radius() + centre.x.abs() + centre.y.abs())
    }

    /// Whether a horizontal ray from `p` may cross a chord between samples
    /// `lo` and `hi`. Every sample between them lies within the arc's length
    /// between their parameters of sample `lo`, so where that reach stays
    /// wholly above or below `p`'s height no chord there straddles it. A
    /// parameter folded across its domain's end inside the block, or a
    /// sample that is not finite, leaves the block to be asked chord by
    /// chord.
    fn may_cross(&self, lo: usize, hi: usize, p: Point2, tol: Tolerances) -> bool {
        let (t0, t1) = (self.parameter(lo), self.parameter(hi));
        #[allow(clippy::cast_precision_loss)]
        let step = (self.to - self.from) * ((hi - lo) as f64 / FINE as f64);
        let unfolded = self
            .wrap
            .is_none_or(|(a, b)| ((t1 - t0) - step).abs() <= 0.5 * (b - a).abs());
        let start = self.point(lo, tol);
        let reach = self.circle.radius() * (t1 - t0).abs() + self.slack();
        if !(unfolded && reach.is_finite() && start.y.is_finite() && p.y.is_finite()) {
            return true;
        }
        start.y - reach <= p.y && start.y + reach >= p.y
    }

    /// The greatest x and the least and greatest y any sample may have,
    /// none drawn: a sample is the centre moved along each axis by the
    /// radius times a cosine or a sine.
    fn extent(&self) -> (f64, f64, f64) {
        let frame = self.circle.frame();
        let (x, y) = (frame.x().vector(), frame.y().vector());
        let r = self.circle.radius();
        let centre = self.circle.centre();
        let across = r * (x.x.abs() + y.x.abs()) + self.slack();
        let up = r * (x.y.abs() + y.y.abs()) + self.slack();
        (centre.x + across, centre.y - up, centre.y + up)
    }
}

impl<'a> FineRings<'a> {
    /// `None` for a face that is not a plane, or a ring with a strand this
    /// cannot resample (a contact or a pole). Where `only` is given, a ring
    /// it marks false is left empty.
    fn of(
        face: &'a GFace,
        rings: &[Vec<Traversal<Tag>>],
        sections: &'a [SectionRec],
        from_a: bool,
        only: Option<&[bool]>,
        tol: Tolerances,
    ) -> Option<Self> {
        if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
            return None;
        }
        let mut lines: Vec<Vec<FinePart<'a>>> = Vec::with_capacity(rings.len());
        for (k, ring) in rings.iter().enumerate() {
            let mut line: Vec<FinePart<'a>> = Vec::new();
            if only.is_some_and(|wanted| !wanted.get(k).copied().unwrap_or(true)) {
                lines.push(line);
                continue;
            }
            for traversal in ring {
                let (curve, a, b, wrap) = fine_strand(face, sections, from_a, &traversal.tag)?;
                if let PlanarCurve::Circle(circle) = curve {
                    line.push(FinePart::Arc(FineArc {
                        curve,
                        circle: circle.circle(),
                        from: a,
                        to: b,
                        wrap,
                        reversed: traversal.reversed,
                        drawn: core::cell::RefCell::new(Vec::new()),
                    }));
                    continue;
                }
                let count = if matches!(curve, PlanarCurve::Line(_)) {
                    1
                } else {
                    FINE
                };
                let mut points: Vec<Point2> = Vec::with_capacity(count + 1);
                for i in 0..=count {
                    #[allow(clippy::cast_precision_loss)]
                    let t = a + (b - a) * (i as f64 / count as f64);
                    let t = wrap.map_or(t, |domain| fold(t, domain));
                    points.push(curve.point_at(t, tol).ok()?);
                }
                if traversal.reversed {
                    points.reverse();
                }
                line.push(FinePart::Points(points));
            }
            lines.push(line);
        }
        Some(Self { rings: lines, tol })
    }

    /// A curved face's rings resampled finely in its chart, as [`Self::of`]
    /// resamples a plane's, every strand drawn in full: for a piece lying on
    /// a partner's surface, where a probe standing outside the piece reads
    /// the other solid's side of the partner, not the piece's. `None` where
    /// a strand cannot be resampled or a ring does not close in the chart
    /// (one winding the chart's period, which the count does not hold).
    fn of_curved(
        face: &'a GFace,
        rings: &[Vec<Traversal<Tag>>],
        sections: &'a [SectionRec],
        from_a: bool,
        tol: Tolerances,
    ) -> Option<Self> {
        use ogeom_geom::Surface as _;
        // A periodic chart's strands may run a whole period (a rim's
        // circle), which no chord between its ends draws.
        if face.surface.is_periodic_u() || face.surface.is_periodic_v() {
            return None;
        }
        let mut lines: Vec<Vec<FinePart<'a>>> = Vec::with_capacity(rings.len());
        for ring in rings {
            let mut line: Vec<FinePart<'a>> = Vec::new();
            let mut drawn: Vec<Point2> = Vec::new();
            for traversal in ring {
                let (curve, a, b, wrap) = fine_strand(face, sections, from_a, &traversal.tag)?;
                let count = if matches!(curve, PlanarCurve::Line(_)) {
                    1
                } else {
                    FINE
                };
                let mut points: Vec<Point2> = Vec::with_capacity(count + 1);
                for i in 0..=count {
                    #[allow(clippy::cast_precision_loss)]
                    let t = a + (b - a) * (i as f64 / count as f64);
                    let t = wrap.map_or(t, |domain| fold(t, domain));
                    points.push(curve.point_at(t, tol).ok()?);
                }
                if traversal.reversed {
                    points.reverse();
                }
                drawn.extend_from_slice(&points);
                line.push(FinePart::Points(points));
            }
            let (Some(first), Some(last)) = (drawn.first(), drawn.last()) else {
                return None;
            };
            let span = drawn.iter().fold(0.0_f64, |m, p| m.max(p.distance(*first)));
            if last.distance(*first) > span * 1e-6 + tol.parametric() {
                return None;
            }
            lines.push(line);
        }
        Some(Self { rings: lines, tol })
    }

    /// Even-odd containment of `p` in the rings, counted as
    /// [`inside_rings`] counts it over the rings drawn in full.
    fn inside(&self, p: Point2) -> bool {
        self.rings
            .iter()
            .fold(false, |acc, ring| acc != self.ring_inside(ring, p))
    }

    /// Whether a ray from `p` crosses `ring` an odd number of times: each
    /// strand's chords, the joins between strands, and the closing join,
    /// whether a ray crosses a chord not depending on the way it runs.
    fn ring_inside(&self, ring: &[FinePart<'_>], p: Point2) -> bool {
        let tol = self.tol;
        let mut inside = false;
        let mut flip = |a: Point2, b: Point2| {
            if arrange::ray_crosses(a, b, p) {
                inside = !inside;
            }
        };
        let mut first: Option<Point2> = None;
        let mut last: Option<Point2> = None;
        for part in ring {
            let (head, tail) = match part {
                FinePart::Points(points) => {
                    for w in points.windows(2) {
                        flip(w[0], w[1]);
                    }
                    (points[0], points[points.len() - 1])
                }
                FinePart::Arc(arc) => {
                    let mut lo = 0;
                    while lo < FINE {
                        let hi = (lo + FINE_BLOCK).min(FINE);
                        if arc.may_cross(lo, hi, p, tol) {
                            for i in lo..hi {
                                flip(arc.point(i, tol), arc.point(i + 1, tol));
                            }
                        }
                        lo = hi;
                    }
                    let ends = (arc.point(0, tol), arc.point(FINE, tol));
                    if arc.reversed { (ends.1, ends.0) } else { ends }
                }
            };
            if let Some(at) = last {
                flip(at, head);
            }
            first.get_or_insert(head);
            last = Some(tail);
        }
        if let (Some(first), Some(last)) = (first, last) {
            flip(last, first);
        }
        inside
    }

    /// Whether a ray from `p` cannot reach `ring`: every sample of it lies
    /// below or above `p`, or to its left, so no chord of it straddles the
    /// ray where the ray runs.
    fn out_of_reach(ring: &[FinePart<'_>], p: Point2) -> bool {
        let (mut right, mut low, mut high) = (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        let mut reach = |x: f64, y: f64| {
            right = right.max(x);
            low = low.min(y);
            high = high.max(y);
            x.is_finite() && y.is_finite()
        };
        for part in ring {
            let finite = match part {
                FinePart::Points(points) => points.iter().all(|q| reach(q.x, q.y)),
                FinePart::Arc(arc) => {
                    let (x, lo, hi) = arc.extent();
                    arc.from.is_finite() && arc.to.is_finite() && reach(x, lo) && reach(x, hi)
                }
            };
            if !finite {
                return false;
            }
        }
        ring.is_empty() || p.y < low || p.y > high || p.x > right
    }

    /// [`FineRings::inside`] for `p` over the rings `reached` marks and any
    /// others it could be asked over, where those others cannot change the
    /// answer: `None` where a ring `reached` leaves out is crossed by the
    /// ray an odd number of times, so the answer depends on whether it is
    /// counted.
    fn inside_beside(&self, p: Point2, reached: &[bool]) -> Option<bool> {
        let mut inside = false;
        for (k, ring) in self.rings.iter().enumerate() {
            let counted = reached.get(k).copied().unwrap_or(true);
            if !counted && Self::out_of_reach(ring, p) {
                continue;
            }
            let crossed = self.ring_inside(ring, p);
            if !counted && crossed {
                return None;
            }
            inside ^= crossed;
        }
        Some(inside)
    }

    /// Every ring's samples, drawn in full, in the ring's direction.
    fn points(&self) -> Vec<Vec<Point2>> {
        self.rings
            .iter()
            .map(|ring| {
                let mut line = Vec::new();
                for part in ring {
                    match part {
                        FinePart::Points(points) => line.extend_from_slice(points),
                        FinePart::Arc(arc) => {
                            let start = line.len();
                            line.extend((0..=FINE).map(|i| arc.point(i, self.tol)));
                            if arc.reversed {
                                line[start..].reverse();
                            }
                        }
                    }
                }
                line
            })
            .collect()
    }
}

/// A pcurve, the range of it a strand covers, and the domain the range
/// folds into where it does.
type FineStrand<'a> = (&'a PlanarCurve, f64, f64, Option<(f64, f64)>);

/// The pcurve a strand of a planar piece is resampled along, its range,
/// and the domain a closed section's parameters fold into. `None` for a
/// contact or a pole, or a tag naming nothing.
fn fine_strand<'a>(
    face: &'a GFace,
    sections: &'a [SectionRec],
    from_a: bool,
    tag: &Tag,
) -> Option<FineStrand<'a>> {
    Some(match *tag {
        Tag::Boundary { edge, range } => {
            let e = face.edges.get(edge)?;
            (
                &e.pcurve,
                rescale(range.0, e.crange, e.prange),
                rescale(range.1, e.crange, e.prange),
                None,
            )
        }
        Tag::Section { section, range } => {
            let record = sections.get(section)?;
            let pcurve = if from_a { &record.pc_a } else { &record.pc_b };
            let wrap = record.closed.then(|| record.curve.domain());
            (pcurve, range.0, range.1, wrap)
        }
        Tag::Contact { .. } | Tag::Pole { .. } => return None,
    })
}

/// Which of a piece's rings a horizontal ray from any of `probes` may
/// cross: those whose outline, widened past what its chords bow, spans the
/// probe's height and reaches to its right. A ray crosses no other, coarse
/// or fine, and an empty ring in its place counts the same.
fn rings_near(outlines: &[Vec<Point2>], probes: &[Point2]) -> Vec<bool> {
    outlines
        .iter()
        .map(|line| {
            let (mut low, mut high) = (
                Point2::new(f64::INFINITY, f64::INFINITY),
                Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            );
            for p in line {
                low = Point2::new(low.x.min(p.x), low.y.min(p.y));
                high = Point2::new(high.x.max(p.x), high.y.max(p.y));
            }
            if low
                .x
                .partial_cmp(&high.x)
                .is_none_or(core::cmp::Ordering::is_gt)
            {
                return true;
            }
            let margin = 1e-3 + 1e-2 * low.distance(high);
            probes
                .iter()
                .any(|p| p.y >= low.y - margin && p.y <= high.y + margin && p.x <= high.x + margin)
        })
        .collect()
}

/// Whether a chart point lies inside a planar face, counted exactly: the
/// crossings of a ray from it with the boundary's own pcurves, not their
/// polylines. A trim's chords bow inside a curved rim by more than the
/// width of a sliver beside it (the annulus a drum a few hundredths of a
/// millimetre wider leaves on a cap), and a point in that sliver reads
/// outside the polyline. `None` for a face that is not a plane or has no
/// curved edge, or where every ray tried meets the boundary at an end, a
/// tangency or the start.
fn inside_plane_exactly(face: &GFace, p: Point2, tol: Tolerances) -> Option<bool> {
    // A boundary of straight edges is its own polyline: nothing to correct.
    if !matches!(face.surface, SurfaceGeometry::Plane(_))
        || face
            .edges
            .iter()
            .all(|e| matches!(e.pcurve, PlanarCurve::Line(_)))
    {
        return None;
    }
    // On the boundary to well within the weld distance (three quarters of
    // it, or an edge's own doubt where wider) is on it: a sliver that thin is
    // welded, no parity speaks for such a point, and the outline's reading
    // stands. One at the weld distance itself is read exactly, as two drums
    // that far apart are sectioned.
    let at = face.surface.point_at(p.x, p.y, tol).ok()?;
    for e in &face.edges {
        let doubt = (e.standoff * 2.0).max(tol.confusion() * 75.0);
        if !e.bound.expanded(doubt).contains(at) {
            continue;
        }
        if distance_to_edge_curve(&e.curve, e.crange, at, tol).ok()? <= doubt {
            return None;
        }
    }
    let reach = face
        .edges
        .iter()
        .map(|e| e.bound.diagonal())
        .sum::<f64>()
        .max(1.0)
        * 4.0;
    let mut curves: Vec<(PlanarCurve, (f64, f64))> = Vec::new();
    for e in &face.edges {
        curves.push((e.pcurve.clone(), e.prange));
        if let Some((other, range)) = &e.other_side {
            curves.push((other.clone(), *range));
        }
    }
    for pole in &face.poles {
        curves.push((pole.pcurve.clone(), pole.prange));
    }
    let options = ogeom_intersect::CurveCurveOptions::default();
    'directions: for (dx, dy) in [(1.0, 0.0137), (-0.31, 1.0), (-1.0, -0.47), (0.53, -1.0)] {
        let length = f64::hypot(dx, dy);
        let far = Point2::new(p.x + dx / length * reach, p.y + dy / length * reach);
        let Ok(ray) = ogeom_geom::Line2d::segment(p, far, tol) else {
            return None;
        };
        let ray: PlanarCurve = ray.into();
        let mut crossings = 0_usize;
        for (curve, (a, b)) in &curves {
            let (lo, hi) = (a.min(*b), a.max(*b));
            if hi - lo <= tol.parametric() {
                continue;
            }
            let Ok(piece) = ogeom_geom::Trimmed2d::new(curve.clone(), lo, hi, tol) else {
                continue 'directions;
            };
            let piece: PlanarCurve = piece.into();
            let Ok(met) = ogeom_intersect::intersect_curves_2d(&ray, &piece, options, tol) else {
                continue 'directions;
            };
            if !met.overlaps.is_empty() {
                continue 'directions;
            }
            let span = hi - lo;
            for hit in &met.crossings {
                let at_start = hit.on_a <= tol.confusion();
                let at_end =
                    (hit.on_b - lo).abs() <= span * 1e-9 || (hi - hit.on_b).abs() <= span * 1e-9;
                if at_start || at_end || hit.reach > 0.0 {
                    continue 'directions;
                }
                crossings += 1;
            }
        }
        return Some(crossings % 2 == 1);
    }
    None
}

// --- the general fuse --------------------------------------------------------

/// One piece of one argument face, classified against the other argument.
struct FacePiece {
    /// Which argument's face list, and which face.
    from_a: bool,
    face: usize,
    rings: Vec<Vec<Traversal<Tag>>>,
    /// The rings as chart polylines, for the coincidence pairing below.
    outlines: Vec<Vec<Point2>>,
    /// The interior point the state was decided at, in space.
    probe: Point,
    state: PieceState,
    /// Whether the *other* argument already contributes this same patch of
    /// this same surface.
    ///
    /// Two faces lying on one surface with their material on the same side
    /// bound the union and the intersection once between them, not twice,
    /// but "once" is a statement about the *patch*, not about which argument
    /// it came from. Where a tool's cap fills a hole its own bore left, the
    /// part has no piece there at all, and dropping the tool's by argument
    /// identity leaves the result with a hole in it. So the duplicate is
    /// identified by containment, and only a piece genuinely stood in for is
    /// dropped.
    covered: bool,
    /// Whether the piece is its whole face, untouched: nothing crosses it
    /// and none of its edges is split.
    whole: bool,
}

struct GeneralFused {
    a: GSolid,
    b: GSolid,
    sections: Vec<SectionRec>,
    contacts: Vec<ContactRec>,
    /// Curves the two boundaries touch along without crossing. They take no
    /// part in the classification below (that is what tangency means) and
    /// are carried for the consumers that want the contact itself.
    tangents: Vec<TangentRec>,
    pieces: Vec<FacePiece>,
    /// Pieces of the two arguments standing for one patch of one surface,
    /// as indices into `pieces`: the first argument's piece, then the
    /// second's.
    coincident: Vec<(usize, usize)>,
    /// Junctions several paves describe, each resolved once for every
    /// strand that ends in it.
    junctions: Vec<Junction>,
}

/// One junction the paving found several times over.
///
/// A tolerant rail meeting a wedge's faces near one corner collects a
/// cluster of crossings inside its own stated radius: each face's section
/// stops at its own crossing with the rail, and the crossings sit a few
/// tenths of a micron to a few hundred apart along it. The rail's strands
/// split once per cluster, at its first pave. The sections' ends still name
/// their own crossings. And the two descriptions of the junction can sit
/// apart by the cluster's whole span, which no single strand's honesty
/// covers. So the junction is one vertex standing at the first pave, owning
/// the span the paves disagree by, and every strand end inside that span
/// names it, whichever pave its own trim stopped at.
#[derive(Debug, Clone, Copy)]
struct Junction {
    /// Where the cluster's first pave sits.
    at: Point,
    /// How far a strand end may sit from `at` and still be this junction:
    /// the cluster's span plus the rail's own reach.
    reach: f64,
    /// Whether `at` is an edge's own end vertex that a section tangent to
    /// the edge was welded onto: the section's chart image is bent onto it
    /// rather than left where the section stopped.
    onto_vertex: bool,
}

/// Paves on each seam of `solid`'s faces where a junction the seam passes
/// through holds section ends of the face on both sides of it.
///
/// The ends of one junction are one point, and a strand end is carried onto
/// a node the junction holds, or onto another end, across the face's chart.
/// Ends either side of a seam stand near its two columns, a period apart in
/// the chart, and no path inside the chart joins them: with no node on the
/// seam they dangle on either side of it. Split where it passes nearest the
/// junction, the seam puts a node on each column, and each end is carried
/// onto the one on its side, as ends a junction holds on one side of the
/// seam are carried onto one another.
fn seams_through_junctions(
    solid: &GSolid,
    from_a: bool,
    (sections, pieces): (&[SectionRec], &[SectionPiece]),
    junctions: &[Junction],
    paves: &mut ogeom_core::FastMap<EdgeKey, Vec<Pave>>,
    tol: Tolerances,
) -> OgeomResult<()> {
    for (fi, face) in solid.faces.iter().enumerate() {
        let seams: Vec<&BoundaryEdge> = face
            .edges
            .iter()
            .filter(|e| e.other_side.is_some())
            .collect();
        if seams.is_empty() {
            continue;
        }
        // Each end of the face's section pieces, in space and in the chart.
        let mut ends: Vec<(Point, Point2)> = Vec::new();
        for sp in pieces {
            let section = &sections[sp.section];
            let (face_of, pcurve, hugs) = if from_a {
                (section.face_a, &section.pc_a, sp.hugs[0])
            } else {
                (section.face_b, &section.pc_b, sp.hugs[1])
            };
            if face_of != fi || hugs || section.closed {
                continue;
            }
            for t in [sp.range.0, sp.range.1] {
                ends.push((section.curve.point_at(t, tol)?, pcurve.point_at(t, tol)?));
            }
        }
        for e in seams {
            let Some((other, orange)) = &e.other_side else {
                continue;
            };
            for junction in junctions {
                if !e.bound.expanded(junction.reach).contains(junction.at) {
                    continue;
                }
                let (gap, t) = nearest_on_edge_curve(&e.curve, e.crange, junction.at, tol)?;
                if gap > junction.reach {
                    continue;
                }
                let paved = paves.get(&e.node).is_some_and(|list| {
                    list.iter().any(|p| {
                        e.curve
                            .point_at(p.t, tol)
                            .is_ok_and(|q| junction.at.distance(q) <= junction.reach)
                    })
                });
                if paved {
                    continue;
                }
                let columns = [
                    e.pcurve.point_at(rescale(t, e.crange, e.prange), tol)?,
                    other.point_at(rescale(t, e.crange, *orange), tol)?,
                ];
                let mut sides = [false; 2];
                for (at, chart) in &ends {
                    if junction.at.distance(*at) <= junction.reach {
                        let nearer =
                            usize::from(chart.distance(columns[1]) < chart.distance(columns[0]));
                        sides[nearer] = true;
                    }
                }
                if sides == [true; 2] {
                    ogeom_core::diag!(
                        Channel::Wire,
                        "PAVE seam e{} at {t:.6}: a junction at {:?} holds ends on both sides",
                        e.node.index(),
                        junction.at
                    );
                    paves.entry(e.node).or_default().push(Pave {
                        t,
                        honesty: gap,
                        across: None,
                        spread: 0.0,
                    });
                }
            }
        }
    }
    Ok(())
}

/// Junctions whose balls overlap, merged transitively into one junction
/// each: at the members' centroid, reaching as far as the farthest member's
/// own reach extends from it.
fn merge_junctions(junctions: Vec<Junction>) -> Vec<Junction> {
    let n = junctions.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    // Swept along x: two balls overlap only if their centres are within
    // both reaches along every axis, so each junction is asked only of those
    // ahead of it by no more than its reach and the widest. The groups are
    // the same as asking every pair, which for the tens of thousands of
    // junctions a converted part carries costs seconds.
    let widest = junctions.iter().map(|j| j.reach).fold(0.0_f64, f64::max);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| junctions[x].at.x.total_cmp(&junctions[y].at.x));
    for (k, &i) in order.iter().enumerate() {
        let limit = junctions[i].at.x + junctions[i].reach + widest;
        for &j in &order[k + 1..] {
            if junctions[j].at.x > limit {
                break;
            }
            if junctions[i].at.distance(junctions[j].at) <= junctions[i].reach + junctions[j].reach
            {
                let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
                if ri != rj {
                    parent[rj] = ri;
                }
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut of_root: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let r = root(&mut parent, i);
        match of_root[r] {
            Some(g) => groups[g].push(i),
            None => {
                of_root[r] = Some(groups.len());
                groups.push(vec![i]);
            }
        }
    }
    // A group is one junction only while it stays a junction's size: a
    // chain of overlapping balls along an edge would otherwise merge into
    // one reaching the chain's length, and the rebuild would weld the edge
    // into a point. A group reaching beyond four of its widest member stays
    // as it was.
    let mut out = Vec::with_capacity(groups.len());
    for members in groups {
        if members.len() == 1 {
            out.push(junctions[members[0]]);
            continue;
        }
        let mut sum = ogeom_math::Vector::ZERO;
        for &m in &members {
            sum += junctions[m].at - Point::ORIGIN;
        }
        #[allow(clippy::cast_precision_loss)]
        let at = Point::ORIGIN + sum / (members.len() as f64);
        let reach = members
            .iter()
            .map(|&m| at.distance(junctions[m].at) + junctions[m].reach)
            .fold(0.0_f64, f64::max);
        let widest = members
            .iter()
            .map(|&m| junctions[m].reach)
            .fold(0.0_f64, f64::max);
        if reach > widest * 4.0 {
            out.extend(members.iter().map(|&m| junctions[m]));
        } else {
            out.push(Junction {
                at,
                reach,
                onto_vertex: members.iter().any(|&m| junctions[m].onto_vertex),
            });
        }
    }
    out
}

/// The junctions the paves describe more than once, or more loosely than
/// the edge itself, per edge in face order.
fn pave_junctions(
    ga: &GSolid,
    gb: &GSolid,
    paves: &ogeom_core::FastMap<EdgeKey, Vec<Pave>>,
    along: &AlongEdges,
    tol: Tolerances,
) -> OgeomResult<Vec<Junction>> {
    let mut seen: ogeom_core::FastSet<EdgeKey> = ogeom_core::FastSet::default();
    let mut junctions = Vec::new();
    for e in ga
        .faces
        .iter()
        .chain(gb.faces.iter())
        .flat_map(|f| f.edges.iter())
    {
        if !seen.insert(e.node) {
            continue;
        }
        let Some(ts) = paves.get(&e.node) else {
            continue;
        };
        // A tolerant edge (a merged facet group's boundary, a fitted rail)
        // is met by a section up to its own radius off it, so the section's
        // end and the edge's split point are one junction that far apart,
        // on either side of the edge.
        let tolerant = e.tolerance > tol.confusion() * 1e2;
        let floor = if tolerant {
            e.tolerance * 2.0
        } else {
            e.tolerance
        }
        .max(tol.confusion() * 10.0);
        let off = |face, p| off_face(ga, gb, face, p, tol);
        let along = along.get(&e.node).map_or(&[][..], Vec::as_slice);
        for cluster in cluster_paves(&e.curve, e.crange, e.tolerance, ts, Some(&off), along, tol)? {
            if cluster.members > 1 || cluster.honesty > tol.confusion() * 1e2 || tolerant {
                junctions.push(Junction {
                    at: cluster.at,
                    reach: cluster.span + cluster.honesty.max(floor),
                    onto_vertex: false,
                });
            }
        }
    }
    Ok(junctions)
}

/// How long an edge of a face left alone must be: a shorter one may be
/// dust some other chart collapses, which only the face's own split finds.
const QUIET_EDGE: f64 = 1e-2;

/// The one piece of a face the other solid leaves alone: its own wires as
/// rings, its state read at the middle of one of its edges, which no
/// section touches and which therefore stands inside the other solid or
/// out of it with the whole face, or [`PieceState::Unread`] where there is no
/// `boundary` to read it against. `None` where that edge reads on the other
/// boundary after all, or a wire's edge is not one the face was gathered
/// with.
fn quiet_piece(
    model: &Model,
    from_a: bool,
    fi: usize,
    face: &GFace,
    boundary: Option<&ogeom_algo::SolidBoundary>,
    tol: Tolerances,
) -> OgeomResult<Option<FacePiece>> {
    let forward = face.face.oriented(ogeom_topo::Orientation::Forward);
    // A face with many holes is looked up by key, not walked once per edge.
    let indexed: Option<ogeom_core::FastMap<EdgeKey, usize>> = (face.edges.len() > 16).then(|| {
        let mut first =
            ogeom_core::FastMap::with_capacity_and_hasher(face.edges.len(), Default::default());
        for (index, e) in face.edges.iter().enumerate() {
            first.entry(e.node).or_insert(index);
        }
        first
    });
    let aside: ogeom_core::FastSet<ogeom_topo::TShapeId> = face
        .holes_aside
        .iter()
        .map(|hole| hole.wire.node())
        .collect();
    let mut rings = Vec::new();
    for wire in model.ordered_children_of(&forward)? {
        if aside.contains(&wire.node()) {
            continue;
        }
        let mut ring = Vec::new();
        for edge in model.ordered_children_of(&wire)? {
            let key = EdgeKey::of(&edge);
            if face.slits.contains(&key) {
                continue;
            }
            let found = match &indexed {
                Some(first) => first.get(&key).copied(),
                None => face.edges.iter().position(|e| e.node == key),
            };
            let Some(index) = found else {
                return Ok(None);
            };
            ring.push(Traversal {
                tag: Tag::Boundary {
                    edge: index,
                    range: face.edges[index].crange,
                },
                reversed: edge.orientation() == ogeom_topo::Orientation::Reversed,
                nodes: None,
            });
        }
        rings.push(ring);
    }
    let Some(e) = face.edges.first() else {
        return Ok(None);
    };
    let probe = e
        .curve
        .point_at(f64::midpoint(e.crange.0, e.crange.1), tol)?;
    let Some(boundary) = boundary else {
        return Ok(Some(FacePiece {
            from_a,
            face: fi,
            rings,
            outlines: Vec::new(),
            probe,
            state: PieceState::Unread,
            covered: false,
            whole: true,
        }));
    };
    let state = match boundary.holds(model, probe, tol) {
        Ok(Containment::In) => PieceState::In,
        Ok(Containment::Out) => PieceState::Out,
        Ok(Containment::On) | Err(ogeom_core::OgeomError::NotDone(_)) => return Ok(None),
        Err(e) => return Err(e),
    };
    Ok(Some(FacePiece {
        from_a,
        face: fi,
        rings,
        outlines: Vec::new(),
        probe,
        state,
        covered: false,
        whole: true,
    }))
}

/// A face's strands with its lone holes taken out: the rest, the place of
/// each among all of them, and the holes, each walked alone.
type Apart = (Vec<Strand<Tag>>, Vec<usize>, Vec<arrange::Lone<Tag>>);

/// The holes of a planar face whose strands stand apart from every other
/// strand, each walked alone (see [`arrange::Lone`]), with the rest of the
/// strands and their places among all of them; `None` where no hole does.
///
/// A hole stands apart when no strand of another wire ends within a few
/// snaps of its box, no section, contact or pole strand's box comes that
/// near it, and its box lies inside the outer wire's.
fn lone_holes(
    model: &Model,
    face: &GFace,
    strands: &[Strand<Tag>],
    snap: f64,
) -> OgeomResult<Option<Apart>> {
    // Fewer strands than this are walked together for no cost worth
    // sparing.
    const FEW: usize = 64;
    if strands.len() < FEW
        || !matches!(face.surface, SurfaceGeometry::Plane(_))
        || !face.poles.is_empty()
    {
        return Ok(None);
    }
    // Each edge's wire, `None` for an edge in no wire or in two.
    let index: ogeom_core::FastMap<EdgeKey, usize> = face
        .edges
        .iter()
        .enumerate()
        .map(|(i, e)| (e.node, i))
        .collect();
    let mut wire_of: Vec<Option<Option<usize>>> = vec![None; face.edges.len()];
    let wires = model.ordered_children_of(&face.face.oriented(ogeom_topo::Orientation::Forward))?;
    for (w, wire) in wires.iter().enumerate() {
        for edge in model.ordered_children_of(wire)? {
            if let Some(&i) = index.get(&EdgeKey::of(&edge)) {
                wire_of[i] = match wire_of[i] {
                    None => Some(Some(w)),
                    Some(_) => Some(None),
                };
            }
        }
    }
    let strand_wire: Vec<Option<usize>> = strands
        .iter()
        .map(|st| match st.tag {
            Tag::Boundary { edge, .. } if st.boundary => {
                wire_of.get(edge).copied().flatten().flatten()
            }
            _ => None,
        })
        .collect();
    let mut by_wire: Vec<Vec<usize>> = vec![Vec::new(); wires.len()];
    for (k, w) in strand_wire.iter().enumerate() {
        if let Some(w) = w {
            by_wire[*w].push(k);
        }
    }
    let boxed = |ks: &[usize]| -> (Point2, Point2) {
        let mut lo = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &k in ks {
            for p in &strands[k].polyline {
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
        (lo, hi)
    };
    let Some(outer) = by_wire.first().filter(|ks| !ks.is_empty()) else {
        return Ok(None);
    };
    let (outer_lo, outer_hi) = boxed(outer);
    // Every strand end, by abscissa; and the box of every strand no wire
    // owns.
    let mut ends: Vec<(f64, f64, usize)> = Vec::with_capacity(strands.len() * 2);
    let mut loose: Vec<(Point2, Point2)> = Vec::new();
    for (k, st) in strands.iter().enumerate() {
        for p in [st.polyline.first(), st.polyline.last()]
            .into_iter()
            .flatten()
        {
            ends.push((p.x, p.y, k));
        }
        if strand_wire[k].is_none() {
            loose.push(boxed(&[k]));
        }
    }
    ends.sort_by(|a, b| a.0.total_cmp(&b.0));
    let margin = snap * 4.0;
    let mut lone = Vec::new();
    let mut taken = vec![false; strands.len()];
    for (w, ks) in by_wire.iter().enumerate().skip(1) {
        if ks.is_empty() {
            continue;
        }
        let (lo, hi) = boxed(ks);
        let (lo, hi) = (
            Point2::new(lo.x - margin, lo.y - margin),
            Point2::new(hi.x + margin, hi.y + margin),
        );
        if !(lo.x > outer_lo.x && lo.y > outer_lo.y && hi.x < outer_hi.x && hi.y < outer_hi.y) {
            continue;
        }
        let first = ends.partition_point(|e| e.0 < lo.x);
        let crowded = ends[first..]
            .iter()
            .take_while(|e| e.0 <= hi.x)
            .any(|&(_, y, k)| y >= lo.y && y <= hi.y && strand_wire[k] != Some(w));
        let crossed = loose
            .iter()
            .any(|(a, b)| a.x <= hi.x && b.x >= lo.x && a.y <= hi.y && b.y >= lo.y);
        if crowded || crossed {
            continue;
        }
        let mine: Vec<Strand<Tag>> = ks.iter().map(|&k| strands[k].clone()).collect();
        let Some(ring) = arrange::lone_ring(&mine, ks, snap) else {
            continue;
        };
        for &k in ks {
            taken[k] = true;
        }
        lone.push(ring);
    }
    if lone.is_empty() {
        return Ok(None);
    }
    let (rest, places): (Vec<Strand<Tag>>, Vec<usize>) = strands
        .iter()
        .enumerate()
        .filter(|(k, _)| !taken[*k])
        .map(|(k, st)| (st.clone(), k))
        .unzip();
    Ok(Some((rest, places, lone)))
}

/// The face's outward normal at a chart point: the surface's, flipped when
/// the face presents its other side.
fn outward_normal(face: &GFace, at: Point2, tol: Tolerances) -> OgeomResult<ogeom_math::Vector> {
    let n = face.surface.normal_at(at.x, at.y, tol)?.vector();
    Ok(
        if face.face.orientation() == ogeom_topo::Orientation::Reversed {
            -n
        } else {
            n
        },
    )
}

/// How far apart two of a face's edges may honestly end in its chart: the
/// loosest edge's tolerance, or the hug width where an edge is a fitted
/// section that once ran along another before parting from it.
fn outline_snap(face: &GFace, tol: Tolerances) -> f64 {
    face.edges
        .iter()
        .fold(tol.confusion() * 1e2, |acc, e| acc.max(e.tolerance * 2.0))
        .max(
            if face.edges.iter().any(|e| e.tolerance > tol.confusion()) {
                tol.confusion() * 1e3
            } else {
                0.0
            },
        )
}

/// Polyline ends sorted by their chart `x`, for finding the ends within a
/// reach of a point without comparing every pair.
struct EndsByX {
    /// Indices into the ends, by ascending `x`.
    order: Vec<usize>,
    xs: Vec<f64>,
}

impl EndsByX {
    fn new(ends: &[(usize, bool, Point2)]) -> Self {
        let mut order: Vec<usize> = (0..ends.len()).collect();
        order.sort_by(|&a, &b| ends[a].2.x.total_cmp(&ends[b].2.x));
        let xs = order.iter().map(|&k| ends[k].2.x).collect();
        Self { order, xs }
    }

    /// The ends whose `x` lies within `reach` of `p`'s, in ascending order
    /// of their index: a superset of the ends within `reach` of `p`.
    fn near(&self, p: Point2, reach: f64) -> Vec<usize> {
        let lo = self.xs.partition_point(|&x| x < p.x - reach);
        let hi = self.xs.partition_point(|&x| x <= p.x + reach);
        let mut out = self.order[lo..hi.max(lo)].to_vec();
        out.sort_unstable();
        out
    }
}

/// Join the outline's ends the weld left apart, each to the unmatched end
/// nearest it where that end is nearest it too and the opening is small
/// against the outline.
///
/// The topology says a face's boundary closes. Where two edges' pcurves end
/// further apart in the chart than their tolerances reach (a fitted blend
/// edge meeting a column on a narrow torus wedge, a pole's row with no edge
/// along it), a ray cast for containment passes through the opening, counts
/// one crossing, and a point level with it reads inside. A short segment
/// across each opening closes the outline as the topology does.
fn bridge_outline_gaps(lines: &mut Vec<Vec<Point2>>) {
    let ends: Vec<(usize, bool, Point2)> = lines
        .iter()
        .enumerate()
        .filter(|(_, line)| line.len() >= 2)
        .flat_map(|(i, line)| [(i, false, line[0]), (i, true, line[line.len() - 1])])
        .collect();
    // Ends that meet within the parametric snap but not exactly are joined
    // where they stand: the weld leaves them, since moving one by rounding
    // noise takes a probe on a seam off it, but a ray level with the sliver
    // between them still slips through.
    let by_x = EndsByX::new(&ends);
    let mut joins: Vec<Vec<Point2>> = Vec::new();
    for (k, &(i, end_i, p)) in ends.iter().enumerate() {
        for m in by_x.near(p, PARAM_SNAP) {
            let (j, end_j, q) = ends[m];
            let d = p.distance(q);
            if m > k && (j, end_j) != (i, end_i) && d > 0.0 && d <= PARAM_SNAP {
                joins.push(vec![p, q]);
            }
        }
    }
    let open: Vec<Point2> = ends
        .iter()
        .filter(|&&(i, end_i, p)| {
            !by_x.near(p, PARAM_SNAP).into_iter().any(|m| {
                let (j, end_j, q) = ends[m];
                (j, end_j) != (i, end_i) && p.distance(q) <= PARAM_SNAP
            })
        })
        .map(|&(_, _, p)| p)
        .collect();
    // Only an opening small against the outline itself: a pole's row
    // spans the chart and is the pole's to bound, not a gap.
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in lines.iter().flatten() {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let reach = lo.distance(hi) * 1e-3;
    let nearest = |k: usize| {
        (0..open.len()).filter(|&m| m != k).min_by(|&a, &b| {
            open[k]
                .distance(open[a])
                .total_cmp(&open[k].distance(open[b]))
        })
    };
    for k in 0..open.len() {
        if let Some(m) = nearest(k)
            && k < m
            && nearest(m) == Some(k)
            && open[k].distance(open[m]) <= reach
        {
            lines.push(vec![open[k], open[m]]);
        }
    }
    lines.extend(joins);
}

/// Close the gaps between a face's outline polylines.
///
/// The edges are polylined one by one, and where two of them meet at a
/// vertex that owns some tolerance (a fitted section's end welded to the
/// edge it hugged, a few dozen microns off), their polylines stop that far
/// apart. A ray cast for containment slips through such a gap, and a probe
/// standing inside the face within a hug's width of the seam reads as
/// outside. Ends within the face's snap ([`outline_snap`]) of another
/// edge's end are made one point, so the outline is closed exactly as the
/// topology says it is.
///
/// The snap is a length in space, and so is every gap it is held against:
/// two ends meet only where their points on the surface do too. A chart
/// stretched along one direction (a wall a micron tall charted over a
/// tenth of a unit) holds a whole edge within a chart distance the snap
/// would read as a gap. Nor does the weld reach as far as half the face's
/// shortest edge, measured the same way: the two ends of that edge would
/// then read as one point, and the outline would close over nothing. Two
/// ends of one edge are never joined.
fn weld_outline_ends(lines: &mut [Vec<Point2>], face: &GFace, tol: Tolerances) {
    let snap = outline_snap(face, tol);
    if snap <= 0.0 {
        return;
    }
    // Which edge drew each line: a seam draws two, one per column.
    let owner: Vec<usize> = face
        .edges
        .iter()
        .enumerate()
        .flat_map(|(i, e)| std::iter::repeat_n(i, 1 + usize::from(e.other_side.is_some())))
        .collect();
    if owner.len() != lines.len() {
        return;
    }
    let lift = |p: Point2| face.surface.point_at(p.x, p.y, tol).ok();
    // How far apart in space two ends may be welded: the snap, and no
    // further than half the face's shortest edge there. A line with no
    // length in space (a pole's or an apex's row) has no ends to tell
    // apart.
    let mut within = snap;
    for line in lines.iter() {
        let points: Vec<Point> = line.iter().filter_map(|p| lift(*p)).collect();
        let length: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
        if points.len() == line.len() && length > tol.confusion() {
            within = within.min(length / 2.0);
        }
    }
    // Every polyline end, as (line, is its last point, where), and its
    // point on the surface.
    let ends: Vec<(usize, bool, Point2)> = lines
        .iter()
        .enumerate()
        .flat_map(|(i, line)| {
            let first = line.first().map(|p| (i, false, *p));
            let last = line.last().map(|p| (i, true, *p));
            first.into_iter().chain(last)
        })
        .collect();
    let lifted: Vec<Option<Point>> = ends.iter().map(|&(_, _, p)| lift(p)).collect();
    // Each end moves onto the nearest other end within reach that sorts
    // before it, so a matched pair lands on one point rather than trading
    // places. Within reach in the chart, as a gap the snap is read across
    // there, and on the surface.
    let by_x = EndsByX::new(&ends);
    let mut moves: Vec<(usize, bool, Point2)> = Vec::new();
    for (k, &(i, end_i, p)) in ends.iter().enumerate() {
        let Some(at) = lifted[k] else {
            continue;
        };
        let nearest = by_x
            .near(p, snap)
            .into_iter()
            .filter_map(|m| {
                let (j, end_j, q) = ends[m];
                let d = p.distance(q);
                // Ends closer than the parametric snap already meet for
                // every purpose here, and moving one by rounding noise
                // takes a probe that stands exactly on a chart's seam off
                // it.
                if (j, end_j) >= (i, end_i) || owner[j] == owner[i] || d <= PARAM_SNAP || d > snap {
                    return None;
                }
                (lifted[m]?.distance(at) <= within).then_some((d, q))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        if let Some((_, q)) = nearest {
            moves.push((i, end_i, q));
        }
    }
    for (i, end_i, q) in moves {
        let slot = if end_i {
            lines[i].last_mut()
        } else {
            lines[i].first_mut()
        };
        if let Some(p) = slot {
            *p = q;
        }
    }
}

/// A curve's length and the diagonal of the box its points span, both from
/// the same 64 samples: a trace lapping one place many times is far longer
/// than its box, and a loop is about a turn round it.
fn length_and_reach(curve: &Curve, tol: Tolerances) -> OgeomResult<(f64, f64)> {
    use ogeom_geom::Curve3d as _;
    let (lo, hi) = curve.domain();
    let mut length = 0.0_f64;
    let mut last: Option<Point> = None;
    let mut span = ogeom_math::Aabb::default();
    for k in 0..=64 {
        let t = if k == 64 {
            hi
        } else {
            lo + (hi - lo) * f64::from(k) / 64.0
        };
        let p = curve.point_at(t, tol)?;
        if let Some(q) = last {
            length += q.distance(p);
        }
        last = Some(p);
        span = span.with_point(p);
    }
    Ok((length, span.diagonal()))
}

/// Whether a point a face's surface holds to within `margin` lies inside
/// the face's trim, or within `margin` of its boundary as the face draws it.
fn within_face(face: &GFace, p: Point, margin: f64, tol: Tolerances) -> OgeomResult<bool> {
    if chart_point_within(face, p, margin.max(tol.confusion() * 1e2), tol).is_some() {
        return Ok(true);
    }
    for e in &face.edges {
        if !e.bound.expanded(margin).contains(p) {
            continue;
        }
        let (drawn, drange, _) = e.drawn();
        if distance_to_edge_curve(drawn, drange, p, tol)? <= margin {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Where a point sits in a face's chart, for a point up to `reach` off the
/// face's surface and inside its trim. `None` elsewhere.
fn chart_point_within(face: &GFace, p: Point, reach: f64, tol: Tolerances) -> Option<Point2> {
    // Closed-form inversion for the analytic surfaces: the same-domain
    // resolution asks "where does this probe sit in the partner's chart", and
    // the partner may be any surface a face melts along: a plane against a
    // plane, but equally a wall band against the cylinder it copies. The
    // reach check keeps the answer honest: a point off the surface has no
    // chart position, whatever the inversion returns.
    use ogeom_math::elementary;
    let raw = match &face.surface {
        SurfaceGeometry::Plane(x) => {
            let local = x.plane().frame().to_local(p);
            if local.z.abs() > reach {
                return None;
            }
            Point2::new(local.x, local.y)
        }
        SurfaceGeometry::Cylinder(x) => {
            let cylinder = x.cylinder();
            if cylinder.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::cylinder_parameters(&cylinder, p, tol).ok()?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Cone(x) => {
            let cone = x.cone();
            if cone.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::cone_parameters(&cone, p, tol).ok()?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Sphere(x) => {
            let sphere = x.sphere();
            if sphere.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::sphere_parameters(&sphere, p, tol).ok()?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Torus(x) => {
            let torus = x.torus();
            if torus.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::torus_parameters(&torus, p, tol).ok()?;
            Point2::new(u, v)
        }
        // No closed form (a fitted patch, a swept or revolved surface)
        // inverts by projection, held to the same reach: a blend's leg on
        // a spline host is a partner like any other.
        other => {
            let foot = ogeom_algo::project_on_surface(other, p, 24, tol).ok()?;
            if foot.distance > reach {
                return None;
            }
            Point2::new(foot.parameters.0, foot.parameters.1)
        }
    };
    let at = fold_point_into_chart(raw, &face.surface, tol);
    let lines = face.outline(tol)?;
    let trim = face.outline_trim(tol)?;
    // The face's boundary polylines are unwrapped (a winding ring may span
    // any one period's window, not necessarily the chart's canonical one,
    // and a wire chained onto one branch of the chart may sit whole periods
    // away from it), so the probe is carried to the outline's own branch
    // first, and asked at the neighbouring images as well.
    let mut shifts = vec![0.0];
    if face.surface.is_periodic_u() {
        let ((ua, ub), _) = face.surface.domain();
        if ub > ua {
            let period = ub - ua;
            let (sum, count) = lines
                .iter()
                .flatten()
                .fold((0.0_f64, 0_usize), |(s, n), q| (s + q.x, n + 1));
            #[allow(clippy::cast_precision_loss)]
            let centre = if count > 0 { sum / count as f64 } else { at.x };
            let home = ((centre - at.x) / period).round() * period;
            shifts = vec![home, home + period, home - period];
        }
    }
    for shift in shifts {
        let shifted = Point2::new(at.x + shift, at.y);
        if trim.inside(shifted) {
            return Some(shifted);
        }
    }
    None
}

/// Whether `p`, on the face's surface, stands inside the face, read at
/// the face's edge nearest it: where the nearest point of the boundary is
/// inside an edge, `p` stands off it square to it, and is inside just
/// where it stands on the side the face keeps its material on, the left
/// of the edge as the face walks it about its outward normal. `None`
/// where the nearest point is an edge's end (a corner reads both sides),
/// a seam's, or within the edge's own radius of it.
fn inside_beside_edge(face: &GFace, p: Point, tol: Tolerances) -> OgeomResult<Option<bool>> {
    let mut nearest: Option<(f64, usize, f64)> = None;
    for (k, e) in face.edges.iter().enumerate() {
        let (off, t) = nearest_on_edge_curve(&e.curve, e.crange, p, tol)?;
        if nearest.is_none_or(|n| off < n.0) {
            nearest = Some((off, k, t));
        }
    }
    let Some((off, k, t)) = nearest else {
        return Ok(None);
    };
    let e = &face.edges[k];
    let (lo, hi) = (e.crange.0.min(e.crange.1), e.crange.0.max(e.crange.1));
    let Some(forward) = e.walked else {
        return Ok(None);
    };
    if off <= e.tolerance.max(tol.confusion())
        || t <= lo + tol.parametric()
        || t >= hi - tol.parametric()
    {
        return Ok(None);
    }
    let at = e.pcurve.point_at(rescale(t, e.crange, e.prange), tol)?;
    let tangent = e.curve.d1_at(t, tol)?;
    let walked = if forward { tangent } else { -tangent };
    let left = outward_normal(face, at, tol)?.cross(walked);
    Ok(Some((p - e.curve.point_at(t, tol)?).dot(left) > 0.0))
}

/// Whether the point of `face`'s chart at `uv` lies on the face, read as a
/// partner's probe is read on it: a probe a piece's coarse outline puts
/// outside its face (beyond the chords' bow) is not the piece's.
fn on_own_face(face: &GFace, uv: Point2, tol: Tolerances) -> OgeomResult<bool> {
    let at = face.surface.point_at(uv.x, uv.y, tol)?;
    Ok(chart_point_within(face, at, partner_reach(face, face, tol), tol).is_some())
}

/// What the strands bounding a piece say of it, where its probes cannot.
enum Thin {
    /// The piece is no thinner than the band, or not planar: its probes
    /// stand.
    Wide,
    /// Thin, but some strand bounding it says nothing.
    Undecided,
    /// Thin, and every strand the other solid draws round it agrees.
    Reads(PieceState),
    /// Thin, and the strands disagree.
    Disagree,
}

/// A thin planar piece read by the other solid's edges it is bounded by.
///
/// A piece whose area is at most the classification band times half its
/// perimeter is on average narrower than the band, and no probe of it can
/// stand outside the band: each reads it on the other solid's boundary. A
/// strand that is the image of the other solid's edge `E`, between its
/// faces `g` on this face's surface and `h`, separates "on `g`" from
/// "beyond `g`". Which side the piece lies on is the arrangement's: its
/// rings keep it on their left about the chart's normal, as `g` keeps its
/// material on the left of its own rings about its outward normal. Beyond
/// `g` the other solid's wedge at `E` answers: past a convex edge (where
/// `h` turns away behind `g`) the plane beyond `g` is outside it, past a
/// concave one inside. Neither answer depends on how thin the piece is.
///
/// Only contact strands speak: a piece a section bounds is left to its
/// probes, as is one beside an edge whose two faces run on tangentially,
/// one it lies on `g`'s side of (a coincidence, which the partner reads),
/// or one whose faces are not both gathered.
fn thin_by_strands(
    face: &GFace,
    rings: &[Vec<Traversal<Tag>>],
    outlines: &[Vec<Point2>],
    contacts: &[ContactRec],
    theirs: &GSolid,
    tol: Tolerances,
) -> OgeomResult<Thin> {
    if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
        return Ok(Thin::Wide);
    }
    // A plane's chart is its own frame, so lengths and areas in it are
    // the plane's.
    let (mut area, mut perimeter) = (0.0, 0.0);
    for (k, ring) in outlines.iter().enumerate() {
        let n = ring.len();
        let mut twice = 0.0;
        for i in 0..n {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            twice += a.x * b.y - b.x * a.y;
            perimeter += (b - a).magnitude();
        }
        area += if k == 0 { twice.abs() } else { -twice.abs() } / 2.0;
    }
    if area > tol.confusion() * 1e4 * perimeter / 2.0 {
        return Ok(Thin::Wide);
    }
    let mut read: Option<PieceState> = None;
    let mut spoke = false;
    for traversal in rings.iter().flatten() {
        let (contact, range) = match &traversal.tag {
            Tag::Contact { contact, range } => (&contacts[*contact], *range),
            Tag::Section { .. } => return Ok(Thin::Undecided),
            Tag::Boundary { .. } | Tag::Pole { .. } => continue,
        };
        let Some(state) = beyond_edge(face, contact, range, traversal.reversed, theirs, tol)?
        else {
            return Ok(Thin::Undecided);
        };
        spoke = true;
        match read {
            Some(said) if said != state => return Ok(Thin::Disagree),
            _ => read = Some(state),
        }
    }
    Ok(match read {
        Some(state) if spoke => Thin::Reads(state),
        _ => Thin::Undecided,
    })
}

/// Where the piece on the left of the contact strand over `range`, walked
/// back where `reversed`, stands against the solid whose edge the contact
/// is: `None` where the strand says nothing (see [`thin_by_strands`]).
fn beyond_edge(
    face: &GFace,
    contact: &ContactRec,
    range: (f64, f64),
    reversed: bool,
    theirs: &GSolid,
    tol: Tolerances,
) -> OgeomResult<Option<PieceState>> {
    let t = f64::midpoint(range.0, range.1);
    let p = contact.curve.point_at(t, tol)?;
    let along = contact.curve.d1_at(t, tol)?;
    let walked = if (range.1 < range.0) != reversed {
        -along
    } else {
        along
    };
    let uv = contact
        .pcurve
        .point_at(rescale(t, contact.crange, contact.prange), tol)?;
    let chart = face.surface.normal_at(uv.x, uv.y, tol)?.vector();
    let Ok(left) = chart.cross(walked).normalized(tol) else {
        return Ok(None);
    };
    // The edge's two faces, each with the way into its material from the
    // edge and its outward normal there.
    let mut sides = Vec::with_capacity(2);
    for other in &theirs.faces {
        for e in other.edges.iter().filter(|e| e.node == contact.node) {
            let Some(forward) = e.walked else {
                return Ok(None);
            };
            let (_, te) = nearest_on_edge_curve(&e.curve, e.crange, p, tol)?;
            let at = e.pcurve.point_at(rescale(te, e.crange, e.prange), tol)?;
            let normal = outward_normal(other, at, tol)?;
            let tangent = e.curve.d1_at(te, tol)?;
            let tangent = if forward { tangent } else { -tangent };
            let Ok(into) = normal.cross(tangent).normalized(tol) else {
                return Ok(None);
            };
            let Ok(normal) = normal.normalized(tol) else {
                return Ok(None);
            };
            sides.push((into, normal));
        }
    }
    let [(into_a, normal_a), (into_b, normal_b)] = sides[..] else {
        return Ok(None);
    };
    // `g` is the face lying along this one: its way in runs along the
    // piece's, one way or the other.
    let along_this = |into: ogeom_math::Vector| into.cross(left).magnitude() < TANGENT_TURN;
    let ((into_g, normal_g), into_h) = match (along_this(into_a), along_this(into_b)) {
        (true, false) => ((into_a, normal_a), into_b),
        (false, true) => ((into_b, normal_b), into_a),
        _ => return Ok(None),
    };
    if into_g.dot(left) > 0.0 {
        return Ok(None);
    }
    let turn = into_h.dot(normal_g);
    if turn.abs() < TANGENT_TURN {
        return Ok(None);
    }
    Ok(Some(if turn < 0.0 {
        PieceState::Out
    } else {
        PieceState::In
    }))
}

/// The sine of the turn under which two curves meet tangentially at a
/// crossing: the parallel-noise gate for pave placement. One degree is far
/// below any deliberate crossing and far above the smear a tangency leaves
/// in the general intersector.
pub(crate) const TANGENT_TURN: f64 = 2e-2;

/// Whether two curves meet tangentially at a crossing: their tangents turn
/// by less than [`TANGENT_TURN`].
fn tangential(a: &Curve, ta: f64, b: &Curve, tb: f64, tol: Tolerances) -> OgeomResult<bool> {
    let da = a.d1_at(ta, tol)?;
    let db = b.d1_at(tb, tol)?;
    let (ma, mb) = (da.magnitude(), db.magnitude());
    if ma <= tol.confusion() || mb <= tol.confusion() {
        return Ok(true);
    }
    Ok(da.cross(db).magnitude() / (ma * mb) < TANGENT_TURN)
}

/// The chart image of an exact curve on a surface, fitted at the curve's
/// own parameters from its points inverted on the surface, with the
/// distance the image honestly sits off the curve, or nothing, when the
/// fit misses its budget by more than a chart's swing beside a pole allows.
fn fitted_image(
    curve: &Curve,
    surface: &SurfaceGeometry,
    budget: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(PlanarCurve, f64)>> {
    use ogeom_geom::Surface as _;
    const SAMPLES: usize = 96;
    let (lo, hi) = curve.domain();
    if !hi.is_finite() || !lo.is_finite() || hi <= lo {
        return Ok(None);
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let periods = (
        if surface.is_periodic_u() {
            ub - ua
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            vb - va
        } else {
            0.0
        },
    );
    // The image sampled where it moves: beside a pole a circle's image
    // swings half a turn within a hair of arc, and a uniform sampling puts
    // one point in the swing. Between any two samples whose images sit
    // further apart than a small step of the chart a sample is added, until
    // the image is walked at that step or the budget of samples is spent.
    let invert = |t: f64, guess: Option<(f64, f64)>| -> OgeomResult<Option<((f64, f64), f64)>> {
        let p = curve.point_at(t, tol)?;
        // A plane's image is its frame's own coordinates, exact and
        // unbounded: the face's chart window is a trim, not a limit on
        // where the curve may be read, and a projection clamped to it puts
        // the circle's far side off the surface by the window's shortfall.
        if let SurfaceGeometry::Plane(plane) = surface {
            let local = plane.plane().frame().to_local(p);
            let back = plane
                .plane()
                .frame()
                .to_world(ogeom_math::Point::new(local.x, local.y, 0.0));
            return Ok(Some(((local.x, local.y), back.distance(p))));
        }
        // A sphere's image is its frame's longitude and latitude, exact
        // where the projection's Newton cannot refine beside a pole and
        // hands back a coarse scan sample instead. The reading is checked
        // round the trip and left to the projection if the chart's
        // convention is not this one.
        if let SurfaceGeometry::Sphere(ball) = surface {
            let frame = ball.sphere().frame();
            let local = frame.to_local(p);
            let flat = local.x.hypot(local.y);
            // Kept a hair inside the poles: a spline fitted through points
            // on the pole itself overshoots it between them, and the chart
            // refuses a latitude past its end.
            let half = core::f64::consts::FRAC_PI_2 - 1e-6;
            let v = local.z.atan2(flat).clamp(-half, half);
            let u = local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU);
            let back = surface.point_at(u, v, tol)?;
            if back.distance(p) <= budget.max(tol.confusion() * 1e2) {
                return Ok(Some(((u, v), back.distance(p))));
            }
        }
        // A cone's image is read in closed form too, past its window's ends:
        // an open conic it cuts is stated over a span that climbs beyond the
        // window, which the face trims, and a projection clamped to the
        // window puts those points off the surface.
        if let SurfaceGeometry::Cone(cone) = surface
            && let Ok((u, v)) = ogeom_math::elementary::cone_parameters(&cone.cone(), p, tol)
        {
            let back = ogeom_math::elementary::cone_at(&cone.cone(), u, v).point;
            if back.distance(p) <= budget.max(tol.confusion() * 1e2) {
                return Ok(Some(((u, v), back.distance(p))));
            }
        }
        let mut foot = match guess {
            Some(g) => ogeom_algo::project_on_surface_from(surface, p, g, tol)?,
            None => ogeom_algo::project_on_surface(surface, p, 24, tol)?,
        };
        if foot.distance > budget.max(tol.confusion() * 1e2) && guess.is_some() {
            foot = ogeom_algo::project_on_surface(surface, p, 48, tol)?;
        }
        if foot.distance > budget.max(tol.confusion() * 1e2) {
            ogeom_core::diag!(
                Channel::Wire,
                "IMAGE FIT: a point sits {:.2e} off the surface; no image",
                foot.distance
            );
            return Ok(None);
        }
        Ok(Some((foot.parameters, foot.distance)))
    };
    let unwrap = |uv: (f64, f64), before: Option<Point2>| -> Point2 {
        let mut uv = Point2::new(uv.0, uv.1);
        if let Some(prev) = before {
            for (coord, period, before) in [
                (&mut uv.x, periods.0, prev.x),
                (&mut uv.y, periods.1, prev.y),
            ] {
                if period > 0.0 {
                    while *coord - before > period / 2.0 {
                        *coord -= period;
                    }
                    while before - *coord > period / 2.0 {
                        *coord += period;
                    }
                }
            }
        }
        uv
    };
    let mut samples: Vec<(f64, Point2)> = Vec::with_capacity(SAMPLES + 1);
    let mut guess: Option<(f64, f64)> = None;
    for k in 0..=SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * (k as f64) / (SAMPLES as f64);
        let Some((uv, _)) = invert(t, guess)? else {
            return Ok(None);
        };
        guess = Some(uv);
        let before = samples.last().map(|(_, p)| *p);
        samples.push((t, unwrap(uv, before)));
    }
    let step = {
        let (su, sv) = (
            if periods.0 > 0.0 { periods.0 } else { ub - ua },
            if periods.1 > 0.0 { periods.1 } else { vb - va },
        );
        su.min(sv) * 0.01
    };
    const CAP: usize = 4096;
    let mut refined = true;
    while refined && samples.len() < CAP {
        refined = false;
        let mut next: Vec<(f64, Point2)> = Vec::with_capacity(samples.len() * 2);
        for pair in samples.windows(2) {
            let ((t0, p0), (t1, p1)) = (pair[0], pair[1]);
            next.push((t0, p0));
            if p0.distance(p1) > step && t1 - t0 > (hi - lo) * 1e-7 && next.len() < CAP {
                let tm = f64::midpoint(t0, t1);
                let Some((uv, _)) = invert(tm, Some((p0.x, p0.y)))? else {
                    return Ok(None);
                };
                next.push((tm, unwrap(uv, Some(p0))));
                refined = true;
            }
        }
        next.push(samples[samples.len() - 1]);
        samples = next;
    }
    let params: Vec<f64> = samples.iter().map(|(t, _)| *t).collect();
    let image: Vec<Point2> = samples.iter().map(|(_, p)| *p).collect();
    let fitted =
        ogeom_geom::fit::fit_points_2d_at(&params, &image, 3, tol.confusion() * 10.0, tol)?;
    let planar: PlanarCurve = fitted.curve.into();
    // The honest error is in space: the surface read through the image
    // against the curve itself, between the samples as well as at them.
    // Checked at the samples and between each pair of them: where the
    // image swings the samples are dense, and a uniform check would step
    // over the swing and read the fit as honest.
    // An image that leaves the chart between its samples (a spline
    // overshooting the pole it was fitted up to) is no image either.
    let mut off = 0.0_f64;
    let mut checks: Vec<f64> = Vec::with_capacity(params.len() * 2);
    for pair in params.windows(2) {
        checks.push(pair[0]);
        checks.push(f64::midpoint(pair[0], pair[1]));
    }
    if let Some(&t) = params.last() {
        checks.push(t);
    }
    // A cone is read past its window's ends, as its points were inverted.
    // So is a cylinder where its window refuses the point: a circle a
    // rounding off square to the axis climbs a hair past the wall's end
    // round half its turn, and the fit through points held at that end
    // follows it there.
    let lift = |u: f64, v: f64| match surface {
        SurfaceGeometry::Cone(cone) => {
            Ok(ogeom_math::elementary::cone_at(&cone.cone(), u, v).point)
        }
        SurfaceGeometry::Cylinder(drum) => surface
            .point_at(u, v, tol)
            .or_else(|_| Ok(ogeom_math::elementary::cylinder_at(&drum.cylinder(), u, v).point)),
        _ => surface.point_at(u, v, tol),
    };
    for t in checks {
        let uv = planar.point_at(t, tol)?;
        let Ok(on_surface) = lift(uv.x, uv.y) else {
            ogeom_core::diag!(
                Channel::Wire,
                "IMAGE FIT: the image leaves the chart at {uv:?}; no image"
            );
            return Ok(None);
        };
        off = off.max(on_surface.distance(curve.point_at(t, tol)?));
    }
    ogeom_core::diag!(
        Channel::Wire,
        "IMAGE FIT: image off the curve by {off:.2e} against {:.2e} (fit error {:.2e} met {})",
        budget.max(tol.confusion() * 1e3),
        fitted.error,
        fitted.met
    );
    if off > budget.max(tol.confusion() * 1e3) {
        return Ok(None);
    }
    Ok(Some((planar, off)))
}

/// The distance from a point to a bounded edge curve, through a sampling
/// fine enough for the along-boundary question it answers.
fn distance_to_edge_curve(
    curve: &Curve,
    crange: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    Ok(nearest_on_edge_curve(curve, crange, p, tol)?.0)
}

/// Absolute distance for the hug's strict radius test. Other sampling
/// queries compare the growth of their sampled distances, or use their
/// feet to seed paving; they retain their own sampling parameterization.
fn hug_distance(curve: &Curve, crange: (f64, f64), p: Point, tol: Tolerances) -> OgeomResult<f64> {
    let (_, seed) = nearest_on_edge_curve(curve, crange, p, tol)?;
    // A chord's foot is only a seed: its distance can be above or below
    // the distance to the curve, even for a point lying exactly on it.
    // Refine on the curve itself before asking the micron-wide hug filters.
    let fine = (crange.1 - crange.0) / (48.0 * 8.0);
    let distance_at = |t| curve.point_at(t, tol).map(|at| at.square_distance(p));
    let refined = ogeom_math::solve::minimize(
        |t| distance_at(t).unwrap_or(f64::INFINITY),
        (seed - fine).max(crange.0),
        (seed + fine).min(crange.1),
        ogeom_math::solve::Criteria {
            residual: 0.0,
            step: tol.parametric(),
            max_iterations: 100,
        },
    )?;
    let mut best = distance_at(seed)?;
    for t in [crange.0, crange.1, refined.value] {
        best = best.min(distance_at(t)?);
    }
    Ok(best.sqrt())
}

/// [`distance_to_edge_curve`], with the parameter of the nearest point.
fn nearest_on_edge_curve(
    curve: &Curve,
    crange: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    // Coarse bracket, then two rounds of local refinement: the answer feeds
    // the hug filters, whose widths are fractions of a millimetre, and a
    // long arc's 48-segment polyline sags by more than that on its own.
    let scan = |lo: f64, hi: f64, steps: u32| -> OgeomResult<(f64, f64)> {
        let mut best = f64::INFINITY;
        let mut best_t = lo;
        let mut previous: Option<(f64, Point)> = None;
        for i in 0..=steps {
            let t = lo + (hi - lo) * f64::from(i) / f64::from(steps);
            let at = curve.point_at(t, tol)?;
            if let Some((t0, last)) = previous {
                let d = at - last;
                let len2 = d.dot(d);
                let s = if len2 > 0.0 {
                    ((p - last).dot(d) / len2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let dist = p.distance(last + d * s);
                if dist < best {
                    best = dist;
                    best_t = (t - t0).mul_add(s, t0);
                }
            }
            previous = Some((t, at));
        }
        Ok((best, best_t))
    };
    let span = crange.1 - crange.0;
    let (_, t1) = scan(crange.0, crange.1, 48)?;
    let step = span / 48.0;
    let (_, t2) = scan((t1 - step).max(crange.0), (t1 + step).min(crange.1), 16)?;
    let fine = span / (48.0 * 8.0);
    scan((t2 - fine).max(crange.0), (t2 + fine).min(crange.1), 16)
}

/// Pair up the pieces two coincident faces contribute for the same patch of
/// one surface, and mark the second argument's copy as stood in for where
/// the material is on the same side.
///
/// Returns each pair, first argument's piece first, the opposed ones
/// included: whichever of the two the operation drops, its face is recorded
/// as having become the one it keeps.
///
/// The substitution is by *region*, not by argument: a piece of `b` is a
/// duplicate exactly where a piece of `a`, on the same side and on a surface
/// they share, already covers the point `b`'s piece was classified at. Where
/// `a` has nothing there (a bore refilled by the cylinder that cut it, whose
/// caps fill holes the part no longer has faces for), nothing stands in, and
/// `b`'s piece is the only description of that patch there is.
fn mark_covered_coincidences(
    ga: &GSolid,
    gb: &GSolid,
    pieces: &mut [FacePiece],
    tol: Tolerances,
) -> Vec<(usize, usize)> {
    let on = |state: PieceState| matches!(state, PieceState::OnAligned | PieceState::OnOpposed);
    // Which pieces of the first argument stand on a shared surface.
    let from_a: Vec<(usize, usize, PieceState)> = pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.from_a && on(p.state))
        .map(|(i, p)| (i, p.face, p.state))
        .collect();
    if from_a.is_empty() {
        return Vec::new();
    }
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for (index, piece) in pieces.iter().enumerate() {
        if piece.from_a || !on(piece.state) {
            continue;
        }
        for &(other, face_a, state) in &from_a {
            if state != piece.state {
                continue;
            }
            let host = &ga.faces[face_a];
            // The point `b`'s piece stands at, read in `a`'s face's chart,
            // as far off it as the two faces say they may stand. A point
            // that is not on that surface at all has no chart position.
            let reach = partner_reach(&gb.faces[piece.face], host, tol);
            let Some(at) = chart_point_within(host, piece.probe, reach, tol) else {
                continue;
            };
            // The host piece's rings, not its boundary strands: a ring is a
            // closed loop whose last point does not repeat its first, so the
            // test has to close it. Asked as though the rings were strands
            // that jointly close, the segment from the ring's end back to its
            // start goes uncounted, and a point the ring plainly encloses
            // comes back outside whenever that missing segment would have
            // been crossed.
            if inside_rings(&pieces[other].outlines, at) {
                pairs.push((other, index));
                break;
            }
        }
    }
    for &(_, index) in &pairs {
        if pieces[index].state == PieceState::OnAligned {
            pieces[index].covered = true;
        }
    }
    pairs
}

thread_local! {
    /// Whether a piece read on the other solid's boundary, with no coincident
    /// partner, is settled by asking just off it on both sides. Off on the
    /// first attempt, so a configuration the nested fallback answers exactly
    /// keeps that answer. On for the second. Read once per boolean on the
    /// calling thread.
    static SETTLE_FROM_BOTH_SIDES: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Whether a piece that reads in or out is read again at the centre of
    /// its area. Off but for the attempt after kept pieces fail to close.
    /// Read once per boolean on the calling thread.
    static PROBE_AT_CENTRE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The farthest apart a contact and a boundary edge it was read along
    /// stand, where that is past what either edge states: the paving takes
    /// the two for one edge and splits neither face along them, so only the
    /// sew can join them. Zero where every such pair lies within its stated
    /// tolerances. Gathered by the attempts since [`settled`] began.
    static ALONG_APART: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
    /// The gap the kept pieces are sewn across: zero but for the attempt
    /// [`settled`] makes last, across the distance [`ALONG_APART`] measured.
    /// Read once per assembly on the calling thread.
    static SEW_ACROSS: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
    /// Whether contacts are read at the radii their curves hold to (see
    /// [`settled`]). Off but for the last attempt. Read once per boolean
    /// on the calling thread.
    static AT_STATED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Which edges a section line also ends where it crosses, for the edges
    /// it carries onto the other face (see [`graze_contacts`]). Off but for
    /// the retries run again once every other has refused a piece read on
    /// the other solid's boundary. Read once per boolean on the calling
    /// thread.
    static FROM_CROSSINGS: std::cell::Cell<Crossings> = const { std::cell::Cell::new(Crossings::Off) };
    /// The farthest apart two faces may stand and still be dropped as a
    /// membrane (see [`without_membranes`]), whatever weld the attempt
    /// itself runs at: during the retries at a widened confusion, the weld
    /// their result is held to, and unbounded otherwise. A shell dropped
    /// whole leaves no edge or vertex behind to state how far apart its
    /// faces stood. Read once per assembly on the calling thread.
    static MEMBRANE_WITHIN: std::cell::Cell<f64> = const { std::cell::Cell::new(f64::INFINITY) };
    /// Whether a plane's section with an extrusion of a B-spline profile
    /// crossing a rim of it is split at the rim in closed form
    /// ([`ogeom_intersect::plane_extrusion_at_rims`]) rather than marched.
    /// Off but for the retries [`settled`] runs once every other has
    /// refused. Read once per assembly on the calling thread.
    static AT_RIMS: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Whether an assembly since [`settled`] began met such a pair with
    /// [`AT_RIMS`] off.
    static RIM_CROSSED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Whether two faces are asked whether they coincide within the
    /// tolerances they state ([`coincide_as_stated`]) where those are the
    /// default ones too. Off but for the attempt [`settled`] makes last.
    /// Read once per assembly on the calling thread.
    static AS_STATED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// An operation run once, and where it refuses, again in three steps. Where
/// the refusal is a piece on the other solid's boundary
/// with no coincident partner, once more with such pieces settled from both
/// sides of them. Where the kept pieces did not close, once more with each
/// piece read at the centre of its area: a face
/// leaning off the other's by microradians, the two meeting along an edge
/// and nothing splitting it, stands within rounding of the other near that
/// edge, and a probe there can read the wrong side. And where it still
/// refuses, at a confusion distance ten
/// and then a hundred times wider: a sliver a few microns thick (two walls
/// a few microns apart, a face tilted a few microradians) is then the kind
/// the boolean welds, and the result is kept only if nothing in it states
/// a tolerance past the weld distance, a hundred confusion distances (and
/// the ten the welding run states its own vertices to). A shell those
/// retries drop as a membrane leaves nothing behind to state that, so it
/// is dropped only where its faces stand within that same weld of each
/// other, and a sliver any thicker is never lost to the retry.
///
/// Then, at the given confusion, with contacts read at the radii their
/// curves hold to rather than at the loose width a fitted curve is
/// measured within: a fitted curve leaving an exact one tangentially
/// leaves a lune between them wider than those radii, which is geometry.
/// Each stretch the two share ends where they part by those radii,
/// measured once for the pair so both sides split at one vertex; a shallow
/// crossing splits both curves, at a vertex of either it cannot be told
/// from; a contact stretch near an edge is read inside or out at that
/// edge; a fitted section's hug width welds only a face it leaves a strand
/// on; and a piece thinner than the classification band is read by the
/// other solid's edges bounding it ([`thin_by_strands`]).
///
/// Last, where the refusal was a piece on the other solid's boundary and
/// every retry refused too, the same retries again with the line edges a
/// section line crosses carried onto the other face from the crossing on
/// ([`graze_contacts`]): a strip of one wall between such a crossing and
/// its edge, standing on another wall a hair's angle off it, is then split
/// off both and read as one contact. Where those refuse too, once more
/// with edges of any curve carried so, the result kept only where nothing
/// in it states a tolerance past the weld distance: an arc's strip is
/// shorter than a line's, and a junction welding its two ends together
/// would weld the faces beside it into one. What the retries without them
/// resolve is resolved as before.
///
/// Where all of that refuses and a plane crossed an extrusion of a
/// B-spline profile across a rim of it, whose marched section stops short
/// of the rim, the run once more from the first attempt with that section
/// split where it crosses the rim in closed form ([`AT_RIMS`]), its result
/// kept only where nothing in it states a tolerance past the weld
/// distance.
///
/// And where the kept pieces still did not close and the first attempt
/// read a contact along a boundary edge it stands further from than either
/// states (a floor edge lying its tolerance off a mirror plane, against its
/// mirror image as far off the other way), once more with the pieces sewn
/// across that distance ([`SEW_ACROSS`]), where it is within the weld
/// distance the retries above hold to: the paving took the two for one
/// edge and split neither face along them, so the edge kept widens to
/// reach where the other ran. The attempt runs at the given confusion, and
/// the sew widens nothing further than that distance beyond what the
/// pieces already state, so its result is not held to the weld distance a
/// second time: a part may itself state wider tolerances.
///
/// After all of that, where the refusal was a piece on the other solid's
/// boundary with no coincident partner and every retry above refused too,
/// once more with every pair of faces asked whether they coincide within
/// the tolerances they state, the default ones included ([`AS_STATED`]). A
/// face leaning a rounding off another it meets along an edge (a pocket's
/// floor sketched a tenth of a micron off the plane it stands on) parts
/// from it by a hair past the confusion distance at its far edge: the
/// closed form reads the two as crossing along the edge they share, which
/// splits neither, and the leaning face reads on the other's boundary
/// everywhere with no partner. Within what the two faces state they are one
/// surface, and each face is split by the other's boundary as a coincident
/// face is.
fn settled(
    model: &mut Model,
    tol: Tolerances,
    mut attempt: impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    RIM_CROSSED.set(false);
    ALONG_APART.set(0.0);
    let refusal = match attempt(model, tol) {
        Ok(built) => return Ok(built),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Err(e) => e,
    };
    let apart = ALONG_APART.get();
    if let Some(built) = resettled(model, tol, &refusal, &mut attempt)? {
        return Ok(built);
    }
    if apart > 0.0
        && apart + tol.confusion() <= weld_limit(tol)
        && refusal.to_string().contains("did not close into a shell")
    {
        struct Plain;
        impl Drop for Plain {
            fn drop(&mut self) {
                SEW_ACROSS.set(0.0);
            }
        }
        SEW_ACROSS.set(apart + tol.confusion());
        let _plain = Plain;
        match attempt(model, tol) {
            Ok(built) => return Ok(built),
            Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
                return Err(e);
            }
            Err(_) => {}
        }
    }
    if refusal.to_string().contains("no coincident partner face") {
        struct Plain;
        impl Drop for Plain {
            fn drop(&mut self) {
                AS_STATED.set(false);
            }
        }
        AS_STATED.set(true);
        let _plain = Plain;
        match attempt(model, tol) {
            Ok(built) => return Ok(built),
            Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
                return Err(e);
            }
            Err(_) => {}
        }
    }
    Err(refusal)
}

/// The retries [`settled`] runs after the first attempt refused with
/// `refusal`, short of the last, and the first result one of them keeps.
fn resettled(
    model: &mut Model,
    tol: Tolerances,
    refusal: &ogeom_core::OgeomError,
    attempt: &mut impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    if let Some(built) = settled_after(model, tol, refusal, attempt)? {
        return Ok(Some(built));
    }
    if !RIM_CROSSED.get() {
        return Ok(None);
    }
    struct Rims;
    impl Drop for Rims {
        fn drop(&mut self) {
            AT_RIMS.set(false);
        }
    }
    AT_RIMS.set(true);
    let _rims = Rims;
    let limit = weld_limit(tol);
    let at_rims = match attempt(model, tol) {
        Ok(built) if welded_within(model, &built.shape, limit)? => return Ok(Some(built)),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Ok(_) => return Ok(None),
        Err(e) => e,
    };
    if let Some(built) = settled_after(model, tol, &at_rims, attempt)?
        && welded_within(model, &built.shape, limit)?
    {
        return Ok(Some(built));
    }
    Ok(None)
}

/// The retries [`settled`] runs after an attempt refused with `refusal`,
/// pieces on the other solid's boundary settled from both sides of them,
/// and the first result one of them keeps.
fn settled_after(
    model: &mut Model,
    tol: Tolerances,
    refusal: &ogeom_core::OgeomError,
    attempt: &mut impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            SETTLE_FROM_BOTH_SIDES.set(false);
        }
    }
    SETTLE_FROM_BOTH_SIDES.set(true);
    let _reset = Reset;
    if let Some(built) = retried(model, tol, refusal, attempt)? {
        return Ok(Some(built));
    }
    if refusal.to_string().contains("no coincident partner face") {
        struct Plain;
        impl Drop for Plain {
            fn drop(&mut self) {
                FROM_CROSSINGS.set(Crossings::Off);
            }
        }
        let _plain = Plain;
        FROM_CROSSINGS.set(Crossings::Lines);
        if let Some(built) = retried(model, tol, refusal, attempt)? {
            return Ok(Some(built));
        }
        FROM_CROSSINGS.set(Crossings::Curves);
        if let Some(built) = retried(model, tol, refusal, attempt)?
            && welded_within(model, &built.shape, weld_limit(tol))?
        {
            return Ok(Some(built));
        }
    }
    Ok(None)
}

/// The retries [`settled`] runs after a first attempt refused with
/// `refusal`, each in turn, and the first result one of them keeps.
fn retried(
    model: &mut Model,
    tol: Tolerances,
    refusal: &ogeom_core::OgeomError,
    attempt: &mut impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    if refusal.to_string().contains("no coincident partner face") {
        match attempt(model, tol) {
            Ok(built) => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            Err(_) => {}
        }
    }
    if refusal.to_string().contains("did not close into a shell") {
        struct Clear;
        impl Drop for Clear {
            fn drop(&mut self) {
                PROBE_AT_CENTRE.set(false);
            }
        }
        PROBE_AT_CENTRE.set(true);
        let _clear = Clear;
        match attempt(model, tol) {
            Ok(built) => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            Err(_) => {}
        }
    }
    let limit = weld_limit(tol);
    struct Restore(f64);
    impl Drop for Restore {
        fn drop(&mut self) {
            MEMBRANE_WITHIN.set(self.0);
        }
    }
    let restore = Restore(MEMBRANE_WITHIN.replace(limit.min(MEMBRANE_WITHIN.get())));
    for widen in [10.0, 100.0] {
        let Ok(wider) = Tolerances::with_scale(tol.scale() / widen) else {
            continue;
        };
        match attempt(model, wider) {
            Ok(built) if welded_within(model, &built.shape, limit)? => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            _ => {}
        }
    }
    drop(restore);
    {
        struct Off;
        impl Drop for Off {
            fn drop(&mut self) {
                AT_STATED.set(false);
            }
        }
        AT_STATED.set(true);
        let _off = Off;
        #[cfg(feature = "test-hooks")]
        work::add(|w| w.stated_retries += 1);
        match attempt(model, tol) {
            Ok(built) => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            Err(_) => {}
        }
    }
    Ok(None)
}

/// A pole's pcurve over `range` as the straight line it runs along in the
/// chart, where it does (a line trimmed to the pole's stretch, as a
/// boolean writes a pole edge back): the poles are cut where strands end
/// on them by reading their line, and a pole read as some other curve was
/// never cut. Any other curve is kept as it is.
fn straight_pole(
    planar: &PlanarCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<(PlanarCurve, (f64, f64))> {
    use ogeom_geom::Curve2d as _;
    if matches!(planar, PlanarCurve::Line(_)) {
        return Ok((planar.clone(), range));
    }
    let (a, b) = (
        planar.point_at(range.0, tol)?,
        planar.point_at(range.1, tol)?,
    );
    let length = a.distance(b);
    if length <= tol.parametric() {
        return Ok((planar.clone(), range));
    }
    for k in 1..8 {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 8.0;
        let (run, off) = (b - a, planar.point_at(t, tol)? - a);
        if (run.x * off.y - run.y * off.x).abs() / length > tol.parametric() {
            return Ok((planar.clone(), range));
        }
    }
    let line = ogeom_geom::Line2d::segment(a, b, tol)?;
    let domain = line.domain();
    Ok((line.into(), domain))
}

/// The widest tolerance a retried result may state: the weld distance, and
/// the wider confusion the welding run itself states its vertices to.
fn weld_limit(tol: Tolerances) -> f64 {
    (tol.confusion() * 1e2 + tol.confusion() * 10.0) * (1.0 + 1e-9)
}

/// Whether every edge and vertex of `shape` states a tolerance within
/// `weld`.
fn welded_within(model: &Model, shape: &Shape, weld: f64) -> OgeomResult<bool> {
    for kind in [ShapeType::Edge, ShapeType::Vertex] {
        for sub in explore_unique(model, shape, kind)? {
            if model.tolerance_of(&sub)?.is_some_and(|t| t.get() > weld) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The strict audit's verdict: the filtered fill and the unfiltered fill
/// keep the same material.
///
/// Pieces are compared as merged parameter intervals per section: the
/// unfiltered run sees more curves, and a dropped pair's dead section can
/// still cross a live one and split its pieces differently, so individual
/// piece boundaries are noise and the kept *union* is the signal. Sections
/// are keyed by their face pair and a sampled point, because indices differ
/// between the runs.
fn audit_fill_equivalence(
    filtered: (&[SectionRec], &[SectionPiece]),
    unfiltered: (&[SectionRec], &[SectionPiece]),
    tol: Tolerances,
) {
    type Key = (usize, usize, [i64; 3]);
    let key_of = |section: &SectionRec, tol: Tolerances| -> Key {
        let (lo, hi) = section.curve.domain();
        let p = section
            .curve
            .point_at(f64::midpoint(lo, hi), tol)
            .unwrap_or(ogeom_math::Point::ORIGIN);
        let grid = (tol.confusion() * 100.0).max(1e-9);
        #[allow(clippy::cast_possible_truncation)]
        let q = |x: f64| (x / grid).round() as i64;
        (section.face_a, section.face_b, [q(p.x), q(p.y), q(p.z)])
    };
    let merged = |sections: &[SectionRec], pieces: &[SectionPiece], tol: Tolerances| {
        let mut kept: ogeom_core::FastMap<Key, Vec<(f64, f64)>> = ogeom_core::FastMap::default();
        for piece in pieces {
            kept.entry(key_of(&sections[piece.section], tol))
                .or_default()
                .push(piece.range);
        }
        let slop = tol.parametric() * 10.0;
        for ranges in kept.values_mut() {
            ranges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
            let mut out: Vec<(f64, f64)> = Vec::with_capacity(ranges.len());
            for &(lo, hi) in ranges.iter() {
                match out.last_mut() {
                    Some(last) if lo <= last.1 + slop => last.1 = last.1.max(hi),
                    _ => out.push((lo, hi)),
                }
            }
            *ranges = out;
        }
        kept
    };
    let a = merged(filtered.0, filtered.1, tol);
    let b = merged(unfiltered.0, unfiltered.1, tol);
    let slop = tol.parametric() * 20.0;
    let matches = |x: &Vec<(f64, f64)>, y: &Vec<(f64, f64)>| {
        x.len() == y.len()
            && x.iter()
                .zip(y)
                .all(|(p, q)| (p.0 - q.0).abs() <= slop && (p.1 - q.1).abs() <= slop)
    };
    for (key, ranges) in &b {
        let held = a.get(key);
        assert!(
            held.is_some_and(|r| matches(r, ranges)),
            "strict bound-filter audit: the unfiltered fill keeps {ranges:?} \
             on section {key:?}, the filtered fill keeps {held:?}; the \
             filter drops material"
        );
    }
    for (key, ranges) in &a {
        assert!(
            b.contains_key(key),
            "strict bound-filter audit: the filtered fill keeps {ranges:?} \
             on section {key:?} the unfiltered fill never made"
        );
    }
}

fn general_fuse(
    model: &Model,
    a: &Shape,
    b: &Shape,
    aside: Option<&[aside::Aside; 2]>,
    tol: Tolerances,
) -> OgeomResult<GeneralFused> {
    general_fuse_as(model, a, b, Operands::Solids, aside, tol)
}

/// [`general_fuse`] for an operation that keeps pieces by their state,
/// refusing a classification that cannot hold: the boundaries meet along
/// sections, every piece of both reads outside the other, and yet the
/// solids overlap.
///
/// Two solids sharing volume share it within a region bounded by pieces of
/// their boundaries, and each such piece lies inside the other solid or on
/// its boundary. So where no piece is in or on, the solids are apart, and a
/// section that leaves them overlapping was not arranged into the faces it
/// splits (a loop through a point where the surfaces touch, read as one
/// ring round both its lobes). Kept as read, the cut would keep the whole
/// solid and the union both solids whole. The overlap is asked of sampled
/// boundaries, only where every piece reads out and sections were found.
fn general_fuse_classified(
    model: &Model,
    a: &Shape,
    b: &Shape,
    aside: Option<&[aside::Aside; 2]>,
    tol: Tolerances,
) -> OgeomResult<GeneralFused> {
    let fused = general_fuse(model, a, b, aside, tol)?;
    #[cfg(feature = "test-hooks")]
    work::add(|w| w.pieces += fused.pieces.len());
    if !fused.sections.is_empty()
        && fused.pieces.iter().all(|p| p.state == PieceState::Out)
        && !apart(model, a, b, tol)?
    {
        ogeom_bail!(
            NotDone,
            "the boundaries cross yet every piece of both reads outside the other, \
             though the solids overlap; a section was not arranged into the faces it \
             splits"
        );
    }
    Ok(fused)
}

/// For each gathered face, the faces of the same solid it shares a vertex
/// with, itself among them, by their places in `solid.faces`.
fn neighbouring_faces(
    model: &Model,
    solid: &GSolid,
) -> OgeomResult<Vec<ogeom_core::FastSet<usize>>> {
    // A vertex is its node where it is placed, as an edge is (see
    // [`EdgeKey`]): one node placed twice is two corners in space.
    let mut by_vertex: ogeom_core::FastMap<EdgeKey, Vec<usize>> = ogeom_core::FastMap::default();
    for (index, face) in solid.faces.iter().enumerate() {
        for vertex in explore_unique(model, &face.face, ShapeType::Vertex)? {
            by_vertex
                .entry(EdgeKey::of(&vertex))
                .or_default()
                .push(index);
        }
    }
    let mut neighbours: Vec<ogeom_core::FastSet<usize>> = (0..solid.faces.len())
        .map(|index| std::iter::once(index).collect())
        .collect();
    for faces in by_vertex.values() {
        for &i in faces {
            neighbours[i].extend(faces.iter().copied());
        }
    }
    Ok(neighbours)
}

/// The general fuse of `a` and `b` taken as `operands` says. Where `a` is a
/// sheet only its pieces are made, and the pieces of `b` are not. The
/// faces `aside` names are left out.
fn general_fuse_as(
    model: &Model,
    a: &Shape,
    b: &Shape,
    operands: Operands,
    aside: Option<&[aside::Aside; 2]>,
    tol: Tolerances,
) -> OgeomResult<GeneralFused> {
    // Read here, on the caller's thread: the pieces are classified on
    // worker threads, which do not share it.
    let settle = SETTLE_FROM_BOTH_SIDES.get();
    let read_at_centre = PROBE_AT_CENTRE.get();
    let read_thin = AT_STATED.get();
    ogeom_core::progress::stage("boolean: gather");
    let mut ga = gather(
        model,
        a,
        matches!(operands, Operands::SheetBySolid | Operands::SheetBySheet),
        aside.map(|aside| &aside[0]),
        tol,
    )?;
    let mut gb = gather(
        model,
        b,
        operands == Operands::SheetBySheet,
        aside.map(|aside| &aside[1]),
        tol,
    )?;
    let (onto_a, onto_b) = (
        drums_onto_planes(&ga, &gb, tol),
        drums_onto_planes(&gb, &ga, tol),
    );
    move_drums(&mut ga, &onto_a, tol)?;
    move_drums(&mut gb, &onto_b, tol)?;
    #[cfg(feature = "test-hooks")]
    work::add(|w| w.faces_gathered += ga.faces.len() + gb.faces.len());
    let neighbours = if operands == Operands::Itself {
        Some(neighbouring_faces(model, &ga)?)
    } else {
        None
    };
    ogeom_core::progress::stage("boolean: intersect");
    let (
        sections,
        section_pieces,
        contacts,
        tangents,
        contact_along,
        mut paves,
        same_a,
        same_b,
        hugs,
    ) = fill(&ga, &gb, false, neighbours.as_deref(), tol)?;
    if operands == Operands::SheetBySheet && same_a.iter().any(|partners| !partners.is_empty()) {
        ogeom_bail!(
            NotDone,
            "the sheets lie on one surface where they overlap; a sheet split \
             by another lying along it is not resolved"
        );
    }
    for section in &sections {
        for (pcurve, face) in [
            (&section.pc_a, &ga.faces[section.face_a]),
            (&section.pc_b, &gb.faces[section.face_b]),
        ] {
            past_a_pole(section, pcurve, face, tol)?;
        }
    }
    let along = crossed_along(&contacts);
    let mut junctions = pave_junctions(&ga, &gb, &paves, &along, tol)?;
    if diag::enabled(Channel::Wire) {
        for j in &junctions {
            ogeom_core::diag!(
                Channel::Wire,
                "JUNCTION from paves at {:?} reach {:.3e}",
                j.at,
                j.reach
            );
        }
        for j in &hugs {
            ogeom_core::diag!(
                Channel::Wire,
                "JUNCTION from hugs at {:?} reach {:.3e}",
                j.at,
                j.reach
            );
        }
    }
    junctions.extend(hugs);
    // A vertex of an edge shared with a face set aside is that vertex in
    // the rebuild, which every strand ending within its span names.
    let mut kept_bins = bins::Bins::new(tol.confusion() * 1e2);
    let kept_points: Vec<Point> = ga
        .shared_vertices
        .iter()
        .chain(&gb.shared_vertices)
        .map(|(at, _)| *at)
        .collect();
    for (index, at) in kept_points.iter().enumerate() {
        kept_bins.insert(*at, index);
    }
    let mut near_kept = Vec::new();
    for face in ga.faces.iter().chain(gb.faces.iter()) {
        // An input vertex that owns a span (the corner an earlier boolean
        // welded three rims into, each ending a fraction of a micron from
        // the others) is a junction here as well: the rebuild's positional
        // weld reaches only a few honesties, and minted separately those
        // ends stand as three vertices with hairlines between them.
        for e in &face.edges {
            for (at, radius) in e.ends {
                let mut kept = || {
                    !kept_points.is_empty()
                        && if kept_bins.near_into(at, radius, &mut near_kept) {
                            near_kept
                                .iter()
                                .any(|&i| kept_points[i].distance(at) <= radius)
                        } else {
                            kept_points.iter().any(|p| p.distance(at) <= radius)
                        }
                };
                if radius > tol.confusion() * 10.0 && !kept() {
                    ogeom_core::diag!(
                        Channel::Wire,
                        "JUNCTION from vertex at {at:?} radius {radius:.3e}"
                    );
                    junctions.push(Junction {
                        at,
                        reach: radius,
                        onto_vertex: false,
                    });
                }
            }
        }
    }
    // Junctions whose balls overlap are one junction: two pave clusters a
    // fraction of a micron apart at a corner where rims meet, each owning
    // a span that reaches into the other's, would weld a strand end to
    // whichever covers it first, and the corner would stand as two
    // vertices with a hairline between them, kept on one face, dropped
    // on another. Merged, the corner is one vertex whose reach covers
    // every member's span, on every face alike.
    junctions = merge_junctions(junctions);
    for (from_a, gathered) in [(true, &ga), (false, &gb)] {
        seams_through_junctions(
            gathered,
            from_a,
            (&sections, &section_pieces),
            &junctions,
            &mut paves,
            tol,
        )?;
    }
    if diag::enabled(Channel::Wire) {
        for j in &junctions {
            ogeom_core::diag!(
                Channel::Wire,
                "JUNCTION at {:?} reach {:.3e}",
                j.at,
                j.reach
            );
        }
        let kind = |sf: &SurfaceGeometry| match sf {
            SurfaceGeometry::Plane(_) => "plane",
            SurfaceGeometry::Cylinder(_) => "cyl",
            SurfaceGeometry::BSpline(_) => "bspline",
            _ => "other",
        };
        for (si, sec) in sections.iter().enumerate() {
            let fa = &ga.faces[sec.face_a];
            let fb = &gb.faces[sec.face_b];
            let d = sec.curve.domain();
            let mut worst_a: f64 = 0.0;
            let mut worst_b: f64 = 0.0;
            let mut len = 0.0;
            let mut prev: Option<Point> = None;
            for i in 0..=32 {
                let t = if i == 32 {
                    d.1
                } else {
                    d.0 + (d.1 - d.0) * f64::from(i) / 32.0
                };
                let at = sec.curve.point_at(t, tol)?;
                if let Some(p) = prev {
                    len += p.distance(at);
                }
                prev = Some(at);
                let near = |f: &GFace| -> OgeomResult<f64> {
                    let mut best = f64::INFINITY;
                    for e in &f.edges {
                        best = best.min(distance_to_edge_curve(&e.curve, e.crange, at, tol)?);
                    }
                    Ok(best)
                };
                worst_a = worst_a.max(near(fa)?);
                worst_b = worst_b.max(near(fb)?);
            }
            let pieces = section_pieces.iter().filter(|p| p.section == si).count();
            ogeom_core::diag!(
                Channel::Wire,
                "SECTION s{si} a{}({}) x b{}({}) tol {:.2e} closed {} len {:.4} pieces {pieces} hug-dist a {:.3e} b {:.3e}",
                sec.face_a,
                kind(&fa.surface),
                sec.face_b,
                kind(&fb.surface),
                sec.tolerance,
                sec.closed,
                len,
                worst_a,
                worst_b
            );
        }
        ogeom_core::diag!(
            Channel::Wire,
            "TANGENTS {}  CONTACTS {}",
            tangents.len(),
            contacts.len()
        );
    }
    // The strict audit: fill again with every pair admitted and demand the
    // same kept material. The production result above is always the
    // filtered run (under audit too), so the audit compares rather than
    // substitutes, and a filtered run *is* replayed exactly. Zero cost with
    // the variable unset.
    if diag::enabled(Channel::BoolAudit) {
        // A pair only the unfiltered fill tries may refuse (a march beside a
        // pole the filter never sends it), which says nothing of the
        // filter. The audit says so and compares nothing.
        match fill(&ga, &gb, true, neighbours.as_deref(), tol) {
            Ok((audit_sections, audit_pieces, ..)) => audit_fill_equivalence(
                (&sections, &section_pieces),
                (&audit_sections, &audit_pieces),
                tol,
            ),
            Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
            Err(e) => ogeom_core::diag!(
                Channel::BoolAudit,
                "bound-filter audit skipped: the unfiltered fill refused: {e}"
            ),
        }
    }

    ogeom_core::progress::stage("boolean: split");
    let mut pieces: Vec<FacePiece> = Vec::new();
    // Every face's strands, before any face is arranged: the dust decision
    // (a piece shorter than its face's snap collapses to one node) must
    // be one decision per shared edge piece, and a sub-piece of one edge
    // lies in several charts, each with its own snap and its own metric.
    // A band's rail hugging a wedge's cylinder for a quarter of a micron
    // is dust in the wall's chart and a strand in the cylinder's, and the
    // sew then finds the cylinder's piece used once.
    let chart_length =
        |line: &[Point2]| -> f64 { line.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>() };
    // Each edge's end vertices where they stand, each with its radius and
    // the edge's own: a contact strand on the other solid's face ends at
    // one of them.
    let mut edge_vertices: ogeom_core::FastMap<EdgeKey, [(Point, f64); 2]> =
        ogeom_core::FastMap::default();
    for e in ga
        .faces
        .iter()
        .chain(gb.faces.iter())
        .flat_map(|f| f.edges.iter())
    {
        edge_vertices.entry(e.node).or_insert([
            (e.vertices[0], e.ends[0].1 + e.tolerance),
            (e.vertices[1], e.ends[1].1 + e.tolerance),
        ]);
    }
    let strands_of = |from_a: bool,
                      fi: usize,
                      face: &GFace|
     -> OgeomResult<(Vec<Strand<Tag>>, f64, Vec<Junction>)> {
        let mut strands: Vec<Strand<Tag>> = Vec::new();
        // How far apart the paves one junction stands for lie. A section
        // meeting the edge all but tangentially (its fitted end sliding
        // along it) paves it a few microns from where an exact section
        // crosses. The first pave speaks for both, and the other's strand
        // must still reach it.
        // Each such junction where it stands in space and in this chart,
        // and how far its members reach: the strands ending at a member are
        // carried onto it below.
        let mut gathered: Vec<(Point, Point2, f64)> = Vec::new();
        for (ei, e) in face.edges.iter().enumerate() {
            let mut stops = vec![e.crange.0];
            if let Some(ts) = paves.get(&e.node) {
                // Paves the edge itself cannot tell apart are one
                // junction, and the cluster's representative speaks for it.
                let off = |face, p| off_face(&ga, &gb, face, p, tol);
                let along = along.get(&e.node).map_or(&[][..], Vec::as_slice);
                for c in cluster_paves(&e.curve, e.crange, e.tolerance, ts, Some(&off), along, tol)?
                {
                    if c.members > 1 {
                        let chart = e.pcurve.point_at(rescale(c.t, e.crange, e.prange), tol)?;
                        gathered.push((c.at, chart, c.span + c.honesty));
                    }
                    stops.push(c.t);
                }
            }
            stops.push(e.crange.1);
            // A closed boundary edge (a cap's full circle) needs two
            // distinct endpoints per strand.
            let closed_edge = e
                .curve
                .point_at(e.crange.0, tol)?
                .distance(e.curve.point_at(e.crange.1, tol)?)
                <= tol.confusion();
            if closed_edge && stops.len() == 2 {
                stops.insert(1, f64::midpoint(e.crange.0, e.crange.1));
            }
            for pair in stops.windows(2) {
                let sub = (pair[0], pair[1]);
                if sub.1 - sub.0 <= tol.parametric() {
                    continue;
                }
                strands.push(Strand {
                    polyline: pcurve_polyline(
                        &e.pcurve,
                        e.prange,
                        e.crange,
                        sub,
                        &face.surface,
                        tol,
                    )?,
                    tag: Tag::Boundary {
                        edge: ei,
                        range: sub,
                    },
                    boundary: true,
                    sag: f64::INFINITY,
                });
                if let Some((other_pc, orange)) = &e.other_side {
                    strands.push(Strand {
                        polyline: pcurve_polyline(
                            other_pc,
                            *orange,
                            e.crange,
                            sub,
                            &face.surface,
                            tol,
                        )?,
                        tag: Tag::Boundary {
                            edge: ei,
                            range: sub,
                        },
                        boundary: true,
                        sag: f64::INFINITY,
                    });
                }
            }
        }
        // The trim the strands are folded inside: the boundary strands as
        // they stand and the poles whole, with their boxes, drawn once for
        // the face. The poles become strands only after the sections that
        // cut them, and a trim without them is open where it reaches one: a
        // point inside reads outside along a ray leaving through the gap.
        let pole_lines = face
            .poles
            .iter()
            .map(|pole| {
                pcurve_polyline(
                    &pole.pcurve,
                    pole.prange,
                    pole.prange,
                    pole.prange,
                    &face.surface,
                    tol,
                )
            })
            .collect::<OgeomResult<Vec<_>>>()?;
        let mut section_trim: Option<arrange::Trim> = None;
        for sp in &section_pieces {
            let section = &sections[sp.section];
            let (belongs, pcurve) = if from_a {
                (section.face_a == fi && !sp.hugs[0], &section.pc_a)
            } else {
                (section.face_b == fi && !sp.hugs[1], &section.pc_b)
            };
            if !belongs {
                continue;
            }
            let domain = section.curve.domain();
            let sub = if section.closed {
                (
                    fold(sp.range.0, domain),
                    sp.range.1 - sp.range.0 + fold(sp.range.0, domain),
                )
            } else {
                sp.range
            };
            // The pcurve shares the curve's parameterization. Sampling
            // uses folded parameters for periodic curves.
            let fold_at = |t: f64| if section.closed { fold(t, domain) } else { t };
            // A section's sampling on a plane holds its chords to a sag,
            // which the arrangement can read it to.
            let (count, sag) = if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
                (32, f64::INFINITY)
            } else if matches!(pcurve, PlanarCurve::Line(_)) {
                (32, 0.0)
            } else {
                match circle_sag_count(pcurve, sub.0, sub.1, 32) {
                    Some(count) => (
                        count,
                        circle_sag(pcurve, sub.0, sub.1, count).unwrap_or(f64::INFINITY),
                    ),
                    None => sag_count(&|t| pcurve.point_at(fold_at(t), tol), sub.0, sub.1, 32)?,
                }
            };
            let mut line = Vec::with_capacity(count + 1);
            for i in 0..=count {
                #[allow(clippy::cast_precision_loss)]
                let t = sub.0 + (sub.1 - sub.0) * i as f64 / count as f64;
                let tf = if section.closed { fold(t, domain) } else { t };
                line.push(pcurve.point_at(tf, tol)?);
            }
            // Folding the parameter can tear the sampled polyline at the
            // period. Unwrap it pointwise, then bring the whole strand
            // into the chart with one shift.
            unwrap_polyline(&mut line, &face.surface, tol);
            let trim = section_trim.get_or_insert_with(|| {
                arrange::Trim::new(
                    strands
                        .iter()
                        .filter(|st| st.boundary)
                        .map(|st| st.polyline.clone())
                        .chain(pole_lines.iter().cloned())
                        .collect(),
                )
            });
            fold_line_inside(&mut line, &face.surface, trim, tol);
            strands.push(Strand {
                polyline: line,
                tag: Tag::Section {
                    section: sp.section,
                    range: sp.range,
                },
                boundary: false,
                sag,
            });
        }
        // A section's end the paving read at a pole (its crossing within
        // the window that reads a crossing as an edge's end vertex, at the
        // pole vertex of a seam) stands beside the pole in space by at most
        // that window, but in the chart that can be anywhere along the
        // pole's line, the chart's free parameter shrinking to nothing
        // there: off the pole row by the window over the radius, and along
        // it by any amount, past the line's end where it crosses the seam
        // beside the pole. It is the pole, and is carried onto the pole's
        // line at its foot, held to the line (every point of which is the
        // pole), where the pole is cut for it below.
        for strand in &mut strands {
            let Tag::Section { section, .. } = strand.tag else {
                continue;
            };
            let s = &sections[section];
            let reach = tol.confusion().max(s.tolerance * 2.0);
            let window = (reach + honest(s.tolerance, tol)).max(tol.confusion() * 1e2);
            let last = strand.polyline.len() - 1;
            for end in [0, last] {
                let at = strand.polyline[end];
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                for pole in &face.poles {
                    let PlanarCurve::Line(line) = &pole.pcurve else {
                        continue;
                    };
                    if q.distance(pole.point) > window {
                        continue;
                    }
                    let axis = line.axis();
                    let along = (at - axis.location)
                        .dot(axis.direction.vector())
                        .clamp(pole.prange.0, pole.prange.1);
                    strand.polyline[end] = axis.point_at(along);
                }
            }
        }
        // A boundary edge ending at a pole (a patch's side running into its
        // collapsed corner) ends there in space, but its pcurve may end off
        // the pole's line by the edge's tolerance over the chart's shrinking
        // scale, and anywhere along it: it is carried onto the line at its
        // foot, as a section's end is above.
        for strand in &mut strands {
            let Tag::Boundary { edge, .. } = strand.tag else {
                continue;
            };
            let window = (face.edges[edge].tolerance * 2.0).max(tol.confusion() * 1e2);
            let last = strand.polyline.len() - 1;
            for end in [0, last] {
                let at = strand.polyline[end];
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                for pole in &face.poles {
                    let PlanarCurve::Line(line) = &pole.pcurve else {
                        continue;
                    };
                    if q.distance(pole.point) > window {
                        continue;
                    }
                    let axis = line.axis();
                    let along = (at - axis.location)
                        .dot(axis.direction.vector())
                        .clamp(pole.prange.0, pole.prange.1);
                    strand.polyline[end] = axis.point_at(along);
                }
            }
        }
        // The poles, after the sections, because a section can end *on*
        // a pole (a plane through a ball's axis cuts it exactly there),
        // and the pole has to be cut where that happens or the two meet
        // at no shared node and the arrangement sees a dangling section.
        // A boundary edge ending on it cuts it too. Where any strand ends
        // on the pole, a piece of it with an end no strand meets lies past
        // the face's own boundary (the face spans only part of the pole's
        // line) and would only dangle: it is left out.
        for (pi, pole) in face.poles.iter().enumerate() {
            let mut stops = vec![pole.prange.0, pole.prange.1];
            let mut met: Vec<f64> = Vec::new();
            if let PlanarCurve::Line(line) = &pole.pcurve {
                let axis = line.axis();
                for strand in &strands {
                    if matches!(strand.tag, Tag::Pole { .. }) {
                        continue;
                    }
                    for end in [strand.polyline.first(), strand.polyline.last()]
                        .into_iter()
                        .flatten()
                    {
                        let along = (*end - axis.location).dot(axis.direction.vector());
                        let foot = axis.point_at(along);
                        if foot.distance(*end) > PARAM_SNAP {
                            continue;
                        }
                        met.push(along);
                        if along > pole.prange.0 + PARAM_SNAP && along < pole.prange.1 - PARAM_SNAP
                        {
                            stops.push(along);
                        }
                    }
                }
            }
            stops.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            stops.dedup_by(|a, b| (*a - *b).abs() <= PARAM_SNAP);
            let meets = |t: f64| met.iter().any(|m| (m - t).abs() <= PARAM_SNAP);
            for pair in stops.windows(2) {
                if !met.is_empty() && !(meets(pair[0]) && meets(pair[1])) {
                    continue;
                }
                let sub = (pair[0], pair[1]);
                strands.push(Strand {
                    polyline: pcurve_polyline(
                        &pole.pcurve,
                        pole.prange,
                        pole.prange,
                        sub,
                        &face.surface,
                        tol,
                    )?,
                    tag: Tag::Pole {
                        pole: pi,
                        range: sub,
                    },
                    boundary: true,
                    sag: f64::INFINITY,
                });
            }
        }
        let mut contact_trim: Option<arrange::Trim> = None;
        for (ci, contact) in contacts.iter().enumerate() {
            if contact.target_from_a != from_a || contact.target_face != fi {
                continue;
            }
            // The owner's paves split its edge. The contact strands split
            // at the same parameters, so the sub-edges rebuilt from both
            // sides are the same edges and sew shared.
            let mut stops = vec![contact.crange.0];
            if let Some(ts) = paves.get(&contact.node) {
                stops.extend(
                    cluster_paves(
                        &contact.curve,
                        contact.crange,
                        contact.tolerance,
                        ts,
                        None,
                        &[],
                        tol,
                    )?
                    .iter()
                    .map(|c| c.t),
                );
            }
            stops.push(contact.crange.1);
            let closed_contact = contact
                .curve
                .point_at(contact.crange.0, tol)?
                .distance(contact.curve.point_at(contact.crange.1, tol)?)
                <= tol.confusion();
            if closed_contact && stops.len() == 2 {
                stops.insert(1, f64::midpoint(contact.crange.0, contact.crange.1));
            }
            for pair in stops.windows(2) {
                let sub = (pair[0], pair[1]);
                if sub.1 - sub.0 <= tol.parametric() {
                    continue;
                }
                let mid_t = f64::midpoint(sub.0, sub.1);
                if contact_along[ci]
                    .iter()
                    .any(|(lo, hi)| mid_t >= *lo && mid_t <= *hi)
                {
                    // Already boundary on both sides.
                    ogeom_core::diag!(
                        Channel::Strands,
                        "CONTACT c{ci} {sub:?} on fi={fi}: along boundary, skipped"
                    );
                    continue;
                }
                // Keep only what lies inside this face's trim. The rest
                // of the owner's boundary splits nothing here.
                let mut line = pcurve_polyline(
                    &contact.pcurve,
                    contact.prange,
                    contact.crange,
                    sub,
                    &face.surface,
                    tol,
                )?;
                unwrap_polyline(&mut line, &face.surface, tol);
                let trim = contact_trim.get_or_insert_with(|| {
                    arrange::Trim::new(
                        strands
                            .iter()
                            .filter(|st| st.boundary)
                            .map(|st| st.polyline.clone())
                            .collect(),
                    )
                });
                fold_line_inside(&mut line, &face.surface, trim, tol);
                let mid = interior_of(&line);
                // On the last attempt, where the stretch stands nearest an
                // edge's inside, it is inside the face just where it stands
                // on the side the face keeps its material on: the edge's
                // own curve says so where the outline's chords sag past it.
                let beside = if read_thin {
                    inside_beside_edge(face, contact.curve.point_at(mid_t, tol)?, tol)?
                } else {
                    None
                };
                if !beside.unwrap_or_else(|| trim.inside_slanted(mid)) {
                    ogeom_core::diag!(
                        Channel::Strands,
                        "CONTACT c{ci} {sub:?} on fi={fi}: outside the trim at {mid:?}, skipped"
                    );
                    continue;
                }
                strands.push(Strand {
                    polyline: line,
                    tag: Tag::Contact {
                        contact: ci,
                        range: sub,
                    },
                    boundary: false,
                    sag: f64::INFINITY,
                });
            }
        }

        // A tolerant contact's chart image meets the boundary it paved
        // only as closely as its own slop allows. The arrangement's node
        // weld reaches that far on this face, or the strand dangles a
        // few microns from the junction it belongs to.
        // A fitted section that hugs one of this face's edges leaves it
        // at the hug's far end by up to the hug's own width. Its strand
        // must still find the edge's node there.
        // A strand ending at another member of a junction than the one that
        // speaks for it ends where its curve all but runs along the edge (a
        // section crossing it at a shallow angle, whose crossing is only
        // known along the edge to within the run). An exact section moved
        // along to the junction stays within tolerance of both curves, and
        // carried there it meets the edge's split exactly, so the face's
        // weld need not widen to the cluster's span, which would swallow a
        // sliver beside the edge whole. A fitted one wanders by more than
        // that, and the weld widens for it.
        // On the last attempt that weld is each such end's alone: the face's
        // node weld also decides which regions are too thin to hold a point,
        // and a junction's run along one edge says nothing of a sliver
        // elsewhere on the face.
        let mut spread = 0.0_f64;
        let mut spread_at: Vec<[f64; 2]> = vec![[0.0; 2]; strands.len()];
        for (si, strand) in strands
            .iter_mut()
            .enumerate()
            .filter(|(_, st)| !st.boundary)
        {
            let exact = matches!(strand.tag, Tag::Section { section, .. }
                if sections.get(section).is_some_and(|s| s.tolerance <= 0.0));
            let last = strand.polyline.len().saturating_sub(1);
            for (side, end) in [0, last].into_iter().enumerate() {
                let Some(&at) = strand.polyline.get(end) else {
                    continue;
                };
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                let landed = gathered
                    .iter()
                    .filter(|(junction, _, reach)| junction.distance(q) <= *reach)
                    .min_by(|x, y| x.0.distance(q).total_cmp(&y.0.distance(q)));
                match landed {
                    // A junction on a seam stands on both sides of it, and
                    // the strand is carried to the side it arrives from.
                    Some(&(_, chart, _)) if exact => {
                        strand.polyline[end] = period_shifts(&face.surface, tol)
                            .into_iter()
                            .map(|(du, dv)| Point2::new(chart.x + du, chart.y + dv))
                            .min_by(|p, r| p.distance(at).total_cmp(&r.distance(at)))
                            .unwrap_or(chart);
                    }
                    Some(&(_, _, reach)) => {
                        spread = spread.max(reach);
                        spread_at[si][side] = reach;
                    }
                    None => {}
                }
            }
        }
        // What the face holds everywhere in its chart: its own tolerance
        // and its holes set aside, and twice the run's confusion, which two
        // ends each within it of one point keep between them (a fuzzy
        // run's confusion is its fuzz). What its edges and vertices add
        // holds only near them, and each strand end is welded at the doubt
        // where it stands, below.
        let everywhere = |face: &GFace| -> f64 {
            face.holes_aside
                .iter()
                .fold(face.tolerance, |acc, hole| acc.max(hole.doubt))
                .min(doubt_of(face))
        };
        let near = contacts
            .iter()
            .filter(|c| c.target_from_a == from_a && c.target_face == fi)
            .fold(PARAM_SNAP, |acc, c| acc.max(c.tolerance * 2.0))
            .max(tol.confusion() * 2.0)
            .max(everywhere(face))
            .max(if read_thin { 0.0 } else { spread })
            .max(
                if !read_thin
                    && sections.iter().any(|s| {
                        s.tolerance > 0.0
                            && if from_a {
                                s.face_a == fi
                            } else {
                                s.face_b == fi
                            }
                    })
                {
                    tol.confusion() * 1e3
                } else {
                    0.0
                },
            );
        // The face a section shares with this one.
        let other = |section: usize| -> &GFace {
            let s = &sections[section];
            if from_a {
                &gb.faces[s.face_b]
            } else {
                &ga.faces[s.face_a]
            }
        };
        let face_snap = near.max(
            (0..sections.len())
                .filter(|&k| {
                    if from_a {
                        sections[k].face_a == fi
                    } else {
                        sections[k].face_b == fi
                    }
                })
                .map(|k| everywhere(other(k)))
                .fold(0.0_f64, f64::max),
        );
        // The weld at a strand's end, in space. An edge's end lands in the
        // chart within its vertex's own recorded doubt plus its own image's:
        // a projected pcurve is honest to the edge's tolerance, and the
        // vertex it ends at was welded to some earlier gap. A section ends
        // where it crosses the other face's boundary, and lands there within
        // that boundary's own doubt: two sections through neighbouring
        // facets of a converted mesh stop on the edge they share a few
        // microns apart, each on its own facet's plane, and meet on this
        // face only if the weld reaches that far. Both doubts are those of
        // the edges and vertices near the end ([`doubt_at`]): a loose edge
        // welds what crosses or ends near it, and a section crossing a
        // tight edge elsewhere on the face keeps its crossing, as it does on
        // the face across that edge.
        let mut ends: Vec<(usize, usize, Point2, Point, f64)> = Vec::new();
        for (si, st) in strands.iter().enumerate() {
            if st.polyline.len() < 2 || matches!(st.tag, Tag::Pole { .. }) {
                continue;
            }
            let last = st.polyline.len() - 1;
            for end in [0, last] {
                let at = st.polyline[end];
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                let mut weld = doubt_at(face, q, tol)?.max(near);
                if let Tag::Section { section, .. } = st.tag {
                    weld = weld.max(doubt_at(other(section), q, tol)?);
                    // An exact section's stop read at an exact edge's end
                    // vertex (see the crossings in [`fill`]) stands beside
                    // the vertex by as far as the section passes it, up to
                    // the weld that read it there, and is that vertex.
                    let s = &sections[section];
                    if s.tolerance <= 0.0 && exact(&s.curve) {
                        let reach = tol.confusion();
                        for e in face
                            .edges
                            .iter()
                            .filter(|e| e.image.is_none() && exact(&e.curve))
                        {
                            let read_at = (reach + e.tolerance.max(e.ends_tolerance))
                                .max(tol.confusion() * 10.0);
                            for end_t in [e.crange.0, e.crange.1] {
                                let gap = e.curve.point_at(end_t, tol)?.distance(q);
                                if gap <= read_at {
                                    weld = weld.max(gap + tol.confusion());
                                }
                            }
                        }
                    }
                }
                // A contact's end is its own edge's end vertex, within that
                // vertex's radius, as an edge of this face ends at its own.
                if let Tag::Contact { contact, .. } = st.tag
                    && let Some(vertices) = edge_vertices.get(&contacts[contact].node)
                {
                    for &(v, radius) in vertices {
                        if q.distance(v) <= radius {
                            weld = weld.max(radius);
                        }
                    }
                }
                // On the last attempt the face weld leaves out a junction's
                // run and a fitted section's hug width, and they hold at the
                // ends they concern: a fitted section's end, hugging an edge
                // within the two curves' radii or crossing it, meets the
                // edge's node within its own reach, twice its tolerance.
                if read_thin {
                    weld = weld.max(spread_at[si][usize::from(end != 0)]);
                    if let Tag::Section { section, .. } = st.tag {
                        weld = weld.max(honest(sections[section].tolerance, tol) * 2.0);
                    }
                }
                ends.push((si, end, at, q, weld));
            }
        }
        // The doubts above are lengths in space, and the weld runs in the
        // chart. Where the chart stretches every direction (a sphere of a
        // few millimetres runs a radian over its whole radius), a gap in
        // space is that many times narrower in the chart, and a weld taken
        // at the space length would swallow whole edges of a small patch.
        // Where some direction shrinks instead (a pole, a thin cylinder's
        // turn), the length stands.
        let stretch = strands
            .iter()
            .flat_map(|st| st.polyline.iter())
            .step_by(8)
            .filter_map(|p| face.surface.d1_at(p.x, p.y, tol).ok())
            .map(|(du, dv)| du.magnitude().min(dv.magnitude()))
            .fold(f64::INFINITY, f64::min);
        let in_chart = |length: f64| -> f64 {
            if stretch.is_finite() && stretch > 1.0 {
                (length / stretch).max(PARAM_SNAP)
            } else {
                length
            }
        };
        let face_snap = in_chart(face_snap);
        let mut reach: Vec<f64> = ends.iter().map(|e| in_chart(e.4)).collect();
        // Where only one direction stretches (a cylinder's turn against its
        // straight length) no single scale holds, and a strand can be short
        // in the chart yet long in space. A strand longer in space than its
        // own doubt is no dust, so no weld near it reaches as far as its
        // chart length, or it collapses and takes its neighbours' ends with
        // it.
        let own_doubt = |si: usize| -> f64 {
            let at_ends = ends
                .iter()
                .filter(|e| e.0 == si)
                .fold(near, |acc, e| acc.max(e.4));
            match strands[si].tag {
                Tag::Section { section, .. } => at_ends.max(everywhere(other(section))),
                Tag::Pole { .. } => f64::INFINITY,
                Tag::Boundary { .. } | Tag::Contact { .. } => at_ends,
            }
        };
        let loosest = reach.iter().copied().fold(face_snap, f64::max);
        let kept: Vec<(Point2, Point2, f64)> = strands
            .iter()
            .enumerate()
            .filter(|(_, st)| st.polyline.len() >= 2 && chart_length(&st.polyline) <= loosest)
            .filter_map(|(si, st)| {
                let points: Vec<Point> = st
                    .polyline
                    .iter()
                    .filter_map(|p| face.surface.point_at(p.x, p.y, tol).ok())
                    .collect();
                let space: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
                (space > own_doubt(si)).then(|| {
                    (
                        st.polyline[0],
                        st.polyline[st.polyline.len() - 1],
                        chart_length(&st.polyline) / 2.0,
                    )
                })
            })
            .collect();
        let face_snap = kept
            .iter()
            .filter(|k| k.2 * 2.0 <= face_snap)
            .fold(face_snap, |acc, k| acc.min(k.2))
            .max(PARAM_SNAP);
        for (e, r) in ends.iter().zip(reach.iter_mut()) {
            for &(a, b, half) in &kept {
                if half * 2.0 <= *r && (e.2.distance(a) <= *r || e.2.distance(b) <= *r) {
                    *r = r.min(half);
                }
            }
            *r = r.max(face_snap);
        }
        // Each end whose own weld reaches past the face's is welded here,
        // onto one node with every end within its reach: a node of the
        // face's boundary where one stands among them, the nearest to a
        // vertex of the face first. A strand whose ends both land on the
        // node and which is no longer than the weld there is dust.
        let mut class: Vec<usize> = (0..ends.len()).collect();
        fn root(class: &mut [usize], mut i: usize) -> usize {
            while class[i] != i {
                class[i] = class[class[i]];
                i = class[i];
            }
            i
        }
        let mut welds: Vec<Junction> = Vec::new();
        let mut welded = false;
        // Two ends that each stand within a vertex's radius (its own and
        // its edge's) of the same vertex are that vertex: the curves a file
        // closes on one vertex can stop either side of it, each within the
        // radius. The vertex is one of this face's edges', or one ending
        // the other solid's edge a contact strand runs along. In the chart they stand as far apart as both their
        // reaches, or further where the chart stretches more than the
        // reaches allow for (a cylinder narrower than the radius turns a
        // quarter within it); there the stretch between them, read as an
        // edge is at seventeen points, lies within the vertex too, on a
        // chart with no pole's line: about a pole the whole chart may lie
        // within one vertex, and the pole's line is what joins the ends
        // there. Two ends at one point in space are no such gap: apart in
        // the chart, they stand on a pole's line, whose ends the pole keeps
        // apart.
        let poled = strands.iter().any(|st| matches!(st.tag, Tag::Pole { .. }));
        let one_vertex = |i: usize, j: usize| -> bool {
            let (a, b) = (ends[i].3, ends[j].3);
            let (from, to) = (ends[i].2, ends[j].2);
            if a.distance(b) <= tol.confusion() {
                return false;
            }
            let close = from.distance(to) <= reach[i] + reach[j];
            let contact_vertices = [i, j]
                .into_iter()
                .filter_map(|k| match strands[ends[k].0].tag {
                    Tag::Contact { contact, .. } => edge_vertices.get(&contacts[contact].node),
                    _ => None,
                })
                .flatten()
                .copied();
            face.edges
                .iter()
                .flat_map(|e| {
                    e.vertices
                        .iter()
                        .zip(e.ends)
                        .map(|(&v, (_, radius))| (v, radius + e.tolerance))
                })
                .chain(contact_vertices)
                .any(|(v, within)| {
                    a.distance(v) <= within
                        && b.distance(v) <= within
                        && (close
                            || !poled
                                && (1..16).all(|k| {
                                    let t = f64::from(k) / 16.0;
                                    let at = Point2::new(
                                        from.x + (to.x - from.x) * t,
                                        from.y + (to.y - from.y) * t,
                                    );
                                    face.surface
                                        .point_at(at.x, at.y, tol)
                                        .is_ok_and(|q| q.distance(v) <= within)
                                }))
                })
        };
        for i in 0..ends.len() {
            if reach[i] <= face_snap {
                continue;
            }
            for j in 0..ends.len() {
                if j != i && (ends[i].2.distance(ends[j].2) <= reach[i] || one_vertex(i, j)) {
                    let (ri, rj) = (root(&mut class, i), root(&mut class, j));
                    if ri != rj {
                        class[ri.max(rj)] = ri.min(rj);
                        welded = true;
                    }
                }
            }
        }
        if welded {
            let off_vertex = |q: Point| -> f64 {
                face.edges
                    .iter()
                    .flat_map(|e| e.ends.iter())
                    .map(|(v, _)| v.distance(q))
                    .fold(f64::INFINITY, f64::min)
            };
            let rank = |m: usize| {
                (
                    u8::from(!strands[ends[m].0].boundary),
                    off_vertex(ends[m].3),
                )
            };
            let mut node: Vec<Option<usize>> = vec![None; ends.len()];
            let mut weld_of = vec![0.0_f64; ends.len()];
            for (i, &weld) in reach.iter().enumerate() {
                let r = root(&mut class, i);
                weld_of[r] = weld_of[r].max(weld);
                let better = node[r].is_none_or(|k| {
                    let (x, y) = (rank(i), rank(k));
                    x.0.cmp(&y.0).then(x.1.total_cmp(&y.1)).is_lt()
                });
                if better {
                    node[r] = Some(i);
                }
            }
            // A node holding ends that stand apart in space is one vertex
            // there, on every face: a junction reaching each of them, unless
            // one the paving found holds them all already.
            let mut spread_of = vec![0.0_f64; ends.len()];
            let mut held = vec![true; ends.len()];
            for i in 0..ends.len() {
                let r = root(&mut class, i);
                if let Some(k) = node[r] {
                    spread_of[r] = spread_of[r].max(ends[k].3.distance(ends[i].3));
                    held[r] &= junctions.iter().any(|j| {
                        j.at.distance(ends[k].3) <= j.reach && j.at.distance(ends[i].3) <= j.reach
                    });
                }
            }
            for (r, k) in node.iter().enumerate() {
                if let Some(k) = *k
                    && spread_of[r] > tol.confusion() * 1e2
                    && !held[r]
                {
                    welds.push(Junction {
                        at: ends[k].3,
                        reach: spread_of[r] + tol.confusion() * 1e2,
                        onto_vertex: false,
                    });
                }
            }
            let lengths: Vec<f64> = strands
                .iter()
                .map(|st| chart_length(&st.polyline))
                .collect();
            let mut landed: Vec<[Option<usize>; 2]> = vec![[None, None]; strands.len()];
            for i in 0..ends.len() {
                let r = root(&mut class, i);
                let Some(k) = node[r] else {
                    continue;
                };
                let (si, end, at, _, _) = ends[i];
                let target = ends[k].2;
                // A node on a seam stands on both sides of it, and the end
                // is carried to the side it arrives from.
                strands[si].polyline[end] = period_shifts(&face.surface, tol)
                    .into_iter()
                    .map(|(du, dv)| Point2::new(target.x + du, target.y + dv))
                    .min_by(|p, q| p.distance(at).total_cmp(&q.distance(at)))
                    .unwrap_or(target);
                landed[si][usize::from(end != 0)] = Some(r);
            }
            for (si, st) in strands.iter_mut().enumerate() {
                if let [Some(a), Some(b)] = landed[si]
                    && a == b
                    && lengths[si] <= weld_of[a]
                {
                    ogeom_core::diag!(
                        Channel::Strands,
                        "WELD from_a={from_a} fi={fi}: strand {si} is dust at its ends' weld {:.3e}",
                        weld_of[a]
                    );
                    st.polyline = vec![st.polyline[0], st.polyline[0]];
                }
            }
        }
        // A section that runs along one of this face's edges into its end
        // stops within the hug's width of the node the edge was split at,
        // and two sections meeting a rail that grazes this face stop on it
        // as far apart as the rail stays within tolerance of the face. The
        // junction there says each pair is one point. Where a chart unit is
        // shorter than a millimetre, or the graze long, that gap is wider in
        // the chart than the weld, and the section dangles beside the node it
        // ends at. Its end moves onto a boundary node the junction holds, or
        // onto the first section end already there.
        let mut anchors: Vec<(Option<usize>, Point2, Point)> = strands
            .iter()
            .filter(|st| st.boundary && st.polyline.len() >= 2)
            .flat_map(|st| [st.polyline[0], st.polyline[st.polyline.len() - 1]])
            .filter_map(|p| {
                face.surface
                    .point_at(p.x, p.y, tol)
                    .ok()
                    .map(|q| (None, p, q))
            })
            .collect();
        for st in strands
            .iter_mut()
            .filter(|st| matches!(st.tag, Tag::Section { .. }) && st.polyline.len() >= 2)
        {
            let last = st.polyline.len() - 1;
            let mut moved_into: Option<usize> = None;
            for at in [0, last] {
                let p = st.polyline[at];
                if anchors
                    .iter()
                    .any(|(held, n, _)| held.is_none() && n.distance(p) <= face_snap)
                {
                    continue;
                }
                let Ok(q) = face.surface.point_at(p.x, p.y, tol) else {
                    continue;
                };
                let Some((ji, junction)) = junctions
                    .iter()
                    .enumerate()
                    .filter(|(_, j)| j.at.distance(q) <= j.reach)
                    .min_by(|(a, x), (b, y)| {
                        x.at.distance(q).total_cmp(&y.at.distance(q)).then(a.cmp(b))
                    })
                else {
                    continue;
                };
                // The nearest in space the junction also holds, and reached
                // across the chart without leaving it: the chart's midpoint
                // of the two lies midway in space too, where a node in
                // another period's copy would not.
                let anchor = anchors
                    .iter()
                    .filter(|(held, _, s)| {
                        held.is_none_or(|h| h == ji) && junction.at.distance(*s) <= junction.reach
                    })
                    .filter(|(_, n, s)| {
                        face.surface
                            .point_at((n.x + p.x) / 2.0, (n.y + p.y) / 2.0, tol)
                            .is_ok_and(|m| m.distance(q.midpoint(*s)) <= junction.reach)
                    })
                    .min_by(|a, b| {
                        a.2.distance(q)
                            .partial_cmp(&b.2.distance(q))
                            .unwrap_or(core::cmp::Ordering::Equal)
                    })
                    .map(|(_, n, _)| *n);
                match anchor {
                    Some(n) => {
                        st.polyline[at] = n;
                        moved_into = Some(ji);
                    }
                    None => anchors.push((Some(ji), p, q)),
                }
            }
            // A piece of section lying wholly inside one junction, its ends
            // carried onto one node, is that junction: two crossings a
            // shallow angle leaves a hair apart along an edge bound a stub
            // of the section between them, and kept as a strand from the
            // node back to itself it encloses a sliver with no inside.
            // Carried to a point, it is dust, and the faces either side of
            // it drop it alike.
            if let Some(ji) = moved_into
                && st.polyline[0].distance(st.polyline[last]) <= face_snap
            {
                let junction = &junctions[ji];
                let mut inside = true;
                for p in &st.polyline {
                    inside &= face
                        .surface
                        .point_at(p.x, p.y, tol)
                        .is_ok_and(|q| junction.at.distance(q) <= junction.reach);
                }
                if inside {
                    st.polyline = vec![st.polyline[0], st.polyline[0]];
                }
            }
        }
        if diag::enabled(Channel::Strands) {
            ogeom_core::diag!(
                Channel::Strands,
                "FACE-SNAP from_a={from_a} fi={fi}: snap {face_snap:.3e} edges {:?}",
                face.edges
                    .iter()
                    .map(|e| (
                        format!("{:.2e}", e.tolerance),
                        format!("{:.2e}", e.ends_tolerance)
                    ))
                    .collect::<Vec<_>>()
            );
            for (si, st) in strands.iter().enumerate() {
                let tag = match st.tag {
                    Tag::Boundary { edge, range } => format!("Boundary e{edge} {range:?}"),
                    Tag::Contact { contact, range } => format!("Contact c{contact} {range:?}"),
                    Tag::Section { section, range } => format!("Section s{section} {range:?}"),
                    Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                };
                let (a, b) = (st.polyline[0], st.polyline[st.polyline.len() - 1]);
                ogeom_core::diag!(
                    Channel::Strands,
                    "STRAND from_a={from_a} fi={fi} {si}: boundary={} {tag} pts={} {a:?} .. {b:?}",
                    st.boundary,
                    st.polyline.len()
                );
            }
        }
        Ok((strands, face_snap, welds))
    };
    // The identity of a piece across faces: a boundary or contact piece by
    // its edge node, a section piece by its section: the same section is
    // split the same way on both faces it cuts, and a piece that is dust
    // in a sphere's chart at the pole must be dust on the plane it also
    // lies in.
    let node_of = |face: &GFace, tag: &Tag| -> Option<(usize, usize, (f64, f64))> {
        match tag {
            Tag::Boundary { edge, range } => {
                Some((0, face.edges[*edge].node.index() as usize, *range))
            }
            Tag::Contact { contact, range } => {
                Some((1, contacts[*contact].node.index() as usize, *range))
            }
            Tag::Section { section, range } => Some((2, *section, *range)),
            Tag::Pole { .. } => None,
        }
    };
    /// A face's strands, its snap and the junctions its welds make,
    /// waiting to be arranged.
    type Prepared = Option<(Vec<Strand<Tag>>, f64, Vec<Junction>)>;
    // A face the other solid leaves alone (no section on it, no contact
    // into it, no edge of it split, no coincident partner, no pole, no
    // edge short enough to be anyone's dust) is one piece, itself, and is
    // neither split nor arranged: its rings are its own wires and its state
    // is read once. On a large solid and a small tool that is nearly every
    // face.
    let mut quiet: [Vec<bool>; 2] = [vec![true; ga.faces.len()], vec![true; gb.faces.len()]];
    for s in &sections {
        quiet[0][s.face_a] = false;
        quiet[1][s.face_b] = false;
    }
    for c in &contacts {
        quiet[usize::from(!c.target_from_a)][c.target_face] = false;
    }
    for t in &tangents {
        quiet[0][t.face_a] = false;
        quiet[1][t.face_b] = false;
    }
    for (side, solid, same) in [(0_usize, &ga, &same_a), (1, &gb, &same_b)] {
        for (fi, face) in solid.faces.iter().enumerate() {
            if !same[fi].is_empty()
                || !face.poles.is_empty()
                || face.edges.iter().any(|e| {
                    paves.contains_key(&e.node)
                        || e.ends[0].0.distance(e.ends[1].0) <= QUIET_EDGE
                        || e.bound.diagonal() <= QUIET_EDGE
                })
            {
                quiet[side][fi] = false;
            }
        }
    }
    let mut prepared: [Vec<Prepared>; 2] = [Vec::new(), Vec::new()];
    let mut dust: Vec<(usize, usize, (f64, f64))> = Vec::new();
    for (side, solid, from_a) in [(0_usize, &ga, true), (1, &gb, false)] {
        // A sheet's tool only cuts: its own faces are not split.
        if side == 1 && operands != Operands::Solids {
            break;
        }
        for (fi, face) in solid.faces.iter().enumerate() {
            ogeom_core::progress::checkpoint()?;
            if quiet[side][fi] {
                prepared[side].push(Some((Vec::new(), 0.0, Vec::new())));
                continue;
            }
            let (strands, snap, welds) = strands_of(from_a, fi, face)?;
            for st in &strands {
                if st.polyline.len() >= 2
                    && chart_length(&st.polyline) <= snap
                    && let Some(key) = node_of(face, &st.tag)
                {
                    dust.push(key);
                }
            }
            prepared[side].push(Some((strands, snap, welds)));
        }
    }
    ogeom_core::progress::stage("boolean: split arrange");
    let same_key = |a: &(usize, usize, (f64, f64)), b: &(usize, usize, (f64, f64))| {
        a.0 == b.0
            && a.1 == b.1
            && (a.2.0 - b.2.0).abs() <= tol.parametric()
            && (a.2.1 - b.2.1).abs() <= tol.parametric()
    };
    // The other solid's boundary, prepared once per side. It is asked once
    // per face piece, and what it costs to prepare (every face's trimming
    // rings, polylined) does not depend on the point being asked about.
    // Rebuilt per question, the preparation costs hundreds of times the
    // question.
    // A sheet bounds nothing, so nothing is classified against one.
    let boundaries = [
        if matches!(operands, Operands::SheetBySheet | Operands::Itself) {
            None
        } else {
            Some(ogeom_algo::SolidBoundary::of_kept(
                model,
                &gb.solid,
                tol.confusion() * 1e4,
                gb.closed,
                tol,
            )?)
        },
        // A solid of one shell whose edges setting faces aside read closed
        // is not asked again.
        if operands == Operands::Solids {
            Some(ogeom_algo::SolidBoundary::of_kept(
                model,
                &ga.solid,
                tol.confusion() * 1e4,
                ga.closed,
                tol,
            )?)
        } else {
            None
        },
    ];
    // Each face is split and its pieces classified on its own: it reads
    // the two solids, the sections, the contacts and the junctions found
    // so far, and writes only its own pieces and the junctions it adds. So
    // the faces run in parallel, and their findings are joined in face
    // order, the order one thread would have made them in.
    type FaceWork = std::sync::Mutex<Prepared>;
    let mut work: Vec<(bool, usize, FaceWork)> = Vec::new();
    for (side, slots) in prepared.iter_mut().enumerate() {
        for (fi, slot) in slots.iter_mut().enumerate() {
            let Some(prepared_face) = slot.take() else {
                ogeom_bail!(Construction, "a face was prepared twice");
            };
            work.push((side == 0, fi, std::sync::Mutex::new(Some(prepared_face))));
        }
    }
    let found = ogeom_core::parallel::map_ordered(
        &work,
        |_, (from_a, fi, prepared_face)| -> OgeomResult<(Vec<FacePiece>, Vec<Junction>)> {
            ogeom_core::progress::checkpoint()?;
            let (from_a, fi) = (*from_a, *fi);
            let (own, other, boundary) = if from_a {
                (&ga, &gb.solid, boundaries[0].as_ref())
            } else {
                (&gb, &ga.solid, boundaries[1].as_ref())
            };
            let face = &own.faces[fi];
            let mut junctions: Vec<Junction> = Vec::new();
            let Some((mut strands, mut face_snap, welds)) = prepared_face
                .lock()
                .map_err(|_| ogeom_core::ogeom_err!(Construction, "a face's work was poisoned"))?
                .take()
            else {
                ogeom_bail!(Construction, "a face was prepared twice");
            };
            junctions.extend(welds);
            if quiet[usize::from(!from_a)][fi] {
                if let Some(piece) = quiet_piece(model, from_a, fi, face, boundary, tol)? {
                    return Ok((vec![piece], Vec::new()));
                }
                // Read on the other solid after all: split it as any face.
                let welds;
                (strands, face_snap, welds) = strands_of(from_a, fi, face)?;
                junctions.extend(welds);
            }
            let mut pieces: Vec<FacePiece> = Vec::new();
            // A strand's ends in space, or none for a pole.
            let space_ends = |face: &GFace, tag: &Tag| -> OgeomResult<Option<(Point, Point)>> {
                Ok(Some(match tag {
                    Tag::Boundary { edge, range } => {
                        let e = &face.edges[*edge];
                        (
                            e.curve.point_at(range.0, tol)?,
                            e.curve.point_at(range.1, tol)?,
                        )
                    }
                    Tag::Contact { contact, range } => {
                        let c = &contacts[*contact];
                        (
                            c.curve.point_at(range.0, tol)?,
                            c.curve.point_at(range.1, tol)?,
                        )
                    }
                    Tag::Section { section, range } => {
                        let sec = &sections[*section];
                        let domain = sec.curve.domain();
                        (
                            sec.curve
                                .point_at(at_param(range.0, domain, sec.closed), tol)?,
                            sec.curve
                                .point_at(at_param(range.1, domain, sec.closed), tol)?,
                        )
                    }
                    Tag::Pole { .. } => return Ok(None),
                }))
            };
            // A strand this face's own weld collapses is one point in space
            // as well: its ends become one junction, or the pieces either
            // side of it end on two vertices the chart says are one, and a
            // loop cut at a seam a few microns from its own start closes in
            // the chart and stays open in space.
            for st in &strands {
                if st.polyline.len() >= 2
                    && chart_length(&st.polyline) <= face_snap
                    && let Some((from, to)) = space_ends(face, &st.tag)?
                    && from.distance(to) > tol.confusion() * 1e2
                {
                    junctions.push(Junction {
                        at: from.midpoint(to),
                        reach: from.distance(to) / 2.0 + tol.confusion() * 1e2,
                        onto_vertex: false,
                    });
                }
            }
            // A piece some other chart already collapsed collapses here too:
            // its ends become one node, every neighbour meeting them moves
            // onto it, and one junction owns the span in space so the
            // rebuilt vertices agree on every face.
            let forced: Vec<usize> = strands
                .iter()
                .enumerate()
                .filter(|(_, st)| {
                    st.polyline.len() >= 2
                        && chart_length(&st.polyline) > face_snap
                        && node_of(face, &st.tag)
                            .is_some_and(|k| dust.iter().any(|d| same_key(d, &k)))
                })
                .map(|(i, _)| i)
                .collect();
            if !forced.is_empty() {
                let mut reps: Vec<(Point2, usize)> = Vec::new();
                let canon = |p: Point2, reps: &mut Vec<(Point2, usize)>| -> usize {
                    match reps.iter().position(|(q, _)| q.distance(p) <= face_snap) {
                        Some(i) => i,
                        None => {
                            reps.push((p, reps.len()));
                            reps.len() - 1
                        }
                    }
                };
                fn root(reps: &mut [(Point2, usize)], mut i: usize) -> usize {
                    while reps[i].1 != i {
                        reps[i].1 = reps[reps[i].1].1;
                        i = reps[i].1;
                    }
                    i
                }
                for &i in &forced {
                    let st = &strands[i];
                    let a = canon(st.polyline[0], &mut reps);
                    let b = canon(st.polyline[st.polyline.len() - 1], &mut reps);
                    let (ra, rb) = (root(&mut reps, a), root(&mut reps, b));
                    if ra != rb {
                        reps[rb].1 = ra;
                    }
                    let Some((from, to)) = space_ends(face, &st.tag)? else {
                        continue;
                    };
                    ogeom_core::diag!(
                        Channel::Strands,
                        "DUST from_a={from_a} fi={fi}: strand {i} collapses with its partners, {from:?} -> {to:?}"
                    );
                    junctions.push(Junction {
                        at: from.midpoint(to),
                        reach: from.distance(to) / 2.0 + tol.confusion() * 1e2,
                        onto_vertex: false,
                    });
                }
                for (i, st) in strands.iter_mut().enumerate() {
                    if forced.contains(&i) {
                        continue;
                    }
                    let last = st.polyline.len() - 1;
                    for end in [0, last] {
                        if let Some(k) = reps
                            .iter()
                            .position(|(q, _)| q.distance(st.polyline[end]) <= face_snap)
                        {
                            let r = root(&mut reps, k);
                            st.polyline[end] = reps[r].0;
                        }
                    }
                }
                let mut index = 0_usize;
                strands.retain(|_| {
                    let keep = !forced.contains(&index);
                    index += 1;
                    keep
                });
            }
            // A hole nothing comes near is walked on its own, and the face's
            // walk is spared it: a face crossed by a small tool can hold
            // hundreds of holes the tool never reaches.
            let lone = lone_holes(model, face, &strands, face_snap)?;
            let every: Vec<usize>;
            let (walked, places, lone): (&[Strand<Tag>], &[usize], &[arrange::Lone<Tag>]) =
                match &lone {
                    Some((rest, places, lone)) => (rest, places, lone),
                    None => {
                        every = (0..strands.len()).collect();
                        (&strands, &every, &[])
                    }
                };
            let split = match arrange_pieces(walked, places, face_snap, lone) {
                Ok(split) => split,
                Err(err) => {
                    if diag::enabled(Channel::Arrange) {
                        ogeom_core::diag!(
                            Channel::Arrange,
                            "ARRANGE from_a={from_a} fi={fi} failed: {err}"
                        );
                        for (si, st) in strands.iter().enumerate() {
                            let (a, b) = (st.polyline[0], st.polyline[st.polyline.len() - 1]);
                            let describe = |curve: &Curve, range: (f64, f64)| -> String {
                                let mid = curve
                                    .point_at(f64::midpoint(range.0, range.1), tol)
                                    .map(|p| format!("({:.4},{:.4},{:.4})", p.x, p.y, p.z))
                                    .unwrap_or_default();
                                match curve {
                                    Curve::Circle(c) => format!(
                                        "circle centre {:?} r {:.4} mid {mid}",
                                        c.circle().centre(),
                                        c.circle().radius()
                                    ),
                                    other => {
                                        format!("{:?} mid {mid}", core::mem::discriminant(other))
                                    }
                                }
                            };
                            let geometry = match st.tag {
                                Tag::Boundary { edge, range } => {
                                    describe(&face.edges[edge].curve, range)
                                }
                                Tag::Contact { contact, range } => {
                                    format!(
                                        "{} window {:?} along {:?} node {}",
                                        describe(&contacts[contact].curve, range),
                                        contacts[contact].crange,
                                        contact_along[contact],
                                        contacts[contact].node.index()
                                    )
                                }
                                _ => String::new(),
                            };
                            ogeom_core::diag!(Channel::Arrange, "    {geometry}");
                            ogeom_core::diag!(
                                Channel::Arrange,
                                "  strand {si} boundary={} {:?} {a:?} .. {b:?}",
                                st.boundary,
                                match st.tag {
                                    Tag::Boundary { edge, range } =>
                                        format!("Boundary e{edge} {range:?}"),
                                    Tag::Contact { contact, range } =>
                                        format!("Contact c{contact} {range:?}"),
                                    Tag::Section { section, range } =>
                                        format!("Section s{section} {range:?}"),
                                    Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                                }
                            );
                        }
                    }
                    return Err(err);
                }
            };
            for mut piece in split {
                let Some(boundary) = boundary else {
                    let probe = face.surface.point_at(piece.first.x, piece.first.y, tol)?;
                    pieces.push(FacePiece {
                        from_a,
                        face: fi,
                        rings: piece.rings,
                        outlines: piece.outlines,
                        probe,
                        state: PieceState::Unread,
                        covered: false,
                        whole: false,
                    });
                    continue;
                };
                // Where a piece stands is asked at its interior probes in
                // turn. The first is the roomiest, and usually the only one
                // needed. The rest are for the piece that merely *touches*
                // the other solid, whose roomiest probe can land on the
                // contact and read neither in nor out.
                //
                // A partner face is asked before the classifier, not after.
                // Where a piece sits on a surface the other solid also
                // carries, the partner *is* the answer, and getting there
                // through the classifier means every ray grazing the shared
                // face, its whole fan of directions exhausted, before it
                // reports the On the partner list already knew, which on a
                // part whose bore is refilled by its own cylinder is orders
                // of magnitude slower.
                //
                // A probe every ray from which is ambiguous (one a hair off
                // the other solid's faces, where they run along this one)
                // says nothing either way, and the next is asked. Only a
                // piece none of whose probes can be read is refused.
                let partners = if from_a { &same_a[fi] } else { &same_b[fi] };
                let mut chosen = None;
                let mut unread = None;
                // What a probe reads: on a coincident partner, or what the
                // other solid's boundary says.
                let ask = |candidate: Point2| -> OgeomResult<(Point, OgeomResult<Containment>)> {
                    let at = face.surface.point_at(candidate.x, candidate.y, tol)?;
                    let shared = !partners.is_empty()
                        && partners.iter().any(|&pi| {
                            let partner = if from_a { &gb.faces[pi] } else { &ga.faces[pi] };
                            chart_point_within(partner, at, partner_reach(face, partner, tol), tol)
                                .is_some()
                        });
                    if shared {
                        return Ok((at, Ok(Containment::On)));
                    }
                    match boundary.holds(model, at, tol) {
                        Err(e @ ogeom_core::OgeomError::NotDone(_)) => Ok((at, Err(e))),
                        Err(e) => Err(e),
                        says => Ok((at, says)),
                    }
                };
                // The roomiest probe settles most pieces on its own, so it is
                // asked before the rest are found, wherever the rest could not
                // change whether it is asked first: on a face that is not a
                // plane every probe is kept, and on a plane bounded by lines
                // and circles its exact rings' count for it does not depend on
                // the rings drawn for the rest. Where it reads neither in nor
                // out, every probe is asked in turn below, as if it had not
                // been.
                let first_kept = match face.surface {
                    SurfaceGeometry::Plane(_) => {
                        let drawn = piece.rings.iter().flatten().all(|t| {
                            fine_strand(face, &sections, from_a, &t.tag).is_some_and(|s| {
                                matches!(s.0, PlanarCurve::Line(_) | PlanarCurve::Circle(_))
                            })
                        });
                        drawn
                            .then(|| {
                                FineRings::of(face, &piece.rings, &sections, from_a, None, tol)
                            })
                            .flatten()
                            .and_then(|every| {
                                let reached = rings_near(&piece.outlines, &[piece.first]);
                                every.inside_beside(piece.first, &reached)
                            })
                    }
                    // On a curved face, a piece on a partner's surface has
                    // its roomiest probe held to its rings drawn finely: a
                    // thin piece's outline chords bow past the piece, and a
                    // probe beyond them reads the other solid off the
                    // partner, as if the piece lay there.
                    // Where the rings cannot be drawn finely (a contact
                    // strand among them), a probe off the face itself is
                    // not the piece's either.
                    _ if !partners.is_empty() => {
                        match FineRings::of_curved(face, &piece.rings, &sections, from_a, tol) {
                            Some(every) => Some(every.inside(piece.first)),
                            None => Some(on_own_face(face, piece.first, tol)?),
                        }
                    }
                    _ => Some(true),
                };
                if first_kept == Some(true)
                    && let (at, Ok(says @ (Containment::In | Containment::Out))) = ask(piece.first)?
                {
                    chosen = Some((piece.first, at, says));
                }
                // A planar piece's probes are held to its exact rings, and
                // where none of them is inside those, the rings offer their
                // own: a sliver narrower than the outline's chords bow has
                // every coarse probe in a neighbour.
                // Only the rings a probe's ray can cross are drawn finely: on a
                // face with hundreds of holes the rest are far from every probe.
                let mut fine = None;
                if chosen.is_none() {
                    let mut interiors: Vec<Point2> = piece.interiors().to_vec();
                    let reached = rings_near(&piece.outlines, &interiors);
                    fine =
                        FineRings::of(face, &piece.rings, &sections, from_a, Some(&reached), tol)
                            .or_else(|| {
                                (!partners.is_empty())
                                    .then(|| {
                                        FineRings::of_curved(
                                            face,
                                            &piece.rings,
                                            &sections,
                                            from_a,
                                            tol,
                                        )
                                    })
                                    .flatten()
                            });
                    if fine.is_none() && !partners.is_empty() {
                        let mut kept = Vec::with_capacity(interiors.len());
                        for p in &interiors {
                            if on_own_face(face, *p, tol)? {
                                kept.push(*p);
                            }
                        }
                        // Where none stands on the face as read, the reading
                        // is not trusted over the piece's own outline.
                        if !kept.is_empty() {
                            interiors = kept;
                        }
                    }
                    if let Some(fine) = &fine {
                        interiors.retain(|p| fine.inside(*p));
                        if interiors.is_empty()
                            && let Some(every) =
                                FineRings::of(face, &piece.rings, &sections, from_a, None, tol)
                                    .or_else(|| {
                                        FineRings::of_curved(
                                            face,
                                            &piece.rings,
                                            &sections,
                                            from_a,
                                            tol,
                                        )
                                    })
                        {
                            interiors =
                                arrange::interior_points_of(&every.points(), tol.parametric());
                        }
                    }
                    for candidate in &interiors {
                        let says = match ask(*candidate)? {
                            (at, Ok(says)) => (at, says),
                            (_, Err(e)) => {
                                unread = Some(e);
                                continue;
                            }
                        };
                        if chosen.is_none() || !matches!(says.1, Containment::On) {
                            chosen = Some((*candidate, says.0, says.1));
                        }
                        if !matches!(says.1, Containment::On) {
                            break;
                        }
                    }
                }
                // On the attempt after kept pieces failed to close, a piece
                // is read at the centre of its area where that lies inside
                // it: the roomiest probe can stand by a corner, and a face
                // leaning microradians off the other solid's, the two
                // sharing an edge, reads within rounding of it there.
                if read_at_centre
                    && matches!(chosen, Some((_, _, Containment::In | Containment::Out)))
                    && let Some(centre) = arrange::centroid_of(&piece.outlines)
                    && arrange::inside_rings(&piece.outlines, centre)
                {
                    if fine.is_none() {
                        let reached = rings_near(&piece.outlines, &[centre]);
                        fine = FineRings::of(
                            face,
                            &piece.rings,
                            &sections,
                            from_a,
                            Some(&reached),
                            tol,
                        );
                    }
                    if fine.as_ref().is_none_or(|f| f.inside(centre)) {
                        let at = face.surface.point_at(centre.x, centre.y, tol)?;
                        if let Ok(says @ (Containment::In | Containment::Out)) =
                            boundary.holds(model, at, tol)
                        {
                            chosen = Some((centre, at, says));
                        }
                    }
                }
                let Some((interior, probe, said)) = chosen else {
                    if let Some(e) = unread {
                        return Err(e);
                    }
                    ogeom_bail!(
                        Construction,
                        "a piece of a face has no interior point to classify at"
                    );
                };
                // On the last attempt, a piece thinner than the band the
                // other solid is read within, every probe of it on that
                // solid's boundary, is read by the strands bounding it.
                let thin = if read_thin && matches!(said, Containment::On) {
                    let theirs = if from_a { &gb } else { &ga };
                    thin_by_strands(face, &piece.rings, &piece.outlines, &contacts, theirs, tol)?
                } else {
                    Thin::Wide
                };
                match thin {
                    Thin::Reads(state) => {
                        pieces.push(FacePiece {
                            from_a,
                            face: fi,
                            rings: piece.rings,
                            outlines: piece.outlines,
                            probe,
                            state,
                            covered: false,
                            whole: false,
                        });
                        continue;
                    }
                    Thin::Disagree => ogeom_bail!(
                        NotDone,
                        "a piece thinner than the boolean reads the other solid \
                         within lies beside that solid's edges on sides that \
                         disagree on whether it is inside: it is dust the \
                         boolean does not resolve"
                    ),
                    Thin::Wide | Thin::Undecided => {}
                }
                let state = match said {
                    Containment::In => PieceState::In,
                    Containment::Out => PieceState::Out,
                    Containment::On => {
                        // On the other boundary: same-domain contact. The
                        // partner face on the shared surface decides whether
                        // the two materials lie on the same side or oppose.
                        let own_normal = outward_normal(face, interior, tol)?;
                        let mut resolved = None;
                        for &pi in partners {
                            let partner = if from_a { &gb.faces[pi] } else { &ga.faces[pi] };
                            let Some(uv) = chart_point_within(
                                partner,
                                probe,
                                partner_reach(face, partner, tol),
                                tol,
                            ) else {
                                continue;
                            };
                            let theirs = outward_normal(partner, uv, tol)?;
                            resolved = Some(if own_normal.dot(theirs) > 0.0 {
                                PieceState::OnAligned
                            } else {
                                PieceState::OnOpposed
                            });
                            break;
                        }
                        let Some(state) = resolved else {
                            // On with no partner containing the probe: the
                            // band's generosity may have read proximity as
                            // coincidence. Ask again at a width where it
                            // cannot. What stays on the boundary there is a
                            // genuine contact, settled below.
                            match ogeom_algo::classify_in_solid_exact_banded(
                                model,
                                other,
                                probe,
                                tol.confusion() * 10.0,
                                tol,
                            )? {
                                Containment::In => {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        outlines: piece.outlines,
                                        probe,
                                        state: PieceState::In,
                                        covered: false,
                                        whole: false,
                                    });
                                    continue;
                                }
                                Containment::Out => {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        outlines: piece.outlines,
                                        probe,
                                        state: PieceState::Out,
                                        covered: false,
                                        whole: false,
                                    });
                                    continue;
                                }
                                Containment::On => {}
                            }
                            // Every probe the piece offered stands on the
                            // other solid's boundary, and none on a coincident
                            // face. Where the two touch only along lines or at
                            // points the piece does not cross the other
                            // boundary there (a crossing would have split it),
                            // so it lies wholly on one side, and any point of
                            // it off the contact says which. Points at golden
                            // fractions of the piece's extent miss a contact
                            // that met its halves and quarters. Where they
                            // read on the boundary too, the contact covers a
                            // region, which is not this case.
                            let mut off_contact = None;
                            for candidate in
                                arrange::off_contact_points(&piece.outlines, tol.parametric())
                            {
                                if !arrange::inside_rings(&piece.outlines, candidate)
                                    || fine.as_ref().is_some_and(|f| !f.inside(candidate))
                                {
                                    continue;
                                }
                                let at = face.surface.point_at(candidate.x, candidate.y, tol)?;
                                match boundary.holds(model, at, tol) {
                                    Ok(Containment::In) => {
                                        off_contact = Some(PieceState::In);
                                        break;
                                    }
                                    Ok(Containment::Out) => {
                                        off_contact = Some(PieceState::Out);
                                        break;
                                    }
                                    Ok(Containment::On)
                                    | Err(ogeom_core::OgeomError::NotDone(_)) => {}
                                    Err(e) => return Err(e),
                                }
                            }
                            if let Some(state) = off_contact {
                                pieces.push(FacePiece {
                                    from_a,
                                    face: fi,
                                    rings: piece.rings,
                                    outlines: piece.outlines,
                                    probe,
                                    state,
                                    covered: false,
                                    whole: false,
                                });
                                continue;
                            }
                            // Still on it, on the boolean's second attempt
                            // (see [`settled`]): a face all but coplanar with
                            // the other's, crossing it at a hair's angle, keeps
                            // a sliver within tolerance of it with no
                            // coincident partner. Just off the piece on both
                            // sides the other solid answers plainly: the
                            // same on both, the piece is that. Inside only
                            // on the piece's own material side, the other's
                            // face runs along it backed the same way, the
                            // aligned contact, which the first solid's piece
                            // speaks for. Inside only beyond it, opposed.
                            if settle {
                                let off = own_normal * (tol.confusion() * 1e3);
                                let side = |q: Point| {
                                    ogeom_algo::classify_in_solid_exact_banded(
                                        model,
                                        other,
                                        q,
                                        tol.confusion() * 10.0,
                                        tol,
                                    )
                                };
                                let (beyond, within) = (side(probe + off)?, side(probe - off)?);
                                let read = match (beyond, within) {
                                    (Containment::In, Containment::In) => {
                                        Some((PieceState::In, false))
                                    }
                                    (Containment::Out, Containment::Out) => {
                                        Some((PieceState::Out, false))
                                    }
                                    (Containment::Out, Containment::In) => {
                                        Some((PieceState::OnAligned, !from_a))
                                    }
                                    (Containment::In, Containment::Out) => {
                                        Some((PieceState::OnOpposed, false))
                                    }
                                    _ => None,
                                };
                                if let Some((state, covered)) = read {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        outlines: piece.outlines,
                                        probe,
                                        state,
                                        covered,
                                        whole: false,
                                    });
                                    continue;
                                }
                            }
                            ogeom_core::diag!(
                                Channel::Wire,
                                "ON-CONTACT piece of {} face {fi} probe {probe:?} partners {partners:?}",
                                if from_a { "A" } else { "B" }
                            );
                            ogeom_bail!(
                                NotDone,
                                "a piece lies on the other solid's boundary \
                             at every point it was asked at, with no \
                             coincident partner face to compare sides \
                             against; a contact over a region of faces not \
                             recognised as coincident is refused; see the \
                             open items in docs/PLAN.md"
                            );
                        };
                        state
                    }
                };
                pieces.push(FacePiece {
                    from_a,
                    face: fi,
                    rings: piece.rings,
                    outlines: piece.outlines,
                    probe,
                    state,
                    covered: false,
                    whole: false,
                });
            }
            // A face nothing crosses, its boundary unsplit, is one piece:
            // itself.
            let untouched = strands.iter().all(|st| match &st.tag {
                Tag::Boundary { edge, range } => {
                    let full = face.edges[*edge].crange;
                    (range.0.min(range.1) - full.0.min(full.1)).abs() <= tol.parametric()
                        && (range.0.max(range.1) - full.0.max(full.1)).abs() <= tol.parametric()
                }
                Tag::Pole { .. } => true,
                Tag::Section { .. } | Tag::Contact { .. } => false,
            });
            if untouched
                && let [piece] = pieces.as_mut_slice()
                && matches!(
                    piece.state,
                    PieceState::In | PieceState::Out | PieceState::Unread
                )
            {
                piece.whole = true;
            }

            Ok((pieces, junctions))
        },
    );
    for face in found {
        let (face_pieces, face_junctions) = face?;
        pieces.extend(face_pieces);
        junctions.extend(face_junctions);
    }
    let coincident = mark_covered_coincidences(&ga, &gb, &mut pieces, tol);
    if diag::enabled(Channel::Arrange) {
        for (fi, partners) in same_a.iter().enumerate() {
            if !partners.is_empty() {
                ogeom_core::diag!(Channel::Arrange, "SAME a{fi} with b{partners:?}");
            }
        }
        for (ci, contact) in contacts.iter().enumerate() {
            ogeom_core::diag!(
                Channel::Arrange,
                "CONTACT c{ci} onto {}{} over {:?} prange {:?} node {} pcurve {:?}",
                if contact.target_from_a { "a" } else { "b" },
                contact.target_face,
                contact.crange,
                contact.prange,
                contact.node.index(),
                contact.pcurve
            );
        }
        for (i, p) in pieces.iter().enumerate() {
            let own = if p.from_a {
                &ga.faces[p.face]
            } else {
                &gb.faces[p.face]
            };
            let kind = match &own.surface {
                SurfaceGeometry::Plane(_) => "plane",
                SurfaceGeometry::Cylinder(_) => "cyl",
                SurfaceGeometry::Torus(_) => "torus",
                SurfaceGeometry::BSpline(_) => "bspline",
                _ => "other",
            };
            ogeom_core::diag!(
                Channel::Arrange,
                "PIECE {i} from_a={} face={} {kind} state={:?} covered={} probe={:?} rings {:?}",
                p.from_a,
                p.face,
                p.state,
                p.covered,
                p.probe,
                p.rings.iter().map(Vec::len).collect::<Vec<_>>()
            );
        }
    }
    ogeom_core::progress::stage("boolean: classified");
    // Merged again: the arrangement adds junctions of its own (a strand
    // that collapsed with its partners becomes one) after the first
    // merge, and two of those at one triple point, each within the other's
    // reach, weld the ends of one band to two vertices with a hairline
    // between them.
    let junctions = merge_junctions(junctions);
    Ok(GeneralFused {
        a: ga,
        b: gb,
        sections,
        contacts,
        tangents,
        pieces,
        coincident,
        junctions,
    })
}

// --- rebuilding --------------------------------------------------------------

/// Everything the rebuild shares across pieces.
struct Rebuild<'m> {
    model: &'m mut Model,
    /// World surface ids, minted once per source face.
    surfaces_a: Vec<Option<ogeom_topo::SurfaceId>>,
    surfaces_b: Vec<Option<ogeom_topo::SurfaceId>>,
    /// Vertices shared by position: a wire's connectivity is checked by node
    /// identity, so two sub-edges meeting at a point must *name* the same
    /// vertex, not merely coincide there.
    vertices: Vec<(Point, Shape)>,
    /// `vertices` binned by position.
    vertex_bins: bins::Bins,
    /// The entries a bin query found, kept between queries.
    near: Vec<usize>,
    /// The widest tolerance any vertex in `vertices` holds: every widening
    /// of one goes through [`Rebuild::vertex`] or [`Rebuild::widen`], so an
    /// end is compared only with the vertices this far and a weld from it.
    widest: f64,
    /// How far two honest descriptions of one junction may sit apart: a
    /// hundred confusions as the floor, widened to three times the loosest
    /// contact edge's or fitted section's own tolerance when one took part:
    /// the strand's own, the crossing's, and the edge's it stopped on.
    weld: f64,
    /// The paving's junctions, each minted as a vertex the first time a
    /// strand end lands inside it.
    junctions: Vec<(Junction, Option<Shape>)>,
    /// `junctions` binned by position.
    junction_bins: bins::Bins,
    /// The widest reach of any junction.
    junction_reach: f64,
    /// Vertices minted from junctions at an edge's own end vertex, which a
    /// tangent section's chart image is bent onto.
    onto_vertex: ogeom_core::FastSet<ogeom_topo::TShapeId>,
    /// Boundary edges rebuilt whole, by the edge they were: the face on the
    /// other side of one finds it built and adds its own pcurve, where a
    /// second build would leave the sew a twin to find.
    whole_edges: ogeom_core::FastMap<EdgeKey, Shape>,
    /// The vertices of edges shared with faces set aside, which the rebuild
    /// names as they are and never changes, each where it stands.
    kept: ogeom_core::FastMap<ogeom_topo::TShapeId, Point>,
    /// Whether a strand end asked a kept vertex to reach further than it
    /// does.
    strained: bool,
}

impl Rebuild<'_> {
    fn vertex(&mut self, p: Point, tol: Tolerances) -> Shape {
        // A junction several paves described is one vertex owning the span
        // they disagree by, and every strand end inside it names that
        // vertex, whichever pave its own trim stopped at.
        // The rim is as generous as the positional weld below: two strand
        // ends a weld apart must not fall on opposite sides of it.
        let slack = self.weld.max(tol.confusion() * 1e2);
        ogeom_core::diag!(Channel::Wire, "VERTEX ask {p:?}");
        let reach = self.junction_reach;
        let inside = |slot: &usize| {
            let j = &self.junctions[*slot].0;
            j.at.distance(p) <= j.reach + slack
        };
        let mut near = core::mem::take(&mut self.near);
        let first = if self.junction_bins.near_into(p, reach + slack, &mut near) {
            near.iter().copied().find(inside)
        } else {
            (0..self.junctions.len()).find(inside)
        };
        if let Some(slot) = first {
            let junction = self.junctions[slot].0;
            let shape = match &self.junctions[slot].1 {
                Some(shape) => shape.clone(),
                None => {
                    let shape = make_vertex(self.model, junction.at).shape;
                    self.junctions[slot].1 = Some(shape.clone());
                    self.remember(junction.at, &shape);
                    shape
                }
            };
            self.near = near;
            if junction.onto_vertex {
                self.onto_vertex.insert(shape.node());
            }
            // The junction owns its span, and an end welded in from the rim
            // widens it by what it actually sat off by, as any weld does.
            let gap = junction.at.distance(p);
            if let Some(node) = self.model.node_mut(&shape)
                && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
            {
                // A confusion of slack past the measured gap: an edge built
                // against this vertex measures the same gap through its own
                // rounding, and a tolerance equal to it fails by an ulp.
                data.tolerance = data
                    .tolerance
                    .widen_to(junction.reach.max(gap) + tol.confusion());
            }
            self.note_tolerance(&shape);
            return shape;
        }
        // The weld reach covers what the inputs may honestly disagree by: a
        // boundary curve that is itself a fitted intersection from an
        // earlier boolean carries a couple of microns of slop, and two
        // descriptions of one junction arrive that far apart. A hundred
        // confusions at millimetre tolerances is ten microns, below any
        // feature this pipeline can resolve, and every weld wider than
        // confusion is recorded on the vertex, not papered over.
        // And as far again as the vertex already owns: a vertex widened by
        // the descriptions it has taken in claims its point within that
        // tolerance, an end arriving within that and the weld of it may be
        // the same point, and a boolean cannot show it is not. Several
        // descriptions of one junction (three bands meeting at an apex)
        // can arrive in turn, each within the last's reach and none within
        // the first's, and measured from the first alone the wire round
        // the band's end has two vertices where it needs one. So every end
        // taken in is remembered where it arrived, and the next end is
        // measured from the description nearest it.
        let floor = self.weld.max(tol.confusion() * 1e2);
        if !self
            .vertex_bins
            .near_into(p, floor + self.widest, &mut near)
        {
            near.extend(0..self.vertices.len());
        }
        let found = near
            .iter()
            .filter_map(|&index| {
                let (q, shape) = &self.vertices[index];
                let own = self
                    .model
                    .node(shape)
                    .and_then(|n| n.data().as_vertex())
                    .map_or(0.0, |d| d.tolerance.get());
                let gap = q.distance(p);
                (gap <= floor.max(own + floor)).then_some((gap, index))
            })
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal))
            .map(|(gap, index)| (gap, self.vertices[index].1.clone()));
        self.near = near;
        if let Some((gap, shape)) = found {
            // Two descriptions of one junction may disagree by a general
            // crossing's residual. The vertex's tolerance is where that
            // disagreement is recorded, so the sub-edges built against
            // either description still reach it honestly.
            let (at, own) = self
                .model
                .node(&shape)
                .and_then(|n| n.data().as_vertex())
                .map_or((p, 0.0), |d| (d.point, d.tolerance.get()));
            let off = at.distance(p);
            if let Some(stands) = self.kept.get(&shape.node()) {
                // A kept vertex is never widened: an end it does not reach
                // as it stands strains the rebuild.
                let off = stands.distance(p);
                self.strained |= off > tol.confusion() && off + tol.confusion() > own;
            } else if off > tol.confusion()
                && let Some(node) = self.model.node_mut(&shape)
                && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
            {
                data.tolerance = data.tolerance.widen_to(off + tol.confusion());
            }
            self.note_tolerance(&shape);
            if gap > tol.confusion() {
                self.remember(p, &shape);
            }
            return shape;
        }
        let shape = make_vertex(self.model, p).shape;
        self.remember(p, &shape);
        shape
    }

    /// Take in `shape` as described at `p`.
    fn remember(&mut self, p: Point, shape: &Shape) {
        self.vertex_bins.insert(p, self.vertices.len());
        self.vertices.push((p, shape.clone()));
        self.note_tolerance(shape);
    }

    /// Keep `widest` covering `shape`'s tolerance.
    fn note_tolerance(&mut self, shape: &Shape) {
        if let Some(data) = self.model.node(shape).and_then(|n| n.data().as_vertex()) {
            self.widest = self.widest.max(data.tolerance.get());
        }
    }

    /// Widen a vertex this rebuild handed out. A kept vertex is not
    /// widened: one that would need it strains the rebuild.
    fn widen(&mut self, vertex: &Shape, to: f64) -> OgeomResult<()> {
        if self.kept.contains_key(&vertex.node()) {
            let own = self
                .model
                .node(vertex)
                .and_then(|n| n.data().as_vertex())
                .map_or(0.0, |d| d.tolerance.get());
            self.strained |= own < to;
            return Ok(());
        }
        self.model.widen(vertex, ogeom_core::Tolerance::new(to)?)?;
        self.note_tolerance(vertex);
        Ok(())
    }

    fn surface_id(
        &mut self,
        fused: &GeneralFused,
        from_a: bool,
        face: usize,
    ) -> ogeom_topo::SurfaceId {
        let slot = if from_a {
            &mut self.surfaces_a[face]
        } else {
            &mut self.surfaces_b[face]
        };
        if let Some(id) = slot {
            return *id;
        }
        let own = if from_a {
            fused.a.faces[face].own_surface
        } else {
            fused.b.faces[face].own_surface
        };
        if let Some(id) = own {
            *slot = Some(id);
            return id;
        }
        let surface = if from_a {
            fused.a.faces[face].surface.clone()
        } else {
            fused.b.faces[face].surface.clone()
        };
        let id = self.model.geometry_mut().add_surface(surface);
        *slot = Some(id);
        id
    }
}

/// The sub-edges built of each strand of a piece, keyed by strand, each
/// with the parameter range it covers, in the order they were built.
type StrandEdges = ogeom_core::FastMap<(usize, u8), Vec<((f64, f64), Shape)>>;

/// Build one piece as a face, orientation matching its source face's side.
fn build_piece(
    rebuild: &mut Rebuild,
    fused: &GeneralFused,
    piece: &FacePiece,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let own = if piece.from_a { &fused.a } else { &fused.b };
    let face = &own.faces[piece.face];
    let surface_id = rebuild.surface_id(fused, piece.from_a, piece.face);

    // Sub-edges cached within the piece, so a seam used from both sides is
    // one edge appearing twice. Keyed by strand, each with the ranges built
    // of it in order.
    let mut cache: StrandEdges = ogeom_core::FastMap::default();
    let mut wires = Vec::new();
    for ring in &piece.rings {
        let mut edges = Vec::with_capacity(ring.len());
        // The edges shared with faces set aside this ring has taken, each
        // the way it runs.
        let mut taken: Vec<(EdgeKey, bool)> = Vec::new();
        for traversal in ring {
            // An edge shared with a face set aside is itself, whole: nothing
            // reaches it to split it, and it carries its pcurve on the
            // face's own surface already. Its strands (a closed edge's two
            // halves) stand for it once each way it is walked.
            if let Tag::Boundary { edge, .. } = &traversal.tag
                && let Some(kept) = own.shared.get(&face.edges[*edge].node)
            {
                if face.own_surface != Some(surface_id) {
                    rebuild.strained = true;
                }
                let key = (face.edges[*edge].node, traversal.reversed);
                if !taken.contains(&key) {
                    taken.push(key);
                    edges.push(if traversal.reversed {
                        kept.reversed()
                    } else {
                        kept.clone()
                    });
                }
                continue;
            }
            let (key_edge, key_kind, range) = match &traversal.tag {
                Tag::Boundary { edge, range } => (*edge, 0_u8, *range),
                Tag::Section { section, range } => (*section, 1, *range),
                Tag::Contact { contact, range } => (*contact, 2, *range),
                Tag::Pole { pole, range } => (*pole, 3, *range),
            };
            let near = |a: (f64, f64), b: (f64, f64)| {
                (a.0 - b.0).abs() <= tol.parametric() && (a.1 - b.1).abs() <= tol.parametric()
            };
            if diag::enabled(Channel::Wire) && (range.1 - range.0).abs() < 1e-4 {
                let kind = match &traversal.tag {
                    Tag::Boundary { .. } => "boundary",
                    Tag::Section { .. } => "section",
                    Tag::Contact { .. } => "contact",
                    Tag::Pole { .. } => "pole",
                };
                ogeom_core::diag!(
                    Channel::Wire,
                    "DUST {kind} {key_edge} range {range:?} in piece from_a={} face={}",
                    piece.from_a,
                    piece.face
                );
            }
            let built = if let Some((_, shape)) = cache
                .get(&(key_edge, key_kind))
                .and_then(|built| built.iter().find(|(r, _)| near(*r, range)))
            {
                shape.clone()
            } else {
                let shape = build_sub_edge(
                    rebuild,
                    fused,
                    piece.from_a,
                    face,
                    surface_id,
                    (&traversal.tag, traversal.nodes),
                    tol,
                )?;
                cache
                    .entry((key_edge, key_kind))
                    .or_default()
                    .push((range, shape.clone()));
                shape
            };
            // A piece whose two ends welded to one vertex and whose whole
            // length lies within the vertex's reach is that vertex's dust
            // (a rim's last fraction of a micron before a pole corner) and
            // no edge of the wire. A closed edge on one vertex, a full
            // circle, reaches far from it and stays.
            if !matches!(traversal.tag, Tag::Pole { .. })
                && let Some((v0, v1)) = ogeom_algo::edge_vertices(rebuild.model, &built)?
                && v0.node() == v1.node()
                && let Some(vd) = rebuild
                    .model
                    .node(&v0)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| (d.point, d.tolerance.get()))
                && let Some(ed) = rebuild.model.node(&built).and_then(|n| n.data().as_edge())
                && let Some(ogeom_topo::EdgeRepr::Curve3d {
                    curve, range: r3, ..
                }) = ed.curve3d()
                && let Some(g) = rebuild.model.geometry().curve(*curve)
                && g.point_at(f64::midpoint(r3.0, r3.1), tol)?.distance(vd.0)
                    <= vd.1.max(tol.confusion() * 1e4)
            {
                ogeom_core::diag!(
                    Channel::Wire,
                    "DUST piece closing on one vertex at {:?} dropped from piece from_a={} face={}",
                    vd.0,
                    piece.from_a,
                    piece.face
                );
                continue;
            }
            edges.push(if traversal.reversed {
                built.reversed()
            } else {
                built
            });
        }
        // A ring every piece of which is dust stands within its vertices'
        // reach of one point: a hole of it is no hole, and a boundary of
        // it bounds no face.
        if edges.is_empty() {
            if wires.is_empty() {
                ogeom_core::diag!(
                    Channel::Wire,
                    "DUST piece from_a={} face={}: its boundary closes on one point",
                    piece.from_a,
                    piece.face
                );
                return Ok(None);
            }
            continue;
        }
        let wire = match make_wire(rebuild.model, &edges, tol) {
            Ok(w) => w.shape,
            Err(e) => {
                if diag::enabled(Channel::Wire) {
                    ogeom_core::diag!(
                        Channel::Wire,
                        "WIRE FAIL piece from_a={} face={}",
                        piece.from_a,
                        piece.face
                    );
                    for built in &edges {
                        if let Some((a, b)) = ogeom_algo::edge_vertices(rebuild.model, built)? {
                            let at = |v: &Shape| {
                                rebuild
                                    .model
                                    .node(v)
                                    .and_then(|n| n.data().as_vertex())
                                    .map(|d| (d.point, d.tolerance.get()))
                            };
                            ogeom_core::diag!(
                                Channel::Wire,
                                "   built {:?}{}: v{} {:?} -> v{} {:?}",
                                built.node(),
                                if built.orientation() == ogeom_topo::Orientation::Reversed {
                                    " rev"
                                } else {
                                    ""
                                },
                                a.node().index(),
                                at(&a),
                                b.node().index(),
                                at(&b)
                            );
                        }
                    }
                    for t in ring {
                        let tag = match &t.tag {
                            Tag::Boundary { edge, range } => format!(
                                "Boundary e{edge} {range:?} tol {:.2e}",
                                face.edges[*edge].tolerance
                            ),
                            Tag::Contact { contact, range } => format!(
                                "Contact c{contact} {range:?} tol {:.2e}",
                                fused.contacts[*contact].tolerance
                            ),
                            Tag::Section { section, range } => format!(
                                "Section s{section} {range:?} tol {:.2e}",
                                fused.sections[*section].tolerance
                            ),
                            Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                        };
                        ogeom_core::diag!(Channel::Wire, "   {tag} reversed={}", t.reversed);
                    }
                }
                return Err(e);
            }
        };
        wires.push(wire);
    }
    let built = match make_face_on(rebuild.model, surface_id, &wires, tol) {
        Ok(b) => b.shape,
        Err(e) => {
            if diag::enabled(Channel::Wire) {
                ogeom_core::diag!(
                    Channel::Wire,
                    "FACE FAIL piece from_a={} face={}: {e}",
                    piece.from_a,
                    piece.face
                );
                for (wi, wire) in wires.iter().enumerate() {
                    for edge in rebuild.model.ordered_children_of(wire)? {
                        if let Some((a, b)) = ogeom_algo::edge_vertices(rebuild.model, &edge)? {
                            let at = |v: &Shape| {
                                rebuild
                                    .model
                                    .node(v)
                                    .and_then(|n| n.data().as_vertex())
                                    .map(|d| d.point)
                            };
                            ogeom_core::diag!(
                                Channel::Wire,
                                "   wire {wi} edge {:?}{}: v{} {:?} -> v{} {:?}",
                                edge.node(),
                                if edge.orientation() == ogeom_topo::Orientation::Reversed {
                                    " rev"
                                } else {
                                    ""
                                },
                                a.node().index(),
                                at(&a),
                                b.node().index(),
                                at(&b)
                            );
                        }
                    }
                }
            }
            return Err(e);
        }
    };
    Ok(Some(
        if face.face.orientation() == ogeom_topo::Orientation::Reversed {
            built.reversed()
        } else {
            built
        },
    ))
}

/// Build the exact sub-edge a tag names, pcurves attached, with the nodes
/// the face's arrangement closed its ring on at the strand's two ends.
fn build_sub_edge(
    rebuild: &mut Rebuild,
    fused: &GeneralFused,
    from_a: bool,
    face: &GFace,
    surface_id: ogeom_topo::SurfaceId,
    (tag, nodes): (&Tag, Option<[Point2; 2]>),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    match tag {
        Tag::Pole { pole, range } => {
            // A pole rebuilds as what it was: one vertex, an edge with no
            // curve bounded by it twice, and the chart line that says where
            // it runs in this face's parameters.
            let p = &face.poles[*pole];
            let at = rebuild.vertex(p.point, tol);
            let model = &mut *rebuild.model;
            let mut data = ogeom_topo::EdgeData::new();
            data.degenerate = true;
            let built = model.add_edge(data, &[at.clone(), at])?;
            ogeom_algo::attach_pcurve(
                model,
                &built,
                p.pcurve.clone(),
                surface_id,
                Location::identity(),
                *range,
            )?;
            Ok(built)
        }
        Tag::Boundary { edge, range } => {
            let e = &face.edges[*edge];
            let whole = e.other_side.is_none()
                && (range.0 - e.crange.0).abs() <= tol.parametric()
                && (range.1 - e.crange.1).abs() <= tol.parametric();
            // On the face's own surface, which other faces may stand on
            // too, an edge its node states as a seam stays one: the faces
            // either side each walk their own column of it, which two
            // plain pcurves on the one surface could not tell apart.
            let seam = match &e.seam {
                Some(pair) if e.other_side.is_none() && face.own_surface == Some(surface_id) => {
                    Some(pair)
                }
                _ => None,
            };
            if whole && let Some(built) = rebuild.whole_edges.get(&e.node).cloned() {
                let attached = rebuild
                    .model
                    .node(&built)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.pcurve_on(surface_id).is_some());
                match seam {
                    // The face on the other side, on the same surface,
                    // attached the seam with the edge's first build.
                    Some(_) if attached => {}
                    Some((forward, reversed)) => ogeom_algo::attach_seam(
                        rebuild.model,
                        &built,
                        forward.clone(),
                        reversed.clone(),
                        surface_id,
                        Location::identity(),
                        e.prange,
                    )?,
                    None => ogeom_algo::attach_pcurve(
                        rebuild.model,
                        &built,
                        e.pcurve.clone(),
                        surface_id,
                        Location::identity(),
                        e.prange,
                    )?,
                }
                return Ok(built);
            }
            let from = e.curve.point_at(range.0, tol)?;
            let to = e.curve.point_at(range.1, tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, (*e.curve).clone(), *range, &v0, &v1, tol)?.shape;
            // A piece of a tolerant edge is the same curve with the same
            // honest radius: a fitted rail's stated slop must survive the
            // split, or the next boolean over this solid measures the rail
            // against a tolerance it never had. The ends own it too: a
            // junction on a tolerant curve is a junction to the curve's own
            // resolution, and the wires rebuilt through it meet within that.
            // A fuzzy boolean's confusion is its fuzz, which is no floor
            // here: a piece is born at the smallest tolerance an edge
            // states, and keeps whatever its edge stated above that.
            if e.tolerance > ogeom_core::Tolerance::MIN.get() {
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(e.tolerance);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, e.tolerance)?;
                }
            }
            let model = &mut *rebuild.model;
            let sub_p = (
                rescale(range.0, e.crange, e.prange),
                rescale(range.1, e.crange, e.prange),
            );
            match (&e.other_side, seam) {
                (None, Some((forward, reversed))) => ogeom_algo::attach_seam(
                    model,
                    &built,
                    forward.clone(),
                    reversed.clone(),
                    surface_id,
                    Location::identity(),
                    sub_p,
                )?,
                (None, None) => ogeom_algo::attach_pcurve(
                    model,
                    &built,
                    e.pcurve.clone(),
                    surface_id,
                    Location::identity(),
                    sub_p,
                )?,
                (Some((other, _)), _) => {
                    // A seam: both sides attach over the one range a seam
                    // carries, and an occurrence picks its side by its
                    // orientation.
                    ogeom_algo::attach_seam(
                        model,
                        &built,
                        e.pcurve.clone(),
                        other.clone(),
                        surface_id,
                        Location::identity(),
                        sub_p,
                    )?;
                }
            }
            // An edge's pcurve need only lie on its curve within the
            // tolerance, not keep pace with it: lifted, it may run a little
            // ahead or behind. Cut short, the piece's end is read at one
            // parameter on both, and that lag stands between the piece's
            // vertex and its image in the chart. The piece states it.
            if !whole
                && let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &built, tol)?
                && let Some(data) = model.node(&built).and_then(|n| n.data().as_edge())
                && gap > data.tolerance.get()
            {
                let widened = ogeom_core::Tolerance::new(gap * (1.0 + 1e-6))?;
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen(widened);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, widened.get())?;
                }
            }
            if whole {
                rebuild.whole_edges.insert(e.node, built.clone());
            }
            Ok(built)
        }
        Tag::Contact { contact, range } => {
            let c = &fused.contacts[*contact];
            let from = c.curve.point_at(range.0, tol)?;
            let to = c.curve.point_at(range.1, tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, c.curve.clone(), *range, &v0, &v1, tol)?.shape;
            // A tolerant contact's pieces and ends own its stated slop, as a
            // boundary's and a section's do: three kinds of strand close one
            // ring, and every junction meets at the strands' own honesty.
            if c.tolerance > ogeom_core::Tolerance::MIN.get() {
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(c.tolerance);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, c.tolerance)?;
                }
            }
            let model = &mut *rebuild.model;
            // The stored image keeps its own window. The attached copy names
            // the sub-window this piece covers under the proportional map.
            let sub_p = (
                rescale(range.0, c.crange, c.prange),
                rescale(range.1, c.crange, c.prange),
            );
            let mid = c.pcurve.point_at(f64::midpoint(sub_p.0, sub_p.1), tol)?;
            let trim: Vec<&[Point2]> = face.trim_lines(tol).iter().map(Vec::as_slice).collect();
            let folded = fold_inside(mid, &face.surface, &trim, tol);
            let shifted = c
                .pcurve
                .transformed(&ogeom_math::Transform2::translation(folded - mid), tol)?;
            ogeom_algo::attach_pcurve(
                model,
                &built,
                shifted,
                surface_id,
                Location::identity(),
                sub_p,
            )?;
            // A contact's curve lies on the target's surface only to the
            // weld or the fuzz (a rim on a drum a micron off its own), and
            // its image there, lifted, stands off the curve by as much. The
            // piece states that, and its ends hold it.
            if let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &built, tol)?
                && let Some(data) = model.node(&built).and_then(|n| n.data().as_edge())
                && gap > data.tolerance.get()
            {
                let widened = ogeom_core::Tolerance::new(gap * (1.0 + 1e-6))?;
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen(widened);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, widened.get())?;
                }
            }
            Ok(built)
        }
        Tag::Section { section, range } => {
            let s = &fused.sections[*section];
            let domain = s.curve.domain();
            let (f0, f1) = folded_range(*range, domain, s.closed);
            let from = s.curve.point_at(at_param(f0, domain, s.closed), tol)?;
            let to = s.curve.point_at(at_param(f1, domain, s.closed), tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, s.curve.clone(), (f0, f1), &v0, &v1, tol)?.shape;
            // A piece welded into a junction some way off its own end (a
            // section tangent to an edge at the edge's vertex, stopped
            // where it began to hug the edge) owns that gap as the vertex
            // does: the face it bounds meets its neighbour there only that
            // closely, in its chart as in space.
            let at = |v: &Shape| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| d.point)
            };
            let ends = [at(&v0), at(&v1)];
            let welded = [(ends[0], from), (ends[1], to)]
                .iter()
                .filter_map(|(v, p)| v.map(|v| v.distance(*p)))
                .fold(0.0, f64::max);
            let own = s.tolerance.max(welded);
            // A fitted section's pieces and their ends own the section's
            // stated slop, exactly as a tolerant boundary's do: the ring
            // they close alternates between the two, and both sides must
            // meet at the junction's own resolution.
            if own > ogeom_core::Tolerance::MIN.get() {
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(own);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, own)?;
                }
            }
            let model = &mut *rebuild.model;
            // The section's pcurve is unwrapped across any seam. The face's
            // triangulator lives in one chart, so the attached copy is folded
            // home by the same period shift the arrangement gave this
            // strand's polyline (to the side of the face's seam its trim is
            // on), decided by the sub-range's midpoint, so an
            // endpoint sitting exactly on the chart's edge stays on the side
            // the arc's body is.
            let pcurve = if from_a { &s.pc_a } else { &s.pc_b };
            let mid = pcurve.point_at(f64::midpoint(f0, f1), tol)?;
            let trim: Vec<&[Point2]> = face.trim_lines(tol).iter().map(Vec::as_slice).collect();
            let folded = fold_inside(mid, &face.surface, &trim, tol);
            let shifted =
                pcurve.transformed(&ogeom_math::Transform2::translation(folded - mid), tol)?;
            // The chart image ends where the vertex is, as the edge does in
            // space to its tolerance: a welded end left where the section
            // stopped leaves the face's outline open by the weld in its
            // chart, wider there than any later arrangement's snap.
            let loose = s.tolerance.max(tol.confusion() * 1e2);
            let bent = [
                rebuild.onto_vertex.contains(&v0.node()),
                rebuild.onto_vertex.contains(&v1.node()),
            ];
            let targets = [
                ends[0].filter(|v| bent[0] && v.distance(from) > loose),
                ends[1].filter(|v| bent[1] && v.distance(to) > loose),
            ];
            // Elsewhere the image ends on the node the arrangement closed
            // the ring on, where its own end stands further off it than any
            // later arrangement's snap and the node lifts nearer the piece's
            // end in space than the image's end does. A fitted image can lag
            // its curve by a few microns. Beside a pole, where a micron spans
            // a hundred times as much of the chart, that leaves it a
            // ten-thousandth of a radian off the seam column its neighbour
            // ends on, and the face's ring open there to every later
            // boolean. A piece whose own end stops short of its node (a loop
            // whose last hair was welded away) keeps the image of that end.
            let mut onto = [None, None];
            for (k, ((target, t), at)) in targets.iter().zip([f0, f1]).zip([from, to]).enumerate() {
                let end = shifted.point_at(t, tol)?;
                if let Some(target) = target {
                    onto[k] = Some(chart_foot(&face.surface, *target, end, tol)?);
                    continue;
                }
                let Some(node) = nodes.map(|n| n[k]) else {
                    continue;
                };
                let node = period_shifts(&face.surface, tol)
                    .into_iter()
                    .map(|(du, dv)| Point2::new(node.x + du, node.y + dv))
                    .min_by(|p, q| p.distance(end).total_cmp(&q.distance(end)))
                    .unwrap_or(node);
                let vertex = if k == 0 { &v0 } else { &v1 };
                let radius = model.tolerance_of(vertex)?.map_or(0.0, |r| r.get());
                if node.distance(end) > PARAM_SNAP
                    && !on_collapsed_row(&face.surface, node, radius, tol)
                    && face.surface.point_at(node.x, node.y, tol)?.distance(at)
                        < face.surface.point_at(end.x, end.y, tol)?.distance(at)
                {
                    onto[k] = Some(node);
                }
            }
            let shifted = if onto.iter().any(Option::is_some) {
                pcurve_onto_ends(&shifted, (f0, f1), &face.surface, onto, tol)?
            } else {
                shifted
            };
            ogeom_algo::attach_pcurve(
                model,
                &built,
                shifted,
                surface_id,
                Location::identity(),
                (f0, f1),
            )?;
            // The section's image was fitted through samples and may have
            // been bent onto its welded ends; lifted, it stands off the
            // curve by what is measured along the whole piece, and the
            // piece and its ends state that.
            if let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &built, tol)?
                && let Some(data) = model.node(&built).and_then(|n| n.data().as_edge())
                && gap > data.tolerance.get()
            {
                let widened = ogeom_core::Tolerance::new(gap * (1.0 + 1e-6))?;
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen(widened);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, widened.get())?;
                }
            }
            Ok(built)
        }
    }
}

/// The point of `surface`'s chart nearest `end` that lifts to the foot of
/// `target` on the surface: the foot may land a turn away from `end` on a
/// periodic chart, whether or not the surface's wrapper says it is
/// periodic, and a shifted candidate counts where the surface agrees it is
/// the foot.
fn chart_foot(
    surface: &SurfaceGeometry,
    target: Point,
    end: Point2,
    tol: Tolerances,
) -> OgeomResult<Point2> {
    let foot = ogeom_algo::project_on_surface_from(surface, target, (end.x, end.y), tol)?;
    let raw = Point2::new(foot.parameters.0, foot.parameters.1);
    let turn = core::f64::consts::TAU;
    let mut shifts = period_shifts(surface, tol);
    for du in [-turn, 0.0, turn] {
        for dv in [-turn, 0.0, turn] {
            shifts.push((du, dv));
        }
    }
    Ok(shifts
        .into_iter()
        .map(|(du, dv)| Point2::new(raw.x + du, raw.y + dv))
        .filter(|c| {
            surface
                .point_at(c.x, c.y, tol)
                .is_ok_and(|q| q.distance(foot.point) <= tol.confusion() * 1e2)
        })
        .min_by(|a, b| {
            a.distance(end)
                .partial_cmp(&b.distance(end))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .unwrap_or(raw))
}

/// Whether the chart row or column through `at` lies, across the chart's
/// whole width that way, within `radius` of one point in space: a pole's
/// row or an apex's, along which a point moves nowhere, so a pcurve's end
/// is placed on it by the way it arrives and not by a node there.
fn on_collapsed_row(surface: &SurfaceGeometry, at: Point2, radius: f64, tol: Tolerances) -> bool {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let Ok((du, dv)) = surface.d1_at(at.x, at.y, tol) else {
        return true;
    };
    du.magnitude() * (u1 - u0).abs() <= radius || dv.magnitude() * (v1 - v0).abs() <= radius
}

/// A pcurve over `range` bent at its ends onto the chart points `targets`,
/// where given, and refitted at the same parameters, at every span of the
/// image, so the image stays same-parameter with its edge to the
/// correction's size.
///
/// Each end's correction fades linearly to nothing at the other end, and,
/// read on the surface's first derivatives, nowhere moves the image
/// further in space than it moves the end: where a chart step costs more
/// in space than it does at the end (away from a pole the end stands
/// beside), the step is scaled down by as much. A ten-thousandth of a
/// radian round a ball's axis is a few microns beside its pole and a
/// ten-thousandth of its radius at its equator.
fn pcurve_onto_ends(
    pcurve: &PlanarCurve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    targets: [Option<Point2>; 2],
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    const SAMPLES: usize = 32;
    let ends = [
        pcurve.point_at(range.0, tol)?,
        pcurve.point_at(range.1, tol)?,
    ];
    // The length in space of a chart step `by` taken at `at`.
    let cost = |at: Point2, by: ogeom_math::Vector2| -> Option<f64> {
        let (du, dv) = surface.d1_at(at.x, at.y, tol).ok()?;
        Some((du * by.x + dv * by.y).magnitude())
    };
    let mut deltas = [ogeom_math::Vector2::new(0.0, 0.0); 2];
    let mut at_end = [None; 2];
    for (k, target) in targets.iter().enumerate() {
        let Some(target) = target else {
            continue;
        };
        deltas[k] = *target - ends[k];
        at_end[k] = cost(ends[k], deltas[k]);
    }
    // How much of an end's correction a point of the image takes, its
    // linear share scaled so the point moves no further in space than the
    // end does.
    let share = |k: usize, p: Point2, linear: f64| -> f64 {
        let (Some(end), Some(here)) = (at_end[k], cost(p, deltas[k])) else {
            return linear;
        };
        if here > end {
            linear * end / here
        } else {
            linear
        }
    };
    // Refitted at every span of the image and held to the target between
    // the samples as well as at them; the edge states what the bend and
    // the fit leave once the image is attached.
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let stations = ogeom_algo::traced::stations_2d(pcurve, (lo, hi), SAMPLES);
    let fitted = ogeom_algo::traced::fit_traced_2d(
        |t| {
            let along = if hi > lo {
                (t - range.0) / (range.1 - range.0)
            } else {
                0.0
            };
            let p = pcurve.point_at(t, tol)?;
            Ok(p + deltas[0] * share(0, p, 1.0 - along) + deltas[1] * share(1, p, along))
        },
        &stations,
        3,
        tol.confusion() * 10.0,
        tol,
    )?;
    Ok(PlanarCurve::BSpline(fitted.curve))
}

/// Sew kept pieces, demand closure, and nest shells into solids and voids.
///
/// The faces each solid set aside come through as they are where
/// `keep_aside` keeps that solid's outside, and are deleted where it does
/// not.
fn assemble_result(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    keep_aside: [bool; 2],
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_core::progress::stage("boolean: assemble");
    let local = !fused.a.aside.is_empty() || !fused.b.aside.is_empty();
    let mut history = History::new();
    let source_face = |piece: &FacePiece| -> Shape {
        if piece.from_a {
            fused.a.faces[piece.face].face.clone()
        } else {
            fused.b.faces[piece.face].face.clone()
        }
    };

    if kept.is_empty() {
        if (keep_aside[0] && !fused.a.aside.is_empty())
            || (keep_aside[1] && !fused.b.aside.is_empty())
        {
            ogeom_bail!(NotDone, "no piece is kept beside the faces set aside");
        }
        // A legitimate answer: cutting a solid away entirely leaves nothing.
        let empty = model.add_compound(&[])?;
        for piece in &fused.pieces {
            history.delete(&source_face(piece));
        }
        for (face, _) in fused.a.aside.iter().chain(&fused.b.aside) {
            history.delete(face);
        }
        history.modify(a, empty.clone());
        history.modify(b, empty.clone());
        return Ok(Built::new(empty, history));
    }

    let Rebuilt {
        faces,
        settled,
        weld: floor,
        holes,
        kept: built_from,
    } = rebuilt_pieces(model, fused, kept, &mut history, tol)?;
    let sewn = sew_around_within(model, &faces, &settled, SEW_ACROSS.get(), tol)?;
    // Sewing rebuilds the pieces onto shared edges: the result's faces are
    // its faces, reached from the inputs through both steps.
    history = history.followed_by(&sewn.history);
    // What the sew and the passes after it made of each rebuilt face, where
    // faces set aside are put back beside them.
    let mut after = if local {
        sewn.history.clone()
    } else {
        History::new()
    };
    let mut sewn = sewn;
    let dropped = without_membranes(model, &mut sewn.shells, floor, tol)?;
    history = history.followed_by(&dropped);
    let held: ogeom_core::FastSet<ogeom_topo::TShapeId> = fused
        .a
        .shared
        .values()
        .chain(fused.b.shared.values())
        .map(Shape::node)
        .collect();
    let joined = seam_join::join_across_seams(model, &mut sewn.shells, &held, tol)?;
    history = history.followed_by(&joined).without_repeated_images();
    let mut closed = None;
    let mut made = None;
    if local {
        after = after.followed_by(&dropped).followed_by(&joined);
        let (shells, shut, new) = with_aside(
            model,
            fused,
            &built_from,
            &faces,
            &holes,
            &after,
            &sewn.shells,
            keep_aside,
            &mut history,
        )?;
        sewn.shells = shells;
        closed = Some(shut);
        made = Some(new);
    }
    if sewn.shells.is_empty() {
        // Membranes all through: what was kept encloses nothing.
        let empty = model.add_compound(&[])?;
        history.modify(a, empty.clone());
        history.modify(b, empty.clone());
        return Ok(Built::new(empty, history));
    }
    for (i, shell) in sewn.shells.iter().enumerate() {
        let shut = match &closed {
            Some(closed) => closed[i],
            None => is_shell_closed(model, shell)?,
        };
        if !shut {
            // Env-gated forensics: the open shell's unshared edges, the
            // question every failure here starts from.
            if diag::enabled(Channel::Arrange) {
                use ogeom_geom::Curve3d as _;
                for edge in ogeom_topo::explore_unique(model, shell, ShapeType::Edge)? {
                    let users = ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))?
                        .iter()
                        .filter(|f| {
                            ogeom_topo::explore_unique(model, f, ShapeType::Edge)
                                .map(|es| es.iter().any(|e2| e2.node() == edge.node()))
                                .unwrap_or(false)
                        })
                        .count();
                    let mut occurrences = 0_usize;
                    for f in ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))? {
                        for wire in model.children_of(&f)? {
                            for e2 in model.children_of(&wire)? {
                                if e2.node() == edge.node() {
                                    occurrences += 1;
                                }
                            }
                        }
                    }
                    if occurrences % 2 == 1
                        && let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge())
                        && let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) =
                            data.curve3d()
                        && let Some(g) = model.geometry().curve(*curve)
                    {
                        let a = g.point_at(range.0, tol)?;
                        let b = g.point_at(range.1, tol)?;
                        ogeom_core::diag!(
                            Channel::Arrange,
                            "  open edge uses={occurrences} faces={users} {:?} range={range:?}: {a:?} -> {b:?}",
                            core::mem::discriminant(g)
                        );
                        // Twins: other edges of the shell with the same ends,
                        // which the sew should have found to be this edge.
                        for other in ogeom_topo::explore_unique(model, shell, ShapeType::Edge)? {
                            if other.node() == edge.node() {
                                continue;
                            }
                            let Some(od) = model.node(&other).and_then(|n| n.data().as_edge())
                            else {
                                continue;
                            };
                            let Some(ogeom_topo::EdgeRepr::Curve3d {
                                curve: oc,
                                range: orange,
                                ..
                            }) = od.curve3d()
                            else {
                                continue;
                            };
                            let Some(og) = model.geometry().curve(*oc) else {
                                continue;
                            };
                            let oa = og.point_at(orange.0, tol)?;
                            let ob = og.point_at(orange.1, tol)?;
                            let same = (oa.distance(a) <= 1e-3 && ob.distance(b) <= 1e-3)
                                || (oa.distance(b) <= 1e-3 && ob.distance(a) <= 1e-3);
                            if same {
                                let om = og.point_at(f64::midpoint(orange.0, orange.1), tol)?;
                                let m = g.point_at(f64::midpoint(range.0, range.1), tol)?;
                                ogeom_core::diag!(
                                    Channel::Arrange,
                                    "    twin {:?} range={orange:?} mid gap {:.2e} tol {:.2e} vs {:.2e}",
                                    core::mem::discriminant(og),
                                    om.distance(m),
                                    od.tolerance.get(),
                                    data.tolerance.get()
                                );
                            }
                        }
                        for f in ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))?
                        {
                            let uses_it = ogeom_topo::explore_unique(model, &f, ShapeType::Edge)
                                .map(|es| es.iter().any(|e2| e2.node() == edge.node()))
                                .unwrap_or(false);
                            if !uses_it {
                                continue;
                            }
                            if let Some(ogeom_topo::NodeData::Face(fd)) =
                                model.node(&f).map(|n| n.data())
                                && let Some(sg) = model.geometry().surface(fd.surface)
                            {
                                let bound = shape_bounds(model, &f, tol)?;
                                ogeom_core::diag!(
                                    Channel::Arrange,
                                    "    used by face on {:?} bound {:?}",
                                    core::mem::discriminant(sg),
                                    bound
                                );
                            }
                        }
                    }
                }
            }
            if diag::enabled(Channel::Arrange) {
                use ogeom_geom::Curve3d as _;
                for f in ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))? {
                    let Some(ogeom_topo::NodeData::Face(fd)) = model.node(&f).map(|n| n.data())
                    else {
                        continue;
                    };
                    let Some(sg) = model.geometry().surface(fd.surface) else {
                        continue;
                    };
                    let edges = ogeom_topo::explore_unique(model, &f, ShapeType::Edge)?;
                    ogeom_core::diag!(
                        Channel::Arrange,
                        "  face {:?} with {} edges",
                        core::mem::discriminant(sg),
                        edges.len()
                    );
                    if matches!(sg, SurfaceGeometry::Sphere(_)) {
                        for e in &edges {
                            if let Some(ed) = model.node(e).and_then(|n| n.data().as_edge())
                                && let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) =
                                    ed.curve3d()
                                && let Some(g) = model.geometry().curve(*curve)
                            {
                                let a = g.point_at(range.0, tol)?;
                                let b = g.point_at(range.1, tol)?;
                                ogeom_core::diag!(
                                    Channel::Arrange,
                                    "    sphere edge {:?} ({:.4},{:.4},{:.4}) -> ({:.4},{:.4},{:.4})",
                                    core::mem::discriminant(g),
                                    a.x,
                                    a.y,
                                    a.z,
                                    b.x,
                                    b.y,
                                    b.z
                                );
                            }
                        }
                    }
                }
            }
            ogeom_bail!(
                NotDone,
                "the kept pieces did not close into a shell; the configuration \
                 is beyond what the boolean currently resolves"
            );
        }
    }

    // Nest: a shell inside another is that solid's void, and a shell inside
    // a void is a solid of its own. The kept pieces are already oriented
    // as the material's boundary, so a void's faces face into it as they
    // are.
    let mut solids = Vec::new();
    for (outer, voids) in nest_shells(model, &sewn.shells, tol)? {
        let mut group = vec![sewn.shells[outer].clone()];
        group.extend(voids.iter().map(|&j| sewn.shells[j].clone()));
        solids.push(model.add_solid(&group)?);
    }
    let result = if solids.len() == 1 {
        solids.remove(0)
    } else {
        model.add_compound(&solids)?
    };
    // A face made here whose boundary does not close in its surface's chart
    // (a section stopped short of where its neighbour starts, at a contact
    // the marcher did not follow through) holds a strand that dangles, and
    // the result would check broken. It is refused instead, for the retries
    // to take another way.
    // Set aside, the faces made are known; otherwise they are the result's
    // faces that are not the operands'.
    let built_faces = match &made {
        Some(made) => made.clone(),
        None => {
            let given: ogeom_core::FastSet<ogeom_topo::TShapeId> = [a, b]
                .into_iter()
                .map(|operand| explore_unique(model, operand, ShapeType::Face))
                .collect::<OgeomResult<Vec<_>>>()?
                .into_iter()
                .flatten()
                .map(|face| face.node())
                .collect();
            explore_unique(model, &result, ShapeType::Face)?
                .into_iter()
                .filter(|f| !given.contains(&f.node()))
                .collect()
        }
    };
    // A wire of a face the boolean gathered that comes through as the same
    // wire (a hole clear of the other solid, put back into the piece that
    // holds it) closes as it did.
    let given_wires: ogeom_core::FastSet<ogeom_topo::TShapeId> = fused
        .a
        .faces
        .iter()
        .chain(&fused.b.faces)
        .map(|gathered| model.ordered_children_of(&gathered.face))
        .collect::<OgeomResult<Vec<_>>>()?
        .into_iter()
        .flatten()
        .map(|wire| wire.node())
        .collect();
    let taken = |wire: &Shape| given_wires.contains(&wire.node());
    for face in &built_faces {
        if !ogeom_algo::chart_breaks(model, face, &taken, tol)?.is_empty() {
            ogeom_bail!(
                NotDone,
                "a face the boolean made does not close in its surface's chart; the \
                 configuration is beyond what the boolean currently resolves"
            );
        }
    }
    // The faces set aside keep the boxes they had: only the faces made
    // are looked at.
    match &made {
        Some(made) => keep_face_boxes(model, made)?,
        None => keep_face_boxes(model, &explore_unique(model, &result, ShapeType::Face)?)?,
    }
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Built::new(result, history))
}

/// Each of `faces` on a ruled surface (a plane, a drum, a cone) keeps its
/// box, found from its edges alone: the next boolean on the result reads
/// every face's box to find the faces it can set aside.
fn keep_face_boxes(model: &Model, faces: &[Shape]) -> OgeomResult<()> {
    for face in faces {
        if model.face_bounds(face).is_some() {
            continue;
        }
        let ruled = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| {
                matches!(
                    s,
                    SurfaceGeometry::Plane(_)
                        | SurfaceGeometry::Cylinder(_)
                        | SurfaceGeometry::Cone(_)
                )
            });
        if ruled {
            ogeom_algo::face_bounds(model, face)?;
        }
    }
    Ok(())
}

/// The kept pieces of a sheet sewn along the edges they share: one shell,
/// or a compound of shells where they fall apart. Nothing is asked of the
/// shells' closure. A sheet that is one face has its history in its
/// pieces; a shell or a compound also becomes the result.
fn assemble_sheet(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[usize],
    sheet: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_core::progress::stage("boolean: assemble");
    let mut history = History::new();
    let whole = model.kind_of(sheet)? != ShapeType::Face;
    if kept.is_empty() {
        let empty = model.add_compound(&[])?;
        for piece in &fused.pieces {
            history.delete(&fused.a.faces[piece.face].face);
        }
        if whole {
            history.modify(sheet, empty.clone());
        }
        return Ok(Built::new(empty, history));
    }
    let kept: Vec<(usize, bool)> = kept.iter().map(|&i| (i, false)).collect();
    let Rebuilt { faces, settled, .. } = rebuilt_pieces(model, fused, &kept, &mut history, tol)?;
    let mut sewn = sew_around(model, &faces, &settled, tol)?;
    history = history.followed_by(&sewn.history);
    let joined = seam_join::join_across_seams(
        model,
        &mut sewn.shells,
        &ogeom_core::FastSet::default(),
        tol,
    )?;
    history = history.followed_by(&joined).without_repeated_images();
    let result = if sewn.shells.len() == 1 {
        sewn.shells.remove(0)
    } else {
        model.add_compound(&sewn.shells)?
    };
    if whole {
        history.modify(sheet, result.clone());
    }
    Ok(Built::new(result, history))
}

/// The kept pieces rebuilt as faces, each recorded in `history` against the
/// face it came from and every face no piece is kept of deleted, whether
/// each is a copy already sewn to all its neighbours, and the distance the
/// rebuild welds their ends within.
fn rebuilt_pieces(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    history: &mut History,
    tol: Tolerances,
) -> OgeomResult<Rebuilt> {
    let source_face = |piece: &FacePiece| -> Shape {
        if piece.from_a {
            fused.a.faces[piece.face].face.clone()
        } else {
            fused.b.faces[piece.face].face.clone()
        }
    };
    let weld = fused
        .contacts
        .iter()
        .map(|c| c.tolerance)
        .chain(fused.sections.iter().map(|s| s.tolerance))
        .fold(0.0_f64, |acc, t| acc.max(honest(t, tol) * 3.0));
    let floor = weld.max(tol.confusion() * 1e2);
    let junction_reach = fused
        .junctions
        .iter()
        .fold(0.0_f64, |acc, j| acc.max(j.reach));
    let mut junction_bins = bins::Bins::new(junction_reach + floor);
    for (index, j) in fused.junctions.iter().enumerate() {
        junction_bins.insert(j.at, index);
    }
    // The vertices' cells are sized to the loosest tolerance a strand
    // brings to its ends, so an end is compared across a few cells, not
    // every vertex.
    let loosest = fused
        .a
        .faces
        .iter()
        .chain(fused.b.faces.iter())
        .flat_map(|f| f.edges.iter().map(|e| e.tolerance))
        .chain(fused.contacts.iter().map(|c| c.tolerance))
        .chain(fused.sections.iter().map(|s| s.tolerance))
        .fold(junction_reach, f64::max);
    let mut rebuild = Rebuild {
        model,
        surfaces_a: vec![None; fused.a.faces.len()],
        surfaces_b: vec![None; fused.b.faces.len()],
        vertices: Vec::new(),
        vertex_bins: bins::Bins::new(floor + loosest),
        near: Vec::new(),
        widest: 0.0,
        weld,
        junctions: fused.junctions.iter().map(|j| (*j, None)).collect(),
        junction_bins,
        junction_reach,
        onto_vertex: ogeom_core::FastSet::default(),
        whole_edges: ogeom_core::FastMap::default(),
        kept: ogeom_core::FastMap::default(),
        strained: false,
    };
    for (at, vertex) in fused
        .a
        .shared_vertices
        .iter()
        .chain(&fused.b.shared_vertices)
    {
        rebuild.remember(*at, vertex);
        rebuild.kept.insert(vertex.node(), *at);
    }
    // A face the operation never touched (nothing crosses it, no edge of it
    // is split, no junction lands at its corners) is rebuilt as an exact
    // copy of itself. Its neighbours across its unsplit edges build those
    // edges from the same curves, so whichever twin the sew keeps is its
    // own boundary.
    let clean = |piece: &FacePiece, model: &Model| -> OgeomResult<bool> {
        if !piece.whole {
            return Ok(false);
        }
        let face = &source_face(piece);
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            let at = vertex.transform(model.datums())?.apply(data.point);
            if fused
                .junctions
                .iter()
                .any(|j| j.at.distance(at) <= j.reach + floor)
            {
                return Ok(false);
            }
        }
        Ok(true)
    };
    // Each kept piece built first: one whose boundary closes on a point is
    // dust, and none of what follows sees it.
    let mut pieces_built: Vec<(usize, bool, Shape)> = Vec::with_capacity(kept.len());
    for &(index, flip) in kept {
        if let Some(built) = build_piece(&mut rebuild, fused, &fused.pieces[index], tol)? {
            pieces_built.push((index, flip, built));
        }
    }
    let kept: Vec<(usize, bool)> = pieces_built.iter().map(|&(i, f, _)| (i, f)).collect();
    let kept = kept.as_slice();
    // One kept turned over (a tool's face lining the cavity a cut leaves)
    // faces the other way, which no copy does.
    let mut copies = Vec::with_capacity(kept.len());
    for &(index, flip) in kept {
        copies.push(!flip && clean(&fused.pieces[index], rebuild.model)?);
    }
    let mut faces = Vec::new();
    let mut kept_sources: ogeom_core::FastSet<Shape> = ogeom_core::FastSet::default();
    let mut built_of: ogeom_core::FastMap<usize, usize> = ogeom_core::FastMap::default();
    for (slot, (index, flip, built)) in pieces_built.into_iter().enumerate() {
        built_of.insert(index, slot);
        let piece = &fused.pieces[index];
        let mut built = built;
        if flip {
            built = built.reversed();
        }
        if copies[slot] {
            history.copy(&source_face(piece), built.clone());
        } else {
            history.modify(&source_face(piece), built.clone());
        }
        kept_sources.insert(source_face(piece));
        faces.push(built);
    }
    if rebuild.strained {
        ogeom_bail!(
            NotDone,
            "a strand ends further from a vertex shared with a face set aside \
             than the vertex reaches"
        );
    }
    let holes = holes_put_back(rebuild.model, fused, kept, tol)?;
    // A piece dropped because the other argument's coincident piece is
    // kept for the same patch: its face became that piece's face.
    let mut stood_for: ogeom_core::FastSet<Shape> = ogeom_core::FastSet::default();
    for &(index_a, index_b) in &fused.coincident {
        let (dropped, standing) = match (built_of.get(&index_a), built_of.get(&index_b)) {
            (Some(&slot), None) => (index_b, slot),
            (None, Some(&slot)) => (index_a, slot),
            _ => continue,
        };
        let source = source_face(&fused.pieces[dropped]);
        if !history
            .modified(&source)
            .iter()
            .any(|f| f.is_same(&faces[standing]))
        {
            history.modify(&source, faces[standing].clone());
        }
        stood_for.insert(source);
    }
    for piece in &fused.pieces {
        let source = source_face(piece);
        if !kept_sources.contains(&source) && !stood_for.contains(&source) {
            history.delete(&source);
        }
    }

    // A copy shares each edge with the face of its own solid across it,
    // built once for both: no edge a copy holds has a twin in its own
    // solid. Where it meets a face of the other solid, its edges may have
    // twins only the sew can find, so only a copy all of whose vertices are
    // held by faces of its own solid alone comes through the sew as it
    // stands.
    let class = |slot: usize| -> u8 {
        match (copies[slot], fused.pieces[kept[slot].0].from_a) {
            (true, true) => 1,
            (true, false) => 2,
            (false, true) => 4,
            (false, false) => 8,
        }
    };
    let mut held: ogeom_core::FastMap<ogeom_topo::TShapeId, u8> = ogeom_core::FastMap::default();
    let mut corners = Vec::with_capacity(faces.len());
    for (slot, face) in faces.iter().enumerate() {
        let vertices = explore_unique(rebuild.model, face, ShapeType::Vertex)?;
        for vertex in &vertices {
            *held.entry(vertex.node()).or_default() |= class(slot);
        }
        corners.push(vertices);
    }
    let settled = corners
        .iter()
        .enumerate()
        .map(|(slot, vertices)| {
            let own = if fused.pieces[kept[slot].0].from_a {
                1 | 4
            } else {
                2 | 8
            };
            copies[slot]
                && vertices
                    .iter()
                    .all(|v| held.get(&v.node()).is_some_and(|by| by & !own == 0))
        })
        .collect();

    Ok(Rebuilt {
        faces,
        settled,
        weld: floor,
        holes,
        kept: kept.to_vec(),
    })
}

/// The shells the result is bounded by, with the faces set aside put back:
/// each face taking back holes left out of its arrangement is made again
/// with them, the faces set aside join the faces sewn where `keep_aside`
/// keeps them, and the shells are found again from the edges they all
/// share. `history` records the faces set aside and the faces made again.
#[allow(
    clippy::too_many_arguments,
    reason = "the assembly's own state, read once"
)]
fn with_aside(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    faces: &[Shape],
    holes: &[Vec<usize>],
    after: &History,
    shells: &[Shape],
    keep_aside: [bool; 2],
    history: &mut History,
) -> OgeomResult<(Vec<Shape>, Vec<bool>, Vec<Shape>)> {
    if kept.len() != faces.len() {
        ogeom_bail!(
            Construction,
            "{} kept pieces for {} rebuilt faces; each face is read with its own piece",
            kept.len(),
            faces.len()
        );
    }
    // Each face made again with its holes, with the wires it walks itself
    // and the faces set aside beside its holes.
    let mut remade: ogeom_core::FastMap<ogeom_topo::TShapeId, Remade> =
        ogeom_core::FastMap::default();
    let mut put_back = History::new();
    // Each face's place: its solid, the place there of the face it came
    // from, and its piece's.
    let mut rank: ogeom_core::FastMap<ogeom_topo::TShapeId, (bool, usize, usize)> =
        ogeom_core::FastMap::default();
    for (slot, built) in faces.iter().enumerate() {
        let piece = &fused.pieces[kept[slot].0];
        let own = if piece.from_a { &fused.a } else { &fused.b };
        let source = &own.faces[piece.face].face;
        let place = (
            !piece.from_a,
            own.order
                .get(&ogeom_topo::SameKey(source.clone()))
                .copied()
                .unwrap_or(usize::MAX),
            slot,
        );
        rank.insert(built.node(), place);
        for image in after.modified(built) {
            rank.insert(image.node(), place);
        }
    }
    for (slot, (built, holes)) in faces.iter().zip(holes).enumerate() {
        if holes.is_empty() {
            continue;
        }
        let image = match after.modified(built) {
            [] if !after.is_deleted(built) => built.clone(),
            [image] => image.clone(),
            _ => ogeom_bail!(NotDone, "a face taking back its holes was split or dropped"),
        };
        let piece = &fused.pieces[kept[slot].0];
        let own = if piece.from_a { &fused.a } else { &fused.b };
        let source = &own.faces[piece.face];
        let side = usize::from(!piece.from_a);
        let wires: Vec<Shape> = holes
            .iter()
            .map(|&h| source.holes_aside[h].wire.clone())
            .collect();
        let joins: Vec<(usize, u32)> = holes
            .iter()
            .flat_map(|&h| source.holes_aside[h].beside.iter().map(move |&j| (side, j)))
            .collect();
        let walk = model.node(&image).map_or(0, |n| n.children().len());
        let made = with_holes(model, &image, &wires)?;
        // A plane's holes lie inside its outer boundary, which alone sets
        // its box: the face with its holes back keeps the box of the face
        // without them.
        ogeom_algo::face_bounds(model, &image)?;
        let found = model.kept_face_box(&image, |_| {
            Err(ogeom_core::ogeom_err!(
                NotDone,
                "the face's box was not kept"
            ))
        });
        if let Ok(found) = found {
            model.kept_face_box(&made, |_| Ok(found))?;
        }
        let copied = after.copy_of(built).is_some() || after.modified(built).is_empty();
        if copied {
            put_back.copy(&image, made.clone());
        } else {
            put_back.modify(&image, made.clone());
        }
        remade.insert(image.node(), (made, walk, joins));
    }
    // Each face with its place: its solid, the place there of the face it
    // came from, and its piece's.
    let mut all = Vec::new();
    for shell in shells {
        for face in model.children_of(shell)? {
            let place = rank
                .get(&face.node())
                .copied()
                .unwrap_or((true, usize::MAX, usize::MAX));
            all.push(match remade.remove(&face.node()) {
                Some((made, walk, joins)) => (
                    place,
                    aside::Member {
                        face: made,
                        aside: None,
                        walk,
                        joins,
                    },
                ),
                None => (
                    place,
                    aside::Member {
                        face,
                        aside: None,
                        walk: usize::MAX,
                        joins: Vec::new(),
                    },
                ),
            });
        }
    }
    *history = core::mem::take(history).followed_by(&put_back);
    // The faces set aside are no image of anything the history holds so
    // far, so they are recorded as they stand.
    for (g, keep, from_b) in [
        (&fused.a, keep_aside[0], false),
        (&fused.b, keep_aside[1], true),
    ] {
        for (face, place) in &g.aside {
            if keep {
                history.copy(face, face.clone());
                all.push((
                    (from_b, *place, 0),
                    aside::Member {
                        face: face.clone(),
                        aside: Some((usize::from(from_b), *place)),
                        walk: 0,
                        joins: Vec::new(),
                    },
                ));
            } else {
                history.delete(face);
            }
        }
    }
    // In the order the solids list the faces they came from, as the
    // general fuse of every face lists them.
    all.sort_by_key(|(place, _)| *place);
    let all: Vec<aside::Member> = all.into_iter().map(|(_, member)| member).collect();
    if all.is_empty() {
        return Ok((Vec::new(), Vec::new(), Vec::new()));
    }
    let made = all
        .iter()
        .filter(|member| member.aside.is_none())
        .map(|member| member.face.clone())
        .collect();
    let (shells, closed) = aside::shells_of(model, &all, [&fused.a.beside, &fused.b.beside])?;
    Ok((shells, closed, made))
}

/// A face made again with its holes, how many of its wires it walks
/// itself, and the faces set aside beside its holes.
type Remade = (Shape, usize, Vec<(usize, u32)>);

/// What [`rebuilt_pieces`] makes of the kept pieces.
struct Rebuilt {
    /// One face per kept piece, in order.
    faces: Vec<Shape>,
    /// Whether each comes through the sew as it stands.
    settled: Vec<bool>,
    /// The distance the rebuild welds ends within.
    weld: f64,
    /// The holes left out of the arrangement each face takes back, by
    /// their place among its source face's.
    holes: Vec<Vec<usize>>,
    /// The kept pieces each face is built from, in order: those given, less
    /// any whose boundary closed on one point.
    kept: Vec<(usize, bool)>,
}

/// The holes left out of each face's arrangement, given to the kept piece
/// holding each.
///
/// A hole left out stands clear of the other solid, so the piece holding it
/// lies outside the other solid. Where the operation keeps what lies
/// outside, the face's one kept piece holds it, or of several the one whose
/// outline it lies within. Where it keeps what lies inside, no kept piece
/// holds one.
fn holes_put_back(
    model: &Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<usize>>> {
    let mut holes = vec![Vec::new(); kept.len()];
    let mut by_face: std::collections::BTreeMap<(bool, usize), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (slot, &(index, _)) in kept.iter().enumerate() {
        let piece = &fused.pieces[index];
        let own = if piece.from_a { &fused.a } else { &fused.b };
        if !own.faces[piece.face].holes_aside.is_empty() {
            by_face
                .entry((piece.from_a, piece.face))
                .or_default()
                .push(slot);
        }
    }
    for ((from_a, fi), slots) in by_face {
        let face = if from_a {
            &fused.a.faces[fi]
        } else {
            &fused.b.faces[fi]
        };
        let outside = slots
            .iter()
            .any(|&slot| fused.pieces[kept[slot].0].state == PieceState::Out);
        if let [slot] = slots.as_slice()
            && outside
        {
            holes[*slot] = (0..face.holes_aside.len()).collect();
            continue;
        }
        for (h, hole) in face.holes_aside.iter().enumerate() {
            let Some(at) = hole_point(model, face, &hole.wire, tol)? else {
                ogeom_bail!(NotDone, "a hole set aside has no point in its face's chart");
            };
            let mut within = slots
                .iter()
                .filter(|&&slot| inside_rings(&fused.pieces[kept[slot].0].outlines, at));
            match (within.next(), within.next()) {
                (Some(&slot), None) if fused.pieces[kept[slot].0].state == PieceState::Out => {
                    holes[slot].push(h);
                }
                (None, _) if !outside => {}
                _ => ogeom_bail!(
                    NotDone,
                    "a hole set aside lies in no one piece standing outside the other solid"
                ),
            }
        }
    }
    Ok(holes)
}

/// A point of a hole's first edge in its face's chart.
fn hole_point(
    model: &Model,
    face: &GFace,
    hole: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Point2>> {
    let Some(surface) = face.own_surface else {
        return Ok(None);
    };
    let Some(edge) = model.children_of(hole)?.into_iter().next() else {
        return Ok(None);
    };
    let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let (pc, range) = match data.pcurve_for(surface, edge.location()) {
        Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
        Some(EdgeRepr::Seam { forward, range, .. }) => (*forward, *range),
        _ => return Ok(None),
    };
    let Some(pcurve) = model.geometry().pcurve(pc) else {
        return Ok(None);
    };
    Ok(Some(pcurve.point_at(f64::midpoint(range.0, range.1), tol)?))
}

/// `image` with the holes `holes` added: a new face on the same surface,
/// facing the same way.
fn with_holes(model: &mut Model, image: &Shape, holes: &[Shape]) -> OgeomResult<Shape> {
    if !image.location().is_identity() {
        ogeom_bail!(
            NotDone,
            "a face taking back its holes stands under a placement"
        );
    }
    let Some(NodeData::Face(data)) = model.node(image).map(|n| n.data().clone()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let mut data = *data;
    data.triangulation = None;
    // The holes' edges are the solid's own and are not widened: a face
    // looser than one of them takes none back. No edge is tighter than the
    // least tolerance.
    if data.tolerance.get() > ogeom_core::Tolerance::MIN.get() {
        for hole in holes {
            for edge in model.children_of(hole)? {
                if model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .is_none_or(|e| e.tolerance.get() < data.tolerance.get())
                {
                    ogeom_bail!(NotDone, "a face is looser than a hole it takes back");
                }
            }
        }
    }
    let mut wires = model.children_of(&Shape::of(image.node()))?;
    wires.extend(holes.iter().cloned());
    Ok(model.add_face(data, &wires)?.oriented(image.orientation()))
}

/// Each shell without its membranes, and the history of their removal.
///
/// Two faces standing apart by less than the weld (two boxes stacked a few
/// microns apart) have their corners welded into one vertex each, and their
/// edges sewn into one: the two faces then bound one shell along the same
/// edges, lying on each other and facing apart, a wall of no thickness with
/// material on both sides. That is the opposed pair the boolean drops where
/// the contact is exact, and it is dropped here the same way. A pair bounded
/// by the same edges but standing apart (a lens's two caps on their one
/// rim) is a solid's own boundary and stays. A shell that is membranes
/// through (a sliver thinner than the weld) encloses nothing and goes.
/// The faces of a pair stand within [`MEMBRANE_WITHIN`] of each other too.
fn without_membranes(
    model: &mut Model,
    shells: &mut Vec<Shape>,
    weld: f64,
    tol: Tolerances,
) -> OgeomResult<History> {
    use ogeom_geom::Transformable as _;
    let mut dropped = History::new();
    let within = MEMBRANE_WITHIN.get();
    let bounding = |model: &Model, face: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
        let mut edges: Vec<ogeom_topo::TShapeId> = explore_unique(model, face, ShapeType::Edge)?
            .iter()
            .filter(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| !d.degenerate)
            })
            .map(Shape::node)
            .collect();
        edges.sort_unstable();
        Ok(edges)
    };
    // Whether `second` lies on `first`, facing away from it.
    let membrane = |model: &Model, first: &Shape, second: &Shape| -> OgeomResult<bool> {
        let (at, n) = ogeom_algo::face_normal(model, first, tol)?;
        let (_, m) = ogeom_algo::face_normal(model, second, tol)?;
        if n.dot(m) > -0.99 * n.magnitude() * m.magnitude() {
            return Ok(false);
        }
        // The first's inner points, every one on the second's surface: a
        // thin lens (a tube's tip poking through a cap, cut off by it) is
        // two faces on the same loop facing apart too, as thin as a weld
        // only near its rim.
        let mesh =
            ogeom_mesh::triangulate_face(model, first, ogeom_mesh::Deflection::default(), tol)?;
        let mut probes: Vec<Point> = mesh
            .triangles
            .iter()
            .map(|t| {
                let [p, q, r] = t.map(|i| mesh.positions[i as usize]);
                Point::from_vector((p.to_vector() + q.to_vector() + r.to_vector()) / 3.0)
            })
            .collect();
        if probes.is_empty() {
            probes.push(at);
        }
        let Some(NodeData::Face(data)) = model.node(second).map(|n| n.data()) else {
            return Ok(false);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(false);
        };
        let placed = surface.transformed(&second.transform(model.datums())?, tol)?;
        // As far as the weld that joined their corners reached: its floor
        // past what the corners themselves own.
        let owned = explore_unique(model, first, ShapeType::Vertex)?
            .iter()
            .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
            .fold(0.0_f64, |acc, d| acc.max(d.tolerance.get()));
        let reach = (weld + owned).min(within);
        for probe in probes {
            let foot = ogeom_algo::project_on_surface(&placed, probe, 16, tol)?;
            if foot.point.distance(probe) > reach {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let mut emptied = Vec::new();
    for (index, shell) in shells.iter_mut().enumerate() {
        let faces = ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))?;
        let edges = faces
            .iter()
            .map(|f| bounding(model, f))
            .collect::<OgeomResult<Vec<_>>>()?;
        // Only faces bounded by the same edges can be a membrane's pair, so
        // the pairs are taken within each such group, in face order.
        let mut alike: ogeom_core::FastMap<&[ogeom_topo::TShapeId], Vec<usize>> =
            ogeom_core::FastMap::default();
        for (i, bounded_by) in edges.iter().enumerate() {
            if !bounded_by.is_empty() {
                alike.entry(bounded_by.as_slice()).or_default().push(i);
            }
        }
        let mut groups: Vec<Vec<usize>> = alike.into_values().filter(|g| g.len() > 1).collect();
        groups.sort_unstable();
        let mut gone = vec![false; faces.len()];
        for group in &groups {
            for (k, &i) in group.iter().enumerate() {
                for &j in &group[k + 1..] {
                    if gone[i] || gone[j] {
                        continue;
                    }
                    if membrane(model, &faces[i], &faces[j])? {
                        gone[i] = true;
                        gone[j] = true;
                    }
                }
            }
        }
        if !gone.contains(&true) {
            continue;
        }
        let mut kept = Vec::new();
        for (face, gone) in faces.iter().zip(&gone) {
            if *gone {
                dropped.delete(face);
            } else {
                kept.push(face.clone());
            }
        }
        if kept.is_empty() {
            emptied.push(index);
        } else {
            *shell = model.add_shell(&kept)?;
        }
    }
    for index in emptied.into_iter().rev() {
        shells.remove(index);
    }
    Ok(dropped)
}

/// Closed shells grouped into solids: each outer shell with the voids
/// directly inside it, as indices into `shells`.
///
/// A shell's depth is how many others enclose it. Even depth is a solid's
/// outer boundary (a free shell, or an island standing in another's
/// cavity), odd depth a void, belonging to the enclosing shell one level
/// out. Enclosure is decided by classifying the inner shell's vertices,
/// not by bounding boxes: a part standing in another's notch sits inside
/// its box and outside its material. A pair whose vertices all lie on the
/// other (shells touching everywhere they could be tested) falls back to
/// the boxes.
fn nest_shells(
    model: &Model,
    shells: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<(usize, Vec<usize>)>> {
    if shells.len() == 1 {
        return Ok(vec![(0, Vec::new())]);
    }
    let mut bounds = Vec::with_capacity(shells.len());
    for shell in shells {
        bounds.push(shape_bounds(model, shell, tol)?);
    }
    // A container's box widened by the tolerance a touch is read within: a
    // void meeting its outer shell at a point (a ball resting against the
    // inside of a wall) may reach a hair past it, and is still within it.
    let reach = tol.confusion() * 1e2;
    let mut boundaries: Vec<Option<ogeom_algo::SolidBoundary>> =
        (0..shells.len()).map(|_| None).collect();
    let mut encloses = vec![vec![false; shells.len()]; shells.len()];
    for i in 0..shells.len() {
        for j in 0..shells.len() {
            if i == j || !bounds[i].expanded(reach).contains_box(&bounds[j]) {
                continue;
            }
            if boundaries[i].is_none() {
                boundaries[i] = Some(ogeom_algo::SolidBoundary::of(
                    model,
                    &shells[i],
                    tol.confusion() * 1e4,
                    tol,
                )?);
            }
            let Some(boundary) = boundaries[i].as_ref() else {
                continue;
            };
            let mut verdict = None;
            for vertex in explore_unique(model, &shells[j], ShapeType::Vertex)? {
                let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
                    continue;
                };
                let at = vertex.transform(model.datums())?.apply(data.point);
                match boundary.holds(model, at, tol)? {
                    Containment::In => verdict = Some(true),
                    Containment::Out => verdict = Some(false),
                    Containment::On => continue,
                }
                break;
            }
            encloses[i][j] = verdict.unwrap_or(true);
        }
    }
    let depth: Vec<usize> = (0..shells.len())
        .map(|j| (0..shells.len()).filter(|&i| encloses[i][j]).count())
        .collect();
    let mut out: Vec<(usize, Vec<usize>)> = (0..shells.len())
        .filter(|&i| depth[i].is_multiple_of(2))
        .map(|i| (i, Vec::new()))
        .collect();
    for j in (0..shells.len()).filter(|&j| !depth[j].is_multiple_of(2)) {
        let parent = (0..shells.len()).find(|&i| encloses[i][j] && depth[i] + 1 == depth[j]);
        if let Some(parent) = parent
            && let Some(slot) = out.iter_mut().find(|(i, _)| *i == parent)
        {
            slot.1.push(j);
        }
    }
    Ok(out)
}

/// Solids from an unordered soup of faces: sew, demand closure, nest
/// shells into solids and voids: the pipeline's own final stages offered
/// as a builder.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// shell fails to close; an open soup encloses no volume, and saying so
/// beats guessing.
pub fn make_volume(model: &mut Model, faces: &[Shape], tol: Tolerances) -> OgeomResult<Built> {
    let sewn = ogeom_algo::sew(model, faces, tol)?;
    for shell in &sewn.shells {
        if !ogeom_algo::is_shell_closed(model, shell)? {
            ogeom_bail!(
                Construction,
                "the faces do not close into shells; an open soup encloses no volume"
            );
        }
    }
    let mut solids = Vec::new();
    for (outer, voids) in nest_shells(model, &sewn.shells, tol)? {
        let mut group = vec![sewn.shells[outer].clone()];
        // A void bounds its solid from inside: material lies outside it,
        // so the sewn outward orientation reverses.
        group.extend(voids.iter().map(|&j| sewn.shells[j].reversed()));
        solids.push(model.add_solid(&group)?);
    }
    let mut history = History::new();
    let result = if solids.len() == 1 {
        solids.remove(0)
    } else {
        model.add_compound(&solids)?
    };
    for face in faces {
        history.modify(face, result.clone());
    }
    Ok(Built::new(result, history))
}

/// The three cells two solids cut space into, each one boolean's answer:
/// what is only in `a`, what is only in `b`, and what is in both. Arbitrary
/// set expressions compose by fusing a selection of these.
#[derive(Debug)]
pub struct Cells {
    /// `a` with `b` removed.
    pub a_not_b: Built,
    /// `b` with `a` removed.
    pub b_not_a: Built,
    /// The overlap.
    pub common: Built,
}

/// Split two solids into their three cells.
///
/// # Errors
///
/// As the operations themselves.
pub fn cells(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Cells> {
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    // One arrangement answers all three: each cell is a different choice
    // of the same classified pieces, the choices `cut` and `common` make,
    // and the cut the other way round with the arguments' roles swapped.
    let fused = general_fuse_classified(model, a, b, None, tol)?;
    let keep = |choose: &dyn Fn(&FacePiece) -> Option<bool>| -> Vec<(usize, bool)> {
        fused
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| choose(p).map(|flip| (i, flip)))
            .collect()
    };
    let cut_from = |first: bool| {
        move |p: &FacePiece| match (p.from_a == first, p.state) {
            (true, PieceState::Out | PieceState::OnOpposed) => Some(false),
            (false, PieceState::In) => Some(true),
            _ => None,
        }
    };
    let a_only = keep(&cut_from(true));
    let b_only = keep(&cut_from(false));
    let both = keep(&|p: &FacePiece| {
        (p.state == PieceState::In || (p.state == PieceState::OnAligned && !p.covered))
            .then_some(false)
    });
    Ok(Cells {
        a_not_b: assemble_result(model, &fused, &a_only, [true, false], a, b, tol)?,
        b_not_a: assemble_result(model, &fused, &b_only, [false, true], a, b, tol)?,
        common: assemble_result(model, &fused, &both, [false, false], a, b, tol)?,
    })
}

/// A tolerance whose confusion *is* the stated fuzz: every gap, pave and
/// weld decision inherits it coherently, which is what a fuzzy boolean
/// means.
fn fuzzed(fuzz: f64, tol: Tolerances) -> OgeomResult<Tolerances> {
    if !fuzz.is_finite() || fuzz <= 0.0 {
        ogeom_bail!(Construction, "a fuzz of {fuzz} is not a distance");
    }
    if fuzz <= tol.confusion() {
        return Ok(tol);
    }
    Tolerances::with_scale(ogeom_core::tolerance::CONFUSION / fuzz)
}

/// [`fuse`] with geometry within `fuzz` of touching counted as touching.
///
/// The fuzz loosens every decision the boolean makes, its welds and probes
/// included, and those reach well past the fuzz itself: at a fuzz near a
/// hundredth of the operands' size they can swallow a face. A result is
/// kept only where its volume agrees, within the fuzz over the operands'
/// area, with what the exact operation finds, or with what any union of
/// the two can hold where the exact one is refused; otherwise the fuzz is
/// taken a tenth smaller, down to the exact operation itself.
///
/// # Errors
///
/// As [`fuse`], plus a non-positive fuzz.
pub fn fuse_fuzzy(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fuzz: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    at_a_fuzz(model, a, b, fuzz, tol, fuse, |va, vb| (va.max(vb), va + vb))
}

/// [`cut`] at a stated fuzz.
///
/// # Errors
///
/// As [`fuse_fuzzy`].
pub fn cut_fuzzy(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fuzz: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    at_a_fuzz(model, a, b, fuzz, tol, cut, |va, vb| {
        ((va - vb).max(0.0), va)
    })
}

/// [`common`] at a stated fuzz.
///
/// # Errors
///
/// As [`fuse_fuzzy`].
pub fn common_fuzzy(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fuzz: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    at_a_fuzz(model, a, b, fuzz, tol, common, |va, vb| (0.0, va.min(vb)))
}

/// One boolean at a fuzz, as [`fuse_fuzzy`] says: the largest fuzz up to
/// `fuzz`, a tenth smaller at a time, whose result holds the volume the
/// exact `operation` finds, or lies within `span` of the operands'
/// volumes where the exact one is refused.
fn at_a_fuzz(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fuzz: f64,
    tol: Tolerances,
    operation: fn(&mut Model, &Shape, &Shape, Tolerances) -> OgeomResult<Built>,
    span: fn(f64, f64) -> (f64, f64),
) -> OgeomResult<Built> {
    let loosened = fuzzed(fuzz, tol)?;
    if loosened.confusion() <= tol.confusion() {
        return operation(model, a, b, tol);
    }
    let measure = |model: &Model, shape: &Shape| -> Option<(f64, f64)> {
        let found =
            ogeom_algo::volume_properties(model, shape, ogeom_mesh::Deflection::default(), tol)
                .ok()?;
        Some((found.mass, found.deflection))
    };
    // Where the operands cannot be measured there is nothing to hold a
    // result to, and the fuzz is taken as asked.
    let (Some((va, da)), Some((vb, db))) = (measure(model, a), measure(model, b)) else {
        return operation(model, a, b, loosened);
    };
    let area = box_area(model, a, tol) + box_area(model, b, tol);
    let exact = operation(model, a, b, tol);
    let held = match &exact {
        Ok(built) => measure(model, &built.shape).map(|(v, _)| (v, v)),
        Err(_) => Some(span(va, vb)),
    };
    let Some((low, high)) = held else {
        return exact;
    };
    let mut first_error = None;
    let mut f = fuzz;
    while f > tol.confusion() {
        let at = fuzzed(f, tol)?;
        match operation(model, a, b, at) {
            Ok(built) => {
                // A meshed measure stands off by its chord over the area.
                let slack = (f + da.max(db)) * area + (va + vb) * 1e-9;
                if let Some((v, dv)) = measure(model, &built.shape)
                    && v >= low - slack - dv * area
                    && v <= high + slack + dv * area
                {
                    return Ok(built);
                }
            }
            Err(e) => {
                first_error.get_or_insert(e);
            }
        }
        f /= 10.0;
    }
    match (exact, first_error) {
        (Ok(built), _) => Ok(built),
        (Err(_), Some(e)) => Err(e),
        (Err(e), None) => Err(e),
    }
}

/// The area of a shape's bounding box: what the fuzz moves its faces over,
/// at most, for a shape that is not a sponge.
fn box_area(model: &Model, shape: &Shape, tol: Tolerances) -> f64 {
    let Ok(bound) = ogeom_algo::tight_bounds(model, shape, tol) else {
        return 0.0;
    };
    let (Some(lo), Some(hi)) = (bound.low(), bound.high()) else {
        return 0.0;
    };
    let (x, y, z) = (hi.x - lo.x, hi.y - lo.y, hi.z - lo.z);
    2.0 * (x * y + y * z + z * x)
}

/// A shape repeated `count` times along a direction at a period and fused
/// into one: the periodic pattern as a composition of what exists.
///
/// # Errors
///
/// As [`fuse`], plus an unusable count or period.
pub fn make_periodic(
    model: &mut Model,
    shape: &Shape,
    step: ogeom_math::Vector,
    count: usize,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if count == 0 {
        ogeom_bail!(Construction, "a pattern of zero copies is nothing");
    }
    if step.magnitude() <= tol.confusion() {
        ogeom_bail!(Construction, "a zero period stacks every copy on the first");
    }
    // Each face and edge of the input traces to its own image in the
    // pattern and generates its copies' images, through every fuse.
    let mut parts: Vec<Shape> = Vec::new();
    for kind in [ShapeType::Face, ShapeType::Edge] {
        parts.extend(ogeom_topo::explore_unique(model, shape, kind)?);
    }
    let mut history = History::new();
    // The pattern so far as disjoint solids: a copy is fused into each piece
    // it joins, and stands as a piece of its own where a gap parts it from
    // the rest.
    let mut pieces: Vec<Shape> = vec![shape.clone()];
    for i in 1..count {
        #[allow(clippy::cast_precision_loss, reason = "pattern counts are small")]
        let offset = step * i as f64;
        let moved =
            ogeom_algo::transformed(model, shape, ogeom_math::Transform::translation(offset))?
                .shape;
        // A placement moves every part with the whole: the copy's parts are
        // the input's, placed as the copy is.
        let placement = moved.location().then(&shape.location().inverted());
        let mut images: Vec<(usize, Vec<Shape>)> = parts
            .iter()
            .enumerate()
            .map(|(k, part)| (k, vec![part.moved(&placement)]))
            .collect();
        let mut current = moved;
        let mut apart: Vec<Shape> = Vec::new();
        for piece in pieces {
            let joined = fuse(model, &piece, &current, tol)?;
            if model.kind_of(&joined.shape)? != ShapeType::Solid {
                apart.push(piece);
                continue;
            }
            history = history.then(&joined.history);
            for (_, traced) in &mut images {
                *traced = traced
                    .iter()
                    .flat_map(|image| joined.history.trace(image).to_vec())
                    .collect();
            }
            current = joined.shape;
        }
        for (k, traced) in images {
            for image in traced {
                history.generate(&parts[k], image);
            }
        }
        apart.push(current);
        pieces = apart;
    }
    let result = if let [one] = pieces.as_slice() {
        one.clone()
    } else {
        ogeom_algo::build::make_compound(model, &pieces)?.shape
    };
    history.generate(shape, result.clone());
    Ok(Built::new(result, history))
}

// --- the operations ----------------------------------------------------------

/// A shape whose placements carry scale, rebuilt with the scale baked
/// into its geometry before the pipeline runs.
///
/// A scale changes a surface's parameterization out from under its
/// pcurves (the refusal the gather makes), so the bake is the whole-shape
/// conversion:
/// surfaces restated in world space, edges moved exactly, pcurves
/// re-derived against the new parameterizations. Unscaled shapes pass
/// through untouched.
fn baked_if_scaled(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Shape> {
    let mut restate = false;
    for face in ogeom_topo::explore(model, shape, Filter::OfType(ShapeType::Face))? {
        if face.location().is_identity() {
            continue;
        }
        let placement = face.transform(model.datums())?;
        // A scale changes lengths the melt compares. A reflection flips
        // every chart's natural normal against its face's flag. Either way
        // the operand is restated in world coordinates first, where both
        // effects are already folded in.
        if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 || !placement.preserves_handedness()
        {
            restate = true;
            break;
        }
    }
    if !restate {
        return Ok(shape.clone());
    }
    Ok(ogeom_algo::baked_shape(model, shape, tol)?.shape)
}

/// Whether a shape is a half space: one shell of one face, open by
/// construction or closed inside out.
fn is_half_space(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    Ok(half_space::half_space_face(model, shape, tol)?.is_some())
}

/// A half space resolved into the solid the operation can act on.
///
/// A planar boundary becomes a box filling the material side, sized past
/// the other argument's whole reach. The box's plane-side face is
/// *coplanar with the boundary itself*, so the cut the caller sees is the
/// exact plane. Its far faces stand outside everything the other shape
/// reaches and never appear in the result. A curved boundary is resolved
/// by [`half_space::resolved`]. A shape that is not a half space passes
/// through untouched.
fn resolved_half_space(
    model: &mut Model,
    shape: &Shape,
    other: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let Some(face) = half_space::half_space_face(model, shape, tol)? else {
        return Ok(shape.clone());
    };
    let planar = model
        .node(&face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
    if !planar {
        return half_space::resolved(model, &face, other, tol);
    }
    // The boundary's outward normal points away from the material.
    let (at, outward) = ogeom_algo::face_normal(model, &face, tol)?;
    let bound = ogeom_algo::shape_bounds(model, other, tol)?;
    let Some(centre) = bound.centre() else {
        ogeom_bail!(
            Construction,
            "the other argument has no bound to fill against"
        );
    };
    let reach = bound.diagonal().max(tol.confusion() * 1e3) * 2.0;
    let into = ogeom_math::Direction::new(-outward, tol)?;
    let foot = centre - outward * outward.dot(centre - at);
    let seed = if into.vector().x.abs() < 0.9 {
        ogeom_math::Vector::new(1.0, 0.0, 0.0)
    } else {
        ogeom_math::Vector::new(0.0, 1.0, 0.0)
    };
    let frame_x = ogeom_math::Direction::from_cross(into.vector(), seed, tol)?;
    let oriented = ogeom_math::Frame::new(foot, into, frame_x, tol)?;
    let corner =
        foot - oriented.x().vector() * (reach / 2.0) - oriented.y().vector() * (reach / 2.0);
    let placed = ogeom_math::Frame::new(corner, into, frame_x, tol)?;
    // Deep enough to pass the other shape's far side however far the
    // plane stands from it: a stand-in whose far face fell inside the
    // shape would cut it there, where the half space cuts nothing.
    let depth = reach + (centre - at).dot(into.vector()).max(0.0);
    Ok(ogeom_algo::make_box(model, placed, (reach, reach, depth), tol)?.shape)
}

/// The plane bounding a half space, where `shape` is one and its boundary
/// is planar.
fn half_space_plane(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_math::Plane>> {
    let Some(face) = half_space::half_space_face(model, shape, tol)? else {
        return Ok(None);
    };
    let Some(surface) = model
        .node(&face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .cloned()
    else {
        return Ok(None);
    };
    use ogeom_geom::Transformable as _;
    match surface.transformed(&face.transform(model.datums())?, tol)? {
        SurfaceGeometry::Plane(p) => Ok(Some(p.plane())),
        _ => Ok(None),
    }
}

/// The union of two solids.
///
/// # Errors
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) for configurations the
/// boolean refuses: tangential or same-domain contact, scaled placements;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for arguments
/// that are not closed solids.
pub fn fuse(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if is_sheet(model, a)? || is_sheet(model, b)? {
        ogeom_bail!(
            Construction,
            "a union takes solids; common and cut trim a sheet, and split_sheet \
             splits one"
        );
    }
    let built = noting_held(model, |model| {
        ball_chart::with_poles_clear(model, a, b, tol, |model, a, b| {
            settled(model, tol, |model, tol| fuse_once(model, a, b, tol))
        })
    })?;
    sub_shapes::with_sub_shapes(model, [a, b], built, tol)
}

/// What `run` makes, with the nodes it passed through from the shapes it
/// was given noted as held by those shapes as well
/// ([`Model::note_held`]): an operation editing the result in place
/// copies them first ([`Model::unshare`]).
fn noting_held(
    model: &mut Model,
    run: impl FnOnce(&mut Model) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    let since = model.node_count();
    let built = run(model)?;
    model.note_held(&built.shape, since)?;
    Ok(built)
}

fn fuse_once(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if is_half_space(model, a, tol)? || is_half_space(model, b, tol)? {
        ogeom_bail!(
            Construction,
            "the union with a half space is unbounded; a half space serves cut, \
             common and section"
        );
    }
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    // Outward pieces bound the union. A same-domain pair with aligned
    // material keeps one copy, and one with opposed material is interior to
    // the union and vanishes.
    let built = fused_and_assembled(model, a, b, [true, true], tol, |p| {
        (p.state == PieceState::Out || (p.state == PieceState::OnAligned && !p.covered))
            .then_some(false)
    });
    // A union of two solids is never empty.
    let built = match built {
        Ok(built) if holds_nothing(model, &built.shape)? => {
            Err(ogeom_core::ogeom_err!(NotDone, "the union came out empty"))
        }
        other => other,
    };
    // The union of a solid and one lying within it is the outer one,
    // whatever contact their boundaries make.
    or_nested(model, built, a, b, true, tol)
}

/// Which of two solids lies within the other, where one does: every sample
/// of the inner one's boundary in or on the outer, none of the outer's
/// strictly inside the inner. `Some(true)` for `b` within `a`, `Some(false)`
/// for `a` within `b`.
///
/// Asked only once the general boolean has refused: boundaries that run
/// along each other within their tolerances (a pad sketched on a converted
/// part's own outline, its sides chords of the part's walls) defeat the
/// arrangement, and a solid within the other is what they usually are.
fn nested(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Option<bool>> {
    sampled_relation(model, a, b, Relation::Nested, tol)
}

/// Whether two solids share no volume: no sample of either's boundary
/// strictly inside the other. Boxes apart settle it without sampling.
fn apart(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    let (box_a, box_b) = (
        ogeom_algo::shape_bounds(model, a, tol)?,
        ogeom_algo::shape_bounds(model, b, tol)?,
    );
    if !box_a.expanded(tol.confusion()).intersects(&box_b) {
        return Ok(true);
    }
    Ok(sampled_relation(model, a, b, Relation::Apart, tol)? == Some(true))
}

/// What [`sampled_relation`] reads of two boundaries.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Relation {
    /// Which solid lies within the other, as [`nested`] answers.
    Nested,
    /// Whether neither boundary enters the other solid: `Some(true)` when
    /// none does.
    Apart,
    /// Whether the first solid's boundary stands nowhere outside the
    /// second: `Some(true)` when no sample of it reads out.
    Covered,
}

/// Two solids' boundaries sampled and read against each other.
fn sampled_relation(
    model: &Model,
    a: &Shape,
    b: &Shape,
    relation: Relation,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    // The boundary's samples, each on its face's exact surface: the face
    // mesh's vertices, and its triangles' middles taken in the chart and
    // lifted, since a middle taken in space stands off a curved face by the
    // mesh's own sag. A triangle along a trim bending into the face reaches
    // past it by its chord's sag, and so may its middle: a middle is kept
    // only where it reads inside the trim drawn at `fine`.
    let samples = |shape: &Shape, fine: f64| -> OgeomResult<Vec<Point>> {
        let mut out = Vec::new();
        for face in explore_unique(model, shape, ShapeType::Face)? {
            let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
                continue;
            };
            let Some(surface) = model.geometry().surface(data.surface) else {
                ogeom_bail!(Dangling, "surface is not in this model");
            };
            let placement = face.transform(model.datums())?;
            let mesh =
                ogeom_mesh::triangulate_face(model, &face, ogeom_mesh::Deflection::default(), tol)?;
            out.extend(mesh.positions.iter().copied());
            if !mesh.has_parameters() {
                continue;
            }
            let trim = ogeom_algo::classify::FaceTrim::of(
                model,
                &face,
                ogeom_mesh::Deflection {
                    chord: fine,
                    ..ogeom_mesh::Deflection::default()
                },
                tol,
            )?;
            for [p, q, r] in (0..mesh.triangles.len()).filter_map(|i| mesh.triangle_parameters(i)) {
                let (u, v) = ((p.0 + q.0 + r.0) / 3.0, (p.1 + q.1 + r.1) / 3.0);
                if trim.holds(model, (u, v), 0.0, tol)? == Containment::In {
                    out.push(placement.apply(surface.point_at(u, v, tol)?));
                }
            }
        }
        Ok(out)
    };
    // Every sample on the right side: `allowed` names the containments that
    // fit. An ambiguous sample says nothing. A boundary most of whose
    // samples say nothing is not settled.
    let all_read = |points: &[Point],
                    boundary: &ogeom_algo::SolidBoundary,
                    allowed: &dyn Fn(Containment) -> bool|
     -> OgeomResult<bool> {
        let mut unread = 0_usize;
        for p in points {
            ogeom_core::progress::checkpoint()?;
            match boundary.holds(model, *p, tol) {
                Ok(says) if allowed(says) => {}
                Ok(_) => return Ok(false),
                Err(ogeom_core::OgeomError::NotDone(_)) => unread += 1,
                Err(e) => return Err(e),
            }
        }
        Ok(unread * 10 <= points.len())
    };
    // Read at the looseness the solids state: a converted part's walls are
    // fits standing off its mesh's points, and a solid drawn on those
    // points lies on the walls only to that.
    let mut loosest = tol.confusion();
    for shape in [a, b] {
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in explore_unique(model, shape, kind)? {
                let stated = match model.node(&sub).map(|n| n.data()) {
                    Some(NodeData::Edge(d)) => d.tolerance.get(),
                    Some(NodeData::Vertex(d)) => d.tolerance.get(),
                    _ => continue,
                };
                loosest = loosest.max(stated);
            }
        }
    }
    // The rings are polylined at the band the caller's tolerance sets, or
    // at the solids' own looseness where that is coarser. Taken from the
    // loosened tolerance instead, a few microns of slop would polyline the
    // rings at tens of millimetres, and every ray would land within that of
    // some trim and read nothing.
    let band = (tol.confusion() * 1e4).max(loosest);
    let tol = if loosest > tol.confusion() {
        fuzzed(loosest, tol)?
    } else {
        tol
    };
    let (of_a, of_b) = (
        ogeom_algo::SolidBoundary::of(model, a, band, tol)?,
        ogeom_algo::SolidBoundary::of(model, b, band, tol)?,
    );
    let (on_a, on_b) = (samples(a, band)?, samples(b, band)?);
    if relation == Relation::Apart {
        let clear = all_read(&on_a, &of_b, &|c| c != Containment::In)?
            && all_read(&on_b, &of_a, &|c| c != Containment::In)?;
        return Ok(Some(clear));
    }
    if relation == Relation::Covered {
        return Ok(Some(all_read(&on_a, &of_b, &|c| c != Containment::Out)?));
    }
    let within = |inner: &[Point],
                  outer_boundary: &ogeom_algo::SolidBoundary,
                  outer: &[Point],
                  inner_boundary: &ogeom_algo::SolidBoundary|
     -> OgeomResult<bool> {
        Ok(all_read(inner, outer_boundary, &|c| c != Containment::Out)?
            && all_read(outer, inner_boundary, &|c| c != Containment::In)?)
    };
    if within(&on_b, &of_a, &on_a, &of_b)? && !pokes_out(model, b, a, &of_a, tol)? {
        return Ok(Some(true));
    }
    if within(&on_a, &of_b, &on_b, &of_a)? && !pokes_out(model, a, b, &of_b, tol)? {
        return Ok(Some(false));
    }
    Ok(None)
}

/// Whether `inner`'s boundary pokes out through `outer`'s between the
/// samples [`nested`] reads.
///
/// A drum a few microns wider than the box it stands in bulges out of each
/// wall along a strip a fraction of a millimetre wide, and samples a mesh's
/// chord apart step over it. So each inner face climbs its own chart from
/// the samples standing nearest each outer face, toward that face's outside,
/// and the furthest point reached is asked of the outer boundary: out there,
/// and on the inner face, the inner solid is not within the outer.
fn pokes_out(
    model: &Model,
    inner: &Shape,
    outer: &Shape,
    outer_boundary: &ogeom_algo::SolidBoundary,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;
    use ogeom_geom::Transformable as _;
    // A start this close to a wall, inside it, may sit beside a bulge the
    // samples straddle: a chord's sag, a few times over.
    let deflection = ogeom_mesh::Deflection::default();
    let reach = deflection.chord * 4.0;
    // The starts climbed per inner face and wall: the nearest few, which
    // stand beside the bulges if there are any.
    const STARTS: usize = 6;
    let placed = |face: &Shape| -> OgeomResult<SurfaceGeometry> {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "surface is not in this model");
        };
        surface.transformed(&face.transform(model.datums())?, tol)
    };
    // Each outer face: its surface placed, which way its outside lies, and
    // the box a near start stands in.
    let mut walls = Vec::new();
    for face in explore_unique(model, outer, ShapeType::Face)? {
        let outward = if face.orientation() == ogeom_topo::Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        let bound = ogeom_algo::shape_bounds(model, &face, tol)?.expanded(reach);
        walls.push((placed(&face)?, outward, bound));
    }
    // How far a point stands outside a wall's surface, along its normal.
    let outside = |p: Point, wall: &(SurfaceGeometry, f64, ogeom_math::Aabb)| -> Option<f64> {
        let found = ogeom_algo::project_on_surface(&wall.0, p, 8, tol).ok()?;
        let (u, v) = found.parameters;
        let n = wall.0.normal_at(u, v, tol).ok()?;
        Some((p - found.point).dot(n.vector()) * wall.1)
    };
    for face in explore_unique(model, inner, ShapeType::Face)? {
        ogeom_core::progress::checkpoint()?;
        let surface = placed(&face)?;
        let mesh = ogeom_mesh::triangulate_face(model, &face, deflection, tol)?;
        if !mesh.has_parameters() {
            continue;
        }
        let corners: Vec<[(f64, f64); 3]> = (0..mesh.triangles.len())
            .filter_map(|i| mesh.triangle_parameters(i))
            .collect();
        let Some(&[first, ..]) = corners.first() else {
            continue;
        };
        let (mut lo, mut hi) = (first, first);
        for &(u, v) in corners.iter().flatten() {
            lo = (lo.0.min(u), lo.1.min(v));
            hi = (hi.0.max(u), hi.1.max(v));
        }
        let span = (hi.0 - lo.0).max(hi.1 - lo.1);
        if span <= 0.0 {
            continue;
        }
        let lift = |at: (f64, f64)| surface.point_at(at.0, at.1, tol);
        let mut starts: Vec<(f64, f64)> = mesh.vertex_parameters().into_iter().flatten().collect();
        for &[p, q, r] in &corners {
            starts.push(((p.0 + q.0 + r.0) / 3.0, (p.1 + q.1 + r.1) / 3.0));
        }
        for wall in &walls {
            let mut near: Vec<((f64, f64), f64)> = Vec::new();
            for &at in &starts {
                let p = lift(at)?;
                if !wall.2.contains(p) {
                    continue;
                }
                if let Some(d) = outside(p, wall)
                    && d >= -reach
                {
                    near.push((at, d));
                }
            }
            near.sort_by(|x, y| y.1.total_cmp(&x.1));
            near.truncate(STARTS);
            for (start, d) in near {
                // A pattern search in the chart, held to the face's own
                // parameter box: a step each way along u and v, halved
                // where none climbs.
                let (mut at, mut best, mut step) = (start, d, span * 0.05);
                let mut rounds = 0;
                while step > span * 1e-7 && rounds < 400 {
                    rounds += 1;
                    let mut climbed = false;
                    for (du, dv) in [(step, 0.0), (-step, 0.0), (0.0, step), (0.0, -step)] {
                        let next = ((at.0 + du).clamp(lo.0, hi.0), (at.1 + dv).clamp(lo.1, hi.1));
                        if let Ok(q) = lift(next)
                            && let Some(g) = outside(q, wall)
                            && g > best
                        {
                            (at, best) = (next, g);
                            climbed = true;
                        }
                    }
                    if !climbed {
                        step *= 0.5;
                    }
                }
                if best <= tol.confusion() {
                    continue;
                }
                let q = lift(at)?;
                if ogeom_algo::classify_on_face(model, &face, q, deflection, tol)?
                    != Containment::Out
                    && outer_boundary
                        .holds(model, q, tol)
                        .is_ok_and(|c| c == Containment::Out)
                {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// A boolean's result, or where the general boolean refused and one solid
/// lies within the other, the answer that nesting gives: the outer one for
/// a union (`outer` true), the inner one for an intersection. The solid
/// returned is that argument itself, the other's faces deleted.
fn or_nested(
    model: &mut Model,
    built: OgeomResult<Built>,
    a: &Shape,
    b: &Shape,
    outer: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let refusal = match built {
        Ok(built) => return Ok(built),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Err(e) => e,
    };
    let Some(b_in_a) = nested(model, a, b, tol)? else {
        return Err(refusal);
    };
    // For the union the container stands. For the intersection, what it
    // contains.
    let (kept, gone) = if b_in_a == outer { (a, b) } else { (b, a) };
    let mut history = ogeom_algo::History::new();
    for face in explore_unique(model, gone, ShapeType::Face)? {
        history.delete(&face);
    }
    Ok(Built::new(kept.clone(), history))
}

/// The intersection of two solids.
///
/// Either argument may instead be a sheet (a face, a shell or a compound of
/// those) when the other is a solid or a half space: the result is then the
/// part of the sheet inside the other. Each face of the sheet is split where
/// the other's boundary crosses it, and the pieces inside are sewn along the
/// edges they share into one shell, or a compound of shells where they fall
/// apart, or an empty compound where none is inside. Each piece lies on its
/// face's surface, and the history records each against the face it came
/// from and every face nothing is kept of as deleted.
///
/// # Errors
///
/// As [`fuse`].
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for two
/// sheets, which share no volume.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a piece of
/// a sheet lies on the other's boundary, which keeps it on neither side.
pub fn common(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    match (is_sheet(model, a)?, is_sheet(model, b)?) {
        (true, true) => ogeom_bail!(
            Construction,
            "two sheets share no volume to keep; split_sheet splits one by the other"
        ),
        (true, false) => trimmed_sheet(model, a, b, true, tol),
        (false, true) => trimmed_sheet(model, b, a, true, tol),
        (false, false) => {
            let built = noting_held(model, |model| {
                ball_chart::with_poles_clear(model, a, b, tol, |model, a, b| {
                    settled(model, tol, |model, tol| common_once(model, a, b, tol))
                })
            })?;
            sub_shapes::with_sub_shapes(model, [a, b], built, tol)
        }
    }
}

fn common_once(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    // Inward pieces bound the intersection. An aligned same-domain pair
    // bounds it too, once. An opposed pair encloses no volume between them.
    let built = fused_and_assembled(model, a, b, [false, false], tol, |p| {
        (p.state == PieceState::In || (p.state == PieceState::OnAligned && !p.covered))
            .then_some(false)
    });
    // Nothing in common is an answer only where the two share no volume.
    let built = match built {
        Ok(built) if holds_nothing(model, &built.shape)? && !apart(model, a, b, tol)? => {
            Err(ogeom_core::ogeom_err!(
                NotDone,
                "the intersection came out empty though the solids overlap"
            ))
        }
        other => other,
    };
    // The intersection of a solid and one lying within it is the inner one.
    or_nested(model, built, a, b, false, tol)
}

/// Two solids fused generally and the pieces `keep` takes assembled, each
/// with whether it is turned over; `keep_aside` says, per solid, whether
/// its faces outside the other are kept.
///
/// The faces each solid has clear of the other are set aside first and
/// come through untouched. Where that cannot be assembled (an edge the
/// faces set aside share would have to move), every face is taken.
fn fused_and_assembled(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    keep_aside: [bool; 2],
    tol: Tolerances,
    keep: impl Fn(&FacePiece) -> Option<bool>,
) -> OgeomResult<Built> {
    let attempt = |model: &mut Model, aside: Option<&[aside::Aside; 2]>| {
        let fused = general_fuse_classified(model, a, b, aside, tol)?;
        let kept: Vec<(usize, bool)> = fused
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| keep(p).map(|flip| (i, flip)))
            .collect();
        assemble_result(model, &fused, &kept, keep_aside, a, b, tol)
    };
    if let Ok(Some(aside)) = aside::set_aside(model, a, b, tol) {
        // The faces set aside keep the solids' own edges and vertices, and
        // what the attempt widens of those is put back where it fails: the
        // solids are given back as they came, and taking every face starts
        // from them.
        let mark = model.note_widened();
        match attempt(model, Some(&aside)) {
            Ok(built) => {
                model.keep_widened(mark);
                return Ok(built);
            }
            Err(e) => {
                model.undo_widened(mark);
                if matches!(e, ogeom_core::OgeomError::Cancelled) {
                    return Err(e);
                }
            }
        }
    }
    attempt(model, None)
}

/// Whether a boolean's result has no face: an empty compound, or one of
/// empty shells.
fn holds_nothing(model: &Model, shape: &Shape) -> OgeomResult<bool> {
    // Down from the shape until the first face, which settles it.
    let mut stack = vec![shape.clone()];
    while let Some(at) = stack.pop() {
        let Some(node) = model.node(&at) else {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        };
        match node.kind() {
            ShapeType::Face => return Ok(false),
            ShapeType::Compound | ShapeType::CompSolid | ShapeType::Solid | ShapeType::Shell => {
                stack.extend(node.children().iter().cloned());
            }
            _ => {}
        }
    }
    Ok(true)
}

/// The first solid with the second removed.
///
/// The pieces of `b` that close the cut into `a`'s material bound the removed
/// volume from `b`'s side, so they join the result with their material side
/// flipped.
///
/// `a` may instead be a sheet (a face, a shell or a compound of those) when
/// `b` is a solid or a half space: the result is then the part of the sheet
/// outside `b`, made as [`common`] makes the part inside.
///
/// # Errors
///
/// As [`common`].
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// `b` is a sheet, which bounds no volume to remove.
pub fn cut(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if is_sheet(model, b)? {
        ogeom_bail!(
            Construction,
            "a sheet bounds no volume to cut away; split_sheet splits by a sheet"
        );
    }
    if is_sheet(model, a)? {
        return trimmed_sheet(model, a, b, false, tol);
    }
    let built = noting_held(model, |model| {
        ball_chart::with_poles_clear(model, a, b, tol, |model, a, b| {
            settled(model, tol, |model, tol| cut_once(model, a, b, tol))
        })
    })?;
    sub_shapes::with_sub_shapes(model, [a, b], built, tol)
}

/// [`cut`], where `b` is a blend's wedge, or a compound of wedges standing
/// apart, that sits flush on `a`: each wedge bounded by planar seats lying
/// on faces of `a` and one strip of a cylinder, the blend; or each a rim's
/// wedge, an annulus of a plane and a band of a cylinder about one axis
/// lying on faces of `a`, and the quarter of a torus's tube between them;
/// or a piece of a rim's wedge between two planar ends through its axis,
/// at most half a turn, each seat and end lying on a face of `a`.
/// Where every seat lies within the face it is on, facing its way, and no
/// other face of `a` reaches into the wedge, the faces the seats lie on are
/// trimmed along the wedge's own edges and the blends close what is left:
/// the general fuse is not run, and the result and its history are those
/// it gives. A straight band's wedge whose end runs, at a block's corner,
/// into a blend of `a` of its radius rounding the face it rounds is trimmed
/// so as its part inside `a`, cut back at the mitre where the two blends
/// cross. Where `b` mixes rims' wedges with other lumps, each rim's box
/// standing further from every other lump's than the boolean welds, the
/// rims are trimmed out so and the other lumps cut from what is left, to
/// the solid and history the general fuse gives cutting all of `b`.
/// Otherwise, and where `a` or `b` is anything else, this is [`cut`].
///
/// # Errors
///
/// As [`cut`].
pub fn cut_seated(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if is_sheet(model, a)? || is_sheet(model, b)? {
        return cut(model, a, b, tol);
    }
    let built = noting_held(model, |model| {
        ball_chart::with_poles_clear(model, a, b, tol, |model, a, b| {
            match seated::seated_cut(model, a, b, false, tol) {
                Ok(Some(built)) => return Ok(built),
                Err(
                    e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_)),
                ) => {
                    return Err(e);
                }
                Ok(None) | Err(_) => {}
            }
            let rest = |model: &mut Model, a: &Shape, b: &Shape| {
                settled(model, tol, |model, tol| cut_once(model, a, b, tol))
            };
            match seated::rings_first(model, a, b, false, tol, rest) {
                Ok(Some(built)) => return Ok(built),
                Err(
                    e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_)),
                ) => {
                    return Err(e);
                }
                Ok(None) | Err(_) => {}
            }
            settled(model, tol, |model, tol| cut_once(model, a, b, tol))
        })
    })?;
    sub_shapes::with_sub_shapes(model, [a, b], built, tol)
}

/// [`fuse`], where `b` is a concave rim's blend wedge, or a compound of
/// such wedges standing apart, that sits flush on `a`: each an annulus of a
/// plane and a band of a cylinder about one axis lying on faces of `a`,
/// facing against them, and the quarter of a torus's tube between them; or
/// a concave straight band's wedge or a piece of a rim's, as [`cut_seated`]
/// reads them, its seats lying on faces of `a` facing against them.
/// Where every seat lies within the face it is on and no other face of `a`
/// reaches into the region the wedge lies in, the faces the seats lie on are
/// trimmed along the wedge's own edges and the blends close what is left:
/// the general fuse is not run, and the result and its history are those
/// it gives. Where `b` mixes such wedges with other lumps, each rim's box
/// standing further from every other lump's than the boolean welds, the
/// rims are filled so and the other lumps fused on to what is left, to the
/// solid and history the general fuse gives fusing all of `b`. Otherwise,
/// and where `a` or `b` is anything else, this is [`fuse`].
///
/// # Errors
///
/// As [`fuse`].
pub fn fuse_seated(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if is_sheet(model, a)? || is_sheet(model, b)? {
        return fuse(model, a, b, tol);
    }
    let built = noting_held(model, |model| {
        ball_chart::with_poles_clear(model, a, b, tol, |model, a, b| {
            match seated::seated_cut(model, a, b, true, tol) {
                Ok(Some(built)) => return Ok(built),
                Err(
                    e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_)),
                ) => {
                    return Err(e);
                }
                Ok(None) | Err(_) => {}
            }
            let rest = |model: &mut Model, a: &Shape, b: &Shape| {
                settled(model, tol, |model, tol| fuse_once(model, a, b, tol))
            };
            match seated::rings_first(model, a, b, true, tol, rest) {
                Ok(Some(built)) => return Ok(built),
                Err(
                    e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_)),
                ) => {
                    return Err(e);
                }
                Ok(None) | Err(_) => {}
            }
            settled(model, tol, |model, tol| fuse_once(model, a, b, tol))
        })
    })?;
    sub_shapes::with_sub_shapes(model, [a, b], built, tol)
}

/// What of a sheet lies inside `tool` (`inside` true) or outside it, as a
/// sheet: [`common`] and [`cut`] with a sheet.
fn trimmed_sheet(
    model: &mut Model,
    sheet: &Shape,
    tool: &Shape,
    inside: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !is_sheet(model, sheet)? {
        ogeom_bail!(
            Construction,
            "the sheet to trim is a face, a shell or a compound of those"
        );
    }
    trim_sheet_edges(model, sheet, tol)?;
    let built = noting_held(model, |model| {
        settled(model, tol, |model, tol| {
            let tool = &resolved_half_space(model, tool, sheet, tol)?;
            let (sheet, tool) = (
                &baked_if_scaled(model, sheet, tol)?,
                &baked_if_scaled(model, tool, tol)?,
            );
            let fused = general_fuse_as(model, sheet, tool, Operands::SheetBySolid, None, tol)?;
            let mut kept = Vec::new();
            for (index, piece) in fused.pieces.iter().enumerate() {
                match piece.state {
                    PieceState::In if inside => kept.push(index),
                    PieceState::Out if !inside => kept.push(index),
                    PieceState::In | PieceState::Out => {}
                    PieceState::OnAligned | PieceState::OnOpposed => ogeom_bail!(
                        NotDone,
                        "a piece of the sheet lies on the tool's boundary, which \
                     keeps it on neither side"
                    ),
                    PieceState::Unread => ogeom_bail!(
                        Construction,
                        "a piece of the sheet was not classified against the tool"
                    ),
                }
            }
            assemble_sheet(model, &fused, &kept, sheet, tol)
        })
    })?;
    sub_shapes::with_sub_shapes(model, [sheet, tool], built, tol)
}

/// A sheet split where another sheet crosses it: the pieces, with history.
///
/// `sheet` and `by` are each a face, a shell or a compound of those. Each
/// face of `sheet` is split along the curves where `by`'s faces cross it,
/// and along `by`'s edges where they lie in it. Only `sheet` is split: two
/// sheets are split by each other in two calls, each sheet by the other.
/// The pieces are sewn along the edges they share into one shell, or a
/// compound of shells where they fall apart. Each piece lies on its face's
/// surface, and the history records each against the face it came from;
/// a face nothing crosses is its own copy. A piece reaching across the seam
/// of a face on a surface of revolution is one face, on that surface turned
/// about its axis so the chart's join falls outside the piece.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// either argument is not a sheet; a solid or a half space trims a sheet
/// through [`common`] and [`cut`].
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where the two
/// sheets overlap on one surface, and as [`fuse`] for configurations the
/// splitting refuses.
pub fn split_sheet(
    model: &mut Model,
    sheet: &Shape,
    by: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for shape in [sheet, by] {
        if !is_sheet(model, shape)? {
            ogeom_bail!(
                Construction,
                "split_sheet splits a sheet by a sheet, each a face, a shell or \
                 a compound of those; common and cut trim a sheet by a solid"
            );
        }
        trim_sheet_edges(model, shape, tol)?;
    }
    let built = noting_held(model, |model| {
        settled(model, tol, |model, tol| {
            let (sheet, by) = (
                &baked_if_scaled(model, sheet, tol)?,
                &baked_if_scaled(model, by, tol)?,
            );
            let fused = general_fuse_as(model, sheet, by, Operands::SheetBySheet, None, tol)?;
            let kept: Vec<usize> = (0..fused.pieces.len()).collect();
            assemble_sheet(model, &fused, &kept, sheet, tol)
        })
    })?;
    sub_shapes::with_sub_shapes(model, [sheet, by], built, tol)
}

fn cut_once(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    // The first argument's outward pieces stay. The tool's inward pieces
    // close the cut with their material side flipped. On the shared
    // surface: an opposed pair means the tool's material is entirely on the
    // other side, so the first argument's face survives untouched; an
    // aligned pair means the tool's material backs the same wall, which the
    // cut removes.
    let built = fused_and_assembled(model, a, b, [true, false], tol, |p| {
        match (p.from_a, p.state) {
            (true, PieceState::Out | PieceState::OnOpposed) => Some(false),
            (false, PieceState::In) => Some(true),
            _ => None,
        }
    });
    let refusal = match built {
        // Nothing left is an answer only where no part of the solid's
        // boundary stands outside the tool.
        Ok(built)
            if !holds_nothing(model, &built.shape)?
                || sampled_relation(model, a, b, Relation::Covered, tol)? == Some(true) =>
        {
            return Ok(built);
        }
        Ok(_) => ogeom_core::ogeom_err!(
            NotDone,
            "the cut left nothing of a solid that does not lie within the tool"
        ),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Err(e) => e,
    };
    // A solid lying within the tool, whatever contact their boundaries
    // make (the same solid twice among them), is cut away entirely.
    if nested(model, b, a, tol)? != Some(true) {
        return Err(refusal);
    }
    let empty = model.add_compound(&[])?;
    let mut history = ogeom_algo::History::new();
    for face in explore_unique(model, a, ShapeType::Face)? {
        history.delete(&face);
    }
    history.modify(a, empty.clone());
    Ok(Built::new(empty, history))
}

/// Refused where `section`'s image on `face`, a sphere face, leaves the
/// chart at a pole: past the pole row by more than the section's
/// tolerance, or through one of the face's poles between its ends. An
/// image past the pole row stands for the point half a turn round on the
/// near side, and a curve through a pole leaves it half a turn round the
/// chart from where it came, which no one image can follow; the
/// arrangement would read either where the section is not. Each is read
/// at its worst among evenly spaced samples, refined by golden section
/// between their neighbours.
fn past_a_pole(
    section: &SectionRec,
    pcurve: &PlanarCurve,
    face: &GFace,
    tol: Tolerances,
) -> OgeomResult<()> {
    const SAMPLES: usize = 256;
    let SurfaceGeometry::Sphere(sphere) = &face.surface else {
        return Ok(());
    };
    let radius = sphere.sphere().radius();
    let reach = section.tolerance.max(tol.confusion());
    let range = section.curve.domain();
    let worst = |measure: &dyn Fn(f64) -> OgeomResult<f64>| -> OgeomResult<f64> {
        #[allow(clippy::cast_precision_loss)]
        let step = |k: usize| range.0 + (range.1 - range.0) * k as f64 / SAMPLES as f64;
        let mut top = (f64::NEG_INFINITY, 0);
        for k in 0..=SAMPLES {
            let m = measure(step(k))?;
            if m > top.0 {
                top = (m, k);
            }
        }
        let (mut lo, mut hi) = (
            step(top.1.saturating_sub(1)),
            step((top.1 + 1).min(SAMPLES)),
        );
        let golden = (5.0_f64.sqrt() - 1.0) / 2.0;
        for _ in 0..100 {
            let (a, b) = (hi - golden * (hi - lo), lo + golden * (hi - lo));
            if measure(a)? > measure(b)? {
                hi = b;
            } else {
                lo = a;
            }
        }
        Ok(measure(f64::midpoint(lo, hi))?.max(top.0))
    };
    let past =
        worst(&|t| Ok((pcurve.point_at(t, tol)?.y.abs() - core::f64::consts::FRAC_PI_2) * radius))?;
    if past > reach {
        ogeom_bail!(
            NotDone,
            "a section's image on a sphere runs {past} past the chart's pole, beyond \
             its tolerance of {reach}; a section through a sphere chart's pole is not \
             resolved; see docs/CAPABILITIES.md, bool.booleans"
        );
    }
    let (start, end) = (
        section.curve.point_at(range.0, tol)?,
        section.curve.point_at(range.1, tol)?,
    );
    for pole in &face.poles {
        if start.distance(pole.point) <= reach || end.distance(pole.point) <= reach {
            continue;
        }
        let near = -worst(&|t| Ok(-section.curve.point_at(t, tol)?.distance(pole.point)))?;
        if near <= reach {
            ogeom_bail!(
                NotDone,
                "a section passes within {near} of a sphere chart's pole between its \
                 ends, inside its tolerance of {reach}; a section through a sphere \
                 chart's pole is not resolved; see docs/CAPABILITIES.md, bool.booleans"
            );
        }
    }
    Ok(())
}

/// The edges where the two solids' boundaries cross.
///
/// Either argument may instead be a face, a shell or a compound of faces
/// with no solid behind it, when the other is a half space bounded by a
/// plane: the section is then [`section_face`](fn@section_face)'s, the curves where its
/// faces cross the plane.
///
/// # Errors
///
/// As [`fuse`].
pub fn section(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    // A face, shell or compound with no solid behind it, against a half
    // space bounded by a plane: the curves where its faces cross the plane.
    for (sheet, other) in [(a, b), (b, a)] {
        if model.kind_of(sheet)? != ShapeType::Solid
            && let Some(plane) = half_space_plane(model, other, tol)?
        {
            return section_face(model, sheet, &plane, tol);
        }
    }
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    let fused = general_fuse(model, a, b, None, tol)?;
    // Every section sub-edge that survived into some piece's ring, built
    // once per distinct sub-range.
    let mut wanted: Vec<(usize, (f64, f64))> = Vec::new();
    for piece in &fused.pieces {
        for ring in &piece.rings {
            for traversal in ring {
                if let Tag::Section { section, range } = &traversal.tag {
                    let near = |a: (f64, f64), b: (f64, f64)| {
                        (a.0 - b.0).abs() <= tol.parametric()
                            && (a.1 - b.1).abs() <= tol.parametric()
                    };
                    if !wanted.iter().any(|(s, r)| s == section && near(*r, *range)) {
                        wanted.push((*section, *range));
                    }
                }
            }
        }
    }
    let mut history = History::new();
    let mut edges = Vec::new();
    // Consecutive pieces of the section meet end to end, and there they
    // share one vertex: the section is a wire's worth of edges, not a
    // scatter of segments each with vertices of its own.
    let mut vertices: Vec<(Point, Shape)> = Vec::new();
    for (si, range) in wanted {
        let s = &fused.sections[si];
        let domain = s.curve.domain();
        let (f0, f1) = folded_range(range, domain, s.closed);
        let from = s.curve.point_at(at_param(f0, domain, s.closed), tol)?;
        let to = s.curve.point_at(at_param(f1, domain, s.closed), tol)?;
        let v0 = section_face::vertex(model, &mut vertices, from, tol)?;
        let v1 = section_face::vertex(model, &mut vertices, to, tol)?;
        edges.push(make_edge_between(model, s.curve.clone(), (f0, f1), &v0, &v1, tol)?.shape);
    }
    // Contacts are not crossings, so no piece's ring carries them and the
    // loop above cannot see them, but a section through a tangency has a
    // curve in it, and this is where it comes from.
    for contact in &fused.tangents {
        for (lo, hi) in contact_intervals(&fused, contact, tol)? {
            let from = contact.curve.point_at(lo, tol)?;
            let to = contact.curve.point_at(hi, tol)?;
            let v0 = section_face::vertex(model, &mut vertices, from, tol)?;
            let v1 = section_face::vertex(model, &mut vertices, to, tol)?;
            edges.push(
                make_edge_between(model, contact.curve.clone(), (lo, hi), &v0, &v1, tol)?.shape,
            );
        }
    }
    let result = model.add_compound(&edges)?;
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Built::new(result, history))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_algo::{check, make_box, make_cylinder, volume_properties};
    use ogeom_math::{Direction, Frame};
    use ogeom_mesh::Deflection;

    const T: Tolerances = Tolerances::millimetres();
    const PI: f64 = core::f64::consts::PI;

    #[test]
    fn a_hug_measures_the_curve_rather_than_its_sampling_chords() {
        let circle = ogeom_math::Circle::new(Frame::WORLD, 100.0, T).unwrap();
        let curve = Curve::Circle(ogeom_geom::CircleCurve::new(circle));
        for range in [(0.0, 2.0 * PI), (-2.0 * PI, 2.0 * PI), (2.0 * PI, 4.0 * PI)] {
            for angle in [
                0.0_f64,
                PI,
                1.5 * PI,
                2.0 * PI,
                2.5 * PI,
                4.0 * PI,
                -PI,
                1.123456,
                -1e-8,
                1e-8,
                2.0 * PI - 1e-8,
                2.0 * PI + 1e-8,
            ] {
                for offset in [-0.123, 0.0, 0.123] {
                    let p = Point::new(
                        (100.0 + offset) * angle.cos(),
                        (100.0 + offset) * angle.sin(),
                        0.0,
                    );
                    let distance = hug_distance(&curve, range, p, T).unwrap();
                    assert!(
                        (distance - offset.abs()).abs() < 1e-8,
                        "{range:?} {angle} {offset}: {distance}"
                    );
                }
            }
        }
        let p = curve.point_at(0.8, T).unwrap();
        let expected = p.distance(curve.point_at(0.7, T).unwrap());
        let distance = hug_distance(&curve, (0.4, 0.7), p, T).unwrap();
        assert!((distance - expected).abs() < 1e-10);
    }

    /// One cylinder, one cone and one torus, each spelled on frames turned
    /// about the axis, slid along it (the cylinder's and the cone's), with
    /// the axis reversed and on
    /// either handedness: every pair's chart map carries a point of the
    /// owner's chart to the target's point at the same place in space.
    /// A frame of the other handedness turns its angle the other way
    /// round, which the map must mirror rather than turn.
    #[test]
    fn a_coaxial_chart_map_lands_on_the_same_point_across_handedness() {
        let frames = |slide: f64| -> Vec<Frame> {
            let mut out = Vec::new();
            for up in [Direction::Z, -Direction::Z] {
                let x = Direction::new(ogeom_math::Vector::new(0.6, 0.8, 0.0), T).unwrap();
                let frame = Frame::new(Point::new(0.0, 0.0, slide), up, x, T).unwrap();
                out.extend([frame, frame.mirrored()]);
            }
            out
        };
        let cylinder = |f: Frame| -> SurfaceGeometry {
            ogeom_geom::CylinderSurface::new(
                ogeom_math::Cylinder::new(f, 2.0, T).unwrap(),
                (-5.0, 5.0),
            )
            .unwrap()
            .into()
        };
        let torus = |f: Frame| -> SurfaceGeometry {
            ogeom_geom::TorusSurface::new(ogeom_math::Torus::new(f, 2.0, 0.5, T).unwrap()).into()
        };
        // One cone, its radius 2 + z / 4 at a height z on the world's axis,
        // spelled from each frame's origin with the half angle signed by
        // which way the frame's axis runs.
        let cone = |f: Frame| -> SurfaceGeometry {
            let up = f.z().vector().z;
            let height = f.origin().z;
            ogeom_geom::ConeSurface::new(
                ogeom_math::Cone::new(f, 0.25_f64.mul_add(height, 2.0), up * 0.25_f64.atan(), T)
                    .unwrap(),
                (-5.0, 5.0),
            )
            .unwrap()
            .into()
        };
        let kinds: [(&dyn Fn(Frame) -> SurfaceGeometry, f64); 3] =
            [(&cylinder, 1.5), (&cone, 1.5), (&torus, 0.0)];
        for (make, slide) in kinds {
            for target_frame in [Frame::WORLD, Frame::WORLD.mirrored()] {
                let target = make(target_frame);
                for owner_frame in frames(slide) {
                    let owner = make(owner_frame);
                    let map = coaxial_chart_map(&owner, &target, 0.0, T).unwrap();
                    for (u, v) in [(0.3, 0.2), (2.0, -1.0), (5.0, 1.1)] {
                        let at = owner.point_at(u, v, T).unwrap();
                        let uv = map.apply(Point2::new(u, v));
                        let tau = core::f64::consts::TAU;
                        let lifted = target
                            .point_at(uv.x.rem_euclid(tau), uv.y, T)
                            .or_else(|_| {
                                target.point_at(uv.x.rem_euclid(tau), uv.y.rem_euclid(tau), T)
                            })
                            .unwrap();
                        // Two closed forms of one point: rounding only.
                        assert!(
                            at.distance(lifted) < 1e-12,
                            "{:?} owner on {:?} target, ({u}, {v}): {at:?} against {lifted:?}",
                            owner_frame.handedness(),
                            target_frame.handedness()
                        );
                    }
                }
            }
        }
    }

    /// An L-shaped prism across `y` whose step's underside lies on
    /// `z = 0` over `1 <= x <= 2`: past its edge at `x = 1` the plane runs
    /// into the leg below (a concave edge), past its edge at `x = 2` out of
    /// the solid (a convex one). A thin piece of a face on that plane,
    /// lying past each edge, reads inside past the first, outside past the
    /// second, and between both it is dust whose strands disagree.
    #[test]
    fn a_thin_piece_reads_the_other_solid_s_wedge_at_its_edges() {
        let mut model = Model::new();
        let profile = [
            (0.0, -1.0),
            (1.0, -1.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
        ]
        .map(|(x, z)| Point::new(x, 0.0, z));
        let wire = ogeom_algo::make_polygon(&mut model, &profile, true, T)
            .unwrap()
            .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            ogeom_math::Plane::new(
                Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap(),
            ),
            (-5.0, 5.0),
            (-5.0, 5.0),
        )
        .unwrap();
        let section = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let step = ogeom_algo::make_prism(
            &mut model,
            &section,
            ogeom_math::Vector::new(0.0, 1.0, 0.0),
            T,
        )
        .unwrap()
        .shape;
        let floor = make_box(
            &mut model,
            Frame::new(Point::new(-1.0, -1.0, -3.0), Direction::Z, Direction::X, T).unwrap(),
            (4.0, 3.0, 3.0),
            T,
        )
        .unwrap()
        .shape;
        let theirs = gather(&model, &step, false, None, T).unwrap();
        let ours = gather(&model, &floor, false, None, T).unwrap();
        let top = ours
            .faces
            .iter()
            .find(|f| {
                matches!(f.surface, SurfaceGeometry::Plane(_))
                    && f.edges
                        .iter()
                        .all(|e| e.curve.point_at(e.crange.0, T).unwrap().z.abs() < 1e-9)
            })
            .unwrap();
        // The step's edge along `y` at `x`, on `z = 0`, as a contact on
        // the floor's top, with the way the ring walks it so the piece lies
        // on the side `x` grows towards when `beyond` is positive.
        let contact = |x: f64, beyond: f64| -> (ContactRec, (f64, f64), bool) {
            let e = theirs
                .faces
                .iter()
                .flat_map(|f| f.edges.iter())
                .find(|e| {
                    [e.crange.0, e.crange.1].iter().all(|&t| {
                        let p = e.curve.point_at(t, T).unwrap();
                        (p.x - x).abs() < 1e-9 && p.z.abs() < 1e-9
                    })
                })
                .unwrap();
            let ends = [e.crange.0, e.crange.1].map(|t| {
                chart_point_within(top, e.curve.point_at(t, T).unwrap(), 1e-6, T).unwrap()
            });
            let line = ogeom_geom::Line2d::segment(ends[0], ends[1], T).unwrap();
            let prange = ogeom_geom::Curve2d::domain(&line);
            let rec = ContactRec {
                curve: (*e.curve).clone(),
                crange: e.crange,
                pcurve: line.into(),
                prange,
                node: e.node,
                tolerance: e.tolerance,
                target_from_a: true,
                target_face: 0,
                bound: ogeom_math::Aabb::EMPTY,
                crossed: false,
            };
            // The piece is on the left of its ring about the chart's
            // normal, which is `+z` or `-z`.
            let t = f64::midpoint(e.crange.0, e.crange.1);
            let along = e.curve.d1_at(t, T).unwrap();
            let uv = ends[0];
            let normal = top.surface.normal_at(uv.x, uv.y, T).unwrap().vector();
            let left = normal.cross(along);
            (rec, e.crange, left.x * beyond < 0.0)
        };
        let thin = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1e-6),
            Point2::new(0.0, 1e-6),
        ];
        let read = |strands: &[(f64, f64)]| -> Thin {
            let mut contacts = Vec::new();
            let mut ring = Vec::new();
            for (k, &(x, beyond)) in strands.iter().enumerate() {
                let (rec, range, reversed) = contact(x, beyond);
                contacts.push(rec);
                ring.push(Traversal {
                    tag: Tag::Contact { contact: k, range },
                    reversed,
                    nodes: None,
                });
            }
            thin_by_strands(
                top,
                &[ring],
                std::slice::from_ref(&thin),
                &contacts,
                &theirs,
                T,
            )
            .unwrap()
        };
        assert!(matches!(read(&[(1.0, -1.0)]), Thin::Reads(PieceState::In)));
        assert!(matches!(read(&[(2.0, 1.0)]), Thin::Reads(PieceState::Out)));
        assert!(matches!(read(&[(1.0, -1.0), (2.0, 1.0)]), Thin::Disagree));
        // On the step's side of either edge the piece lies on its face: a
        // coincidence the partner reads, and the strand says nothing.
        assert!(matches!(read(&[(2.0, -1.0)]), Thin::Undecided));
    }

    /// A quarter of a cylinder wall of radius 0.002, its edges and vertices
    /// widened to 0.003, its seam turned away from the quarter, one of its
    /// rulings restated a period off. Where that ruling meets the arcs, its
    /// ends stand a period from theirs, and that whole turn lifts within
    /// the edges' radii, yet the surface moves along it: the joint is a
    /// jump of branch, not one across an apex. Gathered, every pcurve of
    /// the wall stands where it stands when the ruling is not restated.
    #[test]
    fn a_period_jump_on_a_wall_smaller_than_its_radii_is_folded_back() {
        let wall = |restate: f64| -> (Model, Vec<(EdgeKey, Point2)>) {
            let mut model = Model::new();
            let x = Direction::new(ogeom_math::Vector::new(-1.0, -1.0, 0.0), T).unwrap();
            let whole = make_cylinder(
                &mut model,
                Frame::new(Point::ORIGIN, Direction::Z, x, T).unwrap(),
                0.002,
                0.01,
                T,
            )
            .unwrap()
            .shape;
            let quadrant = make_box(
                &mut model,
                Frame::new(Point::new(0.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
                (1.0, 1.0, 2.0),
                T,
            )
            .unwrap()
            .shape;
            let part = common(&mut model, &whole, &quadrant, T).unwrap().shape;
            let face = explore_unique(&model, &part, ShapeType::Face).unwrap()[0].clone();
            let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
            let mut restated = false;
            for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
                model
                    .widen(&edge, ogeom_core::Tolerance::new(0.003).unwrap())
                    .unwrap();
                let Some(EdgeRepr::PCurve { curve, range, .. }) = model
                    .node(&edge)
                    .unwrap()
                    .data()
                    .as_edge()
                    .unwrap()
                    .pcurve_for(surface, edge.location())
                    .cloned()
                else {
                    continue;
                };
                let stated = model.geometry().pcurve(curve).unwrap().clone();
                let (a, b) = (
                    stated.point_at(range.0, T).unwrap(),
                    stated.point_at(range.1, T).unwrap(),
                );
                if restate == 0.0 || restated || (a.x - b.x).abs() > 1e-9 {
                    continue;
                }
                let moved = stated
                    .transformed(
                        &ogeom_math::Transform2::translation(ogeom_math::Vector2::new(
                            restate, 0.0,
                        )),
                        T,
                    )
                    .unwrap();
                let moved = model.geometry_mut().add_pcurve(moved);
                if let Some(node) = model.node_mut(&edge)
                    && let NodeData::Edge(data) = node.data_mut()
                {
                    for repr in &mut data.representations {
                        if let EdgeRepr::PCurve {
                            curve, surface: s, ..
                        } = repr
                            && *s == surface
                        {
                            *curve = moved;
                        }
                    }
                }
                restated = true;
            }
            assert!(restate == 0.0 || restated);
            let gathered = gather(&model, &part, false, None, T).unwrap();
            let starts = gathered
                .faces
                .iter()
                .find(|f| matches!(f.surface, SurfaceGeometry::Cylinder(_)))
                .unwrap()
                .edges
                .iter()
                .map(|e| (e.node, e.pcurve.point_at(e.prange.0, T).unwrap()))
                .collect();
            (model, starts)
        };
        let (_, want) = wall(0.0);
        for by in [core::f64::consts::TAU, -core::f64::consts::TAU] {
            let (_, got) = wall(by);
            assert_eq!(got.len(), want.len());
            for ((_, g), (_, w)) in got.iter().zip(&want) {
                assert!(
                    g.distance(*w) < 1e-9,
                    "moved {by}: {got:?} against {want:?}"
                );
            }
        }
    }

    /// An edge passing under a drill of radius 23.6 at right angles,
    /// tangent to it a few microns off at the middle of the edge, is
    /// crossed by the sections on the drill either side of the touch
    /// 0.03 apart: the stretch between stays within five microns of the
    /// drill. The two crossings are one junction at the touch, whether the
    /// sections state a fit error of a micron's quarter or of a micron.
    #[test]
    fn a_graze_s_crossings_are_one_junction_however_tightly_fitted() {
        let radius = 23.6;
        let curve = Curve::Line(
            ogeom_geom::LineCurve::segment(
                Point::new(-10.0, 0.0, 0.0),
                Point::new(10.0, 0.0, 0.0),
                T,
            )
            .unwrap(),
        );
        let crange = ogeom_geom::Curve3d::domain(&curve);
        let off_drill = |_: (bool, usize), p: Point| -> OgeomResult<f64> {
            Ok(p.x.hypot(p.z - radius) - radius + 4.6e-6)
        };
        let half = 0.0148;
        for honesty in [2.6e-4, 1e-3] {
            let paves: Vec<Pave> = [-half, half]
                .into_iter()
                .map(|x| Pave {
                    t: x + 10.0,
                    honesty,
                    across: Some((false, 0)),
                    spread: 0.0,
                })
                .collect();
            let clusters =
                cluster_paves(&curve, crange, 1e-7, &paves, Some(&off_drill), &[], T).unwrap();
            assert_eq!(clusters.len(), 1, "honesty {honesty}: {clusters:?}");
            assert!(clusters[0].at.x.abs() < half, "{clusters:?}");
        }
    }

    /// A box whose front top edge, y = 0 at z = 10, is drawn on the top
    /// face `shift` inside where it runs, and states `loose`: the edge's
    /// tolerance owed to the top face alone, as a chord beside a drum owes
    /// its sag to the drum.
    fn box_with_an_edge_drawn_off_its_top(model: &mut Model, shift: f64, loose: f64) -> Shape {
        let solid = make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let top = explore_unique(model, &solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let b = shape_bounds(model, f, T).unwrap();
                b.low().unwrap().z > 9.9
            })
            .unwrap();
        let top_surface = model.node(&top).unwrap().data().as_face().unwrap().surface;
        let edge = explore_unique(model, &top, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let b = shape_bounds(model, e, T).unwrap();
                b.high().unwrap().y.abs() < 1e-6
            })
            .unwrap();
        let (old, range) = match model
            .node(&edge)
            .unwrap()
            .data()
            .as_edge()
            .unwrap()
            .pcurve_on(top_surface)
        {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
            _ => unreachable!("a box's edge has a plain pcurve on its top"),
        };
        let pcurve = model.geometry().pcurve(old).unwrap().clone();
        let (from, to) = (
            pcurve.point_at(range.0, T).unwrap(),
            pcurve.point_at(range.1, T).unwrap(),
        );
        // Square to the edge in the chart, toward the face's middle.
        let along = (to - from) / (to - from).magnitude();
        let mut across = ogeom_math::Vector2::new(-along.y, along.x);
        let middle = Point2::new(f64::midpoint(from.x, to.x), f64::midpoint(from.y, to.y));
        let lift = |p: Point2| {
            let SurfaceGeometry::Plane(plane) = model.geometry().surface(top_surface).unwrap()
            else {
                unreachable!("a box's top is a plane")
            };
            plane.point_at(p.x, p.y, T).unwrap()
        };
        if lift(middle + across).y < 0.0 {
            across = -across;
        }
        let line =
            ogeom_geom::Line2d::segment(from + across * shift, to + across * shift, T).unwrap();
        let domain = (0.0, from.distance(to));
        let id = model.geometry_mut().add_pcurve(PlanarCurve::Line(line));
        let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) else {
            unreachable!("an edge")
        };
        for repr in &mut data.representations {
            if let EdgeRepr::PCurve {
                curve,
                surface,
                range,
                ..
            } = repr
                && *surface == top_surface
            {
                *curve = id;
                *range = domain;
            }
        }
        data.tolerance = data.tolerance.widen_to(loose);
        solid
    }

    /// Each face reads an edge as far as it stands off that face: the top
    /// face draws its edge as its own pcurve on its plane, the front face
    /// reads the same edge exactly. An edge stated loose that stands off
    /// no face of it is read to its tolerance on both.
    #[test]
    fn a_face_reads_an_edge_as_far_as_it_stands_off_that_face() {
        let mut model = Model::new();
        let solid = box_with_an_edge_drawn_off_its_top(&mut model, 0.004, 0.005);
        let gathered = gather(&model, &solid, false, None, T).unwrap();
        let mut read = 0;
        for face in &gathered.faces {
            for e in &face.edges {
                if e.tolerance < 1e-3 {
                    continue;
                }
                read += 1;
                let normal = face.surface.normal_at(0.0, 0.0, T).unwrap().vector();
                let (drawn, range, off) = e.drawn();
                let middle = drawn.point_at(f64::midpoint(range.0, range.1), T).unwrap();
                if normal.z.abs() > 0.5 {
                    assert!((e.standoff - 0.004).abs() < 1e-6, "top: {}", e.standoff);
                    assert!(e.image.is_some() && off == 0.0);
                    assert!((middle.y - 0.004).abs() < 1e-6, "top drawn at {middle:?}");
                } else {
                    assert!(e.standoff < 1e-9, "front: {}", e.standoff);
                    assert!(e.image.is_none());
                    assert!(middle.y.abs() < 1e-9, "front drawn at {middle:?}");
                }
            }
        }
        assert_eq!(read, 2);

        let mut model = Model::new();
        let solid = box_with_an_edge_drawn_off_its_top(&mut model, 1e-4, 0.005);
        let gathered = gather(&model, &solid, false, None, T).unwrap();
        for e in gathered.faces.iter().flat_map(|f| f.edges.iter()) {
            if e.tolerance >= 1e-3 {
                assert_eq!(e.standoff, e.tolerance);
                assert!(e.image.is_none());
            }
        }
    }

    /// The crossings either face of a loose edge finds, a lean apart along
    /// it, are one cluster; tight, they stay two. A run of such crossings,
    /// each within its spread of the last, makes no cluster longer than
    /// two spreads.
    #[test]
    fn a_loose_edge_s_crossings_join_across_their_spread_and_no_further() {
        let curve = Curve::Line(
            ogeom_geom::LineCurve::segment(Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), T)
                .unwrap(),
        );
        let crange = ogeom_geom::Curve3d::domain(&curve);
        let pave = |t: f64, spread: f64| Pave {
            t,
            honesty: 0.0,
            across: None,
            spread,
        };
        let two = [pave(1.0, 0.03), pave(1.05, 0.03)];
        let clusters = cluster_paves(&curve, crange, 0.0096, &two, None, &[], T).unwrap();
        assert_eq!(clusters.len(), 1, "{clusters:?}");
        assert!((clusters[0].span - 0.05).abs() < 1e-9, "{clusters:?}");
        let tight = [pave(1.0, 0.0), pave(1.05, 0.0)];
        assert_eq!(
            cluster_paves(&curve, crange, 0.0096, &tight, None, &[], T)
                .unwrap()
                .len(),
            2
        );
        let run: Vec<Pave> = (0..8)
            .map(|k| pave(0.5 + 0.05 * f64::from(k), 0.03))
            .collect();
        let clusters = cluster_paves(&curve, crange, 0.0096, &run, None, &[], T).unwrap();
        assert!(clusters.len() > 1, "{clusters:?}");
        for c in &clusters {
            assert!(c.span <= 0.06 + 1e-9, "{clusters:?}");
        }
    }

    /// A line touching a cylinder of radius 2 parts from it as `s^2 / 4`
    /// along it. A crossing on the line found at the touch, by a section a
    /// micron off the cylinder, may stand `sqrt(2 * 1e-6 * 2)` from the
    /// true one either way; found `d` = 1e-3 along the line, the valley
    /// reaches past the touch to `d + sqrt(d^2 + 4e-6)`.
    /// A line crossing the cylinder at an angle, and an exact section,
    /// spread nothing.
    #[test]
    fn a_crossing_where_the_edge_touches_the_other_surface_spreads_across_its_valley() {
        let cylinder: SurfaceGeometry = ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(Frame::WORLD, 2.0, T).unwrap(),
            (-5.0, 5.0),
        )
        .unwrap()
        .into();
        let line = |x: f64| {
            Curve::Line(
                ogeom_geom::LineCurve::segment(
                    Point::new(x, -3.0, 0.0),
                    Point::new(x, 3.0, 0.0),
                    T,
                )
                .unwrap(),
            )
        };
        let touching = line(2.0);
        let spread_at = |curve: &Curve, y: f64, honesty: f64| {
            let t = y + 3.0;
            let at = curve.point_at(t, T).unwrap();
            touch_spread((curve, t), &cylinder, at, honesty, T)
        };
        let at_touch = spread_at(&touching, 0.0, 1e-6);
        assert!((at_touch - 2e-3).abs() < 1e-9, "{at_touch}");
        let d = 1e-3;
        let beside = spread_at(&touching, d, 1e-6);
        let want = d + d.mul_add(d, 4e-6).sqrt();
        assert!((beside - want).abs() < 1e-6, "{beside} against {want}");
        assert_eq!(spread_at(&touching, 0.0, 0.0), 0.0);
        assert_eq!(spread_at(&line(1.0), 3.0_f64.sqrt(), 1e-6), 0.0);
    }

    /// A crossing of a loose edge spreads along it as its tolerance over
    /// the sine of its lean out of the crossing face, held to ten
    /// tolerances; a tight edge's does not; and where the face's normal
    /// cannot be read, at a cone's apex, the spread is none.
    #[test]
    fn a_crossing_spreads_along_a_loose_edge_by_its_lean() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let mut gathered = gather(&model, &solid, false, None, T).unwrap();
        let e = gathered
            .faces
            .iter_mut()
            .flat_map(|f| f.edges.iter_mut())
            .find(|e| matches!(&*e.curve, Curve::Line(_)))
            .unwrap();
        let t = f64::midpoint(e.crange.0, e.crange.1);
        let at = e.curve.point_at(t, T).unwrap();
        let along = e.curve.d1_at(t, T).unwrap();
        let along = along / along.magnitude();
        let square = if along.z.abs() < 0.9 {
            along.cross(ogeom_math::Vector::Z)
        } else {
            along.cross(ogeom_math::Vector::X)
        };
        let square = square / square.magnitude();
        let plane_leaning = |lean: f64| -> SurfaceGeometry {
            let normal = square * lean.cos() + along * lean.sin();
            let frame = Frame::new(
                at,
                Direction::new(normal, T).unwrap(),
                Direction::new(normal.cross(square).cross(normal), T)
                    .or_else(|_| Direction::new(along, T))
                    .unwrap(),
                T,
            )
            .unwrap();
            ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)).into()
        };
        assert_eq!(lean_spread(e, t, &plane_leaning(0.3), at, T), 0.0);
        e.tolerance = 0.01;
        for lean in [0.3_f64, 1.2] {
            let spread = lean_spread(e, t, &plane_leaning(lean), at, T);
            assert!(
                (spread - 0.01 / lean.sin()).abs() < 1e-9,
                "{lean}: {spread}"
            );
        }
        let spread = lean_spread(e, t, &plane_leaning(0.01), at, T);
        assert!((spread - 0.1).abs() < 1e-12, "{spread}");
        let apex = ogeom_geom::ConeSurface::new(
            ogeom_math::Cone::new(
                Frame::new(
                    at,
                    Direction::new(square, T).unwrap(),
                    Direction::new(along, T).unwrap(),
                    T,
                )
                .unwrap(),
                0.0,
                0.5,
                T,
            )
            .unwrap(),
            (0.0, 5.0),
        )
        .unwrap();
        let spread = lean_spread(e, t, &SurfaceGeometry::Cone(apex), at, T);
        assert!(spread.is_finite() && spread >= 0.0, "{spread}");
    }

    /// A point on a face's surface inside its trim lies within it; one
    /// beyond its boundary lies within it only to the margin asked.
    #[test]
    fn a_point_beside_a_face_lies_within_it_only_to_the_margin() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let gathered = gather(&model, &solid, false, None, T).unwrap();
        let top = gathered
            .faces
            .iter()
            .find(|f| f.bound.low().unwrap().z > 9.9)
            .unwrap();
        assert!(within_face(top, Point::new(10.0, 10.0, 10.0), 1e-7, T).unwrap());
        let beside = Point::new(10.0, -0.05, 10.0);
        assert!(!within_face(top, beside, 1e-4, T).unwrap());
        assert!(within_face(top, beside, 0.06, T).unwrap());
    }

    /// A circle's count from its sag in closed form stands its chords within
    /// the strand sag, short of the cap both counts share, and never above
    /// the count measuring the chords reaches.
    #[test]
    fn a_circle_is_sampled_to_the_strand_sag_in_closed_form() {
        for radius in [0.01, 1.0, 40.0] {
            let circle = ogeom_math::Circle2::new(ogeom_math::Frame2::WORLD, radius, T).unwrap();
            let pcurve = PlanarCurve::Circle(ogeom_geom::Circle2d::new(circle));
            for (lo, hi) in [(0.0, 2.0 * PI), (1.0, 2.5)] {
                let count = circle_sag_count(&pcurve, lo, hi, 8).unwrap();
                let step = (hi - lo) / f64::from(u32::try_from(count).unwrap());
                assert!(
                    count == 512
                        || radius * (1.0 - (step / 2.0).cos()) <= STRAND_SAG * (1.0 + 1e-9)
                );
                let (measured, _) = sag_count(&|t| pcurve.point_at(t, T), lo, hi, 8).unwrap();
                assert!(count <= measured, "{count} against {measured}");
            }
        }
    }

    /// A narrow wedge whose last edge stops short of the column it meets,
    /// once by more than the weld and once by less than the snap: closed,
    /// a point level with the opening and outside the wedge reads outside,
    /// and one inside still reads inside.
    #[test]
    fn an_outline_closes_across_the_openings_its_ends_leave() {
        let top = PI / 2.0;
        for short in [4e-4, 2.6e-8] {
            let mut lines = vec![
                vec![Point2::new(0.0, 0.0), Point2::new(0.06, 0.0)],
                vec![Point2::new(0.06, 0.0), Point2::new(0.06, top)],
                vec![Point2::new(0.06, top - short), Point2::new(0.0, 0.0)],
            ];
            bridge_outline_gaps(&mut lines);
            let lines: Vec<&[Point2]> = lines.iter().map(Vec::as_slice).collect();
            // Level with the opening, a turn round from the wedge, so its
            // ray runs across the wedge.
            let level = Point2::new(-5.298, top - short / 2.0);
            assert!(!arrange::inside_many(&lines, level), "{short}");
            assert!(
                arrange::inside_many(&lines, Point2::new(0.05, 0.5)),
                "{short}"
            );
        }
    }

    fn frame_at(origin: Point) -> Frame {
        Frame::new(origin, Direction::Z, Direction::X, T).unwrap()
    }

    fn boxes(model: &mut Model) -> (Shape, Shape) {
        // Overlapping in the corner cube [1,2]^3.
        let a = make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            model,
            frame_at(Point::new(1.0, 1.0, 1.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        (a.shape, b.shape)
    }

    fn volume(model: &Model, shape: &Shape) -> f64 {
        let fine = Deflection {
            chord: 1e-3,
            ..Deflection::default()
        };
        volume_properties(model, shape, fine, T).unwrap().mass
    }

    fn assert_valid(model: &Model, shape: &Shape) {
        let diagnosis = check(model, shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "the result fails validity: {:?}",
            diagnosis.problems
        );
    }

    #[test]
    fn fuse_of_overlapping_boxes_has_the_inclusion_exclusion_volume() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let fused = fuse(&mut model, &a, &b, T).unwrap();
        assert_valid(&model, &fused.shape);
        assert!((volume(&model, &fused.shape) - 15.0).abs() < 1e-9);
        assert_eq!(
            fused.history.modified(&a),
            std::slice::from_ref(&fused.shape)
        );
    }

    #[test]
    fn common_of_overlapping_boxes_is_the_overlap_cube() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = common(&mut model, &a, &b, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn cut_removes_the_overlap_from_the_first_argument() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = cut(&mut model, &a, &b, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 7.0).abs() < 1e-9);
    }

    #[test]
    fn section_of_overlapping_boxes_is_the_six_segment_seam() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = section(&mut model, &a, &b, T).unwrap();
        let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
        assert_eq!(edges.len(), 6);
        // One closed loop: six edges end to end on six shared vertices.
        let vertices =
            ogeom_topo::explore_unique(&model, &result.shape, ShapeType::Vertex).unwrap();
        assert_eq!(vertices.len(), 6);
    }

    #[test]
    fn fuse_of_disjoint_boxes_is_a_compound_of_both() {
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(5.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let fused = fuse(&mut model, &a.shape, &b.shape, T).unwrap();
        assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Compound);
        assert!((volume(&model, &fused.shape) - 16.0).abs() < 1e-9);
    }

    #[test]
    fn cutting_a_through_post_leaves_a_slab_with_a_hole() {
        let mut model = Model::new();
        let slab = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let post = make_box(
            &mut model,
            frame_at(Point::new(1.5, 1.5, -1.0)),
            (1.0, 1.0, 3.0),
            T,
        )
        .unwrap();
        let result = cut(&mut model, &slab.shape, &post.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 15.0).abs() < 1e-9);
    }

    #[test]
    fn common_with_a_contained_box_is_that_box() {
        let mut model = Model::new();
        let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
        let inner = make_box(
            &mut model,
            frame_at(Point::new(2.0, 2.0, 2.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let result = common(&mut model, &outer.shape, &inner.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 8.0).abs() < 1e-9);
    }

    #[test]
    fn cutting_everything_away_leaves_an_empty_compound() {
        let mut model = Model::new();
        let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
        let inner = make_box(
            &mut model,
            frame_at(Point::new(2.0, 2.0, 2.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let result = cut(&mut model, &inner.shape, &outer.shape, T).unwrap();
        assert_eq!(model.kind_of(&result.shape).unwrap(), ShapeType::Compound);
        assert!(
            explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn drilling_a_box_leaves_a_cylindrical_hole() {
        // Box minus a through-post cylinder: the box's top and bottom faces
        // come back with *circular* holes, the hole's wall is the cylinder's
        // own surface with its material side flipped, and the cylinder's
        // seam and both section circles all split and sew for the shell to
        // close.
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let drill = make_cylinder(
            &mut model,
            frame_at(Point::new(2.0, 2.0, -1.0)),
            0.5,
            3.0,
            T,
        )
        .unwrap();
        let result = cut(&mut model, &block.shape, &drill.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        let exact = 16.0 - PI * 0.25;
        let got = volume(&model, &result.shape);
        assert!(
            (got - exact).abs() / exact < 2e-3,
            "volume {got} against {exact}"
        );
    }

    /// A long rod of radius 5 at right angles to a long drill of radius 56,
    /// the axes 46 apart, runs in through the drill's wall and out again:
    /// its side meets the drill in two loops 64 apart. Both cut the rod, and
    /// what the drill holds of it is the integral of the rod's slices below
    /// the drill's wall.
    #[test]
    fn a_rod_through_a_long_drill_s_wall_is_cut_on_both_sides() {
        let (big, small, offset) = (56.0, 5.0, 46.0);
        let mut model = Model::new();
        let rod = make_cylinder(
            &mut model,
            Frame::new(
                Point::new(5.0, -1000.0, offset),
                Direction::Y,
                Direction::X,
                T,
            )
            .unwrap(),
            small,
            2000.0,
            T,
        )
        .unwrap();
        let drill = make_cylinder(
            &mut model,
            Frame::new(Point::new(-900.0, 0.0, 0.0), Direction::X, Direction::Z, T).unwrap(),
            big,
            1800.0,
            T,
        )
        .unwrap();
        // The slice of the rod at `y` lies in the drill below z = h(y): a
        // disc's segment, integrated by Simpson's rule.
        let segment = |y: f64| -> f64 {
            let d = (big * big - y * y).max(0.0).sqrt() - offset;
            if d <= -small {
                0.0
            } else if d >= small {
                PI * small * small
            } else {
                small * small * (-d / small).acos() + d * (small * small - d * d).sqrt()
            }
        };
        let reach = (big * big - (offset - small).powi(2)).sqrt();
        let steps = 20_000_u32;
        let h = 2.0 * reach / f64::from(steps);
        let held = (0..=steps)
            .map(|i| {
                let w = if i == 0 || i == steps {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                w * segment(-reach + h * f64::from(i))
            })
            .sum::<f64>()
            * h
            / 3.0;
        let whole = PI * small * small * 2000.0;
        let kept = cut(&mut model, &rod.shape, &drill.shape, T).unwrap();
        assert_valid(&model, &kept.shape);
        let inside = common(&mut model, &rod.shape, &drill.shape, T).unwrap();
        assert_valid(&model, &inside.shape);
        let got = volume(&model, &inside.shape);
        assert!(
            (got - held).abs() / held < 2e-3,
            "common {got} against {held}"
        );
        let got = volume(&model, &kept.shape);
        let exact = whole - held;
        assert!(
            (got - exact).abs() / exact < 2e-3,
            "cut {got} against {exact}"
        );
    }

    #[test]
    fn common_of_a_box_and_a_cylinder_is_the_post_inside_it() {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let post = make_cylinder(
            &mut model,
            frame_at(Point::new(2.0, 2.0, -1.0)),
            0.5,
            3.0,
            T,
        )
        .unwrap();
        let result = common(&mut model, &block.shape, &post.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        let exact = PI * 0.25;
        let got = volume(&model, &result.shape);
        assert!(
            (got - exact).abs() / exact < 2e-3,
            "volume {got} against {exact}"
        );
    }

    #[test]
    fn fusing_a_post_onto_a_slab_adds_what_stands_proud() {
        let mut model = Model::new();
        let slab = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let post = make_cylinder(
            &mut model,
            frame_at(Point::new(2.0, 2.0, -1.0)),
            0.5,
            3.0,
            T,
        )
        .unwrap();
        let result = fuse(&mut model, &slab.shape, &post.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        let exact = PI.mul_add(0.25 * 3.0, 16.0) - PI * 0.25;
        let got = volume(&model, &result.shape);
        assert!(
            (got - exact).abs() / exact < 2e-3,
            "volume {got} against {exact}"
        );
    }

    #[test]
    fn crossed_cylinders_run_through_the_marched_sections() {
        // No closed form exists for cylinder/cylinder: the sections are
        // marched and fitted, and every stage downstream (paves against
        // seams on both charts, splitting, classification, rebuilding with
        // fitted pcurves, sewing, meshing) has to work within the fit's
        // stated budget. Unequal radii keep the tangential branch points
        // away. The volumes have no easy closed form, so the operations are
        // held to each other: fuse = A + B - common and cut = A - common are
        // identities whatever the shapes.
        let mut model = Model::new();
        let upright = make_cylinder(&mut model, Frame::WORLD, 1.0, 4.0, T).unwrap();
        let across_frame =
            Frame::new(Point::new(-2.0, 0.0, 2.0), Direction::X, Direction::Y, T).unwrap();
        let across = make_cylinder(&mut model, across_frame, 0.6, 4.0, T).unwrap();

        let both = fuse(&mut model, &upright.shape, &across.shape, T).unwrap();
        assert_valid(&model, &both.shape);
        let shared = common(&mut model, &upright.shape, &across.shape, T).unwrap();
        assert_valid(&model, &shared.shape);
        let pierced = cut(&mut model, &upright.shape, &across.shape, T).unwrap();
        assert_valid(&model, &pierced.shape);

        let va = volume(&model, &upright.shape);
        let vb = volume(&model, &across.shape);
        let vf = volume(&model, &both.shape);
        let vc = volume(&model, &shared.shape);
        let vx = volume(&model, &pierced.shape);

        assert!(vc > 0.0 && vc < vb, "the overlap is real and partial: {vc}");
        assert!(
            (vf - (va + vb - vc)).abs() / vf < 2e-3,
            "fuse {vf} against A + B - common {}",
            va + vb - vc
        );
        assert!(
            (vx - (va - vc)).abs() / vx < 2e-3,
            "cut {vx} against A - common {}",
            va - vc
        );
    }

    #[test]
    fn stacked_boxes_fuse_into_one_solid_and_the_contact_vanishes() {
        // Same-domain contact with opposed materials: the shared rectangle is
        // interior to the union and no face of the result may carry it.
        let mut model = Model::new();
        let lower = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 1.0), T).unwrap();
        let upper = make_box(
            &mut model,
            frame_at(Point::new(0.0, 0.0, 1.0)),
            (2.0, 2.0, 1.0),
            T,
        )
        .unwrap();
        let fused = fuse(&mut model, &lower.shape, &upper.shape, T).unwrap();
        assert_valid(&model, &fused.shape);
        assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
        assert!((volume(&model, &fused.shape) - 8.0).abs() < 1e-9);
    }

    #[test]
    fn a_small_box_on_a_big_one_fuses_with_a_partial_contact() {
        // The big top face splits into the contact rectangle (swallowed)
        // and the surround, which stays and must sew to the small box's
        // walls along the contact's edges.
        let mut model = Model::new();
        let big = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let small = make_box(
            &mut model,
            frame_at(Point::new(1.0, 1.0, 1.0)),
            (2.0, 2.0, 1.0),
            T,
        )
        .unwrap();
        let fused = fuse(&mut model, &big.shape, &small.shape, T).unwrap();
        assert_valid(&model, &fused.shape);
        assert!((volume(&model, &fused.shape) - 20.0).abs() < 1e-9);

        // Cutting the same pair removes nothing but the measure-zero contact:
        // the big box survives whole, its top face's contact piece intact.
        let mut model = Model::new();
        let big = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let small = make_box(
            &mut model,
            frame_at(Point::new(1.0, 1.0, 1.0)),
            (2.0, 2.0, 1.0),
            T,
        )
        .unwrap();
        let result = cut(&mut model, &big.shape, &small.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 16.0).abs() < 1e-9);
    }

    #[test]
    fn flush_walls_fuse_cut_and_meet_with_aligned_contact() {
        // Overlapping boxes sharing flush walls: same-domain contact with
        // *aligned* materials. A = [0,2]^3, B = [1,3]x[0,2]x[0,2]: the y and
        // z walls of the overlap are coplanar with aligned outward normals.
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(1.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let fused = fuse(&mut model, &a.shape, &b.shape, T).unwrap();
        assert_valid(&model, &fused.shape);
        assert!((volume(&model, &fused.shape) - 12.0).abs() < 1e-9);

        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(1.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let shared = common(&mut model, &a.shape, &b.shape, T).unwrap();
        assert_valid(&model, &shared.shape);
        assert!((volume(&model, &shared.shape) - 4.0).abs() < 1e-9);

        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(1.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let cut_result = cut(&mut model, &a.shape, &b.shape, T).unwrap();
        assert_valid(&model, &cut_result.shape);
        assert!((volume(&model, &cut_result.shape) - 4.0).abs() < 1e-9);
    }

    #[test]
    fn every_source_face_is_accounted_for_in_the_history() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = cut(&mut model, &a, &b, T).unwrap();
        for face in explore(&model, &a, Filter::OfType(ShapeType::Face)).unwrap() {
            assert!(
                result.history.is_affected(&face),
                "a face of the first argument vanished from the history"
            );
        }
    }

    /// A plate 10 thick, `rows` by `rows` holes of radius 1.5 drilled
    /// through it 10 apart and 10 from its sides.
    fn drilled_plate(model: &mut Model, rows: u32) -> Shape {
        let size = 10.0 * f64::from(rows) + 10.0;
        let block = make_box(model, Frame::WORLD, (size, size, 10.0), T)
            .unwrap()
            .shape;
        let mut pins = Vec::new();
        for i in 0..rows {
            for j in 0..rows {
                let at = Point::new(10.0 + 10.0 * f64::from(i), 10.0 + 10.0 * f64::from(j), -5.0);
                pins.push(
                    make_cylinder(model, frame_at(at), 1.5, 20.0, T)
                        .unwrap()
                        .shape,
                );
            }
        }
        let pins = model.add_compound(&pins).unwrap();
        cut(model, &block, &pins, T).unwrap().shape
    }

    /// The operation with every face of both solids taken into the fuse.
    fn with_every_face(
        model: &mut Model,
        a: &Shape,
        b: &Shape,
        keep_aside: [bool; 2],
        keep: impl Fn(&FacePiece) -> Option<bool>,
    ) -> Built {
        let fused = general_fuse_classified(model, a, b, None, T).unwrap();
        let kept: Vec<(usize, bool)> = fused
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| keep(p).map(|flip| (i, flip)))
            .collect();
        assemble_result(model, &fused, &kept, keep_aside, a, b, T).unwrap()
    }

    fn cut_keeps(p: &FacePiece) -> Option<bool> {
        match (p.from_a, p.state) {
            (true, PieceState::Out | PieceState::OnOpposed) => Some(false),
            (false, PieceState::In) => Some(true),
            _ => None,
        }
    }

    /// Whether `face`'s box stands clear of each of `tool`'s lumps' by a
    /// millimetre.
    fn clear_of(model: &Model, face: &Shape, tool: &Shape) -> bool {
        let face = ogeom_algo::face_bounds(model, face).unwrap().expanded(1.0);
        explore_unique(model, tool, ShapeType::Solid)
            .unwrap()
            .iter()
            .all(|lump| !face.intersects(&ogeom_algo::tight_bounds(model, lump, T).unwrap()))
    }

    /// The result with faces set aside against the one with every face
    /// taken: the same volume and face count, valid, and every face of `a`
    /// clear of `b` passed through as itself.
    fn same_as_with_every_face(model: &Model, a: &Shape, b: &Shape, local: &Built, whole: &Built) {
        assert_valid(model, &local.shape);
        assert_valid(model, &whole.shape);
        let (v_local, v_whole) = (volume(model, &local.shape), volume(model, &whole.shape));
        assert!(
            (v_local - v_whole).abs() <= 1e-9 * v_whole.abs(),
            "{v_local} against {v_whole}"
        );
        let count = |shape: &Shape| explore_unique(model, shape, ShapeType::Face).unwrap().len();
        assert_eq!(count(&local.shape), count(&whole.shape));
        let mut passed = 0;
        for face in explore_unique(model, a, ShapeType::Face).unwrap() {
            if clear_of(model, &face, b) {
                let copy = local.history.copy_of(&face);
                assert!(
                    copy.is_some_and(|c| c.is_same(&face)),
                    "a face clear of the tool was not passed through"
                );
                passed += 1;
            }
        }
        assert!(passed > 0);
    }

    /// A tool of one lump along two edges of a drilled plate, its box
    /// spanning the plate while its faces reach only the edges: the hole
    /// walls miss every face of it and stand outside it, so they are set
    /// aside and passed through as they were, the result the one with every
    /// face taken.
    #[test]
    fn a_one_lump_tool_spanning_a_plate_leaves_the_holes_it_misses_as_they_were() {
        let mut model = Model::new();
        let plate = drilled_plate(&mut model, 3);
        let along_y = make_box(
            &mut model,
            frame_at(Point::new(-1.0, -1.0, 8.0)),
            (5.0, 42.0, 3.0),
            T,
        )
        .unwrap()
        .shape;
        let along_x = make_box(
            &mut model,
            frame_at(Point::new(-1.0, -1.0, 9.0)),
            (42.0, 5.0, 2.0),
            T,
        )
        .unwrap()
        .shape;
        let corner = fuse(&mut model, &along_y, &along_x, T).unwrap().shape;
        assert_eq!(
            explore_unique(&model, &corner, ShapeType::Solid)
                .unwrap()
                .len(),
            1
        );
        let local = cut(&mut model, &plate, &corner, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &corner, [true, false], cut_keeps);
        assert_valid(&model, &local.shape);
        let (v_local, v_whole) = (volume(&model, &local.shape), volume(&model, &whole.shape));
        assert!(
            (v_local - v_whole).abs() <= 1e-9 * v_whole,
            "{v_local} against {v_whole}"
        );
        let walls: Vec<Shape> = explore_unique(&model, &plate, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                let data = model.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    model.geometry().surface(data.surface),
                    Some(SurfaceGeometry::Cylinder(_))
                )
            })
            .collect();
        assert_eq!(walls.len(), 9);
        for wall in &walls {
            let copy = local.history.copy_of(wall);
            assert!(
                copy.is_some_and(|c| c.is_same(wall)),
                "a hole wall clear of the tool was not passed through"
            );
        }
    }

    #[test]
    fn a_hole_drilled_at_a_plates_corner_leaves_the_rest_as_it_was() {
        let mut model = Model::new();
        let plate = drilled_plate(&mut model, 3);
        let corner = make_cylinder(
            &mut model,
            frame_at(Point::new(5.0, 5.0, -5.0)),
            1.5,
            20.0,
            T,
        )
        .unwrap()
        .shape;
        assert!(
            aside::set_aside(&model, &plate, &corner, T)
                .unwrap()
                .is_some()
        );
        let local = cut(&mut model, &plate, &corner, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &corner, [true, false], cut_keeps);
        same_as_with_every_face(&model, &plate, &corner, &local, &whole);
        let holed = 40.0 * 40.0 * 10.0 - 10.0 * PI * 1.5 * 1.5 * 10.0;
        assert!((volume(&model, &local.shape) - holed).abs() < 1e-6 * holed);
    }

    /// Two pins at opposite corners, cut as one compound whose box spans
    /// the plate: the faces clear of both pins are set aside all the same.
    #[test]
    fn holes_at_opposite_corners_leave_the_faces_between_them_as_they_were() {
        let mut model = Model::new();
        let plate = drilled_plate(&mut model, 3);
        let mut pins = Vec::new();
        for at in [5.0, 35.0] {
            pins.push(
                make_cylinder(&mut model, frame_at(Point::new(at, at, -5.0)), 1.5, 20.0, T)
                    .unwrap()
                    .shape,
            );
        }
        let pins = model.add_compound(&pins).unwrap();
        let local = cut(&mut model, &plate, &pins, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &pins, [true, false], cut_keeps);
        same_as_with_every_face(&model, &plate, &pins, &local, &whole);
        let holed = 40.0 * 40.0 * 10.0 - 11.0 * PI * 1.5 * 1.5 * 10.0;
        assert!((volume(&model, &local.shape) - holed).abs() < 1e-6 * holed);
    }

    /// Four bars along a box's top edges, fused into a frame whose top
    /// stays split where the bars met: each edge between two of its pieces
    /// lies on the box's top and is laid there once, and the cut leaves the
    /// box less the frame.
    #[test]
    fn a_frame_fused_from_bars_cuts_the_rim_off_a_box() {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let bars = [
            ((0.0, 0.0), (10.0, 1.0)),
            ((0.0, 9.0), (10.0, 1.0)),
            ((0.0, 0.0), (1.0, 10.0)),
            ((9.0, 0.0), (1.0, 10.0)),
        ];
        let mut frame: Option<Shape> = None;
        for ((x, y), (dx, dy)) in bars {
            let bar = make_box(
                &mut model,
                frame_at(Point::new(x, y, 9.0)),
                (dx, dy, 1.0),
                T,
            )
            .unwrap()
            .shape;
            frame = Some(match frame {
                None => bar,
                Some(held) => fuse(&mut model, &held, &bar, T).unwrap().shape,
            });
        }
        let frame = frame.unwrap();
        let rimmed = cut(&mut model, &block, &frame, T).unwrap();
        assert_valid(&model, &rimmed.shape);
        let want = 1000.0 - (100.0 - 64.0);
        let got = volume(&model, &rimmed.shape);
        assert!((got - want).abs() < 1e-9 * want, "{got} against {want}");
    }

    #[test]
    fn a_notch_through_a_plates_edge_meets_the_faces_left_as_they_were() {
        // A box notching the corner through two sides, the top and the
        // bottom: the outer rings of all four are split, and the edges
        // they share with the faces set aside stay as they were.
        let mut model = Model::new();
        let plate = drilled_plate(&mut model, 3);
        let notch = make_box(
            &mut model,
            frame_at(Point::new(-1.0, -1.0, -1.0)),
            (5.0, 4.0, 12.0),
            T,
        )
        .unwrap()
        .shape;
        let local = cut(&mut model, &plate, &notch, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &notch, [true, false], cut_keeps);
        same_as_with_every_face(&model, &plate, &notch, &local, &whole);
    }

    #[test]
    fn a_boss_fused_and_kept_in_common_beside_faces_set_aside() {
        let mut model = Model::new();
        let plate = drilled_plate(&mut model, 3);
        let boss = make_cylinder(
            &mut model,
            frame_at(Point::new(5.0, 5.0, 5.0)),
            3.0,
            10.0,
            T,
        )
        .unwrap()
        .shape;
        let local = fuse(&mut model, &plate, &boss, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &boss, [true, true], |p| {
            (p.state == PieceState::Out || (p.state == PieceState::OnAligned && !p.covered))
                .then_some(false)
        });
        same_as_with_every_face(&model, &plate, &boss, &local, &whole);

        let local = common(&mut model, &plate, &boss, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &boss, [false, false], |p| {
            (p.state == PieceState::In || (p.state == PieceState::OnAligned && !p.covered))
                .then_some(false)
        });
        assert_valid(&model, &local.shape);
        let (v_local, v_whole) = (volume(&model, &local.shape), volume(&model, &whole.shape));
        assert!(
            (v_local - v_whole).abs() <= 1e-9 * v_whole,
            "{v_local} against {v_whole}"
        );
        for face in explore_unique(&model, &plate, ShapeType::Face).unwrap() {
            if clear_of(&model, &face, &boss) {
                assert!(local.history.is_deleted(&face));
            }
        }
    }

    /// The corner hole into a plate of 582 faces, with every face but the
    /// top, the bottom and the tool's set aside, comes out the solid the fuse
    /// of every face makes, and passes every face it leaves alone through as
    /// itself.
    #[test]
    #[ignore = "heavy"]
    fn a_corner_hole_in_a_large_plate_is_the_solid_every_face_makes() {
        let mut model = Model::new();
        let plate = drilled_plate(&mut model, 24);
        assert_eq!(
            explore_unique(&model, &plate, ShapeType::Face)
                .unwrap()
                .len(),
            582
        );
        let corner = make_cylinder(
            &mut model,
            frame_at(Point::new(5.0, 5.0, -5.0)),
            1.5,
            20.0,
            T,
        )
        .unwrap()
        .shape;
        let local = cut(&mut model, &plate, &corner, T).unwrap();
        let whole = with_every_face(&mut model, &plate, &corner, [true, false], cut_keeps);
        same_as_with_every_face(&model, &plate, &corner, &local, &whole);
        // The top and the bottom take the new hole; every other face is
        // untouched.
        let copied = explore_unique(&model, &plate, ShapeType::Face)
            .unwrap()
            .iter()
            .filter(|face| local.history.copy_of(face).is_some())
            .count();
        assert_eq!(copied, 580);
    }

    /// The face a primitive builder gave `role`.
    fn face_with_role(model: &Model, solid: &Shape, role: ogeom_core::Role) -> Shape {
        explore_unique(model, solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                matches!(
                    model.provenance_of(f),
                    Some(
                        ogeom_core::Provenance::Primitive { role: r, .. }
                            | ogeom_core::Provenance::Derived { role: r, .. }
                    ) if *r == role
                )
            })
            .unwrap()
    }

    /// The faces of `result` that `history` traces from `source`.
    fn traced_into(model: &Model, result: &Shape, history: &History, source: &Shape) -> Vec<Shape> {
        explore_unique(model, result, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|r| history.trace(source).iter().any(|t| t.is_same(r)))
            .collect()
    }

    /// A cylinder's caps lying in a cube's faces: the common keeps the
    /// cube's piece of each, and each cap traces to that face; the fuse
    /// keeps one face per cap plane too, which both trace to. The
    /// cylinder's side, half outside the cube, is not lost in the common.
    #[test]
    fn a_coincident_piece_dropped_for_the_others_traces_to_it() {
        use ogeom_algo::primitive::roles::{FACE_MAX_Z, FACE_MIN_Z};
        let mut model = Model::new();
        let cube = make_box(&mut model, Frame::WORLD, (100.0, 100.0, 100.0), T)
            .unwrap()
            .shape;
        let cylinder = make_cylinder(&mut model, Frame::WORLD, 50.0, 100.0, T)
            .unwrap()
            .shape;
        for keep_in in [true, false] {
            let built = if keep_in {
                common(&mut model, &cube, &cylinder, T).unwrap()
            } else {
                fuse(&mut model, &cube, &cylinder, T).unwrap()
            };
            let expected = if keep_in {
                PI * 50.0 * 50.0 * 100.0 / 4.0
            } else {
                1e6 + PI * 50.0 * 50.0 * 100.0 * 3.0 / 4.0
            };
            assert_valid(&model, &built.shape);
            let v = volume(&model, &built.shape);
            assert!(
                (v - expected).abs() <= 1e-3 * expected,
                "{v} against {expected}"
            );
            for role in [FACE_MIN_Z, FACE_MAX_Z] {
                let (side, cap) = (
                    face_with_role(&model, &cube, role),
                    face_with_role(&model, &cylinder, role),
                );
                assert!(!built.history.is_deleted(&cap));
                let from_cap = traced_into(&model, &built.shape, &built.history, &cap);
                let from_side = traced_into(&model, &built.shape, &built.history, &side);
                if keep_in {
                    // The common keeps one face per cap plane, the one both
                    // sources trace to and nothing else.
                    assert_eq!(from_cap.len(), 1, "{role:?}: the cap");
                    assert_eq!(from_side.len(), 1, "{role:?}: the cube's face");
                }
                let shared = from_cap
                    .iter()
                    .filter(|f| from_side.iter().any(|g| g.is_same(f)))
                    .count();
                assert_eq!(shared, 1, "common {keep_in}");
            }
        }
    }

    /// A slot cut into a step already cut, its floor partly on the step's
    /// floor: the cut keeps the step's floor there and the slot's floor
    /// traces to it as well as to its own piece; the slot's ends, which lie
    /// in the block's sides with the material on the same side, are cut
    /// away and stay deleted.
    #[test]
    fn a_tool_face_on_a_floor_a_cut_left_traces_to_that_floor() {
        use ogeom_algo::primitive::roles::{FACE_MAX_Z, FACE_MIN_X, FACE_MIN_Z};
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (100.0, 100.0, 50.0), T)
            .unwrap()
            .shape;
        let step = make_box(
            &mut model,
            frame_at(Point::new(50.0, -10.0, 30.0)),
            (60.0, 120.0, 20.0),
            T,
        )
        .unwrap()
        .shape;
        let slot = make_box(
            &mut model,
            frame_at(Point::new(-10.0, 40.0, 30.0)),
            (120.0, 20.0, 30.0),
            T,
        )
        .unwrap()
        .shape;
        let (step_floor, slot_floor) = (
            face_with_role(&model, &step, FACE_MIN_Z),
            face_with_role(&model, &slot, FACE_MIN_Z),
        );
        let stepped = cut(&mut model, &block, &step, T).unwrap();
        let slotted = cut(&mut model, &stepped.shape, &slot, T).unwrap();
        assert_valid(&model, &slotted.shape);
        let history = stepped.history.then(&slotted.history);
        let v = volume(&model, &slotted.shape);
        assert!((v - 380_000.0).abs() <= 1e-6 * 380_000.0, "{v}");
        let from_slot = traced_into(&model, &slotted.shape, &history, &slot_floor);
        let from_step = traced_into(&model, &slotted.shape, &history, &step_floor);
        assert_eq!(from_slot.len(), 2);
        let shared = from_slot
            .iter()
            .filter(|f| from_step.iter().any(|g| g.is_same(f)))
            .count();
        assert_eq!(shared, 1);
        // Faces whose material is gone stay gone: the slot's top and its
        // end outside the block.
        for role in [FACE_MAX_Z, FACE_MIN_X] {
            assert!(history.is_deleted(&face_with_role(&model, &slot, role)));
        }
    }

    /// A cut followed by a move: each face of the cut's inputs that is not
    /// deleted traces to faces of the moved result, and to nothing else.
    #[test]
    fn a_cut_then_a_move_traces_its_inputs_to_the_moved_result() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let cut_out = cut(&mut model, &a, &b, T).unwrap();
        let moved = ogeom_algo::transformed(
            &mut model,
            &cut_out.shape,
            ogeom_math::Transform::translation(ogeom_math::Vector::new(100.0, 0.0, 0.0)),
        )
        .unwrap();
        let history = cut_out.history.then(&moved.history);
        let result = explore_unique(&model, &moved.shape, ShapeType::Face).unwrap();
        let mut traced = 0;
        for face in explore_unique(&model, &a, ShapeType::Face)
            .unwrap()
            .into_iter()
            .chain(explore_unique(&model, &b, ShapeType::Face).unwrap())
        {
            if history.is_deleted(&face) {
                continue;
            }
            let images = history.trace(&face);
            assert!(!images.is_empty());
            for image in images {
                assert!(result.iter().any(|r| r.is_same(image)));
            }
            traced += 1;
        }
        assert_eq!(traced, 9);
    }

    /// A cut followed by a shear: each face, edge and vertex of the cut's
    /// inputs that the cut carries into its result traces through the
    /// composed history to as many sub-shapes of its kind in the sheared
    /// result, and to nothing else.
    #[test]
    fn a_cut_then_a_shear_traces_its_inputs_to_the_sheared_result() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let cut_out = cut(&mut model, &a, &b, T).unwrap();
        let shear = ogeom_math::GeneralTransform::new(
            ogeom_math::Matrix3 {
                rows: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
            ogeom_math::Vector::ZERO,
        );
        let sheared =
            ogeom_algo::general_transformed_shape(&mut model, &cut_out.shape, &shear, T).unwrap();
        let history = cut_out.history.then(&sheared.history);
        for (kind, expected) in [
            (ShapeType::Face, 9),
            (ShapeType::Edge, 15),
            (ShapeType::Vertex, 8),
        ] {
            let between = explore_unique(&model, &cut_out.shape, kind).unwrap();
            let result = explore_unique(&model, &sheared.shape, kind).unwrap();
            let mut traced = 0;
            for part in explore_unique(&model, &a, kind)
                .unwrap()
                .into_iter()
                .chain(explore_unique(&model, &b, kind).unwrap())
            {
                let carried = cut_out.history.trace(&part);
                if carried.is_empty()
                    || !carried.iter().all(|c| between.iter().any(|r| r.is_same(c)))
                {
                    continue;
                }
                let images = history.trace(&part);
                assert_eq!(images.len(), carried.len(), "{kind:?}");
                assert!(history.is_affected(&part), "{kind:?}");
                for image in images {
                    assert!(result.iter().any(|r| r.is_same(image)), "{kind:?}");
                }
                traced += 1;
            }
            assert_eq!(traced, expected, "{kind:?}");
        }
    }

    /// Each edge and vertex of `argument` traced through `history`: deleted,
    /// or traced to shapes of its kind that `result` holds. The counts of
    /// those deleted, unchanged (held by the result as they stand, or
    /// copied) and modified, by kind.
    fn sub_shapes_traced(
        model: &Model,
        argument: &Shape,
        result: &Shape,
        history: &History,
    ) -> [(usize, usize, usize); 2] {
        let mut counts = [(0, 0, 0); 2];
        for (slot, kind) in [ShapeType::Edge, ShapeType::Vertex].into_iter().enumerate() {
            let held = explore_unique(model, result, kind).unwrap();
            for part in explore_unique(model, argument, kind).unwrap() {
                if history.is_deleted(&part) {
                    counts[slot].0 += 1;
                    continue;
                }
                let images = history.trace(&part);
                assert!(!images.is_empty(), "{kind:?}");
                for image in images {
                    assert!(
                        held.iter().any(|h| h.is_same(image)),
                        "a {kind:?} traces to a shape the result does not hold"
                    );
                }
                if history.copy_of(&part).is_some() || !history.is_affected(&part) {
                    counts[slot].1 += 1;
                } else {
                    counts[slot].2 += 1;
                }
            }
        }
        counts
    }

    fn vertex_point(model: &Model, vertex: &Shape) -> Point {
        let data = model.node(vertex).unwrap().data().as_vertex().unwrap();
        vertex.transform(model.datums()).unwrap().apply(data.point)
    }

    /// A drill through a block whose top line, its seam, stands a hair
    /// above the block's top face: the face is cut along two lines either
    /// side of the seam, closer than the weld, which sew into one line the
    /// drill's face meets on both sides. That line is the drill's seam, and
    /// the result is valid with the volume the drill exactly tangent to the
    /// top cuts.
    #[test]
    fn a_drill_whose_seam_grazes_a_face_keeps_it_a_seam() {
        let chord = Deflection::with_chord(0.01).unwrap();
        let volume_at = |rise: f64| {
            let mut model = Model::new();
            let block = make_box(
                &mut model,
                frame_at(Point::new(-20.0, 0.0, 0.0)),
                (40.0, 30.0, 10.0),
                T,
            )
            .unwrap()
            .shape;
            let up = ogeom_math::Direction::Y;
            let seam = ogeom_math::Direction::Z;
            let drill = make_cylinder(
                &mut model,
                Frame::new(Point::new(0.0, -5.0, 7.0 + rise), up, seam, T).unwrap(),
                3.0,
                40.0,
                T,
            )
            .unwrap()
            .shape;
            let out = cut(&mut model, &block, &drill, T).unwrap();
            let d = check(&model, &out.shape, T).unwrap();
            assert!(d.is_valid(), "rise {rise}: {d}");
            volume_properties(&model, &out.shape, chord, T)
                .unwrap()
                .mass
        };
        let tangent = volume_at(0.0);
        for rise in [1e-9, 1e-8, 1e-7, 2e-7] {
            let v = volume_at(rise);
            assert!(
                (v - tangent).abs() < 1e-6 * tangent,
                "rise {rise}: {v} against {tangent}"
            );
        }
    }

    /// A drill tangent to a wall along the edge where a concave round
    /// leaves it, the drill flatter than the round: the drill touches the
    /// wall from the empty side and cuts into the material behind the
    /// round, so its wall lines the cut from the tangent line to where it
    /// leaves the solid. The material is a slab with a notch whose inner
    /// corner is rounded, and the volume removed is the drill's disc over
    /// the material's section, in closed form.
    #[test]
    fn a_drill_tangent_to_a_wall_along_a_rounds_edge_cuts_the_sliver_behind_it() {
        let mut model = Model::new();
        let slab = make_box(&mut model, Frame::WORLD, (40.0, 40.0, 10.0), T)
            .unwrap()
            .shape;
        let notch = make_box(
            &mut model,
            frame_at(Point::new(20.0, 20.0, -1.0)),
            (21.0, 21.0, 12.0),
            T,
        )
        .unwrap()
        .shape;
        let notched = cut(&mut model, &slab, &notch, T).unwrap().shape;
        // The round: the corner square less a cylinder of radius 2 about
        // (22, 22), tangent to the walls y = 20 and x = 20 at (22, 20) and
        // (20, 22).
        let square = make_box(
            &mut model,
            frame_at(Point::new(20.0, 20.0, 0.0)),
            (2.0, 2.0, 10.0),
            T,
        )
        .unwrap()
        .shape;
        let bore = make_cylinder(
            &mut model,
            frame_at(Point::new(22.0, 22.0, -1.0)),
            2.0,
            12.0,
            T,
        )
        .unwrap()
        .shape;
        let round = cut(&mut model, &square, &bore, T).unwrap().shape;
        let part = fuse(&mut model, &notched, &round, T).unwrap().shape;
        let chord = Deflection::with_chord(0.01).unwrap();
        let before = volume_properties(&model, &part, chord, T).unwrap().mass;
        for seam in [
            Direction::X,
            Direction::new(ogeom_math::Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
        ] {
            let drill = make_cylinder(
                &mut model,
                Frame::new(Point::new(22.0, 50.0, -1.0), Direction::Z, seam, T).unwrap(),
                30.0,
                12.0,
                T,
            )
            .unwrap()
            .shape;
            let out = cut(&mut model, &part, &drill, T).unwrap();
            check(&model, &out.shape, T).unwrap();
            // The integral of sqrt(r^2 - u^2).
            let seg = |u: f64, r: f64| 0.5 * (u * (r * r - u * u).sqrt() + r * r * (u / r).asin());
            // Over the block x in [0, 20]: between y = 40 and the drill's
            // lower arc, y = 50 - sqrt(900 - (x - 22)^2).
            let block = (seg(-2.0, 30.0) - seg(-22.0, 30.0)) - 10.0 * 20.0;
            // Over the corner square x in [20, 22]: the square less its
            // quarter disc, less what lies below the drill's arc.
            let corner = (4.0 - (seg(0.0, 2.0) - seg(-2.0, 2.0)))
                - (60.0 - (seg(0.0, 30.0) - seg(-2.0, 30.0)));
            let expected = before - 10.0 * (block + corner);
            let after = volume_properties(&model, &out.shape, chord, T)
                .unwrap()
                .mass;
            assert!(
                (after - expected).abs() < 1e-6 * expected,
                "seam {seam:?}: {after} against {expected}"
            );
        }
    }

    /// A box cut by one overlapping its corner: the target's three edges
    /// through the corner are shortened and the corner deleted; the tool's
    /// three edges through its corner inside the target line the cavity,
    /// cut short, and the rest of the tool is deleted.
    #[test]
    fn a_cut_records_every_edge_and_vertex_of_both_boxes() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let out = cut(&mut model, &a, &b, T).unwrap();
        let h = &out.history;
        // Edges (deleted, copied, modified), then vertices.
        assert_eq!(
            sub_shapes_traced(&model, &a, &out.shape, h),
            [(0, 9, 3), (1, 7, 0)]
        );
        assert_eq!(
            sub_shapes_traced(&model, &b, &out.shape, h),
            [(9, 0, 3), (7, 0, 1)]
        );
        for vertex in explore_unique(&model, &a, ShapeType::Vertex).unwrap() {
            let at = vertex_point(&model, &vertex);
            assert_eq!(
                h.is_deleted(&vertex),
                at.distance(Point::new(2.0, 2.0, 2.0)) < 1e-9
            );
        }
        // Each shortened edge is one edge of length one.
        for edge in explore_unique(&model, &a, ShapeType::Edge).unwrap() {
            if h.copy_of(&edge).is_some() || !h.is_affected(&edge) {
                continue;
            }
            let [image] = h.modified(&edge) else {
                panic!("a shortened edge is one edge");
            };
            let ends: Vec<Point> = explore_unique(&model, image, ShapeType::Vertex)
                .unwrap()
                .iter()
                .map(|v| vertex_point(&model, v))
                .collect();
            assert_eq!(ends.len(), 2);
            assert!((ends[0].distance(ends[1]) - 1.0).abs() < 1e-9);
        }
    }

    /// A plate drilled through: its edges and vertices come through as
    /// copies, the drill's side seam is cut down to the plate's thickness,
    /// and its rims and their vertices, outside the plate, are deleted.
    #[test]
    fn a_drill_records_every_edge_and_vertex_of_plate_and_drill() {
        let mut model = Model::new();
        let plate = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let drill = make_cylinder(
            &mut model,
            frame_at(Point::new(10.0, 10.0, -1.0)),
            3.0,
            12.0,
            T,
        )
        .unwrap()
        .shape;
        let out = cut(&mut model, &plate, &drill, T).unwrap();
        let h = &out.history;
        assert_eq!(
            sub_shapes_traced(&model, &plate, &out.shape, h),
            [(0, 12, 0), (0, 8, 0)]
        );
        assert_eq!(
            sub_shapes_traced(&model, &drill, &out.shape, h),
            [(2, 0, 1), (2, 0, 0)]
        );
    }

    /// A cut composed with a move: what the cut records of each edge and
    /// vertex is carried to the moved result, and a deletion stays one.
    #[test]
    fn a_cut_then_a_move_traces_every_edge_and_vertex() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let out = cut(&mut model, &a, &b, T).unwrap();
        let moved = ogeom_algo::transformed(
            &mut model,
            &out.shape,
            ogeom_math::Transform::translation(ogeom_math::Vector::new(5.0, 0.0, 0.0)),
        )
        .unwrap();
        let composed = out.history.then(&moved.history);
        // A move modifies what it moves: the cut's copies are modified
        // in the composed history, and the cut's deletions stay deleted.
        let carried = |counts: [(usize, usize, usize); 2]| counts.map(|(d, c, m)| (d, c + m));
        for argument in [&a, &b] {
            let composed_counts = sub_shapes_traced(&model, argument, &moved.shape, &composed);
            assert!(composed_counts.iter().all(|&(_, copied, _)| copied == 0));
            assert_eq!(
                carried(composed_counts),
                carried(sub_shapes_traced(
                    &model,
                    argument,
                    &out.shape,
                    &out.history
                )),
            );
            for kind in [ShapeType::Edge, ShapeType::Vertex] {
                for part in explore_unique(&model, argument, kind).unwrap() {
                    assert_eq!(composed.trace(&part).len(), out.history.trace(&part).len());
                }
            }
        }
    }
}
