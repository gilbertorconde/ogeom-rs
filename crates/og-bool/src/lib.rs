//! Boolean operations.
//!
//! *Elsewhere:* `BOPDS`, `BOPAlgo`, `BOPTools`, `IntTools` and `BRepAlgoAPI`.
//!
//! The architectural insight worth preserving: **fuse, common, cut, section and
//! split are one algorithm — general fuse — plus a selection predicate.** The
//! pipeline is
//!
//! 1. a data structure of indexed sub-shapes, per-type interference lists
//!    (V/V, V/E, V/F, E/E, E/F, F/F), pave blocks and common blocks;
//! 2. a pave filler running in strictly increasing dimension, ending in face/face
//!    intersection, section-edge construction and pcurve generation;
//! 3. a builder that splits faces in 2D parametric space, unifies same-domain
//!    faces, rebuilds solids from face sets, and repairs tolerances;
//! 4. filters that select from the general-fuse result.
//!
//! # One pipeline, any analytic face
//!
//! The pipeline runs on solids whose faces are any surface §7's intersectors
//! answer for — planes, cylinders, spheres, cones, tori — with the planar
//! case as nothing more than the case where every curve is a line. Sections
//! come from [`intersect_surfaces`](og_intersect::intersect_surfaces), exact
//! with same-parameter pcurves where the projection has a closed form and
//! *marched and fitted to a stated tolerance* where it does not. Paves — the
//! points where sections meet boundary edges and each other — come from the
//! exact curve/curve intersection. Every face is then split in its own
//! parameter space by an arrangement of *strands*: polyline scaffolding, each
//! naming the exact sub-curve it stands for, so the combinatorics run on
//! polylines and the rebuilt result is exact curves, pcurves attached both
//! sides, sewn back into shared topology by the sewing whose flipped-carry
//! bug this crate found and §9 fixed.
//!
//! Pieces are classified against the other solid by the exact ray classifier
//! and the filters select: fuse keeps what is outside, common what is inside,
//! cut flips the tool's contribution. Same-domain contact — a piece lying
//! *on* the other boundary — and tangential touching are refused with an
//! error naming the deferred entry, never silently mishandled.

mod arrange;

use og_algo::{
    Built, Containment, History, classify_in_solid_exact, is_shell_closed, make_edge_between,
    make_face_on, make_vertex, make_wire, sew, shape_bounds,
};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve2d as _;
use og_geom::Curve3d as _;
use og_geom::Surface as _;
use og_geom::{Curve, PlanarCurve, SurfaceGeometry, Transformable};
use og_math::{Point, Point2};
use og_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, explore, explore_unique,
};

use arrange::{Strand, Traversal, assemble as arrange_pieces, inside_many};

/// Parameter-space chord for the polyline scaffolding.
const SCAFFOLD_CHORD: f64 = 1e-3;

/// Parameter-space node snap for the arrangement.
const PARAM_SNAP: f64 = 1e-6;

// --- gathering ---------------------------------------------------------------

/// One boundary strand source: an edge as one face uses it, with the pcurve
/// that side of it.
struct BoundaryEdge {
    /// The edge's node identity, for sharing paves across faces.
    node: og_topo::TShapeId,
    /// The curve in world space.
    curve: Curve,
    /// The portion the edge covers, in the curve's parameter.
    crange: (f64, f64),
    /// The pcurve on this face's surface.
    pcurve: PlanarCurve,
    /// The portion of the pcurve, mapped proportionally from `crange`.
    prange: (f64, f64),
    /// The other side of a seam, where this edge is one.
    other_side: Option<(PlanarCurve, (f64, f64))>,
}

/// One face, gathered and vetted.
struct GFace {
    face: Shape,
    /// The surface in world space, placement applied.
    surface: SurfaceGeometry,
    /// A generous world bound: the boundary edges' sampled extent, expanded
    /// by most of its own diagonal so a face bulging past its boundary — a
    /// dome past its equator — stays covered. Used only to *gate refusals*:
    /// two faces on one surface, or touching tangentially, are only a
    /// conflict where the faces could actually meet.
    bound: og_math::Aabb,
    edges: Vec<BoundaryEdge>,
}

/// An argument solid.
struct GSolid {
    solid: Shape,
    faces: Vec<GFace>,
}

fn gather(model: &Model, solid: &Shape, tol: Tolerances) -> OgResult<GSolid> {
    if model.kind_of(solid)? != ShapeType::Solid {
        og_bail!(Construction, "boolean arguments are solids");
    }
    for shell in explore_unique(model, solid, ShapeType::Shell)? {
        if !is_shell_closed(model, &shell)? {
            og_bail!(Construction, "an open shell bounds no volume to operate on");
        }
    }

    let mut faces = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface) else {
            og_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
            og_bail!(
                NotDone,
                "a scaled placement changes a surface's parameterization out \
                 from under its pcurves; bake the scale before a boolean"
            );
        }
        let surface = stored.transformed(&placement, tol)?;
        let surface_id = data.surface;

        let mut edges = Vec::new();
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            let Some(edge_node) = model.node(&edge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let NodeData::Edge(edge_data) = edge_node.data() else {
                og_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
                // A degenerate edge — a cone's apex — has no extent to split.
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            let world = geometry.transformed(&edge.transform(model.datums())?, tol)?;
            let (pcurve, prange, other_side) =
                match edge_data.pcurve_for(surface_id, edge.location()) {
                    Some(EdgeRepr::PCurve {
                        curve: pc, range, ..
                    }) => {
                        let Some(planar) = model.geometry().pcurve(*pc) else {
                            og_bail!(Dangling, "pcurve is not in this model");
                        };
                        (planar.clone(), *range, None)
                    }
                    Some(EdgeRepr::Seam {
                        forward,
                        reversed,
                        range,
                        ..
                    }) => {
                        let (Some(f), Some(r)) = (
                            model.geometry().pcurve(*forward),
                            model.geometry().pcurve(*reversed),
                        ) else {
                            og_bail!(Dangling, "seam pcurve is not in this model");
                        };
                        (f.clone(), *range, Some((r.clone(), *range)))
                    }
                    _ => og_bail!(
                        Construction,
                        "an edge with no pcurve on its face cannot be split in \
                         that face's parameter space"
                    ),
                };
            edges.push(BoundaryEdge {
                node: edge.node(),
                curve: world,
                crange: *range,
                pcurve,
                prange,
                other_side,
            });
        }
        if edges.is_empty() {
            og_bail!(Construction, "a face with no boundary bounds nothing");
        }
        let mut bound = og_math::Aabb::EMPTY;
        for e in &edges {
            for i in 0..=16 {
                #[allow(clippy::cast_precision_loss)]
                let t = e.crange.0 + (e.crange.1 - e.crange.0) * f64::from(i) / 16.0;
                bound = bound.with_point(e.curve.point_at(t, tol)?);
            }
        }
        // A plane never bulges past its boundary; anything curved may — a
        // dome past its equator — and gets most of its own diagonal as
        // allowance.
        let bulge = match &surface {
            SurfaceGeometry::Plane(_) => 0.0,
            _ => bound.diagonal() * 0.75,
        };
        let bound = bound.expanded(bulge + tol.confusion() * 1e2);
        faces.push(GFace {
            face,
            surface,
            bound,
            edges,
        });
    }
    if faces.is_empty() {
        og_bail!(Construction, "a solid with no faces bounds nothing");
    }
    Ok(GSolid {
        solid: solid.clone(),
        faces,
    })
}

/// A parameter carried proportionally from one range to another.
fn rescale(t: f64, from: (f64, f64), to: (f64, f64)) -> f64 {
    let span = from.1 - from.0;
    if span.abs() <= f64::MIN_POSITIVE {
        return to.0;
    }
    to.0 + (to.1 - to.0) * (t - from.0) / span
}

/// The pcurve's course over a sub-range of the *curve's* parameters.
fn pcurve_polyline(
    pcurve: &PlanarCurve,
    prange: (f64, f64),
    crange: (f64, f64),
    sub: (f64, f64),
    tol: Tolerances,
) -> OgResult<Vec<Point2>> {
    let lo = rescale(sub.0, crange, prange);
    let hi = rescale(sub.1, crange, prange);
    // Enough samples that the first step approximates the tangent and the
    // scanline interior test has a faithful outline. Straight pcurves get
    // two points; everything else a fixed fine sampling.
    let count = match pcurve {
        PlanarCurve::Line(_) => 1,
        _ => {
            let span = (hi - lo).abs().max(SCAFFOLD_CHORD);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let n = (span / SCAFFOLD_CHORD).sqrt().ceil() as usize;
            n.clamp(8, 256)
        }
    };
    let mut out = Vec::with_capacity(count + 1);
    for i in 0..=count {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / count as f64;
        out.push(pcurve.point_at(t, tol)?);
    }
    Ok(out)
}

// --- the filler --------------------------------------------------------------

/// One section curve between a face of `a` and a face of `b`.
struct SectionRec {
    curve: Curve,
    /// The pcurve on the A face's surface, sharing the curve's parameter.
    pc_a: PlanarCurve,
    /// The same on the B face.
    pc_b: PlanarCurve,
    face_a: usize,
    face_b: usize,
    closed: bool,
}

/// One kept sub-range of one section.
#[derive(Clone)]
struct SectionPiece {
    section: usize,
    range: (f64, f64),
}

/// What a face's arrangement strand stands for.
#[derive(Clone)]
enum Tag {
    /// A sub-range of boundary edge `edge` (index into the face's edges).
    Boundary { edge: usize, range: (f64, f64) },
    /// A sub-range of a global section.
    Section { section: usize, range: (f64, f64) },
}

/// Fold a param-space polyline into a periodic surface's chart, by one
/// constant offset per axis.
///
/// A section's pcurve is *unwrapped* across a seam — continuous, and allowed
/// to leave the stated domain, because that is what crossing a seam is. The
/// arrangement lives in one chart, and the filler has already split every
/// section at its seam crossings, so each strand spans at most one chart
/// width and a single period shift per axis brings it home. The shift is
/// chosen by the strand's midpoint, so endpoints sitting exactly on the
/// chart's edge stay on whichever side the strand's body is.
fn fold_into_chart(line: &mut [Point2], surface: &SurfaceGeometry) {
    let ((ua, ub), (va, vb)) = surface.domain();
    if line.is_empty() {
        return;
    }
    let mid = line[line.len() / 2];
    if surface.is_periodic_u() {
        let span = ub - ua;
        if span > 0.0 {
            let shift = (ua + (mid.x - ua).rem_euclid(span)) - mid.x;
            for p in line.iter_mut() {
                p.x += shift;
            }
        }
    }
    if surface.is_periodic_v() {
        let span = vb - va;
        if span > 0.0 {
            let shift = (va + (mid.y - va).rem_euclid(span)) - mid.y;
            for p in line.iter_mut() {
                p.y += shift;
            }
        }
    }
}

/// Remove period tears from a sampled polyline, axis by axis.
fn unwrap_polyline(line: &mut [Point2], surface: &SurfaceGeometry) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let spans = (
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
    for i in 1..line.len() {
        if spans.0 > 0.0 {
            while line[i].x - line[i - 1].x > spans.0 * 0.5 {
                line[i].x -= spans.0;
            }
            while line[i].x - line[i - 1].x < -spans.0 * 0.5 {
                line[i].x += spans.0;
            }
        }
        if spans.1 > 0.0 {
            while line[i].y - line[i - 1].y > spans.1 * 0.5 {
                line[i].y -= spans.1;
            }
            while line[i].y - line[i - 1].y < -spans.1 * 0.5 {
                line[i].y += spans.1;
            }
        }
    }
}

/// Fold a point into a periodic surface's chart.
fn fold_point_into_chart(p: Point2, surface: &SurfaceGeometry) -> Point2 {
    let mut one = [p];
    fold_into_chart(&mut one, surface);
    one[0]
}

/// A sub-range of a closed curve, brought into its domain.
///
/// The filler split every wrap interval at the domain end, so a piece fits
/// within one period; the fold of its start may still land the end a hair
/// past the domain, which clamps.
fn folded_range(range: (f64, f64), domain: (f64, f64), closed: bool) -> (f64, f64) {
    if !closed {
        return range;
    }
    let f0 = fold(range.0, domain);
    let f1 = (f0 + (range.1 - range.0)).min(domain.1);
    (f0, f1)
}

/// Fold a parameter into a closed curve's domain.
fn fold(t: f64, domain: (f64, f64)) -> f64 {
    let span = domain.1 - domain.0;
    if span <= 0.0 {
        return domain.0;
    }
    domain.0 + (t - domain.0).rem_euclid(span)
}

/// The sections between two gathered solids, with the paves they put on
/// boundary edges.
#[allow(clippy::type_complexity)]
fn fill(
    ga: &GSolid,
    gb: &GSolid,
    tol: Tolerances,
) -> OgResult<(
    Vec<SectionRec>,
    Vec<SectionPiece>,
    std::collections::HashMap<og_topo::TShapeId, Vec<f64>>,
)> {
    use og_intersect::{
        CurveCurveOptions, IntersectOptions, SurfaceIntersection, intersect_curves,
        intersect_surfaces,
    };
    let mut sections: Vec<SectionRec> = Vec::new();
    for (ia, fa) in ga.faces.iter().enumerate() {
        for (ib, fb) in gb.faces.iter().enumerate() {
            if !fa.bound.intersects(&fb.bound) {
                // The faces cannot meet, whatever their surfaces do.
                continue;
            }
            let met =
                intersect_surfaces(&fa.surface, &fb.surface, IntersectOptions::default(), tol)?;
            match met {
                SurfaceIntersection::Apart => {}
                SurfaceIntersection::Same => og_bail!(
                    NotDone,
                    "two faces lie on the same surface; same-domain contact is \
                     refused rather than resolved — see the deferred table in \
                     docs/SCOPE.md"
                ),
                SurfaceIntersection::Touching(points)
                    if points
                        .iter()
                        .any(|p| fa.bound.contains(*p) && fb.bound.contains(*p)) =>
                {
                    og_bail!(
                        NotDone,
                        "two surfaces touch tangentially where both faces \
                         reach; tangential contact is refused rather than \
                         resolved — see the deferred table in docs/SCOPE.md"
                    )
                }
                SurfaceIntersection::Touching(_) => {}
                SurfaceIntersection::Along(curves) => {
                    for sc in curves {
                        match (sc.on_a, sc.on_b) {
                            (Some(pa), Some(pb)) => sections.push(SectionRec {
                                curve: sc.curve,
                                pc_a: pa,
                                pc_b: pb,
                                face_a: ia,
                                face_b: ib,
                                closed: sc.closed,
                            }),
                            _ => {
                                // An exact curve whose projection has no
                                // closed form: march the pair instead, so
                                // curve and pcurves are fitted *together*.
                                for fitted in march_pair(&fa.surface, &fb.surface, tol)? {
                                    sections.push(SectionRec {
                                        closed: fitted.curve.is_closed(tol),
                                        curve: fitted.curve.into(),
                                        pc_a: fitted.on_a.into(),
                                        pc_b: fitted.on_b.into(),
                                        face_a: ia,
                                        face_b: ib,
                                    });
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // Boundary polylines per face, for the trim tests.
    let outline = |face: &GFace| -> OgResult<Vec<Vec<Point2>>> {
        let mut lines = Vec::new();
        for e in &face.edges {
            lines.push(pcurve_polyline(
                &e.pcurve, e.prange, e.crange, e.crange, tol,
            )?);
            if let Some((other, orange)) = &e.other_side {
                lines.push(pcurve_polyline(other, *orange, e.crange, e.crange, tol)?);
            }
        }
        Ok(lines)
    };
    let mut outlines_a = Vec::new();
    for f in &ga.faces {
        outlines_a.push(outline(f)?);
    }
    let mut outlines_b = Vec::new();
    for f in &gb.faces {
        outlines_b.push(outline(f)?);
    }

    // Crossings of each section with the boundary edges of both its faces,
    // and with every other section sharing a face.
    let mut paves: std::collections::HashMap<og_topo::TShapeId, Vec<f64>> =
        std::collections::HashMap::new();
    let options = CurveCurveOptions::default();
    let mut pieces: Vec<SectionPiece> = Vec::new();
    for (si, section) in sections.iter().enumerate() {
        let domain = section.curve.domain();
        let mut trim_ts: Vec<f64> = Vec::new();
        let mut edge_hits: Vec<(og_topo::TShapeId, f64, f64)> = Vec::new();
        for (face_index, own) in [
            (section.face_a, &ga.faces[section.face_a]),
            (section.face_b, &gb.faces[section.face_b]),
        ] {
            let _ = face_index;
            for e in &own.edges {
                let found = intersect_curves(&section.curve, &e.curve, options, tol)?;
                for crossing in &found.crossings {
                    if crossing.gap > tol.confusion() {
                        continue;
                    }
                    trim_ts.push(crossing.on_a);
                    edge_hits.push((e.node, crossing.on_b, crossing.on_a));
                }
                if !found.overlaps.is_empty() {
                    og_bail!(
                        NotDone,
                        "a section runs along a boundary edge; same-domain \
                         contact is refused rather than resolved — see the \
                         deferred table in docs/SCOPE.md"
                    );
                }
            }
        }
        let mut cross_ts: Vec<f64> = Vec::new();
        for (sj, other) in sections.iter().enumerate() {
            if sj == si {
                continue;
            }
            if other.face_a != section.face_a && other.face_b != section.face_b {
                continue;
            }
            let found = intersect_curves(&section.curve, &other.curve, options, tol)?;
            for crossing in &found.crossings {
                if crossing.gap <= tol.confusion() {
                    cross_ts.push(crossing.on_a);
                }
            }
        }

        // Candidate intervals between trim crossings, kept where the middle
        // sits inside both faces' trims.
        trim_ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
        trim_ts.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());
        let mut candidates: Vec<(f64, f64)> = Vec::new();
        if section.closed {
            let period = domain.1 - domain.0;
            if trim_ts.is_empty() {
                candidates.push(domain);
            } else {
                for i in 0..trim_ts.len() {
                    let lo = trim_ts[i];
                    let hi = if i + 1 < trim_ts.len() {
                        trim_ts[i + 1]
                    } else {
                        trim_ts[0] + period
                    };
                    candidates.push((lo, hi));
                }
            }
        } else {
            let mut stops = vec![domain.0];
            stops.extend(
                trim_ts
                    .iter()
                    .copied()
                    .filter(|t| *t > domain.0 && *t < domain.1),
            );
            stops.push(domain.1);
            for pair in stops.windows(2) {
                candidates.push((pair[0], pair[1]));
            }
        }

        let inside_both = |t: f64| -> OgResult<bool> {
            let tf = if section.closed { fold(t, domain) } else { t };
            let ua = fold_point_into_chart(
                section.pc_a.point_at(tf, tol)?,
                &ga.faces[section.face_a].surface,
            );
            let ub = fold_point_into_chart(
                section.pc_b.point_at(tf, tol)?,
                &gb.faces[section.face_b].surface,
            );
            let la: Vec<&[Point2]> = outlines_a[section.face_a]
                .iter()
                .map(Vec::as_slice)
                .collect();
            let lb: Vec<&[Point2]> = outlines_b[section.face_b]
                .iter()
                .map(Vec::as_slice)
                .collect();
            Ok(inside_many(&la, ua) && inside_many(&lb, ub))
        };

        for (lo, hi) in candidates {
            if hi - lo <= tol.parametric() {
                continue;
            }
            if !inside_both(f64::midpoint(lo, hi))? {
                continue;
            }
            // Keep the paves that end a kept interval: those are where edges
            // genuinely split.
            for (node, on_edge, on_section) in &edge_hits {
                let s = *on_section;
                let near = |x: f64| {
                    (s - x).abs() <= tol.parametric()
                        || (section.closed
                            && ((s + (domain.1 - domain.0)) - x).abs() <= tol.parametric())
                };
                if near(lo) || near(hi) {
                    paves.entry(*node).or_default().push(*on_edge);
                }
            }
            // Split at section/section crossings inside the kept interval,
            // so every face sees the same subdivision.
            let mut cuts = vec![lo];
            // A wrap interval also splits at the curve's own domain end, so
            // every piece lives within one period and evaluates in-domain
            // after a single fold.
            if section.closed
                && domain.1 > lo + tol.parametric()
                && domain.1 < hi - tol.parametric()
            {
                cuts.push(domain.1);
            }
            for &c in &cross_ts {
                let c2 = if section.closed && c < lo {
                    c + (domain.1 - domain.0)
                } else {
                    c
                };
                if c2 > lo + tol.parametric() && c2 < hi - tol.parametric() {
                    cuts.push(c2);
                }
            }
            cuts.push(hi);
            cuts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            for pair in cuts.windows(2) {
                let (lo2, hi2) = (pair[0], pair[1]);
                let from = section.curve.point_at(fold(lo2, domain), tol)?;
                let to = section.curve.point_at(fold(hi2, domain), tol)?;
                if from.distance(to) <= tol.confusion() {
                    // A full loop: two arcs, so every strand has two
                    // distinct endpoints.
                    let mid = f64::midpoint(lo2, hi2);
                    pieces.push(SectionPiece {
                        section: si,
                        range: (lo2, mid),
                    });
                    pieces.push(SectionPiece {
                        section: si,
                        range: (mid, hi2),
                    });
                } else {
                    pieces.push(SectionPiece {
                        section: si,
                        range: (lo2, hi2),
                    });
                }
            }
        }
    }
    Ok((sections, pieces, paves))
}

/// March a pair whose exact section has no closed-form pcurve.
fn march_pair(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> OgResult<Vec<og_intersect::IntersectionCurve>> {
    use og_intersect::{Marching, approximate_branch, branches};
    let traced = branches(a, b, Marching::default(), tol)?;
    let mut out = Vec::new();
    for branch in &traced {
        out.push(approximate_branch(a, b, branch, tol.confusion(), tol)?);
    }
    if out.is_empty() {
        og_bail!(
            NotDone,
            "an exact section has no closed-form pcurve and marching resolved \
             no branch; the configuration is beyond the intersector's current \
             reach"
        );
    }
    Ok(out)
}

// --- the general fuse --------------------------------------------------------

/// One piece of one argument face, classified against the other argument.
struct FacePiece {
    /// Which argument's face list, and which face.
    from_a: bool,
    face: usize,
    rings: Vec<Vec<Traversal<Tag>>>,
    state: Containment,
}

struct GeneralFused {
    a: GSolid,
    b: GSolid,
    sections: Vec<SectionRec>,
    pieces: Vec<FacePiece>,
}

fn general_fuse(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<GeneralFused> {
    let ga = gather(model, a, tol)?;
    let gb = gather(model, b, tol)?;
    let (sections, section_pieces, paves) = fill(&ga, &gb, tol)?;

    let mut pieces: Vec<FacePiece> = Vec::new();
    for (from_a, own, other) in [(true, &ga, &gb.solid), (false, &gb, &ga.solid)] {
        for (fi, face) in own.faces.iter().enumerate() {
            let mut strands: Vec<Strand<Tag>> = Vec::new();
            for (ei, e) in face.edges.iter().enumerate() {
                let mut stops = vec![e.crange.0];
                if let Some(ts) = paves.get(&e.node) {
                    let mut ts: Vec<f64> = ts
                        .iter()
                        .copied()
                        .filter(|t| {
                            *t > e.crange.0 + tol.parametric() && *t < e.crange.1 - tol.parametric()
                        })
                        .collect();
                    ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
                    ts.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
                    stops.extend(ts);
                }
                stops.push(e.crange.1);
                // A closed boundary edge — a cap's full circle — needs two
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
                        polyline: pcurve_polyline(&e.pcurve, e.prange, e.crange, sub, tol)?,
                        tag: Tag::Boundary {
                            edge: ei,
                            range: sub,
                        },
                        boundary: true,
                    });
                    if let Some((other_pc, orange)) = &e.other_side {
                        strands.push(Strand {
                            polyline: pcurve_polyline(other_pc, *orange, e.crange, sub, tol)?,
                            tag: Tag::Boundary {
                                edge: ei,
                                range: sub,
                            },
                            boundary: true,
                        });
                    }
                }
            }
            for sp in &section_pieces {
                let section = &sections[sp.section];
                let (belongs, pcurve) = if from_a {
                    (section.face_a == fi, &section.pc_a)
                } else {
                    (section.face_b == fi, &section.pc_b)
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
                // The pcurve shares the curve's parameterization; sampling
                // uses folded parameters for periodic curves.
                let count = 32;
                let mut line = Vec::with_capacity(count + 1);
                for i in 0..=count {
                    #[allow(clippy::cast_precision_loss)]
                    let t = sub.0 + (sub.1 - sub.0) * i as f64 / count as f64;
                    let tf = if section.closed { fold(t, domain) } else { t };
                    line.push(pcurve.point_at(tf, tol)?);
                }
                // Folding the parameter can tear the sampled polyline at the
                // period; unwrap it pointwise, then bring the whole strand
                // into the chart with one shift.
                unwrap_polyline(&mut line, &face.surface);
                fold_into_chart(&mut line, &face.surface);
                strands.push(Strand {
                    polyline: line,
                    tag: Tag::Section {
                        section: sp.section,
                        range: sp.range,
                    },
                    boundary: false,
                });
            }

            let split = arrange_pieces(&strands, PARAM_SNAP)?;
            for piece in split {
                let probe = face
                    .surface
                    .point_at(piece.interior.x, piece.interior.y, tol)?;
                let state = classify_in_solid_exact(model, other, probe, tol)?;
                if state == Containment::On {
                    og_bail!(
                        NotDone,
                        "a piece of a face lies on the other solid's boundary; \
                         same-domain contact is refused rather than resolved — \
                         see the deferred table in docs/SCOPE.md"
                    );
                }
                pieces.push(FacePiece {
                    from_a,
                    face: fi,
                    rings: piece.rings,
                    state,
                });
            }
        }
    }
    Ok(GeneralFused {
        a: ga,
        b: gb,
        sections,
        pieces,
    })
}

// --- rebuilding --------------------------------------------------------------

/// Everything the rebuild shares across pieces.
struct Rebuild<'m> {
    model: &'m mut Model,
    /// World surface ids, minted once per source face.
    surfaces_a: Vec<Option<og_topo::SurfaceId>>,
    surfaces_b: Vec<Option<og_topo::SurfaceId>>,
    /// Vertices shared by position: a wire's connectivity is checked by node
    /// identity, so two sub-edges meeting at a point must *name* the same
    /// vertex, not merely coincide there.
    vertices: Vec<(Point, Shape)>,
}

impl Rebuild<'_> {
    fn vertex(&mut self, p: Point, tol: Tolerances) -> Shape {
        if let Some((_, shape)) = self
            .vertices
            .iter()
            .find(|(q, _)| q.distance(p) <= tol.confusion() * 10.0)
        {
            return shape.clone();
        }
        let shape = make_vertex(self.model, p).shape;
        self.vertices.push((p, shape.clone()));
        shape
    }

    fn surface_id(
        &mut self,
        fused: &GeneralFused,
        from_a: bool,
        face: usize,
    ) -> og_topo::SurfaceId {
        let slot = if from_a {
            &mut self.surfaces_a[face]
        } else {
            &mut self.surfaces_b[face]
        };
        if let Some(id) = slot {
            return *id;
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

/// Build one piece as a face, orientation matching its source face's side.
fn build_piece(
    rebuild: &mut Rebuild,
    fused: &GeneralFused,
    piece: &FacePiece,
    tol: Tolerances,
) -> OgResult<Shape> {
    let own = if piece.from_a { &fused.a } else { &fused.b };
    let face = &own.faces[piece.face];
    let surface_id = rebuild.surface_id(fused, piece.from_a, piece.face);

    // Sub-edges cached within the piece, so a seam used from both sides is
    // one edge appearing twice.
    let mut cache: Vec<(usize, bool, (f64, f64), Shape)> = Vec::new();
    let mut wires = Vec::new();
    for ring in &piece.rings {
        let mut edges = Vec::with_capacity(ring.len());
        for traversal in ring {
            let (key_edge, key_section, range) = match &traversal.tag {
                Tag::Boundary { edge, range } => (*edge, false, *range),
                Tag::Section { section, range } => (*section, true, *range),
            };
            let near = |a: (f64, f64), b: (f64, f64)| {
                (a.0 - b.0).abs() <= tol.parametric() && (a.1 - b.1).abs() <= tol.parametric()
            };
            let built = if let Some((.., shape)) = cache
                .iter()
                .find(|(k, s, r, _)| *k == key_edge && *s == key_section && near(*r, range))
            {
                shape.clone()
            } else {
                let shape = build_sub_edge(
                    rebuild,
                    fused,
                    piece.from_a,
                    face,
                    surface_id,
                    &traversal.tag,
                    tol,
                )?;
                cache.push((key_edge, key_section, range, shape.clone()));
                shape
            };
            edges.push(if traversal.reversed {
                built.reversed()
            } else {
                built
            });
        }
        if let Err(e) = make_wire(rebuild.model, &edges, tol) {
            for (i, edge) in edges.iter().enumerate() {
                if let Ok(Some((a, b))) = og_algo::edge_vertices(rebuild.model, edge) {
                    let pa = rebuild
                        .model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|v| v.point));
                    let pb = rebuild
                        .model
                        .node(&b)
                        .and_then(|n| n.data().as_vertex().map(|v| v.point));
                    eprintln!(
                        "DBG edge {i} rev={} {:?} -> {:?}",
                        edge.orientation() == og_topo::Orientation::Reversed,
                        pa,
                        pb
                    );
                }
            }
            return Err(e);
        }
        wires.push(make_wire(rebuild.model, &edges, tol)?.shape);
    }
    let built = make_face_on(rebuild.model, surface_id, &wires, tol)?.shape;
    Ok(
        if face.face.orientation() == og_topo::Orientation::Reversed {
            built.reversed()
        } else {
            built
        },
    )
}

/// Build the exact sub-edge a tag names, pcurves attached.
fn build_sub_edge(
    rebuild: &mut Rebuild,
    fused: &GeneralFused,
    from_a: bool,
    face: &GFace,
    surface_id: og_topo::SurfaceId,
    tag: &Tag,
    tol: Tolerances,
) -> OgResult<Shape> {
    match tag {
        Tag::Boundary { edge, range } => {
            let e = &face.edges[*edge];
            let from = e.curve.point_at(range.0, tol)?;
            let to = e.curve.point_at(range.1, tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, e.curve.clone(), *range, &v0, &v1, tol)?.shape;
            let sub_p = (
                rescale(range.0, e.crange, e.prange),
                rescale(range.1, e.crange, e.prange),
            );
            match &e.other_side {
                None => og_algo::attach_pcurve(
                    model,
                    &built,
                    e.pcurve.clone(),
                    surface_id,
                    Location::identity(),
                    sub_p,
                )?,
                Some((other, orange)) => {
                    // A seam: both sides attach, and an occurrence picks its
                    // side by its orientation.
                    let sub_o = (
                        rescale(range.0, e.crange, *orange),
                        rescale(range.1, e.crange, *orange),
                    );
                    let _ = sub_o;
                    og_algo::attach_seam(
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
            Ok(built)
        }
        Tag::Section { section, range } => {
            let s = &fused.sections[*section];
            let domain = s.curve.domain();
            let (f0, f1) = folded_range(*range, domain, s.closed);
            let from = s.curve.point_at(fold(f0, domain), tol)?;
            let to = s.curve.point_at(fold(f1, domain), tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, s.curve.clone(), (f0, f1), &v0, &v1, tol)?.shape;
            // The section's pcurve is unwrapped across any seam; the face's
            // triangulator lives in one chart, so the attached copy is folded
            // home by the same period shift the arrangement gave this
            // strand's polyline — decided by the sub-range's midpoint, so an
            // endpoint sitting exactly on the chart's edge stays on the side
            // the arc's body is.
            let pcurve = if from_a { &s.pc_a } else { &s.pc_b };
            let mid = pcurve.point_at(f64::midpoint(f0, f1), tol)?;
            let folded = fold_point_into_chart(mid, &face.surface);
            let shifted =
                pcurve.transformed(&og_math::Transform2::translation(folded - mid), tol)?;
            og_algo::attach_pcurve(
                model,
                &built,
                shifted,
                surface_id,
                Location::identity(),
                (f0, f1),
            )?;
            Ok(built)
        }
    }
}

/// Sew kept pieces, demand closure, and nest shells into solids and voids.
fn assemble_result(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgResult<Built> {
    let mut history = History::new();
    let source_face = |piece: &FacePiece| -> Shape {
        if piece.from_a {
            fused.a.faces[piece.face].face.clone()
        } else {
            fused.b.faces[piece.face].face.clone()
        }
    };

    if kept.is_empty() {
        // A legitimate answer: cutting a solid away entirely leaves nothing.
        let empty = model.add_compound(&[])?;
        for piece in &fused.pieces {
            history.delete(&source_face(piece));
        }
        history.modify(a, empty.clone());
        history.modify(b, empty.clone());
        return Ok(Built::new(empty, history));
    }

    let mut rebuild = Rebuild {
        model,
        surfaces_a: vec![None; fused.a.faces.len()],
        surfaces_b: vec![None; fused.b.faces.len()],
        vertices: Vec::new(),
    };
    let mut faces = Vec::new();
    let mut kept_sources: Vec<Shape> = Vec::new();
    for &(index, flip) in kept {
        let piece = &fused.pieces[index];
        let mut built = build_piece(&mut rebuild, fused, piece, tol)?;
        if flip {
            built = built.reversed();
        }
        history.modify(&source_face(piece), built.clone());
        kept_sources.push(source_face(piece));
        faces.push(built);
    }
    for piece in &fused.pieces {
        let source = source_face(piece);
        if !kept_sources.iter().any(|s| s == &source) {
            history.delete(&source);
        }
    }

    let model = rebuild.model;
    let sewn = sew(model, &faces, tol)?;
    for shell in &sewn.shells {
        if !is_shell_closed(model, shell)? {
            og_bail!(
                NotDone,
                "the kept pieces did not close into a shell; the configuration \
                 is beyond what the boolean currently resolves"
            );
        }
    }

    // Nest: a shell whose bound sits inside another's is that solid's void.
    let mut bounds = Vec::new();
    for shell in &sewn.shells {
        bounds.push(shape_bounds(model, shell, tol)?);
    }
    let mut solids = Vec::new();
    for (i, shell) in sewn.shells.iter().enumerate() {
        let contained = bounds
            .iter()
            .enumerate()
            .any(|(j, other)| j != i && other.contains_box(&bounds[i]));
        if contained {
            continue;
        }
        let mut group = vec![shell.clone()];
        for (j, candidate) in sewn.shells.iter().enumerate() {
            if j != i && bounds[i].contains_box(&bounds[j]) {
                group.push(candidate.clone());
            }
        }
        solids.push(model.add_solid(&group)?);
    }
    let result = if solids.len() == 1 {
        solids.remove(0)
    } else {
        model.add_compound(&solids)?
    };
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Built::new(result, history))
}

// --- the operations ----------------------------------------------------------

/// The union of two solids.
///
/// # Errors
///
/// [`OgError::NotDone`](og_core::OgError::NotDone) for configurations the
/// boolean refuses — tangential or same-domain contact, scaled placements;
/// [`OgError::Construction`](og_core::OgError::Construction) for arguments
/// that are not closed solids.
pub fn fuse(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.state == Containment::Out)
        .map(|(i, _)| (i, false))
        .collect();
    assemble_result(model, &fused, &kept, a, b, tol)
}

/// The intersection of two solids.
///
/// # Errors
///
/// As [`fuse`].
pub fn common(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.state == Containment::In)
        .map(|(i, _)| (i, false))
        .collect();
    assemble_result(model, &fused, &kept, a, b, tol)
}

/// The first solid with the second removed.
///
/// The pieces of `b` that close the cut into `a`'s material bound the removed
/// volume from `b`'s side, so they join the result with their material side
/// flipped.
///
/// # Errors
///
/// As [`fuse`].
pub fn cut(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| match (p.from_a, p.state) {
            (true, Containment::Out) => Some((i, false)),
            (false, Containment::In) => Some((i, true)),
            _ => None,
        })
        .collect();
    assemble_result(model, &fused, &kept, a, b, tol)
}

/// The edges where the two solids' boundaries cross.
///
/// # Errors
///
/// As [`fuse`].
pub fn section(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
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
    for (si, range) in wanted {
        let s = &fused.sections[si];
        let domain = s.curve.domain();
        let (f0, f1) = folded_range(range, domain, s.closed);
        let from = s.curve.point_at(fold(f0, domain), tol)?;
        let to = s.curve.point_at(fold(f1, domain), tol)?;
        let v0 = make_vertex(model, from).shape;
        let v1 = make_vertex(model, to).shape;
        edges.push(make_edge_between(model, s.curve.clone(), (f0, f1), &v0, &v1, tol)?.shape);
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
    use og_algo::{check, make_box, make_cylinder, volume_properties};
    use og_math::{Direction, Frame};
    use og_mesh::Deflection;

    const T: Tolerances = Tolerances::millimetres();
    const PI: f64 = core::f64::consts::PI;

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
        // The curved milestone: box minus a through-post cylinder. The box's
        // top and bottom faces come back with *circular* holes, the hole's
        // wall is the cylinder's own surface with its material side flipped,
        // and the cylinder's seam and both section circles all had to split
        // and sew for the shell to close.
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
}
