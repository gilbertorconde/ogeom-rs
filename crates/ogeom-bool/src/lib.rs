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
//! come from [`intersect_surfaces`](ogeom_intersect::intersect_surfaces), exact
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

use ogeom_algo::{
    Built, Containment, History, classify_in_solid_exact, is_shell_closed, make_edge_between,
    make_face_on, make_vertex, make_wire, sew, shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve2d as _;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry, Transformable};
use ogeom_math::{Point, Point2};
use ogeom_topo::{
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
    node: ogeom_topo::TShapeId,
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
    bound: ogeom_math::Aabb,
    edges: Vec<BoundaryEdge>,
}

/// An argument solid.
struct GSolid {
    solid: Shape,
    faces: Vec<GFace>,
}

fn gather(model: &Model, solid: &Shape, tol: Tolerances) -> OgeomResult<GSolid> {
    if model.kind_of(solid)? != ShapeType::Solid {
        ogeom_bail!(Construction, "boolean arguments are solids");
    }
    for shell in explore_unique(model, solid, ShapeType::Shell)? {
        if !is_shell_closed(model, &shell)? {
            ogeom_bail!(Construction, "an open shell bounds no volume to operate on");
        }
    }

    let mut faces = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
            ogeom_bail!(
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
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let NodeData::Edge(edge_data) = edge_node.data() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
                // A degenerate edge — a cone's apex — has no extent to split.
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            let world = geometry.transformed(&edge.transform(model.datums())?, tol)?;
            let (pcurve, prange, other_side) =
                match edge_data.pcurve_for(surface_id, edge.location()) {
                    Some(EdgeRepr::PCurve {
                        curve: pc, range, ..
                    }) => {
                        let Some(planar) = model.geometry().pcurve(*pc) else {
                            ogeom_bail!(Dangling, "pcurve is not in this model");
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
                            ogeom_bail!(Dangling, "seam pcurve is not in this model");
                        };
                        (f.clone(), *range, Some((r.clone(), *range)))
                    }
                    _ => ogeom_bail!(
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
            ogeom_bail!(Construction, "a face with no boundary bounds nothing");
        }
        let mut bound = ogeom_math::Aabb::EMPTY;
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
        ogeom_bail!(Construction, "a solid with no faces bounds nothing");
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
) -> OgeomResult<Vec<Point2>> {
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
    /// How far the curve may sit from the true intersection: zero for an
    /// exact section, the trace-plus-fit budget for a marched one. Crossing
    /// filters widen their acceptance by this, or a fitted section would
    /// never register against the edges it genuinely meets.
    tolerance: f64,
}

/// One kept sub-range of one section.
#[derive(Clone)]
struct SectionPiece {
    section: usize,
    range: (f64, f64),
}

/// One boundary edge of one argument's face, lying in a face of the other
/// argument's own surface — the splitting curve same-domain contact
/// contributes, since coincident surfaces have no section curve to offer.
struct ContactRec {
    /// The owner's world curve and the sub-range its edge covers.
    curve: Curve,
    crange: (f64, f64),
    /// The curve spoken in the *target* face's chart.
    pcurve: PlanarCurve,
    /// The owner edge's node, whose paves this record shares.
    node: ogeom_topo::TShapeId,
    /// The target: which argument's face list, and which face.
    target_from_a: bool,
    target_face: usize,
}

/// What a face's arrangement strand stands for.
#[derive(Clone)]
enum Tag {
    /// A sub-range of boundary edge `edge` (index into the face's edges).
    Boundary { edge: usize, range: (f64, f64) },
    /// A sub-range of a global section.
    Section { section: usize, range: (f64, f64) },
    /// A sub-range of a global contact edge — another face's boundary edge
    /// lying in this face's own surface, splitting it.
    Contact { contact: usize, range: (f64, f64) },
}

/// Where a piece stands relative to the other solid.
///
/// `In` and `Out` are the classifier's words. The two `On` states split the
/// case the classifier cannot decide alone: a piece lying on the other
/// boundary bounds material on one side here and one side there, and whether
/// those sides agree — outward normals aligned — or oppose is what every
/// operation's filter turns on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PieceState {
    In,
    Out,
    /// On the other boundary, outward normals aligned: both solids' material
    /// on the same side. One copy of the piece bounds the union and the
    /// intersection alike.
    OnAligned,
    /// On the other boundary, outward normals opposed: material on both
    /// sides. The contact is interior to the union and vanishes from it.
    OnOpposed,
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
) -> OgeomResult<(
    Vec<SectionRec>,
    Vec<SectionPiece>,
    Vec<ContactRec>,
    Vec<Vec<(f64, f64)>>,
    std::collections::HashMap<ogeom_topo::TShapeId, Vec<f64>>,
    Vec<Vec<usize>>,
    Vec<Vec<usize>>,
)> {
    use ogeom_intersect::{
        CurveCurveOptions, IntersectOptions, SurfaceIntersection, intersect_curves,
        intersect_surfaces,
    };
    let mut sections: Vec<SectionRec> = Vec::new();
    let mut contacts: Vec<ContactRec> = Vec::new();
    let mut same_pairs: Vec<(usize, usize)> = Vec::new();
    // Marched sections are fitted; the fit is driven below the confusion
    // tolerance so a fitted curve meets edges, vertices and the mesh welder
    // on the same terms as an exact one. The budget each still carries is
    // recorded per section and widens the crossing filters.
    let options = IntersectOptions {
        tolerance: tol.confusion() * 0.5,
        marching: ogeom_intersect::Marching {
            chord: tol.confusion() * 0.5,
            ..ogeom_intersect::Marching::default()
        },
    };
    for (ia, fa) in ga.faces.iter().enumerate() {
        for (ib, fb) in gb.faces.iter().enumerate() {
            if !fa.bound.intersects(&fb.bound) {
                // The faces cannot meet, whatever their surfaces do.
                continue;
            }
            let met = intersect_surfaces(&fa.surface, &fb.surface, options, tol)?;
            match met {
                SurfaceIntersection::Apart => {}
                SurfaceIntersection::Same => {
                    // Coincident surfaces offer no section curve; what splits
                    // each face is the *other* face's boundary. Exact pcurve
                    // projection carries an edge into the other chart, and
                    // planes always have one; a curved same-domain pair whose
                    // edges do not project in closed form is still refused.
                    same_pairs.push((ia, ib));
                    for (owner_from_a, owner, target_from_a, target, target_face) in
                        [(false, fb, true, fa, ia), (true, fa, false, fb, ib)]
                    {
                        let _ = owner_from_a;
                        for e in &owner.edges {
                            let Some(pcurve) =
                                ogeom_intersect::exact_pcurve_of(&e.curve, &target.surface, tol)
                            else {
                                ogeom_bail!(
                                    NotDone,
                                    "same-domain contact whose edges have no \
                                     closed-form projection into the shared \
                                     surface's chart is refused — see the \
                                     deferred table in docs/SCOPE.md"
                                );
                            };
                            contacts.push(ContactRec {
                                curve: e.curve.clone(),
                                crange: e.crange,
                                pcurve,
                                node: e.node,
                                target_from_a,
                                target_face,
                            });
                        }
                    }
                }
                SurfaceIntersection::Touching(points)
                    if points
                        .iter()
                        .any(|p| fa.bound.contains(*p) && fb.bound.contains(*p)) =>
                {
                    ogeom_bail!(
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
                                tolerance: sc.tolerance,
                            }),
                            _ => {
                                // An exact curve whose projection has no
                                // closed form: march the pair instead, so
                                // curve and pcurves are fitted *together*.
                                for fitted in march_pair(&fa.surface, &fb.surface, &options, tol)? {
                                    sections.push(SectionRec {
                                        closed: fitted.curve.is_closed(tol),
                                        tolerance: options.marching.chord + fitted.fit_error,
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
    let outline = |face: &GFace| -> OgeomResult<Vec<Vec<Point2>>> {
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
    let mut paves: std::collections::HashMap<ogeom_topo::TShapeId, Vec<f64>> =
        std::collections::HashMap::new();
    let mut pieces: Vec<SectionPiece> = Vec::new();
    for (si, section) in sections.iter().enumerate() {
        // A fitted section meets an edge within its own budget, not within
        // rounding.
        let reach = tol.confusion().max(section.tolerance * 2.0);
        let cc = CurveCurveOptions {
            gap: reach.max(CurveCurveOptions::default().gap),
            ..CurveCurveOptions::default()
        };
        let domain = section.curve.domain();
        let mut trim_ts: Vec<f64> = Vec::new();
        let mut edge_hits: Vec<(ogeom_topo::TShapeId, f64, f64)> = Vec::new();
        // Spans of the section running *along* a boundary edge. The split
        // such a span would make already exists as boundary — stacked boxes'
        // perpendicular side planes meet exactly at the boxes' own edges —
        // so the span is excluded from the kept intervals rather than
        // refused or duplicated.
        let mut along: Vec<(f64, f64)> = Vec::new();
        for (face_index, own) in [
            (section.face_a, &ga.faces[section.face_a]),
            (section.face_b, &gb.faces[section.face_b]),
        ] {
            let _ = face_index;
            for e in &own.edges {
                let found = intersect_curves(&section.curve, &e.curve, cc, tol)?;
                for crossing in &found.crossings {
                    if crossing.gap > reach {
                        continue;
                    }
                    trim_ts.push(crossing.on_a);
                    edge_hits.push((e.node, crossing.on_b, crossing.on_a));
                }
                for overlap in &found.overlaps {
                    let (lo, hi) = if overlap.on_a.0 <= overlap.on_a.1 {
                        overlap.on_a
                    } else {
                        (overlap.on_a.1, overlap.on_a.0)
                    };
                    trim_ts.push(lo);
                    trim_ts.push(hi);
                    along.push((lo, hi));
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
            let both = reach.max(tol.confusion().max(other.tolerance * 2.0));
            let cc2 = CurveCurveOptions {
                gap: both.max(CurveCurveOptions::default().gap),
                ..CurveCurveOptions::default()
            };
            let found = intersect_curves(&section.curve, &other.curve, cc2, tol)?;
            for crossing in &found.crossings {
                if crossing.gap <= both {
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

        let inside_both = |t: f64| -> OgeomResult<bool> {
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
            let mid = f64::midpoint(lo, hi);
            let mid_folded = if section.closed {
                fold(mid, domain)
            } else {
                mid
            };
            if along
                .iter()
                .any(|(alo, ahi)| mid_folded >= *alo && mid_folded <= *ahi)
            {
                continue;
            }
            if !inside_both(mid)? {
                continue;
            }
            // A section that runs along a boundary edge of either face
            // splits nothing: the split already exists as boundary. The
            // analytic overlap detection above catches the same-support
            // cases; this catches the rest — a fitted section tracing a
            // boundary curve, a surface meeting another exactly at its own
            // trim — by measurement rather than by recognising supports.
            let hugs_boundary = {
                let mut all_near = true;
                for i in 0..=4 {
                    let t = lo + (hi - lo) * f64::from(i) / 4.0;
                    let tf = if section.closed { fold(t, domain) } else { t };
                    let at = section.curve.point_at(tf, tol)?;
                    let mut near = false;
                    'owners: for own in [&ga.faces[section.face_a], &gb.faces[section.face_b]] {
                        for e in &own.edges {
                            // Wider than the crossing filters on purpose: a
                            // tangentially-traced curve wobbles about the
                            // boundary it hugs by far more than a fit budget,
                            // and a genuine section keeps a distance of
                            // feature scale, not microns.
                            if distance_to_edge_curve(&e.curve, e.crange, at, tol)?
                                <= reach.max(tol.confusion() * 1e3)
                            {
                                near = true;
                                break 'owners;
                            }
                        }
                    }
                    if !near {
                        all_near = false;
                        break;
                    }
                }
                all_near
            };
            if hugs_boundary {
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
    // A contact edge splits where it crosses the target face's boundary —
    // and that crossing is a pave on *both* edges, so the owner's own face
    // splits its boundary consistently and the pieces sew back shared.
    let cc = CurveCurveOptions::default();
    let mut contact_along: Vec<Vec<(f64, f64)>> = vec![Vec::new(); contacts.len()];
    for (ci, contact) in contacts.iter().enumerate() {
        let target = if contact.target_from_a {
            &ga.faces[contact.target_face]
        } else {
            &gb.faces[contact.target_face]
        };
        for e in &target.edges {
            let found = intersect_curves(&contact.curve, &e.curve, cc, tol)?;
            for crossing in &found.crossings {
                if crossing.gap > tol.confusion() {
                    continue;
                }
                if crossing.on_a < contact.crange.0 + tol.parametric()
                    || crossing.on_a > contact.crange.1 - tol.parametric()
                {
                    continue;
                }
                // Touching is not crossing, at curve level as at surface
                // level: where the two curves run nearly parallel — a blend
                // arc meeting the edge it is tangent to — the "crossing" is
                // numerical noise smeared along the contact, and paving it
                // would plant a vertex a hair off both curves.
                if tangential(&contact.curve, crossing.on_a, &e.curve, crossing.on_b, tol)? {
                    continue;
                }
                paves.entry(contact.node).or_default().push(crossing.on_a);
                if crossing.on_b > e.crange.0 + tol.parametric()
                    && crossing.on_b < e.crange.1 - tol.parametric()
                {
                    paves.entry(e.node).or_default().push(crossing.on_b);
                }
            }
            // A span of the contact running along a target boundary edge
            // splits nothing: it is already boundary on both sides —
            // identically stacked boxes are all such spans — and duplicating
            // it as a strand would cancel the boundary it copies.
            for overlap in &found.overlaps {
                let (lo, hi) = if overlap.on_a.0 <= overlap.on_a.1 {
                    overlap.on_a
                } else {
                    (overlap.on_a.1, overlap.on_a.0)
                };
                paves.entry(contact.node).or_default().push(lo);
                paves.entry(contact.node).or_default().push(hi);
                contact_along[ci].push((lo, hi));
                // The *target* edge splits where the shared stretch ends,
                // exactly as the contact does. Without this, the face across
                // the overlap keeps one long boundary edge where its new
                // neighbours carry two short ones, and sew — which matches
                // edges whole — can pair it with neither.
                let (blo, bhi) = if overlap.on_b.0 <= overlap.on_b.1 {
                    overlap.on_b
                } else {
                    (overlap.on_b.1, overlap.on_b.0)
                };
                for t in [blo, bhi] {
                    if t > e.crange.0 + tol.parametric() && t < e.crange.1 - tol.parametric() {
                        paves.entry(e.node).or_default().push(t);
                    }
                }
            }
        }
    }

    let mut same_a: Vec<Vec<usize>> = vec![Vec::new(); ga.faces.len()];
    let mut same_b: Vec<Vec<usize>> = vec![Vec::new(); gb.faces.len()];
    for (ia, ib) in same_pairs {
        same_a[ia].push(ib);
        same_b[ib].push(ia);
    }
    Ok((
        sections,
        pieces,
        contacts,
        contact_along,
        paves,
        same_a,
        same_b,
    ))
}

/// March a pair whose exact section has no closed-form pcurve.
fn march_pair(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: &ogeom_intersect::IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<Vec<ogeom_intersect::IntersectionCurve>> {
    use ogeom_intersect::{approximate_branch, branches};
    let traced = branches(a, b, options.marching, tol)?;
    let mut out = Vec::new();
    for branch in &traced {
        out.push(approximate_branch(a, b, branch, options.tolerance, tol)?);
    }
    if out.is_empty() {
        // No branch found is two different stories. Surfaces that measurably
        // stand apart over their stated extents — a blend's cylinder and a
        // far-off bevel plane whose ellipse of intersection lies beyond
        // both — simply do not interact, and an empty section is the true
        // answer. Only a pair that comes close and still resolves nothing
        // is beyond the intersector.
        if surfaces_stand_apart(a, b, tol) {
            return Ok(out);
        }
        ogeom_bail!(
            NotDone,
            "an exact section has no closed-form pcurve and marching resolved \
             no branch; the configuration is beyond the intersector's current \
             reach"
        );
    }
    Ok(out)
}

/// Whether two surfaces measurably keep their distance over their stated
/// extents.
///
/// A conservative grid measurement: sample the smaller-extent surface and
/// project each sample onto the other; apart means every sample clears a
/// margin scaled to the extents. Surfaces with unbounded or enormous stated
/// domains — an imported plane's billion units — are never called apart this
/// way, because a grid over them samples nothing.
fn surfaces_stand_apart(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: Tolerances) -> bool {
    use ogeom_geom::Surface as _;
    let extent_of = |s: &SurfaceGeometry| -> f64 {
        let ((ua, ub), (va, vb)) = s.domain();
        (ub - ua).abs().max((vb - va).abs())
    };
    let (sample, against) = if extent_of(a) <= extent_of(b) {
        (a, b)
    } else {
        (b, a)
    };
    let span = extent_of(sample);
    if !span.is_finite() || span > 1e4 {
        return false;
    }
    let ((ua, ub), (va, vb)) = sample.domain();
    const GRID: usize = 9;
    let mut clearance = f64::INFINITY;
    for i in 0..=GRID {
        for j in 0..=GRID {
            #[allow(clippy::cast_precision_loss)]
            let u = ua + (ub - ua) * i as f64 / GRID as f64;
            #[allow(clippy::cast_precision_loss)]
            let v = va + (vb - va) * j as f64 / GRID as f64;
            let Ok(p) = sample.point_at(u, v, tol) else {
                return false;
            };
            let Ok(projection) = ogeom_algo::project_on_surface(against, p, 16, tol) else {
                return false;
            };
            clearance = clearance.min(projection.distance);
            if clearance <= tol.confusion() * 1e3 {
                return false;
            }
        }
    }
    clearance > tol.confusion() * 1e3
}

// --- the general fuse --------------------------------------------------------

/// One piece of one argument face, classified against the other argument.
struct FacePiece {
    /// Which argument's face list, and which face.
    from_a: bool,
    face: usize,
    rings: Vec<Vec<Traversal<Tag>>>,
    state: PieceState,
}

struct GeneralFused {
    a: GSolid,
    b: GSolid,
    sections: Vec<SectionRec>,
    contacts: Vec<ContactRec>,
    pieces: Vec<FacePiece>,
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

/// The chart point of a world point lying on a planar face, if it lands
/// inside the face's trim.
fn chart_point_of(face: &GFace, p: Point, tol: Tolerances) -> Option<Point2> {
    // Closed-form inversion for the analytic surfaces: the same-domain
    // resolution asks "where does this probe sit in the partner's chart", and
    // the partner may be any surface a face melts along — a plane against a
    // plane, but equally a wall band against the cylinder it copies. The
    // reach check keeps the answer honest: a point off the surface has no
    // chart position, whatever the inversion returns.
    use ogeom_math::elementary;
    let reach = tol.confusion() * 10.0;
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
        _ => return None,
    };
    let at = fold_point_into_chart(raw, &face.surface);
    let mut lines: Vec<Vec<Point2>> = Vec::new();
    for e in &face.edges {
        lines.push(pcurve_polyline(&e.pcurve, e.prange, e.crange, e.crange, tol).ok()?);
        // A seam bounds the chart twice — once per column — and a trim test
        // that sees only one side reads half the band as outside.
        if let Some((other, orange)) = &e.other_side {
            lines.push(pcurve_polyline(other, *orange, e.crange, e.crange, tol).ok()?);
        }
    }
    let borrowed: Vec<&[Point2]> = lines.iter().map(Vec::as_slice).collect();
    if !inside_many(&borrowed, at) {
        return None;
    }
    Some(at)
}

/// Whether two curves meet tangentially at a crossing: the parallel-noise
/// gate for pave placement. One degree is far below any deliberate crossing
/// and far above the smear a tangency leaves in the general intersector.
fn tangential(a: &Curve, ta: f64, b: &Curve, tb: f64, tol: Tolerances) -> OgeomResult<bool> {
    let da = a.d1_at(ta, tol)?;
    let db = b.d1_at(tb, tol)?;
    let (ma, mb) = (da.magnitude(), db.magnitude());
    if ma <= tol.confusion() || mb <= tol.confusion() {
        return Ok(true);
    }
    Ok(da.cross(db).magnitude() / (ma * mb) < 2e-2)
}

/// The distance from a point to a bounded edge curve, through a sampling
/// fine enough for the along-boundary question it answers.
fn distance_to_edge_curve(
    curve: &Curve,
    crange: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let mut best = f64::INFINITY;
    let mut previous: Option<Point> = None;
    for i in 0..=48 {
        let t = crange.0 + (crange.1 - crange.0) * f64::from(i) / 48.0;
        let at = curve.point_at(t, tol)?;
        if let Some(last) = previous {
            let d = at - last;
            let len2 = d.dot(d);
            let s = if len2 > 0.0 {
                ((p - last).dot(d) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            best = best.min(p.distance(last + d * s));
        }
        previous = Some(at);
    }
    Ok(best)
}

fn general_fuse(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<GeneralFused> {
    ogeom_core::progress::stage("boolean: gather");
    let ga = gather(model, a, tol)?;
    let gb = gather(model, b, tol)?;
    ogeom_core::progress::stage("boolean: intersect");
    let (sections, section_pieces, contacts, contact_along, paves, same_a, same_b) =
        fill(&ga, &gb, tol)?;

    ogeom_core::progress::stage("boolean: split");
    let mut pieces: Vec<FacePiece> = Vec::new();
    for (from_a, own, other) in [(true, &ga, &gb.solid), (false, &gb, &ga.solid)] {
        for (fi, face) in own.faces.iter().enumerate() {
            ogeom_core::progress::checkpoint()?;
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
            for (ci, contact) in contacts.iter().enumerate() {
                if contact.target_from_a != from_a || contact.target_face != fi {
                    continue;
                }
                // The owner's paves split its edge; the contact strands split
                // at the same parameters, so the sub-edges rebuilt from both
                // sides are the same edges and sew shared.
                let mut stops = vec![contact.crange.0];
                if let Some(ts) = paves.get(&contact.node) {
                    let mut ts: Vec<f64> = ts
                        .iter()
                        .copied()
                        .filter(|t| {
                            *t > contact.crange.0 + tol.parametric()
                                && *t < contact.crange.1 - tol.parametric()
                        })
                        .collect();
                    ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
                    ts.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
                    stops.extend(ts);
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
                        continue;
                    }
                    // Keep only what lies inside this face's trim; the rest
                    // of the owner's boundary splits nothing here.
                    let count = 16;
                    let mut line = Vec::with_capacity(count + 1);
                    for i in 0..=count {
                        #[allow(clippy::cast_precision_loss)]
                        let t = sub.0 + (sub.1 - sub.0) * i as f64 / count as f64;
                        line.push(contact.pcurve.point_at(t, tol)?);
                    }
                    unwrap_polyline(&mut line, &face.surface);
                    fold_into_chart(&mut line, &face.surface);
                    let mid = line[line.len() / 2];
                    let boundary_lines: Vec<&[Point2]> = strands
                        .iter()
                        .filter(|st| st.boundary)
                        .map(|st| st.polyline.as_slice())
                        .collect();
                    if !inside_many(&boundary_lines, mid) {
                        continue;
                    }
                    strands.push(Strand {
                        polyline: line,
                        tag: Tag::Contact {
                            contact: ci,
                            range: sub,
                        },
                        boundary: false,
                    });
                }
            }

            let split = arrange_pieces(&strands, PARAM_SNAP)?;
            for piece in split {
                let probe = face
                    .surface
                    .point_at(piece.interior.x, piece.interior.y, tol)?;
                let state = match classify_in_solid_exact(model, other, probe, tol)? {
                    Containment::In => PieceState::In,
                    Containment::Out => PieceState::Out,
                    Containment::On => {
                        // On the other boundary: same-domain contact. The
                        // partner face on the shared surface decides whether
                        // the two materials lie on the same side or oppose.
                        let partners = if from_a { &same_a[fi] } else { &same_b[fi] };
                        let own_normal = outward_normal(face, piece.interior, tol)?;
                        let mut resolved = None;
                        for &pi in partners {
                            let partner = if from_a { &gb.faces[pi] } else { &ga.faces[pi] };
                            let Some(uv) = chart_point_of(partner, probe, tol) else {
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
                            // band's generosity read proximity as
                            // coincidence. Ask again at a width where it
                            // cannot, and only a genuine edge contact
                            // remains refused.
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
                                        state: PieceState::In,
                                    });
                                    continue;
                                }
                                Containment::Out => {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        state: PieceState::Out,
                                    });
                                    continue;
                                }
                                Containment::On => {}
                            }
                            ogeom_bail!(
                                NotDone,
                                "a piece lies on the other solid's boundary \
                                 with no coincident partner face to compare \
                                 sides against — edge or vertex contact is \
                                 refused rather than resolved — see the \
                                 deferred table in docs/SCOPE.md"
                            );
                        };
                        state
                    }
                };
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
        contacts,
        pieces,
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
}

impl Rebuild<'_> {
    fn vertex(&mut self, p: Point, tol: Tolerances) -> Shape {
        let found = self
            .vertices
            .iter()
            .find(|(q, _)| q.distance(p) <= tol.confusion() * 10.0)
            .map(|(q, shape)| (q.distance(p), shape.clone()));
        if let Some((gap, shape)) = found {
            // Two descriptions of one junction may disagree by a general
            // crossing's residual; the vertex's tolerance is where that
            // disagreement is recorded, so the sub-edges built against
            // either description still reach it honestly.
            if gap > tol.confusion()
                && let Some(node) = self.model.node_mut(&shape)
                && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
            {
                data.tolerance = data.tolerance.widen_to(gap + tol.confusion());
            }
            return shape;
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
    ) -> ogeom_topo::SurfaceId {
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
) -> OgeomResult<Shape> {
    let own = if piece.from_a { &fused.a } else { &fused.b };
    let face = &own.faces[piece.face];
    let surface_id = rebuild.surface_id(fused, piece.from_a, piece.face);

    // Sub-edges cached within the piece, so a seam used from both sides is
    // one edge appearing twice.
    let mut cache: Vec<(usize, u8, (f64, f64), Shape)> = Vec::new();
    let mut wires = Vec::new();
    for ring in &piece.rings {
        let mut edges = Vec::with_capacity(ring.len());
        for traversal in ring {
            let (key_edge, key_kind, range) = match &traversal.tag {
                Tag::Boundary { edge, range } => (*edge, 0_u8, *range),
                Tag::Section { section, range } => (*section, 1, *range),
                Tag::Contact { contact, range } => (*contact, 2, *range),
            };
            let near = |a: (f64, f64), b: (f64, f64)| {
                (a.0 - b.0).abs() <= tol.parametric() && (a.1 - b.1).abs() <= tol.parametric()
            };
            let built = if let Some((.., shape)) = cache
                .iter()
                .find(|(k, s, r, _)| *k == key_edge && *s == key_kind && near(*r, range))
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
                cache.push((key_edge, key_kind, range, shape.clone()));
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
                if let Ok(Some((a, b))) = ogeom_algo::edge_vertices(rebuild.model, edge) {
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
                        edge.orientation() == ogeom_topo::Orientation::Reversed,
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
        if face.face.orientation() == ogeom_topo::Orientation::Reversed {
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
    surface_id: ogeom_topo::SurfaceId,
    tag: &Tag,
    tol: Tolerances,
) -> OgeomResult<Shape> {
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
                None => ogeom_algo::attach_pcurve(
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
            let mid = c.pcurve.point_at(f64::midpoint(range.0, range.1), tol)?;
            let folded = fold_point_into_chart(mid, &face.surface);
            let shifted = c
                .pcurve
                .transformed(&ogeom_math::Transform2::translation(folded - mid), tol)?;
            ogeom_algo::attach_pcurve(
                model,
                &built,
                shifted,
                surface_id,
                Location::identity(),
                *range,
            )?;
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
                pcurve.transformed(&ogeom_math::Transform2::translation(folded - mid), tol)?;
            ogeom_algo::attach_pcurve(
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
) -> OgeomResult<Built> {
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
            ogeom_bail!(
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

/// A shape whose placements carry scale, rebuilt with the scale baked
/// into its geometry before the pipeline runs.
///
/// A scale changes a surface's parameterization out from under its
/// pcurves — the old refusal — so the bake is the whole-shape conversion:
/// surfaces restated in world space, edges moved exactly, pcurves
/// re-derived against the new parameterizations. Unscaled shapes pass
/// through untouched.
fn baked_if_scaled(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Shape> {
    let mut scaled = false;
    for face in ogeom_topo::explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let placement = face.transform(model.datums())?;
        if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
            scaled = true;
            break;
        }
    }
    if !scaled {
        return Ok(shape.clone());
    }
    Ok(ogeom_algo::baked_shape(model, shape, tol)?.shape)
}

/// Whether a shape is a half space: one shell of one planar face, open by
/// construction.
fn is_half_space(model: &Model, shape: &Shape) -> OgeomResult<bool> {
    if model.kind_of(shape)? != ShapeType::Solid {
        return Ok(false);
    }
    let shells = ogeom_topo::explore_unique(model, shape, ShapeType::Shell)?;
    let faces = ogeom_topo::explore_unique(model, shape, ShapeType::Face)?;
    if shells.len() != 1 || faces.len() != 1 {
        return Ok(false);
    }
    let Some(data) = model.node(&faces[0]).and_then(|n| n.data().as_face()) else {
        return Ok(false);
    };
    let planar = matches!(
        model.geometry().surface(data.surface),
        Some(SurfaceGeometry::Plane(_))
    );
    Ok(planar && !ogeom_algo::is_shell_closed(model, &shells[0])?)
}

/// A half space resolved into the solid the operation can act on: a box
/// filling the material side of the boundary plane, sized past the other
/// argument's whole reach.
///
/// The box's plane-side face is *coplanar with the boundary itself*, so the
/// cut the caller sees is the exact plane; the box's far faces stand
/// outside everything the other shape reaches and never appear in the
/// result. A shape that is not a half space passes through untouched.
fn resolved_half_space(
    model: &mut Model,
    shape: &Shape,
    other: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if !is_half_space(model, shape)? {
        return Ok(shape.clone());
    }
    let Some(face) = ogeom_topo::explore_unique(model, shape, ShapeType::Face)?
        .first()
        .cloned()
    else {
        ogeom_bail!(Construction, "the half space lost its face between checks");
    };
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
    Ok(ogeom_algo::make_box(model, placed, (reach, reach, reach), tol)?.shape)
}

/// The union of two solids.
///
/// # Errors
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) for configurations the
/// boolean refuses — tangential or same-domain contact, scaled placements;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for arguments
/// that are not closed solids.
pub fn fuse(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if is_half_space(model, a)? || is_half_space(model, b)? {
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
    let fused = general_fuse(model, a, b, tol)?;
    // Outward pieces bound the union; a same-domain pair with aligned
    // material keeps one copy, and one with opposed material is interior to
    // the union and vanishes.
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            p.state == PieceState::Out || (p.state == PieceState::OnAligned && p.from_a)
        })
        .map(|(i, _)| (i, false))
        .collect();
    assemble_result(model, &fused, &kept, a, b, tol)
}

/// The intersection of two solids.
///
/// # Errors
///
/// As [`fuse`].
pub fn common(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    let fused = general_fuse(model, a, b, tol)?;
    // Inward pieces bound the intersection; an aligned same-domain pair
    // bounds it too, once. An opposed pair encloses no volume between them.
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            p.state == PieceState::In || (p.state == PieceState::OnAligned && p.from_a)
        })
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
pub fn cut(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
    let fused = general_fuse(model, a, b, tol)?;
    // The first argument's outward pieces stay; the tool's inward pieces
    // close the cut with their material side flipped. On the shared surface:
    // an opposed pair means the tool's material is entirely on the other
    // side, so the first argument's face survives untouched; an aligned pair
    // means the tool's material backs the same wall, which the cut removes.
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| match (p.from_a, p.state) {
            (true, PieceState::Out) => Some((i, false)),
            (true, PieceState::OnOpposed) => Some((i, false)),
            (false, PieceState::In) => Some((i, true)),
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
pub fn section(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let (a, b) = (
        &resolved_half_space(model, a, b, tol)?,
        &resolved_half_space(model, b, a, tol)?,
    );
    let (a, b) = (
        &baked_if_scaled(model, a, tol)?,
        &baked_if_scaled(model, b, tol)?,
    );
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
    use ogeom_algo::{check, make_box, make_cylinder, volume_properties};
    use ogeom_math::{Direction, Frame};
    use ogeom_mesh::Deflection;

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
    fn crossed_cylinders_run_through_the_marched_sections() {
        // No closed form exists for cylinder/cylinder: the sections are
        // marched and fitted, and every stage downstream — paves against
        // seams on both charts, splitting, classification, rebuilding with
        // fitted pcurves, sewing, meshing — has to work within the fit's
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
        // The big top face splits into the contact rectangle — swallowed —
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
}
