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
mod defeature;

pub use defeature::remove_faces;

use ogeom_algo::{
    Built, Containment, History, is_shell_closed, make_edge_between, make_face_on, make_vertex,
    make_wire, sew, shape_bounds,
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

use arrange::{
    Strand, Traversal, assemble as arrange_pieces, inside_many, inside_many_slanted, inside_rings,
};

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
    /// The radius within which the edge honestly lies — a fitted rail can
    /// carry a few dozen microns of construction slop, and every filter that
    /// compares this edge's curve against exact geometry widens by it.
    tolerance: f64,
}

/// A pole: an edge that bounds a face in parameter space and collapses to
/// a point in space.
///
/// A sphere's poles and a cone's apex have no curve to split, but they are
/// half the chart's boundary — leave them out and the face's outline does
/// not close, and nothing can be arranged inside it.
struct PoleEdge {
    /// Where the pole sits in space.
    point: Point,
    /// The chart line the pole runs along.
    pcurve: PlanarCurve,
    prange: (f64, f64),
}

/// One face, gathered and vetted.
struct GFace {
    face: Shape,
    /// The surface in world space, placement applied.
    surface: SurfaceGeometry,
    /// A conservative world bound: the boundary edges' sampled extent plus
    /// the surface's own allowance — nothing for a plane, measured sampling
    /// slack for the ruled kinds whose rulings pin them to their boundary's
    /// hull, most of the diagonal for anything that can genuinely bulge — a
    /// dome past its equator. The poles join after. Gates the pair filter
    /// and refusals; `OGEOM_BOOL_AUDIT_BOUNDS` audits its conservatism.
    bound: ogeom_math::Aabb,
    /// The scale the marching chord derives from — deliberately *not* the
    /// filter box's diagonal, so the filter can tighten without silently
    /// tightening the marcher.
    chord_scale: f64,
    edges: Vec<BoundaryEdge>,
    poles: Vec<PoleEdge>,
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
        let mut poles = Vec::new();
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            let Some(edge_node) = model.node(&edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let NodeData::Edge(edge_data) = edge_node.data() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
                // A degenerate edge — a sphere's pole, a cone's apex — has
                // no extent to split, but it *does* bound the chart, and a
                // chart whose top is missing bounds nothing.
                if let Some(EdgeRepr::PCurve {
                    curve: pc, range, ..
                }) = edge_data.pcurve_for(surface_id, edge.location())
                {
                    let Some(planar) = model.geometry().pcurve(*pc) else {
                        ogeom_bail!(Dangling, "pcurve is not in this model");
                    };
                    let at = explore_unique(model, &edge, ShapeType::Vertex)?;
                    let Some(point) = at
                        .first()
                        .and_then(|v| model.node(v))
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                    else {
                        ogeom_bail!(Construction, "a pole with no vertex is nowhere");
                    };
                    let placed = edge.transform(model.datums())?.apply(point);
                    poles.push(PoleEdge {
                        point: placed,
                        pcurve: planar.clone(),
                        prange: *range,
                    });
                }
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
                tolerance: edge_data.tolerance.get(),
            });
        }
        if edges.is_empty() {
            ogeom_bail!(Construction, "a face with no boundary bounds nothing");
        }

        // A seam edge bounds the chart twice *only when this face wraps* —
        // when its other boundary reaches both columns, as a full drum's
        // rims do. A face that merely sits against the meridian — a sphere
        // octant whose boundary happens to be the seam — uses one column,
        // and feeding the far copy into the arrangement leaves a strand
        // nothing connects to. The decision is made here, once, by chart
        // connectivity, and everything downstream — the arrangement, the
        // trim tests, the rebuild — inherits it.
        let seam_indices: Vec<usize> = edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.other_side.is_some())
            .map(|(i, _)| i)
            .collect();
        for i in seam_indices {
            let mut pool: Vec<Point2> = Vec::new();
            for (j, e) in edges.iter().enumerate() {
                if j == i {
                    continue;
                }
                for t in [e.prange.0, e.prange.1] {
                    pool.push(e.pcurve.point_at(t, tol)?);
                }
                if let Some((other, orange)) = &e.other_side {
                    for t in [orange.0, orange.1] {
                        pool.push(other.point_at(t, tol)?);
                    }
                }
            }
            for p in &poles {
                for t in [p.prange.0, p.prange.1] {
                    pool.push(p.pcurve.point_at(t, tol)?);
                }
            }
            let connected = |pc: &PlanarCurve, range: (f64, f64)| -> OgeomResult<bool> {
                for t in [range.0, range.1] {
                    let at = pc.point_at(t, tol)?;
                    if !pool.iter().any(|q| q.distance(at) <= PARAM_SNAP) {
                        return Ok(false);
                    }
                }
                Ok(true)
            };
            let e = &edges[i];
            let primary_connects = connected(&e.pcurve, e.prange)?;
            let (other_pc, orange) = e
                .other_side
                .clone()
                .unwrap_or_else(|| (e.pcurve.clone(), e.prange));
            let other_connects = connected(&other_pc, orange)?;
            match (primary_connects, other_connects) {
                // Both columns meet the rest of the boundary: the face
                // wraps, and both belong.
                (true, true) => {}
                // One column is this face's; the far copy would dangle.
                (true, false) => edges[i].other_side = None,
                (false, true) => {
                    let e = &mut edges[i];
                    e.pcurve = other_pc;
                    e.prange = orange;
                    e.other_side = None;
                }
                // Neither connects — leave both, and the arrangement's own
                // refusal names the failure as it always did.
                (false, false) => {}
            }
        }
        let mut bound = ogeom_math::Aabb::EMPTY;
        // For a ruled surface the box will be trusted to the boundary's own
        // hull, so the boundary's sampling slack must be measured: how far
        // the true edge sags from the 16-chord polyline, read at each
        // chord's midpoint and doubled for the sag's asymmetry. Nothing else
        // uses the measurement, and the extra evaluations are priced on
        // spline edges, so nothing else pays for it: an unruled face keeps
        // the exact box it always had, and the exact admit set with it.
        let ruled = matches!(
            &surface,
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)
        );
        let mut slack = 0.0_f64;
        for e in &edges {
            let mut previous: Option<Point> = None;
            for i in 0..=16 {
                #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
                let t = e.crange.0 + (e.crange.1 - e.crange.0) * f64::from(i) / 16.0;
                let p = e.curve.point_at(t, tol)?;
                if ruled && let Some(q) = previous {
                    let step = (e.crange.1 - e.crange.0) / 32.0;
                    let mid = e.curve.point_at(t - step, tol)?;
                    slack = slack.max(mid.distance(Point::midpoint(q, p)) * 2.0);
                    bound = bound.with_point(mid);
                }
                bound = bound.with_point(p);
                previous = Some(p);
            }
        }
        // A plane never bulges past its boundary. A ruled surface — cylinder,
        // cone — cannot either: every surface point lies on a straight ruling
        // whose ends are on the boundary, so the face sits inside its
        // boundary's hull and only the boundary's own sampling slack is owed.
        // Anything else may genuinely bulge — a dome past its equator — and
        // keeps most of its own diagonal as allowance. The audit behind
        // OGEOM_BOOL_AUDIT_BOUNDS holds every arm to conservatism.
        let bulge = match &surface {
            SurfaceGeometry::Plane(_) => 0.0,
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => slack,
            _ => bound.diagonal() * 0.75,
        };
        // The scale the marching chord is derived from, decoupled from the
        // filter box. It reproduces exactly what the heuristic was tuned
        // against — the diagonal as the blanket three-quarter bulge left it —
        // because the chord is a tolerance, not a bound: tightening the
        // filter must not silently tighten the marcher, which is exactly
        // what happened when the box fed both (issue #12).
        let margin = tol.confusion() * 1e2;
        let chord_scale = match &surface {
            SurfaceGeometry::Plane(_) => bound.expanded(margin).diagonal(),
            _ => bound.expanded(bound.diagonal() * 0.75 + margin).diagonal(),
        };
        // A pole bounds the chart with no edge to sample: a cone drilled to
        // its apex reaches the apex, and a bound that omits it would let the
        // filter drop a pair the apex genuinely meets. It joins *after* the
        // bulge is taken from the boundary's own diagonal — the pole is an
        // exact point, and letting it stretch the diagonal would inflate a
        // dome's allowance by its own height over again.
        let mut bound = bound.expanded(bulge);
        for pole in &poles {
            bound = bound.with_point(pole.point);
        }
        let bound = bound.expanded(tol.confusion() * 1e2);
        faces.push(GFace {
            poles,
            face,
            surface,
            bound,
            chord_scale,
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

/// A parameter brought onto the turn its edge actually covers.
///
/// A periodic curve's own domain and an edge's range over it need not agree
/// on which turn to count from: a sphere's seam runs its half meridian over
/// `[-π/2, π/2]`, while every crossing found on it comes back in `[0, 2π)`.
/// Left alone, a crossing at latitude `-0.96` arrives as `5.32`, sits outside
/// the range, and is discarded — so the seam never splits where a section
/// genuinely meets it, and the arrangement finds the chain hanging.
fn onto_range(t: f64, curve: &Curve, range: (f64, f64), tol: Tolerances) -> f64 {
    if !curve.is_periodic() {
        return t;
    }
    let (a, b) = curve.domain();
    let period = b - a;
    if period <= 0.0 {
        return t;
    }
    for k in [0.0, -1.0, 1.0, -2.0, 2.0] {
        let shifted = period.mul_add(k, t);
        if shifted >= range.0 - tol.parametric() && shifted <= range.1 + tol.parametric() {
            return shifted;
        }
    }
    t
}

/// The part of an overlap that falls within the second curve's own bounded
/// range, stated in the *first* curve's parameter.
///
/// An overlap is between two curves; an edge covers only part of its curve.
/// The correspondence the overlap states is affine, so the edge's range
/// carries across as an interval — and on a periodic curve it carries across
/// up to whole turns, so the shift that meets the overlap is the one meant.
/// `None` where the edge's own stretch and the overlap do not meet.
fn overlap_within(
    overlap: &ogeom_intersect::Overlap,
    range: (f64, f64),
    curve: &Curve,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let ordered = |r: (f64, f64)| if r.0 <= r.1 { r } else { (r.1, r.0) };
    let span_a = overlap.on_a.1 - overlap.on_a.0;
    let span_b = overlap.on_b.1 - overlap.on_b.0;
    if span_b.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let to_a = |t: f64| overlap.on_a.0 + span_a * (t - overlap.on_b.0) / span_b;
    let (wlo, whi) = ordered((to_a(range.0), to_a(range.1)));
    let (lo, hi) = ordered(overlap.on_a);
    let period = if curve.is_periodic() {
        let (a, b) = curve.domain();
        // The turn measured in the *first* curve's parameter, which is where
        // the intersection is being taken.
        (b - a) * (span_a / span_b).abs()
    } else {
        0.0
    };
    let mut best: Option<(f64, f64)> = None;
    for k in [0.0, 1.0, -1.0, 2.0, -2.0] {
        let candidate = (
            lo.max(period.mul_add(k, wlo)),
            hi.min(period.mul_add(k, whi)),
        );
        if candidate.1 - candidate.0 > best.map_or(tol.parametric(), |(x, y)| y - x) {
            best = Some(candidate);
        }
        if period == 0.0 {
            break;
        }
    }
    best
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

/// The parameter intervals of a contact curve over which *both* faces
/// actually reach it.
///
/// Two surfaces touch along the whole of their contact; two faces touch
/// along whatever part of it their trims both hold. That part is found by
/// sampling — sixty-four stations along the curve, each asked of both
/// charts — rather than by intersecting the contact with the boundary
/// edges, because a contact meets those boundaries tangentially too and the
/// crossing finder is the wrong instrument for it. The cost of sampling is
/// the usual one: a stretch shorter than a station can be missed, and an
/// endpoint is placed within a station of the truth.
fn contact_intervals(
    fused: &GeneralFused,
    contact: &TangentRec,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, f64)>> {
    let outline = |face: &GFace| -> OgeomResult<Vec<Vec<Point2>>> {
        let mut lines = Vec::new();
        for e in &face.edges {
            lines.push(pcurve_polyline(
                &e.pcurve, e.prange, e.crange, e.crange, tol,
            )?);
        }
        Ok(lines)
    };
    let rings_a = outline(&fused.a.faces[contact.face_a])?;
    let rings_b = outline(&fused.b.faces[contact.face_b])?;
    let refs_a: Vec<&[Point2]> = rings_a.iter().map(Vec::as_slice).collect();
    let refs_b: Vec<&[Point2]> = rings_b.iter().map(Vec::as_slice).collect();

    const STATIONS: usize = 64;
    let domain = contact.curve.domain();
    let span = domain.1 - domain.0;
    let mut runs: Vec<(f64, f64)> = Vec::new();
    let mut open: Option<f64> = None;
    for k in 0..=STATIONS {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = span.mul_add(k as f64 / STATIONS as f64, domain.0);
        // A periodic chart's parameters run out past its own window; the
        // trim test is only meaningful once they are folded back into it.
        let held = matches!(
            (
                contact.pc_a.point_at(t, tol),
                contact.pc_b.point_at(t, tol)
            ),
            (Ok(pa), Ok(pb))
                if arrange::inside_many(
                    &refs_a,
                    fold_point_into_chart(pa, &fused.a.faces[contact.face_a].surface),
                ) && arrange::inside_many(
                    &refs_b,
                    fold_point_into_chart(pb, &fused.b.faces[contact.face_b].surface),
                )
        );
        match (held, open) {
            (true, None) => open = Some(t),
            (false, Some(from)) => {
                if t - from > tol.parametric() {
                    runs.push((from, t));
                }
                open = None;
            }
            _ => {}
        }
    }
    if let Some(from) = open
        && domain.1 - from > tol.parametric()
    {
        runs.push((from, domain.1));
    }
    Ok(runs)
}

/// One curve along which two faces *touch* without crossing.
///
/// A contact carries no boundary parity — neither face passes through the
/// other — so it takes no part in the arrangement or the classification.
/// It is still a curve that exists on the result, and a section view
/// through a tangency has to show it, so it is carried alongside.
struct TangentRec {
    curve: Curve,
    pc_a: PlanarCurve,
    pc_b: PlanarCurve,
    face_a: usize,
    face_b: usize,
}

/// One kept sub-range of one section.
#[derive(Clone)]
struct SectionPiece {
    section: usize,
    range: (f64, f64),
    /// Whether this sub-range already *is* boundary on the A face, the B
    /// face, or neither. A piece that is boundary on one side is still the
    /// splitting curve on the other, and only the face it duplicates leaves
    /// it out.
    hugs: [bool; 2],
}

/// One boundary edge of one argument's face, lying in a face of the other
/// argument's own surface — the splitting curve same-domain contact
/// contributes, since coincident surfaces have no section curve to offer.
struct ContactRec {
    /// The owner's world curve and the sub-range its edge covers.
    curve: Curve,
    crange: (f64, f64),
    /// The curve spoken in the *target* face's chart, over its own window,
    /// mapped proportionally from `crange` exactly as an edge's own stored
    /// pcurve is.
    pcurve: PlanarCurve,
    prange: (f64, f64),
    /// The owner edge's node, whose paves this record shares.
    node: ogeom_topo::TShapeId,
    /// The owner edge's own tolerance: how far its curve may honestly sit
    /// from the exact geometry it meets, which is how far a crossing filter
    /// must reach to see a fitted rail cross an exact seam.
    tolerance: f64,
    /// The target: which argument's face list, and which face.
    target_from_a: bool,
    target_face: usize,
}

/// Whether two surfaces are the *identical chart* — the same
/// parameterization, frame and all, not merely the same point set.
///
/// The five analytics only: a fitted surface never qualifies, because two
/// independent fits of one geometry agree nowhere in parameter space. This
/// is what lets a stored pcurve stand in for a closed-form projection in the
/// same-domain melt — on the identical chart, the owner's pcurve already is
/// the projection.
fn same_chart(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: Tolerances) -> bool {
    use SurfaceGeometry as S;
    let frames = |fa: ogeom_math::Frame, fb: ogeom_math::Frame| -> bool {
        fa.origin().distance(fb.origin()) <= tol.confusion()
            && fa.z().vector().dot(fb.z().vector()) >= 1.0 - tol.angular()
            && fa.x().vector().dot(fb.x().vector()) >= 1.0 - tol.angular()
    };
    match (a, b) {
        (S::Plane(x), S::Plane(y)) => frames(x.plane().frame(), y.plane().frame()),
        (S::Cylinder(x), S::Cylinder(y)) => {
            frames(x.cylinder().frame(), y.cylinder().frame())
                && (x.cylinder().radius() - y.cylinder().radius()).abs() <= tol.confusion()
        }
        (S::Cone(x), S::Cone(y)) => {
            frames(x.cone().frame(), y.cone().frame())
                && (x.cone().reference_radius() - y.cone().reference_radius()).abs()
                    <= tol.confusion()
                && (x.cone().half_angle() - y.cone().half_angle()).abs() <= tol.angular()
        }
        (S::Sphere(x), S::Sphere(y)) => {
            frames(x.sphere().frame(), y.sphere().frame())
                && (x.sphere().radius() - y.sphere().radius()).abs() <= tol.confusion()
        }
        (S::Torus(x), S::Torus(y)) => {
            frames(x.torus().frame(), y.torus().frame())
                && (x.torus().major_radius() - y.torus().major_radius()).abs() <= tol.confusion()
                && (x.torus().minor_radius() - y.torus().minor_radius()).abs() <= tol.confusion()
        }
        _ => false,
    }
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
    /// A sub-range of a pole of the face's own chart. The edge is a point in
    /// space whatever the range says; the range is where it runs in the
    /// chart, which is the only place it has length.
    Pole { pole: usize, range: (f64, f64) },
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
/// A point in a polyline's *interior*: its half-way member, except that a
/// straight strand is two points and its half-way member is an endpoint —
/// which may sit exactly on a seam or a trim. The chord midpoint is on the
/// curve for a straight image and strictly inside either way.
fn interior_of(line: &[Point2]) -> Point2 {
    if line.len() == 2 {
        Point2::new(
            f64::midpoint(line[0].x, line[1].x),
            f64::midpoint(line[0].y, line[1].y),
        )
    } else {
        line[line.len() / 2]
    }
}

fn fold_into_chart(line: &mut [Point2], surface: &SurfaceGeometry) {
    let ((ua, ub), (va, vb)) = surface.domain();
    if line.is_empty() {
        return;
    }
    let mid = interior_of(line);
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

/// The parameter to evaluate a section at: folded where the curve is closed
/// and its pieces may run past the domain end, left alone where it is not.
///
/// An open section's own end *is* the domain end, and folding it lands on the
/// domain start instead — the far end of the curve. On a plane's section
/// through a ball's poles, where the marcher hands back two open half circles
/// each running pole to pole, that is the difference between an edge bounded
/// by the pole it reaches and one bounded by the pole on the other side.
fn at_param(t: f64, domain: (f64, f64), closed: bool) -> f64 {
    if closed { fold(t, domain) } else { t }
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
    Vec<TangentRec>,
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
    let mut tangents: Vec<TangentRec> = Vec::new();
    let mut same_pairs: Vec<(usize, usize)> = Vec::new();
    // Marched sections are fitted; the fit is driven below the confusion
    // tolerance so a fitted curve meets edges, vertices and the mesh welder
    // on the same terms as an exact one. The budget each still carries is
    // recorded per section and widens the crossing filters.
    for (ia, fa) in ga.faces.iter().enumerate() {
        for (ib, fb) in gb.faces.iter().enumerate() {
            let admitted = fa.bound.intersects(&fb.bound);
            if !admitted && !*AUDIT_BOUNDS {
                // The faces cannot meet, whatever their surfaces do.
                continue;
            }
            let scale = fa.chord_scale.min(fb.chord_scale);
            let chord = (scale * 1e-7).max(tol.confusion() * 0.5);
            let options = IntersectOptions {
                // The fit budget scales with the chord: a section against a
                // fitted surface cannot honestly land closer than the
                // geometry it cuts, and whatever it carries is stated on
                // the record and widens every filter downstream.
                tolerance: chord,
                marching: ogeom_intersect::Marching {
                    chord,
                    ..ogeom_intersect::Marching::default()
                },
            };
            // Coincidence is asked before the intersector is, but only where
            // the closed forms have already declined the pair. The marcher
            // documents that it is not the one to answer it — it seeds on
            // sign changes, and a pair that never separates has none — so
            // what it traces over a coincident pair is noise wearing a
            // section's name, and it costs seconds to produce.
            let met = if ogeom_intersect::surface_surface(&fa.surface, &fb.surface, tol).is_err()
                && surfaces_coincide(&fa.surface, &fb.surface, options.tolerance, tol)
            {
                SurfaceIntersection::Same
            } else {
                intersect_surfaces(&fa.surface, &fb.surface, options, tol)?
            };
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
                            let (pcurve, prange) = match ogeom_intersect::exact_pcurve_of(
                                &e.curve,
                                &target.surface,
                                tol,
                            ) {
                                Some(exact) => (exact, e.crange),
                                // A fitted edge has no closed-form projection
                                // — but when the two faces sit on the
                                // *identical chart*, which is exactly the
                                // situation `Same` names for the analytics,
                                // the owner's own stored pcurve already is
                                // the projection, attached at construction,
                                // and it travels with its own window the way
                                // every stored pcurve does. A chart that
                                // merely coincides as a point set is still
                                // refused.
                                None if same_chart(&owner.surface, &target.surface, tol) => {
                                    (e.pcurve.clone(), e.prange)
                                }
                                None => ogeom_bail!(
                                    NotDone,
                                    "same-domain contact whose edges have no \
                                     closed-form projection into the shared \
                                     surface's chart is refused — see the \
                                     remaining work in docs/PLAN.md"
                                ),
                            };
                            contacts.push(ContactRec {
                                curve: e.curve.clone(),
                                crange: e.crange,
                                pcurve,
                                prange,
                                node: e.node,
                                tolerance: e.tolerance,
                                target_from_a,
                                target_face,
                            });
                        }
                    }
                }
                // A touch at a point bounds nothing: no curve to split a
                // face along, no side that is inside on one hand and
                // outside on the other. It contributes to the arrangement
                // exactly what a tangential contact curve does, which is
                // nothing, and the result carries the touch as the
                // non-manifold contact it is.
                SurfaceIntersection::Touching(_) => {}
                SurfaceIntersection::Along(curves) => {
                    for sc in curves {
                        // A tangential curve is contact, not crossing: the
                        // two faces meet along it and neither passes
                        // through the other, so it is set aside from the
                        // arrangement entirely and kept for the consumers
                        // that draw contact rather than classify by it.
                        if sc.tangential {
                            if let (Some(pa), Some(pb)) = (sc.on_a, sc.on_b) {
                                tangents.push(TangentRec {
                                    curve: sc.curve,
                                    pc_a: pa,
                                    pc_b: pb,
                                    face_a: ia,
                                    face_b: ib,
                                });
                            }
                            continue;
                        }
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
                                // A section running through a chart
                                // degeneracy — a plane cutting a ball on its
                                // own axis meets it at both poles — has no
                                // single chart image, because the longitude
                                // jumps half a turn there. Each piece
                                // *between* the poles does have one, and it
                                // is exact. Split first; march only if that
                                // fails.
                                if let Some(split) = split_at_degeneracies(&sc.curve, fa, fb, tol)?
                                {
                                    for (curve, pa, pb) in split {
                                        sections.push(SectionRec {
                                            curve,
                                            pc_a: pa,
                                            pc_b: pb,
                                            face_a: ia,
                                            face_b: ib,
                                            closed: false,
                                            tolerance: 0.0,
                                        });
                                    }
                                    continue;
                                }
                                // An exact curve whose projection has no
                                // closed form: march the pair instead, so
                                // curve and pcurves are fitted *together*.
                                let shared = if admitted {
                                    fa.bound.intersection(&fb.bound)
                                } else {
                                    // Audit only: disjoint bounds have no
                                    // window, and an empty window would mask
                                    // the very miss being hunted.
                                    fa.bound.union(&fb.bound)
                                };
                                for fitted in march_pair(
                                    &windowed_to(&fa.surface, &shared),
                                    &windowed_to(&fb.surface, &shared),
                                    &options,
                                    tol,
                                )? {
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
    // Each section's paving depends only on the sections and the two
    // gathered solids, all read-only here, and writes nothing the next
    // section reads. So the measuring runs in parallel and the accumulating
    // runs afterwards in section order — the same split `tessellate` uses,
    // and the same reason: nothing about scheduling can reach the answer.
    type SectionWork = (Vec<(ogeom_topo::TShapeId, f64)>, Vec<SectionPiece>);
    let paved: Vec<OgeomResult<SectionWork>> =
        ogeom_core::parallel::map_ordered(&sections, |si, section: &SectionRec| {
            ogeom_core::progress::checkpoint()?;
            let mut paves: Vec<(ogeom_topo::TShapeId, f64)> = Vec::new();
            let mut pieces: Vec<SectionPiece> = Vec::new();
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
            // so the span is excluded from that face's strands rather than
            // refused or duplicated.
            //
            // Which face, though, is the whole point. A plane through a bore's
            // axis meets the wall along two rulings, and one of them is the
            // wall's own seam: on the wall that ruling is boundary already, and
            // on the plane it is the curve that separates the two halves the
            // section leaves. Dropped from both, the plane keeps one region where
            // it has two, and the result does not close. So the exclusion is
            // recorded per face.
            let mut along: [Vec<(f64, f64)>; 2] = [Vec::new(), Vec::new()];
            for (side, own) in [
                (0_usize, &ga.faces[section.face_a]),
                (1, &gb.faces[section.face_b]),
            ] {
                for e in &own.edges {
                    let found = intersect_curves(&section.curve, &e.curve, cc, tol)?;
                    for crossing in &found.crossings {
                        if crossing.gap > reach {
                            continue;
                        }
                        let on_b = onto_range(crossing.on_b, &e.curve, e.crange, tol);
                        // A crossing at a boundary edge's own end *is* that end's
                        // vertex, exactly. The stop the section keeps must be the
                        // vertex's parameter on the section — the meet of two
                        // curves a fit tolerance apart lands a couple of microns
                        // off, the boundary side keeps its exact vertex, and the
                        // rebuilt wire gapes by the difference.
                        let mut on_a = crossing.on_a;
                        // The window is the fit-slop scale, not the section's
                        // own: the boundary curve may be a fitted intersection
                        // from an earlier boolean carrying a couple of microns
                        // of wobble, and a stop that misses the vertex by that
                        // much gapes the rebuilt wire by the same.
                        let weld = reach.max(tol.confusion() * 1e2);
                        for end in [e.crange.0, e.crange.1] {
                            let vertex = e.curve.point_at(end, tol)?;
                            if vertex.distance(crossing.point) <= weld {
                                let snapped =
                                    ogeom_algo::project_on_curve(&section.curve, vertex, 64, tol)?;
                                if snapped.distance <= weld {
                                    on_a = snapped.parameter;
                                }
                                break;
                            }
                        }
                        trim_ts.push(on_a);
                        edge_hits.push((e.node, on_b, on_a));
                    }
                    for overlap in &found.overlaps {
                        // The curves overlap; what is *boundary* is the stretch
                        // the edge actually covers. A sphere's seam and the far
                        // half of the same great circle lie on one curve, and
                        // reading the whole curve as boundary makes the meridian
                        // opposite the seam disappear — which is the octant's own
                        // edge, on a ball cut at its corner.
                        let Some((lo, hi)) = overlap_within(overlap, e.crange, &e.curve, tol)
                        else {
                            continue;
                        };
                        trim_ts.push(lo);
                        trim_ts.push(hi);
                        along[side].push((lo, hi));
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
                if !inside_both(mid)? {
                    continue;
                }
                // A section that runs along a boundary edge of a face splits
                // nothing *there*: the split already exists as boundary. The
                // analytic overlap detection above catches the same-support
                // cases; this catches the rest — a fitted section tracing a
                // boundary curve, a surface meeting another exactly at its own
                // trim — by measurement rather than by recognising supports.
                let mut hugs = [false; 2];
                for (side, own, side_from_a, side_face) in [
                    (0_usize, &ga.faces[section.face_a], true, section.face_a),
                    (1, &gb.faces[section.face_b], false, section.face_b),
                ] {
                    if along[side]
                        .iter()
                        .any(|(alo, ahi)| mid_folded >= *alo && mid_folded <= *ahi)
                    {
                        hugs[side] = true;
                        continue;
                    }
                    let mut all_near = true;
                    for i in 0..=4 {
                        let t = lo + (hi - lo) * f64::from(i) / 4.0;
                        let tf = if section.closed { fold(t, domain) } else { t };
                        let at = section.curve.point_at(tf, tol)?;
                        // Wider than the crossing filters on purpose: a
                        // tangentially-traced curve wobbles about the boundary it
                        // hugs by far more than a fit budget, and a genuine
                        // section keeps a distance of feature scale, not microns.
                        let width = reach.max(tol.confusion() * 1e3);
                        let mut near = false;
                        for e in &own.edges {
                            if distance_to_edge_curve(&e.curve, e.crange, at, tol)? <= width {
                                near = true;
                                break;
                            }
                        }
                        // A contact edge is a strand on this face too — the other
                        // solid's boundary, carried into a chart they share — so a
                        // section tracing one would be the same curve twice, and
                        // the arrangement cannot walk a line it meets from both
                        // sides at once. Two boxes side by side put the low one's
                        // lid exactly there.
                        if !near {
                            for c in &contacts {
                                if c.target_from_a != side_from_a || c.target_face != side_face {
                                    continue;
                                }
                                if distance_to_edge_curve(&c.curve, c.crange, at, tol)? <= width {
                                    near = true;
                                    break;
                                }
                            }
                        }
                        if !near {
                            all_near = false;
                            break;
                        }
                    }
                    hugs[side] = all_near;
                }
                if hugs[0] && hugs[1] {
                    // Boundary on both sides: the split exists twice over and
                    // adding it a third time would cancel what it copies.
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
                        paves.push((*node, *on_edge));
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
                            hugs,
                        });
                        pieces.push(SectionPiece {
                            section: si,
                            range: (mid, hi2),
                            hugs,
                        });
                    } else {
                        pieces.push(SectionPiece {
                            section: si,
                            range: (lo2, hi2),
                            hugs,
                        });
                    }
                }
            }
            Ok((paves, pieces))
        });
    for work in paved {
        let (found, made) = work?;
        for (node, at) in found {
            paves.entry(node).or_default().push(at);
        }
        pieces.extend(made);
    }
    // The audit's verdict. A dropped pair whose surfaces intersect is not a
    // filter bug — planes meet along an infinite line the paving then trims
    // to the faces, usually to nothing. A dropped pair whose section
    // *survives paving* is: some of that curve lies inside both faces, so
    // the faces genuinely meet and the filter's boxes failed to. Contacts
    // need no check — a contact is an owner edge lying in the target face,
    // which forces the boxes to overlap where the edge does.
    //
    // Evidence, not proof: the audit paves with every pair admitted, and
    // paving is not compositional — a section's kept intervals see the other
    // sections' paves — so a filtered run is not replayed exactly. It has
    // already earned its keep the other way, by clearing a suspected filter
    // miss and pointing the hunt at the real coupling (the marching chord
    // fed from the filter box).
    if *AUDIT_BOUNDS {
        for piece in &pieces {
            let section = &sections[piece.section];
            let (fa, fb) = (&ga.faces[section.face_a], &gb.faces[section.face_b]);
            assert!(
                fa.bound.intersects(&fb.bound),
                "bound filter audit: faces {}/{} were dropped by the bound \
                 filter, yet their section paved a surviving piece over \
                 {:?}; the filter under-approximates",
                section.face_a,
                section.face_b,
                piece.range,
            );
        }
    }
    // A contact edge splits where it crosses the target face's boundary —
    // and that crossing is a pave on *both* edges, so the owner's own face
    // splits its boundary consistently and the pieces sew back shared.
    let mut contact_along: Vec<Vec<(f64, f64)>> = vec![Vec::new(); contacts.len()];
    for (ci, contact) in contacts.iter().enumerate() {
        let target = if contact.target_from_a {
            &ga.faces[contact.target_face]
        } else {
            &gb.faces[contact.target_face]
        };
        // A fitted contact sits off exact geometry by its own tolerance, and
        // the crossings it genuinely makes gape by the same — the filter and
        // the finder both widen, or a fitted rail never registers against
        // the exact seam it crosses.
        let reach = tol.confusion().max(contact.tolerance * 2.0);
        let cc = CurveCurveOptions {
            gap: reach.max(CurveCurveOptions::default().gap),
            ..CurveCurveOptions::default()
        };
        for e in &target.edges {
            let found = intersect_curves(&contact.curve, &e.curve, cc, tol)?;
            for crossing in &found.crossings {
                if crossing.gap > reach {
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
                let on_b = onto_range(crossing.on_b, &e.curve, e.crange, tol);
                if on_b > e.crange.0 + tol.parametric() && on_b < e.crange.1 - tol.parametric() {
                    paves.entry(e.node).or_default().push(on_b);
                }
            }
            // A span of the contact running along a target boundary edge
            // splits nothing: it is already boundary on both sides —
            // identically stacked boxes are all such spans — and duplicating
            // it as a strand would cancel the boundary it copies.
            for overlap in &found.overlaps {
                let ordered = |r: (f64, f64)| if r.0 <= r.1 { r } else { (r.1, r.0) };
                let (lo, hi) = ordered(overlap.on_a);
                // The overlap is between the two *curves*; what interferes is
                // the stretch both *edges* actually cover. A hole's arc and
                // the disc that fills it lie on one circle, so the curves
                // overlap over the whole turn while the arc covers three
                // quarters of it — and paving at the turn's ends says nothing,
                // where paving at the arc's ends is exactly the split the
                // other side needs to sew against.
                let (lo, hi) = (lo.max(contact.crange.0), hi.min(contact.crange.1));
                if hi - lo <= tol.parametric() {
                    continue;
                }
                paves.entry(contact.node).or_default().push(lo);
                paves.entry(contact.node).or_default().push(hi);
                contact_along[ci].push((lo, hi));
                // The *target* edge splits where the shared stretch ends,
                // exactly as the contact does. Without this, the face across
                // the overlap keeps one long boundary edge where its new
                // neighbours carry two short ones, and sew — which matches
                // edges whole — can pair it with neither. The clamped ends are
                // carried across by the correspondence the overlap itself
                // states, which is affine over the shared stretch.
                let span = overlap.on_a.1 - overlap.on_a.0;
                let carry = |t: f64| -> f64 {
                    if span.abs() <= f64::MIN_POSITIVE {
                        overlap.on_b.0
                    } else {
                        overlap.on_b.0
                            + (overlap.on_b.1 - overlap.on_b.0) * (t - overlap.on_a.0) / span
                    }
                };
                let (mut tlo, mut thi) = (carry(lo), carry(hi));
                if tlo > thi {
                    core::mem::swap(&mut tlo, &mut thi);
                }
                let target_domain = e.curve.domain();
                let periodic = e.curve.is_periodic();
                for t in [tlo, thi] {
                    // A correspondence across a full turn can run the
                    // parameter past the domain; the pave belongs where the
                    // edge actually is.
                    let t = if periodic { fold(t, target_domain) } else { t };
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
        tangents,
        contact_along,
        paves,
        same_a,
        same_b,
    ))
}

/// An exact section cut at the chart degeneracies it runs through, each
/// piece carrying exact pcurves on both faces.
///
/// The degeneracies are not guessed from the surfaces: they are the faces'
/// own *pole edges*, the degenerate edges the topology already carries, so a
/// face whose chart has a pole says so and one that has none costs nothing.
///
/// `None` means the split is no use here — no degeneracy on the curve, or a
/// piece whose projection still has no closed form — and the caller falls
/// back to marching, which is the honest answer rather than a fitted pcurve
/// pretending to be exact.
fn split_at_degeneracies(
    curve: &Curve,
    fa: &GFace,
    fb: &GFace,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<(Curve, PlanarCurve, PlanarCurve)>>> {
    let mut stops: Vec<f64> = Vec::new();
    let domain = curve.domain();
    // The split points are the *surfaces'* chart degeneracies, not merely the
    // faces' pole edges: a meridian section runs through both of a sphere's
    // poles, and a face that owns only the north one still cannot chart an
    // arc that wraps through the south. The surface knows where its chart
    // collapses whether or not the face's trim reaches there.
    let mut candidates: Vec<Point> = Vec::new();
    for face in [fa, fb] {
        for pole in &face.poles {
            candidates.push(pole.point);
        }
        match &face.surface {
            SurfaceGeometry::Sphere(s) => {
                let sphere = s.sphere();
                let axis = sphere.frame().z().vector();
                candidates.push(sphere.centre() + axis * sphere.radius());
                candidates.push(sphere.centre() - axis * sphere.radius());
            }
            SurfaceGeometry::Cone(c) => {
                candidates.push(c.cone().apex());
            }
            _ => {}
        }
    }
    for point in candidates {
        let found = ogeom_algo::project_on_curve(curve, point, 256, tol)?;
        if found.distance <= tol.confusion() {
            stops.push(found.parameter);
        }
    }
    if stops.is_empty() {
        return Ok(None);
    }
    stops.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    stops.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());

    // A closed curve is cut into the arcs between successive stops, the last
    // wrapping past the domain end; an open one keeps its own ends as stops.
    let closed = curve.is_closed(tol) || curve.is_periodic();
    let mut arcs: Vec<(f64, f64)> = Vec::new();
    if closed {
        if stops.len() < 2 {
            // One stop on a closed curve leaves one arc, from the stop right
            // round to itself — which still has the pole at both ends.
            let period = domain.1 - domain.0;
            arcs.push((stops[0], stops[0] + period));
        } else {
            let period = domain.1 - domain.0;
            for i in 0..stops.len() {
                let lo = stops[i];
                let hi = if i + 1 < stops.len() {
                    stops[i + 1]
                } else {
                    stops[0] + period
                };
                arcs.push((lo, hi));
            }
        }
    } else {
        let mut cuts = vec![domain.0];
        cuts.extend(
            stops
                .iter()
                .copied()
                .filter(|t| *t > domain.0 + tol.parametric() && *t < domain.1 - tol.parametric()),
        );
        cuts.push(domain.1);
        for pair in cuts.windows(2) {
            arcs.push((pair[0], pair[1]));
        }
    }

    let mut out = Vec::with_capacity(arcs.len());
    for (lo, hi) in arcs {
        if hi - lo <= tol.parametric() {
            continue;
        }
        let (Some(pa), Some(pb)) = (
            ogeom_intersect::exact_pcurve_over(curve, (lo, hi), &fa.surface, tol),
            ogeom_intersect::exact_pcurve_over(curve, (lo, hi), &fb.surface, tol),
        ) else {
            return Ok(None);
        };
        let Ok(piece) = ogeom_geom::TrimmedCurve::new(curve.clone(), lo, hi, tol) else {
            return Ok(None);
        };
        out.push((piece.into(), pa, pb));
    }
    if out.is_empty() {
        return Ok(None);
    }
    Ok(Some(out))
}

/// A surface narrowed to the reach of a bound, for the marcher's benefit.
///
/// Seeding samples a surface's *parameter box*, so a plane stored over
/// ±10^9 — which is how an unbounded carrier reaches this code — is sampled
/// at a spacing a hundred million times anything it could be meeting, and
/// no seed ever lands near the curve. Narrowing to where the two faces
/// actually are is what lets the seeding see it. The chart is untouched:
/// only the window the marcher walks changes, so pcurves fitted in it mean
/// the same thing on the surface as stored.
fn windowed_to(surface: &SurfaceGeometry, bound: &ogeom_math::Aabb) -> SurfaceGeometry {
    let corners = bound.corners();
    if corners.is_empty() {
        return surface.clone();
    }
    match surface {
        SurfaceGeometry::Plane(p) => {
            let frame = p.plane().frame();
            let (mut u0, mut u1, mut v0, mut v1) = (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            );
            for c in &corners {
                let local = frame.to_local(*c);
                u0 = u0.min(local.x);
                u1 = u1.max(local.x);
                v0 = v0.min(local.y);
                v1 = v1.max(local.y);
            }
            let margin = ((u1 - u0) + (v1 - v0)).mul_add(0.25, 1.0);
            let (want_u, want_v) = ((u0 - margin, u1 + margin), (v0 - margin, v1 + margin));
            let (have_u, have_v) = ogeom_geom::Surface::domain(p);
            if have_u.1 - have_u.0 <= want_u.1 - want_u.0
                && have_v.1 - have_v.0 <= want_v.1 - want_v.0
            {
                return surface.clone();
            }
            ogeom_geom::PlaneSurface::over(p.plane(), want_u, want_v)
                .map_or_else(|_| surface.clone(), Into::into)
        }
        SurfaceGeometry::Cylinder(c) => {
            let frame = c.cylinder().frame();
            let (mut h0, mut h1) = (f64::INFINITY, f64::NEG_INFINITY);
            for corner in &corners {
                let h = (*corner - frame.origin()).dot(frame.z().vector());
                h0 = h0.min(h);
                h1 = h1.max(h);
            }
            let margin = (h1 - h0).mul_add(0.25, 1.0);
            let want = (h0 - margin, h1 + margin);
            let have = ogeom_geom::Surface::domain(c).1;
            if have.1 - have.0 <= want.1 - want.0 {
                return surface.clone();
            }
            ogeom_geom::CylinderSurface::new(c.cylinder(), want)
                .map_or_else(|_| surface.clone(), Into::into)
        }
        other => other.clone(),
    }
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
    // One seeding grid over the far surface, asked a hundred times: the
    // same seeds and the same Newton the per-call projection would use, so
    // the verdict is bit-identical, at hundreds of evaluations instead of
    // tens of thousands — the price issue #26 named for every genuine
    // near-miss.
    let Ok(seeds) = ogeom_algo::SurfaceSeeds::over(against, 16, tol) else {
        return false;
    };
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
            let Ok(projection) = seeds.project(against, p, tol) else {
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

/// Whether two surfaces are one surface wherever they overlap, measured.
///
/// The closed forms answer this for the pairs they know. Where there is no
/// closed form it does not stop being a fair question — two patches restated
/// from one plane are the same surface, and nothing in their control points
/// says so — but it stops being answerable exactly, so it is measured here
/// and only for the pairs the analytic layer has already declined.
///
/// Sampled on the *smaller* window, because the answer is about the region
/// the two share and a stated window is not that region: a plane's own
/// extends for a billion units either way, and a grid over it samples
/// nothing. A sample whose foot lands on the rim of the other is skipped
/// rather than counted against — the other patch simply does not reach that
/// far, and a distance measured to its rim is about the window, not the
/// surface.
///
/// One-sided by construction: a pair that crosses puts interior samples well
/// off the other, so it cannot pass, and a pair this cannot resolve marches
/// exactly as it did before.
fn surfaces_coincide(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    reach: f64,
    tol: Tolerances,
) -> bool {
    use ogeom_geom::Surface as _;
    /// Samples per direction over the window, and how many must land inside
    /// the other before agreement means anything.
    const GRID: usize = 6;
    const EVIDENCE: usize = 4;

    let span = |s: &SurfaceGeometry| -> f64 {
        let ((ua, ub), (va, vb)) = s.domain();
        (ub - ua).abs().max((vb - va).abs())
    };
    let (sampled, against) = if span(a) <= span(b) { (a, b) } else { (b, a) };
    let ((ua, ub), (va, vb)) = sampled.domain();
    if !(ua.is_finite() && ub.is_finite() && va.is_finite() && vb.is_finite()) {
        return false;
    }
    let ((wu0, wu1), (wv0, wv1)) = against.domain();
    // A twentieth of the window in from each rim: enough that a foot the
    // search pinned to the rim is not read as one the surface truly reaches.
    let (mu, mv) = ((wu1 - wu0) * 0.05, (wv1 - wv0) * 0.05);

    let mut evidence = 0_usize;
    for i in 0..=GRID {
        for j in 0..=GRID {
            #[allow(clippy::cast_precision_loss)]
            let u = ua + (ub - ua) * (i as f64 / GRID as f64);
            #[allow(clippy::cast_precision_loss)]
            let v = va + (vb - va) * (j as f64 / GRID as f64);
            let Ok(p) = sampled.point_at(u, v, tol) else {
                return false;
            };
            let Ok(foot) = ogeom_algo::project_on_surface(against, p, 16, tol) else {
                return false;
            };
            let (fu, fv) = foot.parameters;
            if fu <= wu0 + mu || fu >= wu1 - mu || fv <= wv0 + mv || fv >= wv1 - mv {
                continue;
            }
            if foot.distance > reach {
                return false;
            }
            evidence += 1;
        }
    }
    evidence >= EVIDENCE
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
    /// bound the union and the intersection once between them, not twice —
    /// but "once" is a statement about the *patch*, not about which argument
    /// it came from. Where a tool's cap fills a hole its own bore left, the
    /// part has no piece there at all, and dropping the tool's by argument
    /// identity leaves the result with a hole in it. So the duplicate is
    /// identified by containment, and only a piece genuinely stood in for is
    /// dropped.
    covered: bool,
}

struct GeneralFused {
    a: GSolid,
    b: GSolid,
    sections: Vec<SectionRec>,
    contacts: Vec<ContactRec>,
    /// Curves the two boundaries touch along without crossing. They take no
    /// part in the classification below — that is what tangency means — and
    /// are carried for the consumers that want the contact itself.
    tangents: Vec<TangentRec>,
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
    // The face's boundary polylines are unwrapped — a winding ring may span
    // any one period's window, not necessarily the chart's canonical one —
    // so the probe is asked at every period image that could land inside.
    let mut shifts = vec![0.0];
    if face.surface.is_periodic_u() {
        let ((ua, ub), _) = face.surface.domain();
        if ub > ua {
            shifts.push(ub - ua);
            shifts.push(ua - ub);
        }
    }
    for shift in shifts {
        let shifted = Point2::new(at.x + shift, at.y);
        if inside_many(&borrowed, shifted) {
            return Some(shifted);
        }
    }
    None
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
    let (best, _) = scan((t2 - fine).max(crange.0), (t2 + fine).min(crange.1), 16)?;
    Ok(best)
}

/// Pair up the pieces two coincident faces contribute for the same patch of
/// one surface, and mark the second argument's copy as stood in for.
///
/// The substitution is by *region*, not by argument: a piece of `b` is a
/// duplicate exactly where a piece of `a`, on the same side and on a surface
/// they share, already covers the point `b`'s piece was classified at. Where
/// `a` has nothing there — a bore refilled by the cylinder that cut it, whose
/// caps fill holes the part no longer has faces for — nothing stands in, and
/// `b`'s piece is the only description of that patch there is.
fn mark_covered_coincidences(ga: &GSolid, pieces: &mut [FacePiece], tol: Tolerances) {
    // Which pieces of the first argument stand on a shared surface.
    let from_a: Vec<(usize, usize)> = pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.from_a && p.state == PieceState::OnAligned)
        .map(|(i, p)| (i, p.face))
        .collect();
    if from_a.is_empty() {
        return;
    }
    let mut covered: Vec<usize> = Vec::new();
    for (index, piece) in pieces.iter().enumerate() {
        if piece.from_a || piece.state != PieceState::OnAligned {
            continue;
        }
        for &(other, face_a) in &from_a {
            let host = &ga.faces[face_a];
            // The point `b`'s piece stands at, read in `a`'s face's chart. A
            // point that is not on that surface at all has no chart position,
            // and `chart_point_of` says so.
            let Some(at) = chart_point_of(host, piece.probe, tol) else {
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
                covered.push(index);
                break;
            }
        }
    }
    for index in covered {
        pieces[index].covered = true;
    }
}

/// The debug dumps, read once rather than once per face.
///
/// `env::var` takes a process-wide lock and allocates; the strand dump asked
/// it inside the per-face loop, where a large model asks thousands of times
/// to be told no.
/// The face-bound filter's audit: admit every pair, and name any the filter
/// would have dropped that then produces a record. A conservative filter is
/// a correctness precondition — a pair wrongly dropped is absorbed by the
/// empty-result fallback today, and would be a wrong solid if that fallback
/// ever came up empty too — and this is the check that makes the
/// precondition falsifiable (issue #12). Costs one branch per pair when off.
static AUDIT_BOUNDS: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("OGEOM_BOOL_AUDIT_BOUNDS").is_ok());
static DEBUG_STRANDS: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("OGEOM_DEBUG_STRANDS").is_ok());
static ARRANGE_DEBUG: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("OGEOM_ARRANGE_DEBUG").is_ok());

fn general_fuse(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<GeneralFused> {
    ogeom_core::progress::stage("boolean: gather");
    let ga = gather(model, a, tol)?;
    let gb = gather(model, b, tol)?;
    ogeom_core::progress::stage("boolean: intersect");
    let (sections, section_pieces, contacts, tangents, contact_along, paves, same_a, same_b) =
        fill(&ga, &gb, tol)?;

    ogeom_core::progress::stage("boolean: split");
    let mut pieces: Vec<FacePiece> = Vec::new();
    for (from_a, own, other) in [(true, &ga, &gb.solid), (false, &gb, &ga.solid)] {
        // The other solid's boundary, prepared once for the whole side. It is
        // asked once per face piece, and what it costs to prepare — every
        // face's trimming rings, polylined — does not depend on the point
        // being asked about. Rebuilt per question it dwarfed the question:
        // 3.5 ms of preparation against 5.6 µs of ray casting.
        let boundary = ogeom_algo::SolidBoundary::of(model, other, tol.confusion() * 1e4, tol)?;
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
            // The poles, after the sections, because a section can end *on*
            // a pole — a plane through a ball's axis cuts it exactly there —
            // and the pole has to be cut where that happens or the two meet
            // at no shared node and the arrangement sees a dangling section.
            for (pi, pole) in face.poles.iter().enumerate() {
                let mut stops = vec![pole.prange.0, pole.prange.1];
                if let PlanarCurve::Line(line) = &pole.pcurve {
                    let axis = line.axis();
                    for strand in &strands {
                        if strand.boundary {
                            continue;
                        }
                        for end in [strand.polyline.first(), strand.polyline.last()]
                            .into_iter()
                            .flatten()
                        {
                            let along = (*end - axis.location).dot(axis.direction.vector());
                            let foot = axis.point_at(along);
                            if foot.distance(*end) <= PARAM_SNAP
                                && along > pole.prange.0 + PARAM_SNAP
                                && along < pole.prange.1 - PARAM_SNAP
                            {
                                stops.push(along);
                            }
                        }
                    }
                }
                stops.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
                stops.dedup_by(|a, b| (*a - *b).abs() <= PARAM_SNAP);
                for pair in stops.windows(2) {
                    let sub = (pair[0], pair[1]);
                    strands.push(Strand {
                        polyline: pcurve_polyline(
                            &pole.pcurve,
                            pole.prange,
                            pole.prange,
                            sub,
                            tol,
                        )?,
                        tag: Tag::Pole {
                            pole: pi,
                            range: sub,
                        },
                        boundary: true,
                    });
                }
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
                    let mut line =
                        pcurve_polyline(&contact.pcurve, contact.prange, contact.crange, sub, tol)?;
                    unwrap_polyline(&mut line, &face.surface);
                    fold_into_chart(&mut line, &face.surface);
                    let mid = interior_of(&line);
                    let boundary_lines: Vec<&[Point2]> = strands
                        .iter()
                        .filter(|st| st.boundary)
                        .map(|st| st.polyline.as_slice())
                        .collect();
                    if !inside_many_slanted(&boundary_lines, mid) {
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

            // A tolerant contact's chart image meets the boundary it paved
            // only as closely as its own slop allows; the arrangement's node
            // weld reaches that far on this face, or the strand dangles a
            // few microns from the junction it belongs to.
            let face_snap = contacts
                .iter()
                .filter(|c| c.target_from_a == from_a && c.target_face == fi)
                .fold(PARAM_SNAP, |acc, c| acc.max(c.tolerance * 2.0))
                .max(
                    face.edges
                        .iter()
                        .fold(0.0_f64, |acc, e| acc.max(e.tolerance * 2.0)),
                );
            if *DEBUG_STRANDS {
                for (si, st) in strands.iter().enumerate() {
                    let tag = match st.tag {
                        Tag::Boundary { edge, range } => format!("Boundary e{edge} {range:?}"),
                        Tag::Contact { contact, range } => format!("Contact c{contact} {range:?}"),
                        Tag::Section { section, range } => format!("Section s{section} {range:?}"),
                        Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                    };
                    let (a, b) = (st.polyline[0], st.polyline[st.polyline.len() - 1]);
                    eprintln!(
                        "STRAND from_a={from_a} fi={fi} {si}: boundary={} {tag} pts={} {a:?} .. {b:?}",
                        st.boundary,
                        st.polyline.len()
                    );
                }
            }
            let split = arrange_pieces(&strands, face_snap)?;
            for piece in split {
                // Where a piece stands is asked at its interior probes in
                // turn. The first is the roomiest, and usually the only one
                // needed; the rest are for the piece that merely *touches*
                // the other solid, whose roomiest probe can land on the
                // contact and read neither in nor out.
                //
                // A partner face is asked before the classifier, not after.
                // Where a piece sits on a surface the other solid also
                // carries, the partner *is* the answer — and getting there
                // through the classifier means every ray grazing the shared
                // face, its whole fan of directions exhausted, before it
                // reports the On the partner list already knew. On a part
                // whose bore is refilled by its own cylinder that is the
                // difference between a tenth of a second and a minute.
                let partners = if from_a { &same_a[fi] } else { &same_b[fi] };
                let mut chosen = None;
                for candidate in &piece.interiors {
                    let at = face.surface.point_at(candidate.x, candidate.y, tol)?;
                    let shared = !partners.is_empty()
                        && partners.iter().any(|&pi| {
                            let partner = if from_a { &gb.faces[pi] } else { &ga.faces[pi] };
                            chart_point_of(partner, at, tol).is_some()
                        });
                    let says = if shared {
                        Containment::On
                    } else {
                        boundary.holds(model, at, tol)?
                    };
                    if chosen.is_none() || !matches!(says, Containment::On) {
                        chosen = Some((*candidate, at, says));
                    }
                    if !matches!(says, Containment::On) {
                        break;
                    }
                }
                let Some((interior, probe, said)) = chosen else {
                    ogeom_bail!(
                        Construction,
                        "a piece of a face has no interior point to classify at"
                    );
                };
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
                                        outlines: piece.outlines,
                                        probe,
                                        state: PieceState::In,
                                        covered: false,
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
                                 remaining work in docs/PLAN.md"
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
                });
            }
        }
    }
    mark_covered_coincidences(&ga, &mut pieces, tol);
    if *ARRANGE_DEBUG {
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
            eprintln!(
                "PIECE {i} from_a={} face={} {kind} state={:?} covered={} probe={:?}",
                p.from_a, p.face, p.state, p.covered, p.probe
            );
        }
    }
    ogeom_core::progress::stage("boolean: classified");
    Ok(GeneralFused {
        a: ga,
        b: gb,
        sections,
        contacts,
        tangents,
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
    /// How far two honest descriptions of one junction may sit apart: a
    /// hundred confusions as the floor, widened to twice the loosest contact
    /// edge's own tolerance when a fitted rail took part in the melt.
    weld: f64,
}

impl Rebuild<'_> {
    fn vertex(&mut self, p: Point, tol: Tolerances) -> Shape {
        // The weld reach covers what the inputs may honestly disagree by: a
        // boundary curve that is itself a fitted intersection from an
        // earlier boolean carries a couple of microns of slop, and two
        // descriptions of one junction arrive that far apart. A hundred
        // confusions at millimetre tolerances is ten microns — below any
        // feature this pipeline can resolve, and every weld wider than
        // confusion is recorded on the vertex, not papered over.
        let found = self
            .vertices
            .iter()
            .find(|(q, _)| q.distance(p) <= self.weld.max(tol.confusion() * 1e2))
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
                Tag::Pole { pole, range } => (*pole, 3, *range),
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
            // The stored image keeps its own window; the attached copy names
            // the sub-window this piece covers under the proportional map.
            let sub_p = (
                rescale(range.0, c.crange, c.prange),
                rescale(range.1, c.crange, c.prange),
            );
            let mid = c.pcurve.point_at(f64::midpoint(sub_p.0, sub_p.1), tol)?;
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
                sub_p,
            )?;
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
    ogeom_core::progress::stage("boolean: assemble");
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
        weld: fused
            .contacts
            .iter()
            .fold(0.0_f64, |acc, c| acc.max(c.tolerance * 2.0)),
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
            // Env-gated forensics: the open shell's unshared edges, the
            // question every failure here starts from.
            if *ARRANGE_DEBUG {
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
                        eprintln!(
                            "  open edge uses={occurrences} faces={users} {:?} range={range:?}: {a:?} -> {b:?}",
                            core::mem::discriminant(g)
                        );
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

/// Solids from an unordered soup of faces: sew, demand closure, nest
/// shells into solids and voids — the pipeline's own final stages offered
/// as a builder.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// shell fails to close — an open soup encloses no volume, and saying so
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
    let mut bounds = Vec::new();
    for shell in &sewn.shells {
        bounds.push(ogeom_algo::shape_bounds(model, shell, tol)?);
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
                // A void bounds its solid from inside: material lies outside
                // it, so the sewn outward orientation reverses.
                group.push(candidate.reversed());
            }
        }
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
    Ok(Cells {
        a_not_b: cut(model, a, b, tol)?,
        b_not_a: cut(model, b, a, tol)?,
        common: common(model, a, b, tol)?,
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
    let loosened = fuzzed(fuzz, tol)?;
    fuse(model, a, b, loosened)
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
    let loosened = fuzzed(fuzz, tol)?;
    cut(model, a, b, loosened)
}

/// A shape repeated `count` times along a direction at a period and fused
/// into one — the periodic pattern as a composition of what exists.
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
    let mut history = History::new();
    let mut result = shape.clone();
    for i in 1..count {
        #[allow(clippy::cast_precision_loss, reason = "pattern counts are small")]
        let offset = step * i as f64;
        let moved =
            ogeom_algo::transformed(model, shape, ogeom_math::Transform::translation(offset))?
                .shape;
        let joined = fuse(model, &result, &moved, tol)?;
        history.modify(&moved, joined.shape.clone());
        result = joined.shape;
    }
    history.generate(shape, result.clone());
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
    let mut restate = false;
    for face in ogeom_topo::explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let placement = face.transform(model.datums())?;
        // A scale changes lengths the melt compares; a reflection flips
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
            p.state == PieceState::Out || (p.state == PieceState::OnAligned && !p.covered)
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
            p.state == PieceState::In || (p.state == PieceState::OnAligned && !p.covered)
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
        let from = s.curve.point_at(at_param(f0, domain, s.closed), tol)?;
        let to = s.curve.point_at(at_param(f1, domain, s.closed), tol)?;
        let v0 = make_vertex(model, from).shape;
        let v1 = make_vertex(model, to).shape;
        edges.push(make_edge_between(model, s.curve.clone(), (f0, f1), &v0, &v1, tol)?.shape);
    }
    // Contacts are not crossings, so no piece's ring carries them and the
    // loop above cannot see them — but a section through a tangency has a
    // curve in it, and this is where it comes from.
    for contact in &fused.tangents {
        for (lo, hi) in contact_intervals(&fused, contact, tol)? {
            let from = contact.curve.point_at(lo, tol)?;
            let to = contact.curve.point_at(hi, tol)?;
            let v0 = make_vertex(model, from).shape;
            let v1 = if from.distance(to) <= tol.confusion() {
                v0.clone()
            } else {
                make_vertex(model, to).shape
            };
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
    fn coincidence_is_measured_over_the_overlap_and_nowhere_else() {
        // Patches restated from planes: the geometry no longer says "plane",
        // which is the whole reason this measurement exists.
        let patch = |plane: ogeom_math::Plane, u: (f64, f64), v: (f64, f64)| {
            let surface: SurfaceGeometry =
                ogeom_geom::PlaneSurface::over(plane, u, v).unwrap().into();
            SurfaceGeometry::from(surface.to_bspline(T).unwrap())
        };
        let reach = T.confusion() * 1e2;

        // Two windows on one plane, overlapping over a quarter of each. They
        // are the same surface exactly where they meet, which is the claim.
        let here = patch(ogeom_math::Plane::XY, (0.0, 10.0), (0.0, 10.0));
        let over = patch(ogeom_math::Plane::XY, (5.0, 15.0), (5.0, 15.0));
        assert!(surfaces_coincide(&here, &over, reach, T));

        // The same plane lifted clear of itself is not the same surface, and
        // a plane square to it crosses rather than coincides — the case that
        // must keep marching, since a crossing has a section to find.
        let above = patch(
            ogeom_math::Plane::new(
                Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap(),
            ),
            (0.0, 10.0),
            (0.0, 10.0),
        );
        assert!(!surfaces_coincide(&here, &above, reach, T));
        let across = patch(
            ogeom_math::Plane::new(
                Frame::new(Point::new(5.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
            ),
            (0.0, 10.0),
            (0.0, 10.0),
        );
        assert!(!surfaces_coincide(&here, &across, reach, T));
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
