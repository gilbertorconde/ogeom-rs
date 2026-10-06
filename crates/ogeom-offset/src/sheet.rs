//! Sheets: ruled surfaces, lofts and sweeps that bound no volume.
//!
//! Sections, profiles and rails are edges or wires, open or closed, planar
//! or not. Each section edge is restated as its exact B-spline (a conic as
//! its rational form), the edges the sections pair up are made compatible
//! (one degree, one knot vector), and the skin interpolates their control
//! points in homogeneous coordinates. Every section therefore lies on the
//! skin exactly, not to a fit's tolerance. A ruled skin is degree one
//! across, so every ruling is a straight line. A loft with guide curves is
//! a Gordon surface over the sections and the guides (see `guided`).
//!
//! A sweep places copies of its profile along the path (rigid copies under
//! a frame law, similar copies between two rails) and skins them the same
//! way. The copies lie on the skin exactly; between them the skin departs
//! from the swept profile, and that departure is measured against copies
//! placed halfway, by projection, and held to ten confusions, the copies
//! doubled until it is.
//!
//! A sheet has one face per section edge (per span between neighbouring
//! sections for a ruled loft, per spine edge for a sweep), the faces meeting on shared edges. Sections
//! with different edge counts are matched by arc length: each section's
//! edges are split where the others' breaks fall, as fractions of its
//! length, the pieces keeping their exact form. A section edge that bounds
//! the sheet, belongs to the model and was not split is the caller's own
//! edge, so the sheet sews to what it was built from. Each face's normal
//! is its chart's: `u` runs along the sections, `v` across them.

use ogeom_algo::{
    Built, History, attach_pcurve, attach_seam, edge_vertices, make_edge_between, make_face_on,
    make_face_with_pcurves, make_shell, make_vertex, make_wire,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail, ogeom_err};
use ogeom_geom::Curve2d as _;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{
    BSpline2d, BSplineCurve, BSplineSurface, Curve, Line2d, LineCurve, PlanarCurve, PlaneSurface,
    SurfaceGeometry,
};
use ogeom_math::Blend as _;
use ogeom_math::{
    Axis2, ControlGrid, Direction, Direction2, Frame, KnotVector, Plane, Point, Point2, Transform,
    TransformKind, Vector, Weighted,
};
use ogeom_topo::{Location, Model, Orientation, Shape, ShapeType};

use crate::sweep::{PipeLaw, SpineStation, law_normals, spine_curve_of, station_frame};

mod guided;

/// Knots closer than this, on the unit domain every section is restated
/// over, are one knot.
const KNOT_SAME: f64 = 1e-12;

/// The most sections a sweep's skin is built through.
const MOST_SWEEP_SECTIONS: usize = 1025;

/// The ruled surface between two curves: each point of `a` joined by a
/// straight line to the point of `b` at the same fraction of its parameter.
///
/// `a` and `b` are edges or wires, open or closed, taken in their own
/// traversal sense; wires pair edge for edge, one face per pair, and wires
/// of different edge counts are first split to match by arc length (see
/// [`make_loft_surface`]). Between two straight segments the face is the
/// plane when the four corners share one and the bilinear patch otherwise;
/// between curves it is the exact rational B-spline of degree one across. The long edges are `a`'s and
/// `b`'s own edges where they were not split, and the rulings at the ends
/// are straight segments.
///
/// # Errors
///
/// As [`make_loft_surface`] with two sections and `ruled`.
pub fn make_ruled(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    make_loft_surface(model, &[a.clone(), b.clone()], false, &[], true, tol)
}

/// A sheet lofted through sections, in order.
///
/// The sections are edges or wires, open or closed, planar or not, each
/// taken in its own traversal sense from its own start (aligning those is
/// the caller's authorship). Wires pair edge for edge, and the sheet has
/// one face per edge (one per edge and span between neighbouring sections
/// when `ruled`). Sections with different edge counts are matched by arc
/// length first: every section is split at the fractions of its length
/// where any section has a break (breaks closer than ten confusions along
/// the longest section are one), each piece the exact restriction of the
/// edge it came from, so the sections pair edge for edge. The sheet
/// passes through every section exactly. A smooth loft is the B-spline
/// interpolating the sections across, cubic from four sections up,
/// parameterized by the mean distance between their control points; a
/// ruled loft is degree one between each pair of neighbours.
///
/// With `closed` the sheet runs from the last section back to the first:
/// smoothly, with no seam in its slope, or ruled, with a last span. A
/// closed smooth loft spaces its sections evenly in its parameter.
///
/// The first and last sections (every section, when ruled) bound the
/// sheet, and where they are edges of the model, not split to match, they
/// are those edges.
///
/// With `guides` the sheet also follows each guide, an edge or a wire
/// crossing every section once, in the same order along every section and
/// in the sections' order along itself. The guided sheet is one face (the
/// sections are single edges) built as a Gordon surface: each section is
/// paced so every guide crosses it at one parameter, each guide so it
/// crosses every section at one parameter, and the skin across the
/// sections is corrected along each guide by the guide's departure from
/// it. The sheet passes through every section exactly and through every
/// guide to within the distance by which the guide misses the sections;
/// both are measured on the built skin and held to ten confusions. Through
/// closed sections the seam runs along the first guide.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// with fewer than two sections, or three for a closed loft; if a section
/// is not an edge or a wire, or has a curve with no exact B-spline form (a
/// helix, an offset); if the sections differ in being closed; if they
/// differ in edge count and one has no length or their breaks fall too
/// close together to match; if two neighbouring sections coincide; or if
/// neighbouring edges of the sections carry weights at their shared corner
/// that would part their faces. With guides, also if the loft is ruled or closed, a section
/// has more than one edge, a guide misses a section by more than ten
/// confusions, crosses the sections out of their order, or crosses some
/// sections at an end and others inside, or two guides cross the sections
/// in different orders.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the guided
/// skin strays more than ten confusions from a section or a guide.
pub fn make_loft_surface(
    model: &mut Model,
    sections: &[Shape],
    closed: bool,
    guides: &[Shape],
    ruled: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !guides.is_empty() {
        return guided::guided_loft(model, sections, closed, guides, ruled, tol);
    }
    let least = if closed { 3 } else { 2 };
    if sections.len() < least {
        ogeom_bail!(
            Construction,
            "a {}loft surface needs at least {least} sections, given {}",
            if closed { "closed " } else { "" },
            sections.len()
        );
    }
    let mut read: Vec<Section> = Vec::with_capacity(sections.len());
    for shape in sections {
        read.push(read_section(model, shape, "loft section", tol)?);
    }
    let closed_u = read[0].closed;
    for (k, s) in read.iter().enumerate() {
        if s.closed != closed_u {
            ogeom_bail!(
                Construction,
                "loft sections are all closed or all open; section {k} is {} and section 0 \
                 is not",
                if s.closed { "closed" } else { "open" }
            );
        }
    }
    if read.iter().any(|s| s.edges.len() != read[0].edges.len()) {
        matched_by_length(&mut read, tol)?;
    }
    let count = read[0].edges.len();
    // One degree and one knot vector per edge the sections pair up.
    for e in 0..count {
        let curves: Vec<BSplineCurve> = read.iter().map(|s| s.edges[e].curve.clone()).collect();
        let matched = made_compatible(&curves, tol)?;
        for (s, c) in read.iter_mut().zip(matched) {
            s.edges[e].curve = c;
        }
    }
    let n = read.len();
    let mut chords = Vec::with_capacity(n);
    for k in 0..n - 1 {
        chords.push(section_chord(&read[k], &read[k + 1], k, tol)?);
    }
    if closed {
        chords.push(section_chord(&read[n - 1], &read[0], n - 1, tol)?);
    }

    let (spans, surfaces, seam_v) = if ruled {
        let mut spans: Vec<(usize, usize)> = (0..n - 1).map(|k| (k, k + 1)).collect();
        if closed {
            spans.push((n - 1, 0));
        }
        let mut surfaces = Vec::with_capacity(count);
        for e in 0..count {
            let mut per_span = Vec::with_capacity(spans.len());
            for &(lo, hi) in &spans {
                let rows = vec![
                    read[lo].edges[e].curve.control_points().to_vec(),
                    read[hi].edges[e].curve.control_points().to_vec(),
                ];
                let v_knots = KnotVector::clamped_uniform(1, 2)?;
                per_span.push(surface_of(
                    read[lo].edges[e].curve.knots(),
                    v_knots,
                    &rows,
                    tol,
                )?);
            }
            surfaces.push(per_span);
        }
        (spans, surfaces, false)
    } else if closed {
        let mut surfaces = Vec::with_capacity(count);
        for e in 0..count {
            let rows: Vec<Vec<Weighted<Point>>> = read
                .iter()
                .map(|s| s.edges[e].curve.control_points().to_vec())
                .collect();
            let (v_knots, net) = interpolated_closed(&rows, tol)?;
            surfaces.push(vec![surface_of(
                read[0].edges[e].curve.knots(),
                v_knots,
                &net,
                tol,
            )?]);
        }
        // The seam is the first section's curve, built fresh: one edge
        // bounding the face on both sides of its chart.
        for edge in &mut read[0].edges {
            edge.edge = None;
        }
        (vec![(0, 0)], surfaces, true)
    } else {
        let params = unit_params(&chords);
        let mut surfaces = Vec::with_capacity(count);
        for e in 0..count {
            let rows: Vec<Vec<Weighted<Point>>> = read
                .iter()
                .map(|s| s.edges[e].curve.control_points().to_vec())
                .collect();
            let (v_knots, net) = interpolated(&rows, &params, tol)?;
            surfaces.push(vec![surface_of(
                read[0].edges[e].curve.knots(),
                v_knots,
                &net,
                tol,
            )?]);
        }
        (vec![(0, n - 1)], surfaces, false)
    };
    let shape = sheet(model, &read, &surfaces, &spans, seam_v, ruled, tol)?;
    let mut history = History::new();
    for section in sections {
        history.generate(section, shape.clone());
    }
    Ok(Built { shape, history })
}

/// Sweep a profile along a spine into a sheet.
///
/// The profile is an edge or a wire, open or closed, planar or not, swept
/// as it stands: it is the sheet's first section, and `law` turns it about
/// the spine the way it turns a pipe's section ([`PipeLaw::Fixed`] carries
/// it by translation). The sheet has one face per profile edge and spine
/// edge, each spine edge's faces skinned through their own copies and
/// meeting the next edge's on the copy at the edges' shared vertex, and no
/// caps; its far end is the profile where the law carries it to the
/// spine's end. The skin passes through rigid copies of the profile placed
/// along the spine, exactly, and keeps within ten confusions of the swept
/// profile between them, measured halfway between each pair of copies.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the profile is not an edge or a wire or has a curve with no exact
/// B-spline form; if the spine is closed (a sweep surface along a closed
/// spine is not built), turns a sharp corner, or stands still; and as the
/// law refuses its spine (a Frenet frame along a spine that never bends,
/// a guide the stations' planes do not cross).
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the skin
/// cannot reach the target with the most sections it is built through.
pub fn make_sweep_surface(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let section = read_section(model, profile, "sweep profile", tol)?;
    let target = tol.confusion() * 10.0;
    let motions = |model: &Model, density: usize| -> OgeomResult<Placements> {
        let stations = sweep_stations(model, spine, density, tol)?;
        let breaks = stations
            .windows(2)
            .enumerate()
            .filter(|(_, pair)| pair[0].edge != pair[1].edge)
            .map(|(k, _)| k)
            .collect();
        if matches!(law, PipeLaw::Fixed) {
            let start = stations[0].at;
            return Ok(Placements {
                motions: stations
                    .iter()
                    .map(|s| Transform::translation(s.at - start))
                    .collect(),
                breaks,
            });
        }
        let normals = law_normals(model, &stations, law, target, tol)?;
        let start = station_frame(&stations[0], normals[0], tol)?;
        let mut out = Vec::with_capacity(stations.len());
        for (s, n) in stations.iter().zip(&normals) {
            let frame = station_frame(s, *n, tol)?;
            out.push(Transform::from_frame(&frame) * Transform::to_frame(&start));
        }
        Ok(Placements {
            motions: out,
            breaks,
        })
    };
    let shape = swept_sheet(model, &section, motions, target, tol)?;
    let mut history = History::new();
    history.generate(profile, shape.clone());
    history.generate(spine, shape.clone());
    if let PipeLaw::Auxiliary { guide } = law {
        history.generate(guide, shape.clone());
    }
    Ok(Built { shape, history })
}

/// Sweep a profile between two rails into a sheet.
///
/// The profile is an edge or a wire, open or closed, planar or not, whose
/// start sits on the start of one rail and whose end on the start of the
/// other. At each fraction of the rails' lengths the profile is placed by
/// the similarity carrying its ends onto the rails' points there: scaled by
/// the distance between those points, turned so the chord between its ends
/// follows the chord between them, and turned about that chord with the
/// rails' mean direction. The profile is the sheet's first section; the
/// sheet has one face per profile edge, passes through the placed copies
/// exactly and keeps within ten confusions of the swept profile and the
/// rails between them, measured halfway between each pair of copies. The
/// long edges are the skin's own, built on the rails' line.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the profile or a rail is not an edge or a wire, a rail is closed, the
/// profile's ends do not sit on the rails' starts, the rails meet, or the
/// rails' mean direction runs along the chord between them somewhere.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the skin
/// cannot reach the target with the most sections it is built through.
pub fn make_sweep_two_rails(
    model: &mut Model,
    profile: &Shape,
    rail_a: &Shape,
    rail_b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let section = read_section(model, profile, "sweep profile", tol)?;
    if section.closed {
        ogeom_bail!(
            Construction,
            "a two-rail profile runs from one rail to the other; a closed profile has one end"
        );
    }
    let start = section.edges[0].curve.point_at(0.0, tol)?;
    let end = section.edges[section.edges.len() - 1]
        .curve
        .point_at(1.0, tol)?;
    let mut first = Rail::read(model, rail_a, tol)?;
    let mut second = Rail::read(model, rail_b, tol)?;
    let (a0, b0) = (first.at(0.0, tol)?.0, second.at(0.0, tol)?.0);
    let near = tol.confusion() * 10.0;
    if start.distance(b0) <= near && end.distance(a0) <= near {
        core::mem::swap(&mut first, &mut second);
    }
    let (a0, b0) = (first.at(0.0, tol)?.0, second.at(0.0, tol)?.0);
    if start.distance(a0) > near || end.distance(b0) > near {
        ogeom_bail!(
            Construction,
            "the profile's ends must sit on the rails' starts; they are {:.3e} and {:.3e} away",
            start.distance(a0),
            end.distance(b0)
        );
    }
    let frame_at = |f: f64| -> OgeomResult<(Frame, f64)> {
        let (a, ta) = first.at(f, tol)?;
        let (b, tb) = second.at(f, tol)?;
        let chord = b - a;
        let width = chord.magnitude();
        if width <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "the rails meet at {a:?}; a profile between them has no width there"
            );
        }
        let x = chord / width;
        let mean = ta + tb;
        let along = mean - x * mean.dot(x);
        if along.magnitude() <= 1e-6 * mean.magnitude().max(tol.confusion()) {
            ogeom_bail!(
                Construction,
                "the rails' mean direction runs along the chord between them at {a:?}"
            );
        }
        Ok((
            Frame::new(a, Direction::new(along, tol)?, Direction::new(x, tol)?, tol)?,
            width,
        ))
    };
    let (start_frame, start_width) = frame_at(0.0)?;
    let target = tol.confusion() * 10.0;
    let motions = |_: &Model, density: usize| -> OgeomResult<Placements> {
        let count = 32 * density;
        let mut out = Vec::with_capacity(count + 1);
        for i in 0..=count {
            #[allow(clippy::cast_precision_loss)]
            let f = i as f64 / count as f64;
            let (frame, width) = frame_at(f)?;
            let scale = Transform::scaling(Point::ORIGIN, width / start_width, tol)?;
            out.push(Transform::from_frame(&frame) * scale * Transform::to_frame(&start_frame));
        }
        Ok(Placements {
            motions: out,
            breaks: Vec::new(),
        })
    };
    let shape = swept_sheet(model, &section, motions, target, tol)?;
    let mut history = History::new();
    for input in [profile, rail_a, rail_b] {
        history.generate(input, shape.clone());
    }
    Ok(Built { shape, history })
}

/// One edge of a section: the caller's edge where the sheet may bound
/// itself with it, and its exact B-spline.
#[derive(Clone)]
struct SectionEdge {
    /// The model edge, oriented the way the section runs.
    edge: Option<Shape>,
    /// The edge's curve in the section's sense, over `[0, 1]`, its first
    /// weight one unless it is a piece cut from an edge.
    curve: BSplineCurve,
    /// Whether `curve`'s parameter is the edge's own, mapped affinely.
    paced: bool,
}

/// A section as the skin reads it: its edges in traversal order.
#[derive(Clone)]
struct Section {
    edges: Vec<SectionEdge>,
    closed: bool,
}

/// Read an edge or a wire as a section. Its edges are kept for the sheet to
/// bound itself with when they stand unplaced and chain vertex to vertex.
fn read_section(model: &Model, shape: &Shape, what: &str, tol: Tolerances) -> OgeomResult<Section> {
    let edges = match model.kind_of(shape)? {
        ShapeType::Edge => vec![shape.clone()],
        ShapeType::Wire => model.ordered_children_of(shape)?,
        other => ogeom_bail!(
            Construction,
            "a {what} is an edge or a wire, not a {other:?}"
        ),
    };
    if edges.is_empty() {
        ogeom_bail!(Construction, "the {what} has no edges");
    }
    let mut adoptable = true;
    let mut out = Vec::with_capacity(edges.len());
    for edge in &edges {
        let (curve, range) = spine_curve_of(model, edge)?;
        let placement = edge.transform(model.datums())?;
        if placement.kind() != TransformKind::Identity {
            adoptable = false;
        }
        let exact = curve.to_bspline_over(range, tol)?;
        let exact = moved(&exact, &placement, tol)?;
        let exact = if edge.orientation() == Orientation::Reversed {
            let (knots, control) =
                ogeom_math::bspline::reverse(exact.knots(), exact.control_points());
            BSplineCurve::rational(knots, control)?
        } else {
            exact
        };
        out.push(SectionEdge {
            edge: Some(edge.clone()),
            curve: standard(&exact)?,
            paced: true,
        });
    }
    let head = out[0].curve.point_at(0.0, tol)?;
    let tail = out[out.len() - 1].curve.point_at(1.0, tol)?;
    let closed = head.distance(tail) <= tol.confusion();
    if adoptable {
        let mut ends = Vec::with_capacity(edges.len());
        for edge in &edges {
            match edge_vertices(model, edge)? {
                Some(pair) => ends.push(pair),
                None => adoptable = false,
            }
        }
        if adoptable {
            adoptable = ends.windows(2).all(|w| w[0].1.is_partner(&w[1].0))
                && (!closed || ends[ends.len() - 1].1.is_partner(&ends[0].0));
        }
    }
    if !adoptable {
        for e in &mut out {
            e.edge = None;
        }
    }
    Ok(Section { edges: out, closed })
}

/// A B-spline moved by a similarity, exactly: its control points carried,
/// its weights kept.
fn moved(curve: &BSplineCurve, motion: &Transform, tol: Tolerances) -> OgeomResult<BSplineCurve> {
    if motion.kind() == TransformKind::Identity {
        return Ok(curve.clone());
    }
    let control = curve
        .control_points()
        .iter()
        .map(|w| Weighted::new(motion.apply(w.point()), w.weight, tol))
        .collect::<OgeomResult<Vec<_>>>()?;
    BSplineCurve::rational(curve.knots().clone(), control)
}

/// The same curve over `[0, 1]` with its first weight one: scaling every
/// weight alike leaves a rational curve where it is.
fn standard(curve: &BSplineCurve) -> OgeomResult<BSplineCurve> {
    let knots = curve.knots().reparameterized(0.0, 1.0)?;
    let first = curve.control_points()[0].weight;
    let control = curve
        .control_points()
        .iter()
        .map(|w| w.scale(1.0 / first))
        .collect();
    BSplineCurve::rational(knots, control)
}

/// Sections split to pair edge for edge by arc length: each section is cut
/// at the fractions of its length where any section has a break, so every
/// section ends with the same breaks, at the same fractions. Breaks closer
/// than ten confusions along the longest section are one, and a section
/// whose own break is among them keeps it. A piece is the exact
/// restriction of the edge it is cut from. A section that is cut bounds the
/// sheet with fresh edges rather than its own.
fn matched_by_length(sections: &mut [Section], tol: Tolerances) -> OgeomResult<()> {
    let mut lengths: Vec<Vec<f64>> = Vec::with_capacity(sections.len());
    for (k, s) in sections.iter().enumerate() {
        let each = s
            .edges
            .iter()
            .map(|e| ogeom_algo::curve_length(&Curve::BSpline(e.curve.clone()), (0.0, 1.0), tol))
            .collect::<OgeomResult<Vec<f64>>>()?;
        if each.iter().sum::<f64>() <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "section {k} has no length to match the other sections' edges along"
            );
        }
        lengths.push(each);
    }
    let longest = lengths
        .iter()
        .map(|l| l.iter().sum::<f64>())
        .fold(0.0_f64, f64::max);
    let same = tol.confusion() * 10.0 / longest;
    // Each section's own breaks, as fractions of its length.
    let own: Vec<Vec<f64>> = lengths
        .iter()
        .map(|l| {
            let total: f64 = l.iter().sum();
            let mut run = 0.0;
            l[..l.len() - 1]
                .iter()
                .map(|x| {
                    run += x;
                    run / total
                })
                .collect()
        })
        .collect();
    let mut union: Vec<f64> = own.iter().flatten().copied().collect();
    union.sort_by(f64::total_cmp);
    let mut breaks: Vec<f64> = Vec::with_capacity(union.len());
    for f in union {
        if f <= same || f >= 1.0 - same {
            continue;
        }
        if breaks.last().is_none_or(|b| f - b > same) {
            breaks.push(f);
        }
    }
    for ((section, l), mine) in sections.iter_mut().zip(&lengths).zip(&own) {
        let total: f64 = l.iter().sum();
        let mut pieces: Vec<SectionEdge> = Vec::with_capacity(breaks.len() + 1);
        let mut cut = false;
        let mut start = 0.0;
        for (edge, length) in section.edges.iter().zip(l) {
            let (lo, hi) = (start / total, (start + length) / total);
            // The fractions inside this edge none of its own breaks stands
            // for, as lengths from its start.
            let inside: Vec<f64> = breaks
                .iter()
                .filter(|f| **f > lo + same && **f < hi - same)
                .filter(|f| mine.iter().all(|m| (*m - **f).abs() > same))
                .map(|f| f * total - start)
                .collect();
            start += length;
            if inside.is_empty() {
                pieces.push(edge.clone());
                continue;
            }
            cut = true;
            let whole = Curve::BSpline(edge.curve.clone());
            let mut at = Vec::with_capacity(inside.len());
            for along in inside {
                at.push(ogeom_algo::parameter_at_length(
                    &whole,
                    (0.0, 1.0),
                    along,
                    tol,
                )?);
            }
            let mut rest = (
                edge.curve.knots().clone(),
                edge.curve.control_points().to_vec(),
            );
            // The rest keeps the edge's parameter, so each cut is made at
            // the parameter found on the whole edge. The pieces keep their
            // weights as cut, so where two meet the weights agree and the
            // corner keeps the ratio a section with a vertex there has.
            let piece = |(knots, control): (KnotVector, Vec<Weighted<Point>>)| {
                Ok::<_, ogeom_core::OgeomError>(SectionEdge {
                    edge: None,
                    curve: BSplineCurve::rational(knots.reparameterized(0.0, 1.0)?, control)?,
                    paced: false,
                })
            };
            for t in at {
                let (before, after) = ogeom_math::bspline::split(&rest.0, &rest.1, t, tol)?;
                pieces.push(piece(before)?);
                rest = after;
            }
            pieces.push(piece(rest)?);
        }
        if cut {
            for piece in &mut pieces {
                piece.edge = None;
            }
        }
        section.edges = pieces;
    }
    let count = sections[0].edges.len();
    if let Some(k) = sections.iter().position(|s| s.edges.len() != count) {
        ogeom_bail!(
            Construction,
            "the sections' breaks fall too close together to match by length: section 0 \
             splits into {count} edges and section {k} into {}",
            sections[k].edges.len()
        );
    }
    Ok(())
}

/// Curves raised to one degree and refined to one knot vector, each
/// unchanged as a curve.
fn made_compatible(curves: &[BSplineCurve], tol: Tolerances) -> OgeomResult<Vec<BSplineCurve>> {
    let degree = curves.iter().map(BSplineCurve::degree).max().unwrap_or(1);
    let mut raised = Vec::with_capacity(curves.len());
    for c in curves {
        let mut c = c.clone();
        while c.degree() < degree {
            c = c.elevated(tol)?;
        }
        raised.push(c);
    }
    let interior = |c: &BSplineCurve| -> Vec<(f64, usize)> {
        c.knots()
            .distinct()
            .into_iter()
            .filter(|(v, _)| *v > KNOT_SAME && *v < 1.0 - KNOT_SAME)
            .collect()
    };
    let mut union: Vec<(f64, usize)> = Vec::new();
    for c in &raised {
        for (value, mult) in interior(c) {
            match union
                .iter_mut()
                .find(|(v, _)| (*v - value).abs() <= KNOT_SAME)
            {
                Some(entry) => entry.1 = entry.1.max(mult),
                None => union.push((value, mult)),
            }
        }
    }
    for c in &mut raised {
        for &(value, mult) in &union {
            let have = interior(c)
                .iter()
                .find(|(v, _)| (*v - value).abs() <= KNOT_SAME)
                .map_or(0, |entry| entry.1);
            if have < mult {
                *c = c.with_knot_inserted(value, mult - have, tol)?;
            }
        }
    }
    let knots = raised[0].knots().clone();
    for c in &raised {
        let same = c.knots().knots().len() == knots.knots().len()
            && c.knots()
                .knots()
                .iter()
                .zip(knots.knots())
                .all(|(a, b)| (a - b).abs() <= KNOT_SAME * 10.0);
        if !same {
            ogeom_bail!(Construction, "the sections' knots could not be matched");
        }
    }
    raised
        .iter()
        .map(|c| BSplineCurve::rational(knots.clone(), c.control_points().to_vec()))
        .collect()
}

/// The mean distance between two compatible sections' control points.
fn section_chord(a: &Section, b: &Section, index: usize, tol: Tolerances) -> OgeomResult<f64> {
    let mut sum = 0.0;
    let mut count = 0.0;
    for (ea, eb) in a.edges.iter().zip(&b.edges) {
        for (p, q) in ea
            .curve
            .control_points()
            .iter()
            .zip(eb.curve.control_points())
        {
            sum += p.point().distance(q.point());
            count += 1.0;
        }
    }
    let chord = sum / count;
    if chord <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "sections {index} and {} coincide; a skin between them has no extent",
            index + 1
        );
    }
    Ok(chord)
}

/// Cumulative chords over `[0, 1]`.
fn unit_params(chords: &[f64]) -> Vec<f64> {
    let total: f64 = chords.iter().sum();
    let mut params = Vec::with_capacity(chords.len() + 1);
    let mut run = 0.0;
    params.push(0.0);
    for c in chords {
        run += c;
        params.push(run / total);
    }
    let last = params.len() - 1;
    params[last] = 1.0;
    params
}

/// The rows of control points (one per section, homogeneous) interpolated
/// across at `params`: cubic from four rows up, over the averaged knots.
/// Returns the knots across and the interpolating rows.
fn interpolated(
    rows: &[Vec<Weighted<Point>>],
    params: &[f64],
    tol: Tolerances,
) -> OgeomResult<(KnotVector, Vec<Vec<Weighted<Point>>>)> {
    let n = rows.len();
    let degree = (n - 1).min(3);
    let knots = KnotVector::averaged(degree, params)?;
    let mut matrix = vec![vec![0.0; n]; n];
    for (k, v) in params.iter().enumerate() {
        let span = knots.span(*v, tol)?;
        for (b, j) in knots.basis(span, *v).iter().zip(span - degree..=span) {
            matrix[k][j] = *b;
        }
    }
    let inverse = inverted(matrix)?;
    Ok((knots, combined(&inverse, rows)))
}

/// The rows interpolated round a loop: a periodic B-spline through every
/// row and back to the first, the rows evenly spaced in its parameter,
/// cubic from four rows up and quadratic for three. Cut open at the first
/// row and clamped there, so its first and last rows are that row.
fn interpolated_closed(
    rows: &[Vec<Weighted<Point>>],
    tol: Tolerances,
) -> OgeomResult<(KnotVector, Vec<Vec<Weighted<Point>>>)> {
    let n = rows.len();
    let degree = (n - 1).min(3);
    // An even degree is interpolated mid-span, where its basis is
    // dominated by its own control; an odd one at the knots.
    let shift = if degree.is_multiple_of(2) { 0.5 } else { 0.0 };
    // The ring's controls Q[0..n) wrapped one before and degree + 1 after,
    // over uniform knots: the periodic curve over a domain wide enough to
    // cut a whole turn from inside it.
    let count = n + degree + 2;
    #[allow(clippy::cast_precision_loss)]
    let knots = KnotVector::new((0..=count + degree).map(|i| i as f64).collect(), degree)?;
    let ring = |j: usize| (j + n - 1) % n;
    #[allow(clippy::cast_precision_loss)]
    let at = |k: usize| (degree + 1 + k) as f64 + shift;
    let mut matrix = vec![vec![0.0; n]; n];
    for (k, row) in matrix.iter_mut().enumerate() {
        let v = at(k);
        let span = knots.span(v, tol)?;
        for (b, j) in knots.basis(span, v).iter().zip(span - degree..=span) {
            row[ring(j)] += *b;
        }
    }
    let inverse = inverted(matrix)?;
    let solved = combined(&inverse, rows);
    let wrapped: Vec<Vec<Weighted<Point>>> = (0..count).map(|j| solved[ring(j)].clone()).collect();
    #[allow(clippy::cast_precision_loss)]
    let (from, to) = (at(0), at(n));
    let width = rows[0].len();
    let mut columns: Vec<Vec<Weighted<Point>>> = Vec::with_capacity(width);
    let mut cut_knots: Option<KnotVector> = None;
    for i in 0..width {
        let column: Vec<Weighted<Point>> = wrapped.iter().map(|r| r[i]).collect();
        let (_, (k1, c1)) = ogeom_math::bspline::split(&knots, &column, from, tol)?;
        let ((k2, c2), _) = ogeom_math::bspline::split(&k1, &c1, to, tol)?;
        cut_knots = Some(k2);
        columns.push(c2);
    }
    let Some(cut_knots) = cut_knots else {
        ogeom_bail!(Construction, "a section has no control points");
    };
    let l = columns[0].len();
    let mut net: Vec<Vec<Weighted<Point>>> = (0..l)
        .map(|j| columns.iter().map(|c| c[j]).collect())
        .collect();
    // The two ends are the first row, to rounding; one row, exactly.
    net[l - 1] = net[0].clone();
    Ok((cut_knots.reparameterized(0.0, 1.0)?, net))
}

/// `inverse` applied to the rows, row by row in homogeneous coordinates.
fn combined(inverse: &[Vec<f64>], rows: &[Vec<Weighted<Point>>]) -> Vec<Vec<Weighted<Point>>> {
    let width = rows[0].len();
    inverse
        .iter()
        .map(|line| {
            (0..width)
                .map(|i| {
                    line.iter()
                        .zip(rows)
                        .fold(Weighted::<Point>::zero(), |acc, (a, row)| {
                            acc.add(row[i].scale(*a))
                        })
                })
                .collect()
        })
        .collect()
}

/// The inverse of a small dense matrix, by Gauss-Jordan elimination with
/// partial pivoting.
fn inverted(mut a: Vec<Vec<f64>>) -> OgeomResult<Vec<Vec<f64>>> {
    let n = a.len();
    let mut inv: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))
            .unwrap_or(col);
        if a[pivot][col].abs() <= 1e-12 {
            ogeom_bail!(Numeric, "the skin's interpolation system is singular");
        }
        a.swap(col, pivot);
        inv.swap(col, pivot);
        let d = a[col][col];
        for j in 0..n {
            a[col][j] /= d;
            inv[col][j] /= d;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = a[r][col];
            if f == 0.0 {
                continue;
            }
            for j in 0..n {
                a[r][j] -= f * a[col][j];
                inv[r][j] -= f * inv[col][j];
            }
        }
    }
    Ok(inv)
}

/// The patch whose `v`-rows are `rows`, over `u_knots` along them.
fn surface_of(
    u_knots: &KnotVector,
    v_knots: KnotVector,
    rows: &[Vec<Weighted<Point>>],
    tol: Tolerances,
) -> OgeomResult<BSplineSurface> {
    let (k, l) = (rows[0].len(), rows.len());
    let mut points = Vec::with_capacity(k * l);
    for i in 0..k {
        for row in rows {
            let w = row[i];
            if !w.weight.is_finite() || w.weight <= tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the skin's weights fall to {} between the sections; sections this unlike \
                     cannot be skinned exactly",
                    w.weight
                );
            }
            points.push(w);
        }
    }
    BSplineSurface::rational(u_knots.clone(), v_knots, ControlGrid::new(points, k, l)?)
}

/// One boundary side of a sheet face: the edge as the face runs it (bottom
/// along `u`, rails along `v`), its image in the chart at its own
/// parameter, and that parameter's range.
struct Side {
    edge: Shape,
    image: PlanarCurve,
    range: (f64, f64),
}

/// The faces of a sheet over its surfaces, `surfaces[e][s]` the skin of
/// section edge `e` over span `s`, bounded by section edges and rails, the
/// faces chained into a shell where there are several.
///
/// `seam_v` marks a skin closed across, whose single span's first section
/// bounds the face on both sides of its chart.
#[allow(clippy::too_many_lines, reason = "one assembly, spelled out")]
fn sheet(
    model: &mut Model,
    sections: &[Section],
    surfaces: &[Vec<BSplineSurface>],
    spans: &[(usize, usize)],
    seam_v: bool,
    ruled: bool,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let count = sections[0].edges.len();
    let closed_u = sections[0].closed;
    let joints = if closed_u { count } else { count + 1 };
    let end_joint = |e: usize| if closed_u { (e + 1) % count } else { e + 1 };

    // Neighbouring edges' skins share the rail at their corner only where
    // the corner's weights keep one ratio through every section.
    for j in 0..joints {
        let (before, after) = if closed_u {
            ((j + count - 1) % count, j)
        } else if j == 0 || j == count {
            continue;
        } else {
            (j - 1, j)
        };
        let ratio = |s: &Section| {
            let end = s.edges[before].curve.control_points();
            let start = s.edges[after].curve.control_points();
            end[end.len() - 1].weight / start[0].weight
        };
        let first = ratio(&sections[0]);
        if sections
            .iter()
            .any(|s| (ratio(s) - first).abs() > 1e-9 * first.abs())
        {
            ogeom_bail!(
                Construction,
                "neighbouring section edges carry weights at their shared corner {j} that \
                 change from section to section; their skins would part"
            );
        }
    }

    // The sections that bound faces: their joint vertices and edges.
    let mut bounding: Vec<usize> = spans.iter().flat_map(|&(a, b)| [a, b]).collect();
    bounding.sort_unstable();
    bounding.dedup();
    let mut vertices: Vec<Option<Vec<Shape>>> = vec![None; sections.len()];
    let mut borders: Vec<Option<Vec<(Shape, bool)>>> = vec![None; sections.len()];
    for &k in &bounding {
        let section = &sections[k];
        let adopted = section.edges.iter().all(|e| e.edge.is_some());
        let mut joint_vertices = Vec::with_capacity(joints);
        let mut edges = Vec::with_capacity(count);
        if adopted {
            for e in &section.edges {
                let Some(edge) = &e.edge else {
                    ogeom_bail!(Construction, "an adopted section lost an edge");
                };
                let Some((start, _)) = edge_vertices(model, edge)? else {
                    ogeom_bail!(Construction, "a section edge has no vertices");
                };
                joint_vertices.push(start);
                edges.push((edge.clone(), true));
            }
            if !closed_u {
                let Some(last) = &section.edges[count - 1].edge else {
                    ogeom_bail!(Construction, "an adopted section lost an edge");
                };
                let Some((_, end)) = edge_vertices(model, last)? else {
                    ogeom_bail!(Construction, "a section edge has no vertices");
                };
                joint_vertices.push(end);
            }
        } else {
            for e in &section.edges {
                let at = e.curve.point_at(0.0, tol)?;
                joint_vertices.push(make_vertex(model, at).shape);
            }
            if !closed_u {
                let at = section.edges[count - 1].curve.point_at(1.0, tol)?;
                joint_vertices.push(make_vertex(model, at).shape);
            }
            for (e, piece) in section.edges.iter().enumerate() {
                let edge = make_edge_between(
                    model,
                    Curve::BSpline(piece.curve.clone()),
                    (0.0, 1.0),
                    &joint_vertices[e],
                    &joint_vertices[end_joint(e)],
                    tol,
                )?
                .shape;
                edges.push((edge, false));
            }
        }
        vertices[k] = Some(joint_vertices);
        borders[k] = Some(edges);
    }
    let vertex = |k: usize, j: usize| -> OgeomResult<Shape> {
        vertices[k]
            .as_ref()
            .map(|v| v[j].clone())
            .ok_or_else(|| ogeom_err!(Construction, "section {k} bounds no face"))
    };

    // Rails: one per span and joint, shared by the faces either side.
    let mut rails: Vec<Vec<(Shape, bool)>> = Vec::with_capacity(spans.len());
    for (s, &(lo, hi)) in spans.iter().enumerate() {
        let mut per_joint = Vec::with_capacity(joints);
        for j in 0..joints {
            let (e, u) = if j < count {
                (j, 0.0)
            } else {
                (count - 1, 1.0)
            };
            let surface = &surfaces[e][s];
            let control = sections[lo].edges[e].curve.control_points();
            let index = if u == 0.0 { 0 } else { control.len() - 1 };
            let w_lo = control[index].weight;
            let w_hi = sections[hi].edges[e].curve.control_points()[index].weight;
            let (from, to) = (vertex(lo, j)?, vertex(hi, j)?);
            // A ruled rail between equal weights runs at an even pace: the
            // straight segment, its chart image the matching column.
            if ruled && (w_lo - w_hi).abs() <= 1e-12 * w_lo {
                let a = surface.point_at(u, 0.0, tol)?;
                let b = surface.point_at(u, 1.0, tol)?;
                let line: Curve = LineCurve::segment(a, b, tol)?.into();
                let range = line.domain();
                let edge = make_edge_between(model, line, range, &from, &to, tol)?.shape;
                per_joint.push((edge, true));
            } else {
                let curve = Curve::BSpline(surface.iso_u_curve(u, tol)?);
                let edge = make_edge_between(model, curve, (0.0, 1.0), &from, &to, tol)?.shape;
                per_joint.push((edge, false));
            }
        }
        rails.push(per_joint);
    }

    let column = |u: f64, straight: bool, range: (f64, f64)| -> OgeomResult<PlanarCurve> {
        if straight {
            let knots = KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?;
            Ok(BSpline2d::new(knots, vec![Point2::new(u, 0.0), Point2::new(u, 1.0)], tol)?.into())
        } else {
            Ok(Line2d::over(Axis2::new(Point2::new(u, 0.0), Direction2::Y), -1.0, 2.0)?.into())
        }
    };
    let rail_range = |model: &Model, edge: &Shape| -> OgeomResult<(f64, f64)> {
        Ok(spine_curve_of(model, edge)?.1)
    };

    let mut faces = Vec::with_capacity(count * spans.len());
    for e in 0..count {
        for (s, &(lo, hi)) in spans.iter().enumerate() {
            let surface = &surfaces[e][s];
            let geometry: SurfaceGeometry = surface.clone().into();
            let border = |k: usize| -> OgeomResult<(Shape, bool)> {
                borders[k]
                    .as_ref()
                    .map(|b| b[e].clone())
                    .ok_or_else(|| ogeom_err!(Construction, "section {k} bounds no face"))
            };
            let (bottom, bottom_adopted) = border(lo)?;
            let (top, top_adopted) = border(hi)?;
            let (rail0, straight0) = rails[s][e].clone();
            let (rail1, straight1) = rails[s][end_joint(e)].clone();

            if ruled
                && let Some(plane) =
                    ruled_plane(&sections[lo].edges[e], &sections[hi].edges[e], tol)?
            {
                let reach = [0.0, 1.0]
                    .iter()
                    .flat_map(|u| [(*u, 0.0), (*u, 1.0)])
                    .map(|(u, v)| {
                        surface
                            .point_at(u, v, tol)
                            .map(|p| p.distance(plane.origin()))
                    })
                    .collect::<OgeomResult<Vec<f64>>>()?
                    .into_iter()
                    .fold(1.0_f64, f64::max)
                    * 2.0;
                let flat: SurfaceGeometry =
                    PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
                let face = make_face_with_pcurves(
                    model,
                    flat,
                    &[vec![bottom, rail1, top.reversed(), rail0.reversed()]],
                    tol,
                )?
                .shape;
                faces.push(face);
                continue;
            }

            let row = |model: &mut Model,
                       edge: &Shape,
                       adopted: bool,
                       v: f64|
             -> OgeomResult<Side> {
                let (image, range) = if adopted {
                    let section = &sections[if v == 0.0 { lo } else { hi }].edges[e];
                    adopted_row(
                        model,
                        edge,
                        &section.curve,
                        section.paced,
                        v,
                        &geometry,
                        tol,
                    )?
                } else {
                    (
                        Line2d::over(Axis2::new(Point2::new(0.0, v), Direction2::X), -1.0, 2.0)?
                            .into(),
                        (0.0, 1.0),
                    )
                };
                Ok(Side {
                    edge: edge.clone(),
                    image,
                    range,
                })
            };
            let bottom_side = row(model, &bottom, bottom_adopted, 0.0)?;
            let top_side = row(model, &top, top_adopted, 1.0)?;
            let range0 = rail_range(model, &rail0)?;
            let range1 = rail_range(model, &rail1)?;
            let rail0_side = Side {
                image: column(0.0, straight0, range0)?,
                edge: rail0,
                range: range0,
            };
            let rail1_side = Side {
                image: column(1.0, straight1, range1)?,
                edge: rail1,
                range: range1,
            };
            debug_assert!(!seam_v || bottom_side.edge.is_partner(&top_side.edge));
            faces.push(sheet_face(
                model,
                geometry,
                bottom_side,
                rail1_side,
                top_side,
                rail0_side,
                tol,
            )?);
        }
    }
    if faces.len() == 1 {
        return Ok(faces.swap_remove(0));
    }
    Ok(make_shell(model, &faces)?.shape)
}

/// The plane a ruled span between two straight segments lies in, framed so
/// its normal is the chart's (along the segments, then across), if their
/// four ends share one.
fn ruled_plane(a: &SectionEdge, b: &SectionEdge, tol: Tolerances) -> OgeomResult<Option<Plane>> {
    let straight = |c: &BSplineCurve| c.degree() == 1 && c.control_points().len() == 2;
    if !straight(&a.curve) || !straight(&b.curve) {
        return Ok(None);
    }
    let (a0, a1) = (a.curve.point_at(0.0, tol)?, a.curve.point_at(1.0, tol)?);
    let (b0, b1) = (b.curve.point_at(0.0, tol)?, b.curve.point_at(1.0, tol)?);
    let normal = [
        (a1 - a0).cross(b0 - a0),
        (a1 - a0).cross(b1 - a1),
        (b1 - b0).cross(b0 - a0),
    ]
    .into_iter()
    .max_by(|x, y| x.magnitude().total_cmp(&y.magnitude()))
    .unwrap_or(Vector::new(0.0, 0.0, 0.0));
    if normal.magnitude() <= tol.confusion() * (a1 - a0).magnitude().max(1.0) {
        return Ok(None);
    }
    let plane = Plane::through(a0, Direction::new(normal, tol)?);
    Ok([a1, b0, b1]
        .iter()
        .all(|p| plane.distance_to(*p) <= tol.confusion())
        .then_some(plane))
}

/// A face on `surface` bounded by its four sides, each side's image
/// attached; a side that bounds the face twice is a seam.
fn sheet_face(
    model: &mut Model,
    surface: SurfaceGeometry,
    bottom: Side,
    rail1: Side,
    top: Side,
    rail0: Side,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let id = model.geometry_mut().add_surface(surface);
    let here = Location::identity;
    if bottom.edge.is_partner(&top.edge) {
        attach_seam(
            model,
            &bottom.edge,
            bottom.image,
            top.image,
            id,
            here(),
            bottom.range,
        )?;
    } else {
        attach_pcurve(model, &bottom.edge, bottom.image, id, here(), bottom.range)?;
        attach_pcurve(model, &top.edge, top.image, id, here(), top.range)?;
    }
    if rail0.edge.is_partner(&rail1.edge) {
        attach_seam(
            model,
            &rail1.edge,
            rail1.image,
            rail0.image,
            id,
            here(),
            rail1.range,
        )?;
    } else {
        attach_pcurve(model, &rail1.edge, rail1.image, id, here(), rail1.range)?;
        attach_pcurve(model, &rail0.edge, rail0.image, id, here(), rail0.range)?;
    }
    let wire = make_wire(
        model,
        &[
            bottom.edge,
            rail1.edge,
            top.edge.reversed(),
            rail0.edge.reversed(),
        ],
        tol,
    )?
    .shape;
    Ok(make_face_on(model, id, &[wire], tol)?.shape)
}

/// The image of a caller's section edge on the row `v` of a skin, at the
/// edge's own parameter. A segment's or a spline's parameter maps onto the
/// row evenly, exactly, where `paced` says the section keeps it; a conic's
/// does not (its rational form runs at another pace), nor does a section
/// paced anew, and its image is fitted through the row at the edge's
/// parameters, the fit measured on the surface against the edge and the
/// edge widened to what it reached.
fn adopted_row(
    model: &mut Model,
    edge: &Shape,
    section: &BSplineCurve,
    paced: bool,
    v: f64,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<(PlanarCurve, (f64, f64))> {
    let (curve, range) = spine_curve_of(model, edge)?;
    let reversed = edge.orientation() == Orientation::Reversed;
    let even = paced
        && match &curve {
            Curve::Line(_) => true,
            Curve::BSpline(b) => b.knots().is_clamped(),
            _ => false,
        };
    if even {
        let (a, b) = if reversed { (1.0, 0.0) } else { (0.0, 1.0) };
        let knots = KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?;
        let image = BSpline2d::new(knots, vec![Point2::new(a, v), Point2::new(b, v)], tol)?;
        return Ok((image.into(), range));
    }
    // The section's rational spans meet with a jump in how fast the
    // angle turns into `u`: each span is fitted on its own and the pieces
    // joined at the spans' ends, where both agree exactly.
    let u_at = |t: f64, guess: f64| -> OgeomResult<f64> {
        foot_on(section, curve.point_at(t, tol)?, guess, tol)
    };
    // `u` along the edge's own parameter: rising, or falling where the
    // edge runs against the section.
    let ends = if reversed { (1.0, 0.0) } else { (0.0, 1.0) };
    let mut breaks: Vec<(f64, f64)> = vec![(range.0, ends.0)];
    let mut inner: Vec<f64> = section
        .knots()
        .distinct()
        .into_iter()
        .map(|(k, _)| k)
        .filter(|k| *k > KNOT_SAME && *k < 1.0 - KNOT_SAME)
        .collect();
    if reversed {
        inner.reverse();
    }
    for knot in inner {
        let (mut lo, mut hi) = (breaks[breaks.len() - 1].0, range.1);
        for _ in 0..80 {
            let mid = f64::midpoint(lo, hi);
            let below = u_at(mid, knot)? < knot;
            if below != reversed {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        breaks.push((f64::midpoint(lo, hi), knot));
    }
    breaks.push((range.1, ends.1));
    // The chart's `u` runs the section once: a step in it moves the point
    // at most the control polygon's length.
    let reach: f64 = section
        .control_points()
        .windows(2)
        .map(|w| w[0].point().distance(w[1].point()))
        .sum();
    let target = tol.confusion() * 0.01 / reach.max(tol.confusion());
    const SAMPLES: u32 = 64;
    let mut joined: Option<(KnotVector, Vec<Weighted<Point2>>)> = None;
    for w in breaks.windows(2) {
        let ((ta, ua), (tb, ub)) = (w[0], w[1]);
        let mut params = Vec::with_capacity(SAMPLES as usize + 1);
        let mut image = Vec::with_capacity(SAMPLES as usize + 1);
        for i in 0..=SAMPLES {
            let f = f64::from(i) / f64::from(SAMPLES);
            let t = ta + (tb - ta) * f;
            let u = if i == 0 {
                ua
            } else if i == SAMPLES {
                ub
            } else {
                u_at(t, ua + (ub - ua) * f)?
            };
            params.push(t);
            image.push(Point2::new(u, v));
        }
        let fitted = ogeom_geom::fit::fit_points_2d_at(&params, &image, 3, target, tol)?;
        let mut control = fitted.curve.control_points().to_vec();
        // Pinned to the span's ends, which are exact.
        let last = control.len() - 1;
        control[0] = Weighted::new(Point2::new(ua, v), 1.0, tol)?;
        control[last] = Weighted::new(Point2::new(ub, v), 1.0, tol)?;
        let piece = (fitted.curve.knots().clone(), control);
        joined = Some(match joined {
            None => piece,
            Some(before) => ogeom_math::bspline::join(&before, &piece)?,
        });
    }
    let Some((knots, control)) = joined else {
        ogeom_bail!(Construction, "a section edge has no extent");
    };
    let pcurve: PlanarCurve = BSpline2d::rational(knots, control)?.into();
    let mut worst: f64 = 0.0;
    for i in 0..=4096 {
        let t = range.0 + (range.1 - range.0) * f64::from(i) / 4096.0;
        let q = pcurve.point_at(t, tol)?;
        let on = surface.point_at(q.x, q.y, tol)?;
        worst = worst.max(on.distance(curve.point_at(t, tol)?));
    }
    if worst > tol.confusion() * 100.0 {
        ogeom_bail!(
            NotDone,
            "a section edge's image on the skin strays {worst:.3e} from the edge"
        );
    }
    if worst > tol.confusion() {
        let widened = ogeom_core::Tolerance::new(worst + tol.confusion())?;
        model.widen(edge, widened)?;
        if let Some((a, b)) = edge_vertices(model, edge)? {
            model.widen(&a, widened)?;
            model.widen(&b, widened)?;
        }
    }
    Ok((pcurve, range))
}

/// The parameter of `section` at `p`, which lies on it: Newton's method
/// from `guess`.
fn foot_on(section: &BSplineCurve, p: Point, guess: f64, tol: Tolerances) -> OgeomResult<f64> {
    let mut u = guess;
    for _ in 0..64 {
        let c = section.point_at(u, tol)?;
        let d = section.d1_at(u, tol)?;
        let speed = d.dot(d);
        if speed <= f64::MIN_POSITIVE {
            break;
        }
        let next = (u + (p - c).dot(d) / speed).clamp(0.0, 1.0);
        let moved = (next - u).abs();
        u = next;
        if moved <= 1e-15 {
            break;
        }
    }
    let off = section.point_at(u, tol)?.distance(p);
    if off > tol.confusion() * 10.0 {
        ogeom_bail!(
            NotDone,
            "a section edge's point {p:?} was not found on its exact form ({off:.3e} away)"
        );
    }
    Ok(u)
}

/// Where a sweep places its profile at one density: a run of motions and
/// the motions at which the path passes from one piece to the next.
struct Placements {
    /// An odd run of motions, the first the identity.
    motions: Vec<Transform>,
    /// Indices into `motions`, each even, where one piece of the path
    /// ends and the next begins.
    breaks: Vec<usize>,
}

/// A sheet swept from copies of `profile` placed by `motions`: for a
/// density, an odd run of placements whose even members the skin passes
/// through and whose odd members (halfway between) it is measured against.
/// The skin is one face per profile edge and piece of the path, each
/// piece's skin interpolating its own copies, so a path that is smooth
/// only to its tangent at a break is not smoothed across it. The density
/// doubles until the skin keeps within `target` of the halfway copies.
fn swept_sheet(
    model: &mut Model,
    profile: &Section,
    motions: impl Fn(&Model, usize) -> OgeomResult<Placements>,
    target: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let count = profile.edges.len();
    let mut density = 1;
    let mut reached = (f64::INFINITY, 0);
    loop {
        let Placements {
            motions: placed,
            breaks,
        } = motions(model, density)?;
        let skin_count = placed.len().div_ceil(2);
        if skin_count > MOST_SWEEP_SECTIONS {
            ogeom_bail!(
                NotDone,
                "the sweep's skin reached {:.3e} against a target of {target:.3e} through {} \
                 sections",
                reached.0,
                reached.1
            );
        }
        let mut sections: Vec<Section> = Vec::with_capacity(skin_count);
        for (k, motion) in placed.iter().step_by(2).enumerate() {
            let mut edges = Vec::with_capacity(count);
            for piece in &profile.edges {
                edges.push(SectionEdge {
                    edge: if k == 0 { piece.edge.clone() } else { None },
                    curve: moved(&piece.curve, motion, tol)?,
                    paced: piece.paced,
                });
            }
            sections.push(Section {
                edges,
                closed: profile.closed,
            });
        }
        let mut bounds = vec![0];
        bounds.extend(
            breaks
                .iter()
                .filter(|b| **b % 2 == 0)
                .map(|b| b / 2)
                .filter(|k| (1..skin_count - 1).contains(k)),
        );
        bounds.push(skin_count - 1);
        let spans: Vec<(usize, usize)> = bounds.windows(2).map(|w| (w[0], w[1])).collect();
        let mut surfaces: Vec<Vec<BSplineSurface>> = vec![Vec::with_capacity(spans.len()); count];
        let mut span_params = Vec::with_capacity(spans.len());
        for &(lo, hi) in &spans {
            let mut chords = Vec::with_capacity(hi - lo);
            for k in lo..hi {
                chords.push(section_chord(&sections[k], &sections[k + 1], k, tol)?);
            }
            let params = unit_params(&chords);
            for (e, per_span) in surfaces.iter_mut().enumerate() {
                let rows: Vec<Vec<Weighted<Point>>> = sections[lo..=hi]
                    .iter()
                    .map(|s| s.edges[e].curve.control_points().to_vec())
                    .collect();
                let (v_knots, net) = interpolated(&rows, &params, tol)?;
                per_span.push(surface_of(
                    profile.edges[e].curve.knots(),
                    v_knots,
                    &net,
                    tol,
                )?);
            }
            span_params.push(params);
        }
        // Halfway copies, projected onto the skin of their span from where
        // they should land.
        let mut worst: f64 = 0.0;
        for (s, &(lo, hi)) in spans.iter().enumerate() {
            let params = &span_params[s];
            for (e, piece) in profile.edges.iter().enumerate() {
                let geometry: SurfaceGeometry = surfaces[e][s].clone().into();
                for k in lo..hi {
                    let motion = &placed[2 * k + 1];
                    let guess_v = f64::midpoint(params[k - lo], params[k - lo + 1]);
                    for i in 0..=16 {
                        let u = f64::from(i) / 16.0;
                        let p = motion.apply(piece.curve.point_at(u, tol)?);
                        let foot =
                            ogeom_algo::project_on_surface_from(&geometry, p, (u, guess_v), tol)?;
                        worst = worst.max(foot.distance);
                    }
                }
            }
        }
        reached = (worst, skin_count);
        if worst <= target {
            return sheet(model, &sections, &surfaces, &spans, false, false, tol);
        }
        density *= 2;
    }
}

/// Stations along an open spine with no sharp corner, evenly by parameter
/// on each edge, an even number per edge by turning (at most a
/// sixty-fourth of a turn apart at density one), `density` times as many.
fn sweep_stations(
    model: &Model,
    spine: &Shape,
    density: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<SpineStation>> {
    let edges: Vec<Shape> = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        ShapeType::Wire => model.ordered_children_of(spine)?,
        other => ogeom_bail!(
            Construction,
            "a sweep runs along an edge or a wire, not a {other:?}"
        ),
    };
    if edges.is_empty() {
        ogeom_bail!(Construction, "the spine has no edge to run along");
    }
    let mut out: Vec<SpineStation> = Vec::new();
    for (ei, edge) in edges.iter().enumerate() {
        let (curve, range) = spine_curve_of(model, edge)?;
        let curve = curve.transformed(&edge.transform(model.datums())?, tol)?;
        let reversed = edge.orientation() == Orientation::Reversed;
        let mut turning = 0.0_f64;
        let mut last: Option<Vector> = None;
        for i in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(i) / 16.0;
            let d = curve.d1_at(t, tol)?;
            if d.magnitude() <= tol.confusion() {
                continue;
            }
            let unit = d / d.magnitude();
            if let Some(prev) = last {
                turning += prev.dot(unit).clamp(-1.0, 1.0).acos();
            }
            last = Some(unit);
        }
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a station count, bounded"
        )]
        let count =
            (turning / (core::f64::consts::TAU / 64.0)).ceil().max(8.0) as usize * 2 * density;
        for i in 0..=count {
            #[allow(clippy::cast_precision_loss)]
            let f = i as f64 / count as f64;
            let t = if reversed {
                range.1 - (range.1 - range.0) * f
            } else {
                range.0 + (range.1 - range.0) * f
            };
            let p = curve.point_at(t, tol)?;
            let d = curve.d1_at(t, tol)?;
            if d.magnitude() <= tol.confusion() {
                ogeom_bail!(Construction, "the spine stands still at {p:?}");
            }
            let tangent = (if reversed { -d } else { d }) / d.magnitude();
            if let Some(prev) = out.last()
                && prev.at.distance(p) <= tol.confusion()
            {
                if i == 0 && prev.tangent.dot(tangent) >= 1.0 - 1e-12 {
                    continue;
                }
                ogeom_bail!(
                    Construction,
                    "the spine turns a sharp corner at {p:?}; a sweep surface follows a \
                     spine whose direction is continuous"
                );
            }
            out.push(SpineStation {
                at: p,
                tangent,
                edge: ei,
                t,
            });
        }
    }
    if out[0].at.distance(out[out.len() - 1].at) <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a sweep surface along a closed spine is not built; sweep along an open one"
        );
    }
    Ok(out)
}

/// A rail as a run of curves measured by length.
struct Rail {
    /// Each edge's curve where it stands, its range, whether it runs
    /// backwards, and its length.
    pieces: Vec<(Curve, (f64, f64), bool, f64)>,
    total: f64,
}

impl Rail {
    fn read(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Self> {
        let edges = match model.kind_of(shape)? {
            ShapeType::Edge => vec![shape.clone()],
            ShapeType::Wire => model.ordered_children_of(shape)?,
            other => ogeom_bail!(Construction, "a rail is an edge or a wire, not a {other:?}"),
        };
        let mut pieces = Vec::with_capacity(edges.len());
        let mut total = 0.0;
        for edge in &edges {
            let (curve, range) = spine_curve_of(model, edge)?;
            let curve = curve.transformed(&edge.transform(model.datums())?, tol)?;
            let length = ogeom_algo::curve_length(&curve, range, tol)?;
            total += length;
            pieces.push((
                curve,
                range,
                edge.orientation() == Orientation::Reversed,
                length,
            ));
        }
        if total <= tol.confusion() {
            ogeom_bail!(Construction, "a rail has no length");
        }
        let rail = Self { pieces, total };
        let (head, tail) = (rail.at(0.0, tol)?.0, rail.at(1.0, tol)?.0);
        if head.distance(tail) <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "a two-rail sweep runs along open rails; a rail is closed"
            );
        }
        Ok(rail)
    }

    /// The point and unit tangent at fraction `f` of the rail's length.
    fn at(&self, f: f64, tol: Tolerances) -> OgeomResult<(Point, Vector)> {
        let mut along = self.total * f.clamp(0.0, 1.0);
        let last = self.pieces.len() - 1;
        for (i, (curve, range, reversed, length)) in self.pieces.iter().enumerate() {
            if along > *length && i < last {
                along -= length;
                continue;
            }
            let into = along.min(*length);
            let t = if *reversed {
                ogeom_algo::parameter_at_length(curve, *range, length - into, tol)?
            } else {
                ogeom_algo::parameter_at_length(curve, *range, into, tol)?
            };
            let d = curve.d1_at(t, tol)?;
            if d.magnitude() <= tol.confusion() {
                ogeom_bail!(Construction, "a rail stands still at {t}");
            }
            let tangent = (if *reversed { -d } else { d }) / d.magnitude();
            return Ok((curve.point_at(t, tol)?, tangent));
        }
        Err(ogeom_err!(Construction, "a rail has no edges"))
    }
}
