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
mod fine;
mod fuse_all;
mod gather;
mod general_fuse;
mod half_space;
mod mesh_region;
mod pieces;
mod rebuild;
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
use fine::*;
use gather::*;
use general_fuse::*;
use pieces::*;
use rebuild::*;
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
