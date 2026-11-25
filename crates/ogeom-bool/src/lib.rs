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
mod tests;
