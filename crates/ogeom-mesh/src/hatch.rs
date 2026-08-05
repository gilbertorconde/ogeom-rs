//! Hatching in parametric space: parallel lines clipped to a face's trim.
//!
//! The rings are the same chart boundary every other consumer walks —
//! [`face_boundary`] — so the hatch and
//! the triangulation cannot disagree about where the face is. The lines run
//! at an angle in the chart, spaced evenly, and each is cut to the inside
//! intervals by even-odd crossing counting; holes split segments the same
//! way they split everything else.
//!
//! Scanlines sit at half-spacing offsets, so a boundary lying exactly on a
//! round coordinate — every axis-aligned face — is not grazed.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::Point2;
use ogeom_topo::{Model, Shape};

use crate::discretize::Deflection;
use crate::triangulate::face_boundary;

/// Hatch a face's chart with parallel lines.
///
/// `angle` is the line direction in the chart, radians from the `u` axis;
/// `spacing` the perpendicular distance between lines. Returns the clipped
/// segments in chart coordinates, each as its two endpoints — lifting them
/// to space is the surface's own `point_at`, and which chart step is fine
/// enough for that lift is the caller's deflection question, not this one.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// spacing is not finite and positive, plus whatever the boundary walk
/// refuses.
pub fn hatch_face(
    model: &Model,
    face: &Shape,
    angle: f64,
    spacing: f64,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<[Point2; 2]>> {
    if !spacing.is_finite() || spacing <= tol.confusion() {
        ogeom_bail!(Construction, "a hatch spacing of {spacing} draws nothing");
    }
    if !angle.is_finite() {
        ogeom_bail!(Construction, "a hatch angle of {angle} is not a direction");
    }
    let rings = face_boundary(model, face, deflection, tol)?;

    // Rotate the rings so the hatch runs horizontal, scan, rotate back.
    let (sin, cos) = angle.sin_cos();
    let turn = |p: Point2| Point2::new(p.x.mul_add(cos, p.y * sin), p.y.mul_add(cos, -(p.x * sin)));
    let back = |p: Point2| Point2::new(p.x.mul_add(cos, -(p.y * sin)), p.y.mul_add(cos, p.x * sin));
    let turned: Vec<Vec<Point2>> = rings
        .iter()
        .map(|ring| ring.iter().map(|p| turn(*p)).collect())
        .collect();

    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for ring in &turned {
        for p in ring {
            lo = lo.min(p.y);
            hi = hi.max(p.y);
        }
    }
    if !lo.is_finite() || hi <= lo {
        return Ok(Vec::new());
    }

    let mut out = Vec::new();
    let first = (lo / spacing).floor();
    let mut k = first;
    loop {
        let level = (k + 0.5) * spacing;
        k += 1.0;
        if level >= hi {
            break;
        }
        if level <= lo {
            continue;
        }
        // Every crossing of every ring, in order: even-odd pairs are the
        // inside intervals, holes included by the same counting.
        let mut crossings: Vec<f64> = Vec::new();
        for ring in &turned {
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                if (a.y > level) != (b.y > level) {
                    crossings.push((b.x - a.x).mul_add((level - a.y) / (b.y - a.y), a.x));
                }
            }
        }
        crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
        for pair in crossings.chunks_exact(2) {
            if pair[1] - pair[0] > tol.parametric() {
                out.push([
                    back(Point2::new(pair[0], level)),
                    back(Point2::new(pair[1], level)),
                ]);
            }
        }
    }
    Ok(out)
}
