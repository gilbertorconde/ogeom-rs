//! Where a curve first (or last) meets a shape: the face it meets, the
//! point and the curve's parameter there.
//!
//! The meeting is solved against each face's true surface and kept where
//! it falls inside the face's trim. It is not read off a tessellation, so a
//! hit is as exact as the curve/surface intersection, not off by a chord.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, LineCurve, Transformable as _, TrimmedCurve};
use ogeom_math::Point;
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};

use crate::classify::{Containment, FaceTrim};

/// Which meeting along the curve [`curve_hit`] answers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hit {
    /// The meeting at the smallest parameter.
    First,
    /// The meeting at the largest parameter.
    Last,
}

/// Where a curve meets a face of a shape.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveHit {
    /// The face met, as the shape holds it.
    pub face: Shape,
    /// Where, on the curve.
    pub point: Point,
    /// The curve's parameter there.
    pub parameter: f64,
}

/// Where `curve` over `range` first (or last, by `which`) meets a face of
/// `shape`: a face, or any shape whose faces all count.
///
/// A face is met where the curve crosses or touches its surface at a point
/// inside its trim, the trim's own boundary included: a meeting on the
/// carrier surface outside the trim is no meeting. A stretch of the curve
/// lying in a face's surface does not meet that face: the curve runs along
/// it and crosses nothing, and the faces where the stretch leaves the
/// surface are met instead. Where two faces are met at one parameter (an
/// edge they share), the first of them in the shape's order is named.
///
/// Each meeting is solved against the face's surface, placed, in closed form
/// for a line against a plane, a sphere or a cylinder and by Newton on the
/// curve/surface system otherwise; the trim is then asked of the point's
/// parameters on that surface. Placements are taken as rigid motions, under
/// which the curve's parameter is kept.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `range` is
/// empty or leaves the curve's domain;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a face's
/// surface is not in the model; as the intersection and the trim otherwise.
pub fn curve_hit(
    model: &Model,
    shape: &Shape,
    curve: &Curve,
    range: (f64, f64),
    which: Hit,
    tol: Tolerances,
) -> OgeomResult<Option<CurveHit>> {
    let piece = restricted(curve, range, tol)?;
    let mut best: Option<CurveHit> = None;
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        if let Some(hit) = face_hit(model, &face, &piece, which, tol)? {
            let better = best.as_ref().is_none_or(|held| match which {
                Hit::First => hit.parameter < held.parameter,
                Hit::Last => hit.parameter > held.parameter,
            });
            if better {
                best = Some(hit);
            }
        }
    }
    Ok(best)
}

/// The meeting of `curve`, already restricted to its range, with one face.
fn face_hit(
    model: &Model,
    face: &Shape,
    curve: &Curve,
    which: Hit,
    tol: Tolerances,
) -> OgeomResult<Option<CurveHit>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "a face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "a face refers to a surface not in this model");
    };
    // The trim lives in the surface's own chart, so the curve is taken
    // into the face's frame rather than the surface out of it.
    let placement = face.transform(model.datums())?;
    let local = if placement.is_equal(&ogeom_math::Transform::IDENTITY, tol) {
        curve.clone()
    } else {
        curve.transformed(&placement.inverse()?, tol)?
    };
    let found = ogeom_intersect::intersect_curve_surface(
        &local,
        surface,
        ogeom_intersect::CurveSurfaceOptions::for_tolerances(tol),
        tol,
    )?;
    if found.crossings.is_empty() {
        return Ok(None);
    }
    let reach = tol.confusion().max(data.tolerance.get());
    let trim = FaceTrim::of(model, face, trim_deflection(model, face, tol)?, tol)?;
    let mut order: Vec<_> = found.crossings.iter().collect();
    if which == Hit::Last {
        order.reverse();
    }
    for crossing in order {
        if trim.holds(model, crossing.on_surface, reach, tol)? != Containment::Out {
            return Ok(Some(CurveHit {
                face: face.clone(),
                point: placement.apply(crossing.point),
                parameter: crossing.on_curve,
            }));
        }
    }
    Ok(None)
}

/// The deflection a face's trim is drawn at for asking where a hit falls:
/// a hundred-thousandth of the face's size, never under the confusion
/// distance.
fn trim_deflection(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Deflection> {
    let size = crate::measure::shape_bounds(model, face, tol)?.diagonal();
    Ok(Deflection {
        chord: (size * 1e-5).max(tol.confusion()),
        angular: 0.05,
        ..Deflection::for_tolerances(tol)
    })
}

/// `curve` over `range` alone: a line kept a line, so the closed forms
/// answer it.
fn restricted(curve: &Curve, range: (f64, f64), tol: Tolerances) -> OgeomResult<Curve> {
    let (a, b) = range;
    let (lo, hi) = curve.domain();
    if !a.is_finite() || !b.is_finite() || b <= a + tol.parametric() {
        ogeom_bail!(Domain, "the range [{a}, {b}] is empty");
    }
    if !curve.is_periodic() && (a < lo - tol.parametric() || b > hi + tol.parametric()) {
        ogeom_bail!(
            Domain,
            "the range [{a}, {b}] leaves the curve's domain [{lo}, {hi}]"
        );
    }
    Ok(match curve {
        Curve::Line(line) => LineCurve::over(line.axis(), a, b)?.into(),
        other => TrimmedCurve::new(other.clone(), a, b, tol)?.into(),
    })
}
