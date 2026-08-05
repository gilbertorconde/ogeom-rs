//! Surface filling: a face fitted over the region four edges bound.
//!
//! The construction is the transfinite Coons blend of the four boundary
//! curves — which interpolates them exactly — sampled and fitted through
//! the grid machinery, error reported. What the caller gets is a *natural*
//! face over the fitted patch: the patch's own chart rectangle is the trim,
//! and the patch boundary stands within the stated fit tolerance of the
//! edges it was asked to fill.

use ogeom_algo::{Built, History, make_natural_face};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, Reversible as _, SurfaceGeometry, TrimmedCurve};
use ogeom_topo::{EdgeRepr, Model, Shape, ShapeType};

/// Fill the loop `edges` bound with a fitted patch face.
///
/// The four edges must chain head to tail into a closed loop, in order —
/// the first runs along the patch's `u` direction. `samples` controls the
/// Coons sampling per direction and `tolerance` the fit target; the fit
/// that cannot meet it refuses with the error it reached.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// edges do not chain into a loop, an edge carries no curve, or the fit
/// misses the tolerance.
pub fn make_filling(
    model: &mut Model,
    edges: &[Shape; 4],
    samples: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut curves: Vec<Curve> = Vec::with_capacity(4);
    for edge in edges {
        if model.kind_of(edge)? != ShapeType::Edge {
            ogeom_bail!(Construction, "a filling is bounded by edges");
        }
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "edge holds no edge data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "a filling edge needs a 3D curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Construction, "edge refers to a curve not in this model");
        };
        curves.push(Curve::Trimmed(Box::new(TrimmedCurve::new(
            geometry.clone(),
            range.0,
            range.1,
            tol,
        )?)));
    }

    // Chain head to tail, reversing edges whose stored direction runs
    // against the loop.
    let slack = tol.confusion() * 1e3;
    let start = curves[0].start(tol)?;
    let mut cursor = curves[0].end(tol)?;
    for curve in curves.iter_mut().skip(1) {
        if curve.start(tol)?.distance(cursor) > slack {
            if curve.end(tol)?.distance(cursor) > slack {
                ogeom_bail!(
                    Construction,
                    "the edges do not chain into a loop; a gap of {} remains",
                    curve
                        .start(tol)?
                        .distance(cursor)
                        .min(curve.end(tol)?.distance(cursor))
                );
            }
            *curve = curve.reversed();
        }
        cursor = curve.end(tol)?;
    }
    if cursor.distance(start) > slack {
        ogeom_bail!(
            Construction,
            "the loop does not close; a gap of {} remains",
            cursor.distance(start)
        );
    }

    // Loop order to Coons orientation: bottom with u, right with v, top
    // and left reversed back into the same directions.
    let bottom = curves[0].clone();
    let right = curves[1].clone();
    let top = curves[2].reversed();
    let left = curves[3].reversed();
    let fitted =
        ogeom_geom::fit::fill_boundary(&bottom, &top, &left, &right, samples, tolerance, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the filling reached {} against a target of {tolerance}",
            fitted.error
        );
    }

    let built = make_natural_face(model, SurfaceGeometry::BSpline(fitted.curve))?;
    let mut history = History::new();
    for edge in edges {
        history.modify(edge, built.shape.clone());
    }
    Ok(Built::new(built.shape, history))
}
