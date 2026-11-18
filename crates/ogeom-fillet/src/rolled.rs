//! The old band of a concave edge's round whose ball rolls on the far
//! edge of the narrow face beside it, built from the run
//! [`crate::engine::rolled`] reads.

use crate::engine::rolled::{Rolled, Run, rolled_run};
use crate::marched::ExactBand;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::Point;
use ogeom_topo::{Model, Shape};

/// Round `edge` with the ball resting on the far edge of the narrow face
/// beside it, where [`rolled_run`] finds that run; `None` where the edge
/// is not of that kind.
///
/// # Errors
///
/// As [`rolled_run`], and [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where the band's fit
/// stands further off the balls' envelope than its target.
pub(crate) fn rolled_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    match rolled_run(model, solid, edge, radius, tol)? {
        None => Ok(None),
        Some(Run::Rolled(run)) => rolled_band(model, solid, edge, &run, radius, tol).map(Some),
        // A ball on a shelf, a rim, carried across a change of support or
        // down a staircase is the engine's or refused by name.
        Some(Run::Shelf | Run::Rim | Run::Carried | Run::Stepped) => Err(crate::fillet::declined()),
    }
}

/// The fill of a run found by [`rolled_run`], applied to the solid.
pub(crate) fn rolled_band(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    run: &Rolled,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let path = run.path;
    let (sa, sb) = run.span;
    let n = run.blend.len();
    // The distance from a point to the envelope: from its nearest centre on
    // the run, found by golden section about the nearest station.
    let off = |p: Point| -> f64 {
        let nearest = (0..n)
            .min_by(|&i, &j| {
                p.distance(run.blend.spine[i])
                    .total_cmp(&p.distance(run.blend.spine[j]))
            })
            .unwrap_or(0);
        #[allow(clippy::cast_precision_loss)]
        let at = |i: f64| sa + (sb - sa) * i / ((n - 1) as f64);
        #[allow(clippy::cast_precision_loss)]
        let (mut a, mut b) = (at(nearest as f64 - 1.0), at(nearest as f64 + 1.0));
        let f = |s: f64| p.distance(path.centre(s));
        let g = 0.5 * (5.0_f64.sqrt() - 1.0);
        let (mut x1, mut x2) = (b - g * (b - a), a + g * (b - a));
        let (mut f1, mut f2) = (f(x1), f(x2));
        for _ in 0..80 {
            if f1 < f2 {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - g * (b - a);
                f1 = f(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + g * (b - a);
                f2 = f(x2);
            }
        }
        (f1.min(f2) - radius).abs()
    };
    let exact = ExactBand {
        rails: run.rails,
        off: &off,
        strips: &run.strips,
    };
    let [(first, first_sign), (second, second_sign)] = &run.hosts;
    crate::marched::build_held_band(
        model,
        solid,
        Some(edge),
        &run.blend,
        &run.guide,
        [(first, *first_sign), (second, *second_sign)],
        radius,
        false,
        [false, false],
        Some(&exact),
        tol,
    )
}
