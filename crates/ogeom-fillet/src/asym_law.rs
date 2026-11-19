//! Asymmetric rounds whose distances vary along a straight edge.
//!
//! On a straight edge between two planes the section at a fraction `s` of
//! the edge is the conic arc of [`crate::asym`] with control points
//! `A(s) + d1(s) f1`, `A(s)`, `A(s) + d2(s) f2`: the crease `A` runs
//! linearly along the edge, the legs `f1`, `f2` and the weight stay as they
//! are, and each distance is its law's cubic Hermite piece. Every control
//! point of the section, and of the pieces [`crate::asym::conic_pieces`]
//! splits it into, is affine in the two distances, so over each stretch
//! where both laws are one cubic the band is the tensor product of their
//! Bezier forms with the section's: a rational B-spline, cubic along the
//! edge, exact. The contact lines are the laws' own cubics pushed along
//! each face, so the band meets each face with its tangent plane all
//! along.

use crate::asym::{RoundProfile, conic_pieces};
use crate::law::{Profile, RadiusLaw};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{BSplineSurface, SurfaceGeometry};
use ogeom_math::bspline::Weighted;
use ogeom_math::{ControlGrid, KnotVector, Point, Vector};
use ogeom_topo::{Model, Shape};

/// Round a straight edge between two planar faces with a band that meets
/// `face` at the distance `on_face` gives and the edge's other face at the
/// distance `on_other` gives, each a law read by the fraction of the edge
/// from its start, its section `profile` between the two contacts.
///
/// The band is exact: over each stretch where both laws are one cubic it
/// is the section's rational arc swept with its control points moving
/// along the laws' own cubics, and it meets each face with its tangent
/// plane along a contact curve. Naming the other face with the laws
/// swapped builds the same solid, to the bit, and two laws holding one
/// distance each are [`fillet_edge_distances`](crate::fillet_edge_distances)
/// by the same call. The band's ends are square caps: each end of the edge
/// must meet only faces square to it.
///
/// What `solid` shares with other shapes is copied before it is worked on
/// ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where a law states no usable distance, `face` is not one of the edge's
/// two faces, a conic's `rho` is not strictly between 0 and 1, the edge is
/// not straight between two planes, an end of the edge meets a face not
/// square to it, or a contact sets back past its face's side; and as
/// [`fillet_edge_distances`](crate::fillet_edge_distances) where both laws
/// hold one distance.
#[allow(
    clippy::too_many_arguments,
    reason = "fillet_edge_distances' own, with laws for the distances"
)]
pub fn fillet_edge_distances_variable(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: &RadiusLaw,
    on_other: &RadiusLaw,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if let Some(built) = crate::compound::each_solid(
        model,
        solid,
        std::slice::from_ref(edge),
        |edge| edge,
        |model, part, _| {
            fillet_edge_distances_variable(model, part, edge, face, on_face, on_other, profile, tol)
        },
    )? {
        return Ok(built);
    }
    let laws = [
        Profile::new(on_face, false, tol)?,
        Profile::new(on_other, false, tol)?,
    ];
    if laws.iter().all(|law| law.level(tol)) {
        let [a, b] = [laws[0].at(0.0).0, laws[1].at(0.0).0];
        return crate::fillet_edge_distances(model, solid, edge, face, a, b, profile, tol);
    }
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &[edge.clone(), face.clone()],
        |model, solid, parts| {
            let mut built = crate::engine::varying_distances(
                model, solid, &parts[0], &parts[1], &laws, profile, tol,
            )?;
            crate::support::credit_new_faces(model, solid, &parts[0], &mut built)?;
            Ok(built)
        },
    )
}

/// The band of a varying asymmetric round on a straight seat between two
/// planes: the faces in the order the section runs, the net's rows (one
/// per control value along the edge, each the section's pieces across),
/// the pieces' middle weight, both knot vectors, and the surface.
pub(crate) struct VaryingNet {
    /// The seat's face the section starts on.
    pub first: usize,
    /// The control points, row by row along the edge.
    pub rows: Vec<Vec<Point>>,
    /// The knots along the edge, on the edge's fraction.
    pub u_knots: Vec<f64>,
    /// The band.
    pub surface: SurfaceGeometry,
}

/// [`VaryingNet`] of `seat` whose distance on each face, in the seat's
/// face order, follows `laws` by the fraction of the edge, its legs `f`
/// and its section's weight `w`.
pub(crate) fn varying_net(
    seat: &crate::support::Seat,
    laws: [&Profile; 2],
    f: [Vector; 2],
    w: f64,
    tol: Tolerances,
) -> OgeomResult<VaryingNet> {
    let length = seat.start.distance(seat.end);
    let crease = |s: f64| seat.start + seat.along * (length * s);
    // The faces in the order that turns positively about the crease, so the
    // section runs from the first face's contact to the second's.
    let (first, second) = if seat
        .along
        .dot(seat.normals[0].cross(seat.normals[1]))
        .is_sign_positive()
    {
        (0, 1)
    } else {
        (1, 0)
    };

    // Along the edge: the union of both laws' pieces, each law one cubic
    // over each stretch, in Bezier form.
    let pieces = [laws[first].pieces(tol), laws[second].pieces(tol)];
    let mut breaks: Vec<f64> = pieces
        .iter()
        .flat_map(|p| p.iter().flat_map(|(a, b, _)| [*a, *b]))
        .collect();
    breaks.sort_by(f64::total_cmp);
    breaks.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());
    let windowed = |which: usize, a: f64, b: f64| -> [f64; 4] {
        let middle = 0.5 * (a + b);
        let &(x0, x1, bezier) = pieces[which]
            .iter()
            .find(|(x0, x1, _)| *x0 <= middle && middle <= *x1)
            .unwrap_or(&pieces[which][0]);
        bezier_window(bezier, (x0, x1), (a, b))
    };
    // The distances' control values along the edge, end to end, and the
    // crease's at the same places.
    let mut along: Vec<(f64, f64, Point)> = Vec::new();
    for (k, span) in breaks.windows(2).enumerate() {
        let (a, b) = (span[0], span[1]);
        let (da, db) = (windowed(0, a, b), windowed(1, a, b));
        for i in usize::from(k > 0)..4 {
            #[allow(clippy::cast_precision_loss, reason = "a Bezier index")]
            let t = a + (b - a) * i as f64 / 3.0;
            along.push((da[i], db[i], crease(t)));
        }
    }
    let mut u_knots = vec![0.0; 4];
    for &k in &breaks[1..breaks.len() - 1] {
        u_knots.extend([k; 3]);
    }
    u_knots.extend([1.0; 4]);

    // Across: the section's pieces, one row of the net per control value
    // along the edge.
    let row = |(da, db, apex): (f64, f64, Point)| {
        conic_pieces([apex + f[first] * da, apex, apex + f[second] * db], w)
    };
    let (_, weight, v_knots) = row(along[0]);
    let rows: Vec<Vec<Point>> = along.iter().map(|&c| row(c).0).collect();
    let (nu, nv) = (rows.len(), rows[0].len());
    let mut net = Vec::with_capacity(nu * nv);
    for r in &rows {
        for (j, q) in r.iter().enumerate() {
            net.push(Weighted::new(
                *q,
                if j % 2 == 1 { weight } else { 1.0 },
                tol,
            )?);
        }
    }
    let surface: SurfaceGeometry = BSplineSurface::rational(
        KnotVector::new(u_knots.clone(), 3)?,
        KnotVector::new(v_knots.clone(), 2)?,
        ControlGrid::new(net, nu, nv)?,
    )?
    .into();
    Ok(VaryingNet {
        first,
        rows,
        u_knots,
        surface,
    })
}

/// The cubic Bezier `b` over `(x0, x1)` restricted to `(a, b)` within it,
/// in Bezier form there: its blossom at the window's ends.
fn bezier_window(bezier: [f64; 4], (x0, x1): (f64, f64), (a, b): (f64, f64)) -> [f64; 4] {
    let h = x1 - x0;
    let (ta, tb) = ((a - x0) / h, (b - x0) / h);
    // The blossom of the cubic at three parameters, by de Casteljau with a
    // different parameter at each level.
    let blossom = |t: [f64; 3]| -> f64 {
        let mut p = bezier;
        for (level, &x) in t.iter().enumerate() {
            for i in 0..3 - level {
                p[i] = p[i] * (1.0 - x) + p[i + 1] * x;
            }
        }
        p[0]
    };
    [
        blossom([ta, ta, ta]),
        blossom([ta, ta, tb]),
        blossom([ta, tb, tb]),
        blossom([tb, tb, tb]),
    ]
}
