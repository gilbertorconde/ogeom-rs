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

use crate::asym::{RoundProfile, conic_pieces, legs, straight_seat};
use crate::law::{Profile, RadiusLaw};
use crate::support::{face_from_edges, segment_between};
use ogeom_algo::{Built, make_edge_between, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{BSplineCurve, BSplineSurface, Curve, PlaneSurface, SurfaceGeometry};
use ogeom_math::bspline::Weighted;
use ogeom_math::{ControlGrid, Direction, KnotVector, Plane, Point, Point2, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

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
            let mut built = varying(model, solid, &parts[0], &parts[1], &laws, profile, tol)?;
            crate::support::credit_new_faces(model, solid, &parts[0], &mut built)?;
            Ok(built)
        },
    )
}

/// How many stations along the edge a varying band's contacts are held to
/// their faces at, ends included.
const FIT_STATIONS: u32 = 16;

/// [`fillet_edge_distances_variable`] on a shape every node below which is
/// its own, the laws in the order the request names the faces.
fn varying(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    laws: &[Profile; 2],
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (seat, named) = straight_seat(model, solid, edge, face, tol)?;
    // The laws in the seat's face order.
    let law = |i: usize| if i == named { &laws[0] } else { &laws[1] };
    let f = legs(&seat, tol)?;
    let w = profile.weight(f[0].dot(f[1]).clamp(-1.0, 1.0))?;
    square_ends(model, solid, edge, &seat, tol)?;
    let length = seat.start.distance(seat.end);
    let crease = |s: f64| seat.start + seat.along * (length * s);
    for k in 0..=FIT_STATIONS {
        let s = f64::from(k) / f64::from(FIT_STATIONS);
        for (i, (leg, host)) in f.iter().zip(&seat.faces).enumerate() {
            let distance = law(i).at(s).0;
            let contact = crease(s) + *leg * distance;
            if !crate::support::on_face_or_boundary(model, solid, host, contact, tol)? {
                ogeom_bail!(
                    Construction,
                    "an asymmetric round sets back {distance} along a face of the edge, at \
                     {s:.3} of it, past that face's side; see docs/CAPABILITIES.md, \
                     fillet.asymmetric-rounds"
                );
            }
        }
    }

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
    let sign = if seat.convex { 1.0 } else { -1.0 };

    // Along the edge: the union of both laws' pieces, each law one cubic
    // over each stretch, in Bezier form.
    let pieces = [law(first).pieces(tol), law(second).pieces(tol)];
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
    let surface_id = model.geometry_mut().add_surface(surface.clone());

    let apex = [seat.start, seat.end];
    let contact = |end: usize, side: usize| {
        rows[if end == 0 { 0 } else { nu - 1 }][if side == 0 { 0 } else { nv - 1 }]
    };
    let ca = [contact(0, 0), contact(1, 0)];
    let cb = [contact(0, 1), contact(1, 1)];
    let va = [
        make_vertex(model, ca[0]).shape,
        make_vertex(model, ca[1]).shape,
    ];
    let vb = [
        make_vertex(model, cb[0]).shape,
        make_vertex(model, cb[1]).shape,
    ];
    let vapex = [
        make_vertex(model, apex[0]).shape,
        make_vertex(model, apex[1]).shape,
    ];

    // The contact curves: each face's column of the net, on the band's
    // rows `v = 0` and `v = 1`.
    let rail =
        |model: &mut Model, column: usize, ends: &[Shape; 2], row_v: f64| -> OgeomResult<Shape> {
            let control = rows.iter().map(|r| r[column]).collect::<Vec<_>>();
            let curve = Curve::BSpline(BSplineCurve::new(
                KnotVector::new(u_knots.clone(), 3)?,
                control,
                tol,
            )?);
            let built = make_edge_between(model, curve, (0.0, 1.0), &ends[0], &ends[1], tol)?.shape;
            let pcurve = ogeom_geom::BSpline2d::new(
                KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?,
                vec![Point2::new(0.0, row_v), Point2::new(1.0, row_v)],
                tol,
            )?;
            ogeom_algo::attach_pcurve(
                model,
                &built,
                pcurve.into(),
                surface_id,
                ogeom_topo::Location::identity(),
                (0.0, 1.0),
            )?;
            Ok(built)
        };
    let rail_a = rail(model, 0, &va, 0.0)?;
    let rail_b = rail(model, nv - 1, &vb, 1.0)?;

    // The end sections: the net's first and last rows, on the band's
    // columns `u = 0` and `u = 1`.
    let arc = |model: &mut Model, end: usize| -> OgeomResult<Shape> {
        let r = &rows[if end == 0 { 0 } else { nu - 1 }];
        let control = r
            .iter()
            .enumerate()
            .map(|(j, q)| Weighted::new(*q, if j % 2 == 1 { weight } else { 1.0 }, tol))
            .collect::<OgeomResult<Vec<_>>>()?;
        let curve = Curve::BSpline(BSplineCurve::rational(
            KnotVector::new(v_knots.clone(), 2)?,
            control,
        )?);
        let built = make_edge_between(model, curve, (0.0, 1.0), &va[end], &vb[end], tol)?.shape;
        #[allow(clippy::cast_precision_loss, reason = "0 or 1")]
        let column = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(end as f64, 0.0), ogeom_math::Direction2::Y),
            0.0,
            1.0,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &built,
            column.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            (0.0, 1.0),
        )?;
        Ok(built)
    };
    let arc0 = arc(model, 0)?;
    let arc1 = arc(model, 1)?;

    // The band, facing away from the crease, which the wedge holds.
    let band = {
        use ogeom_geom::Surface as _;
        let wire = ogeom_algo::make_wire(
            model,
            &[
                rail_a.clone(),
                arc1.clone(),
                rail_b.reversed(),
                arc0.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
        let middle = surface.point_at(0.5, 0.5, tol)?;
        let (du, dv) = surface.d1_at(0.5, 0.5, tol)?;
        if du.cross(dv).dot(middle - crease(0.5)) >= 0.0 {
            face
        } else {
            face.reversed()
        }
    };

    let crease_edge = segment_between(model, (&vapex[0], apex[0]), (&vapex[1], apex[1]), tol)?;
    // A leg on each face, between the crease and its contact curve, facing
    // as the face does on a convex edge and against it on a concave one.
    let leg = |model: &mut Model,
               side: usize,
               rail: &Shape,
               contacts: &[Point; 2],
               ends: &[Shape; 2]|
     -> OgeomResult<Shape> {
        let outward = seat.normals[side] * sign;
        let plane = Plane::through(apex[0], Direction::new(outward, tol)?);
        let reach =
            (length + contacts[0].distance(apex[0]) + contacts[1].distance(apex[1]) + 1.0) * 2.0;
        let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
        let at0 = segment_between(model, (&vapex[0], apex[0]), (&ends[0], contacts[0]), tol)?;
        let at1 = segment_between(model, (&ends[1], contacts[1]), (&vapex[1], apex[1]), tol)?;
        let edges = vec![at0, rail.clone(), at1, crease_edge.reversed()];
        // Crease start, contact, contact, crease end: walked back where
        // that turns against `outward`.
        let turn = (contacts[0] - apex[0]).cross(apex[1] - apex[0]);
        let edges: Vec<Shape> = if turn.dot(outward) < 0.0 {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        };
        face_from_edges(model, surface.into(), &edges, tol)
    };
    let leg_a = leg(model, first, &rail_a, &ca, &va)?;
    let leg_b = leg(model, second, &rail_b, &cb, &vb)?;

    let cap = |model: &mut Model, end: usize, arc: &Shape, outward: Vector| -> OgeomResult<Shape> {
        let plane = Plane::through(apex[end], Direction::new(outward, tol)?);
        let reach = (apex[end].distance(ca[end]) + apex[end].distance(cb[end]) + 1.0) * 2.0;
        let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
        let edges = vec![
            segment_between(model, (&vapex[end], apex[end]), (&va[end], ca[end]), tol)?,
            arc.clone(),
            segment_between(model, (&vb[end], cb[end]), (&vapex[end], apex[end]), tol)?,
        ];
        let turn = (ca[end] - apex[end]).cross(cb[end] - apex[end]);
        let edges: Vec<Shape> = if turn.dot(outward) < 0.0 {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        };
        face_from_edges(model, surface.into(), &edges, tol)
    };
    let cap0 = cap(model, 0, &arc0, -seat.along)?;
    let cap1 = cap(model, 1, &arc1, seat.along)?;

    let faces = [leg_a, leg_b, cap0, cap1, band];
    let wedge = crate::support::wedge_solid(model, &faces, false, tol)?;
    crate::support::apply_wedge_solid(model, solid, Some(edge), &wedge, !seat.convex, tol)
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

/// Refuses by name where an end of the straight seat's edge meets a face
/// of `solid`, other than the edge's own two, whose plane is not square to
/// the edge: a band capped square there would leave or cut what that face
/// bounds.
fn square_ends(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &crate::support::Seat,
    tol: Tolerances,
) -> OgeomResult<()> {
    let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(());
    };
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if seat.faces.iter().any(|f| f.is_same(&face)) {
            continue;
        }
        let vertices = explore_unique(model, &face, ShapeType::Vertex)?;
        let meets = vertices.iter().any(|v| {
            crate::support::same_occurrence(model, v, &a, tol)
                || crate::support::same_occurrence(model, v, &b, tol)
        });
        if !meets {
            continue;
        }
        let square = match face_normal(model, &face, tol)? {
            Some(normal) => normal.cross(seat.along).magnitude() <= tol.angular(),
            None => false,
        };
        if !square {
            ogeom_bail!(
                Construction,
                "an asymmetric round whose distances vary is capped square at its edge's \
                 ends, and an end of this edge meets a face not square to it; see \
                 docs/CAPABILITIES.md, fillet.asymmetric-rounds"
            );
        }
    }
    Ok(())
}

/// A planar face's unit normal, placed; `None` for a face on any other
/// surface.
fn face_normal(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    use ogeom_geom::Transformable as _;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placed = surface
        .clone()
        .transformed(&face.transform(model.datums())?, tol)?;
    Ok(match placed {
        SurfaceGeometry::Plane(p) => Some(p.plane().normal().vector()),
        _ => None,
    })
}
