//! A round that takes a strip beside its edge off whole (ADR 0001,
//! amendment "marched seats", a strip taken off whole): a bevel or a
//! smaller round on a crease, or a ring on a rim, narrower than the round
//! asked of the edge between it and one of its neighbours.
//!
//! No ball touches both of the edge's own faces there. The round is the
//! plain crease's: the strip taken off a copy of the solid
//! ([`ogeom_bool::remove_faces`], its two neighbours grown back to meet
//! along the crease the strip cut), and that crease rounded by the engine
//! at the asked radius. The copy rounded is the result where the strip lies
//! wholly in what the round takes (cut away at a convex crease, filled at
//! a concave one): every point of the strip then stands on or past the
//! rounded copy's boundary on the round's side. Elsewhere part of the strip
//! would stand, and the round is not this one.

use ogeom_algo::{Built, Containment};
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

thread_local! {
    /// Whether a strip is being taken off: the crease it restores is
    /// rounded by the engine, which must not take a strip off that crease
    /// in turn.
    static TAKING: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
}

/// The faces of `solid` holding `edge`.
fn faces_at(model: &Model, solid: &Shape, edge: &Shape) -> OgeomResult<Vec<Shape>> {
    Ok(
        ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)?
            .into_iter()
            .fold(Vec::new(), |mut found: Vec<Shape>, f| {
                if !found.iter().any(|g| g.is_same(&f)) {
                    found.push(f);
                }
                found
            }),
    )
}

/// An edge's length, summed over chords.
fn edge_length(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(0.0);
    };
    let mut length = 0.0;
    let mut last = curve.point_at(a, tol)?;
    for k in 1..=64 {
        let p = curve.point_at((b - a).mul_add(f64::from(k) / 64.0, a), tol)?;
        length += p.distance(last);
        last = p;
    }
    Ok(length)
}

/// The face across `strip` from `outer`: every edge of the strip but its
/// seams is shared with `outer`, with that face, or with a face at one of
/// the strip's ends, and the face across is the one sharing the most of
/// the strip's border after `outer`. Both sides hold some.
fn beyond(
    model: &Model,
    solid: &Shape,
    strip: &Shape,
    outer: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let mut shares: Vec<(Shape, f64)> = Vec::new();
    let mut on_outer = false;
    for e in explore_unique(model, strip, ShapeType::Edge)? {
        let holders = faces_at(model, solid, &e)?;
        let others: Vec<&Shape> = holders.iter().filter(|f| !f.is_same(strip)).collect();
        match others.as_slice() {
            [] => {}
            [f] if f.is_same(outer) => on_outer = true,
            [f] => {
                let length = edge_length(model, &e, tol)?;
                match shares.iter_mut().find(|(g, _)| g.is_same(f)) {
                    Some((_, total)) => *total += length,
                    None => shares.push(((*f).clone(), length)),
                }
            }
            _ => return Ok(None),
        }
    }
    if !on_outer {
        return Ok(None);
    }
    Ok(shares
        .into_iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(f, _)| f))
}

/// The faces of `built` that `face` of the solid before it became.
fn images(built: &Built, face: &Shape) -> Vec<Shape> {
    let mut found: Vec<Shape> = built.history.modified(face).to_vec();
    if found.is_empty() && !built.history.is_deleted(face) {
        found.push(face.clone());
    }
    found
}

/// The edges of `copy` that two faces of `one` and `two` share.
fn shared_edges(model: &Model, one: &[Shape], two: &[Shape]) -> OgeomResult<Vec<Shape>> {
    let mut found: Vec<Shape> = Vec::new();
    for f in one {
        let edges = explore_unique(model, f, ShapeType::Edge)?;
        for g in two {
            for e in explore_unique(model, g, ShapeType::Edge)? {
                if edges.iter().any(|x| x.is_same(&e)) && !found.iter().any(|x| x.is_same(&e)) {
                    found.push(e);
                }
            }
        }
    }
    Ok(found)
}

/// Points spread over `strip`, off its border: its mesh's triangles'
/// centroids.
fn spread(model: &Model, strip: &Shape, tol: Tolerances) -> OgeomResult<Vec<ogeom_math::Point>> {
    let bounds = ogeom_algo::shape_bounds(model, strip, tol)?;
    let deflection = ogeom_mesh::Deflection {
        chord: bounds.diagonal() * 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let mesh = ogeom_mesh::triangulate_face(model, strip, deflection, tol)?;
    Ok(mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            ogeom_math::Point::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            )
        })
        .collect())
}

/// `solid`'s one `edge` rounded at `radius` with a strip beside it taken
/// off whole: the copy with the strip removed, its restored crease rounded
/// by the engine, where the strip lies wholly in what that round takes.
/// `None` where neither face beside the edge is such a strip, the strip
/// does not come off, the crease does not round, or part of the strip
/// would stand.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] where the work is cancelled.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if TAKING.with(core::cell::Cell::get) {
        return Ok(None);
    }
    struct Leave;
    impl Drop for Leave {
        fn drop(&mut self) {
            TAKING.with(|t| t.set(false));
        }
    }
    TAKING.with(|t| t.set(true));
    let _leave = Leave;
    let beside = faces_at(model, solid, edge)?;
    let [one, two] = beside.as_slice() else {
        return Ok(None);
    };
    for (strip, outer) in [(one, two), (two, one)] {
        let mut scratch = model.clone();
        match taken(&mut scratch, solid, edge, strip, outer, radius, tol) {
            Ok(Some(built)) => {
                *model = scratch;
                return Ok(Some(built));
            }
            Ok(None) => {}
            Err(e @ OgeomError::Cancelled) => return Err(e),
            Err(_) => {}
        }
    }
    Ok(None)
}

/// [`rounds`] with `strip` the face taken off, `outer` the edge's other
/// face.
fn taken(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    strip: &Shape,
    outer: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(across) = beyond(model, solid, strip, outer, tol)? else {
        return Ok(None);
    };
    let points = spread(model, strip, tol)?;
    if points.is_empty() {
        return Ok(None);
    }
    let removed = ogeom_bool::remove_faces(model, solid, std::slice::from_ref(strip), tol)?;
    // The crease the strip cut, restored: the edges the two faces share in
    // the copy that the solid did not hold (they may share others, such as
    // a bore's far rim).
    let before = explore_unique(model, solid, ShapeType::Edge)?;
    let creases: Vec<Shape> =
        shared_edges(model, &images(&removed, outer), &images(&removed, &across))?
            .into_iter()
            .filter(|e| !before.iter().any(|b| b.is_same(e)))
            .collect();
    if creases.is_empty() {
        return Ok(None);
    }
    let asked: Vec<(Shape, f64)> = creases.iter().map(|c| (c.clone(), radius)).collect();
    let rounded = super::marched_seats::outside_seats(|| {
        crate::fillet::rounds_by_engine(model, &removed.shape, &asked, tol)
    })?;
    let Some(rounded) = rounded else {
        return Ok(None);
    };
    // The strip cut material away where it lies inside the copy, and
    // filled where it lies outside: what the round takes must take it
    // all, every point of it on or past the rounded copy's boundary.
    let convex = matches!(
        ogeom_algo::classify_in_solid_exact(model, &removed.shape, points[0], tol)?,
        Containment::In
    );
    let kept = if convex {
        Containment::In
    } else {
        Containment::Out
    };
    for &p in &points {
        if ogeom_algo::classify_in_solid_exact(model, &rounded.shape, p, tol)? == kept {
            return Ok(None);
        }
    }
    if !ogeom_algo::check(model, &rounded.shape, tol)?.is_valid() {
        return Ok(None);
    }
    let mut history = removed.history.then(&rounded.history);
    history.delete(edge);
    let mut built = Built::new(rounded.shape, history);
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(Some(built))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
    use ogeom_math::{Direction, Frame, Point};

    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    /// A drum of radius 10 and 30 tall, bored across at z = 15 by a bore of
    /// radius 3 along y.
    fn holed_drum(m: &mut Model) -> Shape {
        let drum = ogeom_algo::make_cylinder(m, Frame::WORLD, 10.0, 30.0, T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(0.0, -20.0, 15.0), Direction::Y, Direction::X, T).unwrap();
        let bore = ogeom_algo::make_cylinder(m, at, 3.0, 40.0, T)
            .unwrap()
            .shape;
        ogeom_bool::cut(m, &drum, &bore, T).unwrap().shape
    }

    /// The edge of `s` passing nearest `q`.
    fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
        let off = |e: &Shape| {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return f64::INFINITY;
            };
            (0..=64)
                .map(|i| {
                    c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                        .unwrap()
                        .distance(q)
                })
                .fold(f64::INFINITY, f64::min)
        };
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .min_by(|a, b| off(a).total_cmp(&off(b)))
            .unwrap()
    }

    fn volume(m: &Model, s: &Shape) -> f64 {
        ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::with_chord(1e-3).unwrap(), T)
            .unwrap()
            .mass
    }

    /// The holed drum with its bore's rim at y > 0 bevelled 0.5, the
    /// bevel's face, and the bevel's edge on the bore.
    fn bevelled(m: &mut Model) -> (Shape, Shape, Shape) {
        let plain = holed_drum(m);
        let rim = edge_near(m, &plain, Point::new(0.0, 9.5, 18.0));
        let bevelled = crate::chamfer::chamfer_edges(m, &plain, &[rim], 0.5, T)
            .unwrap()
            .shape;
        let edge = edge_near(m, &bevelled, Point::new(0.0, 9.0, 18.0));
        let bevel = faces_at(m, &bevelled, &edge)
            .unwrap()
            .into_iter()
            .find(|f| {
                matches!(
                    m.node(f)
                        .and_then(|n| n.data().as_face())
                        .and_then(|d| m.geometry().surface(d.surface)),
                    Some(ogeom_geom::SurfaceGeometry::BSpline(_))
                )
            })
            .expect("the bevel beside the edge");
        (bevelled, bevel, edge)
    }

    /// The bevel taken off: the copy is valid, the plain holed drum's
    /// volume, and its drum and bore meet along the restored rim (a closed
    /// clamped spline, cut where the faces' seams reach it) and along the
    /// bore's far rim, which the solid already held.
    #[test]
    fn a_bevel_taken_off_a_bore_rim_restores_the_plain_drum() {
        let mut m = Model::new();
        let plain = holed_drum(&mut m);
        let (bevelled, bevel, _) = bevelled(&mut m);
        let copy = ogeom_bool::remove_faces(&mut m, &bevelled, &[bevel], T).unwrap();
        assert!(ogeom_algo::check(&m, &copy.shape, T).unwrap().is_valid());
        let (got, want) = (volume(&m, &copy.shape), volume(&m, &plain));
        assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
    }

    /// The bevelled rim rounded at 2 (the ball sets back past the bevel)
    /// with the bevel taken off whole: the plain rim's round.
    #[test]
    fn a_bevelled_bore_rim_rounded_past_its_bevel_is_the_plain_rim_s_round() {
        let mut m = Model::new();
        let plain = holed_drum(&mut m);
        let rim = edge_near(&m, &plain, Point::new(0.0, 9.5, 18.0));
        let round = crate::fillet::fillet_edges(&mut m, &plain, &[rim], 2.0, T)
            .unwrap()
            .shape;
        let want = volume(&m, &round);
        let (bevelled, _, edge) = bevelled(&mut m);
        let built = rounds(&mut m, &bevelled, &edge, 2.0, T)
            .unwrap()
            .expect("the bevel taken off and the restored rim rounded");
        assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
        let got = volume(&m, &built.shape);
        assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
    }

    /// The boolean alone: the restored rim's marched band at 2, as the
    /// engine's seats build it from either of the rim's two edges, cut
    /// from the copy the bevel left. The copy's regrown drum and bore each
    /// walk their seam both ways, and the cut needs both columns of it.
    #[test]
    fn the_restored_rim_s_band_cuts_the_copy_the_bevel_left() {
        let mut m = Model::new();
        let (bevelled, bevel, _) = bevelled(&mut m);
        let before = explore_unique(&m, &bevelled, ShapeType::Edge).unwrap();
        let copy = ogeom_bool::remove_faces(&mut m, &bevelled, &[bevel], T)
            .unwrap()
            .shape;
        // The restored rim: the copy's spline edges the bevelled solid did
        // not hold (a regrown face's seam is a line).
        let restored: Vec<Shape> = explore_unique(&m, &copy, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| !before.iter().any(|b| b.is_same(e)))
            .filter(|e| {
                matches!(
                    ogeom_algo::placed_edge_curve(&m, e, T).unwrap(),
                    Some((ogeom_geom::Curve::BSpline(_), _))
                )
            })
            .collect();
        assert_eq!(restored.len(), 2, "the rim in two edges");
        for start in &restored {
            let wedges =
                super::super::marched_seats::wedges(&mut m, &copy, &[(start.clone(), 2.0)], T)
                    .unwrap()
                    .expect("the march's band");
            for wedge in wedges {
                let cut = ogeom_bool::cut(&mut m, &copy, &wedge.solid, T)
                    .unwrap_or_else(|e| panic!("the band's cut: {e}"));
                assert!(ogeom_algo::check(&m, &cut.shape, T).unwrap().is_valid());
            }
        }
    }
}
