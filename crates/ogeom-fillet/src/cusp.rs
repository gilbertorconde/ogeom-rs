//! A cusp corner: where a flat side runs tangent into a cylindrical wall,
//! the flat top between them narrows to a cusp at the tangent point, and a
//! request rounds the top's edge with the side (convex) and its edge with
//! the wall (concave) together.
//!
//! Within `w = 2r` of the cusp the two bands no longer fit side by side on
//! the top, and the sound round there is an S: off the side in a convex
//! turn, onto the wall in a concave one, flattening into the walls' common
//! tangent plane at the tangent point. The bands are built, a cell is cut
//! round the tangent point out to where the cusp is `4r` wide, and the cell
//! is closed again: the side and the wall each down to (or up from) a trim
//! along the S's own tangent line on it, and on each side of the tangent
//! point one patch, G1 to both bands' ends, the top, the side and the wall,
//! held inside the solid and the wall's concave fill.
//!
//! On the side the trim is the S's own foot, section by section along the
//! side. On the wall it leaves the tangent point up the wall's line, along
//! which the wall's normal does not turn, and reaches the concave band's
//! rail along it: at the tangent point the side's normal holds still along
//! the side's trim, and a patch tangent to both there can only be one whose
//! normal also holds still along the wall's.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Continuity, Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Direction, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::fill_n::{FillBoundary, make_filling_n_region};
use crate::setback::point_of;

/// A cusp corner of a request, read off the solid it was asked on.
#[derive(Debug, Clone)]
pub(crate) struct CuspCorner {
    /// The tangent point, a vertex of the solid.
    pub vertex: Shape,
    corner: Point,
    /// The top's outward normal.
    up: Vector,
    /// The side's outward normal.
    out: Vector,
    /// The wall's radius; its axis runs along `up` at `radius` behind the
    /// side, through `corner - out * radius`.
    radius: f64,
    /// The directions along the top's edge with the side, away from the
    /// tangent point, along which a convex and a concave band leave it
    /// together.
    sides: Vec<Vector>,
    /// The bands' radius.
    blend: f64,
}

/// A face's surface where it stands.
fn placed(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<SurfaceGeometry>> {
    use ogeom_geom::Transformable as _;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    Ok(Some(
        surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?,
    ))
}

/// An edge's unit tangent at `at`, an end of it, pointing away from it.
fn leaving(model: &Model, edge: &Shape, at: Point, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(None);
    };
    let (a, b) = (curve.point_at(lo, tol)?, curve.point_at(hi, tol)?);
    let (t, sense) = if a.distance(at) <= b.distance(at) {
        (lo, 1.0)
    } else {
        (hi, -1.0)
    };
    Ok(Some(curve.d1_at(t, tol)?.normalized(tol)? * sense))
}

/// The cusp corner at `vertex`, where `at` (each request edge there with
/// its radius and kind) meet, or `None` where it is not one: the vertex's
/// faces a planar top, a planar side square to it and a cylindrical wall
/// about an axis square to the top, the side tangent to the wall at the
/// vertex; each convex edge a line between the top and the side, each
/// concave one an arc between the top and the wall, every one leaving the
/// vertex along with one of the other kind; one radius for all.
pub(crate) fn cusp_corner(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    at: &[(Shape, f64, bool)],
    tol: Tolerances,
) -> OgeomResult<Option<CuspCorner>> {
    let Some(blend) = at.first().map(|e| e.1) else {
        return Ok(None);
    };
    if at.iter().any(|e| (e.1 - blend).abs() > tol.confusion()) {
        return Ok(None);
    }
    let corner = point_of(model, vertex)?;
    let mut planes: Vec<(Shape, Vector)> = Vec::new();
    let mut walls: Vec<(Shape, ogeom_math::Cylinder)> = Vec::new();
    for f in explore_unique(model, solid, ShapeType::Face)? {
        if !explore_unique(model, &f, ShapeType::Vertex)?
            .iter()
            .any(|v| v.is_same(vertex))
        {
            continue;
        }
        match placed(model, &f, tol)? {
            Some(SurfaceGeometry::Plane(_)) => {
                let (_, outward) = ogeom_algo::face_normal(model, &f, tol)?;
                planes.push((f, outward));
            }
            Some(SurfaceGeometry::Cylinder(c)) => walls.push((f, c.cylinder())),
            _ => return Ok(None),
        }
    }
    let ([(top_a, na), (top_b, nb)], [(wall, cylinder)]) = (planes.as_slice(), walls.as_slice())
    else {
        return Ok(None);
    };
    // The top is the plane a concave edge of the request runs on.
    let on = |face: &Shape, edge: &Shape| -> OgeomResult<bool> {
        Ok(explore_unique(model, face, ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(edge)))
    };
    let Some(concave) = at.iter().find(|e| !e.2) else {
        return Ok(None);
    };
    let (up, out) = if on(top_a, &concave.0)? {
        (*na, *nb)
    } else if on(top_b, &concave.0)? {
        (*nb, *na)
    } else {
        return Ok(None);
    };
    let axis = cylinder.axis();
    let radius = cylinder.radius();
    let square = tol.angular() * 10.0;
    if up.dot(out).abs() > square || axis.direction.vector().cross(up).magnitude() > square {
        return Ok(None);
    }
    // The side tangent to the wall at the vertex: the axis `radius` behind
    // it, square to it, through the vertex.
    let behind = corner - out * radius;
    let off = (behind - axis.location) - up * (behind - axis.location).dot(up);
    if off.magnitude() > tol.confusion() * 100.0 {
        return Ok(None);
    }
    let _ = wall;
    // Each convex edge a line on the side, leaving along with a concave arc.
    let mut sides: Vec<Vector> = Vec::new();
    let mut paired = vec![false; at.len()];
    for (i, (edge, _, convex)) in at.iter().enumerate() {
        if !convex {
            continue;
        }
        let Some((curve, _)) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            return Ok(None);
        };
        if !matches!(curve, Curve::Line(_)) {
            return Ok(None);
        }
        let Some(t) = leaving(model, edge, corner, tol)? else {
            return Ok(None);
        };
        let mut found = false;
        for (j, (other, _, other_convex)) in at.iter().enumerate() {
            if *other_convex || paired[j] {
                continue;
            }
            if leaving(model, other, corner, tol)?.is_some_and(|u| u.dot(t) > 1.0 - square) {
                paired[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            return Ok(None);
        }
        paired[i] = true;
        sides.push(t);
    }
    if sides.is_empty() || paired.iter().any(|p| !p) {
        return Ok(None);
    }
    Ok(Some(CuspCorner {
        vertex: vertex.clone(),
        corner,
        up,
        out,
        radius,
        sides,
        blend,
    }))
}

impl CuspCorner {
    /// How wide the top stands between the side and the wall a distance
    /// `s` along the side from the tangent point.
    fn width(&self, s: f64) -> f64 {
        self.radius - (self.radius * self.radius - s * s).max(0.0).sqrt()
    }

    /// The S's turn at `s`: the angle each of its arcs turns through, a
    /// right angle where the top holds both bands.
    fn turn(&self, s: f64) -> f64 {
        (1.0 - self.width(s) / (2.0 * self.blend))
            .clamp(0.0, 1.0)
            .acos()
    }

    /// Where the S leaves the side at `s` along `side`: its convex arc's
    /// foot, below the top by the arc's rise.
    fn foot(&self, side: Vector, s: f64) -> Point {
        self.corner + side * s - self.up * (self.blend * self.turn(s).sin())
    }
}

/// The S's section through a top `w = 2 r (1 - cos theta)` wide, as a
/// distance in from the side and a rise from its foot, at `u` along it
/// (0 at its foot on the side, 1 at its contact with the wall). Its tangent
/// turns from upright by `theta (1 - cos 2 pi u) / 2`, a raised cosine of
/// arc length: the curvature is continuous, nought at both ends and at the
/// inflection, so a section across it asks no kink of a patch meeting it.
/// It is point-symmetric about its middle, `w` wide and `2 r sin theta`
/// high, as the two arcs of radius `r` turning `theta` are, its length
/// `2 r sin theta / (cos(theta/2) J0(theta/2))`. Read by Gauss-Legendre
/// quadrature, eight points on each of eight panels.
fn smooth_s(theta: f64, r: f64, u: f64) -> (f64, f64) {
    const NODES: [(f64, f64); 4] = [
        (0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    ];
    let turn = |t: f64| theta * (1.0 - (std::f64::consts::TAU * t).cos()) / 2.0;
    let integral = |to: f64| -> (f64, f64) {
        let panels = 8;
        let step = to / f64::from(panels);
        let (mut across, mut up) = (0.0, 0.0);
        for k in 0..panels {
            let mid = step * (f64::from(k) + 0.5);
            for (x, weight) in NODES {
                for t in [mid - 0.5 * step * x, mid + 0.5 * step * x] {
                    across += weight * turn(t).sin();
                    up += weight * turn(t).cos();
                }
            }
        }
        (across * 0.5 * step, up * 0.5 * step)
    };
    let height = 2.0 * r * theta.sin();
    let (_, whole) = integral(1.0);
    let length = height / whole;
    let (across, up) = integral(u);
    (across * length, up * length)
}

/// Close the cusp corner `c` of `banded`, the request's bands built on
/// `solid`. The convex band is stopped on a plane square to its edge, the
/// concave band on a plane through the wall's axis, each past where their
/// rails cross on the top, and the cell short of those planes round the
/// tangent point cut away. Put back: the top beyond a conic from the one
/// band's rail end to the other's, the side below its trims and the wall
/// above them, and on each side of the tangent point one patch G1 to both
/// band ends, the top, the side and the wall.
#[allow(clippy::too_many_lines, reason = "one corner, assembled end to end")]
pub(crate) fn close_cusp(
    model: &mut Model,
    solid: &Shape,
    banded: Built,
    c: &CuspCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let r = c.blend;
    let big = c.radius;
    if 4.0 * r >= big {
        ogeom_bail!(
            Construction,
            "a cusp where a flat side runs tangent into a wall of radius {big} closes \
             with bands of radius below a quarter of it, not {r}; see \
             docs/CAPABILITIES.md, fillet.edge-blends"
        );
    }
    let [along, ..] = c.sides.as_slice() else {
        ogeom_bail!(Construction, "a cusp corner has no side");
    };
    let along = *along;
    // The signs of the sides rounded along `along`: both ways from the
    // tangent point, or one way only, the cell then ending on the plane
    // through the tangent point and the wall's axis, square to the side.
    let signs: &[f64] = match c.sides.len() {
        1 => &[1.0],
        2 if c.sides[1].dot(along) < -0.5 => &[1.0, -1.0],
        _ => ogeom_bail!(
            Construction,
            "a cusp corner's sides do not run apart along one line"
        ),
    };
    let both = signs.len() == 2;
    // The rails cross on the top `2 sqrt(R r)` along the side from the
    // tangent point. The convex band stops past that on a plane square to
    // the side, the concave band on a plane through the axis at an angle
    // its whole section stands inside of, and the convex band's section
    // stands inside the concave one's.
    let cross = 2.0 * (big * r).sqrt();
    let stop = cross + 0.7 * (big - cross);
    let lowest = (stop / (big - r)).atan();
    let highest = (stop / (big + r)).asin();
    if lowest >= highest || lowest.is_nan() || highest.is_nan() {
        ogeom_bail!(
            Construction,
            "a cusp corner's bands leave no room to stop them apart"
        );
    }
    let angle = 0.5 * (lowest + highest);
    let axis = c.corner - c.out * big;
    let radial = |phi: f64| c.out * phi.cos() + along * phi.sin();
    let wedge = |sign: f64| -> (Vector, Point) {
        // Square to the plane through the axis at `sign * angle`, pointing
        // towards the tangent point.
        let u = radial(sign * angle);
        let n = c.up.cross(u);
        let n = if n.dot(c.corner - axis) >= 0.0 { n } else { -n };
        (n, axis)
    };
    let mut walls: Vec<(Vector, Point)> = vec![(-along, c.corner + along * stop), wedge(1.0)];
    if both {
        walls.extend([(along, c.corner - along * stop), wedge(-1.0)]);
    } else {
        walls.push((along, c.corner));
    }
    walls.extend([
        (-c.out, c.corner + c.out * r),
        (c.up, c.corner - c.up * 2.0 * r),
        (-c.up, c.corner + c.up * 2.0 * r),
    ]);
    let cell = crate::corner::convex_block(model, &walls, tol)?;
    let cut = ogeom_bool::cut(model, &banded.shape, &cell, tol)?;
    let planes = [(c.up, c.corner), (c.out, c.corner)];
    crate::setback::guard_cell(
        model,
        solid,
        &banded.shape,
        &cut.shape,
        &c.vertex,
        &planes,
        &walls,
        tol,
    )?;
    let history = banded.history.then(&cut.history);
    let slack = tol.confusion() * 100.0;
    let on_plane = |p: Point, (n, at): (Vector, Point)| n.dot(p - at).abs() <= slack;
    let in_cell = |p: Point| walls.iter().all(|(n, at)| n.dot(p - *at) >= -slack);
    // What the cut leaves, less the faces it leaves on the cell's walls.
    let mut kept: Vec<Shape> = Vec::new();
    for face in explore_unique(model, &cut.shape, ShapeType::Face)? {
        let on_wall = matches!(placed(model, &face, tol)?, Some(SurfaceGeometry::Plane(_)))
            && {
                let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
                walls.iter().any(|w| {
                    outward.cross(w.0).magnitude() <= tol.angular() * 10.0 && on_plane(origin, *w)
                })
            }
            && explore_unique(model, &face, ShapeType::Vertex)?
                .iter()
                .all(|v| point_of(model, v).is_ok_and(in_cell));
        if !on_wall {
            kept.push(face);
        }
    }
    // A band's end on a wall: the edge of a kept face of the kind asked
    // lying on that wall inside the cell, with its face.
    let end_on = |model: &Model,
                  wall: (Vector, Point),
                  kind: &dyn Fn(&SurfaceGeometry) -> bool|
     -> OgeomResult<Option<(Shape, Shape)>> {
        for face in &kept {
            let Some(surface) = placed(model, face, tol)? else {
                continue;
            };
            if !kind(&surface) {
                continue;
            }
            for e in explore_unique(model, face, ShapeType::Edge)? {
                let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                    continue;
                };
                let lies = [lo, f64::midpoint(lo, hi), hi]
                    .iter()
                    .map(|t| curve.point_at(*t, tol))
                    .collect::<OgeomResult<Vec<Point>>>()?
                    .into_iter()
                    .all(|p| on_plane(p, wall) && in_cell(p));
                if lies && !matches!(curve, Curve::Line(_)) {
                    return Ok(Some((e, face.clone())));
                }
            }
        }
        Ok(None)
    };
    let convex_kind = |s: &SurfaceGeometry| match s {
        SurfaceGeometry::Cylinder(cy) => {
            cy.cylinder()
                .axis()
                .direction
                .vector()
                .cross(along)
                .magnitude()
                <= tol.angular() * 10.0
        }
        _ => false,
    };
    let concave_kind = |s: &SurfaceGeometry| {
        !matches!(s, SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_))
    };
    let top_face = {
        let mut found = None;
        for face in &kept {
            if matches!(placed(model, face, tol)?, Some(SurfaceGeometry::Plane(_))) {
                let (origin, n) = ogeom_algo::face_normal(model, face, tol)?;
                if n.dot(c.up) > 1.0 - tol.angular() * 10.0 && on_plane(origin, (c.up, c.corner)) {
                    found = Some(face.clone());
                    break;
                }
            }
        }
        found
    };
    let Some(_top_face) = top_face else {
        ogeom_bail!(Construction, "the cusp corner's top is not kept");
    };

    let fit = crate::marched::band_fit_target(tol);
    let ts: Vec<f64> = (0..=256).map(|i| f64::from(i) / 256.0).collect();
    let tip = ogeom_algo::make_vertex(model, c.corner).shape;
    let segment = |model: &mut Model,
                   (va, a): (&Shape, Point),
                   (vb, b): (&Shape, Point)|
     -> OgeomResult<Shape> {
        crate::support::segment_between(model, (va, a), (vb, b), tol)
    };
    let end_vertex = |model: &Model, edge: &Shape, near: Point| -> OgeomResult<(Shape, Point)> {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            ogeom_bail!(Construction, "a band's end has no vertices");
        };
        let (pa, pb) = (point_of(model, &a)?, point_of(model, &b)?);
        Ok(if pa.distance(near) <= pb.distance(near) {
            (a, pa)
        } else {
            (b, pb)
        })
    };

    struct Side {
        convex: (Shape, Shape),
        concave: (Shape, Shape),
        trim: Shape,
        top: Shape,
        /// The trim on the side, in three pieces out from the tangent
        /// point: to the strip, along it, past it to the convex rail.
        low: [Shape; 3],
        /// The trim on the wall, in three pieces likewise.
        high: [Shape; 3],
        low_at: (Shape, Point),
        high_at: (Shape, Point),
        /// The S's strip and its two section edges, nearer and farther.
        strip: Shape,
        near: Shape,
        far: Shape,
    }
    // The S's own strip runs between these shares of its reach: patches
    // close it to the tangent point and on to the bands.
    let full = 2.0 * (r * (big - r)).sqrt();
    let (s1, s2) = (
        full * std::env::var("CUSPS1")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.35),
        full * std::env::var("CUSPS2")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0.85),
    );
    let mut sides: Vec<Side> = Vec::new();
    for &sign in signs {
        let side = along * sign;
        let slab = (-side, c.corner + side * stop);
        let fan = wedge(sign);
        let (Some(convex), Some(concave)) = (
            end_on(model, slab, &convex_kind)?,
            end_on(model, fan, &concave_kind)?,
        ) else {
            ogeom_bail!(
                Construction,
                "a band of the cusp corner does not end on its plane"
            );
        };
        // The rail ends: on the top, on the side, on the wall.
        let top_convex = c.corner + side * stop - c.out * r;
        let side_convex = c.corner + side * stop - c.up * r;
        let u = radial(sign * angle);
        let top_concave = axis + u * (big + r);
        let wall_concave = axis + u * big + c.up * r;
        let (v_tc, p_tc) = end_vertex(model, &convex.0, top_convex)?;
        let (v_sc, p_sc) = end_vertex(model, &convex.0, side_convex)?;
        let (v_tw, p_tw) = end_vertex(model, &concave.0, top_concave)?;
        let (v_ww, p_ww) = end_vertex(model, &concave.0, wall_concave)?;
        // The top's trim: a conic from the convex rail's end, along it
        // towards the tangent point, to the concave rail's, along it.
        let towards = -(c.up.cross(u));
        let towards = if towards.dot(c.corner - top_concave) >= 0.0 {
            towards
        } else {
            -towards
        };
        let (conic, apex) = crate::setback::trim_conic((p_tc, -side), (p_tw, towards), tol)?;
        let trim =
            ogeom_algo::make_edge_between(model, conic, (0.0, 1.0), &v_tc, &v_tw, tol)?.shape;
        let Some(q) = crate::corner::planes_meet(&(c.up, c.corner), &slab, &fan) else {
            ogeom_bail!(
                Construction,
                "the cusp corner's stopping planes do not meet on the top"
            );
        };
        let top = crate::setback::filler_inside(
            model,
            vec![trim.clone()],
            apex,
            (&v_tc, p_tc),
            (&v_tw, p_tw),
            q,
            (c.up, c.corner),
            tol,
        )?;
        // The S at `at` along this side, `u` from its foot on the side
        // (0) through its convex turn, its inflection (1/2) and its
        // concave turn to its contact with the wall (1).
        let s_point = |at: f64, u: f64| -> Point {
            let (inward, rise) = smooth_s(c.turn(at), r, u);
            c.corner + side * at - c.up * (r * c.turn(at).sin()) - c.out * inward + c.up * rise
        };
        let fit_curve = |f: &dyn Fn(f64) -> Point| -> OgeomResult<Curve> {
            let fitted = ogeom_geom::fit::fit_curve_sampled(
                |t| Ok(f(t)),
                &ts,
                false,
                3,
                tol.confusion() * 0.5,
                tol,
            )?;
            Ok(Curve::BSpline(fitted.curve))
        };
        let edge = |model: &mut Model, curve: Curve, a: &Shape, b: &Shape| -> OgeomResult<Shape> {
            let (t0, t1) = curve.domain();
            Ok(ogeom_algo::make_edge_between(model, curve, (t0, t1), a, b, tol)?.shape)
        };
        // The pieces' ends on the side (the S's foot) and on the wall (its
        // contact) at the strip's two sections.
        let (f1, f2) = (s_point(s1, 0.0), s_point(s2, 0.0));
        let (w1, w2) = (s_point(s1, 1.0), s_point(s2, 1.0));
        let [vf1, vf2, vw1, vw2] =
            [f1, f2, w1, w2].map(|p| ogeom_algo::make_vertex(model, p).shape);
        // On the side the trim is the S's own foot all the way out.
        let low = [
            edge(model, fit_curve(&|t| c.foot(side, s1 * t))?, &tip, &vf1)?,
            edge(
                model,
                fit_curve(&|t| c.foot(side, s1 + (s2 - s1) * t))?,
                &vf1,
                &vf2,
            )?,
            edge(
                model,
                fit_curve(&|t| c.foot(side, s2 + (stop - s2) * t))?,
                &vf2,
                &v_sc,
            )?,
        ];
        // On the wall the trim leaves the tangent point up the wall's own
        // line, where its normal does not turn, and bends onto the S's
        // contact with the wall by the strip, standing at most about
        // 0.03 r above it. Past the strip it is the S's contact, its
        // sections stretched along the wall where the wedge ends short of
        // the S's full turn, reaching the concave band's rail along it.
        let bump: f64 = std::env::var("CUSPA").ok().and_then(|v| v.parse().ok()).unwrap_or(0.3);
        let wall_at = |phi: f64, rise: f64| axis + radial(sign * phi) * big + c.up * rise;
        let (phi1, phi2) = ((s1 / big).asin(), (s2 / big).asin());
        let stretch = ((full - s2) / (big * angle.sin() - s2) - 1.0).max(0.0);
        let high = [
            edge(
                model,
                fit_curve(&|u| {
                    let phi = phi1 * u * u;
                    let rise = r * c.turn(big * phi.sin()).sin() + bump * r * u * (1.0 - u).powi(3);
                    wall_at(phi, rise)
                })?,
                &tip,
                &vw1,
            )?,
            edge(
                model,
                fit_curve(&|t| s_point(s1 + (s2 - s1) * t, 1.0))?,
                &vw1,
                &vw2,
            )?,
            edge(
                model,
                fit_curve(&|t| {
                    let phi = phi2 + (angle - phi2) * t;
                    let at = s2 + (big * phi.sin() - s2) * stretch.mul_add(t * t, 1.0);
                    wall_at(phi, r * c.turn(at).sin())
                })?,
                &vw2,
                &v_ww,
            )?,
        ];
        // The strip: the S's sections from the nearer to the farther, a
        // surface tangent to the side along its foot and to the wall along
        // its contact, as every section is.
        let fitted = ogeom_geom::fit::fit_surface_sampled(
            |u, v| Ok(s_point(s1 + (s2 - s1) * v, u)),
            &ogeom_geom::fit::Sampling::even(((0.0, 1.0), (0.0, 1.0)), 16),
            3,
            tol.confusion() * 0.5,
            tol,
        )?;
        let near = edge(model, fit_curve(&|u| s_point(s1, u))?, &vf1, &vw1)?;
        let far = edge(model, fit_curve(&|u| s_point(s2, u))?, &vf2, &vw2)?;
        // The strip's loop runs counter-clockwise about its surface's own
        // normal, and the face is turned to face out of the solid, as the
        // side does along its foot.
        let normal = {
            use ogeom_geom::Surface as _;
            fitted.curve.normal_at(0.0, 0.5, tol)?.vector()
        };
        let ring = vec![low[1].clone(), far.clone(), high[1].reversed(), near.reversed()];
        let ring = if side.cross(c.up).dot(normal) > 0.0 {
            ring
        } else {
            ring.iter().rev().map(Shape::reversed).collect()
        };
        let strip_wire = ogeom_algo::make_wire(model, &ring, tol)?.shape;
        let strip = ogeom_algo::make_face(
            model,
            SurfaceGeometry::BSpline(fitted.curve),
            &[strip_wire],
            tol,
        )?
        .shape;
        ogeom_algo::trimmed_or_fitted_where_bare(model, &strip, tol)?;
        let strip = if normal.dot(c.out) > 0.0 {
            strip
        } else {
            strip.reversed()
        };
        sides.push(Side {
            convex,
            concave,
            trim,
            top,
            low,
            high,
            low_at: (v_sc, p_sc),
            high_at: (v_ww, p_ww),
            strip,
            near,
            far,
        });
    }

    // The side put back below its trims, out to the stopping planes.
    let bottom = c.corner - c.up * 2.0 * r;
    let floor: Vec<(Shape, Point)> = [1.0, -1.0]
        .iter()
        .map(|sign| {
            let p = bottom + along * (stop * sign);
            (ogeom_algo::make_vertex(model, p).shape, p)
        })
        .collect();
    let side_ring = if both {
        let mut ring: Vec<Shape> = sides[0].low.to_vec();
        ring.extend([
            segment(
                model,
                (&sides[0].low_at.0, sides[0].low_at.1),
                (&floor[0].0, floor[0].1),
            )?,
            segment(model, (&floor[0].0, floor[0].1), (&floor[1].0, floor[1].1))?,
            segment(
                model,
                (&floor[1].0, floor[1].1),
                (&sides[1].low_at.0, sides[1].low_at.1),
            )?,
        ]);
        ring.extend(sides[1].low.iter().rev().map(Shape::reversed));
        ring
    } else {
        // Down the line the side and the wall share, from the tangent point.
        let under = ogeom_algo::make_vertex(model, bottom).shape;
        let mut ring: Vec<Shape> = sides[0].low.to_vec();
        ring.extend([
            segment(
                model,
                (&sides[0].low_at.0, sides[0].low_at.1),
                (&floor[0].0, floor[0].1),
            )?,
            segment(model, (&floor[0].0, floor[0].1), (&under, bottom))?,
            segment(model, (&under, bottom), (&tip, c.corner))?,
        ]);
        ring
    };
    let side_surface: SurfaceGeometry = ogeom_geom::PlaneSurface::over(
        Plane::through(c.corner, Direction::new(c.out, tol)?),
        (-4.0 * stop, 4.0 * stop),
        (-4.0 * stop, 4.0 * stop),
    )?
    .into();
    let side_wire = ogeom_algo::make_wire(model, &side_ring, tol)?.shape;
    let side_face = ogeom_algo::make_face(model, side_surface, &[side_wire], tol)?.shape;
    ogeom_algo::trimmed_or_fitted_where_bare(model, &side_face, tol)?;

    // The wall put back above its trims, up to the cell's top.
    let crown = axis + c.up * 2.0 * r;
    let circle = ogeom_math::Circle::new(
        ogeom_math::Frame::new(
            crown,
            Direction::new(c.up, tol)?,
            Direction::new(c.out, tol)?,
            tol,
        )?,
        big,
        tol,
    )?;
    let ring_curve = Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    let angle_of = |p: Point| -> f64 {
        let d = p - crown;
        d.dot(c.up.cross(c.out)).atan2(d.dot(c.out))
    };
    // The crown's ends: over each side's wedge, or over the tangent point
    // where the cell ends there.
    let tops: Vec<(Shape, Point)> = [angle, if both { -angle } else { 0.0 }]
        .iter()
        .map(|phi| {
            let p = crown + radial(*phi) * big;
            (ogeom_algo::make_vertex(model, p).shape, p)
        })
        .collect();
    let (ta, tb) = (angle_of(tops[1].1), angle_of(tops[0].1));
    let (t0, t1, v0, v1, minus_first) = if ta <= tb {
        (ta, tb, &tops[1].0, &tops[0].0, true)
    } else {
        (tb, ta, &tops[0].0, &tops[1].0, false)
    };
    let crown_edge = ogeom_algo::make_edge_between(model, ring_curve, (t0, t1), v0, v1, tol)?.shape;
    let crown_edge = if minus_first {
        crown_edge
    } else {
        crown_edge.reversed()
    };
    let mut wall_ring = if both {
        let mut ring: Vec<Shape> = sides[1].high.to_vec();
        ring.push(segment(
            model,
            (&sides[1].high_at.0, sides[1].high_at.1),
            (&tops[1].0, tops[1].1),
        )?);
        ring
    } else {
        // Up the wall's line through the tangent point.
        vec![segment(model, (&tip, c.corner), (&tops[1].0, tops[1].1))?]
    };
    wall_ring.extend([
        crown_edge,
        segment(
            model,
            (&tops[0].0, tops[0].1),
            (&sides[0].high_at.0, sides[0].high_at.1),
        )?,
    ]);
    wall_ring.extend(sides[0].high.iter().rev().map(Shape::reversed));
    let wall_surface = {
        let mut found = None;
        for face in &kept {
            if let Some(SurfaceGeometry::Cylinder(cy)) = placed(model, face, tol)?
                && (cy.cylinder().radius() - big).abs() <= tol.confusion() * 100.0
                && cy
                    .cylinder()
                    .axis()
                    .direction
                    .vector()
                    .cross(c.up)
                    .magnitude()
                    <= tol.angular() * 10.0
            {
                found = Some(SurfaceGeometry::Cylinder(cy));
                break;
            }
        }
        let Some(surface) = found else {
            ogeom_bail!(Construction, "the cusp corner's wall is not kept");
        };
        surface
    };
    let wall_wire = ogeom_algo::make_wire(model, &wall_ring, tol)?.shape;
    let wall_face = ogeom_algo::make_face(model, wall_surface, &[wall_wire], tol)?.shape;
    ogeom_algo::trimmed_or_fitted_where_bare(model, &wall_face, tol)?;

    // One patch on each side, held inside the solid and the wall's fill:
    // nothing past the side; above the top, nothing inside the wall and
    // nothing outside the concave band's own fill, the corner of the top
    // and the wall less the band's tube.
    let region = |p: Point| -> (f64, Vector) {
        let front = (p - c.corner).dot(c.out);
        let high = (p - c.corner).dot(c.up);
        if high <= 0.0 {
            return (front, c.out);
        }
        let d = p - axis;
        let flat = d - c.up * d.dot(c.up);
        let rho = flat.magnitude();
        let outward = if rho > 0.0 { flat / rho } else { c.out };
        let mut worst = (front, c.out);
        let mut take = |by: f64, n: Vector| {
            if by > worst.0 {
                worst = (by, n);
            }
        };
        take(big - rho, -outward);
        take(high - r, c.up);
        let across = rho - (big + r);
        if across >= 0.0 {
            take(high, c.up);
        } else if high < r {
            let reach = across.hypot(high - r);
            if reach < r && reach > 0.0 {
                take(r - reach, (-outward * across - c.up * (high - r)) / reach);
            }
        }
        worst
    };
    let mut patches: Vec<Shape> = Vec::with_capacity(4);
    let mut strips: Vec<Shape> = Vec::with_capacity(2);
    let mut tops_kept: Vec<Shape> = Vec::with_capacity(2);
    for s in &sides {
        let g1 = |edge: &Shape, support: &Shape| FillBoundary {
            edge: edge.clone(),
            support: Some(support.clone()),
            continuity: Continuity::G1,
        };
        // To the tangent point: the side's and the wall's trims there and
        // the strip's nearer section.
        let to_tip = vec![
            g1(&s.low[0], &side_face),
            g1(&s.near, &s.strip),
            g1(&s.high[0], &wall_face),
        ];
        let filled = make_filling_n_region(model, &to_tip, &[], fit, &region, tol)?;
        if std::env::var_os("DOWNWHERE").is_some() {
            eprintln!("TOTIP kind {:?} faces {}", model.kind_of(&filled.built.shape)?, explore_unique(model, &filled.built.shape, ShapeType::Face)?.len());
        }
        patches.push(filled.built.shape);
        // On to the bands: their ends, the top's trim, the trims past the
        // strip and its farther section.
        let to_bands = vec![
            g1(&s.convex.0, &s.convex.1),
            g1(&s.trim, &s.top),
            g1(&s.concave.0, &s.concave.1),
            g1(&s.high[2], &wall_face),
            g1(&s.far, &s.strip),
            g1(&s.low[2], &side_face),
        ];
        let filled = make_filling_n_region(model, &to_bands, &[], fit, &region, tol)?;
        patches.push(filled.built.shape);
        strips.push(s.strip.clone());
        tops_kept.push(s.top.clone());
    }
    let mut faces = kept;
    faces.extend(tops_kept);
    faces.push(side_face);
    faces.push(wall_face);
    faces.extend(patches.iter().cloned());
    faces.extend(strips.iter().cloned());
    let faces_given = faces.clone();
    let sew_at: f64 = std::env::var("CUSPSEW").ok().and_then(|v| v.parse().ok()).unwrap_or(fit);
    let sewn = ogeom_algo::sew_within(model, &faces, sew_at, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(
            Construction,
            "the cusp corner's faces did not close: {} shells",
            sewn.shells.len()
        );
    }
    let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
    // What the bands ran past the side with inside the cell is rebuilt.
    crate::support::settle_cut_offs_within(&walls, slack);
    if std::env::var_os("DOWNWHERE").is_some() {
        let count = |sh: &Shape| -> OgeomResult<usize> {
            let mut n = 0;
            for f in explore_unique(model, sh, ShapeType::Face)? {
                if matches!(placed(model, &f, tol)?, Some(SurfaceGeometry::BSpline(_))) {
                    n += 1;
                }
            }
            Ok(n)
        };
        let mut given = 0;
        for f in &faces_given {
            if matches!(placed(model, f, tol)?, Some(SurfaceGeometry::BSpline(_))) {
                given += 1;
            }
        }
        eprintln!("SEWCOUNT given {given} sewn {}", count(&rounded)?);
    }
    let mut history: History = history;
    for patch in patches.into_iter().chain(strips) {
        history.generate(&c.vertex, patch);
    }
    Ok(Built::new(rounded, history.then(&sewn.history)))
}
