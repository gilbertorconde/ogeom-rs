//! Setback openings between round, conic and bevel tools.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Continuity, Curve, Curve3d as _};
use ogeom_math::{Direction, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::fill_n::{FillBoundary, make_filling_n_within};
use crate::setback::{CornerEnd, guard_cell, in_cell, point_of, trim_conic};

/// The bands of `banded` meeting at `vertex`, whose corner reads as
/// `frame`, stopped at `setback` along their edges and the opening between
/// them filled by one patch at `continuity`, each band ending as `ends` says.
#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    reason = "one corner, assembled end to end"
)]
pub(crate) fn close(
    model: &mut Model,
    solid: &Shape,
    banded: Built,
    vertex: &Shape,
    frame: &crate::engine::vertex_read::CornerFrame,
    setback: &[f64],
    ends: &[CornerEnd],
    continuity: Continuity,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = frame.d.len();
    // The cell the patch replaces: in from every plane through the vertex,
    // and short of each band's section.
    let mut walls: Vec<(Vector, Point)> = frame.m.iter().map(|mk| (*mk, frame.corner)).collect();
    let sections: Vec<(Vector, Point)> = (0..n)
        .map(|k| {
            let dk = frame.d[k].vector();
            (-dk, frame.corner + dk * setback[k])
        })
        .collect();
    walls.extend(sections.iter().copied());
    let cell = crate::engine::vertex_geometry::convex_block(model, &walls, tol)?;
    let cut = ogeom_bool::cut(model, &banded.shape, &cell, tol)?;
    let planes: Vec<(Vector, Point)> = frame.m.iter().map(|mk| (*mk, frame.corner)).collect();
    guard_cell(
        model,
        solid,
        &banded.shape,
        &cut.shape,
        vertex,
        &planes,
        &walls,
        tol,
    )?;
    let history = banded.history.then(&cut.history);

    let on_plane = |p: Point, (normal, at): (Vector, Point)| -> bool {
        normal.dot(p - at).abs() <= tol.confusion() * 100.0
    };
    let section_of = |face: &Shape, model: &Model| -> OgeomResult<Option<usize>> {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            return Ok(None);
        };
        if !matches!(
            model.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        ) {
            return Ok(None);
        }
        let (origin, outward) = ogeom_algo::face_normal(model, face, tol)?;
        let Some(k) = (0..n).find(|&k| {
            outward.cross(sections[k].0).magnitude() <= tol.angular() * 10.0
                && on_plane(origin, sections[k])
        }) else {
            return Ok(None);
        };
        for v in explore_unique(model, face, ShapeType::Vertex)? {
            if !in_cell(&walls, point_of(model, &v)?, tol) {
                return Ok(None);
            }
        }
        Ok(Some(k))
    };

    // Each band's end arc, on its section plane, and the band face it ends.
    let faces = explore_unique(model, &cut.shape, ShapeType::Face)?;
    let mut kept: Vec<Shape> = Vec::with_capacity(faces.len());
    for face in &faces {
        if section_of(face, model)?.is_none() {
            kept.push(face.clone());
        }
    }
    let mut arcs: Vec<(Shape, Shape)> = Vec::with_capacity(n);
    for (k, section) in sections.iter().enumerate() {
        if ends[k] == CornerEnd::Bevel {
            arcs.push(bevel_end(model, &kept, frame, *section, &walls, tol)?);
            continue;
        }
        arcs.push(round_end(model, &kept, k, *section, &walls, tol)?);
    }

    // On host plane `ring[k + 1]`, between band `k` and band `k + 1`: each
    // band's arc end there, and where the two section planes meet it.
    let host = |t: usize| -> (Vector, Point) { (frame.m[frame.ring[t]], frame.corner) };
    let end_on =
        |model: &Model, arc: &Shape, plane: (Vector, Point)| -> OgeomResult<(Shape, Point)> {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, arc)? else {
                ogeom_bail!(Construction, "a band's end arc has no vertices");
            };
            for v in [a, b] {
                let p = point_of(model, &v)?;
                if on_plane(p, plane) {
                    return Ok((v, p));
                }
            }
            ogeom_bail!(Construction, "a band's end arc does not reach its face");
        };
    let fit = crate::marched::band_fit_target(tol);
    let mut boundary: Vec<FillBoundary> = Vec::with_capacity(2 * n);
    let mut fillers: Vec<Shape> = Vec::with_capacity(n);
    for k in 0..n {
        let next = (k + 1) % n;
        let plane = host(next);
        let (v0, p0) = end_on(model, &arcs[k].0, plane)?;
        let (v1, p1) = end_on(model, &arcs[next].0, plane)?;
        let Some(q) =
            crate::engine::vertex_geometry::planes_meet(&plane, &sections[k], &sections[next])
        else {
            ogeom_bail!(
                Construction,
                "two section planes do not meet on the face between them"
            );
        };
        // The rails run along the edges: each conic leaves its band's rail
        // along that edge's direction, towards the vertex.
        let (t0, t1) = (-frame.d[k].vector(), -frame.d[next].vector());
        let (conic, apex) = trim_conic((p0, t0), (p1, t1), tol)?;
        let (trim, filler) =
            trim_and_filler(model, conic, apex, (&v0, p0), (&v1, p1), q, plane, tol)?;
        boundary.push(FillBoundary {
            edge: arcs[k].0.clone(),
            support: Some(arcs[k].1.clone()),
            continuity: match ends[k] {
                CornerEnd::Round => continuity,
                CornerEnd::Bevel => Continuity::C0,
            },
        });
        // A trim from a bevel's end meets its face in position only: the
        // patch runs along the bevel's end there, which leaves the face at
        // the bevel's angle, so it cannot take the face's tangent plane.
        // Between two rounds both ends lie tangent to the face.
        boundary.push(FillBoundary {
            edge: trim,
            support: Some(filler.clone()),
            continuity: if ends[k] == CornerEnd::Round && ends[next] == CornerEnd::Round {
                continuity
            } else {
                Continuity::C0
            },
        });
        fillers.push(filler);
    }

    // The vertex is convex, so what of the patch stands past one of its
    // planes lies outside the solid: the patch is held within them.
    let within: Vec<(Vector, Point)> = frame.m.iter().map(|m| (-*m, frame.corner)).collect();
    let filled = make_filling_n_within(model, &boundary, &[], fit, &within, tol)?;
    let patch = filled.built.shape.clone();
    let mut faces = kept;
    faces.extend(fillers);
    faces.push(patch.clone());
    let sewn = ogeom_algo::sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the setback corner's faces did not close");
    }
    let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
    // Sewing may rebuild a face whose edges it joins: an earlier corner's
    // patch, as well as this one, is followed to the face it became.
    let mut history: History = history.then(&sewn.history);
    for image in sewn.history.trace(&patch) {
        history.generate(vertex, image.clone());
    }
    Ok(Built::new(rounded, history))
}

/// The edge on `conic` from `v0` to `v1`, and the face's filler inside the
/// corner's cell: the planar face on `plane` bounded by the trim and the
/// two section lines from its ends to `q`, where they meet.
#[allow(clippy::too_many_arguments, reason = "one face's trim and its filler")]
pub(crate) fn trim_and_filler(
    model: &mut Model,
    conic: Curve,
    apex: Point,
    (v0, p0): (&Shape, Point),
    (v1, p1): (&Shape, Point),
    q: Point,
    plane: (Vector, Point),
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let trim = ogeom_algo::make_edge_between(model, conic, (0.0, 1.0), v0, v1, tol)?.shape;
    let filler = filler_inside(
        model,
        vec![trim.clone()],
        apex,
        (v0, p0),
        (v1, p1),
        q,
        plane,
        tol,
    )?;
    Ok((trim, filler))
}

/// The face's filler inside the corner's cell: the planar face on `plane`
/// bounded by `trim`, the edges in order from `v0` to `v1` along a trim
/// whose middle control point is `apex`, and the two section lines from
/// its ends to `q`, where they meet.
#[allow(clippy::too_many_arguments, reason = "one face's filler")]
pub(crate) fn filler_inside(
    model: &mut Model,
    trim: Vec<Shape>,
    apex: Point,
    (v0, p0): (&Shape, Point),
    (v1, p1): (&Shape, Point),
    q: Point,
    plane: (Vector, Point),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let vq = ogeom_algo::make_vertex(model, q).shape;
    let mut ring = trim;
    ring.push(crate::support::segment_between(
        model,
        (v1, p1),
        (&vq, q),
        tol,
    )?);
    ring.push(crate::support::segment_between(
        model,
        (&vq, q),
        (v0, p0),
        tol,
    )?);
    let outward = -plane.0;
    let reach = (p0.distance(q) + p1.distance(q) + apex.distance(q)) * 2.0 + 1.0;
    let surface: ogeom_geom::SurfaceGeometry = ogeom_geom::PlaneSurface::over(
        Plane::through(q, Direction::new(outward, tol)?),
        (-reach, reach),
        (-reach, reach),
    )?
    .into();
    let turn = (p0 - q).cross(p1 - q);
    let ring: Vec<Shape> = if turn.dot(outward) > 0.0 {
        ring
    } else {
        ring.iter().rev().map(Shape::reversed).collect()
    };
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &[ring], tol)?.shape)
}

/// Band `k`'s end arc on its section plane `section`, and the band face it
/// ends: the curved edge of a curved face of `kept` lying on the plane
/// inside the cell `walls` bound.
pub(crate) fn round_end(
    model: &Model,
    kept: &[Shape],
    k: usize,
    section: (Vector, Point),
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let mut found: Option<(Shape, Shape)> = None;
    for face in kept {
        let is_plane = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        if is_plane {
            continue;
        }
        for e in explore_unique(model, face, ShapeType::Edge)? {
            let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                continue;
            };
            if matches!(c, Curve::Line(_)) {
                continue;
            }
            let on = [lo, f64::midpoint(lo, hi), hi]
                .iter()
                .map(|t| c.point_at(*t, tol))
                .collect::<OgeomResult<Vec<Point>>>()?
                .into_iter()
                .all(|p| {
                    section.0.dot(p - section.1).abs() <= tol.confusion() * 100.0
                        && in_cell(walls, p, tol)
                });
            if on {
                found = Some((e, face.clone()));
            }
        }
    }
    let Some(found) = found else {
        ogeom_bail!(
            Construction,
            "band {k} does not end in an arc on its section plane"
        );
    };
    Ok(found)
}

/// A bevel's straight end on its section plane `section`, and the bevel's
/// face: the planar face of `kept`, on none of the corner's host planes,
/// with a straight edge lying on the section plane inside the cell `walls`
/// bound.
pub(crate) fn bevel_end(
    model: &Model,
    kept: &[Shape],
    frame: &crate::engine::vertex_read::CornerFrame,
    section: (Vector, Point),
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let on_plane = |p: Point| {
        section.0.dot(p - section.1).abs() <= tol.confusion() * 100.0 && in_cell(walls, p, tol)
    };
    for face in kept {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        if !matches!(
            model.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        ) {
            continue;
        }
        let (_, outward) = ogeom_algo::face_normal(model, face, tol)?;
        let host = frame
            .m
            .iter()
            .any(|mk| mk.cross(outward).magnitude() <= tol.angular() * 10.0);
        if host || outward.cross(section.0).magnitude() <= tol.angular() * 10.0 {
            continue;
        }
        for e in explore_unique(model, face, ShapeType::Edge)? {
            let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                continue;
            };
            if !matches!(c, Curve::Line(_)) {
                continue;
            }
            let on = [lo, hi]
                .iter()
                .map(|t| c.point_at(*t, tol))
                .collect::<OgeomResult<Vec<Point>>>()?
                .into_iter()
                .all(on_plane);
            if on {
                return Ok((e, face.clone()));
            }
        }
    }
    ogeom_bail!(
        Construction,
        "a bevel does not end in a straight edge on its section plane"
    );
}
