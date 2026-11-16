//! Vertex setback requests and the point opening beside an edge left sharp.

use ogeom_algo::Built;
use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Reversible as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::classify::faces_round;
use super::setbacks::{FaceTrim, SetbackCorner, SetbackEnd};
use super::vertex_read::CornerFrame;

/// Refuse foreign faces intersecting the setback region before P2 trims
/// its hosts. This reads geometry only: no preview boolean may mutate the
/// input or appear in an otherwise sheet-built feature's history.
pub(crate) fn guard(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    frame: &CornerFrame,
    setback: &[f64],
    tol: Tolerances,
) -> OgeomResult<()> {
    let near = tol.confusion() * 100.0;
    let mut walls: Vec<_> = frame.m.iter().map(|m| (*m, frame.corner)).collect();
    walls.extend(
        frame
            .d
            .iter()
            .zip(setback)
            .map(|(d, s)| (-d.vector(), frame.corner + d.vector() * *s)),
    );
    let inside = |p: Point| walls.iter().all(|(n, at)| n.dot(p - *at) >= -near);
    let refused = || {
        ogeom_core::ogeom_err!(
            Construction,
            "another feature of the solid lies inside the setback corner's cell"
        )
    };
    let extent = setback.iter().copied().fold(1.0, f64::max) * 2.0;
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let vertices = explore_unique(model, &face, ShapeType::Vertex)?;
        if vertices.iter().any(|v| v.is_same(vertex)) {
            continue;
        }
        let planar = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
        let (at, normal) = ogeom_algo::face_normal(model, &face, tol)?;
        if planar
            && frame.m.iter().any(|m| {
                normal.cross(*m).magnitude() <= tol.angular() * 10.0
                    && m.dot(at - frame.corner).abs() <= near
            })
        {
            continue;
        }
        for v in vertices {
            if inside(crate::setback::point_of(model, &v)?) {
                return Err(refused());
            }
        }
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &edge, tol)? else {
                continue;
            };
            if matches!(curve, ogeom_geom::Curve::Line(_)) {
                let a = curve.point_at(lo, tol)?;
                let b = curve.point_at(hi, tol)?;
                let mut range = (0.0_f64, 1.0_f64);
                for (n, at) in &walls {
                    let (da, db) = (n.dot(a - *at) + near, n.dot(b - *at) + near);
                    if da < 0.0 && db < 0.0 {
                        range = (1.0, 0.0);
                        break;
                    }
                    if da < 0.0 {
                        range.0 = range.0.max(da / (da - db));
                    }
                    if db < 0.0 {
                        range.1 = range.1.min(da / (da - db));
                    }
                }
                if range.0 <= range.1 {
                    return Err(refused());
                }
            } else {
                for (n, at) in &walls {
                    let plane =
                        ogeom_math::Plane::through(*at, ogeom_math::Direction::new(*n, tol)?);
                    let surface = ogeom_geom::PlaneSurface::over(
                        plane,
                        (-extent, extent),
                        (-extent, extent),
                    )?
                    .into();
                    let hit = ogeom_intersect::intersect_curve_surface(
                        &curve,
                        &surface,
                        ogeom_intersect::CurveSurfaceOptions {
                            gap: near,
                            ..Default::default()
                        },
                        tol,
                    )?;
                    if hit
                        .crossings
                        .iter()
                        .any(|p| p.on_curve >= lo && p.on_curve <= hi && inside(p.point))
                    {
                        return Err(refused());
                    }
                }
            }
        }
        // A planar foreign face can cover the cell's section without any
        // of its own boundary entering it. Test the section's vertices
        // against the face's actual trim, including its holes.
        if planar {
            for i in 0..walls.len() {
                for j in i + 1..walls.len() {
                    let normals = [normal, walls[i].0, walls[j].0];
                    let offsets = [
                        normal.dot(at - frame.corner),
                        walls[i].0.dot(walls[i].1 - frame.corner),
                        walls[j].0.dot(walls[j].1 - frame.corner),
                    ];
                    let Some(p) =
                        super::vertex_geometry::solve3(&normals, offsets).map(|p| frame.corner + p)
                    else {
                        continue;
                    };
                    if inside(p)
                        && ogeom_algo::classify_on_face(
                            model,
                            &face,
                            p,
                            ogeom_mesh::Deflection::default(),
                            tol,
                        )?
                        .is_inside_or_on()
                    {
                        return Err(refused());
                    }
                }
            }
        }
    }
    Ok(())
}

/// Refused by name unless the cell cut from `banded` at `vertex` holds
/// nothing of `solid`, the solid the corner was asked on, but the corner
/// itself. What is foreign is a face of `solid` that does not meet the
/// vertex and lies on none of the corner's own planes through it
/// (`planes`, each a normal and a point): the bands, the bevels and the
/// ball's own corner are made by the rounding, on surfaces of their own.
/// Refused where the cut takes a face on a foreign surface whole (every
/// corner of it inside the cell, `walls`, each an inward normal and a
/// point), where a face on a foreign surface that the cut leaves has an
/// edge on a wall inside the cell (the cut passed through it), or where the
/// cut leaves fewer shells. Any of these (a void, a hole, another feature)
/// would be left out of what is sewn back, and the solid would close
/// without it.
#[allow(
    clippy::too_many_arguments,
    reason = "one cell and what it is read against"
)]
pub(crate) fn guard_cell(
    model: &Model,
    solid: &Shape,
    banded: &Shape,
    cut: &Shape,
    vertex: &Shape,
    planes: &[(Vector, Point)],
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<()> {
    let near = tol.confusion() * 100.0;
    let surface_of = |face: &Shape| {
        model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
    };
    let on_a_plane = |face: &Shape| -> OgeomResult<bool> {
        let planar = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        if !planar {
            return Ok(false);
        }
        let (origin, outward) = ogeom_algo::face_normal(model, face, tol)?;
        Ok(planes.iter().any(|(normal, at)| {
            outward.cross(*normal).magnitude() <= tol.angular() * 10.0
                && normal.dot(origin - *at).abs() <= near
        }))
    };
    let mut foreign = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let meets = explore_unique(model, &face, ShapeType::Vertex)?
            .iter()
            .any(|v| v.is_same(vertex));
        if !meets
            && !on_a_plane(&face)?
            && let Some(surface) = surface_of(&face)
        {
            foreign.push(surface);
        }
    }
    let is_foreign = |face: &Shape| surface_of(face).is_some_and(|s| foreign.contains(&s));
    let inside = |p: Point| walls.iter().all(|(n, at)| n.dot(p - *at) >= -near);
    let refused = || {
        ogeom_core::ogeom_err!(
            Construction,
            "another feature of the solid lies inside the setback corner's cell (a face \
             that is neither one of the corner's faces nor made by its bands, or a shell \
             of its own), and cutting the cell would remove it; see \
             docs/CAPABILITIES.md, fillet.setback-corners"
        )
    };
    // A face the cut takes whole.
    for face in explore_unique(model, banded, ShapeType::Face)? {
        if !is_foreign(&face) {
            continue;
        }
        let mut all_in = true;
        for v in explore_unique(model, &face, ShapeType::Vertex)? {
            if !inside(crate::setback::point_of(model, &v)?) {
                all_in = false;
                break;
            }
        }
        if all_in {
            return Err(refused());
        }
    }
    // A face the cut passes through.
    for face in explore_unique(model, cut, ShapeType::Face)? {
        if !is_foreign(&face) {
            continue;
        }
        for e in explore_unique(model, &face, ShapeType::Edge)? {
            let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                continue;
            };
            let mid = curve.point_at(f64::midpoint(lo, hi), tol)?;
            if walls.iter().any(|(n, at)| n.dot(mid - *at).abs() <= near) && inside(mid) {
                return Err(refused());
            }
        }
    }
    let shells = |shape: &Shape| explore_unique(model, shape, ShapeType::Shell).map(|s| s.len());
    if shells(cut)? < shells(banded)? {
        return Err(refused());
    }
    Ok(())
}

/// P0/P1/P2 for a vertex with two circular strips and an edge of zero
/// width. P3 receives the three face trims, two end arcs and their supports.
#[allow(clippy::too_many_arguments, reason = "the validated vertex request")]
pub(crate) fn sharp(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    edges: &[Shape],
    setback: &[f64],
    sharp: usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let baked = ogeom_algo::baked_where_placed(model, solid, tol)?;
    if !baked.shape.is_same(solid) {
        let trace = |part: &Shape| -> OgeomResult<Shape> {
            let [only] = baked.history.trace(part) else {
                ogeom_bail!(Construction, "a placed setback request has no unique image");
            };
            Ok(only.clone())
        };
        let vertex = trace(vertex)?;
        let edges = edges.iter().map(trace).collect::<OgeomResult<Vec<_>>>()?;
        let built = self::sharp(
            model,
            &baked.shape,
            &vertex,
            &edges,
            setback,
            sharp,
            radius,
            tol,
        )?;
        return Ok(Built {
            shape: built.shape,
            history: baked.history.then(&built.history),
        });
    }
    let Some(planar) = super::classify::read_planar(model, solid, tol)? else {
        ogeom_bail!(Construction, "a setback corner needs planar host faces");
    };
    let Some(seats) = planar.seats() else {
        ogeom_bail!(
            Construction,
            "a setback corner's edges do not have two faces"
        );
    };
    let rounds: Vec<_> = edges
        .iter()
        .enumerate()
        .filter(|(k, _)| *k != sharp)
        .map(|(_, e)| (e.clone(), radius))
        .collect();
    let request = super::Request::from_rounds(&rounds, tol)?;
    let sections: FastMap<_, _> = request
        .blends
        .iter()
        .map(|blend| (super::id(&blend.edge), &blend.section))
        .collect();
    let mut distances = FastMap::default();
    for seat in seats.values() {
        if let Some(section) = sections.get(&seat.edge) {
            let d = super::sections::law(section).setbacks(seat, &planar)?;
            for (k, f) in seat.faces.iter().enumerate() {
                distances.insert((super::id(&planar.faces[*f].shape), seat.edge), d[k]);
            }
        }
    }
    let mut steps = FastMap::default();
    let Some(corners) = super::setbacks::setback_corners(&planar, &distances, &mut steps, tol)?
    else {
        ogeom_bail!(Construction, "the setback corner's contact curves cross");
    };
    let v = super::id(vertex);
    let Some(order) = faces_round(&planar.faces, v) else {
        ogeom_bail!(
            Construction,
            "a setback corner's faces do not close round it"
        );
    };
    let at = planar.points[&v];
    let mut ends = Vec::new();
    for k in 0..order.len() {
        let (f, g) = (order[k], order[(k + 1) % order.len()]);
        let Some(seat) = seats.values().find(|s| {
            (s.ends.0 == v || s.ends.1 == v) && s.faces.contains(&f) && s.faces.contains(&g)
        }) else {
            ogeom_bail!(Construction, "a setback corner's faces meet along no edge");
        };
        let Some(index) = edges.iter().position(|e| super::id(e) == seat.edge) else {
            ogeom_bail!(
                Construction,
                "a setback corner names no edge between its faces"
            );
        };
        let far = if seat.ends.0 == v {
            seat.ends.1
        } else {
            seat.ends.0
        };
        let along = (planar.points[&far] - at).normalized(tol)?;
        let on = at + along * setback[index];
        let (r, d) = if index == sharp {
            (0.0, [0.0; 2])
        } else {
            (
                radius,
                super::sections::law(sections[&seat.edge]).setbacks(seat, &planar)?,
            )
        };
        let mut rails = [(0, Point::ORIGIN); 2];
        for j in 0..2 {
            let fi = seat.faces[j];
            let mut inward = planar.faces[fi].normal.cross(along);
            if inward.dot(seat.normals[1 - j]) > 0.0 {
                inward = -inward;
            }
            rails[j] = (fi, on + inward.normalized(tol)? * d[j]);
        }
        let centre = rails[0].1 - seat.normals[0] * r;
        if rails[0].0 != f {
            rails.swap(0, 1);
        }
        ends.push(SetbackEnd {
            edge: seat.edge,
            at: on,
            along,
            centre,
            radius: r,
            rails,
        });
    }
    let mut trims = Vec::new();
    for k in 0..ends.len() {
        let next = (k + 1) % ends.len();
        let face = order[next];
        let (from, to) = (ends[k].rails[1].1, ends[next].rails[0].1);
        let tangent = |i: usize, p: Point, j: usize, q: Point| -> OgeomResult<Vector> {
            if ends[i].radius > 0.0 {
                return Ok(-ends[i].along);
            }
            let distance = (ends[j].at - at).magnitude();
            (q - ends[j].along * (distance / 2.0) - p).normalized(tol)
        };
        let (conic, apex) = crate::setback::trim_conic(
            (from, tangent(k, from, next, to)?),
            (to, tangent(next, to, k, from)?),
            tol,
        )?;
        let normal = planar.faces[face].normal;
        trims.push(FaceTrim {
            face,
            plane: (normal, normal.dot(from.to_vector())),
            from,
            to,
            conic,
            apex,
        });
    }
    let corner = SetbackCorner {
        vertex: v,
        ends,
        trims,
    };
    super::rounds::assemble(
        model,
        &planar,
        &seats,
        &sections,
        &corners,
        &steps,
        &[corner],
        tol,
    )?
    .ok_or_else(|| {
        ogeom_core::ogeom_err!(
            Construction,
            "the blend engine could not close the sharp setback corner"
        )
    })
}

/// The point patch and the remainder of a five-sided opening. Its trims
/// borrow the sheet's host surfaces; only the two trims meeting at the
/// sharp point are split for the point patch.
pub(super) fn opening(
    model: &mut Model,
    sides: &[crate::round::transfinite::HoleSide],
    supports: &[Shape],
    within: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let plane = |i: usize| matches!(sides[i].support, SurfaceGeometry::Plane(_));
    let Some(start) =
        (0..sides.len()).find(|i| plane(*i) && plane((*i + sides.len() - 1) % sides.len()))
    else {
        ogeom_bail!(
            Construction,
            "a sharp corner opening has no point between two face trims"
        );
    };
    let mut trims = Vec::new();
    let mut arcs = Vec::new();
    for k in 0..5 {
        let i = (start + k) % sides.len();
        let side = &sides[i];
        if k % 2 == 1 {
            arcs.push((side.edge.clone(), supports[i].clone()));
            continue;
        }
        let curve = if side.forward {
            side.curve.clone()
        } else {
            side.curve.reversed()
        };
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &side.edge)? else {
            ogeom_bail!(Construction, "a face trim has no endpoints");
        };
        let start = curve.point_at(0.0, tol)?;
        let (a, b) =
            if crate::setback::point_of(model, &a)?.distance(start) <= tol.confusion() * 100.0 {
                (a, b)
            } else {
                (b, a)
            };
        let edges = super::point_corner::trim_edges(model, &curve, (&a, &b), k != 2, tol)?;
        let SurfaceGeometry::Plane(plane) = &side.support else {
            ogeom_bail!(Construction, "a setback trim needs a planar support");
        };
        let from = crate::setback::point_of(model, &a)?;
        let to = crate::setback::point_of(model, &b)?;
        let chord = ogeom_geom::LineCurve::segment(to, from, tol)?;
        let close = ogeom_algo::make_edge_between(
            model,
            chord.into(),
            (0.0, from.distance(to)),
            &b,
            &a,
            tol,
        )?
        .shape;
        let mut boundary = edges.clone();
        boundary.push(close);
        let wire = ogeom_algo::make_wire(model, &boundary, tol)?.shape;
        let filler = ogeom_algo::make_face(model, side.support.clone(), &[wire], tol)?.shape;
        trims.push(super::point_corner::Trim {
            curve,
            edges,
            normal: plane.plane().frame().z().vector(),
            filler,
        });
    }
    let held = super::point_corner::Held {
        within,
        seam: &[],
        region: None,
    };
    let (point, rest) = super::point_corner::close_through_q(
        model,
        [&trims[0], &trims[1], &trims[2]],
        [&arcs[0], &arcs[1]],
        [true; 5],
        &held,
        tol,
    )?;
    Ok(vec![point, rest])
}

/// A vertex API's single face over its setback opening, tangent to every
/// strip and host face and held within the corner's material half-spaces.
pub(super) fn filled_opening(
    model: &mut Model,
    sides: &[crate::round::transfinite::HoleSide],
    supports: &[Shape],
    within: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let boundary: Vec<_> = sides
        .iter()
        .zip(supports)
        .map(|(side, support)| crate::fill_n::FillBoundary {
            edge: side.edge.clone(),
            support: Some(support.clone()),
            continuity: ogeom_geom::Continuity::G1,
        })
        .collect();
    let filled = crate::fill_n::make_filling_n_within(
        model,
        &boundary,
        &[],
        crate::marched::band_fit_target(tol),
        within,
        tol,
    )?;
    Ok(vec![filled.built.shape])
}
