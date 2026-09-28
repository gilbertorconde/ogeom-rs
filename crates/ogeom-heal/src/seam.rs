//! Seams and poles for periodic faces that arrive without them.
//!
//! A STEP face on a cylinder may be a band between two rings that each go
//! round once, with no seam edge: one ring a whole circle, the other a
//! notched run of arcs and rulings. A face on a sphere may be a cap bounded
//! by one circle, or an octant whose two meridians meet at the pole with no
//! edge along the pole's row. Each is a valid solid's boundary and meshes,
//! but in the surface's chart the boundary does not close, and anything
//! that works in the chart (a boolean's arrangement) cannot use it.
//!
//! The band gets a seam at a column where each ring passes once: a vertex
//! of the notched ring, with the whole circle re-anchored there when it has
//! no vertex of its own on that column. The cap gets a pole vertex, a
//! degenerate pole edge and a seam down to it. The octant gets a degenerate
//! edge along the pole's row between the two meridians.

use std::collections::HashMap;

use ogeom_algo::{
    Built, History, is_shell_closed, make_band_of_rings, make_edge_between, make_face_on,
    make_shell, make_solid, make_vertex, make_wire,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Axis2, Direction2, Point, Point2, Vector2};
use ogeom_topo::{
    EdgeData, EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};

use crate::reanchor::{attach_face_pcurves, circle_parameter};

/// Give every seamless periodic face of a solid the seam or pole edge its
/// chart needs to close.
///
/// Three shapes of face are repaired, each only on a cylinder, cone or
/// sphere and only where every edge has a closed-form pcurve: a band between
/// two rings that each go once round, a sphere's cap bounded by one ring,
/// and a sphere's face reaching a pole with no edge along the pole's row.
/// A face that resists is left as it was. Returns the solid and how many
/// faces were repaired; a solid with none comes back as itself with an
/// empty history.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape is not a solid, or the rebuilt shell does not close.
pub fn seam_periodic_faces(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    if model.kind_of(shape)? != ShapeType::Solid {
        ogeom_bail!(Construction, "seaming heals solids");
    }
    let mut substitution: HashMap<TShapeId, Vec<Shape>> = HashMap::new();
    let mut face_map: HashMap<TShapeId, Shape> = HashMap::new();
    // Edges a rebuilt face already holds, which no later repair may replace.
    let mut settled: std::collections::HashSet<TShapeId> = std::collections::HashSet::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some((surface_id, surface)) = seamless_periodic(model, &face)? else {
            continue;
        };
        let rings: Vec<Vec<Shape>> = model
            .ordered_children_of(&face)?
            .iter()
            .map(|w| model.ordered_children_of(w))
            .collect::<OgeomResult<_>>()?;
        // A face whose edges another repair here already replaced takes
        // them as replaced where it can: a cap on a re-anchored circle
        // seams wherever the circle now starts. Anything else is left
        // alone, since its rings would have to agree with that repair's
        // column.
        let resolved: Vec<Vec<Shape>> = rings
            .iter()
            .map(|ring| {
                ring.iter()
                    .flat_map(|e| match substitution.get(&e.node()) {
                        Some(new_edges) if e.orientation() == ogeom_topo::Orientation::Reversed => {
                            new_edges.iter().rev().map(Shape::reversed).collect()
                        }
                        Some(new_edges) => new_edges.clone(),
                        None => vec![e.clone()],
                    })
                    .collect()
            })
            .collect();
        if resolved != rings {
            let rebuilt = match &resolved[..] {
                [ring] if lone_circle(model, ring).is_some() => {
                    cap(model, &face, &surface, ring, tol)
                }
                _ => Ok(None),
            };
            if let Ok(Some(rebuilt)) = rebuilt {
                settled.extend(
                    explore_unique(model, &rebuilt, ShapeType::Edge)?
                        .iter()
                        .map(Shape::node),
                );
                face_map.insert(face.node(), rebuilt);
            }
            continue;
        }
        let walks: Vec<Option<Walk>> = rings
            .iter()
            .map(|ring| walk(model, surface_id, &surface, ring, tol))
            .collect::<OgeomResult<_>>()?;
        let repaired = match (rings.len(), &walks[..]) {
            (2, [Some(a), Some(b)]) if a.wraps && b.wraps => band(
                model,
                &surface,
                [&rings[0], &rings[1]],
                [a, b],
                (&mut substitution, &settled),
                tol,
            ),
            (1, [Some(a)]) if a.wraps => match (0..rings[0].len()).find(|&i| a.passes_once_at(i)) {
                Some(i) => {
                    let mut ring = rings[0].clone();
                    ring.rotate_left(i);
                    cap(model, &face, &surface, &ring, tol)
                }
                None => Ok(None),
            },
            (1, [Some(a)]) => pole_gap(model, surface_id, &surface, &rings[0], a, tol),
            _ => Ok(None),
        };
        // A face the repair cannot rebuild is kept unchanged.
        if let Ok(Some(rebuilt)) = repaired {
            settled.extend(
                explore_unique(model, &face, ShapeType::Edge)?
                    .iter()
                    .map(Shape::node),
            );
            face_map.insert(face.node(), rebuilt);
        }
    }
    if face_map.is_empty() {
        return Ok((Built::new(shape.clone(), History::new()), 0));
    }
    let repaired = face_map.len();

    // Neighbours of a re-anchored ring use its new edge.
    let mut history = History::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        if face_map.contains_key(&face.node()) {
            continue;
        }
        let uses_any = explore_unique(model, &face, ShapeType::Edge)?
            .iter()
            .any(|e| substitution.contains_key(&e.node()));
        if uses_any {
            let rebuilt = rebuild_with(model, &face, &substitution, tol)?;
            face_map.insert(face.node(), rebuilt);
        }
    }
    for face in explore_unique(model, shape, ShapeType::Face)? {
        if let Some(rebuilt) = face_map.get(&face.node()) {
            let oriented = if face.orientation() == ogeom_topo::Orientation::Reversed {
                rebuilt.reversed()
            } else {
                rebuilt.clone()
            };
            history.modify(&face, oriented);
        }
    }
    let mut shells = Vec::new();
    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        let faces: Vec<Shape> = explore(model, &shell, Filter::OfType(ShapeType::Face))?
            .into_iter()
            .map(|f| {
                face_map.get(&f.node()).map_or(f.clone(), |n| {
                    if f.orientation() == ogeom_topo::Orientation::Reversed {
                        n.reversed()
                    } else {
                        n.clone()
                    }
                })
            })
            .collect();
        let rebuilt = make_shell(model, &faces)?.shape;
        if !is_shell_closed(model, &rebuilt)? {
            ogeom_bail!(
                Construction,
                "seaming left the shell open; the shape resists this repair"
            );
        }
        history.modify(&shell, rebuilt.clone());
        shells.push(rebuilt);
    }
    let solid = make_solid(model, &shells)?.shape;
    history.modify(shape, solid.clone());
    Ok((Built::new(solid, history), repaired))
}

/// The face's surface, where it closes in `u` only, is a cylinder, cone or
/// sphere, and the face has no seam on it.
fn seamless_periodic(
    model: &Model,
    face: &Shape,
) -> OgeomResult<Option<(ogeom_topo::SurfaceId, SurfaceGeometry)>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    if !surface.is_periodic_u()
        || surface.is_periodic_v()
        || !matches!(
            surface,
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) | SurfaceGeometry::Sphere(_)
        )
    {
        return Ok(None);
    }
    for wire in model.ordered_children_of(face)? {
        let edges = model.ordered_children_of(&wire)?;
        for edge in &edges {
            if edges.iter().filter(|e| e.node() == edge.node()).count() > 1 {
                return Ok(None);
            }
            let seam = model
                .node(edge)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.pcurve_for(data.surface, edge.location()))
                .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }));
            if seam {
                return Ok(None);
            }
        }
    }
    Ok(Some((data.surface, surface.clone())))
}

/// A ring walked through its stored pcurves, each piece moved by whole
/// periods to continue the one before.
struct Walk {
    /// Each occurrence's samples, first to last.
    samples: Vec<Vec<Point2>>,
    /// Whether the walk goes once round in `u`.
    wraps: bool,
    /// How far round it goes: plus or minus the period, or less.
    turn: f64,
    span: f64,
    /// The surface the pcurves were read on.
    surface: ogeom_topo::SurfaceId,
    /// How far in `u` each occurrence's stored pcurve was moved.
    shifts: Vec<f64>,
}

impl Walk {
    /// Occurrence `i`'s start, in the walk's unwrapped chart.
    fn start(&self, i: usize) -> Point2 {
        self.samples[i][0]
    }

    /// Whether the column through occurrence `i`'s start meets the ring
    /// there only: every other sample stands strictly between the column
    /// and its copy a turn along.
    fn passes_once_at(&self, i: usize) -> bool {
        let from = self.start(i).x;
        let n = self.samples.len();
        let margin = 1e-6 * self.span;
        for k in 0..n {
            let occurrence = (i + k) % n;
            let shift = if i + k >= n { self.turn } else { 0.0 };
            for (j, p) in self.samples[occurrence].iter().enumerate() {
                if k == 0 && j == 0 {
                    continue;
                }
                if k == n - 1 && j == self.samples[occurrence].len() - 1 {
                    continue;
                }
                let along = (p.x + shift - from) / self.turn;
                if along * self.span.abs() <= margin || (1.0 - along) * self.span <= margin {
                    return false;
                }
            }
        }
        true
    }
}

fn walk(
    model: &Model,
    surface_id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    ring: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<Walk>> {
    const SAMPLES: u32 = 16;
    let ((ua, ub), _) = surface.domain();
    let span = ub - ua;
    let mut samples: Vec<Vec<Point2>> = Vec::with_capacity(ring.len());
    let mut shifts = Vec::with_capacity(ring.len());
    let mut at: Option<Point2> = None;
    for occurrence in ring {
        let Some(data) = model.node(occurrence).and_then(|n| n.data().as_edge()) else {
            return Ok(None);
        };
        let Some(EdgeRepr::PCurve { curve, range, .. }) =
            data.pcurve_for(surface_id, occurrence.location())
        else {
            return Ok(None);
        };
        let Some(pcurve) = model.geometry().pcurve(*curve) else {
            return Ok(None);
        };
        let reversed = occurrence.orientation() == ogeom_topo::Orientation::Reversed;
        let mut piece = Vec::with_capacity(SAMPLES as usize + 1);
        for k in 0..=SAMPLES {
            let f = f64::from(k) / f64::from(SAMPLES);
            let t = if reversed {
                range.1 - (range.1 - range.0) * f
            } else {
                range.0 + (range.1 - range.0) * f
            };
            piece.push(pcurve.point_at(t, tol)?);
        }
        let shift = at.map_or(0.0, |at| ((at.x - piece[0].x) / span).round() * span);
        for p in &mut piece {
            p.x += shift;
        }
        shifts.push(shift);
        at = piece.last().copied();
        samples.push(piece);
    }
    let (Some(first), Some(last)) = (samples.first(), at) else {
        return Ok(None);
    };
    let turn = last.x - first[0].x;
    let wraps = (turn.abs() - span).abs() <= 1e-6 * span;
    Ok(Some(Walk {
        samples,
        wraps,
        turn,
        span,
        surface: surface_id,
        shifts,
    }))
}

/// A ring's single closed circle edge, if that is all it is.
fn lone_circle(model: &Model, ring: &[Shape]) -> Option<(Curve, (f64, f64))> {
    let [edge] = ring else {
        return None;
    };
    let data = model.node(edge).and_then(|n| n.data().as_edge())?;
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return None;
    };
    let curve = model.geometry().curve(*curve)?;
    matches!(curve, Curve::Circle(_)).then(|| (curve.clone(), *range))
}

/// A band between two wrapping rings, seamed at a column clear of every
/// vertex where each ring passes once: a lone circle is re-anchored there,
/// any other ring has the edge crossing it split.
fn band(
    model: &mut Model,
    surface: &SurfaceGeometry,
    rings: [&[Shape]; 2],
    walks: [&Walk; 2],
    (substitution, settled): (
        &mut HashMap<TShapeId, Vec<Shape>>,
        &std::collections::HashSet<TShapeId>,
    ),
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let span = walks[0].span;
    // An edge a rebuilt face holds cannot be replaced.
    if rings
        .iter()
        .flat_map(|r| r.iter())
        .any(|e| settled.contains(&e.node()))
    {
        return Ok(None);
    }
    // Candidate columns midway between neighbouring vertex columns, the
    // widest gaps first.
    let mut columns: Vec<f64> = walks
        .iter()
        .flat_map(|w| (0..w.samples.len()).map(|i| w.start(i).x.rem_euclid(span)))
        .collect();
    columns.sort_by(f64::total_cmp);
    columns.dedup_by(|a, b| (*a - *b).abs() <= 1e-9 * span);
    let mut gaps: Vec<(f64, f64)> = (0..columns.len())
        .map(|k| {
            let (a, b) = (
                columns[k],
                columns.get(k + 1).copied().unwrap_or(columns[0] + span),
            );
            (b - a, f64::midpoint(a, b))
        })
        .collect();
    gaps.sort_by(|x, y| y.0.total_cmp(&x.0));
    for (_, column) in gaps {
        let plans = [
            plan(model, rings[0], walks[0], column, tol)?,
            plan(model, rings[1], walks[1], column, tol)?,
        ];
        let [Some(first), Some(second)] = plans else {
            continue;
        };
        let mut built = Vec::with_capacity(2);
        let mut added = Vec::new();
        for (ring, planned) in rings.iter().zip([first, second]) {
            let (ring, replaced) =
                apply(model, surface, walks[0].surface, ring, planned, column, tol)?;
            built.push(ring);
            added.extend(replaced);
        }
        match make_band_of_rings(model, surface, &built[0], &built[1], tol) {
            Ok(face) => {
                substitution.extend(added);
                return Ok(Some(face));
            }
            Err(_) => continue,
        }
    }
    Ok(None)
}

/// How one ring comes to have a vertex on a column.
enum Plan {
    /// A lone circle, anchored anew there.
    Reanchor,
    /// Occurrence `k` crosses there, at its stored pcurve's parameter `t`.
    Split { k: usize, t: f64 },
}

/// How a ring gets a vertex on `column`, where it passes the column once
/// and not at a vertex.
fn plan(
    model: &Model,
    ring: &[Shape],
    walk: &Walk,
    column: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Plan>> {
    if lone_circle(model, ring).is_some() {
        return Ok(Some(Plan::Reanchor));
    }
    let span = walk.span;
    let level = |p: Point2| ((p.x - column) / span).floor();
    let margin = 1e-6;
    let mut crossing = None;
    let mut crossings = 0.0;
    for (k, piece) in walk.samples.iter().enumerate() {
        for pair in piece.windows(2) {
            let near = |p: Point2| {
                let g = (p.x - column) / span;
                (g - g.round()).abs() <= margin
            };
            if near(pair[0]) || near(pair[1]) {
                return Ok(None);
            }
            let steps = (level(pair[1]) - level(pair[0])).abs();
            if steps > 0.0 {
                crossings += steps;
                crossing = Some((k, pair[0], pair[1]));
            }
        }
    }
    let Some((k, _, _)) = crossing else {
        return Ok(None);
    };
    if crossings != 1.0 {
        return Ok(None);
    }
    // Where along the occurrence's stored pcurve the column falls.
    let Some((pcurve, range, shift)) = stored_piece(model, walk, ring, k)? else {
        return Ok(None);
    };
    let target = {
        let at = pcurve.point_at(range.0, tol)?.x + shift;
        column + ((at - column) / span).round() * span
    };
    let u = |t: f64| -> OgeomResult<f64> { Ok(pcurve.point_at(t, tol)?.x + shift) };
    let (mut lo, mut hi) = range;
    let (mut f_lo, f_hi) = (u(lo)? - target, u(hi)? - target);
    if f_lo.signum() == f_hi.signum() {
        return Ok(None);
    }
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        let f_mid = u(mid)? - target;
        if f_mid.signum() == f_lo.signum() {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
    }
    Ok(Some(Plan::Split {
        k,
        t: f64::midpoint(lo, hi),
    }))
}

/// A stored pcurve, its range, and the shift a walk moved it by.
type StoredPiece = (ogeom_geom::PlanarCurve, (f64, f64), f64);

/// Occurrence `k`'s stored pcurve, its range, and the shift the walk moved
/// it by.
fn stored_piece(
    model: &Model,
    walk: &Walk,
    ring: &[Shape],
    k: usize,
) -> OgeomResult<Option<StoredPiece>> {
    let occurrence = &ring[k];
    let Some(data) = model.node(occurrence).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let surface = walk.surface;
    let Some(EdgeRepr::PCurve { curve, range, .. }) =
        data.pcurve_for(surface, occurrence.location())
    else {
        return Ok(None);
    };
    let Some(pcurve) = model.geometry().pcurve(*curve) else {
        return Ok(None);
    };
    Ok(Some((pcurve.clone(), *range, walk.shifts[k])))
}

/// Give a ring its vertex on `column` as planned: the ring's occurrences
/// with the new edges in place, starting at that vertex, and what each
/// replaced edge became, in its stored direction.
#[allow(clippy::type_complexity)]
fn apply(
    model: &mut Model,
    surface: &SurfaceGeometry,
    surface_id: ogeom_topo::SurfaceId,
    ring: &[Shape],
    planned: Plan,
    column: f64,
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, Vec<(TShapeId, Vec<Shape>)>)> {
    let edge_curve = |model: &Model, edge: &Shape| -> OgeomResult<(Curve, (f64, f64))> {
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.curve3d())
        else {
            ogeom_bail!(Construction, "a ring edge has no curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        Ok((geometry.clone(), *range))
    };
    let reversed = |e: &Shape| e.orientation() == ogeom_topo::Orientation::Reversed;
    match planned {
        Plan::Reanchor => {
            let old = &ring[0];
            let (curve, range) = edge_curve(model, old)?;
            let row = {
                let Some((a, _)) = ogeom_algo::edge_vertices(model, old)? else {
                    ogeom_bail!(Construction, "a ring has no vertex");
                };
                let Some(data) = model.node(&a).and_then(|n| n.data().as_vertex()) else {
                    ogeom_bail!(Dangling, "vertex is not in this model");
                };
                ogeom_algo::measure::project_on_surface(surface, data.point, 32, tol)?
                    .parameters
                    .1
            };
            let target = surface.point_at(column, row, tol)?;
            let Some(at) = circle_parameter(&curve, target) else {
                ogeom_bail!(Construction, "an anchor point fell off its circle");
            };
            let period = range.1 - range.0;
            let vertex = make_vertex(model, target).shape;
            let edge =
                make_edge_between(model, curve, (at, at + period), &vertex, &vertex, tol)?.shape;
            let placed = if reversed(old) {
                edge.reversed()
            } else {
                edge.clone()
            };
            Ok((vec![placed], vec![(old.node(), vec![edge])]))
        }
        Plan::Split { k, t } => {
            let old = &ring[k];
            let (curve, range) = edge_curve(model, old)?;
            // The stored pcurve and the edge share their parameter
            // proportionally across their ranges.
            let Some(EdgeRepr::PCurve { range: prange, .. }) = model
                .node(old)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.pcurve_for(surface_id, old.location()))
                .cloned()
            else {
                ogeom_bail!(Construction, "a ring edge has no pcurve");
            };
            let at = range.0 + (t - prange.0) / (prange.1 - prange.0) * (range.1 - range.0);
            let bounds = model.children_of(old)?;
            let (Some(first), Some(last)) = (bounds.first(), bounds.last()) else {
                ogeom_bail!(Construction, "a ring edge has no vertices");
            };
            let point = curve.point_at(at, tol)?;
            let vertex = make_vertex(model, point).shape;
            let head =
                make_edge_between(model, curve.clone(), (range.0, at), first, &vertex, tol)?.shape;
            let tail = make_edge_between(model, curve, (at, range.1), &vertex, last, tol)?.shape;
            // In the ring's direction: the occurrence's two halves, the
            // second starting on the column.
            let halves = if reversed(old) {
                [tail.reversed(), head.reversed()]
            } else {
                [head.clone(), tail.clone()]
            };
            let mut out: Vec<Shape> = Vec::with_capacity(ring.len() + 1);
            out.push(halves[1].clone());
            out.extend(ring[k + 1..].iter().cloned());
            out.extend(ring[..k].iter().cloned());
            out.push(halves[0].clone());
            Ok((out, vec![(old.node(), vec![head, tail])]))
        }
    }
}

/// A sphere's cap bounded by one wrapping ring, given from a vertex the
/// seam may start at: a pole vertex inside it, a degenerate edge there, and
/// a seam down to the ring.
fn cap(
    model: &mut Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    ring: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let SurfaceGeometry::Sphere(_) = surface else {
        return Ok(None);
    };
    let half = core::f64::consts::FRAC_PI_2;
    // The pole the face covers: its mesh reaches that one. The chart does
    // not close, so asking the face itself of a point is unreliable; its
    // triangulation works from the edges and is not.
    let mesh = ogeom_mesh::triangulate(model, face, ogeom_mesh::Deflection::default(), tol)?;
    let reach = |pole: Point| {
        mesh.positions
            .iter()
            .map(|p| p.distance(pole))
            .fold(f64::INFINITY, f64::min)
    };
    let poles = [
        surface.point_at(0.0, half, tol)?,
        surface.point_at(0.0, -half, tol)?,
    ];
    let pole = if mesh.positions.is_empty() {
        None
    } else if reach(poles[0]) < reach(poles[1]) {
        Some(poles[0])
    } else {
        Some(poles[1])
    };
    let Some(pole) = pole else {
        return Ok(None);
    };
    let apex = make_vertex(model, pole).shape;
    let mut data = EdgeData::new();
    data.degenerate = true;
    let degenerate = model.add_edge(data, &[apex.clone(), apex])?;
    Ok(Some(make_band_of_rings(
        model,
        surface,
        ring,
        &[degenerate],
        tol,
    )?))
}

/// A sphere's face whose boundary reaches a pole from one column and leaves
/// it on another: a degenerate edge along the pole's row between them.
fn pole_gap(
    model: &mut Model,
    surface_id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    ring: &[Shape],
    walk: &Walk,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let SurfaceGeometry::Sphere(_) = surface else {
        return Ok(None);
    };
    let half = core::f64::consts::FRAC_PI_2;
    let n = ring.len();
    let mut edges = Vec::with_capacity(n + 1);
    let mut inserted = false;
    for (k, occurrence) in ring.iter().enumerate() {
        edges.push(occurrence.clone());
        let end = *walk.samples[k].last().unwrap_or(&walk.samples[k][0]);
        let next = walk.samples[(k + 1) % n][0];
        // The walk moved each piece to continue the last, so a piece
        // leaving the pole on another column shows as a jump along the row;
        // the ring does not go round, so it closes where it started.
        let gap = next.x - end.x;
        if (end.y.abs() - half).abs() > 1e-9 || gap.abs() <= 1e-9 * walk.span {
            continue;
        }
        let Some((_, vertex)) = ogeom_algo::edge_vertices(model, occurrence)? else {
            return Ok(None);
        };
        let mut data = EdgeData::new();
        data.degenerate = true;
        let degenerate = model.add_edge(data, &[vertex.clone(), vertex])?;
        // The pole's pcurve joins the stored pieces where they stand.
        let stored_end = stored_end(model, surface_id, occurrence, tol)?;
        let Some(stored_end) = stored_end else {
            return Ok(None);
        };
        let pcurve: ogeom_geom::PlanarCurve = ogeom_geom::Line2d::over(
            Axis2::new(
                stored_end,
                Direction2::new(Vector2::new(gap.signum(), 0.0), tol)?,
            ),
            0.0,
            gap.abs(),
        )?
        .into();
        ogeom_algo::attach_pcurve(
            model,
            &degenerate,
            pcurve,
            surface_id,
            Location::identity(),
            (0.0, gap.abs()),
        )?;
        edges.push(degenerate);
        inserted = true;
    }
    if !inserted {
        return Ok(None);
    }
    // Which way the closed loop goes round the chart, read on the pcurves
    // as stored, where the pole's edges join them: clockwise, the face
    // faces against the surface.
    let mut area = 0.0;
    let stored: Vec<Point2> = walk
        .samples
        .iter()
        .zip(&walk.shifts)
        .flat_map(|(piece, shift)| piece.iter().map(move |p| Point2::new(p.x - shift, p.y)))
        .collect();
    for (k, p) in stored.iter().enumerate() {
        let q = stored[(k + 1) % stored.len()];
        area += p.x * q.y - q.x * p.y;
    }
    let wire = make_wire(model, &edges, tol)?.shape;
    let face = make_face_on(model, surface_id, &[wire], tol)?.shape;
    Ok(Some(if area < 0.0 { face.reversed() } else { face }))
}

/// Where an occurrence's stored pcurve ends, as stored.
fn stored_end(
    model: &Model,
    surface_id: ogeom_topo::SurfaceId,
    occurrence: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Point2>> {
    let Some(EdgeRepr::PCurve { curve, range, .. }) = model
        .node(occurrence)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| d.pcurve_for(surface_id, occurrence.location()))
    else {
        return Ok(None);
    };
    let Some(pcurve) = model.geometry().pcurve(*curve) else {
        return Ok(None);
    };
    let t = if occurrence.orientation() == ogeom_topo::Orientation::Reversed {
        range.0
    } else {
        range.1
    };
    Ok(Some(pcurve.point_at(t, tol)?))
}

/// Rebuild a face around replaced edges, each by the edges it became (in
/// its stored direction), pcurves on the face's surface found in closed
/// form for the new ones.
fn rebuild_with(
    model: &mut Model,
    face: &Shape,
    substitution: &HashMap<TShapeId, Vec<Shape>>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let surface_id = data.surface;
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let mut wires = Vec::new();
    let mut replaced = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let mut edges = Vec::new();
        for edge in model.ordered_children_of(&wire)? {
            match substitution.get(&edge.node()) {
                Some(new_edges) => {
                    replaced.extend(new_edges.iter().cloned());
                    if edge.orientation() == ogeom_topo::Orientation::Reversed {
                        edges.extend(new_edges.iter().rev().map(Shape::reversed));
                    } else {
                        edges.extend(new_edges.iter().cloned());
                    }
                }
                None => edges.push(edge),
            }
        }
        wires.push(make_wire(model, &edges, tol)?.shape);
    }
    attach_face_pcurves(model, &surface, surface_id, &replaced, tol)?;
    Ok(make_face_on(model, surface_id, &wires, tol)?.shape)
}
