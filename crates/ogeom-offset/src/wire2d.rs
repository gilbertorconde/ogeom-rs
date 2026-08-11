//! Offsetting a closed planar wire: the sketch-plane operation everything in
//! this crate stands on.
//!
//! Each edge is offset on its own — a segment to the parallel segment, an arc
//! to the concentric arc — and the corners decide the rest. Where the corner
//! turns *away* from the offset side the pieces separate, and the gap is
//! closed by the chosen join: an arc about the old corner, or the extension
//! of both pieces to their meeting. Where it turns *toward* the offset side
//! the pieces overlap, and both are trimmed to their intersection. The result
//! is a new wire with history: every input edge modified into its offset,
//! every join generated from the corner it rounds.
//!
//! The honest limits, refused by name rather than mishandled: edges that are
//! neither straight nor circular, offsets that consume an edge whole, arcs
//! whose concentric offset would have no radius left, and results that
//! self-intersect — the global arrangement that trims a collapsed offset into
//! its valid loops is recorded in docs/PARITY.md (offset.wire-offset).

use ogeom_algo::{
    Built, History, edge_vertices, find_plane, make_edge_between, make_vertex, make_wire,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{CircleCurve, Curve, LineCurve, PlanarCurve};
use ogeom_math::{Circle, Frame, Point, Point2, Vector2};
use ogeom_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, explore};

/// How a gap at a corner is closed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Join {
    /// An arc about the old corner, radius the offset — the rounded corner.
    Arc,
    /// Both pieces extended to their meeting — the sharp corner. Supported
    /// where both sides are straight; extending arcs to a meeting that may
    /// not exist is refused.
    Intersection,
}

/// A piece of the offset wire, in the plane's own coordinates.
#[derive(Debug, Clone)]
enum Piece {
    Seg {
        from: Point2,
        to: Point2,
    },
    /// Angles are unwrapped along traversal: increasing for a
    /// counter-clockwise piece, decreasing for a clockwise one.
    Arc {
        centre: Point2,
        radius: f64,
        start: f64,
        end: f64,
    },
}

impl Piece {
    fn start_point(&self) -> Point2 {
        match self {
            Self::Seg { from, .. } => *from,
            Self::Arc {
                centre,
                radius,
                start,
                ..
            } => at_angle(*centre, *radius, *start),
        }
    }

    fn end_point(&self) -> Point2 {
        match self {
            Self::Seg { to, .. } => *to,
            Self::Arc {
                centre,
                radius,
                end,
                ..
            } => at_angle(*centre, *radius, *end),
        }
    }

    /// Unit tangent along traversal at the start or end.
    fn tangent(&self, at_end: bool) -> Vector2 {
        match self {
            Self::Seg { from, to } => {
                let d = *to - *from;
                d / d.magnitude()
            }
            Self::Arc {
                centre,
                radius,
                start,
                end,
            } => {
                let a = if at_end { *end } else { *start };
                let radial = (at_angle(*centre, *radius, a) - *centre) / *radius;
                let ccw = end > start;
                if ccw {
                    Vector2::new(-radial.y, radial.x)
                } else {
                    Vector2::new(radial.y, -radial.x)
                }
            }
        }
    }
}

fn at_angle(centre: Point2, radius: f64, angle: f64) -> Point2 {
    Point2::new(
        radius.mul_add(angle.cos(), centre.x),
        radius.mul_add(angle.sin(), centre.y),
    )
}

/// Offset a closed planar wire by `offset`: positive moves outward, away from
/// the enclosed region, negative moves inward.
///
/// The wire's edges must be straight or circular. Gaps at corners are closed
/// per `join`; overlaps are trimmed to the pieces' intersection. The history
/// reports every input edge modified into its offset piece, every join
/// generated from the corner vertex it replaces, and the wire modified into
/// the result.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the wire is
/// open or not planar, an edge is neither straight nor circular, the offset
/// consumes an edge or an arc's radius, an `Intersection` join is asked of a
/// curved side, or the offset self-intersects.
pub fn offset_wire(
    model: &mut Model,
    wire: &Shape,
    offset: f64,
    join: Join,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !offset.is_finite() || offset.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "an offset of {offset} moves nothing");
    }
    if model.kind_of(wire)? != ShapeType::Wire {
        ogeom_bail!(Construction, "offsetting starts from a wire");
    }
    let open = !ogeom_algo::is_wire_closed(model, wire, tol)?;
    let plane = match find_plane(model, wire, tol)? {
        Some(plane) => plane,
        None if open => {
            // A straight open path spans no plane of its own; its outline
            // is symmetric about it, so any plane holding the line serves.
            let ends = explore(model, wire, Filter::OfType(ShapeType::Vertex))?;
            let mut points = Vec::new();
            for v in &ends {
                if let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) {
                    points.push(v.transform(model.datums())?.apply(data.point));
                }
            }
            if points.len() < 2 {
                ogeom_bail!(Construction, "the wire is not planar");
            }
            let along = ogeom_math::Direction::new(points[1] - points[0], tol)?;
            let reference = if along.vector().cross(ogeom_math::Vector::Z).magnitude() > 0.5 {
                ogeom_math::Direction::Z
            } else {
                ogeom_math::Direction::X
            };
            let normal = ogeom_math::Direction::new(along.vector().cross(reference.vector()), tol)?;
            ogeom_math::Plane::through(points[0], normal)
        }
        None => ogeom_bail!(Construction, "the wire is not planar"),
    };
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };

    // Every edge as a piece in plane coordinates, traversal order.
    let edges = explore(model, wire, Filter::OfType(ShapeType::Edge))?;
    let mut pieces: Vec<Piece> = Vec::with_capacity(edges.len());
    for edge in &edges {
        let (curve, range) = {
            let Some(node) = model.node(edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "an edge has no curve to offset");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let Some((sv, ev)) = edge_vertices(model, edge)? else {
            ogeom_bail!(Construction, "an edge has no bounding vertices");
        };
        let position = |v: &Shape| -> OgeomResult<Point> {
            let Some(node) = model.node(v) else {
                ogeom_bail!(Dangling, "vertex is not in this model");
            };
            let Some(data) = node.data().as_vertex() else {
                ogeom_bail!(Construction, "vertex node holds no point");
            };
            Ok(v.transform(model.datums())?.apply(data.point))
        };
        let from = flat(position(&sv)?);
        let to = flat(position(&ev)?);
        match &curve {
            Curve::Line(_) => pieces.push(Piece::Seg { from, to }),
            Curve::Circle(c) => {
                let centre = flat(c.circle().centre());
                let radius = c.circle().radius();
                let mid = flat(curve.point_at(f64::midpoint(range.0, range.1), tol)?);
                // The traversal's sweep direction, read from three points of
                // the arc rather than from flags that compose.
                let ccw = (mid - from).cross(to - mid) > 0.0;
                let a0 = (from - centre).y.atan2((from - centre).x);
                let mut a1 = (to - centre).y.atan2((to - centre).x);
                let tau = core::f64::consts::TAU;
                if ccw {
                    while a1 <= a0 + tol.parametric() {
                        a1 += tau;
                    }
                } else {
                    while a1 >= a0 - tol.parametric() {
                        a1 -= tau;
                    }
                }
                pieces.push(Piece::Arc {
                    centre,
                    radius,
                    start: a0,
                    end: a1,
                });
            }
            _ => ogeom_bail!(
                Construction,
                "offsetting an edge that is neither straight nor circular \
                 needs the offset-curve machinery — docs/PARITY.md, offset.wire-offset"
            ),
        }
    }

    // An open wire offsets on both sides and closes over its ends: the
    // outline of the thick path. The return pass is the same pieces walked
    // backwards, and every corner rule below then applies unchanged — the
    // ends become 180-degree corners, which is what the caps close.
    let source_count = pieces.len();
    if open {
        for piece in pieces.clone().iter().rev() {
            pieces.push(match piece {
                Piece::Seg { from, to } => Piece::Seg {
                    from: *to,
                    to: *from,
                },
                Piece::Arc {
                    centre,
                    radius,
                    start,
                    end,
                } => Piece::Arc {
                    centre: *centre,
                    radius: *radius,
                    start: *end,
                    end: *start,
                },
            });
        }
    }

    // The wire's winding decides which side is out. Signed area by shoelace
    // over a sampling fine enough for the decision it makes. An open path's
    // outline encloses nothing yet; its offset surrounds it, so the side is
    // fixed and the amount is the magnitude.
    let offset = if open { offset.abs() } else { offset };
    let winding = if open {
        1.0
    } else {
        let mut area = 0.0;
        let mut samples: Vec<Point2> = Vec::new();
        for piece in &pieces {
            match piece {
                Piece::Seg { from, .. } => samples.push(*from),
                Piece::Arc {
                    centre,
                    radius,
                    start,
                    end,
                } => {
                    for i in 0..32 {
                        let a = start + (end - start) * f64::from(i) / 32.0;
                        samples.push(at_angle(*centre, *radius, a));
                    }
                }
            }
        }
        for i in 0..samples.len() {
            let (p, q) = (samples[i], samples[(i + 1) % samples.len()]);
            area += p.x.mul_add(q.y, -(q.x * p.y));
        }
        if area.abs() <= tol.confusion() {
            ogeom_bail!(Construction, "the wire encloses no area to offset");
        }
        area.signum()
    };
    let _ = source_count;
    // Rotating the traversal tangent by -90° (times the winding) points out
    // of the enclosed region; `offset` moves along it.
    let outward = |tangent: Vector2| Vector2::new(tangent.y, -tangent.x) * winding;

    // Each piece offset on its own support.
    let mut moved: Vec<Piece> = Vec::with_capacity(pieces.len());
    for piece in &pieces {
        match piece {
            Piece::Seg { from, to } => {
                let shift = outward(piece.tangent(false)) * offset;
                moved.push(Piece::Seg {
                    from: *from + shift,
                    to: *to + shift,
                });
            }
            Piece::Arc {
                centre,
                radius,
                start,
                end,
            } => {
                // Whether "outward" is radially out here depends on the arc's
                // own sweep against the wire's winding; the midpoint says.
                let mid = f64::midpoint(*start, *end);
                let radial = (at_angle(*centre, *radius, mid) - *centre) / *radius;
                let tangent_mid = if end > start {
                    Vector2::new(-radial.y, radial.x)
                } else {
                    Vector2::new(radial.y, -radial.x)
                };
                let sign = outward(tangent_mid).dot(radial).signum();
                let grown = radius + offset * sign;
                if grown <= tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "the offset consumes the arc's radius entirely"
                    );
                }
                moved.push(Piece::Arc {
                    centre: *centre,
                    radius: grown,
                    start: *start,
                    end: *end,
                });
            }
        }
    }

    let n = moved.len();
    let mut chain: Vec<(Piece, Provenance)> = Vec::with_capacity(n * 2);
    for (i, piece) in moved.iter().enumerate() {
        chain.push((piece.clone(), Provenance::Offset(i)));
    }
    // Corners, walked over the original geometry: the turn between the
    // incoming and outgoing tangents against the offset side says gap or
    // overlap.
    for i in 0..n {
        let j = (i + 1) % n;
        let turn = pieces[i].tangent(true).cross(pieces[j].tangent(false));
        let at_i = chain
            .iter()
            .position(|(_, p)| *p == Provenance::Offset(i))
            .unwrap_or(0);
        let at_j = chain
            .iter()
            .position(|(_, p)| *p == Provenance::Offset(j))
            .unwrap_or(0);
        let e = chain[at_i].0.end_point();
        let s = chain[at_j].0.start_point();
        if e.distance(s) <= tol.confusion() * 10.0 {
            continue; // Tangent-continuous: nothing to do.
        }
        let corner = pieces[i].end_point();
        // A 180-degree corner — an open path's end — is a gap whatever the
        // turn's vanishing cross product says.
        let is_cap =
            turn.abs() <= 1e-9 && pieces[i].tangent(true).dot(pieces[j].tangent(false)) < 0.0;
        if is_cap && join == Join::Intersection {
            // The square cap: both offsets extended one width past the end,
            // joined across.
            let d = pieces[i].tangent(true);
            let e_ext = e + d * offset.abs();
            let s_ext = s + d * offset.abs();
            chain.insert(
                at_i + 1,
                (Piece::Seg { from: e, to: e_ext }, Provenance::Join(i)),
            );
            chain.insert(
                at_i + 2,
                (
                    Piece::Seg {
                        from: e_ext,
                        to: s_ext,
                    },
                    Provenance::Join(i),
                ),
            );
            chain.insert(
                at_i + 3,
                (Piece::Seg { from: s_ext, to: s }, Provenance::Join(i)),
            );
            continue;
        }
        if is_cap || turn * offset * winding > 0.0 {
            // Gap.
            match join {
                Join::Arc => {
                    let a0 = (e - corner).y.atan2((e - corner).x);
                    let mut a1 = (s - corner).y.atan2((s - corner).x);
                    // The short way round is the way the gap opens — and a
                    // cap's exact semicircle has no short way, so the end
                    // tangent breaks the tie: the cap bulges past the end,
                    // not back through the path.
                    let tau = core::f64::consts::TAU;
                    while a1 - a0 > core::f64::consts::PI {
                        a1 -= tau;
                    }
                    while a0 - a1 > core::f64::consts::PI {
                        a1 += tau;
                    }
                    if ((a1 - a0).abs() - core::f64::consts::PI).abs() < 1e-9 {
                        let mid = at_angle(corner, offset.abs(), f64::midpoint(a0, a1));
                        let ahead = pieces[i].tangent(true);
                        if (mid - corner).dot(ahead) < 0.0 {
                            a1 -= tau * (a1 - a0).signum();
                        }
                    }
                    chain.insert(
                        at_i + 1,
                        (
                            Piece::Arc {
                                centre: corner,
                                radius: offset.abs(),
                                start: a0,
                                end: a1,
                            },
                            Provenance::Join(i),
                        ),
                    );
                }
                Join::Intersection => {
                    let (Piece::Seg { from: f1, to: t1 }, Piece::Seg { from: f2, to: t2 }) =
                        (chain[at_i].0.clone(), chain[at_j].0.clone())
                    else {
                        ogeom_bail!(
                            Construction,
                            "an intersection join between curved sides may \
                             never meet; use the arc join"
                        );
                    };
                    let met = intersect_lines(f1, t1, f2, t2, tol)?;
                    if let Piece::Seg { to, .. } = &mut chain[at_i].0 {
                        *to = met;
                    }
                    if let Piece::Seg { from, .. } = &mut chain[at_j].0 {
                        *from = met;
                    }
                }
            }
        } else {
            // Overlap: trim both to the crossing nearest the corner.
            let met = nearest_crossing(&chain[at_i].0, &chain[at_j].0, corner, tol)?;
            trim_end(&mut chain[at_i].0, met, tol)?;
            trim_start(&mut chain[at_j].0, met, tol)?;
        }
    }
    // Consumed pieces are the collapse announcing itself; drop them and let
    // the distance filter and the reconnection below say what survives.
    chain.retain(|(piece, _)| match piece {
        Piece::Seg { from, to } => from.distance(*to) > tol.confusion(),
        Piece::Arc { start, end, .. } => (end - start).abs() > tol.parametric(),
    });
    if chain.is_empty() {
        ogeom_bail!(Construction, "the offset consumes the wire whole");
    }

    // Where the raw offset crosses itself, split; every sub-piece then
    // stands or falls by the offset's own definition — a point of the true
    // offset boundary is a full offset from the source, and a collapsed
    // sliver is closer.
    let m = chain.len();
    let mut cuts: Vec<Vec<Point2>> = vec![Vec::new(); m];
    for i in 0..m {
        for j in i + 1..m {
            if j == i + 1 || (i == 0 && j == m - 1) {
                continue;
            }
            for p in crossings(&chain[i].0, &chain[j].0, tol)? {
                if within(&chain[i].0, p, tol) && within(&chain[j].0, p, tol) {
                    cuts[i].push(p);
                    cuts[j].push(p);
                }
            }
        }
    }
    let mut resolved: Vec<(Piece, Provenance)> = Vec::new();
    for (k, (piece, provenance)) in chain.iter().enumerate() {
        for sub in split_at(piece, &cuts[k]) {
            resolved.push((sub, *provenance));
        }
    }

    // The source, densely enough to measure against.
    let source: Vec<Point2> = {
        let mut out = Vec::new();
        for piece in pieces.iter().take(source_count) {
            match piece {
                Piece::Seg { from, to } => {
                    out.push(*from);
                    out.push(*to);
                }
                Piece::Arc {
                    centre,
                    radius,
                    start,
                    end,
                } => {
                    for i in 0..=32 {
                        let a = start + (end - start) * f64::from(i) / 32.0;
                        out.push(at_angle(*centre, *radius, a));
                    }
                }
            }
        }
        out
    };
    let source_distance = |p: Point2| -> f64 {
        let mut best = f64::INFINITY;
        for w in source.windows(2) {
            let d = w[1] - w[0];
            let len2 = d.dot(d);
            let t = if len2 > 0.0 {
                ((p - w[0]).dot(d) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            best = best.min(p.distance(w[0] + d * t));
        }
        best
    };
    let keep_beyond = offset.abs() - (tol.confusion() * 1e3).max(offset.abs() * 1e-3);
    let had_cuts = cuts.iter().any(|c| !c.is_empty());
    let survivors: Vec<(Piece, Provenance)> = resolved
        .into_iter()
        .filter(|(piece, _)| {
            let mid = match piece {
                Piece::Seg { from, to } => from.midpoint(*to),
                Piece::Arc {
                    centre,
                    radius,
                    start,
                    end,
                } => at_angle(*centre, *radius, f64::midpoint(*start, *end)),
            };
            source_distance(mid) >= keep_beyond
        })
        .collect();
    if survivors.is_empty() {
        ogeom_bail!(Construction, "the offset consumes the wire whole");
    }

    // Reconnect the survivors into loops by endpoint adjacency.
    let loops: Vec<Vec<(Piece, Provenance)>> = if !had_cuts && survivors.len() == chain.len() {
        vec![survivors]
    } else {
        let eps = tol.confusion() * 100.0;
        let mut pool = survivors;
        let mut out: Vec<Vec<(Piece, Provenance)>> = Vec::new();
        while let Some(first) = pool.pop() {
            let mut current = vec![first];
            loop {
                let tail = current.last().map(|(p, _)| p.end_point());
                let Some(tail) = tail else { break };
                let head = current.first().map(|(p, _)| p.start_point());
                if head.is_some_and(|h| h.distance(tail) <= eps) && current.len() > 1 {
                    out.push(current);
                    break;
                }
                let Some(next) = pool
                    .iter()
                    .position(|(p, _)| p.start_point().distance(tail) <= eps)
                else {
                    // An open fragment closes nothing; the collapse ate its
                    // continuation.
                    break;
                };
                current.push(pool.remove(next));
            }
        }
        if out.is_empty() {
            ogeom_bail!(
                Construction,
                "the offset's survivors close no loop; the collapse consumed \
                 the wire"
            );
        }
        out
    };

    // Lift each loop back to space and rebuild, vertices shared along the
    // way. One loop is a wire; a collapse that split the offset into islands
    // is a compound of them.
    let lift = |p: Point2| frame.origin() + frame.x().vector() * p.x + frame.y().vector() * p.y;
    let normal = frame.z();
    let x_ref = frame.x();
    let mut history = History::new();
    let mut wires: Vec<Shape> = Vec::new();
    for ring in &loops {
        let count = ring.len();
        let mut new_edges: Vec<Shape> = Vec::with_capacity(count);
        let vertices: Vec<Shape> = (0..count)
            .map(|k| make_vertex(model, lift(ring[k].0.start_point())).shape)
            .collect();
        for (k, (piece, provenance)) in ring.iter().enumerate() {
            let from = &vertices[k];
            let to = &vertices[(k + 1) % count];
            let built = match piece {
                Piece::Seg { from: a, to: b } => {
                    let line = LineCurve::segment(lift(*a), lift(*b), tol)?;
                    let curve = Curve::Line(line);
                    let domain = curve.domain();
                    make_edge_between(model, curve, domain, from, to, tol)?.shape
                }
                Piece::Arc {
                    centre,
                    radius,
                    start,
                    end,
                } => {
                    // Always built counter-clockwise about the plane normal;
                    // a clockwise piece enters the wire reversed.
                    let ccw = end > start;
                    let circle =
                        Circle::new(Frame::new(lift(*centre), normal, x_ref, tol)?, *radius, tol)?;
                    let curve = Curve::Circle(CircleCurve::new(circle));
                    let (lo, hi) = if ccw { (*start, *end) } else { (*end, *start) };
                    let (va, vb) = if ccw { (from, to) } else { (to, from) };
                    let edge = make_edge_between(model, curve, (lo, hi), va, vb, tol)?.shape;
                    if ccw { edge } else { edge.reversed() }
                }
            };
            match provenance {
                Provenance::Offset(i) => history.modify(&edges[i % edges.len()], built.clone()),
                Provenance::Join(i) => {
                    // The join stands where the corner vertex stood — for an
                    // open path's cap, the end vertex itself.
                    if let Some((_, corner_vertex)) = edge_vertices(model, &edges[i % edges.len()])?
                    {
                        history.generate(&corner_vertex, built.clone());
                    }
                }
            }
            new_edges.push(built);
        }
        wires.push(make_wire(model, &new_edges, tol)?.shape);
    }
    let shape = if wires.len() == 1 {
        wires.remove(0)
    } else {
        ogeom_algo::build::make_compound(model, &wires)?.shape
    };
    history.modify(wire, shape.clone());
    Ok(Built::new(shape, history))
}

/// Split a piece at the given points on it, in traversal order.
fn split_at(piece: &Piece, points: &[Point2]) -> Vec<Piece> {
    if points.is_empty() {
        return vec![piece.clone()];
    }
    match piece {
        Piece::Seg { from, to } => {
            let d = *to - *from;
            let len = d.magnitude();
            let mut ts: Vec<f64> = points
                .iter()
                .map(|p| (*p - *from).dot(d / len))
                .filter(|t| *t > 1e-12 && *t < len - 1e-12)
                .collect();
            ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            ts.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            let mut out = Vec::with_capacity(ts.len() + 1);
            let mut last = *from;
            for t in ts {
                let p = *from + (d / len) * t;
                out.push(Piece::Seg { from: last, to: p });
                last = p;
            }
            out.push(Piece::Seg {
                from: last,
                to: *to,
            });
            out
        }
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => {
            let toward = (end - start).signum();
            let mut ts: Vec<f64> = points
                .iter()
                .filter_map(|p| {
                    let v = *p - *centre;
                    let mut a = v.y.atan2(v.x);
                    let tau = core::f64::consts::TAU;
                    while (a - start) * toward < 0.0 {
                        a += tau * toward;
                    }
                    while (a - end) * toward > 0.0 {
                        a -= tau * toward;
                    }
                    ((a - start) * toward > 1e-12 && (end - a) * toward > 1e-12).then_some(a)
                })
                .collect();
            ts.sort_by(|a, b| {
                ((a - start) * toward)
                    .partial_cmp(&((b - start) * toward))
                    .unwrap_or(core::cmp::Ordering::Equal)
            });
            ts.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            let mut out = Vec::with_capacity(ts.len() + 1);
            let mut last = *start;
            for t in ts {
                out.push(Piece::Arc {
                    centre: *centre,
                    radius: *radius,
                    start: last,
                    end: t,
                });
                last = t;
            }
            out.push(Piece::Arc {
                centre: *centre,
                radius: *radius,
                start: last,
                end: *end,
            });
            out
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Provenance {
    Offset(usize),
    Join(usize),
}

/// Where two infinite lines through the segments meet.
fn intersect_lines(
    f1: Point2,
    t1: Point2,
    f2: Point2,
    t2: Point2,
    tol: Tolerances,
) -> OgeomResult<Point2> {
    let d1 = t1 - f1;
    let d2 = t2 - f2;
    let cross = d1.cross(d2);
    if cross.abs() <= tol.angular() * d1.magnitude() * d2.magnitude() {
        ogeom_bail!(Construction, "parallel sides never meet at a corner");
    }
    let t = (f2 - f1).cross(d2) / cross;
    Ok(f1 + d1 * t)
}

/// The crossing of two pieces nearest `corner`.
fn nearest_crossing(a: &Piece, b: &Piece, corner: Point2, tol: Tolerances) -> OgeomResult<Point2> {
    let candidates = crossings(a, b, tol)?;
    candidates
        .into_iter()
        .min_by(|p, q| {
            p.distance(corner)
                .partial_cmp(&q.distance(corner))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .ok_or_else(|| {
            ogeom_core::ogeom_err!(
                Construction,
                "overlapping offset pieces never cross; the offset collapses \
                 here"
            )
        })
}

/// All crossings of the two pieces' supports.
fn crossings(a: &Piece, b: &Piece, tol: Tolerances) -> OgeomResult<Vec<Point2>> {
    let support = |p: &Piece| -> OgeomResult<PlanarCurve> {
        Ok(match p {
            Piece::Seg { from, to } => ogeom_geom::Line2d::segment(*from, *to, tol)?.into(),
            Piece::Arc { centre, radius, .. } => {
                ogeom_geom::Circle2d::new(ogeom_math::Circle2::new(
                    ogeom_math::Frame2::new(*centre, ogeom_math::Direction2::X),
                    *radius,
                    tol,
                )?)
                .into()
            }
        })
    };
    let found = ogeom_intersect::intersect_curves_2d(
        &support(a)?,
        &support(b)?,
        ogeom_intersect::CurveCurveOptions::default(),
        tol,
    )?;
    Ok(found.crossings.into_iter().map(|c| c.point).collect())
}

/// Whether a support crossing lands within the piece's own span, endpoints
/// excluded.
fn within(piece: &Piece, p: Point2, tol: Tolerances) -> bool {
    let margin = tol.confusion() * 100.0;
    match piece {
        Piece::Seg { from, to } => {
            let d = *to - *from;
            let len = d.magnitude();
            let t = (p - *from).dot(d / len);
            t > margin && t < len - margin
        }
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => {
            let a = (p - *centre).y.atan2((p - *centre).x);
            let (lo, hi) = if end > start {
                (*start, *end)
            } else {
                (*end, *start)
            };
            let tau = core::f64::consts::TAU;
            let mut folded = a;
            while folded < lo {
                folded += tau;
            }
            folded > lo + margin / radius && folded < hi - margin / radius
        }
    }
}

/// Trim a piece's end back to `at`.
fn trim_end(piece: &mut Piece, at: Point2, tol: Tolerances) -> OgeomResult<()> {
    match piece {
        Piece::Seg { from, to } => {
            let d = (*to - *from).magnitude();
            let kept = (at - *from).dot((*to - *from) / d);
            if kept <= tol.confusion() {
                ogeom_bail!(Construction, "the trim consumes the offset edge whole");
            }
            *to = at;
        }
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => {
            let a = (at - *centre).y.atan2((at - *centre).x);
            *end = align_angle(a, *start, *end, *radius, tol)?;
        }
    }
    Ok(())
}

/// Trim a piece's start forward to `at`.
fn trim_start(piece: &mut Piece, at: Point2, tol: Tolerances) -> OgeomResult<()> {
    match piece {
        Piece::Seg { from, to } => {
            let d = (*to - *from).magnitude();
            let kept = (*to - at).dot((*to - *from) / d);
            if kept <= tol.confusion() {
                ogeom_bail!(Construction, "the trim consumes the offset edge whole");
            }
            *from = at;
        }
        Piece::Arc {
            centre,
            radius,
            start,
            end,
        } => {
            let a = (at - *centre).y.atan2((at - *centre).x);
            *start = align_angle(a, *end, *start, *radius, tol)?;
        }
    }
    Ok(())
}

/// Fold `angle` next to the span so the sweep from `anchor` keeps its sense
/// and stays non-empty.
fn align_angle(
    angle: f64,
    anchor: f64,
    replaced: f64,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let tau = core::f64::consts::TAU;
    let mut a = angle;
    // Choose the representative nearest the value being replaced.
    while a - replaced > core::f64::consts::PI {
        a -= tau;
    }
    while replaced - a > core::f64::consts::PI {
        a += tau;
    }
    if ((a - anchor).abs() * radius) <= tol.confusion() {
        ogeom_bail!(Construction, "the trim consumes the offset arc whole");
    }
    if (a - anchor).signum() != (replaced - anchor).signum() {
        ogeom_bail!(Construction, "the trim runs past the offset arc's start");
    }
    Ok(a)
}
