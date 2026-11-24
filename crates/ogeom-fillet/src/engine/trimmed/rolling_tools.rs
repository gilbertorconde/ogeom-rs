//! Tools for bevels rolling over or reaching past their faces.

use super::*;

/// Chamfers rolling over their faces that no other route builds, by the
/// one boolean (ADR 0001, amendment "M5b"): each chamfer takes the
/// material between its edge and its bevel, its tool bounded by the bevel,
/// its two faces' surfaces carried on past the faces, and at each end of
/// the edge the plane square to it there; where two bevels meet at a
/// vertex their tools meet there, and the union of the tools is cut (or
/// fused, at a concave edge) in one boolean. A wall rising at an edge's end
/// is cut only where the tool reaches it, and a far face only where a tool
/// reaches it. A rim's tool is its section turned over its own arc. Only
/// where some asked bevel's contact leaves its face (a bevel rolling over
/// or reaching past it); `None` otherwise, or where a tool will not close,
/// a boolean refuses or a bevel is left nothing of.
///
/// # Errors
///
/// Refuses by name the two corners no tool closes: three or more such
/// bevels at one vertex (their corner closes with a facet), and two
/// meeting at a reflex corner of the face they share (each tool would run
/// on past the vertex to the other's). As reading the solid and building
/// the tools, other than their refusals.
pub(crate) fn rolling_tools(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if specs.is_empty() {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut strips: Vec<Strip> = Vec::with_capacity(specs.len());
    for (edge, spec) in specs {
        let e = id(edge);
        if !curved.edges.contains_key(&e) {
            return Ok(None);
        }
        let strip = match bevel_of(model, solid, &curved, e, edge, spec, false, tol)? {
            Some(strip) => strip,
            None => match concave_straight(model, &curved, e, spec, tol)? {
                Some(strip) => strip,
                None => return Ok(None),
            },
        };
        strips.push(strip);
    }
    // Some contact leaves its face: a bevel rolling over or past it. Read
    // on the face's trim drawn finer than the bevel's distances, so a
    // contact past the face's far side by less than the default chord's
    // band reads out, not on.
    let mut rolling = false;
    'strips: for strip in &strips {
        for (k, contact) in strip.contacts.iter().enumerate() {
            let (lo, hi) = contact.domain();
            let edge = &curved.edges[&strip.edge];
            let at = edge.curve.point_at(edge.range.0, tol)?;
            let away = ogeom_algo::project_on_curve(contact, at, 64, tol)?.distance;
            let fine = ogeom_mesh::Deflection {
                chord: (away * 1e-3).max(tol.confusion() * 10.0),
                ..ogeom_mesh::Deflection::default()
            };
            for share in [0.25, 0.5, 0.75] {
                let p = contact.point_at(lo + (hi - lo) * share, tol)?;
                let face = &curved.faces[strip.faces[k]].shape;
                if ogeom_algo::classify_on_face(model, face, p, fine, tol)?
                    == ogeom_algo::Containment::Out
                {
                    rolling = true;
                    break 'strips;
                }
            }
        }
    }
    if !rolling {
        return Ok(None);
    }
    // The corners the tools close: at most two bevels at a vertex, meeting
    // where the face they share turns less than half a turn.
    let asked: Vec<Id> = strips.iter().map(|s| s.edge).collect();
    let mut at_vertex: FastMap<Id, Vec<Id>> = FastMap::default();
    for e in &asked {
        let c = &curved.edges[e];
        at_vertex.entry(c.ends.0).or_default().push(*e);
        if c.ends.1 != c.ends.0 {
            at_vertex.entry(c.ends.1).or_default().push(*e);
        }
    }
    let mut reflex: Vec<(Id, Id, Id)> = Vec::new();
    let mut faceted: Vec<Id> = Vec::new();
    let mut vertices: Vec<&Id> = at_vertex.keys().collect();
    vertices.sort_unstable();
    for v in vertices {
        let meeting = &at_vertex[v];
        match meeting[..] {
            [_] => {}
            [a, b] => {
                if reflex_between(model, &curved, *v, (a, b), tol)? {
                    reflex.push((*v, a, b));
                }
            }
            _ => faceted.push(*v),
        }
    }
    let mut tools = Vec::with_capacity(strips.len());
    for ((edge, _), strip) in specs.iter().zip(&strips) {
        let c = &curved.edges[&strip.edge];
        if let Curve::Circle(_) = &c.curve
            && let Some((section, sweep)) = rim_section_of(&curved, strip, tol)?
        {
            let Some((solid, named, sense)) = revolved_tool(model, &curved, section, sweep, tol)?
            else {
                return Ok(None);
            };
            tools.push(super::fallback::Tool {
                solid,
                sense,
                faces: vec![(named, edge.clone())],
            });
            continue;
        }
        if c.ends.0 == c.ends.1 {
            return Ok(None);
        }
        let mut caps = [EndCap::Own, EndCap::Own];
        for (k, t) in [c.range.0, c.range.1].into_iter().enumerate() {
            let h = (c.range.1 - c.range.0) * 1e-6;
            let p = c.curve.point_at(t, tol)?;
            let q = c.curve.point_at(if k == 0 { t + h } else { t - h }, tol)?;
            caps[k] = EndCap::Plane {
                at: p,
                normal: q - p,
            };
        }
        let Some(tool) = close_tool(model, &curved, strip, caps, tol)? else {
            return Ok(None);
        };
        tools.push(tool);
    }
    // At a reflex corner each tool would run on past the vertex until it
    // met the other's: the corner piece both would take there, beyond both
    // caps, outside both bevels and inside every seat face.
    for (v, a, b) in reflex {
        let by_edge = |e: Id| strips.iter().find(|s| s.edge == e);
        let (Some(sa), Some(sb)) = (by_edge(a), by_edge(b)) else {
            return Ok(None);
        };
        // Both edges convex: their tools cut; both concave: fused.
        let sense_of = |e: Id| {
            strips
                .iter()
                .position(|s| s.edge == e)
                .map(|k| tools[k].sense)
        };
        let (Some(sense), Some(other)) = (sense_of(a), sense_of(b)) else {
            return Ok(None);
        };
        let mut piece = None;
        if sense == other {
            if sense == super::fallback::Sense::Cut {
                piece = reflex_corner(model, &curved, v, [sa, sb], tol)?.map(Some);
            }
            if piece.is_none() {
                piece = reflex_piece(model, &curved, v, [sa, sb], tol)?;
            }
        }
        let Some(piece) = piece else {
            ogeom_core::ogeom_bail!(
                Construction,
                "two chamfers rolling over their faces meet at a reflex corner, one \
                 convex and one concave, or where the piece their tools meet in is not \
                 bounded near the vertex; no tool closes that corner"
            );
        };
        if let Some(piece) = piece {
            tools.push(super::fallback::Tool {
                solid: piece,
                sense,
                faces: Vec::new(),
            });
        }
    }
    // Three or more at a vertex close with the facet through their setback
    // corners: the piece inside every bevel on the vertex's side of it. Of
    // mixed convexity (ADR 0001, amendment "mixed rolling corners"), the
    // vertex's side of the facet takes the sense most of the vertex's
    // bevels take (the complement's corner is the complement), and the
    // bevels' tools spare what the facet settles: the cuts after the fuses
    // neither undo the piece nor keep what the facet gives up.
    for v in faceted {
        let meeting = &at_vertex[&v];
        let mut cuts = 0_usize;
        let mut fuses = 0_usize;
        let mut own: Vec<usize> = Vec::new();
        for e in meeting {
            let Some(k) = strips.iter().position(|s| s.edge == *e) else {
                return Ok(None);
            };
            own.push(k);
            match tools[k].sense {
                super::fallback::Sense::Cut => cuts += 1,
                super::fallback::Sense::Fuse => fuses += 1,
            }
        }
        if cuts == fuses {
            ogeom_core::ogeom_bail!(
                Construction,
                "chamfers rolling over their faces meet at one vertex, as many convex as \
                 concave; which side of their corner's facet the material keeps is not \
                 settled by the bevels"
            );
        }
        let sense = if cuts > fuses {
            super::fallback::Sense::Cut
        } else {
            super::fallback::Sense::Fuse
        };
        let at_v: Vec<&Strip> = strips
            .iter()
            .filter(|s| meeting.contains(&s.edge))
            .collect();
        let mixed = cuts > 0 && fuses > 0;
        // Of mixed convexity the facet's side runs on into the solid (or
        // out of it) past the bevels: what changes is what of it lies
        // outside the solid where the piece is fused, inside where it is
        // cut.
        let within = match sense {
            super::fallback::Sense::Fuse => Side::Drop(solid.clone()),
            super::fallback::Sense::Cut => Side::Keep(solid.clone(), None),
        };
        let Some(piece) = facet_piece(
            model,
            &curved,
            v,
            &at_v,
            mixed.then_some(&within),
            true,
            false,
            tol,
        )?
        else {
            ogeom_core::ogeom_bail!(
                Construction,
                "three or more chamfers rolling over their faces meet at one vertex, \
                 and the facet through their setback corners does not close it"
            );
        };
        if let Some(piece) = piece {
            // The facet's two sides inside every bevel, as the box clips
            // them: the vertex's (what the piece's sense keeps there) spared
            // by the tools of the other sense, the far one (what it gives
            // up) by the tools of the piece's sense.
            let cells = if mixed {
                let mut both = Vec::with_capacity(2);
                for far in [false, true] {
                    match facet_piece(model, &curved, v, &at_v, None, false, far, tol)? {
                        Some(cell) => both.push(cell),
                        None => return Ok(None),
                    }
                }
                Some(both)
            } else {
                None
            };
            for &k in &own {
                let Some(cells) = &cells else {
                    continue;
                };
                let Some(cell) = &cells[usize::from(tools[k].sense == sense)] else {
                    continue;
                };
                let Ok(spared) = ogeom_bool::cut(model, &tools[k].solid, cell, tol) else {
                    return Ok(None);
                };
                tools[k].faces = tools[k]
                    .faces
                    .iter()
                    .flat_map(|(face, from)| {
                        spared
                            .history
                            .trace(face)
                            .iter()
                            .map(|f| (f.clone(), from.clone()))
                            .collect::<Vec<_>>()
                    })
                    .collect();
                tools[k].solid = spared.shape;
            }
            tools.push(super::fallback::Tool {
                solid: piece,
                sense,
                faces: Vec::new(),
            });
        }
    }
    // The corner pieces first: each bevel's tool then meets them across a
    // face, not only along the line where two tools touch at a corner.
    let pieces = tools.split_off(strips.len());
    tools.splice(0..0, pieces);
    let edges: Vec<Shape> = specs.iter().map(|(e, _)| e.clone()).collect();
    super::fallback::applied(model, solid, &edges, &tools, tol)
}

/// A solid on one side of a surface over a box round a corner, and whether
/// the corner keeps what lies in it (else what lies outside it). A kept
/// cone's inside stops short of its apex: `Keep`'s bound is the axis
/// (a point on it and its direction) and the height along it the piece
/// must stay short of.
#[derive(Clone)]
enum Side {
    Keep(Shape, Option<(Point, Vector, f64)>),
    Drop(Shape),
}

/// The side of the plane through `through` of normal `n` holding `keep`,
/// as a block reaching well past the box `(centre, half)`. `None` where
/// `keep` lies on the plane.
fn plane_side(
    model: &mut Model,
    (n, through): (Vector, Point),
    keep: Point,
    (centre, half): (Point, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Side>> {
    let n = n / n.magnitude();
    let off = (keep - through).dot(n);
    if off.abs() <= tol.confusion() {
        return Ok(None);
    }
    let up = n * off.signum();
    let reach = 4.0 * half;
    let foot = centre - n * (centre - through).dot(n);
    let z = Direction::new(up, tol)?;
    let x = Direction::new(
        if up.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        },
        tol,
    )?;
    let frame = Frame::new(foot, z, x, tol)?;
    let origin = foot - frame.x().vector() * reach - frame.y().vector() * reach;
    let frame = Frame::new(origin, z, frame.x(), tol)?;
    let block = ogeom_algo::make_box(model, frame, (2.0 * reach, 2.0 * reach, reach), tol)?.shape;
    Ok(Some(Side::Keep(block, None)))
}

/// The side of `surface` holding `keep`, over the box `(centre, half)`: a
/// plane's half-space, or the inside or outside of a cylinder, cone or
/// sphere. `None` for another surface, a cone whose apex falls by the box,
/// or `keep` on the surface.
fn surface_side(
    model: &mut Model,
    surface: &SurfaceGeometry,
    keep: Point,
    (centre, half): (Point, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Side>> {
    let eps = tol.confusion();
    let revolved = |model: &mut Model,
                    frame: Frame,
                    radius: &dyn Fn(f64) -> f64|
     -> OgeomResult<Option<Side>> {
        let (o, a) = (frame.origin(), frame.z().vector());
        let mut heights = Vec::with_capacity(8);
        for k in 0..8 {
            let corner = centre
                + Vector::new(
                    if k & 1 == 0 { -half } else { half },
                    if k & 2 == 0 { -half } else { half },
                    if k & 4 == 0 { -half } else { half },
                );
            heights.push((corner - o).dot(a));
        }
        let mut low = heights.iter().copied().fold(f64::INFINITY, f64::min) - 0.01 * half;
        let mut high = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max) + 0.01 * half;
        let q = keep - o;
        let h = q.dot(a);
        let rho = (q - a * h).magnitude();
        if (rho - radius(h)).abs() <= eps || radius(h) <= eps {
            return Ok(None);
        }
        // A cone's apex among the heights: the nappe holding `keep` only,
        // stopped where its radius is a thousandth of the box's.
        let mut bound = None;
        let slope = radius(h + 1.0) - radius(h);
        if slope.abs() > eps {
            let apex = h - radius(h) / slope;
            let stop = apex + (1e-3 * half) / slope;
            if apex > low && apex < high {
                if apex > h {
                    high = stop.min(high);
                    bound = Some((o, a, high));
                } else {
                    low = stop.max(low);
                    bound = Some((o, -a, -low));
                }
            }
        }
        let (r_low, r_high) = (radius(low), radius(high));
        if !(r_low > eps && r_high > eps && high > low) {
            return Ok(None);
        }
        let base = Frame::new(o + a * low, frame.z(), frame.x(), tol)?;
        let solid = if (r_low - r_high).abs() <= eps {
            ogeom_algo::make_cylinder(model, base, r_low, high - low, tol)?.shape
        } else {
            ogeom_algo::make_cone(model, base, r_low, r_high, high - low, tol)?.shape
        };
        Ok(Some(if rho < radius(h) {
            Side::Keep(solid, bound)
        } else {
            Side::Drop(solid)
        }))
    };
    match surface {
        SurfaceGeometry::Plane(p) => {
            let plane = p.plane();
            plane_side(
                model,
                (plane.normal().vector(), plane.origin()),
                keep,
                (centre, half),
                tol,
            )
        }
        SurfaceGeometry::Cylinder(c) => {
            let r = c.cylinder().radius();
            revolved(model, c.cylinder().frame(), &|_| r)
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            revolved(model, cone.frame(), &|z| cone.radius_at(z))
        }
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let d = keep.distance(sphere.centre());
            if (d - sphere.radius()).abs() <= eps {
                return Ok(None);
            }
            let frame = Frame::new(sphere.centre(), Direction::Z, Direction::X, tol)?;
            let ball = ogeom_algo::make_sphere(model, frame, sphere.radius(), tol)?.shape;
            Ok(Some(if d < sphere.radius() {
                Side::Keep(ball, None)
            } else {
                Side::Drop(ball)
            }))
        }
        _ => Ok(None),
    }
}

/// The box `(centre, half)` taken in common with each side kept and less
/// each side dropped, and where `near` is given only the lumps of it that
/// reach that point: `Some(None)` where nothing is left, `None` where a
/// boolean refuses or what is kept reaches the box's boundary (the piece
/// is not bounded near the corner).
fn corner_piece(
    model: &mut Model,
    (centre, half): (Point, f64),
    sides: Vec<Side>,
    near: Option<Point>,
    tol: Tolerances,
) -> OgeomResult<Option<Option<Shape>>> {
    corner_piece_bounded(model, (centre, half), sides, true, near, None, tol)
}

/// [`corner_piece`], its piece held bounded near the corner only where
/// `bounded` asks (else the piece as the box clips it, for a tool to be
/// spared what the corner's facet keeps), only the lumps reaching
/// `reaching` where given (as [`corner_piece`] keeps them), and of a piece
/// in lumps only the lump nearest `nearest_to`, where given.
fn corner_piece_bounded(
    model: &mut Model,
    (centre, half): (Point, f64),
    sides: Vec<Side>,
    bounded: bool,
    reaching: Option<Point>,
    nearest_to: Option<Point>,
    tol: Tolerances,
) -> OgeomResult<Option<Option<Shape>>> {
    let low = centre - Vector::new(half, half, half);
    let frame = Frame::new(low, Direction::Z, Direction::X, tol)?;
    let mut piece =
        ogeom_algo::make_box(model, frame, (2.0 * half, 2.0 * half, 2.0 * half), tol)?.shape;
    for side in &sides {
        let step = match side {
            Side::Keep(s, _) => ogeom_bool::common(model, &piece, s, tol),
            Side::Drop(s) => ogeom_bool::cut(model, &piece, s, tol),
        };
        let Ok(step) = step else {
            return Ok(None);
        };
        piece = step.shape;
        if ogeom_topo::explore_unique(model, &piece, ogeom_topo::ShapeType::Face)?.is_empty() {
            return Ok(Some(None));
        }
    }
    if !bounded {
        return Ok(Some(Some(piece)));
    }
    // The lumps reaching `reaching` alone: the region past a bevel carried on
    // can part round the other bevel, and only the corner's own lump is
    // the corner's.
    if let Some(at) = reaching {
        let shells = ogeom_topo::explore_unique(model, &piece, ogeom_topo::ShapeType::Shell)?;
        let reach = 1e-6 * half + tol.confusion() * 10.0;
        let mut kept: Vec<Shape> = Vec::new();
        for shell in &shells {
            let b = ogeom_algo::shape_bounds(model, shell, tol)?;
            let (Some(lo), Some(hi)) = (b.low(), b.high()) else {
                continue;
            };
            if at.x >= lo.x - reach
                && at.y >= lo.y - reach
                && at.z >= lo.z - reach
                && at.x <= hi.x + reach
                && at.y <= hi.y + reach
                && at.z <= hi.z + reach
            {
                kept.push(shell.clone());
            }
        }
        if kept.is_empty() {
            return Ok(Some(None));
        }
        if kept.len() < shells.len() {
            let mut lumps = Vec::with_capacity(kept.len());
            for shell in kept {
                lumps.push(ogeom_algo::make_solid(model, &[shell])?.shape);
            }
            let mut joined = lumps[0].clone();
            for lump in &lumps[1..] {
                let Ok(step) = ogeom_bool::fuse(model, &joined, lump, tol) else {
                    return Ok(None);
                };
                joined = step.shape;
            }
            piece = joined;
        }
    }
    if let Some(at) = nearest_to {
        let lumps = ogeom_topo::explore_unique(model, &piece, ogeom_topo::ShapeType::Solid)?;
        if lumps.len() > 1 {
            let mut nearest: Option<(f64, Shape)> = None;
            for lump in lumps {
                let b = ogeom_algo::shape_bounds(model, &lump, tol)?;
                let (Some(l), Some(h)) = (b.low(), b.high()) else {
                    continue;
                };
                let clamp = Point::new(
                    at.x.clamp(l.x, h.x),
                    at.y.clamp(l.y, h.y),
                    at.z.clamp(l.z, h.z),
                );
                let d = clamp.distance(at);
                if nearest.as_ref().is_none_or(|(best, _)| d < *best) {
                    nearest = Some((d, lump));
                }
            }
            let Some((_, lump)) = nearest else {
                return Ok(Some(None));
            };
            piece = lump;
        }
    }
    let bounds = ogeom_algo::shape_bounds(model, &piece, tol)?;
    let (Some(lo), Some(hi)) = (bounds.low(), bounds.high()) else {
        return Ok(Some(None));
    };
    let margin = 1e-6 * half;
    let high = centre + Vector::new(half, half, half);
    if lo.x <= low.x + margin
        || lo.y <= low.y + margin
        || lo.z <= low.z + margin
        || hi.x >= high.x - margin
        || hi.y >= high.y - margin
        || hi.z >= high.z - margin
    {
        return Ok(None);
    }
    // Short of each kept cone's stop: else the stop, not the cone, bounds it.
    for side in &sides {
        if let Side::Keep(_, Some((o, a, stop))) = side {
            for v in ogeom_topo::explore_unique(model, &piece, ogeom_topo::ShapeType::Vertex)? {
                let Some(p) = model
                    .node(&v)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| d.point)
                else {
                    continue;
                };
                if (p - *o).dot(*a) >= stop - margin {
                    return Ok(None);
                }
            }
        }
    }
    Ok(Some(Some(piece)))
}

/// `p` reflected across the plane through `through` of normal `n`.
fn reflected(p: Point, (n, through): (Vector, Point)) -> Point {
    let unit = n / n.magnitude();
    p - unit * (2.0 * (p - through).dot(unit))
}

/// How far a strip's contacts stand from `v`, at most.
fn contact_reach(strip: &Strip, at: Point, tol: Tolerances) -> OgeomResult<f64> {
    let mut far = 0.0_f64;
    for contact in &strip.contacts {
        far = far.max(ogeom_algo::project_on_curve(contact, at, 64, tol)?.distance);
    }
    Ok(far)
}

/// The piece two tools meeting at the reflex corner `v` take past both
/// their caps, for any bevels and faces the sides are read on: beyond each
/// edge's plane square to it at `v`, on the edge's side of each bevel, and
/// on each seat face's side that holds the strip's other contact. `None`
/// where a side is not read, or the piece is not bounded near `v`;
/// `Some(None)` where it is empty.
fn reflex_piece(
    model: &mut Model,
    curved: &Curved,
    v: Id,
    strips: [&Strip; 2],
    tol: Tolerances,
) -> OgeomResult<Option<Option<Shape>>> {
    let at = curved.points[&v];
    let mut reach = 0.0_f64;
    for strip in strips {
        reach = reach.max(contact_reach(strip, at, tol)?);
    }
    let reach = reach.max(tol.confusion() * 1e3);
    // The box grown until the piece stands clear of it.
    for scale in [4.0, 16.0, 64.0] {
        if let Some(piece) = reflex_piece_in(model, curved, v, strips, (at, scale * reach), tol)? {
            return Ok(Some(piece));
        }
    }
    Ok(None)
}

/// [`reflex_piece`] over the box `bx` round the vertex.
fn reflex_piece_in(
    model: &mut Model,
    curved: &Curved,
    v: Id,
    strips: [&Strip; 2],
    bx: (Point, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Option<Shape>>> {
    let (at, half) = bx;
    let reach = half / 4.0;
    let mut sides = Vec::new();
    // The faces both strips stand on (or two on one carrier, a face split
    // along seams).
    let mut shared: Vec<usize> = Vec::new();
    for fa in strips[0].faces {
        for fb in strips[1].faces {
            if fa == fb
                || ogeom_heal::same_domain(
                    model,
                    &curved.faces[fa].shape,
                    &curved.faces[fb].shape,
                    tol,
                )?
            {
                shared.extend([fa, fb]);
            }
        }
    }
    // The half-spaces taken so far, each `n . p <= off` with `n` a unit
    // normal: a plane met twice (a cap on another edge's face) is taken
    // once, as two booleans on one plane would leave nothing.
    let mut seen: Vec<(Vector, f64)> = Vec::new();
    let mut half_space = |n: Vector, through: Point, keep: Point| -> Option<(Vector, f64)> {
        let n = n / n.magnitude();
        let n = if (keep - through).dot(n) > 0.0 { -n } else { n };
        let off = n.dot(through.to_vector());
        if seen
            .iter()
            .any(|(m, o)| m.dot(n) >= 1.0 - tol.angular() && (o - off).abs() <= tol.confusion())
        {
            return None;
        }
        seen.push((n, off));
        Some((n, off))
    };
    for strip in strips {
        let c = &curved.edges[&strip.edge];
        let away = {
            let (t, sign) = if c.ends.0 == v {
                (c.range.0, 1.0)
            } else {
                (c.range.1, -1.0)
            };
            let d = c.curve.d1_at(t, tol)? * sign;
            d / d.magnitude()
        };
        let beyond = at - away * reach;
        if half_space(away, at, beyond).is_some() {
            let Some(cap) = plane_side(model, (away, at), beyond, bx, tol)? else {
                return Ok(None);
            };
            sides.push(cap);
        }
        let Some(bevel) = surface_side(model, &strip.surface, at, bx, tol)? else {
            return Ok(None);
        };
        sides.push(bevel);
        for k in 0..2 {
            // Only the face both edges bound: each bevel carried on past
            // its own other face, which would stop the piece short of the
            // neighbouring bevel (a sliver beside a curved face's edge).
            let face = strip.faces[k];
            if !shared.contains(&face) {
                continue;
            }
            let other = ogeom_algo::project_on_curve(&strip.contacts[1 - k], at, 64, tol)?.point;
            let surface = &curved.faces[strip.faces[k]].surface;
            if let SurfaceGeometry::Plane(p) = surface {
                let plane = p.plane();
                if half_space(plane.normal().vector(), plane.origin(), other).is_none() {
                    continue;
                }
            }
            let Some(face) = surface_side(model, surface, other, bx, tol)? else {
                return Ok(None);
            };
            sides.push(face);
        }
    }
    corner_piece(model, bx, sides, Some(at), tol)
}

/// The facet piece at the vertex `v` of three or more rolling bevels: the
/// facet through the setback corners (each face's corner where the curves
/// bounding it at `v` meet: an asked edge's contact on it, else the edge),
/// a plane where they are coplanar, else a fan of triangles about their
/// centroid (a G0 corner, ADR 0001, "chamfer corners are G0"); the piece is
/// what lies on the vertex's side of it and inside every bevel. `None`
/// where a corner or a side is not read or the piece is not bounded near
/// `v`; `Some(None)` where it is empty.
#[allow(
    clippy::too_many_arguments,
    reason = "one vertex's corner and the three ways it is read"
)]
fn facet_piece(
    model: &mut Model,
    curved: &Curved,
    v: Id,
    strips: &[&Strip],
    within: Option<&Side>,
    bounded: bool,
    far: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Option<Shape>>> {
    let at = curved.points[&v];
    let mut reach = 0.0_f64;
    for strip in strips {
        reach = reach.max(contact_reach(strip, at, tol)?);
    }
    let long = 100.0 * reach.max(1.0);
    // The curve bounding face `f` at `v` along edge `e`: its contact there
    // where asked, else the edge; lines carried on far either way.
    let bounding = |f: usize, e: Id| -> OgeomResult<Option<Curve>> {
        let curve = match strips.iter().find(|s| s.edge == e) {
            Some(s) => match s.faces.iter().position(|x| *x == f) {
                Some(k) => s.contacts[k].clone(),
                None => return Ok(None),
            },
            None => curved.edges[&e].curve.clone(),
        };
        Ok(Some(match &curve {
            Curve::Line(line) => {
                let (lo, hi) = curve.domain();
                let p = curve.point_at(lo, tol)?;
                let q = curve.point_at(hi, tol)?;
                let d = (q - p) / p.distance(q);
                let _ = line;
                Curve::Line(LineCurve::segment(p - d * long, p + d * long, tol)?)
            }
            _ => curve,
        }))
    };
    // Each face round `v` and its two edges there.
    let mut corners: Vec<Point> = Vec::new();
    for (fi, face) in curved.faces.iter().enumerate() {
        let at_v: Vec<Id> = face
            .loops
            .iter()
            .flatten()
            .map(|(e, _)| *e)
            .filter(|e| {
                let c = &curved.edges[e];
                c.ends.0 == v || c.ends.1 == v
            })
            .collect();
        if at_v.is_empty() {
            continue;
        }
        let [e1, e2] = at_v[..] else {
            return Ok(None);
        };
        if !strips.iter().any(|s| s.edge == e1 || s.edge == e2) {
            return Ok(None);
        }
        let (Some(c1), Some(c2)) = (bounding(fi, e1)?, bounding(fi, e2)?) else {
            return Ok(None);
        };
        let found = ogeom_intersect::curves::intersect_curves(
            &c1,
            &c2,
            ogeom_intersect::curves::CurveCurveOptions {
                gap: 1e3 * tol.confusion(),
                ..ogeom_intersect::curves::CurveCurveOptions::default()
            },
            tol,
        )?;
        let Some(nearest) = found
            .crossings
            .iter()
            .min_by(|a, b| a.point.distance(at).total_cmp(&b.point.distance(at)))
        else {
            return Ok(None);
        };
        corners.push(nearest.point);
    }
    if corners.len() < 3 {
        return Ok(None);
    }
    #[expect(clippy::cast_precision_loss, reason = "a vertex's few corners")]
    let n = corners.len() as f64;
    let centroid = Point::new(
        corners.iter().map(|p| p.x).sum::<f64>() / n,
        corners.iter().map(|p| p.y).sum::<f64>() / n,
        corners.iter().map(|p| p.z).sum::<f64>() / n,
    );
    let half = 4.0 * corners.iter().map(|p| p.distance(at)).fold(reach, f64::max);
    let bx = (at, half);
    // Inside every bevel: on the side of each holding the centroid.
    let bevels = |model: &mut Model| -> OgeomResult<Option<Vec<Side>>> {
        let mut sides = Vec::new();
        for strip in strips {
            let Some(side) = surface_side(model, &strip.surface, centroid, bx, tol)? else {
                return Ok(None);
            };
            sides.push(side);
        }
        if let Some(side) = within {
            sides.push(side.clone());
        }
        Ok(Some(sides))
    };
    // The corners in order round the vertex, by their angle about the
    // centroid's direction.
    let axis = centroid - at;
    if axis.magnitude() <= tol.confusion() {
        return Ok(None);
    }
    let axis = axis / axis.magnitude();
    let x = {
        let r = corners[0] - centroid;
        let r = r - axis * r.dot(axis);
        r / r.magnitude()
    };
    let y = axis.cross(x);
    let mut order: Vec<Point> = corners.clone();
    order.sort_by(|a, b| {
        let angle = |p: &Point| {
            let r = *p - centroid;
            r.dot(y).atan2(r.dot(x))
        };
        angle(a).total_cmp(&angle(b))
    });
    // Coplanar corners: one facet.
    let normal = (order[1] - order[0]).cross(order[2] - order[0]);
    if normal.magnitude() <= tol.confusion() {
        return Ok(None);
    }
    let flat = order
        .iter()
        .all(|p| (*p - order[0]).dot(normal / normal.magnitude()).abs() <= tol.confusion());
    if flat {
        let Some(mut sides) = bevels(model)? else {
            return Ok(None);
        };
        let keep = if far {
            reflected(at, (normal, order[0]))
        } else {
            at
        };
        let Some(facet) = plane_side(model, (normal, order[0]), keep, bx, tol)? else {
            return Ok(None);
        };
        sides.push(facet);
        return corner_piece_bounded(model, bx, sides, bounded, None, within.map(|_| at), tol);
    }
    // A fan of triangles about the centroid, each over its sector.
    let mut pieces: Vec<Shape> = Vec::new();
    for i in 0..order.len() {
        let (a, b) = (order[i], order[(i + 1) % order.len()]);
        let Some(mut sides) = bevels(model)? else {
            return Ok(None);
        };
        let tri = (a - centroid).cross(b - centroid);
        let (Some(facet), Some(from), Some(to)) = (
            plane_side(
                model,
                (tri, centroid),
                if far {
                    reflected(at, (tri, centroid))
                } else {
                    at
                },
                bx,
                tol,
            )?,
            plane_side(model, ((a - at).cross(centroid - at), at), b, bx, tol)?,
            plane_side(model, ((b - at).cross(centroid - at), at), a, bx, tol)?,
        ) else {
            return Ok(None);
        };
        sides.extend([facet, from, to]);
        match corner_piece_bounded(model, bx, sides, bounded, None, within.map(|_| at), tol)? {
            None => return Ok(None),
            Some(None) => {}
            Some(Some(piece)) => pieces.push(piece),
        }
    }
    let Some(mut joined) = pieces.first().cloned() else {
        return Ok(Some(None));
    };
    for piece in &pieces[1..] {
        let Ok(step) = ogeom_bool::fuse(model, &joined, piece, tol) else {
            return Ok(None);
        };
        joined = step.shape;
    }
    Ok(Some(Some(joined)))
}

/// The piece two tools meeting at the reflex corner `v` take past both
/// their caps: the polytope beyond each edge's plane square to it at `v`,
/// outside each bevel's plane and inside each seat face's plane. `None`
/// where a bevel or a seat face is not a plane, an edge is not straight,
/// or the piece is not bounded.
fn reflex_corner(
    model: &mut Model,
    curved: &Curved,
    v: Id,
    strips: [&Strip; 2],
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let at = curved.points[&v];
    let mut planes: Vec<(Vector, f64)> = Vec::new();
    // `n . p <= n . through`.
    let mut keep = |n: Vector, through: Point| {
        let n = n / n.magnitude();
        planes.push((n, n.dot(through.to_vector())));
    };
    for strip in strips {
        let c = &curved.edges[&strip.edge];
        let Curve::Line(_) = &c.curve else {
            return Ok(None);
        };
        let SurfaceGeometry::Plane(bevel) = &strip.surface else {
            return Ok(None);
        };
        let away = {
            let (t, sign) = if c.ends.0 == v {
                (c.range.0, 1.0)
            } else {
                (c.range.1, -1.0)
            };
            c.curve.d1_at(t, tol)? * sign
        };
        // Beyond the cap: back from the vertex along the edge.
        keep(away, at);
        for fi in strip.faces {
            let SurfaceGeometry::Plane(_) = &curved.faces[fi].surface else {
                return Ok(None);
            };
            let Some(n) = outward(curved, fi, at) else {
                return Ok(None);
            };
            keep(n, at);
        }
        let plane = bevel.plane();
        let m = plane.normal().vector();
        let side = (at - plane.origin()).dot(m);
        if side.abs() <= tol.confusion() {
            return Ok(None);
        }
        // Outside the bevel: on the edge's side of its plane.
        keep(m * -side.signum(), plane.origin());
    }
    super::fallback::polytope(model, &planes, tol)
}

/// A concave straight edge's strip between two planes, its legs as asked:
/// each contact the edge moved its leg into its face (the side read on the
/// face itself, [`into_face`]), the length of the edge and half again past
/// each end; the bevel the plane through both. `None` for any other edge,
/// a convex one, or a spelling not read.
fn concave_straight(
    model: &Model,
    curved: &Curved,
    edge: Id,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let e = &curved.edges[&edge];
    let Curve::Line(_) = &e.curve else {
        return Ok(None);
    };
    let [f0, f1] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    let faces = [f0, f1];
    if faces
        .iter()
        .any(|f| !matches!(curved.faces[*f].surface, SurfaceGeometry::Plane(_)))
    {
        return Ok(None);
    }
    let Some(d) = distances(curved, faces, spec) else {
        return Ok(None);
    };
    if d.iter().any(|x| !(x.is_finite() && *x > tol.confusion())) {
        return Ok(None);
    }
    let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    let length = a.distance(b);
    if length <= tol.confusion() {
        return Ok(None);
    }
    let u = (b - a) / length;
    let run = u * (length / 2.0);
    let tm = f64::midpoint(e.range.0, e.range.1);
    let mid = e.curve.point_at(tm, tol)?;
    let mut lines: Vec<(Point, Point)> = Vec::with_capacity(2);
    for k in 0..2 {
        let (Some(into), Some(other)) = (
            into_face(model, curved, (faces[k], edge), (mid, tm), tol)?,
            outward(curved, faces[1 - k], mid),
        ) else {
            return Ok(None);
        };
        // A concave edge: into each face goes over the other.
        if into.dot(other) <= 0.0 {
            return Ok(None);
        }
        let shift = into * d[k];
        lines.push((a + shift - run, b + shift + run));
    }
    let contacts = [
        Curve::from(LineCurve::segment(lines[0].0, lines[0].1, tol)?),
        Curve::from(LineCurve::segment(lines[1].0, lines[1].1, tol)?),
    ];
    let across = lines[1].0 - lines[0].0;
    let n = u.cross(across);
    if n.magnitude() <= tol.confusion() {
        return Ok(None);
    }
    let frame = Frame::new(
        lines[0].0,
        Direction::new(n, tol)?,
        Direction::new(u, tol)?,
        tol,
    )?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: PlaneSurface::new(Plane::new(frame)).into(),
        section: Section::Line,
    }))
}

/// Whether the asked edges `a` and `b` meet at `v` at a reflex corner of a
/// face they share (or of two faces on one carrier, as a boolean leaves a
/// face split along seams): `b`, leaving `v`, runs out of the side of `a`
/// its face lies on. Edges sharing no carrier, or running on into each
/// other, are not reflex.
fn reflex_between(
    model: &Model,
    curved: &Curved,
    v: Id,
    (a, b): (Id, Id),
    tol: Tolerances,
) -> OgeomResult<bool> {
    // An edge's direction away from `v`.
    let away = |e: Id| -> OgeomResult<Vector> {
        let c = &curved.edges[&e];
        let (t, sign) = if c.ends.0 == v {
            (c.range.0, 1.0)
        } else {
            (c.range.1, -1.0)
        };
        let d = c.curve.d1_at(t, tol)? * sign;
        Ok(d / d.magnitude())
    };
    let (da, db) = (away(a)?, away(b)?);
    for fa in faces_of(curved, a) {
        let mut shared = false;
        for fb in faces_of(curved, b) {
            if fa == fb
                || ogeom_heal::same_domain(
                    model,
                    &curved.faces[fa].shape,
                    &curved.faces[fb].shape,
                    tol,
                )?
            {
                shared = true;
            }
        }
        if !shared {
            continue;
        }
        // Into the face from `a`, read a tenth of the way along it from `v`.
        let c = &curved.edges[&a];
        let t = if c.ends.0 == v {
            c.range.0 + (c.range.1 - c.range.0) * 0.1
        } else {
            c.range.1 - (c.range.1 - c.range.0) * 0.1
        };
        let near = c.curve.point_at(t, tol)?;
        let Some(into) = into_face(model, curved, (fa, a), (near, t), tol)? else {
            continue;
        };
        if da.cross(db).magnitude() <= tol.angular() {
            return Ok(false);
        }
        return Ok(db.dot(into) < -tol.angular());
    }
    Ok(false)
}

/// The section [`revolved_tool`] turns for a rim's strip, and the rim's
/// own sweep: `None` where its contacts are not circles about its axis.
fn rim_section_of(
    curved: &Curved,
    strip: &Strip,
    tol: Tolerances,
) -> OgeomResult<Option<(RimSection, f64)>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let c = &curved.edges[&strip.edge];
    let Curve::Circle(circle) = &c.curve else {
        return Ok(None);
    };
    let circle = circle.circle();
    let (origin, axis) = (circle.centre(), circle.frame().z().vector());
    let start = c.curve.point_at(c.range.0, tol)?;
    let q = start - origin;
    let radial = q - axis * q.dot(axis);
    let u = radial / radial.magnitude();
    let mut contacts = [Point::ORIGIN; 2];
    for (k, contact) in strip.contacts.iter().enumerate() {
        let Curve::Circle(ring) = contact else {
            return Ok(None);
        };
        let ring = ring.circle();
        let off = ring.centre() - origin;
        if ring.frame().z().vector().cross(axis).magnitude() > ang
            || (off - axis * off.dot(axis)).magnitude() > eps
        {
            return Ok(None);
        }
        contacts[k] = ring.centre() + u * ring.radius();
    }
    let sweep = if c.ends.0 == c.ends.1 {
        core::f64::consts::TAU
    } else {
        c.range.1 - c.range.0
    };
    Ok(Some((
        RimSection {
            origin,
            axis,
            radial: u,
            rim: start,
            faces: strip.faces,
            contacts,
        },
        sweep,
    )))
}

/// Whether another asked edge is a piece of `edge`'s circle (pieces of
/// one rim asked together; alike, they close the turn as one tool, and
/// with different chamfers they are asked apart).
pub(super) fn asked_apart(
    curved: &Curved,
    specs: &[(Shape, Chamfer)],
    edge: Id,
    tol: Tolerances,
) -> bool {
    let Curve::Circle(own) = &curved.edges[&edge].curve else {
        return false;
    };
    let own = own.circle();
    specs.iter().any(|(x, spec)| {
        let x = id(x);
        x != edge
            && curved.edges.get(&x).is_some_and(|piece| {
                matches!(&piece.curve, Curve::Circle(c)
                    if c.circle().centre().distance(own.centre()) <= tol.confusion()
                        && (c.circle().radius() - own.radius()).abs() <= tol.confusion())
            })
            && specs
                .iter()
                .find(|(y, _)| id(y) == edge)
                .is_some_and(|(_, mine)| !crate::chamfer::same_spelling(mine, spec, tol))
    })
}

/// Whether the pieces of `edge`'s circle in `curved`, each between faces
/// on the carriers either side of `edge`, turn once round it together.
pub(super) fn pieces_close(
    model: &Model,
    curved: &Curved,
    edge: Id,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let Curve::Circle(own) = &curved.edges[&edge].curve else {
        return Ok(false);
    };
    let own = own.circle();
    let [a0, a1] = faces_of(curved, edge)[..] else {
        return Ok(false);
    };
    let one = |a: usize, b: usize| -> OgeomResult<bool> {
        ogeom_heal::same_domain(model, &curved.faces[a].shape, &curved.faces[b].shape, tol)
    };
    let mut turned = 0.0;
    for (x, piece) in &curved.edges {
        let Curve::Circle(c) = &piece.curve else {
            continue;
        };
        let c = c.circle();
        if c.centre().distance(own.centre()) > eps
            || (c.radius() - own.radius()).abs() > eps
            || c.frame()
                .z()
                .vector()
                .cross(own.frame().z().vector())
                .magnitude()
                > ang
        {
            continue;
        }
        let [b0, b1] = faces_of(curved, *x)[..] else {
            continue;
        };
        if (one(a0, b0)? && one(a1, b1)?) || (one(a0, b1)? && one(a1, b0)?) {
            turned += (piece.range.1 - piece.range.0).abs();
        }
    }
    Ok((turned - core::f64::consts::TAU).abs() <= ang)
}

/// A rim's section in the half-plane through its axis ([`rims_revolved`]).
pub(super) struct RimSection {
    pub(super) origin: Point,
    pub(super) axis: Vector,
    pub(super) radial: Vector,
    pub(super) rim: Point,
    pub(super) faces: [usize; 2],
    pub(super) contacts: [Point; 2],
}

/// The curve a surface of revolution about `section`'s axis traces in its
/// half-plane from `from` to `to`, and its range: a line on a plane square
/// to the axis or a coaxial cylinder or cone, an arc about the centre on a
/// coaxial sphere or the tube's centre there on a coaxial torus. `None`
/// for another surface or one not about the axis.
fn meridian_trace(
    surface: &SurfaceGeometry,
    section: &RimSection,
    (from, to): (Point, Point),
    tol: Tolerances,
) -> OgeomResult<Option<(Curve, (f64, f64))>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let (o, a) = (section.origin, section.axis);
    let off_axis = |p: Point| {
        let q = p - o;
        (q - a * q.dot(a)).magnitude()
    };
    let parallel = |d: Vector| d.cross(a).magnitude() <= ang;
    let about = match surface {
        SurfaceGeometry::Plane(p) if parallel(p.plane().normal().vector()) => None,
        SurfaceGeometry::Cylinder(c)
            if parallel(c.cylinder().axis().direction.vector())
                && off_axis(c.cylinder().axis().location) <= eps =>
        {
            None
        }
        SurfaceGeometry::Cone(c)
            if parallel(c.cone().axis().direction.vector())
                && off_axis(c.cone().axis().location) <= eps =>
        {
            None
        }
        SurfaceGeometry::Sphere(s) if off_axis(s.sphere().centre()) <= eps => {
            Some((s.sphere().centre(), s.sphere().radius()))
        }
        SurfaceGeometry::Torus(t)
            if parallel(t.torus().frame().z().vector())
                && off_axis(t.torus().frame().origin()) <= eps =>
        {
            let torus = t.torus();
            Some((
                torus.frame().origin() + section.radial * torus.major_radius(),
                torus.minor_radius(),
            ))
        }
        _ => return Ok(None),
    };
    let Some((centre, radius)) = about else {
        let Ok(line) = LineCurve::segment(from, to, tol) else {
            return Ok(None);
        };
        return Ok(Some((Curve::Line(line), (0.0, from.distance(to)))));
    };
    if (from.distance(centre) - radius).abs() > eps || (to.distance(centre) - radius).abs() > eps {
        return Ok(None);
    }
    // The shorter way round, in the half-plane.
    let normal = section.radial.cross(a);
    let x = Direction::new(from - centre, tol)?;
    let mut z = Direction::new(normal, tol)?;
    let mut frame = Frame::new(centre, z, x, tol)?;
    let reach = |frame: &Frame| {
        let w = to - centre;
        w.dot(frame.y().vector())
            .atan2(w.dot(frame.x().vector()))
            .rem_euclid(core::f64::consts::TAU)
    };
    let mut turn = reach(&frame);
    if turn > core::f64::consts::PI {
        z = Direction::new(-normal, tol)?;
        frame = Frame::new(centre, z, x, tol)?;
        turn = reach(&frame);
    }
    if turn <= ang {
        return Ok(None);
    }
    Ok(Some((
        Curve::Circle(CircleCurve::new(Circle::new(frame, radius, tol)?)),
        (0.0, turn),
    )))
}

/// The tool [`rims_revolved`] turns from `section` through `sweep`: the
/// solid, its face on the bevel, and whether it is cut or fused. `None`
/// where a face's trace is not a line or an arc about the axis, or the
/// turn will not build.
pub(super) fn revolved_tool(
    model: &mut Model,
    curved: &Curved,
    section: RimSection,
    sweep: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape, super::fallback::Sense)>> {
    let [c0, c1] = section.contacts;
    let mut curves: Vec<(Curve, (f64, f64))> = Vec::with_capacity(3);
    // Round the section: the rim to the first contact, and the second
    // contact back to the rim.
    for (k, ends) in [(0, (section.rim, c0)), (1, (c1, section.rim))] {
        let surface = &curved.faces[section.faces[k]].surface;
        let Some(trace) = meridian_trace(surface, &section, ends, tol)? else {
            return Ok(None);
        };
        curves.push(trace);
    }
    let Ok(bevel) = LineCurve::segment(c0, c1, tol) else {
        return Ok(None);
    };
    // The rim convex where the far contact stands inside the first face.
    let Some(out0) = outward(curved, section.faces[0], section.rim) else {
        return Ok(None);
    };
    let convex = (c1 - section.rim).dot(out0) < 0.0;
    let built = (|| -> OgeomResult<(Shape, Shape)> {
        let [r, a, b] = [section.rim, c0, c1].map(|p| ogeom_algo::make_vertex(model, p).shape);
        let to_c0 =
            ogeom_algo::make_edge_between(model, curves[0].0.clone(), curves[0].1, &r, &a, tol)?
                .shape;
        let across = ogeom_algo::make_edge_between(
            model,
            Curve::Line(bevel),
            (0.0, c0.distance(c1)),
            &a,
            &b,
            tol,
        )?
        .shape;
        let to_c1 =
            ogeom_algo::make_edge_between(model, curves[1].0.clone(), curves[1].1, &b, &r, tol)?
                .shape;
        let wire = ogeom_algo::make_wire(model, &[to_c0, across.clone(), to_c1], tol)?.shape;
        let normal = Direction::new(section.radial.cross(section.axis), tol)?;
        let plane = Frame::new(
            section.rim,
            normal,
            Direction::new(section.radial, tol)?,
            tol,
        )?;
        let face = ogeom_algo::make_face(
            model,
            PlaneSurface::new(Plane::new(plane)).into(),
            &[wire],
            tol,
        )?
        .shape;
        let turned = ogeom_algo::make_revolution(
            model,
            &face,
            ogeom_math::Axis {
                location: section.origin,
                direction: Direction::new(section.axis, tol)?,
            },
            sweep,
            tol,
        )?;
        // The face's own edge on the bevel's line, by its middle.
        let middle = c0 + (c1 - c0) * 0.5;
        let mut across = across;
        for e in ogeom_topo::explore_unique(model, &face, ogeom_topo::ShapeType::Edge)? {
            if let Some((curve, (t0, t1))) = ogeom_algo::placed_edge_curve(model, &e, tol)?
                && curve.point_at(f64::midpoint(t0, t1), tol)?.distance(middle) <= tol.confusion()
            {
                across = e;
            }
        }
        let swept: Vec<Shape> = turned
            .history
            .generated(&across)
            .iter()
            .filter(|x| model.kind_of(x).ok() == Some(ogeom_topo::ShapeType::Face))
            .cloned()
            .collect();
        let [ref named] = swept[..] else {
            ogeom_core::ogeom_bail!(Construction, "the bevel's turn is not one face");
        };
        let named = named.clone();
        let solid = if ogeom_algo::volume_properties(
            model,
            &turned.shape,
            ogeom_mesh::Deflection::default(),
            tol,
        )?
        .mass
            < 0.0
        {
            turned.shape.reversed()
        } else {
            turned.shape
        };
        Ok((solid, named))
    })();
    let Ok((solid, named)) = built else {
        return Ok(None);
    };
    Ok(Some((
        solid,
        named,
        if convex {
            super::fallback::Sense::Cut
        } else {
            super::fallback::Sense::Fuse
        },
    )))
}

/// How a strip's end is closed in its tool ([`close_tool`]).
#[derive(Debug, Clone, Copy)]
pub(crate) enum EndCap {
    /// The strip's own end, on the face it ends at: for now a plane
    /// square to the edge.
    Own,
    /// The plane through `at` with normal `normal`: a step square to the
    /// edge at a joint.
    Plane { at: Point, normal: Vector },
}

/// The closed tool of `strip` (ADR 0001, M5a, shared with rounds): the
/// region between its edge and the strip, bounded by the strip, its two
/// flaps (each face's piece between the edge and its contact, on the
/// face's own surface) and its end caps. Cut from the solid where the edge
/// is convex, fused into it where concave; the strip's face is named for
/// its edge. The strip is capped on a plane at each end, for now square to
/// the edge there (the strip's own end where it ends square, or a step at
/// a joint, as between the arcs of a rim), each cap's sides on the faces
/// (the arc about a wall's axis square to the cap, else straight). `None` for a strip along a closed edge (no ends to cap), any
/// other end, or where the pieces will not close.
///
/// # Errors
///
/// As building the faces.
pub(crate) fn close_tool(
    model: &mut Model,
    curved: &Curved,
    strip: &Strip,
    ends: [EndCap; 2],
    tol: Tolerances,
) -> OgeomResult<Option<super::fallback::Tool>> {
    let Some(edge) = curved.edges.get(&strip.edge) else {
        return Ok(None);
    };
    if edge.ends.0 == edge.ends.1 {
        return Ok(None);
    }
    open_tool(model, curved, strip, ends, tol)
}

/// Where `curve`, over `range`, crosses the plane through `at` with normal
/// `normal`, the crossing nearest `near`: by sampling, then bisection.
fn crossing(
    curve: &Curve,
    range: (f64, f64),
    at: Point,
    normal: Vector,
    near: Point,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let side = |t: f64| -> OgeomResult<f64> { Ok((curve.point_at(t, tol)? - at).dot(normal)) };
    let (lo, hi) = range;
    let n = 256;
    let mut best: Option<(f64, f64)> = None;
    let mut prev = (lo, side(lo)?);
    for i in 1..=n {
        let t = lo + (hi - lo) * f64::from(i) / f64::from(n);
        let here = (t, side(t)?);
        if prev.1 == 0.0 || prev.1 * here.1 < 0.0 || here.1 == 0.0 {
            let (mut a, mut b) = (prev, here);
            for _ in 0..200 {
                let m = 0.5 * (a.0 + b.0);
                let sm = side(m)?;
                if sm == 0.0 || (b.0 - a.0).abs() <= f64::EPSILON * (1.0 + m.abs()) {
                    a = (m, sm);
                    break;
                }
                if a.1 * sm < 0.0 {
                    b = (m, sm);
                } else {
                    a = (m, sm);
                }
            }
            let t = if a.1.abs() <= b.1.abs() { a.0 } else { b.0 };
            let d = curve.point_at(t, tol)?.distance(near);
            if best.is_none_or(|(_, e)| d < e) {
                best = Some((t, d));
            }
        }
        prev = here;
    }
    Ok(best.map(|(t, _)| t))
}

/// [`close_tool`] for an open strip, capped on a plane square to its edge
/// at each end.
fn open_tool(
    model: &mut Model,
    curved: &Curved,
    strip: &Strip,
    ends: [EndCap; 2],
    tol: Tolerances,
) -> OgeomResult<Option<super::fallback::Tool>> {
    open_tool_over(model, curved, strip, None, ends, tol)
}

/// [`close_tool`] for an open strip along `over` (a curve and its range)
/// in place of its edge where given: a strip running on past its edge's
/// ends, each cap a [`EndCap::Plane`] through an end of `over`.
///
/// # Errors
///
/// As building the faces.
pub(crate) fn open_tool_over(
    model: &mut Model,
    curved: &Curved,
    strip: &Strip,
    over: Option<(&Curve, (f64, f64))>,
    ends: [EndCap; 2],
    tol: Tolerances,
) -> OgeomResult<Option<super::fallback::Tool>> {
    let eps = tol.confusion();
    let edge = &curved.edges[&strip.edge];
    let (line, (e0, e1)) = over.unwrap_or((&edge.curve, edge.range));
    let h = (e1 - e0) * 1e-6;
    // Each end's point and its cap's plane.
    let mut caps = [(Point::ORIGIN, Vector::ZERO); 2];
    for (k, end) in ends.iter().enumerate() {
        let (t, inward) = if k == 0 { (e0, h) } else { (e1, -h) };
        let p = line.point_at(t, tol)?;
        let along = line.point_at(t + inward, tol)? - p;
        let along = along / along.magnitude();
        caps[k] = match end {
            // The strip's own end: on the one other face at its vertex,
            // a plane square to the edge there (any other end declines).
            EndCap::Own => {
                if over.is_some() {
                    return Ok(None);
                }
                let v = if k == 0 { edge.ends.0 } else { edge.ends.1 };
                let others: Vec<usize> = (0..curved.faces.len())
                    .filter(|fi| {
                        !strip.faces.contains(fi)
                            && curved.faces[*fi].loops.iter().flatten().any(|(x, _)| {
                                let c = &curved.edges[x];
                                c.ends.0 == v || c.ends.1 == v
                            })
                    })
                    .collect();
                let [other] = others[..] else {
                    return Ok(None);
                };
                let SurfaceGeometry::Plane(plane) = &curved.faces[other].surface else {
                    return Ok(None);
                };
                let n = plane.plane().normal().vector();
                if n.cross(along).magnitude() > tol.angular() {
                    return Ok(None);
                }
                (p, along)
            }
            EndCap::Plane { at, normal } => {
                if (p - *at).dot(*normal).abs() > eps * normal.magnitude() {
                    return Ok(None);
                }
                (*at, *normal / normal.magnitude())
            }
        };
    }
    // Each contact between its caps, and the section's centre at each.
    let mut contact_ranges = [(0.0, 0.0); 2];
    let mut corners = [[Point::ORIGIN; 2]; 2];
    for (k, contact) in strip.contacts.iter().enumerate() {
        let (lo, hi) = contact.domain();
        let period = hi - lo;
        let closed = contact
            .point_at(lo, tol)?
            .distance(contact.point_at(hi, tol)?)
            <= eps;
        let mut at = [0.0; 2];
        for (j, (cap_at, normal)) in caps.iter().enumerate() {
            let Some(t) = crossing(contact, (lo, hi), *cap_at, *normal, *cap_at, tol)? else {
                return Ok(None);
            };
            at[j] = t;
            corners[k][j] = contact.point_at(t, tol)?;
        }
        // The piece between the crossings beside the edge.
        let middle = line.point_at(0.5 * (e0 + e1), tol)?;
        let (a, b) = (at[0].min(at[1]), at[0].max(at[1]));
        let mut range = (a, b);
        if closed {
            let inner = contact.point_at(0.5 * (a + b), tol)?.distance(middle);
            let outer = contact
                .point_at(0.5 * (b + a + period), tol)?
                .distance(middle);
            if outer < inner {
                range = (b, a + period);
            }
        }
        contact_ranges[k] = range;
    }
    // A round's section is an arc about its spine's point at each cap; a
    // chamfer's, the segment between its contacts' ends there.
    let spine = match &strip.section {
        Section::Arc { spine } => Some(spine),
        Section::Line => None,
    };
    let mut centres = [Point::ORIGIN; 2];
    if let Some(spine) = spine {
        for (j, (cap_at, normal)) in caps.iter().enumerate() {
            let Some(t) = crossing(spine, spine.domain(), *cap_at, *normal, *cap_at, tol)? else {
                return Ok(None);
            };
            centres[j] = spine.point_at(t, tol)?;
        }
    }
    // The sheet's points: the edge's ends, then each contact's at each end.
    let points = vec![
        caps[0].0,
        caps[1].0,
        corners[0][0],
        corners[0][1],
        corners[1][0],
        corners[1][1],
    ];
    let corner_index = |k: usize, j: usize| 2 + 2 * k + j;
    let mut curves: Vec<SheetCurve> = Vec::new();
    let Some((edge_curve, edge_range)) = within_domain(line, (e0, e1), tol)? else {
        return Ok(None);
    };
    curves.push(SheetCurve {
        ends: (0, 1),
        curve: edge_curve,
        range: edge_range,
    });
    for k in 0..2 {
        let (lo, hi) = contact_ranges[k];
        let start = strip.contacts[k].point_at(lo, tol)?;
        let ends = if start.distance(corners[k][0]) <= start.distance(corners[k][1]) {
            (corner_index(k, 0), corner_index(k, 1))
        } else {
            (corner_index(k, 1), corner_index(k, 0))
        };
        let Some((curve, range)) = within_domain(&strip.contacts[k], (lo, hi), tol)? else {
            return Ok(None);
        };
        curves.push(SheetCurve { ends, curve, range });
    }
    // Each cap's sides, from the edge's end to each contact's, on the face:
    // the arc about a wall's axis square to the cap, else straight.
    for j in 0..2 {
        for (k, (corner, fi)) in corners.iter().zip(strip.faces).enumerate() {
            let (from, to) = (points[j], corner[j]);
            let face = &curved.faces[fi];
            let axis = match &face.surface {
                SurfaceGeometry::Cylinder(c) => {
                    let a = c.cylinder().axis();
                    Some((a.location, a.direction.vector()))
                }
                SurfaceGeometry::Cone(c) => {
                    let f = c.cone().frame();
                    Some((f.origin(), f.z().vector()))
                }
                _ => None,
            };
            let (curve, range) = match axis {
                Some((o, a)) if a.cross(caps[j].1).magnitude() <= tol.angular() => {
                    let centre = o + a * (from - o).dot(a);
                    let (x, y) = (from - centre, to - centre);
                    let z = x.cross(y);
                    if z.magnitude() <= eps * x.magnitude() {
                        return Ok(None);
                    }
                    let frame = Frame::new(
                        centre,
                        Direction::new(z, tol)?,
                        Direction::new(x, tol)?,
                        tol,
                    )?;
                    let arc: Curve =
                        CircleCurve::new(Circle::new(frame, x.magnitude(), tol)?).into();
                    (arc, (0.0, x.angle(y, tol)?))
                }
                _ => (
                    Curve::Line(LineCurve::segment(from, to, tol)?),
                    (0.0, from.distance(to)),
                ),
            };
            let middle = curve.point_at(0.5 * (range.0 + range.1), tol)?;
            // On a plane, off its plane: a cap past the end faces stands
            // beyond the face's own patch, on the same plane.
            let off = match &face.surface {
                SurfaceGeometry::Plane(plane) => {
                    let plane = plane.plane();
                    (middle - plane.origin()).dot(plane.normal().vector()).abs()
                }
                other => ogeom_algo::project_on_surface(other, middle, 16, tol)?.distance,
            };
            if off > eps * 10.0 {
                return Ok(None);
            }
            curves.push(SheetCurve {
                ends: (j, corner_index(k, j)),
                curve,
                range,
            });
        }
    }
    // Each end's section, the arc about its centre between the contacts,
    // or a chamfer's segment between them.
    for j in 0..2 {
        let (c0, c1) = (corners[0][j], corners[1][j]);
        if spine.is_none() {
            curves.push(SheetCurve {
                ends: (corner_index(0, j), corner_index(1, j)),
                curve: Curve::Line(LineCurve::segment(c0, c1, tol)?),
                range: (0.0, c0.distance(c1)),
            });
            continue;
        }
        let centre = centres[j];
        let (x, y) = (c0 - centre, c1 - centre);
        let radius = x.magnitude();
        let mut normal = x.cross(y);
        if normal.magnitude() <= eps * radius {
            return Ok(None);
        }
        normal = normal / normal.magnitude();
        let frame = Frame::new(
            centre,
            Direction::new(normal, tol)?,
            Direction::new(x, tol)?,
            tol,
        )?;
        let sweep = x.angle(y, tol)?;
        curves.push(SheetCurve {
            ends: (corner_index(0, j), corner_index(1, j)),
            curve: CircleCurve::new(Circle::new(frame, radius, tol)?).into(),
            range: (0.0, sweep),
        });
    }
    // The faces' loops by curve: flaps, the strip, the caps.
    let side = |j: usize, k: usize| 3 + 2 * j + k;
    let section = |j: usize| 7 + j;
    let faces: Vec<(SurfaceGeometry, Vec<usize>)> = vec![
        (
            curved.faces[strip.faces[0]].surface.clone(),
            vec![0, side(1, 0), 1, side(0, 0)],
        ),
        (
            curved.faces[strip.faces[1]].surface.clone(),
            vec![0, side(1, 1), 2, side(0, 1)],
        ),
        (strip.surface.clone(), vec![1, section(1), 2, section(0)]),
        (
            PlaneSurface::new(Plane::through(caps[0].0, Direction::new(caps[0].1, tol)?)).into(),
            vec![side(0, 0), section(0), side(0, 1)],
        ),
        (
            PlaneSurface::new(Plane::through(caps[1].0, Direction::new(caps[1].1, tol)?)).into(),
            vec![side(1, 0), section(1), side(1, 1)],
        ),
    ];
    // Each loop walked end to end, then turned so the face lies on its left
    // about its surface's own normal.
    let mut rings: Vec<Vec<(usize, bool)>> = Vec::with_capacity(faces.len());
    for (surface, loop_) in &faces {
        let mut walk: Vec<(usize, bool)> = Vec::with_capacity(loop_.len());
        let first = &curves[loop_[0]];
        walk.push((loop_[0], true));
        let mut at = first.ends.1;
        let mut left: Vec<usize> = loop_[1..].to_vec();
        while !left.is_empty() {
            let Some(k) = left
                .iter()
                .position(|c| curves[*c].ends.0 == at || curves[*c].ends.1 == at)
            else {
                return Ok(None);
            };
            let c = left.remove(k);
            let forward = curves[c].ends.0 == at;
            at = if forward {
                curves[c].ends.1
            } else {
                curves[c].ends.0
            };
            walk.push((c, forward));
        }
        if at != first.ends.0 {
            return Ok(None);
        }
        // Across the loop from its first curve's middle, toward its
        // opposite side: the third curve of four, or the corner facing it
        // of three.
        let mid = |c: usize| -> OgeomResult<Point> {
            let (lo, hi) = curves[c].range;
            curves[c].curve.point_at(0.5 * (lo + hi), tol)
        };
        let p = mid(walk[0].0)?;
        let opposite = match walk.len() {
            4 => mid(walk[2].0)?,
            3 => {
                let (c, forward) = walk[1];
                points[if forward {
                    curves[c].ends.1
                } else {
                    curves[c].ends.0
                }]
            }
            _ => return Ok(None),
        };
        let toward = opposite - p;
        let (lo, hi) = curves[walk[0].0].range;
        let dt = (hi - lo) * 1e-4;
        let t = 0.5 * (lo + hi);
        let along = curves[walk[0].0].curve.point_at(t + dt, tol)?
            - curves[walk[0].0].curve.point_at(t - dt, tol)?;
        let n = chart_normal(surface, p, tol)?;
        let toward = toward - n * (toward.dot(n) / n.dot(n).max(f64::MIN_POSITIVE));
        let turn = n.cross(along).dot(toward);
        if n.magnitude() <= eps || turn.abs() <= f64::EPSILON {
            return Ok(None);
        }
        if turn < 0.0 {
            walk = walk.iter().rev().map(|(c, f)| (*c, !f)).collect();
        }
        rings.push(walk);
    }
    // Convex where the section's chord stands inside the first face's
    // material.
    let out0 = outward_at(curved, strip.faces[0], corners[0][0]);
    if out0.magnitude() <= eps {
        return Ok(None);
    }
    let middle = corners[0][0] + (corners[1][0] - corners[0][0]) * 0.5;
    let convex = (middle - corners[0][0]).dot(out0) < 0.0;
    // Each face turned so every curve is walked once each way: the first
    // flap as its face stands where the tool lies in the solid.
    let mut reversed: Vec<Option<bool>> = vec![None; faces.len()];
    reversed[0] = Some(curved.faces[strip.faces[0]].reversed == convex);
    let flag = |f: usize, c: usize| rings[f].iter().find(|(x, _)| *x == c).map(|(_, d)| *d);
    for _ in 0..faces.len() {
        for f in 0..faces.len() {
            let Some(rf) = reversed[f] else {
                continue;
            };
            for g in 0..faces.len() {
                if g == f || reversed[g].is_some() {
                    continue;
                }
                for (c, _) in &rings[g] {
                    if let (Some(df), Some(dg)) = (flag(f, *c), flag(g, *c)) {
                        reversed[g] = Some(dg == (df != rf));
                        break;
                    }
                }
            }
        }
    }
    let Some(reversed): Option<Vec<bool>> = reversed.into_iter().collect() else {
        return Ok(None);
    };
    for c in 0..curves.len() {
        let uses: Vec<bool> = (0..faces.len())
            .filter_map(|f| flag(f, c).map(|d| d != reversed[f]))
            .collect();
        if uses.len() != 2 || uses[0] == uses[1] {
            return Ok(None);
        }
    }
    let sheet_rings: Vec<SheetRings> = faces
        .into_iter()
        .zip(rings)
        .zip(&reversed)
        .map(|(((surface, _), ring), r)| SheetRings {
            surface,
            rings: vec![ring],
            reversed: *r,
        })
        .collect();
    let Ok(sheet) = super::curved::build(model, &points, &curves, &sheet_rings, tol) else {
        return Ok(None);
    };
    let shell = ogeom_algo::make_shell(model, &sheet.faces)?.shape;
    if !ogeom_algo::is_shell_closed(model, &shell)? {
        return Ok(None);
    }
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    let Some(from) = curved.shapes.get(&strip.edge) else {
        return Ok(None);
    };
    Ok(Some(super::fallback::Tool {
        solid,
        sense: if convex {
            super::fallback::Sense::Cut
        } else {
            super::fallback::Sense::Fuse
        },
        faces: vec![(sheet.faces[2].clone(), from.clone())],
    }))
}

/// `curve` over `range` as a curve whose own domain holds the range: a
/// circle's piece running past its domain (a turn re-anchored) on the same
/// circle turned to start there. `None` for another curve past its domain.
fn within_domain(
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Option<(Curve, (f64, f64))>> {
    let (lo, hi) = curve.domain();
    if range.0 >= lo - f64::EPSILON && range.1 <= hi + f64::EPSILON {
        return Ok(Some((curve.clone(), range)));
    }
    let Curve::Circle(c) = curve else {
        return Ok(None);
    };
    let circle = c.circle();
    let frame = circle.frame();
    let start = curve.point_at(range.0, tol)? - frame.origin();
    let turned = Frame::new(frame.origin(), frame.z(), Direction::new(start, tol)?, tol)?;
    Ok(Some((
        CircleCurve::new(Circle::new(turned, circle.radius(), tol)?).into(),
        (0.0, range.1 - range.0),
    )))
}

/// `surface`'s own normal at `p`, its chart's `du x dv` at the foot of `p`:
/// the sense the sheet walks a face's rings about.
fn chart_normal(surface: &SurfaceGeometry, p: Point, tol: Tolerances) -> OgeomResult<Vector> {
    let foot = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?;
    let (u, v) = foot.parameters;
    let (du, dv) = ogeom_geom::Surface::d1_at(surface, u, v, tol)?;
    Ok(du.cross(dv))
}
