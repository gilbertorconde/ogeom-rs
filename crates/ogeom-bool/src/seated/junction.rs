//! A straight band's wedge that runs, at an end, into a blend the solid
//! already has, read for the straight trim ([`super::trimmed_cut`]) as the
//! part of it inside the solid.
//!
//! Where two edges of a block meeting at a corner are rounded at one radius
//! one after the other, the second band's wedge runs past the first band's
//! blend into the region the first cut took away. The general fuse cuts it
//! there along the mitre: the two blend cylinders, of one radius with axes
//! that meet, cross in two ellipses in the bisector planes of their axes,
//! and the one between the bands is where the second blend leaves the
//! first.
//!
//! In the corner's frame the three faces at the corner are square to each
//! other: the plane both bands round (`T`), the plane the wedge's other leg
//! lies on (`S`), and the plane its cap at that end lies on, which the first
//! band rounded from `T`. With the corner at the origin, the material at
//! `x, y, z <= 0` behind `S: x = 0`, the cap `y = 0` and `T: z = 0`, and
//! the radius `r`, the wedge's axis runs along `y` at `x = z = -r` and the
//! first band's along `x` at `y = r, z = -r`; the axes meet at
//! `O = (-r, r, -r)`. The conditions read, where one fails that end is not
//! read as a junction:
//!
//! - The wedge is a straight band's: one cylinder strip between two
//!   rulings and two arcs ([`super::blend_cylinder`]), two legs along the
//!   axis square to each other, each at the radius from it, and two caps
//!   square to it; its vertices and arcs lie in the region its four planes
//!   and the chord through the rulings bound ([`super::strip_walls`]).
//! - The solid has a vertex where the wedge's `S` contact meets the cap,
//!   and a face there on a cylinder of the wedge's radius, within the
//!   confusion distance, whose axis runs through `O` along `S`'s normal.
//!   It faces away from its axis for a cut, toward it for a fill. It is a
//!   strip as the wedge's blend is, its rulings on `T` and on the cap's
//!   plane, its arcs one on `S` and one further from `S` than every vertex
//!   of the wedge. That face is the first band's blend, the band from
//!   `(-r, r, 0)` and `(-r, 0, -r)` the quarter turn round toward the
//!   corner, over the whole depth of the wedge.
//! - The edge along the corner keeps more than the boolean's weld reach of
//!   its length once each junction end takes the radius off it.
//!
//! Let `Q` be the region the wedge lies in: the triangular prism behind its
//! legs, the chord and its caps. Every face of the solid but the ones its
//! seats lie on (checked where the trim is planned) stands beyond one of
//! `Q`'s planes, the faces on `T` and `S` and at an end not a junction on
//! the cap lie on them, and at each junction end the first band's face
//! meets the inside of `Q` exactly in the quarter of its cylinder there:
//! its strip spans `Q`'s depth, and the cylinder meets `Q` only on the half
//! turn above its axis's level, of which the quarter on the corner's side is
//! the strip's. Each slice of `Q` square to the wedge's axis is the triangle
//! between `T`, `S` and the chord; at depth `0 < y < r` the first band's
//! face crosses it along `z = -r + sqrt(2ry - y^2)`, and the part above,
//! between the face and the corner, is outside the solid (the face's
//! outward side). The rest of each slice is convex and the slices run on
//! unbroken, so the rest of `Q` is one region, inside the solid (it lies
//! under the seat on `T`). On the wedge's cylinder, `(x + r)^2 + (z + r)^2
//! = r^2`, the first band's is crossed where `(r - y)^2 = (x + r)^2`, that
//! is on `x + y = 0`: the mitre ellipse, centred on `O`, its axes `r` along
//! `T`'s normal and `r sqrt 2` along `S`'s normal plus the cap's.
//!
//! So the wedge's part inside the solid is the wedge less, at each junction
//! end, the region past the first band's face: bounded by the seat on `T`
//! between the bands' contacts there, the seat on `S` down to the first
//! band's arc on it, a patch of the first band's face between that arc, its
//! contact on `T` and the ellipse, and the wedge's blend cut at the ellipse.
//! Its seats lie on faces of the solid as the trim asks, and the solid less
//! it is the solid less the wedge; it is built here as a solid of its own
//! for the trim to read, and the history is recorded against the wedge's
//! faces, as the general fuse records it.

use super::ring;
use super::{Solid, Walked, blend_cylinder, curve_of, plane_of, strip_walls, vertex_at, walk_face};
use ogeom_algo::{make_edge_between, make_wire};
use ogeom_core::{FastMap, OgeomResult, Tolerance, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, NodeData, Orientation, Shape, TShapeId};

/// A wedge cut back where it runs into blends of the solid.
pub(super) struct Junction {
    /// The part of the wedge inside the solid, as a solid.
    pub(super) wedge: Shape,
    /// Its blend, the wedge's cut at the mitres.
    pub(super) blend: Shape,
    /// The faces of the solid it runs into, by their place in its shell.
    pub(super) bands: Vec<usize>,
    /// The planes of the region the whole wedge lies in, each a point of it
    /// and its outward unit normal.
    pub(super) walls: Vec<(Point, Vector)>,
    /// The wedge's blend, which the result's blend is recorded against.
    pub(super) given_blend: Shape,
    /// The wedge's planar faces, which the result deletes.
    pub(super) given: Vec<Shape>,
}

/// One end of the wedge read as a junction.
struct End {
    /// The axes' meeting point.
    centre: Point,
    /// The outward unit normals of `T`, `S` and the cap.
    t: Vector,
    s: Vector,
    cap: Vector,
    /// The wedge's vertices at the cap: on both legs, on `T` and the blend,
    /// on `S` and the blend.
    corner: TShapeId,
    on_t: TShapeId,
    on_s: TShapeId,
    /// The cap's edges: on `T`, on `S`, and its arc.
    line: TShapeId,
    side: TShapeId,
    mitre: TShapeId,
    /// The band's place in the solid's shell.
    band: usize,
}

/// A planar face of the wedge: the face, its loop, a point of its plane and
/// its outward unit normal.
type Planar = (Shape, Vec<Walked>, Point, Vector);

/// The wedge cut back at each end where it runs into a blend of the solid
/// as the module reads one; `None` where no end does, where an end reads
/// two ways, or where the wedge is not a straight band's as it asks.
pub(super) fn cut_back(
    model: &mut Model,
    solid: &Solid,
    wedge: &Shape,
    fill: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Junction>> {
    let shells = model.children_of(wedge)?;
    let [shell] = shells.as_slice() else {
        return Ok(None);
    };
    let faces = model.children_of(shell)?;
    if faces.len() != 5 {
        return Ok(None);
    }
    let mut blend: Option<(Shape, Vec<Walked>)> = None;
    let mut planar: Vec<Planar> = Vec::with_capacity(4);
    for face in &faces {
        let Some(mut rings) = walk_face(model, face)? else {
            return Ok(None);
        };
        if rings.len() != 1 {
            return Ok(None);
        }
        let ring = rings.remove(0);
        if let Some((_, origin, outward)) = plane_of(model, face) {
            planar.push((face.clone(), ring, origin, outward));
        } else if blend.is_none() {
            blend = Some((face.clone(), ring));
        } else {
            return Ok(None);
        }
    }
    let Some((blend, blend_loop)) = blend else {
        return Ok(None);
    };
    let Some(strip) = blend_cylinder(model, &blend, &blend_loop, tol)? else {
        return Ok(None);
    };
    let Some((_, SurfaceGeometry::Cylinder(c))) = ring::surface_of(model, &blend) else {
        return Ok(None);
    };
    let cylinder = c.cylinder();
    let (base, axis, radius) = (
        cylinder.frame().origin(),
        cylinder.frame().z().vector(),
        cylinder.radius(),
    );
    // The wedge's vertices and the widest tolerance it states.
    let mut points: FastMap<TShapeId, (Shape, Point, f64)> = FastMap::default();
    let mut stated = tol.confusion();
    for walked in planar.iter().flat_map(|p| &p.1).chain(&blend_loop) {
        stated = stated.max(model.tolerance_of(&walked.edge)?.map_or(0.0, |t| t.get()));
        let Some((start, end)) = ogeom_algo::edge_vertices(model, &walked.edge)? else {
            return Ok(None);
        };
        for v in [start, end] {
            let Some((p, t)) = vertex_at(model, &v) else {
                return Ok(None);
            };
            stated = stated.max(t);
            points
                .entry(v.node())
                .or_insert((v.oriented(Orientation::Forward), p, t));
        }
    }
    // Legs share a ruling with the blend, caps an arc.
    let (mut legs, mut caps): (Vec<&Planar>, Vec<&Planar>) = (Vec::new(), Vec::new());
    for face in &planar {
        let mut along = face.1.iter().filter_map(|w| {
            blend_loop
                .iter()
                .find(|b| b.edge.node() == w.edge.node())
                .and_then(|b| curve_of(model, &b.edge))
        });
        match (along.next(), along.next()) {
            (Some((Curve::Line(_), _, _)), None) if face.3.dot(axis).abs() <= tol.angular() => {
                legs.push(face);
            }
            (Some((Curve::Circle(_), _, _)), None) if ring::parallel(face.3, axis, tol) => {
                caps.push(face);
            }
            _ => return Ok(None),
        }
    }
    if legs.len() != 2 || caps.len() != 2 {
        return Ok(None);
    }
    let corner_points: Vec<Point> = points.values().map(|(_, p, _)| *p).collect();
    let mut ends: Vec<(usize, End)> = Vec::new();
    for (k, cap) in caps.iter().enumerate() {
        let (_, _, cap_origin, _) = cap;
        let foot = base + axis * (*cap_origin - base).dot(axis);
        let read = junction_end(
            model,
            solid,
            cap,
            &legs,
            (foot, radius),
            &blend_loop,
            &points,
            fill,
            tol,
        )?;
        if let Some(end) = read {
            ends.push((k, end));
        }
    }
    if ends.is_empty() {
        return Ok(None);
    }
    // The wedge's region, which each end's reading assumes.
    let walls_from: Vec<(Point, Vector)> = planar.iter().map(|p| (p.2, p.3)).collect();
    let Some(walls) = strip_walls(&strip, walls_from, &corner_points, stated, tol)? else {
        return Ok(None);
    };
    // Each end's cut takes the radius off the corner's length.
    let along = |cap: &Planar| (cap.2 - base).dot(axis);
    let length = (along(caps[0]) - along(caps[1])).abs();
    #[allow(clippy::cast_precision_loss, reason = "one or two ends")]
    let left = length - radius * ends.len() as f64;
    if left <= crate::sub_shapes::reach(tol) + stated {
        return Ok(None);
    }
    let made = build(
        model,
        (&blend, &blend_loop),
        &planar,
        &caps,
        &ends,
        solid,
        radius,
        stated,
        tol,
    )?;
    let Some((cut, cut_blend)) = made else {
        return Ok(None);
    };
    Ok(Some(Junction {
        wedge: cut,
        blend: cut_blend,
        bands: ends.iter().map(|(_, e)| e.band).collect(),
        walls,
        given_blend: blend,
        given: planar.iter().map(|p| p.0.clone()).collect(),
    }))
}

/// The cap's end read as a junction ([`End`]), `None` where it is not one;
/// `foot` is where the wedge's axis meets the cap's plane.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
fn junction_end(
    model: &Model,
    solid: &Solid,
    cap: &Planar,
    legs: &[&Planar],
    (foot, radius): (Point, f64),
    blend_loop: &[Walked],
    points: &FastMap<TShapeId, (Shape, Point, f64)>,
    fill: bool,
    tol: Tolerances,
) -> OgeomResult<Option<End>> {
    let (_, ring, cap_origin, e) = cap;
    if ring.len() != 3 {
        return Ok(None);
    }
    let Some(arc) = ring
        .iter()
        .find(|w| blend_loop.iter().any(|b| b.edge.node() == w.edge.node()))
    else {
        return Ok(None);
    };
    let Some(corner) = ring
        .iter()
        .map(|w| w.from)
        .find(|v| *v != arc.from && *v != arc.to)
    else {
        return Ok(None);
    };
    let mut found: Option<End> = None;
    for (ti, si) in [(0, 1), (1, 0)] {
        let (t_leg, s_leg) = (legs[ti], legs[si]);
        let (n_t, n_s) = (t_leg.3, s_leg.3);
        // The wedge's own corner frame: legs square to each other, the cap's
        // vertices where the frame and the radius put them.
        if n_t.dot(n_s).abs() > tol.angular() {
            continue;
        }
        let Some((on_t, on_s, reach)) =
            cap_corner(points, corner, arc, foot, (n_t, n_s), radius, tol)
        else {
            continue;
        };
        let centre = foot - *e * radius;
        let joining = |x: TShapeId, y: TShapeId| {
            ring.iter()
                .find(|w| (w.from, w.to) == (x, y) || (w.from, w.to) == (y, x))
                .map(|w| w.edge.node())
        };
        let (Some(line), Some(side)) = (joining(corner, on_t), joining(corner, on_s)) else {
            continue;
        };
        // The solid's vertex where the S contact meets the cap.
        let (_, s_point, s_tol) = &points[&on_s];
        let mut near = solid
            .vertices
            .iter()
            .filter(|(_, q, t)| s_point.distance(*q) <= s_tol + t);
        let (Some((held, _, _)), None) = (near.next(), near.next()) else {
            continue;
        };
        let mut seen: Vec<usize> = Vec::new();
        for edge in solid
            .by_vertex
            .get(&held.node())
            .map_or(&[][..], |v| v.as_slice())
        {
            for &f in solid
                .walkers
                .get(&edge.node())
                .map_or(&[][..], |v| v.as_slice())
            {
                if seen.contains(&f) {
                    continue;
                }
                seen.push(f);
                let end = End {
                    centre,
                    t: n_t,
                    s: n_s,
                    cap: *e,
                    corner,
                    on_t,
                    on_s,
                    line,
                    side,
                    mitre: arc.edge.node(),
                    band: f,
                };
                let planes = (t_leg.2, s_leg.2, *cap_origin);
                if is_band(model, solid, &end, radius, planes, points, reach, fill, tol)? {
                    if found.is_some() {
                        return Ok(None);
                    }
                    found = Some(end);
                }
            }
        }
    }
    Ok(found)
}

/// The wedge's vertices on `T` and on `S` at a cap, and the reach the
/// reading allows, where the cap's corner and the arc's ends stand within
/// their tolerances of where the corner frame puts them: off the axis's
/// foot on the cap by the radius along both legs' normals, along `T`'s, and
/// along `S`'s.
fn cap_corner(
    points: &FastMap<TShapeId, (Shape, Point, f64)>,
    corner: TShapeId,
    arc: &Walked,
    foot: Point,
    (n_t, n_s): (Vector, Vector),
    radius: f64,
    tol: Tolerances,
) -> Option<(TShapeId, TShapeId, f64)> {
    let (_, corner_at, corner_tol) = points.get(&corner)?;
    let (_, a, a_tol) = points.get(&arc.from)?;
    let (_, b, b_tol) = points.get(&arc.to)?;
    let reach = corner_tol.max(*a_tol).max(*b_tol) + tol.confusion();
    let (t_at, s_at) = (foot + n_t * radius, foot + n_s * radius);
    let (on_t, on_s) = if a.distance(t_at) <= reach && b.distance(s_at) <= reach {
        (arc.from, arc.to)
    } else if b.distance(t_at) <= reach && a.distance(s_at) <= reach {
        (arc.to, arc.from)
    } else {
        return None;
    };
    (corner_at.distance(foot + (n_t + n_s) * radius) <= reach).then_some((on_t, on_s, reach))
}

/// Whether the solid's face at `end.band` is the band the wedge runs into
/// at `end`, as the module reads one.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
fn is_band(
    model: &Model,
    solid: &Solid,
    end: &End,
    radius: f64,
    (t_origin, s_origin, cap_origin): (Point, Point, Point),
    points: &FastMap<TShapeId, (Shape, Point, f64)>,
    reach: f64,
    fill: bool,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let face = &solid.faces[end.band];
    let Some((_, SurfaceGeometry::Cylinder(c))) = ring::surface_of(model, face) else {
        return Ok(false);
    };
    let cylinder = c.cylinder();
    let frame = cylinder.frame();
    let b = frame.z().vector();
    let off = end.centre - frame.origin();
    // One radius, the axes meeting, as the boolean reads the mitre.
    if (cylinder.radius() - radius).abs() > tol.confusion()
        || !ring::parallel(b, end.s, tol)
        || (off - b * off.dot(b)).magnitude() > tol.confusion()
        || ring::facing(face, &frame) != if fill { -1.0 } else { 1.0 }
    {
        return Ok(false);
    }
    let [ring] = solid.loops[end.band].as_slice() else {
        return Ok(false);
    };
    let Some(strip) = blend_cylinder(model, face, ring, tol)? else {
        return Ok(false);
    };
    let on = |p: Point, o: Point, n: Vector| (p - o).dot(n).abs() <= reach;
    let [p, q] = strip.rulings;
    let rulings = (on(p, t_origin, end.t) && on(q, cap_origin, end.cap))
        || (on(q, t_origin, end.t) && on(p, cap_origin, end.cap));
    if !rulings {
        return Ok(false);
    }
    let centres: Vec<Point> = strip
        .arcs
        .iter()
        .filter_map(|(curve, _)| match curve {
            Curve::Circle(c) => Some(c.circle().centre()),
            _ => None,
        })
        .collect();
    let [m, n] = centres.as_slice() else {
        return Ok(false);
    };
    let (near, far) = if on(*m, s_origin, end.s) {
        (*m, *n)
    } else {
        (*n, *m)
    };
    if !on(near, s_origin, end.s) {
        return Ok(false);
    }
    // Every vertex of the wedge, and so its region, stands inside the
    // strip's depth by more than the boolean welds.
    let margin = crate::sub_shapes::reach(tol) + reach;
    Ok(points
        .values()
        .all(|(_, v, _)| (*v - far).dot(end.s) > margin))
}

/// A replacement edge of the cut-back wedge, running from one of its
/// vertices to another.
struct Run {
    edge: Shape,
    from: TShapeId,
}

impl Run {
    /// The edge turned to run from `from`.
    fn from(&self, from: TShapeId) -> Shape {
        if from == self.from {
            self.edge.clone()
        } else {
            self.edge.reversed()
        }
    }
}

/// An edge on `curve` over `range`, from `a` to `b` where the range runs
/// up from `a`'s end; turned to run from `a` either way.
fn run_between(
    model: &mut Model,
    curve: Curve,
    (ta, tb): (f64, f64),
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Run> {
    let edge = if ta <= tb {
        make_edge_between(model, curve, (ta, tb), a, b, tol)?.shape
    } else {
        make_edge_between(model, curve, (tb, ta), b, a, tol)?
            .shape
            .reversed()
    };
    Ok(Run {
        edge,
        from: a.node(),
    })
}

/// The parameter of a point on an ellipse curve at no reversal.
fn ellipse_at(curve: &Curve, p: Point) -> Option<f64> {
    let Curve::Ellipse(e) = curve else {
        return None;
    };
    if e.is_reversed() {
        return None;
    }
    let ellipse = e.ellipse();
    let local = ellipse.frame().to_local(p);
    Some((local.y / ellipse.minor_radius()).atan2(local.x / ellipse.major_radius()))
}

/// The cut-back wedge made as a solid, and its blend, `None` where a face
/// of it does not take its ring.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn build(
    model: &mut Model,
    (blend, blend_loop): (&Shape, &[Walked]),
    planar: &[Planar],
    caps: &[&Planar],
    ends: &[(usize, End)],
    solid: &Solid,
    radius: f64,
    stated: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape)>> {
    let tolerance = Tolerance::new(stated.max(tol.confusion()))?;
    // At a junction end the wedge's corner, `T` contact and `S` contact
    // move to the seats' corner, the crossing of the two bands' contacts on
    // `T`, and the `S` contact, each made anew.
    let mut moved: FastMap<TShapeId, Shape> = FastMap::default();
    for (_, end) in ends {
        let at = [
            (end.corner, end.centre + (end.t + end.s) * radius),
            (end.on_t, end.centre + end.t * radius),
            (end.on_s, end.centre + (end.s + end.cap) * radius),
        ];
        for (id, point) in at {
            let v = model.add_vertex(ogeom_topo::VertexData { point, tolerance });
            moved.insert(id, v);
        }
    }
    let renamed = |id: TShapeId| moved.get(&id).map_or(id, Shape::node);
    // The faces carried over: the blend, the legs, and each cap at an end
    // that is not a junction.
    let mut sources: Vec<(Shape, Vec<Walked>)> = Vec::with_capacity(planar.len() + 1);
    sources.push((blend.clone(), blend_loop.to_vec()));
    for face in planar {
        let at = caps.iter().position(|c| c.0.is_same(&face.0));
        if at.is_some_and(|k| ends.iter().any(|(j, _)| *j == k)) {
            continue;
        }
        sources.push((face.0.clone(), face.1.clone()));
    }
    let mut runs: FastMap<TShapeId, Run> = FastMap::default();
    let mut faces: Vec<Shape> = Vec::with_capacity(sources.len() + ends.len());
    let mut cut_blend: Option<Shape> = None;
    for (face, ring) in &sources {
        let mut walked: Vec<Shape> = Vec::with_capacity(ring.len());
        let mut changed = false;
        for w in ring {
            let id = w.edge.node();
            if !runs.contains_key(&id) {
                let Some(run) = remade(model, &w.edge, &moved, ends, solid, blend, radius, tol)?
                else {
                    return Ok(None);
                };
                runs.insert(id, run);
            }
            let run = &runs[&id];
            changed |= run.edge.node() != id;
            walked.push(run.from(renamed(w.from)));
        }
        let made = if changed {
            face_from(model, face, face.orientation(), &walked, tol)?
        } else {
            face.clone()
        };
        if face.is_same(blend) {
            cut_blend = Some(made.clone());
        }
        faces.push(made);
    }
    // At each junction end, the patch of the band between its contact on
    // `T`, its arc on `S` and the mitre, facing away from the band's axis.
    for (_, end) in ends {
        let (Some(line), Some(arc), Some(mitre)) = (
            runs.get(&end.line),
            runs.get(&end.side),
            runs.get(&end.mitre),
        ) else {
            return Ok(None);
        };
        let [on_t, corner, on_s] = [end.on_t, end.corner, end.on_s].map(|id| moved[&id].clone());
        let at = |v: &Shape| vertex_at(model, v).map(|(p, _)| p);
        let (Some(p1), Some(p2), Some(p3)) = (at(&on_t), at(&corner), at(&on_s)) else {
            return Ok(None);
        };
        let ring = if (p3 - p1).cross(p2 - p3).dot(end.t + end.cap) > 0.0 {
            vec![
                mitre.from(on_t.node()),
                arc.from(on_s.node()),
                line.from(corner.node()),
            ]
        } else {
            vec![
                line.from(on_t.node()),
                arc.from(corner.node()),
                mitre.from(on_s.node()),
            ]
        };
        let band = &solid.faces[end.band];
        let Some((_, SurfaceGeometry::Cylinder(c))) = ring::surface_of(model, band) else {
            return Ok(None);
        };
        let away = if c.cylinder().frame().handedness() == ogeom_math::Handedness::Left {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
        faces.push(face_from(model, band, away, &ring, tol)?);
    }
    let Some(cut_blend) = cut_blend else {
        return Ok(None);
    };
    let shell = model.add_shell(&faces)?;
    let cut = model.add_solid(&[shell])?;
    Ok(Some((cut, cut_blend)))
}

/// A face on the surface of `like`, used `orientation`, bounded by `ring`
/// walked as the face walks it.
fn face_from(
    model: &mut Model,
    like: &Shape,
    orientation: Orientation,
    ring: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let Some(NodeData::Face(data)) = model.node(like).map(|n| n.data()) else {
        ogeom_core::ogeom_bail!(Construction, "a face of the wedge holds no face data");
    };
    let mut data = (**data).clone();
    data.triangulation = None;
    data.natural_restriction = false;
    let edges: Vec<Shape> = if orientation == Orientation::Reversed {
        ring.iter().rev().map(Shape::reversed).collect()
    } else {
        ring.to_vec()
    };
    let wire = make_wire(model, &edges, tol)?.shape;
    Ok(model.add_face(data, &[wire])?.oriented(orientation))
}

/// A wedge edge as the cut-back wedge holds it, between its ends as they
/// moved: itself where neither did; at a junction end, the contacts'
/// crossing to the corner along the band's contact, the band's arc on `S`
/// from the corner to the `S` contact, and the mitre between the contacts;
/// and an edge along the axis on its own line, over the range between its
/// ends.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn remade(
    model: &mut Model,
    edge: &Shape,
    moved: &FastMap<TShapeId, Shape>,
    ends: &[(usize, End)],
    solid: &Solid,
    blend: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Run>> {
    let forward = edge.oriented(Orientation::Forward);
    let Some((start, finish)) = ogeom_algo::edge_vertices(model, &forward)? else {
        return Ok(None);
    };
    let ids = (start.node(), finish.node());
    let place = |v: &Shape| {
        moved
            .get(&v.node())
            .cloned()
            .unwrap_or_else(|| v.oriented(Orientation::Forward))
    };
    let (a, b) = (place(&start), place(&finish));
    if a.node() == ids.0 && b.node() == ids.1 {
        return Ok(Some(Run {
            edge: forward,
            from: ids.0,
        }));
    }
    let (Some((pa, _)), Some((pb, _))) = (vertex_at(model, &a), vertex_at(model, &b)) else {
        return Ok(None);
    };
    let joins = |x: TShapeId, y: TShapeId| ids == (x, y) || ids == (y, x);
    for (_, end) in ends {
        if joins(end.corner, end.on_t) {
            let line = ogeom_geom::LineCurve::segment(pa, pb, tol)?;
            let curve = Curve::Line(line);
            return run_between(model, curve, (0.0, pa.distance(pb)), &a, &b, tol).map(Some);
        }
        if joins(end.corner, end.on_s) {
            let centre = end.centre + end.s * radius;
            let (Ok(x), Ok(y)) = (Direction::new(end.t, tol), Direction::new(end.cap, tol)) else {
                return Ok(None);
            };
            let Ok(z) = Direction::new(end.t.cross(end.cap), tol) else {
                return Ok(None);
            };
            let frame = Frame::from_axes(centre, x, y, z, tol)?;
            let circle = ogeom_math::Circle::new(frame, radius, tol)?;
            let curve = Curve::Circle(ogeom_geom::CircleCurve::new(circle));
            let quarter = core::f64::consts::FRAC_PI_2;
            let range = if ids.0 == end.corner {
                (0.0, quarter)
            } else {
                (quarter, 0.0)
            };
            return run_between(model, curve, range, &a, &b, tol).map(Some);
        }
        if joins(end.on_t, end.on_s) {
            return mitre(
                model,
                &solid.faces[end.band],
                blend,
                (&a, pa),
                (&b, pb),
                tol,
            );
        }
    }
    let Some((curve @ Curve::Line(_), _, _)) = curve_of(model, &forward) else {
        return Ok(None);
    };
    let Curve::Line(line) = &curve else {
        return Ok(None);
    };
    let axis = line.axis();
    let at = |p: Point| (p - axis.location).dot(axis.direction.vector());
    let range = (at(pa), at(pb));
    run_between(model, curve, range, &a, &b, tol).map(Some)
}

/// The mitre between the band and the wedge's blend from `a` to `b`: the
/// ellipse of the two cylinders' crossing the boolean reads, through both,
/// over the shorter way between them.
fn mitre(
    model: &mut Model,
    band: &Shape,
    blend: &Shape,
    (a, pa): (&Shape, Point),
    (b, pb): (&Shape, Point),
    tol: Tolerances,
) -> OgeomResult<Option<Run>> {
    let (Some((_, own)), Some((_, theirs))) = (
        ring::surface_of(model, band),
        ring::surface_of(model, blend),
    ) else {
        return Ok(None);
    };
    let Ok(ogeom_intersect::Meeting::Along(curves)) =
        ogeom_intersect::surface_surface(&own, &theirs, tol)
    else {
        return Ok(None);
    };
    let tau = core::f64::consts::TAU;
    for curve in curves {
        let (Some(ta), Some(tb)) = (ellipse_at(&curve, pa), ellipse_at(&curve, pb)) else {
            continue;
        };
        let on = |t: f64, p: Point| {
            ogeom_geom::Curve3d::point_at(&curve, t, tol)
                .is_ok_and(|q| q.distance(p) <= tol.confusion())
        };
        if !(on(ta, pa) && on(tb, pb)) {
            continue;
        }
        let turn = (tb - ta + core::f64::consts::PI).rem_euclid(tau) - core::f64::consts::PI;
        let start = ta.rem_euclid(tau);
        let range = (start, start + turn);
        let shift = (range.0.min(range.1) / tau).floor() * tau;
        return run_between(model, curve, (range.0 - shift, range.1 - shift), a, b, tol).map(Some);
    }
    Ok(None)
}
