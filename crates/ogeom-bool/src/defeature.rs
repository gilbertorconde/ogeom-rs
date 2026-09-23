//! Defeaturing by face removal: delete faces, close the wound from the
//! neighbours' own surfaces.
//!
//! The input is a set of faces; what those faces *mean* is the caller's
//! business, and the operation works on a solid whose history is gone.
//! Three wounds exist, and they close differently.
//!
//! A feature whose rim is an **inner loop** of a surviving face (a bore in a
//! lid, a boss on a base, a pocket in the middle of a top) leaves survivors
//! whose boundary is already right except for that loop. The cure is wire
//! surgery: the surviving face is rebuilt without the rim wire, edges,
//! pcurves and all, and nothing is re-intersected because nothing new meets.
//!
//! A feature that **interrupts** its neighbours' outer boundaries (a fillet
//! band or a chamfer along an edge) leaves a gap no surviving boundary
//! closes. The cure is the neighbours themselves: the two side faces'
//! surfaces are re-intersected to recover the edge the blend replaced, the
//! end faces' edges are extended along their own curves to the recovered
//! corners, and the faces are rebuilt on the result. Extension here is the
//! surfaces' and curves' own unbounded carriers: no new geometry is
//! invented, only wider windows of what is already there.
//!
//! A feature that takes a **whole ring** out of a neighbour (a rim
//! blend, round a drum's top, a bore's mouth or a boss's seat) looks
//! like the first wound and closes like the second. The neighbours' own
//! surfaces tell the two apart: a bore's two mouths sit in faces that
//! never meet, so the rings are dropped and the faces grow over them,
//! while a rim blend's cap and wall meet along the very circle it
//! replaced, in the wound's own room. There the ring is replaced rather
//! than dropped, and a neighbour's outer boundary may be the ring: a
//! drum's cap grows back to its own rim. The wall's chart has a seam, and
//! the seam reaches the recovered circle: that is where the circle is
//! cut, and the seam extends to meet it, exactly as a band's end faces
//! extend to their corners. One corner leaves the rim one closed edge,
//! re-anchored so the whole turn stands in the curve's own domain: the
//! shape the rim had before the feature was cut.
//!
//! Several bands close together. Each removed band recovers its own
//! crease; where two creases meet (two blends that met at a corner, or
//! one blend's flush cap standing against another's band, the cap named
//! with its band), the corner is where one crease pierces the other's
//! side, and it is one vertex for both.
//!
//! Every rebuilt wire is spliced in the face's own order rather than
//! re-chained from a bag of edges: a chart's seam stands in its wire
//! twice, and a bag cannot say so. Where a gap leaves and arrives at one
//! vertex, the rim it replaces says which way round it goes; nothing in
//! the topology notices a face inside out along its own rim, and the
//! mesher finds it as a boundary that will not close.
//!
//! What this does not yet close is refused by name: a wound whose sides
//! do not meet in a curve, a removal that would leave a face with no
//! boundary and no edge to grow to, and a gap the recovered edges do not
//! bridge.

use crate::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_algo::{Built, History, make_edge_between, make_solid, make_vertex, sew};
use ogeom_core::ogeom_err;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_intersect::{
    CurveSurfaceOptions, IntersectOptions, SurfaceIntersection, intersect_curve_surface,
    intersect_surfaces,
};
use ogeom_math::Point;
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, TShapeId, explore};
use std::collections::{HashMap, HashSet};

/// Remove `faces` from `solid` and close the openings from the neighbours'
/// own geometry.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, when the removal is not one this operation closes: no face named,
/// every face named, a named shape that is not a face of the solid, a wound
/// whose side surfaces do not meet in a single curve, more than one band, or
/// geometry whose pcurves have no closed form to rebuild with.
pub fn remove_faces(
    model: &mut Model,
    solid: &Shape,
    faces: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Built> {
    if faces.is_empty() {
        ogeom_bail!(Construction, "no faces named; there is nothing to remove");
    }
    // Separate features remove separately. Two bores named in one call are
    // two wounds; classifying their ring edges together declares the two
    // longest interrupted faces "the sides" across both and recovers a
    // nonsense edge. Faces group into features by shared edges, and each
    // feature runs the whole machinery on the previous feature's result,
    // sequential exactly as a caller would have called it, so one call
    // means what N calls mean, in the order given.
    let groups = feature_groups(model, faces)?;
    if groups.len() > 1 {
        let mut current = Built::from_nothing(solid.clone());
        for group in &groups {
            // A later group's faces survive the earlier surgeries untouched
            // (different regions), but the solid they belong to is new.
            let step = remove_faces(model, &current.shape, group, tol)?;
            current = Built {
                shape: step.shape,
                history: current.history.then(&step.history),
            };
        }
        return Ok(current);
    }
    let all_faces = explore(model, solid, Filter::OfType(ShapeType::Face))?;
    let removed: HashSet<TShapeId> = faces.iter().map(Shape::node).collect();
    for face in faces {
        if !all_faces.iter().any(|f| f.node() == face.node()) {
            ogeom_bail!(
                Construction,
                "a face named for removal is not a face of this solid"
            );
        }
    }
    let survivors: Vec<Shape> = all_faces
        .iter()
        .filter(|f| !removed.contains(&f.node()))
        .cloned()
        .collect();
    if survivors.is_empty() {
        ogeom_bail!(
            Construction,
            "every face was named for removal; nothing remains to close"
        );
    }

    // Which edges the removed set shares with the world: an edge is a ring
    // edge when a removed face and a surviving face both use it.
    let mut users: HashMap<TShapeId, Vec<Shape>> = HashMap::new();
    for face in &all_faces {
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            users.entry(edge.node()).or_default().push(face.clone());
        }
    }
    let is_ring = |edge: &Shape| -> bool {
        users.get(&edge.node()).is_some_and(|fs| {
            fs.iter().any(|f| removed.contains(&f.node()))
                && fs.iter().any(|f| !removed.contains(&f.node()))
        })
    };

    // Sort survivors: untouched, whole-ring, and interrupted: a wire with
    // no rim edge, a wire that is all rim, a wire that is part rim.
    struct Touched {
        face: Shape,
        /// The wires with no rim edge at all, which stand either way.
        kept: Vec<Shape>,
        /// Some wire of this face is all rim.
        whole: bool,
        /// Some wire of this face is part rim.
        partial: bool,
        /// And the outer wire is one of the whole ones.
        outer: bool,
    }
    let mut untouched: Vec<Shape> = Vec::new();
    let mut touched: Vec<Touched> = Vec::new();
    for face in &survivors {
        let wires = model.ordered_children_of(face)?;
        let mut kept = Vec::new();
        let (mut whole, mut partial, mut outer) = (false, false, false);
        for (index, wire) in wires.iter().enumerate() {
            let edges = model.ordered_children_of(wire)?;
            let ring_count = edges.iter().filter(|e| is_ring(e)).count();
            if ring_count == 0 {
                kept.push(wire.clone());
            } else if ring_count == edges.len() {
                whole = true;
                outer |= index == 0;
            } else {
                partial = true;
            }
        }
        if whole || partial {
            touched.push(Touched {
                face: face.clone(),
                kept,
                whole,
                partial,
                outer,
            });
        } else {
            untouched.push(face.clone());
        }
    }

    // A neighbour that lost a whole ring is closed one of two ways, and the
    // neighbours' surfaces say which: where they meet in the wound's own
    // room the ring is replaced by the edge they meet along, and where they
    // do not meet at all it is simply dropped and the face grows over what
    // the feature stood in.
    let sides: Vec<Shape> = touched.iter().map(|t| t.face.clone()).collect();
    let recovers = touched.iter().any(|t| t.whole) && wound_recovers(model, faces, &sides, tol)?;
    let mut rim_surgery: Vec<(Shape, Vec<Shape>)> = Vec::new(); // face, kept wires
    let mut interrupted: Vec<Shape> = Vec::new();
    for entry in touched {
        if entry.partial || recovers {
            interrupted.push(entry.face);
        } else {
            // Dropping the outer boundary would leave a face with nothing
            // to stand on; replacing it, where the neighbours meet, is the
            // branch above.
            if entry.outer {
                ogeom_bail!(
                    Construction,
                    "removing these faces erases a neighbour's whole outer \
                     boundary; that face has nothing left to stand on"
                );
            }
            rim_surgery.push((entry.face, entry.kept));
        }
    }

    let mut history = History::new();
    for face in faces {
        history.delete(face);
    }

    let mut rebuilt: Vec<Shape> = untouched;
    for (face, kept_wires) in rim_surgery {
        let new_face = {
            let Some(data) = model.node(&face).and_then(|n| match n.data() {
                NodeData::Face(d) => Some(d.clone()),
                _ => None,
            }) else {
                ogeom_bail!(Construction, "a surviving face holds no face data");
            };
            // The kept wires carry their edges, and the edges their pcurves
            // for this very surface: nothing to recompute.
            let built = ogeom_algo::make_face_on(model, data.surface, &kept_wires, tol)?.shape;
            orient_like(&face, built)
        };
        history.modify(&face, new_face.clone());
        rebuilt.push(new_face);
    }

    if !interrupted.is_empty() {
        let band = close_wound(model, faces, &interrupted, &removed, &users, &is_ring, tol)?;
        for (old, new) in band {
            history.modify(&old, new.clone());
            rebuilt.push(new);
        }
    }

    let sewn = sew(model, &rebuilt, tol)?;
    let [shell] = sewn.shells.as_slice() else {
        ogeom_bail!(
            Construction,
            "closing the wound left {} shells where one solid's worth was \
             expected; the removal disconnected the boundary",
            sewn.shells.len()
        );
    };
    if !ogeom_algo::is_shell_closed(model, shell)? {
        ogeom_bail!(
            Construction,
            "the boundary does not close after removal; the wound needs a \
             closure this operation does not construct yet"
        );
    }
    let built = make_solid(model, std::slice::from_ref(shell))?;
    let mut solid_history = history;
    solid_history.modify(solid, built.shape.clone());
    Ok(Built::new(built.shape, solid_history))
}

/// A face's surface, carried into space by the face's own placement.
fn placed_surface(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let placement = face.transform(model.datums())?;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face().cloned()) else {
        ogeom_bail!(Construction, "a band face holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Construction, "a band face's surface is not in this model");
    };
    surface.clone().transformed(&placement, tol)
}

/// Whether the wound's neighbours meet each other where it sat.
///
/// A whole ring taken out of a neighbour is two different wounds, and only
/// the neighbours' own surfaces tell them apart. A bore's wall leaves its
/// two mouths as whole inner wires, and the faces holding them (a block's
/// top and bottom) never meet: dropping the wires is the closure, and the
/// block comes back whole. A rim blend leaves a whole ring too (the
/// annulus a mouth fillet takes out of the top, the circle a boss's seat
/// takes out of the wall), but there the cap and the wall meet along the
/// very circle the blend replaced, and dropping would leave the boundary
/// open where that circle belongs.
///
/// Meeting *somewhere* is not enough: two faces of any solid meet if their
/// surfaces are carried far enough, and a bore through a wedge would
/// recover the line where the wedge closes. The meeting must stand in the
/// wound's own room (the removed faces' bounds), which is where the edge
/// the feature replaced stood.
fn wound_recovers(
    model: &Model,
    removed_faces: &[Shape],
    candidates: &[Shape],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mut room = ogeom_math::Aabb::default();
    for face in removed_faces {
        room = room.union(&ogeom_algo::shape_bounds(model, face, tol)?);
    }
    let room = room.expanded(tol.confusion() * 1e3);
    let Some(centre) = room.centre() else {
        return Ok(false);
    };
    for (i, first) in candidates.iter().enumerate() {
        let sa = placed_surface(model, first, tol)?;
        for second in &candidates[i + 1..] {
            let sb = placed_surface(model, second, tol)?;
            let Ok(SurfaceIntersection::Along(sections)) =
                intersect_surfaces(&sa, &sb, IntersectOptions::default(), tol)
            else {
                continue;
            };
            for section in sections {
                let foot = ogeom_algo::project_on_curve(&section.curve, centre, 64, tol)?;
                if room.contains(foot.point) {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// One crease the wound recovers: the edge a removed band replaced, from
/// its two side faces' own surfaces.
struct Crease {
    /// The side faces, by node, in a fixed order.
    sides: [Shape; 2],
    /// The recovered curve, the branch nearest the removed faces.
    curve: Curve,
    /// The removed faces' extent along the curve, unwrapped about the
    /// anchor on a periodic curve.
    extent: (f64, f64),
    /// The parameter nearest the removed faces' centre: what a periodic
    /// curve's parameters are unwrapped about, so a band straddling the
    /// curve's seam reads as one run and not its complement.
    anchor: f64,
}

/// A periodic curve's parameter brought within half a period of `about`;
/// any other curve's parameter as it is.
fn unwrapped(curve: &Curve, t: f64, about: f64) -> f64 {
    if !curve.is_periodic() {
        return t;
    }
    let (lo, hi) = curve.domain();
    let period = hi - lo;
    if period <= 0.0 {
        return t;
    }
    about + (t - about + period / 2.0).rem_euclid(period) - period / 2.0
}

/// Close a wound: each removed band's two side faces re-intersected into
/// the crease it replaced, the creases' ends placed where they pierce the
/// other interrupted faces (or one another's sides, which is where two
/// bands meeting at a corner share their corner), every dangling edge
/// extended along its own curve to the corner standing on it, and every
/// interrupted face rebuilt with the creases it borders.
///
/// A band's sides are the two survivors it shares the most ring length
/// with; a wedge's cap named alongside its band shares the band's sides
/// and folds into the same crease. A crease's ends are the nearest
/// piercings just past the removed faces' own extent along it, so a
/// survivor the curve merely runs through far away is not mistaken for an
/// end. Corners are one vertex wherever two creases place them within
/// tolerance of each other.
#[allow(clippy::too_many_lines, reason = "one wound, one narrative")]
fn close_wound(
    model: &mut Model,
    removed_faces: &[Shape],
    interrupted: &[Shape],
    removed: &HashSet<TShapeId>,
    users: &HashMap<TShapeId, Vec<Shape>>,
    is_ring: &dyn Fn(&Shape) -> bool,
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Shape)>> {
    let surface_of = |model: &Model, face: &Shape| placed_surface(model, face, tol);
    let vertices_of = |model: &Model, face: &Shape| -> OgeomResult<Vec<Point>> {
        let mut out = Vec::new();
        for vertex in explore(model, face, Filter::OfType(ShapeType::Vertex))? {
            let placement = vertex.transform(model.datums())?;
            if let Some(d) = model.node(&vertex).and_then(|nd| nd.data().as_vertex()) {
                out.push(placement.apply(d.point));
            }
        }
        Ok(out)
    };
    let interrupted_by_node: HashMap<TShapeId, Shape> =
        interrupted.iter().map(|f| (f.node(), f.clone())).collect();

    // Each removed face's sides: the two interrupted survivors it shares
    // the most ring length with. Creases are keyed by the side pair. A
    // removed face with fewer than two such neighbours (a wedge's cap
    // standing against another blend's band, bordering one wall and two
    // removed faces) joins the crease of a removed neighbour it shares an
    // edge with, once that neighbour has one.
    let mut creases: Vec<(TShapeId, TShapeId, [Shape; 2], Vec<Point>)> = Vec::new();
    let mut crease_of: HashMap<TShapeId, usize> = HashMap::new();
    let mut leftovers: Vec<Shape> = Vec::new();
    for face in removed_faces {
        let mut shared: HashMap<TShapeId, f64> = HashMap::new();
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            if !is_ring(&edge) {
                continue;
            }
            let length = edge_length(model, &edge, tol)?;
            for user in users.get(&edge.node()).into_iter().flatten() {
                if !removed.contains(&user.node()) && interrupted_by_node.contains_key(&user.node())
                {
                    *shared.entry(user.node()).or_default() += length;
                }
            }
        }
        let mut ranked: Vec<(TShapeId, f64)> = shared.into_iter().collect();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.index().cmp(&b.0.index())));
        let [(a, _), (b, _), ..] = ranked.as_slice() else {
            leftovers.push(face.clone());
            continue;
        };
        let (lo, hi) = if a.index() <= b.index() {
            (*a, *b)
        } else {
            (*b, *a)
        };
        let points = vertices_of(model, face)?;
        let index = match creases.iter().position(|c| c.0 == lo && c.1 == hi) {
            Some(i) => {
                creases[i].3.extend(points);
                i
            }
            None => {
                creases.push((
                    lo,
                    hi,
                    [
                        interrupted_by_node[&lo].clone(),
                        interrupted_by_node[&hi].clone(),
                    ],
                    points,
                ));
                creases.len() - 1
            }
        };
        crease_of.insert(face.node(), index);
    }
    for face in leftovers {
        let mut joined = None;
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            for user in users.get(&edge.node()).into_iter().flatten() {
                if let Some(&index) = crease_of.get(&user.node()) {
                    joined = Some(index);
                }
            }
        }
        let Some(index) = joined else {
            ogeom_bail!(
                Construction,
                "a removed face shares ring edges with fewer than two \
                 interrupted neighbours and borders no removed face with a \
                 crease; closing it needs a neighbour to meet itself, which \
                 is not constructed yet"
            );
        };
        let points = vertices_of(model, &face)?;
        creases[index].3.extend(points);
    }

    // The recovered curve of each crease, and the removed faces' extent on it.
    let mut recovered: Vec<Crease> = Vec::new();
    for (_, _, sides, points) in creases {
        let sa = surface_of(model, &sides[0])?;
        let sb = surface_of(model, &sides[1])?;
        let meeting = intersect_surfaces(&sa, &sb, IntersectOptions::default(), tol)?;
        let SurfaceIntersection::Along(sections) = meeting else {
            ogeom_bail!(
                Construction,
                "the band's side surfaces do not meet along a curve; the edge \
                 the feature replaced cannot be recovered from them"
            );
        };
        let anchor = {
            let mut sum = ogeom_math::Vector::ZERO;
            for p in &points {
                sum += p.to_vector();
            }
            #[allow(clippy::cast_precision_loss)]
            let n = points.len().max(1) as f64;
            Point::ORIGIN + sum * (1.0 / n)
        };
        let section = sections
            .into_iter()
            .min_by(|p, q| {
                nearest_distance(&p.curve, anchor, tol)
                    .total_cmp(&nearest_distance(&q.curve, anchor, tol))
            })
            .ok_or_else(|| ogeom_err!(Construction, "the side surfaces meet along no branch"))?;
        let curve = section.curve;
        let anchor = parameter_near(&curve, anchor, tol)?;
        let mut extent = (f64::INFINITY, f64::NEG_INFINITY);
        for p in &points {
            let t = unwrapped(&curve, parameter_near(&curve, *p, tol)?, anchor);
            extent = (extent.0.min(t), extent.1.max(t));
        }
        recovered.push(Crease {
            sides,
            curve,
            extent,
            anchor,
        });
    }

    // Corners: where each crease pierces an interrupted face that is not
    // one of its sides, the nearest piercing past each end of its extent.
    // A shared corner is one vertex.
    let mut corner_vertices: Vec<(Point, Shape)> = Vec::new();
    let mut vertex_at = |model: &mut Model, p: Point| -> Shape {
        if let Some((_, v)) = corner_vertices
            .iter()
            .find(|(q, _)| q.distance(p) <= tol.confusion() * 1e3)
        {
            return v.clone();
        }
        let v = make_vertex(model, p).shape;
        corner_vertices.push((p, v.clone()));
        v
    };
    let mut new_edges: Vec<(Shape, [TShapeId; 2])> = Vec::new(); // edge, its sides
    let mut corners: Vec<Shape> = Vec::new();
    for (index, crease) in recovered.iter().enumerate() {
        let mut piercings: Vec<(f64, Point)> = Vec::new();
        for face in interrupted {
            if crease.sides.iter().any(|s| s.node() == face.node()) {
                continue;
            }
            let se = surface_of(model, face)?;
            let hit =
                intersect_curve_surface(&crease.curve, &se, CurveSurfaceOptions::default(), tol)?;
            for c in &hit.crossings {
                piercings.push((unwrapped(&crease.curve, c.on_curve, crease.anchor), c.point));
            }
        }
        // A tangent junction: two bands of one chain meeting flush (a
        // stadium's straight run into its semicircular end) share the
        // cross-section edge where they meet, and their creases touch
        // there without either piercing the other's side; a wall the
        // crease merely grazes yields no piercing, and a touch found as a
        // closest approach sits anywhere in a valley the width of the
        // slop. The shared edge says exactly where: the cross-section
        // stands in the plane normal to the rim at the junction, so the
        // junction is the foot of that edge on either crease: a
        // transversal projection, exact to the last bit. Taken only where
        // the two creases are tangent there; bands meeting at a corner
        // place theirs by piercing.
        for face in removed_faces {
            if crease_of.get(&face.node()) != Some(&index) {
                continue;
            }
            for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
                let Some(other) = users
                    .get(&edge.node())
                    .into_iter()
                    .flatten()
                    .filter(|user| user.node() != face.node())
                    .filter_map(|user| crease_of.get(&user.node()).copied())
                    .find(|&other| other != index)
                else {
                    continue;
                };
                let samples = sample_edge(model, &edge, tol)?;
                let Some(&middle) = samples.get(samples.len() / 2) else {
                    continue;
                };
                let here = ogeom_algo::project_on_curve(&crease.curve, middle, 64, tol)?;
                let there = ogeom_algo::project_on_curve(&recovered[other].curve, middle, 64, tol)?;
                if here.point.distance(there.point) > tol.confusion() * 1e3 {
                    continue;
                }
                let ta = crease.curve.d1_at(here.parameter, tol)?;
                let tb = recovered[other].curve.d1_at(there.parameter, tol)?;
                if ta.cross(tb).magnitude() > ta.magnitude() * tb.magnitude() * 1e-3 {
                    continue;
                }
                piercings.push((
                    unwrapped(&crease.curve, here.parameter, crease.anchor),
                    here.point,
                ));
            }
        }
        let slack = tol.parametric().max(1e-6);
        let below = piercings
            .iter()
            .filter(|(t, _)| *t <= crease.extent.0 + slack)
            .max_by(|a, b| a.0.total_cmp(&b.0));
        let above = piercings
            .iter()
            .filter(|(t, _)| *t >= crease.extent.1 - slack)
            .min_by(|a, b| a.0.total_cmp(&b.0));
        match (below, above) {
            (None, None) if piercings.is_empty() => {
                // No ends: the band wraps, and the recovered edge closes on
                // itself, but not always as one edge. A chart's seam
                // reaches the recovered curve too: a cylinder wall's wire
                // runs up its seam, round the rim and back down, and a
                // closed edge carrying a vertex of its own leaves that wire
                // two chains that never meet. So each side's own dangling
                // boundary names a corner where it reaches the curve, the
                // curve is cut there, and the dangling edge then extends to
                // it like any other.
                let (lo, hi) = crease.curve.domain();
                if !crease.curve.is_periodic() {
                    ogeom_bail!(
                        Construction,
                        "a wrapping band recovered an open curve; the closure is \
                         not constructible from it"
                    );
                }
                let period = hi - lo;
                let mut stops: Vec<f64> = Vec::new();
                for side in &crease.sides {
                    let mut rim_vertices: HashSet<TShapeId> = HashSet::new();
                    for edge in explore(model, side, Filter::OfType(ShapeType::Edge))? {
                        if is_ring(&edge) {
                            for v in model.ordered_children_of(&edge)? {
                                rim_vertices.insert(v.node());
                            }
                        }
                    }
                    for edge in explore(model, side, Filter::OfType(ShapeType::Edge))? {
                        if is_ring(&edge) {
                            continue;
                        }
                        let free: Vec<Point> = model
                            .ordered_children_of(&edge)?
                            .iter()
                            .filter(|v| rim_vertices.contains(&v.node()))
                            .filter_map(|v| {
                                model
                                    .node(v)
                                    .and_then(|n| n.data().as_vertex())
                                    .map(|d| d.point)
                            })
                            .collect();
                        if free.is_empty() {
                            continue;
                        }
                        let Some(geometry) = edge_geometry(model, &edge, tol)? else {
                            continue;
                        };
                        for start in free {
                            // The two curves' closest approach, from the
                            // dangling end: each in turn answers where the
                            // other's nearest point is, and an edge that
                            // genuinely reaches the curve settles on it.
                            let mut point = start;
                            for _ in 0..8 {
                                let on_curve =
                                    ogeom_algo::project_on_curve(&crease.curve, point, 64, tol)?;
                                let t = parameter_near(&geometry, on_curve.point, tol)?;
                                let Ok(on_edge) = geometry.point_at(t, tol) else {
                                    break;
                                };
                                if on_edge.distance(on_curve.point) <= tol.confusion() * 10.0 {
                                    stops.push(lo + (on_curve.parameter - lo).rem_euclid(period));
                                    break;
                                }
                                point = on_edge;
                            }
                        }
                    }
                }
                // Re-anchored at the first corner, so the whole turn stands
                // inside the curve's own domain: a rim written from a corner
                // right round to itself would otherwise end a turn past the
                // end of it. One corner then leaves the rim one closed edge,
                // which is the shape it had before the feature was cut,
                // and the shape the exact volume integrator reads as a disc.
                let (curve, stops) = match (&crease.curve, stops.first().copied()) {
                    (Curve::Circle(circle), Some(first)) => {
                        let at = crease.curve.point_at(first, tol)?;
                        let held = circle.circle();
                        let anchored = ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
                            ogeom_math::Frame::new(
                                held.centre(),
                                held.frame().z(),
                                ogeom_math::Direction::new(at - held.centre(), tol)?,
                                tol,
                            )?,
                            held.radius(),
                            tol,
                        )?);
                        let mut anchored = Curve::Circle(anchored);
                        // The same circle, and the same way round it: the
                        // re-anchoring moves where the parameter starts and
                        // must not turn the rim over.
                        if anchored
                            .d1_at(0.0, tol)?
                            .dot(crease.curve.d1_at(first, tol)?)
                            < 0.0
                        {
                            anchored = ogeom_geom::Reversible::reversed(&anchored);
                        }
                        let shifted = stops
                            .iter()
                            .map(|t| (t - first).rem_euclid(period))
                            .collect::<Vec<f64>>();
                        (anchored, shifted)
                    }
                    // A rim that is not a circle cannot be re-anchored, so
                    // it is cut at the chart's start as well and the turn
                    // stays inside the domain in pieces instead.
                    _ => {
                        let mut kept = stops;
                        kept.push(lo);
                        (crease.curve.clone(), kept)
                    }
                };
                let mut cuts = stops;
                cuts.sort_by(f64::total_cmp);
                cuts.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric().max(1e-9));
                if cuts.is_empty() {
                    cuts.push(lo);
                }
                let mut cut: Vec<Shape> = Vec::with_capacity(cuts.len());
                for t in &cuts {
                    let at = curve.point_at(*t, tol)?;
                    cut.push(vertex_at(model, at));
                }
                for (index, t) in cuts.iter().enumerate() {
                    let next = if index + 1 == cuts.len() {
                        cuts[0] + period
                    } else {
                        cuts[index + 1]
                    };
                    let edge = make_edge_between(
                        model,
                        curve.clone(),
                        (*t, next),
                        &cut[index],
                        &cut[(index + 1) % cuts.len()],
                        tol,
                    )?
                    .shape;
                    new_edges.push((edge, [crease.sides[0].node(), crease.sides[1].node()]));
                }
                corners.extend(cut);
            }
            (Some(&(t0, p0)), Some(&(t1, p1))) => {
                let v0 = vertex_at(model, p0);
                let v1 = vertex_at(model, p1);
                // Unwrapped about the anchor, a window can start before a
                // periodic curve's domain; slid by whole turns to start
                // inside it, it is the same run, and may end a turn past
                // the end as any run across the seam does.
                let window = if crease.curve.is_periodic() {
                    let (lo, hi) = crease.curve.domain();
                    let turns = ((t0 - lo) / (hi - lo)).floor() * (hi - lo);
                    (t0 - turns, t1 - turns)
                } else {
                    (t0, t1)
                };
                let edge =
                    make_edge_between(model, crease.curve.clone(), window, &v0, &v1, tol)?.shape;
                corners.push(v0);
                corners.push(v1);
                new_edges.push((edge, [crease.sides[0].node(), crease.sides[1].node()]));
            }
            _ => ogeom_bail!(
                Construction,
                "an end face's surface never meets the recovered edge; the \
                 corner cannot be placed"
            ),
        }
    }

    let mut out = Vec::new();
    let mut extended: HashMap<TShapeId, Shape> = HashMap::new();
    for face in interrupted {
        let rims: Vec<Shape> = explore(model, face, Filter::OfType(ShapeType::Edge))?
            .into_iter()
            .filter(|e| is_ring(e))
            .collect();
        let borders: Vec<Shape> = new_edges
            .iter()
            .filter(|(_, sides)| sides.contains(&face.node()))
            .map(|(e, _)| e.clone())
            .collect();
        let new_face =
            rebuild_interrupted(model, face, &rims, &borders, &corners, &mut extended, tol)?;
        out.push((face.clone(), new_face));
    }
    Ok(out)
}

/// The connected components of the removal set: faces joined by shared
/// edges belong to one feature and close as one wound.
fn feature_groups(model: &Model, faces: &[Shape]) -> OgeomResult<Vec<Vec<Shape>>> {
    let mut edge_sets: Vec<HashSet<TShapeId>> = Vec::with_capacity(faces.len());
    for face in faces {
        edge_sets.push(
            explore(model, face, Filter::OfType(ShapeType::Edge))?
                .iter()
                .map(Shape::node)
                .collect(),
        );
    }
    let mut group_of: Vec<usize> = (0..faces.len()).collect();
    // Union by scan: small sets, clarity over asymptotics.
    fn root(group_of: &mut [usize], mut i: usize) -> usize {
        while group_of[i] != i {
            group_of[i] = group_of[group_of[i]];
            i = group_of[i];
        }
        i
    }
    for i in 0..faces.len() {
        for j in i + 1..faces.len() {
            if edge_sets[i].intersection(&edge_sets[j]).next().is_some() {
                let (a, b) = (root(&mut group_of, i), root(&mut group_of, j));
                group_of[a.max(b)] = a.min(b);
            }
        }
    }
    let mut groups: HashMap<usize, Vec<Shape>> = HashMap::new();
    for (i, face) in faces.iter().enumerate() {
        groups
            .entry(root(&mut group_of, i))
            .or_default()
            .push(face.clone());
    }
    let mut out: Vec<Vec<Shape>> = groups.into_values().collect();
    // Deterministic order: by each group's smallest node index.
    out.sort_by_key(|g| g.iter().map(|f| f.node().index()).min());
    Ok(out)
}

/// Rebuild one interrupted face: drop its ring edges, extend the edges that
/// now dangle to the corner vertex standing on their own curve, add the
/// recovered edges this face borders, and rechain.
fn rebuild_interrupted(
    model: &mut Model,
    face: &Shape,
    rims: &[Shape],
    borders: &[Shape],
    corners: &[Shape],
    extended: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let placement = face.transform(model.datums())?;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face().cloned()) else {
        ogeom_bail!(Construction, "an interrupted face holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(
            Construction,
            "an interrupted face's surface is not in this model"
        );
    };
    let surface = surface.transformed(&placement, tol)?;
    let rim_nodes: HashSet<TShapeId> = rims.iter().map(Shape::node).collect();

    let mut wires: Vec<Vec<Shape>> = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let edges = model.ordered_children_of(&wire)?;
        let touched = edges.iter().any(|e| rim_nodes.contains(&e.node()));
        if !touched {
            wires.push(edges);
            continue;
        }
        // Which vertices the dropped rim owned: an edge that shared one now
        // dangles there and must reach a corner instead.
        let mut rim_vertices: HashSet<TShapeId> = HashSet::new();
        for edge in &edges {
            if rim_nodes.contains(&edge.node()) {
                for v in model.ordered_children_of(edge)? {
                    rim_vertices.insert(v.node());
                }
            }
        }
        // Substituted in the wire's own order rather than re-chained from
        // a bag of edges. A chart's seam stands in its wire *twice* (up
        // one column and down the other), and a bag cannot say so: four
        // edge ends meet at each of the seam's vertices, which reads as a
        // branching network and not a wire. The order the face already has
        // is the answer the bag was being asked to guess.
        let mut ring: Vec<(Shape, Option<Shape>)> = Vec::with_capacity(edges.len());
        for edge in &edges {
            if rim_nodes.contains(&edge.node()) {
                ring.push((edge.clone(), None));
                continue;
            }
            // A face that has already extended this edge decided for
            // everyone; sewing rejoins on the shared node.
            let replaced = match extended.get(&edge.node()) {
                Some(found) => Some(found.clone()),
                None => {
                    let dangles = model
                        .ordered_children_of(edge)?
                        .iter()
                        .any(|v| rim_vertices.contains(&v.node()));
                    match (dangles, corner_on_edge(model, edge, corners, tol)?) {
                        (true, Some(corner)) => {
                            Some(extend_to_corner(model, edge, &corner, extended, tol)?)
                        }
                        _ => None,
                    }
                }
            };
            // A replacement is built forward; the wire's own use decides
            // which way it runs here.
            ring.push((
                edge.clone(),
                Some(match replaced {
                    Some(fresh) if edge.orientation() == ogeom_topo::Orientation::Reversed => {
                        fresh.reversed()
                    }
                    Some(fresh) => fresh,
                    None => edge.clone(),
                }),
            ));
        }
        if std::env::var_os("OGEOM_DEBUG_DEFEATURE").is_some() {
            let point = |v: &Shape| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| d.point)
            };
            for (edge, spliced) in &ring {
                let (a, b) = edge_ends(model, spliced.as_ref().unwrap_or(edge))?;
                eprintln!(
                    "DEFEATURE wire edge {} {} {:?} .. {:?}",
                    edge.node().index(),
                    if spliced.is_none() { "RIM" } else { "kept" },
                    point(&a),
                    point(&b)
                );
            }
        }
        let mut pool: Vec<Shape> = borders.to_vec();
        let mut chained: Vec<Shape> = Vec::new();
        if ring.iter().all(|(_, spliced)| spliced.is_none()) {
            // The whole wire was the wound's rim: the recovered edges are
            // the wire, closing on themselves.
            let Some(start) = pool.first().cloned() else {
                ogeom_bail!(
                    Construction,
                    "a neighbour lost a whole ring and no recovered edge \
                     borders it; the wound needs a closure this operation \
                     does not construct yet"
                );
            };
            let (from, _) = edge_ends(model, &start)?;
            chained = bridge_gap(model, &mut pool, &from, &from)?;
            let mut was = Vec::new();
            for (edge, _) in &ring {
                was.extend(sample_edge(model, edge, tol)?);
            }
            wind_like(model, &mut chained, &was, tol)?;
        } else {
            // Rotated so the wire begins on an edge that survived, which
            // puts every run of rim edges between two of them.
            let first = ring
                .iter()
                .position(|(_, spliced)| spliced.is_some())
                .unwrap_or(0);
            ring.rotate_left(first);
            let mut index = 0;
            while index < ring.len() {
                if let Some(edge) = ring[index].1.clone() {
                    chained.push(edge);
                    index += 1;
                    continue;
                }
                let run_end = ring[index..]
                    .iter()
                    .position(|(_, spliced)| spliced.is_some())
                    .map_or(ring.len(), |k| index + k);
                let Some(previous) = chained.last() else {
                    ogeom_bail!(Construction, "a wound's rim opens a wire that has no start");
                };
                let (_, from) = edge_ends(model, previous)?;
                let to = match ring.get(run_end).and_then(|(_, spliced)| spliced.as_ref()) {
                    Some(next) => edge_ends(model, next)?.0,
                    // The run closes the ring: it comes back to the start.
                    None => edge_ends(model, &chained[0])?.0,
                };
                let mut bridge = bridge_gap(model, &mut pool, &from, &to)?;
                // A gap that leaves and arrives at one vertex could be
                // walked either way round; the rim it replaces says which.
                if from.node() == to.node() {
                    let mut was = Vec::new();
                    for (edge, _) in &ring[index..run_end] {
                        was.extend(sample_edge(model, edge, tol)?);
                    }
                    wind_like(model, &mut bridge, &was, tol)?;
                }
                chained.extend(bridge);
                index = run_end;
            }
        }
        if !pool.is_empty() {
            ogeom_bail!(
                Construction,
                "{} recovered edges border this face and its wound's rim has \
                 nowhere to put them",
                pool.len()
            );
        }
        wires.push(chained);
    }
    let built = ogeom_algo::make_face_with_pcurves(model, surface, &wires, tol)?.shape;
    // The face's own side of its surface, carried over. A rebuilt face is
    // born forward, and a bore's wall is not: the mesher reads the wires'
    // winding and forgives it, but the exact integrator reads the flag and
    // hands back the bore as material.
    Ok(orient_like(face, built))
}

/// A rebuilt face put back on the side of its surface the old one was on.
fn orient_like(was: &Shape, built: Shape) -> Shape {
    if was.orientation() == ogeom_topo::Orientation::Reversed {
        built.reversed()
    } else {
        built
    }
}

/// Points along an edge, in the direction this use of it runs.
fn sample_edge(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    const STATIONS: usize = 12;
    let Some(geometry) = edge_geometry(model, edge, tol)? else {
        return Ok(Vec::new());
    };
    let Some(range) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { range, .. } => Some(*range),
            _ => None,
        })
    else {
        return Ok(Vec::new());
    };
    let backwards = edge.orientation() == ogeom_topo::Orientation::Reversed;
    let mut out = Vec::with_capacity(STATIONS + 1);
    for i in 0..=STATIONS {
        #[allow(clippy::cast_precision_loss)]
        let f = i as f64 / STATIONS as f64;
        let f = if backwards { 1.0 - f } else { f };
        out.push(geometry.point_at(range.0 + (range.1 - range.0) * f, tol)?);
    }
    Ok(out)
}

/// Twice the area a closed run of points sweeps about its own centre, as a
/// vector: which way round the run goes, in the only terms two runs of
/// different shapes can be compared in.
fn swept_area(points: &[Point]) -> ogeom_math::Vector {
    if points.len() < 3 {
        return ogeom_math::Vector::ZERO;
    }
    let mut sum = ogeom_math::Vector::ZERO;
    for p in points {
        sum += p.to_vector();
    }
    #[allow(clippy::cast_precision_loss)]
    let centre = Point::ORIGIN + sum / points.len() as f64;
    let mut area = ogeom_math::Vector::ZERO;
    for pair in points.windows(2) {
        area += (pair[0] - centre).cross(pair[1] - centre);
    }
    area
}

/// An edge's curve, carried into space by the edge's own placement.
///
/// The curve, not the trim: a segment cut short by the feature still
/// carries the line the whole edge was cut from, which is what an
/// extension runs along.
fn edge_geometry(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<Curve>> {
    let placement = edge.transform(model.datums())?;
    let Some(curve) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, .. } => Some(*curve),
            _ => None,
        })
        .and_then(|id| model.geometry().curve(id).cloned())
    else {
        return Ok(None);
    };
    Ok(Some(curve.transformed(&placement, tol)?))
}

/// An edge's vertices as this use of it runs: start first.
fn edge_ends(model: &Model, edge: &Shape) -> OgeomResult<(Shape, Shape)> {
    ogeom_algo::edge_vertices(model, edge)?
        .ok_or_else(|| ogeom_err!(Construction, "an edge of a rebuilt wire has no vertices"))
}

/// Turn a bridging chain to run the way the rim it replaces ran.
///
/// A chain that leaves and arrives at one vertex closes either way round,
/// and the walk that built it took whichever direction its first edge
/// happened to be stored in. The rim the wound took out went one way round
/// its face, and the recovered one must go the same way or the face is
/// inside out along it, which nothing in the topology notices, and the
/// mesher finds as a boundary that will not close.
fn wind_like(
    model: &Model,
    chain: &mut [Shape],
    was: &[Point],
    tol: Tolerances,
) -> OgeomResult<()> {
    let mut now = Vec::new();
    for edge in chain.iter() {
        now.extend(sample_edge(model, edge, tol)?);
    }
    if swept_area(&now).dot(swept_area(was)) >= 0.0 {
        return Ok(());
    }
    chain.reverse();
    for edge in chain.iter_mut() {
        *edge = edge.clone().reversed();
    }
    Ok(())
}

/// Walk `pool` from `from` to `to`, orienting each edge to run the way the
/// walk goes and consuming what it uses.
///
/// The gap a wound's rim leaves in a wire is bridged by the recovered
/// edges, and which of them and which way round is decided by their own
/// vertices rather than by any ordering they arrive in. Some gaps need no
/// edge at all: an end face's rim was the band's cap, and once the two
/// edges either side of it reach the corner they meet there themselves. So
/// the walk steps only when the pool offers a step, and arriving with
/// nothing taken is an answer.
fn bridge_gap(
    model: &Model,
    pool: &mut Vec<Shape>,
    from: &Shape,
    to: &Shape,
) -> OgeomResult<Vec<Shape>> {
    let mut chain = Vec::new();
    let mut here = from.node();
    loop {
        let mut found = None;
        for (index, edge) in pool.iter().enumerate() {
            let (a, b) = edge_ends(model, edge)?;
            if a.node() == here {
                found = Some((index, false, b));
                break;
            }
            if b.node() == here {
                found = Some((index, true, a));
                break;
            }
        }
        let Some((index, backwards, next)) = found else {
            if here == to.node() {
                // Nothing to bridge: the wire's own edges already meet
                // where the rim used to run.
                return Ok(chain);
            }
            if std::env::var_os("OGEOM_DEBUG_DEFEATURE").is_some() {
                let point = |v: &Shape| {
                    model
                        .node(v)
                        .and_then(|n| n.data().as_vertex())
                        .map(|d| d.point)
                };
                eprintln!(
                    "DEFEATURE bridge stuck at {:?} heading for {:?}",
                    point(from),
                    point(to)
                );
                for edge in pool.iter() {
                    let (a, b) = edge_ends(model, edge)?;
                    eprintln!("  border {:?} .. {:?}", point(&a), point(&b));
                }
            }
            ogeom_bail!(
                Construction,
                "no recovered edge bridges the wound's rim; the closure is \
                 not constructible from what the neighbours meet along"
            );
        };
        let edge = pool.remove(index);
        chain.push(if backwards { edge.reversed() } else { edge });
        here = next.node();
        if here == to.node() {
            return Ok(chain);
        }
        if pool.is_empty() {
            ogeom_bail!(
                Construction,
                "the recovered edges do not reach across the wound's rim; \
                 the closure is not constructible from them"
            );
        }
    }
}

/// The corner vertex standing on an edge's own curve, nearest the edge,
/// when one does.
fn corner_on_edge(
    model: &Model,
    edge: &Shape,
    corners: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let placement = edge.transform(model.datums())?;
    let Some((curve, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        return Ok(None);
    };
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        return Ok(None);
    };
    let geometry = geometry.transformed(&placement, tol)?;
    let head = geometry.point_at(range.0, tol)?;
    let tail = geometry.point_at(range.1, tol)?;
    let mut best: Option<(f64, Shape)> = None;
    for corner in corners {
        let Some(p) = model
            .node(corner)
            .and_then(|n| n.data().as_vertex())
            .map(|d| d.point)
        else {
            continue;
        };
        let t = parameter_near(&geometry, p, tol)?;
        if geometry.point_at(t, tol)?.distance(p) > tol.confusion() * 1e3 {
            continue;
        }
        let gap = head.distance(p).min(tail.distance(p));
        if best.as_ref().is_none_or(|(g, _)| gap < *g) {
            best = Some((gap, corner.clone()));
        }
    }
    Ok(best.map(|(_, c)| c))
}

/// The edge, extended along its own curve so its dangling end reaches the
/// corner vertex, shared across the faces that use it, so sewing rejoins
/// them on one node.
fn extend_to_corner(
    model: &mut Model,
    edge: &Shape,
    corner: &Shape,
    extended: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if let Some(found) = extended.get(&edge.node()) {
        return Ok(found.clone());
    }
    let placement = edge.transform(model.datums())?;
    let Some((curve, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        ogeom_bail!(Construction, "a dangling edge has no curve to extend");
    };
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        ogeom_bail!(Construction, "a dangling edge's curve is not in this model");
    };
    let geometry = geometry.transformed(&placement, tol)?;
    let corner_point = model
        .node(corner)
        .and_then(|n| n.data().as_vertex())
        .map(|d| d.point)
        .ok_or_else(|| ogeom_err!(Construction, "a corner vertex holds no point"))?;

    // Which end dangles: the one nearer the corner. The corner's parameter
    // on this curve comes from the geometry the curve already has.
    let head = geometry.point_at(range.0, tol)?;
    let tail = geometry.point_at(range.1, tol)?;
    let t_corner = parameter_near(&geometry, corner_point, tol)?;
    let (vertices, new_range, dangle_head) = {
        // Storage order, deliberately: the range is the stored curve's, and
        // `edge_vertices` would swap the pair for a reversed use.
        let bounds = model.children_of(edge)?;
        let (Some(va), Some(vb)) = (bounds.first().cloned(), bounds.last().cloned()) else {
            ogeom_bail!(Construction, "a dangling edge has no vertices");
        };
        if head.distance(corner_point) <= tail.distance(corner_point) {
            ((corner.clone(), vb), (t_corner, range.1), true)
        } else {
            ((va, corner.clone()), (range.0, t_corner), false)
        }
    };
    let _ = dangle_head;
    if new_range.1 <= new_range.0 {
        ogeom_bail!(
            Construction,
            "extending an edge to its corner inverted its range; the corner \
             sits on the wrong side of the edge"
        );
    }
    // A segment's stored domain ends at its own vertices; the extension is
    // the same line over a wider window.
    let geometry = match geometry {
        Curve::Line(line) => {
            let (lo, hi) = ogeom_geom::Curve3d::domain(&line);
            Curve::Line(ogeom_geom::LineCurve::over(
                line.axis(),
                lo.min(new_range.0),
                hi.max(new_range.1),
            )?)
        }
        other => other,
    };
    let built = make_edge_between(model, geometry, new_range, &vertices.0, &vertices.1, tol)?;
    extended.insert(edge.node(), built.shape.clone());
    Ok(built.shape)
}

/// The corner's parameter on a curve, by closed form where one exists and by
/// projection where not.
fn parameter_near(curve: &Curve, p: Point, tol: Tolerances) -> OgeomResult<f64> {
    match curve {
        Curve::Line(line) => {
            let axis = line.axis();
            Ok((p - axis.location).dot(axis.direction.vector()))
        }
        Curve::Circle(c) => {
            let local = c.circle().frame().to_local(p);
            Ok(local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU))
        }
        _ => Ok(ogeom_algo::project_on_curve(curve, p, 64, tol)?.parameter),
    }
}

fn nearest_distance(curve: &Curve, p: Point, tol: Tolerances) -> f64 {
    ogeom_algo::project_on_curve(curve, p, 32, tol).map_or(f64::INFINITY, |pr| pr.distance)
}

fn edge_length(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let Some((curve, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        return Ok(0.0);
    };
    let Some(geometry) = model.geometry().curve(curve) else {
        return Ok(0.0);
    };
    ogeom_algo::curve_length(geometry, range, tol)
}
