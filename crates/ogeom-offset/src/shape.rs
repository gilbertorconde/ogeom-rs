//! Offsetting a solid, and the shelling built on it.
//!
//! The topology-preserving offset: every face's surface moves along its own
//! outward normal (a plane translates, a cylinder's radius grows or
//! shrinks), and the topology is rebuilt one-for-one on the moved surfaces.
//! Vertices re-solve where their planes now meet, edges re-derive on the
//! moved supports with their directions and parameterizations preserved, and
//! band faces rebuild through [`make_revolution_band`] so seams stay seams.
//! Corners stay sharp: this is the parallel solid of the intersection join,
//! not the rounded Minkowski body.
//!
//! Shelling is the offset pointed inward and the boolean pointed at the
//! result: the cavity is the inward offset with the *removed* faces left
//! exactly where they were, so it reaches the boundary at the openings,
//! and the cut's same-domain resolution melts the flush faces away, which is
//! what opens the shell.
//!
//! The honest limits, refused by name: faces whose surfaces are not among
//! the five analytics (a spline, revolution, extrusion, trimmed or offset
//! surface has no same-family parallel to move to, though a face something
//! *replaces* (a draft's turned wall) rides through on the replacement and
//! a face moved by nothing keeps its own surface whatever the family),
//! vertices whose seats leave them under-determined, and offsets that
//! collapse the solid. Edges between moved supports re-derive exactly where
//! a line or circle exists; anywhere else the pair's own intersection is
//! marched and fitted, with its stated slop widening the edge; an edge
//! that sits unmoved on both supports rebuilds on its own curve.

use crate::wire2d::Join;
use ogeom_algo::{
    Built, History, edge_vertices, make_edge, make_edge_between, make_face_with_pcurves,
    make_revolution_band, make_vertex, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, CylinderSurface, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Cylinder, Frame, Plane, Point, Vector};
use ogeom_topo::{
    EdgeData, EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};

use ogeom_core::FastMap;

/// The displacement constraint one face puts on a point of itself.
type Displacement<'a> = dyn Fn(&Model, usize, Point) -> OgeomResult<Option<(Vector, f64)>> + 'a;

/// Canonicalize a solid whose topology is *placed*: a node placed twice (a
/// prism's far cap reusing the profile's nodes under the travel), or any
/// occurrence standing away from its node under a location (a solid moved
/// as a whole).
///
/// The rebuild below resolves everything by node, in the node's own frame,
/// which is one name for two places on an instanced solid and the wrong
/// place on a moved one. Baking restates every occurrence as its own node
/// in world coordinates, and the caller's face handles ride the bake's
/// history. A solid whose nodes each stand where they are placed passes
/// through untouched.
pub(crate) fn canonical_input(
    model: &mut Model,
    solid: &Shape,
    handles: &[Shape],
    tol: Tolerances,
) -> OgeomResult<(Shape, Vec<Shape>, Option<ogeom_algo::History>)> {
    let probe = Point::new(0.123_456_789, 9.87, -3.21);
    let mut seen: FastMap<TShapeId, Point> = FastMap::default();
    let mut instanced = false;
    'outer: for kind in [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face] {
        for occurrence in explore(model, solid, Filter::OfType(kind))? {
            let at = occurrence.transform(model.datums())?.apply(probe);
            if at.distance(probe) > tol.confusion() {
                instanced = true;
                break 'outer;
            }
            match seen.entry(occurrence.node()) {
                ogeom_core::collections::hash_map::Entry::Occupied(held) => {
                    if held.get().distance(at) > tol.confusion() {
                        instanced = true;
                        break 'outer;
                    }
                }
                ogeom_core::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(at);
                }
            }
        }
    }
    if !instanced {
        return Ok((solid.clone(), handles.to_vec(), None));
    }
    let baked = ogeom_algo::baked_shape(model, solid, tol)?;
    let mapped = handles
        .iter()
        .map(|h| match baked.history.trace(h) {
            [one] => Ok(one.clone()),
            traced => ogeom_bail!(
                Construction,
                "a face handle resolved to {} faces through the canonical \
                 rebuild; the reference is ambiguous",
                traced.len()
            ),
        })
        .collect::<OgeomResult<Vec<Shape>>>()?;
    Ok((baked.shape, mapped, Some(baked.history)))
}

/// Offset a solid by `offset`: positive grows it, negative shrinks it, and
/// the topology is preserved one-for-one.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a face,
/// edge or vertex falls outside the analytic vocabulary this rebuild speaks
/// (see the module documentation), or the offset collapses the solid.
pub fn offset_shape(
    model: &mut Model,
    solid: &Shape,
    offset: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !offset.is_finite() || offset.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "an offset of {offset} moves nothing");
    }
    let (canonical, _, prefix) = canonical_input(model, solid, &[], tol)?;
    if let Some(prefix) = prefix {
        let mut out = offset_shape(model, &canonical, offset, tol)?;
        out.history = prefix.then(&out.history);
        return Ok(out);
    }
    let built = rebuilt(model, solid, &|_| offset, &|_| None, tol)?;
    right_side_out(model, built, tol)
}

/// The offset solid, unless the offset ran its faces through each other.
///
/// Driven inward past half the solid's thickness, opposite faces cross and
/// the rebuild closes a body turned inside out there, which the boolean and
/// every measure downstream would read as material.
fn right_side_out(model: &Model, built: Built, tol: Tolerances) -> OgeomResult<Built> {
    if !ogeom_algo::inside_out_faces(model, &built.shape, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the offset runs the solid's faces through each other and \
             collapses it"
        );
    }
    Ok(built)
}

/// Offset some faces of a solid by `distance` along their outward normals
/// (negative moves them into the material), the faces around them
/// following: a plane moves parallel, a cylinder, cone, sphere or torus
/// changes radius about its own axis or centre, and each neighbour stays on
/// its own surface, its edges re-derived where the moved faces now meet it.
/// The topology is kept one for one, so the history maps every face to the
/// face it became.
///
/// # Errors
///
/// As [`offset_shape`], and additionally if `distance` moves nothing or a
/// face is not a face of `solid`. A move that would need a face the solid
/// does not have (a step where a moved face runs past a neighbour) or would
/// collapse one fails by name, as the rebuild's vertices and edges do.
pub fn offset_faces(
    model: &mut Model,
    solid: &Shape,
    faces: &[Shape],
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !distance.is_finite() || distance.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "an offset of {distance} moves nothing");
    }
    let (canonical, mapped, prefix) = canonical_input(model, solid, faces, tol)?;
    if let Some(prefix) = prefix {
        let mut out = offset_faces(model, &canonical, &mapped, distance, tol)?;
        out.history = prefix.then(&out.history);
        return Ok(out);
    }
    let chosen = chosen_faces(model, solid, faces)?;
    let built = rebuilt(
        model,
        solid,
        &|face| {
            if chosen.contains(&face.node()) {
                distance
            } else {
                0.0
            }
        },
        &|_| None,
        tol,
    )?;
    still_sound(model, solid, &chosen, built, tol)
}

/// Move some faces of a solid rigidly by `transform` (a translation or a
/// rotation), the faces around them following as in [`offset_faces`]:
/// each moved face keeps its surface, carried by the transform, and each
/// neighbour stays on its own.
///
/// # Errors
///
/// As [`offset_faces`], and if `transform` is not a rigid motion: a scale
/// or a reflection changes the faces themselves, not where they stand.
pub fn move_faces(
    model: &mut Model,
    solid: &Shape,
    faces: &[Shape],
    transform: &ogeom_math::Transform,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Transformable as _;
    use ogeom_math::TransformKind;
    if !matches!(
        transform.kind(),
        TransformKind::Identity | TransformKind::Translation | TransformKind::Rotation
    ) {
        ogeom_bail!(
            Construction,
            "moving faces takes a translation or a rotation; a scale or a \
             reflection reshapes them"
        );
    }
    if transform.kind() == TransformKind::Identity {
        ogeom_bail!(Construction, "the identity moves nothing");
    }
    let (canonical, mapped, prefix) = canonical_input(model, solid, faces, tol)?;
    if let Some(prefix) = prefix {
        let mut out = move_faces(model, &canonical, &mapped, transform, tol)?;
        out.history = prefix.then(&out.history);
        return Ok(out);
    }
    let chosen = chosen_faces(model, solid, faces)?;
    // Each surface is stated in its face's own frame: the motion is taken
    // into that frame, so the placed surface moves as the world one would.
    let mut moved: FastMap<TShapeId, SurfaceGeometry> = FastMap::default();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        if !chosen.contains(&face.node()) || moved.contains_key(&face.node()) {
            continue;
        }
        let Some(NodeData::Face(data)) = model.node(&face).map(ogeom_topo::TShape::data) else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        let local = placement.inverse()? * *transform * placement;
        moved.insert(face.node(), surface.transformed(&local, tol)?);
    }
    let built = rebuilt(
        model,
        solid,
        &|_| 0.0,
        &|face| moved.get(&face.node()).cloned(),
        tol,
    )?;
    still_sound(model, solid, &chosen, built, tol)
}

/// The edited solid, unless the edit ran a face through the rest of it.
///
/// Only the faces the edit touched can have gone wrong: the ones it moved
/// and every face sharing a vertex with them, whose boundaries followed.
/// Driven past the faces across from it, a moved face turns the solid
/// inside out there; driven into them, it crosses them.
fn still_sound(
    model: &Model,
    solid: &Shape,
    chosen: &ogeom_core::FastSet<TShapeId>,
    built: Built,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let mut corners: ogeom_core::FastSet<TShapeId> = ogeom_core::FastSet::default();
    for face in faces.iter().filter(|f| chosen.contains(&f.node())) {
        for v in explore_unique(model, face, ShapeType::Vertex)? {
            corners.insert(v.node());
        }
    }
    let mut touched: Vec<Shape> = Vec::new();
    for face in &faces {
        let near = chosen.contains(&face.node())
            || explore_unique(model, face, ShapeType::Vertex)?
                .iter()
                .any(|v| corners.contains(&v.node()));
        if near {
            touched.extend(built.history.trace(face).iter().cloned());
        }
    }
    if !ogeom_algo::inside_out_faces(model, &built.shape, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the moved faces run past the faces across from them and turn \
             the solid inside out"
        );
    }
    if !ogeom_algo::check_self_intersection_near(model, &built.shape, &touched, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the moved faces run into the rest of the solid; the edit would \
             make it cross itself"
        );
    }
    Ok(built)
}

/// The faces named, each checked to be one of the solid's.
fn chosen_faces(
    model: &Model,
    solid: &Shape,
    faces: &[Shape],
) -> OgeomResult<ogeom_core::FastSet<TShapeId>> {
    if faces.is_empty() {
        ogeom_bail!(Construction, "no face was named to move");
    }
    let own: ogeom_core::FastSet<TShapeId> =
        explore(model, solid, Filter::OfType(ShapeType::Face))?
            .iter()
            .map(Shape::node)
            .collect();
    for face in faces {
        if !own.contains(&face.node()) {
            ogeom_bail!(Construction, "a named face is not a face of the solid");
        }
    }
    Ok(faces.iter().map(Shape::node).collect())
}

/// Hollow a solid into a shell of the given wall `thickness`, opening it at
/// the `removed` faces.
///
/// Two constructions, chosen by the opening's neighbours. When a removed
/// face meets every neighbour across a corner, the cavity is the inward
/// offset of every kept face with the removed faces left in place,
/// subtracted through the boolean; the flush faces melt away, which is
/// what opens the shell. When a removed face has a *tangent* neighbour (a
/// blend melting into the face it rounds), leaving it in place would tear
/// the shared vertices, so instead the whole solid offsets inward and each
/// removed face's cavity image extrudes back out through the opening; the
/// rim a tangent opening leaves is the tapering strip a true
/// constant-thickness wall has there, which is correct rather than a
/// defect.
///
/// # Errors
///
/// As [`offset_shape`], and additionally if `thickness` is not a usable
/// length, a removed face is not a face of `solid`, or a tangent opening is
/// not planar.
pub fn make_thick_solid(
    model: &mut Model,
    solid: &Shape,
    removed: &[Shape],
    thickness: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    make_thick_solid_with(model, solid, removed, thickness, Join::Intersection, tol)
}

/// [`make_thick_solid`] with a choice of how the walls meet across an edge.
///
/// [`Join::Intersection`] extends the walls until they meet: sharp corners,
/// what [`make_thick_solid`] builds. [`Join::Arc`] rounds them about each
/// edge that is convex on the side the walls grow toward, a cylinder of the
/// wall thickness about the edge and a ball about a corner where several
/// meet: the rolling ball's parallel body. Where the growing side is
/// concave the two joins agree, and the opening faces stay flush either
/// way.
///
/// # Errors
///
/// As [`make_thick_solid`], and for [`Join::Arc`] where an opening meets a
/// tangent neighbour or the rounding fails.
pub fn make_thick_solid_with(
    model: &mut Model,
    solid: &Shape,
    removed: &[Shape],
    thickness: f64,
    join: Join,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !thickness.is_finite() || thickness.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "a wall of {thickness} holds nothing");
    }
    let (canonical, mapped, prefix) = canonical_input(model, solid, removed, tol)?;
    if let Some(prefix) = prefix {
        let mut out = make_thick_solid_with(model, &canonical, &mapped, thickness, join, tol)?;
        out.history = prefix.then(&out.history);
        return Ok(out);
    }
    // The sign is the side: positive hollows inward, negative builds the
    // walls outward around the solid, which becomes the cavity itself.
    let outward_walls = thickness < 0.0;
    let reach = thickness.abs();
    let own: Vec<TShapeId> = explore(model, solid, Filter::OfType(ShapeType::Face))?
        .iter()
        .map(Shape::node)
        .collect();
    for face in removed {
        if !own.contains(&face.node()) {
            ogeom_bail!(Construction, "a removed face is not a face of the solid");
        }
    }

    let mut tangent_opening = false;
    for face in removed {
        if has_tangent_neighbour(model, solid, face, tol)? {
            tangent_opening = true;
            break;
        }
    }
    if !tangent_opening {
        let skip: Vec<TShapeId> = removed.iter().map(Shape::node).collect();
        let moved = rebuilt(
            model,
            solid,
            &|face| {
                if skip.contains(&face.node()) {
                    0.0
                } else if outward_walls {
                    reach
                } else {
                    -reach
                }
            },
            &|_| None,
            tol,
        )?;
        let moved = right_side_out(model, moved, tol)?;
        // The arc join rounds the moved copy about every edge between two
        // moved faces that is convex on the growing side: a ball of the
        // wall's thickness touching both moved walls stands on the old edge.
        let moved = if join == Join::Arc {
            let held: Vec<TShapeId> = removed
                .iter()
                .flat_map(|f| {
                    moved
                        .history
                        .modified(f)
                        .iter()
                        .map(Shape::node)
                        .collect::<Vec<_>>()
                })
                .collect();
            let edges = growing_edges(model, &moved.shape, &held, outward_walls, tol)?;
            if edges.is_empty() {
                moved
            } else {
                ogeom_fillet::fillet_edges(model, &moved.shape, &edges, reach, tol)?
            }
        } else {
            moved
        };
        // Inward, the moved copy is the cavity carved from the solid;
        // outward, the solid is the cavity carved from the moved copy. The
        // held-in-place opening faces coincide either way, and the melt is
        // what leaves them open.
        let mut result = if outward_walls {
            ogeom_bool::cut(model, &moved.shape, solid, tol)?
        } else {
            ogeom_bool::cut(model, solid, &moved.shape, tol)?
        };
        for face in removed {
            result.history.delete(face);
        }
        return Ok(result);
    }

    if join == Join::Arc {
        ogeom_bail!(
            Construction,
            "the arc join is built where every opening meets its neighbours \
             across a corner; an opening with a tangent neighbour is not yet"
        );
    }
    // The tangent construction: everything moves together (which is what
    // keeps the tangencies intact), and each opening is drilled back out by
    // extruding its opening image through where the wall now stands.
    let displaced = if outward_walls { reach } else { -reach };
    let moved = rebuilt(model, solid, &|_| displaced, &|_| None, tol)?;
    let moved = right_side_out(model, moved, tol)?;
    let opening_normal = |model: &Model, face: &Shape| -> OgeomResult<Vector> {
        let Some(NodeData::Face(data)) = model.node(face).map(ogeom_topo::TShape::data) else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
            ogeom_bail!(
                Construction,
                "a tangent opening must be planar; a curved opening needs \
                 the general rebuild; see docs/PARITY.md, offset.shell-thicken"
            );
        };
        let mut normal = p.plane().normal().vector();
        if face.orientation() == Orientation::Reversed {
            normal = -normal;
        }
        Ok(normal)
    };
    let mut result = if outward_walls {
        // The solid itself is the cavity; the openings drill outward from
        // its own faces through the new walls.
        let mut tool = solid.clone();
        for face in removed {
            let outward = opening_normal(model, face)?;
            let punch = ogeom_algo::make_prism(model, &face.clone(), outward * (2.0 * reach), tol)?;
            tool = ogeom_bool::fuse(model, &tool, &punch.shape, tol)?.shape;
        }
        ogeom_bool::cut(model, &moved.shape, &tool, tol)?
    } else {
        let mut tool = moved.shape.clone();
        for face in removed {
            let outward = opening_normal(model, face)?;
            let [image] = moved.history.modified(face) else {
                ogeom_bail!(Construction, "a removed face has no single cavity image");
            };
            let punch =
                ogeom_algo::make_prism(model, &image.clone(), outward * (2.0 * reach), tol)?;
            tool = ogeom_bool::fuse(model, &tool, &punch.shape, tol)?.shape;
        }
        ogeom_bool::cut(model, solid, &tool, tol)?
    };
    for face in removed {
        result.history.delete(face);
    }
    Ok(result)
}

/// Whether any neighbour meets `face` tangentially along a shared edge.
fn has_tangent_neighbour(
    model: &Model,
    solid: &Shape,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;

    let own_edges: Vec<TShapeId> = explore(model, face, Filter::OfType(ShapeType::Edge))?
        .iter()
        .map(Shape::node)
        .collect();
    let normal_at = |model: &Model, face: &Shape, at: Point| -> OgeomResult<Option<Vector>> {
        let Some(NodeData::Face(data)) = model.node(face).map(ogeom_topo::TShape::data) else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let projection = ogeom_algo::project_on_surface(surface, at, 32, tol)?;
        if projection.distance > tol.confusion() * 100.0 {
            return Ok(None);
        }
        let (u, v) = projection.parameters;
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m <= tol.confusion() {
            return Ok(None);
        }
        Ok(Some(n / m))
    };
    for other in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        if other.node() == face.node() {
            continue;
        }
        for edge in explore(model, &other, Filter::OfType(ShapeType::Edge))? {
            if !own_edges.contains(&edge.node()) {
                continue;
            }
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                continue;
            };
            let mid = geometry.point_at(f64::midpoint(range.0, range.1), tol)?;
            let (Some(a), Some(b)) = (normal_at(model, face, mid)?, normal_at(model, &other, mid)?)
            else {
                continue;
            };
            if a.cross(b).magnitude() <= 1e-6 {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// A face prepared for the rebuild.
struct Prepared {
    shape: Shape,
    /// The moved surface.
    surface: SurfaceGeometry,
    /// The outward normal amount this face moved.
    amount: f64,
    /// The sign relating the face's outward side to the surface's own
    /// normal: `+1` for a Forward face.
    sign: f64,
    /// For a full revolution band (seam and two closed rings), the rings.
    rings: Option<[Shape; 2]>,
}

/// Rebuild a solid's topology on moved supports: the rebuild under both
/// entry points, every face offset by its own amount and the topology
/// re-derived on the moved surfaces.
///
/// One rule serves every element. A surface moves along its own normal (a
/// plane translates, a revolution surface's radius grows), which makes the
/// *displacement* constraint at any point of it exactly planar: normal
/// there, offset amount along it. Vertices solve those constraints in the
/// least-squares sense and then Newton-polish onto the moved surfaces
/// themselves, edges re-derive from the moved pair (a line from its planes'
/// constraints, a circle from the pair's analytic intersection re-framed on
/// its old axes so parameters and orientations carry), and faces rebuild
/// wire by wire with exact pcurves, or wholesale through
/// [`make_revolution_band`] where a seam says the face wraps.
///
/// `amount_of` says how far each face travels along its own outward normal;
/// `instead_of` may hand back a surface to use *in place* of that move,
/// which is how an operation that turns a face rather than translating it
/// (a draft) rides the same rebuild. The two are exclusive per face: a
/// surface supplied by `instead_of` is taken as it stands.
pub(crate) fn rebuilt(
    model: &mut Model,
    solid: &Shape,
    amount_of: &dyn Fn(&Shape) -> f64,
    instead_of: &dyn Fn(&Shape) -> Option<SurfaceGeometry>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Surface as _;
    let faces = explore(model, solid, Filter::OfType(ShapeType::Face))?;

    // Move every surface.
    let span = ogeom_algo::shape_bounds(model, solid, tol)?.diagonal();
    let mut prepared: Vec<Prepared> = Vec::with_capacity(faces.len());
    for face in &faces {
        let amount = amount_of(face);
        let Some(node) = model.node(face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        let edges = explore(model, face, Filter::OfType(ShapeType::Edge))?;
        let mut counts: FastMap<TShapeId, usize> = FastMap::default();
        for e in &edges {
            *counts.entry(e.node()).or_insert(0) += 1;
        }
        let has_seam = counts.values().any(|c| *c >= 2);
        let closed_rings: Vec<Shape> = edges
            .iter()
            .filter(|e| {
                edge_vertices(model, e)
                    .ok()
                    .flatten()
                    .is_some_and(|(a, b)| a.node() == b.node())
            })
            .cloned()
            .collect();

        let both_poles = closed_rings.len() == 2
            && closed_rings.iter().all(|ring| {
                model
                    .node(ring)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
            });
        let grow = amount.abs() * 4.0 + 1.0;
        let replacement = instead_of(face);
        let moved: SurfaceGeometry = if let Some(given) = replacement {
            given
        } else if amount == 0.0 {
            // A face a draft or a partial offset leaves alone stays on its
            // own surface, whatever family that is: moving by nothing is
            // identity, not a construction the family has to support. Its
            // window opens by the solid's own size, so a neighbour moved
            // beyond the face's edge still meets it; the face's trim, not
            // the window, says what is kept.
            match surface {
                SurfaceGeometry::Plane(p) => {
                    let ((u0, u1), (v0, v1)) = surface.domain();
                    PlaneSurface::over(p.plane(), (u0 - span, u1 + span), (v0 - span, v1 + span))?
                        .into()
                }
                SurfaceGeometry::Cylinder(c) => {
                    let (_, (v0, v1)) = surface.domain();
                    CylinderSurface::new(c.cylinder(), (v0 - span, v1 + span))?.into()
                }
                SurfaceGeometry::Cone(co) => {
                    // Up to the apex and no farther: past it is the other
                    // nappe.
                    let (_, (v0, v1)) = surface.domain();
                    let apex = co.apex_height();
                    let (lo, hi) = if apex <= v0 {
                        ((v0 - span).max(apex), v1 + span)
                    } else {
                        (v0 - span, (v1 + span).min(apex))
                    };
                    ogeom_geom::ConeSurface::new(co.cone(), (lo, hi))?.into()
                }
                _ => surface.clone(),
            }
        } else {
            match surface {
                SurfaceGeometry::Plane(p) => {
                    let plane = p.plane();
                    let ((u0, u1), (v0, v1)) = surface.domain();
                    let shifted = Plane::new(Frame::new(
                        plane.origin() + plane.normal().vector() * (sign * amount),
                        plane.normal(),
                        plane.frame().x(),
                        tol,
                    )?);
                    PlaneSurface::over(shifted, (u0 - grow, u1 + grow), (v0 - grow, v1 + grow))?
                        .into()
                }
                SurfaceGeometry::Cylinder(c) => {
                    let cylinder = c.cylinder();
                    let grown = sign.mul_add(amount, cylinder.radius());
                    if grown <= tol.confusion() {
                        ogeom_bail!(Construction, "the offset consumes the cylinder's radius");
                    }
                    let (_, (v0, v1)) = surface.domain();
                    CylinderSurface::new(
                        Cylinder::new(cylinder.frame(), grown, tol)?,
                        (v0 - grow, v1 + grow),
                    )?
                    .into()
                }
                SurfaceGeometry::Sphere(sp) => {
                    let sphere = sp.sphere();
                    let grown = sign.mul_add(amount, sphere.radius());
                    if grown <= tol.confusion() {
                        ogeom_bail!(Construction, "the offset consumes the sphere's radius");
                    }
                    // Concentric in the same frame: the chart carries over,
                    // so the face's seam and poles stay where its trim had
                    // them rather than landing on its boundary.
                    ogeom_geom::SphereSurface::new(ogeom_math::Sphere::new(
                        sphere.frame(),
                        grown,
                        tol,
                    )?)
                    .into()
                }
                SurfaceGeometry::Torus(t) => {
                    let torus = t.torus();
                    let grown = sign.mul_add(amount, torus.minor_radius());
                    if grown <= tol.confusion() {
                        ogeom_bail!(Construction, "the offset consumes the torus's tube");
                    }
                    ogeom_geom::TorusSurface::new(ogeom_math::Torus::new(
                        torus.frame(),
                        torus.major_radius(),
                        grown,
                        tol,
                    )?)
                    .into()
                }
                SurfaceGeometry::Cone(co) => {
                    let cone = co.cone();
                    // The parallel cone: same axis and half-angle, the reference
                    // radius moved by the offset over the slant's cosine.
                    let grown = (sign * amount / cone.half_angle().cos())
                        .mul_add(1.0, cone.reference_radius());
                    if grown <= tol.confusion() {
                        ogeom_bail!(Construction, "the offset consumes the cone's throat");
                    }
                    let (_, (v0, v1)) = surface.domain();
                    ogeom_geom::ConeSurface::new(
                        ogeom_math::Cone::new(cone.frame(), grown, cone.half_angle(), tol)?,
                        (v0 - grow, v1 + grow),
                    )?
                    .into()
                }
                _ => ogeom_bail!(
                    Construction,
                    "offsetting a face on this surface needs a construction \
                     the rebuild does not yet speak; see docs/PARITY.md, offset.shell-thicken"
                ),
            }
        };
        // A band rebuilds wholesale only on a surface of revolution; a
        // drafted wall on a fitted support is a band the wire path
        // assembles, seam and all.
        let fitted_support = matches!(moved, SurfaceGeometry::BSpline(_));
        prepared.push(Prepared {
            shape: face.clone(),
            surface: moved,
            amount,
            sign,
            // A band needs a ring with an angle to anchor its chart; a whole
            // sphere, bounded by its two poles alone, is assembled wire by
            // wire like any other seamed face.
            rings: if has_seam && closed_rings.len() == 2 && !fitted_support && !both_poles {
                Some([closed_rings[0].clone(), closed_rings[1].clone()])
            } else {
                None
            },
        });
    }

    // Which faces meet each edge, seams excluded by their double use.
    let mut edge_faces: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    for (fi, face) in faces.iter().enumerate() {
        for e in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            let entry = edge_faces.entry(e.node()).or_default();
            if !entry.contains(&fi) {
                entry.push(fi);
            }
        }
    }

    // The displacement constraint each face puts on a point of itself: the
    // surface normal there, moved its amount along it. Exact, because a
    // normal offset moves every point of a surface along its own normal.
    let constraint = |model: &Model, fi: usize, at: Point| -> OgeomResult<Option<(Vector, f64)>> {
        let face = &faces[fi];
        let Some(node) = model.node(face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let projection = ogeom_algo::project_on_surface(surface, at, 32, tol)?;
        if projection.distance > tol.confusion() * 100.0 {
            return Ok(None);
        }
        let (u, v) = projection.parameters;
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m <= tol.confusion() {
            return Ok(None);
        }
        let outward = n / m * prepared[fi].sign;
        Ok(Some((outward, prepared[fi].amount)))
    };

    // New vertices: the linear constraint solve seeds a Newton polish onto
    // the moved surfaces themselves; the tangent-plane answer is exact for
    // planes and off by the surfaces' own curvature otherwise.
    let mut new_vertices: FastMap<TShapeId, (Shape, Point)> = FastMap::default();
    for vertex in explore_unique(model, solid, ShapeType::Vertex)? {
        let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            continue;
        };
        let at = vertex.transform(model.datums())?.apply(data.point);
        let mut seats: Vec<usize> = Vec::new();
        for (fi, face) in faces.iter().enumerate() {
            for v in explore(model, face, Filter::OfType(ShapeType::Vertex))? {
                if v.node() == vertex.node() && !seats.contains(&fi) {
                    seats.push(fi);
                }
            }
        }
        if seats.is_empty() {
            continue;
        }
        // Independent constraints only: tangent faces share their normal and
        // must agree on the displacement, or the vertex tears.
        let mut normals: Vec<Vector> = Vec::new();
        let mut amounts: Vec<f64> = Vec::new();
        let mut kept: Vec<usize> = Vec::new();
        for fi in &seats {
            let Some((n, w)) = constraint(model, *fi, at)? else {
                continue;
            };
            if let Some(k) = normals
                .iter()
                .position(|m| m.cross(n).magnitude() <= tol.angular().max(1e-6))
            {
                if (amounts[k] - w).abs() > tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "two tangent faces move a shared vertex by different \
                         amounts; the offset tears it"
                    );
                }
                continue;
            }
            normals.push(n);
            amounts.push(w);
            kept.push(*fi);
        }
        if normals.is_empty() {
            // A cone's apex has no normal to offer (the projection there is
            // degenerate), but the parallel cone knows exactly where its own
            // apex went.
            let mut apex: Option<Point> = None;
            for fi in &seats {
                let Some(node) = model.node(&faces[*fi]) else {
                    continue;
                };
                let NodeData::Face(data) = node.data() else {
                    continue;
                };
                let Some(SurfaceGeometry::Cone(old)) = model.geometry().surface(data.surface)
                else {
                    continue;
                };
                if old.cone().apex().distance(at) > tol.confusion() * 100.0 {
                    continue;
                }
                if let SurfaceGeometry::Cone(moved_cone) = &prepared[*fi].surface {
                    apex = Some(moved_cone.cone().apex());
                    break;
                }
            }
            // A sphere's pole has no normal to offer either, and a moved
            // sphere keeps its chart (concentric under an offset, carried
            // along under a rigid move): the pole is where the moved surface
            // stands at the old one's parameters there.
            if apex.is_none() {
                for fi in &seats {
                    let Some(NodeData::Face(data)) = model.node(&faces[*fi]).map(|n| n.data())
                    else {
                        continue;
                    };
                    let Some(old @ SurfaceGeometry::Sphere(_)) =
                        model.geometry().surface(data.surface)
                    else {
                        continue;
                    };
                    if !matches!(prepared[*fi].surface, SurfaceGeometry::Sphere(_)) {
                        continue;
                    }
                    let local = faces[*fi].transform(model.datums())?;
                    let found =
                        ogeom_algo::project_on_surface(old, local.inverse()?.apply(at), 32, tol)?;
                    if found.distance > tol.confusion() * 100.0 {
                        continue;
                    }
                    let (u, v) = found.parameters;
                    apex = Some(local.apply(prepared[*fi].surface.point_at(u, v, tol)?));
                    break;
                }
            }
            let Some(moved) = apex else {
                ogeom_bail!(
                    Construction,
                    "a vertex with no seat the rebuild can read cannot be \
                     re-solved"
                );
            };
            new_vertices.insert(vertex.node(), (make_vertex(model, moved).shape, moved));
            continue;
        }
        if normals.len() == 1 {
            // Every seat is tangent to the rest, and the dedup above made
            // them agree on the amount. A normal offset moves each point of a
            // surface along its own normal, so the shared normal is the exact
            // answer: no corner to solve, nothing to polish.
            let moved = at + normals[0] * amounts[0];
            corner_met(&prepared, &kept, moved, tol)?;
            new_vertices.insert(vertex.node(), (make_vertex(model, moved).shape, moved));
            continue;
        }
        // A vertex where a seam ends has two seats, not three, and the
        // third constraint is the seam itself: the vertex is where the moved
        // support's seam column meets the other seat. Solved as that
        // crossing where the seam has an iso-curve to offer; the nearest
        // point two seats agree on is somewhere along their whole edge.
        if kept.len() == 2
            && let Some(moved) = seam_end(model, &faces, &prepared, &vertex, &kept, at, tol)?
        {
            new_vertices.insert(vertex.node(), (make_vertex(model, moved).shape, moved));
            continue;
        }
        let mut moved = at + solve_corner(&normals, &amounts, tol)?;
        // Newton onto the moved surfaces: residuals are the signed
        // distances, gradients the normals, and the same least-squares
        // machinery takes the step.
        for _ in 0..8 {
            let mut ns: Vec<Vector> = Vec::new();
            let mut rs: Vec<f64> = Vec::new();
            for fi in &kept {
                let projection =
                    ogeom_algo::project_on_surface(&prepared[*fi].surface, moved, 32, tol)?;
                let (u, v) = projection.parameters;
                let (du, dv) = prepared[*fi].surface.d1_at(u, v, tol)?;
                let n = du.cross(dv);
                let m = n.magnitude();
                if m <= tol.confusion() {
                    continue;
                }
                let n = n / m;
                let foot = prepared[*fi].surface.point_at(u, v, tol)?;
                ns.push(n);
                rs.push((moved - foot).dot(n));
            }
            if ns.len() < 2 {
                break;
            }
            let worst = rs.iter().fold(0.0_f64, |a, r| a.max(r.abs()));
            if worst <= tol.confusion() * 0.1 {
                break;
            }
            let step: Vec<f64> = rs.iter().map(|r| -r).collect();
            moved += solve_corner(&ns, &step, tol)?;
        }
        corner_met(&prepared, &kept, moved, tol)?;
        new_vertices.insert(vertex.node(), (make_vertex(model, moved).shape, moved));
    }

    // How many times each edge occurs across all faces; a seam is one face
    // using an edge twice, which face-deduplicated sides cannot see.
    let mut edge_uses: FastMap<TShapeId, usize> = FastMap::default();
    for face in &faces {
        for e in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            *edge_uses.entry(e.node()).or_insert(0) += 1;
        }
    }

    // New edges on the moved supports.
    let mut new_edges: FastMap<TShapeId, Shape> = FastMap::default();
    let mut history = History::new();
    for edge in explore_unique(model, solid, ShapeType::Edge)? {
        let sides = edge_faces.get(&edge.node()).cloned().unwrap_or_default();
        if sides.len() != 2 {
            if edge_uses.get(&edge.node()).copied().unwrap_or(0) >= 2 {
                // A seam. A band face rebuilds its own; a face assembled wire
                // by wire (a band a boolean split into arc rings) needs the
                // moved seam here: the same iso-column on the moved surface,
                // which chart preservation makes exact.
                if let [fi] = sides.as_slice()
                    && let Some(built) =
                        rebuilt_seam_edge(model, &edge, &prepared[*fi], &new_vertices, tol)?
                {
                    history.modify(&edge, built.clone());
                    new_edges.insert(edge.node(), built);
                }
                continue;
            }
            // A genuinely single-sided edge: the ring a boolean left
            // coincident with a neighbour's twin, or a cone's apex.
            let Some(built) =
                rebuilt_lone_edge(model, &edge, &sides, &constraint, &new_vertices, tol)?
            else {
                ogeom_bail!(
                    Construction,
                    "an edge with one face is neither a ring nor an apex; \
                     the offset cannot re-derive it"
                );
            };
            history.modify(&edge, built.clone());
            new_edges.insert(edge.node(), built);
            continue;
        }
        let (curve, range) = {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
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
        let forward = if edge.orientation() == Orientation::Reversed {
            edge.reversed()
        } else {
            edge.clone()
        };
        let built = match &curve {
            Curve::Line(_) => {
                // A straight edge is the line through its own re-solved
                // ends. That is true whether the supports were translated
                // or turned (an offset leaves the direction alone and this
                // reproduces it, a draft does not and this follows it),
                // whereas a line anchored where the old one sat misses its
                // own vertices the moment either end moves sideways.
                let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                    ogeom_bail!(Construction, "a straight edge has no vertices");
                };
                let (Some((v_from, p_from)), Some((v_to, p_to))) = (
                    new_vertices.get(&sv.node()).cloned(),
                    new_vertices.get(&ev.node()).cloned(),
                ) else {
                    ogeom_bail!(Construction, "an edge end has no re-solved vertex");
                };
                if p_to.distance(p_from) <= tol.parametric() {
                    ogeom_bail!(Construction, "the offset collapses an edge");
                }
                let segment = LineCurve::segment(p_from, p_to, tol)?;
                let (t0, t1) = segment.domain();
                let moved: Curve = segment.into();
                make_edge_between(model, moved, (t0, t1), &v_from, &v_to, tol)?.shape
            }
            Curve::Circle(c)
                if !matches!(prepared[sides[0]].surface, SurfaceGeometry::BSpline(_))
                    && !matches!(prepared[sides[1]].surface, SurfaceGeometry::BSpline(_)) =>
            {
                // The moved pair's own analytic intersection, taken in the
                // circle's old frame so parameters and orientations carry.
                // Between analytic supports a circle stays a circle; against
                // a fitted support it is whatever the march finds, below.
                let circle = c.circle();
                let found = ogeom_intersect::intersect_surfaces(
                    &prepared[sides[0]].surface,
                    &prepared[sides[1]].surface,
                    ogeom_intersect::IntersectOptions::default(),
                    tol,
                )?;
                let ogeom_intersect::SurfaceIntersection::Along(candidates) = found else {
                    ogeom_bail!(
                        Construction,
                        "the moved faces no longer meet along the edge they \
                         shared; the offset collapses it"
                    );
                };
                let mut best: Option<(ogeom_math::Circle, f64)> = None;
                for section in &candidates {
                    let Curve::Circle(cc) = &section.curve else {
                        continue;
                    };
                    let candidate = cc.circle();
                    let score = candidate.centre().distance(circle.centre())
                        + (candidate.radius() - circle.radius()).abs();
                    if best.as_ref().is_none_or(|(_, held)| score < *held) {
                        best = Some((candidate, score));
                    }
                }
                let Some((candidate, _)) = best else {
                    ogeom_bail!(
                        Construction,
                        "the moved faces meet along nothing circular where a \
                         circle was; the offset needs the general rebuild"
                    );
                };
                let reframed = ogeom_math::Circle::new(
                    Frame::new(
                        candidate.centre(),
                        circle.frame().z(),
                        circle.frame().x(),
                        tol,
                    )?,
                    candidate.radius(),
                    tol,
                )?;
                let moved: Curve = ogeom_geom::CircleCurve::new(reframed).into();
                let closed = {
                    let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                        ogeom_bail!(Construction, "a ring has no vertex");
                    };
                    sv.node() == ev.node()
                };
                if closed {
                    make_edge(model, moved, range, tol)?.shape
                } else {
                    let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                        ogeom_bail!(Construction, "an arc has no vertices");
                    };
                    let (Some((v_from, p_from)), Some((v_to, p_to))) = (
                        new_vertices.get(&sv.node()).cloned(),
                        new_vertices.get(&ev.node()).cloned(),
                    ) else {
                        ogeom_bail!(Construction, "an arc end has no re-solved vertex");
                    };
                    let angle_of = |p: Point| {
                        let l = reframed.frame().to_local(p);
                        l.y.atan2(l.x)
                    };
                    let tau = core::f64::consts::TAU;
                    let mut t0 = angle_of(p_from);
                    let mut t1 = angle_of(p_to);
                    // Keep the new range in the old one's winding and span.
                    while t0 < range.0 - core::f64::consts::PI {
                        t0 += tau;
                    }
                    while t0 > range.0 + core::f64::consts::PI {
                        t0 -= tau;
                    }
                    while t1 <= t0 + tol.parametric() {
                        t1 += tau;
                    }
                    if (t1 - t0) - (range.1 - range.0) > core::f64::consts::PI {
                        t1 -= tau;
                    }
                    if t1 <= t0 + tol.parametric() {
                        ogeom_bail!(Construction, "the offset collapses an arc");
                    }
                    make_edge_between(model, moved, (t0, t1), &v_from, &v_to, tol)?.shape
                }
            }
            _ => {
                // The general edge. First the still question: a hinge edge
                // (a draft's neutral crossing) sits on both moved supports
                // exactly where it always was, and an edge that did not move
                // rebuilds on its own curve rather than on a march of it.
                let unmoved = {
                    let mut worst = 0.0_f64;
                    'probe: for i in 0..9 {
                        #[allow(clippy::cast_precision_loss)]
                        let t = range.0 + (range.1 - range.0) * (i as f64) / 8.0;
                        let p = curve.point_at(t, tol)?;
                        for side in [sides[0], sides[1]] {
                            let Ok(near) =
                                ogeom_algo::project_on_surface(&prepared[side].surface, p, 17, tol)
                            else {
                                worst = f64::INFINITY;
                                break 'probe;
                            };
                            worst = worst.max(near.distance);
                        }
                    }
                    // Within the moved supports' own stated accuracy: a
                    // fitted support holds its points only to the fit
                    // target, and the hinge is exactly on it by less.
                    (worst <= (tol.confusion() * 1e3).max(1e-4)).then_some(worst)
                };
                if let Some(worst) = unmoved {
                    let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                        ogeom_bail!(Construction, "an edge has no vertices");
                    };
                    let closed = sv.node() == ev.node();
                    let built = if closed {
                        // On the vertex the rest of the rebuild uses (a
                        // seam starts from it), not one of the curve's own.
                        match new_vertices.get(&sv.node()).cloned() {
                            Some((v_at, p_at)) => {
                                let gap = curve.point_at(range.0, tol)?.distance(p_at);
                                if gap > tol.confusion() {
                                    model.widen(&v_at, ogeom_core::Tolerance::new(gap * 2.0)?)?;
                                }
                                make_edge_between(model, curve.clone(), range, &v_at, &v_at, tol)?
                                    .shape
                            }
                            None => make_edge(model, curve.clone(), range, tol)?.shape,
                        }
                    } else {
                        let (Some((v_from, p_from)), Some((v_to, p_to))) = (
                            new_vertices.get(&sv.node()).cloned(),
                            new_vertices.get(&ev.node()).cloned(),
                        ) else {
                            ogeom_bail!(Construction, "an edge end has no re-solved vertex");
                        };
                        // The ends re-solved against a fitted support land a
                        // fit's breadth from the curve that did not move; the
                        // vertices own that breadth.
                        let gap = curve
                            .point_at(range.0, tol)?
                            .distance(p_from)
                            .min(curve.point_at(range.0, tol)?.distance(p_to))
                            .max(
                                curve
                                    .point_at(range.1, tol)?
                                    .distance(p_to)
                                    .min(curve.point_at(range.1, tol)?.distance(p_from)),
                            );
                        if gap > tol.confusion() {
                            for v in [&v_from, &v_to] {
                                model.widen(v, ogeom_core::Tolerance::new(gap * 2.0)?)?;
                            }
                        }
                        make_edge_between(model, curve.clone(), range, &v_from, &v_to, tol)?.shape
                    };
                    if worst > tol.confusion() {
                        model.widen(&built, ogeom_core::Tolerance::new(worst)?)?;
                    }
                    history.modify(&edge, built.clone());
                    new_edges.insert(edge.node(), built);
                    continue;
                }
                // Otherwise the moved pair's own intersection, marched where
                // no closed form exists (a drafted spline wall re-meeting
                // its cap plane), with the candidate nearest the old edge
                // kept and trimmed between the re-solved ends. The section's
                // stated slop widens the edge; nothing pretends the fit is
                // exact.
                let mid = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
                let found = ogeom_intersect::intersect_surfaces(
                    &prepared[sides[0]].surface,
                    &prepared[sides[1]].surface,
                    ogeom_intersect::IntersectOptions::default(),
                    tol,
                )?;
                let ogeom_intersect::SurfaceIntersection::Along(candidates) = found else {
                    ogeom_bail!(
                        Construction,
                        "the moved faces no longer meet along the edge they \
                         shared; the offset collapses it"
                    );
                };
                let mut best: Option<(Curve, f64, f64)> = None;
                for section in candidates {
                    let Ok(projected) = ogeom_algo::project_on_curve(&section.curve, mid, 64, tol)
                    else {
                        continue;
                    };
                    if best
                        .as_ref()
                        .is_none_or(|(_, _, held)| projected.distance < *held)
                    {
                        best = Some((section.curve, section.tolerance, projected.distance));
                    }
                }
                let Some((moved, slop, _)) = best else {
                    ogeom_bail!(
                        Construction,
                        "the moved faces meet along nothing where the edge \
                         was; the offset collapses it"
                    );
                };
                let closed = {
                    let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                        ogeom_bail!(Construction, "an edge has no vertices");
                    };
                    sv.node() == ev.node()
                };
                let built = if closed {
                    // A ring's one vertex is a corner the neighbours' seams
                    // start from, re-solved like any other: the marched
                    // section is re-seamed to begin there, so the ring and
                    // the seam meet at one vertex rather than at two a
                    // section's start apart.
                    let Some((sv, _)) = edge_vertices(model, &forward)? else {
                        ogeom_bail!(Construction, "a ring has no vertex");
                    };
                    match (new_vertices.get(&sv.node()).cloned(), &moved) {
                        (Some((v_at, p_at)), Curve::BSpline(spline)) => {
                            let t = ogeom_algo::project_on_curve(&moved, p_at, 64, tol)?;
                            let (lo, hi) = moved.domain();
                            let seamed: Curve = if t.parameter > lo + tol.parametric()
                                && t.parameter < hi - tol.parametric()
                            {
                                Curve::BSpline(spline.reseamed_at(t.parameter, tol)?)
                            } else {
                                moved.clone()
                            };
                            // Run the way the old ring ran: the wire uses the
                            // rebuilt edge with the old orientation, and a
                            // march has no opinion about direction.
                            let seamed = {
                                use ogeom_geom::Reversible as _;
                                let (a, _) = seamed.domain();
                                let old = curve.d1_at(range.0, tol)?;
                                if seamed.d1_at(a, tol)?.dot(old) < 0.0 {
                                    seamed.reversed()
                                } else {
                                    seamed
                                }
                            };
                            let miss = t.distance.max(slop);
                            if miss > tol.confusion() {
                                model.widen(&v_at, ogeom_core::Tolerance::new(miss * 2.0)?)?;
                            }
                            let window = seamed.domain();
                            make_edge_between(model, seamed, window, &v_at, &v_at, tol)?.shape
                        }
                        _ => {
                            let window = moved.domain();
                            make_edge(model, moved, window, tol)?.shape
                        }
                    }
                } else {
                    let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                        ogeom_bail!(Construction, "an edge has no vertices");
                    };
                    let (Some((v_from, p_from)), Some((v_to, p_to))) = (
                        new_vertices.get(&sv.node()).cloned(),
                        new_vertices.get(&ev.node()).cloned(),
                    ) else {
                        ogeom_bail!(Construction, "an edge end has no re-solved vertex");
                    };
                    // The fitted section lands within its stated slop of the
                    // re-solved ends; the vertices own that slop.
                    if slop > tol.confusion() {
                        for v in [&v_from, &v_to] {
                            model.widen(v, ogeom_core::Tolerance::new(slop * 2.0)?)?;
                        }
                    }
                    let ta = ogeom_algo::project_on_curve(&moved, p_from, 64, tol)?.parameter;
                    let tb = ogeom_algo::project_on_curve(&moved, p_to, 64, tol)?.parameter;
                    if (tb - ta).abs() <= tol.parametric() {
                        ogeom_bail!(Construction, "the offset collapses an edge");
                    }
                    // The ends run with the curve or against it; a run
                    // against builds on the reversed parameterization so
                    // the edge still leaves `v_from` first. On a periodic
                    // section the order of the two parameters says nothing
                    // (either arc joins them), so the old edge's direction
                    // picks the run and the far end is taken one period on
                    // where it falls behind the near one.
                    let (moved, ta, tb) = if moved.is_periodic() {
                        let along = moved.d1_at(ta, tol)?.dot(curve.d1_at(range.0, tol)?) >= 0.0;
                        let (moved, ta, tb) = if along {
                            (moved, ta, tb)
                        } else {
                            use ogeom_geom::Reversible as _;
                            let (lo, hi) = moved.domain();
                            (moved.reversed(), lo + hi - ta, lo + hi - tb)
                        };
                        let (lo, hi) = moved.domain();
                        let tb = ta + (tb - ta).rem_euclid(hi - lo);
                        if tb - ta <= tol.parametric() {
                            ogeom_bail!(Construction, "the offset collapses an edge");
                        }
                        (moved, ta, tb)
                    } else if ta <= tb {
                        (moved, ta, tb)
                    } else {
                        use ogeom_geom::Reversible as _;
                        let (lo, hi) = moved.domain();
                        (moved.reversed(), lo + hi - ta, lo + hi - tb)
                    };
                    make_edge_between(model, moved, (ta, tb), &v_from, &v_to, tol)?.shape
                };
                if slop > tol.confusion() {
                    model.widen(&built, ogeom_core::Tolerance::new(slop)?)?;
                }
                built
            }
        };
        history.modify(&edge, built.clone());
        new_edges.insert(edge.node(), built);
    }

    // Faces: bands wholesale, everything else wire by wire with exact
    // pcurves on the moved surface.
    let mut rebuilt_faces: Vec<Shape> = Vec::with_capacity(prepared.len());
    for prep in &prepared {
        let built = if let Some(rings) = &prep.rings {
            let (Some(lo), Some(hi)) = (
                new_edges.get(&rings[0].node()),
                new_edges.get(&rings[1].node()),
            ) else {
                ogeom_bail!(Construction, "a band's ring was not rebuilt");
            };
            let band = make_revolution_band(model, &prep.surface, lo, hi, tol)?;
            if prep.shape.orientation() == Orientation::Reversed {
                band.reversed()
            } else {
                band
            }
        } else {
            let mut wires: Vec<Vec<Shape>> = Vec::new();
            let mut face_uses: FastMap<TShapeId, usize> = FastMap::default();
            // The wires as the face stores them: the rebuilt face takes the
            // old one's sense below.
            let stored = prep.shape.oriented(Orientation::Forward);
            for wire in explore(model, &stored, Filter::OfType(ShapeType::Wire))? {
                let mut edges: Vec<Shape> = Vec::new();
                // The wire's own order, not the walker's: a rebuilt wire is
                // re-chained edge to edge, and the walk order is not a chain.
                for used in model.ordered_children_of(&wire)? {
                    *face_uses.entry(used.node()).or_insert(0) += 1;
                    let Some(fresh) = new_edges.get(&used.node()) else {
                        ogeom_bail!(Construction, "a face edge was not rebuilt");
                    };
                    edges.push(if used.orientation() == Orientation::Reversed {
                        fresh.reversed()
                    } else {
                        fresh.clone()
                    });
                }
                wires.push(edges);
            }
            let face = if face_uses.values().any(|c| *c >= 2) {
                // A seam in a wire-assembled face: a band a boolean split
                // into arc rings. Every ordinary pcurve is recomputed on the
                // moved surface; the seam's columns carry over, which the
                // seam rebuild already validated against the re-solved ends.
                assembled_with_seam(model, prep, &wires, &new_edges, tol)?
            } else {
                make_face_with_pcurves(model, prep.surface.clone(), &wires, tol)?.shape
            };
            if prep.shape.orientation() == Orientation::Reversed {
                face.reversed()
            } else {
                face
            }
        };
        history.modify(&prep.shape, built.clone());
        rebuilt_faces.push(built);
    }

    let sewn = sew(model, &rebuilt_faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the offset solid did not close");
    }
    // The faces carried their use-orientations through; the *shell* has one
    // too, and a solid whose outer shell was used reversed reads inside out
    // if the rebuilt shell forgets it.
    let outer = {
        let old_reversed = model
            .children_of(solid)?
            .first()
            .is_some_and(|s| s.orientation() == Orientation::Reversed);
        if old_reversed {
            sewn.shells[0].reversed()
        } else {
            sewn.shells[0].clone()
        }
    };
    // Put together raw: `make_solid` would turn an inside-out shell to
    // face out, and the guard below reads it.
    let offset = model.add_solid(std::slice::from_ref(&outer))?;

    // The one global guard the local checks cannot give: an offset that
    // moved faces past each other builds a shell that is closed and inside
    // out. Its measured volume is the tell.
    // The guard meshes at the default deflection, and a thin tangential
    // cusp (a small blend meeting its face) can defeat that resolution
    // without anything being wrong. One finer retry separates a mesh that
    // cannot see the cusp from a solid that is genuinely inside out.
    // The finer retry is a fraction of the part's own size, not a fixed
    // length: a metre-sized part meshed at a tenth of a micron is millions
    // of triangles for a sign.
    let size = ogeom_algo::shape_bounds(model, &offset, tol)?.diagonal();
    let fine = (size * 1e-5).clamp(
        tol.confusion() * 1e2,
        ogeom_mesh::Deflection::default().chord,
    );
    let mut mass = None;
    let mut first_error = None;
    for chord in [ogeom_mesh::Deflection::default().chord, fine] {
        let deflection = ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        };
        match ogeom_algo::volume_properties(model, &offset, deflection, tol) {
            Ok(props) => {
                mass = Some(props.mass);
                break;
            }
            Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
                return Err(e);
            }
            Err(e) => {
                first_error.get_or_insert(e);
            }
        }
    }
    let Some(mass) = mass else {
        ogeom_bail!(
            Construction,
            "the offset solid's mesh does not close at any tried resolution{}",
            first_error.map_or_else(String::new, |e| format!(": {e}"))
        );
    };
    if !mass.is_finite() || mass <= tol.confusion() {
        ogeom_bail!(Construction, "the offset collapses the solid");
    }

    history.modify(solid, offset.clone());
    Ok(Built::new(offset, history))
}

/// Rebuild a seam for a face assembled wire by wire: the same iso-column on
/// the moved surface, over the same rows.
///
/// A same-family move preserves the chart (every point travels along its
/// own normal without changing its parameters), so the moved seam sits at
/// the column the old one's own pcurves state, between the re-solved end
/// vertices. `None` when the old edge carries no seam representation on this
/// face's surface.
fn rebuilt_seam_edge(
    model: &mut Model,
    edge: &Shape,
    prep: &Prepared,
    new_vertices: &FastMap<TShapeId, (Shape, Point)>,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    use ogeom_geom::Curve2d as _;

    let old_surface = {
        let Some(NodeData::Face(data)) = model.node(&prep.shape).map(ogeom_topo::TShape::data)
        else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        data.surface
    };
    let found = {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let mut found = None;
        for repr in &data.representations {
            if let EdgeRepr::Seam {
                forward,
                surface,
                range,
                ..
            } = repr
                && *surface == old_surface
            {
                let Some(pcurve) = model.geometry().pcurve(*forward) else {
                    ogeom_bail!(Dangling, "a seam pcurve is not in this model");
                };
                // The pcurve states the column the seam sits at. The rows
                // cannot come from it: the seam's ends are corner vertices,
                // moved by the corner solve rather than by this face alone.
                found = Some(pcurve.point_at(range.0, tol)?.x);
                break;
            }
        }
        found
    };
    let Some(column) = found else {
        return Ok(None);
    };
    // Rebuilt edges are kept in their node's own direction, each wire
    // turning its use of them as it did the old one's, so the ends are read
    // off the node's forward occurrence, however the seam was reached.
    let forward = if edge.orientation() == Orientation::Reversed {
        edge.reversed()
    } else {
        edge.clone()
    };
    let Some((sv, ev)) = edge_vertices(model, &forward)? else {
        ogeom_bail!(Construction, "a seam has no vertices");
    };
    let (Some((v_from, p_from)), Some((v_to, p_to))) = (
        new_vertices.get(&sv.node()).cloned(),
        new_vertices.get(&ev.node()).cloned(),
    ) else {
        ogeom_bail!(Construction, "a seam end has no re-solved vertex");
    };
    let Some(curve) = ogeom_algo::surface_iso_u_curve(&prep.surface, column, tol) else {
        ogeom_bail!(
            Construction,
            "the moved surface's iso-curve has no closed form; no seam can \
             be rebuilt"
        );
    };
    // The parameters the re-solved ends land at, by the iso-curve's own
    // closed form; the ends were Newton-polished onto this very surface, so
    // they lie on the curve exactly.
    let along = |p: Point| -> OgeomResult<f64> {
        match &curve {
            Curve::Line(l) => Ok((p - l.axis().location).dot(l.axis().direction.vector())),
            Curve::Circle(c) => {
                let local = c.circle().frame().to_local(p);
                let mut angle = local.y.atan2(local.x);
                if angle < 0.0 {
                    angle += core::f64::consts::TAU;
                }
                Ok(angle)
            }
            // A fitted support's iso-curve is a B-spline: the parameter is
            // found by projection, and the check below says whether the end
            // lies on it.
            _ => Ok(ogeom_algo::project_on_curve(&curve, p, 64, tol)?.parameter),
        }
    };
    let (t_start, t_end) = (along(p_from)?, along(p_to)?);
    // Self-validation instead of trusting the move: a turned support only
    // keeps its column when the turn was built to; the re-solved ends say
    // whether it was.
    // A fitted support holds its column only to the fit's target, and the
    // ends were polished onto the surface, not the column: the slack is the
    // fit's, on a fitted support, and a hundred confusions elsewhere.
    // A drafted support is two fits, the hinge's and its rulings', and at
    // the far end of a ruling their errors add; a few targets' worth is
    // the support's own honesty, not a wrong column.
    let slack = if matches!(prep.surface, SurfaceGeometry::BSpline(_)) {
        (tol.confusion() * 1e3).max(1e-4) * 4.0
    } else {
        tol.confusion() * 100.0
    };
    for (t, p, v) in [(t_start, p_from, &v_from), (t_end, p_to, &v_to)] {
        let off = curve.point_at(t, tol)?.distance(p);
        if off > slack {
            return Ok(None);
        }
        // The end sits on the column to the fit's breadth, and the vertex
        // owns that breadth.
        if off > tol.confusion() {
            model.widen(v, ogeom_core::Tolerance::new(off * 2.0)?)?;
        }
    }
    Ok(Some(if t_start <= t_end {
        make_edge_between(model, curve, (t_start, t_end), &v_from, &v_to, tol)?.shape
    } else {
        // The old seam ran against the iso-curve's own direction: build it
        // the way the curve runs, then hand back the reversed occurrence so
        // the wire's stored orientations still compose.
        make_edge_between(model, curve, (t_end, t_start), &v_to, &v_from, tol)?
            .shape
            .reversed()
    }))
}

/// Assemble a moved face whose wires contain a seam.
///
/// Every ordinary edge gets its exact pcurve recomputed on the moved
/// surface. The seam is the one edge no closed-form projection can answer
/// (it needs a column per side), so its columns carry over from the old
/// face's own seam representation (a same-family move leaves the columns
/// where they were), rebuilt over the rows the moved seam actually spans.
fn assembled_with_seam(
    model: &mut Model,
    prep: &Prepared,
    wires: &[Vec<Shape>],
    new_edges: &FastMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut rings: Vec<Shape> = Vec::with_capacity(wires.len());
    for edges in wires {
        rings.push(ogeom_algo::make_wire(model, edges, tol)?.shape);
    }
    let face = ogeom_algo::make_face(model, prep.surface.clone(), &rings, tol)?.shape;
    let new_surface = {
        let Some(NodeData::Face(data)) = model.node(&face).map(ogeom_topo::TShape::data) else {
            ogeom_bail!(Construction, "the face just built holds no face data");
        };
        data.surface
    };
    let old_surface = {
        let Some(NodeData::Face(data)) = model.node(&prep.shape).map(ogeom_topo::TShape::data)
        else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        data.surface
    };

    let mut done: Vec<TShapeId> = Vec::new();
    for used in explore(model, &prep.shape, Filter::OfType(ShapeType::Edge))? {
        if done.contains(&used.node()) {
            continue;
        }
        done.push(used.node());
        let Some(fresh) = new_edges.get(&used.node()).cloned() else {
            ogeom_bail!(Construction, "a face edge was not rebuilt");
        };
        // A pole has no curve, only its row, and the moved surface keeps
        // its chart: the old row carries over as it is.
        if model
            .node(&fresh)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|d| d.degenerate)
        {
            let Some(EdgeRepr::PCurve { curve, range, .. }) = model
                .node(&used)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.pcurve_for(old_surface, used.location()))
                .cloned()
            else {
                ogeom_bail!(Construction, "a pole has no row on its face");
            };
            let Some(row) = model.geometry().pcurve(curve).cloned() else {
                ogeom_bail!(Dangling, "a pole's row is not in this model");
            };
            ogeom_algo::attach_pcurve(
                model,
                &fresh,
                row,
                new_surface,
                ogeom_topo::Location::identity(),
                range,
            )?;
            continue;
        }
        let (fresh_curve, fresh_range) = {
            let Some(data) = model.node(&fresh).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "a rebuilt edge holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "a rebuilt edge has no curve");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let columns = {
            let Some(data) = model.node(&used).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let mut columns = None;
            for repr in &data.representations {
                if let EdgeRepr::Seam {
                    forward,
                    reversed,
                    surface,
                    range,
                    ..
                } = repr
                    && *surface == old_surface
                {
                    use ogeom_geom::Curve2d as _;
                    let (Some(f), Some(r)) = (
                        model.geometry().pcurve(*forward),
                        model.geometry().pcurve(*reversed),
                    ) else {
                        ogeom_bail!(Dangling, "a seam pcurve is not in this model");
                    };
                    columns = Some((f.point_at(range.0, tol)?.x, r.point_at(range.0, tol)?.x));
                    break;
                }
            }
            columns
        };
        if let Some((forward_col, reversed_col)) = columns {
            // The rows the moved seam spans, from its own rebuilt range,
            // identical to the curve range except a cone's slant rescale.
            let rows = match &prep.surface {
                SurfaceGeometry::Cone(c) => {
                    let cos = c.cone().half_angle().cos();
                    (fresh_range.0 * cos, fresh_range.1 * cos)
                }
                _ => fresh_range,
            };
            let column = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
                Ok(ogeom_geom::Line2d::over(
                    ogeom_math::Axis2::new(
                        ogeom_math::Point2::new(u, 0.0),
                        ogeom_math::Direction2::new(ogeom_math::Vector2::new(0.0, 1.0), tol)?,
                    ),
                    rows.0 - 1.0,
                    rows.1 + 1.0,
                )?
                .into())
            };
            ogeom_algo::attach_seam(
                model,
                &fresh,
                column(forward_col)?,
                column(reversed_col)?,
                new_surface,
                ogeom_topo::Location::identity(),
                rows,
            )?;
        } else {
            // A closed form where one exists; on a fitted support, the
            // projected fit the face builder trusts, the measured offset
            // widening the edge.
            let pcurve = match ogeom_intersect::exact_pcurve_over(
                &fresh_curve,
                fresh_range,
                &prep.surface,
                tol,
            ) {
                Some(exact) => exact,
                None => {
                    let (fitted, _, _, worst_off, _) =
                        ogeom_algo::pcurve_fit::fit_projected_pcurve(
                            &fresh_curve,
                            fresh_range,
                            &prep.surface,
                            tol,
                        )?;
                    if worst_off > tol.confusion() {
                        // The edge owns the offset, and so must the vertices
                        // that bound it: a bound no looser than what it
                        // bounds is the containment rule.
                        let widened = ogeom_core::Tolerance::new(worst_off + tol.confusion())?;
                        model.widen(&fresh, widened)?;
                        if let Some((a, b)) = edge_vertices(model, &fresh)? {
                            model.widen(&a, widened)?;
                            model.widen(&b, widened)?;
                        }
                    }
                    fitted
                }
            };
            ogeom_algo::attach_pcurve(
                model,
                &fresh,
                pcurve,
                new_surface,
                ogeom_topo::Location::identity(),
                fresh_range,
            )?;
        }
    }
    Ok(face)
}

/// Where a seam ending at `vertex` meets the other seat, on the moved
/// supports: the seam's column as an iso-curve on its face's moved surface,
/// pierced through the other face's; `None` where no seam ends here or the
/// column has no curve.
fn seam_end(
    model: &Model,
    faces: &[Shape],
    prepared: &[Prepared],
    vertex: &Shape,
    kept: &[usize],
    at: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    use ogeom_geom::Curve2d as _;
    for (slot, &fi) in kept.iter().enumerate() {
        let other = kept[1 - slot];
        let face = &faces[fi];
        let Some(NodeData::Face(data)) = model.node(face).map(ogeom_topo::TShape::data) else {
            continue;
        };
        let old_surface = data.surface;
        let mut uses: FastMap<TShapeId, usize> = FastMap::default();
        for e in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            *uses.entry(e.node()).or_insert(0) += 1;
        }
        for e in explore_unique(model, face, ShapeType::Edge)? {
            if uses.get(&e.node()).copied().unwrap_or(0) < 2 {
                continue;
            }
            let Some((a, b)) = edge_vertices(model, &e)? else {
                continue;
            };
            if a.node() != vertex.node() && b.node() != vertex.node() {
                continue;
            }
            let Some(edata) = model.node(&e).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let mut column = None;
            for repr in &edata.representations {
                if let EdgeRepr::Seam {
                    forward,
                    surface,
                    range,
                    ..
                } = repr
                    && *surface == old_surface
                    && let Some(pcurve) = model.geometry().pcurve(*forward)
                {
                    column = Some(pcurve.point_at(range.0, tol)?.x);
                    break;
                }
            }
            let Some(column) = column else {
                continue;
            };
            let Some(iso) = ogeom_algo::surface_iso_u_curve(&prepared[fi].surface, column, tol)
            else {
                continue;
            };
            let found = ogeom_intersect::intersect_curve_surface(
                &iso,
                &prepared[other].surface,
                ogeom_intersect::CurveSurfaceOptions::default(),
                tol,
            )?;
            let nearest = found
                .crossings
                .iter()
                .map(|hit| hit.point)
                .min_by(|p, q| p.distance(at).total_cmp(&q.distance(at)));
            if let Some(p) = nearest {
                return Ok(Some(p));
            }
        }
    }
    Ok(None)
}

/// Rebuild an edge only one face owns: the ring a boolean left coincident
/// with a neighbour's twin, or a cone's apex.
///
/// A normal offset moves every point of a face along the face's own normal,
/// so three displaced samples of a ring pin the moved ring exactly; no
/// second face required. `None` when the edge is neither shape.
fn rebuilt_lone_edge(
    model: &mut Model,
    edge: &Shape,
    sides: &[usize],
    constraint: &Displacement<'_>,
    new_vertices: &FastMap<TShapeId, (Shape, Point)>,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    use ogeom_geom::Curve3d as _;

    let (degenerate, curve) = {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let curve = data.curve3d().and_then(|repr| {
            let EdgeRepr::Curve3d { curve, range, .. } = repr else {
                return None;
            };
            model.geometry().curve(*curve).cloned().map(|c| (c, *range))
        });
        (data.degenerate, curve)
    };
    let Some((start, end)) = edge_vertices(model, edge)? else {
        ogeom_bail!(Construction, "a lone edge has no vertices");
    };
    if degenerate {
        // An apex: a rim of no length at the re-solved vertex.
        let Some((vertex, _)) = new_vertices.get(&start.node()) else {
            ogeom_bail!(Construction, "an apex has no re-solved vertex");
        };
        let mut data = EdgeData::new();
        data.degenerate = true;
        return Ok(Some(
            model.add_edge(data, &[vertex.clone(), vertex.clone()])?,
        ));
    }
    let (Some((Curve::Circle(c), range)), true, &[fi]) = (curve, start.node() == end.node(), sides)
    else {
        return Ok(None);
    };
    let circle = c.circle();
    let mut moved_points = Vec::with_capacity(3);
    for k in 0..3 {
        #[allow(clippy::cast_precision_loss, reason = "k is 0..3")]
        let t = (range.1 - range.0).mul_add(k as f64 / 3.0, range.0);
        let p = Curve::Circle(c).point_at(t, tol)?;
        let Some((n, w)) = constraint(model, fi, p)? else {
            return Ok(None);
        };
        moved_points.push(p + n * w);
    }
    // Equally spaced samples average to the centre; the displacement is
    // rotationally symmetric about the ring's own axis, so the moved ring is
    // concentric on it.
    let centre = Point::from_vector(
        moved_points
            .iter()
            .fold(Vector::new(0.0, 0.0, 0.0), |a, p| a + p.to_vector())
            / 3.0,
    );
    let radius = centre.distance(moved_points[0]);
    let reframed = ogeom_math::Circle::new(
        Frame::new(centre, circle.frame().z(), circle.frame().x(), tol)?,
        radius,
        tol,
    )?;
    let moved: Curve = ogeom_geom::CircleCurve::new(reframed).into();
    Ok(Some(make_edge(model, moved, range, tol)?.shape))
}

/// The displacement that puts a point back on every moved plane: solve
/// `x · nᵢ = wᵢ` for the corner's normals, exactly for three, in the least
/// squares sense beyond.
/// Refuse a re-solved corner that does not lie on every moved support it
/// sits on: a plane moved clear of the drum it was tangent to meets it
/// nowhere, and the solve stops on neither. Within the supports' own
/// stated accuracy, as the hinge test allows a fitted support.
fn corner_met(
    prepared: &[Prepared],
    seats: &[usize],
    at: Point,
    tol: Tolerances,
) -> OgeomResult<()> {
    for fi in seats {
        // An analytic support's distance is its own; a projection is
        // clamped to the surface's parameter window and measures the
        // window's edge where the point lies past it.
        let off = match &prepared[*fi].surface {
            SurfaceGeometry::Plane(s) => s.plane().distance_to(at),
            SurfaceGeometry::Cylinder(s) => s.cylinder().distance_to(at),
            SurfaceGeometry::Cone(s) => s.cone().distance_to(at),
            SurfaceGeometry::Sphere(s) => s.sphere().distance_to(at),
            SurfaceGeometry::Torus(s) => s.torus().distance_to(at),
            other => ogeom_algo::project_on_surface(other, at, 32, tol)?.distance,
        };
        if off > (tol.confusion() * 1e3).max(1e-4) {
            ogeom_bail!(
                Construction,
                "the moved faces no longer meet at a corner they shared; the \
                 edit pulls them apart there"
            );
        }
    }
    Ok(())
}

fn solve_corner(normals: &[Vector], amounts: &[f64], tol: Tolerances) -> OgeomResult<Vector> {
    // Normal equations: (NᵀN) x = Nᵀw, 3×3 whatever the seat count.
    let mut a = [[0.0_f64; 3]; 3];
    let mut b = [0.0_f64; 3];
    for (n, w) in normals.iter().zip(amounts) {
        let row = [n.x, n.y, n.z];
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += row[i] * row[j];
            }
            b[i] += row[i] * w;
        }
    }
    // For an edge between two planes the system is rank two; regularize
    // along the null direction (the edge itself), where the displacement is
    // rightly zero.
    if normals.len() == 2 {
        let along = normals[0].cross(normals[1]);
        let m = along.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "an edge between parallel faces has no corner");
        }
        let d = along / m;
        let row = [d.x, d.y, d.z];
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += row[i] * row[j];
            }
        }
    }
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() <= tol.angular() * tol.angular() {
        ogeom_bail!(
            Construction,
            "a corner's faces are too nearly parallel to re-solve"
        );
    }
    let inv = |r: usize, c: usize| -> f64 {
        let (r1, r2) = ((r + 1) % 3, (r + 2) % 3);
        let (c1, c2) = ((c + 1) % 3, (c + 2) % 3);
        (a[c1][r1] * a[c2][r2] - a[c1][r2] * a[c2][r1]) / det
    };
    let mut x = [0.0_f64; 3];
    for (i, xi) in x.iter_mut().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            *xi += inv(i, j) * bj;
        }
    }
    Ok(Vector::new(x[0], x[1], x[2]))
}

/// The edges of `body` between two faces neither of which is `held`, that
/// are convex (`convex` true) or concave: where the walls, grown toward
/// that side, round about the old edge. An edge is convex where the face
/// across it falls away below its neighbour's outward normal.
fn growing_edges(
    model: &Model,
    body: &Shape,
    held: &[TShapeId],
    convex: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    use ogeom_geom::{Curve3d as _, Surface as _};
    // Which faces hold each edge, gathered in one walk over the faces, in
    // the faces' order: asking per edge walked every face each time.
    let mut holders: ogeom_core::FastMap<TShapeId, Vec<(Shape, Shape)>> =
        ogeom_core::FastMap::default();
    for face in ogeom_topo::explore(model, body, ogeom_topo::Filter::OfType(ShapeType::Face))? {
        for e in ogeom_topo::explore(model, &face, ogeom_topo::Filter::OfType(ShapeType::Edge))? {
            holders.entry(e.node()).or_default().push((e, face.clone()));
        }
    }
    let mut out = Vec::new();
    for edge in ogeom_topo::explore_unique(model, body, ShapeType::Edge)? {
        let mut faces: Vec<Shape> = Vec::new();
        for (e, face) in holders.get(&edge.node()).map_or(&[][..], Vec::as_slice) {
            if e.is_same(&edge) && !faces.iter().any(|f| f.is_same(face)) {
                faces.push(face.clone());
            }
        }
        let mut distinct: Vec<Shape> = Vec::new();
        for f in faces {
            if !distinct.iter().any(|d| d.node() == f.node()) {
                distinct.push(f);
            }
        }
        let [f1, f2] = distinct.as_slice() else {
            continue;
        };
        if held.contains(&f1.node()) || held.contains(&f2.node()) {
            continue;
        }
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let p = edge
            .transform(model.datums())?
            .apply(geometry.point_at(f64::midpoint(range.0, range.1), tol)?);
        // The first face's outward normal at the edge.
        let Some(NodeData::Face(face_data)) = model.node(f1).map(ogeom_topo::TShape::data) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(face_data.surface) else {
            continue;
        };
        use ogeom_geom::Transformable as _;
        let placed = surface
            .clone()
            .transformed(&f1.transform(model.datums())?, tol)?;
        let (u, v) = ogeom_algo::project_on_surface(&placed, p, 16, tol)?.parameters;
        let mut n1 = placed.normal_at(u, v, tol)?.vector();
        if f1.orientation() == Orientation::Reversed {
            n1 = -n1;
        }
        // Which side of the first face's plane the other face lies: its
        // mesh point standing furthest off that plane says.
        let mesh = ogeom_mesh::triangulate_face(model, f2, ogeom_mesh::Deflection::default(), tol)?;
        let below = mesh
            .positions
            .iter()
            .map(|q| (*q - p).dot(n1))
            .max_by(|a, b| a.abs().total_cmp(&b.abs()))
            .unwrap_or(0.0);
        if (convex && below < -tol.confusion()) || (!convex && below > tol.confusion()) {
            out.push(edge);
        }
    }
    Ok(out)
}
