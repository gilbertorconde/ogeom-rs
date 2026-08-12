//! Offsetting a solid, and the shelling built on it.
//!
//! The topology-preserving offset: every face's surface moves along its own
//! outward normal — a plane translates, a cylinder's radius grows or
//! shrinks — and the topology is rebuilt one-for-one on the moved surfaces.
//! Vertices re-solve where their planes now meet, edges re-derive on the
//! moved supports with their directions and parameterizations preserved, and
//! band faces rebuild through [`make_revolution_band`] so seams stay seams.
//! Corners stay sharp: this is the parallel solid of the intersection join,
//! not the rounded Minkowski body.
//!
//! Shelling is the offset pointed inward and the boolean pointed at the
//! result: the cavity is the inward offset with the *removed* faces left
//! exactly where they were, so it reaches the boundary at the openings —
//! and the cut's same-domain resolution melts the flush faces away, which is
//! what opens the shell.
//!
//! The honest limits, refused by name: faces whose surfaces are not among
//! the five analytics — a spline, revolution, extrusion, trimmed or offset
//! surface has no same-family parallel to move to — edges that are neither
//! straight nor circular, vertices whose seats leave them under-determined,
//! and offsets that collapse the solid. Partial bands, tori and cones all
//! move; what stops the rebuild is the edge between two moved surfaces that
//! no longer meets in a line or a circle.

use ogeom_algo::{
    Built, History, edge_vertices, make_edge, make_edge_between, make_face_with_pcurves,
    make_revolution_band, make_solid, make_vertex, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, CylinderSurface, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Cylinder, Frame, Plane, Point, Vector};
use ogeom_topo::{
    EdgeData, EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};

use std::collections::HashMap;

/// The displacement constraint one face puts on a point of itself.
type Displacement<'a> = dyn Fn(&Model, usize, Point) -> OgeomResult<Option<(Vector, f64)>> + 'a;

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
    rebuilt(model, solid, &|_| offset, &|_| None, tol)
}

/// Hollow a solid into a shell of the given wall `thickness`, opening it at
/// the `removed` faces.
///
/// Two constructions, chosen by the opening's neighbours. When a removed
/// face meets every neighbour across a corner, the cavity is the inward
/// offset of every kept face with the removed faces left in place,
/// subtracted through the boolean — the flush faces melt away, which is
/// what opens the shell. When a removed face has a *tangent* neighbour — a
/// blend melting into the face it rounds — leaving it in place would tear
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
    if !thickness.is_finite() || thickness <= tol.confusion() {
        ogeom_bail!(Construction, "a wall of {thickness} holds nothing");
    }
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
        let cavity = rebuilt(
            model,
            solid,
            &|face| {
                if skip.contains(&face.node()) {
                    0.0
                } else {
                    -thickness
                }
            },
            &|_| None,
            tol,
        )?;
        let mut result = ogeom_bool::cut(model, solid, &cavity.shape, tol)?;
        for face in removed {
            result.history.delete(face);
        }
        return Ok(result);
    }

    // The tangent construction: everything moves together — which is what
    // keeps the tangencies intact — and each opening is drilled back out by
    // extruding its cavity image through where the wall used to be.
    let cavity = rebuilt(model, solid, &|_| -thickness, &|_| None, tol)?;
    let mut tool = cavity.shape.clone();
    for face in removed {
        let outward = {
            let Some(NodeData::Face(data)) = model.node(face).map(ogeom_topo::TShape::data) else {
                ogeom_bail!(Construction, "face node holds no face data");
            };
            let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
                ogeom_bail!(
                    Construction,
                    "a tangent opening must be planar; a curved opening needs \
                     the general rebuild — docs/PARITY.md, offset.shell-thicken"
                );
            };
            let mut normal = p.plane().normal().vector();
            if face.orientation() == Orientation::Reversed {
                normal = -normal;
            }
            normal
        };
        let [image] = cavity.history.modified(face) else {
            ogeom_bail!(Construction, "a removed face has no single cavity image");
        };
        let punch =
            ogeom_algo::make_prism(model, &image.clone(), outward * (2.0 * thickness), tol)?;
        tool = ogeom_bool::fuse(model, &tool, &punch.shape, tol)?.shape;
    }
    let mut result = ogeom_bool::cut(model, solid, &tool, tol)?;
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
    /// For a full revolution band — seam and two closed rings — the rings.
    rings: Option<[Shape; 2]>,
    /// Whether `instead_of` supplied the surface: a *turned* move, under
    /// which the chart is not preserved and pcurves cannot be carried.
    turned: bool,
}

/// The rebuild under both entry points: every face offset by its own amount,
/// the topology re-derived on the moved surfaces.
///
/// One rule serves every element. A surface moves along its own normal — a
/// plane translates, a revolution surface's radius grows — which makes the
/// *displacement* constraint at any point of it exactly planar: normal
/// there, offset amount along it. Vertices solve those constraints in the
/// least-squares sense and then Newton-polish onto the moved surfaces
/// themselves, edges re-derive from the moved pair — a line from its planes'
/// constraints, a circle from the pair's analytic intersection re-framed on
/// its old axes so parameters and orientations carry — and faces rebuild
/// wire by wire with exact pcurves, or wholesale through
/// [`make_revolution_band`] where a seam says the face wraps.
/// Rebuild a solid's topology on moved supports.
///
/// `amount_of` says how far each face travels along its own outward normal;
/// `instead_of` may hand back a surface to use *in place* of that move,
/// which is how an operation that turns a face rather than translating it —
/// a draft — rides the same rebuild. The two are exclusive per face: a
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
        let mut counts: HashMap<TShapeId, usize> = HashMap::new();
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

        let grow = amount.abs() * 4.0 + 1.0;
        let replacement = instead_of(face);
        let turned = replacement.is_some();
        let moved: SurfaceGeometry = if let Some(given) = replacement {
            given
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
                    ogeom_geom::SphereSurface::new(ogeom_math::Sphere::centred(
                        sphere.centre(),
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
                     the rebuild does not yet speak — docs/PARITY.md, offset.shell-thicken"
                ),
            }
        };
        prepared.push(Prepared {
            shape: face.clone(),
            surface: moved,
            amount,
            sign,
            rings: if has_seam && closed_rings.len() == 2 {
                Some([closed_rings[0].clone(), closed_rings[1].clone()])
            } else {
                None
            },
            turned,
        });
    }

    // Which faces meet each edge, seams excluded by their double use.
    let mut edge_faces: HashMap<TShapeId, Vec<usize>> = HashMap::new();
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
    // the moved surfaces themselves — the tangent-plane answer is exact for
    // planes and off by the surfaces' own curvature otherwise.
    let mut new_vertices: HashMap<TShapeId, (Shape, Point)> = HashMap::new();
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
            // A cone's apex has no normal to offer — the projection there is
            // degenerate — but the parallel cone knows exactly where its own
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
            // answer — no corner to solve, nothing to polish.
            let moved = at + normals[0] * amounts[0];
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
        new_vertices.insert(vertex.node(), (make_vertex(model, moved).shape, moved));
    }

    // How many times each edge occurs across all faces — a seam is one face
    // using an edge twice, which face-deduplicated sides cannot see.
    let mut edge_uses: HashMap<TShapeId, usize> = HashMap::new();
    for face in &faces {
        for e in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            *edge_uses.entry(e.node()).or_insert(0) += 1;
        }
    }

    // New edges on the moved supports.
    let mut new_edges: HashMap<TShapeId, Shape> = HashMap::new();
    let mut history = History::new();
    for edge in explore_unique(model, solid, ShapeType::Edge)? {
        let sides = edge_faces.get(&edge.node()).cloned().unwrap_or_default();
        if sides.len() != 2 {
            if edge_uses.get(&edge.node()).copied().unwrap_or(0) >= 2 {
                // A seam. A band face rebuilds its own; a face assembled wire
                // by wire — a band a boolean split into arc rings — needs the
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
                // or turned — an offset leaves the direction alone and this
                // reproduces it, a draft does not and this follows it —
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
            Curve::Circle(c) => {
                // The moved pair's own analytic intersection, taken in the
                // circle's old frame so parameters and orientations carry.
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
            _ => ogeom_bail!(
                Construction,
                "offsetting an edge that is neither straight nor circular \
                 needs the general rebuild — docs/PARITY.md, offset.shell-thicken"
            ),
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
            let mut face_uses: HashMap<TShapeId, usize> = HashMap::new();
            for wire in explore(model, &prep.shape, Filter::OfType(ShapeType::Wire))? {
                let mut edges: Vec<Shape> = Vec::new();
                for used in explore(model, &wire, Filter::OfType(ShapeType::Edge))? {
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
                // A seam in a wire-assembled face — a band a boolean split
                // into arc rings. The chart survives a same-family move
                // untouched, so every pcurve carries over verbatim; a turned
                // face cannot make that claim.
                if prep.turned {
                    ogeom_bail!(
                        Construction,
                        "a turned face with a seam cannot carry its pcurves; \
                         the rebuild has no chart to preserve"
                    );
                }
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
    let built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;

    // The one global guard the local checks cannot give: an offset that
    // moved faces past each other builds a shell that is closed and inside
    // out. Its measured volume is the tell.
    // The guard meshes at the default deflection, and a thin tangential
    // cusp — a small blend meeting its face — can defeat that resolution
    // without anything being wrong. One finer retry separates a mesh that
    // cannot see the cusp from a solid that is genuinely inside out.
    let mut mass = None;
    for chord in [ogeom_mesh::Deflection::default().chord, 1e-4] {
        let deflection = ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        };
        if let Ok(props) = ogeom_algo::volume_properties(model, &built.shape, deflection, tol) {
            mass = Some(props.mass);
            break;
        }
    }
    let Some(mass) = mass else {
        ogeom_bail!(
            Construction,
            "the offset solid's mesh does not close at any tried resolution"
        );
    };
    if !mass.is_finite() || mass <= tol.confusion() {
        ogeom_bail!(Construction, "the offset collapses the solid");
    }

    history.modify(solid, built.shape.clone());
    Ok(Built::new(built.shape, history))
}

/// The displacement that puts a point back on every moved plane: solve
/// `x · nᵢ = wᵢ` for the corner's normals, exactly for three, in the least
/// squares sense beyond.
/// Rebuild a seam for a face assembled wire by wire: the same iso-column on
/// the moved surface, over the same rows.
///
/// A same-family move preserves the chart — every point travels along its
/// own normal without changing its parameters — so the moved seam sits at
/// the column the old one's own pcurves state, between the re-solved end
/// vertices. `None` when the old edge carries no seam representation on this
/// face's surface.
fn rebuilt_seam_edge(
    model: &mut Model,
    edge: &Shape,
    prep: &Prepared,
    new_vertices: &HashMap<TShapeId, (Shape, Point)>,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    use ogeom_geom::Curve2d as _;

    if prep.turned {
        return Ok(None);
    }
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
    let Some((sv, ev)) = edge_vertices(model, edge)? else {
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
    // closed form — the ends were Newton-polished onto this very surface, so
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
            _ => ogeom_bail!(
                Construction,
                "the moved seam's iso-curve is neither straight nor circular"
            ),
        }
    };
    let (t_start, t_end) = (along(p_from)?, along(p_to)?);
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
/// surface. The seam is the one edge no closed-form projection can answer —
/// it needs a column per side — so its columns carry over from the old
/// face's own seam representation (a same-family move leaves the columns
/// where they were), rebuilt over the rows the moved seam actually spans.
fn assembled_with_seam(
    model: &mut Model,
    prep: &Prepared,
    wires: &[Vec<Shape>],
    new_edges: &HashMap<TShapeId, Shape>,
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
            // The rows the moved seam spans, from its own rebuilt range —
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
            let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&fresh_curve, &prep.surface, tol)
            else {
                ogeom_bail!(
                    Construction,
                    "a moved face edge has no closed-form pcurve on its surface"
                );
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

/// Rebuild an edge only one face owns: the ring a boolean left coincident
/// with a neighbour's twin, or a cone's apex.
///
/// A normal offset moves every point of a face along the face's own normal,
/// so three displaced samples of a ring pin the moved ring exactly — no
/// second face required. `None` when the edge is neither shape.
fn rebuilt_lone_edge(
    model: &mut Model,
    edge: &Shape,
    sides: &[usize],
    constraint: &Displacement<'_>,
    new_vertices: &HashMap<TShapeId, (Shape, Point)>,
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

fn solve_corner(normals: &[Vector], amounts: &[f64], tol: Tolerances) -> OgeomResult<Vector> {
    // Normal equations: (NᵀN) x = Nᵀw — 3×3 whatever the seat count.
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
