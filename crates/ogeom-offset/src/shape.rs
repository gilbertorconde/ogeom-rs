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
//! The honest limits, refused by name: faces that are not planes or full
//! cylindrical bands, edges that are not straight or axis-normal circles,
//! vertices not seated on three planes, and offsets that collapse the solid.

use ogeom_algo::{
    Built, History, edge_vertices, make_edge, make_edge_between, make_face_with_pcurves,
    make_revolution_band, make_solid, make_vertex, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, CylinderSurface, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Cylinder, Frame, Plane, Point, Vector};
use ogeom_topo::{
    EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};

use std::collections::HashMap;

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
    rebuilt(model, solid, &|_| offset, tol)
}

/// Hollow a solid into a shell of the given wall `thickness`, opening it at
/// the `removed` faces.
///
/// The cavity is the inward offset of every kept face with the removed faces
/// left in place, subtracted through the boolean — so the history reads as
/// the cut it is, with the removed faces deleted.
///
/// # Errors
///
/// As [`offset_shape`], and additionally if `thickness` is not a usable
/// length or a removed face is not a face of `solid`.
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
        tol,
    )?;
    let mut result = ogeom_bool::cut(model, solid, &cavity.shape, tol)?;
    for face in removed {
        result.history.delete(face);
    }
    Ok(result)
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
fn rebuilt(
    model: &mut Model,
    solid: &Shape,
    amount_of: &dyn Fn(&Shape) -> f64,
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
        let moved: SurfaceGeometry = match surface {
            SurfaceGeometry::Plane(p) => {
                let plane = p.plane();
                let ((u0, u1), (v0, v1)) = surface.domain();
                let shifted = Plane::new(Frame::new(
                    plane.origin() + plane.normal().vector() * (sign * amount),
                    plane.normal(),
                    plane.frame().x(),
                    tol,
                )?);
                PlaneSurface::over(shifted, (u0 - grow, u1 + grow), (v0 - grow, v1 + grow))?.into()
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
                let grown =
                    (sign * amount / cone.half_angle().cos()).mul_add(1.0, cone.reference_radius());
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
                "offsetting a face on this surface needs a construction the \
                 rebuild does not yet speak — see the deferred table"
            ),
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
        if seats.len() < 2 {
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
        if normals.len() < 2 {
            ogeom_bail!(
                Construction,
                "a vertex with fewer than two independent seats cannot be \
                 re-solved"
            );
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

    // New edges on the moved supports.
    let mut new_edges: HashMap<TShapeId, Shape> = HashMap::new();
    let mut history = History::new();
    for edge in explore_unique(model, solid, ShapeType::Edge)? {
        let sides = edge_faces.get(&edge.node()).cloned().unwrap_or_default();
        if sides.len() != 2 {
            // A seam: its band rebuilds it.
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
            Curve::Line(line) => {
                // The direction survives; the anchor solves the two faces'
                // displacement constraints — one constraint where the faces
                // are tangent along the line, two where they cross there.
                let mid = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
                let mut normals: Vec<Vector> = Vec::new();
                let mut amounts: Vec<f64> = Vec::new();
                for fi in &sides {
                    let Some((n, w)) = constraint(model, *fi, mid)? else {
                        continue;
                    };
                    if let Some(k) = normals
                        .iter()
                        .position(|m| m.cross(n).magnitude() <= tol.angular().max(1e-6))
                    {
                        if (amounts[k] - w).abs() > tol.confusion() {
                            ogeom_bail!(
                                Construction,
                                "tangent faces move a shared edge by different \
                                 amounts; the offset tears it"
                            );
                        }
                        continue;
                    }
                    normals.push(n);
                    amounts.push(w);
                }
                let shift = if normals.len() == 1 {
                    normals[0] * amounts[0]
                } else {
                    solve_corner(&normals, &amounts, tol)?
                };
                let anchor = line.axis().location + shift;
                let moved: Curve = LineCurve::new(Axis::new(anchor, line.axis().direction)).into();
                let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                    ogeom_bail!(Construction, "a straight edge has no vertices");
                };
                let (Some((v_from, p_from)), Some((v_to, p_to))) = (
                    new_vertices.get(&sv.node()).cloned(),
                    new_vertices.get(&ev.node()).cloned(),
                ) else {
                    ogeom_bail!(Construction, "an edge end has no re-solved vertex");
                };
                let d = line.axis().direction.vector();
                let (t0, t1) = ((p_from - anchor).dot(d), (p_to - anchor).dot(d));
                if t1 <= t0 + tol.parametric() {
                    ogeom_bail!(Construction, "the offset collapses an edge");
                }
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
                 needs the general rebuild — see the deferred table"
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
            for wire in explore(model, &prep.shape, Filter::OfType(ShapeType::Wire))? {
                let mut edges: Vec<Shape> = Vec::new();
                for used in explore(model, &wire, Filter::OfType(ShapeType::Edge))? {
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
            let face = make_face_with_pcurves(model, prep.surface.clone(), &wires, tol)?.shape;
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
