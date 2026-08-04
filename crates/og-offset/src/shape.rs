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

use og_algo::{
    Built, History, edge_vertices, make_edge, make_edge_between, make_face_with_pcurves,
    make_revolution_band, make_solid, make_vertex, sew,
};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve, CylinderSurface, LineCurve, PlaneSurface, SurfaceGeometry};
use og_math::{Axis, Circle, Cylinder, Frame, Plane, Point, Vector};
use og_topo::{
    EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};

use std::collections::HashMap;

/// Offset a solid by `offset`: positive grows it, negative shrinks it, and
/// the topology is preserved one-for-one.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a face,
/// edge or vertex falls outside the analytic vocabulary this rebuild speaks
/// (see the module documentation), or the offset collapses the solid.
pub fn offset_shape(
    model: &mut Model,
    solid: &Shape,
    offset: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    if !offset.is_finite() || offset.abs() <= tol.confusion() {
        og_bail!(Construction, "an offset of {offset} moves nothing");
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
) -> OgResult<Built> {
    if !thickness.is_finite() || thickness <= tol.confusion() {
        og_bail!(Construction, "a wall of {thickness} holds nothing");
    }
    let own: Vec<TShapeId> = explore(model, solid, Filter::OfType(ShapeType::Face))?
        .iter()
        .map(Shape::node)
        .collect();
    for face in removed {
        if !own.contains(&face.node()) {
            og_bail!(Construction, "a removed face is not a face of the solid");
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
    let mut result = og_bool::cut(model, solid, &cavity.shape, tol)?;
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
    /// For a planar face, its outward unit normal.
    outward: Option<Vector>,
    /// For a full revolution band, its two ring edges.
    rings: Option<[Shape; 2]>,
}

/// The rebuild under both entry points: every face offset by its own amount,
/// the topology re-derived on the moved surfaces.
fn rebuilt(
    model: &mut Model,
    solid: &Shape,
    amount_of: &dyn Fn(&Shape) -> f64,
    tol: Tolerances,
) -> OgResult<Built> {
    let faces = explore(model, solid, Filter::OfType(ShapeType::Face))?;

    // Move every surface, and note what each face is.
    let mut prepared: Vec<Prepared> = Vec::with_capacity(faces.len());
    for face in &faces {
        let amount = amount_of(face);
        let Some(node) = model.node(face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            og_bail!(Dangling, "face refers to a surface not in this model");
        };
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        // A face's edge list names its seam twice; a face with a seam and two
        // closed rings is a band, rebuilt wholesale.
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

        match surface {
            SurfaceGeometry::Plane(p) => {
                let plane = p.plane();
                let outward = plane.normal().vector() * sign;
                let ((u0, u1), (v0, v1)) = {
                    use og_geom::Surface as _;
                    surface.domain()
                };
                let grow = amount.abs() * 4.0 + 1.0;
                let moved = Plane::new(Frame::new(
                    plane.origin() + outward * amount,
                    plane.normal(),
                    plane.frame().x(),
                    tol,
                )?);
                prepared.push(Prepared {
                    shape: face.clone(),
                    surface: PlaneSurface::over(
                        moved,
                        (u0 - grow, u1 + grow),
                        (v0 - grow, v1 + grow),
                    )?
                    .into(),
                    amount,
                    outward: Some(outward),
                    rings: None,
                });
            }
            SurfaceGeometry::Cylinder(c) => {
                if !has_seam || closed_rings.len() != 2 {
                    og_bail!(
                        Construction,
                        "offsetting a partial cylindrical face needs the \
                         general offset-surface rebuild — see the deferred \
                         table"
                    );
                }
                let cylinder = c.cylinder();
                let grown = cylinder.radius() + amount * sign;
                if grown <= tol.confusion() {
                    og_bail!(Construction, "the offset consumes the cylinder's radius");
                }
                let ((_, _), (v0, v1)) = {
                    use og_geom::Surface as _;
                    surface.domain()
                };
                let reach = amount.abs() + 1.0;
                prepared.push(Prepared {
                    shape: face.clone(),
                    surface: CylinderSurface::new(
                        Cylinder::new(cylinder.frame(), grown, tol)?,
                        (v0 - reach, v1 + reach),
                    )?
                    .into(),
                    amount,
                    outward: None,
                    rings: Some([closed_rings[0].clone(), closed_rings[1].clone()]),
                });
            }
            _ => og_bail!(
                Construction,
                "offsetting a face that is neither a plane nor a cylindrical \
                 band needs the general offset-surface rebuild — see the \
                 deferred table"
            ),
        }
    }

    // Which faces share each edge, seams excluded by their double use.
    let mut edge_faces: HashMap<TShapeId, Vec<usize>> = HashMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for e in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            let entry = edge_faces.entry(e.node()).or_default();
            if !entry.contains(&fi) {
                entry.push(fi);
            }
        }
    }

    // New vertices: each old vertex re-solved where its planes now meet.
    let mut new_vertices: HashMap<TShapeId, (Shape, Point)> = HashMap::new();
    for vertex in explore_unique(model, solid, ShapeType::Vertex)? {
        let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            continue;
        };
        let at = vertex.transform(model.datums())?.apply(data.point);
        // The faces this vertex sits on.
        let mut seats: Vec<usize> = Vec::new();
        for (fi, face) in faces.iter().enumerate() {
            for v in explore(model, face, Filter::OfType(ShapeType::Vertex))? {
                if v.node() == vertex.node() && !seats.contains(&fi) {
                    seats.push(fi);
                }
            }
        }
        if seats.len() < 3 {
            // The anchor of a closed ring edge: it has no corner to re-solve
            // and the rebuilt ring makes its own.
            continue;
        }
        let mut normals: Vec<Vector> = Vec::new();
        let mut amounts: Vec<f64> = Vec::new();
        for fi in &seats {
            let Some(outward) = prepared[*fi].outward else {
                og_bail!(
                    Construction,
                    "a vertex seated on a curved face needs the general \
                     offset rebuild — see the deferred table"
                );
            };
            normals.push(outward);
            amounts.push(prepared[*fi].amount);
        }
        let moved = at + solve_corner(&normals, &amounts, tol)?;
        new_vertices.insert(vertex.node(), (make_vertex(model, moved).shape, moved));
    }

    // New edges on the moved supports, directions preserved.
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
                og_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                og_bail!(Construction, "an edge has no curve to offset");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let built = match &curve {
            Curve::Line(line) => {
                // Between two planes: the direction survives, the anchor
                // re-solves against both moved planes.
                let (Some(na), Some(nb)) = (prepared[sides[0]].outward, prepared[sides[1]].outward)
                else {
                    og_bail!(
                        Construction,
                        "a straight edge on a curved face needs the general \
                         offset rebuild — see the deferred table"
                    );
                };
                let shift = solve_corner(
                    &[na, nb],
                    &[prepared[sides[0]].amount, prepared[sides[1]].amount],
                    tol,
                )?;
                let anchor = line.axis().location + shift;
                let moved: Curve = LineCurve::new(Axis::new(anchor, line.axis().direction)).into();
                // The stored order, whatever this occurrence's orientation:
                // the new edge keeps the old curve's direction, so old wire
                // orientation flags carry over unchanged.
                let forward = if edge.orientation() == Orientation::Reversed {
                    edge.reversed()
                } else {
                    edge.clone()
                };
                let Some((sv, ev)) = edge_vertices(model, &forward)? else {
                    og_bail!(Construction, "a straight edge has no vertices");
                };
                let (Some((v_from, p_from)), Some((v_to, p_to))) = (
                    new_vertices.get(&sv.node()).cloned(),
                    new_vertices.get(&ev.node()).cloned(),
                ) else {
                    og_bail!(Construction, "an edge end has no re-solved vertex");
                };
                let d = line.axis().direction.vector();
                let (t0, t1) = ((p_from - anchor).dot(d), (p_to - anchor).dot(d));
                if t1 <= t0 + tol.parametric() {
                    og_bail!(Construction, "the offset collapses an edge");
                }
                make_edge_between(model, moved, (t0, t1), &v_from, &v_to, tol)?.shape
            }
            Curve::Circle(c) => {
                // An axis-normal ring between a cap plane and a cylindrical
                // band: the centre rides the cap's own motion, the radius is
                // the band's new radius.
                let circle = c.circle();
                let z = circle.frame().z().vector();
                let mut new_radius: Option<f64> = None;
                let mut lift = Vector::new(0.0, 0.0, 0.0);
                for fi in &sides {
                    match (&prepared[*fi].surface, prepared[*fi].outward) {
                        (_, Some(outward)) => {
                            if outward.cross(z).magnitude() > tol.angular() {
                                og_bail!(
                                    Construction,
                                    "a ring's cap is not perpendicular to its \
                                     axis; the general rebuild is deferred"
                                );
                            }
                            lift = outward * prepared[*fi].amount;
                        }
                        (SurfaceGeometry::Cylinder(cy), None) => {
                            let axis = cy.cylinder().axis();
                            if axis.direction.vector().cross(z).magnitude() > tol.angular()
                                || axis.distance_to(circle.centre()) > tol.confusion() * 10.0
                            {
                                og_bail!(
                                    Construction,
                                    "a ring off its wall's axis needs the \
                                     general rebuild — see the deferred table"
                                );
                            }
                            new_radius = Some(cy.cylinder().radius());
                        }
                        _ => og_bail!(
                            Construction,
                            "a ring between these faces needs the general \
                             rebuild — see the deferred table"
                        ),
                    }
                }
                let Some(radius) = new_radius else {
                    og_bail!(
                        Construction,
                        "a circular edge not seated on a cylindrical band \
                         needs the general rebuild — see the deferred table"
                    );
                };
                let moved_circle = Circle::new(
                    Frame::new(
                        circle.centre() + lift,
                        circle.frame().z(),
                        circle.frame().x(),
                        tol,
                    )?,
                    radius,
                    tol,
                )?;
                let moved: Curve = og_geom::CircleCurve::new(moved_circle).into();
                let closed = {
                    let Some((sv, ev)) = edge_vertices(model, &edge)? else {
                        og_bail!(Construction, "a ring has no vertex");
                    };
                    sv.node() == ev.node()
                };
                if !closed {
                    og_bail!(
                        Construction,
                        "a partial circular edge needs the general rebuild — \
                         see the deferred table"
                    );
                }
                make_edge(model, moved, range, tol)?.shape
            }
            _ => og_bail!(
                Construction,
                "offsetting an edge that is neither straight nor circular \
                 needs the general rebuild — see the deferred table"
            ),
        };
        history.modify(&edge, built.clone());
        new_edges.insert(edge.node(), built);
    }

    // Faces: bands wholesale, planar faces wire by wire.
    let mut rebuilt_faces: Vec<Shape> = Vec::with_capacity(prepared.len());
    for prep in &prepared {
        let built = if let Some(rings) = &prep.rings {
            let (Some(lo), Some(hi)) = (
                new_edges.get(&rings[0].node()),
                new_edges.get(&rings[1].node()),
            ) else {
                og_bail!(Construction, "a band's ring was not rebuilt");
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
                        og_bail!(Construction, "a face edge was not rebuilt");
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
    if sewn.shells.len() != 1 || !og_algo::is_shell_closed(model, &sewn.shells[0])? {
        og_bail!(Construction, "the offset solid did not close");
    }
    let built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;

    // The one global guard the local checks cannot give: an offset that
    // moved faces past each other builds a shell that is closed and inside
    // out. Its measured volume is the tell.
    let mass =
        og_algo::volume_properties(model, &built.shape, og_mesh::Deflection::default(), tol)?.mass;
    if !mass.is_finite() || mass <= tol.confusion() {
        og_bail!(Construction, "the offset collapses the solid");
    }

    history.modify(solid, built.shape.clone());
    Ok(Built::new(built.shape, history))
}

/// The displacement that puts a point back on every moved plane: solve
/// `x · nᵢ = wᵢ` for the corner's normals, exactly for three, in the least
/// squares sense beyond.
fn solve_corner(normals: &[Vector], amounts: &[f64], tol: Tolerances) -> OgResult<Vector> {
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
            og_bail!(Construction, "an edge between parallel faces has no corner");
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
        og_bail!(
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
