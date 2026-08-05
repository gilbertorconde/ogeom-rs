//! The recognizers, and the shared measurements they stand on.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Surface as _, SurfaceGeometry};
use ogeom_math::{Axis, Point, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, TShapeId, explore};

/// Whether a hole goes through or stops.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HoleKind {
    /// Open at both ends.
    Through,
    /// Capped by a floor or a drill tip.
    Blind,
}

/// A recognized hole: a chain of coaxial concave revolution faces.
#[derive(Debug, Clone)]
pub struct Hole {
    /// The bore walls, cylinders and cones alike, in chain order.
    pub faces: Vec<Shape>,
    /// The common axis.
    pub axis: Axis,
    /// The smallest bore radius in the chain.
    pub radius: f64,
    /// Through or blind.
    pub kind: HoleKind,
    /// Whether the chain steps between cylinder radii — a counterbore.
    pub counterbored: bool,
    /// Whether a cone opens the chain at an entry — a countersink.
    pub countersunk: bool,
}

/// A recognized blend: a face tangent to its neighbours on both sides.
#[derive(Debug, Clone)]
pub struct Fillet {
    /// The blend face — a partial cylinder along a straight edge, or a
    /// torus band around a circular one.
    pub face: Shape,
    /// The rolling-ball radius.
    pub radius: f64,
    /// `true` for a concave blend filling an inside corner; `false` for a
    /// convex round easing an outside one.
    pub concave: bool,
}

/// A recognized chamfer: a planar bevel tangent to neither neighbour.
#[derive(Debug, Clone)]
pub struct Chamfer {
    /// The bevel face.
    pub face: Shape,
}

/// A recognized depression: a floor whose whole boundary folds inward.
#[derive(Debug, Clone)]
pub struct Pocket {
    /// The planar floor.
    pub floor: Shape,
    /// The faces rising from its boundary.
    pub walls: Vec<Shape>,
    /// Whether the floor is an obround — two parallel lines closed by two
    /// arcs — which is what a slot leaves behind.
    pub slot: bool,
}

/// One recognized feature.
#[derive(Debug, Clone)]
pub enum Feature {
    /// A hole.
    Hole(Hole),
    /// A fillet or round.
    Fillet(Fillet),
    /// A chamfer.
    Chamfer(Chamfer),
    /// A pocket or slot.
    Pocket(Pocket),
}

/// Tangency threshold between face normals at a shared edge, radians —
/// loose enough for a real file's slop (the corpus carries hundredths of a
/// degree), still fifty times under any deliberate chamfer angle.
const TANGENT_ANGLE: f64 = 1e-2;

/// Recognize features on a solid.
///
/// Holes claim their faces first, then fillets, chamfers and pockets over
/// what remains, so a bore is a hole rather than a curious pocket with a
/// cylindrical wall. Faces already claimed by a feature are not offered to
/// the later recognizers as seeds, though they may still appear as another
/// feature's walls.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `shape`
/// is not a solid, or its tessellation — the convexity oracle — cannot be
/// built.
pub fn recognize(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<Feature>> {
    if model.kind_of(shape)? != ShapeType::Solid {
        ogeom_bail!(Construction, "features are recognized on solids");
    }
    let scene = Scene::build(model, shape, tol)?;
    let mut features = Vec::new();
    let mut claimed: Vec<TShapeId> = Vec::new();

    holes(&scene, &mut features, &mut claimed)?;
    fillets(&scene, &mut features, &mut claimed)?;
    chamfers(&scene, &mut features, &mut claimed)?;
    pockets(&scene, &mut features, &mut claimed)?;

    Ok(features)
}

// --- the shared measurements ------------------------------------------------

/// The solid, measured once: faces with their surfaces, adjacency over
/// edges, and the parity oracle for convexity.
struct Scene<'m> {
    model: &'m Model,
    tol: Tolerances,
    faces: Vec<FaceInfo>,
    /// Edge node to the indices of the faces using it.
    adjacency: std::collections::HashMap<TShapeId, Vec<usize>>,
}

struct FaceInfo {
    shape: Shape,
    surface: SurfaceGeometry,
    /// Edge nodes per wire, first wire outer.
    wires: Vec<Vec<Shape>>,
}

impl<'m> Scene<'m> {
    fn build(model: &'m Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Self> {
        let mut faces = Vec::new();
        let mut adjacency: std::collections::HashMap<TShapeId, Vec<usize>> =
            std::collections::HashMap::new();
        for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "face is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "face node holds no face data");
            };
            let Some(surface) = model.geometry().surface(data.surface) else {
                ogeom_bail!(Dangling, "face refers to a surface not in this model");
            };
            let surface = surface.clone();
            let mut wires = Vec::new();
            for wire in model.ordered_children_of(&face)? {
                wires.push(model.ordered_children_of(&wire)?);
            }
            let index = faces.len();
            for wire in &wires {
                for edge in wire {
                    let uses = adjacency.entry(edge.node()).or_default();
                    if !uses.contains(&index) {
                        uses.push(index);
                    }
                }
            }
            faces.push(FaceInfo {
                shape: face,
                surface,
                wires,
            });
        }

        Ok(Self {
            model,
            tol,
            faces,
            adjacency,
        })
    }

    /// The outward normal of a face at a point on it.
    fn outward_normal(&self, face: usize, at: Point) -> Option<Vector> {
        let info = &self.faces[face];
        let placement = info.shape.transform(self.model.datums()).ok()?;
        let local = placement.inverse().ok()?.apply(at);
        let mut projection =
            ogeom_algo::project_on_surface(&info.surface, local, 16, self.tol).ok()?;
        if projection.distance > self.tol.confusion() * 1e3 {
            // The point is *on* the face; a projection that says otherwise
            // fell into the wrong basin of a torus or sphere. Seed denser
            // before trusting the normal it implies.
            let denser = ogeom_algo::project_on_surface(&info.surface, local, 96, self.tol).ok()?;
            if denser.distance < projection.distance {
                projection = denser;
            }
        }
        let (u, v) = projection.parameters;
        let normal = info.surface.normal_at(u, v, self.tol).ok()?;
        let world = placement.apply_vector(normal.vector());
        Some(
            if info.shape.orientation() == ogeom_topo::Orientation::Reversed {
                -world
            } else {
                world
            },
        )
    }

    /// The midpoint of an edge's curve, placed.
    fn edge_midpoint(&self, edge: &Shape) -> Option<Point> {
        use ogeom_geom::Curve3d as _;
        let data = self.model.node(edge)?.data().as_edge()?;
        let EdgeRepr::Curve3d { curve, range, .. } = data.curve3d()? else {
            return None;
        };
        let geometry = self.model.geometry().curve(*curve)?;
        let p = geometry
            .point_at(f64::midpoint(range.0, range.1), self.tol)
            .ok()?;
        let placement = edge.transform(self.model.datums()).ok()?;
        Some(placement.apply(p))
    }

    /// How the solid folds at an edge between two faces.
    ///
    /// The classical reading: `t` is the edge as `f1` traverses it —
    /// orientation composition hands it over directly — and the sign of
    /// `(n1 × n2) · t` says which way the surface turns. Positive is an
    /// outside edge, negative an inside one; parallel normals are a smooth
    /// join and have no fold to name.
    fn edge_fold(&self, edge: &Shape, f1: usize, f2: usize) -> Option<Fold> {
        let p = self.edge_midpoint(edge)?;
        let n1 = self.outward_normal(f1, p)?;
        let n2 = self.outward_normal(f2, p)?;
        if n1.cross(n2).magnitude() <= TANGENT_ANGLE {
            return Some(Fold::Smooth);
        }
        let occurrence = self.faces[f1]
            .wires
            .iter()
            .flatten()
            .find(|e| e.node() == edge.node())?;
        let tangent = self.edge_tangent(edge)?;
        let t = if occurrence.orientation() == ogeom_topo::Orientation::Reversed {
            -tangent
        } else {
            tangent
        };
        Some(if n1.cross(n2).dot(t) > 0.0 {
            Fold::Convex
        } else {
            Fold::Concave
        })
    }

    /// The direction of an edge's curve at its middle.
    fn edge_tangent(&self, edge: &Shape) -> Option<Vector> {
        use ogeom_geom::Curve3d as _;
        let data = self.model.node(edge)?.data().as_edge()?;
        let EdgeRepr::Curve3d { curve, range, .. } = data.curve3d()? else {
            return None;
        };
        let geometry = self.model.geometry().curve(*curve)?;
        let mid = f64::midpoint(range.0, range.1);
        let h = (range.1 - range.0).abs().max(1e-9) * 1e-4;
        let a = geometry.point_at(mid - h, self.tol).ok()?;
        let b = geometry.point_at(mid + h, self.tol).ok()?;
        let placement = edge.transform(self.model.datums()).ok()?;
        let d = placement.apply(b) - placement.apply(a);
        let magnitude = d.magnitude();
        (magnitude > f64::MIN_POSITIVE).then(|| d / magnitude)
    }

    /// The other face across an edge from `face`, where the edge is
    /// manifold.
    fn other_face(&self, edge: &Shape, face: usize) -> Option<usize> {
        let uses = self.adjacency.get(&edge.node())?;
        uses.iter().copied().find(|&f| f != face)
    }

    /// Whether an edge appears twice within one face's wires — a seam,
    /// which is what a full revolution face carries.
    fn is_seam_in(&self, face: usize, edge: &Shape) -> bool {
        self.faces[face]
            .wires
            .iter()
            .flatten()
            .filter(|e| e.node() == edge.node())
            .count()
            == 2
    }

    /// Whether a face wraps its surface's full period: some edge of it is
    /// a seam, or a wire of it is a single closed ring.
    fn is_full_revolution(&self, face: usize) -> bool {
        let info = &self.faces[face];
        let seam = info
            .wires
            .iter()
            .flatten()
            .any(|e| self.is_seam_in(face, e));
        let ring_wires = info
            .wires
            .iter()
            .all(|wire| wire.len() == 1 || wire.iter().any(|e| self.is_seam_in(face, e)));
        seam || (info.wires.len() >= 2 && ring_wires)
    }

    /// The line-curve edges of a face — a cylinder's rulings.
    fn straight_edges(&self, face: usize) -> Vec<Shape> {
        self.faces[face]
            .wires
            .iter()
            .flatten()
            .filter(|edge| {
                self.model
                    .node(edge)
                    .and_then(|n| n.data().as_edge())
                    .and_then(|d| d.curve3d())
                    .and_then(|r| match r {
                        EdgeRepr::Curve3d { curve, .. } => self.model.geometry().curve(*curve),
                        _ => None,
                    })
                    .is_some_and(|c| matches!(c, Curve::Line(_)))
            })
            .cloned()
            .collect()
    }
}

/// How a solid folds along an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Fold {
    /// The faces meet tangentially.
    Smooth,
    /// Material angle under a half turn — an outside edge.
    Convex,
    /// Material angle over a half turn — an inside edge.
    Concave,
}

// --- holes -------------------------------------------------------------------

/// The revolution axis of a cylinder or cone surface.
fn revolution_axis(surface: &SurfaceGeometry) -> Option<Axis> {
    match surface {
        SurfaceGeometry::Cylinder(c) => Some(c.cylinder().frame().axis()),
        SurfaceGeometry::Cone(c) => Some(c.cone().frame().axis()),
        _ => None,
    }
}

/// Whether a full revolution face is a bore: its outward normal points
/// toward its own axis.
fn is_bore(scene: &Scene, face: usize) -> bool {
    let info = &scene.faces[face];
    let Some(axis) = revolution_axis(&info.surface) else {
        return false;
    };
    // Sample at any boundary edge midpoint.
    let Some(p) = info
        .wires
        .iter()
        .flatten()
        .find_map(|e| scene.edge_midpoint(e))
    else {
        return false;
    };
    let Some(n) = scene.outward_normal(face, p) else {
        return false;
    };
    let to_axis = axis.location
        + axis.direction.vector() * (p - axis.location).dot(axis.direction.vector())
        - p;
    let magnitude = to_axis.magnitude();
    magnitude > scene.tol.confusion() && n.dot(to_axis / magnitude) > 0.5
}

fn coaxial(a: &Axis, b: &Axis, tol: Tolerances) -> bool {
    a.direction.vector().cross(b.direction.vector()).magnitude() < 1e-6
        && (b.location - a.location)
            .cross(a.direction.vector())
            .magnitude()
            < tol.confusion() * 1e3
}

/// How many wires of a planar face are closed rings centred on `axis` —
/// zero when the face is not such a plane at all.
fn ring_plane_wires(scene: &Scene, face: usize, axis: &Axis) -> usize {
    let info = &scene.faces[face];
    let SurfaceGeometry::Plane(p) = &info.surface else {
        return 0;
    };
    if p.plane()
        .frame()
        .z()
        .vector()
        .cross(axis.direction.vector())
        .magnitude()
        > 1e-6
    {
        return 0;
    }
    // Every edge an on-axis circle — whole rings or the arcs a boolean
    // may have split them into.
    let all_rings = info.wires.iter().all(|wire| {
        !wire.is_empty()
            && wire.iter().all(|edge| {
                scene
                    .model
                    .node(edge)
                    .and_then(|n| n.data().as_edge())
                    .and_then(|d| d.curve3d())
                    .and_then(|r| match r {
                        EdgeRepr::Curve3d { curve, .. } => scene.model.geometry().curve(*curve),
                        _ => None,
                    })
                    .is_some_and(|c| match c {
                        Curve::Circle(circle) => {
                            (circle.circle().centre() - axis.location)
                                .cross(axis.direction.vector())
                                .magnitude()
                                < scene.tol.confusion() * 1e3
                        }
                        _ => false,
                    })
            })
    });
    if all_rings { info.wires.len() } else { 0 }
}

fn holes(
    scene: &Scene,
    features: &mut Vec<Feature>,
    claimed: &mut Vec<TShapeId>,
) -> OgeomResult<()> {
    // Seeds: full-revolution concave cylinders.
    let seeds: Vec<usize> = (0..scene.faces.len())
        .filter(|&f| {
            matches!(scene.faces[f].surface, SurfaceGeometry::Cylinder(_))
                && scene.is_full_revolution(f)
                && is_bore(scene, f)
                && !claimed.contains(&scene.faces[f].shape.node())
        })
        .collect();

    let mut grouped = vec![false; scene.faces.len()];
    for seed in seeds {
        if grouped[seed] {
            continue;
        }
        let Some(seed_axis) = revolution_axis(&scene.faces[seed].surface) else {
            continue;
        };
        // Flood along shared edges over coaxial concave revolution faces —
        // and across the annular shoulders between them, which is what a
        // counterbore's step is: a plane whose every boundary is a ring on
        // the same axis, with more hole on the far side.
        let mut chain = vec![seed];
        grouped[seed] = true;
        let mut cursor = 0;
        while cursor < chain.len() {
            let here = chain[cursor];
            cursor += 1;
            for edge in scene.faces[here].wires.iter().flatten() {
                let Some(next) = scene.other_face(edge, here) else {
                    continue;
                };
                if grouped[next] || chain.contains(&next) {
                    continue;
                }
                if let Some(axis) = revolution_axis(&scene.faces[next].surface) {
                    if coaxial(&seed_axis, &axis, scene.tol)
                        && scene.is_full_revolution(next)
                        && is_bore(scene, next)
                    {
                        grouped[next] = true;
                        chain.push(next);
                    }
                } else if ring_plane_wires(scene, next, &seed_axis) >= 2 {
                    // A shoulder, not a floor: an annulus bridges deeper.
                    grouped[next] = true;
                    chain.push(next);
                }
            }
        }

        // The chain's open ends: ring edges adjoining faces outside it.
        let mut entries = 0;
        let mut caps = 0;
        for &f in &chain {
            for edge in scene.faces[f].wires.iter().flatten() {
                if scene.is_seam_in(f, edge) {
                    continue;
                }
                let Some(outside) = scene.other_face(edge, f) else {
                    continue;
                };
                if chain.contains(&outside) {
                    continue;
                }
                let outside_info = &scene.faces[outside];
                if matches!(outside_info.surface, SurfaceGeometry::Plane(_)) {
                    // Entry if the shared ring is a hole loop of the plane —
                    // any wire but its first — and a cap if it is the
                    // plane's own outer boundary.
                    let in_outer = outside_info
                        .wires
                        .first()
                        .is_some_and(|w| w.iter().any(|e| e.node() == edge.node()));
                    if in_outer {
                        caps += 1;
                    } else {
                        entries += 1;
                    }
                } else {
                    caps += 1;
                }
            }
        }

        let cylinders: Vec<f64> = chain
            .iter()
            .filter_map(|&f| match &scene.faces[f].surface {
                SurfaceGeometry::Cylinder(c) => Some(c.cylinder().radius()),
                _ => None,
            })
            .collect();
        let cones = chain
            .iter()
            .any(|&f| matches!(scene.faces[f].surface, SurfaceGeometry::Cone(_)));
        let Some(radius) = cylinders.iter().copied().reduce(f64::min) else {
            continue;
        };
        let mut radii = cylinders.clone();
        radii.sort_by(f64::total_cmp);
        radii.dedup_by(|a, b| (*a - *b).abs() <= scene.tol.confusion() * 10.0);

        let kind = if entries >= 2 && caps == 0 {
            HoleKind::Through
        } else {
            HoleKind::Blind
        };
        for &f in &chain {
            claimed.push(scene.faces[f].shape.node());
        }
        features.push(Feature::Hole(Hole {
            faces: chain
                .iter()
                .map(|&f| scene.faces[f].shape.clone())
                .collect(),
            axis: seed_axis,
            radius,
            kind,
            counterbored: radii.len() > 1,
            countersunk: cones,
        }));
    }
    Ok(())
}

// --- fillets -----------------------------------------------------------------

fn fillets(
    scene: &Scene,
    features: &mut Vec<Feature>,
    claimed: &mut Vec<TShapeId>,
) -> OgeomResult<()> {
    for f in 0..scene.faces.len() {
        if claimed.contains(&scene.faces[f].shape.node()) {
            continue;
        }
        let (radius, tangent_edges): (f64, Vec<Shape>) = match &scene.faces[f].surface {
            // A straight-edge blend: a partial cylinder, tangent across its
            // rulings.
            SurfaceGeometry::Cylinder(c) if !scene.is_full_revolution(f) => {
                (c.cylinder().radius(), scene.straight_edges(f))
            }
            // A circular-edge blend: a torus band, tangent across its rings.
            SurfaceGeometry::Torus(t) => {
                let rings: Vec<Shape> = scene.faces[f]
                    .wires
                    .iter()
                    .flatten()
                    .filter(|e| !scene.is_seam_in(f, e))
                    .cloned()
                    .collect();
                (t.torus().minor_radius(), rings)
            }
            _ => continue,
        };
        if tangent_edges.len() < 2 {
            continue;
        }
        let mut smooth = 0;
        let mut fold_probe: Option<Fold> = None;
        for edge in &tangent_edges {
            let Some(other) = scene.other_face(edge, f) else {
                continue;
            };
            match scene.edge_fold(edge, f, other) {
                Some(Fold::Smooth) => smooth += 1,
                fold => fold_probe = fold.or(fold_probe),
            }
        }
        if smooth < 2 {
            continue;
        }
        // Concavity of the blend itself: nudge off the blend's own middle
        // along its outward normal's *reverse* — inside means the blend
        // bulges outward (a convex round), outside means it hollows an
        // inside corner (a concave fillet).
        let _ = fold_probe;
        let Some(sample) = scene.faces[f]
            .wires
            .iter()
            .flatten()
            .find_map(|e| scene.edge_midpoint(e))
        else {
            continue;
        };
        let Some(normal) = scene.outward_normal(f, sample) else {
            continue;
        };
        let concave = {
            let axis_side = match &scene.faces[f].surface {
                SurfaceGeometry::Cylinder(c) => {
                    let axis = c.cylinder().frame().axis();
                    let foot = axis.location
                        + axis.direction.vector()
                            * (sample - axis.location).dot(axis.direction.vector());
                    (sample - foot).dot(normal)
                }
                SurfaceGeometry::Torus(t) => {
                    let frame = t.torus().frame();
                    let local = frame.to_local(sample);
                    let radial = ogeom_math::Vector::new(local.x, local.y, 0.0);
                    let spine_local = if radial.magnitude() > scene.tol.confusion() {
                        radial / radial.magnitude() * t.torus().major_radius()
                    } else {
                        radial
                    };
                    let spine = frame.origin()
                        + frame.x().vector() * spine_local.x
                        + frame.y().vector() * spine_local.y;
                    (sample - spine).dot(normal)
                }
                _ => 1.0,
            };
            axis_side < 0.0
        };
        claimed.push(scene.faces[f].shape.node());
        features.push(Feature::Fillet(Fillet {
            face: scene.faces[f].shape.clone(),
            radius,
            concave,
        }));
    }
    Ok(())
}

// --- chamfers ----------------------------------------------------------------

fn chamfers(
    scene: &Scene,
    features: &mut Vec<Feature>,
    claimed: &mut Vec<TShapeId>,
) -> OgeomResult<()> {
    for f in 0..scene.faces.len() {
        if claimed.contains(&scene.faces[f].shape.node()) {
            continue;
        }
        if !matches!(scene.faces[f].surface, SurfaceGeometry::Plane(_)) {
            continue;
        }
        let info = &scene.faces[f];
        if info.wires.len() != 1 || info.wires[0].len() != 4 {
            continue;
        }
        // A bevel meets its neighbours at a deliberate slant: the normals
        // disagree by something strictly between along and across.
        let mut slanted = 0;
        let mut smooth = 0;
        for edge in &info.wires[0] {
            let Some(other) = scene.other_face(edge, f) else {
                continue;
            };
            let Some(p) = scene.edge_midpoint(edge) else {
                continue;
            };
            let (Some(n1), Some(n2)) = (scene.outward_normal(f, p), scene.outward_normal(other, p))
            else {
                continue;
            };
            let angle = n1.dot(n2).clamp(-1.0, 1.0).acos();
            if angle <= TANGENT_ANGLE {
                smooth += 1;
            } else if (0.26..=1.31).contains(&angle) {
                // 15 to 75 degrees.
                slanted += 1;
            }
        }
        if slanted >= 2 && smooth == 0 {
            claimed.push(info.shape.node());
            features.push(Feature::Chamfer(Chamfer {
                face: info.shape.clone(),
            }));
        }
    }
    Ok(())
}

// --- pockets -----------------------------------------------------------------

fn pockets(
    scene: &Scene,
    features: &mut Vec<Feature>,
    claimed: &mut Vec<TShapeId>,
) -> OgeomResult<()> {
    for f in 0..scene.faces.len() {
        if claimed.contains(&scene.faces[f].shape.node()) {
            continue;
        }
        if !matches!(scene.faces[f].surface, SurfaceGeometry::Plane(_)) {
            continue;
        }
        let info = &scene.faces[f];
        let Some(outer) = info.wires.first() else {
            continue;
        };
        if outer.is_empty() {
            continue;
        }
        let mut walls = Vec::new();
        let mut all_concave = true;
        for edge in outer {
            let Some(other) = scene.other_face(edge, f) else {
                all_concave = false;
                break;
            };
            match scene.edge_fold(edge, f, other) {
                Some(Fold::Concave) => {
                    let shape = scene.faces[other].shape.clone();
                    if !walls.iter().any(|held: &Shape| held.is_same(&shape)) {
                        walls.push(shape);
                    }
                }
                _ => {
                    all_concave = false;
                    break;
                }
            }
        }
        if !all_concave || walls.is_empty() {
            continue;
        }
        // An obround floor — two parallel lines closed by two arcs — is
        // what a slot leaves behind.
        let slot = {
            let mut lines: Vec<Vector> = Vec::new();
            let mut arcs = 0;
            for edge in outer {
                let curve = scene
                    .model
                    .node(edge)
                    .and_then(|n| n.data().as_edge())
                    .and_then(|d| d.curve3d())
                    .and_then(|r| match r {
                        EdgeRepr::Curve3d { curve, .. } => scene.model.geometry().curve(*curve),
                        _ => None,
                    });
                match curve {
                    Some(Curve::Line(l)) => lines.push(l.axis().direction.vector()),
                    Some(Curve::Circle(_)) => arcs += 1,
                    _ => {}
                }
            }
            lines.len() == 2 && arcs == 2 && lines[0].cross(lines[1]).magnitude() < 1e-6
        };
        claimed.push(info.shape.node());
        features.push(Feature::Pocket(Pocket {
            floor: info.shape.clone(),
            walls,
            slot,
        }));
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_math::Frame;

    const T: Tolerances = Tolerances::millimetres();

    /// The box edge whose midpoint is nearest `at`.
    fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
        use ogeom_geom::Curve3d as _;
        let mut best: Option<(f64, Shape)> = None;
        for edge in ogeom_topo::explore_unique(model, solid, ShapeType::Edge).unwrap() {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                continue;
            };
            let p = geometry
                .point_at(f64::midpoint(range.0, range.1), T)
                .unwrap();
            let d = p.distance(at);
            if best.as_ref().is_none_or(|(held, _)| d < *held) {
                best = Some((d, edge));
            }
        }
        best.unwrap().1
    }

    fn drilled(depth_through: bool) -> (Model, Shape) {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        // A drill down the middle: through, or stopping short.
        let (z0, h) = if depth_through {
            (-1.0, 12.0)
        } else {
            (4.0, 7.0)
        };
        let frame = Frame::new(
            Point::new(10.0, 10.0, z0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let drill = ogeom_algo::make_cylinder(&mut model, frame, 3.0, h, T)
            .unwrap()
            .shape;
        let cut = ogeom_bool::cut(&mut model, &block, &drill, T)
            .unwrap()
            .shape;
        (model, cut)
    }

    #[test]
    fn a_plain_box_has_no_features() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let features = recognize(&model, &block, T).unwrap();
        assert!(features.is_empty(), "found {features:?}");
    }

    #[test]
    fn a_through_hole_is_a_through_hole() {
        let (model, cut) = drilled(true);
        let features = recognize(&model, &cut, T).unwrap();
        let holes: Vec<&Hole> = features
            .iter()
            .filter_map(|f| match f {
                Feature::Hole(h) => Some(h),
                _ => None,
            })
            .collect();
        assert_eq!(holes.len(), 1, "features: {features:?}");
        assert_eq!(holes[0].kind, HoleKind::Through);
        assert_relative_eq!(holes[0].radius, 3.0, epsilon = 1e-9);
        assert!(!holes[0].counterbored);
        assert!(!holes[0].countersunk);
        assert_eq!(features.len(), 1, "nothing else: {features:?}");
    }

    #[test]
    fn a_stopped_drill_is_a_blind_hole() {
        let (model, cut) = drilled(false);
        let features = recognize(&model, &cut, T).unwrap();
        let holes: Vec<&Hole> = features
            .iter()
            .filter_map(|f| match f {
                Feature::Hole(h) => Some(h),
                _ => None,
            })
            .collect();
        assert_eq!(holes.len(), 1, "features: {features:?}");
        assert_eq!(holes[0].kind, HoleKind::Blind);
    }

    #[test]
    fn a_counterbore_reads_as_one_stepped_hole() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let narrow = ogeom_algo::make_cylinder(
            &mut model,
            Frame::new(
                Point::new(10.0, 10.0, -1.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            2.0,
            12.0,
            T,
        )
        .unwrap()
        .shape;
        let wide = ogeom_algo::make_cylinder(
            &mut model,
            Frame::new(
                Point::new(10.0, 10.0, 6.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            4.0,
            5.0,
            T,
        )
        .unwrap()
        .shape;
        let cut = ogeom_bool::cut(&mut model, &block, &narrow, T)
            .unwrap()
            .shape;
        let cut = ogeom_bool::cut(&mut model, &cut, &wide, T).unwrap().shape;
        let features = recognize(&model, &cut, T).unwrap();
        let holes: Vec<&Hole> = features
            .iter()
            .filter_map(|f| match f {
                Feature::Hole(h) => Some(h),
                _ => None,
            })
            .collect();
        assert_eq!(holes.len(), 1, "one chained hole: {features:?}");
        assert!(holes[0].counterbored);
        assert_relative_eq!(holes[0].radius, 2.0, epsilon = 1e-9);
    }

    #[test]
    fn a_rounded_edge_is_a_convex_fillet() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let edge = edge_near(&model, &block, Point::new(5.0, 0.0, 10.0));
        let rounded = ogeom_fillet::fillet_edge(&mut model, &block, &edge, 2.0, T)
            .unwrap()
            .shape;
        let features = recognize(&model, &rounded, T).unwrap();
        let fillets: Vec<&Fillet> = features
            .iter()
            .filter_map(|f| match f {
                Feature::Fillet(x) => Some(x),
                _ => None,
            })
            .collect();
        assert_eq!(fillets.len(), 1, "features: {features:?}");
        assert_relative_eq!(fillets[0].radius, 2.0, epsilon = 1e-9);
        assert!(!fillets[0].concave, "an outside round is convex");
    }

    #[test]
    fn a_bevelled_edge_is_a_chamfer() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let edge = edge_near(&model, &block, Point::new(5.0, 0.0, 10.0));
        let bevelled = ogeom_fillet::chamfer_edge(&mut model, &block, &edge, 1.5, T)
            .unwrap()
            .shape;
        let features = recognize(&model, &bevelled, T).unwrap();
        assert!(
            features.iter().any(|f| matches!(f, Feature::Chamfer(_))),
            "features: {features:?}"
        );
    }

    #[test]
    fn a_milled_recess_is_a_pocket_with_its_walls() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let mill = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(6.0, 6.0, 6.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (8.0, 8.0, 6.0),
            T,
        )
        .unwrap()
        .shape;
        let cut = ogeom_bool::cut(&mut model, &block, &mill, T).unwrap().shape;
        let features = recognize(&model, &cut, T).unwrap();
        let pockets: Vec<&Pocket> = features
            .iter()
            .filter_map(|f| match f {
                Feature::Pocket(p) => Some(p),
                _ => None,
            })
            .collect();
        assert_eq!(pockets.len(), 1, "features: {features:?}");
        assert_eq!(pockets[0].walls.len(), 4);
        assert!(!pockets[0].slot);
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod slot_tests {
    use super::*;
    const T: Tolerances = Tolerances::millimetres();

    /// An obround mill — a prism over two parallel lines closed by two
    /// arcs — leaves a slot: a pocket whose floor keeps that outline.
    #[test]
    fn a_slot_mill_leaves_a_slot() {
        use ogeom_geom::{CircleCurve, LineCurve};
        use ogeom_math::{Circle, Direction, Frame};

        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (24.0, 20.0, 10.0), T)
            .unwrap()
            .shape;

        // The obround profile at z = 6: line, arc, line, arc.
        let a = ogeom_algo::make_vertex(&mut model, Point::new(8.0, 8.0, 6.0)).shape;
        let b = ogeom_algo::make_vertex(&mut model, Point::new(16.0, 8.0, 6.0)).shape;
        let c = ogeom_algo::make_vertex(&mut model, Point::new(16.0, 12.0, 6.0)).shape;
        let d = ogeom_algo::make_vertex(&mut model, Point::new(8.0, 12.0, 6.0)).shape;
        let line = |model: &mut Model, from: Point, to: Point, va: &Shape, vb: &Shape| {
            let curve: ogeom_geom::Curve = LineCurve::segment(from, to, T).unwrap().into();
            ogeom_algo::make_edge_between(model, curve, (0.0, from.distance(to)), va, vb, T)
                .unwrap()
                .shape
        };
        let arc = |model: &mut Model, centre: Point, from: &Shape, to: &Shape| {
            // A half circle from angle -pi/2 to +pi/2 in the frame whose x
            // points +y... simplest: frame x toward the `from` vertex.
            let start = model
                .node(from)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .unwrap();
            let x = Direction::new(start - centre, T).unwrap();
            let frame = Frame::new(centre, Direction::Z, x, T).unwrap();
            let circle = Circle::new(frame, 2.0, T).unwrap();
            let curve: ogeom_geom::Curve = CircleCurve::new(circle).into();
            ogeom_algo::make_edge_between(model, curve, (0.0, core::f64::consts::PI), from, to, T)
                .unwrap()
                .shape
        };
        let bottom = line(
            &mut model,
            Point::new(8.0, 8.0, 6.0),
            Point::new(16.0, 8.0, 6.0),
            &a,
            &b,
        );
        let right_arc = arc(&mut model, Point::new(16.0, 10.0, 6.0), &b, &c);
        let top = line(
            &mut model,
            Point::new(16.0, 12.0, 6.0),
            Point::new(8.0, 12.0, 6.0),
            &c,
            &d,
        );
        let left_arc = arc(&mut model, Point::new(8.0, 10.0, 6.0), &d, &a);
        let edges = vec![bottom, right_arc, top, left_arc];
        let probe = ogeom_algo::make_wire(&mut model, &edges, T).unwrap().shape;
        let plane = ogeom_algo::find_plane(&model, &probe, T).unwrap().unwrap();
        let surface: SurfaceGeometry =
            ogeom_geom::PlaneSurface::over(plane, (-30.0, 30.0), (-30.0, 30.0))
                .unwrap()
                .into();
        let profile = ogeom_algo::make_face_with_pcurves(&mut model, surface, &[edges], T)
            .unwrap()
            .shape;
        let mill = ogeom_algo::make_prism(&mut model, &profile, Vector::new(0.0, 0.0, 6.0), T)
            .unwrap()
            .shape;
        let cut = ogeom_bool::cut(&mut model, &block, &mill, T).unwrap().shape;

        let features = recognize(&model, &cut, T).unwrap();
        let pockets: Vec<&Pocket> = features
            .iter()
            .filter_map(|f| match f {
                Feature::Pocket(p) => Some(p),
                _ => None,
            })
            .collect();
        assert_eq!(pockets.len(), 1, "features: {features:?}");
        assert!(pockets[0].slot, "an obround pocket is a slot");
    }
}
