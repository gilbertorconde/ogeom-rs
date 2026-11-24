//! Sweeping a profile of several faces that share edges: a sketch's adjacent
//! regions.
//!
//! Each face sweeps into a solid of its own, and the solids glue into one
//! compsolid along the walls their shared edges sweep. A shared edge sweeps
//! one wall, bounding the solid on each side of it, once from each side;
//! faces touching nowhere, or only at a vertex, stay separate solids in a
//! compound. The outer boundary of the result is the single solid the
//! regions make together ([`outer_solid`](crate::outer_solid)).
//!
//! The pieces may not overlap: overlapping regions would sweep solids
//! through each other, so the profile is refused before anything is built
//! if its faces are not in one plane, if two of their boundaries cross or
//! run along each other without sharing an edge, if one face holds part of
//! another, or (for a revolution) if they lie on both sides of the axis.

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d, Surface as _, SurfaceKind};
use ogeom_intersect::{CurveCurveOptions, intersect_curves};
use ogeom_math::{Axis, Point, Vector};
use ogeom_topo::{Filter, Model, SameKey, Shape, ShapeType, explore};

use crate::build::{edge_vertices, make_compound, make_compsolid, placed_edge_curve};
use crate::history::{Built, History};

/// The faces a profile of several faces holds, in order: the faces of a
/// shell, or those of a compound of faces and shells.
///
/// # Errors
///
/// Construction if the profile holds anything but faces, holds no face, or
/// holds one face twice.
pub(super) fn profile_faces(model: &Model, profile: &Shape) -> OgeomResult<Vec<Shape>> {
    let mut faces: Vec<Shape> = Vec::new();
    let mut stack = vec![profile.clone()];
    while let Some(shape) = stack.pop() {
        match model.kind_of(&shape)? {
            ShapeType::Face => {
                if faces.iter().any(|f| f.is_same(&shape)) {
                    ogeom_bail!(
                        Construction,
                        "the profile holds one face twice; each region of a \
                         profile of several faces is swept once"
                    );
                }
                faces.push(shape);
            }
            ShapeType::Shell | ShapeType::Compound => {
                stack.extend(model.children_of(&shape)?.into_iter().rev());
            }
            other => ogeom_bail!(
                Construction,
                "a profile of several faces holds faces only; this one holds a \
                 {other:?}"
            ),
        }
    }
    if faces.is_empty() {
        ogeom_bail!(Construction, "the profile holds no face to sweep");
    }
    Ok(faces)
}

/// The plane every face lies in, as a point on it and its unit normal.
///
/// # Errors
///
/// Construction if a face is not planar or the faces are not in one plane.
pub(super) fn one_plane(
    model: &Model,
    faces: &[Shape],
    tol: Tolerances,
) -> OgeomResult<(Point, Vector)> {
    let (origin, normal) = crate::measure::face_normal(model, &faces[0], tol)?;
    let normal = normal / normal.magnitude();
    for face in faces {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            ogeom_bail!(Dangling, "a face of the profile is not in this model");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "a face's surface is not in this model");
        };
        if surface.kind() != SurfaceKind::Plane {
            ogeom_bail!(
                Construction,
                "a profile of several faces sweeps planar faces; one face here \
                 is a {:?}",
                surface.kind()
            );
        }
        let (_, own) = crate::measure::face_normal(model, face, tol)?;
        if own.cross(normal).magnitude() > own.magnitude() * tol.angular() {
            ogeom_bail!(
                Construction,
                "the faces of the profile are not in one plane; regions swept \
                 from different planes can pass through each other"
            );
        }
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            for point in edge_samples(model, &edge, tol)? {
                if (point - origin).dot(normal).abs() > tol.confusion() * 10.0 {
                    ogeom_bail!(
                        Construction,
                        "the faces of the profile are not in one plane; regions \
                         swept from different planes can pass through each other"
                    );
                }
            }
        }
    }
    Ok((origin, normal))
}

/// Refuse faces whose regions overlap.
///
/// Two planar regions overlap only if their boundaries cross, run along each
/// other, or one region holds a point of the other's boundary. Edges both
/// faces share are their common boundary and are skipped; a contact at a
/// vertex both share is where the boundaries meet, not a crossing.
///
/// # Errors
///
/// Construction naming how the faces overlap.
pub(super) fn apart(model: &Model, faces: &[Shape], tol: Tolerances) -> OgeomResult<()> {
    struct Bound {
        edge: Shape,
        curve: ogeom_geom::Curve,
        ends: Vec<(Shape, Point)>,
    }
    let near = tol.confusion() * 10.0;
    // Faces whose boxes miss each other cannot overlap.
    let boxes: Vec<ogeom_math::Aabb> = faces
        .iter()
        .map(|f| Ok(crate::face_bounds(model, f)?.expanded(near)))
        .collect::<OgeomResult<_>>()?;
    let mut bounds: Vec<Vec<Bound>> = Vec::with_capacity(faces.len());
    for face in faces {
        let mut own = Vec::new();
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            if own.iter().any(|b: &Bound| b.edge.is_same(&edge)) {
                continue;
            }
            let Some((curve, (lo, hi))) = placed_edge_curve(model, &edge, tol)? else {
                ogeom_bail!(Construction, "an edge of the profile has no curve");
            };
            let curve: ogeom_geom::Curve =
                ogeom_geom::TrimmedCurve::new(curve, lo, hi, tol)?.into();
            let mut ends = Vec::new();
            if let Some((start, end)) = edge_vertices(model, &edge)? {
                for vertex in [start, end] {
                    let point = vertex_point(model, &vertex)?;
                    ends.push((vertex, point));
                }
            }
            own.push(Bound { edge, curve, ends });
        }
        bounds.push(own);
    }

    let options = CurveCurveOptions {
        gap: tol.confusion(),
        ..CurveCurveOptions::default()
    };
    for (i, first) in bounds.iter().enumerate() {
        for (j, second) in bounds.iter().enumerate().skip(i + 1) {
            if !boxes[i].intersects(&boxes[j]) {
                continue;
            }
            for a in first {
                for b in second {
                    if a.edge.is_same(&b.edge) {
                        continue;
                    }
                    let met = intersect_curves(&a.curve, &b.curve, options, tol)?;
                    if !met.overlaps.is_empty() {
                        ogeom_bail!(
                            Construction,
                            "two faces of the profile have edges running along \
                             each other without sharing them; split the longer \
                             edge where the other ends so the faces share it"
                        );
                    }
                    for crossing in &met.crossings {
                        let shared = a.ends.iter().any(|(va, pa)| {
                            (crossing.point - *pa).magnitude() <= near
                                && b.ends.iter().any(|(vb, _)| vb.is_same(va))
                        });
                        if !shared {
                            ogeom_bail!(
                                Construction,
                                "the boundaries of two faces of the profile meet \
                                 at {:?} where they share no vertex: the faces \
                                 overlap, or one edge ends on another that is \
                                 not split there",
                                crossing.point
                            );
                        }
                    }
                }
            }
        }
    }

    // The boundaries meet only where they share topology, so each face
    // either holds the other's unshared edges whole or holds none of them.
    let deflection = ogeom_mesh::Deflection {
        chord: tol.confusion() * 100.0,
        ..ogeom_mesh::Deflection::default()
    };
    for (i, face) in faces.iter().enumerate() {
        for (j, other) in bounds.iter().enumerate() {
            if i == j || !boxes[i].intersects(&boxes[j]) {
                continue;
            }
            for b in other {
                if bounds[i].iter().any(|a| a.edge.is_same(&b.edge)) {
                    continue;
                }
                let (lo, hi) = b.curve.domain();
                let mid = b.curve.point_at(f64::midpoint(lo, hi), tol)?;
                if crate::classify_on_face(model, face, mid, deflection, tol)?
                    != crate::Containment::Out
                {
                    ogeom_bail!(
                        Construction,
                        "one face of the profile holds part of another: the \
                         regions overlap, and their sweeps would pass through \
                         each other"
                    );
                }
            }
        }
    }
    Ok(())
}

/// Refuse a profile of several faces that a revolution about `axis` would
/// sweep through itself: the faces must lie in one plane through the axis,
/// all on one side of it.
///
/// # Errors
///
/// Construction if the axis is not in the profile's plane or the faces lie
/// on both sides of it.
pub(super) fn beside_axis(
    model: &Model,
    faces: &[Shape],
    (origin, normal): (Point, Vector),
    axis: Axis,
    tol: Tolerances,
) -> OgeomResult<()> {
    let along = axis.direction.vector();
    if along.dot(normal).abs() > tol.angular()
        || (axis.location - origin).dot(normal).abs() > tol.confusion() * 10.0
    {
        ogeom_bail!(
            Construction,
            "a profile of several faces revolves about an axis in its own plane; \
             off that plane the regions' turns can pass through each other"
        );
    }
    // Across the axis within the plane: which side a point is on.
    let across = along.cross(normal);
    let mut sides = (false, false);
    for face in faces {
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            for point in edge_samples(model, &edge, tol)? {
                let side = (point - axis.location).dot(across);
                if side > tol.confusion() * 10.0 {
                    sides.0 = true;
                } else if side < -tol.confusion() * 10.0 {
                    sides.1 = true;
                }
            }
        }
    }
    if sides.0 && sides.1 {
        ogeom_bail!(
            Construction,
            "the faces of the profile lie on both sides of the axis, so their \
             turns pass through each other"
        );
    }
    Ok(())
}

/// The walls built so far in one sweep of a profile of several faces, keyed
/// by the profile edge each rose from, so an edge two faces share sweeps one
/// wall that both solids hold.
#[derive(Default)]
pub(super) struct Walls {
    /// The face of the profile being swept.
    pub(super) region: usize,
    /// Each edge's wall: the region that built it, the edge as that region
    /// walked it, and the wall as that region's solid holds it.
    built: FastMap<SameKey, (usize, Shape, Option<Shape>)>,
    /// The edges two regions have shared.
    shared: ogeom_core::FastSet<SameKey>,
    /// The pairs of regions that share a wall.
    links: Vec<(usize, usize)>,
}

impl Walls {
    /// The wall `edge` sweeps, built by `build` unless another face of the
    /// profile already swept it; then that wall, as the solid on this side
    /// holds it, and no history: the first sweep recorded it.
    pub(super) fn sweep(
        &mut self,
        model: &mut Model,
        edge: &Shape,
        build: impl FnOnce(&mut Model) -> OgeomResult<(Option<Shape>, History)>,
    ) -> OgeomResult<(Option<Shape>, History)> {
        let key = SameKey(edge.clone());
        if let Some((region, first, wall)) = self.built.get(&key)
            && *region != self.region
        {
            if first.orientation() == edge.orientation() {
                ogeom_bail!(
                    Construction,
                    "two faces of the profile walk a shared edge the same way, \
                     so they lie on one side of it and overlap"
                );
            }
            if !self.shared.insert(key) {
                ogeom_bail!(
                    Construction,
                    "an edge of the profile bounds more than two of its faces"
                );
            }
            self.links.push((*region, self.region));
            return Ok((wall.as_ref().map(Shape::reversed), History::new()));
        }
        let (wall, history) = build(model)?;
        self.built
            .entry(key)
            .or_insert_with(|| (self.region, edge.clone(), wall.clone()));
        Ok((wall, history))
    }
}

/// Sweep each face of `faces` with `one`, the faces' shared walls through
/// one [`Walls`], and glue the solids: a compsolid for each group of solids
/// joined through shared walls, a lone solid for a face sharing no wall, and
/// a compound of those where there is more than one.
///
/// The history joins each face's and adds the profile generating the result.
pub(super) fn swept_regions(
    model: &mut Model,
    profile: &Shape,
    faces: &[Shape],
    mut one: impl FnMut(&mut Model, &mut Walls, &Shape) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    let mut walls = Walls::default();
    let mut history = History::new();
    let mut solids = Vec::with_capacity(faces.len());
    for (region, face) in faces.iter().enumerate() {
        walls.region = region;
        let built = one(model, &mut walls, face)?;
        history = history.joined(built.history);
        solids.push(built.shape);
    }

    // Groups of regions joined through shared walls, by union-find.
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    fn root(parent: &mut [usize], mut at: usize) -> usize {
        while parent[at] != at {
            parent[at] = parent[parent[at]];
            at = parent[at];
        }
        at
    }
    for &(a, b) in &walls.links {
        let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
        parent[ra.max(rb)] = ra.min(rb);
    }
    let mut groups: Vec<(usize, Vec<Shape>)> = Vec::new();
    for (region, solid) in solids.into_iter().enumerate() {
        let r = root(&mut parent, region);
        match groups.iter_mut().find(|(g, _)| *g == r) {
            Some((_, members)) => members.push(solid),
            None => groups.push((r, vec![solid])),
        }
    }
    let mut pieces = Vec::with_capacity(groups.len());
    for (_, members) in groups {
        if let [solid] = members.as_slice() {
            pieces.push(solid.clone());
        } else {
            let glued = make_compsolid(model, &members)?;
            history = history.joined(glued.history);
            pieces.push(glued.shape);
        }
    }
    let shape = if let [piece] = pieces.as_slice() {
        piece.clone()
    } else {
        let gathered = make_compound(model, &pieces)?;
        history = history.joined(gathered.history);
        gathered.shape
    };
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// Where a vertex occurrence stands.
fn vertex_point(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Dangling, "a vertex of the profile is not in this model");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// A few points along an edge where it stands: its ends and three between.
fn edge_samples(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let Some((curve, (lo, hi))) = placed_edge_curve(model, edge, tol)? else {
        return Ok(Vec::new());
    };
    (0..=4)
        .map(|k| curve.point_at(lo + (hi - lo) * f64::from(k) / 4.0, tol))
        .collect()
}
