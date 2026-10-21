//! A solid whose boundary passes through itself, resolved into the union
//! of the material it encloses.
//!
//! A boundary that crosses itself encloses some regions twice: the winding
//! number, how many more times a ray from a point leaves through the
//! boundary than enters, is two there, one in the material enclosed once
//! and zero outside. The material is where it is one or more, and its
//! boundary is made of the patches of faces with none on their outer side
//! and some on their inner side.
//!
//! So every face is split by every other face it meets (the general fuse
//! of the solid with itself, each face against all but itself and the
//! faces it shares an edge or a corner with), and each piece is read on
//! its two sides, a step off it along its normal that meets no other face
//! on the way. A piece with winding zero beyond it and one or more behind
//! it bounds the union and is kept; every other piece lies inside the
//! union or outside it, and goes. Where faces lie on one surface over a
//! patch, facing the same way, the patch is kept once, from the first of
//! them. The kept pieces are sewn and closed as a boolean's are.

use ogeom_algo::{Built, History, SolidBoundary};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Surface as _;
use ogeom_math::{Point, Point2, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::{
    FacePiece, GFace, GeneralFused, Operands, assemble_result, baked_if_scaled, general_fuse_as,
    outward_normal,
};

/// The union of the material `solid` encloses, where its boundary passes
/// through itself: the regions it encloses twice made one, and its faces
/// cut where they cross and trimmed back to the union's boundary.
///
/// A solid whose boundary does not pass through itself (what
/// [`self_crossings`](ogeom_algo::self_crossings) finds) comes back as it
/// is, with an empty history. A region the boundary winds round negatively,
/// enclosed inside out, is not material and is left out.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `solid` is not a solid;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a face
/// passes through itself or through a face it shares an edge or a corner
/// with (only faces apart are split against each other), where a piece's
/// two sides cannot be read, or where the kept pieces do not close or the
/// union still crosses itself; and as the boolean's operations.
pub fn self_union(model: &mut Model, solid: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    if model.kind_of(solid)? != ShapeType::Solid {
        ogeom_bail!(Construction, "a self-union resolves a solid");
    }
    let crossings = ogeom_algo::self_crossings(model, solid, tol)?;
    if crossings.is_empty() {
        return Ok(Built::new(solid.clone(), History::new()));
    }
    // Only faces apart are split against each other: a face meets those it
    // shares an edge or a corner with there already.
    let corners = |face: &Shape| -> OgeomResult<ogeom_core::FastSet<super::EdgeKey>> {
        Ok(explore_unique(model, face, ShapeType::Vertex)?
            .iter()
            .map(super::EdgeKey::of)
            .collect())
    };
    for crossing in &crossings {
        if crossing.face.is_same(&crossing.other)
            || !corners(&crossing.face)?.is_disjoint(&corners(&crossing.other)?)
        {
            ogeom_bail!(
                NotDone,
                "a face passes through itself or through a face it shares an edge or a \
                 corner with, which the self-union does not split"
            );
        }
    }
    let solid = &baked_if_scaled(model, solid, tol)?;
    let fused = general_fuse_as(model, solid, solid, Operands::Itself, None, tol)?;
    let boundary = SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
    let mut keep: Vec<Option<bool>> = vec![None; fused.pieces.len()];
    for (index, piece) in fused.pieces.iter().enumerate() {
        if piece.from_a && !piece.outlines.is_empty() {
            let Some(bounds) = bounds_the_union(model, &fused, piece, &boundary, tol)? else {
                ogeom_bail!(
                    NotDone,
                    "a piece of a face reads its two sides nowhere a step off it meets no \
                     other face"
                );
            };
            keep[index] = Some(bounds);
        }
    }
    for group in quiet_groups(&fused) {
        let mut read = None;
        for &index in &group {
            let Some(this) = bounds_the_union(model, &fused, &fused.pieces[index], &boundary, tol)?
            else {
                continue;
            };
            if read.is_some_and(|read| read != this) {
                ogeom_bail!(
                    NotDone,
                    "a run of faces nothing crosses bounds the union in part, which only a \
                     face meeting its edges makes"
                );
            }
            read = Some(this);
        }
        let Some(bounds) = read else {
            ogeom_bail!(
                NotDone,
                "a run of faces nothing crosses reads its two sides nowhere a step off it \
                 meets no other face"
            );
        };
        for index in group {
            keep[index] = Some(bounds);
        }
    }
    let kept: Vec<(usize, bool)> = keep
        .iter()
        .enumerate()
        .filter(|(_, keep)| **keep == Some(true))
        .map(|(index, _)| (index, false))
        .collect();
    let built = assemble_result(model, &fused, &kept, [false, false], solid, solid, tol)?;
    for made in explore_unique(model, &built.shape, ShapeType::Solid)? {
        if !ogeom_algo::self_crossings(model, &made, tol)?.is_empty() {
            ogeom_bail!(NotDone, "the union still passes through itself");
        }
    }
    Ok(built)
}

/// Whether `piece` bounds the union: winding zero just beyond it, one or
/// more just behind it, and, where faces lie on its surface there facing
/// the same way, it is the first of them. `None` where the piece has no
/// interior point drawn, no normal there, or no step off it reads.
fn bounds_the_union(
    model: &Model,
    fused: &GeneralFused,
    piece: &FacePiece,
    boundary: &SolidBoundary,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    let face = &fused.a.faces[piece.face];
    let Some((at, uv)) = interior_point(model, face, piece, tol)? else {
        return Ok(None);
    };
    let normal = outward_normal(face, uv, tol)?;
    let length = normal.magnitude();
    if length <= tol.confusion() {
        return Ok(None);
    }
    let normal = normal / length;
    // The faces lying where the piece does, itself among them, and how
    // many of them face the way it does net of those facing against it:
    // what stepping through them all changes the winding by.
    let on = boundary.faces_at(model, at, tol)?;
    let mut net = 0_i64;
    let mut first = true;
    for other in &on {
        let Some(gathered) = fused.a.faces.iter().position(|g| g.face.is_same(other)) else {
            ogeom_bail!(NotDone, "a face of the solid was not gathered");
        };
        let theirs = &fused.a.faces[gathered];
        let foot = ogeom_algo::project_on_surface(&theirs.surface, at, 16, tol)?;
        let (u, v) = foot.parameters;
        let along = outward_normal(theirs, Point2::new(u, v), tol)?.dot(normal);
        if along > 0.0 {
            net += 1;
            if gathered < piece.face {
                first = false;
            }
        } else {
            net -= 1;
        }
    }
    let Some((beyond, behind)) = sides(model, face, at, normal, &on, boundary, tol)? else {
        return Ok(None);
    };
    if behind - beyond != net {
        ogeom_bail!(
            NotDone,
            "a piece of a face reads windings {beyond} and {behind} on its two sides, \
             which the {net} faces lying there do not step between"
        );
    }
    if beyond >= 1 && behind <= 0 {
        ogeom_bail!(
            NotDone,
            "a piece of a face bounds the union facing into it, which no face lying there \
             is turned to"
        );
    }
    Ok(Some(beyond == 0 && behind >= 1 && first))
}

/// The winding numbers just beyond `at` along `normal` and just behind it,
/// read a step off each way along a segment that meets no face past those
/// it starts on, `on`: from a twentieth of the face's size down by tenths to twenty
/// times the band a point counts as on it within. `None` where no step
/// reads.
fn sides(
    model: &Model,
    face: &GFace,
    at: Point,
    normal: Vector,
    on: &[Shape],
    boundary: &SolidBoundary,
    tol: Tolerances,
) -> OgeomResult<Option<(i64, i64)>> {
    let shortest = face.tolerance.max(tol.confusion() * 1e3) * 20.0;
    let mut step = (face.filter.diagonal() * 0.05).max(shortest);
    loop {
        ogeom_core::progress::checkpoint()?;
        let (out_at, in_at) = (at + normal * step, at - normal * step);
        if boundary.clear_from(model, at, out_at, on, tol)?
            && boundary.clear_from(model, at, in_at, on, tol)?
            && let (Some(beyond), Some(behind)) = (
                boundary.winding(model, out_at, tol)?,
                boundary.winding(model, in_at, tol)?,
            )
        {
            return Ok(Some((beyond, behind)));
        }
        if step <= shortest {
            return Ok(None);
        }
        step = (step * 0.1).max(shortest);
    }
}

/// A point inside `piece`, in space and in its face's chart: the probe the
/// arrangement chose for a piece it split off, and for a face nothing split
/// the chart middle of a triangle of its mesh, the largest first, that lies
/// inside the face itself, clear of its trim. `None` where none does.
fn interior_point(
    model: &Model,
    face: &GFace,
    piece: &FacePiece,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Point2)>> {
    if !piece.outlines.is_empty() {
        let foot = ogeom_algo::project_on_surface(&face.surface, piece.probe, 16, tol)?;
        let (u, v) = foot.parameters;
        return Ok(Some((piece.probe, Point2::new(u, v))));
    }
    let mesh = ogeom_mesh::triangulate_face(model, &face.face, Deflection::default(), tol)?;
    if mesh.parameters.len() != mesh.positions.len() {
        return Ok(None);
    }
    let area = |t: &[u32; 3]| {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        (b - a).cross(c - a).magnitude()
    };
    let mut triangles: Vec<&[u32; 3]> = mesh.triangles.iter().collect();
    triangles.sort_by(|x, y| area(y).total_cmp(&area(x)));
    for triangle in triangles {
        let corners = triangle.map(|i| mesh.parameters[i as usize]);
        let uv = Point2::new(
            (corners[0].0 + corners[1].0 + corners[2].0) / 3.0,
            (corners[0].1 + corners[1].1 + corners[2].1) / 3.0,
        );
        let at = face.surface.point_at(uv.x, uv.y, tol)?;
        if ogeom_algo::classify_on_face(model, &face.face, at, Deflection::default(), tol)?
            == ogeom_algo::Containment::In
        {
            return Ok(Some((at, uv)));
        }
    }
    Ok(None)
}

/// The pieces that are whole faces nothing crosses, in runs joined across
/// the edges only two faces of the solid meet at, each run in face order.
///
/// Nothing crosses such a face or comes near its edges, so beside one of
/// its edges the faces either side of it bound the same two regions:
/// a whole run reads alike. Each face of it is read where it can be, and
/// the readings must agree.
fn quiet_groups(fused: &GeneralFused) -> Vec<Vec<usize>> {
    let quiet: Vec<usize> = fused
        .pieces
        .iter()
        .enumerate()
        .filter(|(_, piece)| piece.from_a && piece.outlines.is_empty())
        .map(|(index, _)| index)
        .collect();
    let mut parent: Vec<usize> = (0..quiet.len()).collect();
    let find = |parent: &mut Vec<usize>, mut i: usize| {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    };
    // An edge more faces than two meet at joins none of them: a face beside
    // it may lie between the others.
    let mut uses: ogeom_core::FastMap<super::EdgeKey, usize> = ogeom_core::FastMap::default();
    for face in &fused.a.faces {
        for edge in &face.edges {
            *uses.entry(edge.node).or_default() += 1;
        }
    }
    let mut by_edge: ogeom_core::FastMap<super::EdgeKey, usize> = ogeom_core::FastMap::default();
    for (slot, &index) in quiet.iter().enumerate() {
        for edge in &fused.a.faces[fused.pieces[index].face].edges {
            let key = edge.node;
            if uses.get(&key) != Some(&2) {
                continue;
            }
            if let Some(&other) = by_edge.get(&key) {
                let (a, b) = (find(&mut parent, slot), find(&mut parent, other));
                parent[a.max(b)] = a.min(b);
            } else {
                by_edge.insert(key, slot);
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut root_group: ogeom_core::FastMap<usize, usize> = ogeom_core::FastMap::default();
    for (slot, &index) in quiet.iter().enumerate() {
        let root = find(&mut parent, slot);
        let group = *root_group.entry(root).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[group].push(index);
    }
    groups
}
