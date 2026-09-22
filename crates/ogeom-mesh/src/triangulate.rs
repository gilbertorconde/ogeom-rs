//! Triangulating a face.
//!
//! The second half of tessellation. A face is a trimmed region of a surface, so
//! the triangulation is built in the surface's `(u, v)` parameter space — where
//! the region is an ordinary polygon with holes — and then lifted back into
//! space by evaluating the surface at each vertex.
//!
//! # Why parameter space
//!
//! Triangulating in 3D would mean deciding which side of a curved boundary a
//! point falls on, in space, which is the point-in-solid problem. In parameter
//! space the boundary is a closed 2D polygon and the question is a winding
//! count. The surface does the rest.
//!
//! The cost is that parameter space is distorted: equal steps in `(u, v)` cover
//! very different distances near a sphere's pole than near its equator. So the
//! interior points are chosen by measuring deflection *in space* and the
//! triangulation is done in parameter space — measuring where the answer
//! matters, connecting where it is easy.
//!
//! # Watertightness
//!
//! A face's boundary points come from discretizing the *edge's* 3D curve and
//! evaluating the pcurve at those same parameters. Two faces sharing an edge
//! therefore place their boundary vertices at identical spatial positions, and
//! the join has no gap. Discretizing each face's pcurve independently would
//! give each face its own idea of where the edge runs, and the seams would show.

use ogeom_core::{Exact, OgeomResult, Predicates, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve2d, Surface, SurfaceGeometry};
use ogeom_math::{Direction, Point, Point2, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, Triangulation};
use spade::{
    ConstrainedDelaunayTriangulation, Point2 as SpadePoint, Triangulation as _, mitigate_underflow,
};

use crate::discretize::{Deflection, discretize};

/// Triangulate one face.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `face` is not a
/// face or its geometry is missing;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the boundary
/// cannot be triangulated.
pub fn triangulate_face(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Triangulation> {
    // One pass at the caller's chord, and a second only where the first
    // came back short — the same shape as the whole-shape path, so a caller
    // meshing face by face pays for one triangulation per face, not two.
    let nothing = EdgeChords::new();
    let (mesh, verdict) = triangulate_reporting(model, face, deflection, Some(&nothing), tol)?;
    if verdict != Verdict::Short {
        return Ok(mesh);
    }
    triangulate_with(model, face, deflection, None, tol)
}

/// What one pass over a face found.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Verdict {
    /// The boundary enclosed a region and it was drawn.
    Whole,
    /// The boundary crossed itself: what came back is fragments, or nothing.
    /// The edges want drawing finer.
    Short,
}

/// One face, with the finer edge chords the whole shape agreed on.
///
/// `None` means this face is being meshed on its own and may work out its
/// own: alone it has no neighbour to disagree with.
fn triangulate_with(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    finer: Option<&EdgeChords>,
    tol: Tolerances,
) -> OgeomResult<Triangulation> {
    let (mesh, verdict) = triangulate_reporting(model, face, deflection, finer, tol)?;
    if verdict == Verdict::Short && mesh.triangles.is_empty() {
        ogeom_bail!(
            NotDone,
            "the face's boundary enclosed no triangulable region"
        );
    }
    Ok(mesh)
}

/// The surface's unit normal at `(u, v)`, or where the surface is
/// degenerate there — a cone's apex, a sphere's pole, a patch's collapsed
/// corner — the normal a step inside along the vertex's own column.
///
/// A degenerate point has no normal of its own, and `None` for it left
/// the vertex shading black and dragged every normal it was welded with
/// towards nothing. It has a *limit* normal along any line approaching
/// it: the apex of a cone seen up one ruling is that ruling's normal, and
/// a mesh vertex at the apex carries the `u` of the ruling it closes. The
/// step is a millionth of the domain towards its middle, first along
/// `v`, then `u`, then both; `None` only where all three are degenerate.
fn limit_normal(surface: &SurfaceGeometry, u: f64, v: f64, tol: Tolerances) -> Option<Vector> {
    if let Ok(n) = surface.normal_at(u, v, tol) {
        return Some(n.vector());
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let su = if u < f64::midpoint(ua, ub) { 1.0 } else { -1.0 };
    let sv = if v < f64::midpoint(va, vb) { 1.0 } else { -1.0 };
    // Widening steps: a patch whose corner collapses with its tangent —
    // the control rows drawn together and the next row too — is degenerate
    // to first order for a stretch, and a millionth of the domain is still
    // inside it; on a sphere a millimetre across the tangents a millionth
    // in from the pole cross to less than a direction resolves. A
    // hundredth of the domain is past both, and on a patch that small the
    // normal there is the corner's for every purpose.
    for scale in [1e-6, 1e-4, 1e-2] {
        let du = (ub - ua).abs().max(f64::EPSILON) * scale * su;
        let dv = (vb - va).abs().max(f64::EPSILON) * scale * sv;
        let found = [(u, v + dv), (u + du, v), (u + du, v + dv)]
            .into_iter()
            .find_map(|(nu, nv)| surface.normal_at(nu, nv, tol).ok().map(|n| n.vector()));
        if found.is_some() {
            return found;
        }
    }
    None
}

/// One face, and whether the rings it was built from crossed themselves.
///
/// The flag is how the shape-wide pass learns which faces need their edges
/// drawn finer without building every face's rings twice: the rings are
/// already in hand here, and the sweep over them is the only extra cost a
/// face that does not cross ever pays.
fn triangulate_reporting(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    finer: Option<&EdgeChords>,
    tol: Tolerances,
) -> OgeomResult<(Triangulation, Verdict)> {
    triangulate_reporting_from(model, face, deflection, finer, None, tol)
}

/// A face's rings walked ahead of drawing, and the chord its edges want
/// if it is narrower than a few of the caller's.
type Walked = (Trimming, Option<f64>);

/// [`triangulate_reporting`] with the rings already walked, where the
/// caller has them and none of the face's edges were told to draw finer
/// since — the rings depend on nothing else.
fn triangulate_reporting_from(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    finer: Option<&EdgeChords>,
    prepared: Option<Trimming>,
    tol: Tolerances,
) -> OgeomResult<(Triangulation, Verdict)> {
    deflection.validate()?;
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "expected a face");
    }
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;

    let own;
    let finer = match finer {
        Some(shared) => shared,
        None => {
            own = face_chords(model, face, data.surface, surface, deflection, tol)?;
            &own
        }
    };
    let phase = std::time::Instant::now();
    let Trimming {
        rings: uv,
        anchors,
        met,
    } = match prepared {
        Some(trim) => trim,
        None => trimming_rings(model, face, data.surface, surface, deflection, finer, tol)?,
    };
    let rings_ms = phase.elapsed().as_secs_f64() * 1e3;
    let phase = std::time::Instant::now();
    let planar = triangulate_region(&uv, surface, deflection, tol)?;
    let region_ms = phase.elapsed().as_secs_f64() * 1e3;
    let phase = std::time::Instant::now();

    // Whether the triangulator was handed a region at all, asked of what it
    // returned rather than of what it was given. A well-formed triangulation
    // over `b` boundary points in `w` rings has at least `b + 2w - 4`
    // triangles — exactly that where no interior point is added, more where
    // refinement adds them. Fewer is not a coarse answer, it is a different
    // shape: the boundary crossed itself and what came back is disconnected
    // fragments with holes between them.
    //
    // Counting is why it is asked this way round. Sweeping the rings for a
    // crossing is the direct question and costs a sort and an active list
    // per face; measured over a hundred thousand faces it was the whole of
    // an eighteen per cent regression, to catch four bodies. The count is
    // already in hand and exact for the failure that matters.
    //
    // That count is exact only before interior points go in; over a face
    // that takes hundreds of them, a crossing that costs a handful of
    // triangles is lost in the total. So the region also counts its own
    // triangles the moment the boundary is in and nothing else — where the
    // number is exactly `b + 2w - 4` for a boundary that encloses a region
    // and anything else for one that crosses.
    //
    // And a face narrower than a few chords is short too, whatever its
    // count: drawn at the caller's chord its boundary sags by more than
    // the face is wide, and every triangle across the width stands off the
    // surface by that sag. Its edges want a chord under the width, which
    // is [`face_chords`]' first move; it is asked here so the whole-shape
    // pass tells the neighbours to draw those edges finer too.
    let boundary: usize = uv.iter().map(Vec::len).sum();
    let narrow = narrow_chord(surface, &uv, deflection, tol)
        .is_some_and(|want| !edges_already_at(model, face, finer, want));
    let crossed = narrow || planar.crossed || planar.triangles.len() + 4 < boundary + 2 * uv.len();
    if *MESH_DEBUG_REFINE {
        eprintln!(
            "SHORT {} triangles against {boundary} boundary points in {} rings: crossed {crossed} (boundary pass {})",
            planar.triangles.len(),
            uv.len(),
            planar.crossed
        );
    }

    // Boundary vertices take their positions from their edges' own curves —
    // the shared authority — keyed by their exact parameter-space bits.
    let mut anchored: std::collections::HashMap<(u64, u64), Point> =
        std::collections::HashMap::new();
    for (ring, ring_anchor) in uv.iter().zip(&anchors) {
        for (p, a) in ring.iter().zip(ring_anchor) {
            if let Some(point) = a {
                anchored.insert((p.x.to_bits(), p.y.to_bits()), *point);
            }
        }
    }

    // Lift into space. The normal follows the face's orientation, not the
    // surface's: a reversed face presents the other side, and a renderer or a
    // volume computation that ignored that would have the solid inside out.
    // A reflecting placement flips it once more — the mirrored chart's
    // natural normal points the other way through the same flag.
    let flip = (face.orientation() == Orientation::Reversed)
        != !face.location().preserves_handedness(model.datums())?;
    let mut mesh = Triangulation::new();
    mesh.deflection_met = met;
    let mut hits = 0_usize;
    for (u, v) in planar.parameters {
        // Anchors are already in world coordinates — their edges' own
        // placements applied — where surface lifts still need the face's.
        let point = match anchored.get(&(u.to_bits(), v.to_bits())) {
            Some(anchor) => {
                hits += 1;
                *anchor
            }
            None => placement.apply(surface.point_at(u, v, tol)?),
        };
        let normal =
            limit_normal(surface, u, v, tol).map_or(Vector::ZERO, |n| placement.apply_vector(n));
        mesh.positions.push(point);
        mesh.normals.push(if flip { -normal } else { normal });
        mesh.parameters.push((u, v));
    }
    if *MESH_DEBUG_REFINE && (rings_ms + region_ms) > 50.0 {
        eprintln!(
            "PHASE rings {rings_ms:.0}ms region {region_ms:.0}ms lift {:.0}ms  ({} ring points, {} tris)",
            phase.elapsed().as_secs_f64() * 1e3,
            uv.iter().map(Vec::len).sum::<usize>(),
            planar.triangles.len()
        );
    }
    if *MESH_DEBUG {
        eprintln!(
            "DBG anchors map={} hits={} verts={}",
            anchored.len(),
            hits,
            mesh.positions.len()
        );
    }
    mesh.triangles = planar
        .triangles
        .into_iter()
        .map(|t| if flip { [t[0], t[2], t[1]] } else { t })
        .collect();
    Ok((
        mesh,
        if crossed {
            Verdict::Short
        } else {
            Verdict::Whole
        },
    ))
}

/// Triangulate every face below a shape, welded into one mesh.
///
/// # Errors
///
/// As [`triangulate_face`].
pub fn triangulate(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Triangulation> {
    // Faces in two phases, as `tessellate` does: each face is meshed from a
    // model nothing is writing to, in parallel; the pieces are then appended
    // sequentially in face order. The split is what keeps the answer
    // bit-identical at any thread count — scheduling decides only who does
    // which face, never where its triangles land.
    let faces: Vec<Shape> =
        ogeom_topo::explore(model, shape, ogeom_topo::Filter::OfType(ShapeType::Face))?;
    let read_model: &Model = model;

    // Meshed once at the caller's deflection, each face saying whether the
    // rings it was drawn from crossed themselves. Nearly none do, and those
    // faces are finished: the only cost they carry is one sweep over a ring
    // that had to be built anyway.
    // `Some(&nothing)`, not `None`: a face left to itself refines its own
    // edges, which is right when it is meshed alone and wrong here, where
    // its neighbours must be told to refine the same ones. Phase one draws
    // every face at exactly what the caller asked and reports what crossed.
    //
    // Before that, every face's rings are walked at the caller's chord and
    // its width read off them: a face narrower than a few chords wants its
    // edges drawn finer, and so do the faces across those edges. Known
    // before anything is drawn, those chords go into the first pass, and
    // the big faces round a thousand small fillets are drawn once rather
    // than once and again. The rings are kept for the faces they still
    // describe — every face none of whose edges the map names.
    let nothing = EdgeChords::new();
    let prepared: Vec<OgeomResult<Option<Walked>>> =
        ogeom_core::parallel::map_ordered(&faces, |_, face| {
            ogeom_core::progress::checkpoint()?;
            let Some(node) = read_model.node(face) else {
                return Ok(None);
            };
            let NodeData::Face(data) = node.data() else {
                return Ok(None);
            };
            let Some(surface) = read_model.geometry().surface(data.surface) else {
                return Ok(None);
            };
            let trim = trimming_rings(
                read_model,
                face,
                data.surface,
                surface,
                deflection,
                &nothing,
                tol,
            )?;
            let narrow = narrow_chord(surface, &trim.rings, deflection, tol);
            Ok(Some((trim, narrow)))
        });
    let mut finer = EdgeChords::new();
    let mut kept: Vec<Option<Trimming>> = Vec::with_capacity(faces.len());
    for (face, one) in faces.iter().zip(prepared) {
        match one {
            Ok(Some((prep, narrow))) => {
                if let Some(chord) = narrow {
                    for edge in ogeom_topo::explore(
                        read_model,
                        face,
                        ogeom_topo::Filter::OfType(ShapeType::Edge),
                    )? {
                        let held = finer.entry(edge.node().index()).or_insert(chord);
                        *held = held.min(chord);
                    }
                }
                kept.push(Some(prep));
            }
            Ok(None) => kept.push(None),
            Err(_) => kept.push(None),
        }
    }
    let touched = |face: &Shape| -> bool {
        ogeom_topo::explore(
            read_model,
            face,
            ogeom_topo::Filter::OfType(ShapeType::Edge),
        )
        .is_ok_and(|es| es.iter().any(|e| finer.contains_key(&e.node().index())))
    };
    // The rings are handed over by the job that draws the face; a shared
    // slice cannot give them away, so each sits behind a lock it is taken
    // from once.
    let jobs: Vec<(&Shape, std::sync::Mutex<Option<Trimming>>)> = faces
        .iter()
        .zip(kept)
        .map(|(face, prep)| {
            let trim = prep.filter(|_| !touched(face));
            (face, std::sync::Mutex::new(trim))
        })
        .collect();
    let first: Vec<OgeomResult<(Triangulation, Verdict)>> =
        ogeom_core::parallel::map_ordered(&jobs, |_, (face, slot)| {
            ogeom_core::progress::checkpoint()?;
            let trim = slot.lock().ok().and_then(|mut held| held.take());
            triangulate_reporting_from(read_model, face, deflection, Some(&finer), trim, tol)
        });
    let mut computed: Vec<OgeomResult<Triangulation>> = Vec::with_capacity(faces.len());
    let mut crossed: Vec<usize> = Vec::new();
    for (index, one) in first.into_iter().enumerate() {
        match one {
            Ok((mesh, Verdict::Short)) => {
                crossed.push(index);
                computed.push(Ok(mesh));
            }
            Ok((mesh, Verdict::Whole)) => computed.push(Ok(mesh)),
            Err(e) => computed.push(Err(e)),
        }
    }

    // A face whose boundary crossed itself needs its edges drawn finer —
    // and so does every face that shares one of them, or the two sides of
    // that edge arrive with a different number of points, which is a worse
    // crack than the sliver the refinement was for. Only those faces are
    // drawn again.
    if *MESH_DEBUG_REFINE && !crossed.is_empty() {
        eprintln!(
            "REFINE {} of {} faces came up short",
            crossed.len(),
            faces.len()
        );
    }
    if !crossed.is_empty() {
        // On top of the first pass's map, not instead of it: a neighbour
        // drawn again here must still draw the edges the first pass held
        // finer at that chord, or the two sides of one of them disagree.
        // Only the faces touching an edge whose chord *changed* are drawn
        // again.
        let mut changed: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for &index in &crossed {
            let face = &faces[index];
            let Some(node) = read_model.node(face) else {
                continue;
            };
            let NodeData::Face(data) = node.data() else {
                continue;
            };
            let Some(surface) = read_model.geometry().surface(data.surface) else {
                continue;
            };
            for (edge, chord) in
                face_chords(read_model, face, data.surface, surface, deflection, tol)?
            {
                let held = finer.entry(edge).or_insert(f64::INFINITY);
                if chord < *held {
                    *held = chord;
                    changed.insert(edge);
                }
            }
        }
        let again: Vec<usize> = (0..faces.len())
            .filter(|&i| {
                ogeom_topo::explore(
                    model,
                    &faces[i],
                    ogeom_topo::Filter::OfType(ShapeType::Edge),
                )
                .is_ok_and(|es| es.iter().any(|e| changed.contains(&e.node().index())))
            })
            .collect();
        let redone: Vec<OgeomResult<Triangulation>> =
            ogeom_core::parallel::map_ordered(&again, |_, &index| {
                ogeom_core::progress::checkpoint()?;
                triangulate_with(read_model, &faces[index], deflection, Some(&finer), tol)
            });
        for (index, one) in again.into_iter().zip(redone) {
            computed[index] = one;
        }
    }

    let mut mesh = Triangulation::new();
    let mut pieces: Vec<(usize, usize)> = Vec::with_capacity(faces.len());
    for piece in computed {
        let piece = piece?;
        let t0 = mesh.triangles.len();
        mesh.append(&piece);
        pieces.push((t0, mesh.triangles.len()));
    }
    orient_pieces(&mut mesh, &pieces);
    let mesh = mesh.welded(tol);

    // A second, border-only pass at the tolerance the model itself recorded:
    // imported edges carry the file's slop in their widened tolerances, and
    // two faces lifting the same edge through disagreeing geometry land that
    // far apart. Interior edges are already manifold and are not touched.
    let mut reach = 0.0_f64;
    for kind in [ShapeType::Edge, ShapeType::Vertex] {
        for shape in ogeom_topo::explore(model, shape, ogeom_topo::Filter::OfType(kind))? {
            let recorded = model.node(&shape).map_or(0.0, |n| match n.data() {
                NodeData::Edge(d) => d.tolerance.get(),
                NodeData::Vertex(d) => d.tolerance.get(),
                _ => 0.0,
            });
            reach = reach.max(recorded);
        }
    }
    // Floored at a tenth of the chord the caller asked for: a border gap
    // smaller than that is below the resolution of the mesh they accepted,
    // whether or not the model recorded the slop that caused it.
    let reach = reach.max(deflection.chord * 0.1);
    if reach > tol.confusion() {
        let reach = reach + tol.confusion();
        Ok(mesh.border_welded(reach).border_stitched(reach))
    } else {
        Ok(mesh)
    }
}

/// Make the appended face meshes traverse their shared boundaries in
/// opposite directions, flipping as few faces as possible.
///
/// A closed oriented surface walks each interior edge once each way. A file
/// whose face orientation flags disagree with each other still *pairs* every
/// edge — closed by counting, two-sided nowhere — and every integral over the
/// result is reference-dependent garbage. The face flags cannot be judged
/// from the model's wires (a synthesised seam's occurrence directions are
/// chart bookkeeping, not 3D traversal), but the meshes tell the truth:
/// shared boundary vertices are anchored to their edges' own curves, so two
/// faces meeting along an edge carry bitwise-identical positions, and the
/// direction each walks them is right there in the triangles. Faces are
/// flood-filled across those shared runs, each constrained to oppose its
/// neighbour, and each connected component keeps the polarity that flips the
/// minority — flipping means reversing the piece's windings and negating its
/// normals. Conflicts are left standing: this repairs orientation, not
/// topology.
fn orient_pieces(mesh: &mut Triangulation, pieces: &[(usize, usize)]) {
    use std::collections::HashMap;
    type Key = (u64, u64, u64);
    let key = |p: &Point| -> Key { (p.x.to_bits(), p.y.to_bits(), p.z.to_bits()) };

    // Directed boundary edges per piece, keyed by position: only each
    // piece's *border* edges (used once within the piece) face other pieces.
    let mut owners: HashMap<(Key, Key), Vec<(usize, bool)>> = HashMap::new();
    for (i, &(t0, t1)) in pieces.iter().enumerate() {
        let mut inside: HashMap<(u32, u32), u32> = HashMap::new();
        for t in &mesh.triangles[t0..t1] {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *inside.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        for t in &mesh.triangles[t0..t1] {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if inside.get(&(a.min(b), a.max(b))).copied() != Some(1) {
                    continue;
                }
                let (ka, kb) = (
                    key(&mesh.positions[a as usize]),
                    key(&mesh.positions[b as usize]),
                );
                // One undirected key, direction recorded: `true` where this
                // piece walks the smaller key first.
                let (lo, hi, forward) = if ka <= kb {
                    (ka, kb, true)
                } else {
                    (kb, ka, false)
                };
                owners.entry((lo, hi)).or_default().push((i, forward));
            }
        }
    }

    // Manifold constraints: exactly two pieces sharing a run. Same recorded
    // direction means exactly one must flip. A pair of faces shares many
    // runs and slop can corrupt a few, so each pair's constraint is the
    // majority of its runs — and ties abstain rather than guess. Everything
    // is aggregated in sorted order, so the answer does not depend on hash
    // iteration.
    let mut votes: std::collections::BTreeMap<(usize, usize), (u32, u32)> =
        std::collections::BTreeMap::new();
    for list in owners.values() {
        if let [(a, fa), (b, fb)] = list[..]
            && a != b
        {
            let pair = (a.min(b), a.max(b));
            let entry = votes.entry(pair).or_insert((0, 0));
            if fa == fb {
                entry.0 += 1;
            } else {
                entry.1 += 1;
            }
        }
    }
    // Strongest constraints first: a pair vouched for by many runs outranks
    // one hanging on a sliver, so when the graph carries a contradiction —
    // an odd cycle born of slop — the weakest link is the one disbelieved.
    // Union-find with parity keeps the whole resolution order-independent.
    // A pair whose runs split evenly still says the faces touch: it joins
    // the components with the consistent-shell prior, at zero strength, so
    // a fragment cannot drift off and choose its polarity alone.
    let mut constraints: Vec<(u32, usize, usize, bool)> = votes
        .iter()
        .map(|(&(a, b), &(same, opposite))| (same.abs_diff(opposite), a, b, same > opposite))
        .collect();
    constraints.sort_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));

    let mut parent: Vec<usize> = (0..pieces.len()).collect();
    // Whether a piece is flipped relative to its parent.
    let mut parity: Vec<bool> = vec![false; pieces.len()];
    fn find(parent: &mut [usize], parity: &mut [bool], i: usize) -> (usize, bool) {
        if parent[i] == i {
            return (i, false);
        }
        let (root, above) = find(parent, parity, parent[i]);
        parent[i] = root;
        parity[i] ^= above;
        (root, parity[i])
    }
    for &(_, a, b, same_direction) in &constraints {
        let (ra, pa) = find(&mut parent, &mut parity, a);
        let (rb, pb) = find(&mut parent, &mut parity, b);
        // One of a same-direction pair flips: their parities must differ.
        let need = same_direction;
        if ra == rb {
            // Agreeing or contradicting, the die is cast; a contradiction
            // here lost to stronger evidence.
            continue;
        }
        parent[rb] = ra;
        parity[rb] = (pa != pb) != need;
    }

    let mut flip: Vec<bool> = vec![false; pieces.len()];
    let mut components: std::collections::BTreeMap<usize, Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, f) in flip.iter_mut().enumerate() {
        let (root, p) = find(&mut parent, &mut parity, i);
        *f = p;
        components.entry(root).or_default().push(i);
    }
    for members in components.values() {
        let flipped = members.iter().filter(|&&i| flip[i]).count();
        if flipped * 2 > members.len() {
            for &i in members {
                flip[i] = !flip[i];
            }
        }
    }

    for (i, &(t0, t1)) in pieces.iter().enumerate() {
        if !flip[i] {
            continue;
        }
        let mut flipped_vertices: Vec<u32> = Vec::new();
        for t in &mut mesh.triangles[t0..t1] {
            t.swap(1, 2);
            flipped_vertices.extend_from_slice(&t[..]);
        }
        flipped_vertices.sort_unstable();
        flipped_vertices.dedup();
        for v in flipped_vertices {
            mesh.normals[v as usize] = -mesh.normals[v as usize];
        }
    }
}

/// A face's trimming boundary, in its surface's parameter space.
///
/// The outer wire first, then any holes, each as a closed ring of `(u, v)`
/// points with no repeated closing point. A face with no wires covers its
/// surface's whole domain, and gets that rectangle as its boundary.
///
/// Public because parameter-space trimming is not only the triangulator's
/// concern: classifying a point against a face, splitting a face in a boolean,
/// and hidden-line removal all ask the same question of the same rings, and
/// each deriving them separately would be three chances to disagree.
///
/// # Errors
///
/// As [`triangulate_face`].
pub fn face_boundary(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<Point2>>> {
    deflection.validate()?;
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "expected a face");
    }
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };

    Ok(trimming_rings(
        model,
        face,
        data.surface,
        surface,
        deflection,
        &EdgeChords::new(),
        tol,
    )?
    .rings)
}

/// The rings bounding a face in parameter space, and whether every boundary
/// edge met its deflection.
/// Boundary rings with, per ring vertex, the 3D anchor its edge's own curve
/// provides — `None` where an edge has no 3D curve to defer to.
/// The chord each edge must be drawn with, where the caller's is too coarse.
///
/// Keyed by the edge's node, so both faces bounding it look the same value
/// up and sample it identically. Absent means the caller's own chord.
type EdgeChords = std::collections::HashMap<u32, f64>;

/// How many times a face's boundary may be redrawn finer before its
/// crossing is taken to be something the chord cannot fix.
///
/// Six halvings is a chord sixty-four times tighter than asked. A sliver
/// still crossing itself there is degenerate in a way refinement does not
/// reach — two boundaries genuinely on top of one another — and drawing it
/// a seventh time only spends longer to say so.
const REFINEMENTS: usize = 6;

/// What one face needs its own edges drawn with.
///
/// A face narrower than the chord error its boundary is drawn with has a
/// boundary that crosses *itself*. One body of a real assembly carries a
/// quarter-arc sliver forty-five millimetres long and eighteen microns
/// wide: at a tenth of a millimetre the sagitta of each bounding arc is
/// twenty-nine microns, so the inner polyline bulges straight through the
/// outer one, and what reaches the triangulator is not a region at all. It
/// answered with sixteen triangles in fifteen disconnected pieces, and the
/// holes between them were what stopped the body meshing closed.
///
/// The deflection a caller asks for bounds how far the mesh may sit from
/// the surface; it is not a licence to hand the triangulator a polygon that
/// crosses itself. Refining *lowers* that error, so this never breaks the
/// caller's bound.
///
/// Only this face's own edges are named, and only when they need it.
fn face_chords(
    model: &Model,
    face: &Shape,
    id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<EdgeChords> {
    let mut finer = EdgeChords::new();
    let mut chord = deflection.chord;
    // A face narrower than a few chords first: its edges drawn to a
    // fraction of its width, so the boundary's sag is small against it.
    let first = trimming_rings(model, face, id, surface, deflection, &finer, tol)?.rings;
    if let Some(want) = narrow_chord(surface, &first, deflection, tol) {
        chord = want;
        for edge in ogeom_topo::explore(model, face, ogeom_topo::Filter::OfType(ShapeType::Edge))? {
            finer.insert(edge.node().index(), chord);
        }
    }
    for _ in 0..=REFINEMENTS {
        let rings = trimming_rings(model, face, id, surface, deflection, &finer, tol)?.rings;
        let planar = triangulate_region(&rings, surface, deflection, tol)?;
        let boundary: usize = rings.iter().map(Vec::len).sum();
        if !planar.crossed && planar.triangles.len() + 4 >= boundary + 2 * rings.len() {
            break;
        }
        chord *= 0.5;
        for edge in ogeom_topo::explore(model, face, ogeom_topo::Filter::OfType(ShapeType::Edge))? {
            finer.insert(edge.node().index(), chord);
        }
    }
    Ok(finer)
}

/// How many chords wide a face must be for its edges to be drawn at the
/// caller's chord; narrower, they are drawn at the width over this.
///
/// The boundary of a face drawn at chord `c` sags up to `c` between its
/// points, and a triangle from that boundary across the face to the
/// other side stands off the surface by the sag. On a face `w` wide that
/// is `c / w` of the way to standing on end; held to a quarter, the
/// triangles lean at most fourteen degrees.
const NARROW: f64 = 4.0;

/// The chord a face narrower than [`NARROW`] chords wants its edges drawn
/// to, `None` for a face wide enough at the caller's.
///
/// Width is the rings' chart extent each way scaled by the mean tangent
/// length that way, the smaller of the two: right for a strip along a
/// chart axis, an over-estimate for one across the chart, which then
/// keeps the caller's chord.
fn narrow_chord(
    surface: &SurfaceGeometry,
    rings: &[Vec<Point2>],
    deflection: Deflection,
    tol: Tolerances,
) -> Option<f64> {
    let (lo, hi) = chart_extent(rings);
    if !(lo.x.is_finite() && hi.x.is_finite() && lo.y.is_finite() && hi.y.is_finite()) {
        return None;
    }
    let (mut du, mut dv, mut n) = (0.0_f64, 0.0_f64, 0usize);
    for i in 0..=2 {
        for j in 0..=2 {
            let u = lo.x + (hi.x - lo.x) * (0.25 + 0.25 * f64::from(i));
            let v = lo.y + (hi.y - lo.y) * (0.25 + 0.25 * f64::from(j));
            if let Ok((a, b)) = surface.d1_at(u, v, tol) {
                du += a.magnitude();
                dv += b.magnitude();
                n += 1;
            }
        }
    }
    if n == 0 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    let width = ((hi.x - lo.x) * du / n as f64).min((hi.y - lo.y) * dv / n as f64);
    if !width.is_finite() || width <= tol.confusion() || width >= deflection.chord * NARROW {
        return None;
    }
    Some(width / NARROW)
}

/// Whether every edge of `face` is already held to `want` or finer.
fn edges_already_at(model: &Model, face: &Shape, finer: &EdgeChords, want: f64) -> bool {
    ogeom_topo::explore(model, face, ogeom_topo::Filter::OfType(ShapeType::Edge)).is_ok_and(
        |edges| {
            edges
                .iter()
                .all(|e| finer.get(&e.node().index()).is_some_and(|&c| c <= want))
        },
    )
}

/// A face's trimming rings, walked, folded and cleaned.
struct Trimming {
    /// The outer ring first, then the holes, each closed without a repeated
    /// closing point.
    rings: Vec<Vec<Point2>>,
    /// Each ring point's position in space where its edge's own curve put
    /// it, `None` where it was made up.
    anchors: Vec<Vec<Option<Point>>>,
    /// Whether every edge's polyline honoured the deflection.
    met: bool,
}

/// One walked ring: chart points, anchors, deflection honesty, the ambiguous
/// whole-period folds taken, and the half-period ties left undecided.
type WalkedRing = (
    Vec<Point2>,
    Vec<Option<Point>>,
    bool,
    Vec<(usize, f64)>,
    Vec<usize>,
);

fn trimming_rings(
    model: &Model,
    face: &Shape,
    id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    deflection: Deflection,
    finer: &EdgeChords,
    tol: Tolerances,
) -> OgeomResult<Trimming> {
    let mut rings = Vec::new();
    let mut ring_anchors = Vec::new();
    let mut ring_folds: Vec<Vec<(usize, f64)>> = Vec::new();
    let mut ring_ties: Vec<Vec<usize>> = Vec::new();
    let mut met = true;
    for wire in model.ordered_children_of(face)? {
        let (ring, anchors, ring_met, folds, ties) =
            boundary_ring(model, &wire, id, deflection, finer, tol)?;
        met &= ring_met;
        if ring.len() >= 3 {
            rings.push(ring);
            ring_anchors.push(anchors);
            ring_folds.push(folds);
            ring_ties.push(ties);
        }
    }

    // A band between two *wound* rings: each boundary winds the periodic
    // direction once — a bore wall whose ends are a full rim and a staircase
    // of arcs — and neither ring closes on its own. In the unrolled chart
    // the face is the strip between the two chains, so the pair is merged
    // into one ring by two joining runs standing exactly one period apart:
    // their lifted points are the same 3D points, and the weld closes them
    // the way it closes a seam.
    if geometry_winds(surface) {
        use ogeom_geom::Surface as _;
        let ((ua, ub), _) = surface.domain();
        let period = ub - ua;
        let du_of = |ring: &[Point2]| -> f64 { ring.last().map_or(0.0, |l| l.x - ring[0].x) };
        let open_wound: Vec<usize> = rings
            .iter()
            .enumerate()
            .filter(|(_, ring)| {
                ring.len() >= 3
                    && (du_of(ring).abs() - period).abs() <= period * 1e-3
                    && ring
                        .last()
                        .is_some_and(|l| ring[0].distance(*l) > period * 0.5)
            })
            .map(|(i, _)| i)
            .collect();
        if let [i, j] = open_wound[..]
            && du_of(&rings[i]).signum() != du_of(&rings[j]).signum()
        {
            let sign = du_of(&rings[i]).signum();
            let b = rings.remove(j);
            let b_anchors = ring_anchors.remove(j);
            ring_folds.remove(j);
            ring_ties.remove(j);
            let a = &mut rings[i];
            let a_anchors = &mut ring_anchors[i];
            let a_last = a.last().copied().unwrap_or(a[0]);
            // Whole periods only: the second chain slides along the unrolled
            // chart until its start stands nearest the first chain's end.
            let shift = ((a_last.x - b[0].x) / period).round() * period;
            let b_first = Point2::new(b[0].x + shift, b[0].y);
            let steps = 8;
            for k in 1..steps {
                let f = f64::from(k) / f64::from(steps);
                a.push(Point2::new(
                    a_last.x + (b_first.x - a_last.x) * f,
                    a_last.y + (b_first.y - a_last.y) * f,
                ));
                a_anchors.push(None);
            }
            for (p, anchor) in b.iter().zip(&b_anchors) {
                a.push(Point2::new(p.x + shift, p.y));
                a_anchors.push(*anchor);
            }
            // The way back: the same run, one period over, walked the other
            // way — the two runs lift to identical points.
            for k in (1..steps).rev() {
                let f = f64::from(k) / f64::from(steps);
                a.push(Point2::new(
                    a_last.x + (b_first.x - a_last.x) * f - sign * period,
                    a_last.y + (b_first.y - a_last.y) * f,
                ));
                a_anchors.push(None);
            }
        }
    }

    // A single ring still winding after the pairing had no partner, and if
    // it cannot close it was mis-folded: a jump of exactly one period is the
    // same 3D point, and the greedy fold that zeroed it may have wound the
    // ring instead. Undoing one such ambiguous fold — translating the walk
    // from there on by a period — closes the ring; the join it reopens
    // still lifts to a single point.
    if geometry_winds(surface) {
        use ogeom_geom::Surface as _;
        let ((ua, ub), _) = surface.domain();
        let period = ub - ua;
        for ((ring, folds), ties) in rings.iter_mut().zip(&ring_folds).zip(&ring_ties) {
            let Some(last) = ring.last().copied() else {
                continue;
            };
            let du = last.x - ring[0].x;
            let k = (du / period).round();
            if k == 0.0
                || (du - k * period).abs() > period * 1e-3
                || ring[0].distance(last) <= period * 1e-3
            {
                continue;
            }
            // The last ambiguous fold whose applied shift matches the
            // winding is the one to undo.
            if let Some(&(fold_start, _)) = folds
                .iter()
                .rev()
                .find(|(_, shift)| (shift - k * period).abs() <= period * 1e-6)
                && fold_start < ring.len()
            {
                // Translating from the last undecided tie before the fold —
                // where the walk first guessed — keeps the period jump on
                // the degenerate row, where it lifts to nothing.
                let start = ties
                    .iter()
                    .rev()
                    .find(|&&t| t < fold_start)
                    .copied()
                    .unwrap_or(fold_start);
                for p in &mut ring[start..] {
                    p.x -= k * period;
                }
            }
        }
    }

    // Points closer than the triangulator's own resolution make it refuse
    // the constraint; a real file's chart can carry them. The consecutive
    // near-duplicates collapse, anchors staying aligned.
    for (ring, anchors) in rings.iter_mut().zip(ring_anchors.iter_mut()) {
        if ring.len() < 3 {
            continue;
        }
        let mut extent = 0.0_f64;
        for pair in ring.windows(2) {
            extent = extent.max(pair[0].distance(pair[1]));
        }
        let eps = extent.mul_add(1e-9, 1e-12);
        let mut kept_ring = Vec::with_capacity(ring.len());
        let mut kept_anchors = Vec::with_capacity(anchors.len());
        for (p, a) in ring.iter().zip(anchors.iter()) {
            if kept_ring
                .last()
                .is_some_and(|held: &Point2| held.distance(*p) <= eps)
            {
                continue;
            }
            kept_ring.push(*p);
            kept_anchors.push(*a);
        }
        if kept_ring.len() > 2
            && let (Some(first), Some(last)) = (kept_ring.first(), kept_ring.last())
            && first.distance(*last) <= eps
        {
            kept_ring.pop();
            kept_anchors.pop();
        }
        *ring = kept_ring;
        *anchors = kept_anchors;
    }
    // Folded across a join, a ring on a *closed* surface can come to rest a
    // whole period outside the domain. A periodic surface would not mind —
    // it wraps — but a surface that merely closes on itself evaluates only
    // where its knots are, and refuses everywhere else; three bodies of one
    // assembly stopped meshing that way, every point of one ring a turn
    // past the end. The fold kept the ring continuous, which is the part
    // that matters, and a rigid slide by whole periods keeps it so: the
    // same points on the surface, named inside the chart.
    {
        use ogeom_geom::Surface as _;
        let ((ua, ub), (va, vb)) = surface.domain();
        let slides = [
            (!surface.is_periodic_u() && surface.is_closed_u(tol), ua, ub),
            (!surface.is_periodic_v() && surface.is_closed_v(tol), va, vb),
        ];
        for (across, (closed, lo, hi)) in [true, false].into_iter().zip(slides) {
            if !closed || hi <= lo {
                continue;
            }
            let span = hi - lo;
            for ring in &mut rings {
                let read = |p: &Point2| if across { p.x } else { p.y };
                let (least, most) = ring
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), p| {
                        (a.min(read(p)), b.max(read(p)))
                    });
                if !(least.is_finite() && most.is_finite()) || most - least > span * (1.0 + 1e-9) {
                    continue;
                }
                // A hair past the end is fit noise, not a period: the knots
                // take it, and rounding it up to a whole turn would carry the
                // ring a period the wrong way — which is exactly what it did
                // to the face this was written for, before the slack.
                let slack = span * 1e-6;
                let turns = if least < lo - slack {
                    ((lo - least) / span).ceil()
                } else if most > hi + slack {
                    -((most - hi) / span).ceil()
                } else {
                    0.0
                };
                if turns == 0.0 {
                    continue;
                }
                for p in ring.iter_mut() {
                    if across {
                        p.x += turns * span;
                    } else {
                        p.y += turns * span;
                    }
                }
            }
        }
    }
    if *MESH_DEBUG {
        for (i, (ring, anchors)) in rings.iter().zip(&ring_anchors).enumerate() {
            eprintln!("DBG ring {i}: {} points", ring.len());
            for (p, a) in ring.iter().zip(anchors) {
                eprintln!(
                    "DBG   uv({:.5},{:.5}) anchor {}",
                    p.x,
                    p.y,
                    a.map_or("-".to_string(), |q| format!(
                        "({:.4},{:.4},{:.4})",
                        q.x, q.y, q.z
                    ))
                );
            }
        }
    }
    // A run out along an edge and straight back along it — a ring that
    // reads `p, q, p` — is a spike into the region that bounds nothing:
    // the file's way of drawing a slit of no width at all on the face's
    // own boundary. It triangulates to two hairs and one vertex too many,
    // which is one triangle more than a boundary that encloses a region
    // has, and it is not a crossing. Off it comes, out to in, until the
    // ring reverses nowhere.
    // And two consecutive points a millionth of the ring's size apart are
    // one point: the end of the last edge and the start of the first, each
    // where its own curve put the shared vertex, a file's slop apart. Kept
    // both, the second sits on the first's next segment to the last bit,
    // and a constraint through a vertex is one the triangulation refuses.
    for (ring, anchors) in rings.iter_mut().zip(ring_anchors.iter_mut()) {
        let reach = chart_reach(ring);
        merge_near_duplicates(ring, anchors, reach);
        remove_spikes(ring, anchors, reach);
    }
    rings.retain(|r| r.len() >= 3);
    ring_anchors.retain(|a| a.len() >= 3);

    // An inner ring thinner than a micron is a slit, not a hole. A real file
    // draws one by running out along two arcs and back along two splines
    // fitted to the same arcs: a loop three millimetres long and a fifth of
    // a micron wide, enclosing nothing, which the triangulator can only
    // read as a tangle — one face carrying five of them drew with twelve
    // holes it does not have. Measured in space through the ring's own
    // anchors, so a chart's units do not enter into it; a ring not anchored
    // end to end is left alone, and so is the outer ring, whatever its
    // width, since a face that is itself a slit is a different question.
    // The thinnest real feature in the assembly that showed this is twenty
    // microns across, twenty times the cutoff.
    if rings.len() > 1 {
        let chart_area = |ring: &[Point2]| -> f64 {
            let mut a = 0.0;
            for i in 0..ring.len() {
                let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
                a += p.x * q.y - q.x * p.y;
            }
            a.abs()
        };
        let outer = (0..rings.len())
            .max_by(|&i, &j| chart_area(&rings[i]).total_cmp(&chart_area(&rings[j])))
            .unwrap_or(0);
        let width = |anchors: &[Option<Point>]| -> Option<f64> {
            let pts: Option<Vec<Point>> = anchors.iter().copied().collect();
            let pts = pts?;
            let mut normal = Vector::ZERO;
            let mut perimeter = 0.0;
            for i in 0..pts.len() {
                let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                normal += a.to_vector().cross(b.to_vector());
                perimeter += a.distance(b);
            }
            (perimeter > 0.0).then(|| normal.magnitude() * 0.5 / perimeter)
        };
        let keep: Vec<bool> = (0..rings.len())
            .map(|i| {
                i == outer || width(&ring_anchors[i]).is_none_or(|w| w >= tol.confusion() * 1e4)
            })
            .collect();
        let mut it = keep.iter();
        rings.retain(|_| *it.next().unwrap_or(&true));
        let mut it = keep.iter();
        ring_anchors.retain(|_| *it.next().unwrap_or(&true));
    }

    if rings.is_empty() {
        // A face with no wires covers its surface's whole domain, so the domain
        // rectangle is the boundary.
        let ring = domain_ring(surface, deflection, tol);
        ring_anchors.push(vec![None; ring.len()]);
        rings.push(ring);
    }
    Ok(Trimming {
        rings,
        anchors: ring_anchors,
        met,
    })
}

/// Within this of each other, two chart points of a ring are one point: a
/// millionth of the ring's extent, well under any feature and well over
/// the slop two edges leave at the vertex they share.
fn chart_reach(ring: &[Point2]) -> f64 {
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in ring {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let extent = (hi.x - lo.x).max(hi.y - lo.y);
    if extent.is_finite() && extent > 0.0 {
        extent * 1e-6
    } else {
        0.0
    }
}

/// Merge consecutive ring points within `reach` of each other,
/// cyclically; the earlier point and its anchor stay.
fn merge_near_duplicates(ring: &mut Vec<Point2>, anchors: &mut Vec<Option<Point>>, reach: f64) {
    if ring.len() < 2 || anchors.len() != ring.len() || reach <= 0.0 {
        return;
    }
    let mut i = 0;
    while i < ring.len() && ring.len() >= 2 {
        let next = (i + 1) % ring.len();
        let (p, q) = (ring[i], ring[next]);
        if (p.x - q.x).hypot(p.y - q.y) <= reach {
            // The later point goes — the last one when the ring's end
            // repeats its start — and the earlier is looked at again, in
            // case it now sits next to another near-duplicate.
            if next == 0 {
                ring.remove(i);
                anchors.remove(i);
                break;
            }
            ring.remove(next);
            anchors.remove(next);
        } else {
            i += 1;
        }
    }
}

/// Strip every `p, q, p` from a ring — a point stepped out to and straight
/// back from, the two `p` within `reach` of each other — with its anchors,
/// until none is left. Cyclic: the ring's last point is its first's
/// neighbour.
fn remove_spikes(ring: &mut Vec<Point2>, anchors: &mut Vec<Option<Point>>, reach: f64) {
    loop {
        let n = ring.len();
        if n < 3 || anchors.len() != n {
            return;
        }
        let Some(tip) = (0..n).find(|&i| {
            let (before, after) = (ring[(i + n - 1) % n], ring[(i + 1) % n]);
            (before.x - after.x).hypot(before.y - after.y) <= reach
        }) else {
            return;
        };
        // The tip and one of its two identical neighbours go.
        let neighbour = (tip + 1) % n;
        let (first, second) = if tip < neighbour {
            (neighbour, tip)
        } else {
            (tip, neighbour)
        };
        ring.remove(first);
        anchors.remove(first);
        ring.remove(second);
        anchors.remove(second);
    }
}

/// Whether a surface is periodic in `u` alone — the charts on which a wound
/// ring cannot close against a translate of itself and needs a partner.
fn geometry_winds(surface: &SurfaceGeometry) -> bool {
    use ogeom_geom::Surface as _;
    surface.is_periodic_u() && !surface.is_periodic_v()
}

/// Whether a point in parameter space lies inside the region `rings` bound.
///
/// Even-odd winding: inside the outer ring and outside every hole. The rings
/// come from wires whose direction already encodes outer from inner, but
/// counting crossings does not depend on that being right, which makes it
/// robust to a wire that was built the wrong way round.
///
/// Says nothing about a point *on* a ring — the crossing count of a boundary
/// point is whichever side rounding puts it. A caller that cares has to measure
/// its distance to the boundary and decide, which is what classification does.
///
/// Decided with exact predicates. See [`inside_boundary_with`].
#[must_use]
pub fn inside_boundary(rings: &[Vec<Point2>], p: Point2) -> bool {
    inside_boundary_with::<Exact>(rings, p)
}

/// As [`inside_boundary`], with the predicate implementation named.
///
/// This is the seam `docs/DATA_MODEL.md` §9 describes, and it is here rather
/// than anywhere else because this is where the *combinatorial* decision is.
/// Whether a point is inside a ring is not a measurement that can be a little
/// wrong: it decides whether a triangle is kept or dropped, so an error near a
/// boundary is a hole in the mesh rather than a slightly misplaced one.
///
/// The naive form divides to find where an edge crosses the sampling ray, and
/// that division cancels catastrophically for a point nearly on the edge.
/// `orient2d` answers the same question with no division at all, and
/// [`Exact`] answers it correctly however close the point is.
#[must_use]
pub fn inside_boundary_with<P: Predicates>(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for ring in rings {
        if crosses_odd_times::<P>(ring, p) {
            inside = !inside;
        }
    }
    inside
}

/// The boundary of one wire, in the face's parameter space.
///
/// Each edge is discretized in *space* and its pcurve evaluated at the resulting
/// parameters, so two faces sharing the edge agree on where its points are.
fn boundary_ring(
    model: &Model,
    wire: &Shape,
    surface: ogeom_topo::SurfaceId,
    deflection: Deflection,
    finer: &EdgeChords,
    tol: Tolerances,
) -> OgeomResult<WalkedRing> {
    let mut ring: Vec<Point2> = Vec::new();
    let mut anchors: Vec<Option<Point>> = Vec::new();
    // Ambiguous folds: where an edge was translated a whole period to
    // continue the walk from a start that already coincided modulo the
    // period — the fold was one of two defensible choices, recorded so a
    // mis-wound ring can be unwound.
    let mut folds: Vec<(usize, f64)> = Vec::new();
    // Half-period jumps whose side could not be decided when walked.
    let mut ties: Vec<usize> = Vec::new();
    let mut met = true;
    // Whether each chart direction comes back on itself — periodic, or
    // closed without repeating. Asked once here: closure on a spline is a
    // walk down a control column, and asking it at every edge of every face
    // of a hundred-thousand-face assembly was a tenth of the meshing time.
    let (wraps_u, wraps_v) = model
        .geometry()
        .surface(surface)
        .map_or((false, false), |g| {
            use ogeom_geom::Surface as _;
            (
                g.is_periodic_u() || g.is_closed_u(tol),
                g.is_periodic_v() || g.is_closed_v(tol),
            )
        });

    // Start the walk off a seam if the wire allows it: a seam's side is
    // chosen by continuity with the point already walked to, and continuity
    // needs something to continue from. The ring is cyclic, so rotating the
    // walk changes nothing it reports.
    let mut children = model.ordered_children_of(wire)?;
    let is_seam = |model: &Model, e: &Shape| -> bool {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface, e.location()))
            .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }))
    };
    if let Some(start) = children.iter().position(|e| !is_seam(model, e)) {
        children.rotate_left(start);
    }

    // The column each seam edge's first traversal effectively walked — after
    // folding — and whether its two sides differ in u or in v. A seam bounds
    // its face twice, and the two traversals must bracket the ring exactly
    // one period apart; the walk checks the second against this record.
    let mut seam_walked: std::collections::HashMap<ogeom_topo::TShapeId, Point2> =
        std::collections::HashMap::new();

    for edge in children {
        let Some(node) = model.node(&edge) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        // Whether this edge is a seam, and if so whether its two sides
        // differ in u (true) or in v.
        let mut seam: Option<bool> = None;
        let (pcurve_id, pcurve_range) = match data.pcurve_for(surface, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
            // A seam edge runs along a closed surface's join and bounds its
            // face twice — up one side of the parameter rectangle and down
            // the other. Which side this occurrence takes is decided by the
            // ring itself: the side whose oriented start continues the point
            // already walked to. Orientation flags cannot answer it — a
            // reversed face flips every occurrence while the chart columns
            // stay where they were built — but the chart can.
            Some(EdgeRepr::Seam {
                forward,
                reversed,
                range,
                ..
            }) => {
                let (f, r) = (*forward, *reversed);
                let side_start = |id: ogeom_topo::PCurveId| -> Option<Point2> {
                    model.geometry().pcurve(id)?.point_at(range.0, tol).ok()
                };
                seam = Some(match (side_start(f), side_start(r)) {
                    (Some(a), Some(b)) => (a.x - b.x).abs() >= (a.y - b.y).abs(),
                    _ => true,
                });
                let picked = if let Some(last) = ring.last().copied() {
                    let start_of = |id: ogeom_topo::PCurveId| -> Option<Point2> {
                        let pc = model.geometry().pcurve(id)?;
                        let t = if edge.orientation() == Orientation::Reversed {
                            range.1
                        } else {
                            range.0
                        };
                        pc.point_at(t, tol).ok()
                    };
                    match (start_of(f), start_of(r)) {
                        (Some(a), Some(b)) => {
                            if last.distance(a) <= last.distance(b) {
                                f
                            } else {
                                r
                            }
                        }
                        _ => f,
                    }
                } else if edge.orientation() == Orientation::Reversed {
                    r
                } else {
                    f
                };
                (picked, *range)
            }
            _ => ogeom_bail!(
                Construction,
                "edge has no pcurve on this face, so the face cannot be \
                 triangulated in its own parameter space"
            ),
        };
        let Some(pcurve) = model.geometry().pcurve(pcurve_id) else {
            ogeom_bail!(Dangling, "pcurve is not in this model");
        };

        // Sample where the *3D* curve says to, so an adjacent face lands on
        // the same points — and *anchor* the boundary vertices to that curve
        // too: two faces sharing an edge lift the same parameters through
        // different surfaces, and on an imported file those surfaces
        // disagree by the file's own slop. The edge is the shared authority,
        // so its points are the positions both faces use, and the weld is a
        // matter of identity rather than luck.
        let mut edge_anchors: Vec<Option<Point>> = Vec::new();
        // An edge drawn finer is drawn finer for *every* face that bounds
        // it, which is the whole point: the two sides must agree point for
        // point or the weld has nothing to join.
        let along = match finer.get(&edge.node().index()) {
            Some(chord) => Deflection {
                chord: *chord,
                ..deflection
            },
            None => deflection,
        };
        let samples = match sample_parameters(model, data, along, tol)? {
            Some((parameters, edge_met)) => {
                met &= edge_met;
                if let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d()
                    && let Some(geometry) = model.geometry().curve(*curve)
                    && let Ok(edge_placement) = edge.transform(model.datums())
                {
                    for t in &parameters {
                        edge_anchors.push(
                            geometry
                                .point_at(*t, tol)
                                .ok()
                                .map(|p| edge_placement.apply(p)),
                        );
                    }
                }
                map_to_pcurve(&parameters, data, pcurve_range)
            }
            // No 3D curve to defer to — the pcurve's own shape, measured in
            // space through the surface, so the chord tolerance means the
            // same thing it means everywhere else.
            None => {
                let Some(geometry) = model.geometry().surface(surface) else {
                    ogeom_bail!(Dangling, "face refers to a surface not in this model");
                };
                let (_, parameters) = crate::discretize::discretize_on_surface(
                    pcurve,
                    pcurve_range,
                    geometry,
                    deflection,
                    tol,
                )?;
                parameters
            }
        };

        let mut points: Vec<Point2> = samples
            .iter()
            .map(|u| pcurve.point_at(*u, tol))
            .collect::<OgeomResult<_>>()?;
        if edge.orientation() == Orientation::Reversed {
            points.reverse();
            edge_anchors.reverse();
        }
        if edge_anchors.len() != points.len() {
            edge_anchors = vec![None; points.len()];
        }
        // The ends of an edge belong to its *vertices* — the one authority
        // every face and every neighbouring edge shares — but only within the
        // tolerance the vertex itself records. An imported curve ends within
        // the vertex's widened tolerance of it, and lifting the ends through
        // the curve alone would leave each corner split as many ways as there
        // are curves meeting there; a vertex that sits *beyond* its stated
        // tolerance from the curve is not describing the curve's end at all,
        // and the curve stays the authority.
        if !points.is_empty()
            && let Ok(vs) = model.children_of(&edge)
            && vs.len() >= 2
            && let Ok(edge_placement) = edge.transform(model.datums())
        {
            let point_of = |v: &Shape| -> Option<(Point, f64)> {
                let data = model.node(v)?.data().as_vertex()?;
                Some((edge_placement.apply(data.point), data.tolerance.get()))
            };
            let (from, to) = if edge.orientation() == Orientation::Reversed {
                (&vs[vs.len() - 1], &vs[0])
            } else {
                (&vs[0], &vs[vs.len() - 1])
            };
            if let Some((p, within)) = point_of(from)
                && let Some(a) = edge_anchors.first_mut()
                && a.is_none_or(|end| end.distance(p) <= within + tol.confusion())
            {
                *a = Some(p);
            }
            if let Some((p, within)) = point_of(to)
                && let Some(a) = edge_anchors.last_mut()
                && a.is_none_or(|end| end.distance(p) <= within + tol.confusion())
            {
                *a = Some(p);
            }
        }
        // Fold onto the branch that continues the ring. Two faces sharing an
        // edge share its pcurve, and on a periodic surface the pcurve sits in
        // *one* face's window: a cylinder split into two halves has a ruling
        // at u = 0 that the other half needs at u = 2pi. The chart cannot
        // store both; continuity with the ring being walked recovers the
        // right branch, exactly as the seam sides are chosen.
        if let Some(last) = ring.last().copied()
            && let Some(first) = points.first().copied()
            && let Some(geometry) = model.geometry().surface(surface)
        {
            use ogeom_geom::Surface as _;
            let ((ua, ub), (va, vb)) = geometry.domain();
            // Nearest whole period, ties broken toward *not moving*: a jump
            // of exactly half a period is what a boundary crossing a
            // degenerate row looks like — two rulings into an apex stand
            // half a turn apart, and the connecting run along the apex row
            // lifts to nothing — and folding it would drag the edge a full
            // period from the column its own projection put it on.
            let whole_periods = |gap: f64, span: f64| -> f64 {
                let r = gap / span;
                if (r.fract().abs() - 0.5).abs() <= 1e-9 {
                    r.trunc() * span
                } else {
                    r.round() * span
                }
            };
            let mut shift = Point2::new(0.0, 0.0);
            // A repeated seam edge folds like any other — but never onto its
            // own first traversal. The two traversals bound the face up one
            // side of the chart and down the other, one period apart, and at
            // a degenerate row continuity cannot say so: a cone walked to
            // its apex reaches a corner that maps to the whole row, both
            // sides continue it equally, and folding by nearness closes the
            // ring over nothing. The record decides instead: land exactly a
            // period from the first walk, on the side the ring occupies.
            let prior = seam.and_then(|_| seam_walked.get(&edge.node())).copied();
            if wraps_u && (ub - ua) > 0.0 {
                let span = ub - ua;
                let gap = last.x - first.x;
                shift.x = whole_periods(gap, span);
                let mut bracketed = false;
                if seam == Some(true)
                    && let Some(prior) = prior
                    && (first.x + shift.x - prior.x).abs() < span * 0.5
                {
                    let (lo, hi) = ring
                        .iter()
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                            (lo.min(p.x), hi.max(p.x))
                        });
                    // Bracketing is for the face that wraps the period — a
                    // band whose rims run the whole way round, its two seam
                    // columns a period apart. A *slit* uses one edge twice
                    // without wrapping: the ring stays in a fraction of the
                    // chart, both traversals stand on one column, and
                    // forcing them apart winds the ring, invites a pole row
                    // it never touches, and meshes the complement of the
                    // face — the issue #37 screw head. The ring's own reach
                    // says which face this is.
                    if hi - lo >= span * 0.5 {
                        let side = if f64::midpoint(lo, hi) >= prior.x {
                            1.0
                        } else {
                            -1.0
                        };
                        shift.x = prior.x + side * span - first.x;
                        bracketed = true;
                    }
                }
                if !bracketed {
                    if shift.x != 0.0 && (gap - shift.x).abs() <= span * 1e-6 {
                        // The start already stood a whole period from the walk —
                        // the same 3D point — so this fold is a choice, not a
                        // repair; recorded so a mis-wound ring can be unwound.
                        folds.push((ring.len(), shift.x));
                    } else if ((gap / span).fract().abs() - 0.5).abs() <= 1e-9 {
                        // A half-period tie: either side of the degenerate row
                        // was defensible, and if the ring comes out wound the
                        // unwinding starts here rather than at the later fold.
                        ties.push(ring.len());
                    }
                }
            }
            if wraps_v && (vb - va) > 0.0 {
                let span = vb - va;
                shift.y = whole_periods(last.y - first.y, span);
                if seam == Some(false)
                    && let Some(prior) = prior
                    && (first.y + shift.y - prior.y).abs() < span * 0.5
                {
                    let (lo, hi) = ring
                        .iter()
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                            (lo.min(p.y), hi.max(p.y))
                        });
                    let side = if f64::midpoint(lo, hi) >= prior.y {
                        1.0
                    } else {
                        -1.0
                    };
                    shift.y = prior.y + side * span - first.y;
                }
            }
            if shift.x != 0.0 || shift.y != 0.0 {
                for p in &mut points {
                    p.x += shift.x;
                    p.y += shift.y;
                }
            }
        }
        if seam.is_some()
            && let Some(first) = points.first().copied()
        {
            seam_walked.entry(edge.node()).or_insert(first);
        }
        // The previous edge already contributed the shared vertex — but only
        // where the chart agrees it is shared. Two rulings meeting at an
        // apex share the *vertex* while standing apart in the chart, and the
        // run between them along the degenerate row is boundary the ring
        // needs: dropping its start would cut the corner straight through
        // the face's interior.
        //
        // Two things have to hold for a gap to be such a run, and neither
        // alone is enough. It has to be wide: two ends that disagree by the
        // file's slop stand a micron over a radius apart, and a thousandth
        // of the period is three orders above that. And what lies between
        // has to be degenerate: a run along an apex row lifts to one point
        // the whole way, so its chart midpoint lands on the shared vertex,
        // where a wide gap on a live row lifts to somewhere the width of the
        // gap away. Width alone kept slop on fitted splines; the lift alone
        // kept every gap too small for its midpoint to land anywhere else.
        // And no fraction of the period alone is right at all — one real
        // part's rulings stand exactly a quarter turn apart, which a
        // quarter-period test read as not apart, and the face lost the
        // triangle at its apex.
        let keep_gap = if let (Some(last), Some(first)) = (ring.last(), points.first()) {
            model.geometry().surface(surface).is_some_and(|geometry| {
                use ogeom_geom::Surface as _;
                let ((ua, ub), (va, vb)) = geometry.domain();
                let wide = (wraps_u && (last.x - first.x).abs() > (ub - ua) * 1e-3)
                    || (wraps_v && (last.y - first.y).abs() > (vb - va) * 1e-3);
                if !wide {
                    return false;
                }
                let mid = Point2::new(
                    f64::midpoint(last.x, first.x),
                    f64::midpoint(last.y, first.y),
                );
                let reach = tol.confusion() * 1e4;
                match (
                    geometry.point_at(last.x, last.y, tol),
                    geometry.point_at(mid.x, mid.y, tol),
                    geometry.point_at(first.x, first.y, tol),
                ) {
                    (Ok(a), Ok(m), Ok(b)) => a.distance(m) <= reach && m.distance(b) <= reach,
                    _ => false,
                }
            })
        } else {
            false
        };
        if !ring.is_empty() && !points.is_empty() && !keep_gap {
            points.remove(0);
            edge_anchors.remove(0);
        }
        ring.extend(points);
        anchors.extend(edge_anchors);
    }

    // A closed ring repeats its first point at the end; the triangulator wants
    // it named once.
    if ring.len() > 2
        && let (Some(first), Some(last)) = (ring.first().copied(), ring.last().copied())
    {
        // Equal in the chart, or the same vertex in space: the last edge's
        // curve ends where the first edge's begins to within the file's
        // slop — up to ten microns in a real assembly, recorded on the
        // vertex as its widened tolerance. Kept as two points, the ring
        // closes with a fold back over its own first segment, a crossing
        // a fraction of a micron deep that the triangulation refuses as a
        // constraint and the face is then drawn six times finer for.
        // The same vertex in space is not enough on its own: a ring that
        // winds a periodic chart ends a whole period from where it began
        // and lifts to the same point, and that closing is a seam, not
        // slop. Close in the chart too — within a thousandth of the ring's
        // own extent — or the ring is left to the winding rule below.
        let extent = ring.iter().fold(
            (
                Point2::new(f64::INFINITY, f64::INFINITY),
                Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            ),
            |(lo, hi), p| {
                (
                    Point2::new(lo.x.min(p.x), lo.y.min(p.y)),
                    Point2::new(hi.x.max(p.x), hi.y.max(p.y)),
                )
            },
        );
        let extent = (extent.1.x - extent.0.x).max(extent.1.y - extent.0.y);
        let near_in_chart = (first.x - last.x).hypot(first.y - last.y) <= extent * 1e-3;
        let same_vertex = near_in_chart
            && match (anchors.first(), anchors.last()) {
                (Some(Some(a)), Some(Some(b))) => a.distance(*b) <= tol.confusion() * 1e5,
                _ => false,
            };
        if first.is_equal(last, tol) || same_vertex {
            ring.pop();
            anchors.pop();
        } else if let Some(geometry) = model.geometry().surface(surface) {
            // A ring that winds one periodic direction of a doubly-periodic
            // surface — a diagonal loop on a torus. The folded walk ends a
            // whole period from where it began, and the face is the band
            // between the chain and its own translate one period over in the
            // *other* periodic direction, joined at the ends by columns that
            // lift to one 3D circle. The translate's anchors are the same 3D
            // points, and the joining columns' two copies lift identically,
            // so the weld closes them exactly as it closes a seam.
            use ogeom_geom::Surface as _;
            let ((ua, ub), (va, vb)) = geometry.domain();
            let du = last.x - first.x;
            let dv = last.y - first.y;
            let winds_u = geometry.is_periodic_u()
                && (du.abs() - (ub - ua)).abs() <= (ub - ua) * 1e-3
                && dv.abs() <= (vb - va).max(1.0) * 1e-3;
            let winds_v = geometry.is_periodic_u()
                && geometry.is_periodic_v()
                && (dv.abs() - (vb - va)).abs() <= (vb - va) * 1e-3
                && du.abs() <= (ub - ua).max(1.0) * 1e-3;
            // Where does a u-winding ring close against? On a doubly
            // periodic surface, its own translate one v-period over. On a
            // cone or sphere, the row where the surface collapses to a point
            // — the apex or the pole — which every u reaches: the closure
            // costs no area error because the row has none.
            let degenerate_row = |v: f64| -> bool {
                let (Ok(p), Ok(q), Ok(r)) = (
                    geometry.point_at(ua, v, tol),
                    geometry.point_at(f64::midpoint(ua, ub), v, tol),
                    geometry.point_at(ub, v, tol),
                ) else {
                    return false;
                };
                p.distance(q) <= tol.confusion() * 10.0 && p.distance(r) <= tol.confusion() * 10.0
            };
            let target_v = if winds_u && !geometry.is_periodic_v() {
                // The nearer degenerate row, if either end has one.
                let mid_v = f64::midpoint(first.y, last.y);
                if degenerate_row(va) && (mid_v - va).abs() <= (mid_v - vb).abs() {
                    Some(va)
                } else if degenerate_row(vb) {
                    Some(vb)
                } else if degenerate_row(va) {
                    Some(va)
                } else {
                    None
                }
            } else {
                None
            };
            let column_steps = 8;
            if let Some(v_apex) = target_v {
                // Down the seam column to the apex row, across it, and back
                // up: the row has no length in space, so the closure adds no
                // area and its lifted points weld to the one apex.
                let row_steps = ring.len().max(8);
                for k in 1..=column_steps {
                    let f = f64::from(k) / f64::from(column_steps);
                    ring.push(Point2::new(last.x, last.y + (v_apex - last.y) * f));
                    anchors.push(None);
                }
                for k in 1..row_steps {
                    #[allow(clippy::cast_precision_loss)]
                    let f = k as f64 / row_steps as f64;
                    ring.push(Point2::new(last.x + (first.x - last.x) * f, v_apex));
                    anchors.push(None);
                }
                for k in 0..column_steps {
                    let f = f64::from(column_steps - k) / f64::from(column_steps);
                    ring.push(Point2::new(first.x, first.y + (v_apex - first.y) * f));
                    anchors.push(None);
                }
            } else if (winds_u && geometry.is_periodic_v()) || winds_v {
                let shift = if winds_u {
                    Point2::new(0.0, -(vb - va))
                } else {
                    Point2::new(-(ub - ua), 0.0)
                };
                let chain: Vec<Point2> = ring.clone();
                let chain_anchors = anchors.clone();
                // Down from the chain's end to its translate's end.
                for k in 1..=column_steps {
                    let f = f64::from(k) / f64::from(column_steps);
                    ring.push(Point2::new(last.x + shift.x * f, last.y + shift.y * f));
                    anchors.push(None);
                }
                // The translate, walked back.
                for (p, a) in chain.iter().rev().zip(chain_anchors.iter().rev()).skip(1) {
                    ring.push(Point2::new(p.x + shift.x, p.y + shift.y));
                    anchors.push(*a);
                }
                // Up from the translate's start back to the chain's start,
                // stopping one step short of closing.
                for k in 1..column_steps {
                    let f = f64::from(column_steps - k) / f64::from(column_steps);
                    ring.push(Point2::new(first.x + shift.x * f, first.y + shift.y * f));
                    anchors.push(None);
                }
            }
        }
    }
    Ok((ring, anchors, met, folds, ties))
}

/// Parameters at which to sample an edge, taken from its 3D curve.
fn sample_parameters(
    model: &Model,
    data: &ogeom_topo::EdgeData,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<f64>, bool)>> {
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(None);
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let line = discretize(geometry, *range, deflection, tol)?;
    Ok(Some((line.parameters, line.deflection_met)))
}

/// Map parameters on the 3D curve onto the pcurve's own range.
///
/// The two are parameterized over their own intervals; `same_parameter` means
/// they agree *proportionally*, which is what this converts.
fn map_to_pcurve(
    parameters: &[f64],
    data: &ogeom_topo::EdgeData,
    pcurve_range: (f64, f64),
) -> Vec<f64> {
    let Some(EdgeRepr::Curve3d { range, .. }) = data.curve3d() else {
        return parameters.to_vec();
    };
    let (ca, cb) = *range;
    let (pa, pb) = pcurve_range;
    if (cb - ca).abs() <= f64::MIN_POSITIVE {
        return parameters.to_vec();
    }
    parameters
        .iter()
        .map(|u| pa + (pb - pa) * (u - ca) / (cb - ca))
        .collect()
}

/// The edge of a surface's domain, as a boundary ring.
///
/// Refined, not just the four corners. A triangulation only ever connects the
/// points it is given, so a boundary named by its corners alone forces long
/// triangles reaching right across the domain to find one — on a sphere, a
/// sliver from the equator to the pole. The interior refinement cannot fix that;
/// the missing points are on the boundary.
fn domain_ring(surface: &SurfaceGeometry, deflection: Deflection, tol: Tolerances) -> Vec<Point2> {
    let ((ua, ub), (va, vb)) = surface.domain();
    if ![ua, ub, va, vb].iter().all(|x| x.is_finite()) {
        return Vec::new();
    }

    let along_u = |v: f64| {
        refine_direction(ua, ub, deflection.chord, |a, b| {
            cell_error(surface, (a, v), (b, v), deflection, tol)
        })
    };
    let along_v = |u: f64| {
        refine_direction(va, vb, deflection.chord, |a, b| {
            cell_error(surface, (u, a), (u, b), deflection, tol)
        })
    };

    // Counter-clockwise around the rectangle. Each side drops its final point,
    // which the next side contributes: a repeated vertex would be a
    // zero-length boundary edge, and a constraint of zero length is not one.
    let (bottom, top) = (along_u(va), along_u(vb));
    let (left, right) = (along_v(ua), along_v(ub));
    let mut ring = Vec::new();
    ring.extend(
        bottom[..bottom.len() - 1]
            .iter()
            .map(|u| Point2::new(*u, va)),
    );
    ring.extend(right[..right.len() - 1].iter().map(|v| Point2::new(ub, *v)));
    ring.extend(top[1..].iter().rev().map(|u| Point2::new(*u, vb)));
    ring.extend(left[1..].iter().rev().map(|v| Point2::new(ua, *v)));
    ring
}

/// A triangulation still in parameter space, before it is lifted onto the
/// surface.
struct PlanarMesh {
    /// The `(u, v)` of each vertex.
    parameters: Vec<(f64, f64)>,
    /// Triangles as indices into `parameters`.
    triangles: Vec<[u32; 3]>,
    /// The boundary alone did not triangulate to the count a boundary that
    /// encloses a region gives: it crosses itself somewhere.
    crossed: bool,
}

/// Triangulate a region in parameter space, given its boundary rings.
///
/// The first ring is the outer boundary; the rest are holes.
fn triangulate_region(
    rings: &[Vec<Point2>],
    surface: &SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<PlanarMesh> {
    // A chart mangled enough — trims fitted through millimetres of boundary
    // error — can drive the triangulation library past its own asserts, and
    // a panic in a dependency is a crash in every consumer. The rings are
    // screened for what provably breaks it, and whatever still slips
    // through is caught at this boundary and spoken as the refusal it is:
    // the kernel's contract is refusal by name, never a crash on bad input.
    for ring in rings {
        let mut extent = 0.0_f64;
        for p in ring {
            if !p.x.is_finite() || !p.y.is_finite() {
                ogeom_bail!(
                    NotDone,
                    "a boundary ring carries a non-finite chart coordinate; \
                     the face's trim does not describe a region"
                );
            }
            extent = extent.max(p.x.abs()).max(p.y.abs());
        }
        if extent > 1e12 {
            ogeom_bail!(
                NotDone,
                "a boundary ring reaches {extent:.1e} in the chart; a trim \
                 that far out describes no face"
            );
        }
    }
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        triangulate_region_inner(rings, surface, deflection, tol)
    })) {
        Ok(result) => result,
        Err(_) => ogeom_bail!(
            NotDone,
            "the face's boundary broke the triangulation; the chart is too \
             degenerate to mesh"
        ),
    }
}

/// The chart stretched per axis to the surface's own metric, so that the
/// triangulation — Delaunay in the chart — sees distances as space does.
///
/// A fitted strip a fifth of a millimetre wide and a centimetre long may
/// carry `u` over a fiftieth of a unit and `v` over one: in the chart the
/// long way is the short way, sixty times over, and Delaunay, which
/// connects nearest neighbours *in the chart*, joins points along the
/// strip across several columns rather than to the row beside them. The
/// triangles it makes are slivers in space that lift folded — a flat
/// triangle spanning a bend the surface takes in between, its normal
/// pointing where none of its vertices' do — and the face shades as a
/// quilt of creases. Scaled by the mean tangent length each way, the
/// chart is the surface to first order, and Delaunay in it is Delaunay
/// on the surface.
#[derive(Debug, Clone, Copy)]
struct ChartScale {
    su: f64,
    sv: f64,
}

impl ChartScale {
    /// The mean tangent length each way over the rings' extent, the longer
    /// normalized to one; `(1, 1)` where the surface will not say.
    fn of(surface: &SurfaceGeometry, rings: &[Vec<Point2>], tol: Tolerances) -> Self {
        let (lo, hi) = chart_extent(rings);
        if !(lo.x.is_finite() && hi.x.is_finite() && lo.y.is_finite() && hi.y.is_finite()) {
            return Self { su: 1.0, sv: 1.0 };
        }
        let (mut du, mut dv, mut n) = (0.0_f64, 0.0_f64, 0usize);
        for i in 0..=2 {
            for j in 0..=2 {
                let u = lo.x + (hi.x - lo.x) * (0.25 + 0.25 * f64::from(i));
                let v = lo.y + (hi.y - lo.y) * (0.25 + 0.25 * f64::from(j));
                if let Ok((a, b)) = surface.d1_at(u, v, tol) {
                    du += a.magnitude();
                    dv += b.magnitude();
                    n += 1;
                }
            }
        }
        if n == 0 || !(du > 0.0 && dv > 0.0) || !du.is_finite() || !dv.is_finite() {
            return Self { su: 1.0, sv: 1.0 };
        }
        let longer = du.max(dv);
        Self {
            su: du / longer,
            sv: dv / longer,
        }
    }

    /// A chart point into the scaled chart.
    fn to(self, u: f64, v: f64) -> SpadePoint<f64> {
        SpadePoint::new(u * self.su, v * self.sv)
    }

    /// A scaled-chart point back into the chart.
    fn from(self, x: f64, y: f64) -> (f64, f64) {
        (x / self.su, y / self.sv)
    }
}

/// The rings' bounding box in the chart.
fn chart_extent(rings: &[Vec<Point2>]) -> (Point2, Point2) {
    rings.iter().flatten().fold(
        (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        |(lo, hi), p| {
            (
                Point2::new(lo.x.min(p.x), lo.y.min(p.y)),
                Point2::new(hi.x.max(p.x), hi.y.max(p.y)),
            )
        },
    )
}

fn triangulate_region_inner(
    rings: &[Vec<Point2>],
    surface: &SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<PlanarMesh> {
    let mut cdt: ConstrainedDelaunayTriangulation<SpadePoint<f64>> =
        ConstrainedDelaunayTriangulation::new();
    let sub = std::time::Instant::now();
    let mut refused_total = 0usize;

    // Everything the triangulation sees is in the scaled chart; the rings'
    // own chart points are kept by their scaled bits so a boundary vertex
    // comes back with exactly the parameters its anchor was keyed by.
    let scale = ChartScale::of(surface, rings, tol);
    let scaled: Vec<Vec<Point2>> = rings
        .iter()
        .map(|ring| {
            ring.iter()
                .map(|p| {
                    let q = scale.to(p.x, p.y);
                    Point2::new(q.x, q.y)
                })
                .collect()
        })
        .collect();
    let mut exact: std::collections::HashMap<(u64, u64), (f64, f64)> =
        std::collections::HashMap::new();

    // The boundary edges are constraints, so the triangulation respects the
    // trimming rather than spanning across a hole.
    for (ring, chart) in scaled.iter().zip(rings) {
        let mut ring_handles = Vec::with_capacity(ring.len());
        for (p, uv) in ring.iter().zip(chart) {
            exact.insert((p.x.to_bits(), p.y.to_bits()), (uv.x, uv.y));
            // A chart coordinate can come out subnormal-tiny — the sine of
            // a fold angle, the residue of an exact cancellation — and the
            // triangulation refuses what is, for every purpose, zero.
            let handle = cdt
                .insert(mitigate_underflow(SpadePoint::new(p.x, p.y)))
                .map_err(|e| ogeom_core::ogeom_err!(NotDone, "boundary insertion failed: {e}"))?;
            ring_handles.push(handle);
        }
        let mut refused = 0usize;
        let mut same = 0usize;
        for i in 0..ring_handles.len() {
            let (a, b) = (ring_handles[i], ring_handles[(i + 1) % ring_handles.len()]);
            if a == b {
                same += 1;
            } else if cdt.can_add_constraint(a, b) {
                cdt.add_constraint(a, b);
            } else {
                refused += 1;
                if *MESH_DEBUG_REFINE {
                    let n = ring.len();
                    let partners: Vec<usize> = (0..n)
                        .filter(|&j| j != i && j != (i + 1) % n && j != (i + n - 1) % n)
                        .filter(|&j| {
                            segments_cross(ring[i], ring[(i + 1) % n], ring[j], ring[(j + 1) % n])
                        })
                        .collect();
                    let len = ring[i].distance(ring[(i + 1) % n]);
                    eprintln!(
                        "REFUSED segment {i} of {n} (chart length {len:.3e}) crosses {partners:?}"
                    );
                    for j in [(i + n - 1) % n, i, (i + 1) % n, (i + 2) % n, (i + 3) % n] {
                        eprintln!("   ring[{j}] = ({:.12}, {:.12})", ring[j].x, ring[j].y);
                    }
                }
            }
        }
        if *MESH_DEBUG_REFINE && (refused > 0 || same > 0) {
            eprintln!(
                "CONSTRAINTS ring of {}: {refused} refused, {same} zero-length",
                ring.len()
            );
        }
        refused_total += refused;
    }

    // With the boundary in and nothing else, a boundary that encloses a
    // region triangulates to exactly `b + 2w - 4` triangles inside it, `b`
    // its distinct vertices and `w` its rings — the count any triangulation
    // of a polygon with holes has. One that crosses itself gives another
    // number: a segment refused as a constraint, a lobe wound the wrong
    // way. Asked here, before interior points bury the difference.
    let bands = RingBands::over(&scaled);
    // A segment the triangulation refused as a constraint crossed one
    // already there; that alone is the answer.
    let crossed = refused_total > 0
        || inside_by_parity(&cdt)
            .is_none_or(|inside| inside + 4 != cdt.num_vertices() + 2 * rings.len());
    if *MESH_DEBUG_REFINE && crossed {
        let points: usize = rings.iter().map(Vec::len).sum();
        eprintln!(
            "PARITY inside {:?} vertices {} points {points} rings {} inner faces {}",
            inside_by_parity(&cdt),
            cdt.num_vertices(),
            rings.len(),
            cdt.num_inner_faces()
        );
    }

    // Interior points where the surface bends away from the flat triangle. A
    // planar face needs none, which is why this is driven by measured
    // deflection rather than by a fixed grid.
    let boundary_ms = sub.elapsed().as_secs_f64() * 1e3;
    let sub = std::time::Instant::now();
    add_interior_points(&mut cdt, rings, surface, deflection, scale, tol)?;
    let interior_ms = sub.elapsed().as_secs_f64() * 1e3;
    let interior_points = cdt.num_vertices();
    let sub = std::time::Instant::now();
    let mut rounds_run = 0usize;

    // The scale a degenerate chart triangle is measured against: the
    // region's own span, the longer way. Its *position* is not its size —
    // a face on a cylinder whose axis point sits half a metre away has
    // `v` near −500 000 and a span of twenty, and a scale taken from where
    // the ring sits rather than how far it reaches would call every honest
    // cell a hair.
    let (lo, hi) = chart_extent(&scaled);
    let extent = (hi.x - lo.x).max(hi.y - lo.y);
    let degenerate_area = extent.max(1.0).powi(2) * 1e-12;

    // The grid rows guarantee the deflection along their own lines, but a
    // hole in the face punches a gap through a row, and where the surface is
    // flat in one direction — a cylinder along its axis — there may be no
    // other row for the mesher to reach. The band around the hole then fans
    // from the rim to the far side of the gap, in triangles that sag through
    // the solid by far more than the deflection while every one of their
    // vertices sits exactly on the surface. The boolean caught this as a
    // fused solid whose faces all had the right area and the wrong volume.
    //
    // The repair measures the truth: any kept triangle whose midpoints sag
    // beyond the chord gets its centre inserted, and the loop runs until the
    // mesh is honest or the cap says the surface is being unreasonable.
    // The rings do not change while the mesh is refined, so the containment
    // test they answer is indexed once and reused by every round below and by
    // the output pass.
    if !matches!(surface.kind(), ogeom_geom::SurfaceKind::Plane) {
        for _ in 0..REFINEMENT_ROUNDS {
            rounds_run += 1;
            let before = cdt.num_vertices();
            let mut worst: Vec<SpadePoint<f64>> = Vec::new();
            for triangle in cdt.inner_faces() {
                let vertices = triangle.vertices();
                let centre = triangle.center();
                let at = Point2::new(centre.x, centre.y);
                if !bands.holds(at) {
                    continue;
                }
                let corners: [(f64, f64); 3] = [
                    (vertices[0].position().x, vertices[0].position().y),
                    (vertices[1].position().x, vertices[1].position().y),
                    (vertices[2].position().x, vertices[2].position().y),
                ];
                // A chart-degenerate hair is dropped from the mesh, not
                // refined: its 3D chord can sag enormously, and feeding its
                // centre back in only breeds more hairs along the same line.
                let area = ((corners[1].0 - corners[0].0) * (corners[2].1 - corners[0].1)
                    - (corners[1].1 - corners[0].1) * (corners[2].0 - corners[0].0))
                    .abs()
                    / 2.0;
                if area < degenerate_area {
                    continue;
                }
                // The grid already bounds sag along rows and columns, and a
                // grid triangle's diagonal spanning one cell each way may
                // legitimately sag up to the sum — twice the chord — which
                // was the guarantee before this loop existed. The threshold
                // sits clear above that band so the repair fires only on the
                // fan triangles it exists for, which sag through a hole's
                // gap by tens of chords, and an honest grid — including a
                // perfectly symmetric one, whose mesh must stay symmetric —
                // is left untouched.
                let sagged = (0..3).any(|i| {
                    let (a, b) = (corners[i], corners[(i + 1) % 3]);
                    sag_between(surface, scale.from(a.0, a.1), scale.from(b.0, b.1), tol)
                        > deflection.chord * 3.0
                });
                if sagged {
                    worst.push(SpadePoint::new(centre.x, centre.y));
                }
            }
            if worst.is_empty() {
                break;
            }
            let mut inserted = 0usize;
            for point in worst {
                // A centre that lands on a vertex already there is not a new
                // point. A sliver whose apex sits on its own base — three
                // grid points on a diagonal, the middle one a rounding off
                // the line — has its centre at that apex to the last bits,
                // and inserting it breeds a hair a few ulps wide, whose
                // centre is the same point again: round after round, a
                // stack of hairs the degenerate filter then drops, and a
                // hole in the face where they were.
                if lands_on_the_mesh(&cdt, point, extent * 1e-9, (1.0, 1.0)) {
                    continue;
                }
                cdt.insert(mitigate_underflow(point)).map_err(|e| {
                    ogeom_core::ogeom_err!(NotDone, "refinement insertion failed: {e}")
                })?;
                inserted += 1;
            }
            if inserted == 0 {
                // Everything that sagged was a hair on a vertex; another
                // round would find the same hairs.
                break;
            }
            if *MESH_DEBUG_REFINE {
                eprintln!(
                    "ROUND {rounds_run}: +{} vertices",
                    cdt.num_vertices() - before
                );
            }
        }
    }

    let refine_ms = sub.elapsed().as_secs_f64() * 1e3;
    if *MESH_DEBUG_REFINE && boundary_ms + interior_ms + refine_ms > 50.0 {
        eprintln!(
            "SUB boundary {boundary_ms:.0}ms interior {interior_ms:.0}ms ({interior_points} verts) refine {refine_ms:.0}ms ({rounds_run} rounds, {} verts)",
            cdt.num_vertices()
        );
    }
    let mut parameters = Vec::new();
    let mut index_of = std::collections::HashMap::new();
    for (i, vertex) in cdt.vertices().enumerate() {
        let p = vertex.position();
        index_of.insert(vertex.fix(), i);
        parameters.push(
            exact
                .get(&(p.x.to_bits(), p.y.to_bits()))
                .copied()
                .unwrap_or_else(|| scale.from(p.x, p.y)),
        );
    }

    let mut triangles = Vec::new();
    let (mut dbg_total, mut dbg_outside, mut dbg_degenerate) = (0usize, 0usize, 0usize);
    for triangle in cdt.inner_faces() {
        dbg_total += 1;
        let vertices = triangle.vertices();
        let centre = triangle.center();
        // A constrained Delaunay covers the convex hull of its input, so
        // triangles outside the trimmed region — across a concavity, or inside
        // a hole — have to be discarded. Winding tells them apart.
        if !bands.holds(Point2::new(centre.x, centre.y)) {
            dbg_outside += 1;
            continue;
        }
        // A boundary run whose points differ by last-bit noise — a chart row
        // whose corner came off a different pcurve than its interior — lets
        // the triangulation weave a hair of a triangle along it: chart area
        // measured in ulps, a centroid *on* the boundary that even-odd
        // counting places wherever rounding falls, and a lifted sliver that
        // spans the run in one spurious stroke. It bounds nothing; drop it.
        let area = {
            let [a, b, c] = [
                vertices[0].position(),
                vertices[1].position(),
                vertices[2].position(),
            ];
            ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).abs() / 2.0
        };
        if area < degenerate_area {
            dbg_degenerate += 1;
            continue;
        }
        #[allow(clippy::cast_possible_truncation)]
        let indices = [
            index_of[&vertices[0].fix()] as u32,
            index_of[&vertices[1].fix()] as u32,
            index_of[&vertices[2].fix()] as u32,
        ];
        triangles.push(indices);
    }

    if *MESH_DEBUG_REFINE {
        eprintln!(
            "FILTER {dbg_total} triangles: {dbg_outside} outside, {dbg_degenerate} degenerate (area < {degenerate_area:.3e}), {} kept",
            triangles.len()
        );
    }
    // No triangles is not refused here: a boundary drawn coarsely enough
    // to cross itself can enclose nothing at all — an annulus narrower
    // than the sag of its rims' polygons, two arcs a hair apart — and the
    // caller's answer to a crossing is to draw the edges finer and ask
    // again. Empty counts as short; only a face still empty after that is
    // refused, by [`triangulate_with`].
    Ok(PlanarMesh {
        parameters,
        triangles,
        crossed,
    })
}

/// Add interior points wherever the surface deviates from flat by more than the
/// deflection allows.
///
/// The two parameter directions are refined *independently*, and that is not an
/// optimization — it is the difference between converging and not. A uniform
/// grid on a sphere puts as many meridians through the pole as through the
/// equator, so the triangles there become arbitrarily thin slivers in space.
/// The summed area of such a mesh does not approach the sphere's; it grows
/// without bound as the grid tightens. (Schwarz's lantern is the standard
/// example: an inscribed polyhedron whose area diverges under refinement.)
///
/// Refining each direction by its own measured sag fixes it at the source. Near
/// a pole the circle of latitude has almost no radius, so a chord right across
/// it sags by almost nothing and the direction stops subdividing after one or
/// two steps — while the meridian direction, whose curvature does not change,
/// keeps refining. The mesh degenerates into a fan, which is the right shape.
fn add_interior_points(
    cdt: &mut ConstrainedDelaunayTriangulation<SpadePoint<f64>>,
    rings: &[Vec<Point2>],
    surface: &SurfaceGeometry,
    deflection: Deflection,
    scale: ChartScale,
    tol: Tolerances,
) -> OgeomResult<()> {
    // A plane is flat everywhere; sampling it would add points that buy nothing.
    if matches!(surface.kind(), ogeom_geom::SurfaceKind::Plane) {
        return Ok(());
    }

    let bound = rings
        .iter()
        .flatten()
        .fold(ogeom_math::Aabb::EMPTY, |acc, p| {
            acc.with_point(Point::new(p.x, p.y, 0.0))
        });
    let (Some(low), Some(high)) = (bound.low(), bound.high()) else {
        return Ok(());
    };
    // Within this of a boundary vertex or segment is on it: a hair's width
    // at the scaled chart's scale, the same reach the repair pass keeps.
    let reach = ((high.x - low.x) * scale.su)
        .max((high.y - low.y) * scale.sv)
        .max(1.0)
        * 1e-9;

    // The v resolution has to hold everywhere the region reaches, so its sag is
    // the worst over a spread of u probes rather than the sag along one line.
    // A surface of revolution is the same at every u and a lofted one is not.
    #[allow(clippy::cast_precision_loss)]
    let probes: Vec<f64> = (0..=U_PROBES)
        .map(|i| low.x + (high.x - low.x) * i as f64 / U_PROBES as f64)
        .collect();
    let rows = refine_direction(low.y, high.y, deflection.chord, |a, b| {
        probes
            .iter()
            .map(|u| cell_error(surface, (*u, a), (*u, b), deflection, tol))
            .fold(0.0_f64, f64::max)
    });

    // Sag alone leaves a cylinder one row: it is straight along its axis, so
    // nothing along `v` ever sags. But the triangulation is Delaunay in the
    // chart, and a bore four hundred millimetres long with one row in the
    // middle hands it two-hundred-millimetre spans from each rim to that
    // row. Delaunay bridges those however it likes, and the repair below
    // fires only at three chords; a triangle a quarter turn wide on a two
    // millimetre bore sags less than that, so it stayed, and the bore drew
    // as a square between its holes. Cells are held to a bounded aspect
    // instead: rows close enough, measured in space through the surface,
    // that no triangle between two rows can reach across more than a few
    // columns. Rows are added, never removed, and spread evenly, so a grid
    // that was symmetric stays symmetric. The same the other way round.
    // How far apart two chart points are in space, where the surface says.
    let span = |p: (f64, f64), q: (f64, f64)| -> f64 {
        use ogeom_geom::Surface as _;
        match (
            surface.point_at(p.0, p.1, tol),
            surface.point_at(q.0, q.1, tol),
        ) {
            (Ok(a), Ok(b)) => a.distance(b),
            _ => 0.0,
        }
    };
    let rows = spread_to_aspect(
        rows,
        low.x,
        high.x,
        |v| {
            let columns = refine_direction(low.x, high.x, deflection.chord, |a, b| {
                cell_error(surface, (a, v), (b, v), deflection, tol)
            });
            columns.len()
        },
        |a, b| {
            probes
                .iter()
                .map(|&u| span((u, a), (u, b)))
                .fold(0.0_f64, f64::max)
        },
        |a, b, v| span((a, v), (b, v)),
    );

    if *MESH_DEBUG_REFINE {
        let v = f64::midpoint(low.y, high.y);
        let columns = refine_direction(low.x, high.x, deflection.chord, |a, b| {
            cell_error(surface, (a, v), (b, v), deflection, tol)
        });
        eprintln!(
            "GRID u [{:.3},{:.3}] v [{:.3},{:.3}]: {} rows after aspect, {} columns at the middle row, chord {}",
            low.x,
            high.x,
            low.y,
            high.y,
            rows.len(),
            columns.len(),
            deflection.chord
        );
    }
    for (row, &v) in rows
        .iter()
        .enumerate()
        .take(rows.len().saturating_sub(1))
        .skip(1)
    {
        // The row gap either side of this row, the finer of the two: the
        // chart scale a keep-out band is measured against along `v`.
        let dv = (v - rows[row - 1]).abs().min((rows[row + 1] - v).abs());
        // Each row gets its own u resolution, measured at that row.
        let columns = refine_direction(low.x, high.x, deflection.chord, |a, b| {
            cell_error(surface, (a, v), (b, v), deflection, tol)
        });
        // The same the other way round: a surface straight along `u` gets
        // two columns from sag, and a row a hundred millimetres wide would
        // bridge across the rows as badly as the bore bridged its columns.
        let columns = spread_to_aspect(
            columns,
            low.y,
            high.y,
            |_| rows.len(),
            |a, b| span((a, v), (b, v)),
            |a, b, u| span((u, a), (u, b)),
        );
        for (column, &u) in columns
            .iter()
            .enumerate()
            .take(columns.len().saturating_sub(1))
            .skip(1)
        {
            // Interior points only: the boundary is already constrained, and a
            // point landing just off a constraint would split it.
            if !inside_region(rings, Point2::new(u, v)) {
                continue;
            }
            // Inside, and not *on* the boundary: a grid point can fall
            // exactly on a ring segment that runs diagonally across the
            // chart — the midpoint of two grid corners the ring happens to
            // join — and even-odd counting calls it inside. Inserted, it
            // splits that constraint on this face alone, and the face
            // across the edge is drawn to the unsplit segment.
            let point = scale.to(u, v);
            if lands_on_the_mesh(cdt, point, reach, (1.0, 1.0)) {
                continue;
            }
            // Nor *near* it, measured in cells. A grid point a sliver's
            // width from a boundary chord makes a triangle with that
            // chord's two ends that is thin in the chart and, lifted, is
            // not thin at all: the chord cuts across the curvature by its
            // sag and the point sits on the surface, so the triangle
            // stands off the surface as a fin whose normal is tangent to
            // it and whose sign is whichever way the sliver leaned. A face
            // shades with a crease along every such chord. The point is
            // left out and the boundary's own row of triangles reaches
            // to the next grid line instead.
            let du = (u - columns[column - 1])
                .abs()
                .min((columns[column + 1] - u).abs());
            if du > 0.0
                && dv > 0.0
                && lands_on_the_mesh(
                    cdt,
                    point,
                    KEEP_OUT,
                    (1.0 / (du * scale.su), 1.0 / (dv * scale.sv)),
                )
            {
                continue;
            }
            cdt.insert(mitigate_underflow(point))
                .map_err(|e| ogeom_core::ogeom_err!(NotDone, "interior insertion failed: {e}"))?;
        }
    }
    Ok(())
}

/// Whether segments `a..b` and `c..d` cross properly: at a point interior
/// to both, neither touching the other's end.
fn segments_cross(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let orient =
        |p: Point2, q: Point2, r: Point2| (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
    let (o1, o2) = (orient(a, b, c), orient(a, b, d));
    let (o3, o4) = (orient(c, d, a), orient(c, d, b));
    o1 != 0.0
        && o2 != 0.0
        && o3 != 0.0
        && o4 != 0.0
        && (o1 > 0.0) != (o2 > 0.0)
        && (o3 > 0.0) != (o4 > 0.0)
}

/// How many of the triangulation's faces lie inside its constraints, told
/// by parity rather than by geometry.
///
/// Walking from a face on the convex hull, which is outside, every
/// constraint edge crossed flips inside for outside. Asked of the
/// triangulation of a boundary and nothing else, this is exact where the
/// even-odd test of a triangle's centre is not: a sliver face triangulates
/// to hairs whose centres sit on the boundary to the last bit, and which
/// side rounding puts them is a coin toss. `None` when the walk reaches a
/// face both ways with different answers — the constraints do not enclose
/// consistently, which is a crossing by another name.
fn inside_by_parity(cdt: &ConstrainedDelaunayTriangulation<SpadePoint<f64>>) -> Option<usize> {
    use std::collections::HashMap;
    let mut parity: HashMap<spade::handles::FixedFaceHandle<spade::handles::InnerTag>, bool> =
        HashMap::with_capacity(cdt.num_inner_faces());
    let mut queue = Vec::new();
    for hull in cdt.convex_hull() {
        // The hull edge's far side is the outer face; its near side is a
        // face of the triangulation, outside unless the hull edge itself
        // is a boundary.
        let Some(face) = hull.rev().face().as_inner() else {
            continue;
        };
        let inside = hull.is_constraint_edge();
        match parity.get(&face.fix()) {
            Some(&known) if known != inside => return None,
            Some(_) => {}
            None => {
                parity.insert(face.fix(), inside);
                queue.push((face.fix(), inside));
            }
        }
    }
    while let Some((face, inside)) = queue.pop() {
        for edge in cdt.face(face).adjacent_edges() {
            let Some(next) = edge.rev().face().as_inner() else {
                continue;
            };
            let next_inside = inside != edge.is_constraint_edge();
            match parity.get(&next.fix()) {
                Some(&known) if known != next_inside => return None,
                Some(_) => {}
                None => {
                    parity.insert(next.fix(), next_inside);
                    queue.push((next.fix(), next_inside));
                }
            }
        }
    }
    Some(parity.values().filter(|&&inside| inside).count())
}

/// Whether `point` sits within `reach` of a vertex the triangulation has,
/// or of one of its constraint edges.
///
/// Asked of whatever the point lands on — a vertex, an edge's two ends, a
/// face's three corners and whichever of its sides are constraints —
/// which is where anything that close must be. A point on a vertex is not
/// a new point; a point on a constraint would split it, and a boundary
/// split on one face only is a crack against the face across it.
fn lands_on_the_mesh(
    cdt: &ConstrainedDelaunayTriangulation<SpadePoint<f64>>,
    point: SpadePoint<f64>,
    reach: f64,
    scale: (f64, f64),
) -> bool {
    use spade::PositionInTriangulation as At;
    // Distances in a chart scaled per axis: `scale` is one over the local
    // grid step each way, so `reach` reads in cells, whatever the chart's
    // own units — one face's `u` runs over a fiftieth of its `v`.
    let scaled = |p: SpadePoint<f64>| ((p.x - point.x) * scale.0, (p.y - point.y) * scale.1);
    let near = |v: SpadePoint<f64>| {
        let (x, y) = scaled(v);
        x.hypot(y) <= reach
    };
    let along = |a: SpadePoint<f64>, b: SpadePoint<f64>| {
        // Distance to the segment `a..b`, the point at the origin.
        let (ax, ay) = scaled(a);
        let (bx, by) = scaled(b);
        let (dx, dy) = (bx - ax, by - ay);
        let len2 = dx * dx + dy * dy;
        let t = if len2 > 0.0 {
            ((-ax * dx - ay * dy) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (ax + t * dx).hypot(ay + t * dy) <= reach
    };
    match cdt.locate(point) {
        At::OnVertex(_) => true,
        At::OnEdge(edge) => {
            let edge = cdt.directed_edge(edge);
            edge.is_constraint_edge() || edge.vertices().iter().any(|v| near(v.position()))
        }
        At::OnFace(face) => {
            let face = cdt.face(face);
            face.vertices().iter().any(|v| near(v.position()))
                || face.adjacent_edges().iter().any(|e| {
                    e.is_constraint_edge() && along(e.from().position(), e.to().position())
                })
        }
        At::OutsideOfConvexHull(_) | At::NoTriangulation => false,
    }
}

/// How close to the boundary, in grid cells, an interior point may sit.
///
/// Closer than this the triangle between the point and a boundary chord
/// is a sliver in the chart and a fin in space; at this and beyond the
/// boundary's own row of triangles is at least this tall against the
/// chord, and lifts as a facet on the surface rather than off it.
const KEEP_OUT: f64 = 0.35;

/// How many column widths a grid cell may be tall before rows are added.
///
/// Delaunay in the chart connects nearest neighbours in the chart; held to
/// this aspect, a cell's nearest neighbours across a row gap are the same
/// columns, not columns several away, and the triangles between rows stay
/// as narrow as the columns are.
const CELL_ASPECT: f64 = 6.0;

/// Grid lines in one direction spread so that no cell is longer, in
/// space, than [`CELL_ASPECT`] times its width the other way — extra
/// lines added evenly between the ones sag chose.
///
/// Written for rows against columns and used both ways round. `lines`
/// are the parameters sag chose in this direction; `lo..hi` is the chart
/// range the other way; `crossings_at(t)` counts the lines the other way
/// at parameter `t` of this one; `length(a, b)` is the extent in space
/// between two lines of this direction; `width(a, b, t)` is the extent in
/// space between two parameters of the other direction, along this one's
/// line `t`.
fn spread_to_aspect(
    lines: Vec<f64>,
    lo: f64,
    hi: f64,
    crossings_at: impl Fn(f64) -> usize,
    length: impl Fn(f64, f64) -> f64,
    width: impl Fn(f64, f64, f64) -> f64,
) -> Vec<f64> {
    if lines.len() < 2 {
        return lines;
    }
    let mid = f64::midpoint(lines[0], lines[lines.len() - 1]);
    let crossings = crossings_at(mid);
    if crossings < 4 {
        // Flat the other way as well: a plane in all but name, and nothing
        // to hold an aspect against.
        return lines;
    }
    // One cell's width, not the whole range's: across a closed direction
    // the whole range comes back to its own start and measures nothing.
    #[allow(clippy::cast_precision_loss)]
    let step = (hi - lo) / (crossings - 1) as f64;
    let cell = width(lo, lo + step, mid);
    if cell.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return lines;
    }
    let mut out = Vec::with_capacity(lines.len());
    for pair in lines.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        out.push(a);
        let tall = length(a, b);
        let pieces = (tall / (CELL_ASPECT * cell)).ceil();
        if pieces.is_finite() && pieces > 1.0 {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let n = (pieces as usize).min(MAX_DIRECTION_STEPS);
            #[allow(clippy::cast_precision_loss)]
            for i in 1..n {
                out.push(a + (b - a) * i as f64 / n as f64);
            }
        }
    }
    out.push(lines[lines.len() - 1]);
    out
}

/// How many places across the domain the v resolution is measured at.
const U_PROBES: usize = 8;

/// How many rounds of sag-driven refinement a region may take.
///
/// Each round halves the worst sag roughly; a surface not honest after this
/// many is degenerate, and the cap makes that a coarse mesh rather than an
/// exhausted allocator.
const REFINEMENT_ROUNDS: usize = 12;

/// The most subdivisions one parameter direction may take.
///
/// A surface that has not converged by 512 has a singularity, not a resolution
/// problem, and the ceiling is what makes that a coarse mesh rather than an
/// exhausted allocator.
const MAX_DIRECTION_STEPS: usize = 512;

/// Subdivide `[lo, hi]` until no sub-interval sags further than `chord`.
///
/// The same adaptive bisection [`discretize`] uses on a curve, applied to a
/// line through parameter space. Returns the parameters in increasing order,
/// endpoints included.
fn refine_direction<F: Fn(f64, f64) -> f64>(lo: f64, hi: f64, chord: f64, sag: F) -> Vec<f64> {
    let mut values = vec![lo, f64::midpoint(lo, hi), hi];
    // A cursor rather than a rescan. Splitting an interval cannot change
    // whether an *earlier* one sags — the earlier one's endpoints do not move
    // — so restarting the search at zero re-measures intervals already known
    // to be good, and re-measuring is what costs: each measurement here is
    // several `sag_between` calls and each of those is three surface
    // evaluations. Reaching n points that way costs on the order of n²
    // measurements; walking forward costs n, and splits in the same
    // left-to-right order, so the values come out identical — including where
    // the step cap truncates them.
    let mut i = 0;
    while i + 1 < values.len() && values.len() < MAX_DIRECTION_STEPS {
        if sag(values[i], values[i + 1]) <= chord {
            i += 1;
            continue;
        }
        let mid = f64::midpoint(values[i], values[i + 1]);
        // A split that does not divide the interval means the parameters have
        // reached the resolution of f64, and refining further would loop.
        if mid <= values[i] || mid >= values[i + 1] {
            break;
        }
        values.insert(i + 1, mid);
    }
    values
}

/// How far the surface departs from the chord joining two parameter points.
///
/// Measured in space, which is the only place the number means anything: the
/// same step in `u` covers a metre at a sphere's equator and a millimetre near
/// its pole.
/// How far a grid cell's edge is from honest, as a sag: the chord sag
/// itself, or the normal's turn across it scaled so that a turn of the
/// angular deflection weighs the same as a sag of the chord — whichever
/// is worse.
///
/// The chord alone is what the boundary's edges are *not* drawn to: a
/// curve is discretized to both deflections, so a bore's rims come out
/// round at the angular limit while columns held to the chord alone come
/// out a polygon of far fewer sides, and the bore changes shape a chord's
/// length in from each rim. The interior is held to the same two limits
/// the boundary is.
fn cell_error(
    surface: &SurfaceGeometry,
    from: (f64, f64),
    to: (f64, f64),
    deflection: Deflection,
    tol: Tolerances,
) -> f64 {
    let sag = sag_between(surface, from, to, tol);
    let turn = match (
        surface.normal_at(from.0, from.1, tol),
        surface.normal_at(to.0, to.1, tol),
    ) {
        (Ok(a), Ok(b)) => a.angle(b),
        // A pole or an apex has no normal to compare; the sag still governs.
        _ => 0.0,
    };
    sag.max(turn / deflection.angular * deflection.chord)
}

fn sag_between(
    surface: &SurfaceGeometry,
    from: (f64, f64),
    to: (f64, f64),
    tol: Tolerances,
) -> f64 {
    let mid = (f64::midpoint(from.0, to.0), f64::midpoint(from.1, to.1));
    let (Ok(a), Ok(b), Ok(m)) = (
        surface.point_at(from.0, from.1, tol),
        surface.point_at(to.0, to.1, tol),
        surface.point_at(mid.0, mid.1, tol),
    ) else {
        // Off the surface's domain; nothing to refine towards.
        return 0.0;
    };
    ogeom_math::Axis::through(a, b, tol).map_or_else(|_| a.distance(m), |axis| axis.distance_to(m))
}

/// Whether a point is inside the region the rings bound.
fn inside_region(rings: &[Vec<Point2>], p: Point2) -> bool {
    inside_boundary_with::<Exact>(rings, p)
}

/// Ray-crossing count for one ring, decided by orientation rather than by
/// arithmetic.
///
/// A horizontal ray in `+u`. The half-open `y` comparison counts a vertex lying
/// exactly on the ray once rather than twice or not at all; which side of the
/// edge the point falls on is then an `orient2d` sign.
///
/// Deliberately *not* "solve for where the edge crosses the ray, then compare".
/// That form divides by the edge's `y` extent, which is near zero for a nearly
/// horizontal edge, and subtracts two nearly equal numbers to compare — so for a
/// point close to the boundary it can answer either way. Here there is no
/// division and the comparison is a determinant's sign, which an exact predicate
/// gets right at any separation.
fn crosses_odd_times<P: Predicates>(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    let n = ring.len();
    let at = |q: Point2| [q.x, q.y];
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        if (a.y > p.y) == (b.y > p.y) {
            continue;
        }
        // The edge crosses the ray's line. Whether it crosses the ray *itself*
        // — to the right of `p` — is which side of the directed edge `p` is on,
        // read the right way round for the edge's direction in `y`.
        let side = P::orient2d(at(a), at(b), [p.x, p.y]);
        let rightwards = if b.y > a.y {
            side == ogeom_core::Sign::Positive
        } else {
            side == ogeom_core::Sign::Negative
        };
        if rightwards {
            inside = !inside;
        }
    }
    inside
}

/// Whether the mesh debug dump is on, read once.
///
/// `env::var` takes a process-wide lock and allocates its answer, and this was
/// asked once per face — on an imported assembly, once per face of every part.
/// Whether to report, per shape, how many faces were drawn again finer.
static MESH_DEBUG_REFINE: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("OGEOM_MESH_DEBUG_REFINE").is_ok());

static MESH_DEBUG: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("OGEOM_MESH_DEBUG").is_ok());

/// The ring edges that can cross a horizontal ray, bucketed by height.
///
/// [`crosses_odd_times`] walks every edge of every ring for each point it is
/// asked about. That is fine for a handful of queries and ruinous for the
/// refinement loop, which asks once per triangle per round while the rings
/// themselves never change: a face with 544 boundary points and 41 000
/// triangles pays a quarter of a billion edge visits, nearly all on edges
/// nowhere near the point.
///
/// An edge can only straddle a ray at height `y` if `y` lies within the edge's
/// own `y` span, so bucketing edges by that span and querying one bucket tests
/// a conservative superset of the edges that could contribute. **The answer is
/// therefore identical** — the same straddle test and the same exact predicate
/// decide each candidate; the index only declines to visit edges that could
/// not have counted.
///
/// Parity is taken over all rings at once, which is what
/// [`inside_boundary_with`] computes as an exclusive-or of per-ring parities:
/// the two agree because the parity of the total crossing count is the sum of
/// the rings' parities.
struct RingBands {
    /// Every ring's edges, flattened.
    edges: Vec<(Point2, Point2)>,
    /// Edge indices per band, low `y` first.
    bands: Vec<Vec<u32>>,
    low: f64,
    high: f64,
    /// Band height. Zero when every point shares one `y`, which leaves a
    /// single band holding everything.
    step: f64,
}

impl RingBands {
    /// Index the rings. Cheap enough to build per face and paid back by the
    /// first few hundred queries.
    fn over(rings: &[Vec<Point2>]) -> Self {
        let mut edges = Vec::new();
        for ring in rings {
            for i in 0..ring.len() {
                edges.push((ring[i], ring[(i + 1) % ring.len()]));
            }
        }
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for (a, b) in &edges {
            low = low.min(a.y).min(b.y);
            high = high.max(a.y).max(b.y);
        }
        if edges.is_empty() || !low.is_finite() || !high.is_finite() {
            return Self {
                edges,
                bands: Vec::new(),
                low: 0.0,
                high: 0.0,
                step: 0.0,
            };
        }
        // About four edges to a band: enough to keep the per-query walk short
        // without spreading a long edge across a table of mostly empty bands.
        let count = edges.len().div_ceil(4).clamp(1, 4096);
        #[allow(
            clippy::cast_precision_loss,
            reason = "a band count, far below the integers f64 represents exactly"
        )]
        let step = (high - low) / count as f64;
        let mut bands: Vec<Vec<u32>> = vec![Vec::new(); count];
        for (i, (a, b)) in edges.iter().enumerate() {
            let (lo, hi) = (a.y.min(b.y), a.y.max(b.y));
            let first = Self::band_of(lo, low, step, count);
            let last = Self::band_of(hi, low, step, count);
            for band in &mut bands[first..=last] {
                #[allow(
                    clippy::cast_possible_truncation,
                    reason = "an edge index, bounded by the ring lengths"
                )]
                band.push(i as u32);
            }
        }
        Self {
            edges,
            bands,
            low,
            high,
            step,
        }
    }

    /// Which band a height falls in, clamped to the table.
    fn band_of(y: f64, low: f64, step: f64, count: usize) -> usize {
        if step <= 0.0 {
            return 0;
        }
        let raw = (y - low) / step;
        if raw <= 0.0 {
            return 0;
        }
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "clamped to the band count on the next line"
        )]
        let index = raw as usize;
        index.min(count - 1)
    }

    /// Whether the point lies inside the region the rings bound.
    fn holds(&self, p: Point2) -> bool {
        if self.bands.is_empty() || p.y < self.low || p.y > self.high {
            return false;
        }
        let band = Self::band_of(p.y, self.low, self.step, self.bands.len());
        let at = |q: Point2| [q.x, q.y];
        let mut inside = false;
        for &i in &self.bands[band] {
            let (a, b) = self.edges[i as usize];
            if (a.y > p.y) == (b.y > p.y) {
                continue;
            }
            let side = Exact::orient2d(at(a), at(b), [p.x, p.y]);
            let rightwards = if b.y > a.y {
                side == ogeom_core::Sign::Positive
            } else {
                side == ogeom_core::Sign::Negative
            };
            if rightwards {
                inside = !inside;
            }
        }
        inside
    }
}

/// The unit normal of a triangle, or `None` if it is degenerate.
#[must_use]
pub fn triangle_normal(a: Point, b: Point, c: Point, tol: Tolerances) -> Option<Direction> {
    Direction::from_cross(b - a, c - a, tol).ok()
}

/// Discretize an edge into a polyline in space, for display or coarse queries.
///
/// # Errors
///
/// As [`discretize`].
pub fn polyline_of_edge(
    model: &Model,
    edge: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "expected an edge");
    }
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(Vec::new());
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let placement = edge.transform(model.datums())?;
    let line = discretize(geometry, *range, deflection, tol)?;
    let mut points: Vec<Point> = line.points.iter().map(|p| placement.apply(*p)).collect();
    if edge.orientation() == Orientation::Reversed {
        points.reverse();
    }
    Ok(points)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn a_wild_chart_refuses_instead_of_crashing() {
        let plane: SurfaceGeometry = ogeom_geom::PlaneSurface::over(
            ogeom_math::Plane::through(
                ogeom_math::Point::new(0.0, 0.0, 0.0),
                ogeom_math::Direction::Z,
            ),
            (-10.0, 10.0),
            (-10.0, 10.0),
        )
        .unwrap()
        .into();
        let tol = Tolerances::millimetres();
        // Coordinates wilder than any face: refused by the screen.
        let wild = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(1e15, 0.0),
            Point2::new(0.0, 1e15),
        ]];
        let Err(err) = triangulate_region(&wild, &plane, Deflection::default(), tol) else {
            panic!("a wild chart must refuse");
        };
        assert!(err.to_string().contains("describes no face"), "{err}");
        // A non-finite coordinate: refused by name, not fed to the library.
        let nan = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(f64::NAN, 1.0),
            Point2::new(1.0, 1.0),
        ]];
        let Err(err) = triangulate_region(&nan, &plane, Deflection::default(), tol) else {
            panic!("a NaN chart must refuse");
        };
        assert!(err.to_string().contains("non-finite"), "{err}");
    }

    use approx::assert_relative_eq;
    use ogeom_algo::make_box;
    use ogeom_math::Frame;
    use ogeom_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    /// `refine_direction` as it was written: rescan from zero after every
    /// split. Kept here as the reference the cursor form is held to.
    fn refine_by_rescan<F: Fn(f64, f64) -> f64>(lo: f64, hi: f64, chord: f64, sag: F) -> Vec<f64> {
        let mut values = vec![lo, f64::midpoint(lo, hi), hi];
        while values.len() < MAX_DIRECTION_STEPS {
            let Some(i) = (0..values.len() - 1).find(|&i| sag(values[i], values[i + 1]) > chord)
            else {
                break;
            };
            let mid = f64::midpoint(values[i], values[i + 1]);
            if mid <= values[i] || mid >= values[i + 1] {
                break;
            }
            values.insert(i + 1, mid);
        }
        values
    }

    #[test]
    fn walking_forward_splits_where_rescanning_did() {
        // The cursor is only sound because splitting an interval cannot change
        // whether an earlier one sags. Held to the old form's output exactly,
        // over sag profiles that bite in different places: flat, steep at one
        // end, periodic, and one savage enough to reach the step cap.
        /// A named sag profile to hold both forms to.
        type Profile = (&'static str, Box<dyn Fn(f64, f64) -> f64>);
        let cases: Vec<Profile> = vec![
            ("flat", Box::new(|_a: f64, _b: f64| 0.0)),
            ("width", Box::new(|a: f64, b: f64| (b - a).abs())),
            (
                "steep at the low end",
                Box::new(|a: f64, b: f64| (b - a).abs() / a.abs().max(1e-3)),
            ),
            (
                "periodic",
                Box::new(|a: f64, b: f64| (b - a).abs() * (a * 12.0).sin().abs()),
            ),
            (
                "beyond the cap",
                Box::new(|a: f64, b: f64| (b - a).abs() * 1e6),
            ),
        ];
        for (name, sag) in cases {
            for chord in [1.0, 0.1, 0.01, 1e-3] {
                let walked = refine_direction(0.0, 1.0, chord, &sag);
                let rescanned = refine_by_rescan(0.0, 1.0, chord, &sag);
                assert_eq!(
                    walked, rescanned,
                    "{name} at chord {chord}: the cursor split somewhere the rescan did not"
                );
            }
        }
    }

    #[test]
    fn banded_containment_answers_what_the_full_scan_answers() {
        // The index may only decline to visit edges that could not have
        // counted. Held to the unindexed predicate over a ring with a hole,
        // on a grid that straddles both boundaries and the vertices themselves.
        let outer: Vec<Point2> = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 3.0),
            Point2::new(2.0, 1.5),
            Point2::new(0.0, 3.0),
        ];
        let hole: Vec<Point2> = vec![
            Point2::new(1.0, 0.5),
            Point2::new(1.0, 1.0),
            Point2::new(1.5, 1.0),
            Point2::new(1.5, 0.5),
        ];
        let rings = vec![outer, hole];
        let bands = RingBands::over(&rings);
        for i in 0..=80 {
            for j in 0..=60 {
                #[allow(clippy::cast_precision_loss)]
                let p = Point2::new(f64::from(i) * 0.05 - 0.1, f64::from(j) * 0.05 - 0.1);
                assert_eq!(
                    bands.holds(p),
                    inside_region(&rings, p),
                    "the index disagreed with the full scan at {p:?}"
                );
            }
        }
    }

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    /// A vertex at a cone's apex or a sphere's pole carries the normal the
    /// surface tends to there along its own column, not nothing.
    #[test]
    fn apex_and_pole_vertices_carry_the_limit_normal() {
        let mut model = Model::new();
        let cone = ogeom_algo::make_cone(&mut model, Frame::WORLD, 2.0, 0.0, 3.0, T).unwrap();
        let sphere = ogeom_algo::make_sphere(&mut model, Frame::WORLD, 1.5, T).unwrap();
        for (shape, what) in [(&cone.shape, "cone"), (&sphere.shape, "sphere")] {
            for face in explore_unique(&model, shape, ShapeType::Face).unwrap() {
                let mesh = triangulate_face(&model, &face, fine(), T).unwrap();
                for (i, n) in mesh.normals.iter().enumerate() {
                    assert!(
                        (n.magnitude() - 1.0).abs() < 1e-9,
                        "{what} vertex {i} at {:?} has normal {n:?}",
                        mesh.positions[i]
                    );
                }
            }
        }
        // At the apex the limit normal along a ruling is that ruling's
        // normal: on a cone of half-angle atan(2/3) it leans out by that
        // much from the axis, the same as every other normal in its column.
        let faces = explore_unique(&model, &cone.shape, ShapeType::Face).unwrap();
        let lean = (2.0_f64 / 3.0).atan();
        for face in &faces {
            let mesh = triangulate_face(&model, face, fine(), T).unwrap();
            for (i, p) in mesh.positions.iter().enumerate() {
                if p.distance(Point::new(0.0, 0.0, 3.0)) < 1e-9 {
                    let n = mesh.normals[i];
                    let from_axis = n.z.abs().acos();
                    assert!(
                        ((std::f64::consts::FRAC_PI_2 - from_axis) - lean).abs() < 1e-6,
                        "apex normal {n:?} leans {from_axis} from the axis"
                    );
                }
            }
        }
    }

    #[test]
    fn a_box_face_triangulates_into_two_triangles() {
        // A planar square needs no interior points at all, which is what makes
        // deflection-driven refinement worth having: a fixed grid would add
        // dozens that buy nothing.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();

        for face in &faces {
            let mesh = triangulate_face(&model, face, fine(), T).unwrap();
            assert_eq!(mesh.triangle_count(), 2, "a rectangle is two triangles");
            assert_eq!(mesh.vertex_count(), 4);
            assert!(mesh.deflection_met);
        }
    }

    #[test]
    fn a_boxs_mesh_is_closed_and_reports_the_right_volume() {
        // The end-to-end check: triangulate, weld, and ask the mesh what it
        // encloses. A volume that comes out negative would mean the faces are
        // wound inward; one that is wrong in magnitude would mean the
        // triangulation is not covering the boundary.
        let mut model = Model::new();
        let size = (2.0, 3.0, 4.0);
        let built = make_box(&mut model, Frame::WORLD, size, T).unwrap();
        let mesh = triangulate(&model, &built.shape, fine(), T).unwrap();

        assert_eq!(mesh.triangle_count(), 12, "six faces, two triangles each");
        assert_eq!(mesh.vertex_count(), 8, "welding merged the shared corners");
        assert!(
            mesh.is_closed(),
            "every triangle edge should be shared by two"
        );
        assert_relative_eq!(mesh.volume(), size.0 * size.1 * size.2, epsilon = 1e-9);
        assert_relative_eq!(
            mesh.area(),
            2.0 * (size.0 * size.1 + size.1 * size.2 + size.2 * size.0),
            epsilon = 1e-9
        );
    }

    #[test]
    fn welding_is_what_closes_the_mesh() {
        // Without it each face brings its own copy of every boundary vertex, so
        // no triangle edge is shared and the surface is a pile of loose squares.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();

        let mut loose = Triangulation::new();
        for face in ogeom_topo::explore(
            &model,
            &built.shape,
            ogeom_topo::Filter::OfType(ShapeType::Face),
        )
        .unwrap()
        {
            loose.append(&triangulate_face(&model, &face, fine(), T).unwrap());
        }
        assert_eq!(loose.vertex_count(), 24, "four corners per face, unmerged");
        assert!(!loose.is_closed());

        let welded = loose.welded(T);
        assert_eq!(welded.vertex_count(), 8);
        assert!(welded.is_closed());
    }

    #[test]
    fn a_reversed_face_presents_the_other_side() {
        // A renderer or a volume computation that ignored orientation would
        // have the solid inside out, and nothing about the positions says so.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();

        let forward = triangulate_face(&model, &face, fine(), T).unwrap();
        let backward = triangulate_face(&model, &face.reversed(), fine(), T).unwrap();

        assert_eq!(forward.triangle_count(), backward.triangle_count());
        for (a, b) in forward.normals.iter().zip(&backward.normals) {
            assert!(a.is_equal(-*b, T), "normals did not flip: {a:?} vs {b:?}");
        }
        // And the winding flipped with them, so the two agree.
        let winding = |m: &Triangulation, i: usize| {
            let [a, b, c] = m.triangles[i].map(|k| m.positions[k as usize]);
            (b - a).cross(c - a)
        };
        assert!(winding(&forward, 0).dot(winding(&backward, 0)) < 0.0);
    }

    #[test]
    fn every_triangle_vertex_lies_on_the_surface_it_came_from() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 1.0, 3.0), T).unwrap();

        for face in explore_unique(&model, &built.shape, ShapeType::Face).unwrap() {
            let mesh = triangulate_face(&model, &face, fine(), T).unwrap();
            let data = model.node(&face).unwrap().data().as_face().unwrap().clone();
            let surface = model.geometry().surface(data.surface).unwrap();
            for (position, (u, v)) in mesh.positions.iter().zip(&mesh.parameters) {
                let exact = surface.point_at(*u, *v, T).unwrap();
                assert!(
                    position.is_equal(exact, T),
                    "{position:?} is not on its surface"
                );
            }
        }
    }

    #[test]
    fn a_finer_deflection_never_gives_fewer_triangles() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let coarse = triangulate(&model, &built.shape, Deflection::default(), T).unwrap();
        let detailed = triangulate(&model, &built.shape, fine(), T).unwrap();
        assert!(detailed.triangle_count() >= coarse.triangle_count());
        // Both enclose the same volume, since the faces are flat.
        assert_relative_eq!(coarse.volume(), 1.0, epsilon = 1e-9);
        assert_relative_eq!(detailed.volume(), 1.0, epsilon = 1e-9);
    }

    #[test]
    fn a_sphere_gets_interior_points_and_converges_on_its_true_area() {
        // The whole point of measuring deflection in space rather than in
        // parameter space. A sphere's parameter rectangle is uniform; the
        // surface it maps to is not, and a fixed grid would be dense at the
        // poles and coarse at the equator.
        use ogeom_algo::make_natural_face;
        use ogeom_geom::SphereSurface;
        use ogeom_math::Sphere;

        let radius = 10.0;
        let exact = 4.0 * std::f64::consts::PI * radius * radius;
        let mut previous = 0.0;

        for chord in [1.0_f64, 0.25, 0.05] {
            let mut model = Model::new();
            let surface = SphereSurface::new(Sphere::new(Frame::WORLD, radius, T).unwrap());
            let face = make_natural_face(&mut model, surface.into()).unwrap().shape;
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            let mesh = triangulate_face(&model, &face, deflection, T).unwrap();

            assert!(
                mesh.triangle_count() > 2,
                "a curved face needs interior points, got {} triangles",
                mesh.triangle_count()
            );
            // Every triangle chord-cuts the sphere, so the area comes in under
            // the truth and climbs as the tolerance tightens.
            let area = mesh.area();
            assert!(area < exact, "a chord-cut area cannot exceed the surface's");
            assert!(
                area > previous,
                "tightening the chord from the previous step lost area: \
                 {area} after {previous}"
            );
            previous = area;
        }
        assert!(
            previous > exact * 0.99,
            "at a chord of 0.05 on a radius of 10 the area should be within a \
             percent, got {previous} against {exact}"
        );
    }

    #[test]
    fn every_sphere_vertex_sits_at_the_right_radius() {
        // Lifting through the surface is what makes the mesh curved at all; a
        // vertex left in parameter space, or lifted with the wrong parameters,
        // would land nowhere near the sphere.
        use ogeom_algo::make_natural_face;
        use ogeom_geom::SphereSurface;
        use ogeom_math::Sphere;

        let mut model = Model::new();
        let surface = SphereSurface::new(Sphere::new(Frame::WORLD, 3.0, T).unwrap());
        let face = make_natural_face(&mut model, surface.into()).unwrap().shape;
        let mesh = triangulate_face(&model, &face, Deflection::default(), T).unwrap();

        for p in &mesh.positions {
            assert_relative_eq!(p.to_vector().magnitude(), 3.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn a_curved_domains_boundary_is_refined_not_just_its_corners() {
        // The corners alone would leave the triangulation nothing to connect to
        // along a side, so it reaches right across the domain for one — and a
        // sliver from a sphere's equator to its pole makes the summed area
        // diverge under refinement rather than converge.
        use ogeom_geom::{PlaneSurface, SphereSurface};
        use ogeom_math::{Plane, Sphere};

        let sphere: SurfaceGeometry =
            SphereSurface::new(Sphere::new(Frame::WORLD, 10.0, T).unwrap()).into();
        let coarse = domain_ring(&sphere, Deflection::default(), T);
        let fine_ring = domain_ring(&sphere, fine(), T);
        assert!(coarse.len() > 4, "a sphere's domain edge is curved");
        assert!(
            fine_ring.len() > coarse.len(),
            "a tighter chord should place more boundary points"
        );

        // No side may repeat a corner: a zero-length constraint is not one.
        for w in fine_ring.windows(2) {
            assert!(!w[0].is_equal(w[1], T), "the ring repeats a point");
        }

        // A plane is flat, so its domain needs only what closes the rectangle.
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        assert!(domain_ring(&plane, fine(), T).len() >= 4);
    }

    #[test]
    fn winding_detects_points_inside_and_outside_a_ring() {
        let square = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        assert!(inside_region(
            std::slice::from_ref(&square),
            Point2::new(0.5, 0.5)
        ));
        assert!(!inside_region(
            std::slice::from_ref(&square),
            Point2::new(1.5, 0.5)
        ));
        assert!(!inside_region(
            std::slice::from_ref(&square),
            Point2::new(0.5, -0.5)
        ));

        // With a hole, the middle is outside again.
        let hole = vec![
            Point2::new(0.4, 0.4),
            Point2::new(0.6, 0.4),
            Point2::new(0.6, 0.6),
            Point2::new(0.4, 0.6),
        ];
        let with_hole = vec![square, hole];
        assert!(!inside_region(&with_hole, Point2::new(0.5, 0.5)));
        assert!(inside_region(&with_hole, Point2::new(0.2, 0.2)));
    }

    #[test]
    fn a_polyline_of_an_edge_runs_in_the_edges_direction() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let edge = explore_unique(&model, &built.shape, ShapeType::Edge).unwrap()[0].clone();

        let forward = polyline_of_edge(&model, &edge, fine(), T).unwrap();
        let backward = polyline_of_edge(&model, &edge.reversed(), fine(), T).unwrap();
        assert!(forward.len() >= 2);
        assert!(forward[0].is_equal(backward[backward.len() - 1], T));
        assert!(forward[forward.len() - 1].is_equal(backward[0], T));
    }

    #[test]
    fn triangulating_something_that_is_not_a_face_is_refused() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        assert!(triangulate_face(&model, &built.shape, fine(), T).is_err());

        let vertex = explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap()[0].clone();
        assert!(triangulate_face(&model, &vertex, fine(), T).is_err());
        assert!(polyline_of_edge(&model, &vertex, fine(), T).is_err());
    }

    #[test]
    fn a_triangles_normal_is_none_when_it_is_degenerate() {
        let a = Point::ORIGIN;
        let b = Point::new(1.0, 0.0, 0.0);
        assert!(triangle_normal(a, b, Point::new(0.0, 1.0, 0.0), T).is_some());
        assert!(triangle_normal(a, b, Point::new(2.0, 0.0, 0.0), T).is_none());
        assert!(triangle_normal(a, a, a, T).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod predicate_tests {
    use super::*;
    use ogeom_core::Fast;

    /// A triangle with one very long, very nearly diagonal edge.
    ///
    /// The configuration where a naive crossing test goes wrong. Subtracting
    /// the edge's ends from a query point close to the line cancels almost
    /// every significant digit, and what is left is rounding rather than
    /// geometry.
    fn sliver() -> Vec<Vec<Point2>> {
        vec![vec![
            Point2::new(0.5, 0.5),
            Point2::new(1000.0, 1000.0),
            Point2::new(1000.0, 0.5),
        ]]
    }

    #[test]
    fn the_two_implementations_are_a_real_choice_and_not_a_decoration() {
        // The point of the seam. In a band of near-degenerate queries the two
        // answer differently — and if they never did, routing through the trait
        // would be ceremony rather than a design.
        //
        // The query points straddle the long edge at a spacing far below what
        // the subtraction can resolve, which is exactly the case a mesh hits
        // when a face's boundary passes close to a triangulation vertex.
        let rings = sliver();
        let step = 2.0_f64.powi(-48);
        let mut disagreements = 0;
        for i in 0..256i32 {
            for j in 0..256i32 {
                let p = Point2::new(
                    250.0 + f64::from(i - 128) * step,
                    250.0 + f64::from(j - 128) * step,
                );
                if inside_boundary_with::<Exact>(&rings, p)
                    != inside_boundary_with::<Fast>(&rings, p)
                {
                    disagreements += 1;
                }
            }
        }
        assert!(
            disagreements > 0,
            "the exact and fast predicates never disagreed, so the seam is not \
             carrying anything"
        );
    }

    #[test]
    fn the_exact_predicate_is_the_one_that_is_right_near_the_edge() {
        // Disagreeing is not enough — the exact one has to be the *correct*
        // one, and that needs points whose side is known without asking either
        // implementation.
        //
        // Stepped in units of the last place rather than by a small distance.
        // A distance below the spacing of `f64` at 250 rounds away entirely,
        // leaving two points that are both exactly *on* the diagonal — where
        // either answer is defensible and the test would be asserting nothing.
        let rings = sliver();
        let base = 250.0_f64;
        let mut off = base;
        for k in 1..64 {
            off = off.next_up();
            assert_ne!(off, base, "step {k} did not move the point at all");
            // Inside this triangle is below the diagonal `y = x`: larger x.
            assert!(
                inside_boundary_with::<Exact>(&rings, Point2::new(off, base)),
                "a point {k} ulps below the diagonal was reported outside"
            );
            assert!(
                !inside_boundary_with::<Exact>(&rings, Point2::new(base, off)),
                "a point {k} ulps above the diagonal was reported inside"
            );
        }
    }

    #[test]
    fn the_exact_answer_is_the_one_the_geometry_supports() {
        // Points placed by construction, so the right answer is known without
        // asking either implementation.
        let square = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]];
        assert!(inside_boundary_with::<Exact>(
            &square,
            Point2::new(0.5, 0.5)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &square,
            Point2::new(1.5, 0.5)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &square,
            Point2::new(-0.5, 0.5)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &square,
            Point2::new(0.5, 1.5)
        ));

        // A vertex exactly on the sampling ray is counted once, not twice or
        // not at all — which is what the half-open comparison is for.
        let diamond = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, -1.0),
        ]];
        assert!(inside_boundary_with::<Exact>(
            &diamond,
            Point2::new(1.0, 0.0)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &diamond,
            Point2::new(3.0, 0.0)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &diamond,
            Point2::new(-1.0, 0.0)
        ));
    }

    /// A cylinder's grid has as many rows as its length needs, not as
    /// many as its sag asks for.
    #[test]
    fn rows_are_spread_until_no_cell_is_taller_than_its_aspect() {
        // Sixteen columns a millimetre wide over a region a hundred long:
        // sag left one interior row, the aspect wants cells six tall.
        let sagged = vec![0.0, 50.0, 100.0];
        let spread = spread_to_aspect(
            sagged.clone(),
            0.0,
            16.0,
            |_| 17,
            |a, b| (b - a).abs(),
            |a, b, _| (b - a).abs(),
        );
        assert_eq!(
            spread.len(),
            2 * 9 + 1,
            "each 50 mm half in nine 6 mm pieces: {spread:?}"
        );
        assert_eq!(spread[0], 0.0);
        assert_eq!(spread[9], 50.0, "the rows sag chose stay where they were");
        assert_eq!(spread[18], 100.0);
        assert!(spread.windows(2).all(|w| w[1] > w[0]));

        // Three columns is flat along `u` as well; nothing to hold to.
        let flat = spread_to_aspect(
            sagged.clone(),
            0.0,
            16.0,
            |_| 3,
            |a, b| (b - a).abs(),
            |a, b, _| (b - a).abs(),
        );
        assert_eq!(flat, sagged);

        // Cells already shorter than the aspect are left alone.
        let fine = vec![0.0, 4.0, 8.0];
        let same = spread_to_aspect(
            fine.clone(),
            0.0,
            16.0,
            |_| 17,
            |a, b| (b - a).abs(),
            |a, b, _| (b - a).abs(),
        );
        assert_eq!(same, fine);
    }

    /// A ring's slop-duplicates and spikes come off before it is
    /// triangulated.
    #[test]
    fn a_ring_is_cleaned_of_duplicates_and_spikes() {
        let p = |x: f64, y: f64| Point2::new(x, y);
        // A square whose closing point repeats its first a hair off, and
        // whose right side steps out to a point and straight back.
        let mut ring = vec![
            p(0.0, 0.0),
            p(10.0, 0.0),
            p(10.0, 5.0),
            p(10.2, 5.0),
            p(10.0, 5.0 + 1e-9),
            p(10.0, 10.0),
            p(0.0, 10.0),
            p(1e-9, 1e-9),
        ];
        let mut anchors = vec![None; ring.len()];
        let reach = chart_reach(&ring);
        merge_near_duplicates(&mut ring, &mut anchors, reach);
        assert_eq!(ring.len(), 7, "the closing duplicate is gone: {ring:?}");
        remove_spikes(&mut ring, &mut anchors, reach);
        assert_eq!(
            ring,
            vec![
                p(0.0, 0.0),
                p(10.0, 0.0),
                p(10.0, 5.0),
                p(10.0, 10.0),
                p(0.0, 10.0)
            ],
            "the spike is gone and the ring is the square with a point on one side"
        );
        assert_eq!(anchors.len(), ring.len());
    }

    #[test]
    fn a_hole_is_outside_however_either_ring_is_wound() {
        // Even-odd does not depend on the wires being wound consistently, which
        // is what makes it survive imported geometry.
        let outer = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 4.0),
        ];
        let hole: Vec<Point2> = vec![
            Point2::new(1.0, 1.0),
            Point2::new(3.0, 1.0),
            Point2::new(3.0, 3.0),
            Point2::new(1.0, 3.0),
        ];
        let backwards: Vec<Point2> = hole.iter().rev().copied().collect();
        for inner in [hole, backwards] {
            let rings = vec![outer.clone(), inner];
            assert!(!inside_boundary_with::<Exact>(
                &rings,
                Point2::new(2.0, 2.0)
            ));
            assert!(inside_boundary_with::<Exact>(&rings, Point2::new(0.5, 0.5)));
        }
    }
}
