//! Built faces assembled into shells and solids, a piece nested an odd
//! number of times deep a void of the one around it, and each closed
//! body's volume and turned-in faces checked against its mesh.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::Point;
use ogeom_topo::{Model, Shape};

use super::weld::{Piece, inside};
use super::{Carrier, FallbackReason, Groups, MeshSolidReport};
use crate::recognize::Canonical;

/// One closed body of the result: the solid, and the mesh triangles of its
/// shell and of the voids inside it.
pub(super) struct Body {
    solid: Shape,
    triangles: Vec<u32>,
}

/// Faces into one shell per piece; closed pieces into solids, a piece
/// nested an odd number of times deep being a void of the one around it.
/// The result, and each solid with its triangles.
#[allow(
    clippy::too_many_arguments,
    reason = "the build's own state, passed through"
)]
pub(super) fn assemble(
    model: &mut Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    pieces: &[Piece],
    depth: &[usize],
    all_closed: bool,
    groups: &Groups,
    built: &[Option<Shape>],
) -> OgeomResult<(Shape, Vec<Body>)> {
    let mut shells = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let mut faces: Vec<Shape> = Vec::new();
        let mut taken = vec![false; groups.carriers.len()];
        for &t in &piece.triangles {
            let g = groups.of[t as usize];
            if !taken[g] {
                taken[g] = true;
                if let Some(face) = &built[g] {
                    faces.push(face.clone());
                }
            }
        }
        shells.push(model.add_shell(&faces)?);
    }
    if !all_closed {
        let shape = if shells.len() == 1 {
            shells.swap_remove(0)
        } else {
            model.add_compound(&shells)?
        };
        return Ok((shape, Vec::new()));
    }
    let mut bodies = Vec::new();
    for (i, shell) in shells.iter().enumerate() {
        if depth[i] % 2 == 1 {
            continue;
        }
        let mut members = vec![shell.clone()];
        let mut mine = pieces[i].triangles.clone();
        for (j, void) in shells.iter().enumerate() {
            if depth[j] == depth[i] + 1 && inside(points, triangles, &pieces[i], &pieces[j]) {
                members.push(void.clone());
                mine.extend(&pieces[j].triangles);
            }
        }
        bodies.push(Body {
            solid: model.add_solid(&members)?,
            triangles: mine,
        });
    }
    let shape = if bodies.len() == 1 {
        bodies[0].solid.clone()
    } else {
        let solids: Vec<Shape> = bodies.iter().map(|b| b.solid.clone()).collect();
        model.add_compound(&solids)?
    };
    Ok((shape, bodies))
}

/// The recognized regions to facet so every body of the result is sound: a
/// recognized face turned into the material gives up its region, and a
/// body whose volume is not
/// the mesh's, within what its recognized surfaces may add over their
/// facets, gives up all of its recognized regions (a face closed over the
/// wrong part of its surface passes every local test and is caught only
/// there). Other bodies keep theirs. Returned with which of the two it
/// was.
#[allow(
    clippy::too_many_arguments,
    reason = "the build's own state, passed through"
)]
pub(super) fn body_culprits(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    bodies: &[Body],
    flat: f64,
    kept: &crate::mass::VolumeKept,
    tol: Tolerances,
    report: &mut MeshSolidReport,
) -> OgeomResult<(Vec<usize>, FallbackReason)> {
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    let own = |body: &Body| -> Vec<usize> {
        let mut own: Vec<usize> = body
            .triangles
            .iter()
            .map(|&t| groups.of[t as usize])
            .filter(|&g| curved(g))
            .collect();
        own.sort_unstable();
        own.dedup();
        own
    };
    // Orientation first: a recognized face turned against the triangles it
    // replaces faces into the material, and is facetted whatever the
    // volumes say.
    let culprits = inverted_faces(model, points, triangles, groups, built, tol)?;
    if !culprits.is_empty() {
        return Ok((culprits, FallbackReason::TurnedIn));
    }
    let mut culprits: Vec<usize> = Vec::new();
    for body in bodies {
        let own = own(body);
        if own.is_empty() {
            continue;
        }
        let (mesh_volume, allowance, area) =
            volume_allowance(points, triangles, &body.triangles, groups, flat);
        let diagonal = body_diagonal(points, triangles, &body.triangles);
        // Measured no finer than the allowance needs: at the facets' mean
        // offset from their surfaces, whose error over the area is the
        // allowance again, counted in the slack.
        let chord = (allowance / area.max(f64::MIN_POSITIVE)).max(flat);
        let deflection = ogeom_mesh::Deflection::with_chord(chord)?;
        // A body with a swept face is measured on its tessellation at that
        // chord: the exact integral over a fitted profile is held to far
        // finer than the slack asks, at far greater cost.
        let swept = own
            .iter()
            .any(|&g| matches!(&groups.carriers[g], Carrier::Curved(c) if matches!(c.shape, Canonical::Swept(_))));
        let measured = if swept {
            let (mesh, _, _) = ogeom_mesh::triangulate_with_face_meshes(
                model,
                &body.solid,
                deflection,
                Some(kept.meshes()),
                tol,
            )?;
            mesh.triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| mesh.positions[i as usize].to_vector());
                    a.dot(b.cross(c)) / 6.0
                })
                .sum::<f64>()
        } else {
            crate::mass::volume_as_flagged(model, &body.solid, deflection, Some(kept), tol)?.mass
        };
        // Measured at the facets' corners and edge middles, the allowance
        // misses the surface's rise inside a facet and a fitted boundary's
        // wander between the rows; twice over covers both.
        let slack = allowance * 2.0 + area * chord * 2.0 + diagonal.powi(3) * 1e-12;
        if (measured - mesh_volume).abs() > slack {
            // The sweeps answer first: the least proven of the fits. Only
            // where none is left does the body give up its recognition.
            let sweeps: Vec<usize> = own
                .iter()
                .copied()
                .filter(|&g| matches!(&groups.carriers[g], Carrier::Curved(c) if matches!(c.shape, Canonical::Swept(_))))
                .collect();
            if sweeps.is_empty() {
                report.recognition_withdrawn = true;
                culprits.extend(own);
            } else {
                culprits.extend(sweeps);
            }
        }
    }
    culprits.sort_unstable();
    culprits.dedup();
    Ok((culprits, FallbackReason::VolumeOff))
}

/// The recognized regions whose built face points against the triangles it
/// replaces. The mesh is oriented outward before anything is built, so each
/// triangle's normal is the side the material is not on; the face's own
/// outward normal, at the foot of the triangle's middle on its surface,
/// must agree with it. A few triangles vote, and the face is turned only
/// where most of them say so.
fn inverted_faces(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    use ogeom_geom::Surface as _;
    const VOTES: usize = 9;
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); groups.carriers.len()];
    for (t, &g) in groups.of.iter().enumerate() {
        if matches!(groups.carriers.get(g), Some(Carrier::Curved(_))) {
            members[g].push(t);
        }
    }
    let mut inverted = Vec::new();
    for (g, tris) in members.iter().enumerate() {
        let Some(face) = built.get(g).and_then(Option::as_ref) else {
            continue;
        };
        if tris.is_empty() {
            continue;
        }
        let Some(surface) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|data| model.geometry().surface(data.surface))
        else {
            continue;
        };
        let turned = face.orientation() == ogeom_topo::Orientation::Reversed;
        let votes = VOTES.min(tris.len());
        let (mut against, mut asked) = (0, 0);
        for k in 0..votes {
            let [a, b, c] = triangles[tris[k * tris.len() / votes]].map(|v| points[v as usize]);
            let normal = (b - a).cross(c - a);
            let middle = Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) / 3.0);
            let Ok(foot) = crate::measure::project_on_surface(surface, middle, 8, tol) else {
                continue;
            };
            let Ok(outward) = surface.normal_at(foot.parameters.0, foot.parameters.1, tol) else {
                continue;
            };
            let outward = if turned {
                -outward.vector()
            } else {
                outward.vector()
            };
            asked += 1;
            if outward.dot(normal) < 0.0 {
                against += 1;
            }
        }
        if asked > 0 && against * 2 > asked {
            inverted.push(g);
        }
    }
    Ok(inverted)
}

/// The diagonal of the box around some of the mesh's triangles.
fn body_diagonal(points: &[Point], triangles: &[[u32; 3]], mine: &[u32]) -> f64 {
    let mut lo = Point::new(f64::MAX, f64::MAX, f64::MAX);
    let mut hi = Point::new(f64::MIN, f64::MIN, f64::MIN);
    for &t in mine {
        for &v in &triangles[t as usize] {
            let p = points[v as usize];
            lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    lo.distance(hi)
}

/// Some triangles' own volume, how far the solid built on their groups may
/// differ from it, and their area. Each triangle on a recognized surface may
/// stand off it by its worst corner or edge middle, over its whole area; a
/// triangle on a plane by the coplanar distance.
fn volume_allowance(
    points: &[Point],
    triangles: &[[u32; 3]],
    mine: &[u32],
    groups: &Groups,
    flat: f64,
) -> (f64, f64, f64) {
    let (mut volume, mut allowance, mut area) = (0.0, 0.0, 0.0);
    for &t in mine {
        let t = t as usize;
        let tri = &triangles[t];
        let [a, b, c] = tri.map(|v| points[v as usize]);
        let (va, vb, vc) = (a - Point::ORIGIN, b - Point::ORIGIN, c - Point::ORIGIN);
        volume += va.dot(vb.cross(vc)) / 6.0;
        let size = (b - a).cross(c - a).magnitude() / 2.0;
        area += size;
        let off = match groups.carriers.get(groups.of[t]) {
            Some(Carrier::Curved(curved)) => {
                let middle = |p: Point, q: Point| p + (q - p) * 0.5;
                [
                    a,
                    b,
                    c,
                    middle(a, b),
                    middle(b, c),
                    middle(c, a),
                    middle(a, middle(b, c)),
                ]
                .into_iter()
                .map(|p| curved.shape.distance_to(p))
                .fold(0.0_f64, f64::max)
            }
            _ => flat,
        };
        allowance += size * off;
    }
    (volume, allowance, area)
}
