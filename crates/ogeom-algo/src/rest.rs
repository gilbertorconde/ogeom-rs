//! The ways a part can lie on a flat bed.
//!
//! A part lies on a flat bed on a plane of its convex hull: the bed is a
//! plane every point of the part is on one side of, touching it. Each flat
//! face of the hull is one such pose; a recess the bed cannot reach is
//! none, and a hull face across an opening (the bed spanning the two feet
//! of an L) is one even where the part touches it only along edges.
//!
//! # What a rest reports
//!
//! - The contact area is the area of the part's own faces in the bed's
//!   plane, not the hull face's: an L standing on its two feet touches the
//!   bed with the feet alone, and one standing on two edges has no contact
//!   area at all.
//! - The footprint is the support polygon, the convex hull of everywhere
//!   the part touches the bed. That is the hull face.
//! - The centre of mass is taken at uniform density. The pose is stable
//!   where the centre's foot on the bed falls inside the footprint by more
//!   than the in-plane tolerance; on the boundary or outside, the part tips.
//! - The margin is the static stability margin: the distance from the
//!   centre's foot to the nearest footprint edge, negative outside. The
//!   tip is the least angle the part must be turned about a footprint edge
//!   to carry its centre over that edge, `atan(margin / centre height)`.
//!   It is scale-free, and the part's centre rises by
//!   `centre height × (1 / cos tip − 1)` on the way, the energy it takes.
//!
//! # Curved walls
//!
//! A part does not rest on a curved wall: a cylinder on its side rolls,
//! and a sphere has no pose at all. A mesh of either is a polyhedron, and
//! each of its facets on the wall is a hull face with the centre above it.
//! They are told apart by the hull's turn across a footprint edge, the
//! angle between the hull face and the next one over that edge. Below 30°
//! the part bends smoothly there, as a curved surface's facets do. A hull
//! face is a facet of a curved wall, and no rest, where the hull turns
//! smoothly across two of its edges on opposite sides (their inward
//! directions apart by more than a right angle) and the part tips over one
//! of them before the next hull face comes down. On a wall curved about
//! the centre of mass, the turn to tip and the turn to the next face add up
//! to the hull's turn, both positive, so every facet rolls. The same holds
//! for a rounded edge's facets, which tip off on to the flat beside them.
//!
//! A flat face with rounded edges stays a rest: the part has to be tipped
//! much further than one facet of the round before its centre passes the
//! edge. So does a flat face with one smooth edge and sharp ones, stable or
//! not: the end of an L's foot that its upright tips it off, towards a
//! hull face across the L's opening that turns only 20° from it. A curved
//! bottom that rocks the part back upright (a hemisphere on its dome) is no
//! rest either: the part touches the bed at a point, which holds nothing
//! on a bed. A genuinely faceted part whose faces turn by less than 30°
//! rolls over them, as a curved one would.
//!
//! # Order
//!
//! Stable rests come first. Within each group the larger contact area comes
//! first; areas within the in-plane tolerance across the part are equal,
//! and the lower centre of mass then comes first, and then the rest whose
//! `down` points most nearly along −z, then −y, then −x.
//!
//! # Tolerance
//!
//! A point is in the bed's plane within the larger of the confusion
//! distance and a millionth of the part's bounding diagonal, enough to
//! take a mesh rounded to single precision. Hull faces in one plane within
//! that distance are one rest.

use core::cmp::Ordering;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Point, Point2, Vector};
use ogeom_topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

use crate::hull3::convex_hull;

/// The hull's turn across a footprint edge below which the part is taken to
/// bend smoothly there, as a curved surface's facets do. The crease angle
/// face directions are read with (`orient::Directions`).
const SMOOTH_TURN: f64 = core::f64::consts::FRAC_PI_6;

/// The shape's mesh chord as a fraction of its bounding diagonal.
const SHAPE_CHORD: f64 = 1e-3;

/// The greatest turn across one facet of the shape's mesh, in radians:
/// well under [`SMOOTH_TURN`], so a curved face meshes into facets that
/// roll.
const SHAPE_TURN: f64 = 0.25;

/// A way a part can lie on a flat bed.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Rest {
    /// The outward normal of the plane the part lies on: the way that
    /// faces the bed.
    pub down: Direction,
    /// The area of the part's own faces in that plane; zero where it
    /// stands on edges and corners alone.
    pub area: f64,
    /// The part's height above the bed in this pose.
    pub height: f64,
    /// The height of its centre of mass above the bed.
    pub centre_height: f64,
    /// The bed in this pose: its origin the foot of the centre of mass,
    /// its `z` up, opposite `down`. The footprint is in its `x`, `y`.
    pub bed: Frame,
    /// The support polygon, counter-clockwise about the bed's `z`.
    pub footprint: Vec<Point2>,
    /// The distance from the centre's foot (the bed's origin) to the
    /// nearest footprint edge; negative outside the footprint.
    pub margin: f64,
    /// The least angle, in radians, the part must be turned about a
    /// footprint edge to carry its centre over it; negative where it tips
    /// already.
    pub tip: f64,
    /// Whether the centre's foot is inside the footprint by more than the
    /// in-plane tolerance.
    pub stable: bool,
}

/// The planar supports of a closed mesh: each face of its convex hull the
/// part can lie on, at uniform density. See the [module
/// documentation](self) for what is reported and in what order.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a triangle
/// indexing a vertex the mesh does not have, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for a
/// mesh that encloses no volume (open, flat or wound inward).
pub fn resting_faces(mesh: &Triangulation, tol: Tolerances) -> OgeomResult<Vec<Rest>> {
    crate::orient::check_indices(mesh)?;
    let flat = in_plane(&mesh.positions, tol);
    let Some(centre) = centre_of(&mesh.positions, &mesh.triangles, flat) else {
        ogeom_bail!(
            Construction,
            "the mesh encloses no volume: it is open, flat or wound inward"
        );
    };
    rests_of(&mesh.positions, &mesh.triangles, centre, flat)
}

/// [`resting_faces`] of a solid: its mesh, at a chord of a thousandth of
/// its bounding diagonal and a turn of at most 0.25 rad per facet, with the
/// centre of mass and the contact areas taken from the exact geometry. The
/// contact area is that of the solid's planar faces in the bed's plane, so
/// a curved face touching the bed along a line adds nothing. Heights and
/// footprints are the mesh's, within its chord.
///
/// # Errors
///
/// As [`resting_faces`], and as the meshing and
/// [`volume_properties`](crate::volume_properties) of the shape.
pub fn resting_faces_of(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<Rest>> {
    let mut deflection = crate::relative_deflection(model, shape, SHAPE_CHORD, tol)?;
    deflection.angular = deflection.angular.min(SHAPE_TURN);
    let mesh = ogeom_mesh::triangulate(model, shape, deflection, tol)?;
    let flat = in_plane(&mesh.positions, tol);
    let centre = crate::volume_properties(model, shape, deflection, tol)?.centre;
    let mut rests = rests_of(&mesh.positions, &mesh.triangles, centre, flat)?;

    let mut planar: Vec<(Point, Vector, f64)> = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        if !matches!(
            model.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        ) {
            continue;
        }
        let (at, normal) = crate::face_normal(model, &face, tol)?;
        let area = crate::surface_properties(model, &face, deflection, tol)?.mass;
        planar.push((at, normal, area));
    }
    let diagonal = diagonal(&mesh.positions);
    let parallel = 1.0 - 0.5 * (flat / diagonal).powi(2);
    for rest in &mut rests {
        let down = rest.down.vector();
        rest.area = planar
            .iter()
            .filter(|(at, normal, _)| {
                normal.dot(down) >= parallel * normal.magnitude()
                    && down.dot(*at - rest.bed.origin()).abs() <= flat
            })
            .map(|(.., area)| area)
            .sum();
    }
    sort(&mut rests, flat * diagonal, flat);
    Ok(rests)
}

/// The in-plane tolerance for a part with these points.
fn in_plane(points: &[Point], tol: Tolerances) -> f64 {
    (1e-6 * diagonal(points)).max(tol.confusion())
}

fn diagonal(points: &[Point]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for (k, x) in p.to_array().into_iter().enumerate() {
            lo[k] = lo[k].min(x);
            hi[k] = hi[k].max(x);
        }
    }
    if points.is_empty() {
        return 0.0;
    }
    Vector::new(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]).magnitude()
}

/// The centre of the volume the triangles enclose, by the divergence
/// theorem about the first point; `None` where that volume is not
/// positive beyond what the in-plane tolerance over the part could make.
fn centre_of(points: &[Point], triangles: &[[u32; 3]], flat: f64) -> Option<Point> {
    let origin = *points.first()?;
    let mut volume = 0.0;
    let mut moment = Vector::ZERO;
    for &[a, b, c] in triangles {
        let [a, b, c] = [a, b, c].map(|i| points[i as usize] - origin);
        let v = a.dot(b.cross(c)) / 6.0;
        volume += v;
        moment += (a + b + c) * (v / 4.0);
    }
    let size = diagonal(points);
    (volume > flat * size * size).then(|| origin + moment / volume)
}

/// A hull face grown over its neighbours in one plane.
struct Facet {
    normal: Vector,
    offset: f64,
    corners: Vec<u32>,
    /// Each boundary edge's two corners, in the facet's winding, and the
    /// facet across it, if any.
    edges: Vec<(u32, u32, Option<usize>)>,
}

fn rests_of(
    points: &[Point],
    triangles: &[[u32; 3]],
    centre: Point,
    flat: f64,
) -> OgeomResult<Vec<Rest>> {
    let Some(hull) = convex_hull(points) else {
        ogeom_bail!(
            Construction,
            "the part's convex hull could not be built: its points lie in one plane"
        );
    };
    let at = |i: u32| points[i as usize];
    let normals: Vec<Vector> = hull
        .faces
        .iter()
        .map(|&[a, b, c]| (at(b) - at(a)).cross(at(c) - at(a)))
        .collect();

    // Grown from the largest hull faces first, so a plane's facet is
    // seeded by a face that settles its direction.
    let mut order: Vec<usize> = (0..hull.faces.len()).collect();
    order.sort_by(|&f, &g| {
        normals[g]
            .magnitude()
            .total_cmp(&normals[f].magnitude())
            .then(f.cmp(&g))
    });
    let mut facet_of: Vec<Option<usize>> = vec![None; hull.faces.len()];
    let mut members: Vec<Vec<usize>> = Vec::new();
    for &seed in &order {
        let length = normals[seed].magnitude();
        if facet_of[seed].is_some() || length == 0.0 {
            continue;
        }
        let n = normals[seed] / length;
        let base = at(hull.faces[seed][0]);
        let id = members.len();
        facet_of[seed] = Some(id);
        let mut grown = vec![seed];
        let mut walk = vec![seed];
        while let Some(f) = walk.pop() {
            for &g in &hull.across[f] {
                let in_plane = hull.faces[g]
                    .iter()
                    .all(|&i| n.dot(at(i) - base).abs() <= flat);
                if facet_of[g].is_none() && in_plane && normals[g].dot(n) >= 0.0 {
                    facet_of[g] = Some(id);
                    grown.push(g);
                    walk.push(g);
                }
            }
        }
        members.push(grown);
    }
    let facets: Vec<Facet> = members
        .iter()
        .enumerate()
        .map(|(id, grown)| {
            let sum = grown.iter().fold(Vector::ZERO, |s, &f| s + normals[f]);
            let normal = sum / sum.magnitude();
            let mut corners: Vec<u32> = grown.iter().flat_map(|&f| hull.faces[f]).collect();
            corners.sort_unstable();
            corners.dedup();
            let offset = corners
                .iter()
                .map(|&i| normal.dot(at(i) - Point::ORIGIN))
                .fold(f64::NEG_INFINITY, f64::max);
            let mut edges = Vec::new();
            for &f in grown {
                for k in 0..3 {
                    let g = hull.across[f][k];
                    if facet_of[g] != Some(id) {
                        let face = hull.faces[f];
                        edges.push((face[k], face[(k + 1) % 3], facet_of[g]));
                    }
                }
            }
            Facet {
                normal,
                offset,
                corners,
                edges,
            }
        })
        .collect();

    let mut rests = Vec::new();
    for facet in &facets {
        let n = facet.normal;
        let centre_height = facet.offset - n.dot(centre - Point::ORIGIN);
        let foot = centre + n * centre_height;
        let mut margin = f64::INFINITY;
        // Each edge the hull turns smoothly across: its inward direction
        // in the facet and whether the part tips over it before the next
        // hull face comes down.
        let mut smooth: Vec<(Vector, bool)> = Vec::new();
        for &(u, v, across) in &facet.edges {
            let inward = n.cross(at(v) - at(u));
            let length = inward.magnitude();
            if length == 0.0 {
                continue;
            }
            let inward = inward / length;
            let m = inward.dot(foot - at(u));
            margin = margin.min(m);
            let turn = across.map_or(core::f64::consts::FRAC_PI_2, |g| {
                n.dot(facets[g].normal).clamp(-1.0, 1.0).acos()
            });
            if turn < SMOOTH_TURN {
                smooth.push((inward, m.atan2(centre_height) <= turn));
            }
        }
        let rolls = smooth
            .iter()
            .any(|&(a, tips)| tips && smooth.iter().any(|&(b, _)| a.dot(b) < 0.0));
        if rolls {
            continue;
        }
        let down = Direction::new(n, Tolerances::millimetres())?;
        let up = Direction::new(-n, Tolerances::millimetres())?;
        let bed = Frame::about(foot, up);
        let local: Vec<Point2> = facet
            .corners
            .iter()
            .map(|&i| {
                let p = bed.to_local(at(i));
                Point2::new(p.x, p.y)
            })
            .collect();
        let lowest = points
            .iter()
            .map(|p| n.dot(*p - Point::ORIGIN))
            .fold(f64::INFINITY, f64::min);
        let near = |i: u32| (facet.offset - n.dot(at(i) - Point::ORIGIN)).abs() <= flat;
        let area = triangles
            .iter()
            .filter(|t| t.iter().all(|&i| near(i)))
            .map(|&[a, b, c]| (at(b) - at(a)).cross(at(c) - at(a)))
            .filter(|cross| cross.dot(n) > 0.0)
            .map(|cross| cross.magnitude() / 2.0)
            .sum();
        rests.push(Rest {
            down,
            area,
            height: facet.offset - lowest,
            centre_height,
            bed,
            footprint: crate::fit::segments::hull(&local),
            margin,
            tip: margin.atan2(centre_height),
            stable: margin > flat,
        });
    }
    let size = diagonal(points);
    sort(&mut rests, flat * size, flat);
    Ok(rests)
}

/// The rests in the module's order: stable first, then by contact area
/// within `area_tie`, then by centre height within `height_tie`, then by
/// `down`.
fn sort(rests: &mut [Rest], area_tie: f64, height_tie: f64) {
    let areas = groups(rests.iter().map(|r| -r.area).collect(), area_tie);
    let heights = groups(rests.iter().map(|r| r.centre_height).collect(), height_tie);
    let mut keyed: Vec<(usize, Rest)> = rests.iter().cloned().enumerate().collect();
    keyed.sort_by(|(i, a), (j, b)| {
        let along = |r: &Rest| {
            let d = r.down.vector();
            [d.z, d.y, d.x]
        };
        b.stable
            .cmp(&a.stable)
            .then(areas[*i].cmp(&areas[*j]))
            .then(heights[*i].cmp(&heights[*j]))
            .then_with(|| {
                along(a)
                    .iter()
                    .zip(along(b))
                    .map(|(x, y)| x.total_cmp(&y))
                    .find(|o| *o != Ordering::Equal)
                    .unwrap_or(Ordering::Equal)
            })
    });
    for (slot, (_, rest)) in rests.iter_mut().zip(keyed) {
        *slot = rest;
    }
}

/// Each value's rank among runs of values sorted ascending, a run going on
/// while each next value is within `tie` of the one before.
fn groups(values: Vec<f64>, tie: f64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let mut rank = vec![0; values.len()];
    let mut current = 0;
    for pair in order.windows(2) {
        if values[pair[1]] - values[pair[0]] > tie {
            current += 1;
        }
        rank[pair[1]] = current;
    }
    rank
}
