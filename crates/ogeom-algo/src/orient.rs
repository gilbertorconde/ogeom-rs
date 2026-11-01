//! Directions and frames read off a mesh.
//!
//! Four questions a scan is asked before anything is fitted to it: which
//! ways its flat faces face ([`face_directions`]), which axis its walls turn
//! about ([`axis_of_revolution`]), along which axes its area spreads
//! ([`principal_frame`]), and where it is open ([`open_rims`]).
//!
//! Every answer is a computation over the triangles alone. The same mesh
//! gives the same answer bit for bit on any thread count, and where two
//! candidates tie they are told apart by where they point, never by the
//! order they were found in.
//!
//! # Face directions
//!
//! Each triangle's normal is smoothed over its neighbours, across the edges
//! that are not creases, so a scan's noise averages out within a face and
//! never between two. The smoothed normals are binned on the sphere by area,
//! a direction and its antipode in one bin, since a face and the face
//! opposite it name one direction. Every bin holding more area than its
//! neighbours seeds a search for the direction the area gathers at: the
//! mean of the normals within the resolution, moved to and taken again
//! until it stops moving. Those directions are kept by area, largest first:
//! one within the resolution of a direction already kept is the same
//! direction, and one is only kept where it holds more area within the
//! resolution than in the ring from the resolution to three times it, the
//! area of any direction already kept aside. A cylinder's normals sweep a
//! circle and put twice as much in that ring as within the resolution, so
//! a cylinder gives no face direction, and a plane's put none there.
//!
//! A triangle then belongs to the kept direction its own normal or its
//! smoothed one is nearest, where either is within the resolution. One
//! whose own normal is within the resolution of the direction counts with
//! it; one claimed by its smoothed normal alone counts where its own normal
//! is within the [`Directions::cutoff`] of the face's scatter. A scan's
//! noise scatters a face's own normals past the resolution and the cutoff
//! reaches out with them; an exact face's do not scatter, and a fillet's
//! triangle that smoothing drew towards the face is left out. A face
//! direction is the area-weighted mean of its triangles' own normals,
//! unsmoothed: exact where the faces are exact, and for a scan's noise the
//! mean that noise does not bias.

use core::cmp::Ordering;
use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6, PI};
use std::collections::BTreeSet;

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail, parallel};
use ogeom_math::{Axis, Direction, Frame, Point, Vector};
use ogeom_topo::Triangulation;

use crate::Cutoff;
use crate::recognize::{covariance, smallest_direction};

/// How [`face_directions`] and [`axis_of_revolution`] read a mesh's
/// normals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Directions {
    /// The angle, in radians, within which normals are one direction: face
    /// directions closer than this are one, a face direction is the mean
    /// normal of the triangles within it, and a wall is a triangle whose
    /// normal is within it of perpendicular to an axis. At most π/6, so the
    /// ring out to three times it stays on one side of the perpendicular.
    pub resolution: f64,
    /// The angle, in radians, between two neighbouring triangles' normals
    /// at which smoothing stops: an edge this sharp is a crease, and the
    /// faces either side of it are smoothed apart.
    pub crease: f64,
    /// Rounds of smoothing: each averages every triangle's normal with its
    /// neighbours' across the edges that are not creases, by area.
    pub rounds: u32,
    /// How far past the resolution a face direction reaches for the
    /// triangles only their smoothed normals put within it: to the
    /// cutoff's threshold over the angles of its triangles' own normals
    /// from their mean, weighed by area, and never less than the
    /// resolution. [`Cutoff::All`] takes every one.
    pub cutoff: Cutoff,
}

impl Default for Directions {
    /// A resolution of 5°, a crease of 30°, three rounds and a cutoff at
    /// three RMS: a draft of a degree or two is one direction with its
    /// face, the angle a scan's noise of a tenth of a millimetre opens
    /// between neighbouring triangles a millimetre across stays well inside
    /// the crease, and three rounds average each normal with some thirty
    /// triangles'. A face's noise spreads its normals' angles as a Rayleigh
    /// distribution, so three RMS leave out one of a face's triangles in
    /// some eight thousand (`e⁻⁹`).
    fn default() -> Self {
        Self {
            resolution: 5.0_f64.to_radians(),
            crease: 30.0_f64.to_radians(),
            rounds: 3,
            cutoff: Cutoff::Sigma(3.0),
        }
    }
}

impl Directions {
    /// The options checked: a resolution in (0, π/6], a crease in (0, π]
    /// and the cutoff's number finite and positive.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for any of
    /// them out of its range or not finite.
    pub fn check(&self) -> OgeomResult<()> {
        self.cutoff.check()?;
        if !(self.resolution.is_finite() && self.resolution > 0.0 && self.resolution <= FRAC_PI_6) {
            ogeom_bail!(
                Domain,
                "the resolution {} must be finite, positive and at most π/6",
                self.resolution
            );
        }
        if !(self.crease.is_finite() && self.crease > 0.0 && self.crease <= PI) {
            ogeom_bail!(
                Domain,
                "the crease angle {} must be finite, positive and at most π",
                self.crease
            );
        }
        Ok(())
    }
}

/// A direction a mesh's flat faces face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceDirection {
    /// The area-weighted mean normal of the triangles that face it, a
    /// triangle facing the opposite way counted reversed. Its sense is the
    /// one more of their area faces; where both hold the same area to
    /// within the confusion distance's worth across a face of that area,
    /// the sense whose largest component is positive.
    pub direction: Direction,
    /// The area of those triangles.
    pub area: f64,
    /// The area-weighted RMS angle, in radians, between their normals and
    /// the direction: the noise, and any draft between faces merged into
    /// one direction.
    pub spread: f64,
}

/// The axis a mesh's walls turn about.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Revolution {
    /// The axis, located level with the walls' centroid; its sense is the
    /// one whose largest component is positive.
    pub axis: Axis,
    /// The share of the mesh's area on walls about the axis: triangles
    /// whose normal is within the resolution of pointing straight at it.
    /// Near one for a turned part with small ends, near the resolution's
    /// sine for a prism, whose walls point at an axis only along a strip.
    pub share: f64,
}

/// Where a mesh is open: a chain of edges each of which only one triangle
/// uses.
#[derive(Debug, Clone, PartialEq)]
pub struct Rim {
    /// The vertices in the order the mesh's triangles traverse the rim; the
    /// last joins the first where the rim is closed.
    pub vertices: Vec<u32>,
    /// The vector area of the face that would close the rim, wound the way
    /// the mesh is: the sum of area × normal over any surface spanning it,
    /// which is the same for every one. Its length is the opening's area
    /// where the rim is planar, and it points out of the mesh through the
    /// opening, so it can vote as the face that would close it.
    pub vector_area: Vector,
    /// The rim's centroid, by length.
    pub centre: Point,
    /// The greatest distance of a rim vertex from the plane through the
    /// centre square to the vector area: zero for a planar rim.
    pub flatness: f64,
    /// Whether the chain returns to where it started. A rim through a
    /// vertex where the triangles round it do not meet in one fan breaks
    /// there, and its vector area is the chain's closed by its chord.
    pub closed: bool,
}

/// The directions a mesh's flat faces face, largest area first.
///
/// See the [module documentation](self) for what is reported. A direction
/// and its antipode are one face direction. An empty mesh, or one whose
/// triangles all have no area, has none.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for options out
/// of range, or a triangle indexing a vertex the mesh does not have.
pub fn face_directions(
    mesh: &Triangulation,
    options: &Directions,
    tol: Tolerances,
) -> OgeomResult<Vec<FaceDirection>> {
    options.check()?;
    let read = Read::of(mesh, options)?;
    let kept = read.kept_directions(tol);
    Ok(read.face_directions(&kept, tol))
}

/// The axis a mesh's walls turn about, and the share of its area on them.
///
/// A wall is a triangle whose normal is within the resolution of
/// perpendicular to a direction; the axis direction is the one the walls'
/// normals are most nearly perpendicular to, the least eigenvector of their
/// area-weighted second moments, and the axis line the one their normal
/// lines meet, by least squares: a normal line `(c, n)` meets the axis
/// `(p, a)` where `(c − p)·(a × n) = 0`, which is linear in `p`. The walls
/// are then the triangles whose normal is within the resolution of pointing
/// straight at that line, and the two are found again from them until they
/// stop changing. Walls at any radius serve, so a stepped shaft gives the
/// axis its steps share.
///
/// The search starts from each of the face directions and from the
/// principal directions of the normals, and keeps the axis that gathers the
/// most wall area. A shoulder or an end is not a wall and does not count
/// towards the share, nor does a cone's taper.
///
/// `None` where no direction has walls whose normal lines locate a line:
/// an empty mesh, one plane, a sphere's single point.
///
/// # Errors
///
/// As [`face_directions`].
pub fn axis_of_revolution(
    mesh: &Triangulation,
    options: &Directions,
    tol: Tolerances,
) -> OgeomResult<Option<Revolution>> {
    options.check()?;
    let read = Read::of(mesh, options)?;
    let total: f64 = read.facets.area.iter().sum();
    if total <= 0.0 {
        return Ok(None);
    }
    let mut starts: Vec<Vector> = read.kept_directions(tol);
    let spread = covariance(
        read.smooth
            .iter()
            .zip(&read.facets.area)
            .map(|(s, a)| *s * a.sqrt()),
    );
    let eigen = nalgebra::SymmetricEigen::new(spread);
    for i in 0..3 {
        let v = eigen.eigenvectors.column(i);
        let v = Vector::new(v[0], v[1], v[2]);
        if v.magnitude() > 0.0 {
            starts.push(canonical(v / v.magnitude()));
        }
    }
    let found: Vec<Option<(Vector, Point, f64)>> =
        parallel::map_ordered(&starts, |_, start| read.revolution_from(*start, tol));
    let mut best: Option<(Vector, Point, f64)> = None;
    for (direction, location, area) in found.into_iter().flatten() {
        let better = match best {
            None => true,
            Some((d, _, a)) => match area.total_cmp(&a) {
                Ordering::Greater => true,
                Ordering::Less => false,
                Ordering::Equal => lexical(direction, d) == Ordering::Less,
            },
        };
        if better {
            best = Some((direction, location, area));
        }
    }
    let Some((direction, location, area)) = best else {
        return Ok(None);
    };
    Ok(Some(Revolution {
        axis: Axis::new(location, Direction::new(direction, tol)?),
        share: area / total,
    }))
}

/// The frame a mesh's area spreads along: its origin at the area centroid,
/// `x` along the greatest spread, `y` the next and `z` the least, from the
/// second moments of the surface's area, integrated exactly over each
/// triangle.
///
/// Each of `x` and `y` points the way the area's third moment along it is
/// positive, the way the area reaches further; `z` completes a
/// right-handed frame. A third moment smaller than moving the vertices by
/// the confusion distance could change it does not decide, and the sense
/// whose largest component is positive is taken instead. Spreads equal to
/// within what that move could change do not tell their axes apart either:
/// among them, the frame takes the world axes' projections, X first, the
/// longest of them where it has a choice.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a triangle
/// indexing a vertex the mesh does not have, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for a
/// mesh with no area.
pub fn principal_frame(mesh: &Triangulation, tol: Tolerances) -> OgeomResult<Frame> {
    check_indices(mesh)?;
    let Some(&reference) = mesh.positions.first() else {
        ogeom_bail!(Construction, "an empty mesh has no principal frame");
    };
    let corners = |t: &[u32; 3]| t.map(|i| mesh.positions[i as usize] - reference);
    let mut area = 0.0;
    let mut first = Vector::ZERO;
    for t in &mesh.triangles {
        let [a, b, c] = corners(t);
        let w = (b - a).cross(c - a).magnitude() * 0.5;
        if w.is_finite() {
            area += w;
            first += (a + b + c) * (w / 3.0);
        }
    }
    if area.is_nan() || area <= 0.0 {
        ogeom_bail!(Construction, "a mesh with no area has no principal frame");
    }
    let centre = first / area;
    let mut second = nalgebra::Matrix3::<f64>::zeros();
    for t in &mesh.triangles {
        let [a, b, c] = corners(t).map(|p| p - centre);
        let w = (b - a).cross(c - a).magnitude() * 0.5;
        if !w.is_finite() || w <= 0.0 {
            continue;
        }
        let s = a + b + c;
        let column = |v: Vector| nalgebra::Vector3::new(v.x, v.y, v.z);
        let sum = [a, b, c, s]
            .iter()
            .map(|v| column(*v) * column(*v).transpose())
            .fold(nalgebra::Matrix3::zeros(), |m, x| m + x);
        second += sum * (w / 12.0);
    }
    let second = second / area;
    let eigen = nalgebra::SymmetricEigen::new(second);
    let mut axes: Vec<(f64, Vector)> = (0..3)
        .map(|i| {
            let v = eigen.eigenvectors.column(i);
            (eigen.eigenvalues[i].max(0.0), Vector::new(v[0], v[1], v[2]))
        })
        .collect();
    axes.sort_by(|p, q| q.0.total_cmp(&p.0).then_with(|| lexical(p.1, q.1)));
    let reach = tol.confusion();
    let tied = |i: usize, j: usize| {
        (axes[i].0 - axes[j].0).abs() <= 2.0 * reach * (axes[i].0.sqrt() + axes[j].0.sqrt())
    };
    let third = |e: Vector| -> f64 {
        mesh.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = corners(t).map(|p| p - centre);
                let w = (b - a).cross(c - a).magnitude() * 0.5;
                if !w.is_finite() {
                    return 0.0;
                }
                let [f, g, h] = [a.dot(e), b.dot(e), c.dot(e)];
                let mixed = f * f * (g + h) + g * g * (f + h) + h * h * (f + g) + f * g * h;
                w / 10.0 * (f * f * f + g * g * g + h * h * h + mixed)
            })
            .sum()
    };
    let sensed = |e: Vector, spread: f64| -> Vector {
        let skew = third(e);
        if skew.abs() > 3.0 * reach * area * spread {
            if skew > 0.0 { e } else { -e }
        } else {
            canonical(e)
        }
    };
    let world = [Vector::X, Vector::Y, Vector::Z];
    let (x, y) = match (tied(0, 1), tied(1, 2)) {
        (false, false) => (sensed(axes[0].1, axes[0].0), sensed(axes[1].1, axes[1].0)),
        (true, false) => {
            let z = sensed(axes[2].1, axes[2].0);
            let x = projected_world(&world, z);
            (x, z.cross(x))
        }
        (false, true) => {
            let x = sensed(axes[0].1, axes[0].0);
            (x, projected_world(&world, x))
        }
        (true, true) => (Vector::X, Vector::Y),
    };
    let z = x.cross(y);
    let origin = reference + centre;
    Frame::from_axes(
        origin,
        Direction::new(x, tol)?,
        Direction::new(y, tol)?,
        Direction::new(z, tol)?,
        tol,
    )
}

/// The rims where a mesh is open, the largest vector area first.
///
/// A rim edge is an edge one triangle uses; an edge two use is inside the
/// mesh however they are wound, and one used three times or more is not a
/// rim either. From a rim edge's end the next is found by turning through
/// the triangles round that vertex, across the edges two of them share,
/// until an edge only one uses: the rim the same fan of triangles bounds,
/// so two rims that touch at a vertex stay two.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a triangle
/// indexing a vertex the mesh does not have.
pub fn open_rims(mesh: &Triangulation) -> OgeomResult<Vec<Rim>> {
    check_indices(mesh)?;
    let uses = edge_uses(&mesh.triangles);
    let mut rim: Vec<(u32, u32, u32)> = Vec::new();
    for (t, tri) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if a != b && uses.get(&key(a, b)).is_some_and(|u| u.count == 1) {
                rim.push((a, b, as_index(t)));
            }
        }
    }
    let at: FastMap<(u32, u32), usize> = rim
        .iter()
        .enumerate()
        .map(|(i, &(a, b, _))| ((a, b), i))
        .collect();
    let limit = mesh.triangles.len();
    let next: Vec<Option<usize>> = rim
        .iter()
        .map(|&(a, b, t)| {
            let mut tri = t as usize;
            let mut far = third_vertex(mesh.triangles[tri], a, b)?;
            for _ in 0..=limit {
                let uses = uses.get(&key(b, far))?;
                match uses.count {
                    1 => return at.get(&(b, far)).copied(),
                    2 => {
                        let other = if uses.first as usize == tri {
                            uses.second
                        } else {
                            uses.first
                        } as usize;
                        let turned = third_vertex(mesh.triangles[other], b, far)?;
                        tri = other;
                        far = turned;
                    }
                    _ => return None,
                }
            }
            None
        })
        .collect();
    let mut entered = vec![false; rim.len()];
    for n in next.iter().flatten() {
        entered[*n] = true;
    }
    let mut taken = vec![false; rim.len()];
    let mut chains: Vec<(Vec<usize>, bool)> = Vec::new();
    let walk = |start: usize, taken: &mut Vec<bool>| -> (Vec<usize>, bool) {
        let mut chain = Vec::new();
        let mut at = start;
        loop {
            taken[at] = true;
            chain.push(at);
            match next[at] {
                Some(n) if n == start => return (chain, true),
                Some(n) if !taken[n] => at = n,
                _ => return (chain, false),
            }
        }
    };
    for start in 0..rim.len() {
        if !entered[start] && !taken[start] {
            chains.push(walk(start, &mut taken));
        }
    }
    for start in 0..rim.len() {
        if !taken[start] {
            chains.push(walk(start, &mut taken));
        }
    }
    let position = |v: u32| mesh.positions[v as usize];
    let mut rims: Vec<Rim> = chains
        .into_iter()
        .map(|(mut chain, closed)| {
            if closed {
                let lowest = (0..chain.len())
                    .min_by(|&i, &j| {
                        let (a, b) = (rim[chain[i]].0, rim[chain[j]].0);
                        lexical(position(a).to_vector(), position(b).to_vector()).then(a.cmp(&b))
                    })
                    .unwrap_or(0);
                chain.rotate_left(lowest);
            }
            let mut vertices: Vec<u32> = chain.iter().map(|&e| rim[e].0).collect();
            if !closed && let Some(&last) = chain.last() {
                vertices.push(rim[last].1);
            }
            rim_of(vertices, closed, &position)
        })
        .collect();
    rims.sort_by(|p, q| {
        q.vector_area
            .magnitude()
            .total_cmp(&p.vector_area.magnitude())
            .then_with(|| lexical(p.centre.to_vector(), q.centre.to_vector()))
    });
    Ok(rims)
}

/// A rim's measures from its vertices in order.
fn rim_of(vertices: Vec<u32>, closed: bool, position: &impl Fn(u32) -> Point) -> Rim {
    let count = vertices.len();
    let segments = |k: usize| (position(vertices[k]), position(vertices[(k + 1) % count]));
    let edges = if closed {
        count
    } else {
        count.saturating_sub(1)
    };
    let mut length = 0.0;
    let mut weighted = Vector::ZERO;
    let reference = position(vertices[0]);
    for k in 0..edges {
        let (p, q) = segments(k);
        let l = p.distance(q);
        length += l;
        weighted += ((p - reference) + (q - reference)) * (l * 0.5);
    }
    let centre = if length > 0.0 {
        reference + weighted / length
    } else {
        let sum = vertices
            .iter()
            .fold(Vector::ZERO, |s, v| s + (position(*v) - reference));
        #[allow(
            clippy::cast_precision_loss,
            reason = "vertex counts are far below 2^52"
        )]
        let n = count as f64;
        reference + sum / n
    };
    let mut traversed = Vector::ZERO;
    for k in 0..count {
        let (p, q) = segments(k);
        traversed += (p - centre).cross(q - centre);
    }
    let vector_area = traversed * -0.5;
    let size = vector_area.magnitude();
    let flatness = if size > 0.0 {
        let n = vector_area / size;
        vertices
            .iter()
            .map(|v| (position(*v) - centre).dot(n).abs())
            .fold(0.0, f64::max)
    } else {
        0.0
    };
    Rim {
        vertices,
        vector_area,
        centre,
        flatness,
        closed,
    }
}

/// The vertex of a triangle that is neither `a` nor `b`.
fn third_vertex(tri: [u32; 3], a: u32, b: u32) -> Option<u32> {
    tri.into_iter().find(|&v| v != a && v != b)
}

/// The triangles that use an edge: how many, and the first two.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Uses {
    pub(crate) count: u32,
    pub(crate) first: u32,
    pub(crate) second: u32,
}

pub(crate) fn key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

fn as_index(t: usize) -> u32 {
    #[allow(
        clippy::cast_possible_truncation,
        reason = "triangles are indexed by u32 throughout the mesh"
    )]
    let t = t as u32;
    t
}

pub(crate) fn edge_uses(triangles: &[[u32; 3]]) -> FastMap<(u32, u32), Uses> {
    let mut uses: FastMap<(u32, u32), Uses> = FastMap::default();
    uses.reserve(triangles.len() * 3 / 2 + 3);
    for (t, tri) in triangles.iter().enumerate() {
        let t = as_index(t);
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if a == b {
                continue;
            }
            let entry = uses.entry(key(a, b)).or_insert(Uses {
                count: 0,
                first: t,
                second: t,
            });
            if entry.count == 1 {
                entry.second = t;
            }
            entry.count += 1;
        }
    }
    uses
}

pub(crate) fn check_indices(mesh: &Triangulation) -> OgeomResult<()> {
    let count = mesh.positions.len();
    if let Some(t) = mesh
        .triangles
        .iter()
        .find(|t| t.iter().any(|&i| i as usize >= count))
    {
        ogeom_bail!(
            Domain,
            "a triangle {t:?} indexes a vertex past the mesh's {count}"
        );
    }
    Ok(())
}

/// The sense of `v` whose largest component is positive, the first of
/// equal ones.
fn canonical(v: Vector) -> Vector {
    let major = major_axis(v);
    if v.to_array()[major] < 0.0 { -v } else { v }
}

fn major_axis(v: Vector) -> usize {
    let [x, y, z] = v.to_array().map(f64::abs);
    if x >= y && x >= z {
        0
    } else if y >= z {
        1
    } else {
        2
    }
}

/// Vectors ordered by their components, x first.
fn lexical(a: Vector, b: Vector) -> Ordering {
    a.x.total_cmp(&b.x)
        .then(a.y.total_cmp(&b.y))
        .then(a.z.total_cmp(&b.z))
}

/// The world axis with the longest projection square to `normal`, X first
/// of equal ones, projected and made unit, its sense the axis's.
fn projected_world(world: &[Vector; 3], normal: Vector) -> Vector {
    let mut best = Vector::ZERO;
    for w in world {
        let p = *w - normal * normal.dot(*w);
        if p.magnitude() > best.magnitude() {
            best = p;
        }
    }
    best / best.magnitude()
}

/// Each triangle's unit normal, area and centroid; a triangle without a
/// normal has no area and a zero normal.
#[derive(Debug)]
struct Facets {
    normal: Vec<Vector>,
    area: Vec<f64>,
    centroid: Vec<Point>,
}

fn facets(mesh: &Triangulation) -> Facets {
    let count = mesh.triangles.len();
    let mut normal = Vec::with_capacity(count);
    let mut area = Vec::with_capacity(count);
    let mut centroid = Vec::with_capacity(count);
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let cross = (b - a).cross(c - a);
        let length = cross.magnitude();
        if length.is_finite() && length > 0.0 {
            normal.push(cross / length);
            area.push(length * 0.5);
        } else {
            normal.push(Vector::ZERO);
            area.push(0.0);
        }
        centroid.push(a + ((b - a) + (c - a)) / 3.0);
    }
    Facets {
        normal,
        area,
        centroid,
    }
}

/// Each triangle's normal averaged by area with its neighbours' across the
/// edges that are not creases, `rounds` times over, every round from the
/// last round's normals.
fn smoothed(mesh: &Triangulation, facets: &Facets, options: &Directions) -> Vec<Vector> {
    let mut smooth = facets.normal.clone();
    if options.rounds == 0 {
        return smooth;
    }
    let uses = edge_uses(&mesh.triangles);
    let limit = options.crease.cos();
    let neighbours: Vec<[u32; 3]> = mesh
        .triangles
        .iter()
        .enumerate()
        .map(|(t, tri)| {
            let mut out = [u32::MAX; 3];
            if facets.area[t] <= 0.0 {
                return out;
            }
            for k in 0..3 {
                let Some(u) = uses.get(&key(tri[k], tri[(k + 1) % 3])) else {
                    continue;
                };
                if u.count != 2 {
                    continue;
                }
                let other = if u.first as usize == t {
                    u.second
                } else {
                    u.first
                };
                let o = other as usize;
                if facets.area[o] > 0.0 && facets.normal[t].dot(facets.normal[o]) > limit {
                    out[k] = other;
                }
            }
            out
        })
        .collect();
    for _ in 0..options.rounds {
        let last = smooth;
        smooth = (0..last.len())
            .map(|t| {
                let mut sum = last[t] * facets.area[t];
                for &n in &neighbours[t] {
                    if n != u32::MAX {
                        sum += last[n as usize] * facets.area[n as usize];
                    }
                }
                let m = sum.magnitude();
                if m > 0.0 { sum / m } else { last[t] }
            })
            .collect();
    }
    smooth
}

/// An equiangular cube map over directions with a direction and its
/// antipode in one cell: three faces of `n × n` cells, the face the one of
/// the largest component, made positive.
#[derive(Debug, Clone, Copy)]
struct Grid {
    n: u64,
    width: f64,
}

impl Grid {
    /// Cells half the resolution wide or narrower.
    fn new(resolution: f64) -> Self {
        let cells = (FRAC_PI_2 / (resolution * 0.5)).ceil().max(1.0);
        #[allow(
            clippy::cast_possible_truncation,
            reason = "the resolution is at least what makes cells a few thousand a side"
        )]
        let n = cells as u64;
        #[allow(clippy::cast_precision_loss, reason = "n is small")]
        let width = FRAC_PI_2 / n as f64;
        Self { n, width }
    }

    /// The farthest a point of a cell lies from its centre: half the
    /// cell's diagonal in the map's angles, `w/√2`, stretched by at most
    /// the map's greatest stretch, `2/√3`, which it reaches at the cube's
    /// corners.
    fn reach(self) -> f64 {
        self.width * (2.0_f64 / 3.0).sqrt()
    }

    fn index(self, angle: f64) -> u64 {
        let at = ((angle + FRAC_PI_4) / self.width).floor().max(0.0);
        #[allow(
            clippy::cast_possible_truncation,
            reason = "clamped to the cells a side"
        )]
        let at = at as u64;
        at.min(self.n - 1)
    }

    fn cell(self, d: Vector) -> u64 {
        let d = canonical(d).to_array();
        let m = major_axis(Vector::from_array(d));
        let (a, b) = (d[(m + 1) % 3], d[(m + 2) % 3]);
        let (i, j) = (self.index(a.atan2(d[m])), self.index(b.atan2(d[m])));
        m as u64 * self.n * self.n + i * self.n + j
    }

    fn angles(self, cell: u64) -> (usize, f64, f64) {
        let m = cell / (self.n * self.n);
        let (i, j) = ((cell / self.n) % self.n, cell % self.n);
        #[allow(clippy::cast_precision_loss, reason = "cell indices are small")]
        let at = |k: u64| (k as f64 + 0.5).mul_add(self.width, -FRAC_PI_4);
        #[allow(clippy::cast_possible_truncation, reason = "a face index is 0, 1 or 2")]
        let m = m as usize;
        (m, at(i), at(j))
    }

    fn direction(m: usize, u: f64, v: f64) -> Vector {
        let mut d = [0.0; 3];
        d[m] = 1.0;
        d[(m + 1) % 3] = u.tan();
        d[(m + 2) % 3] = v.tan();
        let d = Vector::from_array(d);
        d / d.magnitude()
    }

    fn centre(self, cell: u64) -> Vector {
        let (m, u, v) = self.angles(cell);
        Self::direction(m, u, v)
    }

    fn neighbours(self, cell: u64) -> impl Iterator<Item = u64> {
        let (m, u, v) = self.angles(cell);
        let w = self.width;
        [-1.0, 0.0, 1.0].into_iter().flat_map(move |du| {
            [-1.0, 0.0, 1.0].into_iter().filter_map(move |dv| {
                if du == 0.0 && dv == 0.0 {
                    return None;
                }
                let next = self.cell(Self::direction(m, u + du * w, v + dv * w));
                (next != cell).then_some(next)
            })
        })
    }
}

/// A mesh read for its normals: the facets, their smoothed normals, and
/// both binned on the cube map.
#[derive(Debug)]
struct Read {
    facets: Facets,
    smooth: Vec<Vector>,
    grid: Grid,
    /// The triangles binned by their smoothed normals.
    smoothed: Binned,
    /// The triangles binned by their own normals.
    own: Binned,
    resolution: f64,
    cutoff: Cutoff,
}

/// Triangles binned by a normal on the cube map: the cells holding one,
/// ascending, where each cell's triangles start, and the triangles by cell
/// and then by index.
#[derive(Debug)]
struct Binned {
    cells: Vec<u64>,
    starts: Vec<usize>,
    members: Vec<u32>,
}

impl Binned {
    fn of(grid: Grid, normals: &[Vector], area: &[f64]) -> Self {
        let mut binned: Vec<(u64, u32)> = normals
            .iter()
            .enumerate()
            .filter(|(t, _)| area[*t] > 0.0)
            .map(|(t, s)| (grid.cell(*s), as_index(t)))
            .collect();
        binned.sort_unstable();
        let mut cells = Vec::new();
        let mut starts = Vec::new();
        for (k, (cell, _)) in binned.iter().enumerate() {
            if cells.last() != Some(cell) {
                cells.push(*cell);
                starts.push(k);
            }
        }
        starts.push(binned.len());
        let members = binned.into_iter().map(|(_, t)| t).collect();
        Self {
            cells,
            starts,
            members,
        }
    }

    fn in_cell(&self, cell: u64) -> &[u32] {
        match self.cells.binary_search(&cell) {
            Ok(k) => &self.members[self.starts[k]..self.starts[k + 1]],
            Err(_) => &[],
        }
    }
}

impl Read {
    fn of(mesh: &Triangulation, options: &Directions) -> OgeomResult<Self> {
        check_indices(mesh)?;
        let facets = facets(mesh);
        let smooth = smoothed(mesh, &facets, options);
        let grid = Grid::new(options.resolution);
        let smoothed = Binned::of(grid, &smooth, &facets.area);
        let own = Binned::of(grid, &facets.normal, &facets.area);
        Ok(Self {
            facets,
            smooth,
            grid,
            smoothed,
            own,
            resolution: options.resolution,
            cutoff: options.cutoff,
        })
    }

    fn in_cell(&self, cell: u64) -> &[u32] {
        self.smoothed.in_cell(cell)
    }

    /// The cells that can hold a direction within `radius` of `d` or of
    /// its antipode, ascending.
    fn cells_near(&self, d: Vector, radius: f64) -> Vec<u64> {
        let limit = (radius + self.grid.reach()).min(FRAC_PI_2).cos();
        let near = |cell: u64| self.grid.centre(cell).dot(d).abs() >= limit;
        let start = self.grid.cell(d);
        let mut seen = BTreeSet::from([start]);
        let mut queue = vec![start];
        while let Some(cell) = queue.pop() {
            for next in self.grid.neighbours(cell) {
                if near(next) && seen.insert(next) {
                    queue.push(next);
                }
            }
        }
        seen.into_iter().collect()
    }

    /// Every triangle whose own normal is within `radius` of `d` or of its
    /// antipode, by cell and then by index, with that normal's dot with
    /// `d`.
    fn near_own(&self, d: Vector, radius: f64) -> Vec<(u32, f64)> {
        let limit = radius.cos();
        let mut out = Vec::new();
        for cell in self.cells_near(d, radius) {
            for &t in self.own.in_cell(cell) {
                let along = self.facets.normal[t as usize].dot(d);
                if along.abs() >= limit {
                    out.push((t, along));
                }
            }
        }
        out
    }

    /// Every triangle whose smoothed normal is within `radius` of `d` or of
    /// its antipode, by cell and then by index, with that normal's dot
    /// with `d`.
    fn near(&self, d: Vector, radius: f64) -> Vec<(u32, f64)> {
        let limit = radius.cos();
        let mut out = Vec::new();
        for cell in self.cells_near(d, radius) {
            for &t in self.in_cell(cell) {
                let along = self.smooth[t as usize].dot(d);
                if along.abs() >= limit {
                    out.push((t, along));
                }
            }
        }
        out
    }

    /// The direction the area gathers at from a cell: the mean of the
    /// smoothed normals within the resolution, aligned, taken again from
    /// the mean until it stops moving; with the area within it.
    fn mode(&self, cell: u64, tol: Tolerances) -> Option<(Vector, f64)> {
        let centre = self.grid.centre(cell);
        let mean = |members: &[(u32, f64)]| -> Option<(Vector, f64)> {
            let mut sum = Vector::ZERO;
            let mut area = 0.0;
            for &(t, along) in members {
                let a = self.facets.area[t as usize];
                let s = self.smooth[t as usize];
                sum += if along < 0.0 { -s } else { s } * a;
                area += a;
            }
            let m = sum.magnitude();
            (m > 0.0).then(|| (canonical(sum / m), area))
        };
        let seed: Vec<(u32, f64)> = self
            .in_cell(cell)
            .iter()
            .map(|&t| (t, self.smooth[t as usize].dot(centre)))
            .collect();
        let (mut d, mut area) = mean(&seed)?;
        // Mean shift with a flat window settles in finitely many steps, the
        // area within the window growing with each; the bound only stops a
        // step that rounding sends back and forth between two sets.
        for _ in 0..64 {
            let (next, gathered) = mean(&self.near(d, self.resolution))?;
            let moved = next.cross(d).magnitude().atan2(next.dot(d).abs());
            d = next;
            area = gathered;
            if moved <= tol.angular() {
                break;
            }
        }
        Some((d, area))
    }

    /// The directions kept from the modes, in the order they were kept.
    fn kept_directions(&self, tol: Tolerances) -> Vec<Vector> {
        let area_of = |cell: u64| -> f64 {
            self.in_cell(cell)
                .iter()
                .map(|&t| self.facets.area[t as usize])
                .sum()
        };
        let areas: Vec<f64> = self.smoothed.cells.iter().map(|&c| area_of(c)).collect();
        let above = |p: (f64, Vector), q: (f64, Vector)| -> bool {
            match p.0.total_cmp(&q.0) {
                Ordering::Greater => true,
                Ordering::Less => false,
                Ordering::Equal => lexical(p.1, q.1) == Ordering::Less,
            }
        };
        let seeds: Vec<u64> = self
            .smoothed
            .cells
            .iter()
            .zip(&areas)
            .filter(|&(&cell, &area)| {
                let own = (area, self.grid.centre(cell));
                self.grid.neighbours(cell).all(|n| {
                    let theirs = match self.smoothed.cells.binary_search(&n) {
                        Ok(k) => areas[k],
                        Err(_) => 0.0,
                    };
                    above(own, (theirs, self.grid.centre(n)))
                })
            })
            .map(|(&cell, _)| cell)
            .collect();
        let mut modes: Vec<(f64, Vector)> =
            parallel::map_ordered(&seeds, |_, &cell| self.mode(cell, tol))
                .into_iter()
                .flatten()
                .map(|(d, area)| (area, d))
                .collect();
        modes.sort_by(|p, q| q.0.total_cmp(&p.0).then_with(|| lexical(p.1, q.1)));
        let one = self.resolution.cos();
        let mut kept: Vec<Vector> = Vec::new();
        for (_, d) in modes {
            if kept.iter().any(|k| k.dot(d).abs() >= one) {
                continue;
            }
            let mut inner = 0.0;
            let mut ring = 0.0;
            for (t, along) in self.near(d, 3.0 * self.resolution) {
                let s = self.smooth[t as usize];
                if kept.iter().any(|k| k.dot(s).abs() >= one) {
                    continue;
                }
                let a = self.facets.area[t as usize];
                if along.abs() >= one {
                    inner += a;
                } else {
                    ring += a;
                }
            }
            if inner > ring {
                kept.push(d);
            }
        }
        kept
    }

    /// The face directions of the kept directions: each triangle whose own
    /// or smoothed normal is within the resolution of one belongs to the
    /// one either is nearest, the first kept of equally near ones; each
    /// keeps those of its triangles the cutoff keeps, and is measured from
    /// their own normals.
    fn face_directions(&self, kept: &[Vector], tol: Tolerances) -> Vec<FaceDirection> {
        let mut owner: Vec<(usize, f64)> = vec![(usize::MAX, -1.0); self.smooth.len()];
        for (k, d) in kept.iter().enumerate() {
            let near = self.near(*d, self.resolution);
            for (t, along) in near.into_iter().chain(self.near_own(*d, self.resolution)) {
                let slot = &mut owner[t as usize];
                if along.abs() > slot.1 {
                    *slot = (k, along.abs());
                }
            }
        }
        let mut owned: Vec<Vec<u32>> = vec![Vec::new(); kept.len()];
        for (t, &(k, _)) in owner.iter().enumerate() {
            if k != usize::MAX {
                owned[k].push(as_index(t));
            }
        }
        let mut out: Vec<FaceDirection> = owned
            .iter()
            .zip(kept)
            .filter_map(|(triangles, &d)| {
                let triangles = self.cut(triangles, d);
                self.measured(&triangles, d, tol)
            })
            .collect();
        out.sort_by(|p, q| {
            q.area.total_cmp(&p.area).then_with(|| {
                lexical(
                    canonical(p.direction.vector()),
                    canonical(q.direction.vector()),
                )
            })
        });
        out
    }

    /// The area-weighted mean of the triangles' own normals, each turned
    /// to agree with `about`.
    fn mean(&self, triangles: &[u32], about: Vector) -> Option<Vector> {
        let mut sum = Vector::ZERO;
        for &t in triangles {
            let n = self.facets.normal[t as usize];
            let a = self.facets.area[t as usize];
            sum += if n.dot(about) < 0.0 { -n } else { n } * a;
        }
        let m = sum.magnitude();
        (m > 0.0).then(|| sum / m)
    }

    /// The triangles the cutoff keeps of those a direction owns. One whose
    /// own normal is within the resolution of the mean is always kept:
    /// directions that near are one. One owned by its smoothed normal alone
    /// is kept where its own normal is within the cutoff's threshold of the
    /// mean, measured over the angles of the triangles kept, weighed by
    /// area: a noisy face's own normals scatter past the resolution and the
    /// threshold reaches out with them, while an exact face's do not
    /// scatter at all, and a fillet's triangle that smoothing drew towards
    /// it is left out. Each round takes the mean of the last round's kept
    /// triangles, and the rounds repeat until the kept set stops changing.
    fn cut(&self, owned: &[u32], about: Vector) -> Vec<u32> {
        let angle_from = |d: Vector| {
            move |t: u32| {
                let n = self.facets.normal[t as usize];
                n.cross(d).magnitude().atan2(n.dot(d).abs())
            }
        };
        let within: Vec<u32> = owned
            .iter()
            .copied()
            .filter(|&t| angle_from(about)(t) <= self.resolution)
            .collect();
        let mut kept = if within.is_empty() {
            owned.to_vec()
        } else {
            within
        };
        // The kept set only changes while the threshold does; the bound
        // stops a set that rounding sends back and forth between two.
        for _ in 0..64 {
            let Some(d) = self.mean(&kept, about) else {
                break;
            };
            let angle = angle_from(d);
            let weighed: Vec<(f64, f64)> = kept
                .iter()
                .map(|&t| (angle(t), self.facets.area[t as usize]))
                .collect();
            let limit = self
                .cutoff
                .threshold_weighted(&weighed)
                .max(self.resolution);
            let next: Vec<u32> = owned
                .iter()
                .copied()
                .filter(|&t| angle(t) <= limit)
                .collect();
            if next.is_empty() || next == kept {
                break;
            }
            kept = next;
        }
        kept
    }

    /// A face direction measured from its triangles' own normals.
    fn measured(&self, triangles: &[u32], about: Vector, tol: Tolerances) -> Option<FaceDirection> {
        let d = canonical(self.mean(triangles, about)?);
        let area: f64 = triangles
            .iter()
            .map(|&t| self.facets.area[t as usize])
            .sum();
        let mut net = 0.0;
        let mut square = 0.0;
        for &t in triangles {
            let n = self.facets.normal[t as usize];
            let a = self.facets.area[t as usize];
            let angle = n.cross(d).magnitude().atan2(n.dot(d).abs());
            net += a * n.dot(d);
            square += a * angle * angle;
        }
        let d = if net < -tol.confusion() * area.sqrt() {
            -d
        } else {
            d
        };
        Some(FaceDirection {
            direction: Direction::new(d, tol).ok()?,
            area,
            spread: (square / area).sqrt(),
        })
    }

    /// The axis found from `start`, with the wall area about it.
    fn revolution_from(&self, start: Vector, tol: Tolerances) -> Option<(Vector, Point, f64)> {
        let sine = self.resolution.sin();
        let walls = |a: Vector| -> Vec<u32> {
            (0..self.smooth.len())
                .filter(|&t| self.facets.area[t] > 0.0 && self.smooth[t].dot(a).abs() <= sine)
                .map(as_index)
                .collect()
        };
        let about = |a: Vector, p: Point| -> Vec<u32> {
            walls(a)
                .into_iter()
                .filter(|&t| {
                    let w = self.facets.centroid[t as usize] - p;
                    let radial = w - a * w.dot(a);
                    let r = radial.magnitude();
                    r > 0.0 && self.smooth[t as usize].dot(a.cross(radial) / r).abs() <= sine
                })
                .collect()
        };
        let mut a = start;
        let mut chosen = walls(a);
        let mut p = self.located(a, &chosen, tol)?;
        // Each round's walls are the triangles pointing at the last round's
        // axis; the bound only stops a set that rounding sends back and
        // forth.
        for _ in 0..64 {
            let next = about(a, p);
            if next.is_empty() {
                return None;
            }
            a = self.least_normal(&next)?;
            p = self.located(a, &next, tol)?;
            if next == chosen {
                break;
            }
            chosen = next;
        }
        let mut area = 0.0;
        let mut weighted = Vector::ZERO;
        for &t in &chosen {
            let w = self.facets.area[t as usize];
            area += w;
            weighted += (self.facets.centroid[t as usize] - p) * w;
        }
        let level = p + a * (weighted / area).dot(a);
        Some((a, level, area))
    }

    /// The direction the triangles' normals are most nearly perpendicular
    /// to.
    fn least_normal(&self, triangles: &[u32]) -> Option<Vector> {
        let m = covariance(
            triangles
                .iter()
                .map(|&t| self.facets.normal[t as usize] * self.facets.area[t as usize].sqrt()),
        );
        let d = smallest_direction(m);
        let length = d.magnitude();
        (length > 0.0).then(|| canonical(d / length))
    }

    /// The point, square to `a`, of the line along `a` the triangles'
    /// normal lines come nearest meeting, by area-weighted least squares;
    /// `None` where their normals do not span the plane square to `a`.
    fn located(&self, a: Vector, triangles: &[u32], tol: Tolerances) -> Option<Point> {
        let reference = self.facets.centroid[*triangles.first()? as usize];
        let e1 = Direction::new(a, tol).ok()?.any_perpendicular().vector();
        let e2 = a.cross(e1);
        let (mut m11, mut m12, mut m22, mut b1, mut b2) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for &t in triangles {
            let w = self.facets.area[t as usize];
            let m = a.cross(self.facets.normal[t as usize]);
            let (x, y) = (m.dot(e1), m.dot(e2));
            let rhs = (self.facets.centroid[t as usize] - reference).dot(m);
            m11 += w * x * x;
            m12 += w * x * y;
            m22 += w * y * y;
            b1 += w * x * rhs;
            b2 += w * y * rhs;
        }
        let det = m11.mul_add(m22, -m12 * m12);
        let trace = m11 + m22;
        let angular = tol.angular();
        if det.is_nan() || det <= angular * angular * trace * trace {
            return None;
        }
        let x = (b1 * m22 - b2 * m12) / det;
        let y = (b2 * m11 - b1 * m12) / det;
        Some(reference + e1 * x + e2 * y)
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn every_direction_lands_in_a_cell_its_antipode_shares_and_its_centre_is_near() {
        let grid = Grid::new(5.0_f64.to_radians());
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut unit = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            #[allow(clippy::cast_precision_loss, reason = "test draws")]
            let u = (seed >> 11) as f64 / (1_u64 << 53) as f64;
            2.0f64.mul_add(u, -1.0)
        };
        for _ in 0..20_000 {
            let v = Vector::new(unit(), unit(), unit());
            if v.magnitude() < 1e-3 {
                continue;
            }
            let d = v / v.magnitude();
            let cell = grid.cell(d);
            assert_eq!(cell, grid.cell(-d));
            let c = grid.centre(cell);
            let angle = c.cross(d).magnitude().atan2(c.dot(d).abs());
            assert!(angle <= grid.reach(), "{angle} past {}", grid.reach());
        }
    }

    #[test]
    fn a_resolution_out_of_range_is_refused() {
        for resolution in [0.0, -0.1, f64::NAN, 0.6] {
            let options = Directions {
                resolution,
                ..Directions::default()
            };
            assert!(options.check().is_err());
        }
        assert!(
            Directions {
                crease: 0.0,
                ..Directions::default()
            }
            .check()
            .is_err()
        );
        let bad = Triangulation {
            positions: vec![Point::ORIGIN],
            triangles: vec![[0, 1, 2]],
            ..Triangulation::new()
        };
        assert!(face_directions(&bad, &Directions::default(), T).is_err());
        assert!(open_rims(&bad).is_err());
    }

    #[test]
    fn an_empty_mesh_has_no_directions_no_axis_and_no_rims() {
        let empty = Triangulation::new();
        assert!(
            face_directions(&empty, &Directions::default(), T)
                .unwrap()
                .is_empty()
        );
        assert!(
            axis_of_revolution(&empty, &Directions::default(), T)
                .unwrap()
                .is_none()
        );
        assert!(open_rims(&empty).unwrap().is_empty());
        assert!(principal_frame(&empty, T).is_err());
    }
}
