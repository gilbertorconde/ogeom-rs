//! Offsetting an editable mesh: every face moved along its outward normal
//! by one distance, the folds that leaves repaired.
//!
//! # Placing a vertex
//!
//! A face's outward normal is the way its corners wind counter-clockwise.
//! About each vertex, the faces of its fan are taken in turn, and an edge
//! between two of them whose normals turn by more than [`OFFSET_EDGE_ANGLE`]
//! (30°) is an *edge* of the offset. The edges, and the rims, part the fan
//! into sectors; a closed fan with fewer than two edges is one sector. Each
//! sector's normal is the mean of its faces' unit normals weighted by the
//! faces' angles at the vertex: the normals smoothed within the sector and
//! never across an edge.
//!
//! The sectors' normals fix the directions the vertex is held in: the
//! heaviest sector's normal, then, of the sectors whose normals turn by
//! more than the edge angle from each normal taken, the one whose normal
//! has the largest part square to the directions taken, that part's
//! direction, up to three. Within the directions taken, the
//! vertex moves by the step `x` that minimises `Σ α (n · x − d)²` over its
//! faces, `n` a face's unit normal and `α` its angle at the vertex: least
//! squares weighted by corner angle, so that each face ends as near
//! `distance` `d` from where it was as the vertex can put it. A vertex on a
//! flat stretch moves by `d` along the normal, one on an edge between faces
//! square to each other by `√2 d`, and the corner of a box by `√3 d`, to
//! where the offset faces' planes meet. No vertex moves farther than
//! [`OFFSET_FARTHEST`] (3) times `|d|`: a longer step is shortened along
//! its line, as computed.
//!
//! # Folds
//!
//! Where the offset faces' planes meet before the faces end, as beside a
//! convex edge offset inward or a concave one offset outward, the vertices
//! within reach of the edge overshoot it and their faces fold over. Folds
//! are repaired in two steps.
//!
//! First, overshooting vertices are pushed back. A vertex moved in one or
//! two directions that ends nearer than `|d|` to the mesh as it was, to a
//! face it stands in front of, on the side the offset goes, and that turns
//! by at least the edge angle from the plane or line the vertex may still
//! slide in, slides in it onto that face's offset plane, and is then held
//! to that plane too. Faces are looked for about where the vertex was
//! placed, the nearest first, and again for a second push along the line
//! the first leaves. A slide that would take the vertex farther than the
//! cap from where it started is not made. A strip that overshot an edge
//! so lands on the edge, and a vertex in reach of two edges on the corner
//! where they meet. Faces of one smooth stretch never push each other, so
//! the residuals least squares leaves on a curved surface stay; nor do
//! faces the vertex stands behind, so the offset of a piece thinner than
//! twice the distance, whose sides pass through each other, is not
//! repaired. These vertices are counted in [`OffsetReport::eased`].
//!
//! Then folded faces are collapsed. A face is *folded* when it turns from
//! the way it faced before the offset by more than the edge angle, over
//! included, or is a *sliver*, its height within the rounding of its
//! corners' coordinates; a face that was a sliver before is passed over.
//! Each folded face, lowest-numbered first, has an edge collapsed: one of
//! its corners merged into another, which keeps its position. The corner
//! merged is the one held in the fewest directions, so a vertex of a flat
//! stretch merges into one on an edge and that into a corner; ties go to
//! the highest-numbered corner. It merges into the corner held in the most
//! directions, ties to the nearest, then the lowest-numbered. A collapse is
//! made only where it keeps the mesh a manifold, the two vertices sharing
//! no neighbour but the corners across the edge, no face about one closing
//! over the edge with a face about the other, and a vertex on a rim
//! merging only along the rim; and only where it turns no face it moves by
//! more than the edge angle short of a sliver. The faces about the vertex
//! merged into are then looked at again. A strip pushed onto an edge so
//! merges into it, and the faces left lie on the offset planes.
//!
//! Each collapse takes a vertex away, so the repair ends. A folded face no
//! collapse may repair is left, and counted in [`OffsetReport::folded`].
//!
//! # A region
//!
//! With a region, the faces whose corners are all in it move and the rest
//! stay. A vertex the region shares with faces outside it is parted in
//! two: the vertex keeps its number and moves with the region's faces,
//! placed from them alone, and a new vertex, with its source, stays with
//! the faces outside. Each edge between a face of the region and one
//! outside becomes a wall of two faces joining the two copies of its ends,
//! so the region moves as a step and a closed mesh stays closed.
//!
//! # What stays
//!
//! Vertex and face numbers stay, but for the vertices merged, which are
//! gone, the faces collapsed, which are removed, and the vertices and faces
//! a region's wall adds, numbered after every one there is. A vertex with
//! no face does not move. A vertex [`Editable::non_manifold`] reports moves
//! with its own fan like any other, so the copies of a point several fans
//! meet part, but is never merged; edges of three or more faces are joined
//! again as they then stand.

use core::f64::consts::FRAC_PI_6;
use std::collections::BTreeSet;

use super::smooth::capped;
use super::{Editable, NONE, Scope, index, prev};
use crate::bvh::Bvh;
use crate::curvature::in_parallel;
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;

/// The angle, in radians, beyond which the normals of two faces meeting at
/// an edge turn for the edge to be one the offset keeps: 30°.
pub const OFFSET_EDGE_ANGLE: f64 = FRAC_PI_6;

/// The farthest a vertex moves, as a multiple of the offset distance.
pub const OFFSET_FARTHEST: f64 = 3.0;

/// A face's height, in units of the rounding of its corners' largest
/// coordinate, at or below which it is a sliver. A vertex is placed and
/// pushed by solving for the meeting of at most three planes, each pair
/// turning by at least the edge angle, so each solve loses at most
/// `1 / sin 30° = 2` of a few units of rounding: sixteen covers three.
const SLIVER: f64 = 16.0;

/// What an offset did besides moving the vertices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OffsetReport {
    /// How many vertices were merged into a neighbour to repair folds.
    pub merged: usize,
    /// How many vertices that overshot an edge were pushed back onto it.
    pub eased: usize,
    /// How many vertices' steps were shortened to [`OFFSET_FARTHEST`] times
    /// the distance.
    pub capped: usize,
    /// How many faces are left folded, no collapse being able to repair
    /// them while keeping the mesh a manifold.
    pub folded: usize,
}

/// Where a vertex goes: its step, and how many directions it is held in.
#[derive(Clone, Copy)]
struct Placed {
    step: Vector,
    rank: u8,
    /// For a vertex moved in one direction, that direction, square to the
    /// plane it may still slide in; in two, the line it may slide along.
    free: Vector,
    capped: bool,
}

/// A sector of a vertex's fan: its faces' unit normals and angles there.
#[derive(Default)]
struct Sector {
    faces: Vec<(Vector, f64)>,
}

impl Sector {
    fn weight(&self) -> f64 {
        self.faces.iter().map(|&(_, a)| a).sum()
    }

    fn normal(&self) -> Option<Vector> {
        let sum = self
            .faces
            .iter()
            .fold(Vector::new(0.0, 0.0, 0.0), |s, &(n, a)| s + n * a);
        unit(sum)
    }
}

fn unit(v: Vector) -> Option<Vector> {
    let length = v.magnitude();
    (length > 0.0 && length.is_finite()).then(|| v / length)
}

impl Editable {
    /// Offset the mesh by `distance` along its outward normals, outward
    /// where positive and inward where negative, and repair the folds
    /// that leaves, as the [module documentation](crate::editable::offset)
    /// describes. With `region`, one entry per vertex number, only the
    /// faces whose corners are all in it move, as a step joined to the
    /// rest by a wall.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if
    /// `distance` is not finite;
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if `region` does
    /// not have one entry per vertex number, or the wall's vertices and
    /// faces would pass `u32` numbering.
    pub fn offset(&mut self, distance: f64, region: Option<&[bool]>) -> OgeomResult<OffsetReport> {
        if !distance.is_finite() {
            ogeom_bail!(Numeric, "an offset of {distance} is not a distance");
        }
        let n = self.positions.len();
        if let Some(region) = region
            && region.len() != n
        {
            ogeom_bail!(
                Range,
                "a region of {} entries for {n} vertices",
                region.len()
            );
        }
        let mut report = OffsetReport {
            merged: 0,
            eased: 0,
            capped: 0,
            folded: 0,
        };
        if distance == 0.0 {
            return Ok(report);
        }
        let faces_before = self.corners.len();
        let area: Vec<Vector> = (0..faces_before).map(|f| self.area_vector(f)).collect();
        let normal: Vec<Option<Vector>> = area
            .iter()
            .zip(&self.live)
            .map(|(&a, &live)| if live { unit(a) } else { None })
            .collect();
        let was_sliver: Vec<bool> = (0..faces_before)
            .map(|f| {
                let [a, b, c] = self.corners[f].map(|v| self.positions[v as usize]);
                shape(a, b, c, area[f]) == Shape::Sliver
            })
            .collect();
        let moves: Vec<bool> = self
            .corners
            .iter()
            .zip(&self.live)
            .map(|(c, &live)| live && region.is_none_or(|r| c.iter().all(|&v| r[v as usize])))
            .collect();

        let candidates: Vec<u32> = (0..n)
            .filter(|&v| self.out[v] != NONE && region.is_none_or(|r| r[v]))
            .map(index)
            .collect();
        let placed = in_parallel(&candidates, |v| self.place(v, distance, &normal, &moves));
        if region.is_some() {
            self.raise_walls(&candidates, &moves)?;
        }

        let before = self.positions.clone();
        let stuck = self.non_manifold();
        let mut rank = vec![0_u8; self.positions.len()];
        let mut movable = vec![false; self.positions.len()];
        let mut sliding: Vec<(u32, Placed)> = Vec::new();
        for (&v, p) in candidates.iter().zip(&placed) {
            let k = v as usize;
            if p.rank > 0 {
                self.positions[k] += p.step;
                rank[k] = p.rank;
                movable[k] = true;
                report.capped += usize::from(p.capped);
                if p.rank < 3 {
                    sliding.push((v, *p));
                }
            }
        }
        for v in stuck {
            movable[v as usize] = false;
        }
        self.rejoin_crowded_edges();

        let pushed = self.push_back(&before, faces_before, &normal, distance, &sliding)?;
        report.eased = pushed.len();
        for (v, to, r) in pushed {
            self.positions[v as usize] = to;
            rank[v as usize] = r;
        }

        let folded = |me: &Self, f: usize| {
            f < faces_before && me.live[f] && !was_sliver[f] && me.folded(f, area[f])
        };
        let mut work: BTreeSet<u32> = (0..faces_before)
            .filter(|&f| folded(self, f))
            .map(index)
            .collect();
        while let Some(f) = work.pop_first() {
            if !folded(self, f as usize) {
                continue;
            }
            if let Some(w) = self.repair(f, &rank, &movable, &area) {
                report.merged += 1;
                work.extend(
                    self.fan(w)
                        .into_iter()
                        .map(|h| h / 3)
                        .filter(|&g| (g as usize) < faces_before),
                );
            }
        }
        report.folded = (0..faces_before).filter(|&f| folded(self, f)).count();
        Ok(report)
    }

    /// Face `f`'s area vector, twice its area along its normal.
    fn area_vector(&self, f: usize) -> Vector {
        let [a, b, c] = self.corners[f].map(|v| self.positions[v as usize]);
        (b - a).cross(c - a)
    }

    /// The half-edges leaving vertex `v`, turning about it from its `out`.
    fn fan(&self, v: u32) -> Vec<u32> {
        let start = self.out[v as usize];
        let mut fan = Vec::new();
        if start == NONE {
            return fan;
        }
        let mut h = start;
        loop {
            fan.push(h);
            h = self.twin[prev(h) as usize];
            if h == NONE || h == start {
                return fan;
            }
        }
    }

    /// Where vertex `v` goes, offset `d` from the faces `moves` marks.
    fn place(&self, v: u32, d: f64, normal: &[Option<Vector>], moves: &[bool]) -> Placed {
        let still = Placed {
            step: Vector::new(0.0, 0.0, 0.0),
            rank: 0,
            free: Vector::new(0.0, 0.0, 0.0),
            capped: false,
        };
        let fan = self.fan(v);
        let closed = fan.first().is_some_and(|&h| self.twin[h as usize] != NONE);
        let p = self.positions[v as usize];
        let cos_edge = OFFSET_EDGE_ANGLE.cos();
        let mut sectors: Vec<Sector> = vec![Sector::default()];
        // Whether the fan was parted where it starts: always for an open
        // fan, and for a closed one by an edge or a face left out there.
        let mut parted_at_start = !closed;
        let mut last: Option<Vector> = None;
        let mut first: Option<Vector> = None;
        for (i, &h) in fan.iter().enumerate() {
            let f = (h / 3) as usize;
            if !moves[f] {
                if i == 0 {
                    parted_at_start = true;
                }
                if !sectors.last().is_some_and(|s| s.faces.is_empty()) {
                    sectors.push(Sector::default());
                }
                last = None;
                continue;
            }
            let Some(nf) = normal[f] else { continue };
            if last.is_some_and(|m| m.dot(nf) < cos_edge)
                && !sectors.last().is_some_and(|s| s.faces.is_empty())
            {
                sectors.push(Sector::default());
            }
            if first.is_none() && sectors.len() == 1 {
                first = Some(nf);
            }
            last = Some(nf);
            let [a, b] = [self.head(h), self.tail(prev(h))].map(|w| self.positions[w as usize] - p);
            let angle = a.cross(b).magnitude().atan2(a.dot(b));
            if let Some(s) = sectors.last_mut() {
                s.faces.push((nf, angle));
            }
        }
        sectors.retain(|s| !s.faces.is_empty());
        if closed
            && !parted_at_start
            && sectors.len() > 1
            && let (Some(m), Some(nf)) = (last, first)
            && m.dot(nf) >= cos_edge
            && let Some(tail) = sectors.pop()
        {
            sectors[0].faces.extend(tail.faces);
        }

        let groups: Vec<(Vector, f64)> = sectors
            .iter()
            .filter_map(|s| s.normal().map(|m| (m, s.weight())))
            .collect();
        let Some(&(heaviest, _)) = groups.iter().reduce(|a, b| if b.1 > a.1 { b } else { a })
        else {
            return still;
        };
        let mut basis = vec![heaviest];
        let mut taken = vec![heaviest];
        while basis.len() < 3 {
            let across = groups
                .iter()
                .filter(|&&(m, _)| taken.iter().all(|&t| m.dot(t) < cos_edge))
                .map(|&(m, _)| (m, basis.iter().fold(m, |r, &e| r - e * r.dot(e))))
                .reduce(|a, b| {
                    if b.1.magnitude() > a.1.magnitude() {
                        b
                    } else {
                        a
                    }
                });
            match across.and_then(|(m, r)| unit(r).map(|e| (m, e))) {
                Some((m, e)) => {
                    taken.push(m);
                    basis.push(e);
                }
                None => break,
            }
        }

        let faces = sectors.iter().flat_map(|s| s.faces.iter());
        let k = basis.len();
        let mut m = [[0.0; 3]; 3];
        let mut r = [0.0; 3];
        for &(nf, a) in faces {
            let mut c = [0.0; 3];
            for (ci, &e) in c.iter_mut().zip(&basis) {
                *ci = nf.dot(e);
            }
            for i in 0..k {
                r[i] += a * d * c[i];
                for j in 0..k {
                    m[i][j] += a * c[i] * c[j];
                }
            }
        }
        let Some(alpha) = solve(m, r, k) else {
            return still;
        };
        let step = basis
            .iter()
            .zip(alpha)
            .fold(Vector::new(0.0, 0.0, 0.0), |s, (&e, t)| s + e * t);
        if !step.magnitude().is_finite() {
            return still;
        }
        let (to, held) = capped(p, p + step, OFFSET_FARTHEST * d.abs());
        let free = match k {
            1 => Some(basis[0]),
            2 => unit(basis[0].cross(basis[1])),
            _ => None,
        };
        Placed {
            step: to - p,
            rank: u8::try_from(k).unwrap_or(u8::MAX),
            free: free.unwrap_or(Vector::new(0.0, 0.0, 0.0)),
            capped: held,
        }
    }

    /// Part each vertex among `candidates` that has faces both moving and
    /// not, as `moves` marks, into itself on the moving faces and a new
    /// vertex on the rest, and join the two by a wall along every edge
    /// between a moving face and one that is not.
    fn raise_walls(&mut self, candidates: &[u32], moves: &[bool]) -> OgeomResult<()> {
        let mut foot = vec![NONE; self.positions.len()];
        let mut parted: Vec<u32> = Vec::new();
        for &v in candidates {
            let fan = self.fan(v);
            let mixed = fan.iter().any(|&h| moves[(h / 3) as usize])
                && fan.iter().any(|&h| !moves[(h / 3) as usize]);
            if mixed {
                parted.push(v);
            }
        }
        if parted.is_empty() {
            return Ok(());
        }
        let walls: Vec<u32> = parted
            .iter()
            .flat_map(|&v| self.fan(v))
            .filter(|&h| {
                let g = self.twin[h as usize];
                moves[(h / 3) as usize] && g != NONE && !moves[(g / 3) as usize]
            })
            .collect();
        let added = parted.len().saturating_add(self.positions.len());
        if added >= NONE as usize || (self.corners.len() + 2 * walls.len()) * 3 >= NONE as usize {
            ogeom_bail!(Range, "a wall past u32 numbering");
        }
        let mut points: Vec<u32> = parted.iter().map(|&v| self.point[v as usize]).collect();
        points.sort_unstable();
        points.dedup();
        let mut faces = self.about(&points);
        for &v in &parted {
            let w = self.add_vertex(self.positions[v as usize]);
            self.source[w as usize] = self.source[v as usize];
            foot[v as usize] = w;
            points.push(w);
        }
        for &f in &faces {
            if !moves[f as usize] {
                for c in &mut self.corners[f as usize] {
                    let w = foot[*c as usize];
                    if w != NONE {
                        *c = w;
                    }
                }
            }
        }
        for h in walls {
            let (a, b) = (self.tail(h), self.head(h));
            let (fa, fb) = (foot[a as usize], foot[b as usize]);
            faces.push(self.add_face([fa, fb, b]));
            faces.push(self.add_face([fa, b, a]));
        }
        faces.sort_unstable();
        points.sort_unstable();
        self.relink(&Scope::Points(points), &faces);
        Ok(())
    }

    /// Push back each vertex of `sliding` that has come nearer than
    /// `|d|` to the mesh as it stood at `before`, its first
    /// `faces_before` faces with their unit normals `normal`: slid along
    /// the plane or line it may still move in onto the offset plane of the
    /// face it came nearest, then held to that plane too, while it stays on
    /// the same side of that plane and within the cap. Returns each vertex
    /// pushed, where it goes, and in how many directions it is then held.
    fn push_back(
        &self,
        before: &[Point],
        faces_before: usize,
        normal: &[Option<Vector>],
        d: f64,
        sliding: &[(u32, Placed)],
    ) -> OgeomResult<Vec<(u32, Point, u8)>> {
        if sliding.is_empty() {
            return Ok(Vec::new());
        }
        let mut faces: Vec<u32> = Vec::new();
        let mut triangles: Vec<[u32; 3]> = Vec::new();
        for (f, n) in normal.iter().enumerate().take(faces_before) {
            if n.is_some() {
                faces.push(index(f));
                triangles.push(self.corners[f]);
            }
        }
        let mesh = Triangulation {
            positions: before.to_vec(),
            triangles,
            ..Triangulation::new()
        };
        let bvh = Bvh::build(&mesh)?;
        let reach = d.abs();
        let least = OFFSET_EDGE_ANGLE.sin();
        let movers: Vec<u32> = sliding.iter().map(|&(v, _)| v).collect();
        let pushed = in_parallel(&movers, |v| {
            let slot = movers.binary_search(&v).ok()?;
            let placed = sliding[slot].1;
            let start = before[v as usize];
            let placed_at = self.positions[v as usize];
            let mut x = placed_at;
            let mut rank = placed.rank;
            let mut free = placed.free;
            let mut moved = false;
            while rank < 3 {
                // The face must turn by the edge angle from the plane or
                // line the vertex slides in, or it is of the same smooth
                // stretch.
                let across = |nf: Vector| {
                    let slant = if rank == 1 {
                        (1.0 - nf.dot(free).powi(2)).max(0.0).sqrt()
                    } else {
                        nf.dot(free).abs()
                    };
                    slant >= least
                };
                // Faces are looked for about where the vertex was placed,
                // and it must stand there on the side of the face it is
                // offset to.
                let side = |f: usize, nf: Vector, at: Point| {
                    nf.dot(at - before[self.corners[f][0] as usize])
                };
                let Some(hit) = bvh.nearest_where(placed_at, reach, |t| {
                    let f = faces[t as usize] as usize;
                    normal[f].is_some_and(|nf| across(nf) && side(f, nf, placed_at) * d > 0.0)
                }) else {
                    break;
                };
                if hit.distance >= reach {
                    break;
                }
                let f = faces[hit.triangle as usize] as usize;
                let Some(nf) = normal[f] else { break };
                let side = side(f, nf, x);
                let along = if rank == 1 {
                    nf - free * nf.dot(free)
                } else {
                    free
                };
                let Some(along) = unit(along) else { break };
                let t = (d - side) / nf.dot(along);
                let to = x + along * t;
                if !t.is_finite() || (to - start).magnitude() > OFFSET_FARTHEST * reach {
                    break;
                }
                x = to;
                moved = true;
                if rank == 1 {
                    let Some(line) = unit(free.cross(nf)) else {
                        break;
                    };
                    free = line;
                }
                rank += 1;
            }
            moved.then_some((v, x, rank))
        });
        Ok(pushed.into_iter().flatten().collect())
    }

    /// Whether face `f`, its area vector `was` before the offset, is
    /// folded.
    fn folded(&self, f: usize, was: Vector) -> bool {
        let [a, b, c] = self.corners[f].map(|v| self.positions[v as usize]);
        shape(a, b, c, was) != Shape::Kept
    }

    /// Repair folded face `f` by merging one of its corners into another,
    /// as the module documentation orders them; the vertex merged into, or
    /// `None` where no collapse keeps the mesh a manifold.
    fn repair(&mut self, f: u32, rank: &[u8], movable: &[bool], area: &[Vector]) -> Option<u32> {
        let corners = self.corners[f as usize];
        let mut from: Vec<u32> = corners
            .iter()
            .copied()
            .filter(|&v| movable[v as usize])
            .collect();
        from.sort_by(|&a, &b| rank[a as usize].cmp(&rank[b as usize]).then(b.cmp(&a)));
        for &u in &from {
            let p = self.positions[u as usize];
            let mut to: Vec<u32> = from.iter().copied().filter(|&w| w != u).collect();
            to.sort_by(|&a, &b| {
                rank[b as usize]
                    .cmp(&rank[a as usize])
                    .then(
                        self.positions[a as usize]
                            .distance(p)
                            .total_cmp(&self.positions[b as usize].distance(p)),
                    )
                    .then(a.cmp(&b))
            });
            for w in to {
                if self.collapse(u, w, area) {
                    return Some(w);
                }
            }
        }
        None
    }

    /// Merge vertex `u` into its neighbour `w`, removing the faces on the
    /// edge between them, if that keeps the mesh a manifold; whether it did.
    fn collapse(&mut self, u: u32, w: u32, area: &[Vector]) -> bool {
        let lone = |me: &Self, v: u32| {
            me.ring(me.point[v as usize])
                .into_iter()
                .all(|c| c == v || me.out[c as usize] == NONE)
        };
        if !lone(self, u) || !lone(self, w) {
            return false;
        }
        let (fan_u, fan_w) = (self.fan(u), self.fan(w));
        let on_rim = |fan: &[u32]| fan.first().is_some_and(|&h| self.twin[h as usize] == NONE);
        let neighbours = |fan: &[u32]| {
            let mut n: Vec<u32> = fan
                .iter()
                .flat_map(|&h| [self.head(h), self.tail(prev(h))])
                .collect();
            n.sort_unstable();
            n.dedup();
            n
        };
        let link = |fan: &[u32]| {
            let mut l: Vec<(u32, u32)> = fan
                .iter()
                .map(|&h| {
                    let (a, b) = (self.head(h), self.tail(prev(h)));
                    (a.min(b), a.max(b))
                })
                .collect();
            l.sort_unstable();
            l
        };
        let shared: Vec<u32> = fan_u
            .iter()
            .copied()
            .filter(|&h| self.corners[(h / 3) as usize].contains(&w))
            .collect();
        if shared.is_empty() || shared.len() > 2 {
            return false;
        }
        if on_rim(&fan_u) && shared.len() != 1 {
            return false;
        }
        let mut across: Vec<u32> = shared
            .iter()
            .flat_map(|&h| self.corners[(h / 3) as usize])
            .filter(|&c| c != u && c != w)
            .collect();
        across.sort_unstable();
        across.dedup();
        let (nu, nw) = (neighbours(&fan_u), neighbours(&fan_w));
        let common: Vec<u32> = nu
            .iter()
            .copied()
            .filter(|c| nw.binary_search(c).is_ok())
            .collect();
        if common != across {
            return false;
        }
        let lw = link(&fan_w);
        if link(&fan_u).iter().any(|e| lw.binary_search(e).is_ok()) {
            return false;
        }

        let p = self.positions[w as usize];
        for &h in &fan_u {
            let f = (h / 3) as usize;
            if f >= area.len() || shared.contains(&h) {
                continue;
            }
            let [a, b] = [self.head(h), self.tail(prev(h))].map(|v| self.positions[v as usize]);
            if shape(p, a, b, area[f]) == Shape::Turned {
                return false;
            }
        }

        let mut points: Vec<u32> = [u, w]
            .into_iter()
            .chain(across.iter().copied())
            .map(|v| self.point[v as usize])
            .collect();
        points.sort_unstable();
        points.dedup();
        let around = self.about(&points);
        for &h in &shared {
            let g = h / 3;
            self.live[g as usize] = false;
            for e in 3 * g..3 * g + 3 {
                self.twin[e as usize] = NONE;
            }
        }
        let kept: Vec<u32> = around
            .into_iter()
            .filter(|&g| self.live[g as usize])
            .collect();
        for &g in &kept {
            for c in &mut self.corners[g as usize] {
                if *c == u {
                    *c = w;
                }
            }
        }
        self.relink(&Scope::Points(points), &kept);
        true
    }
}

/// How a face stands against the way it faced before the offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Shape {
    /// Turned from its normal before by no more than the edge angle.
    Kept,
    /// Turned by more, short of a sliver.
    Turned,
    /// A sliver: its height is within the rounding of its corners.
    Sliver,
}

/// How the face on `a`, `b` and `c`, its area vector `was` before the
/// offset, stands.
fn shape(a: Point, b: Point, c: Point, was: Vector) -> Shape {
    let now = (b - a).cross(c - a);
    let longest = (b - a)
        .magnitude()
        .max((c - b).magnitude())
        .max((a - c).magnitude());
    let size = [a, b, c]
        .iter()
        .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
        .fold(0.0, f64::max);
    if now.magnitude() <= SLIVER * f64::EPSILON * size * longest {
        Shape::Sliver
    } else if now.dot(was) < OFFSET_EDGE_ANGLE.cos() * now.magnitude() * was.magnitude() {
        Shape::Turned
    } else {
        Shape::Kept
    }
}

/// The solution of the first `k` rows and columns of `m x = r`, by
/// elimination with partial pivoting; `None` where a pivot is zero.
fn solve(mut m: [[f64; 3]; 3], mut r: [f64; 3], k: usize) -> Option<[f64; 3]> {
    for col in 0..k {
        let pivot = (col..k).max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))?;
        if m[pivot][col] == 0.0 {
            return None;
        }
        m.swap(col, pivot);
        r.swap(col, pivot);
        for row in col + 1..k {
            let s = m[row][col] / m[col][col];
            let above = m[col];
            for (mj, &aj) in m[row].iter_mut().zip(&above).take(k).skip(col) {
                *mj -= s * aj;
            }
            r[row] -= s * r[col];
        }
    }
    let mut x = [0.0; 3];
    for row in (0..k).rev() {
        let mut s = r[row];
        for j in row + 1..k {
            s -= m[row][j] * x[j];
        }
        x[row] = s / m[row][row];
    }
    Some(x)
}
