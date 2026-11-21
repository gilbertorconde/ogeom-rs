//! Remeshing an editable mesh to an even edge length, and subdividing it.
//!
//! # Remeshing
//!
//! [`Editable::remesh`] makes [`Remesh::passes`] passes of isotropic
//! remeshing. Each pass takes four steps, each over the edges in the order
//! of their half-edge numbers:
//!
//! 1. every edge longer than 4/3 of [`Remesh::edge`] is split at its
//!    midpoint, and each half split again while it is still too long;
//! 2. every edge shorter than 4/5 of it is collapsed, one end moved onto
//!    the other, where the mesh stays a 2-manifold (the two ends share no
//!    neighbour but the corners across the edge, and none of those is left
//!    with fewer than three edges), no edge ends longer than the split
//!    length, and no face turns round;
//! 3. every edge is flipped whose flip brings the edge counts of its four
//!    vertices nearer six, four on a rim, where no face turns round;
//! 4. each vertex free to move goes to the mean of its neighbours, taken
//!    across its tangent plane, and then onto the nearest point of the mesh
//!    as it was before remeshing, among the triangles facing the vertex's
//!    way. A move that would turn a face round is not made.
//!
//! # Sharp edges
//!
//! An edge on a rim is sharp, and with [`Remesh::keep_edges`] so is an edge
//! whose two faces' normals turn by more than that angle. A sharp edge is
//! never flipped, and a collapse along it moves only a vertex on it. A
//! vertex on exactly two sharp edges slides along them where they turn by
//! no more than the crease angle at it ([`RIM_CORNER`] on a rim when there
//! is none): it collapses only along them, and moves to the mean of its two
//! neighbours on them, then onto the nearest point of the sharp edges as
//! they were. Every other vertex on a sharp edge, the end of a crease, a
//! corner where three or more meet, or a sharper turn, never moves.
//!
//! # Region
//!
//! With [`Remesh::region`], a face with a corner outside the region stays
//! as it is, corners and positions, and so does every vertex of such a
//! face. A face [`Editable::non_manifold`] reports a corner of stays alike.
//!
//! # Subdividing
//!
//! [`Editable::subdivide`] cuts every face into four at its edge midpoints,
//! as many times as asked. [`Subdivision::Flat`] puts each new vertex at
//! its edge's midpoint and leaves the others where they are.
//! [`Subdivision::Smooth`] is Loop's scheme: a new vertex is weighted 3/8
//! on its edge's ends and 1/8 on the corners across it, and a vertex there
//! before goes to `(1 - nβ)` itself and `β` each of its `n` neighbours,
//! `β = (5/8 - (3/8 + cos(2π/n)/4)²)/n`. On a sharp edge, as remeshing
//! reads it, a new vertex is the midpoint, and a vertex sliding along two
//! goes to 3/4 itself and 1/8 each of its two neighbours on them; a vertex
//! that would not slide stays. Repeated, the surface tends to Loop's limit
//! surface.
//!
//! # What stays
//!
//! As [`Editable::simplify`], the result is built afresh: vertex `v`
//! afterwards is the `v`th vertex left, and face `f` the `f`th face. A
//! vertex kept keeps its source vertices, a collapse leaves the vertex kept
//! with those of both, and a vertex the edit makes has none. A request
//! refused leaves the mesh as it was.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;

use super::{Editable, NONE, index, next, prev};
use crate::bvh::Bvh;

/// The turn, in radians, beyond which a rim vertex stays where it is when
/// [`Remesh::keep_edges`] or [`Subdivision::Smooth`] gives no crease angle.
pub const RIM_CORNER: f64 = core::f64::consts::FRAC_PI_6;

/// How a mesh is remeshed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Remesh<'a> {
    /// The edge length aimed for.
    pub edge: f64,
    /// Passes of split, collapse, flip and relax; 5 a sensible default.
    pub passes: usize,
    /// The angle, in radians, beyond which an edge is a crease: creases
    /// stay, their vertices sliding along them only. `None`: no creases.
    pub keep_edges: Option<f64>,
    /// Which vertices may change, one entry per vertex number; the faces
    /// with a corner outside stay as they are. `None` for all.
    pub region: Option<&'a [bool]>,
    /// Refused rather than exceeded.
    pub max_faces: Option<usize>,
}

/// What a remeshing did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemeshReport {
    /// How many edges were split.
    pub split: usize,
    /// How many edges were collapsed.
    pub collapsed: usize,
    /// How many edges were flipped.
    pub flipped: usize,
    /// How far the new surface lies from the old, the largest: the
    /// distance to the mesh as it was from each vertex of the result and
    /// each face's centroid.
    pub max: f64,
    /// The root mean square of the same distances.
    pub rms: f64,
}

/// How a subdivision places the points it adds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Subdivision {
    /// Each face into four at its edge midpoints: the surface unchanged.
    Flat,
    /// Loop's scheme, creases (as [`Remesh::keep_edges`]) and rims kept
    /// sharp.
    Smooth {
        /// The crease angle, as [`Remesh::keep_edges`].
        keep_edges: Option<f64>,
    },
}

/// What a subdivision did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SubdivideReport {
    /// How many faces the mesh has after.
    pub faces: usize,
    /// How many vertices were added.
    pub added: usize,
    /// The farthest a vertex there before moved: zero for
    /// [`Subdivision::Flat`].
    pub moved: f64,
}

impl Editable {
    /// Remesh toward edges of [`Remesh::edge`], as the [module
    /// documentation](crate::editable::remesh) describes.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if `edge` is
    /// not positive and finite, or `keep_edges` is not an angle from zero to
    /// `π`; [`OgeomError::Range`](ogeom_core::OgeomError::Range) if `region`
    /// does not have one entry per vertex number, or the result would have
    /// more faces than `max_faces`: judged before any work from the area
    /// to remesh over that of an equilateral triangle of side `edge`, and
    /// after each pass from the faces there are. Nothing changes when it
    /// refuses.
    pub fn remesh(&mut self, remesh: &Remesh<'_>) -> OgeomResult<RemeshReport> {
        let Remesh {
            edge,
            passes,
            keep_edges,
            region,
            max_faces,
        } = *remesh;
        if !(edge.is_finite() && edge > 0.0) {
            ogeom_bail!(Numeric, "an edge length of {edge} is not a length");
        }
        check_angle(keep_edges)?;
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

        let mut w = Work::of(self, region, keep_edges);
        #[cfg(test)]
        w.validate("build");
        let estimate = w.estimate(edge);
        // Splitting longest first keeps triangles near the estimate's
        // shape; eight times it, or the faces there are, bounds a split
        // step that goes wrong.
        #[allow(clippy::cast_possible_truncation, reason = "saturating, and capped")]
        let bound = (8.0 * estimate).min(f64::from(NONE) / 4.0) as usize;
        let mut cap = bound.max(4 * w.faces);
        if let Some(most) = max_faces {
            // A split step overshoots the pass's result, which collapses
            // then thin: it may reach twice the faces allowed.
            cap = cap.min(most.saturating_mul(2));
            if estimate > f64::from(u32::try_from(most).unwrap_or(NONE)) {
                ogeom_bail!(
                    Range,
                    "remeshing to edges of {edge} makes about {estimate:.0} faces, past the {most} allowed"
                );
            }
        }
        let input = w.surface();
        let bvh = Bvh::build(&input)?;
        let facing: Vec<Vector> = input
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|k| input.positions[k as usize]);
                (b - a).cross(c - a)
            })
            .collect();
        let rails = w.rails();
        let rails = if rails.triangles.is_empty() {
            None
        } else {
            Some(Bvh::build(&rails)?)
        };

        let (high, low) = (edge * 4.0 / 3.0, edge * 4.0 / 5.0);
        let (mut split, mut collapsed, mut flipped) = (0, 0, 0);
        for _ in 0..passes {
            let Some(splits) = w.split_long(high, cap) else {
                match max_faces {
                    Some(most) if w.faces > most.saturating_mul(2) => ogeom_bail!(
                        Range,
                        "remeshing to edges of {edge} split to {} faces, past twice the {most} allowed",
                        w.faces
                    ),
                    _ => ogeom_bail!(
                        NotDone,
                        "remeshing to edges of {edge} split past {cap} faces, eight times the estimate"
                    ),
                }
            };
            split += splits;
            #[cfg(test)]
            w.validate("splitting");
            collapsed += w.collapse_short(low, high);
            #[cfg(test)]
            w.validate("collapsing");
            flipped += w.flip_to_valence();
            #[cfg(test)]
            w.validate("flipping");
            w.relax(&bvh, &facing, rails.as_ref());
            if let Some(most) = max_faces
                && w.faces > most
            {
                ogeom_bail!(
                    Range,
                    "remeshing to edges of {edge} reached {} faces, past the {most} allowed",
                    w.faces
                );
            }
        }
        let (max, rms) = w.deviation(&bvh);
        *self = w.finish(self.source_count, self.deflection_met)?;
        Ok(RemeshReport {
            split,
            collapsed,
            flipped,
            max,
            rms,
        })
    }

    /// Subdivide every face into four `levels` times, as the [module
    /// documentation](crate::editable::remesh) describes.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a crease
    /// angle is not from zero to `π`;
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if the result
    /// would have more faces than `max_faces`, or than `u32` numbers,
    /// judged before any work. Nothing changes when it refuses.
    pub fn subdivide(
        &mut self,
        levels: u32,
        scheme: Subdivision,
        max_faces: Option<usize>,
    ) -> OgeomResult<SubdivideReport> {
        let keep_edges = match scheme {
            Subdivision::Flat => None,
            Subdivision::Smooth { keep_edges } => keep_edges,
        };
        check_angle(keep_edges)?;
        let live = self.live.iter().filter(|&&l| l).count();
        let total = 4_usize
            .checked_pow(levels)
            .and_then(|k| k.checked_mul(live))
            .filter(|&t| t.saturating_mul(3) < NONE as usize);
        let Some(total) = total else {
            ogeom_bail!(
                Range,
                "{levels} subdivisions of {live} faces are past u32 numbering"
            );
        };
        if let Some(most) = max_faces
            && total > most
        {
            ogeom_bail!(
                Range,
                "{levels} subdivisions make {total} faces, past the {most} allowed"
            );
        }
        if levels == 0 {
            return Ok(SubdivideReport {
                faces: live,
                added: 0,
                moved: 0.0,
            });
        }
        let smooth = matches!(scheme, Subdivision::Smooth { .. });
        let mut w = Work::of(self, None, keep_edges);
        let first = w.pos.clone();
        for _ in 0..levels {
            w = w.subdivided(smooth);
        }
        let moved = first
            .iter()
            .zip(&w.pos)
            .map(|(p, q)| p.distance(*q))
            .fold(0.0, f64::max);
        let added = w.pos.len() - first.len();
        let faces = w.faces;
        *self = w.finish(self.source_count, self.deflection_met)?;
        Ok(SubdivideReport {
            faces,
            added,
            moved,
        })
    }
}

fn check_angle(keep_edges: Option<f64>) -> OgeomResult<()> {
    if let Some(angle) = keep_edges
        && !(0.0..=core::f64::consts::PI).contains(&angle)
    {
        ogeom_bail!(Numeric, "a crease angle of {angle} is not from 0 to π");
    }
    Ok(())
}

/// How a vertex may move.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// Anywhere on the surface, and collapsed along any edge.
    Free,
    /// Along its two sharp edges only.
    Slide,
    /// Never.
    Fixed,
}

/// A manifold mesh of half-edges as [`Editable`] numbers them, side `k` of
/// face `f` being half-edge `3 f + k`, with every vertex one fan. Faces and
/// vertices removed keep their numbers until [`Work::finish`].
struct Work {
    pos: Vec<Point>,
    sources: Vec<Vec<u32>>,
    kind: Vec<Kind>,
    alive: Vec<bool>,
    /// A half-edge leaving each vertex, the one on its rim where there is
    /// one; `NONE` for a vertex with no face.
    out: Vec<u32>,
    corners: Vec<[u32; 3]>,
    live: Vec<bool>,
    /// Faces that stay as they are.
    frozen: Vec<bool>,
    twin: Vec<u32>,
    sharp: Vec<bool>,
    /// How many faces are live.
    faces: usize,
}

/// A vertex's neighbours, each with whether the edge to it is sharp.
type Ring = Vec<(u32, bool)>;

impl Work {
    /// The vertices left of `e`, in order, and its faces left, joined as
    /// `e` joins them; the vertices outside `region` or not a manifold
    /// fixed with every face at them.
    fn of(e: &Editable, region: Option<&[bool]>, keep_edges: Option<f64>) -> Self {
        let n = e.positions.len();
        let mut map = vec![NONE; n];
        let mut pos = Vec::new();
        let mut sources = Vec::new();
        let mut blocked = Vec::new();
        for v in 0..n {
            if !e.gone[v] {
                map[v] = index(pos.len());
                pos.push(e.positions[v]);
                sources.push(e.source_vertices(index(v)).to_vec());
                blocked.push(region.is_some_and(|r| !r[v]));
            }
        }
        for v in e.non_manifold() {
            blocked[map[v as usize] as usize] = true;
        }
        let mut face_map = vec![NONE; e.corners.len()];
        let mut corners = Vec::new();
        for (f, &live) in e.live.iter().enumerate() {
            if live {
                face_map[f] = index(corners.len());
                corners.push(e.corners[f].map(|v| map[v as usize]));
            }
        }
        let mut twin = Vec::with_capacity(corners.len() * 3);
        for (f, &live) in e.live.iter().enumerate() {
            if live {
                for k in 0..3 {
                    let g = e.twin[3 * f + k];
                    twin.push(if g == NONE {
                        NONE
                    } else {
                        3 * face_map[(g / 3) as usize] + g % 3
                    });
                }
            }
        }
        let frozen: Vec<bool> = corners
            .iter()
            .map(|t| t.iter().any(|&v| blocked[v as usize]))
            .collect();
        let m = pos.len();
        let mut w = Self {
            pos,
            sources,
            kind: vec![Kind::Fixed; m],
            alive: vec![true; m],
            out: vec![NONE; m],
            live: vec![true; corners.len()],
            faces: corners.len(),
            frozen,
            sharp: twin.iter().map(|&g| g == NONE).collect(),
            twin,
            corners,
        };
        w.find_outs();
        if let Some(angle) = keep_edges {
            let cos = angle.cos();
            for h in 0..w.twin.len() {
                let g = w.twin[h];
                if g == NONE || (g as usize) < h {
                    continue;
                }
                let (a, b) = (w.normal(index(h) / 3), w.normal(g / 3));
                let (la, lb) = (a.magnitude(), b.magnitude());
                if la > 0.0 && lb > 0.0 && a.dot(b) / (la * lb) < cos {
                    w.sharp[h] = true;
                    w.sharp[g as usize] = true;
                }
            }
        }
        for (f, &frozen) in w.frozen.iter().enumerate() {
            if frozen {
                for v in w.corners[f] {
                    blocked[v as usize] = true;
                }
            }
        }
        w.classify(&blocked, keep_edges.unwrap_or(RIM_CORNER));
        w
    }

    /// Set each vertex's `out` from the live faces, a rim half-edge where
    /// there is one.
    fn find_outs(&mut self) {
        self.out.iter_mut().for_each(|o| *o = NONE);
        for h in 0..self.twin.len() {
            if !self.live[h / 3] {
                continue;
            }
            let h = index(h);
            let v = self.tail(h) as usize;
            let o = self.out[v];
            if o == NONE || (self.twin[h as usize] == NONE && self.twin[o as usize] != NONE) {
                self.out[v] = h;
            }
        }
    }

    /// Give each vertex its kind: fixed where `blocked`, with no face, or
    /// on sharp edges other than two turning by at most `turn`.
    fn classify(&mut self, blocked: &[bool], turn: f64) {
        let cos = turn.cos();
        let mut ring = Ring::new();
        for v in 0..self.pos.len() {
            self.kind[v] = if blocked[v] || self.out[v] == NONE {
                Kind::Fixed
            } else {
                self.ring(index(v), &mut ring);
                let sharp: Vec<u32> = ring.iter().filter(|r| r.1).map(|r| r.0).collect();
                match *sharp.as_slice() {
                    [] => Kind::Free,
                    [s, t] => {
                        let p = self.pos[v];
                        let (into, from) = (p - self.pos[s as usize], self.pos[t as usize] - p);
                        let lengths = into.magnitude() * from.magnitude();
                        if lengths > 0.0 && into.dot(from) >= cos * lengths {
                            Kind::Slide
                        } else {
                            Kind::Fixed
                        }
                    }
                    _ => Kind::Fixed,
                }
            };
        }
    }

    fn tail(&self, h: u32) -> u32 {
        self.corners[(h / 3) as usize][(h % 3) as usize]
    }

    fn head(&self, h: u32) -> u32 {
        self.tail(next(h))
    }

    fn length2(&self, h: u32) -> f64 {
        let d = self.pos[self.head(h) as usize] - self.pos[self.tail(h) as usize];
        d.dot(d)
    }

    /// Face `f`'s normal, its length twice the face's area.
    fn normal(&self, f: u32) -> Vector {
        let [a, b, c] = self.corners[f as usize].map(|v| self.pos[v as usize]);
        (b - a).cross(c - a)
    }

    fn is_rim(&self, v: u32) -> bool {
        let o = self.out[v as usize];
        o != NONE && self.twin[o as usize] == NONE
    }

    /// The neighbours of `v`, turning from its `out`.
    fn ring(&self, v: u32, ring: &mut Ring) {
        ring.clear();
        let start = self.out[v as usize];
        if start == NONE {
            return;
        }
        let mut h = start;
        // A fan has fewer half-edges than the mesh: the bound only stops a
        // walk round twins that never close.
        for _ in 0..self.twin.len() {
            ring.push((self.head(h), self.sharp[h as usize]));
            let before = prev(h);
            let t = self.twin[before as usize];
            if t == NONE {
                ring.push((self.tail(before), true));
                return;
            }
            h = t;
            if h == start {
                return;
            }
        }
        debug_assert!(false, "the fan of vertex {v} does not close");
    }

    /// The half-edges leaving `v`, turning from its `out`.
    fn leaving(&self, v: u32, leaving: &mut Vec<u32>) {
        leaving.clear();
        let start = self.out[v as usize];
        if start == NONE {
            return;
        }
        let mut h = start;
        for _ in 0..self.twin.len() {
            leaving.push(h);
            h = self.twin[prev(h) as usize];
            if h == NONE || h == start {
                return;
            }
        }
        debug_assert!(false, "the fan of vertex {v} does not close");
    }

    /// Set `v`'s `out` from `h`, a live half-edge leaving it: turned back
    /// to the rim where there is one.
    fn set_out(&mut self, v: u32, h: u32) {
        let mut h = h;
        let start = h;
        for _ in 0..self.twin.len() {
            let t = self.twin[h as usize];
            if t == NONE {
                break;
            }
            let g = next(t);
            if g == start {
                break;
            }
            h = g;
        }
        self.out[v as usize] = h;
    }

    /// Make `x` and `y` twins, either `NONE`; an edge left on a rim is
    /// sharp.
    fn join(&mut self, x: u32, y: u32, sharp: bool) {
        if x != NONE {
            self.twin[x as usize] = y;
            self.sharp[x as usize] = sharp || y == NONE;
        }
        if y != NONE {
            self.twin[y as usize] = x;
            self.sharp[y as usize] = sharp || x == NONE;
        }
    }

    /// The live faces as a mesh, their corners numbered as the vertices.
    fn surface(&self) -> Triangulation {
        Triangulation {
            positions: self.pos.clone(),
            triangles: self
                .corners
                .iter()
                .zip(&self.live)
                .filter(|&(_, &l)| l)
                .map(|(&t, _)| t)
                .collect(),
            ..Triangulation::new()
        }
    }

    /// The sharp edges, each as a triangle with its second end twice.
    fn rails(&self) -> Triangulation {
        let mut triangles = Vec::new();
        for h in 0..self.twin.len() {
            let g = self.twin[h];
            if self.live[h / 3] && self.sharp[h] && (g == NONE || h < g as usize) {
                let h = index(h);
                let (a, b) = (self.tail(h), self.head(h));
                triangles.push([a, b, b]);
            }
        }
        Triangulation {
            positions: self.pos.clone(),
            triangles,
            ..Triangulation::new()
        }
    }

    /// The faces a remeshing to `edge` makes, about: the faces that stay,
    /// and the area of the rest over an equilateral triangle's of side
    /// `edge`.
    fn estimate(&self, edge: f64) -> f64 {
        let mut kept = 0.0;
        let mut area = 0.0;
        for f in 0..self.corners.len() {
            if !self.live[f] {
                continue;
            }
            if self.frozen[f] {
                kept += 1.0;
            } else {
                area += self.normal(index(f)).magnitude() / 2.0;
            }
        }
        kept + area / (3.0_f64.sqrt() / 4.0 * edge * edge)
    }

    fn splittable(&self, h: u32) -> bool {
        let g = self.twin[h as usize];
        self.live[(h / 3) as usize]
            && !self.frozen[(h / 3) as usize]
            && (g == NONE || (g > h && !self.frozen[(g / 3) as usize]))
    }

    /// Split the edges longer than `high`, in sweeps that take each such
    /// edge once, longest first; how many splits, or `None` as soon as the
    /// faces outnumber `cap`.
    ///
    /// A sweep halves every long edge, and an edge it makes is no longer
    /// than the longest it began with, so the longest edge falls below
    /// `high` within `log2(longest / high)` sweeps and a few more for the
    /// edges made across faces. That many sweeps are made at most: beside
    /// an edge that stays, on a region's border, the edges to the corner
    /// across can halve towards a length above `high` without end. What a
    /// pass leaves long, the next pass splits.
    fn split_long(&mut self, high: f64, cap: usize) -> Option<usize> {
        let high2 = high * high;
        let mut count = 0;
        let mut long: Vec<(f64, u32)> = Vec::new();
        let longest = (0..index(self.twin.len()))
            .filter(|&h| self.splittable(h))
            .map(|h| self.length2(h))
            .fold(0.0, f64::max)
            .sqrt();
        let sweeps = if longest > high {
            #[allow(clippy::cast_possible_truncation, reason = "saturating")]
            let halvings = (longest / high).log2().ceil() as usize;
            halvings.saturating_add(3)
        } else {
            0
        };
        for _ in 0..sweeps {
            long.clear();
            for h in 0..index(self.twin.len()) {
                if self.splittable(h) {
                    let l = self.length2(h);
                    if l > high2 {
                        long.push((l, h));
                    }
                }
            }
            if long.is_empty() {
                break;
            }
            long.sort_unstable_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
            for &(_, h) in &long {
                // A split moves the half-edge after the one it splits to a
                // new number, so a number listed may name another edge by
                // now: it is split only if that edge is still long.
                if !self.splittable(h) || self.length2(h) <= high2 {
                    continue;
                }
                self.split(h);
                count += 1;
                if self.faces > cap {
                    return None;
                }
            }
        }
        Some(count)
    }

    /// Split the edge of half-edge `h` at its midpoint, `h` keeping the
    /// half from its tail.
    fn split(&mut self, h: u32) {
        let g = self.twin[h as usize];
        let sharp = self.sharp[h as usize];
        let (a, b) = (self.tail(h), self.head(h));
        let m = index(self.pos.len());
        self.pos
            .push(self.pos[a as usize].midpoint(self.pos[b as usize]));
        self.sources.push(Vec::new());
        self.kind.push(if sharp { Kind::Slide } else { Kind::Free });
        self.alive.push(true);
        self.out.push(NONE);
        let to_b = self.split_face(h, m);
        if g == NONE {
            self.join(h, NONE, true);
            self.join(to_b, NONE, true);
        } else {
            let to_a = self.split_face(g, m);
            self.join(h, to_a, sharp);
            self.join(g, to_b, sharp);
        }
        self.set_out(m, to_b);
    }

    /// Cut the face of `h`, from `a` to `b` with `c` across, into `a m c`,
    /// which keeps the face's number, and a new face `m b c`, whose first
    /// half-edge, from `m` to `b`, is returned unjoined.
    fn split_face(&mut self, h: u32, m: u32) -> u32 {
        let f = (h / 3) as usize;
        let from_b = next(h);
        let (b, c) = (self.head(h), self.tail(prev(h)));
        let (bc_twin, bc_sharp) = (self.twin[from_b as usize], self.sharp[from_b as usize]);
        self.corners[f][(from_b % 3) as usize] = m;
        let nf = index(self.corners.len());
        self.corners.push([m, b, c]);
        self.live.push(true);
        self.frozen.push(false);
        self.twin.extend([NONE; 3]);
        self.sharp.extend([false; 3]);
        let (to_b, bc, cm) = (3 * nf, 3 * nf + 1, 3 * nf + 2);
        self.twin[bc as usize] = bc_twin;
        if bc_twin != NONE {
            self.twin[bc_twin as usize] = bc;
        }
        self.sharp[bc as usize] = bc_sharp;
        self.join(from_b, cm, false);
        if self.out[b as usize] == from_b {
            self.out[b as usize] = bc;
        }
        self.faces += 1;
        to_b
    }

    /// Collapse every edge shorter than `low` that may be, by
    /// [`Work::may_collapse`]; how many collapses.
    fn collapse_short(&mut self, low: f64, high: f64) -> usize {
        let low2 = low * low;
        let mut count = 0;
        let mut scratch = Scratch::default();
        for h in 0..index(self.twin.len()) {
            let g = self.twin[h as usize];
            if !self.live[(h / 3) as usize] || (g != NONE && g < h) || self.length2(h) >= low2 {
                continue;
            }
            // The freer end goes.
            let (a, b) = (self.tail(h), self.head(h));
            let order =
                if self.kind[b as usize] == Kind::Free && self.kind[a as usize] != Kind::Free {
                    [g, h]
                } else {
                    [h, g]
                };
            for k in order {
                if k != NONE && self.may_collapse(k, high, &mut scratch) {
                    self.collapse(k, &mut scratch.leaving);
                    count += 1;
                    break;
                }
            }
        }
        count
    }

    /// Whether `h`'s tail may be collapsed onto its head: the tail moves
    /// along this edge, the mesh stays a 2-manifold, no edge ends longer
    /// than `high`, and no face turns round.
    fn may_collapse(&self, h: u32, high: f64, s: &mut Scratch) -> bool {
        let (a, b) = (self.tail(h), self.head(h));
        let g = self.twin[h as usize];
        match self.kind[a as usize] {
            Kind::Free => {}
            Kind::Slide if self.sharp[h as usize] => {}
            _ => return false,
        }
        if g != NONE && self.is_rim(a) {
            return false;
        }
        let c = self.tail(prev(h));
        let d = if g == NONE { NONE } else { self.tail(prev(g)) };
        self.ring(a, &mut s.ring_a);
        self.ring(b, &mut s.ring_b);
        let mut common = 0;
        for &(x, _) in &s.ring_a {
            if s.ring_b.iter().any(|&(y, _)| y == x) {
                if x != c && x != d {
                    return false;
                }
                common += 1;
            }
        }
        if common != if g == NONE { 1 } else { 2 } {
            return false;
        }
        if s.ring_a.len() + s.ring_b.len() - common - 2 < 3 {
            return false;
        }
        for x in [c, d] {
            if x != NONE {
                self.ring(x, &mut s.ring_c);
                let least = if self.is_rim(x) { 3 } else { 4 };
                if s.ring_c.len() < least {
                    return false;
                }
            }
        }
        let pb = self.pos[b as usize];
        let high2 = high * high;
        for &(x, _) in &s.ring_a {
            if x != b {
                let d = self.pos[x as usize] - pb;
                if d.dot(d) > high2 {
                    return false;
                }
            }
        }
        self.leaving(a, &mut s.leaving);
        for &l in &s.leaving {
            let f = l / 3;
            if f == h / 3 || (g != NONE && f == g / 3) {
                continue;
            }
            let before = self.normal(f);
            let [p, q, r] =
                self.corners[f as usize].map(|v| if v == a { pb } else { self.pos[v as usize] });
            if (q - p).cross(r - p).dot(before) <= 0.0 {
                return false;
            }
        }
        true
    }

    /// Collapse `h`'s tail onto its head. `leaving` is scratch.
    fn collapse(&mut self, h: u32, leaving: &mut Vec<u32>) {
        let (a, b) = (self.tail(h), self.head(h));
        let g = self.twin[h as usize];
        self.leaving(a, leaving);
        let mut ends = Vec::with_capacity(2);
        for k in [h, g] {
            if k == NONE {
                continue;
            }
            // The face's two other edges become one.
            let (from_head, into_tail) = (next(k), prev(k));
            let across = self.tail(into_tail);
            let (x, y) = (self.twin[from_head as usize], self.twin[into_tail as usize]);
            let sharp = self.sharp[from_head as usize] || self.sharp[into_tail as usize];
            self.join(x, y, sharp);
            self.live[(k / 3) as usize] = false;
            self.faces -= 1;
            ends.push((across, x, y));
        }
        for &l in leaving.iter() {
            if self.live[(l / 3) as usize] {
                self.corners[(l / 3) as usize][(l % 3) as usize] = b;
            }
        }
        let gone = core::mem::take(&mut self.sources[a as usize]);
        let kept = &mut self.sources[b as usize];
        kept.extend(gone);
        kept.sort_unstable();
        kept.dedup();
        self.alive[a as usize] = false;
        self.out[a as usize] = NONE;
        // `x` runs from the corner across to `b`, and `y` from `b` to it,
        // so the half-edge after `y` leaves it too.
        for (across, x, y) in ends {
            match (x, y) {
                (NONE, NONE) => self.out[across as usize] = NONE,
                (NONE, y) => self.set_out(across, next(y)),
                (x, _) => self.set_out(across, x),
            }
        }
        let from_b = leaving
            .iter()
            .copied()
            .find(|&l| self.live[(l / 3) as usize])
            .unwrap_or(NONE);
        if from_b == NONE {
            self.out[b as usize] = NONE;
        } else {
            self.set_out(b, from_b);
        }
    }

    /// The valence a vertex is brought towards, or `None` for one that
    /// never moves.
    fn target(&self, v: u32) -> Option<usize> {
        match self.kind[v as usize] {
            Kind::Fixed => None,
            _ if self.is_rim(v) => Some(4),
            _ => Some(6),
        }
    }

    /// Flip every edge whose flip brings its four vertices' valences
    /// nearer their targets; how many flips.
    fn flip_to_valence(&mut self) -> usize {
        let mut count = 0;
        let mut scratch = Scratch::default();
        for h in 0..index(self.twin.len()) {
            if self.may_flip(h, &mut scratch) {
                self.flip(h);
                count += 1;
            }
        }
        count
    }

    fn may_flip(&self, h: u32, s: &mut Scratch) -> bool {
        let g = self.twin[h as usize];
        let (f, fg) = ((h / 3) as usize, (g / 3) as usize);
        if g == NONE
            || g < h
            || !self.live[f]
            || self.sharp[h as usize]
            || self.frozen[f]
            || self.frozen[fg]
        {
            return false;
        }
        let (a, b) = (self.tail(h), self.head(h));
        let (c, d) = (self.tail(prev(h)), self.tail(prev(g)));
        if c == d {
            return false;
        }
        self.ring(c, &mut s.ring_c);
        if s.ring_c.iter().any(|&(x, _)| x == d) {
            return false;
        }
        let mut before = 0;
        let mut after = 0;
        for (v, step) in [(a, -1), (b, -1), (c, 1), (d, 1)] {
            self.ring(v, &mut s.ring_a);
            let valence = s.ring_a.len() as isize;
            if step < 0 && valence - 1 < if self.is_rim(v) { 2 } else { 3 } {
                return false;
            }
            if let Some(t) = self.target(v) {
                let t = t as isize;
                before += (valence - t).abs();
                after += (valence + step - t).abs();
            }
        }
        if after >= before {
            return false;
        }
        let [pa, pb, pc, pd] = [a, b, c, d].map(|v| self.pos[v as usize]);
        let old = [self.normal(index(f)), self.normal(index(fg))];
        let new = [(pd - pa).cross(pc - pa), (pb - pd).cross(pc - pd)];
        new.iter().all(|n| old.iter().all(|o| n.dot(*o) > 0.0))
    }

    /// Flip `h`'s edge, from `a` to `b` with `c` and `d` across, to run
    /// from `d` to `c`: faces `a d c` and `d b c`.
    fn flip(&mut self, h: u32) {
        let g = self.twin[h as usize];
        let (f, fg) = (h / 3, g / 3);
        let (a, b) = (self.tail(h), self.head(h));
        let (c, d) = (self.tail(prev(h)), self.tail(prev(g)));
        let info = |me: &Self, x: u32| (me.twin[x as usize], me.sharp[x as usize]);
        let (bc, ca, ad, db) = (
            info(self, next(h)),
            info(self, prev(h)),
            info(self, next(g)),
            info(self, prev(g)),
        );
        self.corners[f as usize] = [a, d, c];
        self.corners[fg as usize] = [d, b, c];
        let set = |me: &mut Self, x: u32, (t, s): (u32, bool)| {
            me.twin[x as usize] = t;
            if t != NONE {
                me.twin[t as usize] = x;
            }
            me.sharp[x as usize] = s;
        };
        set(self, 3 * f, ad);
        set(self, 3 * f + 2, ca);
        set(self, 3 * fg, db);
        set(self, 3 * fg + 1, bc);
        self.join(3 * f + 1, 3 * fg + 2, false);
        self.set_out(a, 3 * f);
        self.set_out(b, 3 * fg + 1);
        self.set_out(c, 3 * f + 2);
        self.set_out(d, 3 * fg);
    }

    /// Move each vertex free to, as the module documentation says: to its
    /// neighbours' mean across its tangent plane and onto `bvh`, among the
    /// triangles whose normals in `facing` agree with its own, or along its
    /// sharp edges and onto `rails`. Moves that would turn a face round
    /// are taken back.
    fn relax(&mut self, bvh: &Bvh, facing: &[Vector], rails: Option<&Bvh>) {
        let start = self.pos.clone();
        let mut normals = vec![Vector::new(0.0, 0.0, 0.0); self.pos.len()];
        for f in 0..self.corners.len() {
            if self.live[f] {
                let n = self.normal(index(f));
                for v in self.corners[f] {
                    normals[v as usize] += n;
                }
            }
        }
        let mut free: Vec<(u32, Point, Vector)> = Vec::new();
        let mut sliding: Vec<(u32, Point)> = Vec::new();
        let mut ring = Ring::new();
        for v in 0..self.pos.len() {
            if !self.alive[v] || self.out[v] == NONE {
                continue;
            }
            let p = self.pos[v];
            match self.kind[v] {
                Kind::Fixed => {}
                Kind::Free => {
                    let n = normals[v];
                    let length = n.magnitude();
                    if length == 0.0 {
                        continue;
                    }
                    let n = n / length;
                    self.ring(index(v), &mut ring);
                    let q = mean(ring.iter().map(|r| self.pos[r.0 as usize]));
                    free.push((index(v), q + n * n.dot(p - q), n));
                }
                Kind::Slide => {
                    self.ring(index(v), &mut ring);
                    let on: Vec<Point> = ring
                        .iter()
                        .filter(|r| r.1)
                        .map(|r| self.pos[r.0 as usize])
                        .collect();
                    if on.len() == 2 {
                        sliding.push((index(v), on[0].midpoint(on[1])));
                    }
                }
            }
        }
        let points: Vec<Point> = free.iter().map(|f| f.1).collect();
        let hits = bvh.nearest_many_where(&points, f64::INFINITY, |i, t| {
            facing[t as usize].dot(free[i].2) > 0.0
        });
        let mut moved = vec![false; self.pos.len()];
        for (&(v, q, _), hit) in free.iter().zip(hits) {
            let to = hit
                .or_else(|| bvh.nearest(q, f64::INFINITY))
                .map_or(q, |h| h.point);
            self.pos[v as usize] = to;
            moved[v as usize] = true;
        }
        if let Some(rails) = rails {
            let points: Vec<Point> = sliding.iter().map(|s| s.1).collect();
            for (&(v, q), hit) in sliding
                .iter()
                .zip(rails.nearest_many(&points, f64::INFINITY))
            {
                self.pos[v as usize] = hit.map_or(q, |h| h.point);
                moved[v as usize] = true;
            }
        }
        loop {
            let mut undone = false;
            for f in 0..self.corners.len() {
                let t = self.corners[f];
                if !self.live[f] || !t.iter().any(|&v| moved[v as usize]) {
                    continue;
                }
                let [a, b, c] = t.map(|v| start[v as usize]);
                let before = (b - a).cross(c - a);
                let after = self.normal(index(f));
                let turned = after.dot(before) < 0.0
                    || (after.dot(after) == 0.0 && before.dot(before) > 0.0);
                if turned {
                    for v in t {
                        if moved[v as usize] {
                            self.pos[v as usize] = start[v as usize];
                            moved[v as usize] = false;
                        }
                    }
                    undone = true;
                }
            }
            if !undone {
                break;
            }
        }
    }

    /// The largest and root mean square distance to `bvh` from each vertex
    /// with a face and each live face's centroid.
    fn deviation(&self, bvh: &Bvh) -> (f64, f64) {
        let mut points: Vec<Point> = (0..self.pos.len())
            .filter(|&v| self.alive[v] && self.out[v] != NONE)
            .map(|v| self.pos[v])
            .collect();
        for (t, &live) in self.corners.iter().zip(&self.live) {
            if live {
                points.push(mean(t.iter().map(|&v| self.pos[v as usize])));
            }
        }
        let (mut max, mut sum) = (0.0_f64, 0.0);
        for hit in bvh
            .nearest_many(&points, f64::INFINITY)
            .into_iter()
            .flatten()
        {
            max = max.max(hit.distance);
            sum += hit.distance * hit.distance;
        }
        let rms = if points.is_empty() {
            0.0
        } else {
            (sum / f64::from(index(points.len()))).sqrt()
        };
        (max, rms)
    }

    /// Every face cut into four at its edge midpoints, placed as the module
    /// documentation says for `smooth` or not.
    fn subdivided(&self, smooth: bool) -> Self {
        let m = self.pos.len();
        let mut pos = self.pos.clone();
        let mut kind = self.kind.clone();
        let mut ring = Ring::new();
        if smooth {
            for v in 0..m {
                let p = self.pos[v];
                pos[v] = match self.kind[v] {
                    Kind::Fixed => p,
                    Kind::Free => {
                        self.ring(index(v), &mut ring);
                        let n = f64::from(index(ring.len()));
                        let w = (3.0 / 8.0 + (core::f64::consts::TAU / n).cos() / 4.0).powi(2);
                        let beta = (5.0 / 8.0 - w) / n;
                        let sum = ring.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, r| {
                            s + self.pos[r.0 as usize].to_vector()
                        });
                        Point::from_vector(p.to_vector() * (1.0 - n * beta) + sum * beta)
                    }
                    Kind::Slide => {
                        self.ring(index(v), &mut ring);
                        let sum = ring
                            .iter()
                            .filter(|r| r.1)
                            .fold(Vector::new(0.0, 0.0, 0.0), |s, r| {
                                s + self.pos[r.0 as usize].to_vector()
                            });
                        Point::from_vector(p.to_vector() * 0.75 + sum * 0.125)
                    }
                };
            }
        }
        let mut mid = vec![NONE; self.twin.len()];
        for h in 0..self.twin.len() {
            if !self.live[h / 3] {
                continue;
            }
            let g = self.twin[h];
            if g != NONE && (g as usize) < h {
                mid[h] = mid[g as usize];
                continue;
            }
            let h = index(h);
            mid[h as usize] = index(pos.len());
            let (a, b) = (
                self.pos[self.tail(h) as usize],
                self.pos[self.head(h) as usize],
            );
            let sharp = self.sharp[h as usize];
            pos.push(if smooth && !sharp {
                let c = self.pos[self.tail(prev(h)) as usize];
                let d = self.pos[self.tail(prev(g)) as usize];
                Point::from_vector(
                    (a.to_vector() + b.to_vector()) * 0.375
                        + (c.to_vector() + d.to_vector()) * 0.125,
                )
            } else {
                a.midpoint(b)
            });
            kind.push(if sharp { Kind::Slide } else { Kind::Free });
        }
        let mut corners = Vec::with_capacity(4 * self.faces);
        let mut sharp = Vec::with_capacity(12 * self.faces);
        for f in 0..self.corners.len() {
            if !self.live[f] {
                continue;
            }
            let [a, b, c] = self.corners[f];
            let h = 3 * f;
            let [ab, bc, ca] = [mid[h], mid[h + 1], mid[h + 2]];
            let [sab, sbc, sca] = [self.sharp[h], self.sharp[h + 1], self.sharp[h + 2]];
            corners.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
            sharp.extend([
                sab, false, sca, sab, sbc, false, false, sbc, sca, false, false, false,
            ]);
        }
        let n = pos.len();
        let mut sources = self.sources.clone();
        sources.resize(n, Vec::new());
        let faces = corners.len();
        let mut w = Self {
            pos,
            sources,
            kind,
            alive: vec![true; n],
            out: vec![NONE; n],
            live: vec![true; faces],
            frozen: vec![false; faces],
            twin: vec![NONE; 3 * faces],
            sharp,
            corners,
            faces,
        };
        w.link();
        w
    }

    /// Join the half-edges running one edge opposite ways, the only two on
    /// it, and set the `out`s; a rim edge is sharp.
    fn link(&mut self) {
        let mut edges: Vec<(u64, u32)> = Vec::with_capacity(self.twin.len());
        for h in 0..index(self.twin.len()) {
            let (a, b) = (self.tail(h), self.head(h));
            edges.push(((u64::from(a.min(b)) << 32) | u64::from(a.max(b)), h));
        }
        edges.sort_unstable();
        for run in edges.chunk_by(|x, y| x.0 == y.0) {
            if let [(_, x), (_, y)] = *run
                && self.tail(x) != self.tail(y)
            {
                let sharp = self.sharp[x as usize] || self.sharp[y as usize];
                self.join(x, y, sharp);
            } else {
                for &(_, x) in run {
                    self.join(x, NONE, true);
                }
            }
        }
        self.find_outs();
    }

    /// The mesh as an [`Editable`], its vertices and faces left in order,
    /// each vertex with its source vertices.
    fn finish(self, source_count: usize, deflection_met: bool) -> OgeomResult<Editable> {
        let mut map = vec![NONE; self.pos.len()];
        let mut positions = Vec::new();
        let mut sources = Vec::new();
        for (v, set) in self.sources.into_iter().enumerate() {
            if self.alive[v] {
                map[v] = index(positions.len());
                positions.push(self.pos[v]);
                sources.push(set);
            }
        }
        let triangles = self
            .corners
            .iter()
            .zip(&self.live)
            .filter(|&(_, &l)| l)
            .map(|(t, _)| t.map(|v| map[v as usize]))
            .collect();
        let written = Triangulation {
            positions,
            triangles,
            deflection_met,
            ..Triangulation::new()
        };
        let mut next = Editable::from(&written, 0.0)?;
        // Built without a weld, each vertex's source is the written vertex
        // it is or was split from.
        for v in 0..next.positions.len() {
            let r = next.source[v] as usize;
            next.set_sources(index(v), sources[r].clone());
        }
        next.source_count = source_count;
        Ok(next)
    }
}

#[cfg(test)]
impl Work {
    /// The structure's invariants: twins run one edge opposite ways
    /// between live faces, and each vertex's `out` leaves it from its rim
    /// and sweeps every live face at it.
    fn validate(&self, after: &str) {
        let mut at = vec![0_usize; self.pos.len()];
        for f in 0..self.corners.len() {
            if !self.live[f] {
                continue;
            }
            for h in index(3 * f)..index(3 * f + 3) {
                let v = self.tail(h);
                assert!(self.alive[v as usize], "{after}: face {f} names dead {v}");
                at[v as usize] += 1;
                let g = self.twin[h as usize];
                if g != NONE {
                    assert!(self.live[(g / 3) as usize], "{after}: {h} twinned to dead");
                    assert_eq!(self.twin[g as usize], h, "{after}");
                    assert_eq!(
                        (self.tail(g), self.head(g)),
                        (self.head(h), self.tail(h)),
                        "{after}: {h}"
                    );
                }
            }
        }
        let mut leaving = Vec::new();
        for v in 0..self.pos.len() {
            if at[v] == 0 {
                continue;
            }
            let o = self.out[v];
            assert!(
                o != NONE && self.live[(o / 3) as usize],
                "{after}: out of {v}"
            );
            assert_eq!(self.tail(o), index(v), "{after}: out of {v} leaves another");
            self.leaving(index(v), &mut leaving);
            assert_eq!(leaving.len(), at[v], "{after}: fan of {v}");
        }
    }
}

/// Buffers reused from one edge to the next.
#[derive(Default)]
struct Scratch {
    ring_a: Ring,
    ring_b: Ring,
    ring_c: Ring,
    leaving: Vec<u32>,
}

fn mean(points: impl Iterator<Item = Point>) -> Point {
    let (sum, n) = points.fold((Vector::new(0.0, 0.0, 0.0), 0.0), |(s, n), p| {
        (s + p.to_vector(), n + 1.0)
    });
    Point::from_vector(if n > 0.0 { sum / n } else { sum })
}

#[cfg(test)]
mod tests;
