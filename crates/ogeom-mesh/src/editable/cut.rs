//! Cutting an editable mesh by a plane, and capping the cut flat.
//!
//! # The cut
//!
//! Which side of the plane each vertex lies on is [`section`]'s: decided
//! once per vertex by the exact orientation predicate, a vertex in the plane
//! counting as above. A vertex within rounding of the plane, so near that
//! it and its foot on the plane could take the same coordinates in it, or
//! one onto which the crossing of an edge from it rounds, is snapped onto
//! it: it counts as in the plane, so that no new vertex crowds it and no
//! sliver is left between it and the plane.
//!
//! A face whose corners all lie on the kept side stays, one whose corners
//! all lie on the other goes, and one whose corners disagree is clipped: its
//! kept part, a triangle or a quadrilateral, takes its place as new faces.
//! An edge whose ends disagree is cut at one vertex, shared by every face on
//! the edge: its upper end where that lies in the plane, as [`section`]
//! cuts it, and otherwise a new vertex standing where [`section`] puts that
//! point of its chain, to the bit. A face lying in the plane goes with the
//! side its material is on, as it does for [`section`]: kept above where
//! it faces down, its material above it, and not otherwise.
//!
//! # The caps
//!
//! The rims the cut leaves lie in the plane where both ends of every edge
//! do: a cut vertex, a vertex in the plane, or one snapped onto it. A rim
//! lying wholly in the plane is a loop. A rim that leaves the plane and
//! comes back leaves runs in it, and the end of one run is joined to the
//! start of another where they are one point, or by a straight segment, a
//! bridge, where they lie within the cap's `bridge` of each other, nearest
//! pairs first; runs joined round into a loop make one.
//!
//! The loops are capped together, so loops within loops give discs, rings
//! and discs within rings: each loop, run with the kept side's cap on its
//! left seen from the side the normal points to, winds once about the
//! points of the cap and not at all about the rest. The cap is the faces
//! of the constrained Delaunay triangulation of the loops' points, in the
//! plane's coordinates, that the loops wind about once, each wound outward
//! from the kept side. Its corners are the rim vertices themselves, and
//! every loop edge is an edge of it, so the cap shares each rim edge with
//! the face beyond it.
//!
//! A loop is left open where that cannot hold: where two of its points
//! land on one point of the plane, where an edge of it crosses another or
//! runs through a point, or where the cap's side of an edge is not wound
//! about once and the far side not at all, as where two bodies overlap.
//! The other loops are then capped again without it.
//!
//! [`section`]: crate::section()

use super::{Editable, NONE, Scope, index};
use crate::section::{Cutter, Place};
use ogeom_core::{FastMap, OgeomResult, ogeom_bail};
use ogeom_math::{Plane, Point};
use spade::handles::{FixedDirectedEdgeHandle, FixedVertexHandle};
use spade::{ConstrainedDelaunayTriangulation, Triangulation as _};

/// The side of a plane a cut keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Side {
    /// Against the plane's normal.
    Below,
    /// The way the plane's normal points, the plane itself included.
    Above,
}

/// How a cut is capped.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Cap {
    /// The longest straight segment drawn across a gap in the cut's rim to
    /// close a loop. Zero caps only the loops the rim closes itself.
    pub bridge: f64,
}

/// What a cut capped and what it left open.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CutReport {
    /// The loops capped, each closed by the rim alone or with bridges.
    pub capped_loops: usize,
    /// The bridges drawn in the loops capped.
    pub bridged: usize,
    /// The stretches of the rim in the plane left open: each loop the rim
    /// closes that was not capped, and each run of it in the plane not
    /// closed into a loop that was.
    pub open: usize,
}

/// The vertices a cut has made: by the edge cut, its points lower end
/// first, and by position.
#[derive(Default)]
struct Cuts {
    edges: FastMap<(u32, u32), u32>,
    at: FastMap<[u64; 3], u32>,
}

/// A loop to cap: vertices in order, with the cap on the left of each edge
/// to the next seen from the side the normal points to.
struct Outline {
    vertices: Vec<u32>,
    /// The bridges among its edges.
    bridges: usize,
    /// The stretches of rim it joins: one for a loop the rim closes.
    stretches: usize,
}

/// The rim in the plane: loops the rim closes, and runs of it with ends.
/// Each runs with the kept side's cap on its left seen from the side the
/// normal points to.
struct Stretches {
    loops: Vec<Vec<u32>>,
    runs: Vec<Vec<u32>>,
}

impl Editable {
    /// Cut the mesh by `plane`, keeping the faces on side `keep` and the
    /// kept parts of the faces it crosses, and cap the cut flat when `cap`
    /// asks. See the [module documentation](self) for how the cut is made
    /// and capped.
    ///
    /// Kept vertices keep their numbers and their [`source_vertex`]; the
    /// vertices the cut makes, on the edges it crosses, have none. Faces
    /// removed keep their numbers, unused; the kept parts of faces clipped,
    /// then the cap's faces, take numbers after every face there is. A mesh
    /// closed before the cut and with every loop capped is closed after.
    ///
    /// [`source_vertex`]: Editable::source_vertex
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the
    /// plane's origin is not finite or `bridge` is negative or not finite,
    /// and [`OgeomError::Range`](ogeom_core::OgeomError::Range) if the faces
    /// and vertices the cut may make would pass `u32` numbering.
    pub fn cut(&mut self, plane: &Plane, keep: Side, cap: Option<Cap>) -> OgeomResult<CutReport> {
        if !plane.origin().is_finite() {
            ogeom_bail!(Numeric, "the cutting plane's origin is not finite");
        }
        if let Some(Cap { bridge }) = cap
            && !(bridge.is_finite() && bridge >= 0.0)
        {
            ogeom_bail!(Numeric, "a bridge of {bridge} is not a distance");
        }
        // The cut makes at most a vertex per edge and two faces per face;
        // the cap at most two faces per vertex. So `F` faces and `V`
        // vertices become at most `9F + 2V` faces of `27F + 6V` half-edges.
        let (faces, vertices) = (self.corners.len(), self.positions.len());
        if faces
            .saturating_mul(27)
            .saturating_add(vertices.saturating_mul(6))
            >= NONE as usize
        {
            ogeom_bail!(
                Range,
                "a cut of {faces} faces and {vertices} vertices may pass u32 numbering"
            );
        }

        let cutter = Cutter::new(plane, 0.0);
        let level = self.clip(&cutter, keep);
        let Stretches { loops, runs } = self.stretches(&level, keep);
        let Some(Cap { bridge }) = cap else {
            return Ok(CutReport {
                open: loops.len() + runs.len(),
                ..CutReport::default()
            });
        };
        let (mut outlines, unjoined) = self.join_runs(runs, bridge);
        outlines.extend(loops.into_iter().map(|vertices| Outline {
            vertices,
            bridges: 0,
            stretches: 1,
        }));
        let capped = self.cap(&cutter, keep, &outlines);
        let mut report = CutReport {
            open: unjoined,
            ..CutReport::default()
        };
        for (o, &c) in outlines.iter().zip(&capped) {
            if c {
                report.capped_loops += 1;
                report.bridged += o.bridges;
            } else {
                report.open += o.stretches;
            }
        }
        Ok(report)
    }

    /// Clip every face to the kept side. Returns, by point, whether the
    /// point lies on the plane: in it, snapped onto it, or made by the cut.
    fn clip(&mut self, cutter: &Cutter, keep: Side) -> Vec<bool> {
        let places = self.places(cutter);
        let mut level: Vec<bool> = places.iter().map(|&p| p == Place::On).collect();
        let above = keep == Side::Above;
        let first_new = self.positions.len();
        let mut cuts = Cuts::default();
        let mut removed = Vec::new();
        let mut made: Vec<[u32; 3]> = Vec::new();
        let mut polygon: Vec<u32> = Vec::with_capacity(4);
        for f in 0..self.corners.len() {
            if !self.live[f] {
                continue;
            }
            let c = self.corners[f];
            let kept = c.map(|v| places[v as usize].up() == above);
            if kept == [true; 3] && !(above && self.faces_up(cutter, &places, c)) {
                continue;
            }
            if kept == [false; 3] || kept == [true; 3] {
                removed.push(index(f));
                continue;
            }
            polygon.clear();
            for i in 0..3 {
                let j = (i + 1) % 3;
                if kept[i] {
                    push_new(&mut polygon, c[i]);
                }
                if kept[i] != kept[j] {
                    let x = self.crossing(cutter, &places, &mut cuts, c[i], c[j]);
                    level.resize(self.positions.len(), true);
                    level[self.point[x as usize] as usize] = true;
                    push_new(&mut polygon, x);
                }
            }
            if polygon.len() > 1 && polygon.first() == polygon.last() {
                polygon.pop();
            }
            if polygon.len() == 3 && (0..3).any(|r| polygon.iter().eq(c[r..].iter().chain(&c[..r])))
            {
                continue;
            }
            removed.push(index(f));
            match polygon[..] {
                [a, b, c] => made.push([a, b, c]),
                [a, b, c, d] => {
                    // The shorter diagonal of the convex quadrilateral.
                    let at = |v: u32| self.positions[v as usize];
                    if at(a).distance(at(c)) <= at(b).distance(at(d)) {
                        made.extend([[a, b, c], [a, c, d]]);
                    } else {
                        made.extend([[b, c, d], [b, d, a]]);
                    }
                }
                _ => {}
            }
        }
        if removed.is_empty() {
            return level;
        }

        let mut points: Vec<u32> = removed
            .iter()
            .flat_map(|&f| self.corners[f as usize])
            .map(|v| self.point[v as usize])
            .collect();
        points.sort_unstable();
        points.dedup();
        let mut faces = self.about(&points);
        points.extend((first_new..self.positions.len()).map(index));
        for &f in &removed {
            self.live[f as usize] = false;
            for h in 3 * f..3 * f + 3 {
                self.twin[h as usize] = NONE;
            }
        }
        faces.retain(|f| removed.binary_search(f).is_err());
        for t in made {
            faces.push(self.add_face(t));
        }
        self.relink(&Scope::Points(points), &faces);
        // A cut vertex every face about it left degenerate is not used.
        for v in first_new..self.positions.len() {
            if self.out[v] == NONE {
                self.gone[v] = true;
            }
        }
        level
    }

    /// Each vertex's place: [`section`]'s exact side, except that a vertex
    /// the cut snaps onto the plane is in it. That is one within rounding
    /// of the plane ([`Cutter::near`]), or one onto which the crossing of
    /// an edge from it rounds. Copies of a point share its place.
    ///
    /// [`section`]: crate::section()
    fn places(&self, cutter: &Cutter) -> Vec<Place> {
        let mut places: Vec<Place> = self
            .positions
            .iter()
            .map(|&p| {
                let place = cutter.side(p);
                if place != Place::On && cutter.near(p) {
                    Place::On
                } else {
                    place
                }
            })
            .collect();
        let mut snapped = Vec::new();
        for (c, _) in self.corners.iter().zip(&self.live).filter(|&(_, &l)| l) {
            for i in 0..3 {
                let (a, b) = (c[i], c[(i + 1) % 3]);
                let (low, high) = match (places[a as usize].up(), places[b as usize].up()) {
                    (false, true) => (a, b),
                    (true, false) => (b, a),
                    _ => continue,
                };
                if places[high as usize] == Place::On {
                    continue;
                }
                let (p, q) = (self.positions[low as usize], self.positions[high as usize]);
                let x = cutter.crossing(p, q);
                if x == q {
                    snapped.push(self.point[high as usize]);
                } else if x == p {
                    snapped.push(self.point[low as usize]);
                }
            }
        }
        for p in snapped {
            for v in self.ring(p) {
                places[v as usize] = Place::On;
            }
        }
        places
    }

    /// Whether the face with corners `c`, every one in the plane, faces
    /// the way the normal points, its material below it.
    fn faces_up(&self, cutter: &Cutter, places: &[Place], c: [u32; 3]) -> bool {
        let [a, b, d] = c.map(|v| self.positions[v as usize]);
        c.iter().all(|&v| places[v as usize] == Place::On)
            && (b - a).cross(d - a).dot(cutter.normal()) > 0.0
    }

    /// The vertex the edge from `a` to `b`, its ends on opposite sides, is
    /// cut at: its upper end where that lies in the plane, and otherwise a
    /// new vertex at [`Cutter::crossing`], one for the edge and one for
    /// all the crossings at one position.
    fn crossing(
        &mut self,
        cutter: &Cutter,
        places: &[Place],
        cuts: &mut Cuts,
        a: u32,
        b: u32,
    ) -> u32 {
        let (low, high) = if places[a as usize].up() {
            (b, a)
        } else {
            (a, b)
        };
        if places[high as usize] == Place::On {
            return high;
        }
        let key = (self.point[low as usize], self.point[high as usize]);
        if let Some(&v) = cuts.edges.get(&key) {
            return v;
        }
        let x = cutter.crossing(self.positions[low as usize], self.positions[high as usize]);
        let v = *cuts
            .at
            .entry([x.x, x.y, x.z].map(f64::to_bits))
            .or_insert_with(|| self.add_vertex(x));
        cuts.edges.insert(key, v);
        v
    }

    /// A new vertex at `at`, a point of its own with no source.
    fn add_vertex(&mut self, at: Point) -> u32 {
        let v = index(self.positions.len());
        self.positions.push(at);
        self.point.push(v);
        self.copies.push(v);
        self.source.push(NONE);
        self.out.push(NONE);
        self.gone.push(false);
        v
    }

    /// A new face, joined to nothing until relinked.
    fn add_face(&mut self, corners: [u32; 3]) -> u32 {
        let f = index(self.corners.len());
        self.corners.push(corners);
        self.live.push(true);
        self.twin.extend([NONE; 3]);
        f
    }

    /// The rims' stretches in the plane, `level` marking its points.
    fn stretches(&self, level: &[bool], keep: Side) -> Stretches {
        let flat = |h: u32| {
            level[self.point[self.tail(h) as usize] as usize]
                && level[self.point[self.head(h) as usize] as usize]
        };
        let mut loops = Vec::new();
        let mut runs = Vec::new();
        for rim in self.rims() {
            let edges = &rim.edges;
            let Some(j) = edges.iter().position(|&h| !flat(h)) else {
                loops.push(edges.iter().map(|&h| self.tail(h)).collect());
                continue;
            };
            let m = edges.len();
            let mut run: Vec<u32> = Vec::new();
            for k in 1..=m {
                let h = edges[(j + k) % m];
                if flat(h) {
                    if run.is_empty() {
                        run.push(self.tail(h));
                    }
                    run.push(self.head(h));
                } else if !run.is_empty() {
                    runs.push(core::mem::take(&mut run));
                }
            }
        }
        // Rims run with their faces on the left. Seen from the side the
        // normal points to, the cap of a part above lies on the left of its
        // rim, and the cap of a part below on the right.
        if keep == Side::Below {
            for s in loops.iter_mut().chain(runs.iter_mut()) {
                s.reverse();
            }
        }
        Stretches { loops, runs }
    }

    /// Join the runs' ends into loops: the end of one to the start of
    /// another where they are one point or lie within `bridge`, nearest
    /// pairs first and ties to the lower-numbered runs. Returns the loops
    /// and how many runs no loop takes.
    fn join_runs(&self, runs: Vec<Vec<u32>>, bridge: f64) -> (Vec<Outline>, usize) {
        let at = |v: u32| self.positions[v as usize];
        let ends: Vec<u32> = runs.iter().map(|r| r[r.len() - 1]).collect();
        let starts: Vec<u32> = runs.iter().map(|r| r[0]).collect();
        let mut by_x: Vec<(f64, usize)> = starts.iter().map(|&v| (at(v).x, 0)).collect();
        for (j, s) in by_x.iter_mut().enumerate() {
            s.1 = j;
        }
        by_x.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        // Each candidate: its length, whether it is a bridge, and the runs
        // it joins.
        let mut candidates: Vec<(f64, bool, usize, usize)> = Vec::new();
        for (i, &e) in ends.iter().enumerate() {
            let p = at(e);
            let from = by_x.partition_point(|s| s.0 < p.x - bridge);
            for &(x, j) in &by_x[from..] {
                if x > p.x + bridge {
                    break;
                }
                let s = starts[j];
                if self.point[s as usize] == self.point[e as usize] {
                    candidates.push((0.0, false, i, j));
                } else {
                    let d = p.distance(at(s));
                    if d <= bridge {
                        candidates.push((d, true, i, j));
                    }
                }
            }
        }
        candidates.sort_by(|a, b| {
            a.0.total_cmp(&b.0)
                .then(a.1.cmp(&b.1))
                .then(a.2.cmp(&b.2))
                .then(a.3.cmp(&b.3))
        });
        let mut link: Vec<Option<(usize, bool)>> = vec![None; runs.len()];
        let mut taken = vec![false; runs.len()];
        for (_, bridged, i, j) in candidates {
            if link[i].is_none() && !taken[j] {
                link[i] = Some((j, bridged));
                taken[j] = true;
            }
        }

        // Each run has at most one run after it and one before, so the
        // runs fall into paths and cycles; a cycle is a loop.
        let mut seen = vec![false; runs.len()];
        let mut outlines = Vec::new();
        let mut joined = 0;
        for first in 0..runs.len() {
            if seen[first] {
                continue;
            }
            let mut cycle = vec![first];
            seen[first] = true;
            let mut at = first;
            let closed = loop {
                match link[at] {
                    Some((j, _)) if j == first => break true,
                    Some((j, _)) if !seen[j] => {
                        seen[j] = true;
                        cycle.push(j);
                        at = j;
                    }
                    _ => break false,
                }
            };
            if !closed {
                continue;
            }
            let mut vertices = Vec::new();
            let mut bridges = 0;
            for &r in &cycle {
                vertices.extend_from_slice(&runs[r]);
                if link[r].is_some_and(|(_, bridged)| bridged) {
                    bridges += 1;
                } else {
                    // The next run starts at this one's end point.
                    vertices.pop();
                }
            }
            joined += cycle.len();
            outlines.push(Outline {
                vertices,
                bridges,
                stretches: cycle.len(),
            });
        }
        (outlines, runs.len() - joined)
    }

    /// Cap the outlines together, leaving open each one that cannot be
    /// capped soundly, and return which were capped.
    fn cap(&mut self, cutter: &Cutter, keep: Side, outlines: &[Outline]) -> Vec<bool> {
        let mut active: Vec<bool> = outlines
            .iter()
            .map(|o| {
                let mut points: Vec<u32> =
                    o.vertices.iter().map(|&v| self.point[v as usize]).collect();
                points.sort_unstable();
                points.dedup();
                points.len() >= 3
            })
            .collect();
        let faces = loop {
            match self.cap_faces(cutter, outlines, &active) {
                Ok(faces) => break faces,
                Err(bad) => {
                    for k in bad {
                        active[k] = false;
                    }
                }
            }
        };
        if faces.is_empty() {
            return active;
        }
        let mut points: Vec<u32> = faces
            .iter()
            .flatten()
            .map(|&v| self.point[v as usize])
            .collect();
        points.sort_unstable();
        points.dedup();
        let mut around = self.about(&points);
        for t in faces {
            let [a, b, c] = t;
            around.push(self.add_face(if keep == Side::Below {
                [a, b, c]
            } else {
                [a, c, b]
            }));
        }
        self.relink(&Scope::Points(points), &around);
        active
    }

    /// The cap of the active outlines, each face counter-clockwise seen
    /// from the side the normal points to, or the outlines to leave open.
    fn cap_faces(
        &self,
        cutter: &Cutter,
        outlines: &[Outline],
        active: &[bool],
    ) -> Result<Vec<[u32; 3]>, Vec<usize>> {
        // Each point once, with a vertex standing for it.
        let mut points: Vec<(u32, u32)> = outlines
            .iter()
            .zip(active)
            .filter(|&(_, &a)| a)
            .flat_map(|(o, _)| o.vertices.iter().map(|&v| (self.point[v as usize], v)))
            .collect();
        points.sort_unstable();
        points.dedup_by_key(|p| p.0);
        if points.is_empty() {
            return Ok(Vec::new());
        }
        let slot = |v: u32| {
            points
                .binary_search_by_key(&self.point[v as usize], |p| p.0)
                .unwrap_or(0)
        };
        let coordinates: Vec<[f64; 2]> = points
            .iter()
            .map(|&(p, _)| cutter.along(self.positions[p as usize]))
            .collect();

        // Two points landing on one point of the plane, or on coordinates
        // the triangulation cannot take, cannot be told apart.
        let mut bad_point = vec![false; points.len()];
        for (k, c) in coordinates.iter().enumerate() {
            if c.iter().any(|&x| spade::validate_coordinate(x).is_err()) {
                bad_point[k] = true;
            }
        }
        let mut order: Vec<usize> = (0..points.len()).collect();
        order.sort_by(|&a, &b| {
            let (p, q) = (coordinates[a], coordinates[b]);
            p[0].total_cmp(&q[0]).then(p[1].total_cmp(&q[1]))
        });
        for pair in order.windows(2) {
            if coordinates[pair[0]] == coordinates[pair[1]] {
                bad_point[pair[0]] = true;
                bad_point[pair[1]] = true;
            }
        }
        let active_loops = || (0..outlines.len()).filter(|&k| active[k]);
        let bad: Vec<usize> = active_loops()
            .filter(|&k| outlines[k].vertices.iter().any(|&v| bad_point[slot(v)]))
            .collect();
        if !bad.is_empty() {
            return Err(bad);
        }

        let vertices: Vec<spade::Point2<f64>> = coordinates
            .iter()
            .map(|c| spade::Point2::new(c[0], c[1]))
            .collect();
        let Ok(mut cdt) = ConstrainedDelaunayTriangulation::<spade::Point2<f64>>::bulk_load_cdt(
            vertices,
            Vec::new(),
        ) else {
            return Err(active_loops().collect());
        };

        // Each loop edge as the directed constraint edge it became.
        let mut edges: Vec<(usize, FixedDirectedEdgeHandle)> = Vec::new();
        let mut bad = Vec::new();
        for k in active_loops() {
            let vertices = &outlines[k].vertices;
            let n = vertices.len();
            for i in 0..n {
                let (a, b) = (slot(vertices[i]), slot(vertices[(i + 1) % n]));
                if a == b {
                    continue;
                }
                let (from, to) = (
                    FixedVertexHandle::from_index(a),
                    FixedVertexHandle::from_index(b),
                );
                let made = cdt.try_add_constraint(from, to);
                match made[..] {
                    [e] => {
                        let forward = cdt.directed_edge(e).from().fix() == from;
                        edges.push((k, if forward { e } else { e.rev() }));
                    }
                    // Crossing another loop's edge: nothing was added.
                    [] => {
                        bad.push(k);
                        break;
                    }
                    // Through a point, split in the triangulation and no
                    // longer the rim's edge.
                    _ => return Err(vec![k]),
                }
            }
        }
        if !bad.is_empty() {
            return Err(bad);
        }

        let mut weight = vec![0_i64; cdt.num_directed_edges()];
        for &(_, e) in &edges {
            weight[e.index()] += 1;
            weight[e.rev().index()] -= 1;
        }
        // Counts by face, walked from the unbounded face. The face on an
        // edge's left counts its weight more than the face on its right.
        let mut count: Vec<Option<i64>> = vec![None; cdt.num_all_faces()];
        let outer = cdt.outer_face().fix().index();
        count[outer] = Some(0);
        let mut sides = vec![[None; 3]; cdt.num_all_faces()];
        for face in cdt.inner_faces() {
            sides[face.fix().index()] = face.adjacent_edges().map(|e| Some(e.fix()));
        }
        let mut stack = vec![outer];
        let mut around: Vec<FixedDirectedEdgeHandle> = Vec::with_capacity(3);
        while let Some(f) = stack.pop() {
            let Some(here) = count[f] else { continue };
            around.clear();
            if f == outer {
                for e in cdt.convex_hull() {
                    around.push(if e.face().is_outer() {
                        e.fix()
                    } else {
                        e.rev().fix()
                    });
                }
            } else {
                around.extend(sides[f].iter().flatten());
            }
            for &e in &around {
                let g = cdt.directed_edge(e).rev().face().fix().index();
                if count[g].is_none() {
                    count[g] = Some(here - weight[e.index()]);
                    stack.push(g);
                }
            }
        }

        let bad: Vec<usize> = {
            let mut bad: Vec<usize> = edges
                .iter()
                .filter(|&&(_, e)| {
                    let h = cdt.directed_edge(e);
                    count[h.face().fix().index()] != Some(1)
                        || count[h.rev().face().fix().index()] != Some(0)
                })
                .map(|&(k, _)| k)
                .collect();
            bad.dedup();
            bad
        };
        if !bad.is_empty() {
            return Err(bad);
        }
        Ok(cdt
            .inner_faces()
            .filter(|face| count[face.fix().index()] == Some(1))
            .map(|face| face.vertices().map(|v| points[v.fix().index()].1))
            .collect())
    }
}

/// Append `v` unless it repeats the last.
fn push_new(polygon: &mut Vec<u32>, v: u32) {
    if polygon.last() != Some(&v) {
        polygon.push(v);
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]
mod tests {
    use super::super::tests::{check, mesh, torus};
    use super::*;
    use crate::section::section;
    use core::f64::consts::{PI, TAU};
    use ogeom_core::Tolerances;
    use ogeom_math::{Direction, Vector};
    use ogeom_topo::Triangulation;
    use proptest::prelude::*;

    fn plane(origin: Point, normal: Vector) -> Plane {
        Plane::through(
            origin,
            Direction::new(normal, Tolerances::millimetres()).unwrap(),
        )
    }

    fn horizontal(z: f64) -> Plane {
        plane(Point::new(0.0, 0.0, z), Vector::Z)
    }

    /// The area of a regular `n`-gon of circumradius `r`.
    fn polygon_area(n: u32, r: f64) -> f64 {
        f64::from(n) / 2.0 * r * r * (TAU / f64::from(n)).sin()
    }

    /// A closed prism about the `z` axis: a regular `sides`-gon of
    /// circumradius `radius` from `z = 0` to `z = height`, its wall
    /// `rings - 1` rows of quads, its ends fans about the axis. Wound
    /// outward, or inward where not `outward`. Wall row `j`, column `i` is
    /// faces `2 (j sides + i)` and the one after.
    fn prism(radius: f64, height: f64, sides: u32, rings: u32, outward: bool) -> Triangulation {
        let mut points = Vec::new();
        for j in 0..rings {
            let z = height * f64::from(j) / f64::from(rings - 1);
            for i in 0..sides {
                let a = TAU * f64::from(i) / f64::from(sides);
                points.push(Point::new(radius * a.cos(), radius * a.sin(), z));
            }
        }
        let at = |i: u32, j: u32| j * sides + i % sides;
        let mut triangles = Vec::new();
        for j in 0..rings - 1 {
            for i in 0..sides {
                let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                triangles.extend([[a, b, c], [a, c, d]]);
            }
        }
        let (low, high) = (rings * sides, rings * sides + 1);
        points.extend([Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, height)]);
        for i in 0..sides {
            triangles.push([low, at(i + 1, 0), at(i, 0)]);
            triangles.push([high, at(i, rings - 1), at(i + 1, rings - 1)]);
        }
        if !outward {
            for t in &mut triangles {
                t.swap(1, 2);
            }
        }
        mesh(points, triangles)
    }

    /// A closed tube about the `z` axis between two radii, two rings high,
    /// with annular ends.
    fn tube(inner: f64, outer: f64, height: f64, sides: u32) -> Triangulation {
        let mut points = Vec::new();
        for r in [outer, inner] {
            for z in [0.0, height] {
                for i in 0..sides {
                    let a = TAU * f64::from(i) / f64::from(sides);
                    points.push(Point::new(r * a.cos(), r * a.sin(), z));
                }
            }
        }
        // Wall `w` (outer 0, inner 1), end `e` (bottom 0, top 1), column `i`.
        let at = |w: u32, e: u32, i: u32| (2 * w + e) * sides + i % sides;
        let mut triangles = Vec::new();
        for i in 0..sides {
            let (a, b, c, d) = (at(0, 0, i), at(0, 0, i + 1), at(0, 1, i + 1), at(0, 1, i));
            triangles.extend([[a, b, c], [a, c, d]]);
            let (a, b, c, d) = (at(1, 0, i), at(1, 0, i + 1), at(1, 1, i + 1), at(1, 1, i));
            triangles.extend([[a, c, b], [a, d, c]]);
            let (o, p, q, r) = (at(0, 0, i), at(0, 0, i + 1), at(1, 0, i + 1), at(1, 0, i));
            triangles.extend([[o, r, q], [o, q, p]]);
            let (o, p, q, r) = (at(0, 1, i), at(0, 1, i + 1), at(1, 1, i + 1), at(1, 1, i));
            triangles.extend([[o, p, q], [o, q, r]]);
        }
        mesh(points, triangles)
    }

    /// A sphere of latitude rings and longitude sides about poles on `z`,
    /// symmetric across `z = 0` to the bit, with a ring in that plane when
    /// `bands` is even.
    fn sphere(radius: f64, bands: u32, sides: u32) -> Triangulation {
        let mut points = vec![Point::new(0.0, 0.0, radius)];
        for j in 1..bands {
            let mirror = bands - j;
            let polar = PI * f64::from(j.min(mirror)) / f64::from(bands);
            let (s, c) = (radius * polar.sin(), radius * polar.cos());
            let z = match (2 * j).cmp(&bands) {
                core::cmp::Ordering::Less => c,
                core::cmp::Ordering::Equal => 0.0,
                core::cmp::Ordering::Greater => -c,
            };
            for i in 0..sides {
                let a = TAU * f64::from(i) / f64::from(sides);
                points.push(Point::new(s * a.cos(), s * a.sin(), z));
            }
        }
        let south = 1 + (bands - 1) * sides;
        points.push(Point::new(0.0, 0.0, -radius));
        let at = |j: u32, i: u32| 1 + (j - 1) * sides + i % sides;
        let mut triangles = Vec::new();
        for i in 0..sides {
            triangles.push([0, at(1, i), at(1, i + 1)]);
            triangles.push([south, at(bands - 1, i + 1), at(bands - 1, i)]);
        }
        for j in 1..bands - 1 {
            for i in 0..sides {
                let (a, b, c, d) = (at(j, i), at(j + 1, i), at(j + 1, i + 1), at(j, i + 1));
                triangles.extend([[a, b, c], [a, c, d]]);
            }
        }
        mesh(points, triangles)
    }

    /// The volume a closed mesh bounds, by the divergence theorem.
    fn volume(e: &Editable) -> f64 {
        let m = e.to_triangulation();
        m.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| m.positions[v as usize] - Point::ORIGIN);
                a.dot(b.cross(c)) / 6.0
            })
            .sum()
    }

    /// The area of the faces lying in `plane`, every corner on it to
    /// within `1e-12`: cut vertices lie within a few rounding units of it,
    /// at coordinates of a few units, and every other corner of the meshes
    /// here at least `1e-3` off it.
    fn area_in(e: &Editable, plane: &Plane) -> f64 {
        let n = plane.normal().vector();
        let height = |p: Point| (p - plane.origin()).dot(n);
        let m = e.to_triangulation();
        m.triangles
            .iter()
            .map(|t| t.map(|v| m.positions[v as usize]))
            .filter(|t| t.iter().all(|&p| height(p).abs() <= 1e-12))
            .map(|[a, b, c]| (b - a).cross(c - a).magnitude() / 2.0)
            .sum()
    }

    /// Equal to within `1e-12` of the larger: sums of a few hundred terms,
    /// each rounded by a few units in `2^-53`.
    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
    }

    const CAPPED: Option<Cap> = Some(Cap { bridge: 0.0 });

    fn one_loop() -> CutReport {
        CutReport {
            capped_loops: 1,
            ..CutReport::default()
        }
    }

    #[test]
    fn a_closed_prism_cut_across_is_a_closed_prism_of_the_height_kept() {
        let (r, n) = (2.0, 24);
        let whole = Editable::from(&prism(r, 3.0, n, 4, true), 0.0).unwrap();
        let cut = horizontal(1.3);
        for (keep, height) in [(Side::Below, 1.3), (Side::Above, 3.0 - 1.3)] {
            let mut e = whole.clone();
            assert_eq!(e.cut(&cut, keep, CAPPED).unwrap(), one_loop());
            check(&e);
            assert!(e.is_closed(), "{keep:?}");
            assert!(e.non_manifold().is_empty());
            assert!(close(volume(&e), polygon_area(n, r) * height), "{keep:?}");
            assert!(close(area_in(&e, &cut), polygon_area(n, r)));
        }
    }

    #[test]
    fn a_sphere_cut_through_its_equator_ring_is_two_closed_halves() {
        let (r, n) = (1.5, 24);
        let m = sphere(r, 16, n);
        let whole = Editable::from(&m, 0.0).unwrap();
        let all = volume(&whole);
        for keep in [Side::Below, Side::Above] {
            let mut e = whole.clone();
            assert_eq!(e.cut(&horizontal(0.0), keep, CAPPED).unwrap(), one_loop());
            check(&e);
            assert!(e.is_closed());
            // The ring in the plane is the rim: nothing new is made.
            assert_eq!(e.vertex_count(), m.positions.len());
            assert!(close(volume(&e), all / 2.0), "{keep:?}");
            assert!(close(area_in(&e, &horizontal(0.0)), polygon_area(n, r)));
        }
    }

    #[test]
    fn a_sphere_cut_aslant_caps_both_halves_with_the_area_of_its_section() {
        let m = sphere(1.5, 15, 22);
        let whole = Editable::from(&m, 0.0).unwrap();
        let cut = plane(Point::new(0.1, -0.2, 0.3), Vector::new(1.0, 2.0, 3.0));
        let s = section(&m, &cut).unwrap();
        assert_eq!(s.chains.len(), 1);
        let n = cut.normal().vector();
        let points = &s.chains[0].points;
        let enclosed: f64 = points
            .windows(2)
            .skip(1)
            .map(|w| (w[0] - points[0]).cross(w[1] - points[0]).dot(n))
            .sum::<f64>()
            / 2.0;
        let mut halves = 0.0;
        for keep in [Side::Below, Side::Above] {
            let mut e = whole.clone();
            assert_eq!(e.cut(&cut, keep, CAPPED).unwrap(), one_loop());
            check(&e);
            assert!(e.is_closed());
            assert!(close(area_in(&e, &cut), enclosed.abs()), "{keep:?}");
            halves += volume(&e);
        }
        assert!(close(halves, volume(&whole)));
    }

    #[test]
    fn the_cut_rim_is_the_section_to_the_bit() {
        let m = torus(24, 12);
        let cut = plane(Point::new(0.3, 0.1, 0.2), Vector::new(0.3, -0.2, 1.0));
        let bits = |p: Point| [p.x, p.y, p.z].map(f64::to_bits);
        let chains = section(&m, &cut).unwrap().chains;
        let mut chain: Vec<[u64; 3]> = chains
            .iter()
            .flat_map(|c| c.points.iter().map(|&p| bits(p)))
            .collect();
        chain.sort_unstable();
        chain.dedup();
        for keep in [Side::Below, Side::Above] {
            let mut e = Editable::from(&m, 0.0).unwrap();
            let report = e.cut(&cut, keep, None).unwrap();
            assert_eq!(report.open, chains.len());
            let mut rim: Vec<[u64; 3]> = e
                .rims()
                .iter()
                .flat_map(|r| r.edges.iter().map(|&h| e.edge(h).unwrap()[0]))
                .map(|v| bits(e.position(v).unwrap()))
                .collect();
            rim.sort_unstable();
            rim.dedup();
            assert_eq!(rim, chain, "{keep:?}");
        }
    }

    #[test]
    fn kept_vertices_keep_their_numbers_and_sources() {
        let m = torus(24, 12);
        let cut = plane(Point::new(0.3, 0.1, 0.2), Vector::new(0.3, -0.2, 1.0));
        let n = cut.normal().vector();
        let mut e = Editable::from(&m, 0.0).unwrap();
        e.cut(&cut, Side::Above, CAPPED).unwrap();
        check(&e);
        assert!(e.is_closed());
        for (v, &p) in m.positions.iter().enumerate() {
            let v = index(v);
            // Clear of the plane by far more than rounding: on the kept side.
            if (p - cut.origin()).dot(n) > 1e-9 {
                assert_eq!(e.position(v), Some(p));
                assert_eq!(e.source_vertex(v), Some(v));
            }
        }
        for v in 0..index(e.vertex_count()) {
            let Some(p) = e.position(v) else { continue };
            match e.source_vertex(v) {
                Some(s) => assert_eq!(p, m.positions[s as usize]),
                None => {
                    assert!(v as usize >= m.positions.len());
                    assert!((p - cut.origin()).dot(n).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn a_tube_cut_across_is_capped_by_a_ring() {
        let (inner, outer, n) = (1.0, 2.5, 20);
        let whole = Editable::from(&tube(inner, outer, 2.0, n), 0.0).unwrap();
        assert!(whole.is_closed());
        let cut = horizontal(0.7);
        let mut e = whole;
        let report = e.cut(&cut, Side::Below, CAPPED).unwrap();
        assert_eq!(
            report,
            CutReport {
                capped_loops: 2,
                ..CutReport::default()
            }
        );
        check(&e);
        assert!(e.is_closed());
        assert_eq!(e.components().len(), 1);
        let ring = polygon_area(n, outer) - polygon_area(n, inner);
        assert!(close(area_in(&e, &cut), ring));
        assert!(close(volume(&e), ring * 0.7));
    }

    #[test]
    fn nested_shells_give_nested_caps() {
        let n = 18;
        let mut e = Editable::from(&nested(n), 0.0).unwrap();
        let cut = horizontal(1.4);
        let report = e.cut(&cut, Side::Above, CAPPED).unwrap();
        assert_eq!(
            report,
            CutReport {
                capped_loops: 3,
                ..CutReport::default()
            }
        );
        check(&e);
        assert!(e.is_closed());
        assert_eq!(e.components().len(), 2);
        let solid = polygon_area(n, 3.0) - polygon_area(n, 2.0) + polygon_area(n, 1.0);
        assert!(close(area_in(&e, &cut), solid));
        let above =
            polygon_area(n, 3.0) * 1.6 - polygon_area(n, 2.0) * 1.1 + polygon_area(n, 1.0) * 0.6;
        assert!(close(volume(&e), above));
    }

    #[test]
    fn a_gap_in_the_rim_is_bridged_within_reach_and_left_open_beyond() {
        // A prism with one quad of its middle row missing; the cut runs
        // through that row, so its rim in the plane stops either side of
        // the gap, one side of the polygon long.
        let (r, n) = (2.0, 16);
        let mut whole = Editable::from(&prism(r, 3.0, n, 4, true), 0.0).unwrap();
        whole.remove_faces(&[2 * n, 2 * n + 1]);
        let side = 2.0 * r * (PI / f64::from(n)).sin();
        let cut = horizontal(1.5);

        let mut e = whole.clone();
        let report = e
            .cut(
                &cut,
                Side::Below,
                Some(Cap {
                    bridge: side * 1.01,
                }),
            )
            .unwrap();
        assert_eq!(
            report,
            CutReport {
                capped_loops: 1,
                bridged: 1,
                open: 0
            }
        );
        check(&e);
        assert!(close(area_in(&e, &cut), polygon_area(n, r)));
        // What is left open is the gap below the plane, closed by the
        // bridge: four edges.
        let rims = e.rims();
        assert_eq!(rims.len(), 1);
        assert_eq!(rims[0].edges.len(), 4);

        let mut e = whole;
        let report = e
            .cut(
                &cut,
                Side::Below,
                Some(Cap {
                    bridge: side * 0.99,
                }),
            )
            .unwrap();
        assert_eq!(
            report,
            CutReport {
                open: 1,
                ..CutReport::default()
            }
        );
        check(&e);
        assert_eq!(area_in(&e, &cut), 0.0);
    }

    /// A solid drum holding a void, holding a solid drum, each clear of
    /// the next: radii 3, 2 and 1, from `z = 0`, `0.5` and `1` to `z = 3`,
    /// `2.5` and `2`.
    fn nested(sides: u32) -> Triangulation {
        let mut m = prism(3.0, 3.0, sides, 4, true);
        for (r, base, outward) in [(2.0, 0.5, false), (1.0, 1.0, true)] {
            let mut d = prism(r, 3.0 - 2.0 * base, sides, 4, outward);
            for p in &mut d.positions {
                *p += Vector::new(0.0, 0.0, base);
            }
            m.append(&d);
        }
        m
    }

    /// `m` moved by `x` along the `x` axis.
    fn shifted(mut m: Triangulation, x: f64) -> Triangulation {
        for p in &mut m.positions {
            *p += Vector::new(x, 0.0, 0.0);
        }
        m
    }

    #[test]
    fn of_two_overlapping_bodies_one_is_left_open() {
        // Two drums overlapping, and a third apart from both.
        let mut m = prism(1.0, 2.0, 12, 3, true);
        m.append(&shifted(prism(1.0, 2.0, 12, 3, true), 1.2));
        m.append(&shifted(prism(1.0, 2.0, 12, 3, true), 5.0));
        let mut e = Editable::from(&m, 0.0).unwrap();
        let report = e.cut(&horizontal(0.5), Side::Below, CAPPED).unwrap();
        assert_eq!(
            report,
            CutReport {
                capped_loops: 2,
                open: 1,
                ..CutReport::default()
            }
        );
        check(&e);
        assert_eq!(e.rims().len(), 1);
    }

    #[test]
    fn bodies_meeting_at_an_edge_are_capped_through_the_point_they_share() {
        // Two square drums whose walls meet along one upright edge, their
        // corners on the axes exactly. Welded, the edge is one; apart, the
        // two edges are cut at one position, one point.
        let mut m = prism(1.0, 2.0, 4, 3, true);
        m.append(&shifted(prism(1.0, 2.0, 4, 3, true), 2.0));
        for p in &mut m.positions {
            *p = Point::new(p.x.round(), p.y.round(), p.z);
        }
        let cut = horizontal(0.5);
        for weld in [0.0, 1e-9] {
            let mut e = Editable::from(&m, weld).unwrap();
            let report = e.cut(&cut, Side::Below, CAPPED).unwrap();
            assert_eq!(
                report,
                CutReport {
                    capped_loops: 2,
                    ..CutReport::default()
                },
                "{weld}"
            );
            check(&e);
            assert!(e.is_closed());
            assert!(close(volume(&e), 2.0 * 2.0 * 0.5));
        }
    }

    #[test]
    fn a_vertex_a_rounding_above_the_plane_is_cut_at_not_beside() {
        // An octahedron whose waist stands `1e-20` above `z = 0`: every
        // edge from a pole to the waist crosses the plane within rounding
        // of the waist.
        let tiny = 1e-20;
        let points = vec![
            Point::new(1.0, 0.0, tiny),
            Point::new(0.0, 1.0, tiny),
            Point::new(-1.0, 0.0, tiny),
            Point::new(0.0, -1.0, tiny),
            Point::new(0.0, 0.0, 1.0),
            Point::new(0.0, 0.0, -1.0),
        ];
        let mut triangles = Vec::new();
        for i in 0..4 {
            let j = (i + 1) % 4;
            triangles.extend([[i, j, 4], [j, i, 5]]);
        }
        let m = mesh(points, triangles);
        for keep in [Side::Below, Side::Above] {
            let mut e = Editable::from(&m, 0.0).unwrap();
            assert_eq!(e.cut(&horizontal(0.0), keep, CAPPED).unwrap(), one_loop());
            check(&e);
            assert!(e.is_closed());
            assert_eq!(e.vertex_count(), 6, "{keep:?}");
            assert!(close(volume(&e), 2.0 / 3.0), "{keep:?}");
        }
    }

    #[test]
    fn an_open_sheet_cut_across_leaves_its_stretch_open() {
        // A cylinder's wall with no ends.
        let mut e = Editable::from(&prism(1.0, 2.0, 12, 3, true), 0.0).unwrap();
        e.remove_faces(&(48..72).collect::<Vec<_>>());
        e.remove_faces(&[0, 1]);
        let report = e.cut(&horizontal(0.5), Side::Below, CAPPED).unwrap();
        assert_eq!(
            report,
            CutReport {
                open: 1,
                ..CutReport::default()
            }
        );
        check(&e);
    }

    #[test]
    fn a_plane_missing_the_mesh_keeps_all_or_nothing() {
        let m = sphere(1.0, 6, 8);
        let mut e = Editable::from(&m, 0.0).unwrap();
        assert_eq!(
            e.cut(&horizontal(2.0), Side::Below, CAPPED).unwrap(),
            CutReport::default()
        );
        assert_eq!(e.to_triangulation().triangles, m.triangles);
        assert_eq!(
            e.cut(&horizontal(2.0), Side::Above, CAPPED).unwrap(),
            CutReport::default()
        );
        assert!(e.to_triangulation().triangles.is_empty());
        check(&e);
    }

    #[test]
    fn a_plane_through_an_end_keeps_the_end_with_its_material() {
        let (r, n, h) = (1.0, 8, 2.0);
        let m = prism(r, h, n, 3, true);
        let full = polygon_area(n, r) * h;
        let cut = |z: f64, keep: Side| {
            let mut e = Editable::from(&m, 0.0).unwrap();
            let report = e.cut(&horizontal(z), keep, CAPPED).unwrap();
            check(&e);
            (e, report)
        };
        // The top faces up: below keeps the body and caps it afresh, above
        // keeps nothing.
        let (e, report) = cut(h, Side::Below);
        assert_eq!(report, one_loop());
        assert!(e.is_closed());
        assert!(close(volume(&e), full));
        let (e, report) = cut(h, Side::Above);
        assert_eq!(report, CutReport::default());
        assert!(e.to_triangulation().triangles.is_empty());
        // The bottom faces down: above keeps the whole, unchanged, below
        // nothing.
        let (e, report) = cut(0.0, Side::Above);
        assert_eq!(report, CutReport::default());
        assert_eq!(e.to_triangulation().triangles, m.triangles);
        let (e, report) = cut(0.0, Side::Below);
        assert_eq!(report, CutReport::default());
        assert!(e.to_triangulation().triangles.is_empty());
    }

    #[test]
    fn the_cut_is_deterministic() {
        let m = torus(30, 14);
        let cut = plane(Point::new(0.0, 0.5, 0.1), Vector::new(1.0, 1.0, 0.4));
        let mut a = Editable::from(&m, 0.0).unwrap();
        let mut b = a.clone();
        let ra = a.cut(&cut, Side::Below, CAPPED).unwrap();
        let rb = b.cut(&cut, Side::Below, CAPPED).unwrap();
        assert_eq!(ra, rb);
        assert_eq!(a, b);
        assert!(a.is_closed());
    }

    #[test]
    fn a_bad_plane_or_bridge_is_refused() {
        let mut e = Editable::from(&sphere(1.0, 6, 8), 0.0).unwrap();
        let far = plane(Point::new(f64::NAN, 0.0, 0.0), Vector::Z);
        assert!(e.cut(&far, Side::Below, None).is_err());
        for bridge in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(
                e.cut(&horizontal(0.0), Side::Below, Some(Cap { bridge }))
                    .is_err()
            );
        }
        assert!(e.is_closed());
    }

    /// Planes through or a rounding from vertices, each of which once left
    /// a part open: through a ring a rounding off an axis plane, through a
    /// vertex others round onto, and along a drum's flat side.
    #[test]
    fn planes_grazing_vertices_close_both_parts() {
        let shift = 0.001_455_055_119_429_773_2;
        let cases = [
            (torus(12, 8), 0, Vector::new(-1.0, -1.0, 1.0), 0.0),
            (torus(12, 8), 5, Vector::new(0.0, -1.0, 0.0), 0.0),
            (torus(12, 8), 91, Vector::new(-1.0, 0.0, 0.0), 0.0),
            (torus(12, 8), 0, Vector::new(-1.0, 3.0, 0.0), shift),
            (nested(12), 5, Vector::new(0.0, -2.0, 0.0), 0.0),
            (nested(12), 19, Vector::new(3.0, 3.0, 0.0), 0.0),
        ];
        for (m, through, normal, shift) in cases {
            let origin = m.positions[through] + Vector::new(shift, shift / 3.0, shift / 7.0);
            let cut = plane(origin, normal);
            let whole = Editable::from(&m, 0.0).unwrap();
            let mut total = 0.0;
            for keep in [Side::Below, Side::Above] {
                let mut e = whole.clone();
                let report = e.cut(&cut, keep, CAPPED).unwrap();
                check(&e);
                assert_eq!(report.open, 0, "{through} {normal:?} {keep:?}");
                assert!(e.is_closed() || e.to_triangulation().triangles.is_empty());
                total += volume(&e);
            }
            assert!(close(total, volume(&whole)), "{through} {normal:?}");
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        /// Any plane, through a vertex or not, cuts a closed mesh, a torus,
        /// a sphere or drums nested in a drum, into two closed parts, capped
        /// without a loop left open, whose volumes add up to the whole.
        #[test]
        fn a_closed_mesh_cut_anywhere_closes_both_parts(
            shape in 0_usize..3,
            through in 0_usize..96,
            on_vertex in any::<bool>(),
            shift in -0.5_f64..0.5,
            normal in (-3_i32..=3, -3_i32..=3, -3_i32..=3),
        ) {
            let (x, y, z) = normal;
            prop_assume!((x, y, z) != (0, 0, 0));
            let m = match shape {
                0 => torus(12, 8),
                1 => sphere(1.5, 8, 14),
                _ => nested(12),
            };
            let at = m.positions[through % m.positions.len()];
            let origin = if on_vertex {
                at
            } else {
                at + Vector::new(shift, shift / 3.0, shift / 7.0)
            };
            let cut = plane(origin, Vector::new(f64::from(x), f64::from(y), f64::from(z)));
            let whole = Editable::from(&m, 0.0).unwrap();
            let mut total = 0.0;
            for keep in [Side::Below, Side::Above] {
                let mut e = whole.clone();
                let report = e.cut(&cut, keep, CAPPED).unwrap();
                check(&e);
                prop_assert_eq!(report.open, 0);
                prop_assert!(e.is_closed() || e.to_triangulation().triangles.is_empty());
                total += volume(&e);
            }
            let all = volume(&whole);
            // Both parts clip each face at the same vertices and cap the
            // same loops, so they split the whole exactly but for rounding.
            prop_assert!(close(total, all), "{} against {}", total, all);
        }
    }
}
