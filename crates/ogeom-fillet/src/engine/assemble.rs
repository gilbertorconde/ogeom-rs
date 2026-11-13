//! P4: one assembly of the result from its rings.
//!
//! The polyhedral classes hand over planar rings over shared points: every
//! shrunk face, every strip, every vertex facet. Here the rings are cleaned
//! (corners repeated in turn, rings of no area), a point lying inside
//! another ring's edge is put into that edge, rings of one plane sharing an
//! edge are joined (facing opposite ways, the smaller is a corner the larger
//! loses), and the solid is built once over the shared points, each face
//! keeping its own plane. Curved pieces (M2's strips and vertex patches) are
//! built by [`build_faces`](super::curved::build_faces), over edges and
//! surfaces the caller gives.

use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_math::{Point, Vector};

use super::Origin;

/// A ring of the result: its points, where it comes from, and the way it
/// must face where that is known.
pub(crate) type Ring = (Vec<usize>, Vec<Origin>, Option<Vector>);

/// The result's rings, over points shared by index.
pub(crate) struct Rings {
    points: Vec<Point>,
    grid: FastMap<(i64, i64, i64), Vec<usize>>,
    cell: f64,
    rings: Vec<Ring>,
    /// Each ring's holes, wound the other way round.
    holes: Vec<Vec<Vec<usize>>>,
    tol: Tolerances,
}

/// The rings ready for a solid: the points they use, each ring wound out of
/// the solid, its holes wound the other way, and each ring's origins.
pub(crate) struct Assembled {
    pub points: Vec<Point>,
    pub rings: Vec<Vec<usize>>,
    pub holes: Vec<Vec<Vec<usize>>>,
    pub origins: Vec<Vec<Origin>>,
}

impl Rings {
    pub fn new(tol: Tolerances) -> Self {
        Self {
            points: Vec::new(),
            grid: FastMap::default(),
            cell: tol.confusion() * 16.0,
            rings: Vec::new(),
            holes: Vec::new(),
            tol,
        }
    }

    /// The index of `p`, one already held within the confusion distance or
    /// a new one.
    pub fn point(&mut self, p: Point) -> usize {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "a grid cell of a finite coordinate over a positive cell"
        )]
        let key = |x: f64, cell: f64| (x / cell).floor() as i64;
        let k = (
            key(p.x, self.cell),
            key(p.y, self.cell),
            key(p.z, self.cell),
        );
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(held) = self.grid.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        for &i in held {
                            if self.points[i].distance(p) <= self.tol.confusion() {
                                return i;
                            }
                        }
                    }
                }
            }
        }
        self.points.push(p);
        let i = self.points.len() - 1;
        self.grid.entry(k).or_default().push(i);
        i
    }

    /// A ring, where it comes from, and the way it must face where that is
    /// known: a shrunk face as the face did, a bevel's quad along the sum of
    /// its two faces' normals.
    pub fn push(&mut self, ring: Vec<usize>, origin: Origin, facing: Option<Vector>) {
        self.push_holed(ring, Vec::new(), origin, facing);
    }

    /// [`Rings::push`] for a ring with holes, each wound the other way.
    pub fn push_holed(
        &mut self,
        ring: Vec<usize>,
        holes: Vec<Vec<usize>>,
        origin: Origin,
        facing: Option<Vector>,
    ) {
        self.rings.push((ring, vec![origin], facing));
        self.holes.push(holes);
    }

    /// Every ring holding a point of `allowed` clipped by each of `cuts`
    /// but its own, keeping the side each cut's normal points away from:
    /// a ring's part beyond a cut is taken off, the cut's line across it
    /// put in. `false` where a clip would take a corner outside `allowed`,
    /// or one `may_take` refuses to the cut given the ring and the vertex
    /// the corner stands at (points the clips make join `allowed` with no
    /// vertex as they are made). A ring holding a hole is outside the clip.
    pub fn clip_local(
        &mut self,
        cuts: &[super::consumed::Cut],
        allowed: &mut FastMap<usize, Option<super::Id>>,
        may_take: impl Fn(&Origin, Option<super::Id>, &Origin) -> bool,
    ) -> bool {
        let eps = self.tol.confusion();
        let local: Vec<usize> = (0..self.rings.len())
            .filter(|ri| self.rings[*ri].0.iter().any(|p| allowed.contains_key(p)))
            .collect();
        for ri in local {
            if !self.holes[ri].is_empty() {
                return false;
            }
            for (origin, n, off) in cuts {
                if self.rings[ri].1.iter().any(|o| same_origin(o, origin)) {
                    continue;
                }
                let ring = std::mem::take(&mut self.rings[ri].0);
                let side = |me: &Self, i: usize| n.dot(me.points[i].to_vector()) - off;
                let mut out: Vec<usize> = Vec::with_capacity(ring.len() + 2);
                for k in 0..ring.len() {
                    let (here, next) = (ring[k], ring[(k + 1) % ring.len()]);
                    let (sh, sn) = (side(self, here), side(self, next));
                    if sh <= eps {
                        out.push(here);
                    } else {
                        let ring_origin = self.rings[ri].1[0];
                        match allowed.get(&here) {
                            Some(at) if may_take(&ring_origin, *at, origin) => {}
                            _ => return false,
                        }
                    }
                    if (sh > eps && sn < -eps) || (sh < -eps && sn > eps) {
                        let (ph, pn) = (self.points[here], self.points[next]);
                        let x = self.point(ph + (pn - ph) * (sh / (sh - sn)));
                        allowed.entry(x).or_insert(None);
                        out.push(x);
                    }
                }
                out.dedup();
                while out.len() > 1 && out.first() == out.last() {
                    out.pop();
                }
                self.rings[ri].0 = out;
            }
        }
        true
    }

    /// A loop with corners repeated in turn dropped, or `None` where it
    /// keeps no area.
    fn cleaned(&self, ring: Vec<usize>) -> Option<Vec<usize>> {
        let tol = self.tol;
        let mut clean: Vec<usize> = Vec::with_capacity(ring.len());
        for i in ring {
            if clean.last() != Some(&i) {
                clean.push(i);
            }
        }
        while clean.len() > 1 && clean.first() == clean.last() {
            clean.pop();
        }
        if clean.len() < 3 {
            return None;
        }
        let n = self.newell(&clean);
        (n.magnitude() > tol.confusion() * tol.confusion()).then_some(clean)
    }

    /// A ring's normal by Newell's sum, twice its area long.
    fn newell(&self, ring: &[usize]) -> Vector {
        let mut n = Vector::ZERO;
        let a = self.points[ring[0]];
        for k in 1..ring.len().saturating_sub(1) {
            n += (self.points[ring[k]] - a).cross(self.points[ring[k + 1]] - a);
        }
        n
    }

    /// The rings cleaned and checked: corners repeated in turn dropped,
    /// rings of no area dropped, a point lying on another ring's edge put
    /// into that edge, rings of one plane sharing an edge joined, and the
    /// whole refused where a face turned inside out or two faces cross.
    pub fn assemble(mut self) -> OgeomResult<Option<Assembled>> {
        let tol = self.tol;
        let mut rings: Vec<Ring> = Vec::new();
        let mut holes: Vec<Vec<Vec<usize>>> = Vec::new();
        let taken = std::mem::take(&mut self.rings);
        let taken_holes = std::mem::take(&mut self.holes);
        for ((ring, origin, facing), inner) in taken.into_iter().zip(taken_holes) {
            let Some(clean) = self.cleaned(ring) else {
                continue;
            };
            rings.push((clean, origin, facing));
            holes.push(inner.into_iter().filter_map(|h| self.cleaned(h)).collect());
        }
        self.split_on_edges(&mut rings, &mut holes);
        self.join_coplanar(&mut rings, &mut holes);
        self.drop_straight_corners(&mut rings, &mut holes);
        if !super::guards::rings_fit_holed(&self.points, &rings, &holes, tol)? {
            return Ok(None);
        }
        // Only the points the rings use, in first use.
        let mut renumber: FastMap<usize, usize> = FastMap::default();
        let mut points = Vec::new();
        let mut renumbered = |ring: &[usize]| -> Vec<usize> {
            ring.iter()
                .map(|i| {
                    *renumber.entry(*i).or_insert_with(|| {
                        points.push(self.points[*i]);
                        points.len() - 1
                    })
                })
                .collect()
        };
        let mut out_rings = Vec::with_capacity(rings.len());
        let mut out_holes = Vec::with_capacity(rings.len());
        let mut origins = Vec::with_capacity(rings.len());
        for ((ring, origin, _), inner) in rings.into_iter().zip(holes) {
            out_rings.push(renumbered(&ring));
            out_holes.push(inner.iter().map(|h| renumbered(h)).collect());
            origins.push(origin);
        }
        Ok(Some(Assembled {
            points,
            rings: out_rings,
            holes: out_holes,
            origins,
        }))
    }

    /// Every point the rings use that lies inside another ring's edge, put
    /// into that edge in order along it.
    fn split_on_edges(&self, rings: &mut [Ring], holes: &mut [Vec<Vec<usize>>]) {
        let used: Vec<usize> = {
            let mut set: Vec<usize> = rings
                .iter()
                .flat_map(|(r, _, _)| r.iter().copied())
                .chain(holes.iter().flatten().flatten().copied())
                .collect();
            set.sort_unstable();
            set.dedup();
            set
        };
        for (ring, _, _) in rings.iter_mut() {
            *ring = self.split_loop(ring, &used);
        }
        for hole in holes.iter_mut().flatten() {
            *hole = self.split_loop(hole, &used);
        }
    }

    /// One loop with every point of `used` lying inside one of its edges
    /// put into that edge, in order along it.
    fn split_loop(&self, ring: &[usize], used: &[usize]) -> Vec<usize> {
        let tol = self.tol;
        let mut out = Vec::with_capacity(ring.len());
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            out.push(a);
            let (pa, pb) = (self.points[a], self.points[b]);
            let ab = pb - pa;
            let length = ab.magnitude();
            let mut inside: Vec<(f64, usize)> = used
                .iter()
                .copied()
                .filter(|c| *c != a && *c != b)
                .filter_map(|c| {
                    let t = (self.points[c] - pa).dot(ab) / (length * length);
                    let off = (self.points[c] - (pa + ab * t)).magnitude();
                    (t > 0.0 && t < 1.0 && off <= tol.confusion()).then_some((t, c))
                })
                .collect();
            inside.sort_by(|x, y| x.0.total_cmp(&y.0));
            out.extend(inside.into_iter().map(|(_, c)| c));
        }
        out
    }

    /// Rings of one plane that share an edge, joined into one across it,
    /// until no two do. Facing the same way they are one face; facing
    /// opposite ways the smaller is a piece the larger loses (a vertex's
    /// corner a bevel cuts off an end face), and joining takes it out. A
    /// ring joined away to nothing is dropped. A ring meeting a hole of
    /// another is joined into that hole the same way, the face's material
    /// lying to the left of each of its loops alike: the hole gains or
    /// loses the piece, and a hole joined away to nothing is dropped.
    fn join_coplanar(&self, rings: &mut Vec<Ring>, holes: &mut Vec<Vec<Vec<usize>>>) {
        let tol = self.tol;
        loop {
            let mut joined = false;
            'search: for i in 0..rings.len() {
                for j in (i + 1)..rings.len() {
                    let (ni, nj) = (self.newell(&rings[i].0), self.newell(&rings[j].0));
                    let (ui, uj) = (ni / ni.magnitude(), nj / nj.magnitude());
                    if ui.dot(uj).abs() < 1.0 - tol.angular()
                        || (self.points[rings[j].0[0]] - self.points[rings[i].0[0]])
                            .dot(ui)
                            .abs()
                            > tol.confusion() * 10.0
                    {
                        continue;
                    }
                    let mut origin = rings[i].1.clone();
                    origin.extend(rings[j].1.iter().copied());
                    if origin.iter().any(|o| matches!(o, Origin::Face(_))) {
                        origin.retain(|o| matches!(o, Origin::Face(_)));
                    }
                    let facing = rings[i].2.or(rings[j].2);
                    let facing = match (&rings[i].1[0], &rings[j].1[0]) {
                        (Origin::Face(_), _) => rings[i].2,
                        (_, Origin::Face(_)) => rings[j].2,
                        _ => facing,
                    };
                    if let Some(ring) = splice(&rings[i].0, &rings[j].0) {
                        // A face takes in what joins it, which then makes no
                        // face of its own; two faces joined are both the one.
                        let gone = ring.len() < 3
                            || self.newell(&ring).magnitude() <= tol.confusion() * tol.confusion();
                        rings.remove(j);
                        let inner = holes.remove(j);
                        if gone {
                            rings.remove(i);
                            holes.remove(i);
                        } else {
                            rings[i] = (ring, origin, facing);
                            holes[i].extend(inner);
                        }
                        joined = true;
                        break 'search;
                    }
                    // A ring of no holes joined into a hole of the other.
                    for (holed, plain) in [(i, j), (j, i)] {
                        if !holes[plain].is_empty() {
                            continue;
                        }
                        let found = holes[holed]
                            .iter()
                            .enumerate()
                            .find_map(|(k, h)| Some((k, splice(h, &rings[plain].0)?)));
                        let Some((k, hole)) = found else {
                            continue;
                        };
                        match self.cleaned(hole) {
                            Some(hole) => holes[holed][k] = hole,
                            None => {
                                holes[holed].remove(k);
                            }
                        }
                        rings[holed].1 = origin;
                        rings[holed].2 = facing;
                        rings.remove(plain);
                        holes.remove(plain);
                        joined = true;
                        break 'search;
                    }
                }
            }
            if !joined {
                break;
            }
        }
    }

    /// A corner where a ring runs straight on and no other ring turns,
    /// taken out of every ring holding it.
    fn drop_straight_corners(&self, rings: &mut [Ring], holes: &mut [Vec<Vec<usize>>]) {
        let tol = self.tol;
        let straight = |ring: &[usize], k: usize| -> bool {
            let n = ring.len();
            let (a, p, b) = (
                self.points[ring[(k + n - 1) % n]],
                self.points[ring[k]],
                self.points[ring[(k + 1) % n]],
            );
            let (u, w) = (p - a, b - p);
            u.cross(w).magnitude() <= tol.angular() * u.magnitude() * w.magnitude()
                && u.dot(w) > 0.0
        };
        let loops: Vec<&[usize]> = rings
            .iter()
            .map(|(r, _, _)| r.as_slice())
            .chain(holes.iter().flatten().map(Vec::as_slice))
            .collect();
        let mut holders: FastMap<usize, Vec<(usize, usize)>> = FastMap::default();
        for (ri, ring) in loops.iter().enumerate() {
            for (k, p) in ring.iter().enumerate() {
                holders.entry(*p).or_default().push((ri, k));
            }
        }
        let mut drop: FastSet<usize> = FastSet::default();
        for (p, held) in &holders {
            if held.iter().all(|(ri, k)| straight(loops[*ri], *k)) {
                drop.insert(*p);
            }
        }
        for (ring, _, _) in rings.iter_mut() {
            ring.retain(|p| !drop.contains(p));
        }
        for hole in holes.iter_mut().flatten() {
            hole.retain(|p| !drop.contains(p));
        }
    }
}

/// Two rings sharing one edge, walked once each way, joined across it:
/// `None` where they share none.
fn splice(a: &[usize], b: &[usize]) -> Option<Vec<usize>> {
    let (na, nb) = (a.len(), b.len());
    for i in 0..na {
        let (p, q) = (a[i], a[(i + 1) % na]);
        if let Some(j) = (0..nb).find(|j| b[*j] == q && b[(j + 1) % nb] == p) {
            // a up to p, then b from p round to q, then a from q on.
            let mut out: Vec<usize> = Vec::with_capacity(na + nb - 2);
            for k in 0..na {
                out.push(a[(i + 1 + k) % na]);
            }
            // out starts at q and ends at p; continue round b from p to q.
            for k in 2..nb {
                out.push(b[(j + k) % nb]);
            }
            return Some(without_spikes(out));
        }
    }
    None
}

/// Whether two origins name the same element.
fn same_origin(a: &Origin, b: &Origin) -> bool {
    match (a, b) {
        (Origin::Face(x), Origin::Face(y)) => x == y,
        (Origin::Edge(x), Origin::Edge(y)) | (Origin::Vertex(x), Origin::Vertex(y)) => x == y,
        _ => false,
    }
}

/// A ring with every there-and-back (`p, q, p`) taken out.
fn without_spikes(mut ring: Vec<usize>) -> Vec<usize> {
    loop {
        let n = ring.len();
        if n < 3 {
            return ring;
        }
        let spike = (0..n).find(|k| ring[(k + n - 1) % n] == ring[(k + 1) % n]);
        let Some(k) = spike else {
            return ring;
        };
        let next = (k + 1) % n;
        let (first, second) = if k < next { (next, k) } else { (k, next) };
        ring.remove(first);
        ring.remove(second);
    }
}
