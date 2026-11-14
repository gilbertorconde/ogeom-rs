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
    /// put in. Then each ring not one of a group's own loses the funnel
    /// above all of that group's planes ([`Rings::clip_funnel`]). `false`
    /// where a clip would take a corner outside `allowed`, or one
    /// `may_take` refuses to the cut given the ring and the vertex the
    /// corner stands at (to a group, to any of its bevels; points the clips
    /// make join `allowed` with no vertex as they are made). A ring holding
    /// a hole is outside the clip.
    pub fn clip_local(
        &mut self,
        cuts: &[super::consumed::Cut],
        groups: &[super::consumed::Group],
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
            let ring_origin = self.rings[ri].1[0];
            for g in groups {
                if self.rings[ri].1.iter().any(|o| g.owns(o)) {
                    continue;
                }
                let may = |at: Option<super::Id>| {
                    g.edges
                        .iter()
                        .any(|e| may_take(&ring_origin, at, &Origin::Edge(*e)))
                };
                if !self.clip_funnel(ri, &g.planes, allowed, &may) {
                    return false;
                }
            }
        }
        true
    }

    /// Ring `ri` less the funnel above every one of `planes` (each an
    /// outward normal and offset; the funnel is where a point stands past
    /// all of them): the ring keeps what lies under any. Where the ring's
    /// boundary runs into the funnel and out again, the funnel's own
    /// boundary across the ring is put in: along one plane, or along two
    /// through the corner where they and the ring's plane meet. `false`
    /// where a corner the funnel takes is not in `allowed` or `may` refuses
    /// it, where the way round the funnel is not along one plane or two
    /// meeting at a corner of it inside the ring, or where the funnel lies
    /// inside the ring, holing it.
    fn clip_funnel(
        &mut self,
        ri: usize,
        planes: &[(Vector, f64)],
        allowed: &mut FastMap<usize, Option<super::Id>>,
        may: &dyn Fn(Option<super::Id>) -> bool,
    ) -> bool {
        let eps = self.tol.confusion();
        let ring = self.rings[ri].0.clone();
        let n = ring.len();
        if n < 3 || planes.is_empty() {
            return true;
        }
        let past = |p: Point| {
            planes
                .iter()
                .map(|(nv, off)| nv.dot(p.to_vector()) - off)
                .fold(f64::INFINITY, f64::min)
        };
        let pts: Vec<Point> = ring.iter().map(|i| self.points[*i]).collect();
        let inside: Vec<bool> = pts.iter().map(|p| past(*p) > eps).collect();
        let normal = self.newell(&ring);
        if normal.magnitude() <= eps * eps {
            return true;
        }
        let normal = normal / normal.magnitude();
        let on_ring = normal.dot(pts[0].to_vector());
        // Where the ring's plane and planes `i` and `j` meet, if they do at
        // a point of the funnel's closure.
        let corner = |i: usize, j: usize| -> Option<Point> {
            let (a, da) = planes[i];
            let (b, db) = planes[j];
            let det = normal.dot(a.cross(b));
            if det.abs() <= 1e-9 {
                return None;
            }
            let x = (a.cross(b) * on_ring + b.cross(normal) * da + normal.cross(a) * db) / det;
            let x = Point::ORIGIN + x;
            (past(x) >= -eps).then_some(x)
        };
        // A ring point's place in the ring's plane, for inside tests.
        let ux = (pts[1] - pts[0]).cross(normal);
        let ux = ux / ux.magnitude();
        let uy = normal.cross(ux);
        let flat = |p: Point| ((p - pts[0]).dot(ux), (p - pts[0]).dot(uy));
        let outline: Vec<(f64, f64)> = pts.iter().map(|p| flat(*p)).collect();
        let within = |p: Point| {
            let q = flat(p);
            super::guards::strictly_inside(&outline, q, eps)
                || (0..n).any(|k| {
                    let (a, b) = (outline[k], outline[(k + 1) % n]);
                    let ab = (b.0 - a.0, b.1 - a.1);
                    let t = (((q.0 - a.0) * ab.0 + (q.1 - a.1) * ab.1)
                        / ab.0.mul_add(ab.0, ab.1 * ab.1))
                    .clamp(0.0, 1.0);
                    (a.0 + ab.0 * t - q.0).hypot(a.1 + ab.1 * t - q.1) <= eps
                })
        };
        // Each edge's run inside the funnel: where it enters and where it
        // leaves.
        let run = |a: Point, b: Point| -> Option<(f64, f64)> {
            let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
            for (nv, off) in planes {
                let (fa, fb) = (nv.dot(a.to_vector()) - off, nv.dot(b.to_vector()) - off);
                if fa <= eps && fb <= eps {
                    return None;
                }
                if fa > eps && fb > eps {
                    continue;
                }
                let t = fa / (fa - fb);
                if fa <= eps {
                    t0 = t0.max(t);
                } else {
                    t1 = t1.min(t);
                }
            }
            (t1 - t0 > 1e-12).then_some((t0, t1))
        };
        // The planes whose boundary a point of the funnel's closure is on.
        let on = |p: Point| -> Vec<usize> {
            (0..planes.len())
                .filter(|j| {
                    let (nv, off) = planes[*j];
                    (nv.dot(p.to_vector()) - off).abs() <= eps
                })
                .collect()
        };
        // The funnel's section in the ring's plane, counter-clockwise about
        // the ring's normal as the ring is: a square far round the ring
        // cut by each plane in turn.
        let section: Vec<(f64, f64)> = {
            let (lo, hi) = outline.iter().fold(
                ((f64::INFINITY, f64::INFINITY), (f64::NEG_INFINITY, f64::NEG_INFINITY)),
                |(lo, hi), q| ((lo.0.min(q.0), lo.1.min(q.1)), (hi.0.max(q.0), hi.1.max(q.1))),
            );
            let far = 10.0 * (hi.0 - lo.0).max(hi.1 - lo.1) + 1.0;
            let mut poly = vec![
                (lo.0 - far, lo.1 - far),
                (hi.0 + far, lo.1 - far),
                (hi.0 + far, hi.1 + far),
                (lo.0 - far, hi.1 + far),
            ];
            for (nv, off) in planes {
                let (a, b, c) = (
                    nv.dot(ux),
                    nv.dot(uy),
                    nv.dot(pts[0].to_vector()) - off,
                );
                let f = |q: (f64, f64)| a.mul_add(q.0, b.mul_add(q.1, c));
                let mut cut = Vec::with_capacity(poly.len() + 1);
                for k in 0..poly.len() {
                    let (here, next) = (poly[k], poly[(k + 1) % poly.len()]);
                    let (fh, fn_) = (f(here), f(next));
                    if fh >= 0.0 {
                        cut.push(here);
                    }
                    if (fh >= 0.0) != (fn_ >= 0.0) {
                        let t = fh / (fh - fn_);
                        cut.push((here.0 + (next.0 - here.0) * t, here.1 + (next.1 - here.1) * t));
                    }
                }
                poly = cut;
            }
            poly
        };
        // The section's corners passed going clockwise round it from `p` to
        // `q`, both on its boundary.
        let round = |p: Point, q: Point| -> Option<Vec<Point>> {
            let m = section.len();
            if m < 3 {
                return None;
            }
            let side_of = |x: Point| {
                let x = flat(x);
                (0..m)
                    .map(|k| {
                        let (a, b) = (section[k], section[(k + 1) % m]);
                        let ab = (b.0 - a.0, b.1 - a.1);
                        let t = (((x.0 - a.0) * ab.0 + (x.1 - a.1) * ab.1)
                            / ab.0.mul_add(ab.0, ab.1 * ab.1))
                        .clamp(0.0, 1.0);
                        ((a.0 + ab.0 * t - x.0).hypot(a.1 + ab.1 * t - x.1), k, t)
                    })
                    .min_by(|a, b| a.0.total_cmp(&b.0))
            };
            let (dp, mut k, tp) = side_of(p)?;
            let (dq, kq, tq) = side_of(q)?;
            if dp > eps || dq > eps || (k == kq && tq > tp) {
                return None;
            }
            let mut way = Vec::new();
            while k != kq {
                let (a, b) = section[k];
                way.push(pts[0] + ux * a + uy * b);
                k = (k + m - 1) % m;
                if way.len() > m {
                    return None;
                }
            }
            Some(way)
        };
        let Some(start) = (0..n).find(|k| !inside[*k]) else {
            // The whole ring lies in the funnel.
            for i in &ring {
                match allowed.get(i) {
                    Some(at) if may(*at) => {}
                    _ => return false,
                }
            }
            self.rings[ri].0 = Vec::new();
            return true;
        };
        let mut out: Vec<usize> = Vec::with_capacity(n + 4);
        let mut entered: Option<Vec<usize>> = None;
        let mut p_in = pts[0];
        let mut crossed = false;
        for step in 0..n {
            let k = (start + step) % n;
            let (a, b) = (pts[k], pts[(k + 1) % n]);
            if inside[k] {
                match allowed.get(&ring[k]) {
                    Some(at) if may(*at) => {}
                    _ => return false,
                }
            } else {
                out.push(ring[k]);
            }
            let Some((t0, t1)) = run(a, b) else {
                continue;
            };
            crossed = true;
            if !inside[k] {
                let p = a + (b - a) * t0;
                let x = self.point(p);
                allowed.entry(x).or_insert(None);
                out.push(x);
                entered = Some(on(p));
                p_in = p;
            }
            if t1 < 1.0 - 1e-12 || !inside[(k + 1) % n] {
                let Some(from) = entered.take() else {
                    return false;
                };
                let q = a + (b - a) * t1;
                let to = on(q);
                // Along a plane both ends are on; else round the funnel's
                // section the way that keeps the ring's rest on its left,
                // every corner passed inside the ring.
                if !from.iter().any(|j| to.contains(j)) {
                    let Some(way) = round(p_in, q) else {
                        return false;
                    };
                    for x in way {
                        if !within(x) {
                            return false;
                        }
                        let x = self.point(x);
                        allowed.entry(x).or_insert(None);
                        out.push(x);
                    }
                }
                let x = self.point(q);
                allowed.entry(x).or_insert(None);
                out.push(x);
            }
        }
        if entered.is_some() {
            return false;
        }
        if !crossed {
            // The ring's boundary keeps clear of the funnel: refused where
            // the funnel stands inside it.
            for i in 0..planes.len() {
                for j in i + 1..planes.len() {
                    if corner(i, j)
                        .is_some_and(|x| super::guards::strictly_inside(&outline, flat(x), eps))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        out.dedup();
        while out.len() > 1 && out.first() == out.last() {
            out.pop();
        }
        self.rings[ri].0 = out;
        true
    }

    /// Each ring of consumed face `face` split where it passes a corner
    /// twice, each loop kept but those keeping no area and those lying in
    /// one of `holes` (each a loop of the face's corners). The clip cuts a
    /// face whose holes it dropped back along its holes' walls, so what it
    /// leaves of the face can hold the room of a hole too, met at corners.
    /// `false` where a loop turns the other way.
    pub fn drop_in(&mut self, face: usize, normal: Vector, holes: &[Vec<Point>]) -> bool {
        let eps = self.tol.confusion();
        let ux = if normal.x.abs() < 0.9 { Vector::X } else { Vector::Y };
        let ux = ux - normal * ux.dot(normal);
        let ux = ux / ux.magnitude();
        let uy = normal.cross(ux);
        let flat = |p: Point| (p.to_vector().dot(ux), p.to_vector().dot(uy));
        let gaps: Vec<Vec<(f64, f64)>> = holes
            .iter()
            .map(|h| h.iter().map(|p| flat(*p)).collect())
            .collect();
        let mine: Vec<usize> = (0..self.rings.len())
            .filter(|ri| matches!(self.rings[*ri].1.as_slice(), [Origin::Face(f)] if *f == face))
            .collect();
        for ri in mine {
            let ring = std::mem::take(&mut self.rings[ri].0);
            let mut stack: Vec<usize> = Vec::with_capacity(ring.len());
            let mut loops: Vec<Vec<usize>> = Vec::new();
            for i in ring {
                if let Some(j) = stack.iter().position(|x| *x == i) {
                    loops.push(stack.split_off(j));
                }
                stack.push(i);
            }
            loops.push(stack);
            let mut kept: Vec<Vec<usize>> = Vec::new();
            for l in loops {
                if l.len() < 3 {
                    continue;
                }
                let area = self.newell(&l).dot(normal);
                if area.abs() <= eps * eps {
                    continue;
                }
                if area < 0.0 {
                    return false;
                }
                // A point just inside the loop off its longest side.
                let m = l.len();
                let Some(k) = (0..m).max_by(|x, y| {
                    let d = |k: usize| (self.points[l[(k + 1) % m]] - self.points[l[k]]).magnitude();
                    d(*x).total_cmp(&d(*y))
                }) else {
                    continue;
                };
                let (a, b) = (self.points[l[k]], self.points[l[(k + 1) % m]]);
                let inward = normal.cross(b - a);
                let at = a + (b - a) * 0.5 + inward * 1e-3;
                let outline: Vec<(f64, f64)> = l.iter().map(|i| flat(self.points[*i])).collect();
                if !super::guards::strictly_inside(&outline, flat(at), eps) {
                    return false;
                }
                if gaps
                    .iter()
                    .any(|g| super::guards::strictly_inside(g, flat(at), eps))
                {
                    continue;
                }
                kept.push(l);
            }
            let mut kept = kept.into_iter();
            self.rings[ri].0 = kept.next().unwrap_or_default();
            for l in kept {
                let (origin, facing) = (self.rings[ri].1.clone(), self.rings[ri].2);
                self.rings.push((l, origin, facing));
                self.holes.push(Vec::new());
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
        despike(&mut clean);
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
        // A clip's way round a cut can double back along a ring's own
        // edge from a corner on the cut's boundary; split onto that edge,
        // the way back is a spike.
        for (ring, _, _) in &mut rings {
            despike(ring);
        }
        for hole in holes.iter_mut().flatten() {
            despike(hole);
        }
        let mut k = 0;
        while k < rings.len() {
            if rings[k].0.len() < 3 {
                rings.remove(k);
                holes.remove(k);
            } else {
                holes[k].retain(|h| h.len() >= 3);
                k += 1;
            }
        }
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

/// A loop with every spike out to a corner and back the same way taken
/// out, which a clip leaves where a cut's boundary only touches the loop.
fn despike(ring: &mut Vec<usize>) {
    let mut k = 0;
    while ring.len() >= 3 && k < ring.len() {
        let m = ring.len();
        if ring[(k + m - 1) % m] == ring[(k + 1) % m] {
            let (a, b) = (k, (k + 1) % m);
            ring.remove(a.max(b));
            ring.remove(a.min(b));
            k = k.saturating_sub(2);
        } else {
            k += 1;
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
