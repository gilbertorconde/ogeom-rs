//! P4: one assembly of the result from its rings.
//!
//! The polyhedral classes hand over planar rings over shared points: every
//! shrunk face, every strip, every vertex facet. Here the rings are cleaned
//! (corners repeated in turn, rings of no area), a point lying inside
//! another ring's edge is put into that edge, rings of one plane sharing an
//! edge are joined (facing opposite ways, the smaller is a corner the larger
//! loses), and the solid is built once over the shared points, each face
//! keeping its own plane. Curved pieces (M2's strips and vertex patches) are
//! built by [`build_faces`], over edges and surfaces the caller gives.

use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

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
    tol: Tolerances,
}

/// The rings ready for a solid: the points they use, each ring wound out of
/// the solid, and each ring's origins.
pub(crate) struct Assembled {
    pub points: Vec<Point>,
    pub rings: Vec<Vec<usize>>,
    pub origins: Vec<Vec<Origin>>,
}

impl Rings {
    pub fn new(tol: Tolerances) -> Self {
        Self {
            points: Vec::new(),
            grid: FastMap::default(),
            cell: tol.confusion() * 16.0,
            rings: Vec::new(),
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
        self.rings.push((ring, vec![origin], facing));
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
        for (ring, origin, facing) in std::mem::take(&mut self.rings) {
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
                continue;
            }
            let n = self.newell(&clean);
            if n.magnitude() <= tol.confusion() * tol.confusion() {
                continue;
            }
            rings.push((clean, origin, facing));
        }
        self.split_on_edges(&mut rings);
        self.join_coplanar(&mut rings);
        self.drop_straight_corners(&mut rings);
        if !super::guards::rings_fit(&self.points, &rings, tol)? {
            return Ok(None);
        }
        // Only the points the rings use, in first use.
        let mut renumber: FastMap<usize, usize> = FastMap::default();
        let mut points = Vec::new();
        let mut out_rings = Vec::with_capacity(rings.len());
        let mut origins = Vec::with_capacity(rings.len());
        for (ring, origin, _) in rings {
            let ring: Vec<usize> = ring
                .iter()
                .map(|i| {
                    *renumber.entry(*i).or_insert_with(|| {
                        points.push(self.points[*i]);
                        points.len() - 1
                    })
                })
                .collect();
            out_rings.push(ring);
            origins.push(origin);
        }
        Ok(Some(Assembled {
            points,
            rings: out_rings,
            origins,
        }))
    }

    /// Every point the rings use that lies inside another ring's edge, put
    /// into that edge in order along it.
    fn split_on_edges(&self, rings: &mut [Ring]) {
        let tol = self.tol;
        let used: Vec<usize> = {
            let mut set: Vec<usize> = rings
                .iter()
                .flat_map(|(r, _, _)| r.iter().copied())
                .collect();
            set.sort_unstable();
            set.dedup();
            set
        };
        for (ring, _, _) in rings.iter_mut() {
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
            *ring = out;
        }
    }

    /// Rings of one plane that share an edge, joined into one across it,
    /// until no two do. Facing the same way they are one face; facing
    /// opposite ways the smaller is a piece the larger loses (a vertex's
    /// corner a bevel cuts off an end face), and joining takes it out. A
    /// ring joined away to nothing is dropped.
    fn join_coplanar(&self, rings: &mut Vec<Ring>) {
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
                    if let Some(ring) = splice(&rings[i].0, &rings[j].0) {
                        // A face takes in what joins it, which then makes no
                        // face of its own; two faces joined are both the one.
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
                        let gone = ring.len() < 3
                            || self.newell(&ring).magnitude() <= tol.confusion() * tol.confusion();
                        rings.remove(j);
                        if gone {
                            rings.remove(i);
                        } else {
                            rings[i] = (ring, origin, facing);
                        }
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
    fn drop_straight_corners(&self, rings: &mut [Ring]) {
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
        let mut holders: FastMap<usize, Vec<(usize, usize)>> = FastMap::default();
        for (ri, (ring, _, _)) in rings.iter().enumerate() {
            for (k, p) in ring.iter().enumerate() {
                holders.entry(*p).or_default().push((ri, k));
            }
        }
        let mut drop: FastSet<usize> = FastSet::default();
        for (p, held) in &holders {
            if held.iter().all(|(ri, k)| straight(&rings[*ri].0, *k)) {
                drop.insert(*p);
            }
        }
        for (ring, _, _) in rings.iter_mut() {
            ring.retain(|p| !drop.contains(p));
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

/// One edge of a curved sheet: its two ends (indices into the sheet's
/// points) and its curve with the range between them, or `None` for the
/// straight segment.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the curved classes build their strips and vertex blends through it"
    )
)]
#[derive(Debug, Clone)]
pub(crate) struct SheetEdge {
    pub ends: (usize, usize),
    pub curve: Option<(Curve, (f64, f64))>,
}

/// One face of a curved sheet: its surface (`None` for the plane its ring
/// spans), its ring of edges each with the way it is walked, wound
/// counter-clockwise seen from outside the solid.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the curved classes build their strips and vertex blends through it"
    )
)]
#[derive(Debug, Clone)]
pub(crate) struct SheetFace {
    pub surface: Option<SurfaceGeometry>,
    pub ring: Vec<(usize, bool)>,
}

/// The faces of a sheet built over shared vertices and edges, in the
/// sheet's order: each edge made once and shared by the two faces walking
/// it, each face on its surface with the pcurves its edges need. The
/// faces are left for the caller to assemble with the rest.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// an edge's curve does not meet its ends, a face's ring does not close, or
/// a face cannot be trimmed on its surface.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the curved classes build their strips and vertex blends through it"
    )
)]
pub(crate) fn build_faces(
    model: &mut Model,
    points: &[Point],
    edges: &[SheetEdge],
    faces: &[SheetFace],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let mut made: Vec<Shape> = Vec::with_capacity(edges.len());
    for e in edges {
        let (a, b) = e.ends;
        let (curve, range) = match &e.curve {
            Some((curve, range)) => (curve.clone(), *range),
            None => {
                let length = points[a].distance(points[b]);
                let line = ogeom_geom::LineCurve::segment(points[a], points[b], tol)?;
                (Curve::Line(line), (0.0, length))
            }
        };
        made.push(
            ogeom_algo::make_edge_between(model, curve, range, &vertices[a], &vertices[b], tol)?
                .shape,
        );
    }
    let mut out = Vec::with_capacity(faces.len());
    for face in faces {
        let ring: Vec<Shape> = face
            .ring
            .iter()
            .map(|(e, forward)| {
                if *forward {
                    made[*e].clone()
                } else {
                    made[*e].reversed()
                }
            })
            .collect();
        let surface = match &face.surface {
            Some(surface) => surface.clone(),
            None => {
                let corners: Vec<Point> = face
                    .ring
                    .iter()
                    .map(|(e, forward)| {
                        let (a, b) = edges[*e].ends;
                        points[if *forward { a } else { b }]
                    })
                    .collect();
                let mut n = Vector::ZERO;
                for k in 1..corners.len().saturating_sub(1) {
                    n += (corners[k] - corners[0]).cross(corners[k + 1] - corners[0]);
                }
                let normal = ogeom_math::Direction::new(n, tol)?;
                let x = ogeom_math::Direction::new(corners[1] - corners[0], tol)?;
                let plane =
                    ogeom_math::Plane::new(ogeom_math::Frame::new(corners[0], normal, x, tol)?);
                ogeom_geom::PlaneSurface::new(plane).into()
            }
        };
        let wire = ogeom_algo::make_wire(model, &ring, tol)?.shape;
        let made_face = ogeom_algo::make_face(model, surface, &[wire], tol)?.shape;
        ogeom_algo::trimmed_or_fitted_where_bare(model, &made_face, tol)?;
        out.push(made_face);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_math::{Circle, Cylinder, Direction, Frame};

    const T: Tolerances = Tolerances::millimetres();

    /// A 10 cube with its top front edge rounded at 2, built from a sheet:
    /// six planar faces, the cylinder strip between the contact lines and
    /// the two arcs ending it. Valid, and the cube less the round's
    /// section (`(1 - pi/4) r^2`) along its 10.
    #[test]
    fn a_sheet_with_a_cylinder_strip_builds_the_rounded_cube() {
        let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let points = vec![
            p(0.0, 0.0, 0.0),
            p(10.0, 0.0, 0.0),
            p(10.0, 10.0, 0.0),
            p(0.0, 10.0, 0.0),
            p(0.0, 0.0, 8.0),
            p(10.0, 0.0, 8.0),
            p(10.0, 2.0, 10.0),
            p(0.0, 2.0, 10.0),
            p(10.0, 10.0, 10.0),
            p(0.0, 10.0, 10.0),
        ];
        let arc = |x: f64| -> Curve {
            let frame = Frame::new(
                p(x, 2.0, 8.0),
                Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
                Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
                T,
            )
            .unwrap();
            Curve::Circle(ogeom_geom::CircleCurve::new(
                Circle::new(frame, 2.0, T).unwrap(),
            ))
        };
        let line = |a: usize, b: usize| SheetEdge {
            ends: (a, b),
            curve: None,
        };
        let edges = vec![
            line(0, 1), // 0
            line(1, 2), // 1
            line(2, 3), // 2
            line(3, 0), // 3
            line(0, 4), // 4
            line(1, 5), // 5
            line(2, 8), // 6
            line(3, 9), // 7
            line(4, 5), // 8: the front contact
            line(7, 6), // 9: the top contact
            line(6, 8), // 10
            line(8, 9), // 11
            line(9, 7), // 12
            SheetEdge {
                ends: (4, 7),
                curve: Some((arc(0.0), (0.0, std::f64::consts::FRAC_PI_2))),
            }, // 13
            SheetEdge {
                ends: (5, 6),
                curve: Some((arc(10.0), (0.0, std::f64::consts::FRAC_PI_2))),
            }, // 14
        ];
        let plane = |ring: Vec<(usize, bool)>| SheetFace {
            surface: None,
            ring,
        };
        let axis = Frame::new(
            p(0.0, 2.0, 8.0),
            Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
            Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
            T,
        )
        .unwrap();
        let cylinder =
            ogeom_geom::CylinderSurface::new(Cylinder::new(axis, 2.0, T).unwrap(), (-1.0, 11.0))
                .unwrap();
        let faces = vec![
            plane(vec![(3, false), (2, false), (1, false), (0, false)]), // bottom
            plane(vec![(9, true), (10, true), (11, true), (12, true)]),  // top
            plane(vec![(0, true), (5, true), (8, false), (4, false)]),   // front
            plane(vec![(7, true), (11, false), (6, false), (2, true)]),  // back
            plane(vec![
                (4, true),
                (13, true),
                (12, false),
                (7, false),
                (3, true),
            ]), // left
            plane(vec![
                (1, true),
                (6, true),
                (10, false),
                (14, false),
                (5, false),
            ]), // right
            SheetFace {
                surface: Some(SurfaceGeometry::Cylinder(cylinder)),
                ring: vec![(8, true), (14, true), (9, false), (13, false)],
            },
        ];
        let mut m = Model::new();
        let built = build_faces(&mut m, &points, &edges, &faces, T).unwrap();
        let shell = ogeom_algo::make_shell(&mut m, &built).unwrap().shape;
        let solid = ogeom_algo::make_solid(&mut m, &[shell]).unwrap().shape;
        let found = ogeom_algo::check(&m, &solid, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let volume = ogeom_algo::volume_properties(
            &m,
            &solid,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass;
        let want = 1000.0 - (1.0 - std::f64::consts::FRAC_PI_4) * 4.0 * 10.0;
        assert!((volume - want).abs() < 1e-6, "{volume} against {want}");
    }
}
