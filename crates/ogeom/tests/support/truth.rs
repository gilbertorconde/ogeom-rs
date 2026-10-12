//! The exact shape of a step-conformance generated part, from its
//! `truth.json`: solids of revolution and prisms.
//!
//! Ported from step-conformance's own `src/truth.rs` (MIT OR Apache-2.0).
//! A solid of revolution is a profile of segments and arcs in the
//! half-plane through an axis, turned a full turn: the nearest point of its
//! surface to `p` lies in the half-plane through `p`, so distance to it is
//! distance to the profile there. A prism is a flat profile of loops swept
//! straight along an axis: the nearest point of a wall lies over the
//! nearest point of the profile, and of a cap straight below or above `p`
//! when `p` lies over the profile.

use super::json::Json;
use ogeom::math::{Point, Vector};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

type V2 = [f64; 2];

fn sub(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn len(a: V2) -> f64 {
    a[0].hypot(a[1])
}

/// A piece of profile, in `[z, r]` for a solid of revolution and `[x, y]`
/// for a prism.
#[derive(Debug, Clone, Copy)]
enum Element {
    Segment {
        from: V2,
        to: V2,
    },
    /// `center + radius (cos t, sin t)`, `t` from `start` through `start + sweep`.
    Arc {
        center: V2,
        radius: f64,
        start: f64,
        sweep: f64,
    },
}

/// The point at angle `t` on the unit circle, exact at every quarter turn.
fn unit(t: f64) -> V2 {
    let quarters = t / FRAC_PI_2;
    if (quarters - quarters.round()).abs() < 1e-12 {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = quarters.round().rem_euclid(4.0) as usize;
        [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]][k]
    } else {
        [t.cos(), t.sin()]
    }
}

impl Element {
    fn read(json: &Json) -> Option<Self> {
        let (kind, body) = json.variant()?;
        match kind {
            "segment" => Some(Self::Segment {
                from: body.get("from")?.numbers()?,
                to: body.get("to")?.numbers()?,
            }),
            "arc" => Some(Self::Arc {
                center: body.get("center")?.numbers()?,
                radius: body.get("radius")?.number()?,
                start: body.get("start")?.number()?,
                sweep: body.get("sweep")?.number()?,
            }),
            _ => None,
        }
    }

    fn on_axis(&self) -> bool {
        matches!(self, Self::Segment { from, to } if from[1].abs() < 1e-12 && to[1].abs() < 1e-12)
    }

    fn at(&self, t: f64) -> V2 {
        match *self {
            Self::Segment { from, to } => [
                from[0] + (to[0] - from[0]) * t,
                from[1] + (to[1] - from[1]) * t,
            ],
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let u = unit(start + sweep * t);
                [center[0] + radius * u[0], center[1] + radius * u[1]]
            }
        }
    }

    fn length(&self) -> f64 {
        match *self {
            Self::Segment { from, to } => len(sub(to, from)),
            Self::Arc { radius, sweep, .. } => radius * sweep.abs(),
        }
    }

    /// A box holding the piece: an arc's whole circle's, enough to index.
    fn bounds(&self) -> (V2, V2) {
        match *self {
            Self::Segment { from, to } => (
                [from[0].min(to[0]), from[1].min(to[1])],
                [from[0].max(to[0]), from[1].max(to[1])],
            ),
            Self::Arc { center, radius, .. } => (
                [center[0] - radius, center[1] - radius],
                [center[0] + radius, center[1] + radius],
            ),
        }
    }

    fn distance(&self, q: V2) -> f64 {
        match *self {
            Self::Segment { from, to } => {
                let ab = sub(to, from);
                let aq = sub(q, from);
                let t = ((aq[0] * ab[0] + aq[1] * ab[1]) / (ab[0] * ab[0] + ab[1] * ab[1]))
                    .clamp(0.0, 1.0);
                len(sub(q, [from[0] + ab[0] * t, from[1] + ab[1] * t]))
            }
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let d = sub(q, center);
                let angle = d[1].atan2(d[0]);
                let gone = if sweep >= 0.0 {
                    (angle - start).rem_euclid(TAU)
                } else {
                    (start - angle).rem_euclid(TAU)
                };
                if gone <= sweep.abs() {
                    (len(d) - radius).abs()
                } else {
                    len(sub(q, self.at(0.0))).min(len(sub(q, self.at(1.0))))
                }
            }
        }
    }

    /// How many times the piece crosses the ray from `q` along +x.
    fn crossings(&self, q: V2) -> usize {
        match *self {
            Self::Segment { from: a, to: b } => {
                if (a[1] > q[1]) == (b[1] > q[1]) {
                    return 0;
                }
                let x = a[0] + (q[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                usize::from(x > q[0])
            }
            Self::Arc {
                center: c,
                radius,
                start,
                sweep,
            } => {
                // Cut where the circle turns, at its top and bottom, into
                // pieces that each cross a level at most once.
                let (lo, hi) = (start.min(start + sweep), start.max(start + sweep));
                let mut cuts = vec![lo];
                let mut t = ((lo - FRAC_PI_2) / PI).floor().mul_add(PI, FRAC_PI_2);
                while t < hi && cuts.len() < 4 {
                    if t > lo {
                        cuts.push(t);
                    }
                    t += PI;
                }
                cuts.push(hi);
                cuts.windows(2)
                    .filter(|w| {
                        let (ua, ub) = (unit(w[0]), unit(w[1]));
                        let (ay, by) = (c[1] + radius * ua[1], c[1] + radius * ub[1]);
                        if (ay > q[1]) == (by > q[1]) {
                            return false;
                        }
                        let side = f64::midpoint(w[0], w[1]).cos().signum();
                        let x = side.mul_add(
                            radius
                                .mul_add(radius, -(q[1] - c[1]).powi(2))
                                .max(0.0)
                                .sqrt(),
                            c[0],
                        );
                        x > q[0]
                    })
                    .count()
            }
        }
    }
}

/// The pieces of a profile by where they lie, so the nearest to a point,
/// and those a ray from it crosses, are found without trying them all.
struct Index {
    pieces: Vec<Element>,
    min: V2,
    cell: f64,
    nx: usize,
    ny: usize,
    /// For each cell, row by row, the pieces whose box meets it.
    cells: Vec<Vec<usize>>,
    /// For each row, the pieces whose box meets it.
    rows: Vec<Vec<usize>>,
}

impl Index {
    fn new(pieces: Vec<Element>) -> Self {
        let boxes: Vec<(V2, V2)> = pieces.iter().map(Element::bounds).collect();
        let mut min = [f64::INFINITY; 2];
        let mut max = [f64::NEG_INFINITY; 2];
        for (lo, hi) in &boxes {
            min = [min[0].min(lo[0]), min[1].min(lo[1])];
            max = [max[0].max(hi[0]), max[1].max(hi[1])];
        }
        let size = [(max[0] - min[0]).max(1e-9), (max[1] - min[1]).max(1e-9)];
        // About a piece to a cell, and at most 1024 cells a side.
        #[allow(clippy::cast_precision_loss)]
        let cell = (size[0] * size[1] / pieces.len().max(1) as f64)
            .sqrt()
            .max(size[0].max(size[1]) / 1024.0);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (nx, ny) = (
            ((size[0] / cell).ceil() as usize).max(1),
            ((size[1] / cell).ceil() as usize).max(1),
        );
        let mut index = Self {
            pieces: Vec::new(),
            min,
            cell,
            nx,
            ny,
            cells: vec![Vec::new(); nx * ny],
            rows: vec![Vec::new(); ny],
        };
        for (k, (lo, hi)) in boxes.iter().enumerate() {
            let ((i0, j0), (i1, j1)) = (index.cell_of(*lo), index.cell_of(*hi));
            for j in j0..=j1 {
                index.rows[j].push(k);
                for i in i0..=i1 {
                    index.cells[j * nx + i].push(k);
                }
            }
        }
        index.pieces = pieces;
        index
    }

    /// The cell holding `q`, or the one nearest it.
    fn cell_of(&self, q: V2) -> (usize, usize) {
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let at = |x: f64, lo: f64, n: usize| {
            ((x - lo) / self.cell).floor().clamp(0.0, (n - 1) as f64) as usize
        };
        (
            at(q[0], self.min[0], self.nx),
            at(q[1], self.min[1], self.ny),
        )
    }

    /// How far `q` is from the nearest piece, the cells searched in rings
    /// out from the one nearest it.
    #[allow(
        clippy::cast_possible_wrap,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    fn nearest(&self, q: V2) -> f64 {
        if self.pieces.is_empty() {
            return f64::INFINITY;
        }
        let (ci, cj) = self.cell_of(q);
        let (ci, cj) = (ci as isize, cj as isize);
        let (nx, ny) = (self.nx as isize, self.ny as isize);
        let mut best = f64::INFINITY;
        for k in 0..=nx.max(ny) {
            for j in (cj - k).max(0)..=(cj + k).min(ny - 1) {
                let edge = (j - cj).abs() == k;
                let mut i = (ci - k).max(0);
                while i <= (ci + k).min(nx - 1) {
                    for &p in &self.cells[(j * nx + i) as usize] {
                        best = best.min(self.pieces[p].distance(q));
                    }
                    i = if edge || i == ci + k { i + 1 } else { ci + k };
                }
            }
            if best <= k as f64 * self.cell {
                break;
            }
        }
        best
    }

    /// Whether `q` lies inside the profile's loops: a ray along +x from it
    /// crosses them an odd number of times.
    fn inside(&self, q: V2) -> bool {
        let row = ((q[1] - self.min[1]) / self.cell).floor();
        #[allow(clippy::cast_precision_loss)]
        if !(row >= 0.0 && row < self.ny as f64) {
            return false;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let row = &self.rows[row as usize];
        row.iter()
            .map(|&p| self.pieces[p].crossings(q))
            .sum::<usize>()
            % 2
            == 1
    }
}

/// A solid of revolution: its axis and the pieces of profile off it.
struct Solid {
    origin: Point,
    axis: Vector,
    index: Index,
}

/// A prism: its frame, height, and every piece of every loop.
struct Prism {
    origin: Point,
    axis: Vector,
    x: Vector,
    y: Vector,
    height: f64,
    index: Index,
}

/// Every solid of a generated part.
pub struct Truth {
    solids: Vec<Solid>,
    prisms: Vec<Prism>,
}

fn point(json: &Json) -> Option<Point> {
    let [x, y, z] = json.numbers()?;
    Some(Point::new(x, y, z))
}

fn unit_vector(json: &Json) -> Option<Vector> {
    let [x, y, z] = json.numbers()?;
    let v = Vector::new(x, y, z);
    Some(v / v.magnitude())
}

impl Truth {
    /// The truth a `truth.json` states; `None` where it cannot be read.
    pub fn read(text: &str) -> Option<Self> {
        let json = Json::parse(text)?;
        let mut solids = Vec::new();
        for s in json.get("solids").map_or(&[][..], Json::list) {
            let pieces: Option<Vec<Element>> =
                s.get("profile")?.list().iter().map(Element::read).collect();
            solids.push(Solid {
                origin: point(s.get("origin")?)?,
                axis: unit_vector(s.get("axis")?)?,
                index: Index::new(pieces?.into_iter().filter(|e| !e.on_axis()).collect()),
            });
        }
        let mut prisms = Vec::new();
        for p in json.get("prisms").map_or(&[][..], Json::list) {
            let axis = unit_vector(p.get("axis")?)?;
            let x = unit_vector(p.get("x")?)?;
            let mut pieces = Vec::new();
            for ring in p.get("loops")?.list() {
                for e in ring.list() {
                    pieces.push(Element::read(e)?);
                }
            }
            prisms.push(Prism {
                origin: point(p.get("origin")?)?,
                axis,
                x,
                y: axis.cross(x),
                height: p.get("height")?.number()?,
                index: Index::new(pieces),
            });
        }
        Some(Self { solids, prisms })
    }

    /// How far `p` is from the nearest surface of the part.
    pub fn distance(&self, p: Point) -> f64 {
        let mut best = f64::INFINITY;
        for s in &self.solids {
            let d = p - s.origin;
            let z = d.dot(s.axis);
            let r = (d - s.axis * z).magnitude();
            best = best.min(s.index.nearest([z, r]));
        }
        for prism in &self.prisms {
            let d = p - prism.origin;
            let q = [d.dot(prism.x), d.dot(prism.y)];
            let w = d.dot(prism.axis);
            let wall = prism.index.nearest(q).hypot(w - w.clamp(0.0, prism.height));
            let inside = prism.index.inside(q);
            best = best.min(if inside {
                wall.min(w.abs()).min((w - prism.height).abs())
            } else {
                wall
            });
        }
        best
    }

    /// Points spread over every surface: `along` down each piece of a
    /// profile and `around` round an axis, or about an `around`th of a
    /// prism's widest side apart.
    pub fn samples(&self, along: usize, around: usize) -> Vec<Point> {
        let mut out = Vec::new();
        for s in &self.solids {
            let seed = if s.axis.x.abs() < 0.9 {
                Vector::new(1.0, 0.0, 0.0)
            } else {
                Vector::new(0.0, 1.0, 0.0)
            };
            let x = seed.cross(s.axis);
            let x = x / x.magnitude();
            let y = s.axis.cross(x);
            let along = (along * 16 / s.index.pieces.len().max(16)).max(4);
            for e in &s.index.pieces {
                for i in 0..=along {
                    #[allow(clippy::cast_precision_loss)]
                    let zr = e.at(i as f64 / along as f64);
                    let turns = if zr[1] < 1e-12 { 1 } else { around };
                    for j in 0..turns {
                        #[allow(clippy::cast_precision_loss)]
                        let t = TAU * j as f64 / turns as f64;
                        out.push(s.origin + s.axis * zr[0] + (x * t.cos() + y * t.sin()) * zr[1]);
                    }
                }
            }
        }
        for prism in &self.prisms {
            let world =
                |q: V2, w: f64| prism.origin + prism.x * q[0] + prism.y * q[1] + prism.axis * w;
            let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
            for e in &prism.index.pieces {
                for i in 0..=16 {
                    let q = e.at(f64::from(i) / 16.0);
                    lo = [lo[0].min(q[0]), lo[1].min(q[1])];
                    hi = [hi[0].max(q[0]), hi[1].max(q[1])];
                }
            }
            #[allow(clippy::cast_precision_loss)]
            let spacing = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(prism.height) / around as f64;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let up = ((prism.height / spacing).ceil() as usize).max(2);
            for e in &prism.index.pieces {
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let steps = ((e.length() / spacing).ceil() as usize).max(8);
                for i in 0..=steps {
                    #[allow(clippy::cast_precision_loss)]
                    let q = e.at(i as f64 / steps as f64);
                    for j in 0..=up {
                        #[allow(clippy::cast_precision_loss)]
                        out.push(world(q, prism.height * j as f64 / up as f64));
                    }
                }
            }
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let (nx, ny) = (
                ((hi[0] - lo[0]) / spacing).ceil() as usize,
                ((hi[1] - lo[1]) / spacing).ceil() as usize,
            );
            for i in 0..nx {
                for j in 0..ny {
                    #[allow(clippy::cast_precision_loss)]
                    let q = [
                        spacing.mul_add(i as f64 + 0.5, lo[0]),
                        spacing.mul_add(j as f64 + 0.5, lo[1]),
                    ];
                    if prism.index.inside(q) {
                        out.push(world(q, 0.0));
                        out.push(world(q, prism.height));
                    }
                }
            }
        }
        out
    }
}
