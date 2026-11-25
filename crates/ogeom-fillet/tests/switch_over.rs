//! Rounds whose ball changes walls off the chain's vertex (ADR 0001,
//! amendment "a ball changing walls off the vertex").
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

/// The free-ball reference for a concave foot loop round a convex wall,
/// independent of the kernel.
///
/// The solid is a plate (z ≤ `base`) under the opening `S = (E ⊕ B_r)` of
/// a convex body by balls of radius `r`, where `E` is the set of the
/// opening balls' centres. Rounding the foot loop at `r` fills what no
/// ball of radius `r` clear of the solid reaches: the free centres are
/// `F = {z ≥ base + r} ∖ int D`, `D = E ⊕ B_2r`. For a point `p` outside
/// `S` with `base ≤ z < base + r`, `dist(p, F) ≤ r` exactly where one of
/// three holds:
/// - its lift to `z = base + r` lies outside `D`;
/// - `p` moved `2r − d` along the normal from its nearest point of `E`
///   (`d` its distance to `E`, the nearest point of `∂D`) stands at or
///   above `z = base + r`;
/// - `p` lies within `r` of the curve `Γ = ∂D ∩ {z = base + r}`.
///
/// Every other candidate nearest point of `F` is either on the plane `z =
/// base + r` (so at the lift or on `Γ`) or a critical point on `∂D`
/// above it, which for `E` convex is the second case or lies farther than
/// `2r`.
mod reference {
    /// A convex body `E`: its nearest point to `p` (and the distance),
    /// `p` itself when inside.
    pub trait Body {
        fn nearest(&self, p: [f64; 3]) -> (f64, [f64; 3]);
    }

    pub struct Fill<'a, B: Body> {
        pub body: &'a B,
        pub base: f64,
        pub r: f64,
        /// A point of `E` at every height the loop's round spans, about
        /// which each section of `S` and `D` is star-shaped.
        pub centre: [f64; 2],
        gamma: Vec<[f64; 2]>,
    }

    /// Six-point Gauss-Legendre nodes and weights on `[-1, 1]`.
    const GL: [(f64, f64); 6] = [
        (-0.932_469_514_203_152, 0.171_324_492_379_170),
        (-0.661_209_386_466_265, 0.360_761_573_048_139),
        (-0.238_619_186_083_197, 0.467_913_934_572_691),
        (0.238_619_186_083_197, 0.467_913_934_572_691),
        (0.661_209_386_466_265, 0.360_761_573_048_139),
        (0.932_469_514_203_152, 0.171_324_492_379_170),
    ];

    fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    }

    fn seg_dist(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
        let (ux, uy) = (b[0] - a[0], b[1] - a[1]);
        let (vx, vy) = (p[0] - a[0], p[1] - a[1]);
        let l2 = ux * ux + uy * uy;
        let t = if l2 > 0.0 {
            ((ux * vx + uy * vy) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (vx - t * ux).hypot(vy - t * uy)
    }

    /// The first `t` in `[lo, hi]` where `inside` turns false, `inside(lo)`
    /// true and `inside(hi)` false.
    fn bisect(mut lo: f64, mut hi: f64, inside: impl Fn(f64) -> bool) -> f64 {
        for _ in 0..48 {
            let m = 0.5 * (lo + hi);
            if inside(m) {
                lo = m;
            } else {
                hi = m;
            }
        }
        0.5 * (lo + hi)
    }

    impl<'a, B: Body> Fill<'a, B> {
        pub fn new(body: &'a B, base: f64, r: f64, centre: [f64; 2], samples: u32) -> Self {
            let z = base + r;
            let gamma = (0..samples)
                .map(|i| {
                    let phi = core::f64::consts::TAU * f64::from(i) / f64::from(samples);
                    let (c, s) = (phi.cos(), phi.sin());
                    let at = |t: f64| [centre[0] + c * t, centre[1] + s * t, z];
                    let mut hi = 1.0;
                    while body.nearest(at(hi)).0 < 2.0 * r {
                        hi *= 2.0;
                    }
                    let t = bisect(0.0, hi, |t| body.nearest(at(t)).0 < 2.0 * r);
                    [centre[0] + c * t, centre[1] + s * t]
                })
                .collect();
            Self {
                body,
                base,
                r,
                centre,
                gamma,
            }
        }

        fn to_gamma(&self, p: [f64; 2]) -> f64 {
            let n = self.gamma.len();
            let stride = 16;
            let mut best = (f64::INFINITY, 0);
            for i in (0..n).step_by(stride) {
                let g = self.gamma[i];
                let d = (g[0] - p[0]).hypot(g[1] - p[1]);
                if d < best.0 {
                    best = (d, i);
                }
            }
            let mut least = f64::INFINITY;
            for k in 0..4 * stride {
                let i = (best.1 + n + k - 2 * stride) % n;
                least = least.min(seg_dist(p, self.gamma[i], self.gamma[(i + 1) % n]));
            }
            least
        }

        fn in_solid(&self, p: [f64; 3]) -> bool {
            p[2] <= self.base || self.body.nearest(p).0 <= self.r
        }

        /// Whether `p`, outside the solid with `base ≤ z < base + r`, is
        /// filled.
        fn filled(&self, p: [f64; 3]) -> bool {
            let (r, top) = (self.r, self.base + self.r);
            if self.body.nearest([p[0], p[1], top]).0 >= 2.0 * r {
                return false;
            }
            let (d, e) = self.body.nearest(p);
            if d > 0.0 && p[2] + (2.0 * r - d) * (p[2] - e[2]) / d >= top {
                return false;
            }
            let h = top - p[2];
            self.to_gamma([p[0], p[1]]) > (r * r - h * h).max(0.0).sqrt()
        }

        /// The fill's area in the section at height `z` between `phi`
        /// and `phi + dphi`, read on the ray at `phi` (per unit angle).
        fn ray(&self, phi: f64, z: f64) -> f64 {
            let (c, s) = (phi.cos(), phi.sin());
            let at = |t: f64| [self.centre[0] + c * t, self.centre[1] + s * t, z];
            let mut t = 0.0;
            while self.in_solid(at(t)) {
                t += 0.05;
            }
            let wall = bisect(t - 0.05, t, |t| self.in_solid(at(t)));
            let eps = 1e-9;
            if !self.filled(at(wall + eps)) {
                return 0.0;
            }
            let end = bisect(wall + eps, wall + 3.0 * self.r, |t| self.filled(at(t)));
            0.5 * (end * end - wall * wall)
        }

        /// The volume the loop's round adds: `rays` rays by `panels`
        /// six-point Gauss panels in `s`, `z = base + r s²`, on each
        /// stretch between the heights in `breaks` (where the wall's
        /// contact stands, and the integrand kinks).
        pub fn volume(&self, rays: u32, panels: u32, breaks: &[f64]) -> f64 {
            let mut cuts = vec![0.0];
            cuts.extend(breaks.iter().map(|z| ((z - self.base) / self.r).sqrt()));
            cuts.push(1.0);
            let dphi = core::f64::consts::TAU / f64::from(rays);
            let mut sum = 0.0;
            for w in cuts.windows(2) {
                let h = (w[1] - w[0]) / f64::from(panels);
                for k in 0..panels {
                    let mid = w[0] + h * (f64::from(k) + 0.5);
                    for (x, g) in GL {
                        let s = mid + 0.5 * h * x;
                        let z = self.base + self.r * s * s;
                        let jac = 2.0 * self.r * s * g * 0.5 * h;
                        let mut ring = 0.0;
                        for i in 0..rays {
                            ring += self.ray(f64::from(i) * dphi, z);
                        }
                        sum += jac * ring * dphi;
                    }
                }
            }
            sum
        }
    }

    /// The free-ball reference for a convex top loop: the solid near its
    /// top is `T = (E ⊕ B_R) ∩ {z ≤ top}` for a convex `E`, and rounding
    /// the top loop at `r < R` keeps the opening of `T` by balls of radius
    /// `r`, `C ⊕ B_r` with `C = (E ⊕ B_q) ∩ {z ≤ top − r}`, `q = R − r`
    /// (for `E` convex a ball of radius `r` lies in `E ⊕ B_R` exactly where
    /// its centre lies within `q` of `E`). A point of `T` is kept where it
    /// lies within `r` of `C`. `C` is convex: its nearest point to `p` is
    /// the nearest point of `E ⊕ B_q` where that stands at or below
    /// `top − r`, and otherwise lies on that plane, in the slice of
    /// `E ⊕ B_q` there, whose outline `Γ` is sampled.
    pub struct Opening<'a, B: Body> {
        pub body: &'a B,
        pub top: f64,
        /// `R`, the radius `E` is grown by to the solid.
        pub grown: f64,
        pub r: f64,
        /// A point inside the opening at every height the round spans,
        /// about which each section of `T` and of the opening is convex.
        pub centre: [f64; 2],
        gamma: Vec<[f64; 2]>,
    }

    impl<'a, B: Body> Opening<'a, B> {
        pub fn new(
            body: &'a B,
            top: f64,
            grown: f64,
            r: f64,
            centre: [f64; 2],
            samples: u32,
        ) -> Self {
            let (q, z) = (grown - r, top - r);
            let gamma = (0..samples)
                .map(|i| {
                    let phi = core::f64::consts::TAU * f64::from(i) / f64::from(samples);
                    let (c, s) = (phi.cos(), phi.sin());
                    let at = |t: f64| [centre[0] + c * t, centre[1] + s * t, z];
                    let mut hi = 1.0;
                    while body.nearest(at(hi)).0 <= q {
                        hi *= 2.0;
                    }
                    let t = bisect(0.0, hi, |t| body.nearest(at(t)).0 <= q);
                    [centre[0] + c * t, centre[1] + s * t]
                })
                .collect();
            Self {
                body,
                top,
                grown,
                r,
                centre,
                gamma,
            }
        }

        fn to_gamma(&self, p: [f64; 2]) -> f64 {
            let n = self.gamma.len();
            let stride = 16;
            let mut best = (f64::INFINITY, 0);
            for i in (0..n).step_by(stride) {
                let g = self.gamma[i];
                let d = (g[0] - p[0]).hypot(g[1] - p[1]);
                if d < best.0 {
                    best = (d, i);
                }
            }
            let mut least = f64::INFINITY;
            for k in 0..4 * stride {
                let i = (best.1 + n + k - 2 * stride) % n;
                least = least.min(seg_dist(p, self.gamma[i], self.gamma[(i + 1) % n]));
            }
            least
        }

        fn in_solid(&self, p: [f64; 3]) -> bool {
            p[2] <= self.top && self.body.nearest(p).0 <= self.grown
        }

        /// Whether `p`, a point of `T`, lies within `r` of `C`.
        fn kept(&self, p: [f64; 3]) -> bool {
            let (q, zc) = (self.grown - self.r, self.top - self.r);
            let (d, e) = self.body.nearest(p);
            if d <= q && p[2] <= zc {
                return true;
            }
            if d > q {
                let kz = e[2] + (p[2] - e[2]) * q / d;
                if kz <= zc {
                    return d - q <= self.r;
                }
            }
            let inside = self.body.nearest([p[0], p[1], zc]).0 <= q;
            let across = if inside {
                0.0
            } else {
                self.to_gamma([p[0], p[1]])
            };
            (p[2] - zc).hypot(across) <= self.r
        }

        /// The area taken off in the section at height `z` between `phi`
        /// and `phi + dphi`, read on the ray at `phi` (per unit angle).
        fn ray(&self, phi: f64, z: f64) -> f64 {
            let (c, s) = (phi.cos(), phi.sin());
            let at = |t: f64| [self.centre[0] + c * t, self.centre[1] + s * t, z];
            let mut t = 0.0;
            while self.in_solid(at(t)) {
                t += 0.05;
            }
            let wall = bisect(t - 0.05, t, |t| self.in_solid(at(t)));
            let eps = 1e-9;
            if self.kept(at(wall - eps)) {
                return 0.0;
            }
            let from = bisect(0.0, wall - eps, |t| self.kept(at(t)));
            0.5 * (wall * wall - from * from)
        }

        /// The volume the top loop's round takes off: `rays` rays by
        /// `panels` six-point Gauss panels in `s`, `z = top − r s²`, on
        /// each stretch between the heights in `breaks`.
        pub fn volume(&self, rays: u32, panels: u32, breaks: &[f64]) -> f64 {
            let mut cuts = vec![0.0];
            cuts.extend(breaks.iter().map(|z| ((self.top - z) / self.r).sqrt()));
            cuts.push(1.0);
            cuts.sort_by(f64::total_cmp);
            let dphi = core::f64::consts::TAU / f64::from(rays);
            let mut sum = 0.0;
            for w in cuts.windows(2) {
                let h = (w[1] - w[0]) / f64::from(panels);
                for k in 0..panels {
                    let mid = w[0] + h * (f64::from(k) + 0.5);
                    for (x, g) in GL {
                        let s = mid + 0.5 * h * x;
                        let z = self.top - self.r * s * s;
                        let jac = 2.0 * self.r * s * g * 0.5 * h;
                        let mut ring = 0.0;
                        for i in 0..rays {
                            ring += self.ray(f64::from(i) * dphi, z);
                        }
                        sum += jac * ring * dphi;
                    }
                }
            }
            sum
        }
    }

    /// The free-ball reference for a concave chain where two convex bodies
    /// meet, independent of the kernel: the solid near the chain is `A ∪
    /// B`, and rounding the chain at `r` fills what no ball of radius `r`
    /// clear of both reaches. The free centres are `F = {d_A ≥ r} ∩ {d_B ≥
    /// r}`, and a point outside the solid is filled where it stands
    /// farther than `r` from `F`. The nearest point of `F` lies on `∂(A ⊕
    /// B_r)` beyond `B ⊕ B_r` (for a point outside `A`, the one point of
    /// that boundary nearest it), on `∂(B ⊕ B_r)` beyond `A ⊕ B_r` the
    /// same way, or on the curve `Γ` where the two boundaries meet. `A` is
    /// the drum `ρ ≤ rho` about the upright axis through `axis`, and `Γ`
    /// is sampled on the cylinder `ρ = rho + r`.
    pub struct Closing<'a, B: Body> {
        pub axis: [f64; 2],
        pub rho: f64,
        pub other: &'a B,
        pub r: f64,
        gamma: Vec<[f64; 3]>,
    }

    impl<'a, B: Body> Closing<'a, B> {
        /// `Γ` read on rays in the cylinder's chart `(arc length, z)` from
        /// `from`, at `samples` angles over `sweep`.
        pub fn new(
            axis: [f64; 2],
            rho: f64,
            other: &'a B,
            r: f64,
            from: [f64; 2],
            sweep: (f64, f64),
            samples: u32,
        ) -> Self {
            let big = rho + r;
            let on = |s: f64, z: f64| {
                let phi = s / big;
                [axis[0] + big * phi.cos(), axis[1] + big * phi.sin(), z]
            };
            let gamma = (0..=samples)
                .map(|i| {
                    let a = sweep.0 + (sweep.1 - sweep.0) * f64::from(i) / f64::from(samples);
                    let (c, s) = (a.cos(), a.sin());
                    let at = |t: f64| on(from[0] + c * t, from[1] + s * t);
                    let mut hi = 1.0;
                    while other.nearest(at(hi)).0 < r {
                        hi *= 2.0;
                    }
                    let t = bisect(0.0, hi, |t| other.nearest(at(t)).0 < r);
                    at(t)
                })
                .collect();
            Self {
                axis,
                rho,
                other,
                r,
                gamma,
            }
        }

        fn rho_of(&self, p: [f64; 3]) -> f64 {
            (p[0] - self.axis[0]).hypot(p[1] - self.axis[1])
        }

        fn to_gamma(&self, p: [f64; 3]) -> f64 {
            let n = self.gamma.len();
            let stride = 64;
            let mut best = (f64::INFINITY, 0);
            for i in (0..n).step_by(stride) {
                let d = dist(self.gamma[i], p);
                if d < best.0 {
                    best = (d, i);
                }
            }
            let mut least = f64::INFINITY;
            let lo = best.1.saturating_sub(2 * stride);
            let hi = (best.1 + 2 * stride).min(n - 1);
            for i in lo..hi {
                let (a, b) = (self.gamma[i], self.gamma[i + 1]);
                let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let v = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
                let l2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
                let t = if l2 > 0.0 {
                    ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / l2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                least = least.min(dist(p, [a[0] + t * u[0], a[1] + t * u[1], a[2] + t * u[2]]));
            }
            least
        }

        pub fn in_solid(&self, p: [f64; 3]) -> bool {
            self.rho_of(p) <= self.rho || self.other.nearest(p).0 <= 0.0
        }

        /// Whether `p`, outside the solid, is filled.
        pub fn filled(&self, p: [f64; 3]) -> bool {
            let r = self.r;
            let rho = self.rho_of(p);
            let (d_a, d_b) = (rho - self.rho, self.other.nearest(p));
            if d_a >= r && d_b.0 >= r {
                return false;
            }
            // On the drum's offset, radially.
            let k = (self.rho + r) / rho;
            let c = [
                self.axis[0] + (p[0] - self.axis[0]) * k,
                self.axis[1] + (p[1] - self.axis[1]) * k,
                p[2],
            ];
            if self.other.nearest(c).0 >= r && (rho - self.rho - r).abs() <= r {
                return false;
            }
            // On the other body's offset, along its nearest point.
            let (d, e) = d_b;
            if d > 0.0 {
                let c = [
                    e[0] + (p[0] - e[0]) * r / d,
                    e[1] + (p[1] - e[1]) * r / d,
                    e[2] + (p[2] - e[2]) * r / d,
                ];
                if self.rho_of(c) - self.rho >= r && (d - r).abs() <= r {
                    return false;
                }
            }
            self.to_gamma(p) > r
        }

        /// The fill on the ray from the axis at angle `theta` and height
        /// `z` (per unit angle and height): `½(ρ_e² − ρ_s²)` from where it
        /// leaves the solid to where the fill ends.
        pub fn ray(&self, theta: f64, z: f64) -> f64 {
            let (c, s) = (theta.cos(), theta.sin());
            let at = |t: f64| [self.axis[0] + c * t, self.axis[1] + s * t, z];
            let reach = self.rho + 4.0 * self.r;
            let start = if self.in_solid(at(self.rho + 1e-12)) {
                if self.in_solid(at(reach)) {
                    return 0.0;
                }
                bisect(self.rho, reach, |t| self.in_solid(at(t)))
            } else {
                self.rho
            };
            let eps = 1e-9;
            if !self.filled(at(start + eps)) {
                return 0.0;
            }
            let end = bisect(start + eps, reach, |t| self.filled(at(t)));
            0.5 * (end * end - start * start)
        }

        /// The fill over angles between the breaks `thetas` and heights
        /// from `z0` to `z1`, each integral adaptive to `tol` (absolute,
        /// shared out over its stretches): at each angle the heights are
        /// split at `z_breaks` and where the fill starts and ends, the last
        /// stretch read in `z = top − h s²` where the fill ends below `z1`.
        pub fn volume(
            &self,
            thetas: &[f64],
            (z0, z1): (f64, f64),
            z_breaks: &[f64],
            tol: f64,
        ) -> f64 {
            let width = thetas[thetas.len() - 1] - thetas[0];
            let column = |theta: f64| -> f64 {
                // Where the fill starts and ends in height, read on a grid
                // and bisected: it fills one stretch.
                let span = 4.0;
                let grid = 48;
                let step = span / f64::from(grid);
                let mut first = None;
                let mut last = None;
                for i in 0..=grid {
                    let z = z0 + step * f64::from(i);
                    if z > z1 + step {
                        break;
                    }
                    if self.ray(theta, z) > 0.0 {
                        first.get_or_insert(z);
                        last = Some(z);
                    }
                }
                let (Some(z_first), Some(z_last)) = (first, last) else {
                    return 0.0;
                };
                let edge = |mut inside: f64, mut outside: f64| {
                    for _ in 0..48 {
                        let m = 0.5 * (inside + outside);
                        if self.ray(theta, m) > 0.0 {
                            inside = m;
                        } else {
                            outside = m;
                        }
                    }
                    0.5 * (inside + outside)
                };
                let end = edge(z_last, z_last + step);
                let top = end.min(z1);
                let bottom = if z_first > z0 {
                    edge(z_first, z_first - step)
                } else {
                    z0
                };
                if top <= bottom {
                    return 0.0;
                }
                let mut cuts = vec![bottom];
                cuts.extend(z_breaks.iter().copied().filter(|z| *z > bottom && *z < top));
                cuts.push(top);
                let mut sum = 0.0;
                let share = tol / width / f64::from(u32::try_from(cuts.len()).unwrap());
                for (k, w) in cuts.windows(2).enumerate() {
                    let h = w[1] - w[0];
                    sum += if k + 2 == cuts.len() && end < z1 {
                        // z = top − h s²: dz = 2 h s ds.
                        adaptive(
                            &|s| 2.0 * h * s * self.ray(theta, w[1] - h * s * s),
                            (0.0, 1.0),
                            share,
                            0,
                        )
                    } else {
                        adaptive(&|z| self.ray(theta, z), (w[0], w[1]), share, 0)
                    };
                }
                sum
            };
            let share = tol / f64::from(u32::try_from(thetas.len() - 1).unwrap());
            thetas
                .windows(2)
                .map(|w| adaptive(&column, (w[0], w[1]), share, 0))
                .sum()
        }
    }

    fn gauss(f: &dyn Fn(f64) -> f64, (a, b): (f64, f64)) -> f64 {
        let (m, h) = (0.5 * (a + b), 0.5 * (b - a));
        GL.iter().map(|(x, g)| g * h * f(m + h * x)).sum()
    }

    /// `f` over `[a, b]` by six-point Gauss, halved until the halves agree
    /// with the whole within `tol`, twelve halvings at most.
    fn adaptive(f: &dyn Fn(f64) -> f64, (a, b): (f64, f64), tol: f64, depth: u32) -> f64 {
        let m = 0.5 * (a + b);
        let whole = gauss(f, (a, b));
        let halves = gauss(f, (a, m)) + gauss(f, (m, b));
        if (whole - halves).abs() <= tol || depth >= 10 {
            return halves;
        }
        adaptive(f, (a, m), 0.5 * tol, depth + 1) + adaptive(f, (m, b), 0.5 * tol, depth + 1)
    }

    /// A rib `x0 ≤ x ≤ x1`, `|y − cy| ≤ half`, `z ≤ top`, its two top edges
    /// rounded `round`, run on below: the product of an interval and the
    /// section `{|y − cy| ≤ half − round, z ≤ top − round} ⊕ B_round`.
    pub struct Rib {
        pub x: [f64; 2],
        pub cy: f64,
        pub half: f64,
        pub top: f64,
        pub round: f64,
    }

    impl Body for Rib {
        fn nearest(&self, p: [f64; 3]) -> (f64, [f64; 3]) {
            let x = p[0].clamp(self.x[0], self.x[1]);
            let core = self.half - self.round;
            let y0 = p[1].clamp(self.cy - core, self.cy + core);
            let z0 = p[2].min(self.top - self.round);
            let (dy, dz) = (p[1] - y0, p[2] - z0);
            let d = dy.hypot(dz);
            let (y, z) = if d <= self.round {
                (p[1], p[2])
            } else {
                (y0 + dy * self.round / d, z0 + dz * self.round / d)
            };
            let q = [x, y, z];
            (dist(p, q), q)
        }
    }

    /// The half frustum's eroded body: `{y ≥ 10 + r, ρ ≤ R(z) − r / cos α}`
    /// about the axis through (10, 10), `R(z) = r0 − (z − 4) slope`.
    pub struct Frustum {
        pub r0: f64,
        pub slope: f64,
        pub r: f64,
    }

    impl Frustum {
        fn cos_a(&self) -> f64 {
            1.0 / self.slope.hypot(1.0)
        }

        fn eroded(&self, z: f64) -> f64 {
            self.r0 - (z - 4.0) * self.slope - self.r / self.cos_a()
        }

        fn inside(&self, p: [f64; 3]) -> bool {
            p[1] >= 10.0 + self.r - 1e-12
                && (p[0] - 10.0).hypot(p[1] - 10.0) <= self.eroded(p[2]) + 1e-12
        }
    }

    /// A half cylinder's eroded body: `{x0 ≤ x ≤ x1, (y − cy)² + (z −
    /// cz)² ≤ ρ²}`, the product of an interval and a disc, so its nearest
    /// point is each factor's.
    pub struct Drum {
        pub x: [f64; 2],
        pub axis: [f64; 2],
        pub rho: f64,
    }

    impl Body for Drum {
        fn nearest(&self, p: [f64; 3]) -> (f64, [f64; 3]) {
            let x = p[0].clamp(self.x[0], self.x[1]);
            let (dy, dz) = (p[1] - self.axis[0], p[2] - self.axis[1]);
            let d = dy.hypot(dz);
            let (y, z) = if d <= self.rho {
                (p[1], p[2])
            } else {
                (
                    self.axis[0] + dy * self.rho / d,
                    self.axis[1] + dz * self.rho / d,
                )
            };
            let q = [x, y, z];
            (dist(p, q), q)
        }
    }

    impl Body for Frustum {
        fn nearest(&self, p: [f64; 3]) -> (f64, [f64; 3]) {
            if self.inside(p) {
                return (0.0, p);
            }
            let (ca, sa) = (self.cos_a(), self.slope * self.cos_a());
            let y0 = 10.0 + self.r;
            let mut best = (f64::INFINITY, p);
            if p[1] < y0 && self.inside([p[0], y0, p[2]]) {
                best = (y0 - p[1], [p[0], y0, p[2]]);
            }
            let (dx, dy) = (p[0] - 10.0, p[1] - 10.0);
            let rho = dx.hypot(dy);
            let d = (rho - self.eroded(p[2])) * ca;
            if d > 0.0 && rho > 0.0 {
                let (rf, zf) = (rho - d * ca, p[2] - d * sa);
                if rf >= 0.0 && 10.0 + dy / rho * rf >= y0 - 1e-12 && d < best.0 {
                    best = (d, [10.0 + dx / rho * rf, 10.0 + dy / rho * rf, zf]);
                }
            }
            for side in [-1.0, 1.0] {
                let on = |z: f64| -> Option<[f64; 3]> {
                    let e = self.eroded(z);
                    (e >= self.r).then(|| [10.0 + side * (e * e - self.r * self.r).sqrt(), y0, z])
                };
                let far = |z: f64| on(z).map_or(f64::INFINITY, |q| dist(p, q));
                let (mut a, mut b) = (p[2] - 4.0, p[2] + 4.0);
                let (mut at, mut least) = (a, f64::INFINITY);
                for i in 0..=64 {
                    let z = a + (b - a) * f64::from(i) / 64.0;
                    if far(z) < least {
                        (at, least) = (z, far(z));
                    }
                }
                let step = (b - a) / 64.0;
                (a, b) = (at - step, at + step);
                let g = (5.0_f64.sqrt() - 1.0) / 2.0;
                for _ in 0..80 {
                    let (c, d) = (b - g * (b - a), a + g * (b - a));
                    if far(c) < far(d) {
                        b = d;
                    } else {
                        a = c;
                    }
                }
                let z = 0.5 * (a + b);
                if let Some(q) = on(z)
                    && dist(p, q) < best.0
                {
                    best = (dist(p, q), q);
                }
            }
            best
        }
    }
}

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{FilletCorner, fillet_edges, fillet_edges_corner};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry, Transformable as _};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The D of radius `big` on its chord with its two corners rounded `r`:
/// its perimeter.
fn rounded_d_perimeter(big: f64, r: f64) -> f64 {
    let x = ((big - r).powi(2) - r * r).sqrt();
    let phi = (r / x).atan();
    2.0 * x + 2.0 * r * (PI / 2.0 + phi) + big * (PI - 2.0 * phi)
}

/// The reference on an upright D (no slope) gives the closed form: the
/// round of radius `r` round a convex outline of perimeter `p` adds
/// `p r²(1 − π/4) + π r³(5/3 − π/2)`.
#[test]
fn the_reference_rounds_an_upright_d_as_the_closed_form() {
    let body = reference::Frustum {
        r0: 5.0,
        slope: 0.0,
        r: 1.0,
    };
    let fill = reference::Fill::new(&body, 4.0, 1.0, [10.0, 12.0], 16384);
    let got = fill.volume(512, 2, &[]);
    let p = rounded_d_perimeter(5.0, 1.0);
    let want = p * (1.0 - PI / 4.0) + PI * (5.0 / 3.0 - PI / 2.0);
    eprintln!(
        "got {got}, want {want}, error {:.2e}",
        (got - want).abs() / want
    );
    assert!((got - want).abs() <= 1e-5 * want, "{got} against {want}");
}

/// The opening reference on an upright D (no slope) gives the closed
/// form: the round of radius `r` along a convex outline of perimeter `p`
/// takes off `p r²(1 − π/4) − π r³(5/3 − π/2)`.
#[test]
fn the_opening_reference_rounds_an_upright_d_s_top_as_the_closed_form() {
    let body = reference::Frustum {
        r0: 5.0,
        slope: 0.0,
        r: 1.0,
    };
    let fill = reference::Opening::new(&body, 10.0, 1.0, 0.5, [10.0, 11.5], 16384);
    let got = fill.volume(1024, 2, &[]);
    let p = rounded_d_perimeter(5.0, 1.0);
    let r: f64 = 0.5;
    let want = p * r * r * (1.0 - PI / 4.0) - PI * r.powi(3) * (5.0 / 3.0 - PI / 2.0);
    eprintln!(
        "got {got}, want {want}, error {:.2e}",
        (got - want).abs() / want
    );
    assert!((got - want).abs() <= 1e-5 * want, "{got} against {want}");
}

/// The half frustum's foot loop rounded 1 once its uprights are: the
/// reference's fill at 2048 rays by 4 panels a stretch; 1024 rays agree
/// within 3e-8, and 8 panels within 3e-9.
const FRUSTUM_FILL: f64 = 3.438_297_00;

/// Where the ball on the cone's foot touches the cone: the integrand's
/// kink in height.
fn frustum_contact() -> f64 {
    5.0 - (1.0 / 3.0) / (1.0_f64 / 9.0 + 1.0).sqrt()
}

fn frustum_body() -> reference::Frustum {
    reference::Frustum {
        r0: 5.0,
        slope: 1.0 / 3.0,
        r: 1.0,
    }
}

/// [`FRUSTUM_FILL`] is the reference's, converged.
#[test]
#[ignore = "heavy: the reference at full resolution takes over two minutes"]
fn the_reference_gives_the_frustum_fill_the_tests_read() {
    let body = frustum_body();
    let fill = reference::Fill::new(&body, 4.0, 1.0, [10.0, 12.0], 16384);
    let fine = fill.volume(2048, 4, &[frustum_contact()]);
    assert!(
        (fine - FRUSTUM_FILL).abs() <= 1e-8,
        "{fine} against {FRUSTUM_FILL}"
    );
    let coarse = fill.volume(1024, 4, &[frustum_contact()]);
    assert!(
        (coarse - FRUSTUM_FILL).abs() <= 1e-7,
        "{coarse} against {FRUSTUM_FILL}"
    );
}

fn frame_at(p: Point, z: Direction) -> Frame {
    let x = if z == Direction::Z {
        Direction::X
    } else {
        Direction::Z
    };
    Frame::new(p, z, x, T).unwrap()
}

fn only_solid(m: &Model, s: &Shape) -> Shape {
    let solids = explore_unique(m, s, ShapeType::Solid).unwrap();
    assert_eq!(solids.len(), 1);
    solids[0].clone()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The edge of `s` passing nearest `q`, 32 samples along each.
fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=32)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 32.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// A plate 20 × 20 × 4 with a frustum (radii 5 to 3, 6 high) on it, cut
/// through its axis and kept on the side y ≥ 10.
fn half_frustum(m: &mut Model) -> Shape {
    let plate = ogeom_algo::make_box(
        m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (20.0, 20.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let cone = ogeom_algo::make_cone(
        m,
        frame_at(Point::new(10.0, 10.0, 4.0), Direction::Z),
        5.0,
        3.0,
        6.0,
        T,
    )
    .unwrap()
    .shape;
    let cutter = ogeom_algo::make_box(
        m,
        frame_at(Point::new(4.0, 4.0, 3.0), Direction::Z),
        (12.0, 6.0, 8.0),
        T,
    )
    .unwrap()
    .shape;
    let half = ogeom_bool::cut(m, &cone, &cutter, T).unwrap().shape;
    let half = only_solid(m, &half);
    let fused = ogeom_bool::fuse(m, &plate, &half, T).unwrap().shape;
    only_solid(m, &fused)
}

fn frustum_uprights(m: &Model, s: &Shape) -> Vec<Shape> {
    [Point::new(6.0, 10.0, 7.0), Point::new(14.0, 10.0, 7.0)]
        .iter()
        .map(|&p| edge_near(m, s, p))
        .collect()
}

/// The edges of `s` lying on the plane `z = top` away from its outline
/// `[0, size]²`: a boss's foot loop.
fn foot_loop(m: &Model, s: &Shape, top: f64, size: f64) -> Vec<Shape> {
    let mut foot = Vec::new();
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        let on = (0..=8).all(|i| {
            let p = c.point_at(a + (b - a) * f64::from(i) / 8.0, T).unwrap();
            (p.z - top).abs() < 1e-6
                && p.x > 1.0
                && p.x < size - 1.0
                && p.y > 1.0
                && p.y < size - 1.0
        });
        if on {
            foot.push(e);
        }
    }
    foot
}

fn face_surface(m: &Model, f: &Shape) -> SurfaceGeometry {
    m.node(f)
        .and_then(|n| n.data().as_face())
        .and_then(|d| m.geometry().surface(d.surface))
        .unwrap()
        .clone()
        .transformed(&f.transform(m.datums()).unwrap(), T)
        .unwrap()
}

/// The widest turn between the tangent planes of `faces` and those of the
/// faces beside them across each of their edges, read at 15 points along
/// each.
fn widest_turn(m: &Model, s: &Shape, faces: &[Shape]) -> f64 {
    widest_turn_where(m, s, faces, &|_| true)
}

/// [`widest_turn`] read only at the points `keep` keeps.
fn widest_turn_where(m: &Model, s: &Shape, faces: &[Shape], keep: &dyn Fn(Point) -> bool) -> f64 {
    let all = explore_unique(m, s, ShapeType::Face).unwrap();
    let normal = |f: &Shape, p: Point| {
        let surface = face_surface(m, f);
        let (u, v) = ogeom_algo::project_on_surface(&surface, p, 16, T)
            .unwrap()
            .parameters;
        let n = surface.normal_at(u, v, T).unwrap().vector();
        if f.orientation() == ogeom_topo::Orientation::Reversed {
            -n
        } else {
            n
        }
    };
    let mut worst = 0.0_f64;
    for face in faces {
        for edge in explore_unique(m, face, ShapeType::Edge).unwrap() {
            let Some(other) = all.iter().find(|f| {
                !f.is_same(face)
                    && explore_unique(m, f, ShapeType::Edge)
                        .unwrap()
                        .iter()
                        .any(|e| e.is_same(&edge))
            }) else {
                continue;
            };
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, &edge, T).unwrap() else {
                continue;
            };
            for i in 1..16 {
                let p = c.point_at(a + (b - a) * f64::from(i) / 16.0, T).unwrap();
                if !keep(p) {
                    continue;
                }
                // Unsigned: a torus face's normal is read in its chart's
                // own sense, and `check` holds the faces' orientation.
                let turn = normal(face, p).dot(normal(other, p)).abs().min(1.0).acos();
                worst = worst.max(turn);
            }
        }
    }
    worst
}

/// The widest turn a join of a pipe corner's rounds may make: the rounds,
/// the bands they run on and the faces those were built against.
const G1: f64 = 1e-4;

/// The half frustum's uprights rounded 1: each band, marched between the
/// cut's plane and the cone and fitted, meets both within [`G1`]. The
/// pipes of the foot and the top loop are solved off those bands, so their
/// joins hold no better than the bands' own.
#[test]
fn a_half_frustum_s_upright_bands_join_their_walls_g1() {
    let mut m = Model::new();
    let s = half_frustum(&mut m);
    let uprights = frustum_uprights(&m, &s);
    let up = fillet_edges(&mut m, &s, &uprights, 1.0, T).unwrap();
    let found = ogeom_algo::check(&m, &up.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let made: Vec<Shape> = uprights
        .iter()
        .flat_map(|e| up.history.generated(e).to_vec())
        .collect();
    assert_eq!(made.len(), 2);
    // Each band runs past the plate and the top and is clipped there: its
    // edges on them are creases, not joins.
    let turn = widest_turn_where(&m, &up.shape, &made, &|p| {
        p.z > 4.0 + 1e-6 && p.z < 10.0 - 1e-6
    });
    eprintln!("widest turn across an upright band's join {turn:.2e} rad");
    assert!(turn <= G1, "a crease of {turn} rad");
}

/// The half frustum with its uprights rounded 1 (each band leaning into
/// the plate, carried past it and clipped there), then its foot loop
/// rounded 1: the plate's line, the two bands' feet and the cone's arc.
/// Where the ball along the line or the arc reaches a band, its contact
/// crosses the band's rail at its own height, off the vertex where the
/// creases meet (0.35 from it on the line, 0.09 on the arc): between
/// those crossings it touches the plate and the band alone, and the round
/// there is the pipe about its centre's curve. The fill is the reference's
/// within 1e-5, sound, the same in either order, G1 across every join.
#[test]
fn a_half_frustum_s_foot_rounds_through_its_pipe_corners() {
    let mut m = Model::new();
    let s = half_frustum(&mut m);
    let uprights = frustum_uprights(&m, &s);
    let up = fillet_edges(&mut m, &s, &uprights, 1.0, T).unwrap();
    let before = volume(&m, &up.shape);
    let foot = foot_loop(&m, &up.shape, 4.0, 20.0);
    assert_eq!(foot.len(), 4);
    let mut first: Option<f64> = None;
    for order in [foot.clone(), foot.iter().rev().cloned().collect()] {
        let mut rounded = m.clone();
        let built = fillet_edges(&mut rounded, &up.shape, &order, 1.0, T).unwrap();
        let found = ogeom_algo::check(&rounded, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let fill = volume(&rounded, &built.shape) - before;
        let error = (fill - FRUSTUM_FILL).abs() / FRUSTUM_FILL;
        eprintln!("fill {fill}, reference {FRUSTUM_FILL}, relative error {error:.2e}");
        assert!(error <= 1e-5, "fill {fill} against {FRUSTUM_FILL}");
        match first {
            None => first = Some(fill),
            Some(v) => assert!((fill - v).abs() <= 1e-9 * v, "{fill} against {v}"),
        }
        let made: Vec<Shape> = foot
            .iter()
            .flat_map(|e| built.history.generated(e).to_vec())
            .collect();
        assert!(!made.is_empty());
        let turn = widest_turn(&rounded, &built.shape, &made);
        eprintln!("widest turn across a join {turn:.2e} rad");
        assert!(turn <= G1, "a crease of {turn} rad");
    }
    // In one request under `Rolling` the uprights go first and the foot
    // rolls round them: the same solid in either order.
    let after = volume(&m, &up.shape) + first.unwrap();
    let mut asked = vec![
        edge_near(&m, &s, Point::new(10.0, 10.0, 4.0)),
        edge_near(&m, &s, Point::new(10.0, 15.0, 4.0)),
    ];
    asked.extend(uprights);
    for order in [asked.clone(), asked.iter().rev().cloned().collect()] {
        let mut rolled = m.clone();
        let built =
            fillet_edges_corner(&mut rolled, &s, &order, 1.0, FilletCorner::Rolling, T).unwrap();
        let found = ogeom_algo::check(&rolled, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let got = volume(&rolled, &built.shape);
        assert!((got - after).abs() <= 1e-9 * after, "{got} against {after}");
    }
}

/// `edges` of `s` rounded `r`: refused, naming a vertex where the ball
/// changes walls off its section by a distance reading `by`.
fn refused_off_vertex(m: &Model, s: &Shape, edges: &[Shape], r: f64, by: &str) {
    for order in [edges.to_vec(), edges.iter().rev().cloned().collect()] {
        let got = fillet_edges(&mut m.clone(), s, &order, r, T)
            .unwrap_err()
            .to_string();
        assert!(
            got.contains("the ball changes walls off the vertex at (")
                && got.contains(&format!("seam between the walls {by}")),
            "{got}"
        );
    }
}

/// The half frustum's top loop rounded 0.5 once its uprights are rounded
/// 1: the opening reference's volume taken off at 2048 rays by 8 panels a
/// stretch; 1024 rays agree within 1e-8, and 4 panels within 2e-7.
const FRUSTUM_TOP_TAKEN: f64 = 0.453_328_0;

/// Where the ball along the top's arc touches the cone: the integrand's
/// kink in height.
fn frustum_top_contact() -> f64 {
    9.5 + 0.5 * (1.0 / 3.0) / (1.0_f64 / 9.0 + 1.0).sqrt()
}

/// [`FRUSTUM_TOP_TAKEN`] is the opening reference's, converged.
#[test]
#[ignore = "heavy: the reference at full resolution takes minutes"]
fn the_reference_gives_the_frustum_top_the_tests_read() {
    let body = frustum_body();
    let taken = reference::Opening::new(&body, 10.0, 1.0, 0.5, [10.0, 11.5], 16384);
    let fine = taken.volume(2048, 8, &[frustum_top_contact()]);
    assert!(
        (fine - FRUSTUM_TOP_TAKEN).abs() <= 1e-7,
        "{fine} against {FRUSTUM_TOP_TAKEN}"
    );
}

/// The half frustum with its uprights rounded 1, then its top loop
/// rounded 0.5: the plate's line, the two bands' tops and the cone's arc.
/// Convex, the ball within the material along the line and along the arc
/// reaches each band's rail at its own height past the vertex where the
/// creases meet (0.19 from it on the line, 0.08 on the arc): between those
/// crossings it touches the top and the band alone, and the line's and the
/// arc's rounds run on there as pipes about the ball's centre. The cut is
/// the opening reference's within 1e-5, sound, the same in either order
/// and under `Rolling`, G1 across every join.
#[test]
fn a_half_frustum_s_top_loop_rounds_through_its_pipe_corners() {
    let mut m = Model::new();
    let s = half_frustum(&mut m);
    let uprights = frustum_uprights(&m, &s);
    let up = fillet_edges(&mut m, &s, &uprights, 1.0, T).unwrap();
    let before = volume(&m, &up.shape);
    let top = foot_loop(&m, &up.shape, 10.0, 20.0);
    assert_eq!(top.len(), 4);
    let mut first: Option<f64> = None;
    for order in [top.clone(), top.iter().rev().cloned().collect()] {
        for corner in [FilletCorner::Setback, FilletCorner::Rolling] {
            let mut rounded = m.clone();
            let built =
                fillet_edges_corner(&mut rounded, &up.shape, &order, 0.5, corner, T).unwrap();
            let found = ogeom_algo::check(&rounded, &built.shape, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            let taken = before - volume(&rounded, &built.shape);
            let error = (taken - FRUSTUM_TOP_TAKEN).abs() / FRUSTUM_TOP_TAKEN;
            eprintln!("taken {taken}, reference {FRUSTUM_TOP_TAKEN}, relative error {error:.2e}");
            assert!(error <= 1e-5, "taken {taken} against {FRUSTUM_TOP_TAKEN}");
            match first {
                None => first = Some(taken),
                Some(v) => assert!((taken - v).abs() <= 1e-9 * v, "{taken} against {v}"),
            }
            // The pipes are the line's and the arc's rounds, the bands'
            // tops marched on the fitted bands.
            let (piped, marched): (Vec<Shape>, Vec<Shape>) = top.iter().cloned().partition(|e| {
                let (c, _) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
                !matches!(c, ogeom_geom::Curve::BSpline(_))
            });
            assert_eq!(piped.len(), 2);
            let made = |edges: &[Shape]| -> Vec<Shape> {
                edges
                    .iter()
                    .flat_map(|e| built.history.generated(e).to_vec())
                    .collect()
            };
            let pipes = made(&piped);
            assert_eq!(pipes.len(), 2);
            let turn = widest_turn(&rounded, &built.shape, &pipes);
            eprintln!("widest turn across a pipe's join {turn:.2e} rad");
            assert!(turn <= G1, "a crease of {turn} rad");
            // A band marched on a fitted host holds its fit in angle as
            // well as in distance.
            let turn = widest_turn(&rounded, &built.shape, &made(&marched));
            eprintln!("widest turn across a marched top's join {turn:.2e} rad");
            assert!(turn <= G1, "a crease of {turn} rad");
        }
    }
}

/// The half cylinder's foot loop rounded 1 once its end arcs are: the
/// reference's fill at 2048 rays by 4 panels a stretch, split where the
/// cylinder's contact stands (z = 6.8); 1024 rays agree within 1e-7
/// relative, and 8 panels within 5e-8.
const DRUM_FILL: f64 = 6.315_826_7;

/// [`DRUM_FILL`] is the reference's, converged.
#[test]
#[ignore = "heavy: the reference at full resolution takes minutes"]
fn the_reference_gives_the_drum_fill_the_tests_read() {
    let body = reference::Drum {
        x: [6.0, 14.0],
        axis: [10.0, 6.0],
        rho: 3.0,
    };
    let fill = reference::Fill::new(&body, 6.0, 1.0, [10.0, 10.0], 16384);
    let fine = fill.volume(2048, 4, &[6.8]);
    assert!(
        (fine - DRUM_FILL).abs() <= 1e-7,
        "{fine} against {DRUM_FILL}"
    );
}

/// A plate 20 × 20 × 6 with a half cylinder of radius 4 lying on it, its
/// axis on the plate's top along x from 5 to 15, its four end arcs
/// rounded 1, and its foot loop: the solid and the loop's eight edges.
fn half_cylinder_with_rounded_ends(m: &mut Model) -> (Shape, Shape, Vec<Shape>) {
    let plate = ogeom_algo::make_box(
        m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (20.0, 20.0, 6.0),
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom_algo::make_cylinder(
        m,
        frame_at(Point::new(5.0, 10.0, 6.0), Direction::X),
        4.0,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom_bool::fuse(m, &plate, &drum, T).unwrap().shape;
    let s = only_solid(m, &fused);
    let arcs: Vec<Shape> = explore_unique(m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            let p = c.point_at(0.5 * (a + b), T).unwrap();
            ((p.x - 5.0).abs() < 1e-9 || (p.x - 15.0).abs() < 1e-9) && p.z > 6.5
        })
        .collect();
    assert_eq!(arcs.len(), 4);
    let up = fillet_edges(m, &s, &arcs, 1.0, T).unwrap();
    let foot = foot_loop(m, &up.shape, 6.0, 20.0);
    assert_eq!(foot.len(), 8);
    (s, up.shape, foot)
}

/// A half cylinder of radius 4 lying on a plate, its end arcs rounded 1,
/// then its foot loop rounded 1. Each arc's torus band meets the end disc
/// along a circle upright at the plate, but at the ball's height the
/// circle stands 0.17 along the disc's foot from the vertex: the ball
/// changes walls there, and the torus band's foot is rounded as the pipe
/// from that crossing to its other vertex, where the band meets the
/// cylinder along the circle in that vertex's own section. The fill is the
/// reference's within 1e-5, sound, the same in either order, G1 across
/// every join.
#[test]
fn a_half_cylinder_s_foot_rounds_through_its_pipe_corners() {
    let mut m = Model::new();
    let (_, up, foot) = half_cylinder_with_rounded_ends(&mut m);
    let before = volume(&m, &up);
    let mut first: Option<f64> = None;
    for order in [foot.clone(), foot.iter().rev().cloned().collect()] {
        let mut rounded = m.clone();
        let built = fillet_edges(&mut rounded, &up, &order, 1.0, T).unwrap();
        let found = ogeom_algo::check(&rounded, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let fill = volume(&rounded, &built.shape) - before;
        let error = (fill - DRUM_FILL).abs() / DRUM_FILL;
        eprintln!("fill {fill}, reference {DRUM_FILL}, relative error {error:.2e}");
        assert!(error <= 1e-5, "fill {fill} against {DRUM_FILL}");
        match first {
            None => first = Some(fill),
            Some(v) => assert!((fill - v).abs() <= 1e-9 * v, "{fill} against {v}"),
        }
        let made: Vec<Shape> = foot
            .iter()
            .flat_map(|e| built.history.generated(e).to_vec())
            .collect();
        assert!(!made.is_empty());
        let turn = widest_turn(&rounded, &built.shape, &made);
        eprintln!("widest turn across a join {turn:.2e} rad");
        assert!(turn <= G1, "a crease of {turn} rad");
    }
}

/// The half cylinder on its plate in one request under `Rolling`: its four
/// end arcs and its foot's four lines, rounded 1. At each mixed vertex the
/// lone round is the arc, its band a torus about the cylinder's axis: it
/// goes first, and the foot rolls round its band's foot through the
/// pipes. The same solid as the two calls, in either order.
#[test]
fn a_half_cylinder_on_a_plate_rolls_in_one_request() {
    let mut m = Model::new();
    let (s, up, foot) = half_cylinder_with_rounded_ends(&mut m);
    let mut two = m.clone();
    let built = fillet_edges(&mut two, &up, &foot, 1.0, T).unwrap();
    let want = volume(&two, &built.shape);
    let asked: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            let p = c.point_at(0.5 * (a + b), T).unwrap();
            let end = (p.x - 5.0).abs() < 1e-9 || (p.x - 15.0).abs() < 1e-9;
            let foot =
                (p.z - 6.0).abs() < 1e-9 && p.x > 1.0 && p.x < 19.0 && p.y > 1.0 && p.y < 19.0;
            (end && p.z > 6.5) || foot
        })
        .collect();
    assert_eq!(asked.len(), 8);
    for order in [asked.clone(), asked.iter().rev().cloned().collect()] {
        let mut rolled = m.clone();
        let built =
            fillet_edges_corner(&mut rolled, &s, &order, 1.0, FilletCorner::Rolling, T).unwrap();
        let found = ogeom_algo::check(&rolled, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let got = volume(&rolled, &built.shape);
        assert!((got - want).abs() <= 1e-9 * want, "{got} against {want}");
    }
}

/// A plate 30 × 20 × 4, a drum of radius 5 standing on it about (10,
/// 10), 8 high, and a rib 2 wide and 4 high fused against the drum's side
/// from x = 12 to 24, lower than the drum's top.
fn rib_under_drum(m: &mut Model) -> Shape {
    let plate = ogeom_algo::make_box(
        m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (30.0, 20.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom_algo::make_cylinder(
        m,
        frame_at(Point::new(10.0, 10.0, 4.0), Direction::Z),
        5.0,
        8.0,
        T,
    )
    .unwrap()
    .shape;
    let rib = ogeom_algo::make_box(
        m,
        frame_at(Point::new(12.0, 9.0, 4.0), Direction::Z),
        (12.0, 2.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let a = ogeom_bool::fuse(m, &plate, &drum, T).unwrap().shape;
    let a = only_solid(m, &a);
    let fused = ogeom_bool::fuse(m, &a, &rib, T).unwrap().shape;
    only_solid(m, &fused)
}

/// The rib under the drum with its two top edges rounded 0.5, and the
/// concave chain where it meets the drum: its two uprights, the two
/// bands' feet and the two arcs of its top.
fn rib_with_rounded_top(m: &mut Model) -> (Shape, Shape, Vec<Shape>) {
    let s = rib_under_drum(m);
    let tops = vec![
        edge_near(m, &s, Point::new(18.0, 11.0, 8.0)),
        edge_near(m, &s, Point::new(18.0, 9.0, 8.0)),
    ];
    let up = fillet_edges(m, &s, &tops, 0.5, T).unwrap().shape;
    let chain: Vec<Shape> = explore_unique(m, &up, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            (0..=8).all(|i| {
                let p = c.point_at(a + (b - a) * f64::from(i) / 8.0, T).unwrap();
                ((p.x - 10.0).hypot(p.y - 10.0) - 5.0).abs() < 1e-6
                    && (4.0 - 1e-6..=8.0 + 1e-6).contains(&p.z)
                    && p.x > 12.0
                    && (p.y - 10.0).abs() < 1.0 + 1e-6
            })
        })
        .collect();
    assert_eq!(chain.len(), 6);
    (s, up, chain)
}

/// The area of the 2D fill of radius `r` in the corner where a line
/// `y = a` meets the circle of radius `big` about the origin, the material
/// the disc and `y ≤ a`: by Green's theorem round its line, the ball's arc
/// and the circle's arc.
fn corner_area(big: f64, a: f64, r: f64) -> f64 {
    let corner = ((big * big - a * a).sqrt(), a);
    let c = (((big + r).powi(2) - (a + r).powi(2)).sqrt(), a + r);
    let on_line = (c.0, a);
    let on_circle = (c.0 * big / (big + r), c.1 * big / (big + r));
    let line = |p: (f64, f64), q: (f64, f64)| 0.5 * (p.0 * q.1 - q.0 * p.1);
    let arc = |centre: (f64, f64), rho: f64, from: f64, to: f64| {
        0.5 * (rho * rho * (to - from)
            + rho * (centre.0 * (to.sin() - from.sin()) - centre.1 * (to.cos() - from.cos())))
    };
    let ball_from = -core::f64::consts::FRAC_PI_2;
    let ball_to = (on_circle.1 - c.1).atan2(on_circle.0 - c.0);
    let circle_from = on_circle.1.atan2(on_circle.0);
    let circle_to = corner.1.atan2(corner.0);
    (line(corner, on_line)
        + arc(c, r, ball_from, ball_to)
        + arc((0.0, 0.0), big, circle_from, circle_to))
    .abs()
}

/// The closing reference's fill of the rib's chain rounded `r` at heights
/// `z0` to `z1`, each integral adaptive to `tol`.
fn rib_reference(r: f64, tol: f64, (z0, z1): (f64, f64)) -> f64 {
    let rib = reference::Rib {
        x: [12.0, 24.0],
        cy: 10.0,
        half: 1.0,
        top: 8.0,
        round: 0.5,
    };
    let closing = reference::Closing::new(
        [10.0, 10.0],
        5.0,
        &rib,
        r,
        [0.0, 5.0],
        (
            -core::f64::consts::FRAC_PI_4,
            PI + core::f64::consts::FRAC_PI_4,
        ),
        16384,
    );
    let wide = ((1.0 + r) / (5.0 + r)).asin();
    // Where the uprights' balls touch the rib's sides.
    let side = (1.0 / ((5.0 + r).powi(2) - (1.0 + r).powi(2)).sqrt()).atan();
    let thetas = [
        -(wide + 0.02),
        -wide,
        -0.2_f64.asin(),
        -side,
        -0.1_f64.asin(),
        0.0,
        0.1_f64.asin(),
        side,
        0.2_f64.asin(),
        wide,
        wide + 0.02,
    ];
    // Breaks: the uprights' balls' reach and their touches on the drum and
    // on the rib's sides, the rib's sides and its flat top at the drum.
    closing.volume(&thetas, (z0, z1), &[7.5, 8.0], tol)
}

/// The rib's chain fill above z = 5 by the closing reference, adaptive to
/// 1e-7: at r = 0.5, then 0.3. At 1e-6 they agree within 6e-7 relative.
const RIB_FILLS: [(f64, f64); 2] = [(0.5, 0.271_171_36), (0.3, 0.099_305_24)];

/// The closing reference reads the rib's uprights below its top, where
/// the fill is the corner's plane section, as that section's area.
#[test]
fn the_closing_reference_gives_the_rib_s_uprights_as_their_section() {
    for r in [0.5, 0.3] {
        let got = rib_reference(r, 1e-6, (5.0, 7.0));
        let want = 2.0 * 2.0 * corner_area(5.0, 1.0, r);
        assert!((got - want).abs() <= 1e-7 * want, "{got} against {want}");
    }
}

/// [`RIB_FILLS`] are the closing reference's, converged.
#[test]
#[ignore = "heavy: the reference above the uprights takes minutes"]
fn the_closing_reference_gives_the_rib_fills_the_tests_read() {
    for (r, want) in RIB_FILLS {
        let got = rib_reference(r, 1e-6, (5.0, 20.0));
        assert!((got - want).abs() <= 1e-6 * want, "{got} against {want}");
    }
}

/// A rib fused against a drum below its top, its two top edges rounded
/// 0.5, then the concave chain where it meets the drum rounded: its two
/// uprights, the two bands' feet and the two arcs of its top, all on the
/// drum, a curved cap. Where an arc's round reaches a band, its ball's
/// contact on the rib's top crosses the band's rail 0.05 (at 0.5) from the
/// vertex, over the arc; at the band's other end the band meets the
/// upright's plane along a line in that vertex's own section. Each band's
/// foot is rounded as the pipe between. The fill is the closing
/// reference's (above the plate's reach, the uprights' section below),
/// sound, the same in either order, G1 across every join. The uprights'
/// rounds are their exact prisms: a march holds its fit in distance only,
/// and left them 8e-5 short at 0.3, a crease of 3.5e-3 rad on the drum.
#[test]
fn a_rib_under_a_drum_s_top_rounds_through_its_pipe_corners() {
    let mut m = Model::new();
    let (s, up, chain) = rib_with_rounded_top(&mut m);
    // The rib's rounded top is the drum and the rounded rib, which the
    // reference reads: what the two rounds took off, by quadrature.
    let removed = volume(&m, &s) - volume(&m, &up);
    let n = 256;
    let mut model = 0.0;
    for i in 0..n {
        let t = core::f64::consts::FRAC_PI_2 * (f64::from(i) + 0.5) / f64::from(n);
        let (y, h) = (10.5 + 0.5 * t.sin(), 0.5 - 0.5 * t.cos());
        let x = 10.0 + (25.0 - (y - 10.0).powi(2)).sqrt();
        model += 2.0 * (24.0 - x) * h * 0.5 * t.cos() * core::f64::consts::FRAC_PI_2 / f64::from(n);
    }
    assert!(
        (removed - model).abs() <= 2e-5 * model,
        "{removed} against {model}"
    );
    let before = volume(&m, &up);
    for (r, reference) in RIB_FILLS {
        let mut first: Option<f64> = None;
        for order in [chain.clone(), chain.iter().rev().cloned().collect()] {
            let mut rounded = m.clone();
            let built = fillet_edges(&mut rounded, &up, &order, r, T).unwrap();
            let found = ogeom_algo::check(&rounded, &built.shape, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            // Below z = 5 the uprights' fill is their plane section, down
            // to the plate.
            let fill = volume(&rounded, &built.shape) - before - 2.0 * corner_area(5.0, 1.0, r);
            let error = (fill - reference).abs() / reference;
            eprintln!("r {r}: fill {fill}, reference {reference}, relative error {error:.2e}");
            assert!(error <= 1e-5, "fill {fill} against {reference}");
            match first {
                None => first = Some(fill),
                Some(v) => assert!((fill - v).abs() <= 1e-9 * v, "{fill} against {v}"),
            }
            let made: Vec<Shape> = chain
                .iter()
                .flat_map(|e| built.history.generated(e).to_vec())
                .collect();
            assert!(!made.is_empty());
            // The uprights' rounds end square on the plate, a face left
            // sharp, not a join.
            let turn = widest_turn_where(&rounded, &built.shape, &made, &|p| p.z > 4.0 + 1e-6);
            eprintln!("widest turn across a join {turn:.2e} rad");
            assert!(turn <= G1, "a crease of {turn} rad");
        }
    }
}

/// The rib's chain without one of its uprights: the band's foot beside it
/// hands the ball over in its vertex's section to an edge the request
/// leaves sharp, so no pipe takes that band's crossing on the arc, and the
/// request is refused by name.
#[test]
fn a_rib_s_chain_short_of_an_upright_is_refused_by_name() {
    let mut m = Model::new();
    let (_, up, chain) = rib_with_rounded_top(&mut m);
    let short: Vec<Shape> = chain
        .iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
            !((p.y - 9.0).abs() < 1e-6 && (q.y - 9.0).abs() < 1e-6)
        })
        .cloned()
        .collect();
    assert_eq!(short.len(), 5);
    refused_off_vertex(&m, &up, &short, 0.5, "0.050");
}

/// The rib under the drum in one request under `Rolling`: its two top
/// edges, its two uprights and the two arcs of its top, rounded 0.5. At
/// each mixed vertex the lone round is the top edge, on the drum: it goes
/// first, and the chain rolls round its foot through the pipes. The same
/// solid as the two calls, in either order.
#[test]
fn a_rib_under_a_drum_s_top_rolls_in_one_request() {
    let mut m = Model::new();
    let (s, up, chain) = rib_with_rounded_top(&mut m);
    let mut two = m.clone();
    let built = fillet_edges(&mut two, &up, &chain, 0.5, T).unwrap();
    let want = volume(&two, &built.shape);
    let x = 24.0_f64.sqrt() + 10.0;
    let asked: Vec<Shape> = [
        Point::new(18.0, 11.0, 8.0),
        Point::new(18.0, 9.0, 8.0),
        Point::new(x, 11.0, 6.0),
        Point::new(x, 9.0, 6.0),
        Point::new(x + 0.05, 10.5, 8.0),
        Point::new(x + 0.05, 9.5, 8.0),
    ]
    .iter()
    .map(|&p| edge_near(&m, &s, p))
    .collect();
    for order in [asked.clone(), asked.iter().rev().cloned().collect()] {
        let mut rolled = m.clone();
        let built =
            fillet_edges_corner(&mut rolled, &s, &order, 0.5, FilletCorner::Rolling, T).unwrap();
        let found = ogeom_algo::check(&rolled, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let got = volume(&rolled, &built.shape);
        assert!((got - want).abs() <= 1e-9 * want, "{got} against {want}");
    }
}
