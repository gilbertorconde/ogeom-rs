//! A concave circular rim beside a planar cap ending in a straight edge.
//!
//! The ball touches the cylindrical wall throughout. Where its planar contact
//! leaves the cap it rests on the cap's far edge instead: its centre stays on
//! the cylinder of radius `a+r` about the wall and on the cylinder of radius
//! `r` about that straight edge. Their intersection gives its height exactly.
//! The face/edge transitions have the same centre and section on both runs.

use ogeom_algo::Built;
use ogeom_algo::classify::FaceTrim;
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::march::{BlendStop, MarchedBlend, Sides};
use crate::marched::{ExactBand, Rail};
use crate::support::{Contact, edge_curve, revolved_seat};

const STATIONS: u32 = 96;

pub(crate) struct Run {
    guide: Curve,
    range: (f64, f64),
    switches: [f64; 2],
    origin: Point,
    radial_axes: [Vector; 2],
    up: Vector,
    wall_radius: f64,
    radius: f64,
    far: (Point, Vector),
    outward: Vector,
    hosts: [(SurfaceGeometry, f64); 2],
    faces: [Shape; 2],
    across: Shape,
}

impl Run {
    fn station(&self, t: f64) -> [Point; 3] {
        let radial = self.radial_axes[0] * t.cos() + self.radial_axes[1] * t.sin();
        let planar = self.origin + radial * (self.wall_radius + self.radius);
        let over = (planar - self.far.0).dot(self.outward).max(0.0);
        // rim_run proves distance(far line, origin) > wall_radius, hence
        // over < radius at every angle, including between the stations.
        let height = (self.radius * self.radius - over * over).sqrt();
        let centre = planar + self.up * height;
        let wall = self.origin + radial * self.wall_radius + self.up * height;
        let cap = planar - self.outward * over;
        [centre, wall, cap]
    }

    pub(crate) fn contacts(&self) -> Vec<Contact> {
        let mut contacts = Vec::new();
        for i in 0..=STATIONS {
            let t = station(self.range, i);
            let [_, wall, cap] = self.station(t);
            for (face, at) in [(&self.faces[0], wall), (&self.faces[1], cap)] {
                contacts.push(Contact {
                    face: face.clone(),
                    at,
                    setback: self.radius,
                });
            }
            if t >= self.switches[0] && t <= self.switches[1] {
                contacts.push(Contact {
                    face: self.across.clone(),
                    at: cap,
                    setback: self.radius,
                });
            }
        }
        contacts
    }
}

fn station((lo, hi): (f64, f64), i: u32) -> f64 {
    lo + (hi - lo) * f64::from(i) / f64::from(STATIONS)
}

/// Recognise an open circular boss-base rim whose cap contact runs off one
/// straight side in its middle, while both ends still touch the cap face.
/// Other seats continue through the ordinary fillet paths.
pub(crate) fn rim_run(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Run>> {
    let (guide, range) = edge_curve(model, edge, tol)?;
    let Curve::Circle(circle) = &guide else {
        return Ok(None);
    };
    let seat = match revolved_seat(model, solid, edge, circle, tol) {
        Ok(seat) => seat,
        Err(OgeomError::Construction(_)) => return Ok(None),
        Err(other) => return Err(other),
    };
    if seat.tau >= 0.0 || seat.sigma < 0.0 || seat.slope != 0.0 || seat.closes {
        return Ok(None);
    }
    let (Some((first, first_sign)), Some((second, second_sign))) = (
        crate::rolled::host(model, &seat.wall_face, tol)?,
        crate::rolled::host(model, &seat.cap_face, tol)?,
    ) else {
        ogeom_bail!(Dangling, "a rim's support surface is unavailable");
    };
    let cap_foot = |t: f64| -> OgeomResult<Point> {
        let radial = (guide.point_at(t, tol)? - seat.centre) / seat.radius;
        Ok(seat.centre + radial * (seat.radius + radius))
    };
    for side in explore_unique(model, &seat.cap_face, ShapeType::Edge)? {
        let Some([a, b]) = crate::consumed::segment(model, &side, tol)? else {
            continue;
        };
        let along = (b - a) / a.distance(b);
        let mut outward = along.cross(seat.up);
        if outward.dot(a - seat.centre) < 0.0 {
            outward = -outward;
        }
        let distance = (a - seat.centre).dot(outward);
        if distance <= seat.radius + tol.confusion() || distance >= seat.radius + radius {
            continue;
        }
        let over = |t: f64| -> OgeomResult<f64> { Ok((cap_foot(t)? - a).dot(outward)) };
        let middle = f64::midpoint(range.0, range.1);
        if over(range.0)? >= 0.0 || over(range.1)? >= 0.0 || over(middle)? <= 0.0 {
            continue;
        }
        let others: Vec<_> = ogeom_topo::ancestors_of(model, solid, &side, ShapeType::Face)?
            .into_iter()
            .filter(|f| !f.is_same(&seat.cap_face))
            .collect();
        let [across] = others.as_slice() else {
            continue;
        };
        let Some(flat) = crate::rolled::flat(model, across, tol)? else {
            continue;
        };
        if flat.normal.dot(outward) < 1.0 - tol.angular() {
            continue;
        }
        let root = |mut lo: f64, mut hi: f64| -> OgeomResult<f64> {
            let sign = over(lo)? > 0.0;
            for _ in 0..64 {
                let mid = f64::midpoint(lo, hi);
                if (over(mid)? > 0.0) == sign {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            Ok(f64::midpoint(lo, hi))
        };
        let switches = [root(range.0, middle)?, root(middle, range.1)?];
        let run = Run {
            guide: guide.clone(),
            range,
            switches,
            origin: seat.centre,
            radial_axes: [
                circle.circle().frame().x().vector(),
                circle.circle().frame().y().vector(),
            ],
            up: seat.up,
            wall_radius: seat.radius,
            radius,
            far: (a, along),
            outward,
            hosts: [(first.clone(), first_sign), (second.clone(), second_sign)],
            faces: [seat.wall_face.clone(), seat.cap_face.clone()],
            across: across.clone(),
        };
        // Both lines of contact must remain on the actual trimmed supports.
        // The edge contact also stays within that edge, not its infinite line.
        // FaceTrim folds periodic parameters into the actual trimming window;
        // the cylinder's principal-domain seam is not a boundary of the rim.
        let deflection = Deflection {
            chord: tol.confusion() * 10.0,
            ..Default::default()
        };
        let trims = [
            FaceTrim::of(model, &run.faces[0], deflection, tol)?,
            FaceTrim::of(model, &run.faces[1], deflection, tol)?,
        ];
        for span in [
            (range.0, switches[0]),
            (switches[0], switches[1]),
            (switches[1], range.1),
        ] {
            for i in 0..=STATIONS {
                let t = station(span, i);
                let [_, wall, cap] = run.station(t);
                for (j, point) in [(0, wall), (1, cap)] {
                    let (u, v) =
                        ogeom_algo::project_on_surface(&run.hosts[j].0, point, 32, tol)?.parameters;
                    if !trims[j]
                        .holds(model, (u, v), tol.confusion(), tol)?
                        .is_inside_or_on()
                    {
                        ogeom_bail!(
                            Construction,
                            "a rolled rim's contact {point:?} at guide {t} leaves support face {j}"
                        );
                    }
                }
                if t >= switches[0] && t <= switches[1] {
                    let along_edge = (cap - a).dot(along);
                    if along_edge < -tol.confusion() || along_edge > a.distance(b) + tol.confusion()
                    {
                        ogeom_bail!(Construction, "a rolled rim's contact leaves its far edge");
                    }
                }
            }
        }
        return Ok(Some(run));
    }
    Ok(None)
}

pub(crate) fn build(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    run: &Run,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let spans = [
        (run.range.0, run.switches[0]),
        (run.switches[0], run.switches[1]),
        (run.switches[1], run.range.1),
    ];
    let centre_at = |t| run.station(t)[0];
    let centre = |t| Ok(centre_at(t));
    let rolling = crate::between::Rolling {
        radius: run.radius,
        stations: spans
            .iter()
            .flat_map(|&span| (0..STATIONS).map(move |i| station(span, i)))
            .chain(std::iter::once(run.range.1))
            .collect(),
        centre: &centre,
        lines: vec![],
        ends: [run.range.0, run.range.1]
            .into_iter()
            .flat_map(|t| {
                let [_, a, b] = run.station(t);
                [a, b]
            })
            .collect(),
        lunes: None,
    };
    if let Some(at) = crate::between::face_inside_balls(
        model,
        solid,
        [&run.faces[0], &run.faces[1]],
        &rolling,
        tol,
    )? {
        ogeom_bail!(
            Construction,
            "a face at {at:?} stands inside a rolled rim's ball"
        );
    }
    let mut blends = Vec::new();
    for span in spans {
        let mut blend = MarchedBlend {
            spine: vec![],
            on_first: vec![],
            on_second: vec![],
            touch_first: vec![],
            touch_second: vec![],
            along: vec![],
            sides: Sides {
                first: 1,
                second: 1,
            },
            stopped: BlendStop::RanPastTheGuide,
        };
        for i in 0..=STATIONS {
            let t = station(span, i);
            let [c, wall, cap] = run.station(t);
            blend.spine.push(c);
            blend.touch_first.push(wall);
            blend.touch_second.push(cap);
            blend
                .on_first
                .push(ogeom_algo::project_on_surface(&run.hosts[0].0, wall, 32, tol)?.parameters);
            blend
                .on_second
                .push(ogeom_algo::project_on_surface(&run.hosts[1].0, cap, 32, tol)?.parameters);
            blend.along.push(t);
        }
        blends.push(blend);
    }
    let offs: Vec<Box<dyn Fn(Point) -> f64 + '_>> = spans
        .iter()
        .map(|&(lo, hi)| {
            Box::new(move |p: Point| {
                let mut best = 0;
                let mut distance = f64::INFINITY;
                for i in 0..=STATIONS {
                    let d = p.distance(centre_at(station((lo, hi), i)));
                    if d < distance {
                        distance = d;
                        best = i;
                    }
                }
                let mut a = station((lo, hi), best.saturating_sub(1));
                let mut b = station((lo, hi), (best + 1).min(STATIONS));
                for _ in 0..64 {
                    let x = a + (b - a) / 3.0;
                    let y = b - (b - a) / 3.0;
                    if p.distance(centre_at(x)) < p.distance(centre_at(y)) {
                        b = y;
                    } else {
                        a = x;
                    }
                }
                (p.distance(centre_at(f64::midpoint(a, b))) - run.radius).abs()
            }) as Box<dyn Fn(Point) -> f64>
        })
        .collect();
    let exact: Vec<_> = offs
        .iter()
        .enumerate()
        .map(|(i, off)| ExactBand {
            rails: [
                Rail::Fitted,
                if i == 1 {
                    Rail::Line(run.far.0, run.far.1)
                } else {
                    Rail::Plane(run.origin, run.up)
                },
            ],
            off: off.as_ref(),
            strips: &[],
        })
        .collect();
    crate::marched::build_chain_band_at(
        model,
        solid,
        Some(edge),
        &blends,
        &exact,
        &run.guide,
        [
            (&run.hosts[0].0, run.hosts[0].1),
            (&run.hosts[1].0, run.hosts[1].1),
        ],
        run.radius,
        tol.confusion() * 20.0,
        tol,
    )
}
