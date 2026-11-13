//! A convex straight round carried from a narrow cap onto its circular
//! concave rim. The ball keeps the outer wall and rests on the rim where
//! its ordinary contact would be inside the cylinder above the cap.

use crate::march::{BlendStop, MarchedBlend, Sides};
use crate::marched::{ExactBand, Rail};
use crate::support::{Contact, edge_curve};
use ogeom_algo::Built;
use ogeom_algo::classify::FaceTrim;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const STATIONS: u32 = 96;

pub(crate) struct Run {
    guide: Curve,
    range: (f64, f64),
    switches: [f64; 2],
    origin: Point,
    up: Vector,
    out: Vector,
    along: Vector,
    wall: f64,
    drum: f64,
    radius: f64,
    faces: [Shape; 2],
    beyond: Shape,
    rim: (Curve, (f64, f64)),
    hosts: [(SurfaceGeometry, f64); 2],
}

impl Run {
    fn station(&self, t: f64, tol: Tolerances) -> OgeomResult<[Point; 3]> {
        let y = (self.guide.point_at(t, tol)? - self.origin).dot(self.along);
        let d = self.wall - self.radius;
        let radial = self.out * d + self.along * y;
        let rho = radial.magnitude();
        let over = (self.drum - rho).max(0.0);
        let height = (self.radius * self.radius - over * over).sqrt();
        let centre = self.origin + radial - self.up * height;
        let cap = self.origin + radial * (if over > 0.0 { self.drum / rho } else { 1.0 });
        Ok([centre, centre + self.out * self.radius, cap])
    }

    pub(crate) fn contacts(&self, tol: Tolerances) -> OgeomResult<Vec<Contact>> {
        let mut found = Vec::new();
        for span in self.spans() {
            for i in 0..=STATIONS {
                let t = at(span, i);
                let [_, wall, cap] = self.station(t, tol)?;
                let p = self.guide.point_at(t, tol)?;
                for (face, point) in [(&self.faces[0], wall), (&self.faces[1], cap)] {
                    found.push(Contact {
                        face: face.clone(),
                        at: point,
                        setback: point.distance(p),
                    });
                }
                if (self.switches[0]..=self.switches[1]).contains(&t) {
                    found.push(Contact {
                        face: self.beyond.clone(),
                        at: cap,
                        setback: cap.distance(p),
                    });
                }
            }
        }
        Ok(found)
    }

    fn spans(&self) -> [(f64, f64); 3] {
        [
            (self.range.0, self.switches[0]),
            (self.switches[0], self.switches[1]),
            (self.switches[1], self.range.1),
        ]
    }
}

fn at((lo, hi): (f64, f64), i: u32) -> f64 {
    lo + (hi - lo) * f64::from(i) / f64::from(STATIONS)
}

/// Recognise the two planar supports and a concave circular rim across
/// one of them. Other transitions continue to the ordinary guards.
pub(crate) fn run(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Run>> {
    let Some(seat) = crate::consumed::seat_of(model, solid, edge, tol)? else {
        return Ok(None);
    };
    if !seat.convex || seat.normals[0].dot(seat.normals[1]).abs() > tol.angular() {
        return Ok(None);
    }
    let (guide, range) = edge_curve(model, edge, tol)?;
    for narrow in 0..2 {
        let up = seat.normals[narrow];
        let out = seat.normals[1 - narrow];
        let cap = &seat.faces[narrow];
        for far in explore_unique(model, cap, ShapeType::Edge)? {
            let (rim_curve, rim_range) = edge_curve(model, &far, tol)?;
            let Curve::Circle(circle) = &rim_curve else {
                continue;
            };
            let rim = match crate::support::revolved_seat(model, solid, &far, circle, tol) {
                Ok(rim) => rim,
                Err(ogeom_core::OgeomError::Construction(_)) => continue,
                Err(e) => return Err(e),
            };
            if rim.tau >= 0.0
                || rim.slope != 0.0
                || !rim.cap_face.is_same(cap)
                || rim.up.dot(up).abs() < 1.0 - tol.angular()
            {
                continue;
            }
            let wall = (seat.start - rim.centre).dot(out);
            let d = wall - radius;
            if wall <= rim.radius + tol.confusion()
                || d <= 0.0
                || d >= rim.radius
                || rim.radius - d >= radius
            {
                continue;
            }
            let y0 = (seat.start - rim.centre).dot(seat.along);
            let y1 = (seat.end - rim.centre).dot(seat.along);
            let switch = (rim.radius * rim.radius - d * d).sqrt();
            if y0 >= -switch || y1 <= switch {
                continue;
            }
            // In the changing section c=(d,y,-h), h²=r²-(a-rho)².
            // The crease crossing a*y/rho is strictly increasing. The
            // section can fold only at c''·(p-c)=|c'|²; p_z-c_z<=r and
            // |c'|²>=1. This upper bound on positive c_z'' proves every
            // section, including between the fitted stations, clear.
            let hmin = (radius * radius - (rim.radius - d).powi(2)).sqrt();
            let gmax = (rim.radius - d) / d;
            let bend = rim.radius * switch.powi(2) / (d.powi(3) * hmin)
                + gmax.powi(2) * switch.powi(2) / hmin.powi(3);
            if radius * bend >= 1.0 {
                ogeom_bail!(Construction, "a carried shelf's sections may cross");
            }
            let Some(first) = crate::rolled::host(model, &seat.faces[1 - narrow], tol)? else {
                continue;
            };
            let Some(second) = crate::rolled::host(model, cap, tol)? else {
                continue;
            };
            let parameter = |y| range.0 + (range.1 - range.0) * (y - y0) / (y1 - y0);
            let run = Run {
                guide: guide.clone(),
                range,
                switches: [parameter(-switch), parameter(switch)],
                origin: rim.centre,
                up,
                out,
                along: seat.along,
                wall,
                drum: rim.radius,
                radius,
                faces: [seat.faces[1 - narrow].clone(), cap.clone()],
                beyond: rim.wall_face,
                rim: (rim_curve, rim_range),
                hosts: [first, second],
            };
            let detail = ogeom_mesh::Deflection {
                chord: tol.confusion(),
                ..Default::default()
            };
            let trims = [
                FaceTrim::of(model, &run.faces[0], detail, tol)?,
                FaceTrim::of(model, &run.faces[1], detail, tol)?,
            ];
            let Some((beyond_surface, _)) = crate::rolled::host(model, &run.beyond, tol)? else {
                continue;
            };
            let beyond_trim = FaceTrim::of(model, &run.beyond, detail, tol)?;
            for contact in run.contacts(tol)? {
                // A circular rim can span several edges. Its contact must
                // lie on the trimmed cylinder as well as the planar cap.
                if contact.face.is_same(&run.beyond) {
                    let q = ogeom_algo::project_on_surface(&beyond_surface, contact.at, 32, tol)?;
                    if !beyond_trim
                        .holds(model, q.parameters, tol.confusion(), tol)?
                        .is_inside_or_on()
                    {
                        ogeom_bail!(
                            Construction,
                            "a carried shelf contact leaves its circular rim at {:?}",
                            contact.at
                        );
                    }
                    continue;
                }
                let j = usize::from(!contact.face.is_same(&run.faces[0]));
                let (u, v) = ogeom_algo::project_on_surface(&run.hosts[j].0, contact.at, 32, tol)?
                    .parameters;
                if !trims[j]
                    .holds(model, (u, v), tol.confusion(), tol)?
                    .is_inside_or_on()
                {
                    ogeom_bail!(
                        Construction,
                        "a carried shelf contact leaves its support at {:?}",
                        contact.at
                    );
                }
            }
            return Ok(Some(run));
        }
    }
    Ok(None)
}

/// Build the carried wedge on the line and supports validated by `run`.
///
/// # Panics
///
/// If a validated line guide cannot be evaluated within its own run.
#[allow(
    clippy::expect_used,
    reason = "recognition validates the line guide over the entire run"
)]
pub(crate) fn build(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    run: &Run,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let spans = run.spans();
    let stations: Vec<_> = spans
        .iter()
        .flat_map(|&span| (0..STATIONS).map(move |i| at(span, i)))
        .chain(std::iter::once(run.range.1))
        .collect();
    let section = |t| run.station(t, tol);
    let lunes = crate::between::TurningLunes {
        radius: run.radius,
        stations: &stations,
        section: &section,
    };
    let search = |_: &Model, _: &Shape, measure: &ogeom_algo::PointDistance, _: Option<Point>| {
        lunes.face_inside(measure, tol)
    };
    let centre = |t| Ok(run.station(t, tol)?[0]);
    let rolling = crate::between::Rolling {
        radius: run.radius,
        stations: stations.clone(),
        centre: &centre,
        lines: vec![(run.guide.clone(), run.range), run.rim.clone()],
        ends: [run.range.0, run.range.1]
            .into_iter()
            .flat_map(|t| {
                let [_, a, b] = run.station(t, tol).expect("validated shelf station");
                [a, b]
            })
            .collect(),
        lunes: Some(&search),
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
            "a face at {at:?} enters a carried shelf's ball"
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
                first: -1,
                second: -1,
            },
            stopped: BlendStop::RanPastTheGuide,
        };
        for i in 0..=STATIONS {
            let t = at(span, i);
            let [c, wall, cap] = run.station(t, tol)?;
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
        .into_iter()
        .map(|span| {
            Box::new(move |p: Point| {
                let mut best = 0;
                let mut distance = f64::INFINITY;
                for i in 0..=STATIONS {
                    let d = p.distance(
                        run.station(at(span, i), tol)
                            .expect("validated shelf station")[0],
                    );
                    if d < distance {
                        best = i;
                        distance = d;
                    }
                }
                let mut lo = at(span, best.saturating_sub(1));
                let mut hi = at(span, (best + 1).min(STATIONS));
                for _ in 0..64 {
                    let a = lo + (hi - lo) / 3.0;
                    let b = hi - (hi - lo) / 3.0;
                    if p.distance(run.station(a, tol).expect("validated shelf station")[0])
                        < p.distance(run.station(b, tol).expect("validated shelf station")[0])
                    {
                        hi = b;
                    } else {
                        lo = a;
                    }
                }
                (p.distance(
                    run.station(f64::midpoint(lo, hi), tol)
                        .expect("validated shelf station")[0],
                ) - run.radius)
                    .abs()
            }) as Box<dyn Fn(Point) -> f64>
        })
        .collect();
    let exact: Vec<_> = offs
        .iter()
        .map(|off| ExactBand {
            rails: [
                Rail::Plane(run.origin + run.out * run.wall, run.out),
                Rail::Plane(run.origin, run.up),
            ],
            off: off.as_ref(),
            strips: &[],
        })
        .collect();
    crate::marched::build_chain_band_signed(
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
        true,
        tol,
    )
}
