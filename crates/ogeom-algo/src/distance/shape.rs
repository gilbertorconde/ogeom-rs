//! Distances to a shape's own geometry.

use std::sync::OnceLock;

use ogeom_core::{OgeomError, OgeomResult, Tolerances, parallel};
use ogeom_geom::Surface as _;
use ogeom_math::{Aabb, Point, Point2, Vector};
use ogeom_mesh::{Deflection, face_boundary_in_chart_met};
use ogeom_topo::{Model, Orientation, SameKey, Shape, ShapeType, explore_unique};

use super::{Closest, DistanceOptions, ShapeFeature, signed, unit, unit_or_zero};
use crate::classify::{Containment, SolidBoundary};
use crate::measure::{curve_bounds, project_on_curve, project_on_surface};
use crate::proximity::{Element, PointDistance, RING_ANGULAR, RING_CHORD, Trim, trim_side};

/// The signed distance from each point to `shape`'s faces, edges and
/// vertices, measured on their own surfaces, curves and points, with the
/// closest point, the feature it lies on and the normal there; in the
/// points' order, `None` where nothing is within [`DistanceOptions::max`].
///
/// The candidates are [`PointDistance`]'s: each vertex, the foot on each
/// edge, and the foot on each face's surface where it lands inside the
/// face. A foot within the trimming rings' polylining band of the face's
/// boundary is placed again against rings drawn a hundred times finer, and
/// again, down to the confusion distance, for as long as it is nearer than
/// every other candidate and the band cannot tell; one still in the band
/// where the rings can be drawn no finer is left to the edges and vertices
/// round it, as [`PointDistance::nearest`] leaves it. Equally near
/// candidates go to the vertex before the edge before the face, and among
/// those to the first the shape holds. A point that is not finite, or with
/// a facing limit a normal that has no direction, answers `None`. See the
/// [module](super) for how the sign is read.
///
/// A face's foot is its surface's nearest point to the query, as
/// [`project_on_surface`] finds it: on a surface with more than one local
/// nearest point, a nearer one inside the trim than the one found outside
/// it is not a candidate.
///
/// # Errors
///
/// As [`to_mesh`](super::to_mesh) for the options and normals; as
/// [`PointDistance::of`] for the shape, and as the projections and
/// [`SolidBoundary::holds`] report.
pub fn to_shape(
    model: &Model,
    shape: &Shape,
    points: &[Point],
    normals: Option<&[Vector]>,
    options: &DistanceOptions,
    tol: Tolerances,
) -> OgeomResult<Vec<Option<Closest<ShapeFeature>>>> {
    let facing = options.checked(points.len(), normals)?;
    let target = Target::of(model, shape, options.signed, tol)?;
    let answers = parallel::map_ordered(points, |i, &p| {
        let wanted = match (facing, normals) {
            (Some(cos), Some(normals)) => match unit(normals[i]) {
                Some(m) => Some((m, cos)),
                None => return Ok(None),
            },
            _ => None,
        };
        let faces = wanted.map(|(m, cos)| move |n: Vector| n.dot(m) >= cos);
        let faces = faces.as_ref().map(|f| f as &dyn Fn(Vector) -> bool);
        target.closest(model, p, faces, options, tol)
    });
    answers.into_iter().collect()
}

/// A face's trimming rings in its chart.
type Rings = Vec<Vec<Point2>>;

/// A shape's elements gathered once, with what each point asks of them.
pub(crate) struct Target {
    elements: PointDistance,
    /// Each element's box, widened by the confusion distance.
    bounds: Vec<Aabb>,
    /// For an edge or a vertex, the faces of the shape that hold it, by
    /// element index; empty for a face.
    faces_of: Vec<Vec<usize>>,
    /// For a face, its rings drawn at each of `chords`, as first needed;
    /// `None` inside where they cannot be drawn within the chord.
    finer: Vec<Vec<OnceLock<Option<Rings>>>>,
    /// The chords rings are drawn again at, finest last.
    chords: Vec<f64>,
    /// The shape's boundary prepared for classifying, where it is closed
    /// and the distances are signed.
    solid: Option<SolidBoundary>,
}

/// One candidate: its unsigned distance, its point, its element, and the
/// normal its side is read from.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Found {
    pub(crate) distance: f64,
    pub(crate) point: Point,
    element: usize,
    normal: Vector,
}

impl Target {
    /// The shape's elements gathered, with its boundary prepared for
    /// classifying where `sign` asks for signed distances.
    pub(crate) fn of(
        model: &Model,
        shape: &Shape,
        sign: bool,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let elements = PointDistance::of(model, shape, tol)?;
        let list = elements.elements();
        let mut places: ogeom_core::FastMap<SameKey, usize> = ogeom_core::FastMap::default();
        let mut bounds = Vec::with_capacity(list.len());
        let reach = tol.confusion().max(model.tolerances().confusion());
        for (i, element) in list.iter().enumerate() {
            let (held, bound) = match element {
                Element::Vertex(v, p) => (v, Aabb::of_point(*p)),
                Element::Edge(e, curve) => (e, curve_bounds(curve, tol)?),
                Element::Face(f, _) => (f, crate::tight::face_bounds(model, f)?),
            };
            places.insert(SameKey(held.clone()), i);
            bounds.push(bound.expanded(reach));
        }
        let mut faces_of = vec![Vec::new(); list.len()];
        let mut finer = Vec::with_capacity(list.len());
        let mut chords = Vec::new();
        let mut chord = RING_CHORD;
        while chord > tol.confusion() {
            chord = (chord / 100.0).max(tol.confusion());
            chords.push(chord);
        }
        for (i, element) in list.iter().enumerate() {
            let Element::Face(face, _) = element else {
                finer.push(Vec::new());
                continue;
            };
            finer.push(chords.iter().map(|_| OnceLock::new()).collect());
            for kind in [ShapeType::Edge, ShapeType::Vertex] {
                for held in explore_unique(model, face, kind)? {
                    if let Some(&at) = places.get(&SameKey(held)) {
                        faces_of[at].push(i);
                    }
                }
            }
        }
        let solid = if sign {
            let chord = crate::classify::exact_ring_chord(tol);
            match SolidBoundary::of(model, shape, chord, tol) {
                Ok(solid) => Some(solid),
                // An open boundary has no inside to classify against.
                Err(OgeomError::Construction(_)) => None,
                Err(other) => return Err(other),
            }
        } else {
            None
        };
        Ok(Self {
            elements,
            bounds,
            faces_of,
            finer,
            chords,
            solid,
        })
    }

    /// The oriented unit normal of face element `face` at `(u, v)` on its
    /// world surface, or `None` where the surface has none.
    fn normal(&self, face: usize, (u, v): (f64, f64), tol: Tolerances) -> Option<Vector> {
        let Element::Face(shape, prepared) = &self.elements.elements()[face] else {
            return None;
        };
        let n = Vector::from(prepared.world.normal_at(u, v, tol).ok()?);
        Some(if shape.orientation() == Orientation::Reversed {
            -n
        } else {
            n
        })
    }

    /// The unit normals at `q`, on the faces that hold edge or vertex
    /// element `element`, of those that have one there.
    fn normals_at(&self, element: usize, q: Point, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let mut out = Vec::with_capacity(self.faces_of[element].len());
        for &face in &self.faces_of[element] {
            let Element::Face(_, prepared) = &self.elements.elements()[face] else {
                continue;
            };
            let foot = project_on_surface(&prepared.world, q, 32, tol)?;
            out.extend(self.normal(face, foot.parameters, tol));
        }
        Ok(out)
    }

    fn closest(
        &self,
        model: &Model,
        p: Point,
        faces: Option<&dyn Fn(Vector) -> bool>,
        options: &DistanceOptions,
        tol: Tolerances,
    ) -> OgeomResult<Option<Closest<ShapeFeature>>> {
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return Ok(None);
        }
        let Some(found) = self.nearest(model, p, faces, options.max, tol)? else {
            return Ok(None);
        };
        let element = &self.elements.elements()[found.element];
        let mut normal = found.normal;
        let feature = match element {
            Element::Vertex(v, _) => ShapeFeature::Vertex(v.clone()),
            Element::Edge(e, _) => ShapeFeature::Edge(e.clone()),
            Element::Face(f, _) => ShapeFeature::Face(f.clone()),
        };
        let mut distance = signed(found.distance, p, found.point, normal, options.signed);
        if options.signed
            && matches!(feature, ShapeFeature::Vertex(_))
            && let Some(solid) = &self.solid
        {
            match solid.holds(model, p, tol)? {
                Containment::In => distance = -found.distance,
                Containment::Out => distance = found.distance,
                Containment::On => {}
            }
        }
        normal = unit_or_zero(normal);
        Ok(Some(Closest {
            distance,
            point: found.point,
            feature,
            normal,
        }))
    }

    /// The nearest candidate within `max` whose face's unit normal there
    /// `faces` accepts, where there is a filter (for an edge or a vertex,
    /// the normal of any face holding it): elements taken nearest box first
    /// and passed over once their box is farther than the best so far,
    /// which changes only the time.
    pub(crate) fn nearest(
        &self,
        model: &Model,
        p: Point,
        faces: Option<&dyn Fn(Vector) -> bool>,
        max: f64,
        tol: Tolerances,
    ) -> OgeomResult<Option<Found>> {
        let list = self.elements.elements();
        let mut order: Vec<(f64, usize)> = self
            .bounds
            .iter()
            .enumerate()
            .map(|(i, b)| (b.distance_to(p), i))
            .filter(|&(gap, _)| gap <= max)
            .collect();
        order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut best: Option<Found> = None;
        let better = |best: &Option<Found>, distance: f64, element: usize| match best {
            None => distance <= max,
            Some(held) => {
                distance < held.distance || (distance == held.distance && element < held.element)
            }
        };
        let kept = |normals: &[Vector]| faces.is_none_or(|faces| normals.iter().any(|&n| faces(n)));
        let sum = |normals: Vec<Vector>| normals.into_iter().fold(Vector::ZERO, |a, n| a + n);
        let mut banded: Vec<Found> = Vec::new();
        for (gap, i) in order {
            if best.is_some_and(|held| gap > held.distance) {
                break;
            }
            match &list[i] {
                Element::Vertex(_, q) => {
                    let distance = p.distance(*q);
                    if !better(&best, distance, i) {
                        continue;
                    }
                    let normals = self.normals_at(i, *q, tol)?;
                    if kept(&normals) {
                        best = Some(Found {
                            distance,
                            point: *q,
                            element: i,
                            normal: sum(normals),
                        });
                    }
                }
                Element::Edge(_, curve) => {
                    let foot = project_on_curve(curve, p, 64, tol)?;
                    if !better(&best, foot.distance, i) {
                        continue;
                    }
                    let normals = self.normals_at(i, foot.point, tol)?;
                    if kept(&normals) {
                        best = Some(Found {
                            distance: foot.distance,
                            point: foot.point,
                            element: i,
                            normal: sum(normals),
                        });
                    }
                }
                Element::Face(_, prepared) => {
                    let foot = project_on_surface(&prepared.world, p, 32, tol)?;
                    if !better(&best, foot.distance, i) {
                        continue;
                    }
                    // A foot where the surface has no normal (a pole, an
                    // apex) is a point of the edges and vertices there.
                    let Some(normal) = self.normal(i, foot.parameters, tol) else {
                        continue;
                    };
                    if !faces.is_none_or(|faces| faces(normal)) {
                        continue;
                    }
                    let found = Found {
                        distance: foot.distance,
                        point: foot.point,
                        element: i,
                        normal,
                    };
                    match trim_side(prepared, None, RING_CHORD, foot.point, tol)? {
                        Trim::Inside => best = Some(found),
                        Trim::Outside => {}
                        Trim::Band => banded.push(found),
                    }
                }
            }
        }
        banded.sort_by(|a, b| {
            a.distance
                .total_cmp(&b.distance)
                .then(a.element.cmp(&b.element))
        });
        for found in banded {
            if !better(&best, found.distance, found.element) {
                break;
            }
            if self.inside_finer(model, found.element, found.point, tol)? {
                best = Some(found);
                break;
            }
        }
        Ok(best)
    }

    /// Whether a foot the coarse rings put in their band is inside its
    /// face, placed against rings drawn finer until it leaves the band;
    /// `false` where it is outside, or still in the band where the rings
    /// can be drawn no finer.
    fn inside_finer(
        &self,
        model: &Model,
        face: usize,
        q: Point,
        tol: Tolerances,
    ) -> OgeomResult<bool> {
        let Element::Face(shape, prepared) = &self.elements.elements()[face] else {
            return Ok(false);
        };
        for (level, &chord) in self.chords.iter().enumerate() {
            let rings = self.finer[face][level].get_or_init(|| {
                // A curve's segments at a chord grow as the chord's inverse
                // square root, and so does the cap: the coarse rings' cap
                // scaled by as much.
                let coarse = Deflection::default().max_segments;
                #[allow(
                    clippy::cast_possible_truncation,
                    clippy::cast_sign_loss,
                    clippy::cast_precision_loss
                )]
                let cap = (coarse as f64 * (RING_CHORD / chord).sqrt()).ceil() as usize;
                let deflection = Deflection {
                    chord,
                    angular: RING_ANGULAR,
                    max_segments: cap,
                    ..Deflection::default()
                };
                // A boundary that will not draw finer leaves the foot in
                // the band, where the coarse rings left it.
                face_boundary_in_chart_met(model, shape, deflection, tol)
                    .ok()
                    .flatten()
            });
            let Some(rings) = rings else {
                return Ok(false);
            };
            match trim_side(prepared, Some(rings), chord, q, tol)? {
                Trim::Inside => return Ok(true),
                Trim::Outside => return Ok(false),
                Trim::Band => {}
            }
        }
        Ok(false)
    }
}
