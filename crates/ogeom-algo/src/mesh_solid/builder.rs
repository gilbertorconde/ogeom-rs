//! The planned vertices, edges and faces built: planar faces, fans and
//! wedges between facets and curved faces, and curved faces bounded in one
//! chart or wrapped round their axis as bands.

use ogeom_core::{FastMap, OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, LineCurve, PlanarCurve, PlaneSurface};
use ogeom_math::{Direction, Frame, Plane, Point, Point2, Vector};
use ogeom_topo::{EdgeData, EdgeRepr, FaceData, Location, Model, Shape, VertexData};

use super::planner::{Corner, EdgeSpec, Plan, TOLERANCE_MARGIN};
use super::seams::{
    RimStart, centred_rings, choose_seam, hole_polygons, ring_area, round_the_tube, swapped,
    swapped_rings, turns_along,
};
use super::segment::{axis_frame, chart, evaluate, unit_normal, unwrapped, walk_entries};
use super::snap::projected_image;
use super::weld::{Half, from_to};
use super::{Carrier, Curved, Fan, Groups, Layout};
use crate::recognize::Canonical;

/// Builds the planned vertices, edges and faces.
pub(super) struct Builder<'a> {
    pub(super) model: &'a mut Model,
    pub(super) points: &'a [Point],
    pub(super) triangles: &'a [[u32; 3]],
    pub(super) groups: &'a Groups,
    pub(super) plan: &'a Plan,
    pub(super) fans: &'a FastMap<usize, Fan>,
    pub(super) tol: Tolerances,
}

impl Builder<'_> {
    pub(super) fn entry(&self, h: Half) -> (usize, bool) {
        let (a, b) = from_to(self.triangles, h);
        let (edge, along) = self.plan.edge_of[&(a.min(b), a.max(b))];
        (edge, along == (a < b))
    }

    /// The faces, one per live group, by group index, and the groups whose
    /// face could not be built: an edge that cannot be made fails the
    /// groups it bounds, and builds no face.
    pub(super) fn build(mut self) -> (Vec<Option<Shape>>, Vec<usize>) {
        let edges = match self.edges() {
            Ok(edges) => edges,
            Err(edge) => {
                let mut failed: Vec<usize> = (0..self.triangles.len() * 3)
                    .filter(|&h| {
                        let (a, b) = from_to(self.triangles, h);
                        self.plan
                            .edge_of
                            .get(&(a.min(b), a.max(b)))
                            .is_some_and(|&(e, _)| e == edge)
                    })
                    .map(|h| self.groups.of[h / 3])
                    .collect();
                failed.sort_unstable();
                failed.dedup();
                return (Vec::new(), failed);
            }
        };
        let (edges, corners) = edges;
        let mut faces = Vec::with_capacity(self.groups.carriers.len());
        let mut failed = Vec::new();
        for g in 0..self.groups.carriers.len() {
            if let Ok(face) = self.face(g, &edges, &corners) {
                faces.push(face);
            } else {
                failed.push(g);
                faces.push(None);
            }
        }
        (faces, failed)
    }

    /// The planned edges and their corners, or the index of the first edge
    /// that cannot be made.
    fn edges(&mut self) -> Result<(Vec<Shape>, FastMap<Corner, Shape>), usize> {
        let mut corners: FastMap<Corner, Shape> = FastMap::default();
        let mut edges: Vec<Shape> = Vec::with_capacity(self.plan.edges.len());
        for (index, spec) in self.plan.edges.iter().enumerate() {
            for corner in spec.ends {
                corners.entry(corner).or_insert_with(|| {
                    let at = match corner {
                        Corner::Mesh(v) => self.points[v as usize],
                        Corner::Placed(i) => self.plan.placed[i],
                    };
                    self.model.add_vertex(VertexData::new(at))
                });
            }
            let id = self.model.geometry_mut().add_curve(spec.curve.clone());
            let mut data = EdgeData::on_curve(id, Location::identity(), spec.range);
            // A tolerance measured as a gap is held a millionth wider: the
            // checker measures the same gap again, a rounding apart.
            data.tolerance =
                Tolerance::new(spec.tolerance * (1.0 + TOLERANCE_MARGIN)).map_err(|_| index)?;
            let bounds = if spec.ends[0] == spec.ends[1] {
                vec![
                    corners[&spec.ends[0]].clone(),
                    corners[&spec.ends[0]].clone(),
                ]
            } else {
                vec![
                    corners[&spec.ends[0]].clone(),
                    corners[&spec.ends[1]].clone(),
                ]
            };
            edges.push(self.model.add_edge(data, &bounds).map_err(|_| index)?);
        }
        Ok((edges, corners))
    }

    /// Group `g`'s face on the built edges; `None` for a group with no face.
    fn face(
        &mut self,
        g: usize,
        edges: &[Shape],
        corners: &FastMap<Corner, Shape>,
    ) -> OgeomResult<Option<Shape>> {
        let (groups, plan) = (self.groups, self.plan);
        let rings = &plan.loops[g];
        Ok(match &groups.carriers[g] {
            Carrier::Gone => None,
            Carrier::Curved(curved) if plan.layouts[g] == Layout::Whole => {
                Some(self.whole_face(curved, g, 0.0)?)
            }
            _ if rings.is_empty() => None,
            Carrier::Plane(plane) => match self.fans.get(&g) {
                Some(fan) => match self.fan_face(*fan, rings, edges, corners)? {
                    Some(face) => Some(face),
                    None => Some(self.plane_face(*plane, rings, edges)?),
                },
                None => Some(self.plane_face(*plane, rings, edges)?),
            },
            Carrier::Curved(curved) => match plan.layouts[g] {
                Layout::Band { round_tube } => {
                    Some(self.band_face(curved, rings, edges, round_tube)?)
                }
                Layout::Cap => Some(match curved.shape {
                    Canonical::Cone(_) => self.cone_tip_face(curved, rings, edges)?,
                    _ => self.cap_face(curved, rings, edges)?,
                }),
                Layout::HoledCap => Some(self.holed_cap_face(curved, g, rings, edges)?),
                Layout::Wrapped => Some(self.wrapped_face(curved, g, rings, edges)?),
                Layout::Holed => Some(self.holed_face(curved, g, rings, edges)?),
                Layout::Threaded => Some(self.threaded_face(curved, g, rings, edges)?),
                Layout::Open | Layout::Whole => Some(self.curved_face(curved, g, rings, edges)?),
            },
        })
    }

    /// A ring's planned edges in walking order, repeats run together.
    pub(super) fn entries(&self, ring: &[Half]) -> Vec<(usize, bool)> {
        walk_entries(ring, |h| self.entry(h))
    }

    pub(super) fn has_pcurve(&self, edge: &Shape, surface: ogeom_topo::SurfaceId) -> bool {
        self.model.node(edge).and_then(|n| n.data().as_edge()).is_some_and(|d| {
            d.representations.iter().any(|rep| {
                matches!(rep, EdgeRepr::PCurve { surface: s, .. } | EdgeRepr::Seam { surface: s, .. } if *s == surface)
            })
        })
    }

    /// A fan's face: the B-spline surface ruled from the apex (its `v = 0`
    /// row) to the seam's curve (its `v = 1` row), over the seam's own
    /// parameter. The seam's image is its row, each straight side's the
    /// ruling at its end, and the apex an edge with no length along the
    /// `v = 0` row. `None` where the face is not the triangle the fan was
    /// planned from or its seam is no polynomial spline, and the facet is
    /// built flat.
    fn fan_face(
        &mut self,
        fan: Fan,
        rings: &[Vec<Half>],
        edges: &[Shape],
        corners: &FastMap<Corner, Shape>,
    ) -> OgeomResult<Option<Shape>> {
        use ogeom_geom::Surface as _;
        if fan.across.is_some() {
            return self.wedge_face(fan, rings, edges, corners);
        }
        let [ring] = rings else {
            return Ok(None);
        };
        let entries = self.entries(ring);
        if entries.len() != 3 {
            return Ok(None);
        }
        let (a, b) = fan.seam;
        let Some(&(seam, _)) = self.plan.edge_of.get(&(a.min(b), a.max(b))) else {
            return Ok(None);
        };
        let Some(spec) = self.plan.edges.get(seam) else {
            return Ok(None);
        };
        let Curve::BSpline(spline) = &spec.curve else {
            return Ok(None);
        };
        if spline.control_points().iter().any(|c| c.weight != 1.0) {
            return Ok(None);
        }
        let Corner::Mesh(start) = spec.ends[0] else {
            return Ok(None);
        };
        let (t0, t1) = spec.range;
        let apex = self.points[fan.apex as usize];
        let mut net = Vec::with_capacity(spline.control_points().len() * 2);
        for c in spline.control_points() {
            net.push(apex);
            net.push(c.point());
        }
        let count = spline.control_points().len();
        let grid = ogeom_math::ControlGrid::new(net, count, 2)?;
        let along = ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
        let geometry: ogeom_geom::SurfaceGeometry =
            ogeom_geom::BSplineSurface::new(spline.knots().clone(), along, &grid, self.tol)?.into();
        let surface = self.model.geometry_mut().add_surface(geometry.clone());
        // A straight piece of the chart from `p` to `q` over `range`.
        let straight = |p: Point2, q: Point2, range: (f64, f64)| -> OgeomResult<PlanarCurve> {
            Ok(ogeom_geom::BSpline2d::new(
                ogeom_math::KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?,
                vec![p, q],
                self.tol,
            )?
            .into())
        };
        // Where each corner stands on the chart's top row.
        let top = |v: u32| if v == start { t0 } else { t1 };
        let mut wire = Vec::with_capacity(4);
        for (k, &(edge, forward)) in entries.iter().enumerate() {
            let spec = &self.plan.edges[edge];
            let ends = spec.ends.map(|c| match c {
                Corner::Mesh(v) => Some(v),
                Corner::Placed(_) => None,
            });
            let [Some(from), Some(to)] = ends else {
                return Ok(None);
            };
            let pcurve = if edge == seam {
                straight(Point2::new(t0, 1.0), Point2::new(t1, 1.0), spec.range)?
            } else if from == fan.apex && (to == a || to == b) {
                let u = top(to);
                straight(Point2::new(u, 0.0), Point2::new(u, 1.0), spec.range)?
            } else if to == fan.apex && (from == a || from == b) {
                let u = top(from);
                straight(Point2::new(u, 1.0), Point2::new(u, 0.0), spec.range)?
            } else {
                return Ok(None);
            };
            crate::build::attach_pcurve(
                self.model,
                &edges[edge],
                pcurve,
                surface,
                Location::identity(),
                spec.range,
            )?;
            wire.push(oriented(&edges[edge], forward));
            // The apex closes between the side that runs into it and the
            // side that runs out.
            let reaches = if forward { to } else { from };
            if reaches == fan.apex {
                let (next, next_forward) = entries[(k + 1) % entries.len()];
                let next_spec = &self.plan.edges[next];
                let leaves = if next_forward {
                    next_spec.ends[1]
                } else {
                    next_spec.ends[0]
                };
                let Corner::Mesh(leaves) = leaves else {
                    return Ok(None);
                };
                let Some(vertex) = corners.get(&Corner::Mesh(fan.apex)) else {
                    return Ok(None);
                };
                let (p, q) = (
                    Point2::new(top(if forward { from } else { to }), 0.0),
                    Point2::new(top(leaves), 0.0),
                );
                let mut data = EdgeData::new();
                data.degenerate = true;
                let pole = self
                    .model
                    .add_edge(data, &[vertex.clone(), vertex.clone()])?;
                crate::build::attach_pcurve(
                    self.model,
                    &pole,
                    ogeom_geom::Line2d::segment(p, q, self.tol)?.into(),
                    surface,
                    Location::identity(),
                    (0.0, p.distance(q)),
                )?;
                wire.push(pole);
            }
        }
        // Outward where the surface's normal mid-way agrees with the facet's.
        let mid = (f64::midpoint(t0, t1), 0.5);
        let normal = geometry.normal_at(mid.0, mid.1, self.tol)?;
        let outward = normal
            .vector()
            .dot(unit_normal(self.points, [a, b, fan.apex]))
            >= 0.0;
        // The triangle's walk keeps the facet on its left about the facet's
        // normal; a face turned against its surface stores it walked back.
        let wire = if outward { wire } else { walked_back(&wire) };
        let wire = self.model.add_wire(&wire)?;
        let face = self
            .model
            .add_face(FaceData::new(surface, Location::identity()), &[wire])?;
        Ok(Some(if outward { face } else { face.reversed() }))
    }

    /// A wedge's face (see [`Fan`]): the B-spline surface ruled between its
    /// two seams, the second (from the shared corner to the apex) its
    /// `v = 0` row and the first its `v = 1` row, both run from the shared
    /// corner over their common parameter. Each seam's image is its row,
    /// the third side's the ruling at the far end, and the shared corner an
    /// edge with no length along the first column. `None` where the face is
    /// not the triangle the wedge was planned from, or its seams are not
    /// polynomial splines on one knot vector the same both ways, and the
    /// facet is built flat.
    fn wedge_face(
        &mut self,
        fan: Fan,
        rings: &[Vec<Half>],
        edges: &[Shape],
        corners: &FastMap<Corner, Shape>,
    ) -> OgeomResult<Option<Shape>> {
        use ogeom_geom::Surface as _;
        let [ring] = rings else {
            return Ok(None);
        };
        let entries = self.entries(ring);
        if entries.len() != 3 {
            return Ok(None);
        }
        let (a, shared, b) = (fan.seam.0, fan.seam.1, fan.apex);
        let seam_of = |p: u32, q: u32| self.plan.edge_of.get(&(p.min(q), p.max(q))).map(|e| e.0);
        let (Some(first), Some(second)) = (seam_of(shared, a), seam_of(shared, b)) else {
            return Ok(None);
        };
        // A seam's control points run from the shared corner, its knots and
        // range.
        let from_shared =
            |edge: usize| -> Option<(ogeom_math::KnotVector, Vec<Point>, (f64, f64))> {
                let spec = self.plan.edges.get(edge)?;
                let Curve::BSpline(spline) = &spec.curve else {
                    return None;
                };
                if spline.control_points().iter().any(|c| c.weight != 1.0) {
                    return None;
                }
                let Corner::Mesh(start) = spec.ends[0] else {
                    return None;
                };
                let mut net: Vec<Point> =
                    spline.control_points().iter().map(|c| c.point()).collect();
                if start != shared {
                    net.reverse();
                }
                Some((spline.knots().clone(), net, spec.range))
            };
        let (Some((knots, top, range)), Some((other, bottom, other_range))) =
            (from_shared(first), from_shared(second))
        else {
            return Ok(None);
        };
        let (t0, t1) = range;
        let values = knots.knots();
        let even = values
            .iter()
            .zip(values.iter().rev())
            .all(|(x, y)| (x + y - t0 - t1).abs() <= 1e-12 * (t1 - t0).abs().max(1.0));
        if knots != other || range != other_range || top.len() != bottom.len() || !even {
            return Ok(None);
        }
        let count = top.len();
        let mut net = Vec::with_capacity(count * 2);
        for (low, high) in bottom.iter().zip(&top) {
            net.push(*low);
            net.push(*high);
        }
        let grid = ogeom_math::ControlGrid::new(net, count, 2)?;
        let along = ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
        let geometry: ogeom_geom::SurfaceGeometry =
            ogeom_geom::BSplineSurface::new(knots, along, &grid, self.tol)?.into();
        let surface = self.model.geometry_mut().add_surface(geometry.clone());
        // Where each corner stands in the chart, read from the edge `edge`.
        let at = |edge: usize, v: u32| -> Option<Point2> {
            if v == shared {
                (edge == first || edge == second)
                    .then(|| Point2::new(t0, if edge == first { 1.0 } else { 0.0 }))
            } else if v == a {
                Some(Point2::new(t1, 1.0))
            } else if v == b {
                Some(Point2::new(t1, 0.0))
            } else {
                None
            }
        };
        let mut wire = Vec::with_capacity(4);
        for (k, &(edge, forward)) in entries.iter().enumerate() {
            let spec = &self.plan.edges[edge];
            let [Corner::Mesh(from), Corner::Mesh(to)] = spec.ends else {
                return Ok(None);
            };
            let (Some(p), Some(q)) = (at(edge, from), at(edge, to)) else {
                return Ok(None);
            };
            let pcurve: PlanarCurve = ogeom_geom::BSpline2d::new(
                ogeom_math::KnotVector::new(
                    vec![spec.range.0, spec.range.0, spec.range.1, spec.range.1],
                    1,
                )?,
                vec![p, q],
                self.tol,
            )?
            .into();
            crate::build::attach_pcurve(
                self.model,
                &edges[edge],
                pcurve,
                surface,
                Location::identity(),
                spec.range,
            )?;
            wire.push(oriented(&edges[edge], forward));
            // The shared corner closes between the seam that runs into it
            // and the one that runs out.
            let reaches = if forward { to } else { from };
            if reaches == shared {
                let (next, _) = entries[(k + 1) % entries.len()];
                let (Some(p), Some(q), Some(vertex)) = (
                    at(edge, shared),
                    at(next, shared),
                    corners.get(&Corner::Mesh(shared)),
                ) else {
                    return Ok(None);
                };
                let mut data = EdgeData::new();
                data.degenerate = true;
                let pole = self
                    .model
                    .add_edge(data, &[vertex.clone(), vertex.clone()])?;
                crate::build::attach_pcurve(
                    self.model,
                    &pole,
                    ogeom_geom::Line2d::segment(p, q, self.tol)?.into(),
                    surface,
                    Location::identity(),
                    (0.0, p.distance(q)),
                )?;
                wire.push(pole);
            }
        }
        // Outward where the surface's normal mid-way agrees with the facet's.
        let normal = geometry.normal_at(f64::midpoint(t0, t1), 0.5, self.tol)?;
        let outward = normal
            .vector()
            .dot(unit_normal(self.points, [a, shared, b]))
            >= 0.0;
        // The triangle's walk keeps the facet on its left about the facet's
        // normal; a face turned against its surface stores it walked back.
        let wire = if outward { wire } else { walked_back(&wire) };
        let wire = self.model.add_wire(&wire)?;
        let face = self
            .model
            .add_face(FaceData::new(surface, Location::identity()), &[wire])?;
        Ok(Some(if outward { face } else { face.reversed() }))
    }

    fn plane_face(
        &mut self,
        plane: Plane,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let geometry: ogeom_geom::SurfaceGeometry = PlaneSurface::new(plane).into();
        let surface = self.model.geometry_mut().add_surface(geometry.clone());
        let local = |p: Point| {
            let l = plane.frame().to_local(p);
            Point2::new(l.x, l.y)
        };
        let mut wires: Vec<(f64, Shape)> = Vec::with_capacity(rings.len());
        for ring in rings {
            let area = ring_area(ring, self.triangles, |p| Some(local(p)), self.points);
            let mut ring_edges = Vec::new();
            for (edge, forward) in self.entries(ring) {
                let spec = &self.plan.edges[edge];
                if !self.has_pcurve(&edges[edge], surface) {
                    let Some(pcurve) = self.plane_image(plane, &geometry, spec)? else {
                        ogeom_bail!(Construction, "an edge has no image in its face's plane");
                    };
                    crate::build::attach_pcurve(
                        self.model,
                        &edges[edge],
                        pcurve,
                        surface,
                        Location::identity(),
                        spec.range,
                    )?;
                }
                ring_edges.push(oriented(&edges[edge], forward));
            }
            wires.push((area, self.model.add_wire(&ring_edges)?));
        }
        wires.sort_by(|a, b| b.0.total_cmp(&a.0));
        let wires: Vec<Shape> = wires.into_iter().map(|(_, w)| w).collect();
        self.model
            .add_face(FaceData::new(surface, Location::identity()), &wires)
    }

    /// An edge's image in a plane face's chart. The closed-form image keeps
    /// a curve's own shape and assumes it lies in the plane, so it is kept
    /// only where it lands within the edge's tolerance: an arc fitted on a
    /// neighbouring face may cross the plane rather than lie in it, and
    /// its closed-form image is then a circle the edge never runs along.
    /// Otherwise the image is the curve's projection, interpolated at the
    /// edge's own parameters.
    fn plane_image(
        &self,
        plane: Plane,
        geometry: &ogeom_geom::SurfaceGeometry,
        spec: &EdgeSpec,
    ) -> OgeomResult<Option<PlanarCurve>> {
        if let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&spec.curve, geometry, self.tol)
            && crate::pcurve_gap::pcurve_gap(
                (&spec.curve, spec.range),
                (&pcurve, spec.range),
                geometry,
                spec.tolerance,
                self.tol,
            )? <= spec.tolerance
        {
            return Ok(Some(pcurve));
        }
        Ok(projected_image(plane, &spec.curve, spec.range, self.tol))
    }

    /// Whether the region's outward side is the surface's own normal side.
    pub(super) fn outward(&self, curved: &Curved, g: usize) -> bool {
        let mut vote = 0.0;
        for (t, tri) in self.triangles.iter().enumerate() {
            if self.groups.of[t] != g {
                continue;
            }
            let corners = tri.map(|v| self.points[v as usize]);
            let centroid = Point::from_vector(
                (corners[0].to_vector() + corners[1].to_vector() + corners[2].to_vector()) / 3.0,
            );
            let n = unit_normal(self.points, *tri);
            vote += n.dot(self.surface_normal(curved, centroid));
        }
        vote >= 0.0
    }

    /// The surface's own normal (the chart's `du × dv`) near a point.
    fn surface_normal(&self, curved: &Curved, p: Point) -> Vector {
        let Some(at) = chart(&curved.shape, p, self.tol) else {
            return Vector::ZERO;
        };
        let e = 1e-6;
        let o = evaluate(&curved.shape, at);
        let du = evaluate(&curved.shape, (at.0 + e, at.1)) - o;
        let dv = evaluate(&curved.shape, (at.0, at.1 + e)) - o;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m > 0.0 { n / m } else { Vector::ZERO }
    }

    fn curved_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(
                Construction,
                "a curved face was planned without its surface"
            );
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let mut wires: Vec<(f64, Shape)> = Vec::with_capacity(rings.len());
        for ring in rings {
            let area = ring_area(
                ring,
                self.triangles,
                |p| unwrapped(curved, p, self.tol).map(|(u, v)| Point2::new(u, v)),
                self.points,
            );
            let mut ring_edges = Vec::new();
            for (edge, forward) in self.entries(ring) {
                let spec = &self.plan.edges[edge];
                if !self.has_pcurve(&edges[edge], surface) {
                    let Some((pcurve, deviation)) = self.plan.pcurves.get(&(edge, g)).cloned()
                    else {
                        ogeom_bail!(
                            Construction,
                            "an edge was planned without its image on a face"
                        );
                    };
                    self.model.widen(
                        &edges[edge],
                        Tolerance::new(deviation.max(self.tol.confusion()))?,
                    )?;
                    crate::build::attach_pcurve(
                        self.model,
                        &edges[edge],
                        pcurve,
                        surface,
                        Location::identity(),
                        spec.range,
                    )?;
                }
                ring_edges.push(oriented(&edges[edge], forward));
            }
            // The triangles' walk keeps the region on its left about its
            // outward side; a face turned against its surface stores each
            // ring walked back, so it keeps the face on its left in the chart.
            let sign = if outward { 1.0 } else { -1.0 };
            let ring_edges = if outward {
                ring_edges
            } else {
                walked_back(&ring_edges)
            };
            wires.push((area * sign, self.model.add_wire(&ring_edges)?));
        }
        wires.sort_by(|a, b| b.0.total_cmp(&a.0));
        let wires: Vec<Shape> = wires.into_iter().map(|(_, w)| w).collect();
        // Fitted images answer on whichever branch their projection chose.
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A face round its axis between two rims of any shape, with holes.
    ///
    /// The seam joins a vertex of one rim to a vertex of the other along a
    /// straight line in the chart: the pair turning least between them
    /// whose line crosses no hole. It is a ruling where the pair stand at
    /// one angle on a cylinder or a cone, and otherwise the curve that line
    /// traces on the surface, interpolated. The outer wire walks the seam
    /// down, one rim round, the seam up a whole turn over, and the other rim
    /// back; each hole is its own wire.
    fn wrapped_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let tau = core::f64::consts::TAU;
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(
                Construction,
                "a face round its axis was planned without its surface"
            );
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        // Every edge's image, as the plan made it.
        for ring in rings {
            for (edge, _) in self.entries(ring) {
                if self.has_pcurve(&edges[edge], surface) {
                    continue;
                }
                let Some((pcurve, deviation)) = self.plan.pcurves.get(&(edge, g)).cloned() else {
                    ogeom_bail!(
                        Construction,
                        "an edge was planned without its image on a face"
                    );
                };
                self.model.widen(
                    &edges[edge],
                    Tolerance::new(deviation.max(self.tol.confusion()))?,
                )?;
                crate::build::attach_pcurve(
                    self.model,
                    &edges[edge],
                    pcurve,
                    surface,
                    Location::identity(),
                    self.plan.edges[edge].range,
                )?;
            }
        }
        let round_tube = round_the_tube(curved);
        let windings: Vec<i32> = rings
            .iter()
            .map(|ring| {
                turns_along(curved, ring, self.triangles, self.points, self.tol).unwrap_or(0)
            })
            .collect();
        let rims: Vec<usize> = (0..rings.len())
            .filter(|&k| windings[k].abs() == 1)
            .collect();
        let [low, high] = rims[..] else {
            ogeom_bail!(Construction, "a face round its axis has two rims");
        };
        if windings[low] != -windings[high] {
            ogeom_bail!(
                Construction,
                "a face round its axis has its rims turning one way"
            );
        }
        let holes: Vec<usize> = (0..rings.len()).filter(|&k| windings[k] == 0).collect();

        // Each rim's entries, and the vertex each one starts from.
        let starts = |this: &Self, ring: &[Half]| -> OgeomResult<Vec<RimStart>> {
            let mut out = Vec::new();
            for (edge, forward) in this.entries(ring) {
                let ends = this.model.children_of(&edges[edge])?;
                let vertex = if forward { ends.first() } else { ends.last() };
                let Some(vertex) = vertex.cloned() else {
                    ogeom_bail!(Construction, "a rim edge has no vertex");
                };
                let Some(ogeom_topo::NodeData::Vertex(data)) =
                    this.model.node(&vertex).map(|n| n.data())
                else {
                    ogeom_bail!(Construction, "a rim vertex has no position");
                };
                out.push(((edge, forward), vertex, data.point));
            }
            Ok(out)
        };
        let from = starts(self, &rings[low])?;
        let to = starts(self, &rings[high])?;
        let hole_rings = swapped_rings(
            centred_rings(
                curved,
                hole_polygons(
                    &curved.shape,
                    &holes
                        .iter()
                        .map(|&k| rings[k].as_slice())
                        .collect::<Vec<_>>(),
                    self.triangles,
                    self.points,
                    self.tol,
                ),
            ),
            round_tube,
        );
        let from_at: Vec<Point> = from.iter().map(|x| x.2).collect();
        let to_at: Vec<Point> = to.iter().map(|x| x.2).collect();
        let Some((i, j, a, b)) =
            choose_seam(curved, round_tube, &from_at, &to_at, &hole_rings, self.tol)
        else {
            ogeom_bail!(Construction, "no seam joins the rims clear of the holes");
        };
        // Back in the chart's own order: the angle round the axis first.
        let (a, b) = (swapped(a, round_tube), swapped(b, round_tube));
        let (pa, pb) = (from[i].2, to[j].2);
        let straight = (b.0 - a.0).abs() <= 1e-12
            && matches!(curved.shape, Canonical::Cylinder(_) | Canonical::Cone(_));
        let (seam_curve, range, deviation): (Curve, (f64, f64), f64) = if straight {
            let line = LineCurve::segment(pa, pb, self.tol)?;
            (line.into(), (0.0, pa.distance(pb)), self.tol.confusion())
        } else {
            chart_trace(&curved.shape, (a, b), (pa, pb), self.tol)?
        };
        let id = self.model.geometry_mut().add_curve(seam_curve);
        let mut data = EdgeData::on_curve(id, Location::identity(), range);
        data.tolerance = Tolerance::new(deviation)?;
        let seam = self
            .model
            .add_edge(data, &[from[i].1.clone(), to[j].1.clone()])?;
        // Down its near side where the wire leaves the high rim, up its far
        // side a whole turn over, where the low rim's walk comes round to.
        // A face turned against its surface walks its wires back, and so
        // each side of the seam the other way.
        let over = swapped((tau * f64::from(windings[low]), 0.0), round_tube);
        let near = linear(a, b, range, self.tol)?;
        let far = linear(
            (a.0 + over.0, a.1 + over.1),
            (b.0 + over.0, b.1 + over.1),
            range,
            self.tol,
        )?;
        let (forward, back) = if outward { (far, near) } else { (near, far) };
        crate::build::attach_seam(
            self.model,
            &seam,
            forward,
            back,
            surface,
            Location::identity(),
            range,
        )?;
        let rotated = |list: &[RimStart], at: usize| -> Vec<Shape> {
            (0..list.len())
                .map(|k| {
                    let ((edge, forward), _, _) = list[(at + k) % list.len()];
                    oriented(&edges[edge], forward)
                })
                .collect()
        };
        let mut outer = vec![seam.reversed()];
        outer.extend(rotated(&from, i));
        outer.push(seam.clone());
        outer.extend(rotated(&to, j));
        let mut lists = vec![outer];
        for &k in &holes {
            lists.push(
                self.entries(&rings[k])
                    .into_iter()
                    .map(|(edge, forward)| oriented(&edges[edge], forward))
                    .collect(),
            );
        }
        // The triangles' walk keeps the region on its left about its
        // outward side; a face turned against its surface stores it walked
        // back, so it keeps the face on its left in the chart.
        let mut wires = Vec::with_capacity(lists.len());
        for list in lists {
            let list = if outward { list } else { walked_back(&list) };
            wires.push(self.model.add_wire(&list)?);
        }
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A band round the axis: the two rim circles, and a seam at the
    /// surface's angle zero joining their vertices.
    fn band_face(
        &mut self,
        curved: &Curved,
        rings: &[Vec<Half>],
        edges: &[Shape],
        round_tube: bool,
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let g = self.groups.of[rings[0][0] / 3];
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a band was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let Some(frame) = axis_frame(&curved.shape) else {
            ogeom_bail!(Construction, "a band has no axis");
        };
        // Each rim: its edge, its chart height, and whether its parameter
        // runs with the surface's angle.
        let mut rims = Vec::with_capacity(2);
        for ring in rings {
            let (edge, _) = self.entry(ring[0]);
            let spec = &self.plan.edges[edge];
            let Curve::Circle(c) = &spec.curve else {
                ogeom_bail!(Construction, "a band's rim is not a circle");
            };
            let circle = c.circle();
            let start = spec.curve.point_at(0.0, self.tol)?;
            let Some((u, v)) = unwrapped(curved, start, self.tol) else {
                ogeom_bail!(Construction, "a band's rim has no chart position");
            };
            // Whether the rim's own parameter runs with the chart's angle
            // it goes round: the axis's for a band round it, the tube's for
            // a band round the tube.
            let turning = if round_tube {
                let radial = circle.centre() - frame.origin();
                radial.cross(frame.z().vector())
            } else {
                frame.z().vector()
            };
            let with = circle.frame().z().vector().dot(turning) > 0.0;
            rims.push((edge, if round_tube { u } else { v }, with, start));
        }
        rims.sort_by(|a, b| a.1.total_cmp(&b.1));
        let [
            (low, v_low, low_with, low_at),
            (high, v_high, high_with, high_at),
        ] = rims[..]
        else {
            ogeom_bail!(Construction, "a band has two rims");
        };
        for (edge, at, with) in [(low, v_low, low_with), (high, v_high, high_with)] {
            let (a, b) = if with { (0.0, tau) } else { (tau, 0.0) };
            let (from, to) = if round_tube {
                ((at, a), (at, b))
            } else {
                ((a, at), (b, at))
            };
            let pcurve = linear(from, to, (0.0, tau), self.tol)?;
            crate::build::attach_pcurve(
                self.model,
                &edges[edge],
                pcurve,
                surface,
                Location::identity(),
                (0.0, tau),
            )?;
        }
        // The seam, low rim to high, on the surface along angle zero: a
        // ruling, a meridian of a sphere, a circle of a torus's tube, or
        // for a band round the tube, an arc of its outer equator.
        let seam_curve: Curve = match curved.shape {
            Canonical::Torus(t) if round_tube => ogeom_geom::CircleCurve::new(
                ogeom_math::Circle::new(frame, t.major_radius() + t.minor_radius(), self.tol)?,
            )
            .into(),
            Canonical::Sphere(s) => {
                let normal =
                    Direction::new(frame.x().vector().cross(frame.z().vector()), self.tol)?;
                ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
                    Frame::new(s.centre(), normal, frame.x(), self.tol)?,
                    s.radius(),
                    self.tol,
                )?)
                .into()
            }
            Canonical::Torus(t) => {
                let spine = frame.origin() + frame.x().vector() * t.major_radius();
                let normal =
                    Direction::new(frame.x().vector().cross(frame.z().vector()), self.tol)?;
                ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
                    Frame::new(spine, normal, frame.x(), self.tol)?,
                    t.minor_radius(),
                    self.tol,
                )?)
                .into()
            }
            _ => LineCurve::segment(low_at, high_at, self.tol)?.into(),
        };
        let seam_range = match curved.shape {
            Canonical::Torus(_) | Canonical::Sphere(_) => (v_low, v_high),
            _ => (0.0, low_at.distance(high_at)),
        };
        // The seam runs between the rims' own vertices, placed where each
        // circle's parameter starts: the surface's angle zero.
        let placed = |at: Point| {
            self.plan
                .placed
                .iter()
                .any(|p| p.distance(at) <= self.tol.confusion())
        };
        if !placed(low_at) || !placed(high_at) {
            ogeom_bail!(
                Construction,
                "a band's rim vertex is not where its seam starts"
            );
        }
        let vertex = |edge: usize| -> OgeomResult<Shape> {
            match self.model.children_of(&edges[edge])?.first() {
                Some(v) => Ok(v.clone()),
                None => ogeom_bail!(Construction, "a rim has no vertex"),
            }
        };
        let (from_vertex, to_vertex) = (vertex(low)?, vertex(high)?);
        let id = self.model.geometry_mut().add_curve(seam_curve);
        let data = EdgeData::on_curve(id, Location::identity(), seam_range);
        let seam = self.model.add_edge(data, &[from_vertex, to_vertex])?;
        // The seam's two images: where the wire walks it forward, then back.
        let (forward, back) = if round_tube {
            (
                linear((v_low, 0.0), (v_high, 0.0), seam_range, self.tol)?,
                linear((v_low, tau), (v_high, tau), seam_range, self.tol)?,
            )
        } else {
            (
                linear((tau, v_low), (tau, v_high), seam_range, self.tol)?,
                linear((0.0, v_low), (0.0, v_high), seam_range, self.tol)?,
            )
        };
        crate::build::attach_seam(
            self.model,
            &seam,
            forward,
            back,
            surface,
            Location::identity(),
            seam_range,
        )?;
        // Counter-clockwise in the chart, whichever side the face faces:
        // a face turned against its surface keeps its material on the left
        // of its ring about the surface's normal. Round the axis: along
        // the low rim, up the seam's far side, back along the high rim, down
        // its near side. Round the tube: along the seam's near side, up the
        // high rim, back along the seam's far side, down the low rim.
        let ring = if round_tube {
            vec![
                seam.clone(),
                oriented(&edges[high], high_with),
                seam.reversed(),
                oriented(&edges[low], !low_with),
            ]
        } else {
            vec![
                oriented(&edges[low], low_with),
                seam.clone(),
                oriented(&edges[high], !high_with),
                seam.reversed(),
            ]
        };
        let wire = self.model.add_wire(&ring)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, std::slice::from_ref(&wire))?;
        Ok(if outward { face } else { face.reversed() })
    }
}

/// The curve a straight chart line from `a` to `b` traces on a surface,
/// interpolated and starting and ending exactly at `pa` and `pb`: the
/// curve, its range, and how far it strays from the trace.
pub(super) fn chart_trace(
    shape: &Canonical,
    (a, b): ((f64, f64), (f64, f64)),
    (pa, pb): (Point, Point),
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64), f64)> {
    use ogeom_geom::Curve3d as _;
    const SAMPLES: u32 = 96;
    let along = |f: f64| (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
    let mut pts: Vec<Point> = (0..=SAMPLES)
        .map(|k| evaluate(shape, along(f64::from(k) / f64::from(SAMPLES))))
        .collect();
    pts[0] = pa;
    pts[SAMPLES as usize] = pb;
    let curve: Curve = crate::fit::interpolate(&pts, 3, crate::fit::Spacing::Uniform, tol)?.into();
    let range = curve.domain();
    let mut deviation = tol.confusion();
    for k in 0..=(SAMPLES * 4) {
        let f = f64::from(k) / f64::from(SAMPLES * 4);
        let t = range.0 + (range.1 - range.0) * f;
        deviation = deviation.max(curve.point_at(t, tol)?.distance(evaluate(shape, along(f))));
    }
    Ok((curve, range, deviation))
}

/// A ring walked the other way round.
pub(super) fn walked_back(ring: &[Shape]) -> Vec<Shape> {
    ring.iter().rev().map(Shape::reversed).collect()
}

pub(super) fn oriented(edge: &Shape, forward: bool) -> Shape {
    if forward {
        edge.clone()
    } else {
        edge.reversed()
    }
}

/// A straight chart segment from `a` to `b`, linear in the edge's own
/// parameter over `range`.
pub(super) fn linear(
    a: (f64, f64),
    b: (f64, f64),
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    let knots = ogeom_math::KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?;
    Ok(ogeom_geom::BSpline2d::new(
        knots,
        vec![Point2::new(a.0, a.1), Point2::new(b.0, b.1)],
        tol,
    )?
    .into())
}
