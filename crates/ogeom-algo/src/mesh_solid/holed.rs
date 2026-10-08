//! Curved faces whose rings do not bound one region of the chart on their
//! own: sphere caps, cone tips, caps and faces with holes, faces threaded
//! along a seam chain through their holes, and whole closed surfaces.

use ogeom_core::{OgeomResult, Tolerance, ogeom_bail};
use ogeom_geom::{Curve, LineCurve};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{EdgeData, FaceData, Location, Shape, VertexData};

use super::Curved;
use super::builder::{Builder, chart_trace, linear, oriented, walked_back};
use super::seams::{
    cap_rings, cap_seam, centred_rings, hole_polygons, swapped, swapped_rings, thread_pieces,
    torus_seam_v,
};
use super::segment::{evaluate, unwrapped};
use super::weld::Half;
use crate::recognize::Canonical;

impl Builder<'_> {
    /// A sphere's cap: its rim, a meridian seam from the rim to the pole,
    /// and the pole as an edge of no length, the frame's axis pointing into
    /// the cap so the pole is the north one.
    pub(super) fn cap_face(
        &mut self,
        curved: &Curved,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let north = core::f64::consts::FRAC_PI_2;
        let g = self.groups.of[rings[0][0] / 3];
        let Canonical::Sphere(sphere) = curved.shape else {
            ogeom_bail!(Construction, "a cap is a sphere's");
        };
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a cap was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let frame = sphere.frame();
        let (rim, _) = self.entry(rings[0][0]);
        let spec = &self.plan.edges[rim];
        let Curve::Circle(c) = &spec.curve else {
            ogeom_bail!(Construction, "a cap's rim is not a circle");
        };
        let with = c.circle().frame().z().vector().dot(frame.z().vector()) > 0.0;
        let start = spec.curve.point_at(0.0, self.tol)?;
        let Some((_, v_rim)) = unwrapped(curved, start, self.tol) else {
            ogeom_bail!(Construction, "a cap's rim has no chart position");
        };
        let (a, b) = if with { (0.0, tau) } else { (tau, 0.0) };
        crate::build::attach_pcurve(
            self.model,
            &edges[rim],
            linear((a, v_rim), (b, v_rim), (0.0, tau), self.tol)?,
            surface,
            Location::identity(),
            (0.0, tau),
        )?;
        let Some(rim_vertex) = self.model.children_of(&edges[rim])?.first().cloned() else {
            ogeom_bail!(Construction, "a rim has no vertex");
        };
        let pole = self.model.add_vertex(VertexData::new(
            sphere.centre() + frame.z().vector() * sphere.radius(),
        ));
        let normal = Direction::new(frame.x().vector().cross(frame.z().vector()), self.tol)?;
        let meridian: Curve = ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
            Frame::new(sphere.centre(), normal, frame.x(), self.tol)?,
            sphere.radius(),
            self.tol,
        )?)
        .into();
        let id = self.model.geometry_mut().add_curve(meridian);
        let seam_range = (v_rim, north);
        let seam = self.model.add_edge(
            EdgeData::on_curve(id, Location::identity(), seam_range),
            &[rim_vertex, pole.clone()],
        )?;
        crate::build::attach_seam(
            self.model,
            &seam,
            linear((tau, v_rim), (tau, north), seam_range, self.tol)?,
            linear((0.0, v_rim), (0.0, north), seam_range, self.tol)?,
            surface,
            Location::identity(),
            seam_range,
        )?;
        let mut data = EdgeData::new();
        data.degenerate = true;
        let tip = self.model.add_edge(data, &[pole.clone(), pole])?;
        crate::build::attach_pcurve(
            self.model,
            &tip,
            linear((0.0, north), (tau, north), (0.0, tau), self.tol)?,
            surface,
            Location::identity(),
            (0.0, tau),
        )?;
        // Counter-clockwise in the chart: along the rim, up the seam's far
        // side, back along the pole, down the seam's near side.
        let ring = vec![
            oriented(&edges[rim], with),
            seam.clone(),
            tip.reversed(),
            seam.reversed(),
        ];
        let wire = self.model.add_wire(&ring)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, std::slice::from_ref(&wire))?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A cone closing at its apex inside its rim: the rim, a ruling seam
    /// from the apex up to the rim's vertex on the planned surface's angle
    /// zero, and the apex as an edge of no length.
    pub(super) fn cone_tip_face(
        &mut self,
        curved: &Curved,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let g = self.groups.of[rings[0][0] / 3];
        let Canonical::Cone(cone) = curved.shape else {
            ogeom_bail!(Construction, "a cone's tip is a cone's");
        };
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a cone's tip was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let frame = cone.frame();
        let v_apex = -cone.reference_radius() / cone.half_angle().tan();
        // The rim walked so the chart's angle falls from a whole turn to
        // nothing: each arc's image is the parallel's line over its angles,
        // read off the arc's own parameter.
        let mut walk = self.entries(&rings[0]);
        let first = {
            let (edge, forward) = walk[0];
            let Curve::Circle(c) = &self.plan.edges[edge].curve else {
                ogeom_bail!(Construction, "a cone's rim is not a circle");
            };
            let with = c.circle().frame().z().vector().dot(frame.z().vector()) > 0.0;
            with == forward
        };
        if first {
            walk.reverse();
            for entry in &mut walk {
                entry.1 = !entry.1;
            }
        }
        let start = {
            let (edge, forward) = walk[0];
            let spec = &self.plan.edges[edge];
            spec.curve
                .point_at(if forward { spec.range.0 } else { spec.range.1 }, self.tol)?
        };
        let v_rim = (start - frame.origin()).dot(frame.z().vector());
        let mut u = tau;
        let mut rim_ring = Vec::with_capacity(walk.len());
        for (k, &(edge, forward)) in walk.iter().enumerate() {
            let spec = &self.plan.edges[edge];
            let Curve::Circle(c) = &spec.curve else {
                ogeom_bail!(Construction, "a cone's rim is not a circle");
            };
            let with = c.circle().frame().z().vector().dot(frame.z().vector()) > 0.0;
            let (t0, t1) = spec.range;
            let span = (t1 - t0).abs();
            let end = if k + 1 == walk.len() { 0.0 } else { u - span };
            // Walked forward, the edge's start is where the walk is; the
            // angle falls along the walk.
            let (a, b) = if forward { (u, end) } else { (end, u) };
            if with == forward && walk.len() > 1 {
                ogeom_bail!(Construction, "a cone's rim arcs turn against each other");
            }
            crate::build::attach_pcurve(
                self.model,
                &edges[edge],
                linear((a, v_rim), (b, v_rim), (t0, t1), self.tol)?,
                surface,
                Location::identity(),
                (t0, t1),
            )?;
            rim_ring.push(oriented(&edges[edge], forward));
            u = end;
        }
        let rim_vertex = {
            let (edge, forward) = walk[0];
            let ends = self.model.children_of(&edges[edge])?;
            let vertex = if forward { ends.first() } else { ends.last() };
            let Some(vertex) = vertex.cloned() else {
                ogeom_bail!(Construction, "a rim has no vertex");
            };
            vertex
        };
        let at_apex = frame.origin() + frame.z().vector() * v_apex;
        let apex = self.model.add_vertex(VertexData::new(at_apex));
        let seam_range = (0.0, at_apex.distance(start));
        let id = self
            .model
            .geometry_mut()
            .add_curve(LineCurve::segment(at_apex, start, self.tol)?.into());
        let seam = self.model.add_edge(
            EdgeData::on_curve(id, Location::identity(), seam_range),
            &[apex.clone(), rim_vertex],
        )?;
        crate::build::attach_seam(
            self.model,
            &seam,
            linear((tau, v_apex), (tau, v_rim), seam_range, self.tol)?,
            linear((0.0, v_apex), (0.0, v_rim), seam_range, self.tol)?,
            surface,
            Location::identity(),
            seam_range,
        )?;
        let mut data = EdgeData::new();
        data.degenerate = true;
        let tip = self.model.add_edge(data, &[apex.clone(), apex])?;
        crate::build::attach_pcurve(
            self.model,
            &tip,
            linear((0.0, v_apex), (tau, v_apex), (0.0, tau), self.tol)?,
            surface,
            Location::identity(),
            (0.0, tau),
        )?;
        // Counter-clockwise in the chart: along the apex, up the seam's far
        // side, back along the rim, down the seam's near side.
        let mut ring = vec![tip, seam.clone()];
        ring.extend(rim_ring);
        ring.push(seam.reversed());
        let wire = self.model.add_wire(&ring)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, std::slice::from_ref(&wire))?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A sphere's cap with holes: the seam down from the pole, the rim in
    /// the mesh's order from the seam's vertex, the seam back up a whole
    /// turn over, and the pole as an edge of no length; each hole an inner
    /// wire. The seam is the meridian where it runs straight up the chart
    /// from a vertex on the sphere, and a curve traced along its chart line
    /// otherwise. The wires run as the region's triangles do, and the face
    /// flips as one where those face against the surface.
    pub(super) fn holed_cap_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let tau = core::f64::consts::TAU;
        let north = core::f64::consts::FRAC_PI_2;
        let Canonical::Sphere(sphere) = curved.shape else {
            ogeom_bail!(Construction, "a cap is a sphere's");
        };
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a cap was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
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
        let Some((rim, turn, holes)) =
            cap_rings(curved, rings, self.triangles, self.points, self.tol)
        else {
            ogeom_bail!(Construction, "a cap has one rim and holes");
        };
        let entries = self.entries(&rings[rim]);
        let mut starts: Vec<(Shape, Point)> = Vec::with_capacity(entries.len());
        for &(edge, forward) in &entries {
            let ends = self.model.children_of(&edges[edge])?;
            let vertex = if forward { ends.first() } else { ends.last() };
            let Some(vertex) = vertex.cloned() else {
                ogeom_bail!(Construction, "a rim edge has no vertex");
            };
            let Some(ogeom_topo::NodeData::Vertex(data)) =
                self.model.node(&vertex).map(|n| n.data())
            else {
                ogeom_bail!(Construction, "a rim vertex has no position");
            };
            let at = data.point;
            starts.push((vertex, at));
        }
        let hole_rings: Vec<&[Half]> = holes.iter().map(|&k| rings[k].as_slice()).collect();
        let at: Vec<Point> = starts.iter().map(|s| s.1).collect();
        let Some((i, a, b)) = cap_seam(
            curved,
            &rings[rim],
            &hole_rings,
            &at,
            self.triangles,
            self.points,
            self.tol,
        ) else {
            ogeom_bail!(Construction, "no seam reaches the pole clear of the holes");
        };
        let frame = sphere.frame();
        let pole_at = sphere.centre() + frame.z().vector() * sphere.radius();
        let pa = starts[i].1;
        let meridian = (b.0 - a.0).abs() <= 1e-12
            && pa.distance(evaluate(&curved.shape, a)) <= self.tol.confusion();
        let (seam_curve, range, deviation): (Curve, (f64, f64), f64) = if meridian {
            let x = Direction::new(
                frame.x().vector() * a.0.cos() + frame.y().vector() * a.0.sin(),
                self.tol,
            )?;
            let normal = Direction::new(x.vector().cross(frame.z().vector()), self.tol)?;
            let circle = ogeom_math::Circle::new(
                Frame::new(sphere.centre(), normal, x, self.tol)?,
                sphere.radius(),
                self.tol,
            )?;
            (
                ogeom_geom::CircleCurve::new(circle).into(),
                (a.1, north),
                self.tol.confusion(),
            )
        } else {
            chart_trace(&curved.shape, (a, b), (pa, pole_at), self.tol)?
        };
        let pole = self.model.add_vertex(VertexData::new(pole_at));
        let id = self.model.geometry_mut().add_curve(seam_curve);
        let mut data = EdgeData::on_curve(id, Location::identity(), range);
        data.tolerance = Tolerance::new(deviation)?;
        let seam = self
            .model
            .add_edge(data, &[starts[i].0.clone(), pole.clone()])?;
        // Down its near side from the pole, round the rim the way the
        // triangles run, and up its far side a whole turn over. A face
        // turned against its surface walks its wires back, and so each side
        // of the seam the other way.
        let over = tau * f64::from(turn);
        let far = linear((a.0 + over, a.1), (b.0 + over, b.1), range, self.tol)?;
        let near = linear(a, b, range, self.tol)?;
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
        let mut tip_data = EdgeData::new();
        tip_data.degenerate = true;
        let tip = self.model.add_edge(tip_data, &[pole.clone(), pole])?;
        crate::build::attach_pcurve(
            self.model,
            &tip,
            linear((b.0 + over, north), (b.0, north), (0.0, tau), self.tol)?,
            surface,
            Location::identity(),
            (0.0, tau),
        )?;
        let mut outer = vec![seam.reversed()];
        for k in 0..entries.len() {
            let (edge, forward) = entries[(i + k) % entries.len()];
            outer.push(oriented(&edges[edge], forward));
        }
        outer.push(seam);
        outer.push(tip);
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
        let mut face_data = FaceData::new(surface, Location::identity());
        face_data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(face_data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A sphere or torus whole but for holes: the whole surface's face, a
    /// torus's seam round its axis on the parallel the plan placed clear of
    /// the holes, with each ring an inner wire. The rings
    /// are turned to run as the primitive's own wires do, so the face
    /// flips as one where the region faces against its surface.
    pub(super) fn holed_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let outward = self.outward(curved, g);
        // The whole surface's wires as stored, counter-clockwise in the
        // chart whichever side the face faces.
        let whole = self.whole_face(curved, g, torus_seam_v(curved))?;
        let whole = whole.oriented(ogeom_topo::Orientation::Forward);
        let Some(ogeom_topo::NodeData::Face(data)) =
            self.model.node(&whole).map(|n| n.data().clone())
        else {
            ogeom_bail!(Construction, "a whole surface's face has no data");
        };
        let surface = data.surface;
        let mut wires = self.model.ordered_children_of(&whole)?;
        for ring in rings {
            let mut ring_edges = Vec::new();
            for (edge, forward) in self.entries(ring) {
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
                        self.plan.edges[edge].range,
                    )?;
                }
                ring_edges.push(oriented(&edges[edge], forward));
            }
            if !outward {
                ring_edges.reverse();
                ring_edges = ring_edges.iter().map(Shape::reversed).collect();
            }
            wires.push(self.model.add_wire(&ring_edges)?);
        }
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut face_data = FaceData::new(surface, Location::identity());
        face_data.tolerance = data.tolerance;
        let face = self.model.add_face(face_data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A torus whole but for holes that leave no parallel (or no meridian)
    /// free: its seam that way is the plan's [`Thread`], straight seam
    /// edges between the holes it passes through, and its seam the other
    /// way a full circle from the chain's start. In the working chart the
    /// outer wire runs along the chain (over each hole it passes through),
    /// up the circle a whole turn on, back along the chain a whole turn up
    /// (under each hole) and down the circle; every other hole is its own
    /// wire.
    ///
    /// [`Thread`]: super::Thread
    pub(super) fn threaded_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        use core::f64::consts::{PI, TAU};
        let Some(thread) = self.plan.threads.get(&g).cloned() else {
            ogeom_bail!(Construction, "a threaded face was planned without its seam");
        };
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a torus was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
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
        let swap = thread.round_tube;
        let slices: Vec<&[Half]> = rings.iter().map(Vec::as_slice).collect();
        let polygons = swapped_rings(
            centred_rings(
                curved,
                hole_polygons(
                    &curved.shape,
                    &slices,
                    self.triangles,
                    self.points,
                    self.tol,
                ),
            ),
            swap,
        );
        if polygons.len() != rings.len() {
            ogeom_bail!(Construction, "a hole has no place in its face's chart");
        }
        // Each hole's edges, clockwise in the working chart, with the
        // vertex each starts from and where that stands.
        let mut walks: Vec<Vec<(usize, bool, Shape, Point)>> = Vec::with_capacity(rings.len());
        for (ring, polygon) in rings.iter().zip(&polygons) {
            let area: f64 = (0..polygon.len())
                .map(|i| {
                    let (p, q) = (polygon[i], polygon[(i + 1) % polygon.len()]);
                    p.0 * q.1 - q.0 * p.1
                })
                .sum();
            let mut list = self.entries(ring);
            if area > 0.0 {
                list.reverse();
                for entry in &mut list {
                    entry.1 = !entry.1;
                }
            }
            let mut walk = Vec::with_capacity(list.len());
            for (edge, forward) in list {
                let ends = self.model.children_of(&edges[edge])?;
                let vertex = if forward { ends.first() } else { ends.last() };
                let Some(vertex) = vertex.cloned() else {
                    ogeom_bail!(Construction, "a hole's edge has no vertex");
                };
                let Some(ogeom_topo::NodeData::Vertex(data)) =
                    self.model.node(&vertex).map(|n| n.data())
                else {
                    ogeom_bail!(Construction, "a hole's vertex has no position");
                };
                let point = data.point;
                walk.push((edge, forward, vertex, point));
            }
            walks.push(walk);
        }
        // Where along its hole's walk the chain enters and leaves.
        let place = |walk: &[(usize, bool, Shape, Point)], v: u32| -> OgeomResult<usize> {
            let at = self.points[v as usize];
            let nearest = (0..walk.len())
                .min_by(|&i, &j| walk[i].3.distance(at).total_cmp(&walk[j].3.distance(at)));
            match nearest {
                Some(i) => Ok(i),
                None => ogeom_bail!(Construction, "a hole has no edges"),
            }
        };
        let arc = |walk: &[(usize, bool, Shape, Point)], from: usize, to: usize| -> Vec<Shape> {
            let n = walk.len();
            let count = (to + n - from) % n;
            let count = if count == 0 { n } else { count };
            (0..count)
                .map(|k| {
                    let (edge, forward, _, _) = &walk[(from + k) % n];
                    oriented(&edges[*edge], *forward)
                })
                .collect()
        };
        let chart_of = |p: (f64, f64)| swapped(p, swap);
        let start_point = evaluate(&curved.shape, chart_of(thread.start));
        let start = self.model.add_vertex(VertexData::new(start_point));
        // The chain's straight pieces, each a seam edge: walked forward on
        // its first image, back on its image a whole turn up.
        let pieces = thread_pieces(thread.start, &thread.stops);
        let mut chain: Vec<Shape> = Vec::with_capacity(pieces.len());
        let mut ends: Vec<(Shape, Point)> = vec![(start.clone(), start_point)];
        let mut uppers: Vec<Vec<Shape>> = Vec::new();
        let mut lowers: Vec<Vec<Shape>> = Vec::new();
        for stop in &thread.stops {
            let walk = &walks[stop.ring];
            let (enter, leave) = (place(walk, stop.enter)?, place(walk, stop.leave)?);
            ends.push((walk[enter].2.clone(), walk[enter].3));
            ends.push((walk[leave].2.clone(), walk[leave].3));
            uppers.push(arc(walk, enter, leave));
            lowers.push(arc(walk, leave, enter));
        }
        ends.push((start.clone(), start_point));
        let up = chart_of((0.0, TAU));
        for (k, &(a, b)) in pieces.iter().enumerate() {
            let (from, to) = (&ends[2 * k], &ends[2 * k + 1]);
            let (ca, cb) = (chart_of(a), chart_of(b));
            let (curve, range, deviation) =
                chart_trace(&curved.shape, (ca, cb), (from.1, to.1), self.tol)?;
            let id = self.model.geometry_mut().add_curve(curve);
            let mut data = EdgeData::on_curve(id, Location::identity(), range);
            data.tolerance = Tolerance::new(deviation)?;
            let piece = self.model.add_edge(data, &[from.0.clone(), to.0.clone()])?;
            let low = linear(ca, cb, range, self.tol)?;
            let high = linear(
                (ca.0 + up.0, ca.1 + up.1),
                (cb.0 + up.0, cb.1 + up.1),
                range,
                self.tol,
            )?;
            // The wire is built counter-clockwise in the working chart and
            // turned over where that chart is the surface's mirrored.
            let (forward, back) = if swap { (high, low) } else { (low, high) };
            crate::build::attach_seam(
                self.model,
                &piece,
                forward,
                back,
                surface,
                Location::identity(),
                range,
            )?;
            chain.push(piece);
        }
        // The full circle the other way, through the chain's start, its
        // parameter running with the working chart's second angle.
        let across = |t: f64| evaluate(&curved.shape, chart_of((thread.start.0, t)));
        let centre = Point::from_vector(
            (across(thread.start.1).to_vector() + across(thread.start.1 + PI).to_vector()) / 2.0,
        );
        let x = Direction::new(start_point - centre, self.tol)?;
        let onward = across(thread.start.1 + 1e-3) - across(thread.start.1 - 1e-3);
        let normal = Direction::new(x.vector().cross(onward), self.tol)?;
        let circle = ogeom_math::Circle::new(
            Frame::new(centre, normal, x, self.tol)?,
            start_point.distance(centre),
            self.tol,
        )?;
        let id = self
            .model
            .geometry_mut()
            .add_curve(ogeom_geom::CircleCurve::new(circle).into());
        let side = self.model.add_edge(
            EdgeData::on_curve(id, Location::identity(), (0.0, TAU)),
            &[start.clone(), start],
        )?;
        let s0 = thread.start;
        let near = linear(
            chart_of(s0),
            chart_of((s0.0, s0.1 + TAU)),
            (0.0, TAU),
            self.tol,
        )?;
        let far = linear(
            chart_of((s0.0 + TAU, s0.1)),
            chart_of((s0.0 + TAU, s0.1 + TAU)),
            (0.0, TAU),
            self.tol,
        )?;
        let (forward, back) = if swap { (near, far) } else { (far, near) };
        crate::build::attach_seam(
            self.model,
            &side,
            forward,
            back,
            surface,
            Location::identity(),
            (0.0, TAU),
        )?;
        let mut outer: Vec<Shape> = Vec::new();
        for (k, piece) in chain.iter().enumerate() {
            outer.push(piece.clone());
            if let Some(over) = uppers.get(k) {
                outer.extend(over.iter().cloned());
            }
        }
        outer.push(side.clone());
        for (k, piece) in chain.iter().enumerate().rev() {
            if let Some(under) = lowers.get(k) {
                outer.extend(under.iter().cloned());
            }
            outer.push(piece.reversed());
        }
        outer.push(side.reversed());
        let mut wire_lists: Vec<Vec<Shape>> = vec![outer];
        for (k, walk) in walks.iter().enumerate() {
            if thread.stops.iter().any(|stop| stop.ring == k) {
                continue;
            }
            wire_lists.push(
                walk.iter()
                    .map(|(edge, forward, _, _)| oriented(&edges[*edge], *forward))
                    .collect(),
            );
        }
        let mut wires = Vec::with_capacity(wire_lists.len());
        for mut list in wire_lists {
            if swap {
                list.reverse();
                list = list.iter().map(Shape::reversed).collect();
            }
            wires.push(self.model.add_wire(&list)?);
        }
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A whole sphere or torus: the face the primitive builds on the
    /// recognized surface, a torus's seam round its axis on the parallel at
    /// tube angle `seam_v`.
    pub(super) fn whole_face(
        &mut self,
        curved: &Curved,
        g: usize,
        seam_v: f64,
    ) -> OgeomResult<Shape> {
        let outward = self.outward(curved, g);
        let face = match curved.shape {
            Canonical::Sphere(s) => {
                let built =
                    crate::primitive::make_sphere(self.model, s.frame(), s.radius(), self.tol)?;
                let Some(face) = ogeom_topo::explore_unique(
                    self.model,
                    &built.shape,
                    ogeom_topo::ShapeType::Face,
                )?
                .into_iter()
                .next() else {
                    ogeom_bail!(Construction, "a primitive came back with no face");
                };
                face
            }
            Canonical::Torus(t) => crate::primitive::torus_face_seamed_at(
                self.model,
                t.frame(),
                t.major_radius(),
                t.minor_radius(),
                seam_v,
                self.tol,
            )?,
            _ => ogeom_bail!(Construction, "only a sphere or a torus is whole"),
        };
        // The facets stand off the surface by up to the fit's deviation, and
        // so do the seam and the poles built on it: the face's tolerance
        // reaches down to every edge and vertex it holds, as the checker
        // requires of a face looser than its boundary.
        self.model.widen(
            &face,
            Tolerance::new(curved.deviation.max(self.tol.confusion()))?,
        )?;
        Ok(if outward { face } else { face.reversed() })
    }
}
