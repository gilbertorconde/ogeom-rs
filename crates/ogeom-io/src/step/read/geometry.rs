//! Geometry: points, placements, surfaces, curves and pcurves, read from their entities.

use super::*;

impl<'a> Reader<'a> {
    // --- geometry ------------------------------------------------------------

    pub(super) fn point(&mut self, id: u64) -> OgeomResult<Point> {
        let args = self.args(id, "CARTESIAN_POINT")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            ogeom_bail!(Construction, "#{id}: a point without coordinates");
        };
        let scale = self.report.scale_mm;
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0) * scale;
        Ok(Point::new(value(0), value(1), value(2)))
    }

    fn direction(&mut self, id: u64) -> OgeomResult<Direction> {
        let args = self.args(id, "DIRECTION")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            ogeom_bail!(Construction, "#{id}: a direction without components");
        };
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0);
        Direction::new(Vector::new(value(0), value(1), value(2)), self.tol)
    }

    /// An `AXIS2_PLACEMENT_3D` as a frame, with the standard's defaults for
    /// what the file leaves out.
    pub(super) fn frame(&mut self, id: u64) -> OgeomResult<Frame> {
        let args = self.args(id, "AXIS2_PLACEMENT_3D")?;
        let origin = args
            .get(1)
            .and_then(Arg::reference)
            .map(|r| self.point(r))
            .transpose()?
            .unwrap_or(Point::ORIGIN);
        let z = args
            .get(2)
            .and_then(Arg::reference)
            .map(|r| self.direction(r))
            .transpose()?
            .unwrap_or(Direction::Z);
        let x = match args.get(3).and_then(Arg::reference) {
            Some(r) => self.direction(r)?,
            None => Direction::from_cross(z.vector(), Vector::new(0.31, 0.52, 0.8), self.tol)?,
        };
        Frame::new(origin, z, x, self.tol)
    }

    pub(super) fn surface(&mut self, id: u64) -> OgeomResult<Option<SurfaceGeometry>> {
        self.nested(|reader| reader.surface_at(id))
    }

    fn surface_at(&mut self, id: u64) -> OgeomResult<Option<SurfaceGeometry>> {
        let (keyword, args) = {
            let instance = self.instance(id)?;
            (
                instance.keyword(),
                instance
                    .parts
                    .first()
                    .map_or(&[][..], |(_, a)| a.as_slice()),
            )
        };
        let scale = self.report.scale_mm;
        let radius_arg = |args: &[Arg], i: usize| args.get(i).and_then(Arg::number);
        // B-spline surfaces arrive two ways: a simple instance with every
        // attribute in one list, or a complex instance whose parts each
        // carry their own slice, the rational form always the latter.
        {
            let (base, knots_part, weights) = {
                let instance = self.instance(id)?;
                (
                    instance.part("B_SPLINE_SURFACE"),
                    instance.part("B_SPLINE_SURFACE_WITH_KNOTS"),
                    instance.part("RATIONAL_B_SPLINE_SURFACE"),
                )
            };
            if let Some(kp) = knots_part {
                let (degrees, grid_arg, mults_knots) = if let Some(base) = base {
                    (
                        (base.first().cloned(), base.get(1).cloned()),
                        base.get(2).cloned(),
                        (
                            kp.first().cloned(),
                            kp.get(1).cloned(),
                            kp.get(2).cloned(),
                            kp.get(3).cloned(),
                        ),
                    )
                } else {
                    (
                        (kp.get(1).cloned(), kp.get(2).cloned()),
                        kp.get(3).cloned(),
                        (
                            kp.get(8).cloned(),
                            kp.get(9).cloned(),
                            kp.get(10).cloned(),
                            kp.get(11).cloned(),
                        ),
                    )
                };
                return self
                    .bspline_surface(id, degrees, grid_arg, mults_knots, weights, None)
                    .map(Some);
            }
        }
        // A complex instance whose knots are implied by its form, or a
        // simple one of those forms.
        {
            let (base, weights, form) = {
                let instance = self.instance(id)?;
                (
                    instance.part("B_SPLINE_SURFACE"),
                    instance.part("RATIONAL_B_SPLINE_SURFACE"),
                    ["BEZIER_SURFACE", "UNIFORM_SURFACE", "QUASI_UNIFORM_SURFACE"]
                        .into_iter()
                        .find(|k| instance.part(k).is_some()),
                )
            };
            if let (Some(base), Some(form)) = (base, form) {
                return self
                    .bspline_surface(
                        id,
                        (base.first().cloned(), base.get(1).cloned()),
                        base.get(2).cloned(),
                        (None, None, None, None),
                        weights,
                        Some(form),
                    )
                    .map(Some);
            }
            if let Some(form) = ["BEZIER_SURFACE", "UNIFORM_SURFACE", "QUASI_UNIFORM_SURFACE"]
                .into_iter()
                .find(|k| *k == keyword)
            {
                return self
                    .bspline_surface(
                        id,
                        (args.get(1).cloned(), args.get(2).cloned()),
                        args.get(3).cloned(),
                        (None, None, None, None),
                        None,
                        Some(form),
                    )
                    .map(Some);
            }
        }
        let out = match keyword {
            "PLANE" => {
                let frame = self.frame(ref_at(args, 1))?;
                Some(
                    PlaneSurface::over(
                        Plane::new(frame),
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                    )?
                    .into(),
                )
            }
            "CYLINDRICAL_SURFACE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let radius = radius_arg(args, 2).unwrap_or(0.0) * scale;
                Some(
                    CylinderSurface::new(
                        Cylinder::new(frame, radius, self.tol)?,
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                    )?
                    .into(),
                )
            }
            "CONICAL_SURFACE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let radius = radius_arg(args, 2).unwrap_or(0.0) * scale;
                let angle = radius_arg(args, 3).unwrap_or(0.0) * self.angle_scale;
                Some(
                    ConeSurface::new(
                        Cone::new(frame, radius, angle, self.tol)?,
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                    )?
                    .into(),
                )
            }
            "SPHERICAL_SURFACE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let radius = radius_arg(args, 2).unwrap_or(0.0) * scale;
                Some(SphereSurface::new(Sphere::new(frame, radius, self.tol)?).into())
            }
            "TOROIDAL_SURFACE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let major = radius_arg(args, 2).unwrap_or(0.0) * scale;
                let minor = radius_arg(args, 3).unwrap_or(0.0) * scale;
                Some(TorusSurface::new(Torus::new(frame, major, minor, self.tol)?).into())
            }
            "SURFACE_OF_LINEAR_EXTRUSION" => {
                // A curve swept along a vector, unbounded either way. A
                // file writes a drum's wall this way as often as it writes a
                // cylinder (a circle swept along its own axis) and a wall
                // as a line swept: those are read as the cylinder and the
                // plane they are, exact and known everywhere downstream.
                // Anything else sweeps as itself, over a window the face's
                // own edges then widen to fit.
                let Some(curve) = self.curve(ref_at(args, 1))? else {
                    self.report.warnings.push(format!(
                        "#{id}: an extrusion's swept curve is not read; its face is skipped"
                    ));
                    return Ok(None);
                };
                let vector = self.args(ref_at(args, 2), "VECTOR")?;
                let direction = self.direction(vector[1].reference().unwrap_or(0))?;
                // Parallel to the file's own precision in directions: a
                // writer states an axis to nine digits, a whisker off the
                // sweep it meant to be exactly along.
                let parallel =
                    |axis: Direction| axis.vector().cross(direction.vector()).magnitude() <= 1e-6;
                match &curve {
                    Curve::Circle(c) if parallel(c.circle().frame().z()) => Some(
                        CylinderSurface::new(
                            Cylinder::new(c.circle().frame(), c.circle().radius(), self.tol)?,
                            (-SURFACE_EXTENT, SURFACE_EXTENT),
                        )?
                        .into(),
                    ),
                    Curve::Line(l) if !parallel(l.axis().direction) => {
                        let axis = l.axis();
                        let normal = Direction::from_cross(
                            axis.direction.vector(),
                            direction.vector(),
                            self.tol,
                        )?;
                        let frame = Frame::new(axis.location, normal, axis.direction, self.tol)?;
                        Some(
                            PlaneSurface::over(
                                Plane::new(frame),
                                (-SURFACE_EXTENT, SURFACE_EXTENT),
                                (-SURFACE_EXTENT, SURFACE_EXTENT),
                            )?
                            .into(),
                        )
                    }
                    _ => Some(
                        ExtrusionSurface::over(
                            curve,
                            direction,
                            (-SURFACE_EXTENT, SURFACE_EXTENT),
                        )?
                        .into(),
                    ),
                }
            }
            "SURFACE_OF_REVOLUTION" => {
                let Some(curve) = self.curve(ref_at(args, 1))? else {
                    self.report.warnings.push(format!(
                        "#{id}: a revolution's swept curve is not read; its face is skipped"
                    ));
                    return Ok(None);
                };
                let placement = self.args(ref_at(args, 2), "AXIS1_PLACEMENT")?;
                let location = self.point(placement[1].reference().unwrap_or(0))?;
                let direction = match placement.get(2).and_then(Arg::reference) {
                    Some(r) => self.direction(r)?,
                    None => Direction::Z,
                };
                Some(
                    ogeom_geom::RevolutionSurface::new(
                        curve,
                        Axis {
                            location,
                            direction,
                        },
                        core::f64::consts::TAU,
                    )?
                    .into(),
                )
            }
            "OFFSET_SURFACE" => {
                let Some(basis) = self.surface(ref_at(args, 1))? else {
                    return Ok(None);
                };
                let distance = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                if distance == 0.0 {
                    return Ok(Some(basis));
                }
                // An analytic basis offsets to the analytic surface it is,
                // which every downstream path speaks exactly.
                let offset = ogeom_geom::OffsetSurface::new(basis, distance)?;
                Some(match offset.analytic(self.tol)? {
                    Some(analytic) => analytic,
                    None => SurfaceGeometry::Offset(Box::new(offset)),
                })
            }
            // A face's own edges bound it; the window these name is the
            // basis's, restated.
            "RECTANGULAR_TRIMMED_SURFACE" | "CURVE_BOUNDED_SURFACE" => {
                self.surface(ref_at(args, 1))?
            }
            "DEGENERATE_TOROIDAL_SURFACE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let major = radius_arg(args, 2).unwrap_or(0.0) * scale;
                let minor = radius_arg(args, 3).unwrap_or(0.0) * scale;
                Some(TorusSurface::new(Torus::new(frame, major, minor, self.tol)?).into())
            }
            "RECTANGULAR_COMPOSITE_SURFACE" => Some(self.composite_surface(id, args)?),
            "SURFACE_REPLICA" => {
                let Some(parent) = self.surface(ref_at(args, 1))? else {
                    return Ok(None);
                };
                let motion = self.transformation_operator(ref_at(args, 2))?;
                Some(parent.transformed(&motion, self.tol)?)
            }
            other => {
                self.report.warnings.push(format!(
                    "#{id}: surface kind {other} is not read yet; its face is skipped"
                ));
                None
            }
        };
        Ok(out)
    }

    /// Expand STEP's multiplicity-compressed knots.
    fn expand_knots(mults: Option<Arg>, knots: Option<Arg>) -> OgeomResult<Vec<f64>> {
        let (Some(Arg::List(mults)), Some(Arg::List(knots))) = (mults, knots) else {
            return Ok(Vec::new());
        };
        expanded_knots(&mults, &knots)
    }

    #[allow(clippy::type_complexity)]
    fn bspline_surface(
        &mut self,
        id: u64,
        degrees: (Option<Arg>, Option<Arg>),
        grid_arg: Option<Arg>,
        mults_knots: (Option<Arg>, Option<Arg>, Option<Arg>, Option<Arg>),
        weights: Option<&[Arg]>,
        form: Option<&str>,
    ) -> OgeomResult<SurfaceGeometry> {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (u_degree, v_degree) = (
            file_degree(degrees.0.as_ref(), 1.0)?,
            file_degree(degrees.1.as_ref(), 1.0)?,
        );
        let Some(Arg::List(rows)) = grid_arg else {
            ogeom_bail!(
                Construction,
                "#{id}: a b-spline surface without control points"
            );
        };
        let mut points = Vec::new();
        let (mut u_count, mut v_count) = (0, 0);
        for row in &rows {
            let Some(cells) = row.list() else { continue };
            u_count += 1;
            v_count = cells.len();
            for cell in cells {
                points.push(self.point(cell.reference().unwrap_or(0))?);
            }
        }
        let (u_raw, v_raw) = match form {
            Some(form) => (
                implied_knots(form, u_degree, u_count),
                implied_knots(form, v_degree, v_count),
            ),
            None => (
                Self::expand_knots(mults_knots.0, mults_knots.2)?,
                Self::expand_knots(mults_knots.1, mults_knots.3)?,
            ),
        };
        let u_knots = KnotVector::new(u_raw, u_degree)?;
        let v_knots = KnotVector::new(v_raw, v_degree)?;
        let surface = if let Some(weights) = weights {
            let flat: Vec<f64> = weights
                .first()
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::list)
                .flatten()
                .filter_map(Arg::number)
                .collect();
            let weighted: Vec<ogeom_math::Weighted<Point>> = points
                .iter()
                .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                .map(|(p, w)| ogeom_math::Weighted::new(*p, *w, self.tol))
                .collect::<OgeomResult<_>>()?;
            ogeom_geom::BSplineSurface::rational(
                u_knots,
                v_knots,
                ogeom_math::ControlGrid::new(weighted, u_count, v_count)?,
            )?
        } else {
            ogeom_geom::BSplineSurface::new(
                u_knots,
                v_knots,
                &ogeom_math::ControlGrid::new(points, u_count, v_count)?,
                self.tol,
            )?
        };
        Ok(surface.into())
    }

    pub(super) fn curve(&mut self, id: u64) -> OgeomResult<Option<Curve>> {
        self.nested(|reader| reader.curve_at(id))
    }

    fn curve_at(&mut self, id: u64) -> OgeomResult<Option<Curve>> {
        let (keyword, args) = {
            let instance = self.instance(id)?;
            (
                instance.keyword(),
                instance
                    .parts
                    .first()
                    .map_or(&[][..], |(_, a)| a.as_slice()),
            )
        };
        let scale = self.report.scale_mm;
        {
            let (base, kp, weights) = {
                let instance = self.instance(id)?;
                (
                    instance.part("B_SPLINE_CURVE"),
                    instance.part("B_SPLINE_CURVE_WITH_KNOTS"),
                    instance.part("RATIONAL_B_SPLINE_CURVE"),
                )
            };
            // A complex instance: the base part carries degree and control
            // points; the knots come from the knots part or are implied by
            // the form (Bézier, uniform, quasi-uniform); weights from the
            // rational part where there is one.
            let form = {
                let instance = self.instance(id)?;
                ["BEZIER_CURVE", "UNIFORM_CURVE", "QUASI_UNIFORM_CURVE"]
                    .into_iter()
                    .find(|k| instance.part(k).is_some())
            };
            if let Some(base) = base
                && (kp.is_some() || form.is_some())
            {
                let degree = file_degree(base.first(), 1.0)?;
                let control: Vec<Point> = base
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .map(|r| self.point(r))
                    .collect::<OgeomResult<_>>()?;
                let knots = match (&kp, form) {
                    (Some(kp), _) => Self::expand_knots(kp.first().cloned(), kp.get(1).cloned())?,
                    (None, Some(form)) => implied_knots(form, degree, control.len()),
                    (None, None) => unreachable!("guarded above"),
                };
                let knots = KnotVector::new(knots, degree)?;
                let flat: Vec<f64> = weights
                    .as_ref()
                    .and_then(|w| w.first())
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::number)
                    .collect();
                let weighted: Vec<ogeom_math::Weighted<Point>> = control
                    .iter()
                    .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                    .map(|(p, w)| ogeom_math::Weighted::new(*p, *w, self.tol))
                    .collect::<OgeomResult<_>>()?;
                return Ok(Some(BSplineCurve::rational(knots, weighted)?.into()));
            }
        }
        let out: Option<Curve> = match keyword {
            "LINE" => {
                let through = self.point(ref_at(args, 1))?;
                // The vector's magnitude scales STEP's parameter; ranges here
                // are re-derived from vertex geometry, so only the direction
                // matters.
                let vector = self.args(ref_at(args, 2), "VECTOR")?;
                let direction = self.direction(vector[1].reference().unwrap_or(0))?;
                Some(
                    LineCurve::new(Axis {
                        location: through,
                        direction,
                    })
                    .into(),
                )
            }
            "CIRCLE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let radius = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(CircleCurve::new(Circle::new(frame, radius, self.tol)?).into())
            }
            "ELLIPSE" => {
                let frame = self.frame(ref_at(args, 1))?;
                let a = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                let b = args.get(3).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(EllipseCurve::new(Ellipse::new(frame, a, b, self.tol)?).into())
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = file_degree(args.get(1), 0.0)?;
                let control: Vec<Point> = args
                    .get(2)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .map(|r| self.point(r))
                    .collect::<OgeomResult<_>>()?;
                let mults = args.get(6).and_then(Arg::list).unwrap_or(&[]).to_vec();
                let knots = args.get(7).and_then(Arg::list).unwrap_or(&[]).to_vec();
                let expanded = expanded_knots(&mults, &knots)?;
                Some(
                    BSplineCurve::new(KnotVector::new(expanded, degree)?, control, self.tol)?
                        .into(),
                )
            }
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" => {
                // A curve dressed in its surface associations: what many
                // exporters write for every edge. The 3D curve is the first
                // argument; the pcurves are read where a face asks for them.
                self.curve(args.get(1).and_then(Arg::reference).unwrap_or(0))?
            }
            "BEZIER_CURVE" | "UNIFORM_CURVE" | "QUASI_UNIFORM_CURVE" => {
                let degree = file_degree(args.get(1), 1.0)?;
                let control: Vec<Point> = args
                    .get(2)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .map(|r| self.point(r))
                    .collect::<OgeomResult<_>>()?;
                let knots = implied_knots(keyword, degree, control.len());
                Some(BSplineCurve::new(KnotVector::new(knots, degree)?, control, self.tol)?.into())
            }
            "HYPERBOLA" => {
                let frame = self.frame(ref_at(args, 1))?;
                let a = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                let b = args.get(3).and_then(Arg::number).unwrap_or(0.0) * scale;
                // The branch's reach in its own parameter: cosh 20 carries a
                // hyperbola far past any part.
                Some(HyperbolaCurve::new(Hyperbola::new(frame, a, b, self.tol)?, 20.0)?.into())
            }
            "PARABOLA" => {
                let frame = self.frame(ref_at(args, 1))?;
                let focal = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(
                    ParabolaCurve::new(Parabola::new(frame, focal, self.tol)?, SURFACE_EXTENT)?
                        .into(),
                )
            }
            "TRIMMED_CURVE" => {
                // An edge's ends are its vertices; the trims restate them.
                self.curve(args.get(1).and_then(Arg::reference).unwrap_or(0))?
            }
            "OFFSET_CURVE_3D" => {
                let Some(basis) = self.curve(args.get(1).and_then(Arg::reference).unwrap_or(0))?
                else {
                    return Ok(None);
                };
                let distance = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                let reference =
                    self.direction(args.get(4).and_then(Arg::reference).unwrap_or(0))?;
                Some(Curve::Offset(Box::new(ogeom_geom::OffsetCurve::new(
                    basis, distance, reference,
                )?)))
            }
            "POLYLINE" => {
                let points: Vec<Point> = args
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .map(|r| self.point(r))
                    .collect::<OgeomResult<_>>()?;
                if points.len() < 2 {
                    ogeom_bail!(Construction, "#{id}: a polyline of fewer than two points");
                }
                // Degree one through the points, a knot at each.
                let mut knots = vec![0.0];
                #[allow(clippy::cast_precision_loss)]
                knots.extend((0..points.len()).map(|i| i as f64));
                #[allow(clippy::cast_precision_loss)]
                knots.push((points.len() - 1) as f64);
                Some(BSplineCurve::new(KnotVector::new(knots, 1)?, points, self.tol)?.into())
            }
            "COMPOSITE_CURVE" => {
                // Its segments, each over its own trim and in its own
                // sense, as splines joined end to start.
                let mut joined: Option<(KnotVector, Vec<ogeom_math::Weighted<Point>>)> = None;
                let segments: Vec<u64> = args
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .collect();
                for segment in segments {
                    let sargs = self.args(segment, "COMPOSITE_CURVE_SEGMENT")?;
                    let same_sense = !sargs.get(1).is_some_and(|a| a.is_enum("F"));
                    let parent = sargs.get(2).and_then(Arg::reference).unwrap_or(0);
                    let Some((curve, range)) = self.bounded_curve(parent)? else {
                        return Ok(None);
                    };
                    let mut spline = curve.to_bspline_over(range, self.tol)?;
                    if !same_sense {
                        let (knots, control) =
                            ogeom_math::bspline::reverse(spline.knots(), spline.control_points());
                        spline = BSplineCurve::rational(knots, control)?;
                    }
                    joined = Some(match joined {
                        None => (spline.knots().clone(), spline.control_points().to_vec()),
                        Some(held) => join_splines(held, &spline, self.tol)?,
                    });
                }
                let Some((knots, control)) = joined else {
                    ogeom_bail!(Construction, "#{id}: a composite curve with no segments");
                };
                Some(BSplineCurve::rational(knots, control)?.into())
            }
            "CURVE_REPLICA" => {
                let Some(parent) = self.curve(args.get(1).and_then(Arg::reference).unwrap_or(0))?
                else {
                    return Ok(None);
                };
                let motion = self
                    .transformation_operator(args.get(2).and_then(Arg::reference).unwrap_or(0))?;
                Some(parent.transformed(&motion, self.tol)?)
            }
            other => {
                self.report.warnings.push(format!(
                    "#{id}: curve kind {other} is not read yet; its edge is skipped"
                ));
                None
            }
        };
        Ok(out)
    }

    /// A rectangular composite surface: its grid of patches joined into one
    /// B-spline surface, so a face across the grid is a face on one surface
    /// and needs no rebuilding. Each patch in its spline form over its own
    /// bounds, turned where its sense says, raised to the grid's degrees,
    /// its knots unified with its column's and its row's, and joined to its
    /// neighbours on their shared boundary. Neighbours that do not meet
    /// there are refused by name.
    fn composite_surface(&mut self, id: u64, args: &[Arg]) -> OgeomResult<SurfaceGeometry> {
        let rows: Vec<Vec<u64>> = args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .map(|row| {
                row.list()
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .collect()
            })
            .collect();
        if rows.is_empty()
            || rows
                .iter()
                .any(|r| r.len() != rows[0].len() || r.is_empty())
        {
            ogeom_bail!(
                Construction,
                "#{id}: a composite surface's grid is not rectangular"
            );
        }
        let mut grid: Vec<Vec<Patch>> = Vec::with_capacity(rows.len());
        for row in &rows {
            let mut out = Vec::with_capacity(row.len());
            for &patch in row {
                let pargs = self.args(patch, "SURFACE_PATCH")?;
                let parent = pargs.get(1).and_then(Arg::reference).unwrap_or(0);
                let u_sense = !pargs.get(4).is_some_and(|a| a.is_enum("F"));
                let v_sense = !pargs.get(5).is_some_and(|a| a.is_enum("F"));
                let spline = self.bounded_surface_spline(parent)?;
                let mut p = Patch::of(&spline)?;
                if !u_sense {
                    p = p.reversed_u();
                }
                if !v_sense {
                    p = p.transposed().reversed_u().transposed();
                }
                out.push(p);
            }
            grid.push(out);
        }
        join_patches(grid, self.tol).map_err(|e| {
            ogeom_core::ogeom_err!(
                Construction,
                "#{id}: a composite surface's patches do not join: {e}"
            )
        })
    }

    /// A composite's patch surface as a B-spline over its own bounds: a
    /// B-spline as it is; a rectangularly trimmed surface's basis over its
    /// window, in the file's reading of the basis's parameters.
    fn bounded_surface_spline(&mut self, id: u64) -> OgeomResult<ogeom_geom::BSplineSurface> {
        let keyword = self.instance(id)?.keyword().to_owned();
        if keyword == "RECTANGULAR_TRIMMED_SURFACE" {
            let args = self.args(id, "RECTANGULAR_TRIMMED_SURFACE")?;
            let basis_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
            let Some(basis) = self.surface(basis_id)? else {
                ogeom_bail!(Construction, "#{id}: a patch's basis is not read");
            };
            let angular = matches!(
                self.instance(basis_id)?.keyword(),
                "CYLINDRICAL_SURFACE"
                    | "CONICAL_SURFACE"
                    | "SPHERICAL_SURFACE"
                    | "TOROIDAL_SURFACE"
            );
            let n = |i: usize| args.get(i).and_then(Arg::number).unwrap_or(0.0);
            let (u_scale, v_scale) = if angular {
                (self.angle_scale, self.report.scale_mm)
            } else {
                (self.report.scale_mm, self.report.scale_mm)
            };
            let (u1, u2) = (n(2) * u_scale, n(3) * u_scale);
            let (v1, v2) = (n(4) * v_scale, n(5) * v_scale);
            let window = ((u1.min(u2), u1.max(u2)), (v1.min(v2), v1.max(v2)));
            let trimmed = ogeom_geom::TrimmedSurface::new(basis, window.0, window.1, self.tol)?;
            return SurfaceGeometry::Trimmed(Box::new(trimmed)).to_bspline(self.tol);
        }
        let Some(surface) = self.surface(id)? else {
            ogeom_bail!(Construction, "#{id}: a patch's surface is not read");
        };
        surface.to_bspline(self.tol)
    }

    /// A curve with the stretch of it an entity bounds: a trimmed curve's
    /// basis between its trims (a point where the file gives one, else its
    /// parameter in the file's own reading of the basis), in its own sense;
    /// any other curve over its own domain, which must be bounded.
    pub(super) fn bounded_curve(&mut self, id: u64) -> OgeomResult<Option<(Curve, (f64, f64))>> {
        self.nested(|reader| reader.bounded_curve_at(id))
    }

    fn bounded_curve_at(&mut self, id: u64) -> OgeomResult<Option<(Curve, (f64, f64))>> {
        let keyword = self.instance(id)?.keyword().to_owned();
        if keyword != "TRIMMED_CURVE" {
            let Some(curve) = self.curve(id)? else {
                return Ok(None);
            };
            let (lo, hi) = curve.domain();
            if !lo.is_finite() || !hi.is_finite() || hi - lo > SURFACE_EXTENT {
                ogeom_bail!(
                    Construction,
                    "#{id}: an unbounded curve cannot stand as a composite's segment"
                );
            }
            return Ok(Some((curve, (lo, hi))));
        }
        let args = self.args(id, "TRIMMED_CURVE")?;
        let basis_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let Some(basis) = self.curve(basis_id)? else {
            return Ok(None);
        };
        let basis_keyword = self.instance(basis_id)?.keyword().to_owned();
        let forward = !args.get(4).is_some_and(|a| a.is_enum("F"));
        let mut ends = Vec::with_capacity(2);
        for trim in [args.get(2), args.get(3)] {
            let items = trim.and_then(Arg::list).unwrap_or(&[]).to_vec();
            let mut at: Option<Point> = None;
            for item in &items {
                if let Some(r) = item.reference() {
                    at = Some(self.point(r)?);
                }
            }
            if at.is_none() {
                for item in &items {
                    if let Arg::Typed(typed) = item
                        && typed.0 == "PARAMETER_VALUE"
                        && let Some(v) = typed.1.first().and_then(Arg::number)
                    {
                        at =
                            Some(self.step_parameter_point(&basis, &basis_keyword, basis_id, v)?);
                    }
                }
            }
            let Some(point) = at else {
                ogeom_bail!(
                    Construction,
                    "#{id}: a trim names neither point nor parameter"
                );
            };
            let t = match self.parameter_of(&basis, point) {
                Some(t) => t,
                None => project_on_curve(&basis, point, 256, self.tol)?.parameter,
            };
            ends.push(t);
        }
        let (mut a, mut b) = (ends[0], ends[1]);
        if !forward {
            core::mem::swap(&mut a, &mut b);
        }
        if basis.is_periodic() && b <= a {
            let (lo, hi) = basis.domain();
            b += hi - lo;
        }
        let range = if a <= b { (a, b) } else { (b, a) };
        let mut curve = basis;
        if !forward || a > b {
            // Travelled backward: the same points, the parameter turned.
            let spline = curve.to_bspline_over(range, self.tol)?;
            let (knots, control) =
                ogeom_math::bspline::reverse(spline.knots(), spline.control_points());
            curve = BSplineCurve::rational(knots, control)?.into();
            let domain = curve.domain();
            return Ok(Some((curve, domain)));
        }
        Ok(Some((curve, range)))
    }

    /// Where the file's parameter `v` stands on a basis curve, read as the
    /// file reads it: a line's by its own vector's length, a conic's as an
    /// angle in the file's angle unit, a spline's as it is.
    fn step_parameter_point(
        &mut self,
        basis: &Curve,
        keyword: &str,
        basis_id: u64,
        v: f64,
    ) -> OgeomResult<Point> {
        match keyword {
            "LINE" => {
                let args = self.args(basis_id, "LINE")?;
                let through = self.point(ref_at(args, 1))?;
                let vector = self.args(ref_at(args, 2), "VECTOR")?;
                let direction = self.direction(vector[1].reference().unwrap_or(0))?;
                let magnitude = vector.get(2).and_then(Arg::number).unwrap_or(1.0);
                Ok(through + direction.vector() * (v * magnitude * self.report.scale_mm))
            }
            "CIRCLE" | "ELLIPSE" => basis.point_at(v * self.angle_scale, self.tol),
            _ => basis.point_at(v, self.tol),
        }
    }

    /// A Cartesian transformation operator as the motion it states: the
    /// axes it names (completed square to one another), its origin, and
    /// its uniform scale.
    fn transformation_operator(&mut self, id: u64) -> OgeomResult<Transform> {
        let args = self.args(id, "CARTESIAN_TRANSFORMATION_OPERATOR_3D")?;
        let direction_at = |this: &mut Self, i: usize| -> OgeomResult<Option<Vector>> {
            match args.get(i).and_then(Arg::reference) {
                Some(r) => Ok(Some(this.direction(r)?.vector())),
                None => Ok(None),
            }
        };
        let (axis1, axis2, axis3) = (
            direction_at(self, 1)?,
            direction_at(self, 2)?,
            direction_at(self, 5)?,
        );
        let origin = self.point(args.get(3).and_then(Arg::reference).unwrap_or(0))?;
        let scale = args.get(4).and_then(Arg::number).unwrap_or(1.0);
        let z = match (axis3, axis1, axis2) {
            (Some(z), _, _) => z,
            (None, Some(x), Some(y)) => x.cross(y),
            _ => Vector::Z,
        };
        let z = Direction::new(z, self.tol)?.vector();
        let x = match axis1 {
            Some(x) => x - z * x.dot(z),
            None => {
                let seed = if z.x.abs() < 0.9 {
                    Vector::X
                } else {
                    Vector::Y
                };
                seed - z * seed.dot(z)
            }
        };
        let x = Direction::new(x, self.tol)?.vector();
        let y = z.cross(x);
        let linear = Matrix3::new([[x.x, y.x, z.x], [x.y, y.y, z.y], [x.z, y.z, z.z]]);
        Transform::from_parts(linear, scale, origin.to_vector(), self.tol.angular())
    }

    /// The parameter of a point on one of this kernel's curves.
    pub(super) fn parameter_of(&self, curve: &Curve, p: Point) -> Option<f64> {
        crate::inversion::parameter_on(curve, p)
    }

    /// The pcurves the file states for edge `edge_id` on surface entity
    /// `surface_id`: the `PCURVE`s of its `SURFACE_CURVE` or `SEAM_CURVE`
    /// drawn on that surface, as written. Nothing where the file states
    /// none, or states them in a form not read here.
    pub(super) fn file_pcurves(&mut self, edge_id: u64, surface_id: u64) -> Vec<PlanarCurve> {
        let mut out = Vec::new();
        let Ok(args) = self.args(edge_id, "EDGE_CURVE") else {
            return out;
        };
        let curve_id = ref_at(args, 3);
        let Ok(instance) = self.instance(curve_id) else {
            return out;
        };
        let Some(dressed) = instance
            .part("SURFACE_CURVE")
            .or_else(|| instance.part("SEAM_CURVE"))
            .or_else(|| instance.part("INTERSECTION_CURVE"))
        else {
            return out;
        };
        let uses: Vec<u64> = dressed
            .get(2)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        for id in uses {
            let Ok(pargs) = self.args(id, "PCURVE") else {
                continue;
            };
            if pargs.get(1).and_then(Arg::reference) != Some(surface_id) {
                continue;
            }
            let Ok(dargs) = self.args(ref_at(pargs, 2), "DEFINITIONAL_REPRESENTATION") else {
                continue;
            };
            let Some(item) = dargs
                .get(1)
                .and_then(Arg::list)
                .and_then(|items| items.first())
                .and_then(Arg::reference)
            else {
                continue;
            };
            if let Ok(Some(curve)) = self.nested(|reader| reader.planar_curve(item)) {
                out.push(curve);
            }
        }
        // Restated in this kernel's units where the file's are known; a
        // pcurve that has no image here at the same parameters stays as
        // written.
        if out.is_empty() {
            return out;
        }
        if let (Some(chart), Some(pace)) =
            (self.chart_units(surface_id), self.curve_units(curve_id))
        {
            for pcurve in &mut out {
                if let Some(restated) =
                    in_kernel_units(pcurve, chart, pace, self.angle_scale, self.tol)
                {
                    *pcurve = restated;
                }
            }
        } else if let (Some(map), Some(pace)) = (
            self.swept_line_chart(surface_id),
            self.curve_units(curve_id),
        ) {
            for pcurve in &mut out {
                if let Some(restated) = through_map(pcurve, map, pace, self.tol) {
                    *pcurve = restated;
                }
            }
        }
        out
    }

    /// How a line swept into a plane reads the file's chart on the plane it
    /// is read as: the map from the file's parameters, the line's own and
    /// the length along the sweep's vector as ISO 10303-42 parameterizes
    /// the extrusion, to this kernel's, the plane's coordinates along its
    /// frame. The sweep may lean on the line, so the map is affine, not a
    /// stretch along each direction. `None` for any other surface.
    fn swept_line_chart(&mut self, id: u64) -> Option<ChartMap> {
        let args = self.args(id, "SURFACE_OF_LINEAR_EXTRUSION").ok()?.to_vec();
        let Some(SurfaceGeometry::Plane(plane)) = self.surface(id).ok()? else {
            return None;
        };
        let line_id = ref_at(&args, 1);
        let Some(Curve::Line(line)) = self.curve(line_id).ok()? else {
            return None;
        };
        let pace = self.curve_units(line_id)?;
        let vector = self.args(ref_at(&args, 2), "VECTOR").ok()?.to_vec();
        let across = vector.get(2).and_then(Arg::number).unwrap_or(1.0) * self.report.scale_mm;
        let sweep = self.direction(vector.get(1)?.reference()?).ok()?.vector();
        let frame = plane.plane().frame();
        let (x, y) = (frame.x().vector(), frame.y().vector());
        let (along, from) = (
            line.axis().direction.vector(),
            line.axis().location - frame.origin(),
        );
        Some(ChartMap {
            linear: [
                [pace * along.dot(x), across * sweep.dot(x)],
                [pace * along.dot(y), across * sweep.dot(y)],
            ],
            offset: (from.dot(x), from.dot(y)),
        })
    }

    /// How a surface's parameters here follow the file's: this kernel's `u`
    /// and `v` per file unit of each, as ISO 10303-42 parameterizes the
    /// surface. A plane's are lengths; a cylinder's and a cone's `u` is an
    /// angle and `v` a length along the axis; a sphere's and a torus's are
    /// both angles; a spline's are its own; a revolution's are the angle
    /// round its axis and its curve's parameter; an extrusion's are its
    /// curve's parameter and a length along its vector; a replica's are its
    /// parent's, each that is a length scaled by the replica's scale (its
    /// motion a rotation and a shift, never a reflection). `None` for a
    /// surface whose chart here is not the file's stretched along each
    /// direction (a line swept into a plane, whose sweep may lean on it:
    /// [`swept_line_chart`](Self::swept_line_chart) maps its chart).
    fn chart_units(&mut self, id: u64) -> Option<(f64, f64)> {
        let (keyword, args) = {
            let instance = self.instance(id).ok()?;
            if [
                "B_SPLINE_SURFACE",
                "B_SPLINE_SURFACE_WITH_KNOTS",
                "BEZIER_SURFACE",
                "UNIFORM_SURFACE",
                "QUASI_UNIFORM_SURFACE",
            ]
            .into_iter()
            .any(|k| instance.part(k).is_some())
            {
                return Some((1.0, 1.0));
            }
            (
                instance.keyword(),
                instance
                    .parts
                    .first()
                    .map_or(&[][..], |(_, a)| a.as_slice()),
            )
        };
        let (length, angle) = (self.report.scale_mm, self.angle_scale);
        match keyword {
            "PLANE" => Some((length, length)),
            "CYLINDRICAL_SURFACE" | "CONICAL_SURFACE" => Some((angle, length)),
            "SPHERICAL_SURFACE" | "TOROIDAL_SURFACE" | "DEGENERATE_TOROIDAL_SURFACE" => {
                Some((angle, angle))
            }
            "OFFSET_SURFACE" | "RECTANGULAR_TRIMMED_SURFACE" | "CURVE_BOUNDED_SURFACE" => {
                self.chart_units(ref_at(args, 1))
            }
            "SURFACE_OF_REVOLUTION" => Some((angle, self.curve_units(ref_at(args, 1))?)),
            // A replica's parameters are its parent's; carried by a motion
            // and a scaling, each of them that is a length here scales.
            "SURFACE_REPLICA" => {
                let parent = ref_at(args, 1);
                let scale = self
                    .args(ref_at(args, 2), "CARTESIAN_TRANSFORMATION_OPERATOR_3D")
                    .ok()?
                    .get(4)
                    .and_then(Arg::number)
                    .unwrap_or(1.0);
                if !(scale > 0.0 && scale.is_finite()) {
                    return None;
                }
                let (u, v) = self.chart_units(parent)?;
                let (lu, lv) = self.length_parameters(parent)?;
                Some((
                    if lu { u * scale } else { u },
                    if lv { v * scale } else { v },
                ))
            }
            "SURFACE_OF_LINEAR_EXTRUSION" => {
                let along = self.curve_units(ref_at(args, 1))?;
                let vector = self.args(ref_at(args, 2), "VECTOR").ok()?;
                let magnitude = vector.get(2).and_then(Arg::number).unwrap_or(1.0);
                let direction = self.direction(vector[1].reference()?).ok()?;
                let across = magnitude * length;
                // The reader takes a circle swept along its axis for the
                // cylinder it is, whose heights run along the circle's
                // axis, with the sweep or against it.
                match self.surface(id).ok()?? {
                    SurfaceGeometry::Cylinder(c) => {
                        let z = c.cylinder().frame().z();
                        Some((along, across.copysign(direction.vector().dot(z.vector()))))
                    }
                    SurfaceGeometry::Extrusion(_) => Some((along, across)),
                    _ => None,
                }
            }
            _ => None,
        }
    }

    /// Whether each of a surface's parameters here is a length, which a
    /// scaling of the surface scales: a plane's both, a cylinder's and a
    /// cone's height, a revolution's curve parameter and an extrusion's
    /// where it is the length along a line or a parabola's axis normal, an
    /// extrusion's sweep; an angle and a spline's own parameter are not.
    /// `None` for a surface [`chart_units`](Self::chart_units) does not
    /// read.
    fn length_parameters(&mut self, id: u64) -> Option<(bool, bool)> {
        let (keyword, args) = {
            let instance = self.instance(id).ok()?;
            if [
                "B_SPLINE_SURFACE",
                "B_SPLINE_SURFACE_WITH_KNOTS",
                "BEZIER_SURFACE",
                "UNIFORM_SURFACE",
                "QUASI_UNIFORM_SURFACE",
            ]
            .into_iter()
            .any(|k| instance.part(k).is_some())
            {
                return Some((false, false));
            }
            (
                instance.keyword().to_owned(),
                instance
                    .parts
                    .first()
                    .map_or_else(Vec::new, |(_, a)| a.clone()),
            )
        };
        match keyword.as_str() {
            "PLANE" => Some((true, true)),
            "CYLINDRICAL_SURFACE" | "CONICAL_SURFACE" => Some((false, true)),
            "SPHERICAL_SURFACE" | "TOROIDAL_SURFACE" | "DEGENERATE_TOROIDAL_SURFACE" => {
                Some((false, false))
            }
            "OFFSET_SURFACE"
            | "RECTANGULAR_TRIMMED_SURFACE"
            | "CURVE_BOUNDED_SURFACE"
            | "SURFACE_REPLICA" => self.length_parameters(ref_at(&args, 1)),
            "SURFACE_OF_REVOLUTION" => Some((false, self.length_parameter(ref_at(&args, 1))?)),
            "SURFACE_OF_LINEAR_EXTRUSION" => match self.surface(id).ok()?? {
                SurfaceGeometry::Cylinder(_) => Some((false, true)),
                SurfaceGeometry::Extrusion(_) => {
                    Some((self.length_parameter(ref_at(&args, 1))?, true))
                }
                _ => None,
            },
            _ => None,
        }
    }

    /// Whether a curve's parameter here is a length, which a scaling of the
    /// curve scales: a line's and a parabola's; a conic's angle, a
    /// hyperbola's and a spline's own parameter are not. `None` for a curve
    /// [`curve_units`](Self::curve_units) does not read.
    fn length_parameter(&mut self, id: u64) -> Option<bool> {
        let (keyword, inner) = {
            let instance = self.instance(id).ok()?;
            if [
                "B_SPLINE_CURVE",
                "B_SPLINE_CURVE_WITH_KNOTS",
                "BEZIER_CURVE",
                "UNIFORM_CURVE",
                "QUASI_UNIFORM_CURVE",
            ]
            .into_iter()
            .any(|k| instance.part(k).is_some())
            {
                return Some(false);
            }
            let args = instance
                .parts
                .first()
                .map_or(&[][..], |(_, a)| a.as_slice());
            (instance.keyword().to_owned(), ref_at(args, 1))
        };
        match keyword.as_str() {
            "LINE" | "PARABOLA" => Some(true),
            "CIRCLE" | "ELLIPSE" | "HYPERBOLA" => Some(false),
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" | "TRIMMED_CURVE"
            | "OFFSET_CURVE_3D" => self.length_parameter(inner),
            _ => None,
        }
    }

    /// How a curve's parameter here follows the file's: this kernel's per
    /// file unit, as ISO 10303-42 parameterizes the curve. A line's runs
    /// by its vector, a length; a circle's and an ellipse's are angles; a
    /// hyperbola's and a spline's are their own; a parabola's here is the
    /// distance along its axis's normal, twice its focal length per file
    /// unit. `None` for a curve not read that way.
    fn curve_units(&mut self, id: u64) -> Option<f64> {
        let (keyword, args) = {
            let instance = self.instance(id).ok()?;
            if [
                "B_SPLINE_CURVE",
                "B_SPLINE_CURVE_WITH_KNOTS",
                "BEZIER_CURVE",
                "UNIFORM_CURVE",
                "QUASI_UNIFORM_CURVE",
            ]
            .into_iter()
            .any(|k| instance.part(k).is_some())
            {
                return Some(1.0);
            }
            (
                instance.keyword(),
                instance
                    .parts
                    .first()
                    .map_or(&[][..], |(_, a)| a.as_slice()),
            )
        };
        let length = self.report.scale_mm;
        match keyword {
            "LINE" => {
                let vector = self.args(ref_at(args, 2), "VECTOR").ok()?;
                Some(vector.get(2).and_then(Arg::number).unwrap_or(1.0) * length)
            }
            "CIRCLE" | "ELLIPSE" => Some(self.angle_scale),
            "HYPERBOLA" => Some(1.0),
            "PARABOLA" => Some(2.0 * args.get(2).and_then(Arg::number)? * length),
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" | "TRIMMED_CURVE"
            | "OFFSET_CURVE_3D" => self.curve_units(ref_at(args, 1)),
            _ => None,
        }
    }

    /// A curve in a surface's parameters, its coordinates as written: a
    /// line, a circle, an ellipse, a spline, or a trim of one, which is its
    /// basis.
    fn planar_curve(&mut self, id: u64) -> OgeomResult<Option<PlanarCurve>> {
        let instance = self.instance(id)?;
        let (base, knotted, weights) = (
            instance.part("B_SPLINE_CURVE"),
            instance.part("B_SPLINE_CURVE_WITH_KNOTS"),
            instance.part("RATIONAL_B_SPLINE_CURVE"),
        );
        // A spline's degree, control points and knots sit at these places
        // in a simple instance, and split between parts in a complex one.
        let spline = match (base, knotted) {
            (Some(base), Some(knotted)) => {
                Some((base.first(), base.get(1), knotted.first(), knotted.get(1)))
            }
            (None, Some(knotted)) => Some((
                knotted.get(1),
                knotted.get(2),
                knotted.get(6),
                knotted.get(7),
            )),
            _ => None,
        };
        if let Some((degree, control, mults, knots)) = spline {
            let degree = file_degree(degree, 0.0)?;
            let mut points = Vec::new();
            for r in control
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::reference)
            {
                points.push(self.point2(r)?);
            }
            let expanded = expanded_knots(
                mults.and_then(Arg::list).unwrap_or(&[]),
                knots.and_then(Arg::list).unwrap_or(&[]),
            )?;
            let flat: Vec<f64> = weights
                .and_then(|w| w.first())
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::number)
                .collect();
            let control = points
                .into_iter()
                .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                .map(|(p, w)| Weighted::new(p, *w, self.tol))
                .collect::<OgeomResult<Vec<_>>>()?;
            let knots = KnotVector::new(expanded, degree)?;
            return Ok(Some(
                ogeom_geom::BSpline2d::rational(knots, control)?.into(),
            ));
        }
        let args = instance
            .parts
            .first()
            .map_or(&[][..], |(_, a)| a.as_slice());
        let out: PlanarCurve = match instance.keyword() {
            "LINE" => {
                let through = self.point2(ref_at(args, 1))?;
                let vector = self.args(ref_at(args, 2), "VECTOR")?;
                let direction = self.direction2(vector[1].reference().unwrap_or(0))?;
                // The vector's magnitude is the parameter's pace along the
                // line. A unit one is the length a line here runs by; any
                // other is the straight spline of degree one at that pace.
                let magnitude = vector.get(2).and_then(Arg::number).unwrap_or(1.0);
                if magnitude == 1.0 {
                    let axis = ogeom_math::Axis2 {
                        location: through,
                        direction,
                    };
                    ogeom_geom::Line2d::over(axis, -SURFACE_EXTENT, SURFACE_EXTENT)?.into()
                } else {
                    let step = direction.vector() * magnitude;
                    let ends = vec![
                        through - step * SURFACE_EXTENT,
                        through + step * SURFACE_EXTENT,
                    ];
                    let knots = vec![
                        -SURFACE_EXTENT,
                        -SURFACE_EXTENT,
                        SURFACE_EXTENT,
                        SURFACE_EXTENT,
                    ];
                    ogeom_geom::BSpline2d::new(KnotVector::new(knots, 1)?, ends, self.tol)?.into()
                }
            }
            "CIRCLE" => {
                let frame = self.frame2(ref_at(args, 1))?;
                let radius = args.get(2).and_then(Arg::number).unwrap_or(0.0);
                ogeom_geom::Circle2d::new(ogeom_math::Circle2::new(frame, radius, self.tol)?).into()
            }
            "ELLIPSE" => {
                let frame = self.frame2(ref_at(args, 1))?;
                let a = args.get(2).and_then(Arg::number).unwrap_or(0.0);
                let b = args.get(3).and_then(Arg::number).unwrap_or(0.0);
                ogeom_geom::Ellipse2d::new(ogeom_math::Ellipse2::new(frame, a, b, self.tol)?).into()
            }
            "TRIMMED_CURVE" => {
                // The trims restate the edge's ends; the basis carries the
                // parameter, so it is the pcurve.
                return self.planar_curve(ref_at(args, 1));
            }
            _ => return Ok(None),
        };
        Ok(Some(out))
    }

    /// A point of a surface's parameters, as written.
    fn point2(&mut self, id: u64) -> OgeomResult<ogeom_math::Point2> {
        let args = self.args(id, "CARTESIAN_POINT")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            ogeom_bail!(Construction, "#{id}: a point without coordinates");
        };
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0);
        Ok(ogeom_math::Point2::new(value(0), value(1)))
    }

    fn direction2(&mut self, id: u64) -> OgeomResult<ogeom_math::Direction2> {
        let args = self.args(id, "DIRECTION")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            ogeom_bail!(Construction, "#{id}: a direction without components");
        };
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0);
        ogeom_math::Direction2::new(ogeom_math::Vector2::new(value(0), value(1)), self.tol)
    }

    /// An `AXIS2_PLACEMENT_2D` as a frame, its reference direction the
    /// standard's default where the file leaves it out.
    fn frame2(&mut self, id: u64) -> OgeomResult<ogeom_math::Frame2> {
        let args = self.args(id, "AXIS2_PLACEMENT_2D")?;
        let origin = self.point2(ref_at(args, 1))?;
        let x = match args.get(2).and_then(Arg::reference) {
            Some(r) => self.direction2(r)?,
            None => ogeom_math::Direction2::X,
        };
        Ok(ogeom_math::Frame2::new(origin, x))
    }
}
