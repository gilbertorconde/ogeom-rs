//! Screw threads cut into a cylindrical face.
//!
//! A thread is a groove, the space between two teeth, swept along a helix
//! and cut from the solid. The groove's section lies in a plane through the
//! axis: a trapezoid whose flanks lean by the flank angle from the plane
//! square to the axis, with a flat at its floor (the root) and a flat land
//! left on the face between neighbouring turns (the crest). Along the
//! axis the cut is held to the thread's length by a turned clip shape,
//! which also carries the run-out, the stretch over which the groove's
//! depth falls to nothing.
//!
//! Which side the material is on is read from the face: a face whose
//! material lies inside it (a shaft) takes an external thread, cut in
//! towards the axis; one whose material lies outside it (a bore) takes an
//! internal thread, cut out away from it.

use core::f64::consts::TAU;

use ogeom_algo::{
    Built, Containment, History, classify_in_solid_exact, extent_along, face_normal, make_face,
    make_polygon, make_revolution,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, SameKey, Shape, ShapeType, explore_unique};

use crate::sweep::make_helical_sweep;

/// The form of a thread: its section, its hand and starts, and where it
/// runs on the face.
///
/// The section is read in a plane through the axis. `crest` is the flat
/// land left on the face between neighbouring grooves and `root` the flat
/// at the groove's floor, `depth` below the face; the flanks join them,
/// each leaning `flank_angle` from the plane square to the axis. The four
/// close one pitch: `crest + root + 2 depth tan(flank_angle) = pitch`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThreadForm {
    /// The distance along the axis between neighbouring turns.
    pub pitch: f64,
    /// How far the groove's floor stands below the face, radially.
    pub depth: f64,
    /// Each flank's lean from the plane square to the axis, in radians:
    /// half the included thread angle.
    pub flank_angle: f64,
    /// The width along the axis of the land left on the face between two
    /// grooves.
    pub crest: f64,
    /// The width along the axis of the groove's floor.
    pub root: f64,
    /// A left-hand thread: its grooves wind anticlockwise as they advance,
    /// seen looking along the advance; a right-hand thread's wind
    /// clockwise.
    pub left_handed: bool,
    /// How many grooves run side by side, each a pitch from the next; the
    /// lead, the advance per turn, is `starts * pitch`.
    pub starts: u32,
    /// How far along the axis the thread runs; `None` for the face's
    /// whole length.
    pub length: Option<f64>,
    /// Start at the face's end the axis points to rather than the one it
    /// points from.
    pub from_far_end: bool,
    /// The length at the thread's end over which the groove's depth falls
    /// linearly to nothing; zero ends the groove square.
    pub run_out: f64,
    /// Turn a 45 degree lead-in, as deep as the thread, off the end the
    /// thread starts from.
    pub chamfer: bool,
}

impl ThreadForm {
    /// The ISO metric basic form for an external thread of `pitch`: a 60
    /// degree thread angle, a depth of `5/8 H` (with `H = pitch * sqrt(3)/2`
    /// the height of the fundamental triangle), a crest of `pitch / 8` and
    /// a root of `pitch / 4`. Right-hand, one start, the face's whole
    /// length, no run-out and no chamfer.
    #[must_use]
    pub fn metric(pitch: f64) -> Self {
        let h = pitch * 3.0_f64.sqrt() / 2.0;
        Self {
            pitch,
            depth: h * 5.0 / 8.0,
            flank_angle: core::f64::consts::FRAC_PI_6,
            crest: pitch / 8.0,
            root: pitch / 4.0,
            left_handed: false,
            starts: 1,
            length: None,
            from_far_end: false,
            run_out: 0.0,
            chamfer: false,
        }
    }

    /// The ISO metric basic form for an internal thread of `pitch`, cut
    /// into a bore of the basic minor diameter: as [`ThreadForm::metric`],
    /// with the crest of `pitch / 4` on the bore and the root of
    /// `pitch / 8` at the major diameter.
    #[must_use]
    pub fn metric_internal(pitch: f64) -> Self {
        Self {
            crest: pitch / 4.0,
            root: pitch / 8.0,
            ..Self::metric(pitch)
        }
    }

    /// This form cut `depth` deep, its pitch, flank angle and root kept and
    /// its crest what is left of the pitch: the form for a bore drilled to
    /// a given size, or a thread truncated at its crest.
    #[must_use]
    pub fn with_depth(self, depth: f64) -> Self {
        Self {
            depth,
            crest: self.pitch - self.root - 2.0 * depth * self.flank_angle.tan(),
            ..self
        }
    }

    /// Half the groove's width along the axis at `t` above the face
    /// (negative below it, towards the floor).
    fn half_width(&self, t: f64) -> f64 {
        self.root / 2.0 + (t + self.depth) * self.flank_angle.tan()
    }

    fn validate(&self, tol: Tolerances) -> OgeomResult<()> {
        let eps = tol.confusion();
        let finite = [
            self.pitch,
            self.depth,
            self.flank_angle,
            self.crest,
            self.root,
            self.run_out,
        ]
        .iter()
        .all(|v| v.is_finite());
        if !finite {
            ogeom_bail!(
                Construction,
                "a thread form holds a value that is not finite"
            );
        }
        if self.pitch <= eps || self.depth <= eps {
            ogeom_bail!(
                Construction,
                "a thread needs a positive pitch and depth; got {} and {}",
                self.pitch,
                self.depth
            );
        }
        if !(0.0..core::f64::consts::FRAC_PI_2).contains(&self.flank_angle) {
            ogeom_bail!(
                Construction,
                "a flank leans from 0 up to a right angle; {} does not",
                self.flank_angle
            );
        }
        if self.crest <= eps || self.root < 0.0 {
            ogeom_bail!(
                Construction,
                "a thread leaves a land on the face between its grooves and a \
                 floor of no negative width; got a crest of {} and a root of {}",
                self.crest,
                self.root
            );
        }
        let closed = self.crest + self.root + 2.0 * self.depth * self.flank_angle.tan();
        if (closed - self.pitch).abs() > eps {
            ogeom_bail!(
                Construction,
                "a thread's crest, root and flanks span {closed} along the axis, \
                 not its pitch {}",
                self.pitch
            );
        }
        if self.starts == 0 {
            ogeom_bail!(Construction, "a thread has at least one start");
        }
        if let Some(length) = self.length
            && !(length.is_finite() && length > eps)
        {
            ogeom_bail!(Construction, "a thread of length {length} cuts nothing");
        }
        if self.run_out < 0.0 {
            ogeom_bail!(
                Construction,
                "a run-out of {} is not a length",
                self.run_out
            );
        }
        Ok(())
    }
}

/// The axis and radius of a cylindrical face, placed.
fn cylinder_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<(Axis, f64)> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        ogeom_bail!(Construction, "a thread is cut on a face");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let mut at = surface;
    let cylinder = loop {
        match at {
            SurfaceGeometry::Cylinder(c) => break c.cylinder(),
            SurfaceGeometry::Trimmed(t) => at = t.basis(),
            _ => ogeom_bail!(Construction, "a thread is cut on a cylindrical face"),
        }
    };
    let placement = face.transform(model.datums())?;
    let frame = cylinder.frame();
    let axis = Axis {
        location: placement.apply(frame.origin()),
        direction: placement.apply_direction(frame.z(), tol)?,
    };
    Ok((axis, cylinder.radius() * placement.scale_factor().abs()))
}

/// Points in the plane through the axis that holds the thread's section,
/// named by `(t, s)`: `t` above the face (towards the air) and `s` along
/// the axis from where the thread starts.
#[derive(Clone, Copy)]
struct Section {
    start: Point,
    along: Vector,
    out: Vector,
    radius: f64,
    /// `1` where the air is outside the face (a shaft), `-1` inside it.
    air: f64,
}

impl Section {
    fn point(&self, t: f64, s: f64) -> Point {
        self.start + self.along * s + self.out * (self.radius + self.air * t)
    }

    /// The planar face bounded by the polygon through `corners`.
    fn face(
        &self,
        model: &mut Model,
        corners: &[(f64, f64)],
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let points: Vec<Point> = corners.iter().map(|&(t, s)| self.point(t, s)).collect();
        let wire = make_polygon(model, &points, true, tol)?.shape;
        let normal = Direction::new(self.out.cross(self.along), tol)?;
        let x = Direction::new(self.along, tol)?;
        let plane = Plane::new(Frame::new(self.start, normal, x, tol)?);
        Ok(make_face(model, PlaneSurface::new(plane).into(), &[wire], tol)?.shape)
    }
}

/// Cut a thread of `form` on `face`, a cylindrical face of `solid`.
///
/// The side is read from the face's orientation: material inside the face
/// (a shaft) takes an external thread, cut in towards the axis, and
/// material outside it (a bore) an internal one, cut out away from it.
/// The thread starts at one end of the face, where the face meets its
/// least extent along its own axis (its greatest with `from_far_end`),
/// and runs `form.length` from there, the run-out included. Along the
/// thread the groove is the form's section exactly: its floor `depth`
/// from the face, its flanks at the flank angle. Where the thread ends
/// inside the face the groove ends square across the axis, or over the
/// run-out on a cone that takes its depth linearly to nothing at the
/// thread's end.
///
/// Where the solid stops at the face's end, the groove runs out of the
/// material there. Where it carries on past the end (a shoulder), the
/// groove ends in the end's plane, the run-out's cone too.
///
/// The groove's walls are fitted through their exact screw images within
/// a hundred times `tol.confusion()`; its ends, on the clip, are exact.
///
/// The result is the solid less the grooves (and the chamfer), with the
/// history of the cuts, in which the face generates the threaded solid.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// form does not close one pitch ([`ThreadForm`]) or holds a value out of
/// range; if `face` is not a cylindrical face of `solid`; if the thread is
/// longer than the face, the run-out not shorter than the thread, or an
/// external groove reaches the axis; if a chamfer is asked at an end the
/// solid carries on past; and as the sweep and the booleans report.
#[allow(clippy::too_many_lines)]
pub fn make_thread(
    model: &mut Model,
    solid: &Shape,
    face: &Shape,
    form: &ThreadForm,
    tol: Tolerances,
) -> OgeomResult<Built> {
    form.validate(tol)?;
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "a thread is cut on a face");
    }
    let key = SameKey(face.clone());
    if !explore_unique(model, solid, ShapeType::Face)?
        .into_iter()
        .any(|f| SameKey(f) == key)
    {
        ogeom_bail!(
            Construction,
            "the face to thread is not a face of the solid"
        );
    }
    let (axis, radius) = cylinder_of(model, face, tol)?;
    let (low, high) = extent_along(model, face, axis, tol)?;
    let span = high - low;
    let length = form.length.unwrap_or(span);
    if length > span + tol.confusion() {
        ogeom_bail!(
            Construction,
            "a thread of length {length} is longer than its face, {span} along the axis"
        );
    }
    if form.run_out > length - tol.confusion() && form.run_out > 0.0 {
        ogeom_bail!(
            Construction,
            "a run-out of {} does not fit a thread of length {length}",
            form.run_out
        );
    }

    // Which side the material is on: the face's normal points to the air.
    let (on, normal) = face_normal(model, face, tol)?;
    let d = axis.direction.vector();
    let radial = {
        let r = on - axis.location;
        r - d * r.dot(d)
    };
    let air = if normal.dot(radial) > 0.0 { 1.0 } else { -1.0 };
    let depth = form.depth;
    if air > 0.0 && depth >= radial.magnitude() - tol.confusion() {
        ogeom_bail!(
            Construction,
            "a groove {depth} deep reaches the axis of a shaft {radius} in radius"
        );
    }

    let (along, s0) = if form.from_far_end {
        (-d, high)
    } else {
        (d, low)
    };
    let section = Section {
        start: axis.location + d * s0,
        along,
        out: radial / radial.magnitude(),
        radius,
        air,
    };
    let pitch = form.pitch;
    let lead = pitch * f64::from(form.starts);
    let tan = form.flank_angle.tan();
    // The groove's section reaches `e` past the face into the air, so that
    // none of the tool's edges lies on the face; it stays narrower than a
    // pitch there, so one turn clears the next.
    let mut e = depth;
    if tan > 0.0 {
        e = e.min(form.crest / (4.0 * tan));
    }
    if air < 0.0 {
        e = e.min(radius / 4.0);
    }
    // A probe a short way past each end, halfway down the groove: whether
    // the solid stops at that end.
    let margin = (pitch / 4.0).min(depth);
    let free_at = |model: &Model, near: f64, far: f64| -> OgeomResult<bool> {
        for s in [near, far] {
            let probe = section.point(-depth / 2.0, s);
            if classify_in_solid_exact(model, solid, probe, tol)? != Containment::Out {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let start_free = free_at(model, -margin / 2.0, -margin)?;
    let end_free =
        length >= span - tol.confusion() && free_at(model, length + margin / 2.0, length + margin)?;

    let mut steps: Vec<History> = Vec::new();
    let mut current = solid.clone();

    if form.chamfer {
        if !start_free {
            ogeom_bail!(
                Construction,
                "a thread's chamfer is turned off an end where the solid stops; \
                 the solid carries on past this one"
            );
        }
        // A 45 degree cone from the groove's floor at the end out to the
        // face a depth along, and on past both a little into the air.
        let g = margin;
        let corners = [(-depth - g, -g), (g, -g), (g, depth + g)];
        let profile = section.face(model, &corners, tol)?;
        let tool = make_revolution(model, &profile, axis, TAU, tol)?.shape;
        let cut = ogeom_bool::cut(model, &current, &tool, tol)?;
        current = cut.shape;
        steps.push(cut.history);
    }

    // The clip: the stretch of the axis the thread runs over, between
    // radii clear of the groove's floor and its top, with the run-out's
    // cone at the far end.
    let g = (0.1 * depth).min(e);
    let s_start = if start_free { -margin } else { 0.0 };
    let (t_lo, t_hi) = (-depth - g, e + g);
    let mut clip_corners = vec![(t_lo, s_start), (t_hi, s_start)];
    if form.run_out > 0.0 {
        // The cone through the floor where the run-out starts and the face
        // where the thread ends, carried a little past the face; where the
        // solid carries on past the face's end, it stops in the end's plane
        // rather than cut past it.
        let slope = form.run_out / depth;
        let blocked = length >= span - tol.confusion() && !end_free;
        let t_top = if blocked { 0.0 } else { e / 2.0 };
        let s_top = length + t_top * slope;
        clip_corners.push((t_hi, s_top));
        clip_corners.push((t_top, s_top));
        let s_low = length + t_lo * slope;
        if s_low > s_start + tol.confusion() {
            clip_corners.push((t_lo, s_low));
        } else {
            clip_corners[0] = ((s_start - length) / slope, s_start);
        }
    } else {
        let s_top = if end_free { length + margin } else { length };
        clip_corners.push((t_hi, s_top));
        clip_corners.push((t_lo, s_top));
    }
    let clip_profile = section.face(model, &clip_corners, tol)?;
    let clip = make_revolution(model, &clip_profile, axis, TAU, tol)?.shape;

    // The groove's section, a trapezoid (a triangle with no root) from the
    // floor up past the face, swept a lead beyond the clip each way.
    let (w0, w1) = (form.root / 2.0, form.half_width(e));
    let s_lo = s_start - lead - pitch;
    let s_hi = length + margin + lead + pitch;
    let turns = (s_hi - s_lo) / lead;
    for k in 0..form.starts {
        let s = s_lo + pitch * f64::from(k);
        let corners: Vec<(f64, f64)> = if w0 > tol.confusion() {
            vec![(-depth, s - w0), (e, s - w1), (e, s + w1), (-depth, s + w0)]
        } else {
            vec![(-depth, s), (e, s - w1), (e, s + w1)]
        };
        let profile = section.face(model, &corners, tol)?;
        let screw = Axis {
            location: section.start,
            direction: Direction::new(along, tol)?,
        };
        let groove = make_helical_sweep(
            model,
            &profile,
            screw,
            lead,
            turns,
            form.left_handed,
            0.0,
            tol,
        )?
        .shape;
        let held = ogeom_bool::common(model, &groove, &clip, tol)?.shape;
        let mut cut = ogeom_bool::cut(model, &current, &held, tol)?;
        cut.history.generate(face, cut.shape.clone());
        current = cut.shape;
        steps.push(cut.history);
    }
    Ok(Built::new(current, History::chain(&steps)))
}
