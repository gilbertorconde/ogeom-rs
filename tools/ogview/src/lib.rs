//! Rendering a tessellation to an image, so that a wrong result is *visible*.
//!
//! Numeric checks catch what you thought to check. A picture catches the rest —
//! a face wound inside out, a hole where two surfaces failed to meet, a seam
//! that did not weld. Those have all been real defects in this kernel, and each
//! would have been obvious in one frame.
//!
//! # Why there is no window here
//!
//! A real-time viewer needs a GPU, a windowing system and a driver, none of
//! which a headless build has. That would make the viewer the one corner of the
//! repository `tools/check.sh` cannot verify — and a visual check nobody can
//! run in CI is a visual check nobody runs.
//!
//! So this renders in software, to an image, deterministically. It has no
//! dependencies beyond the kernel itself, works anywhere `cargo test` does, and
//! can be *tested*: a sphere's silhouette is a circle, a cube seen down an axis
//! is a square, and both of those are assertions rather than opinions.
//!
//! A windowed viewer can wrap this later. It would add interaction, not
//! correctness.
//!
//! # The format is PPM
//!
//! Uncompressed binary PPM: a short text header and three bytes a pixel. Every
//! image viewer reads it, `convert` and `ffmpeg` transcode it, and it costs no
//! dependency. A PNG would be a fifth the size and an encoder to maintain.

use og_core::{OgResult, Tolerances, og_bail};
use og_math::{Direction, Point, Vector};
use og_topo::Triangulation;

/// Where the camera is and what it can see.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Where the camera sits.
    pub eye: Point,
    /// What it looks at.
    pub target: Point,
    /// Which way is up, as a reference — its component along the view is
    /// removed, so it need not be perpendicular.
    pub up: Vector,
    /// Vertical field of view, in radians.
    pub field_of_view: f64,
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            eye: Point::new(1.0, 1.0, 1.0),
            target: Point::ORIGIN,
            up: Vector::Z,
            field_of_view: core::f64::consts::FRAC_PI_4,
            width: 800,
            height: 600,
        }
    }
}

impl Camera {
    /// A camera framing a shape's bounds from a given direction.
    ///
    /// Placed far enough back that the whole bound is inside the field of view,
    /// with a margin — so the answer to "did it render" is never "it was
    /// off-screen".
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) if the bound
    /// is empty, or the direction has no length.
    pub fn framing(bounds: &og_math::Aabb, from: Vector, tol: Tolerances) -> OgResult<Self> {
        let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
            og_bail!(Construction, "there is nothing to frame");
        };
        let direction = Direction::new(from, tol)?;
        let centre = Point::from_vector((low.to_vector() + high.to_vector()) * 0.5);
        let radius = (high - low).magnitude() * 0.5;
        let field_of_view = core::f64::consts::FRAC_PI_4;
        // Far enough that a sphere of that radius fits, plus a fifth for air.
        let distance = radius / (field_of_view * 0.5).sin() * 1.2;

        // Any up reference not parallel to the view will do; the frame removes
        // its along-view component anyway.
        let up = if direction.vector().cross(Vector::Z).magnitude() > tol.confusion() {
            Vector::Z
        } else {
            Vector::Y
        };
        Ok(Self {
            eye: centre + direction.vector() * distance,
            target: centre,
            up,
            field_of_view,
            ..Self::default()
        })
    }

    /// This camera at a different image size.
    #[must_use]
    pub const fn sized(mut self, width: usize, height: usize) -> Self {
        self.width = width;
        self.height = height;
        self
    }
}

/// An image, three bytes a pixel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Image {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Red, green and blue for each pixel, row-major from the top left.
    pub pixels: Vec<[u8; 3]>,
}

impl Image {
    /// The colour at a pixel, or `None` outside the image.
    #[must_use]
    pub fn at(&self, x: usize, y: usize) -> Option<[u8; 3]> {
        if x >= self.width || y >= self.height {
            return None;
        }
        self.pixels.get(y * self.width + x).copied()
    }

    /// How many pixels are not the background.
    ///
    /// The cheapest useful measure of "did anything render", and the basis of
    /// the silhouette checks.
    #[must_use]
    pub fn covered(&self, background: [u8; 3]) -> usize {
        self.pixels.iter().filter(|p| **p != background).count()
    }

    /// The bounding box of everything that is not the background, as
    /// `(min x, min y, max x, max y)`.
    #[must_use]
    pub fn footprint(&self, background: [u8; 3]) -> Option<(usize, usize, usize, usize)> {
        let mut found: Option<(usize, usize, usize, usize)> = None;
        for y in 0..self.height {
            for x in 0..self.width {
                if self.at(x, y) != Some(background) {
                    found = Some(match found {
                        None => (x, y, x, y),
                        Some((lx, ly, hx, hy)) => (lx.min(x), ly.min(y), hx.max(x), hy.max(y)),
                    });
                }
            }
        }
        found
    }

    /// This image as a binary PPM.
    #[must_use]
    pub fn to_ppm(&self) -> Vec<u8> {
        let mut out = format!("P6\n{} {}\n255\n", self.width, self.height).into_bytes();
        for pixel in &self.pixels {
            out.extend_from_slice(pixel);
        }
        out
    }
}

/// How a rendering is shaded.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Style {
    /// The colour of nothing.
    pub background: [u8; 3],
    /// The colour a surface is at full illumination.
    pub surface: [u8; 3],
    /// The colour a surface facing *away* is drawn in.
    ///
    /// Not the same as the lit colour, and that is the point. A face wound
    /// inside out still renders — it is not culled — and it renders in a colour
    /// that says so, because a hole where a face should be looks the same as a
    /// face that was never built, and an obviously wrong colour does not.
    pub reversed: [u8; 3],
    /// How much light a surface gets when it faces fully away from the light.
    pub ambient: f64,
    /// Which way the light comes from.
    pub light: Vector,
}

impl Default for Style {
    fn default() -> Self {
        Self {
            background: [16, 18, 22],
            surface: [200, 205, 215],
            reversed: [190, 70, 60],
            ambient: 0.25,
            light: Vector::new(-0.4, -0.6, 1.0),
        }
    }
}

/// Render a mesh.
///
/// Painter-free: every triangle is depth-tested per pixel, so the result does
/// not depend on the order the triangles arrive in. That matters for a test —
/// a renderer whose output depends on triangle order cannot be compared against
/// anything.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the image size
/// is zero, the field of view is not a positive angle below half a turn, or the
/// camera sits on its own target.
pub fn render(
    mesh: &Triangulation,
    camera: &Camera,
    style: &Style,
    tol: Tolerances,
) -> OgResult<Image> {
    if camera.width == 0 || camera.height == 0 {
        og_bail!(Construction, "an image needs a positive size");
    }
    if !camera.field_of_view.is_finite()
        || camera.field_of_view <= 0.0
        || camera.field_of_view >= core::f64::consts::PI
    {
        og_bail!(
            Construction,
            "a field of view of {} radians sees nothing",
            camera.field_of_view
        );
    }
    let view = Direction::new(camera.target - camera.eye, tol)?;
    let right = Direction::from_cross(view.vector(), camera.up, tol)?;
    let above = Direction::from_cross(right.vector(), view.vector(), tol)?;
    let light = Direction::new(style.light, tol)?;

    let mut pixels = vec![style.background; camera.width * camera.height];
    let mut depth = vec![f64::MAX; camera.width * camera.height];

    #[allow(clippy::cast_precision_loss)]
    let (w, h) = (camera.width as f64, camera.height as f64);
    let half = (camera.field_of_view * 0.5).tan();
    let aspect = w / h;

    // Camera space: x right, y up, z forward. A point is visible when z > 0.
    let to_camera = |p: Point| {
        let d = p - camera.eye;
        Vector::new(
            d.dot(right.vector()),
            d.dot(above.vector()),
            d.dot(view.vector()),
        )
    };
    let to_screen = |c: Vector| {
        let x = (c.x / (c.z * half * aspect)).mul_add(0.5, 0.5) * w;
        // Screen y runs down, camera y runs up.
        let y = (1.0 - (c.y / (c.z * half)).mul_add(0.5, 0.5)) * h;
        (x, y)
    };

    for triangle in &mesh.triangles {
        let Some(corners) = fetch(mesh, *triangle) else {
            continue;
        };
        let camera_space = corners.map(to_camera);
        // Nothing behind the eye. A triangle straddling the plane needs
        // clipping, which this does not do — it is dropped, and the docs say
        // so, because a partly-drawn triangle is worse than a missing one for
        // a picture whose job is to be trusted.
        if camera_space.iter().any(|c| c.z <= tol.confusion()) {
            continue;
        }
        let screen = camera_space.map(to_screen);

        // The facing test is done in *camera space* on the geometry, not from
        // the stored normals: a mesh assembled by hand may have none, and the
        // winding is what a renderer can always see.
        let normal = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
        let facing = normal.dot(camera.eye - corners[0]) > 0.0;
        let shade = if normal.magnitude() <= tol.confusion() {
            style.ambient
        } else {
            let unit = normal * (1.0 / normal.magnitude());
            let towards = if facing { unit } else { -unit };
            style
                .ambient
                .max(towards.dot(-light.vector()).max(0.0))
                .min(1.0)
        };
        let base = if facing {
            style.surface
        } else {
            style.reversed
        };
        let colour = [
            channel(base[0], shade),
            channel(base[1], shade),
            channel(base[2], shade),
        ];

        fill(
            &mut pixels,
            &mut depth,
            camera.width,
            camera.height,
            screen,
            [camera_space[0].z, camera_space[1].z, camera_space[2].z],
            colour,
        );
    }

    Ok(Image {
        width: camera.width,
        height: camera.height,
        pixels,
    })
}

/// The three corner positions of a triangle, or `None` if it names a vertex the
/// mesh does not have.
fn fetch(mesh: &Triangulation, triangle: [u32; 3]) -> Option<[Point; 3]> {
    Some([
        *mesh.positions.get(triangle[0] as usize)?,
        *mesh.positions.get(triangle[1] as usize)?,
        *mesh.positions.get(triangle[2] as usize)?,
    ])
}

/// One channel at a given illumination.
fn channel(value: u8, shade: f64) -> u8 {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "clamped to 0..=255 immediately before the cast"
    )]
    {
        (f64::from(value) * shade).clamp(0.0, 255.0) as u8
    }
}

/// Fill one triangle, depth-testing each pixel.
fn fill(
    pixels: &mut [[u8; 3]],
    depth: &mut [f64],
    width: usize,
    height: usize,
    screen: [(f64, f64); 3],
    z: [f64; 3],
    colour: [u8; 3],
) {
    let area = edge(screen[0], screen[1], screen[2]);
    if area.abs() < f64::EPSILON {
        return;
    }
    let low_x = screen.iter().map(|p| p.0).fold(f64::MAX, f64::min);
    let high_x = screen.iter().map(|p| p.0).fold(f64::MIN, f64::max);
    let low_y = screen.iter().map(|p| p.1).fold(f64::MAX, f64::min);
    let high_y = screen.iter().map(|p| p.1).fold(f64::MIN, f64::max);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        reason = "clamped into the image before the cast"
    )]
    let (x0, x1, y0, y1) = (
        low_x.floor().max(0.0) as usize,
        (high_x.ceil() as isize).clamp(0, width as isize) as usize,
        low_y.floor().max(0.0) as usize,
        (high_y.ceil() as isize).clamp(0, height as isize) as usize,
    );

    for y in y0..y1 {
        for x in x0.min(width)..x1 {
            #[allow(clippy::cast_precision_loss)]
            let at = (x as f64 + 0.5, y as f64 + 0.5);
            let (a, b, c) = (
                edge(screen[1], screen[2], at) / area,
                edge(screen[2], screen[0], at) / area,
                edge(screen[0], screen[1], at) / area,
            );
            if a < 0.0 || b < 0.0 || c < 0.0 {
                continue;
            }
            // Depth interpolated in screen space is not exact for a perspective
            // projection, but it is monotonic in true depth, which is all a
            // visibility test needs.
            let here = a * z[0] + b * z[1] + c * z[2];
            let index = y * width + x;
            if here < depth[index] {
                depth[index] = here;
                pixels[index] = colour;
            }
        }
    }
}

/// Twice the signed area of a triangle in screen space.
fn edge(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0).mul_add(c.1 - a.1, -((b.1 - a.1) * (c.0 - a.0)))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use og_algo::{make_box, make_sphere, shape_bounds};
    use og_math::Frame;
    use og_mesh::{Deflection, triangulate};
    use og_topo::Model;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 0.01,
            ..Deflection::default()
        }
    }

    fn shot(mesh: &Triangulation, bounds: &og_math::Aabb, from: Vector) -> Image {
        let camera = Camera::framing(bounds, from, T).unwrap().sized(200, 200);
        render(mesh, &camera, &Style::default(), T).unwrap()
    }

    #[test]
    fn a_sphere_renders_as_a_circle() {
        // The cheapest real check there is: a sphere looks the same from every
        // direction, and its silhouette is round. A footprint that is not
        // square-ish means the projection is wrong; one that is empty means
        // nothing rendered at all.
        let mut model = Model::new();
        let built = make_sphere(&mut model, Frame::WORLD, 2.0, T).unwrap();
        let mesh = triangulate(&model, &built.shape, fine(), T).unwrap();
        let bounds = shape_bounds(&model, &built.shape, T).unwrap();

        for from in [Vector::X, Vector::Y, Vector::Z, Vector::new(1.0, 1.0, 1.0)] {
            let image = shot(&mesh, &bounds, from);
            let (lx, ly, hx, hy) = image.footprint(Style::default().background).unwrap();
            let (w, h) = (hx - lx, hy - ly);
            // Framing works from the bounding *box*, whose diagonal
            // circumscribes the sphere, so a sphere fills a little under half
            // the frame rather than most of it. Conservative on purpose: the
            // answer to "did it render" should never be "it was off-screen".
            assert!(w > 50, "the sphere barely rendered: {w} wide");
            let ratio = f64::from(u32::try_from(w).unwrap()) / f64::from(u32::try_from(h).unwrap());
            assert!(
                (ratio - 1.0).abs() < 0.06,
                "a sphere's silhouette should be round, got {ratio} from {from:?}"
            );

            // And it is filled, not an outline: a disc covers about pi/4 of its
            // bounding square.
            let covered = image.covered(Style::default().background);
            #[allow(clippy::cast_precision_loss)]
            let expected = core::f64::consts::PI / 4.0 * (w * h) as f64;
            #[allow(clippy::cast_precision_loss)]
            let found = covered as f64;
            assert!(
                (found / expected - 1.0).abs() < 0.05,
                "a disc should fill about pi/4 of its box, got {found} against {expected}"
            );
        }
    }

    #[test]
    fn a_cube_seen_down_an_axis_is_a_square_that_fills_its_footprint() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let mesh = triangulate(&model, &built.shape, fine(), T).unwrap();
        let bounds = shape_bounds(&model, &built.shape, T).unwrap();

        let image = shot(&mesh, &bounds, Vector::Z);
        let (lx, ly, hx, hy) = image.footprint(Style::default().background).unwrap();
        let (w, h) = (hx - lx, hy - ly);
        assert!(w.abs_diff(h) <= 2, "a cube face should be square: {w}x{h}");

        // Solid, not hollow: every pixel inside the footprint is covered.
        let covered = image.covered(Style::default().background);
        assert!(
            covered >= w * h - 2 * (w + h),
            "the square should be filled, {covered} of {}",
            w * h
        );
    }

    #[test]
    fn the_nearer_surface_wins_whatever_order_the_triangles_arrive_in() {
        // A painter's algorithm would make this depend on triangle order, and a
        // renderer whose output depends on order cannot be compared against
        // anything.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let mesh = triangulate(&model, &built.shape, fine(), T).unwrap();
        let bounds = shape_bounds(&model, &built.shape, T).unwrap();
        let camera = Camera::framing(&bounds, Vector::new(1.0, 1.0, 1.0), T)
            .unwrap()
            .sized(120, 120);

        let straight = render(&mesh, &camera, &Style::default(), T).unwrap();
        let mut shuffled = mesh.clone();
        shuffled.triangles.reverse();
        let other = render(&shuffled, &camera, &Style::default(), T).unwrap();
        assert_eq!(straight, other, "the image depended on triangle order");
    }

    #[test]
    fn a_face_wound_inside_out_renders_in_a_colour_that_says_so() {
        // Not culled. A hole where a face should be looks exactly like a face
        // that was never built; a face in an alarming colour does not.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let mut mesh = triangulate(&model, &built.shape, fine(), T).unwrap();
        let bounds = shape_bounds(&model, &built.shape, T).unwrap();
        for triangle in &mut mesh.triangles {
            triangle.swap(1, 2);
        }

        let image = shot(&mesh, &bounds, Vector::new(1.0, 1.0, 1.0));
        let style = Style::default();
        assert!(
            image.covered(style.background) > 1000,
            "an inside-out solid should still render"
        );
        let alarming = image
            .pixels
            .iter()
            .filter(|p| p[0] > p[1] && p[0] > p[2])
            .count();
        assert!(
            alarming > 1000,
            "an inside-out solid should be obviously wrong, got {alarming} \
             warning-coloured pixels"
        );
    }

    #[test]
    fn a_ppm_carries_its_own_size_and_every_pixel() {
        let image = Image {
            width: 2,
            height: 1,
            pixels: vec![[1, 2, 3], [4, 5, 6]],
        };
        let ppm = image.to_ppm();
        assert!(ppm.starts_with(b"P6\n2 1\n255\n"));
        assert_eq!(ppm.len(), "P6\n2 1\n255\n".len() + 6);
        assert_eq!(&ppm[ppm.len() - 6..], &[1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn a_camera_that_sees_nothing_is_refused() {
        let mesh = Triangulation::new();
        let base = Camera::default();
        for camera in [
            Camera { width: 0, ..base },
            Camera { height: 0, ..base },
            Camera {
                field_of_view: 0.0,
                ..base
            },
            Camera {
                field_of_view: core::f64::consts::PI,
                ..base
            },
            Camera {
                eye: Point::ORIGIN,
                target: Point::ORIGIN,
                ..base
            },
        ] {
            assert!(render(&mesh, &camera, &Style::default(), T).is_err());
        }
        assert!(Camera::framing(&og_math::Aabb::EMPTY, Vector::Z, T).is_err());
    }
}
