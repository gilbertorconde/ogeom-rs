//! Pairs of surfaces known to meet only along an edge one of them is
//! bounded by, which the general fuse does not intersect.
//!
//! A blend's band touches each face it was built against along its rail,
//! tangentially, and the wedge carrying it lies on that face with a leg of
//! the face's own surface, which splits the face along the rail. Asking the
//! band and the face for their meeting finds that contact again, and on a
//! fitted band, which touches its host only within its fit, the contact
//! walk has nothing smooth to follow. So the blend says so, and the pair is
//! read as apart.
//!
//! An entry is held on the calling thread. It names the band by its stored
//! id and by its points at a few parameters, read again from the gathered
//! face, so a surface stored later under the same id is another surface and
//! is intersected as any. The host is named by its points alone: a face
//! rebuilt on a copy of the host's surface stands on the same surface.

use core::cell::RefCell;
use std::collections::VecDeque;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_topo::{Model, SurfaceId};

/// How many pairs a thread keeps: the newest, the older dropped.
const KEPT: usize = 256;

/// A surface's points at a three by three grid of parameters, as bits.
type Mark = [u64; 27];

#[derive(Clone)]
struct Entry {
    band: (SurfaceId, Mark),
    host: Mark,
}

thread_local! {
    static PAIRS: RefCell<VecDeque<Entry>> = const { RefCell::new(VecDeque::new()) };
}

fn mark(surface: &SurfaceGeometry, tol: Tolerances) -> Option<Mark> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let at = |lo: f64, hi: f64, t: f64| {
        if (hi - lo).is_finite() {
            lo + (hi - lo) * t
        } else {
            t
        }
    };
    let mut out = [0_u64; 27];
    for (i, t) in [0.25, 0.5, 0.75].into_iter().enumerate() {
        for (j, s) in [0.25, 0.5, 0.75].into_iter().enumerate() {
            let p = surface.point_at(at(u0, u1, t), at(v0, v1, s), tol).ok()?;
            let k = (i * 3 + j) * 3;
            out[k] = p.x.to_bits();
            out[k + 1] = p.y.to_bits();
            out[k + 2] = p.z.to_bits();
        }
    }
    Some(out)
}

/// Note that the surface stored as `band` meets the one stored as `host`
/// only along an edge of a face on `band`, and that whatever carries that
/// face into a boolean with a face on `host` also carries a face on
/// `host`'s surface bounded by that edge: a blend's band and the face it
/// was built against, the wedge's leg lying on that face. A general fuse
/// on this thread then reads a face on `band` against a face on `host`, or
/// on a surface with `host`'s points, both at their stored placement, as
/// apart. A thread keeps the newest 256 such notes.
///
/// The caller vouches for the claim: the boolean does not check it.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) where either
/// id names no surface of `model`.
pub fn note_meeting_along_edge(
    model: &Model,
    band: SurfaceId,
    host: SurfaceId,
    tol: Tolerances,
) -> OgeomResult<()> {
    let read = |id: SurfaceId| -> OgeomResult<Mark> {
        let Some(surface) = model.geometry().surface(id) else {
            ogeom_bail!(
                Dangling,
                "a surface noted as meeting along an edge is not stored"
            );
        };
        mark(surface, tol).ok_or_else(|| {
            ogeom_core::ogeom_err!(Construction, "a noted surface cannot be evaluated")
        })
    };
    let entry = Entry {
        band: (band, read(band)?),
        host: read(host)?,
    };
    PAIRS.with(|pairs| {
        let mut pairs = pairs.borrow_mut();
        if pairs.len() == KEPT {
            pairs.pop_front();
        }
        pairs.push_back(entry);
    });
    Ok(())
}

/// The pairs noted on this thread, read once by a fuse before it hands its
/// pairs to the workers.
pub(crate) struct Noted(Vec<Entry>);

impl Noted {
    pub(crate) fn read() -> Self {
        Self(PAIRS.with(|pairs| pairs.borrow().iter().cloned().collect()))
    }

    /// Whether the surfaces `a` and `b`, stored as `ids` and gathered at
    /// their stored placement, were noted as meeting only along an edge.
    pub(crate) fn holds(
        &self,
        (a, b): (&SurfaceGeometry, &SurfaceGeometry),
        ids: (SurfaceId, SurfaceId),
        tol: Tolerances,
    ) -> bool {
        if !self
            .0
            .iter()
            .any(|e| e.band.0 == ids.0 || e.band.0 == ids.1)
        {
            return false;
        }
        let (Some(ma), Some(mb)) = (mark(a, tol), mark(b, tol)) else {
            return false;
        };
        self.0.iter().any(|e| {
            (e.band == (ids.0, ma) && e.host == mb) || (e.band == (ids.1, mb) && e.host == ma)
        })
    }
}
