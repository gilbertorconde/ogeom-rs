//! What one conversion keeps across its plans and builds: built faces'
//! areas by their content, the curves snapped to chains with what each read
//! of the faces beside it, and edge images on curved faces by every input
//! that decides them. An entry is reused only where everything it read is
//! unchanged, so a kept answer is the one a fresh solve gives.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::PlanarCurve;
use ogeom_math::Plane;
use ogeom_topo::{Model, Shape};

use super::planner::{EdgeSpec, Planner};
use super::snap::{Images, Snapped, image_on};
use super::{Carrier, Curved};
use crate::recognize::Canonical;

/// Faces' areas kept across the builds of one conversion, by the face's
/// content and the deflection (see [`ogeom_topo::face_content`]): a face's
/// area reads nothing else. A face whose area cannot be measured is asked
/// again.
#[derive(Default)]
pub(super) struct AreaCache {
    held: std::sync::Mutex<FastMap<Vec<u8>, f64>>,
}

impl AreaCache {
    /// The area of `face` as [`crate::surface_properties`] measures it at
    /// `deflection`.
    pub(super) fn area(
        &self,
        model: &Model,
        face: &Shape,
        deflection: ogeom_mesh::Deflection,
        tol: Tolerances,
    ) -> OgeomResult<f64> {
        let measure = || crate::surface_properties(model, face, deflection, tol).map(|m| m.mass);
        let Some(mut key) = ogeom_topo::face_content(model, face) else {
            return measure();
        };
        key.extend_from_slice(format!("{deflection:?}{tol:?}").as_bytes());
        let held = self
            .held
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .get(&key)
            .cloned();
        if let Some(found) = held {
            return Ok(found);
        }
        let found = measure()?;
        self.held
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .insert(key, found);
        Ok(found)
    }
}

impl Planner<'_> {
    pub(super) fn snapped(
        &self,
        chain: &[u32],
        faces: &[usize],
    ) -> Option<(Snapped, bool, Images)> {
        let read = self.snap_read(faces);
        let key = (chain.to_vec(), faces.to_vec());
        if let Some(known) = self.snaps_held().get(&key)
            && let Some((_, found)) = known.iter().find(|(was, _)| *was == read)
        {
            return found.clone();
        }
        let found = self.snap(chain, faces);
        self.snaps_held()
            .entry(key)
            .or_default()
            .push((read, found.clone()));
        found
    }

    /// What a snapped curve between `faces` reads of them.
    fn snap_read(&self, faces: &[usize]) -> Vec<SnapFace> {
        faces
            .iter()
            .map(|&g| SnapFace::of(&self.groups.carriers[g]))
            .collect()
    }

    /// [`image_on`] for an edge on a curved face, asked once per distinct
    /// question across the plans of one conversion.
    pub(super) fn image(
        &self,
        curved: &Curved,
        surface: &ogeom_geom::SurfaceGeometry,
        spec: &EdgeSpec,
        reach: f64,
    ) -> Option<(PlanarCurve, f64)> {
        // Everything the answer reads, written out: a float's debug form
        // reads back to the same bits, so two keys match only where every
        // input does.
        let key = format!(
            "{:?}|{:?}|{:?}|{:?}|{:?}|{:x}|{:x}|{:x}",
            curved.shape,
            curved.centre,
            curved.patch,
            surface,
            spec.curve,
            spec.range.0.to_bits(),
            spec.range.1.to_bits(),
            reach.to_bits(),
        );
        let held = |images: &std::sync::Mutex<ImageCache>| {
            images
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .get(&key)
                .cloned()
        };
        if let Some(found) = held(self.images) {
            return found;
        }
        let found = image_on(curved, surface, &spec.curve, spec.range, reach, self.tol);
        self.images
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .insert(key, found.clone());
        found
    }

    /// The curves snapped so far.
    fn snaps_held(&self) -> std::sync::MutexGuard<'_, SnapCache> {
        self.snaps
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Whether the curve for `chain` between `faces` as they stand is known.
    pub(super) fn snap_known(&self, chain: &[u32], faces: &[usize]) -> bool {
        let read = self.snap_read(faces);
        self.snaps_held()
            .get(&(chain.to_vec(), faces.to_vec()))
            .is_some_and(|known| known.iter().any(|(was, _)| *was == read))
    }

    /// The curves for chains between faces, each solved on its own, kept.
    pub(super) fn snap_all(&self, wanted: Vec<(Vec<u32>, Vec<usize>)>) {
        let found =
            ogeom_core::parallel::map_ordered(&wanted, |_, (chain, faces)| self.snap(chain, faces));
        for ((chain, faces), found) in wanted.into_iter().zip(found) {
            let read = self.snap_read(&faces);
            self.snaps_held()
                .entry((chain, faces))
                .or_default()
                .push((read, found));
        }
    }
}

/// What a snapped curve reads of one face beside its chain: the plane, or
/// the curved surface with the chart point its images unwrap around and
/// whether it wraps.
#[derive(PartialEq)]
pub(super) enum SnapFace {
    Plane(Plane),
    Curved(Canonical, (f64, f64), bool, bool),
    Gone,
}

impl SnapFace {
    fn of(carrier: &Carrier) -> Self {
        match carrier {
            Carrier::Plane(plane) => Self::Plane(*plane),
            Carrier::Curved(c) => Self::Curved(c.shape.clone(), c.centre, c.wraps, c.wraps_v),
            Carrier::Gone => Self::Gone,
        }
    }
}

/// The curves snapped to chains, kept across the plans of one conversion:
/// by chain and faces, each with what it read of the faces. The coplanar
/// distance and the points are those of the whole conversion.
pub(super) type SnapCache =
    FastMap<(Vec<u32>, Vec<usize>), Vec<(Vec<SnapFace>, Option<(Snapped, bool, Images)>)>>;

/// Edge images on curved faces, kept across the plans of one conversion:
/// by every input of [`image_on`] written out, with what it answered.
pub(super) type ImageCache = FastMap<String, Option<(PlanarCurve, f64)>>;
