//! The filler: the sections between two gathered solids, and the section
//! curves it splits, windows and marches.

use super::*;

/// The sections between two gathered solids, with the paves they put on
/// boundary edges.
#[allow(clippy::type_complexity)]
pub(crate) fn fill(
    ga: &GSolid,
    gb: &GSolid,
    admit_all: bool,
    neighbours: Option<&[ogeom_core::FastSet<usize>]>,
    tol: Tolerances,
) -> OgeomResult<(
    Vec<SectionRec>,
    Vec<SectionPiece>,
    Vec<ContactRec>,
    Vec<TangentRec>,
    Vec<Vec<(f64, f64)>>,
    ogeom_core::FastMap<EdgeKey, Vec<Pave>>,
    Vec<Vec<usize>>,
    Vec<Vec<usize>>,
    Vec<Junction>,
)> {
    use ogeom_intersect::{
        CurveCurveOptions, IntersectOptions, SurfaceIntersection, intersect_curves,
        intersect_surfaces_within,
    };
    let mut sections: Vec<SectionRec> = Vec::new();
    let mut contacts: Vec<ContactRec> = Vec::new();
    let mut tangents: Vec<TangentRec> = Vec::new();
    let mut same_pairs: Vec<(usize, usize)> = Vec::new();
    // Marched sections are fitted. The fit is driven below the confusion
    // tolerance so a fitted curve meets edges, vertices and the mesh welder
    // on the same terms as an exact one. The budget each still carries is
    // recorded per section and widens the crossing filters.
    // Every pair of faces is its own question: each reads the two solids
    // and writes only its own findings, so the pairs are asked in parallel
    // and their findings joined in pair order, the order a single thread
    // would have produced them in.
    #[derive(Default)]
    struct PairFound {
        sections: Vec<SectionRec>,
        contacts: Vec<ContactRec>,
        tangents: Vec<TangentRec>,
        same_pairs: Vec<(usize, usize)>,
        rim_crossed: bool,
        #[cfg(feature = "test-hooks")]
        marched: usize,
    }
    let at_rims = AT_RIMS.get();
    let as_stated = AS_STATED.get();
    let along_edges = along_edges::Noted::read();
    let mut pairs: Vec<(usize, usize)> = if admit_all {
        (0..ga.faces.len())
            .flat_map(|ia| (0..gb.faces.len()).map(move |ib| (ia, ib)))
            .collect()
    } else {
        // The second solid's face boxes in a tree, asked once per face of
        // the first: the pairs whose boxes meet, in the order every pair
        // tested in turn would find them.
        let filters: Vec<ogeom_math::Aabb> = gb.faces.iter().map(|f| f.filter).collect();
        let tree = box_tree::BoxTree::new(&filters);
        let mut meeting = Vec::new();
        let mut pairs = Vec::new();
        for (ia, fa) in ga.faces.iter().enumerate() {
            tree.meeting(&fa.filter, &mut meeting);
            pairs.extend(meeting.iter().map(|&ib| (ia, ib)));
        }
        pairs
    };
    // A solid against itself: a face meets its neighbours along the edges
    // and corners they share, which are no crossing, and itself everywhere.
    if let Some(neighbours) = neighbours {
        pairs.retain(|(ia, ib)| !neighbours[*ia].contains(ib));
    }
    #[cfg(feature = "test-hooks")]
    work::add(|w| w.face_pairs += pairs.len());
    // Each surface is marched over its faces' parameter boxes, not over
    // its whole stored domain, and a pair of surfaces is intersected once
    // for every pair of faces standing on them at the same chord.
    let (carriers_a, carriers_b) = (carrier_windows(ga, tol), carrier_windows(gb, tol));
    type Met = std::sync::Arc<std::sync::OnceLock<OgeomResult<SurfaceIntersection>>>;
    let met_once: std::sync::Mutex<ogeom_core::FastMap<(usize, usize, u64), Met>> =
        std::sync::Mutex::default();
    let found = ogeom_core::parallel::map_ordered(
        &pairs,
        |_, &(ia, ib)| -> OgeomResult<PairFound> {
            ogeom_core::progress::checkpoint()?;
            let (fa, fb) = (&ga.faces[ia], &gb.faces[ib]);
            let ((carrier_a, window_a), (carrier_b, window_b)) = (carriers_a[ia], carriers_b[ib]);
            let admitted = fa.filter.intersects(&fb.filter);
            let mut out = PairFound::default();
            let scale = fa.chord_scale.min(fb.chord_scale);
            let chord = (scale * 1e-7).max(tol.confusion() * 0.5);
            // A spline extrusion may carry a long profile while only a tiny
            // part meets this face. Its box is not an accuracy budget for
            // the section: that section can end on an exact planar cap,
            // whose vertex must not drift with the profile's overall size.
            let chord = if matches!(fa.surface, SurfaceGeometry::Extrusion(_))
                || matches!(fb.surface, SurfaceGeometry::Extrusion(_))
            {
                chord.min(tol.confusion())
            } else {
                chord
            };
            let options = IntersectOptions {
                // The fit budget scales with the chord: a section against a
                // fitted surface cannot honestly land closer than the
                // geometry it cuts, and whatever it carries is stated on
                // the record and widens every filter downstream.
                tolerance: chord,
                marching: ogeom_intersect::Marching {
                    chord,
                    ..ogeom_intersect::Marching::default()
                },
            };
            // A plane across a B-spline extrusion's rim, split there in
            // closed form on the run that asks for it, and noted on any
            // other.
            let at_rim = if at_rims {
                ogeom_intersect::plane_extrusion_at_rims(&fa.surface, &fb.surface, tol)
            } else {
                out.rim_crossed =
                    ogeom_intersect::crosses_extrusion_rim(&fa.surface, &fb.surface, tol);
                None
            };
            // The intersector measures a pair without a closed form for
            // coincidence at the chord; the faces' own stated looseness is
            // asked here, since only the faces know it.
            // A blend's band and a face it was built against: they meet
            // along the band's rail, which the wedge's leg on that face
            // splits it along.
            let noted = fa.carrier.1.is_identity()
                && fb.carrier.1.is_identity()
                && along_edges.holds(
                    (&fa.surface, &fb.surface),
                    (fa.carrier.0, fb.carrier.0),
                    tol,
                );
            let met = if noted {
                SurfaceIntersection::Apart
            } else if coincide_as_stated(fa, fb, as_stated, tol) {
                SurfaceIntersection::Same
            } else if let Some(found) = at_rim {
                if found.is_empty() {
                    SurfaceIntersection::Apart
                } else {
                    SurfaceIntersection::Along(found)
                }
            } else {
                let once = met_once
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner)
                    .entry((carrier_a, carrier_b, chord.to_bits()))
                    .or_default()
                    .clone();
                let met = once
                    .get_or_init(|| {
                        intersect_surfaces_within(
                            &fa.surface,
                            &fb.surface,
                            (window_a, window_b),
                            options,
                            tol,
                        )
                    })
                    .clone()?;
                let wide_a = reaching(&fa.surface, &fb.surface, &met, tol)?;
                let wide_b = reaching(&fb.surface, &fa.surface, &met, tol)?;
                if wide_a.is_some() || wide_b.is_some() {
                    // A stretched cylinder is walked whole: its faces'
                    // window would cut the stretch back off.
                    intersect_surfaces_within(
                        wide_a.as_ref().unwrap_or(&fa.surface),
                        wide_b.as_ref().unwrap_or(&fb.surface),
                        (
                            window_a.filter(|_| wide_a.is_none()),
                            window_b.filter(|_| wide_b.is_none()),
                        ),
                        options,
                        tol,
                    )?
                } else {
                    met
                }
            };
            match met {
                SurfaceIntersection::Apart => {}
                SurfaceIntersection::Same => {
                    // Coincident surfaces offer no section curve. What splits
                    // each face is the *other* face's boundary. Exact pcurve
                    // projection carries an edge into the other chart, and
                    // planes always have one. A curved pair's edge without
                    // one is carried across coaxial charts, read from the
                    // owner's own pcurve on the identical chart, or fitted
                    // by projection, and refused where none of those holds.
                    out.same_pairs.push((ia, ib));
                    for (owner_from_a, owner, target_from_a, target, target_face) in
                        [(false, fb, true, fa, ia), (true, fa, false, fb, ib)]
                    {
                        for e in &owner.edges {
                            // An edge off the target's surface by no more
                            // than the two faces' stated tolerances lies on
                            // it, as the contact itself says: its closed-form
                            // image is read at that looseness when the strict
                            // one finds it a hair off (a line an edge's slop
                            // from a cylinder it runs along). The target's
                            // looseness is its edges' widest.
                            let slop = target
                                .edges
                                .iter()
                                .map(|t| t.tolerance)
                                .fold(e.tolerance, f64::max)
                                .max(owner.tolerance + target.tolerance);
                            let exact =
                                ogeom_intersect::exact_pcurve_of(&e.curve, &target.surface, tol)
                                    .or_else(|| {
                                        if slop <= tol.confusion() {
                                            return None;
                                        }
                                        let loose = Tolerances::with_scale(
                                            tol.scale() * tol.confusion() / slop,
                                        )
                                        .ok()?;
                                        ogeom_intersect::exact_pcurve_of(
                                            &e.curve,
                                            &target.surface,
                                            loose,
                                        )
                                    });
                            let mut clipped_contact = if exact.is_none()
                                && !same_chart(&owner.surface, &target.surface, tol)
                                && carried_across(
                                    &e.pcurve,
                                    e.prange,
                                    &owner.surface,
                                    &target.surface,
                                    slop,
                                    tol,
                                )
                                .is_none()
                                && projected_into_shared_chart(
                                    &e.curve,
                                    e.crange,
                                    &owner.surface,
                                    &target.surface,
                                    tol,
                                )
                                .is_err()
                            {
                                let reach = slop
                                    .max(owner.tolerance + target.tolerance)
                                    .max(tol.confusion() * 1e2);
                                let stretches = clipped_into_chart(
                                    &e.curve,
                                    e.crange,
                                    &target.surface,
                                    reach,
                                    tol,
                                )?;
                                // Wholly off the target's window is wholly
                                // off the target.
                                if stretches.is_empty() {
                                    continue;
                                }
                                stretches
                            } else {
                                Vec::new()
                            };
                            let clipped = core::mem::take(&mut clipped_contact);
                            let shared_chart = exact.is_none()
                                && clipped.is_empty()
                                && same_chart(&owner.surface, &target.surface, tol);
                            // One cylinder on two frames whose fit the
                            // target's window cannot hold (a converted
                            // part's short wall against a pad's long one):
                            // the owner's own pcurve carried across the
                            // charts exactly.
                            let carried = if exact.is_none() && clipped.is_empty() && !shared_chart
                            {
                                carried_across(
                                    &e.pcurve,
                                    e.prange,
                                    &owner.surface,
                                    &target.surface,
                                    slop,
                                    tol,
                                )
                                .filter(|_| {
                                    projected_into_shared_chart(
                                        &e.curve,
                                        e.crange,
                                        &owner.surface,
                                        &target.surface,
                                        tol,
                                    )
                                    .is_err()
                                })
                            } else {
                                None
                            };
                            let images: Vec<(PlanarCurve, (f64, f64), (f64, f64))> = match (
                                exact,
                                clipped.is_empty(),
                                carried,
                            ) {
                                (Some(exact), _, _) => vec![(exact, e.crange, e.crange)],
                                // What the other branches cannot image whole
                                // but only runs on past the target's window
                                // is imaged over each stretch it has on the
                                // window, each its own contact.
                                (None, false, _) => clipped
                                    .into_iter()
                                    .map(|(pcurve, range)| (pcurve, range, range))
                                    .collect(),
                                // A fitted edge has no closed-form projection,
                                // but when the two faces sit on the
                                // *identical chart*, which is exactly the
                                // situation `Same` names for the analytics,
                                // the owner's own stored pcurve already is
                                // the projection, attached at construction,
                                // and it travels with its own window the way
                                // every stored pcurve does.
                                (None, true, _) if shared_chart => {
                                    vec![(e.pcurve.clone(), e.prange, e.crange)]
                                }
                                (None, true, Some(carried)) => vec![(carried, e.prange, e.crange)],
                                // Two patches that coincide as point sets
                                // without being one chart (a blend's leg
                                // on a spline host continued past the
                                // face, against the face's own patch) get
                                // the edge fitted by projection into the
                                // target's chart, same-parameter with the
                                // edge, the way every reader derives a
                                // pcurve it was not given.
                                (None, true, None) => {
                                    let Some(pcurve) = projected_into_shared_chart(
                                        &e.curve,
                                        e.crange,
                                        &owner.surface,
                                        &target.surface,
                                        tol,
                                    )?
                                    else {
                                        // Wholly outside the target's window
                                        // is wholly outside the target: a
                                        // blend's run-out past the face it
                                        // melts with splits nothing there.
                                        continue;
                                    };
                                    ogeom_core::diag!(
                                        Channel::Wire,
                                        "SAME owner {} face {} (from_a {owner_from_a}) edge {} of kind {:?}: {:?} .. {:?} against target face {target_face}",
                                        if owner_from_a { "a" } else { "b" },
                                        if owner_from_a { ia } else { ib },
                                        e.node.index(),
                                        core::mem::discriminant(&*e.curve),
                                        e.curve.point_at(e.crange.0, tol).ok(),
                                        e.curve.point_at(e.crange.1, tol).ok()
                                    );
                                    vec![(pcurve, e.crange, e.crange)]
                                }
                            };
                            for (pcurve, prange, crange) in images {
                                out.contacts.push(ContactRec {
                                    curve: (*e.curve).clone(),
                                    crange,
                                    pcurve,
                                    prange,
                                    node: e.node,
                                    // The contact stands off the target as
                                    // far as the two faces say they may stand
                                    // off their own surfaces, not only its
                                    // edge.
                                    tolerance: e.tolerance.max(owner.tolerance + target.tolerance),
                                    target_from_a,
                                    target_face,
                                    bound: e.bound,
                                    crossed: false,
                                });
                            }
                        }
                    }
                }
                // A touch at a point bounds nothing: no curve to split a
                // face along, no side that is inside on one hand and
                // outside on the other. It contributes to the arrangement
                // exactly what a tangential contact curve does, which is
                // nothing, and the result carries the touch as the
                // non-manifold contact it is.
                SurfaceIntersection::Touching(_) => {}
                SurfaceIntersection::Along(curves) => {
                    // Two planes all but parallel cross along a line, and on
                    // either side of it stay within tolerance of each other
                    // for a stretch: a block stacked on another and tilted a
                    // hair. An edge of either lying on the other's plane
                    // there splits the other's face as a coincident pair's
                    // edge would, or the two faces walk different splits of
                    // the stretch they share.
                    if let (SurfaceGeometry::Plane(pa), SurfaceGeometry::Plane(pb)) =
                        (&fa.surface, &fb.surface)
                    {
                        let (na, nb) = (pa.plane().normal(), pb.plane().normal());
                        if na.vector().cross(nb.vector()).magnitude() <= 1e-3 {
                            for (owner, target, target_from_a, target_face) in
                                [(fb, fa, true, ia), (fa, fb, false, ib)]
                            {
                                let SurfaceGeometry::Plane(on) = &target.surface else {
                                    continue;
                                };
                                for e in &owner.edges {
                                    if !matches!(&*e.curve, Curve::Line(_)) {
                                        continue;
                                    }
                                    let lies = [e.crange.0, e.crange.1].iter().all(|&t| {
                                        e.curve.point_at(t, tol).is_ok_and(|p| {
                                            on.plane().distance_to(p) <= tol.confusion()
                                        })
                                    });
                                    if !lies {
                                        continue;
                                    }
                                    let Some(pcurve) = ogeom_intersect::exact_pcurve_of(
                                        &e.curve,
                                        &target.surface,
                                        tol,
                                    ) else {
                                        continue;
                                    };
                                    out.contacts.push(ContactRec {
                                        curve: (*e.curve).clone(),
                                        crange: e.crange,
                                        pcurve,
                                        prange: e.crange,
                                        node: e.node,
                                        tolerance: e
                                            .tolerance
                                            .max(owner.tolerance + target.tolerance),
                                        target_from_a,
                                        target_face,
                                        bound: e.bound,
                                        crossed: false,
                                    });
                                }
                            }
                        }
                    }
                    // Two planes at a small angle that share an edge only to
                    // their points' rounding (a pad's wall standing on a
                    // facet just under the top edge they both hold) cross on
                    // a line the solve places by that rounding over the
                    // angle's sine: microns of slop, tens of microns along.
                    // An edge of either lying on the other's plane to within
                    // the weld distance, and off the solved line by more than
                    // the confusion distance yet no more than that slop
                    // allows, is where they cross: it splits the other face
                    // as a contact, and the solved line, a sliver away from
                    // it, splits nothing. The angle may be well under a
                    // thousandth of a radian (a facet leaning off a wall by
                    // less than its chord sags) where the edge is one both
                    // faces hold and the planes part by more than the weld
                    // distance across the smaller face: the band where they
                    // stand within it of each other runs from the shared
                    // edge, and the face past it lies off the other plane.
                    // Where they never part by that much the faces are one
                    // surface to within the weld, and a contact along one
                    // edge does not settle them.
                    let mut replaced = false;
                    if let (SurfaceGeometry::Plane(pa), SurfaceGeometry::Plane(pb)) =
                        (&fa.surface, &fb.surface)
                    {
                        let sine = pa
                            .plane()
                            .normal()
                            .vector()
                            .cross(pb.plane().normal().vector())
                            .magnitude();
                        let lines: Vec<&Curve> = curves
                            .iter()
                            .filter(|c| !c.tangential)
                            .map(|c| &c.curve)
                            .collect();
                        let weld = tol.confusion() * 1e2;
                        let narrow = sine * fa.bound.diagonal().min(fb.bound.diagonal()) > weld;
                        if (sine > 1e-3 || narrow)
                            && let [solved] = lines.as_slice()
                            && matches!(solved, Curve::Line(_))
                        {
                            for (owner, target, target_from_a, target_face) in
                                [(fb, fa, true, ia), (fa, fb, false, ib)]
                            {
                                let SurfaceGeometry::Plane(on) = &target.surface else {
                                    continue;
                                };
                                for e in &owner.edges {
                                    if !matches!(&*e.curve, Curve::Line(_)) {
                                        continue;
                                    }
                                    let mut ends = Vec::with_capacity(2);
                                    for t in [e.crange.0, e.crange.1] {
                                        ends.push(e.curve.point_at(t, tol)?);
                                    }
                                    let stray = ends
                                        .iter()
                                        .map(|p| on.plane().distance_to(*p))
                                        .fold(0.0_f64, f64::max);
                                    if stray > weld {
                                        continue;
                                    }
                                    let mut off = 0.0_f64;
                                    for p in &ends {
                                        off = off.max(
                                            ogeom_algo::project_on_curve(solved, *p, 64, tol)?
                                                .distance,
                                        );
                                    }
                                    // The solve's own slop, the edge's
                                    // stray magnified by the shallow angle;
                                    // a line no further from the edge than
                                    // the edge from the plane is a crossing
                                    // the solve placed well.
                                    if off <= tol.confusion().max(stray * 2.0) || off > weld / sine
                                    {
                                        continue;
                                    }
                                    let Some(pcurve) = ogeom_intersect::exact_pcurve_of(
                                        &e.curve,
                                        &target.surface,
                                        tol,
                                    ) else {
                                        continue;
                                    };
                                    // Under a thousandth of a radian the
                                    // edge must be one the two faces share:
                                    // an edge of the other face beside it
                                    // end to end, to the weld distance. A
                                    // face whose edge only lies on the
                                    // other's plane there (a block stacked
                                    // on another tilted a hair) truly
                                    // crosses it on the solved line.
                                    if sine <= 1e-3 && !shared_line(&ends, target, weld, tol)? {
                                        continue;
                                    }
                                    replaced = true;
                                    // All but parallel, an edge on the
                                    // other's plane to the confusion
                                    // distance is a contact already.
                                    if sine <= 1e-3 && stray <= tol.confusion() {
                                        continue;
                                    }
                                    out.contacts.push(ContactRec {
                                        curve: (*e.curve).clone(),
                                        crange: e.crange,
                                        pcurve,
                                        prange: e.crange,
                                        node: e.node,
                                        tolerance: e
                                            .tolerance
                                            .max(owner.tolerance + target.tolerance),
                                        target_from_a,
                                        target_face,
                                        bound: e.bound,
                                        crossed: false,
                                    });
                                }
                            }
                        }
                    }
                    // A face whose outline lies on one side of the other's
                    // plane, to rounding, does not cross it, and their
                    // planes' line splits neither face. Where the two stand
                    // all but parallel (a corner facet meeting a wall only at
                    // a vertex and leaning off it by less than its points'
                    // rounding over its height) the solve places that line by
                    // the rounding over the angle's sine, across the facet a
                    // tenth of a millimetre from the vertex. At a wider angle
                    // a face touching the plane only at corners gives a line
                    // along its own edge to within the rounding (a facet's
                    // side edge beside the next wall), a stub the arrangement
                    // cannot tell from the edge. An edge lying in the plane is
                    // where the faces meet: all but parallel the edge is a
                    // contact already (above), and at a wider angle the line
                    // stays.
                    if !replaced
                        && let (SurfaceGeometry::Plane(pa), SurfaceGeometry::Plane(pb)) =
                            (&fa.surface, &fb.surface)
                    {
                        let sine = pa
                            .plane()
                            .normal()
                            .vector()
                            .cross(pb.plane().normal().vector())
                            .magnitude();
                        for (face, plane) in [(fa, pb.plane()), (fb, pa.plane())] {
                            if let Some(edge_in) = one_side_of(face, plane, tol)?
                                && (sine <= 1e-3 || !edge_in)
                            {
                                replaced = true;
                            }
                        }
                    }
                    // Where one curve has no chart image, the pair is
                    // marched whole, every branch at once: the sections its
                    // other curves already gave are replaced by the marched
                    // set, and the curves after it add no more, while its
                    // tangential contacts are still kept.
                    let before = out.sections.len();
                    let mut marched = false;
                    for sc in curves {
                        // A tangential curve is contact, not crossing: the
                        // two faces meet along it and neither passes
                        // through the other, so it is set aside from the
                        // arrangement entirely and kept for the consumers
                        // that draw contact rather than classify by it.
                        if sc.tangential {
                            if let (Some(pa), Some(pb)) = (sc.on_a, sc.on_b) {
                                out.tangents.push(TangentRec {
                                    curve: sc.curve,
                                    pc_a: pa,
                                    pc_b: pb,
                                    face_a: ia,
                                    face_b: ib,
                                });
                            }
                            continue;
                        }
                        if marched || replaced {
                            continue;
                        }
                        // A section can be no longer than a turn round the
                        // faces it cuts or round its own extent: a marched
                        // trace that wanders off beside a chart's pole comes
                        // back many times the length of the circle it stands
                        // for, faithfully fitted. That is no section, and is
                        // refused by name. A whole loop of two surfaces
                        // meeting all the way round, faces cutting only an
                        // arc of it, is one turn round its own extent.
                        if !sc.exact {
                            let (length, reach) = length_and_reach(&sc.curve, tol)?;
                            let turn =
                                4.0 * fa.bound.diagonal().max(fb.bound.diagonal()).max(reach);
                            ogeom_core::diag!(
                                Channel::Wire,
                                "SECTION CHECK faces {ia}/{ib} exact {} closed {} length {length:.4} turn {turn:.4} tol {:.2e}",
                                sc.exact,
                                sc.closed,
                                sc.tolerance
                            );
                            if length > turn {
                                ogeom_bail!(
                                    NotDone,
                                    "a marched section of length {length} runs beyond a turn \
                                 round the faces it cuts ({turn}); its trace wandered \
                                 beside a chart's pole; see docs/CAPABILITIES.md, bool.booleans"
                                );
                            }
                        }
                        match (sc.on_a.clone(), sc.on_b.clone()) {
                            (Some(pa), Some(pb)) => {
                                // A marched section's ends snapped onto a
                                // vertex both surfaces pass through, where
                                // the marcher stopped short of it.
                                let (curve, pa, pb, moved) = if sc.exact || sc.closed {
                                    (sc.curve, pa, pb, 0.0)
                                } else {
                                    snapped_to_vertices(fa, fb, sc.curve, pa, pb, tol)?
                                };
                                out.sections.push(SectionRec {
                                    curve,
                                    pc_a: pa,
                                    pc_b: pb,
                                    face_a: ia,
                                    face_b: ib,
                                    closed: sc.closed,
                                    tolerance: sc.tolerance.max(moved),
                                });
                            }
                            _ => {
                                // A section running through a chart
                                // degeneracy (a plane cutting a ball on its
                                // own axis meets it at both poles) has no
                                // single chart image, because the longitude
                                // jumps half a turn there. Each piece
                                // *between* the poles does have one, and it
                                // is exact. Split first. March only if that
                                // fails.
                                if let Some(split) = split_at_degeneracies(&sc.curve, fa, fb, tol)?
                                {
                                    for (curve, pa, pb) in split {
                                        out.sections.push(SectionRec {
                                            curve,
                                            pc_a: pa,
                                            pc_b: pb,
                                            face_a: ia,
                                            face_b: ib,
                                            closed: false,
                                            tolerance: 0.0,
                                        });
                                    }
                                    continue;
                                }
                                // An exact curve whose projection has no
                                // closed form keeps its exactness and has
                                // its chart image fitted: the curve's own
                                // points inverted on the surface, read in
                                // the chart, fitted at the curve's own
                                // parameters: a plane's circle passing
                                // beside a sphere chart's pole, whose image
                                // swings fast but is a curve all the same.
                                // Marching the pair instead wanders beside
                                // the pole.
                                let image = |surface: &SurfaceGeometry,
                                         have: Option<&PlanarCurve>|
                             -> OgeomResult<Option<(PlanarCurve, f64)>> {
                                if let Some(pc) = have {
                                    return Ok(Some((pc.clone(), 0.0)));
                                }
                                fitted_image(&sc.curve, surface, options.tolerance, tol)
                            };
                                if let (Some((pa, ea)), Some((pb, eb))) = (
                                    image(&fa.surface, sc.on_a.as_ref())?,
                                    image(&fb.surface, sc.on_b.as_ref())?,
                                ) {
                                    out.sections.push(SectionRec {
                                        curve: sc.curve,
                                        pc_a: pa,
                                        pc_b: pb,
                                        face_a: ia,
                                        face_b: ib,
                                        closed: sc.closed,
                                        tolerance: sc.tolerance.max(ea).max(eb),
                                    });
                                    continue;
                                }
                                // No image fits: march the pair instead, so
                                // curve and pcurves are fitted *together*.
                                ogeom_core::diag!(
                                    Channel::Wire,
                                    "MARCH PAIR faces {ia}/{ib} a {:?} b {:?} admitted {admitted}",
                                    core::mem::discriminant(&fa.surface),
                                    core::mem::discriminant(&fb.surface)
                                );
                                out.sections.truncate(before);
                                marched = true;
                                let shared = if admitted {
                                    fa.bound.intersection(&fb.bound)
                                } else {
                                    // Audit only: disjoint bounds have no
                                    // window, and an empty window would mask
                                    // the very miss being hunted.
                                    fa.bound.union(&fb.bound)
                                };
                                let (wa, wb) = (
                                    windowed_to(&fa.surface, &shared),
                                    windowed_to(&fb.surface, &shared),
                                );
                                #[cfg(feature = "test-hooks")]
                                {
                                    out.marched += 1;
                                }
                                let (marched, pole_reach) = march_pair(&wa, &wb, &options, tol)?;
                                for fitted in marched {
                                    // The same two refusals the intersector's
                                    // own marched sections meet: a fit a
                                    // thousand chords off its trace, and a
                                    // trace longer than a turn round the
                                    // faces. A plane's circle passing beside
                                    // a sphere chart's pole marches both
                                    // ways, and stated as data either welds
                                    // the tool into a point.
                                    let budget =
                                        (options.marching.chord * 1e3).max(options.tolerance * 1e3);
                                    if fitted.fit_error > budget {
                                        ogeom_bail!(
                                            NotDone,
                                            "a marched section's fit misses its trace by {} \
                                         against a chord of {}; a branch passing beside a \
                                         chart's pole fits nothing yet; see docs/CAPABILITIES.md, \
                                         boolean.general",
                                            fitted.fit_error,
                                            options.marching.chord
                                        );
                                    }
                                    let curve: Curve = fitted.curve.into();
                                    let (length, reach) = length_and_reach(&curve, tol)?;
                                    let turn = 4.0
                                        * fa.bound.diagonal().max(fb.bound.diagonal()).max(reach);
                                    if length > turn {
                                        ogeom_bail!(
                                            NotDone,
                                            "a marched section of length {length} runs beyond a \
                                         turn round the faces it cuts ({turn}); its trace \
                                         wandered beside a chart's pole; see docs/CAPABILITIES.md, \
                                         boolean.general"
                                        );
                                    }
                                    out.sections.push(SectionRec {
                                        // A branch cut at a pole of either
                                        // chart ends on it twice, its images
                                        // apart there: closed in space, open
                                        // in the charts.
                                        closed: curve.is_closed(tol)
                                            && !on_a_chart_pole(&curve, [&wa, &wb], tol)?,
                                        // Cut at a pole within a reach, the
                                        // trace passes the pole up to that
                                        // reach off the truth.
                                        tolerance: (options.marching.chord + fitted.fit_error)
                                            .max(pole_reach),
                                        curve,
                                        pc_a: fitted.on_a.into(),
                                        pc_b: fitted.on_b.into(),
                                        face_a: ia,
                                        face_b: ib,
                                    });
                                }
                            }
                        }
                    }
                }
            }

            Ok(out)
        },
    );
    // An edge two faces of one solid hold, both lying on the other solid's
    // face (a top split in two by a fuse), is carried onto that face once
    // per holder; the second copy is the same strand, which the
    // arrangement would walk as a sliver between two coincident lines.
    let mut carried: ogeom_core::FastSet<(EdgeKey, bool, usize, [u64; 4])> =
        ogeom_core::FastSet::default();
    for pair in found {
        let pair = pair?;
        sections.extend(pair.sections);
        contacts.extend(pair.contacts.into_iter().filter(|c| {
            carried.insert((
                c.node,
                c.target_from_a,
                c.target_face,
                [
                    c.crange.0.to_bits(),
                    c.crange.1.to_bits(),
                    c.prange.0.to_bits(),
                    c.prange.1.to_bits(),
                ],
            ))
        }));
        tangents.extend(pair.tangents);
        same_pairs.extend(pair.same_pairs);
        if pair.rim_crossed {
            RIM_CROSSED.set(true);
        }
        #[cfg(feature = "test-hooks")]
        work::add(|w| w.marched_pairs += pair.marched);
    }
    for section in &mut sections {
        if let Some(rebased) = seamed_at_a_stop(
            section,
            &ga.faces[section.face_a],
            &gb.faces[section.face_b],
            tol,
        ) {
            *section = rebased;
        }
    }
    #[cfg(feature = "test-hooks")]
    work::add(|w| {
        w.widest_section = sections
            .iter()
            .fold(w.widest_section, |widest, s| widest.max(s.tolerance));
    });
    let grazed = graze_contacts(ga, gb, &sections, &contacts, FROM_CROSSINGS.get(), tol)?;
    contacts.extend(grazed);
    let touched = tangent_contacts(ga, gb, &tangents, &contacts, tol)?;
    contacts.extend(touched);

    // Boundary polylines per face, for the trim tests.
    let outline = |face: &GFace| -> OgeomResult<Vec<Vec<Point2>>> {
        let mut lines = Vec::new();
        for e in &face.edges {
            lines.push(pcurve_polyline_within_sag(
                &e.pcurve,
                e.prange,
                e.crange,
                e.crange,
                &face.surface,
                tol,
            )?);
            if let Some((other, orange)) = &e.other_side {
                lines.push(pcurve_polyline_within_sag(
                    other,
                    *orange,
                    e.crange,
                    e.crange,
                    &face.surface,
                    tol,
                )?);
            }
        }
        lines.extend(face.pole_rows(tol)?);
        weld_outline_ends(&mut lines, face, tol);
        bridge_outline_gaps(&mut lines);
        Ok(lines)
    };
    // Only the faces a section crosses are asked about: on a large solid
    // and a small tool that is a handful, and drawing every face's outline
    // was most of the paving's cost.
    let crossed_a: ogeom_core::FastSet<usize> = sections.iter().map(|s| s.face_a).collect();
    let crossed_b: ogeom_core::FastSet<usize> = sections.iter().map(|s| s.face_b).collect();
    // The face's own outline is drawn once and kept for the later trim
    // tests; where it cannot be drawn, drawing it again says why.
    let mut outlines_a: Vec<std::borrow::Cow<'_, [Vec<Point2>]>> =
        Vec::with_capacity(ga.faces.len());
    let mut outlines_b: Vec<std::borrow::Cow<'_, [Vec<Point2>]>> =
        Vec::with_capacity(gb.faces.len());
    for (faces, crossed, outlines) in [
        (&ga.faces, &crossed_a, &mut outlines_a),
        (&gb.faces, &crossed_b, &mut outlines_b),
    ] {
        for (i, f) in faces.iter().enumerate() {
            outlines.push(if !crossed.contains(&i) {
                std::borrow::Cow::Owned(Vec::new())
            } else if let Some(lines) = f.outline(tol) {
                std::borrow::Cow::Borrowed(lines)
            } else {
                std::borrow::Cow::Owned(outline(f)?)
            });
        }
    }

    // Crossings of each section with the boundary edges of both its faces,
    // and with every other section sharing a face.
    let mut paves: ogeom_core::FastMap<EdgeKey, Vec<Pave>> = ogeom_core::FastMap::default();
    let mut pieces: Vec<SectionPiece> = Vec::new();
    // Each section's paving depends only on the sections and the two
    // gathered solids, all read-only here, and writes nothing the next
    // section reads. So the measuring runs in parallel and the accumulating
    // runs afterwards in section order: the same split `tessellate` uses,
    // and the same reason: nothing about scheduling can reach the answer.
    type SectionWork = (
        Vec<(EdgeKey, Pave)>,
        Vec<SectionPiece>,
        Vec<Junction>,
        Option<SectionRec>,
    );
    let mut hug_junctions: Vec<Junction> = Vec::new();
    // Each section's box, for skipping pairs of sections too far apart to
    // cross. A curve whose box cannot be had stands unbounded and is never
    // skipped.
    let section_bounds: Vec<Option<ogeom_math::Aabb>> = sections
        .iter()
        .map(|s| {
            ogeom_algo::curve_bounds(&s.curve, tol)
                .ok()
                .map(|b| b.expanded(s.tolerance))
        })
        .collect();
    // Read here, on the caller's thread: the sections are paved on worker
    // threads, which do not share it.
    let at_stated = AT_STATED.get();
    let paved: Vec<OgeomResult<SectionWork>> = ogeom_core::parallel::map_ordered(
        &sections,
        |si, section: &SectionRec| {
            ogeom_core::progress::checkpoint()?;
            let mut paves: Vec<(EdgeKey, Pave)> = Vec::new();
            let mut pieces: Vec<SectionPiece> = Vec::new();
            let mut junctions: Vec<Junction> = Vec::new();
            // A fitted section meets an edge within its own budget, not within
            // rounding.
            let reach = tol.confusion().max(section.tolerance * 2.0);
            let cc = CurveCurveOptions {
                gap: reach.max(CurveCurveOptions::default().gap),
                ..CurveCurveOptions::default()
            };
            let domain = section.curve.domain();
            let mut trim_ts: Vec<f64> = Vec::new();
            // Crossings with boundary edges: side, edge, parameter on the
            // edge, parameter on the section, and how honestly the stop sits.
            let mut hits: Vec<(usize, EdgeKey, f64, f64, f64, bool, f64)> = Vec::new();
            // Spans of the section running *along* a boundary edge. The split
            // such a span would make already exists as boundary (stacked boxes'
            // perpendicular side planes meet exactly at the boxes' own edges),
            // so the span is excluded from that face's strands rather than
            // refused or duplicated.
            //
            // Which face, though, is the whole point. A plane through a bore's
            // axis meets the wall along two rulings, and one of them is the
            // wall's own seam: on the wall that ruling is boundary already, and
            // on the plane it is the curve that separates the two halves the
            // section leaves. Dropped from both, the plane keeps one region where
            // it has two, and the result does not close. So the exclusion is
            // recorded per face.
            let mut along: [Vec<(f64, f64, usize)>; 2] = [Vec::new(), Vec::new()];
            for (side, own) in [
                (0_usize, &ga.faces[section.face_a]),
                (1, &gb.faces[section.face_b]),
            ] {
                // The other face the section must also lie in: a crossing
                // with this edge outside it stops nothing that is kept, and a
                // face bounded by thousands of edges would try every one.
                let across = if side == 0 {
                    &gb.faces[section.face_b]
                } else {
                    &ga.faces[section.face_a]
                };
                for (ei, e) in own.edges.iter().enumerate() {
                    if !admit_all && !e.bound.expanded(reach).intersects(&across.bound) {
                        continue;
                    }
                    // An edge within its own radius of the face (the boundary
                    // of a merged group of near-coplanar facets, a ruling of a
                    // fitted cone, a fitted rail between a cylinder and a
                    // helical flank) stands off the surface the section lies
                    // on by as much. The face stops where it draws the edge,
                    // its pcurve on its surface, and a section on the face
                    // meets that image exactly; the edge's curve is met only
                    // where it is the face's own to within what it stands
                    // off, and a near miss by that much is the crossing. A
                    // fitted section states only its own error, so its reach
                    // takes that standoff in as well.
                    let (drawn, drange, off) = e.drawn();
                    let reach = if section.tolerance > 0.0 {
                        reach + off
                    } else {
                        reach
                    };
                    let planar = matches!(own.surface, SurfaceGeometry::Plane(_));
                    let widened = planar || section.tolerance <= 0.0;
                    let edge_reach = if widened { reach + off } else { reach };
                    let cc = CurveCurveOptions {
                        gap: cc.gap.max(edge_reach),
                        ..cc
                    };
                    let found = intersect_curves(&section.curve, drawn, cc, tol)?;
                    // Met at the edge's widened reach, an exact section that
                    // dips under an exact edge's curve by no more than that
                    // reach reads as touching it, the touch spread along the
                    // valley. Where the valley takes in one of the edge's
                    // ends, the touch is read at that end, and the stretch
                    // of the valley on the edge, which the section runs along,
                    // is lost. Met at the section's own reach the two cross
                    // where they truly meet, either side of the valley's
                    // bottom, and those crossings stand for the touch: the
                    // one past the edge's end is beyond the edge and dropped
                    // below, and the stretch from the other to the end runs
                    // along the edge, as it does where the dip is deeper
                    // than the widened reach and the crossings are found
                    // in the first place.
                    let near_an_end = |c: &ogeom_intersect::Crossing<Point>| -> OgeomResult<bool> {
                        for end in [e.crange.0, e.crange.1] {
                            if e.curve.point_at(end, tol)?.distance(c.point) <= c.reach {
                                return Ok(true);
                            }
                        }
                        Ok(false)
                    };
                    let mut crossings = found.crossings.clone();
                    if section.tolerance <= 0.0
                        && e.image.is_none()
                        && exact(&section.curve)
                        && exact(drawn)
                        && edge_reach > reach
                    {
                        let mut touches = Vec::new();
                        for c in &found.crossings {
                            if c.reach > 0.0 && near_an_end(c)? {
                                touches.push(*c);
                            }
                        }
                        if !touches.is_empty() {
                            let narrow = intersect_curves(
                                &section.curve,
                                drawn,
                                CurveCurveOptions {
                                    gap: reach.max(CurveCurveOptions::default().gap),
                                    ..cc
                                },
                                tol,
                            )?;
                            let met: Vec<ogeom_intersect::Crossing<Point>> = narrow
                                .crossings
                                .into_iter()
                                .filter(|n| {
                                    n.reach <= 0.0
                                        && n.gap <= reach
                                        && touches
                                            .iter()
                                            .any(|t| t.point.distance(n.point) <= t.reach)
                                })
                                .collect();
                            if !met.is_empty() {
                                crossings.retain(|c| !touches.contains(c));
                                crossings.extend(met);
                            }
                        }
                    }
                    for crossing in &crossings {
                        if crossing.gap > edge_reach {
                            continue;
                        }
                        // On a curved face a rail can pass within its radius
                        // of a section without crossing it, running beside
                        // it: a near miss past the section's own reach counts
                        // only where the two meet at an angle. A rail
                        // crossing at a helix's shallow lead stays within its
                        // radius for a stretch many radii long and still
                        // crosses.
                        if !planar
                            && crossing.gap > reach
                            && tangential(&section.curve, crossing.on_a, drawn, crossing.on_b, tol)?
                        {
                            continue;
                        }
                        let mut on_b =
                            onto_range(e.onto_curve(crossing.on_b), &e.curve, e.crange, tol);
                        let mut honesty = honest(section.tolerance, tol);
                        // A crossing at a boundary edge's own end *is* that end's
                        // vertex, exactly. The stop the section keeps must be the
                        // vertex's parameter on the section: the meet of two
                        // curves a fit tolerance apart lands a couple of microns
                        // off, the boundary side keeps its exact vertex, and the
                        // rebuilt wire gapes by the difference.
                        let mut on_a = crossing.on_a;
                        // The window is the fit-slop scale, not the section's
                        // own: the boundary curve may be a fitted intersection
                        // from an earlier boolean carrying a couple of microns
                        // of wobble, and a stop that misses the vertex by that
                        // much gapes the rebuilt wire by the same. A fitted
                        // section's stop misses by its own honesty once more:
                        // the trace is off the true crossing by that, and
                        // the crossing found against the edge by that again.
                        // An exact section crossing an exact edge has no such
                        // slop and lands where it crosses: within the edge's
                        // and its vertices' own radii of the end it is the
                        // end, and past them it is not. A drum lowered two
                        // microns under a plate's top crosses the wall it
                        // stands in two microns under the wall's corner, and
                        // the strip between is the wall's.
                        let exact_pair = section.tolerance <= 0.0
                            && e.image.is_none()
                            && exact(&section.curve)
                            && exact(drawn);
                        let weld = if exact_pair {
                            (reach + e.tolerance.max(e.ends_tolerance)).max(tol.confusion() * 10.0)
                        } else {
                            (reach + honest(section.tolerance, tol)).max(tol.confusion() * 1e2)
                        };
                        let mut at_end: Option<(f64, Point)> = None;
                        for end in [e.crange.0, e.crange.1] {
                            let vertex = e.curve.point_at(end, tol)?;
                            if vertex.distance(crossing.point) <= weld + crossing.reach {
                                at_end = Some((end, vertex));
                                let snapped =
                                    ogeom_algo::project_on_curve(&section.curve, vertex, 64, tol)?;
                                if snapped.distance <= weld + crossing.reach {
                                    on_a = snapped.parameter;
                                    // And the edge's own parameter is the
                                    // end's: two rim circles through one
                                    // block corner meet at a shallow angle
                                    // there, the crossing lands a few
                                    // hundredths of a millimetre along the
                                    // rim from the corner, and paved there
                                    // it splits the rim into a sliver.
                                    on_b = end;
                                }
                                break;
                            }
                        }
                        // A touch at the section's own end *is* the end. A
                        // band's section through a wall meets the wall's drum
                        // edge tangentially where the band is tangent to the
                        // drum, and the touch of two fitted curves at a
                        // shallow angle wanders along them by far more than
                        // their honesty. But the section stops there because
                        // it leaves its own face, and that stop is the
                        // junction: the crossing takes the end's parameter
                        // and the edge splits under the end itself.
                        let touch = crossing.reach > 0.0
                            || tangential(
                                &section.curve,
                                crossing.on_a,
                                drawn,
                                crossing.on_b,
                                tol,
                            )?;
                        let mut ended = false;
                        // A closed section has no ends: where it starts is
                        // only where its circle's frame happened to point.
                        if touch && !section.closed {
                            let half = (domain.1 - domain.0) * 0.5;
                            for end in [domain.0, domain.1] {
                                if (crossing.on_a - end).abs() > half {
                                    continue;
                                }
                                let tip = section
                                    .curve
                                    .point_at(at_param(end, domain, section.closed), tol)?;
                                if crossing.reach > 0.0
                                    && tip.distance(crossing.point) > crossing.reach + weld
                                {
                                    continue;
                                }
                                // The tip stands where the marcher left the
                                // other face (on a fitted rail, to the rail's
                                // own honesty), so the edge is asked at the
                                // hug width, not the section's.
                                let foot = ogeom_algo::project_on_curve(drawn, tip, 64, tol)?;
                                let width =
                                    (crossing.reach + honest(weld, tol)).max(tol.confusion() * 1e3);
                                ogeom_core::diag!(
                                    Channel::Wire,
                                    "PAVE s{si}: touch on edge {} at {:.6}: tip {:.3e} off the edge, width {width:.2e}",
                                    e.node.index(),
                                    crossing.on_a,
                                    foot.distance
                                );
                                if foot.distance <= width {
                                    ended = true;
                                    on_a = end;
                                    on_b = onto_range(
                                        e.onto_curve(foot.parameter),
                                        &e.curve,
                                        e.crange,
                                        tol,
                                    );
                                    // A touch at the edge's own end vertex
                                    // is that vertex: a section tangent to
                                    // the edge there runs along it for a
                                    // stretch either side, stops somewhere
                                    // on that stretch, and split where it
                                    // stopped the edge keeps a sliver the
                                    // face across the vertex never has.
                                    if let Some((end_b, vertex)) = at_end
                                        && tip.distance(vertex) <= width.min(tol.confusion() * 1e5)
                                    {
                                        on_b = end_b;
                                        let gap = tip.distance(vertex);
                                        if gap > tol.confusion() * 1e2 {
                                            junctions.push(Junction {
                                                at: vertex,
                                                reach: gap + tol.confusion() * 1e2,
                                                onto_vertex: true,
                                            });
                                        }
                                    }
                                    // Along the edge the touch is known only
                                    // as far as the two curves stay together:
                                    // a rail grazing a bore is met there by
                                    // the sections on the faces either side
                                    // of it a hundredth of a millimetre
                                    // apart, and those are one junction. No
                                    // further: a section touching an edge to
                                    // the fourth order stays beside it for
                                    // tenths of a millimetre, and a junction
                                    // that wide swallows the edge's vertex
                                    // and the section's other end with it.
                                    // An exact section's end is where its
                                    // closed form stops, and wanders nowhere.
                                    let wander = if section.tolerance > 0.0 {
                                        crossing.reach.min(tol.confusion() * 1e5)
                                    } else {
                                        0.0
                                    };
                                    honesty = honesty.max(foot.distance).max(wander);
                                }
                                break;
                            }
                        }
                        // A touch along the section is where the edge
                        // touches the other face's surface: the section
                        // lies in both surfaces, the edge in its own face's.
                        // The meet of the two curves has a double root
                        // there and lands anywhere in its valley, the root
                        // of the section's honesty times the radius they
                        // part at long, so the face across, meeting the
                        // section at the same point, splits it elsewhere.
                        // The touch is solved from the edge and the surface
                        // instead, and taken where it lies in the crossing's
                        // own valley with the section within its reach.
                        if touch
                            && !ended
                            && at_end.is_none()
                            && let Some(touched) = ogeom_intersect::touch_curve_surface(
                                drawn,
                                &across.surface,
                                crossing.on_b,
                                tol,
                            )
                            && touched.point.distance(crossing.point) <= crossing.reach + weld
                        {
                            let foot = ogeom_algo::project_on_curve(
                                &section.curve,
                                touched.point,
                                64,
                                tol,
                            )?;
                            if foot.distance <= reach {
                                on_a = foot.parameter;
                                on_b = onto_range(
                                    e.onto_curve(touched.on_curve),
                                    &e.curve,
                                    e.crange,
                                    tol,
                                );
                            }
                        }
                        // The curve is crossed over its whole length, the edge
                        // only over its range. A crossing past either end and
                        // farther than the weld from that end's vertex is on
                        // the curve beyond the edge: a line under a plane's
                        // edge crosses the circle of a round's end arc tangent
                        // to that edge twice, either side of the tangent point,
                        // and only one of the two is on the arc. A seam is
                        // left as it is: its circle runs on past a pole into
                        // the face it bounds, and a drill beside the pole
                        // crossing it there is stopped where the pole's own
                        // reading meets it.
                        if at_end.is_none()
                            && !ended
                            && e.seam.is_none()
                            && (on_b < e.crange.0 - tol.parametric()
                                || on_b > e.crange.1 + tol.parametric())
                            && distance_to_edge_curve(drawn, drange, crossing.point, tol)?
                                > weld + crossing.reach
                        {
                            continue;
                        }
                        // Where the edge runs all but in the other face's plane
                        // (a facet row's edge between two facets nearly
                        // coplanar with a wall that crosses them), the crossing
                        // stands from where it pierces that plane by its own
                        // distance off the plane over the sine of its lean out
                        // of it. A facet's plane misses the edge it shares with
                        // the next by its points' rounding, so the sections on
                        // the facets either side reach the edge apart by that
                        // rounding over the sine (a few tenths of a micron
                        // where the edge is a rung between two rows of facets
                        // lying all but in a wall, a micron where the
                        // rung leans out of the wall by a couple of
                        // thousandths of a radian), and those are one point;
                        // an exact pair meets with no gap and stays exact.
                        // The doubt seeds a junction that owns it, and the
                        // result's vertex there is widened to it. It runs past
                        // the micron any one tolerance is held to only where
                        // the gap is within the weld distance, rounding
                        // rather than a miss, and stops at two microns:
                        // the vertex a joined crossing leaves claims no more
                        // than two crossings each held to a micron would. A
                        // wider gap is held to the micron; crossings further
                        // apart than the doubt stay apart, and the boolean
                        // refuses rather than weld across a sliver it has
                        // not resolved.
                        if at_end.is_none()
                            && let SurfaceGeometry::Plane(plane) = &across.surface
                        {
                            let along_edge = e.curve.d1_at(on_b, tol)?;
                            let length = along_edge.magnitude();
                            if length > 0.0 {
                                let sine =
                                    along_edge.dot(plane.plane().normal().vector()).abs() / length;
                                if sine > 0.0 {
                                    let off = plane
                                        .plane()
                                        .signed_distance_to(e.curve.point_at(on_b, tol)?)
                                        .abs();
                                    let leaning = off / sine;
                                    honesty =
                                        honesty.max(if crossing.gap <= tol.confusion() * 1e2 {
                                            leaning.min(tol.confusion() * 2e4)
                                        } else {
                                            honest(leaning, tol)
                                        });
                                }
                            }
                        }
                        ogeom_core::diag!(
                            Channel::Wire,
                            "PAVE s{si}: side {side} edge {} crossing at {on_a:.6} (edge {on_b:.6}) gap {:.2e} reach {:.2e} honesty {honesty:.2e} at {:?}",
                            e.node.index(),
                            crossing.gap,
                            crossing.reach,
                            crossing.point
                        );
                        // A loose edge (a chord between a facet and the drum
                        // it was converted beside, standing off the drum by
                        // its sag) is met by the faces either side of it each
                        // where the other solid's face crosses its own image
                        // of it: up to the tolerance apart across the edge,
                        // and so that over the sine of the edge's lean out of
                        // the crossing face apart along it. A drill all but
                        // tangent to a round crosses the round's rim and the
                        // facet's chord beside it five hundredths apart.
                        // Held to ten tolerances, where the lean vanishes.
                        let spread = if at_end.is_none() {
                            lean_spread(e, on_b, &across.surface, crossing.point, tol).max(
                                touch_spread(
                                    (drawn, crossing.on_b),
                                    &across.surface,
                                    crossing.point,
                                    honest(section.tolerance, tol),
                                    tol,
                                ),
                            )
                        } else {
                            0.0
                        };
                        hits.push((side, e.node, on_b, on_a, honesty, at_end.is_none(), spread));
                    }
                    for overlap in &found.overlaps {
                        // The curves overlap. What is *boundary* is the stretch
                        // the edge actually covers. A sphere's seam and the far
                        // half of the same great circle lie on one curve, and
                        // reading the whole curve as boundary makes the meridian
                        // opposite the seam disappear, which is the octant's own
                        // edge, on a ball cut at its corner.
                        for ((lo, hi), ends) in overlap_within(overlap, drange, drawn, tol) {
                            // Where the edge's own end bounds the stretch, the
                            // stretch ends at that end's vertex, as a crossing
                            // at an edge's end does: the overlap states its
                            // correspondence by its own ends alone, and carried
                            // across in proportion the edge's end lands off the
                            // vertex along a fitted section by the curves'
                            // differing pace. A section touching the edge at
                            // its end to the fourth order lands there a few
                            // tenths of a micron off, the other crossings at
                            // the vertex fold into that stop and are dropped
                            // with it, and the edges they split stay whole.
                            // The foot is taken only within the overlap's own
                            // stretch, which the vertex, on the edge inside it,
                            // lies beside.
                            let (olo, ohi) = if overlap.on_a.0 <= overlap.on_a.1 {
                                overlap.on_a
                            } else {
                                (overlap.on_a.1, overlap.on_a.0)
                            };
                            let mut window = [lo, hi];
                            for (t, end) in window.iter_mut().zip(ends) {
                                let Some(end) = end else {
                                    continue;
                                };
                                let vertex = drawn.point_at(end, tol)?;
                                let foot =
                                    ogeom_algo::project_on_curve(&section.curve, vertex, 64, tol)?;
                                let mut f = foot.parameter;
                                if section.closed {
                                    let period = domain.1 - domain.0;
                                    f += period * ((*t - f) / period).round();
                                }
                                if foot.distance <= cc.gap
                                    && f >= olo - tol.parametric()
                                    && f <= ohi + tol.parametric()
                                {
                                    *t = f;
                                }
                            }
                            let [lo, hi] = window;
                            if hi - lo <= tol.parametric() {
                                continue;
                            }
                            ogeom_core::diag!(
                                Channel::Wire,
                                "PAVE s{si}: side {side} overlap ({lo:.6}, {hi:.6})"
                            );
                            // The edge splits where the shared stretch ends, as
                            // the section does: the stretch itself is the edge's
                            // to carry, and the section's next piece must meet
                            // the edge at a vertex the edge actually has. A band
                            // meeting a wall tangentially to the wall's own top
                            // edge hugs it for a few microns from their corner;
                            // without the split at the hug's far end the section
                            // dangles there and is pruned, and the wall never
                            // splits.
                            // The section leaves the edge at the hug's end by
                            // what the hug allowed, which no single strand's
                            // honesty covers: that junction owns the gap.
                            for t in [lo, hi] {
                                let at = section
                                    .curve
                                    .point_at(at_param(t, domain, section.closed), tol)?;
                                let foot = ogeom_algo::project_on_curve(drawn, at, 64, tol)?;
                                if foot.distance <= reach.max(tol.confusion() * 1e3) {
                                    let on_b = onto_range(
                                        e.onto_curve(foot.parameter),
                                        &e.curve,
                                        e.crange,
                                        tol,
                                    );
                                    // Drawn apart from its curve, the edge
                                    // splits there a standoff from the foot.
                                    let (split, gap) = if e.image.is_some() {
                                        let split = e.curve.point_at(on_b, tol)?;
                                        (split, at.distance(split))
                                    } else {
                                        (foot.point, foot.distance)
                                    };
                                    hits.push((
                                        side,
                                        e.node,
                                        on_b,
                                        t,
                                        honest(section.tolerance, tol).max(gap),
                                        false,
                                        0.0,
                                    ));
                                    if gap > tol.confusion() * 1e2 {
                                        junctions.push(Junction {
                                            at: split,
                                            reach: gap + tol.confusion() * 1e2,
                                            onto_vertex: false,
                                        });
                                    }
                                }
                            }
                            trim_ts.push(lo);
                            trim_ts.push(hi);
                            along[side].push((lo, hi, ei));
                        }
                    }
                }
            }
            // Inside a stretch the section runs along one of a face's edges,
            // that face's other edges cannot genuinely cross it: a boundary
            // is a simple loop, and they meet the hugged edge only at its
            // ends, which the stretch's own ends already stop at. What the
            // sampler reports there is the tangential dust of a leg touching
            // the arc it ends on: thirty near-crossings inside a micron,
            // which read as stops would shatter the section and pave the
            // leg at each.
            hits.retain(|(side, _, _, on_a, _, _, _)| {
                !along[*side].iter().any(|(lo, hi, _)| {
                    *on_a > lo + tol.parametric() && *on_a < hi - tol.parametric()
                })
            });
            #[allow(clippy::type_complexity)]
            let edge_hits: Vec<(
                EdgeKey,
                f64,
                f64,
                f64,
                Option<(bool, usize)>,
                f64,
                usize,
            )> = hits
                .iter()
                .map(|(side, node, on_b, on_a, honesty, crossed, spread)| {
                    let across = crossed.then(|| {
                        if *side == 0 {
                            (false, section.face_b)
                        } else {
                            (true, section.face_a)
                        }
                    });
                    (*node, *on_b, *on_a, *honesty, across, *spread, *side)
                })
                .collect();
            trim_ts.extend(hits.iter().map(|(_, _, _, on_a, _, _, _)| *on_a));
            let mut cross_ts: Vec<f64> = Vec::new();
            for (sj, other) in sections.iter().enumerate() {
                if sj == si {
                    continue;
                }
                if other.face_a != section.face_a && other.face_b != section.face_b {
                    continue;
                }
                let both = reach.max(tol.confusion().max(other.tolerance * 2.0));
                if !admit_all
                    && let (Some(x), Some(y)) = (&section_bounds[si], &section_bounds[sj])
                    && !x.expanded(both).intersects(y)
                {
                    continue;
                }
                // Where two sections share one face, a crossing matters only
                // inside the faces they do not share, since a section is kept
                // only where it lies in both its faces: one wall met by a
                // thousand facets carries a thousand sections, and only those
                // whose facets' bounds meet can cross where it counts. The
                // audit's unfiltered pass tries every pair still.
                if !admit_all {
                    let apart = |x: &ogeom_math::Aabb, y: &ogeom_math::Aabb| {
                        !x.expanded(both).intersects(&y.expanded(both))
                    };
                    let unshared_a = other.face_a != section.face_a
                        && apart(
                            &ga.faces[section.face_a].filter,
                            &ga.faces[other.face_a].filter,
                        );
                    let unshared_b = other.face_b != section.face_b
                        && apart(
                            &gb.faces[section.face_b].filter,
                            &gb.faces[other.face_b].filter,
                        );
                    if unshared_a || unshared_b {
                        continue;
                    }
                }
                let cc2 = CurveCurveOptions {
                    gap: both.max(CurveCurveOptions::default().gap),
                    ..CurveCurveOptions::default()
                };
                let found = intersect_curves(&section.curve, &other.curve, cc2, tol)?;
                for crossing in &found.crossings {
                    if crossing.gap > both {
                        continue;
                    }
                    let mut at = crossing.on_a;
                    // A touch at this section's end is the end, as against an
                    // edge: two sections through one wall from a band and
                    // the leg it is tangent to touch where the band's ends.
                    let touch = crossing.reach > 0.0
                        || tangential(
                            &section.curve,
                            crossing.on_a,
                            &other.curve,
                            crossing.on_b,
                            tol,
                        )?;
                    if touch {
                        let half = (domain.1 - domain.0) * 0.5;
                        for end in [domain.0, domain.1] {
                            if (crossing.on_a - end).abs() > half {
                                continue;
                            }
                            let tip = section
                                .curve
                                .point_at(at_param(end, domain, section.closed), tol)?;
                            let foot = ogeom_algo::project_on_curve(&other.curve, tip, 64, tol)?;
                            let width =
                                (crossing.reach + honest(both, tol)).max(tol.confusion() * 1e3);
                            if foot.distance <= width {
                                at = end;
                            }
                            break;
                        }
                    }
                    // Two sections through one face that lie on two different
                    // faces of the other solid cross only on the edge those
                    // two share, where each already stops. Meeting there
                    // tangentially (a wall through a plane and the fillet
                    // tangent to it), their own crossing is ill-conditioned
                    // and lands off the edge, as far off as the two stay
                    // within their doubt of each other. The edge's stop is
                    // the one where they cannot be told apart up to it.
                    let shared = |edges: &[BoundaryEdge], node: EdgeKey| {
                        edges.iter().any(|e| e.node == node)
                    };
                    let apart_on = if other.face_a == section.face_a {
                        Some(&gb.faces[other.face_b].edges)
                    } else if other.face_b == section.face_b {
                        Some(&ga.faces[other.face_a].edges)
                    } else {
                        None
                    };
                    if let Some(edges) = apart_on {
                        let close = both.max(tol.confusion() * 1e2);
                        let mut best: Option<(f64, f64)> = None;
                        for (node, _, on_a, ..) in &edge_hits {
                            if !shared(edges, *node)
                                || best.is_some_and(|(bd, _)| (on_a - at).abs() >= bd)
                            {
                                continue;
                            }
                            let mut together = true;
                            for k in 1..8 {
                                let t = (on_a - at).mul_add(f64::from(k) / 8.0, at);
                                let p = section
                                    .curve
                                    .point_at(at_param(t, domain, section.closed), tol)?;
                                if ogeom_algo::project_on_curve(&other.curve, p, 16, tol)?.distance
                                    > close
                                {
                                    together = false;
                                    break;
                                }
                            }
                            if together {
                                best = Some(((on_a - at).abs(), *on_a));
                            }
                        }
                        if let Some((_, on_a)) = best {
                            at = on_a;
                        } else if let (Some(own), Some(beside)) = if other.face_a == section.face_a
                        {
                            (gb.faces.get(section.face_b), gb.faces.get(other.face_b))
                        } else {
                            (ga.faces.get(section.face_a), ga.faces.get(other.face_a))
                        } {
                            // Off the edge the two faces share, the crossing
                            // is one only where the other section is kept:
                            // inside the other section's face, to within that
                            // face's doubt and a loose shared edge's
                            // tolerance. Beside a loose shared edge (a chord
                            // between a facet and the drum it was converted
                            // beside) the two faces' surfaces meet the
                            // crossing face again a sag's lean off the edge,
                            // on neither face, and a cut there leaves a stub
                            // inside the junction the edge's crossings make.
                            // Beside a tight one, the other section's surface
                            // can run on past its face's end: a round's
                            // cylinder meets a drum lowered under the round's
                            // tangent edge on both sides of that edge, the
                            // far branch a few thousandths off the round, and
                            // a cut there splits the drum's own top line
                            // where nothing else meets it.
                            let loose = beside
                                .edges
                                .iter()
                                .filter(|e| {
                                    e.tolerance > tol.confusion() * DRAWN_APART
                                        && own.edges.iter().any(|o| o.node == e.node)
                                })
                                .map(|e| e.tolerance)
                                .fold(None, |most: Option<f64>, t| {
                                    Some(most.map_or(t, |m| m.max(t)))
                                });
                            let doubt = doubt_at(beside, crossing.point, tol)?
                                .max(both)
                                .max(tol.confusion() * 1e3)
                                + loose.unwrap_or(0.0);
                            if !within_face(beside, crossing.point, doubt, tol)? {
                                ogeom_core::diag!(
                                    Channel::Wire,
                                    "PAVE s{si}: cross s{sj} at {at:.6} off s{sj}'s face, no cut"
                                );
                                continue;
                            }
                        }
                    }
                    ogeom_core::diag!(
                        Channel::Wire,
                        "PAVE s{si}: cross s{sj} at {at:.6} gap {:.2e}",
                        crossing.gap
                    );
                    cross_ts.push(at);
                }
            }

            // Candidate intervals between trim crossings, kept where the middle
            // sits inside both faces' trims.
            trim_ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            trim_ts.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());
            // Stops the section cannot tell apart are one stop: two
            // crossings a few nanometres apart along it would make an
            // interval of dust that one face keeps and the other drops.
            // A stop standing on the kept one (within the confusion
            // distance, a nanometre along the section and yet past the
            // parametric tolerance) is the same point: its crossings split
            // their edges at the kept stop too. One further off is dropped
            // with its crossings, which would split their edges where no
            // strand ends, unless the kept stop crosses nothing of that
            // crossing's face: the strand on that face ends at the kept stop
            // all the same, and the face's edge splits at the stop's foot.
            let mut speaks_for: Vec<(f64, f64, bool)> = Vec::with_capacity(trim_ts.len());
            {
                let floor = (honest(section.tolerance, tol) * 3.0).max(tol.confusion() * 10.0);
                let mut kept: Vec<f64> = Vec::with_capacity(trim_ts.len());
                let mut held: Option<Point> = None;
                for t in &trim_ts {
                    let at = section
                        .curve
                        .point_at(at_param(*t, domain, section.closed), tol)?;
                    if held.is_none_or(|h: Point| h.distance(at) > floor) {
                        kept.push(*t);
                        held = Some(at);
                    }
                    if let (Some(&first), Some(h)) = (kept.last(), held) {
                        speaks_for.push((*t, first, h.distance(at) > tol.confusion()));
                    }
                }
                trim_ts = kept;
            }
            // How far from an edge a stretch of the section may run and
            // still run along it, below the edge's own reach (see the hug
            // below). Wider than the crossing filters on purpose for a
            // fitted section: a tangentially-traced curve wobbles about the
            // boundary it hugs by far more than a fit budget, and a genuine
            // section keeps a distance of feature scale, not microns. An
            // exact section keeps to the edge's own honesty: a plane's
            // ellipse across a band touches the band's rail tangentially,
            // and the stretch within a tenth of a millimetre of the rail is
            // the section, not the rail. Read as a hug it leaves the cap's
            // crescent unable to close. On the last attempt a fitted
            // section, too, runs along an edge only within the sum of its
            // radius and the edge's: its tolerance is how far it may sit
            // from the true section, the shallow angle of the two surfaces
            // included, and a section leaving a rim it touches parts from it
            // by more than that well inside the loose width.
            let hug_floor = if section.tolerance > 0.0 && !at_stated {
                tol.confusion() * 1e3
            } else {
                tol.confusion() * 10.0
            };
            // A vertex of either face within that width of an exact section
            // is a stop of it. Running along a face's boundary, the section
            // passes there from one edge onto the next (a drill's top line a
            // few confusion distances under a plate's edge, on past its end
            // along the arc the edge runs into), and the hug reads a stretch
            // along one edge at a time: a piece along two reads as along
            // neither, and splits the face in a sliver within the width of
            // the face's own boundary. A section crossing a face's edges at
            // a vertex already stops there.
            let mut vertex_stops: Vec<f64> = Vec::new();
            if section.tolerance <= 0.0 {
                let width = hug_floor;
                let near_stop = |t: f64, at: Point| -> OgeomResult<bool> {
                    Ok(section
                        .curve
                        .point_at(at_param(t, domain, section.closed), tol)?
                        .distance(at)
                        <= width)
                };
                for own in [&ga.faces[section.face_a], &gb.faces[section.face_b]] {
                    for e in &own.edges {
                        for end in [e.crange.0, e.crange.1] {
                            let vertex = e.curve.point_at(end, tol)?;
                            if section_bounds[si]
                                .as_ref()
                                .is_some_and(|b| !b.expanded(width).contains(vertex))
                            {
                                continue;
                            }
                            let foot =
                                ogeom_algo::project_on_curve(&section.curve, vertex, 64, tol)?;
                            if foot.distance > width
                                || (!section.closed
                                    && (foot.parameter <= domain.0 + tol.parametric()
                                        || foot.parameter >= domain.1 - tol.parametric()))
                            {
                                continue;
                            }
                            let mut taken = false;
                            for &t in trim_ts.iter().chain(&vertex_stops) {
                                if near_stop(t, foot.point)? {
                                    taken = true;
                                    break;
                                }
                            }
                            if taken {
                                continue;
                            }
                            ogeom_core::diag!(
                                Channel::Wire,
                                "PAVE s{si}: stop at {:.6} beside the vertex {vertex:?}, {:.2e} off",
                                foot.parameter,
                                foot.distance
                            );
                            vertex_stops.push(foot.parameter);
                        }
                    }
                }
                trim_ts.extend_from_slice(&vertex_stops);
                trim_ts.sort_by(f64::total_cmp);
            }
            let mut candidates: Vec<(f64, f64)> = Vec::new();
            if section.closed {
                let period = domain.1 - domain.0;
                if trim_ts.is_empty() {
                    candidates.push(domain);
                } else {
                    for i in 0..trim_ts.len() {
                        let lo = trim_ts[i];
                        let hi = if i + 1 < trim_ts.len() {
                            trim_ts[i + 1]
                        } else {
                            trim_ts[0] + period
                        };
                        candidates.push((lo, hi));
                    }
                }
            } else {
                let mut stops = vec![domain.0];
                stops.extend(
                    trim_ts
                        .iter()
                        .copied()
                        .filter(|t| *t > domain.0 && *t < domain.1),
                );
                stops.push(domain.1);
                for pair in stops.windows(2) {
                    candidates.push((pair[0], pair[1]));
                }
            }

            // Whether the section at `t` lies inside each face's trim.
            // Whether a piece of the section lies inside each face. A planar
            // face with a curved rim is asked exactly at the piece's middle,
            // and where the middle stands too close to the rim to say (a
            // sliver a few microns wide, deepest at one end of the piece),
            // at points either side of it; the polyline answers only where
            // no point along the piece is clear of the rim.
            let inside_over = |lo: f64, hi: f64| -> OgeomResult<[bool; 2]> {
                let la: Vec<&[Point2]> = outlines_a[section.face_a]
                    .iter()
                    .map(Vec::as_slice)
                    .collect();
                let lb: Vec<&[Point2]> = outlines_b[section.face_b]
                    .iter()
                    .map(Vec::as_slice)
                    .collect();
                let face_a = &ga.faces[section.face_a];
                let face_b = &gb.faces[section.face_b];
                let at = |t: f64| -> OgeomResult<(Point2, Point2)> {
                    let tf = if section.closed { fold(t, domain) } else { t };
                    Ok((
                        fold_inside(section.pc_a.point_at(tf, tol)?, &face_a.surface, &la, tol),
                        fold_inside(section.pc_b.point_at(tf, tol)?, &face_b.surface, &lb, tol),
                    ))
                };
                let mid = f64::midpoint(lo, hi);
                let (ua, ub) = at(mid)?;
                let mut held = [
                    inside_plane_exactly(face_a, ua, tol),
                    inside_plane_exactly(face_b, ub, tol),
                ];
                for k in [3.0, 5.0, 2.0, 6.0, 1.0, 7.0] {
                    if held.iter().all(Option::is_some) {
                        break;
                    }
                    let (va, vb) = at(lo + (hi - lo) * k / 8.0)?;
                    if held[0].is_none() {
                        held[0] = inside_plane_exactly(face_a, va, tol);
                    }
                    if held[1].is_none() {
                        held[1] = inside_plane_exactly(face_b, vb, tol);
                    }
                }
                Ok([
                    held[0].unwrap_or_else(|| inside_many(&la, ua)),
                    held[1].unwrap_or_else(|| inside_many(&lb, ub)),
                ])
            };

            // Stretches running along an edge of each face, and the stops
            // where a kept piece ends.
            let mut bridges: Vec<(f64, f64, [Vec<usize>; 2])> = Vec::new();
            let mut kept_stops: Vec<f64> = Vec::new();
            let mut reseamed: Option<SectionRec> = None;
            for (lo, hi) in candidates {
                if hi - lo <= tol.parametric() {
                    continue;
                }
                let mid = f64::midpoint(lo, hi);
                let mid_folded = if section.closed {
                    fold(mid, domain)
                } else {
                    mid
                };
                let held = inside_over(lo, hi)?;
                if diag::enabled(Channel::Wire) {
                    let tf = if section.closed {
                        fold(mid, domain)
                    } else {
                        mid
                    };
                    ogeom_core::diag!(
                        Channel::Wire,
                        "PAVE s{si}: candidate ({lo:.6}, {hi:.6}) inside {held:?} at a {:?} b {:?}",
                        section
                            .pc_a
                            .point_at(tf, tol)
                            .ok()
                            .map(|q| fold_point_into_chart(
                                q,
                                &ga.faces[section.face_a].surface,
                                tol
                            )),
                        section
                            .pc_b
                            .point_at(tf, tol)
                            .ok()
                            .map(|q| fold_point_into_chart(
                                q,
                                &gb.faces[section.face_b].surface,
                                tol
                            ))
                    );
                }
                // A section that runs along a boundary edge of a face splits
                // nothing *there*: the split already exists as boundary. The
                // analytic overlap detection above catches the same-support
                // cases. This catches the rest (a fitted section tracing a
                // boundary curve, a surface meeting another exactly at its own
                // trim) by measurement rather than by recognising supports.
                // The hug is asked before the trim: a section on a face's own
                // edge reads inside or outside that face by a hair, and the
                // corner tool's block face meets a band exactly along the arc
                // that bounds it. On the hugging face it is boundary already;
                // what matters is whether the *other* face holds it.
                let mut hugs = [false; 2];
                // The edges each side hugs (a rim already in pieces is
                // several): they must split wherever the section does, and
                // the section wherever they end, or the two faces walk
                // different subdivisions of one curve.
                let mut hugged: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
                for (side, own, side_from_a, side_face) in [
                    (0_usize, &ga.faces[section.face_a], true, section.face_a),
                    (1, &gb.faces[section.face_b], false, section.face_b),
                ] {
                    if let Some((_, _, ei)) = along[side]
                        .iter()
                        .find(|(alo, ahi, _)| mid_folded >= *alo && mid_folded <= *ahi)
                    {
                        hugs[side] = true;
                        hugged[side].push(*ei);
                        continue;
                    }
                    let mut all_near = true;
                    let mut votes: Vec<usize> = vec![0; own.edges.len()];
                    let mut contact_votes: Vec<usize> = vec![0; contacts.len()];
                    let samples = (0..=4)
                        .map(|i| {
                            let t = lo + (hi - lo) * f64::from(i) / 4.0;
                            let tf = if section.closed { fold(t, domain) } else { t };
                            section.curve.point_at(tf, tol)
                        })
                        .collect::<OgeomResult<Vec<Point>>>()?;
                    // A piece shorter than the width is within it of every
                    // edge through either of its ends: a stub between two
                    // crossings near a corner, leaving one edge along the
                    // other, reads as running along both. Along an edge
                    // such a piece keeps within a quarter of its own
                    // length of it (it leans off it by less than fifteen
                    // degrees), and the width is held to that.
                    let stub: f64 = samples.windows(2).map(|w| w[0].distance(w[1])).sum();
                    let held_to = |width: f64, own: f64| -> f64 {
                        if stub < width {
                            (stub * 0.25).max(own).max(tol.confusion() * 10.0)
                        } else {
                            width
                        }
                    };
                    for (i, &at) in samples.iter().enumerate() {
                        let mut near = false;
                        for (ei, e) in own.edges.iter().enumerate() {
                            let (drawn, drange, off) = e.drawn();
                            let radii = if at_stated && section.tolerance > 0.0 {
                                section.tolerance + e.tolerance.max(e.residual)
                            } else {
                                0.0
                            };
                            let width =
                                held_to(reach.max(hug_floor).max(radii).max(off * 2.0), off * 2.0);
                            // The edge's box holds its curve: a point outside
                            // the box grown by the width is not near it.
                            if !e.bound.expanded(width).contains(at) {
                                continue;
                            }
                            let d = hug_distance(drawn, drange, at, tol)?;
                            if d <= width {
                                ogeom_core::diag!(
                                    Channel::Wire,
                                    "PAVE s{si}: side {side} sample {i} hugs edge {} at {d:.2e} (width {width:.2e}) range {:?} point {at:?}",
                                    e.node.index(),
                                    e.crange
                                );
                                // Every edge within reach is hugged, not the
                                // first found: where three rims converge on
                                // a pole corner all three come within a
                                // micron of the section, and the one it
                                // actually runs along must not lose its
                                // vote to a neighbour listed before it.
                                near = true;
                                votes[ei] += 1;
                            }
                        }
                        let width = reach.max(hug_floor);
                        // A contact edge is a strand on this face too (the other
                        // solid's boundary, carried into a chart they share), so a
                        // section tracing one would be the same curve twice, and
                        // the arrangement cannot walk a line it meets from both
                        // sides at once. Two boxes side by side put the low one's
                        // lid exactly there.
                        for (ci, c) in contacts.iter().enumerate() {
                            if c.target_from_a != side_from_a || c.target_face != side_face {
                                continue;
                            }
                            let held = held_to(width.max(c.tolerance * 2.0), c.tolerance * 2.0);
                            if c.bound.expanded(held).contains(at)
                                && hug_distance(&c.curve, c.crange, at, tol)? <= held
                            {
                                near = true;
                                contact_votes[ci] += 1;
                            }
                        }
                        if !near {
                            all_near = false;
                            break;
                        }
                    }
                    // And along one of them all the way: a section winding
                    // round a drum passes the drum's seam once a turn, and a
                    // stretch spanning whole turns puts every sample on the
                    // seam without running along it anywhere between.
                    let one_line = votes.iter().chain(&contact_votes).any(|&v| v == 5);
                    hugs[side] = all_near && one_line;
                    if all_near {
                        hugged[side].extend(
                            votes
                                .iter()
                                .enumerate()
                                .filter(|(_, c)| **c > 0)
                                .map(|(ei, _)| ei),
                        );
                    }
                }
                ogeom_core::diag!(
                    Channel::Wire,
                    "PAVE s{si}: candidate ({lo:.6}, {hi:.6}) hugs {hugs:?}"
                );
                // Inside both faces the piece is a section of both. Hugging
                // one face's own edge while that face's trim reads it out
                // (by a hair, on the edge), it is still the other face's
                // split, and is admitted there, unless that edge is already
                // carried onto the other face as a contact: then the split
                // is laid down once already, and a section would lay it
                // twice.
                let mut hug_key: Option<(EdgeKey, bool, usize)> = None;
                let admitted = if held[0] && held[1] {
                    // Inside both by a hair while running along one face's
                    // own edge, it is that edge's split of the other face,
                    // as a hug-admitted piece is: another face meeting the
                    // other along the same edge (a chamfer beside the face
                    // it bevels) lays the same split, and the two must be
                    // told apart below, or the other face is cut twice
                    // along one line.
                    if hugs[0] != hugs[1] {
                        let hugging = usize::from(hugs[1]);
                        let (own, other_from_a, other_face) = if hugging == 0 {
                            (&ga.faces[section.face_a], false, section.face_b)
                        } else {
                            (&gb.faces[section.face_b], true, section.face_a)
                        };
                        hug_key = hugged[hugging]
                            .first()
                            .map(|&ei| (own.edges[ei].node, other_from_a, other_face));
                    }
                    true
                } else if hugs[0] != hugs[1] {
                    let hugging = usize::from(hugs[1]);
                    let other = 1 - hugging;
                    let (own, other_from_a, other_face) = if hugging == 0 {
                        (&ga.faces[section.face_a], false, section.face_b)
                    } else {
                        (&gb.faces[section.face_b], true, section.face_a)
                    };
                    let carried = hugged[hugging].iter().any(|&ei| {
                        contacts.iter().any(|c| {
                            c.node == own.edges[ei].node
                                && c.target_from_a == other_from_a
                                && c.target_face == other_face
                        })
                    });
                    if held[other] && !carried {
                        hug_key = hugged[hugging]
                            .first()
                            .map(|&ei| (own.edges[ei].node, other_from_a, other_face));
                        true
                    } else {
                        false
                    }
                } else {
                    false
                };
                if hugs[0] && hugs[1] {
                    bridges.push((lo, hi, hugged.clone()));
                }
                if !admitted {
                    continue;
                }
                if hugs[0] && hugs[1] {
                    // Boundary on both sides: the split exists twice over and
                    // adding it a third time would cancel what it copies.
                    continue;
                }
                kept_stops.extend([lo, hi]);
                // Keep the paves that end a kept interval: those are where edges
                // genuinely split.
                let speaker = |on_section: f64| {
                    speaks_for
                        .iter()
                        .find(|(t, ..)| (t - on_section).abs() <= tol.parametric())
                        .map_or((on_section, false), |&(_, first, moved)| (first, moved))
                };
                for (node, on_edge, on_section, honesty, across, spread, side) in &edge_hits {
                    let (s, moved) = speaker(*on_section);
                    let near = |x: f64| {
                        (s - x).abs() <= tol.parametric()
                            || (section.closed
                                && ((s + (domain.1 - domain.0)) - x).abs() <= tol.parametric())
                    };
                    if !(near(lo) || near(hi)) {
                        continue;
                    }
                    if !moved {
                        paves.push((
                            *node,
                            Pave {
                                t: *on_edge,
                                honesty: *honesty,
                                across: *across,
                                spread: *spread,
                            },
                        ));
                        continue;
                    }
                    // The kept stop's own crossings on this face, or those
                    // standing on it, end the strand here already.
                    if edge_hits.iter().any(|(_, _, other, .., other_side)| {
                        other_side == side && speaker(*other) == (s, false)
                    }) {
                        continue;
                    }
                    let own = if *side == 0 {
                        &ga.faces[section.face_a]
                    } else {
                        &gb.faces[section.face_b]
                    };
                    let Some(e) = own.edges.iter().find(|e| e.node == *node) else {
                        continue;
                    };
                    let at = section
                        .curve
                        .point_at(at_param(s, domain, section.closed), tol)?;
                    let foot = ogeom_algo::project_on_curve(e.drawn().0, at, 64, tol)?;
                    ogeom_core::diag!(
                        Channel::Wire,
                        "PAVE s{si}: side {side} edge {} crossing at {on_section:.6} handed to the stop at {s:.6}, {:.2e} off",
                        node.index(),
                        foot.distance
                    );
                    paves.push((
                        *node,
                        Pave {
                            t: onto_range(e.onto_curve(foot.parameter), &e.curve, e.crange, tol),
                            honesty: honesty.max(foot.distance),
                            across: *across,
                            spread: *spread,
                        },
                    ));
                }
                // Split at section/section crossings inside the kept interval,
                // so every face sees the same subdivision.
                let mut cuts = vec![lo];
                // A wrap interval also splits at the curve's own domain end, so
                // every piece lives within one period and evaluates in-domain
                // after a single fold.
                if section.closed
                    && domain.1 > lo + tol.parametric()
                    && domain.1 < hi - tol.parametric()
                {
                    cuts.push(domain.1);
                }
                for &c in &cross_ts {
                    let c2 = if section.closed && c < lo {
                        c + (domain.1 - domain.0)
                    } else {
                        c
                    };
                    if c2 > lo + tol.parametric() && c2 < hi - tol.parametric() {
                        cuts.push(c2);
                    }
                }
                cuts.push(hi);
                // A hugged edge's own ends are cuts of the section too: the
                // edge may already be several pieces (a sphere's rim split
                // at its chart seam), and the section must walk the same
                // pieces, or the two faces never sew along it.
                for side in 0..2 {
                    if !hugs[side] {
                        continue;
                    }
                    let own = if side == 0 {
                        &ga.faces[section.face_a]
                    } else {
                        &gb.faces[section.face_b]
                    };
                    for &ei in &hugged[side] {
                        let e = &own.edges[ei];
                        for end in [e.crange.0, e.crange.1] {
                            let p = e.curve.point_at(end, tol)?;
                            let foot = ogeom_algo::project_on_curve(&section.curve, p, 64, tol)?;
                            ogeom_core::diag!(
                                Channel::Wire,
                                "PAVE s{si}: side {side} hugged edge {} end {end:.6} lands at {:.6} off {:.2e} against ({lo:.6}, {hi:.6})",
                                e.node.index(),
                                foot.parameter,
                                foot.distance
                            );
                            if foot.distance > reach.max(tol.confusion() * 1e3) {
                                continue;
                            }
                            let mut t = foot.parameter;
                            if section.closed {
                                let period = domain.1 - domain.0;
                                while t < lo {
                                    t += period;
                                }
                                while t - period >= lo {
                                    t -= period;
                                }
                            }
                            if t > lo + tol.parametric() && t < hi - tol.parametric() {
                                ogeom_core::diag!(
                                    Channel::Wire,
                                    "PAVE s{si}: side {side} cut at {t:.6} where hugged edge {} ends",
                                    e.node.index()
                                );
                                cuts.push(t);
                            }
                        }
                    }
                }
                cuts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
                // Cuts the section cannot tell apart are one cut: two
                // crossings a few nanometres apart along it would leave a
                // piece of dust one face keeps and the other drops. The
                // interval's own ends stay.
                {
                    let floor = (honest(section.tolerance, tol) * 3.0).max(tol.confusion() * 10.0);
                    let mut kept: Vec<f64> = Vec::with_capacity(cuts.len());
                    let mut held: Option<Point> = None;
                    for (index, c) in cuts.iter().enumerate() {
                        let at = section
                            .curve
                            .point_at(at_param(*c, domain, section.closed), tol)?;
                        let last = index + 1 == cuts.len();
                        if index == 0 || last || held.is_none_or(|h: Point| h.distance(at) > floor)
                        {
                            kept.push(*c);
                            held = Some(at);
                        }
                    }
                    cuts = kept;
                }
                // A marched section's ends weld at a thousand confusions, and
                // a piece no longer than that between the domain-end split
                // and its neighbouring cut is dust the arrangement collapses
                // into one node: the loop then closes in space within that
                // weld and its image stops the piece's length short in the
                // chart. Where nothing moves but the seam, the curve's seam
                // moves to the interval's start instead, and nothing is
                // split there.
                if section.tolerance > 0.0
                    && let Some(at) = cuts.iter().position(|c| *c == domain.1)
                    && at > 0
                    && at + 1 < cuts.len()
                {
                    let span = |from: f64, to: f64| -> OgeomResult<f64> {
                        let mut length = 0.0;
                        let mut last = section.curve.point_at(at_param(from, domain, true), tol)?;
                        for k in 1..=16 {
                            let t = from + (to - from) * f64::from(k) / 16.0;
                            let p = section.curve.point_at(at_param(t, domain, true), tol)?;
                            length += last.distance(p);
                            last = p;
                        }
                        Ok(length)
                    };
                    let shortest = span(cuts[at - 1], cuts[at])?.min(span(cuts[at], cuts[at + 1])?);
                    if shortest <= tol.confusion() * 1e3 {
                        reseamed = reseamed_section(
                            section,
                            lo,
                            [
                                &ga.faces[section.face_a].surface,
                                &gb.faces[section.face_b].surface,
                            ],
                            tol,
                        );
                        ogeom_core::diag!(
                            Channel::Wire,
                            "PAVE s{si}: the domain-end split of ({lo:.6}, {hi:.6}) leaves {shortest:.3e}; {}",
                            if reseamed.is_some() {
                                "the seam moves to the interval's start"
                            } else {
                                "it stays"
                            }
                        );
                        if reseamed.is_some() {
                            cuts.remove(at);
                        }
                    }
                }
                // A hugged edge is boundary where this piece is section: it
                // must split at the piece's interior cuts too, or the face
                // that keeps the edge walks one long strand where its
                // neighbour walks three, and the sew pairs none of them.
                for side in 0..2 {
                    if !hugs[side] {
                        continue;
                    }
                    let own = if side == 0 {
                        &ga.faces[section.face_a]
                    } else {
                        &gb.faces[section.face_b]
                    };
                    for &ei in &hugged[side] {
                        let e = &own.edges[ei];
                        // The piece's own ends as well as its interior cuts: a
                        // stop against another face's edge is a split of the
                        // hugged edge too, or the section's neighbour walks two
                        // strands where the hugged edge's face walks one.
                        for c in &cuts {
                            let at = section
                                .curve
                                .point_at(at_param(*c, domain, section.closed), tol)?;
                            let foot = ogeom_algo::project_on_curve(e.drawn().0, at, 64, tol)?;
                            // A cut off the edge (past the stretch the
                            // section hugs, on a neighbouring edge) is no
                            // split of it. Stated as one with its distance
                            // for honesty it seeds a junction a feature wide.
                            if foot.distance > reach.max(tol.confusion() * 1e3) {
                                continue;
                            }
                            let on_e =
                                onto_range(e.onto_curve(foot.parameter), &e.curve, e.crange, tol);
                            if on_e <= e.crange.0.min(e.crange.1) + tol.parametric()
                                || on_e >= e.crange.0.max(e.crange.1) - tol.parametric()
                            {
                                continue;
                            }
                            ogeom_core::diag!(
                                Channel::Wire,
                                "PAVE s{si}: side {side} hugged edge {} paved at {on_e:.6} for the cut at {c:.6}",
                                e.node.index()
                            );
                            paves.push((
                                e.node,
                                Pave {
                                    t: on_e,
                                    honesty: honest(section.tolerance, tol).max(foot.distance),
                                    across: None,
                                    spread: 0.0,
                                },
                            ));
                        }
                    }
                }
                for pair in cuts.windows(2) {
                    let (lo2, hi2) = (pair[0], pair[1]);
                    let from = section
                        .curve
                        .point_at(at_param(lo2, domain, section.closed), tol)?;
                    let to = section
                        .curve
                        .point_at(at_param(hi2, domain, section.closed), tol)?;
                    if from.distance(to) <= tol.confusion() {
                        // A full loop: two arcs, so every strand has two
                        // distinct endpoints.
                        let mid = f64::midpoint(lo2, hi2);
                        pieces.push(SectionPiece {
                            section: si,
                            range: (lo2, mid),
                            hugs,
                            hug_key,
                        });
                        pieces.push(SectionPiece {
                            section: si,
                            range: (mid, hi2),
                            hugs,
                            hug_key,
                        });
                    } else {
                        pieces.push(SectionPiece {
                            section: si,
                            range: (lo2, hi2),
                            hugs,
                            hug_key,
                        });
                    }
                }
            }
            // A stretch running along an edge of each face is those edges'
            // to carry, and they split where it ends, as the section does,
            // and where it passes a vertex of either face, from one edge
            // onto the next.
            // Two faces' edges all but in line (a merged facet column's
            // edge leaning off a pad's corner edge by a hundred-thousandth)
            // are crossed by a section lying all but along both a micron
            // apart along it and a nanometre across, and the piece kept on
            // one side of that stretch ends at one face's crossing: the
            // other face's edge splits under that end too, or its piece of
            // the section ends a nanometre inside it and the face is never
            // split.
            for (lo, hi, hugged) in &bridges {
                for stop in [*lo, *hi] {
                    if !kept_stops
                        .iter()
                        .chain(&vertex_stops)
                        .any(|k| (k - stop).abs() <= tol.parametric())
                    {
                        continue;
                    }
                    let at = section
                        .curve
                        .point_at(at_param(stop, domain, section.closed), tol)?;
                    for (side, own) in [&ga.faces[section.face_a], &gb.faces[section.face_b]]
                        .into_iter()
                        .enumerate()
                    {
                        for &ei in &hugged[side] {
                            let e = &own.edges[ei];
                            if edge_hits.iter().any(|(node, _, on_section, ..)| {
                                *node == e.node && (on_section - stop).abs() <= tol.parametric()
                            }) {
                                continue;
                            }
                            let foot = ogeom_algo::project_on_curve(e.drawn().0, at, 64, tol)?;
                            if foot.distance > reach.max(tol.confusion() * 1e3) {
                                continue;
                            }
                            let on_e =
                                onto_range(e.onto_curve(foot.parameter), &e.curve, e.crange, tol);
                            if on_e <= e.crange.0.min(e.crange.1) + tol.parametric()
                                || on_e >= e.crange.0.max(e.crange.1) - tol.parametric()
                            {
                                continue;
                            }
                            ogeom_core::diag!(
                                Channel::Wire,
                                "PAVE s{si}: side {side} edge {} along ({lo:.6}, {hi:.6}) paved at {on_e:.6} for the stop at {stop:.6}, {:.2e} off",
                                e.node.index(),
                                foot.distance
                            );
                            paves.push((
                                e.node,
                                Pave {
                                    t: on_e,
                                    honesty: honest(section.tolerance, tol).max(foot.distance),
                                    across: None,
                                    spread: 0.0,
                                },
                            ));
                        }
                    }
                }
            }
            Ok((paves, pieces, junctions, reseamed))
        },
    );
    let mut reseams: Vec<(usize, SectionRec)> = Vec::new();
    for (si, work) in paved.into_iter().enumerate() {
        let (found, made, hugged, reseamed) = work?;
        hug_junctions.extend(hugged);
        for (node, at) in found {
            paves.entry(node).or_default().push(at);
        }
        pieces.extend(made);
        if let Some(section) = reseamed {
            reseams.push((si, section));
        }
    }
    // A section whose seam moved names every point by the parameter it had,
    // a period on where it stood before the new start: each piece's range is
    // carried into the new domain, so nothing evaluates it out of domain.
    for (si, section) in reseams {
        let domain = section.curve.domain();
        for piece in pieces.iter_mut().filter(|p| p.section == si) {
            let start = fold(piece.range.0, domain);
            piece.range = (start, start + (piece.range.1 - piece.range.0));
        }
        sections[si] = section;
    }
    share_hugs(&mut pieces, &sections, [ga, gb], tol)?;
    // A piece running along a face's own edge and that edge are one curve,
    // split twice: the piece where its section's cuts fell, the edge where
    // every section that met it paved it, each worked out apart. A loop
    // hugging a rim is split at its middle, the rim wherever another
    // section crossed it. The two faces then walk different pieces of one
    // circle and never sew. Once every pave is in, each such piece is
    // cut wherever its edge is split, and the edge split wherever the piece
    // ends, twice over so a split one piece adds reaches the others.
    for _ in 0..2 {
        let mut settled: Vec<SectionPiece> = Vec::with_capacity(pieces.len());
        for piece in pieces.drain(..) {
            if !(piece.hugs[0] || piece.hugs[1]) {
                settled.push(piece);
                continue;
            }
            let section = &sections[piece.section];
            let domain = section.curve.domain();
            let period = domain.1 - domain.0;
            let (lo, hi) = piece.range;
            let at = |t: f64| {
                section
                    .curve
                    .point_at(at_param(t, domain, section.closed), tol)
            };
            let width = (tol.confusion() * 1e3).max(section.tolerance * 3.0);
            let mid = at(f64::midpoint(lo, hi))?;
            let mut hugged: Vec<&BoundaryEdge> = Vec::new();
            for (side, &hugs) in piece.hugs.iter().enumerate().take(2) {
                if !hugs {
                    continue;
                }
                let face = if side == 0 {
                    &ga.faces[section.face_a]
                } else {
                    &gb.faces[section.face_b]
                };
                for e in &face.edges {
                    let (drawn, drange, off) = e.drawn();
                    let near = width.max(off * 2.0);
                    if e.bound.expanded(near).contains(mid)
                        && distance_to_edge_curve(drawn, drange, mid, tol)? <= near
                    {
                        hugged.push(e);
                    }
                }
            }
            let mut cuts: Vec<f64> = Vec::new();
            for e in &hugged {
                let mut on_edge: Vec<f64> = vec![e.crange.0, e.crange.1];
                if let Some(list) = paves.get(&e.node) {
                    on_edge.extend(list.iter().map(|pave| pave.t));
                }
                let (drawn, _, off) = e.drawn();
                for t in on_edge {
                    let q = drawn.point_at(e.onto_drawn(t), tol)?;
                    let foot = ogeom_algo::project_on_curve(&section.curve, q, 64, tol)?;
                    if foot.distance > width.max(off * 2.0) {
                        continue;
                    }
                    let mut f = foot.parameter;
                    if section.closed && period > 0.0 {
                        while f < lo {
                            f += period;
                        }
                        while f - period >= lo {
                            f -= period;
                        }
                    }
                    if f > lo + tol.parametric() && f < hi - tol.parametric() {
                        cuts.push(f);
                    }
                }
            }
            cuts.push(lo);
            cuts.push(hi);
            cuts.sort_by(f64::total_cmp);
            cuts.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());
            // The edge split wherever the piece now ends.
            for e in &hugged {
                let (drawn, _, off) = e.drawn();
                for &c in &cuts {
                    let q = at(c)?;
                    let foot = ogeom_algo::project_on_curve(drawn, q, 64, tol)?;
                    if foot.distance > width.max(off * 2.0) {
                        continue;
                    }
                    let on_e = onto_range(e.onto_curve(foot.parameter), &e.curve, e.crange, tol);
                    if on_e <= e.crange.0.min(e.crange.1) + tol.parametric()
                        || on_e >= e.crange.0.max(e.crange.1) - tol.parametric()
                    {
                        continue;
                    }
                    let list = paves.entry(e.node).or_default();
                    if list
                        .iter()
                        .all(|pave| (pave.t - on_e).abs() > tol.parametric())
                    {
                        list.push(Pave {
                            t: on_e,
                            honesty: honest(section.tolerance, tol).max(foot.distance),
                            across: None,
                            spread: 0.0,
                        });
                    }
                }
            }
            for pair in cuts.windows(2) {
                settled.push(SectionPiece {
                    section: piece.section,
                    range: (pair[0], pair[1]),
                    hugs: piece.hugs,
                    hug_key: piece.hug_key,
                });
            }
        }
        pieces = settled;
    }
    // Two sections hugging onto one face along one line lay the same split
    // down twice (a top face's trace and its band's, both along the rail
    // they share, which need not be one edge node once a shell has been
    // built), and a line met from both sides at once breaks the
    // arrangement. Pieces of different sections admitted by a hug onto the
    // same face, each lying on the other's curve, are one split, and the
    // section covering the most of it speaks for it.
    {
        let width = tol.confusion() * 1e3;
        let mut drop: Vec<usize> = Vec::new();
        let hugged_onto: Vec<usize> = pieces
            .iter()
            .enumerate()
            .filter(|(_, p)| p.hug_key.is_some())
            .map(|(i, _)| i)
            .collect();
        let mid_of = |p: &SectionPiece| -> OgeomResult<Point> {
            let section = &sections[p.section];
            let domain = section.curve.domain();
            section.curve.point_at(
                at_param(f64::midpoint(p.range.0, p.range.1), domain, section.closed),
                tol,
            )
        };
        let length_of = |p: &SectionPiece| -> OgeomResult<f64> {
            let section = &sections[p.section];
            let domain = section.curve.domain();
            let a = section
                .curve
                .point_at(at_param(p.range.0, domain, section.closed), tol)?;
            let b = section
                .curve
                .point_at(at_param(p.range.1, domain, section.closed), tol)?;
            Ok(a.distance(b))
        };
        // A hug-admitted piece against every other piece kept on the same
        // face: a section admitted the ordinary way, inside both faces,
        // already speaks for the split, and the hug-admitted one yields
        // to it. Two hug-admitted ones yield to the longer.
        for &x in &hugged_onto {
            let px = &pieces[x];
            let Some((_, target_from_a, target_face)) = px.hug_key else {
                continue;
            };
            for (y, py) in pieces.iter().enumerate() {
                if y == x || py.section == px.section || drop.contains(&x) || drop.contains(&y) {
                    continue;
                }
                let sy = &sections[py.section];
                let on_target = if target_from_a {
                    sy.face_a == target_face && !py.hugs[0]
                } else {
                    sy.face_b == target_face && !py.hugs[1]
                };
                if !on_target {
                    continue;
                }
                let sx = &sections[px.section];
                // A closed section's pieces are ranged past its domain end
                // and folded on use. The distance helper samples the range
                // it is given, so it is given the folded one.
                let stretch = |section: &SectionRec, range: (f64, f64)| -> (f64, f64) {
                    folded_range(range, section.curve.domain(), section.closed)
                };
                let on_y =
                    distance_to_edge_curve(&sy.curve, stretch(sy, py.range), mid_of(px)?, tol)?
                        <= width;
                let on_x =
                    distance_to_edge_curve(&sx.curve, stretch(sx, px.range), mid_of(py)?, tol)?
                        <= width;
                if !(on_x && on_y) {
                    continue;
                }
                let loser = if py.hug_key.is_some() && length_of(px)? > length_of(py)? {
                    y
                } else {
                    x
                };
                ogeom_core::diag!(
                    Channel::Wire,
                    "PAVE: s{} and s{} split one face along one line; s{} speaks for it",
                    px.section,
                    py.section,
                    pieces[if loser == x { y } else { x }].section
                );
                drop.push(loser);
                if loser == x {
                    break;
                }
            }
        }
        let mut index = 0_usize;
        pieces.retain(|_| {
            let keep = !drop.contains(&index);
            index += 1;
            keep
        });
    }
    // The audit's verdict. A dropped pair whose surfaces intersect is not a
    // filter bug: planes meet along an infinite line the paving then trims
    // to the faces, usually to nothing. A dropped pair whose section
    // *survives paving* is: some of that curve lies inside both faces, so
    // the faces genuinely meet and the filter's boxes failed to. Contacts
    // need no check: a contact is an owner edge lying in the target face,
    // which forces the boxes to overlap where the edge does.
    //
    // This membership check is the audit's first tooth. The second is the
    // strict replay in general_fuse, which runs the *filtered* fill for the
    // production result and diffs this unfiltered one against it, so the
    // non-compositionality of paving (a section's kept intervals see the
    // other sections' paves) is caught rather than stated as a limit.
    if admit_all {
        for piece in &pieces {
            let section = &sections[piece.section];
            let (fa, fb) = (&ga.faces[section.face_a], &gb.faces[section.face_b]);
            assert!(
                fa.filter.intersects(&fb.filter),
                "bound filter audit: faces {}/{} were dropped by the bound \
                 filter, yet their section paved a surviving piece over \
                 {:?}; the filter under-approximates",
                section.face_a,
                section.face_b,
                piece.range,
            );
        }
    }
    // A contact edge splits where it crosses the target face's boundary,
    // and that crossing is a pave on *both* edges, so the owner's own face
    // splits its boundary consistently and the pieces sew back shared.
    let mut contact_along: Vec<Vec<(f64, f64)>> = vec![Vec::new(); contacts.len()];
    // Where stretches end at the radii the curves hold to, each edge holds
    // to the widest it is measured to on any face of it, and each edge's
    // own record is kept to measure a stretch from its side.
    let mut held: ogeom_core::FastMap<EdgeKey, (f64, &BoundaryEdge)> =
        ogeom_core::FastMap::default();
    if at_stated {
        for e in ga
            .faces
            .iter()
            .chain(&gb.faces)
            .flat_map(|f| f.edges.iter())
        {
            held.entry(e.node)
                .and_modify(|(most, _)| *most = most.max(e.residual))
                .or_insert((e.residual, e));
        }
    }
    // A stretch two edges share, as measured once for the pair: the
    // contact's side of it on `on_a`.
    let mut shared: ogeom_core::FastMap<(EdgeKey, EdgeKey), Vec<ogeom_intersect::Overlap>> =
        ogeom_core::FastMap::default();
    for (ci, contact) in contacts.iter().enumerate() {
        let target = if contact.target_from_a {
            &ga.faces[contact.target_face]
        } else {
            &gb.faces[contact.target_face]
        };
        // A fitted contact sits off exact geometry by its own tolerance, and
        // the crossings it genuinely makes gape by the same: the filter and
        // the finder both widen, or a fitted rail never registers against
        // the exact seam it crosses.
        let reach = tol.confusion().max(contact.tolerance * 2.0);
        let cc = CurveCurveOptions {
            gap: reach.max(CurveCurveOptions::default().gap),
            ..CurveCurveOptions::default()
        };
        for e in &target.edges {
            if !contact.bound.expanded(reach).intersects(&e.bound) {
                continue;
            }
            let found = intersect_curves(&contact.curve, &e.curve, cc, tol)?;
            for crossing in &found.crossings {
                if crossing.gap > reach {
                    continue;
                }
                // The target face's boundary is what its edges cover, not
                // their whole curves. A blend's rail lying on a wedge's end
                // face meets the circle of the end's arc twice: at the
                // arc's end, and a diameter on, outside the face, where a
                // pave leaves the rail a vertex no face needs. A crossing
                // past the edge's ends is on the boundary only within the
                // edge's radius, its ends' and the crossing's own.
                let on_b = onto_range(crossing.on_b, &e.curve, e.crange, tol);
                if (on_b < e.crange.0 - tol.parametric() || on_b > e.crange.1 + tol.parametric())
                    && distance_to_edge_curve(&e.curve, e.crange, crossing.point, tol)?
                        > reach + crossing.reach + e.tolerance.max(e.ends_tolerance)
                {
                    continue;
                }
                // A contact ending on a target edge needs no split of its
                // own there, but the edge does: the contact's strand ends
                // on it, and an edge left whole leaves the strand dangling.
                let at_end = crossing.on_a < contact.crange.0 + tol.parametric()
                    || crossing.on_a > contact.crange.1 - tol.parametric();
                if at_end
                    && (crossing.on_a < contact.crange.0 - tol.parametric()
                        || crossing.on_a > contact.crange.1 + tol.parametric())
                {
                    continue;
                }
                // Touching is not crossing, at curve level as at surface
                // level: where the two curves run nearly parallel (a blend
                // arc meeting the edge it is tangent to), the "crossing" is
                // numerical noise smeared along the contact, and paving it
                // would plant a vertex a hair off both curves. A crossing at
                // a shallow turn is no touch: the two part in proportion
                // either side of it, and each side of it is a strand of its
                // own. The two cross where they are found to, as at any
                // turn, not anywhere along a touching run.
                let mut reach = crossing.reach;
                let mut onto_end = None;
                if tangential(&contact.curve, crossing.on_a, &e.curve, crossing.on_b, tol)? {
                    let shallow = at_stated
                        && crosses_shallowly(
                            &contact.curve,
                            contact.crange,
                            contact.tolerance,
                            e,
                            (crossing.on_a, crossing.on_b),
                            tol,
                        )?;
                    if !shallow {
                        // A touch at the contact's end is the end's vertex
                        // lying on the edge: the contact's strand ends
                        // there, and the edge splits at the vertex's foot,
                        // which is known where the touch is not. A rounded
                        // box's arc between a corner ball and a round ends
                        // on a drill's seam that runs along the round's
                        // tangent line, tangent to it.
                        if at_end {
                            let end = if crossing.on_a
                                < f64::midpoint(contact.crange.0, contact.crange.1)
                            {
                                contact.crange.0
                            } else {
                                contact.crange.1
                            };
                            let vertex = contact.curve.point_at(end, tol)?;
                            let (gap, on_edge) =
                                nearest_on_edge_curve(&e.curve, e.crange, vertex, tol)?;
                            let width = held_radius(&held, contact.node, contact.tolerance)
                                + held_radius(&held, e.node, e.tolerance);
                            if gap <= width
                                && on_edge > e.crange.0 + tol.parametric()
                                && on_edge < e.crange.1 - tol.parametric()
                            {
                                paves.entry(e.node).or_default().push(Pave {
                                    t: on_edge,
                                    honesty: honest(contact.tolerance, tol).max(gap),
                                    across: None,
                                    spread: 0.0,
                                });
                            }
                        }
                        continue;
                    }
                    reach = 0.0;
                    let width = held_radius(&held, contact.node, contact.tolerance)
                        + held_radius(&held, e.node, e.tolerance);
                    onto_end =
                        shallow_at_end(contact, e, (crossing.on_a, crossing.on_b), width, tol)?;
                }
                let honesty = honest(contact.tolerance, tol).max(reach);
                match onto_end {
                    Some(End::OfEdge { on_contact, gap }) => {
                        if on_contact > contact.crange.0 + tol.parametric()
                            && on_contact < contact.crange.1 - tol.parametric()
                        {
                            paves.entry(contact.node).or_default().push(Pave {
                                t: on_contact,
                                honesty: honesty.max(gap),
                                across: None,
                                spread: 0.0,
                            });
                        }
                        continue;
                    }
                    Some(End::OfContact { on_edge, gap }) => {
                        if on_edge > e.crange.0 + tol.parametric()
                            && on_edge < e.crange.1 - tol.parametric()
                        {
                            paves.entry(e.node).or_default().push(Pave {
                                t: on_edge,
                                honesty: honesty.max(gap),
                                across: None,
                                spread: 0.0,
                            });
                        }
                        continue;
                    }
                    None => {}
                }
                if !at_end {
                    paves.entry(contact.node).or_default().push(Pave {
                        t: crossing.on_a,
                        honesty,
                        across: None,
                        spread: 0.0,
                    });
                }
                if on_b > e.crange.0 + tol.parametric() && on_b < e.crange.1 - tol.parametric() {
                    paves.entry(e.node).or_default().push(Pave {
                        t: on_b,
                        honesty,
                        across: None,
                        spread: 0.0,
                    });
                }
            }
            // A span of the contact running along a target boundary edge
            // splits nothing: it is already boundary on both sides
            // (identically stacked boxes are all such spans), and duplicating
            // it as a strand would cancel the boundary it copies.
            // The intersector answers coincidence in closed form for the
            // analytic pairs only. A fitted curve lying on an exact one
            // (a marched sphere's rim on a wedge cap's arc) comes back as
            // crossings or nothing. That span is measured instead: the
            // stretch of the contact within the pair's honesty of the edge.
            // An analytic overlap survives only where it lands inside the
            // contact's own window. A trimmed circle written a turn up
            // clips to nothing, and is then measured like a fitted one.
            let survives = |overlap: &ogeom_intersect::Overlap| -> bool {
                let (lo, hi) = if overlap.on_a.0 <= overlap.on_a.1 {
                    overlap.on_a
                } else {
                    (overlap.on_a.1, overlap.on_a.0)
                };
                let carried = if contact.curve.is_periodic() {
                    let (dlo, dhi) = contact.curve.domain();
                    let period = dhi - dlo;
                    let turn = ((contact.crange.0 - lo) / period).floor();
                    [turn, turn + 1.0, turn - 1.0]
                        .into_iter()
                        .map(|k| (k.mul_add(period, lo), k.mul_add(period, hi)))
                        .map(|(a, b)| b.min(contact.crange.1) - a.max(contact.crange.0))
                        .fold(f64::NEG_INFINITY, f64::max)
                } else {
                    hi.min(contact.crange.1) - lo.max(contact.crange.0)
                };
                carried > tol.parametric()
            };
            // The closed-form overlap is between the two *curves*. The
            // stretch that is boundary is what the *edge* covers of it. A
            // box cut from an L-bracket flush with the bracket's wall puts
            // the box's wall-side edge on the line of the end face's own
            // edge along the wall, a length below it: read as along that
            // edge over the whole curve, the strip's side is never paved
            // and the end face keeps the strip. Clipped through the
            // overlap's own correspondence so the carry below stays affine.
            let clipped: Vec<ogeom_intersect::Overlap> = found
                .overlaps
                .iter()
                .flat_map(|overlap| {
                    let span = overlap.on_a.1 - overlap.on_a.0;
                    let parts = if span.abs() <= f64::MIN_POSITIVE {
                        Vec::new()
                    } else {
                        overlap_within(overlap, e.crange, &contact.curve, tol)
                            .into_iter()
                            .map(|(part, _)| part)
                            .collect()
                    };
                    let to_b = move |t: f64| {
                        overlap.on_b.0
                            + (overlap.on_b.1 - overlap.on_b.0) * (t - overlap.on_a.0) / span
                    };
                    parts
                        .into_iter()
                        .map(move |(lo, hi)| ogeom_intersect::Overlap {
                            on_a: (lo, hi),
                            on_b: (to_b(lo), to_b(hi)),
                        })
                })
                .collect();
            // A line and a curved conic share no stretch, however loosely
            // either is held: a straight edge tangent to an arc stays within
            // the measuring width of it for a short run either side of the
            // touch, and read as a span there it splits the arc at a vertex
            // the face across the arc never has.
            let conic = |c: &Curve| {
                matches!(
                    c,
                    Curve::Circle(_) | Curve::Ellipse(_) | Curve::Hyperbola(_) | Curve::Parabola(_)
                )
            };
            let straight_on_conic = (matches!(contact.curve, Curve::Line(_)) && conic(&e.curve))
                || (conic(&contact.curve) && matches!(*e.curve, Curve::Line(_)));
            // Where the two are held to their own radii, a fitted pair is
            // measured, as the intersector's coincidence for it ends at its
            // own width and from one side only.
            let fitted = !(exact(&contact.curve) && exact(&e.curve));
            let clipped = if at_stated && fitted {
                Vec::new()
            } else {
                clipped
            };
            let measured: Vec<ogeom_intersect::Overlap> =
                if clipped.iter().any(survives) || straight_on_conic {
                    Vec::new()
                } else if at_stated {
                    shared_overlaps(contact, e, &held, &mut shared, tol)?
                } else {
                    measured_overlaps(
                        &contact.curve,
                        contact.crange,
                        contact.tolerance,
                        e,
                        None,
                        tol,
                    )?
                };
            for overlap in clipped.iter().chain(measured.iter()) {
                let ordered = |r: (f64, f64)| if r.0 <= r.1 { r } else { (r.1, r.0) };
                let (lo, hi) = ordered(overlap.on_a);
                // The overlap is between the two *curves*. What interferes is
                // the stretch both *edges* actually cover. A hole's arc and
                // the disc that fills it lie on one circle, so the curves
                // overlap over the whole turn while the arc covers three
                // quarters of it, and paving at the turn's ends says nothing,
                // where paving at the arc's ends is exactly the split the
                // other side needs to sew against.
                // A periodic curve's overlap is answered on its base turn;
                // an edge written a turn up (a wedge cap's arc at 2π..2.5π)
                // covers it only once the turn is carried across. A stretch
                // and a window that both run near a full turn can meet in
                // two turns, and each meeting is a stretch along the edge.
                let turns: Vec<f64> = if contact.curve.is_periodic() {
                    let (dlo, dhi) = contact.curve.domain();
                    let period = dhi - dlo;
                    let turn = ((contact.crange.0 - lo) / period).floor();
                    [turn - 1.0, turn, turn + 1.0]
                        .into_iter()
                        .map(|k| k * period)
                        .collect()
                } else {
                    vec![0.0]
                };
                for shift in turns {
                    let (lo, hi) = (
                        (lo + shift).max(contact.crange.0),
                        (hi + shift).min(contact.crange.1),
                    );
                    if hi - lo <= tol.parametric() {
                        continue;
                    }
                    let span = overlap.on_a.1 - overlap.on_a.0;
                    let carry = |t: f64| -> f64 {
                        if span.abs() <= f64::MIN_POSITIVE {
                            overlap.on_b.0
                        } else {
                            overlap.on_b.0
                                + (overlap.on_b.1 - overlap.on_b.0) * (t - overlap.on_a.0) / span
                        }
                    };
                    // A stretch ended where the two curves part by the radii
                    // they hold to has one end on each curve, that far apart
                    // or less: one vertex holds both, and the split on
                    // either curve owns the distance there.
                    let honesty = |t: f64| -> OgeomResult<f64> {
                        let own = honest(contact.tolerance, tol);
                        if !at_stated {
                            return Ok(own);
                        }
                        let gap = contact
                            .curve
                            .point_at(t, tol)?
                            .distance(e.curve.point_at(carry(t), tol)?);
                        Ok(own.max(gap))
                    };
                    for t in [lo, hi] {
                        paves.entry(contact.node).or_default().push(Pave {
                            t,
                            honesty: honesty(t)?,
                            across: None,
                            spread: 0.0,
                        });
                    }
                    ogeom_core::diag!(
                        Channel::Wire,
                        "CONTACT c{ci} along edge {} over ({lo:.6}, {hi:.6}) of {:?}: {:?} .. {:?}",
                        e.node.index(),
                        contact.crange,
                        contact.curve.point_at(lo, tol).ok(),
                        contact.curve.point_at(hi, tol).ok()
                    );
                    contact_along[ci].push((lo, hi));
                    let mut apart = 0.0_f64;
                    for t in [lo, f64::midpoint(lo, hi), hi] {
                        apart = apart.max(
                            contact
                                .curve
                                .point_at(t, tol)?
                                .distance(e.curve.point_at(carry(t), tol)?),
                        );
                    }
                    if apart > contact.tolerance.max(e.tolerance) {
                        ALONG_APART.set(ALONG_APART.get().max(apart));
                    }
                    // The *target* edge splits where the shared stretch ends,
                    // exactly as the contact does. Without this, the face across
                    // the overlap keeps one long boundary edge where its new
                    // neighbours carry two short ones, and sew (which matches
                    // edges whole) can pair it with neither. The clamped ends are
                    // carried across by the correspondence the overlap itself
                    // states, which is affine over the shared stretch.
                    let (mut tlo, mut thi) = ((carry(lo), honesty(lo)?), (carry(hi), honesty(hi)?));
                    if tlo.0 > thi.0 {
                        core::mem::swap(&mut tlo, &mut thi);
                    }
                    let target_domain = e.curve.domain();
                    let periodic = e.curve.is_periodic();
                    for (t, honesty) in [tlo, thi] {
                        // A correspondence across a full turn can run the
                        // parameter past the domain. The pave belongs where the
                        // edge actually is.
                        let t = if periodic { fold(t, target_domain) } else { t };
                        if t > e.crange.0 + tol.parametric() && t < e.crange.1 - tol.parametric() {
                            paves.entry(e.node).or_default().push(Pave {
                                t,
                                honesty,
                                across: None,
                                spread: 0.0,
                            });
                        }
                    }
                }
            }
        }
    }

    let mut same_a: Vec<Vec<usize>> = vec![Vec::new(); ga.faces.len()];
    let mut same_b: Vec<Vec<usize>> = vec![Vec::new(); gb.faces.len()];
    for (ia, ib) in same_pairs {
        same_a[ia].push(ib);
        same_b[ib].push(ia);
    }
    Ok((
        sections,
        pieces,
        contacts,
        tangents,
        contact_along,
        paves,
        same_a,
        same_b,
        hug_junctions,
    ))
}

/// An exact section cut at the chart degeneracies it runs through, each
/// piece carrying exact pcurves on both faces.
///
/// The degeneracies are not guessed from the surfaces: they are the faces'
/// own *pole edges*, the degenerate edges the topology already carries, so a
/// face whose chart has a pole says so and one that has none costs nothing.
///
/// `None` means the split is no use here (no degeneracy on the curve, or a
/// piece whose projection still has no closed form), and the caller falls
/// back to marching, which is the honest answer rather than a fitted pcurve
/// pretending to be exact.
pub(crate) fn split_at_degeneracies(
    curve: &Curve,
    fa: &GFace,
    fb: &GFace,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<(Curve, PlanarCurve, PlanarCurve)>>> {
    let mut stops: Vec<f64> = Vec::new();
    let domain = curve.domain();
    // The split points are the *surfaces'* chart degeneracies, not merely the
    // faces' pole edges: a meridian section runs through both of a sphere's
    // poles, and a face that owns only the north one still cannot chart an
    // arc that wraps through the south. The surface knows where its chart
    // collapses whether or not the face's trim reaches there.
    let mut candidates: Vec<Point> = Vec::new();
    for face in [fa, fb] {
        for pole in &face.poles {
            candidates.push(pole.point);
        }
        match &face.surface {
            SurfaceGeometry::Sphere(s) => {
                let sphere = s.sphere();
                let axis = sphere.frame().z().vector();
                candidates.push(sphere.centre() + axis * sphere.radius());
                candidates.push(sphere.centre() - axis * sphere.radius());
            }
            SurfaceGeometry::Cone(c) => {
                candidates.push(c.cone().apex());
            }
            _ => {}
        }
    }
    for point in candidates {
        let found = ogeom_algo::project_on_curve(curve, point, 256, tol)?;
        if found.distance <= tol.confusion() {
            stops.push(found.parameter);
        }
    }
    if stops.is_empty() {
        return Ok(None);
    }
    stops.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    stops.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());

    // A closed curve is cut into the arcs between successive stops, the last
    // wrapping past the domain end. An open one keeps its own ends as stops.
    let closed = curve.is_closed(tol) || curve.is_periodic();
    let mut arcs: Vec<(f64, f64)> = Vec::new();
    if closed {
        if stops.len() < 2 {
            // One stop on a closed curve leaves one arc, from the stop right
            // round to itself, which still has the pole at both ends.
            let period = domain.1 - domain.0;
            arcs.push((stops[0], stops[0] + period));
        } else {
            let period = domain.1 - domain.0;
            for i in 0..stops.len() {
                let lo = stops[i];
                let hi = if i + 1 < stops.len() {
                    stops[i + 1]
                } else {
                    stops[0] + period
                };
                arcs.push((lo, hi));
            }
        }
    } else {
        let mut cuts = vec![domain.0];
        cuts.extend(
            stops
                .iter()
                .copied()
                .filter(|t| *t > domain.0 + tol.parametric() && *t < domain.1 - tol.parametric()),
        );
        cuts.push(domain.1);
        for pair in cuts.windows(2) {
            arcs.push((pair[0], pair[1]));
        }
    }

    let mut out = Vec::with_capacity(arcs.len());
    for (lo, hi) in arcs {
        if hi - lo <= tol.parametric() {
            continue;
        }
        let (Some(pa), Some(pb)) = (
            ogeom_intersect::exact_pcurve_over(curve, (lo, hi), &fa.surface, tol),
            ogeom_intersect::exact_pcurve_over(curve, (lo, hi), &fb.surface, tol),
        ) else {
            return Ok(None);
        };
        let Ok(piece) = ogeom_geom::TrimmedCurve::new(curve.clone(), lo, hi, tol) else {
            return Ok(None);
        };
        out.push((piece.into(), pa, pb));
    }
    if out.is_empty() {
        return Ok(None);
    }
    Ok(Some(out))
}

/// Whether `curve` starts on a pole of either chart.
pub(crate) fn on_a_chart_pole(
    curve: &Curve,
    surfaces: [&SurfaceGeometry; 2],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let start = curve.point_at(curve.domain().0, tol)?;
    Ok(surfaces.iter().any(|surface| {
        ogeom_intersect::march::chart_poles(surface, tol)
            .iter()
            .any(|&(_, _, pole)| pole.distance(start) <= tol.confusion())
    }))
}

/// Each face's carrier, numbered within the solid, with the window of the
/// carrier's parameters its faces' trims reach.
///
/// The window is the union of the parameter boxes of the faces on that
/// surface, read from their sampled trims and poles, grown a quarter of
/// its span on every side so a section crossing a trim is traced past it
/// and cut there by the trim. A face whose trim does not sample leaves its
/// carrier whole.
pub(crate) fn carrier_windows(
    g: &GSolid,
    tol: Tolerances,
) -> Vec<(usize, Option<ogeom_intersect::ParameterWindow>)> {
    let mut numbered: ogeom_core::FastMap<&(ogeom_topo::SurfaceId, Location), usize> =
        ogeom_core::FastMap::default();
    // Per carrier the box reached so far, low corner then high, or `None`
    // once a face on it fails to sample.
    let empty = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    let mut reach: Vec<Option<(Point2, Point2)>> = Vec::new();
    let mut of_face = Vec::with_capacity(g.faces.len());
    for face in &g.faces {
        let next = numbered.len();
        let at = *numbered.entry(&face.carrier).or_insert(next);
        if at == reach.len() {
            reach.push(Some(empty));
        }
        of_face.push(at);
        let Some((low, high)) = reach[at].as_mut() else {
            continue;
        };
        let mut points: Vec<Point2> = face.trim_lines(tol).iter().flatten().copied().collect();
        for pole in &face.poles {
            for k in 0..=4 {
                let t = pole.prange.0 + (pole.prange.1 - pole.prange.0) * f64::from(k) / 4.0;
                if let Ok(p) = pole.pcurve.point_at(t, tol) {
                    points.push(p);
                }
            }
        }
        if points.is_empty() || points.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
            reach[at] = None;
            continue;
        }
        for p in points {
            *low = Point2::new(low.x.min(p.x), low.y.min(p.y));
            *high = Point2::new(high.x.max(p.x), high.y.max(p.y));
        }
    }
    let grow = |lo: f64, hi: f64| {
        let margin = (hi - lo) * 0.25;
        (lo - margin, hi + margin)
    };
    let windows: Vec<Option<ogeom_intersect::ParameterWindow>> = reach
        .into_iter()
        .map(|bound| bound.map(|(low, high)| (grow(low.x, high.x), grow(low.y, high.y))))
        .collect();
    of_face.into_iter().map(|at| (at, windows[at])).collect()
}

/// A surface narrowed to the reach of a bound, for the marcher's benefit.
///
/// Seeding samples a surface's *parameter box*, so a plane stored over
/// ±10^9 (which is how an unbounded carrier reaches this code) is sampled
/// at a spacing a hundred million times anything it could be meeting, and
/// no seed ever lands near the curve. Narrowing to where the two faces
/// actually are is what lets the seeding see it. The chart is untouched:
/// only the window the marcher walks changes, so pcurves fitted in it mean
/// the same thing on the surface as stored.
pub(crate) fn windowed_to(surface: &SurfaceGeometry, bound: &ogeom_math::Aabb) -> SurfaceGeometry {
    let corners = bound.corners();
    if corners.is_empty() {
        return surface.clone();
    }
    match surface {
        SurfaceGeometry::Plane(p) => {
            let frame = p.plane().frame();
            let (mut u0, mut u1, mut v0, mut v1) = (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            );
            for c in &corners {
                let local = frame.to_local(*c);
                u0 = u0.min(local.x);
                u1 = u1.max(local.x);
                v0 = v0.min(local.y);
                v1 = v1.max(local.y);
            }
            let margin = ((u1 - u0) + (v1 - v0)).mul_add(0.25, 1.0);
            let (want_u, want_v) = ((u0 - margin, u1 + margin), (v0 - margin, v1 + margin));
            let (have_u, have_v) = ogeom_geom::Surface::domain(p);
            if have_u.1 - have_u.0 <= want_u.1 - want_u.0
                && have_v.1 - have_v.0 <= want_v.1 - want_v.0
            {
                return surface.clone();
            }
            ogeom_geom::PlaneSurface::over(p.plane(), want_u, want_v)
                .map_or_else(|_| surface.clone(), Into::into)
        }
        SurfaceGeometry::Cylinder(c) => {
            let frame = c.cylinder().frame();
            let (mut h0, mut h1) = (f64::INFINITY, f64::NEG_INFINITY);
            for corner in &corners {
                let h = (*corner - frame.origin()).dot(frame.z().vector());
                h0 = h0.min(h);
                h1 = h1.max(h);
            }
            let margin = (h1 - h0).mul_add(0.25, 1.0);
            let want = (h0 - margin, h1 + margin);
            let have = ogeom_geom::Surface::domain(c).1;
            if have.1 - have.0 <= want.1 - want.0 {
                return surface.clone();
            }
            ogeom_geom::CylinderSurface::new(c.cylinder(), want)
                .map_or_else(|_| surface.clone(), Into::into)
        }
        other => other.clone(),
    }
}

/// A cylinder whose marched sections with a patch stop short of an end of
/// its window, the window stretched along its axis to carry them past it.
///
/// A primitive's wall is stored over exactly its own height, so a section
/// crossing the rim is walked to the window's edge and stops there, a step
/// short of the rim, and the rim is never crossed. Stretched a fiftieth of
/// its height at that end, the walk carries on past the rim, crosses it like
/// any other boundary edge, and the face's trim cuts it there. The chart is
/// untouched, only the window grows. `None` where every section ends on the
/// window's ends or clear of them, for any other surface, and against an
/// analytic surface, whose sections are closed forms already exact to the
/// window.
pub(crate) fn reaching(
    surface: &SurfaceGeometry,
    other: &SurfaceGeometry,
    met: &ogeom_intersect::SurfaceIntersection,
    tol: Tolerances,
) -> OgeomResult<Option<SurfaceGeometry>> {
    let SurfaceGeometry::Cylinder(c) = surface else {
        return Ok(None);
    };
    if !matches!(other, SurfaceGeometry::BSpline(_)) {
        return Ok(None);
    }
    let ogeom_intersect::SurfaceIntersection::Along(curves) = met else {
        return Ok(None);
    };
    let frame = c.cylinder().frame();
    let have = ogeom_geom::Surface::domain(c).1;
    let past = (have.1 - have.0) * 0.02;
    let short = |gap: f64| gap > tol.confusion() && gap <= past;
    let (mut low, mut high) = (false, false);
    for section in curves {
        let (t0, t1) = section.curve.domain();
        for t in [t0, t1] {
            let at = section.curve.point_at(t, tol)?;
            let h = (at - frame.origin()).dot(frame.z().vector());
            low |= short(h - have.0);
            high |= short(have.1 - h);
        }
    }
    if !low && !high {
        return Ok(None);
    }
    let want = (
        if low { have.0 - past } else { have.0 },
        if high { have.1 + past } else { have.1 },
    );
    Ok(ogeom_geom::CylinderSurface::new(c.cylinder(), want)
        .ok()
        .map(Into::into))
}

/// March a pair whose exact section has no closed-form pcurve: the fitted
/// branches, and the reach they were cut at poles within, which their
/// tolerance must cover.
pub(crate) fn march_pair(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: &ogeom_intersect::IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<(Vec<ogeom_intersect::IntersectionCurve>, f64)> {
    use ogeom_intersect::{approximate_branch, branches_cut_at_poles, passes_a_pole};
    // Fitted whole past a pole of either chart within its stated
    // tolerance, a section passes through the pole as far as the fit can
    // say: the pair is traced again, cut there as within the chord.
    let mut pole_reach = 0.0;
    let out = loop {
        let traced = branches_cut_at_poles(a, b, options.marching, pole_reach, tol)?;
        let mut out = Vec::new();
        let mut again: Option<f64> = None;
        for branch in &traced {
            let fitted = approximate_branch(a, b, branch, options.tolerance, tol)?;
            let stated = options.marching.chord + fitted.fit_error;
            if pole_reach == 0.0 && passes_a_pole(&fitted.curve.clone().into(), [a, b], stated, tol)
            {
                again = Some(again.map_or(stated, |r: f64| r.max(stated)));
            }
            out.push(fitted);
        }
        match again {
            Some(reach) => pole_reach = reach,
            None => break out,
        }
    };
    if out.is_empty() {
        // No branch found is two different stories. Surfaces that measurably
        // stand apart over their stated extents (a blend's cylinder and a
        // far-off bevel plane whose ellipse of intersection lies beyond
        // both) simply do not interact, and an empty section is the true
        // answer. Only a pair that comes close and still resolves nothing
        // is beyond the intersector.
        if surfaces_stand_apart(a, b, tol) {
            return Ok((out, pole_reach));
        }
        ogeom_bail!(
            NotDone,
            "an exact section has no closed-form pcurve and marching resolved \
             no branch; the configuration is beyond the intersector's current \
             reach"
        );
    }
    Ok((out, pole_reach))
}

/// Whether two surfaces measurably keep their distance over their stated
/// extents.
///
/// A conservative grid measurement: sample the smaller-extent surface and
/// project each sample onto the other. Apart means every sample clears a
/// margin scaled to the extents. Surfaces with unbounded or enormous stated
/// domains (an imported plane's billion units) are never called apart this
/// way, because a grid over them samples nothing.
pub(crate) fn surfaces_stand_apart(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    use ogeom_geom::Surface as _;
    let extent_of = |s: &SurfaceGeometry| -> f64 {
        let ((ua, ub), (va, vb)) = s.domain();
        (ub - ua).abs().max((vb - va).abs())
    };
    let (sample, against) = if extent_of(a) <= extent_of(b) {
        (a, b)
    } else {
        (b, a)
    };
    let span = extent_of(sample);
    if !span.is_finite() || span > 1e4 {
        return false;
    }
    let ((ua, ub), (va, vb)) = sample.domain();
    const GRID: usize = 9;
    // One seeding grid over the far surface, asked a hundred times: the
    // same seeds and the same Newton a per-call projection uses, at
    // hundreds of evaluations instead of tens of thousands for every
    // genuine near-miss.
    let Ok(seeds) = ogeom_algo::SurfaceSeeds::over(against, 16, tol) else {
        return false;
    };
    let mut clearance = f64::INFINITY;
    for i in 0..=GRID {
        for j in 0..=GRID {
            #[allow(clippy::cast_precision_loss)]
            let u = ua + (ub - ua) * i as f64 / GRID as f64;
            #[allow(clippy::cast_precision_loss)]
            let v = va + (vb - va) * j as f64 / GRID as f64;
            let Ok(p) = sample.point_at(u, v, tol) else {
                return false;
            };
            let Ok(projection) = seeds.project(against, p, tol) else {
                return false;
            };
            clearance = clearance.min(projection.distance);
            if clearance <= tol.confusion() * 1e3 {
                return false;
            }
        }
    }
    clearance > tol.confusion() * 1e3
}

/// Whether two faces lie on one surface where they overlap, to within what
/// they state they may stand off their own: a face fitted to a mesh, and a
/// face swept from a sketch of its outline whose points a float's rounding
/// moved a few microns, are one wall by the fit's own measure, though the
/// arc through the rounded points has another centre and radius and the
/// closed form, held to the confusion distance, sees two.
///
/// Asked only where the faces, or a fuzzy boolean's widened confusion,
/// allow more than the default confusion distance's few, or with
/// `as_stated` at any tolerance (see [`AS_STATED`]). The smaller face
/// is sampled on a grid over its chart, inside its
/// trim. Every sample inside the other face's trim must lie within the
/// stated tolerance of the other's surface, and a few must, for the answer
/// to mean anything. A pair that crosses fails at its first sample off the
/// other.
pub(crate) fn coincide_as_stated(fa: &GFace, fb: &GFace, as_stated: bool, tol: Tolerances) -> bool {
    const GRID: usize = 8;
    const EVIDENCE: usize = 4;
    // What the faces state, or the run's own confusion where a fuzzy
    // boolean has widened it past theirs: a fuzz is a statement of how far
    // the two solids' geometry may disagree and still be one.
    let reach = (fa.tolerance + fb.tolerance).max(tol.confusion() * 2.0);
    if reach <= ogeom_core::tolerance::CONFUSION * 10.0 && !as_stated {
        return false;
    }
    // Two closed forms of different kinds share no patch: a plane is no
    // cylinder anywhere, however close they pass.
    let kind = |s: &SurfaceGeometry| match s {
        SurfaceGeometry::Plane(_) => Some(0),
        SurfaceGeometry::Cylinder(_) => Some(1),
        SurfaceGeometry::Cone(_) => Some(2),
        SurfaceGeometry::Sphere(_) => Some(3),
        SurfaceGeometry::Torus(_) => Some(4),
        _ => None,
    };
    if let (Some(a), Some(b)) = (kind(&fa.surface), kind(&fb.surface))
        && a != b
    {
        return false;
    }
    let (sampled, against) = if fa.bound.diagonal() <= fb.bound.diagonal() {
        (fa, fb)
    } else {
        (fb, fa)
    };
    let Some(lines) = sampled.outline(tol) else {
        return false;
    };
    let borrowed: Vec<&[Point2]> = lines.iter().map(Vec::as_slice).collect();
    let (lo, hi) = lines.iter().flatten().fold(
        (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        |(lo, hi), p| {
            (
                Point2::new(lo.x.min(p.x), lo.y.min(p.y)),
                Point2::new(hi.x.max(p.x), hi.y.max(p.y)),
            )
        },
    );
    if !(lo.x.is_finite() && hi.x.is_finite() && lo.y.is_finite() && hi.y.is_finite()) {
        return false;
    }
    // Wide enough to find the other face's trim about a sample standing
    // off its surface, so a sample off it is caught rather than skipped.
    let search = reach.max(sampled.bound.diagonal().max(against.bound.diagonal()) * 1e-3);
    let mut evidence = 0_usize;
    for i in 0..GRID {
        for j in 0..GRID {
            #[allow(clippy::cast_precision_loss)]
            let at = Point2::new(
                lo.x + (hi.x - lo.x) * ((i as f64 + 0.5) / GRID as f64),
                lo.y + (hi.y - lo.y) * ((j as f64 + 0.5) / GRID as f64),
            );
            if !arrange::inside_many(&borrowed, at) {
                continue;
            }
            let Ok(p) = sampled.surface.point_at(at.x, at.y, tol) else {
                return false;
            };
            if !against.bound.expanded(search).contains(p) {
                continue;
            }
            let Some(uv) = chart_point_within(against, p, search, tol) else {
                continue;
            };
            let Ok(foot) = against.surface.point_at(uv.x, uv.y, tol) else {
                return false;
            };
            if foot.distance(p) > reach {
                return false;
            }
            evidence += 1;
        }
    }
    evidence >= EVIDENCE
}

/// How far off a partner's surface a probe of `face` may stand and still
/// be read on it: what the two faces state, or the confusion distance's
/// few.
pub(crate) fn partner_reach(face: &GFace, partner: &GFace, tol: Tolerances) -> f64 {
    (face.tolerance + partner.tolerance).max(tol.confusion() * 10.0)
}
