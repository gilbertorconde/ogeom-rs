# Plan

What is left to build, why each piece is not built yet, and how each one gets
done. This replaces the scope map, which had become a record of finished work
with the unfinished part scattered through it.

What is built: topology and the geometry vocabulary, the intersector, the
boolean, healing, blending, offsets and sweeps, tessellation, drawings, the
document layer, and STEP, STL, DXF, glTF, OBJ, PLY, VRML and 3MF exchange.

What follows is the remainder, and it is now *derived rather than
remembered*: `docs/PARITY.md` audits every public header of the reference's
modelling modules against 97 named capabilities, the gate in
`tools/check.sh` holds the audit to its evidence, and every `absent` and
`partial` row there anchors back into a section here. This document keeps
the *how* — what each gap needs and why it is not built yet; the ledger
keeps the *what*. An item missing from both is a gate failure, not an
oversight.

## How this project works

These are the rules the work is held to. They are not aspirations; every one
of them has rejected a patch.

- **Measured, not asserted.** A claim about geometry is made by a test that
  measures it against a closed form, a published value, or an independent
  computation. "It looks right" is not a result.
- **An approximation that says so is allowed; one that does not is not.**
  Sampling, fitting and stated tolerances are fine, *stated*. A number that
  hides its own error is not.
- **Refusals are by name, and pinned.** Where the kernel cannot do
  something, it says which thing and why, in the error text, and a test
  holds it to that. A silent wrong answer is the only unacceptable outcome.
- **No stubs.** Either it is implemented, or it is in this document.
- **Scope is parity with the reference kernel's modelling modules, and nothing
  else.** `docs/SCOPE.md` is normative and says how a case is decided. Parity is
  a claim about capability rather than structure; usage data sequences the work
  and never sets its bounds.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against,
  bundles or imports another kernel, and the design is arrived at here rather
  than mirrored from one. The field's vocabulary is used throughout, because
  that is how the field talks about itself.

## The remaining work

### A. The boolean's interference table — **A1–A5 closed**

Five otherwise unrelated gaps, held to be one family. Each is pinned by a
measured test in `crates/ogeom/tests/boolean_interference.rs`, written before
the change and asserting volumes against closed forms, not merely that no
error came back. All five pass.

**The gaps, and what each turned out to be.**

*A1, a tool flush with the part* — refill a bore with the cylinder that cut
it, ends flush with the faces it broke. Two things: an overlap between two
*curves* was being read as an overlap between the two *edges*, so a hole's
arc and the disc filling it paved each other at the turn's ends rather than
the arc's and the two sides sewed against different subdivisions; and the
rule for which copy of a coincident pair to keep was stated by argument
identity rather than by region, so the tool's cap — filling a hole the part
has no face for at all — was dropped and the result had a hole in it. The
substitution is now by containment: a piece is dropped only where a piece of
the other argument genuinely stands where it does.

*A2, a section through a chart pole* — a plane cutting a ball on its own
axis. Two things again. The section's ends land exactly on the poles, and an
*open* section's own end is its domain end, which the rebuild was folding
round to the domain start — the far pole, eight radii away. And the section
itself was being marched, because a meridian's chart image has no closed
form: its longitude jumps half a turn at each pole. Each *half* does have
one, and it is a straight line — writing the sphere's axis as
`Z = cos α·X + sin α·Y` in the circle's own frame, the latitude is
`asin(cos(t − α))`, which on `t − α ∈ [0, π]` is exactly `π/2 − (t − α)`,
affine in the circle's own parameter. So the section is cut at the poles the
face's own degenerate edges name, and each half comes back exact rather than
fitted. The plane's half ball is now measured against `2πr³/3`, not against
a marched approximation of it.

*A3, a ball on a box's corner vertex* — three placements: tangent at the
vertex from outside, centred on it (where the closed form is an octant), and
with the vertex exactly on the sphere so that all three of the box's faces
cut the ball through one point. The last is a triple point, and what it
needed was the seam: a crossing found on a periodic curve comes back in
`[0, 2π)` while the seam edge covers `[-π/2, π/2]`, so every crossing on it
was outside its range, discarded, and the seam never split where the chain
of section arcs actually met it.

*A4, a section through tangential contact* — a plane through a bore's axis
meets the wall along its rulings, and one of those rulings *is* the wall's
own seam. The exclusion for "this section is already boundary" was global,
so the ruling vanished from the plane too — where it is not boundary at all
but the curve separating the two halves the section leaves. It is now
recorded per face, and a section is dropped only from the face that already
carries it, either as a boundary edge or as a contact strand. This buys the
drawing feature D2 as well.

*A5, a shell around a three-cylinder tip* — the corner block less the ball,
which is the tool the three-blend corner needs. The block's three far planes
are *tangent* to the ball, at exactly the three vertices adjacent to the
corner; the near three cut it in a meridian, a meridian along the seam, and
the equator. What it needed was the same lesson as A1 one level up: an
overlap between a section and a boundary edge must be clipped to the range
the *edge* covers. A sphere's seam and the far half of the same great circle
lie on one curve, and reading the whole curve as boundary made the meridian
opposite the seam disappear — which is the octant's own edge.

**What this says about the shape.** The plan proposed an `Interference`
table built before any splitting and a build phase with no classifier in it.
That structure was *not* built, and it is worth being exact about why. Four
of the table's five stated properties turned out to be the load-bearing
ones, and each is now enforced where the question is actually asked:
coincidence identified rather than rediscovered (the overlap correspondence
between two descriptions of one circle is now stated, and clipped to the
ranges the edges cover); an interference range that excludes what its edge
does not reach; degeneracies consumed before the stage that cannot handle
them, rather than after; and substitution by region rather than by name. The
fifth — a per-face state cache that lets the build phase read a
classification instead of probing for it — is not built, and the ray
classifier is still there. It was not what any of A1–A5 needed. It remains
the right thing to do for *speed*, and the case that would force it for
correctness has not been produced.

**A caution worth keeping.** Mature implementations of this design still
state that a face's classification can depend on which point of the face is
chosen, and still keep an open-ended list of configurations that defeat
them. Nothing here changes that. What closed is five named configurations,
each with a test that will say so if it reopens.

> **A7 — four faults behind one family of refusals.** **Closed.** A
> downstream front end brought four failures that all read as boolean gaps
> and were not all boolean gaps. Each is now pinned:
> `crates/ogeom/tests/mirrored.rs`,
> `crates/ogeom/tests/boolean_contact.rs` and `crates/ogeom/tests/groove.rs`.
>
> The two that were never the boolean's. A *copy* applied every stored
> placement and orientation twice — it recursed on occurrences that already
> carried the parent's, stored those, and composed the parent's on again. A
> doubled placement only misplaces; a doubled orientation takes a wire apart,
> since reversing one reverses the walk as well as each edge and only the
> second half of that survives a re-composition. A box cannot show it, having
> nothing stored reversed or displaced; a prism can, its near cap being the
> profile reversed and its far cap that same node moved. And
> `general_transformed_shape` converted every surface to a B-spline for any
> general transform, though a similarity is a placement and the kernel's own
> types say which is which. The cost was not the conversion but the silence
> after it: a restated plane is still that plane and no longer *says* so,
> `surface_surface` has no closed form for a plane against a patch, and the
> coincidence went unrecognized — so the marcher, which documents that it is
> deliberately not a coincidence detector, was handed a coincident pair, and
> a face lying on the other solid's boundary reached the classifier with no
> partner to compare sides against. Both surfaced as boolean refusals two
> stages downstream of themselves.
>
> The two that were. A revolved wall *parallel* to the turn's axis is a
> cylinder, and `revolution_over_edge` named only the plane its perpendicular
> case makes. A revolution is a surface nothing has a closed form for, so a
> ring could never melt against the cylinder it is and every plane meeting it
> was marched into a fitted curve where an exact ruling stood — which is what
> left a groove cut through a block unable to close. It is named now, on a
> frame whose chart is the revolution's own, so only `v` moves and it moves
> affinely; a profile running against the axis puts the chart's normal
> opposite the revolution's, and the face's flag carries that, which a
> rectangular profile exercises in a single ring. The oblique line and the
> meridian circle — the cone and the torus — are not named yet, and the cone
> would want its apex emitted as a pole edge.
>
> Last, the interior probes a piece offers. They varied the scanline's
> *height* and always took each interval's midpoint, so a piece symmetric
> about a chart-vertical line — a cylinder band, a revolved wall, a chart
> rectangle — offered the same column at every height, and a solid touching it
> down that column met every probe at once. The probes now vary along the
> scanline too, and are ranked so that a different column comes before a
> further height: the roomiest candidate in each distinct column first, then
> the rest by room. It is the same rule the quarter *heights* already stated,
> applied to the width.
>
> One thing this pass corrected rather than closed: which coincident piece
> stands in for which is decided by asking whether one piece's probe lies
> inside another piece's outline, and that was asked with the test for
> *strands*, which close only jointly. A piece's outline is a single ring that
> does not repeat its first point, so the segment closing it went uncounted
> and a point the ring plainly encloses came back outside whenever the ray
> crossed exactly there. `common` then kept both descriptions of one disk
> while `cut`, which never asks, closed — which is the shape the report had.
>
> **What remains, and what the refusal means.** A piece that lies on the other
> solid's boundary with no coincident partner face whose trim accepts its probe
> is still refused by name, and the four configurations above no longer reach
> it: each had a cause of its own, and none was that case. Two things do reach
> it. The one the text describes — genuine edge or vertex contact, confined to
> a line or a point, where there is no shared *region* for a partner to be
> found in and the two normals it would compare do not exist. And one it
> describes wrongly: a same-domain pair that is genuinely *not* analytic, a
> sheared copy sharing a plane with its original, where the surfaces are
> B-splines by necessity rather than by oversight and no closed form answers
> `Same` for them. That is the marcher being handed a coincident pair it says
> it must never be handed, and it reports as edge contact something that is
> nothing of the kind, after some seconds of marching. It wants a coincidence
> gate ahead of the marched fallback so the honest same-domain refusal already
> written in the filler is what comes back.


### B. Blending

**B1 — the marching blend.** **Done.** `march_blend` solves the section's two
endpoints directly, as the entry called for: unknowns `(u₁, v₁)` and
`(u₂, v₂)`, three equations saying the ball's centre is the same point from
either side, and a fourth tying the section to a guide. The tangency curves
come back **in the supports' own parameters**, which is the whole point — a
pcurve fitted through them is a pcurve of the curve rather than of a
projection of it, and the test holds that to `1e-12` rather than to a
projection tolerance.

Two things the entry did not say, both from building it.

*The guide's parameter is a fifth unknown, not a loop counter.* Four equations
in five unknowns is a curve, which is exactly what the shared walker follows —
so the step control, the stall reporting and the closure test are the
intersector's own, and the step is set by the sag of the tangency curve being
walked rather than by a guess at how finely to sample the guide. A guide that
closes on itself has its parameter wrapped and the *geometry* decides when the
march is back where it started.

*The seat is tried, not assumed.* Which side of each support the ball rolls on
is one sign apiece, and normals cannot tell a step from a slot; all four
combinations are solved and the one seating a ball that touches two distinct
points wins. A radius the corner cannot hold seats none of them, and that is
what the refusal says.

*The states,* which the formulation gives for free and which are now the stop
reasons a caller acts on: closed, left the first support, the second, both at
once — a corner rather than a run-out — ran past the guide, the section
collapsed because the radius is too large for the local geometry, stalled, and
ran out of steps. What remains is *unhooked*: carrying a blend that leaves one
support's boundary onto the face next door, which is a topological continuation
rather than a state of the solver, and is B2's business.

Measured on a cylinder square on a plane against the torus's own arithmetic,
and on one meeting a plane at twenty degrees — which has no closed form, which
is why it needs marching — against the ball's own definition at every station.

**B2 — corners where blends meet.** **Done for the three-blend vertex**, the
family's centre of gravity: three sequential fillets at a box corner —
which chain cleanly, an assumed impossibility that wasn't — take the A5
tool, and `b2_three_fillets_and_the_corner_tool_round_the_vertex` measures
the rounded vertex against a closed form derived by inclusion–exclusion:
within the corner cube every fillet prism's removal lies inside the
spike's, so V = 10³ − 3(1−π/4)r²(10−r) − r³ + πr³/6. What remains of the
family is the two-blend meeting and the N>3 setback vertex, neither of
which the corner construction covers and neither of which is blocked by
anything below it any more.

*The original entry, kept for its analysis:* The family is right: one blend running out at a
face's boundary, two meeting — resolved by intersecting them or by extending
both to a shared end — three, which is the setback vertex blend, and more,
which no single construction covers.

The three-blend equal-radius tool is right too, and A5 measures it: the corner
block running from the ball's centre out past the corner, *minus the ball*, is
exactly the leftover spike, because anything in that block further from the
centre than the ball lies outside one of the three fillet cylinders and was cut
away when that edge was rounded.

What blocks it is **not** A5, which is closed. It is that **a corner blend is
tangent to everything it rounds, by construction.** The ball sits at the
radius from each of the three faces, so it touches each of them — and it is
inscribed in each of the three fillet cylinders, so it touches those along
whole circles. Applying the tool therefore asks the boolean for tangential
contact in its hardest position: not a tangency in the interior of a face,
which the boolean carries today and a ball seated in a bore pins, but one at a
*vertex* of the tool's own spherical patch, where the octant's three corners
are exactly the three points at which it touches the three planes. That
surfaces as a dangling boundary strand rather than as a refusal by name, which
is worse than either answering or refusing.

So the remaining work is one named thing, and it belongs to the boolean rather
than to blending:

> **A6 — tangency at a face's own corner.** **Closed.** A face touching
> another at a point that is a *vertex* of its own boundary — what every
> corner blend asks for — now cuts, and
> `a6_the_corner_tool_cuts_through_its_own_tangencies` measures the rounded
> corner against its closed form. Three defects stacked under the one
> symptom, and each fix is a rule rather than a patch. A seam edge bounds a
> chart twice *only when the face wraps*: a sphere octant whose boundary is
> the seam meridian uses one column, and the far copy — fed in
> unconditionally before — dangled by construction; the decision is now
> made once at gather, by chart connectivity, and the arrangement, the trim
> tests and the rebuild all inherit it. A tangency was missed for want of a
> pcurve: the tangent circle between the corner ball and a fillet cylinder
> is a meridian through the sphere's poles, whose chart image has no closed
> form — so `touching_along` now inverts sample points through the
> surface's own closed forms when a pcurve is missing, sampling beside the
> degeneracies rather than on them. And the degeneracy splitter stopped
> only at the *face's* pole edges: a meridian section runs through both of
> a sphere's poles, and a face owning only the north one still cannot
> chart an arc that wraps through the south — the split points are now the
> *surfaces'* chart degeneracies, wherever the trim reaches.

The closure that finished it, after A6 and the two stop-welds: sewing and
validation now honour stated tolerances end to end. Edge fingerprints
carry the widest tolerance their edge and vertices state and compare
within it; vertex merging welds within stated reach and the survivor
widens to answer for what it absorbed; and the validity check accepts a
curve end within its *vertex's* tolerance, because that is the same
acceptance construction applies — a checker stricter than the builder
condemns what the builder rightly admitted and honestly recorded. The
principle underneath all three is the data model's: tolerances only
widen, are recorded where the disagreement was measured, and are then
*believed*.

Two things were built on the way and stand on their own. `march_blend` (B1)
is above. And `exact_pcurve` now sees through a **trim**: a trim says *where*
on a curve, not what it is, and since it shares its basis's parameter, its
pcurve is the basis's own pcurve trimmed the same way — on every surface, since
the answer does not depend on the surface. A fillet's end cap is a plane whose
bounding edges are trimmed curves, and the boolean had been refusing that
coincidence for want of putting a trimmed curve into a chart it plainly lies
in.

### C. Sweeping

**C1 — evolved shapes.** **Done.** `make_evolved` sweeps a profile along a
planar spine given as a wire or as a face: straight spine edges extrude the
profile as prisms, arcs turn it about their own axes as revolutions, and each
corner turns it about the corner through exactly the angle the spine turns
there — the join the 2D offset makes, for the same reason. Every piece is
exact; nothing is fitted.

Two things the entry did not say. The assembly is a *union*, not a sew:
consecutive pieces meet on the same placed profile, which is a coincident
face, and identifying that is the boolean's own work rather than something to
re-do here. And the face spine buys a volume by closing the *profile*, not by
capping the sweep afterwards — an open profile whose two ends reach the spine
face's plane is closed against it and swept as a section, so what comes back
is a solid by construction. An open profile along a wire spine has no plane to
close against, and is refused with the spine that would.

Measured against closed forms: a square spine is four runs and four
quarter-annulus wedges, and a bend is Pappus on its own annulus.

### D. Drawings

**D1 — marched silhouettes.** **Done.** A surface with no closed-form
silhouette — a torus, a spline — is *walked* rather than refused: the whole
content of a silhouette is one equation on the surface's own chart,

> `n(u, v) · d = 0`

and one equation in two unknowns is a curve, which is what the shared walker
follows. So it needed no new machinery, exactly as this entry said; the
condition is thirty lines and everything else is inherited.

Two things it did need, both found by measuring rather than by reading.

*The residual must be dimensionless.* Stated with the unnormalized `Sᵤ × Sᵥ`,
the condition's residual carries the surface's own scale, so on a torus of
radius eight the correction had to drive it below a *length* tolerance —
a demand on the angle eight times tighter than anything asked for — and the
walk answered by halving its step until it crawled. The unit normal, with its
own exact derivative, fixes it.

*The step control wants a length, and the face cannot supply one.* A full
torus's face is bounded by a seam and a single vertex, so its vertices'
bounding box is a point; fed that as a scale, the walk stepped round the whole
ring a ten-thousandth at a time and ran out. The extent now comes from the
surface itself, sampled over its own chart.

Measured against the torus's own equators seen down its axis — radii
`major ± minor`, which is arithmetic — and, from an oblique direction that has
no closed form at all, against the defining property itself: the surface's
normal is square to the view at every point of what comes back.

**D2 — the on-axis half-section.** Unblocked: A4 is closed, and the boolean
now cuts a bore on its own axis and reports the two rulings. What is left is
the drawing side — the half-section as a view, with its own hatching
convention.

### E. Documents

**E2 — datum targets and presentation PMI.** **Done.** `DatumTarget` with its
four kinds — point, line, rectangle, circle — placed and sized, tied to the
datum it establishes; and `Callout`, which is the drawing: the plane an
annotation is drawn in, the polylines that make its frame, leader and text,
and which semantic annotation it is a picture of.

Both directions in STEP, and the presentation half is read from NIST's own
annotated part rather than only from what this writer emits — which is what
found the shape of it: a callout's geometry is a *set* of tessellated curve
sets, nested and repositioned by a placement of its own, over one-based
indices into a coordinates list, and the link to the semantic annotation is
made by instance identity rather than by matching a name two annotations may
share. Twenty-three callouts come back with their planes and their 800-odd
drawn points; fourteen link to an annotation, and the ones that do not are
the file's own — a text note has nothing semantic behind it.

No style is written, which is a statement and not an omission: a style is
about rendering, and this kernel keeps no draughting style model to have one
from.

**E3 — saved views and standalone notes.** **Done.** `ogeom_doc::View` —
a name, a camera frame, an optional clipping plane, and indices into the
document's callouts, so restyling a callout restyles it in every view that
shows it — and `ogeom_doc::Note`, text with an author attached to a product
or the document. Both take part in undo like every other attribute, both
persist in the native document format (callouts now persist there too,
which the views made necessary), and views round-trip through STEP as the
named draughting models they are there: camera placement in, callout
membership by identity. A view without PMI is a camera bookmark and lives
in the native format; STEP carries views alongside the PMI they present.

### F. Exchange

**F1 — the `.brep` interchange text format.** **Done.** Read and written from
the format's own published specification: placements, the elementary and
spline curves and surfaces in both dimensions, the trimmed and offset forms
over them, and the whole topology encoding with its backward subshape
numbering, orientations and placement references. Cached meshes are parsed
and skipped; the per-record bookkeeping flags are read and dropped, since
they describe the writing session rather than the shape.

One thing the reader does not take on faith is whether an edge's
representations agree on parameterization — that claim is what everything
downstream relies on, so it is measured rather than believed, and a file can
come back with it *established* where its writer never made it. A drilled
block round-trips to a fixed point, byte for byte; a drum and a ball come
back as a cylinder, a sphere and two degenerate poles; and a file written by
hand against the specification reads as the square it describes.

**F2 — IGES, both directions.** **Built**, from the published record layout,
and held to eight measured round trips: planes, a periodic cylinder wall with
its seam, a doubly periodic torus, a seam-only sphere whose poles the format
cannot spell, a boolean result, a spline-walled prism through the rational
B-spline entities, inch-unit scaling, and a refusal by name. The reader takes
both kinds of file — manifold solid B-rep objects bottom-up, and the older
surface files as trimmed faces sewn into shells, solids where they close —
re-deriving edge ranges on this kernel's own parameterizations exactly as the
STEP reader does. Twenty-nine entity types translate; the live figure and the
refused remainder are the parity ledger's `io.iges` row.

Two findings came out of building it that reach past IGES. The exchange
writers paired both of an edge's vertices with the edge's own placement,
which quietly welds an instanced vertex — a prism's top corner is its bottom
corner, moved — to its other placement; both writers now resolve each
vertex's composed placement. And the fitted-pcurve machinery the STEP reader
grew is now `ogeom-io`'s shared `pcurves` module, because the second reader
needed exactly the first one's policy.

**F3 — reading glTF.** **Done.** `read_glb` and `read_gltf`, with the whole
indirection honoured because a writer chooses it and a reader does not get to
assume: accessors over buffer views over buffers, byte strides, byte offsets
at both levels, all six component types, `normalized` integers scaled into
fractions, and the sparse block applied over the base it sits on. The scene's
node hierarchy is walked and composed, stated as a matrix or as
translation–rotation–scale, and normals come through the inverse transpose so
an uneven scale leaves them normal. A `.gltf` document's data uris are decoded;
an external file reference is refused, as are a non-triangle primitive mode,
a Draco payload and a node that is its own descendant, each by name.

It needed a JSON parser, which is now `ogeom_io::json` — the grammar and
nothing else, no dependency, written because glTF's structure is JSON and
`cargo build` still needs a Rust toolchain and nothing more.

**F5 — a closed spline wall through exchange.** **Done**, and the fix was
one distinction: the fitted-pcurve unwrap engaged on *periodicity*, and a
skinned loft's wall is a clamped B-spline that closes on itself without
being periodic — projections near the joining column land in either copy,
both right pointwise, and only continuity chooses. Unwrapping now engages
on *closure*, in the shared module, so both exchange readers learned it at
once — which `f5_a_closed_spline_wall_survives_both_formats` pins by
demanding the two formats agree with each other a million times tighter
than either must agree with the original.

**F4 — SAT, X\_T, JT.** Refused for want of public documentation. These are
proprietary formats whose specifications are not published; implementing them
would mean reverse-engineering files rather than reading a standard. If a
specification becomes available the refusal lifts, and until then the honest
answer is this row.

### H. Defeaturing

**H1 — removing a set of faces.** **Built**, as
`ogeom_bool::remove_faces(model, solid, &[Shape], tol)`, and the two wounds
close differently. A feature whose rim is an inner loop of a surviving face —
a bore, a boss, a mid-face pocket — is wire surgery: the survivor is rebuilt
without the rim wire, nothing re-intersected, and the drilled block comes
back to its exact volume with **no overshoot at all**, because the face-set
approach fills nothing — the boundary is resewn, and the sliver-band problem
the tool-based approach paid ten microns to avoid never arises. A band
feature — a fillet or chamfer along an edge — is the re-intersection case:
the two side surfaces recover the edge the blend replaced, the end faces'
own edges extend along their own curves to the corners the recovered edge
pierces, and both the chamfered and the filleted box come back to the sharp
box *exactly*, six faces, twelve edges, eight vertices.

One ordering lesson worth keeping: the ends rebuild before the sides,
because extending a cap's dangling edge decides that edge for the side that
shares it, and sewing rejoins them on the one node. Multiple simultaneous
bands, bands meeting at corners (the B2/A6 family), and spline-surfaced
neighbours are refused by name — the parity ledger's `bool.defeaturing`
restriction is the live list.

What follows is the original entry, kept because the overshoot finding is a
real measurement about the tool-based road not taken.

The operation was to be
`remove_faces(model, solid, &[Shape], tol)` in **`ogeom-bool`**: given faces to
delete, extend the neighbours that surrounded them, re-intersect the extensions
against each other, and sew the result back into the shell. Face extension plus
re-intersection is the fuse machinery already there; what is missing is the
driver that decides which neighbours to extend and how far.

It belongs to the boolean and not to a recognizer. A caller hands over faces —
what those faces *mean* is the caller's business, and the operation must work on
a solid whose history is gone. There was a `remove_feature(model, solid,
&Feature, tol)` here that dispatched on a recognized feature and rebuilt the
volume that feature described. That is a different operation with a different
input, and it left with the recognizer; salvaging its code would have preserved
the half that does not generalise.

One finding from it is worth keeping, because a face-set implementation will hit
it again the moment it builds a tool that meets the solid at an opening:

> A filling tool flush with the faces it meets is a coincidence at every opening
> at once, and the boolean does not assemble it. So the tool overshoots. But the
> overshoot cannot be small, and this is the part that is not obvious: a margin
> leaves a sliver band standing past the opening, and that band's interior probes
> must be *decisively* outside the part, or the exact classifier finds every ray
> from them grazing the face they sit against, exhausts its whole fan of
> directions and answers `On` the slow way. A micron of overshoot is inside the
> band the classifier reads as "on the boundary" and costs fifty seconds on a
> part that takes a fifth of one otherwise; ten microns is outside it and costs
> nothing. A hundred thousand confusion tolerances — ten microns at millimetre
> tolerances — was the working figure, and the restored solid is larger than the
> original by that times the openings' area.

Any caller whose tolerance is tighter than the overshoot is the case to watch.

### I. Canonical simplification

**I1 — recognizing that exact geometry is secretly analytic.** **Done.**
`ogeom_heal::canonical_simplify` samples each free-form surface on its own
chart with its own normals, proposes plane, sphere, cylinder or cone from
the classical estimators — the mean normal, the least-squares meeting of
normal lines, the direction the normals avoid, the linear taper of radius
against height — and accepts only when *every* sample verifies at the
caller's tolerance, the certificate being the worst deviation actually
measured. Free-form curves get the same treatment on the way, because a
rim spelt as a B-spline that is exactly a circle must become the circle
before the analytic surface has anything to project in closed form. A
nurbsed drum comes back a cylinder at 1.8e-15 with its volume unchanged to
the last bit; a skinned loft stays what it is. The reference's set —
plane, cylinder, cone, sphere — is matched exactly; a torus candidate is
out on both sides.

### J. The medial axis

**J1 — the medial axis of a planar region.** **Built for the convex
polygonal case**, exactly: the shrinking-polygon construction — every edge
inward at unit speed, every vertex riding its angular bisector, each event
retiring an edge and starting a branch — with convexity as the honesty
condition, because it is what excludes the split events the construction
does not have. Held to closed forms: a rectangle's four diagonals and roof
line to 1e-9, and a 3–4–5 triangle's branches meeting at the incenter with
the inradius (a+b−c)/2 = 1 as the deepest clearance. Holes, reflex corners
and arcs are refused by name; the parity row's restriction is the live
worklist should a caller ever need the general region.

## Decisions, not gaps

These are settled. They are here so nobody reopens them by accident.

- **A pcurve with no closed form is `None`, not a fit.** An exact curve
  carrying a fitted pcurve would be a curve whose two descriptions disagree
  by an amount nothing on it records. The consumer that needs one marches the
  pair instead.
- **A closed exact section partly outside a surface's extent is kept whole.**
  The restriction that matters is the face's trim, which is the boolean's own
  2D stage; the surface extent is only a parameterization window, and cutting
  there would split a curve where no boundary exists.
- **Scaled placements in the boolean are refused, with instructions.** A
  scale changes a surface's parameterization out from under its pcurves. Bake
  it first — `baked_shape` does exactly that, and the boolean calls it.
- **The crossing walker refuses tangential contact.** It is a crossing
  walker; the tangential walker owns that case and the section pipeline
  routes to it. The refusal stays pinned.
- **Bi-tangent construction is subsumed** — the 2D repertoire in 2D, the
  blend family's own envelope in 3D.
- **Glue is subsumed** by the boolean's same-domain unification, which
  already skips nothing it needs and unifies what glue would.

## Order

1. ~~**A** — the interference table.~~ **Done.** Five named failures closed,
   B2 and D2 unblocked.
2. **C**, **E2**, **F3** — contained pieces, any order, each a stone.
3. ~~**The walker abstraction.**~~ **Done.** `ogeom_intersect::walk`: a
   `Condition` in `n` unknowns and `n − 1` equations, and one walk over it.
   The shortfall is the point — the solution set of `n − 1` equations in `n`
   unknowns *is* a curve — and the walker supplies the missing equation
   itself, a plane across the direction of travel. The direction comes free
   as the Jacobian's null vector, so a condition need not know its own
   tangent formula. Step control, stall reporting and closure are written
   once. The intersector's own walk now goes through it, which is what says
   it is general rather than merely present.

   One thing the abstraction had to be taught. A null vector's sign is
   whatever the arithmetic gave it, so the walker keeps its own heading — but
   the intersector's tangent is the cross product of two normals, whose sign
   is the surfaces' own, and whose *flip* at a tangency is what stops the
   march. Turned quietly back round, two thin curves through two touching
   points came back as one confident loop lying on neither. So a condition
   declares whether its tangent's sign is its own, and that declaration is
   load-bearing.
4. ~~**B1** — the marching blend.~~ **Done.** **B2**'s corner family is owed
   still, and on **A6** — tangency at a face's own corner — rather than on A5.
5. ~~**D1** — marched silhouettes.~~ **Done**, and they were indeed a second
   condition for a walker that already existed.
6. **F2** — IGES.
7. **H1** — removing a set of faces.

**F4** has no scheduled slot: it needs a document nobody has published.

**A debt paid.** Fifteen places in `crates/` used to refer the reader to
"the deferred table", which no longer existed. All fifteen now point at
`docs/PARITY.md` rows by id, or state plainly that an earlier plan owed the
thing and it was delivered. Six of them are inside error strings, so a
refusal now names a row a reader can actually open.
