# Plan

This file lists what is left to build, why each piece is not built yet, and how
each one gets done. It replaces the old scope map, which had become a record of
finished work with the unfinished parts scattered through it.

**Already built:** topology and the geometry vocabulary, the intersector, the
boolean, healing, blending, offsets and sweeps, tessellation, drawings, the
document layer, and exchange for STEP, STL, DXF, glTF, OBJ, PLY, VRML and 3MF.

**How the remainder is tracked.** The remaining work is *derived*, not
remembered:

- `docs/PARITY.md` audits every public header of the reference kernel's
  modelling modules against 98 named capabilities.
- The gate in `tools/check.sh` holds the audit to its evidence.
- Every `absent` and `partial` row in the ledger points back to a section here.

The ledger records *what* is missing. This file records *how*: what each gap
needs and why it is not built yet. An item missing from both is a gate failure,
not an oversight.

## How this project works

These rules are enforced. Every one of them has caused a patch to be rejected.

- **Measured, not asserted.** A claim about geometry is backed by a test that
  measures it against a closed form, a published value, or an independent
  computation. "It looks right" is not a result.
- **An approximation must say so.** Sampling, fitting and tolerances are fine
  when they are *stated*. A number that hides its own error is not.
- **Refusals are by name, and pinned.** When the kernel cannot do something,
  the error text says which thing and why, and a test holds it to that. A silent
  wrong answer is the only unacceptable outcome.
- **No stubs.** Either a thing is implemented, or it is in this document.
- **Scope is parity with the reference kernel's modelling modules, and nothing
  else.** `docs/SCOPE.md` is normative and says how to decide a case. Parity is
  about capability, not structure. Usage data sets the order of work, never its
  bounds.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against, bundles
  or imports another kernel, and the design is worked out here, not mirrored.
  The field's vocabulary is used throughout, because that is how the field
  talks about itself.

## The remaining work

### A. The boolean's interference table (A1 to A5 closed)

A1 to A5 were five unrelated-looking gaps, treated as one family. Each is
pinned by a measured test in `crates/ogeom/tests/boolean_interference.rs`. The
tests were written before the fix and check volumes against closed forms, not
just that no error came back. All five pass.

**The gaps, and what each turned out to be.**

*A1: a tool flush with the part.* Refill a bore with the cylinder that cut it,
with ends flush with the faces it broke. Two causes:

- An overlap between two *curves* was read as an overlap between the two
  *edges*. A hole's arc and the disc filling it matched each other at the
  circle's ends instead of at the arc's ends, so the two sides were sewn
  against different subdivisions.
- The rule for which copy of a coincident pair to keep was based on which
  argument it came from, not on region. So the tool's cap, which fills a hole
  the part has no face for, was dropped, and the result had a hole in it.

Fix: substitution is now by containment. A piece is dropped only where a piece
of the other argument actually occupies the same place.

*A2: a section through a chart pole.* A plane cuts a ball through its own axis.
Two causes:

- The section's ends land exactly on the poles. An *open* section's end is its
  domain end, but the rebuild was wrapping it round to the domain start: the
  far pole, eight radii away.
- The section was being marched, because a meridian's chart image has no
  closed form (its longitude jumps half a turn at each pole). Each *half* does
  have one, and it is a straight line. Write the sphere's axis as
  `Z = cos α·X + sin α·Y` in the circle's own frame. The latitude is then
  `asin(cos(t − α))`, which for `t − α ∈ [0, π]` is exactly `π/2 − (t − α)`:
  affine in the circle's own parameter.

Fix: the section is cut at the poles named by the face's own degenerate edges,
and each half comes back exact, not fitted. The half ball is now measured
against `2πr³/3`, not against a marched approximation of it.

*A3: a ball on a box's corner vertex.* Three placements:

- tangent at the vertex from outside;
- centred on the vertex (the closed form is an octant);
- with the vertex exactly on the sphere, so that all three of the box's faces
  cut the ball through one point.

The last is a triple point, and the problem was the seam. A crossing found on a
periodic curve comes back in `[0, 2π)`, but the seam edge covers
`[-π/2, π/2]`. Every crossing on the seam was therefore out of range and
discarded, and the seam was never split where the chain of section arcs met it.

*A4: a section through tangential contact.* A plane through a bore's axis meets
the wall along its rulings, and one of those rulings *is* the wall's own seam.
The rule "this section is already boundary" was applied globally, so the
ruling also vanished from the plane. On the plane it is not boundary at all: it
is the curve separating the two halves the section leaves.

Fix: the rule is now recorded per face. A section is dropped only from the face
that already carries it, as a boundary edge or as a contact strand. This also
enables drawing feature D2.

*A5: a shell around a three-cylinder tip.* The corner block minus the ball,
which is the tool the three-blend corner needs.

- The block's three far planes are *tangent* to the ball, exactly at the three
  vertices next to the corner.
- The three near planes cut it in a meridian, a meridian along the seam, and
  the equator.

The fix is the A1 lesson one level up: an overlap between a section and a
boundary edge must be clipped to the range the *edge* covers. A sphere's seam
and the far half of the same great circle lie on one curve. Treating the whole
curve as boundary made the meridian opposite the seam disappear, and that
meridian is the octant's own edge.

**What this says about the design.** The plan proposed an `Interference` table
built before any splitting, and a build phase with no classifier. That
structure was *not* built. The reason:

- Four of the table's five stated properties turned out to be the essential
  ones, and each is now enforced where the question is actually asked:
  - coincidence is identified, not rediscovered (the overlap correspondence
    between two descriptions of one circle is now stated, and clipped to the
    ranges the edges cover);
  - an interference range excludes what its edge does not reach;
  - degeneracies are handled before the stage that cannot handle them, not
    after;
  - substitution is by region, not by name.
- The fifth, a per-face state cache that lets the build phase read a
  classification instead of probing for it, is not built. The ray classifier is
  still there. None of A1 to A5 needed it. It is still the right thing to do
  for *speed*, but no case has yet been found that needs it for correctness.

**A caution.** Mature implementations of this design still say that a face's
classification can depend on which point of the face is chosen, and still keep
an open-ended list of configurations that defeat them. Nothing here changes
that. What closed is five named configurations, each with a test that will
catch it if it reopens.

> **A7: four faults behind one family of refusals.** **Closed.** A downstream
> front end reported four failures. All looked like boolean gaps; not all were.
> Each is now pinned in `crates/ogeom/tests/mirrored.rs`,
> `crates/ogeom/tests/boolean_contact.rs` and `crates/ogeom/tests/groove.rs`.
>
> **Two that were not the boolean's fault.**
>
> - *A copy applied every stored placement and orientation twice.* It recursed
>   on occurrences that already carried the parent's placement, stored those,
>   and then composed the parent's on again. A doubled placement only misplaces
>   a shape. A doubled orientation breaks a wire apart: reversing a wire
>   reverses the walk as well as each edge, and only the second half of that
>   survives re-composition. A box cannot show the bug, because nothing in it is
>   stored reversed or displaced. A prism can: its near cap is the profile
>   reversed, and its far cap is the same node moved.
> - *`general_transformed_shape` converted every surface to a B-spline for any
>   general transform*, even though a similarity is a placement and the
>   kernel's own types say which is which. The conversion itself was not the
>   cost. The cost was that a restated plane is still that plane but no longer
>   *says* so. `surface_surface` has no closed form for a plane against a patch,
>   so the coincidence went unrecognized. The marcher (which is documented as
>   deliberately not a coincidence detector) was handed a coincident pair, and
>   a face lying on the other solid's boundary reached the classifier with no
>   partner to compare sides against.
>
> Both showed up as boolean refusals two stages after the real fault.
>
> **Two that were the boolean's.**
>
> - *Revolved surfaces were not recognized as the analytic surfaces they are.*
>   A revolved wall *parallel* to the axis is a cylinder, but
>   `revolution_over_edge` only named the plane produced by the perpendicular
>   case. Nothing has a closed form against a general surface of revolution, so
>   a ring could never merge with the cylinder it is, and every plane meeting it
>   was marched into a fitted curve where an exact ruling existed. That is why a
>   groove cut through a block could not close. Now:
>   - The cylinder is named, on a frame whose chart is the revolution's own, so
>     only `v` changes and it changes affinely. A profile running against the
>     axis puts the chart's normal opposite the revolution's; the face's flag
>     records that. A rectangular profile tests this in a single ring.
>   - The oblique line gives a cone. It is measured from the profile's own low
>     end, so its reference radius is one the profile states. Its window stops
>     at the apex, because past it the radius goes negative and that is the
>     *other* nappe.
>   - The meridian circle gives a torus whose `v` is the angle round the tube,
>     which the circle's own parameter gives up to a full turn and a sign.
>   - A revolved cone and a revolved torus now use the same surfaces the
>     primitives build them from, and say so.
> - *Interior probes.* They varied the scanline's *height* but always took each
>   interval's midpoint. So a piece symmetric about a vertical line in its chart
>   (a cylinder band, a revolved wall, a chart rectangle) offered the same column
>   at every height, and a solid touching it along that column hit every probe
>   at once. The probes now also vary along the scanline, ranked so a different
>   column comes before another height: the roomiest candidate in each distinct
>   column first, then the rest by room. This is the rule the quarter *heights*
>   already used, applied to the width.
>
> **One correction (not a closure).** Which coincident piece stands in for which
> is decided by asking whether one piece's probe lies inside another piece's
> outline. That was asked with the test meant for *strands*, which close only
> jointly. A piece's outline is a single ring that does not repeat its first
> point, so the closing segment was not counted, and a point clearly inside the
> ring came back outside whenever the ray crossed exactly there. `common` then
> kept both descriptions of one disk, while `cut`, which never asks, closed.
> That matches the reported shape.
>
> **Coincidence where there is no closed form.** A sheared copy sharing a plane
> with its original used to reach the same refusal, but it is a different case.
> The surfaces are B-splines by necessity, not by oversight, and nothing answers
> `Same` for two patches. The marcher (documented as not a coincidence detector,
> and with no crossing to trace over a coincident pair) was handed the pair
> anyway. It returned a section made of noise, and the boolean refused for edge
> or vertex contact some seconds later.
>
> Coincidence is now measured, with limits:
>
> - only for pairs the closed forms have already declined;
> - only over the region the two surfaces actually share, sampled on the
>   smaller window (a plane's own window runs for a billion units, so a grid
>   over it samples nothing);
> - a sample whose foot lands on the rim of the other surface is skipped, not
>   counted against, because a distance measured to a rim is about the window,
>   not the surface.
>
> The test can only fail one way. A crossing pair puts interior samples well
> off the other surface, so it cannot pass, and a pair the test cannot resolve
> is marched as before. The sheared copy now gets the refusal the filler already
> had written for it, in a few hundred microseconds instead of forty seconds.
>
> Resolving that case, instead of refusing it, needs a same-domain contact
> whose edges project into a *fitted* chart. That is the pcurve-fitting
> question, not this one.
>
> **What remains, and what the refusal means.** A piece that lies on the other
> solid's boundary, with no coincident partner face whose trim accepts its
> probe, is still refused by name. What reaches that refusal now is the
> configuration its text describes: genuine edge or vertex contact, confined to
> a line or a point. There is no shared *region* in which to find a partner, and
> the two normals it would compare do not exist.
>
> **Code in the wrong crate.** The measurement above belongs in the
> intersector, which owns `Same`. But the projection it needs lives in
> `ogeom-algo`, which depends on the intersector, not the other way round. So it
> is done in the boolean, next to the apartness test it mirrors. Moving
> `project_on_surface` down into `ogeom-geom` would let the intersector answer
> for itself, and every caller of `intersect_surfaces` would benefit, not only
> this one.


### B. Blending

**B1: the marching blend.** **Done.** `march_blend` solves directly for the
section's two endpoints, as the entry specified:

- unknowns `(u₁, v₁)` and `(u₂, v₂)`;
- three equations saying the ball's centre is the same point from either side;
- a fourth tying the section to a guide.

The tangency curves come back **in the supports' own parameters**, which is
the point. A pcurve fitted through them is a pcurve of the actual curve, not of
a projection of it, and the test holds that to `1e-12`, not to a projection
tolerance.

Three things the entry did not say, learned while building it:

- *The guide's parameter is a fifth unknown, not a loop counter.* Four
  equations in five unknowns define a curve, which is exactly what the shared
  walker follows. So step control, stall reporting and the closure test are
  the intersector's own. The step size is set by the sag of the tangency curve
  being walked, not by a guess at how finely to sample the guide. For a closed
  guide, the parameter wraps and the *geometry* decides when the march is back
  where it started.
- *The seat is tried, not assumed.* Which side of each support the ball rolls
  on is one sign per support, and normals cannot tell a step from a slot. All
  four combinations are solved, and the one that seats a ball touching two
  distinct points wins. If the corner cannot hold the radius, none of them
  seats, and the refusal says so.
- *The stop states* come free with the formulation and are now the stop reasons
  a caller acts on: closed; left the first support; left the second; left both
  at once (a corner, not a run-out); ran past the guide; section collapsed
  (radius too large for the local geometry); stalled; ran out of steps. Still
  missing is *unhooked*: carrying a blend that leaves one support's boundary
  onto the neighbouring face. That is a topological continuation, not a solver
  state, and it belongs to B2.

Measured on a cylinder square to a plane, against the torus's own arithmetic;
and on a cylinder meeting a plane at twenty degrees (no closed form, which is
why it needs marching), against the ball's own definition at every station.

**B2: corners where blends meet.** **Done for the three-blend vertex**, the
core of the family.

- Three sequential fillets at a box corner chain cleanly (this had been assumed
  impossible, and is not).
- They then take the A5 tool.
- `b2_three_fillets_and_the_corner_tool_round_the_vertex` measures the rounded
  vertex against a closed form derived by inclusion and exclusion. Within the
  corner cube, each fillet prism's removal lies inside the spike's, so
  V = 10³ − 3(1−π/4)r²(10−r) − r³ + πr³/6.

Still to do in this family: the two-blend meeting and the N>3 setback vertex.
The corner construction covers neither, and nothing below them blocks either
any more. Tracked as issue #18. The defeaturing corner family, which shares
this machinery, is #19.

*The original entry, kept for its analysis.*

The family is: one blend running out at a face's boundary; two blends meeting
(resolved by intersecting them, or by extending both to a shared end); three,
which is the setback vertex blend; and more, which no single construction
covers.

The three-blend equal-radius tool is also right, and A5 measures it. The corner
block running from the ball's centre out past the corner, *minus the ball*, is
exactly the leftover spike. Anything in that block further from the centre than
the ball lies outside one of the three fillet cylinders, and was cut away when
that edge was rounded.

What blocked it was **not** A5, which is closed. It is that **a corner blend is
tangent to everything it rounds, by construction.**

- The ball sits at the radius from each of the three faces, so it touches each
  of them.
- It is inscribed in each of the three fillet cylinders, so it touches them
  along whole circles.

Applying the tool therefore asks the boolean for tangential contact in its
hardest position. Not a tangency inside a face (the boolean handles that today,
and a ball seated in a bore pins it), but one at a *vertex* of the tool's own
spherical patch: the octant's three corners are exactly the three points where
it touches the three planes. This showed up as a dangling boundary strand, not
as a refusal by name, which is worse than either answering or refusing.

So the remaining work was one named thing, and it belonged to the boolean, not
to blending:

> **A6: tangency at a face's own corner.** **Closed.** A face touching another
> at a point that is a *vertex* of its own boundary (which every corner blend
> needs) now cuts. `a6_the_corner_tool_cuts_through_its_own_tangencies`
> measures the rounded corner against its closed form. Three defects sat under
> the one symptom, and each fix is a rule, not a patch:
>
> - *A seam edge bounds a chart twice only when the face wraps.* A sphere octant
>   whose boundary is the seam meridian uses one column. The far copy was
>   previously added unconditionally, and so dangled by construction. The
>   decision is now made once, at gather time, by chart connectivity, and the
>   arrangement, the trim tests and the rebuild all inherit it.
> - *A tangency was missed for lack of a pcurve.* The tangent circle between the
>   corner ball and a fillet cylinder is a meridian through the sphere's poles,
>   whose chart image has no closed form. `touching_along` now inverts sample
>   points through the surface's own closed forms when a pcurve is missing,
>   sampling beside the degeneracies, not on them.
> - *The degeneracy splitter stopped only at the face's pole edges.* A meridian
>   section runs through both of a sphere's poles, and a face that owns only the
>   north pole still cannot chart an arc that wraps through the south. The split
>   points are now the *surfaces'* chart degeneracies, wherever the trim
>   reaches.

The final step, after A6 and the two stop-welds: sewing and validation now
respect stated tolerances end to end.

- Edge fingerprints carry the widest tolerance their edge and vertices state,
  and compare within it.
- Vertex merging welds within the stated reach, and the surviving vertex's
  tolerance grows to cover what it absorbed.
- The validity check accepts a curve end within its *vertex's* tolerance,
  because construction applies the same acceptance. A checker stricter than the
  builder would condemn what the builder correctly accepted and recorded.

The principle behind all three comes from the data model: tolerances only
grow, are recorded where the disagreement was measured, and are then *trusted*.

Two things built along the way stand on their own:

- `march_blend` (B1, above).
- `exact_pcurve` now sees through a **trim**. A trim says *where* on a curve,
  not what the curve is, and it shares its basis's parameter. So its pcurve is
  the basis's own pcurve, trimmed the same way, on every surface. A fillet's end
  cap is a plane whose bounding edges are trimmed curves, and the boolean had
  been refusing that coincidence because it could not put a trimmed curve into
  a chart it obviously lies in.

### C. Sweeping

**C1: evolved shapes.** **Done.** `make_evolved` sweeps a profile along a
planar spine, given as a wire or as a face:

- straight spine edges extrude the profile as prisms;
- arcs turn it about their own axes as revolutions;
- each corner turns it about the corner by exactly the angle the spine turns
  there (the same join the 2D offset uses, for the same reason).

Every piece is exact; nothing is fitted.

Two things the entry did not say:

- *The pieces are joined by a union, not by sewing.* Consecutive pieces meet on
  the same placed profile, which is a coincident face, and identifying that is
  the boolean's job; it is not redone here.
- *A face spine gives a volume by closing the profile, not by capping the sweep
  afterwards.* An open profile whose two ends reach the spine face's plane is
  closed against that plane and swept as a section, so the result is a solid by
  construction. An open profile along a wire spine has no plane to close
  against, and is refused, naming the spine that would be needed.

Measured against closed forms: a square spine gives four runs and four
quarter-annulus wedges, and a bend obeys Pappus on its own annulus.

### D. Drawings

**D1: marched silhouettes.** **Done.** A surface with no closed-form silhouette
(a torus, a spline) is now *walked* instead of refused. A silhouette is one
equation on the surface's own chart:

> `n(u, v) · d = 0`

One equation in two unknowns is a curve, which is what the shared walker
follows. So, as this entry predicted, it needed no new machinery: the condition
is thirty lines and the rest is inherited.

It did need two things, both found by measurement:

- *The residual must be dimensionless.* Written with the unnormalized
  `Sᵤ × Sᵥ`, the residual carries the surface's own scale. On a torus of radius
  eight, the correction had to push it below a *length* tolerance, which asks
  for an angle eight times tighter than intended, and the walk kept halving its
  step until it crawled. Using the unit normal, with its own exact derivative,
  fixes this.
- *The step control needs a length, and the face cannot supply one.* A full
  torus face is bounded by a seam and a single vertex, so the bounding box of
  its vertices is a point. Given that as a scale, the walk went round the ring a
  ten-thousandth at a time and ran out of steps. The extent now comes from the
  surface itself, sampled over its own chart.

Measured against the torus's equators seen down its axis (radii
`major ± minor`, which is plain arithmetic). From an oblique direction, which
has no closed form, it is measured against the defining property: the surface
normal is perpendicular to the view at every returned point.

**D2: the on-axis half-section.** **Done.**

- The boolean side is A4: it cuts a bore on its own axis and reports the two
  rulings.
- The drawing side is `ogeom_hlr::half_section`. The plane's frame states the
  whole convention: `+z` is removed over the `+x` half, and the frame's `y`
  axis is the split line.
- `ogeom_hlr::hatch` clips the hatching strokes to the section loops by the
  even-odd rule, so holes interrupt the hatching exactly as they interrupt the
  material.

The bored drum's half-section measures its wall's half-area, and every stroke
lands in material. This was issue #32.

### E. Documents

**E2: datum targets and presentation PMI.** **Done.**

- `DatumTarget`, with its four kinds (point, line, rectangle, circle), placed
  and sized, and tied to the datum it establishes.
- `Callout`, which is the drawing: the plane an annotation is drawn in, the
  polylines that form its frame, leader and text, and which semantic
  annotation it depicts.

Both are read and written in STEP. The presentation half is tested on NIST's
own annotated part, not only on this writer's output, and that is how its
structure was found:

- a callout's geometry is a *set* of tessellated curve sets, nested and
  repositioned by its own placement, over one-based indices into a coordinates
  list;
- the link to the semantic annotation is by instance identity, not by matching
  a name that two annotations may share.

Twenty-three callouts come back with their planes and their 800 or so drawn
points. Fourteen link to an annotation. The ones that do not are the file's
own: a text note has nothing semantic behind it.

No style is written. This is a decision, not an omission: a style is about
rendering, and this kernel has no draughting style model to take one from.

**E3: saved views and standalone notes.** **Done.**

- `ogeom_doc::View`: a name, a camera frame, an optional clipping plane, and
  indices into the document's callouts. Restyling a callout therefore restyles
  it in every view that shows it.
- `ogeom_doc::Note`: text with an author, attached to a product or to the
  document.

Both take part in undo like every other attribute, and both persist in the
native document format. (Callouts now persist there too; views required it.)
Views round-trip through STEP as the named draughting models they are in that
format: camera placement is read, and callout membership is matched by
identity. A view without PMI is a camera bookmark and lives in the native
format. STEP carries views alongside the PMI they present.

### F. Exchange

**F1: the `.brep` interchange text format.** **Done.** Read and written from
the format's published specification. Covered:

- placements;
- elementary and spline curves and surfaces, in both dimensions;
- the trimmed and offset forms over them;
- the whole topology encoding, with its backward subshape numbering,
  orientations and placement references.

Cached meshes are parsed and skipped. The per-record bookkeeping flags are read
and dropped, since they describe the writing session, not the shape.

The reader does not trust a file's claim that an edge's representations agree
on parameterization. Everything downstream relies on that claim, so it is
measured, and a file can come back with it *established* even where its writer
never set it.

Results: a drilled block round-trips to a fixed point, byte for byte; a drum
and a ball come back as a cylinder, a sphere and two degenerate poles; and a
file written by hand against the specification reads as the square it
describes.

**F2: IGES, both directions.** **Built** from the published record layout, and
held to eight measured round trips:

1. planes;
2. a periodic cylinder wall with its seam;
3. a doubly periodic torus;
4. a seam-only sphere whose poles the format cannot express;
5. a boolean result;
6. a spline-walled prism, through the rational B-spline entities;
7. inch-unit scaling;
8. a refusal by name.

The reader handles both kinds of file: manifold solid B-rep objects, read
bottom-up; and the older surface files, read as trimmed faces sewn into shells,
and into solids where they close. It re-derives edge ranges on this kernel's
own parameterizations, exactly as the STEP reader does. Twenty-nine entity
types translate. The live figure and the refused remainder are in the parity
ledger's `io.iges` row.

Coverage grows on demand, not on speculation. When a real file contains an
entity outside the set, an issue names the entity and the file, and the
translation is built against that evidence. Issue #27 records this policy. The
table in the module is the single source; `tools/parity.py exchange`
regenerates the figure from it.

Two findings from building it apply beyond IGES:

- The exchange writers paired both of an edge's vertices with the edge's own
  placement. This quietly welds an instanced vertex to its other placement (a
  prism's top corner is its bottom corner, moved). Both writers now resolve
  each vertex's composed placement.
- The fitted-pcurve machinery the STEP reader had grown is now `ogeom-io`'s
  shared `pcurves` module, because the second reader needed exactly the same
  policy.

**F3: reading glTF.** **Done.** `read_glb` and `read_gltf` follow the whole
indirection chain, because the writer chooses it and the reader cannot assume
anything:

- accessors over buffer views over buffers;
- byte strides, and byte offsets at both levels;
- all six component types;
- `normalized` integers scaled into fractions;
- the sparse block applied over its base.

The scene's node hierarchy is walked and composed, whether a node states a
matrix or translation, rotation and scale. Normals go through the inverse
transpose, so they stay normal under non-uniform scale. A `.gltf` document's
data URIs are decoded. The following are refused, each by name: an external
file reference, a non-triangle primitive mode, a Draco payload, and a node that
is its own descendant.

It needed a JSON parser, which is now `ogeom_io::json`: the grammar and nothing
else, with no dependency. It was written because glTF's structure is JSON and
`cargo build` must still need only a Rust toolchain.

**F5: a closed spline wall through exchange.** **Done.** The fix was one
distinction. The fitted-pcurve unwrap was triggered by *periodicity*. A skinned
loft's wall is a clamped B-spline that closes on itself without being periodic.
Projections near the joining column land in either copy; both are correct
pointwise, and only continuity decides between them. Unwrapping is now
triggered by *closure*, in the shared module, so both exchange readers gained
it at once. `f5_a_closed_spline_wall_survives_both_formats` pins this by
requiring the two formats to agree with each other a million times more
tightly than either must agree with the original.

**F4: SAT, X\_T, JT.** Refused for lack of public documentation. These are
proprietary formats with unpublished specifications. Implementing them would
mean reverse-engineering files instead of reading a standard. If a
specification becomes available, the refusal lifts. Until then, this row is the
honest answer.

### H. Defeaturing

**H1: removing a set of faces.** **Built** as
`ogeom_bool::remove_faces(model, solid, &[Shape], tol)`. The two kinds of gap
close differently.

- *A feature whose rim is an inner loop of a surviving face* (a bore, a boss, a
  mid-face pocket) is wire surgery. The surviving face is rebuilt without the
  rim wire, and nothing is re-intersected. The drilled block comes back to its
  exact volume with **no overshoot at all**. The face-set approach fills
  nothing: the boundary is resewn, so the sliver-band problem that the
  tool-based approach paid ten microns to avoid never arises.
- *A band feature* (a fillet or chamfer along an edge) is the re-intersection
  case. The two side surfaces recover the edge the blend replaced, and the end
  faces' own edges extend along their own curves to the corners that the
  recovered edge pierces. Both the chamfered and the filleted box come back to
  the sharp box *exactly*: six faces, twelve edges, eight vertices.

One ordering lesson: the ends are rebuilt before the sides. Extending a cap's
dangling edge decides that edge for the side that shares it, and sewing rejoins
them on the one node.

Refused by name: several bands at once, bands meeting at corners (the B2/A6
family), and spline-surfaced neighbours. The parity ledger's
`bool.defeaturing` restriction is the live list, and issue #19 is its plan.

*The original entry, kept because the overshoot finding is a real measurement
of the tool-based approach that was not taken.*

The operation was to be `remove_faces(model, solid, &[Shape], tol)` in
**`ogeom-bool`**: given faces to delete, extend the neighbouring faces,
re-intersect the extensions with each other, and sew the result back into the
shell. Face extension plus re-intersection is the existing fuse machinery. What
was missing is the driver that decides which neighbours to extend, and how far.

It belongs to the boolean, not to a recognizer. The caller supplies faces; what
those faces *mean* is the caller's business, and the operation must work on a
solid whose history is gone. There used to be a `remove_feature(model, solid,
&Feature, tol)` here that dispatched on a recognized feature and rebuilt the
volume that feature described. That is a different operation with a different
input, and it left with the recognizer. Reusing its code would have kept the
half that does not generalise.

One finding from it is worth keeping, because a face-set implementation will
hit it as soon as it builds a tool that meets the solid at an opening:

> A filling tool flush with the faces it meets is a coincidence at every
> opening at once, and the boolean does not assemble it. So the tool must
> overshoot. The non-obvious part is that the overshoot cannot be small:
>
> - A margin leaves a sliver band standing past the opening. That band's
>   interior probes must be *clearly* outside the part. Otherwise the exact
>   classifier finds every ray from them grazing the face they sit against,
>   exhausts its whole fan of directions, and answers `On` the slow way.
> - A micron of overshoot is inside the band the classifier reads as "on the
>   boundary", and costs fifty seconds on a part that otherwise takes a fifth of
>   a second. Ten microns is outside it and costs nothing.
> - The working figure was a hundred thousand confusion tolerances (ten microns
>   at millimetre tolerances). The restored solid is larger than the original by
>   that amount times the openings' area.

Watch any caller whose tolerance is tighter than the overshoot.

### I. Canonical simplification

**I1: recognizing that exact geometry is secretly analytic.** **Done.**
`ogeom_heal::canonical_simplify` works as follows:

- It samples each free-form surface on its own chart, with its own normals.
- It proposes a plane, sphere, cylinder or cone using the classical
  estimators: the mean normal, the least-squares meeting point of the normal
  lines, the direction the normals avoid, and the linear taper of radius
  against height.
- It accepts only when *every* sample verifies at the caller's tolerance. The
  certificate is the worst deviation actually measured.

Free-form curves are simplified on the way too. A rim written as a B-spline
that is exactly a circle must become the circle before the analytic surface
can project it in closed form.

Results: a drum converted to NURBS comes back as a cylinder at 1.8e-15, with its
volume unchanged to the last bit; a skinned loft stays as it is. The reference
kernel's set (plane, cylinder, cone, sphere) is matched exactly. Neither side
proposes a torus.

### J. The medial axis

**J1: the medial axis of a planar region.** **Built for the convex polygonal
case**, exactly. It uses the shrinking-polygon construction:

- every edge moves inward at unit speed;
- every vertex rides its angular bisector;
- each event retires an edge and starts a branch.

Convexity is the validity condition, because it excludes the split events this
construction does not handle. Held to closed forms:

- a rectangle's four diagonals and roof line, to 1e-9;
- a 3-4-5 triangle's branches meeting at the incenter, with the inradius
  (a+b−c)/2 = 1 as the deepest clearance.

Holes, reflex corners and arcs are refused by name. The parity row's
restriction is the live worklist, if a caller ever needs the general region.

## Decisions, not gaps

These are settled. They are listed so nobody reopens them by accident.

- **A surface returns its point and derivatives in one evaluation, and states
  that the result agrees with the separate accessors only to rounding.**
  - A foot-point solve needs all six values at the same place. A
    tensor-product patch answering three separate accessors locates its spans,
    builds its basis functions and sums its control grid three times.
  - [`Surface::jet_at`] answers from one order-two table. This saves 3 to 25%
    of a spline-rich STEP read, which matters on files that take a minute.
  - The price is stated, not hidden. A patch sums its point by de Boor and its
    derivatives by basis functions. These reassociate differently, and
    `basis_derivatives` does not return identical lower-order rows at different
    requested orders. So a jet may differ from the separate accessors in the
    last ulp.
  - Consistency within a jet is what a Newton step needs, and what the type
    guarantees. A caller must not mix a jet with the accessors at the same
    parameters and expect identical bits.
- **A pcurve with no closed form is `None`, not a fit.** An exact curve with a
  fitted pcurve would have two descriptions that disagree by an amount nothing
  records. The consumer that needs one marches the pair instead.
- **A closed exact section partly outside a surface's extent is kept whole.**
  The restriction that matters is the face's trim, which is the boolean's own
  2D stage. The surface extent is only a parameterization window, and cutting
  there would split a curve where there is no boundary.
- **Scaled placements in the boolean are refused, with instructions.** A scale
  changes a surface's parameterization underneath its pcurves. Bake it first:
  `baked_shape` does exactly that, and the boolean calls it.
- **The crossing walker refuses tangential contact.** It is a crossing walker.
  The tangential walker owns that case, and the section pipeline routes to it.
  The refusal stays pinned.
- **Bi-tangent construction is subsumed**: by the 2D repertoire in 2D, and by
  the blend family's own envelope in 3D.
- **Glue is subsumed** by the boolean's same-domain unification, which already
  skips nothing it needs and unifies what glue would.

## Order

1. ~~**A**: the interference table.~~ **Done.** Five named failures closed; B2
   and D2 unblocked.
2. ~~**C**, **E2**, **F3**.~~ **Done**, each in its own section. What they left
   owed is recorded with the entries: D2's drawing side, B2's two-blend meeting
   and N>3 setback vertex, and H1's refused-by-name list. All of it is filed,
   with the kernel's other known debts, as issues #17 to #32 (features #17 to
   #22, performance #23 to #26, the tail #27 to #32). The tracker holds the
   state; this file holds the narrative.
3. ~~**The walker abstraction.**~~ **Done.** `ogeom_intersect::walk`: a
   `Condition` in `n` unknowns with `n − 1` equations, and one walk over it.
   - The missing equation is the point: the solution set of `n − 1` equations
     in `n` unknowns *is* a curve. The walker supplies the last equation itself,
     as a plane across the direction of travel.
   - The direction comes free as the Jacobian's null vector, so a condition
     does not need to know its own tangent formula.
   - Step control, stall reporting and closure are written once.
   - The intersector's own walk now goes through it, which shows it is general,
     not just present.

   One thing the abstraction had to learn. A null vector's sign is arbitrary,
   so the walker keeps its own heading. But the intersector's tangent is the
   cross product of two normals, so its sign comes from the surfaces, and its
   *flip* at a tangency is what stops the march. When the walker silently
   turned it back round, two thin curves through two touching points came back
   as one confident loop lying on neither. So a condition declares whether its
   tangent's sign is meaningful, and that declaration matters.
4. ~~**B1**: the marching blend.~~ **Done.** **B2**'s corner family is down to
   the curved-edged corner. Any convex planar vertex rounds, as the exact
   envelope of the rolling ball where no single ball touches all its faces.
   What is left depends on **A6** (tangency at a face's own corner), not on A5.
5. ~~**D1**: marched silhouettes.~~ **Done**. They were indeed a second
   condition for a walker that already existed.
6. **F2**: IGES.
7. **H1**: removing a set of faces.

**F4** has no scheduled slot: it needs a specification nobody has published.

**A debt paid.** Fifteen places in `crates/` used to point the reader to "the
deferred table", which no longer existed. All fifteen now point to
`docs/PARITY.md` rows by id, or state plainly that an earlier plan owed the
thing and it was delivered. Six of them are in error strings, so a refusal now
names a row the reader can actually open.
