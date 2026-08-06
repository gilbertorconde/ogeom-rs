# Plan

What is left to build, why each piece is not built yet, and how each one gets
done. This replaces the scope map, which had become a record of finished work
with the unfinished part scattered through it.

The kernel is substantially complete: topology and the geometry vocabulary,
the intersector, the boolean, healing, blending, offsets and sweeps,
tessellation, drawings, sketching, selection, recognition, the document
layer, and STEP, STL, DXF, glTF, OBJ, PLY, VRML and 3MF exchange. What
follows is the remainder, and it is small enough to enumerate exactly.

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
- **Scope is set by the capability, not by any application's demand for it.**
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

**B2 — corners where blends meet.** *Reframed.* This is not one construction
but a family, classified by **how many blends meet at the corner**: one, two,
three, or more. Our owed "setback vertex blend" is the three-blend case; the
one- and two-blend cases are separate constructions we had not enumerated at
all — a blend running out at a face boundary, and two blends meeting, which
is resolved either by intersecting them with each other or by extending both
to a shared end.

For the three-blend equal-radius corner the tool is already known: after the
three edge fillets, the corner block running from the ball's centre out past
the corner, *minus the ball*, is exactly the leftover spike — anything in
that block further from the centre than the ball lies outside one of the
three fillet cylinders and was cut away when that edge was rounded. The
removal now works (A5 is closed), so this is ordinary work rather than
blocked.

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

**D1 — marched silhouettes.** A torus or a spline silhouette has no closed
form. The useful observation is that it needs no new machinery either: the
silhouette is the zero set of

> `n(u, v) · d = 0`

for view or light direction `d`, which is an implicit condition on the
surface's own chart — exactly the kind of thing our marcher already follows.
And where the contour meets the face's trim, the boundary crossing is found
by a separate solve against the boundary curve, the same shape as the
blend's inverted system in B1.

That is a cross-cutting note worth stating once: **one walker, several
conditions.** Surface intersection tracks "on both surfaces"; a silhouette
tracks "normal perpendicular to the view"; a blend tracks the four-equation
contact system. If the marcher is abstracted over the condition it follows —
value, derivatives, and a boundary-crossing solve — then B1 and D1 share it
and the intersector's hard-won step control, stall detection and branch
handling are inherited rather than rewritten twice.

**D2 — the on-axis half-section.** Unblocked: A4 is closed, and the boolean
now cuts a bore on its own axis and reports the two rulings. What is left is
the drawing side — the half-section as a view, with its own hatching
convention.

### E. Interaction and documents

**E1 — a multi-resolution pick hierarchy.** **Done.** `PickHierarchy` builds
one scene per deflection over one shape, coarsest first, and answers a pick by
descending: each coarse level rules faces out, and only what survives is put
to the finest one. Face indices are the same at every level by construction,
which is what lets an answer travel between them.

The property that matters is that *nothing changes*. Each level's mesh stands
within its own stated chord of the true surface, so widening the coarse
level's boxes by the two chords cannot rule out a face the fine level would
hit — and the test does not argue this, it checks it: every face the fine
level strikes is admitted by the coarse one, ray by ray, and the descent's
hits equal the finest level's own, sub-shape, kind, depth and refined
position. What is measured about the *performance* is a count, not a timing:
the coarse level admits well under half the scene.

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

**F2 — IGES, both directions.** The record layout is a published standard and
the entity-to-topology mapping is documented. A week-scale job, staged the
way STEP was: geometry first, then trimmed surfaces, then assemblies.

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

**F4 — SAT, X\_T, JT.** Refused for want of public documentation. These are
proprietary formats whose specifications are not published; implementing them
would mean reverse-engineering files rather than reading a standard. If a
specification becomes available the refusal lifts, and until then the honest
answer is this row.

### G. Healing residuals

**G1 — tessellating a rebuilt full-period face.** A mesh-to-B-rep drum comes
back as two caps and one recognized cylinder wall at its exact radius, closed
as a shell — but the wall spans a full period, and a full-period face does
not tessellate closed, so the result's volume cannot be measured. The
topology is right; the instrument is what fails.

**G2 — segmentation across tangent-smooth junctions.** Where a fillet band
meets its wall the two share a tolerance-ambiguous strip, and the
recognition-driven segmentation cannot place the boundary from curvature
alone.

## Decisions, not gaps

These are settled. They are here so nobody reopens them by accident.

- **Analytic per-constraint Jacobians in the sketch — declined.** Central
  differences over the parameters each constraint names already agree with
  the analytic value to a part in ten billion, at one extra residual
  evaluation apiece. Taking them exactly means rewriting every residual over
  a scalar trait: one more place for a residual and its derivative to drift
  apart, bought for speed the solver does not need at sketch scale. If
  profiling ever says otherwise, the change is mechanical and this is where
  to start.
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
2. **C**, **E1**, **E2**, **F3** — contained pieces, any order, each a stone.
3. **G1** — the tessellation instrument, which unblocks measuring a class of
   rebuilt shape.
4. ~~**The walker abstraction.**~~ **Done.** `ogeom_intersect::walk`: a
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
5. **B1** — the marching blend, one seat at a time, starting with a cylinder
   meeting a plane at an angle. Then **B2**'s corner family, on top of A5.
6. **D1** — marched silhouettes, which by then are a second condition for a
   walker that already exists.
7. **F2** — IGES.

**G2** and **F4** have no scheduled slot: the first needs an idea nobody has
had yet, the second needs a document nobody has published.
