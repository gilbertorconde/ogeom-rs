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

### A. The boolean's interference table

One family, and the highest-leverage item in this document: five otherwise
unrelated gaps are all the same thing, and studying how the problem is
usually structured says the fix is not five patches but one change of shape.

**The gaps.** *A1, a tool flush with the part*: refill a bore with the
cylinder it was cut by, ends flush with the faces it broke, and the kept
pieces do not close into a shell. *A2, a section through a chart pole*: a
plane cutting a ball on its own axis ends its section exactly at both poles,
where the sphere's derivatives degenerate. *A3, a ball on a box's corner
vertex*: fails in the rebuild with a geometry error rather than by name.
*A4, sections through tangential contact*: a plane through a bore's axis
meets the wall along its rulings — the textbook half-section, so this one
buys a drawing feature too (D2). *A5, a shell around a three-cylinder tip*:
blocks the vertex blend (B2); the tool is known, the removal fails.

**Why they are one gap.** Our boolean computes face-against-face sections
and then decides where each resulting piece stands by *asking* — probing an
interior point and casting rays. Every failure above is that question being
put somewhere it cannot be answered: on a coincidence, at a degeneracy, in a
cusp, on a tangency. The multi-probe retry and the widened bands are
workarounds for asking a question that should never have been asked.

The established shape for this is an **interference table**: a pass that
computes, once and up front, every way the two shapes meet — and a build
phase that only *reads* it. Four properties of that shape matter here.

*Interferences are computed in dimension order.* Vertex against vertex,
then vertex against edge, edge against edge, vertex against face, edge
against face, and only then face against face — each stage consuming what
the earlier ones established. A3 is a vertex-against-face interference, a
first-class thing computed early; today it has nowhere to be recorded and
surfaces as an anomaly during face splitting.

*Coincidence is identified, not rediscovered.* Where two edges lie on each
other, or an edge lies in a face, they are recorded as one shared piece
carrying the set of faces it lies on, and a substitution map says which
entity stands in for which. A1 is exactly this: the flush tool's faces
coincide with the part's, and today nothing says so until the classifier
trips over it.

*An edge's interference range excludes its own ends.* A piece of an edge is
bounded by two vertices, each with a tolerance; the part of it that can
genuinely interfere with anything is the stretch *outside* both tolerance
spheres. Testing on that stretch rather than the nominal one is the
systematic version of what our probe-ranking does by luck.

*Each face carries a state cache.* During the interference pass, every face
accumulates which boundary pieces and vertices are inside it, which are on
it, and which came from section curves. The build phase reads that cache
instead of probing: a piece bounded by pieces known to be *on* the other
face is on it, and no ray is cast at all.

*Degenerate edges get their own late pass.* After pcurves exist, each
degenerate edge is split by finding the boundary pieces that run through its
vertex and pairing against them **in the chart**, where the pole is a line
rather than a point. That is what we already do for pole strands; what we
lack is the state cache that would let the pieces either side be classified
without a probe. A2 is that.

**The construction.** An `Interference` table built before any splitting,
holding: vertex identifications, edge pieces with their interference ranges,
shared pieces with their substitution map, the per-face state cache, and the
section curves. Then a build phase with no classifier in it. This is a
larger change than "fix the arrangement", and it is worth it: it removes an
entire class of failure rather than another instance of it.

**A caution worth recording.** Mature implementations of this design still
state that a face's classification can depend on which point of the face is
chosen, and still keep an open-ended list of configurations that defeat
them. The table does not make the problem disappear; it makes most of the
cases stop needing the fragile question. Expect the campaign to close A1–A5
and to leave a shorter, better-named list behind.

**Verification.** Each of A1–A5 becomes a test written before the change,
asserting the measured result — volume, face count, shell closure — not
merely that no error came back.

### B. Blending

**B1 — the marching blend.** *The plan here has changed.* The previous
entry proposed intersecting the two supports' offset surfaces to get a
spine, projecting back onto each support for the tangency points, and
skinning the arcs between them. That works on paper and is the wrong shape:
the tangency curves arrive by projection, so the legs' pcurves are *fitted*,
and fitted pcurves on a support are what the boolean cannot treat as
same-domain.

The formulation that avoids this solves for the section's two endpoints
directly. **Unknowns: four** — `(u₁, v₁)` on the first support and
`(u₂, v₂)` on the second, the two points where the rolling ball touches.
**Equations: four.** Three say the ball's centre is the same point computed
from either side,

> `P₁ + r·n₁ = P₂ + r·n₂`

where `Pᵢ` and `nᵢ` are the point and unit normal of support *i* at its
parameters. The fourth ties the section to a guide — the edge being blended,
or any curve running along the seat — by requiring the section to lie in the
plane through the guide point normal to the guide's tangent.

March the guide parameter, solving the four-by-four system at each step.
What comes out is worth the change: the tangency curves emerge **in the
supports' own parameters**, so the legs' pcurves are exact by construction
rather than fitted, and the blend's own surface is the skinned sections as
before. This is why the reformulation is not a detail — it is what makes the
result something the boolean can consume.

*A second solver, for the boundaries.* When a section endpoint runs off the
edge of its support, the answer is not to clip afterwards but to switch to an
inverted system: unknowns `(t, w, u, v)` where `t` runs along the support's
own boundary curve, `w` along the guide, and `(u, v)` are the partner's
parameters. Solving that finds exactly where the blend crosses the support's
boundary, so the blend is trimmed on the boundary rather than near it.

*The states to implement,* which are a case checklist obtained for free:
step too large, step too small, the march reversed, the two section
endpoints collapsed onto each other (the radius is too large for the local
geometry), the section reached the boundary of the first support, of the
second, or of both — and *unhooked*, where the blend leaves a support's
boundary and must continue onto the adjacent face.

*Where to start:* a cylinder meeting a plane at an angle. One leg is planar
and the other's tangency curve lies on a cylinder, both of which our pcurve
machinery already handles exactly, so the first seat exercises the solver
without also exercising the fitting.

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
three fillet cylinders and was cut away when that edge was rounded. What
fails is the removal (A5), so this is blocked rather than unknown.

### C. Sweeping

**C1 — evolved shapes.** A profile wire swept along a spine that may be a
wire *or a planar face* — the face case is the one we had not accounted for,
and it is what makes the operation a volume rather than a shell. The
construction is a composition of what §11 already has: the spine's straight
edges sweep the profile as prisms, its arcs sweep it as revolutions about
their own axes, and the convex corners between them take the profile
revolved about the corner — the same join the 2D offset already decides. The
pieces exist; the assembly and the corner bookkeeping do not.

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

**D2 — the on-axis half-section.** Blocked on A4.

### E. Interaction and documents

**E1 — a multi-resolution pick hierarchy.** One structure serving several
deflections instead of independent scenes each built at one. A performance
shape: every answer picking gives today is already exact, including the
refinement onto the true surface.

**E2 — datum targets and presentation PMI.** The annotation polylines and
planes, and the targets a datum feature carries. Data plus their STEP
carriage; the semantic PMI core they attach to is already read and written
both ways.

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

**F3 — reading glTF.** Written as GLB already. Reading means the accessor and
buffer-view indirection, every component type and stride a writer might have
chosen, and the sparse-accessor form. Ordinary work, not yet done.

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

1. **A** — the interference table. It unblocks B2 and D2, closes five named
   failures, and is the difference between a boolean that is usually right
   and one that is dependably right. It is also the largest structural
   change left, so it goes first while there is appetite for it.
2. **C**, **E1**, **E2**, **F3** — contained pieces, any order, each a stone.
   Good work to interleave when the appetite for A runs out.
3. **G1** — the tessellation instrument, which unblocks measuring a class of
   rebuilt shape.
4. **The walker abstraction** — generalize the marcher over the condition it
   follows. Small on its own, and it is the shared foundation of B1 and D1,
   so it comes before either.
5. **B1** — the marching blend, one seat at a time, starting with a cylinder
   meeting a plane at an angle. Then **B2**'s corner family, on top of A5.
6. **D1** — marched silhouettes, which by then are a second condition for a
   walker that already exists.
7. **F2** — IGES.

**G2** and **F4** have no scheduled slot: the first needs an idea nobody has
had yet, the second needs a document nobody has published.
