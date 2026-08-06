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

### A. Degenerate configurations in the boolean

One family, and the highest-leverage item in this document: several
otherwise-unrelated gaps are all the arrangement failing where two shapes
touch instead of crossing. Fixing the family closes five entries.

**A1 — a tool flush with the part.** Refill a bore with the cylinder it was
cut by, ends flush with the faces it broke, and the kept pieces do not close
into a shell. That is the reproduction; it is exact and it is fast to run.

Standing the tool a little past those faces works, and *how far* matters in a
way worth keeping written down. The sliver band an overshoot leaves has its
own interior probes, and if those sit inside the band the exact classifier
reads as on-the-boundary (`confusion * 1e4`), every ray from them grazes the
face beneath, the classifier exhausts its whole fan of directions, and one
fuse takes fifty seconds where it should take a fifth of one. Ten microns of
overshoot is decisively outside that band and costs nothing. `remove_feature`
overshoots by `confusion * 1e5` for exactly this reason.

**A2 — a section running through a chart pole.** A plane cutting a ball on
its own axis ends its section curve exactly at both poles, where the sphere's
derivatives degenerate. The pole strands are paved where sections land on
them; the pieces still do not close.

**A3 — a ball seated on a box's corner vertex.** Fails in the rebuild with a
geometry error rather than by name, which is the worse of the two failures.

**A4 — sections through tangential contact.** A section plane exactly through
a bore's axis meets the wall along its rulings. The textbook half-section is
this cut, so this one buys a drawing feature too (see D2).

**A5 — a shell around a three-cylinder tip.** Blocks the vertex blend (B2).
The tool is known and simple; it is the removal that fails.

*Construction.* These want one campaign, not five patches: a pass over the
arrangement's handling of coincident and tangential strands, with each
reproduction above as a test written first. The classifier's on-boundary band
and its direction fan are the two mechanisms implicated in all of them.

*Verification.* Each reproduction becomes a test asserting the measured
result — volume, face count, shell closure — not merely that no error was
returned.

### B. Blending

**B1 — the marching blend.** Edges whose rolling-ball envelope has no closed
form. The spine is the intersection of the two supports' offset surfaces,
which the intersector already computes; the sections are circular arcs
between the two tangency points, which are the spine point's projections onto
each support; the blend face is those arcs skinned, which the grid fit
already does. The legs must be built on the *supports' own surfaces* with
fitted pcurves, or the boolean sees two coincident-but-distinct surfaces
rather than a same-domain pair.

This is the largest single item in this document. For scale: a mature
implementation of this family runs to tens of thousands of lines. It is a
milestone, not a stone, and it should be approached one seat at a time — a
cylinder meeting a plane at an angle is the first, since one leg is planar
and the other's tangency curve projects onto a cylinder, both of which the
pcurve machinery already handles.

**B2 — the setback vertex blend.** The corner patch where three fillets meet.
The tool is known: after the three edge fillets, the corner block running
from the ball's centre out past the corner, *minus the ball*, is exactly the
leftover spike — anything in that block further from the centre than the ball
lies outside one of the three fillet cylinders and was already cut away when
that edge was rounded. What fails is the removal (A5). Blocked, not unknown.

### C. Sweeping

**C1 — evolved shapes.** A planar spine's straight edges sweep the profile as
prisms, its arcs sweep it as revolutions about their own axes, and the convex
corners between them are the profile revolved about the corner — the same
join the 2D offset already decides. The pieces exist; the assembly and the
corner bookkeeping do not.

### D. Drawings

**D1 — marched silhouettes.** A torus or a spline silhouette has no closed
form and needs the same marching the intersector does. The exact path refuses
these by name today and the polygonal path draws them, which is a working
answer at the mesh's accuracy rather than a missing one.

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

1. **A** — the degenerate-configuration campaign. It unblocks B2 and D2, and
   it is the difference between a boolean that is usually right and one that
   is dependably right.
2. **C**, **E1**, **E2**, **F3** — contained pieces, any order, each a stone.
3. **G1** — the tessellation instrument, which unblocks measuring a class of
   rebuilt shape.
4. **F2** — IGES.
5. **B1** — the marching blend, then **B2** on top of A5.
6. **D1** — marched silhouettes, which share machinery with B1.

**G2** and **F4** have no scheduled slot: the first needs an idea nobody has
had yet, the second needs a document nobody has published.
