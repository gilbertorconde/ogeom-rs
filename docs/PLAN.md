# Plan

What is left to build, and the rules the work follows.

## Where the remaining work is tracked

- **GitHub issues** hold bugs and requests, one per issue.
- **The parity ledger** (`docs/PARITY.md`, generated from
  `docs/parity/parity.toml`) audits the reference kernel's modelling modules
  against named capabilities. The restriction on each `partial` row is the
  list of what that capability does not do yet. `tools/check.sh` holds the
  ledger to its evidence.
- **This file** holds the open items that are neither: known limits of the
  code, refusals by name, design items and speed work.

Finished work is recorded in `CHANGELOG.md` and in the tests that pin it, not
here.

## Open

Within each group the items are ordered by risk: a possible wrong answer
first, then refusals, then cost.

### Checking and measurement

**Slits between a fitted edge and its pcurves.** The exact integral takes the
strip between a lifted pcurve and its edge's curve as a ruled surface where
the curve has a closed form and the pcurve stands off it by more than a
hundred confusions. A fitted curve's pcurves within its stated tolerance are
taken as lying along it. A converted face bounded by fitted edges stating
2e-2 to 5e-2 then leaves slits: converted tori measure 0.2% to 0.7% under a
fine mesh, and the volume moves with the reference point. The torus commons
of the stress run carry up to 1.5e-9 relative for the same reason. Taking a
strip on every edge that stands off at all closes both (the commons within
1e-13), but makes the stress run 5.6 times slower, so it is not taken.

**Fits measured only at their samples.** A surface fit (`to_nurbs_within`, a
patch's degree restriction) states the worst miss at its samples and span
quarter points, which can sit a little under the true miss; a fit whose
surface cannot be evaluated somewhere (an offset sphere's poles) refuses. In
the boolean, `projected_into_shared_chart` ignores the owner fit's `met` (the
contact piece states the dense gap after). In `ogeom-heal`, the joint fit in
`split.rs` and the iso fit in `divide.rs` are measured at their samples only.
A blended pipe section fitted to a tenth of the sweep's tolerance stands 0.8%
past that tenth between the eighths it is checked at.

**The mesher tells a slit from a hole by width alone.** An inner ring whose
mean width (twice its area over its perimeter) is under a micron is a slit
(`is_slit` in `ogeom-mesh/src/triangulate.rs`). Boolean slits measure about
six tenths of a micron, real holes from a micron and a third, so a wall poked
through more thinly loses its hole. Telling them apart by what fills the ring
would hold at any width.

### Booleans

**Nearly coplanar faces are not paired as coincident.** A piece that reads on
the other solid's boundary at every point asked, with no coincident partner
face, is refused by the boolean's piece classification ("a contact over a
region of faces not recognised as coincident is refused"). The case is a
faceted wall leaning off a pad's wall by under a ten-thousandth of a radian:
the band where the two stand within the weld distance covers most of the
facet. A slab drafted 0.00001 under a pad refuses about a third of its
row-count sets; where a corner facet meets the wall only at a vertex, the
planes' solved line crosses it diagonally a tenth of a millimetre off and the
kept pieces do not close. The design: pair the band as same-domain and split
it where the planes part by the weld distance. At drafts of 0.00005 and
0.00002 the remaining refusals are a lower corner row's facets crossing the
pad's walls, or a ray meeting a tangency.

**Section loops pinned at a sphere's pole.** Where a curved wall passes
through a pole of a sphere face, the marched section stalls and doubles back
at the pole, and the arrangement can read one piece round both loops. Before
the boolean, `ball_chart.rs` restates a whole ball or a trimmed sphere face on
a chart whose poles stand clear. It leaves as it is (and the pole stays): a
face whose trim the new chart would have to split (two loops round the axis,
as a band between two planes has, or a loop the seam meridian crosses at no
vertex), a face at a placement or on a left-handed sphere, and a face with an
edge that has no curve in space.

**A section whose fit misses widens its junctions as far.** A marched section
keeps its fit's measured error as its tolerance, and the paving's reach and
the junctions along an edge grow with it, unbounded. Two sources remain. A
section through or beside a sphere's pole fits its samples, but its sphere
pcurve strays up to 0.09 between them; split at the pole sample each piece
measures under 1e-3, but the boolean then refuses a drill it passes whole, so
the split waits on the pole item above. A branch nothing fits (residuals of
0.1 to 10 near a tangency), and a near-tangent crossing whose capped angle
term stands at 1e-3 or more, carry their miss.

**A drill lying in a face that runs into a round.** A drill whose lowest line
lies in a plate's underside meets the round at each end within a few microns
of tangency. The section is a figure eight crossing the round's tangent edge
twice at its double point, or two loops crossing it a hundredth of a
millimetre apart. With the drill's seam toward -y, the marcher returns the
section as six open pieces whose fits miss by up to 9e-3, and the kept pieces
do not close; the same placement moved by a few ulps leaves the drill less
the part unclosed. With the seam toward z it closes.

**A drill's seam along the other part's edge.** A drill whose wall runs along
a round's tangent line with its own seam on that line: the section hugs an
edge on both sides and lays no split, and the boolean refuses. A drill tangent
to a plane along a face's edge, seam there, is the same. Carrying the seam
onto that face as a contact closes it, but the same contact breaks a drill
whose seam crosses a face's interior lying in its plane. Turned away, the
seam closes.

**A section tangent to an edge at its end.** A drill whose circle on a rounded
box's base touches the round's tangent line at the corner vertex, where the
corner ball touches the base, is refused. At seeds 1 to 3 the stress drill
scenario refuses one, six and one drills (on the rounded box, a converted
plate, the frustum, two corpus parts and the shaved cube); the three items
above cover the diagnosed ones.

**A tilted half ball drilled on the diagonal chart.** A half ball charted
about the cube diagonal (1, 1, 1) with its seam toward z, on a drum and drilled
through, gives fuse + common and cut + common 1.2e-8 relative off the
operands, where other charts hold 5e-9. The rim and section are fitted on the
tilted chart (held to 7e-6 and 2e-5), and the measured volumes move by that
much.

### Fillets and blends

**One corner, two results by order.** A cube shaved by a drum: a bottom edge
and the upright edge at its end, where the drum takes over from a side face,
rounded bottom first, take off the union of the two blends, the bottom's flush
cap leaving a sliver standing against the drum. Rounded upright first, the
bottom blend runs on through the upright's band and trims that sliver too
(0.8% more at radius 0.18, 4% at 1). Both are valid solids; one corner should
give one result. `edge_chains.rs` pins the first order.

**Curved face blends refused by name.** `blend_faces` and `fillet_faces`
round curved faces that share no edge where their surfaces share a direction
or an axis, and otherwise march the ball round where the surfaces cross. They
refuse (`curved_face_blends_refuse_by_name` in `tests/blends.rs`):

- faces the ball leaves at different places round its seat and never both at
  once, so neither face gives the round both its ends;
- a marched round closing on itself between separate faces, and one whose
  line of contact crosses a face more than once;
- a round crossing a face that stands between the two (a ball crossing a
  chamfer between the faces it rounds), where the solid beside the round is
  open where the seat says material, or the reverse.

**A pole leg whose edge winds round the pole and whose rail does not.** A
marched band whose rail passes over a sphere's pole has a leg from the
winding loop to the pole with the other loop cut out as a hole. When the
winding loop is the rail, `pole_leg` in `marched.rs` splits it clear of the
hole and starts the seam there. When it is the edge, the seam leaves from the
edge's own vertex and can cross the hole. No measured part builds one.

### Sweeps and drafts

**A draft across a seam closed only to position.** `general_draft` on a closed
face whose seam is closed to position but not to tangent (a skinned loft's
wall): the normal turns across the seam and the exact rulings with it, so the
drafted wall closes on the mean ruling and blends to each side's within a
thirty-second of the hinge. There it stands off the exact draft by half the
turn times the reach (6.5e-4 on a loft of circles of radius 10 skinned at
1e-3, reach 9); elsewhere within 1e-4. Such a draft rebuilds in about 1.3 s,
the time in the exact volume over the wall.

### Construction and exchange

**A planar face whose ring winds against its plane.** `make_face` given a
plane and a ring walked clockwise about the plane's normal builds a face that
is ambiguous. Sewing trusts the ring and turns it with its neighbours, so a
ruled sheet closed by such a cap can come back with the cap facing in; where
the cap has no pcurves the shell cannot be meshed to settle it. `make_face`
could refuse such a ring, or turn the plane to agree with it.

**STEP pcurves the file cannot state exactly.** The reader keeps a file's
pcurve where its lifted gap is within twice the projection's. A pcurve the
file cannot state in the curve's parameter (a reversed or left-handed conic,
a trig or offset trace, a curve written as a spline conversion) is derived on
reading. A foreign file's pcurves are read without its length or angle unit,
so one in inches or degrees fails the gap test and is derived.

**Edges through a pole on reading.** Both readers cut an edge through a pole
of a face it bounds (a sphere's pole, a cone's apex, a patch's collapsed side)
before any face is built. Not cut: IGES trimmed surfaces (144), whose faces
sew whole edge to whole edge, and pcurves other producers build through a
pole. The cut edge stays in the model unused, and a STEP style or shape
aspect naming its file id lands on every piece.

### Mesh conversion

`solid_from_mesh` rebuilds planes, the four canonical surfaces, extrusions,
surfaces of revolution and fitted B-spline patches, and leaves every other
region faceted; `MeshSolidReport::fallbacks` names each curved region it
faceted and why. Each change is held to the converter's standard (verified
against every sample, edges placed on both surfaces, facets as the fallback)
and measured on `mesh_corpus` and the truth bench before and after.

**Seams between fitted regions.** Each region is fitted on its own and each
seam solved point by point along the mesh boundary, or threaded as a chord
through the boundary's vertices. Folds and faces that mesh open are caught
after the build and their curved faces withdrawn to facets, so every
measured part tessellates closed at the cost of those faces. A chord where no
point carries onto both surfaces (a torus and a cylinder fitted to pieces of
one spline blend) can stand 0.2 off, stated as its tolerance. Design items,
in order: absorb single facets at junctions without costing an exact seam
(taking every facet within reach breaks exact rims); trace seams between
solved corners as one branch of the intersection guided by the mesh path;
junction rules (Benière et al. 2012: a junction stays only if its faces are
pairwise adjacent, dangling edges dropped).

**Fillets fitted apart from their supports.** Derived from their supports:
rounds between two planes, a torus between a plane and a coaxial cylinder or
cone, corner balls, and a sphere put on the torus a round turns on beside a
plane. Not derived: a torus between two cylinders or a cylinder and a sphere,
variable radius, and fillets between curved faces off a common axis, which
keep their fitted surfaces and meet their supports at a near tangency. A
torus corner a few facets round between two fillet cylinders and a floor fits
nothing (every sample reaches into the cylinders) and stays facets.

**Sweeps.** A whole turn is recognized, its axis settled on the samples. Not
done: the band layout for a sweep between two circles (its seam the profile
itself), and a closed profile (an extrusion of a closed curve, a revolution
of one). The helical line complex would tell a thread from a surface of
revolution; its threshold has to sit above the twist a mesh's own normals
carry (a pitch of about 1e-3 on a turned wave meshed at 0.01). A thread goes
to the patch.

**Patches.** A smooth region nothing else fits becomes a fitted B-spline
patch where it is one disk with one loop and verifies
(`MeshSolidOptions::patches`). Limits:

- a tangent line inside a row of triangles (a scan, or a mesh without
  vertices along its edges): the curvature's jump falls inside one span, the
  fit runs out at 7e-5 to 4e-4 against 1.5e-5, and a near-flat foot gathers
  into narrow slanted planes;
- a run-out steeper than the fit can follow stays faceted (a hill on a round
  bar at 0.8 high, the fit at 1.24e-5 against 1.23e-5), and a mesh too coarse
  to leave a hill enough vertices keeps the sphere round the hill's facets;
- run-outs into cones and tori are threaded as into cylinders and spheres,
  but not measured;
- between vertices a fit is no closer to the true surface than its knots
  allow: on coarse meshes up to 4 times the distance in the row at a free
  edge, under 1 inside.

**A torus pierced with no circle free either way.** Holes in a whole torus
are inner wires, its seams placed on a free parallel and meridian or threaded
through the holes. Holes that leave no parallel and no meridian free, or that
overlap along the chain at every level tried, stay facets.

**Region edits.** `MeshRegions` `fit` has no sweep kind, and the tangent
passes do not revisit a neighbour an edit changes (a ball beside a merged
round keeps its fit). The corner and boundary graph of the seam work above is
what an application would edit next.

**Inflections on a coarse open mesh.** Where the turns change sign, the
coplanar distance estimate asks them to change linearly, so a coarse mesh of
a surface whose curvature varies fast can read its inflections as scatter.

### Speed

**A per-face state cache.** The boolean's build phase could read a piece's
classification instead of probing for it. No case needs it for correctness.

**The rest of the local boolean.** A face the tool leaves alone is neither
split nor arranged, and the sew compares only rebuilt pieces and the copies
beside the other solid (`sew_around`); `boolean_local` in `tools/ogeom-bench`
times it. Still whole-solid: gathering both solids (each edge sampled at 33
points for its box), rebuilding every kept face as new nodes, the history
over all of them, the seam join and closure passes, and the classifiers'
face bounds; a crossed face's interior probes are sought over all its holes.
Profiled on a 900-hole plate, no phase is above a fifth, so no one change
wins much. A drill into a face of 900 holes spends 45% in that face's
arrangement (the scanline probe search over 200,000 outline points, each hole
walked on its own).

**Sharing unchanged faces' nodes.** `History::copy` records an exact copy on
new nodes. Sharing the nodes would remove up to a third of the rebuild, and
needs: the rebuilt neighbours to name the shared edges and vertices (pcurves
on the input's own surface ids), no weld or tolerance widening on a shared
vertex (`Rebuild::vertex` and `Model::add_face` widen in place), the sew to
leave a shared face's derivation alone, and every operation that edits its
result in place (heal, fillet, offset, sew) to copy on write first.

**Fillets in fewer passes.** A chain is taken in rounds, no two edges of a
round sharing a vertex, each round in one boolean each way. Each edge still
builds its blend against what the earlier rounds left.

## How this project works

These rules are enforced. Every one of them has caused a patch to be rejected.

- **Measured, not asserted.** A claim about geometry is backed by a test that
  measures it against a closed form, a published value, or an independent
  computation.
- **An approximation must say so.** Sampling, fitting and tolerances are fine
  when they are stated. A number that hides its own error is not.
- **Refusals are by name, and pinned.** When the kernel cannot do something,
  the error text says which thing and why, and a test holds it to that. A
  silent wrong answer is the only unacceptable outcome.
- **No stubs.** A thing is implemented, or it is refused by name and tracked
  as above.
- **Scope is parity with the reference kernel's modelling modules, plus what
  `docs/SCOPE.md` admits.** `docs/SCOPE.md` is normative and says how to
  decide a case. Parity is about capability, not structure.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against,
  bundles or imports another kernel, and the design is worked out here, not
  mirrored. The field's vocabulary is used throughout.

## Decisions, not gaps

These are settled, and listed so nobody reopens them by accident.

- **A surface returns its point and derivatives in one evaluation, agreeing
  with the separate accessors only to rounding.** A foot-point solve needs
  all six values at one place, and a patch answering three accessors locates
  its spans, builds its basis and sums its grid three times.
  [`Surface::jet_at`] answers from one order-two table (3 to 25% of a
  spline-rich STEP read). A patch sums its point by de Boor and its
  derivatives by basis functions, which can differ in the last ulp.
  Consistency within a jet is what a Newton step needs and what the type
  guarantees; a caller must not mix a jet with the accessors at the same
  parameters and expect identical bits.
- **A pcurve with no closed form is `None`, not a fit.** An exact curve with a
  fitted pcurve would have two descriptions that disagree by an amount nothing
  records. The consumer that needs one marches the pair instead.
- **A closed exact section partly outside a surface's extent is kept whole.**
  The restriction that matters is the face's trim, the boolean's own 2D
  stage; the surface extent is only a parameterization window.
- **Scaled placements in the boolean are baked first.** A scale changes a
  surface's parameterization under its pcurves; `baked_shape` rebuilds the
  geometry, and the boolean calls it.
- **The crossing walker refuses tangential contact.** The tangential walker
  owns that case, and the section pipeline routes to it. The refusal stays
  pinned.
- **Bi-tangent construction is subsumed**: by the 2D repertoire in 2D, and by
  the blend family's own envelope in 3D.
- **Glue is subsumed** by the boolean's same-domain unification, which skips
  nothing it needs and unifies what glue would.
- **A reversed face walks its wires backward, and its outer wire is stored
  first.** `ordered_children_of` reverses a reversed parent's children, which
  a wire's walk needs and which lists a reversed face's holes before its
  outer wire. Code that wants the outer wire reads `Model::outer_wire`; code
  that wants every wire under the face's sense, outer first, reads
  `children_of`. Most readers pick the outer ring by area.
- **A fillet refuses an edge between tangent faces.** There is no corner to
  round ("the edge's faces are tangent"). The stress fillet scenario's
  refusals are all of this kind or an edge splitting two coplanar faces, and
  its crease list takes both.
- **A marched blend's run-out is capped in the ball's section.** Where a round
  runs off a face, the cap stands in the ball's section through the crease's
  end, which leans with the faces; a cap in a side's plane would end the band
  off a section.
- **Face bounds are not cached in the model.** `shape_bounds` and
  `tight_bounds` take about a millisecond on a four-hundred-face part, which a
  cache and its invalidation would not improve on.
- **The converter uses no learned model.** Published learned reconstruction
  (2024 to 2026) gives a valid solid for 70 to 76% of parts under a hundred
  faces, every face a B-spline, and places its own bottleneck in rebuilding
  the b-rep, which is the converter's verified build. A learned proposal is
  also a result nothing here measures. Global selection by an integer program
  is out for the same reason: its time has no bound.
- **SAT, X\_T and JT are not read or written.** Their specifications are
  unpublished, and implementing them would mean reverse engineering files
  instead of reading a standard. If a specification is published, this lifts.
