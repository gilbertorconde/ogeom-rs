# Changelog

The fifteen library crates are one kernel and move together: they share a
version, and an entry here covers all of them. `tools/` is not published and
is not recorded here.

Dates are the release date. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/), which before 1.0 means a minor
bump may break the API and a patch bump may not.

## [Unreleased]

### Added

- **A closed wire turns sharp corners run by run.** (#394)
  - Extended and round corners, and any corner under the auxiliary,
    binormal or evolving laws, swept only faces. A closed planar wire now
    sweeps there as the walls of the face it bounds, the end caps taken
    off. An open wire is still refused there by name.
  - Measured: a 2 by 2 square wire up and over a right-angled L of 10 by
    10 legs has walls of area 160 extended or mitred and
    `160 - 6 + 3π/2` rounded, to 1e-12, under the rotation-minimizing
    frame and the binormal.

- **IGES external references resolve against the files beside the
  deck.** (#404)
  - `read_iges_file` reads a deck from disk and resolves its external
    references (416) in the same directory. A file is looked for under
    the last component of the identifier the reference gives, and must
    name itself so in its own global section (IGES 5.3, section 3.6.4).
  - A reference to a definition (form 0) takes the entity the referenced
    file's index (402, form 12) lists under its symbolic name. A
    reference to a whole file (form 1) takes every shape that file reads
    as. Either is built at that file's own coordinates and units, for the
    instance naming the reference to place. Files may refer on to files.
  - Refused by name, with the instance placing it: a missing file, a
    file that names itself otherwise, a symbol the index does not list, a
    logical reference (form 2), a definition kept natively by the
    receiving system (forms 3 and 4), and a file that instances itself.
    `read_iges`, handed text, has no directory to look in and refuses
    every reference so.
  - Measured on synthetic decks: an inch library's unit patch placed at
    (100, 0, 0) reads area 645.16 about (112.7, 12.7, 0), and a whole
    file's 2-cube and point placed 50 up read volume 8 about (1, 1, 51)
    and the point at (3, 4, 55).
- **IGES reads and writes stand-alone points.** (#404)
  - An independent point (116) that is geometry was counted as skipped.
    It now reads as a vertex in `IgesImport::points`, carried by its
    matrix, and a point a subfigure definition holds is placed with it.
    A point that is annotation, or that another entity owns, is not one.
  - The writer refused a part holding a free vertex. It now writes it as
    an independent point. STEP still refuses it.
  - Measured on a synthetic deck: the three points stand exactly where
    the file, its matrix and its instance put them, and again after a
    round trip.
- **`fix_shape` strips tails, marks wireless faces natural, and chains
  free edges into contours.** (#397)
  - A wire that ran out along an edge and straight back along the same
    curve kept the pair: a vertex bounding no corner, which every boolean
    through the face meets as a strand and the mesher strips on its own.
    The pair now goes where both edges lie within the wider of their
    tolerances of each other, the chart walks back to the same copy (two
    meridians a period apart bound the whole sphere and stay), neither is
    a seam, and no other face walks either. A hole made of such pairs
    alone goes; an outer wire is never stripped to nothing.
    `FixReport::tails_stripped` counts the pairs.
  - A face with no wires that had lost its natural-restriction mark is
    marked again on a sphere, a torus or a patch, whose domain is the
    whole surface; on a plane or another windowed surface it stays
    reported. `FixReport::natural_bounds_marked` counts them.
  - `ogeom_heal::free_bounds` chains a shape's free edges head to tail
    into closed contours and open chains, and `FixReport::free_bounds`
    holds them for the shape as mended. Nothing is sewn across a gap
    wider than the edges state: on DoomCube's 24 open sheets, sewing at a
    millimetre joins hundreds of edge pairs that are not twins.
  - `ogeom_algo::edges_apart` measures how far two edges' curves stand
    apart.
  - Measured: the corpus spike (`a_spike_on_the_boundary.step`) loses its
    two-edge tail, 13 free edges to 11, its area unchanged to 1e-9 and
    its rim one closed contour. The five slit loops of
    `slit_loops_are_not_holes.step` stand 0.2 to 0.6 µm apart against
    edge tolerances of 1e-7 and stay, as the tolerances say. DoomCube
    and Voron: no solid changed (0 tails, 0 wireless faces), and the 24
    open sheets' free edges chain into 72 closed contours, none open.
- **`make_pipe` along a wire, and profiles brought onto their spine.**
  (#394)
  - `make_pipe` read one spine edge. It now takes a wire as well, swept
    exactly as the pipe shell of the disc square to its start: an open
    wire of lines and arcs meeting tangent (cylinders and torus segments),
    or lines with sharp corners, open or closed (cylinders mitred). A
    wire with arcs that closes or turns a sharp corner is refused by name
    for `make_pipe_shell`, which takes a tolerance. Measured: three
    mitred legs 27 long and a line, quarter arc and line 10 + 3π/2 long
    hold `π r²` times their length to 1e-12, every face a plane,
    cylinder or torus.
  - `place_on_spine` moves a planar profile rigidly onto a spine's start
    as `PipePlacement` asks: `correction` turns it about its centroid by
    the least rotation standing it square to the start tangent,
    `contact` moves its centroid onto the start point. Measured: a square
    drawn leaning off the spine, placed both ways, sweeps the prism of
    volume 16 with its centroid at (0, 0, 2) to 1e-12.
- **Pipe shells: scale and twist laws beyond linear, and evolving
  sections round sharp corners.** (#394)
  - `make_pipe_shell_laws` takes a `PipeLaws`: a `SpineLaw` for the
    scale and one for the twist, each `Linear` to an end value or
    `Through` given `(share, value)` points along the monotone cubic
    Hermite interpolant, which keeps between neighbouring values. A scale
    law through values must stay above naught; shares must rise to one.
    Linear laws build as `make_pipe_shell_evolving` does, which now
    delegates to it.
  - Round a smooth closed spine the twist may follow any law ending on a
    whole turn, or on a part turn the section's symmetry closes.
  - Round sharp corners an evolving face, which was refused, now sweeps
    run by run with mitred, extended or round corners: each run under the
    laws over its share of the length, from the section it starts with.
    Where the laws stand still at a run's ends, that stretch is the exact
    prism or pipe shell of the section.
  - Measured against closed forms: a disc scaled 1, 2, 1 down a line
    lands within 1.6e-3 of `π r² h · 274/105`; a square twisted out to
    0.3 rad and back holds the prism's volume to 2.3e-4 and reaches
    `cos 0.3 + sin 0.3` to 6.5e-5; a square scaled linearly to a half
    round an L lands to 1e-9 of its mitred and rounded closed forms; a
    square twisted 0.3 rad before an L's corner lands within the skin's
    bound of its mitred and rounded closed forms.
  - A skinned run whose corner pieces agree only to the skin's tolerance
    (a section with a corner on the round corner's axis) is beyond the
    boolean and stays refused.
- **`fix_shape` mends a face's boundary that crosses itself in its
  surface's chart, where the mend is the one the boundary says.** (#397)
  - A wire crossing itself was named by nothing and mended by nothing: a
    10 mm bow tie checks valid, and measures 1.3e-13 mm² of area.
  - `ogeom_algo::chart_crossings` finds where a face's boundary crosses
    itself in its chart: two edges' pcurves crossing at an angle, not at
    a vertex of both (within ten times the joint's tolerance, as for a
    gap in the chart), the walk carried across a closed surface's seam
    and each pair also asked a period apart. Two consecutive edges whose
    loop from the crossing to their shared vertex stays within twice
    their tolerances together close a sliver the images of two curves
    meeting there may close, and are no crossing. `check` reports each
    crossing as broken on the face, where its edges cross.
  - Where a face's boundary crosses itself once, within one wire and one
    copy of the chart: two consecutive edges crossing short of the vertex
    they meet at, the crossing and each edge's run past it within the
    trim cap, are cut back to the crossing and the vertex moves there,
    widened to reach where it stood plus what it stated, so every other
    curve ending within it still does. Any other crossing splits the
    face into the figure eight's two loops, each a face facing as the
    face did and wound as an outer boundary; a face with holes, or a
    face handed over alone, is not split. Every face a mend changes must
    check clean and cross itself nowhere after, or nothing is changed.
  - `FixReport::crossings_trimmed` and
    `FixReport::faces_split_at_crossings` count them, and
    `FixReport::crossings_left` counts the faces still crossing
    themselves.
  - Measured: a 10 mm bow tie in a compound becomes two triangles of
    25 mm² each, both facing up; a 50 mm sliver triangle whose edge out
    of its corner starts 1 µm off it, crossing the edge in 0.05 mm from
    the corner, is cut back to the crossing and measures the closed
    form of the triangle left, to 1e-9.
  - A shape whose boundaries do not cross themselves takes the same path
    and comes out the same.
- **`fix_shape` makes the edge a wire lacks where a straight one provably
  fills the gap.** (#397)
  - A gap between two edges of a wire wider than the trim cap, within
    which the two vertices merge, stayed reported, and no edge was made.
  - Such a gap now gets a straight edge, with a pcurve on each face it
    bounds, where that edge is beyond doubt the one missing: the gap is
    shorter than every edge of each face it is in, no other vertex of
    those faces stands as near either end as the other end does, the
    line lies on each face's surface within the face's tolerance along
    the chart segment between the edges either side, and each face,
    filled, checks clean. A gap two faces share is filled by one edge on
    both or not at all. Any other gap stays reported.
  - `FixReport::lacking_edges_made` counts them.
  - Measured: a 10 mm box gaping 2 mm along an edge closes and measures
    1000 again; a 10 mm square missing a corner across √2 measures
    99.5. A chord 2 mm across a 10 mm cylinder's rim, 0.05 mm off the
    surface (its sagitta), is refused, as is a gap as long as an edge of
    its face.
  - A shape with no such gap takes the same path and comes out the same.
- **Pipe shells: part-turn twists that the section's symmetry closes
  round a closed spine, and round corners for holed sweeps under the
  auxiliary and binormal laws.** (#394)
  - `make_pipe_shell_evolving` refused any twist round a smooth closed
    spine that was not a whole number of turns. A part turn that brings
    each loop of the section home onto itself is now built: a square a
    quarter turn, a circle about its centre any angle. A faceted loop's
    edges each sweep a strip open at both ends, from the edge at the start
    to the edge it comes home on, and the strips close on one another. A
    smooth loop is read from a point that slides round it as it turns, so
    each column of the skin comes home to itself. A part turn the symmetry
    does not close (a rectangle a quarter turn) is still refused by name.
  - A round corner checks that the next run starts on the last run's end
    section turned about the corner. That check compared within ten
    confusions, but a skinned run is held to the skin tolerance, and a
    holed section's end under the auxiliary and binormal laws stood off
    the law's placement by up to that. The check now compares within what
    either run is held to: the skin tolerance for a skinned run, ten
    confusions between prisms.
  - Measured against Pappus: a 2×2 square twisted a quarter turn round a
    circle of radius 10 (4.3e-6 relative), three quarters back (1.9e-7),
    and with a centred 0.8 square hole under the level binormal (6.9e-7);
    each reaches the half diagonal √2 in `z` to 5e-7, as only a turned
    square does. A circle twisted a radian is the torus (2.8e-6). The
    holed square down an arc into a quarter turn, rounded under the level
    binormal, is valid and measures within 5.2e-5 of the rotation-
    minimizing sweep (bound 1.6).
  - Every sweep built before takes the same path and builds the same
    solid.
- **IGES network subfigures, and a group as an array's base.** (#404)
  - A network subfigure definition (320) and its instances (420) were
    skipped. A network definition now reads as a subfigure definition
    does, its members built once, and each physical instance places it
    by its scale, matrix and translation in the order IGES 5.3, section
    3.6.2 gives. Connect points (132) are not geometry and are not read.
  - Refused by name: an instance flagged logical, by itself or, where it
    leaves its type unspecified, by its definition (a schematic's, not
    the model's), and an instance whose axes scale unlike, which would
    distort what it places.
  - A group (402, forms 1, 7, 14 and 15) as the base of a rectangular or
    circular array placed nothing. Each of its members is now placed at
    every position the array keeps.
  - Measured on a synthetic deck: a unit patch instanced at three times
    its size reads area 9 centred at (101.5, 1.5); an array of a group of
    a unit patch and a 2-cube places both at each of its two positions.
  - Decks without these entities read as before.

- **Pipe shells under the auxiliary and binormal laws round smooth closed
  spines, with holes along open ones, and twisted whole turns round
  closed spines.** (#394)
  - `make_pipe_shell_with` refused an auxiliary or binormal law round a
    smooth closed spine unless its frames were the rotation-minimizing
    ones turned by a fixed angle. Those laws set a frame at every point of
    the spine, single-valued round it, so the skin now closes on the
    law's own frames and is held to the tolerance against the law there,
    each hole a void tunnel. A circle and a stadium of arcs and lines are
    both swept.
  - Those laws refused a face with holes along an open spine. Each hole
    now sweeps its own solid under the same frames and is cut from the
    outer one.
  - `make_pipe_shell_evolving` refused every closed spine. Unscaled and
    twisted whole turns, the section comes home to itself, and the skin
    is closed round the spine. A scale, or a twist of part of a turn, is
    refused there by name.
  - Measured against Pappus's `A L`: a rectangle round a circle under a
    leaning binormal (1e-7 relative) and toward an off-centre auxiliary
    circle (1e-8), round a stadium (2e-8), and twisted one and two whole
    turns (8e-8, 2e-6). Each reaches its half diagonal in `z`, as only a
    turned section does.
  - Every sweep built before takes the same path and builds the same
    solid.

- **`remove_faces` removes a slot cut right across a face, and a rebate
  round a drum's rim.**
  - Removing a slot across a box's top was refused: the top's two halves
    lie on one plane and meet in no curve. They now come back as one face
    through both halves' surviving edges, and the top's front and back
    edges run on across the slot's mouths.
  - Removing a rebate round a drum's rim was refused, because the cap and
    wall meet along a circle no other face stops. The circle is now taken
    whole, cut where the wall's seam reaches it, and the cap grows back to
    the rim.
  - Both give back the original exactly: the box (3000 to 1e-12, six
    faces on its own planes, twelve edges once merged) and the drum (its
    volume to 1e-12, three faces on its own surfaces), each valid.
  - Every removal that worked before takes the same path. Same-surface
    merging is tried only after both existing closures refuse.

- **A solid whose boundary passes through itself resolves to the union of
  its material: `ogeom_bool::self_union`.** (#401)
  - Before, such a solid had no mend: `check` reported its faces passing
    through each other, and orientation findings the crossing causes.
    `fix_shape` left it as it came.
  - Now the solid is fused with itself. Each face is split by every face
    it shares no edge or corner with. A piece is kept where the winding
    number just beyond it is zero and just behind it one or more. It is
    read a step off each side, along a segment that meets no other face.
    Faces lying over one patch facing the same way keep it once, and a
    wall lying twice facing either way goes.
  - Opt-in. A solid `self_crossings` finds nothing in comes back as it
    is. Refused by name: a face crossing a face it shares an edge or
    corner with, or crossing itself, and a piece whose sides no step reads.
  - Two boxes joined as one solid come out as their union, crossing and
    sharing a floor and roof, to 1e-9 of 10.5 and 10. So does a strip
    wound back through itself in one shell, to 1e-9 of 9.8, with or
    without its sewn joint wall left inside twice.
  - Voron 2.4 solid #907 checks valid after, with no problems found; it
    had 15. Its volume falls by 0.014676 from 168.018836. Sampling its
    winding number on a 2-micron grid finds the region enclosed twice to
    be 0.014672.
  - `SolidBoundary::winding` reads a point's winding number by signed ray
    crossings. `self_crossings` returns what `check` reports as the
    boundary crossing itself.

- **The IGES writer keeps surfaces of revolution, extrusions and offsets
  in their own forms.**
  - These used to be written as B-spline patches (128), exact in geometry
    but losing their kind. They now go out as 120, 122 and 140, and the
    reader brings each back as the same kind of surface.
  - An offset with an analytic form is written as that form. One over
    another basis is still written as a spline patch.
  - A revolved spline profile and a spline-walled prism round-trip with
    their surface kinds kept. An offset sheet comes back as an offset, the
    same way up.
  - A trim on a surface of revolution of an ellipse whose axes lie along
    its plane's now lifts through the ellipse's angle.
  - A solid's whole name, past the eight-column label, is written and read
    as its Name property (406 form 15) (#404).
- **IGES reads and writes wireframe, and places array subfigures.** (#404)
  - An independent curve (100, 102, 104, 106, 110, 112, 126, 130) that is
    geometry was counted as skipped, and a deck of curves only was
    refused. Each now reads as a wire in `IgesImport::wires`. Its pieces
    are chained, and a chain that comes back to its start closes. A
    parameter-space curve, a construction curve and a drawing's own curve
    are not wireframe.
  - The writer refused a part holding a free wire or edge. It now writes
    each as an independent curve, a wire's edges as a composite curve
    (102). A free vertex is still refused, and STEP still refuses
    wireframe.
  - Rectangular (412) and circular (414) array instances place their base
    at each position their pattern keeps, with the DO and DON'T lists.
  - A subfigure definition's untrimmed surfaces, constructive solids and
    curves are placed with it. Before, only its manifold solids, trimmed
    surfaces and nested instances were.
  - A composite curve skips the points and connect points among its
    pieces, which it used to refuse.
  - Measured: wire lengths 10, 7 and 7 to 1e-12 and a circle's `6π` within
    its mesh bound, the same after a round trip. Five patches of area 4
    and three cubes of volume 8, each at its closed-form centre to 1e-9.

- **Guided lofts take sections of several edges.** `make_loft_surface`
  and `make_loft_solid` with guide curves refused any section of more
  than one edge, so a guided loft through rectangles was not built. Now
  the sections share one edge count and a guide runs through every corner
  where two edges meet. Each edge's wall is the Gordon patch through that
  edge of every section, following the guides through its two corners and
  any guide crossing it inside; neighbouring walls meet along their
  corner's guide. A rectangle along four bowed corner arcs, through two
  sections or three, measures its closed-form volume to 1.2e-15, and
  along straight corner guides it is the pyramid's frustum. A corner with
  no guide, two guides through one corner, a guide changing corners, and
  sections of different edge counts are refused by name (#394).

- **The Frenet frame turns sharp corners extended or rounded.** It was
  refused round any sharp corner other than mitred. Now each run between
  corners places the section by its own Frenet frame and the runs are
  joined as under the binormal and auxiliary laws. A square down two arcs
  bending the same way, meeting at a sharp quarter turn, rounds and
  extends to its closed forms to 6e-16. A straight run between sharp
  corners has no Frenet frame and is refused by name (#394).

- **Pipe shells turn sharp corners against curved legs, round closed
  spines, and under the binormal and auxiliary laws.**
  - Before, extended and round corners were built only for a face down an
    open spine of straight legs with the rotation-minimizing frame. The
    binormal and auxiliary laws swept only an open spine with no sharp
    corner.
  - Now the spine is cut at its sharp corners and each smooth run is swept
    on its own. The runs are joined at each corner: mitred on the halving
    plane, extended, or rounded.
  - A rectangle under a tilted binormal, down an L and round a square
    ring, mitres to its legs' volume to 1e-15. An arc into a sharp turn
    rounds and extends to its closed forms to 1e-15, and so does a square
    ring.
  - A binormal along a circle's axis gives the torus by Pappus.
  - Refused by name:
    - a tilted binormal round a circle;
    - a round corner where the law sets a different next section (#394).

- **`make_pipe_shell_evolving` scales and twists a pipe's section along
  its spine.** A `PipeEvolution` gives the scale and the twist reached at the
  spine's end, each growing linearly with the length run. Scaled without a
  twist down one straight edge, the result is the exact ruled loft to the
  scaled image. A square scaled to 0.4 is the pyramid's frustum, and its
  volume is exact to 1e-15. A disc gives the cone's frustum, and a scale of
  zero gives the pyramid or the cone. Twisted or curved sweeps are skinned
  and held to the tolerance against the law. A square twisted a quarter
  turn holds the prism's volume to 2.4e-6. A square scaled round a quarter
  arc holds Pappus's volume to 9e-6. Holes are swept as their own solids
  and cut out (#394).

- **`make_pipe_sections` takes sections with holes.** A face with a hole
  was refused. Each hole is now swept through the sections like the outer
  loop and cut from it; every section must have the same number of holes.
  Two holed squares, the second half the size, give the pyramid's frustum
  minus the cone's, to 1.2e-6 (#394).
- **IGES reads untrimmed surfaces, copious-data paths, and trims on a
  revolved circular arc.**
  - An independent B-spline, ruled, revolved, tabulated, parametric spline
    or offset surface with no trim was dropped and counted as skipped. It
    now reads as the face its own extent bounds, sewn with the rest.
  - Copious data drawn as a path (106, forms 11, 12, 13 and 63) reads as
    the exact polyline through its points, as a model-space curve or as a
    trim.
  - A trim given only in the parameters of a surface of revolution whose
    generatrix is a circular arc was refused. It now lifts through the
    arc's angle, measured from the start angle IGES specifies.
  - The writer writes a face with no boundary as the trimmed surface its
    whole surface bounds, where it refused before.
  - A torus patch trimmed in its parameters, an untrimmed patch and drum,
    and a copious-data triangle each read to their closed-form areas and
    survive a round trip through the writer (#404).

- **`fix_shape` gives a boundary that crosses a pole its degenerate edge
  back.** Take a sphere whose wire walks its seam up and down with no edge
  at either pole, or a cone with none at its apex. In space the boundary
  closes, so `check` passes it, but its chart crosses the pole's row by a
  whole turn. The shortest-way reading of that crossing is wrong: the
  sphere meshed to nothing and measured no volume, and the cone measured
  12.532 where πr²h/3 is 12.566.

  A joint is at a pole where its two pcurve ends stand apart along one
  parameter direction and the whole row through them, sampled across the
  surface's domain, lifts onto the joint vertex within that vertex's
  tolerance. At each such joint the crossing may be the shortest way or
  a turn either way more. The wire is solved for the one choice under
  which its unrolled chart closes and encloses area the right way round.
  Where that choice is not the shortest way, each pole gets a degenerate
  edge spanning its crossing. The step runs after short edges collapse,
  so it also closes the gap a collapsed edge at a pole leaves.

  Both synthetic shapes then measure their closed forms. In the DoomCube
  350 assembly, three spherical caps imported this way (solids 709 and
  711) measured about 3% short: 13.115 against 2πR²(1 − sin v₀) =
  13.525, and 12.188 against πR² = 12.566. Each now measures its closed
  form. The report counts the edges (`pole_edges_made`). This is #397's
  degenerated-edges item.

- **`fix_shape` closes a gap in a wire where two edges meet at two
  vertices.** One edge ending at a vertex and the next beginning at
  another, as edges built apart carry them, is a gap: `check` reports it
  Broken however close the two stand, and `fix_shape` left it. A wire in
  order whose vertices only coincided made `fix_shape` fail outright,
  since putting it in order rebuilt it through `make_wire`, which refuses
  the gap.

  Now the two vertices become one where they stand no farther apart than
  a fitted trim may stand off its surface, the reach within which a
  curve ending away from its vertex is stated as slop. The vertex kept
  widens to reach the one it absorbs, plus that one's own tolerance. A
  class of vertices joined this way is left as it is where merging would
  close an edge on itself or pinch a face's boundary at a point, and a
  gap past the reach stays reported: no edge is made to fill it. A box
  with an edge cut a hundredth of a millimetre short at a vertex of its
  own, and one whose edge begins at a coincident vertex of its own, both
  check clean after and measure 1000. Neither real-world assembly has
  such a gap after the read, since STEP shares vertices by reference;
  both are unchanged. The report counts the vertices merged
  (`gap_vertices_merged`). This is #397's wire-gaps item.

- **`make_loft_solid` lofts a solid through sections of different edge
  counts, and along guide curves.** `make_loft` refused a square to a
  circle, a slot to a hexagon, or any section mixing lines with arcs, and
  `make_loft_skinned` could only fit such sections, approximately and to a
  tolerance. `make_loft_solid` closes the loft sheet (`make_loft_surface`)
  through N closed sections with the first and last sections' planes. It is
  ruled or smooth, the sections split to match by arc length, every piece
  exact. A ruled square to a circle measures its closed-form volume to
  1e-9, and a square of four edges lofted to one of eight is the prism.
  With guide curves the walls are the guided Gordon sheet, capped: circles
  on a cone along two of its rulings are the frustum (#394).
- **`solids_from_shells` builds solids from closed shells, and tells
  outer shells from voids by which shell holds which.** Shells handed
  over in any order are each classified against every other, at all of
  their vertices that lie off it. A shell that no other holds, or that an
  even number of shells hold, is a part's outer boundary. A shell held by
  an odd number is a void of the deepest shell that holds it. For
  example, a block, a cavity in it, an island in the cavity and a block
  apart give three solids: the block with its void, the island, and the
  block apart. Each solid faces out, and its void faces into the cavity.
  Shells whose vertices fall on both sides of another cross each other
  and are refused (#397).

- **`fix_shape` states the slop an imported boundary carries.** An edge
  whose pcurve stands off its curve farther than the edge states, and a
  vertex a curve ends away from farther than either states, stayed
  broken: `check` reported "has a pcurve that leaves its curve" and
  "curve stops … from the vertex it should meet" after the fix as
  before. These are the two most frequent defects the STEP reader meets
  (some 27 000 and 23 000 in the two real-world assemblies), and the
  reader records them by widening; `fix_shape` now does the same for a
  shape no reader built. The edge widens to its pcurves' gap, measured
  at least as densely as `check` samples, and the vertex to the widest
  gap of the curves ending at it, each only where the gap is within the
  reach a fitted trim is allowed (a millimetre at unit scale). A wider
  gap is a missing piece, not slop, and stays reported. The report counts
  both (`pcurve_gaps_stated`, `curve_ends_stated`). The remaining
  ShapeFix gaps are ranked in #397.

- **`apply_draft` drafts a face the neutral plane does not reach.** A
  fitted patch, or a wall of revolution about an oblique neutral plane,
  that the plane missed was refused: "the neutral plane does not cross the
  drafted face". It now turns about where its surface, run on beyond the
  face, meets the plane, as a planar wall, an axial wall of revolution and
  an extruded wall already did: a drum's or a cone's rulings run on to the
  plane, a patch's crossing within its own extent. A drum and a cone
  below a neutral plane tilted 0.35 rad off their axis now draft, every
  ruling leaning the draft off the pull and meeting the plane on the wall's
  own surface. A patch that stops short of the plane is refused by name
  ("meets neither the drafted face nor its surface"), and a face the plane
  crosses in two runs stays refused, its reason now stated: a region beside
  both runs would turn about two hinges (#352).

- **`remove_faces` removes a step, a notch or a pocket cut at an edge.**
  A feature of several faces each bordering only one of the neighbours
  that must meet again was refused: a step along an edge and a notch cut
  from a corner by "an end face's surface never meets the recovered
  edge", a pocket open at an edge by "a wrapping band recovered an open
  curve". Such a wound now closes by extension from the neighbours' own
  surfaces: a step's lost crease comes back between the faces it pierces,
  a notch's three cut-short edges run on to where its three neighbours
  meet, and the edge a pocket opens on joins across its mouth. Each comes
  back as the box it was cut from, its volume exact, its six faces on the
  box's planes and twelve edges once merged. A slot right across a top,
  whose two halves must merge, stays refused (#360).

### Changed

- **A straight round beside a drum runs out across it.** (#342)
  - This is a convention change. A straight band ran out past its edge's
    end only across planes. Where a drum's cylinder followed the crease
    (either end of a D-shaped block's chord, an edge of a cube shaved by a
    drum), the band was capped square at the vertex. What of its section
    the drum still held past the vertex was left standing. On a drafted
    chord wall the ball touches the wall outside the drum at the vertex,
    and the cut refused: "the kept pieces did not close into a shell".
  - The band now runs on along its crease until its whole section has left
    the drum, as a band does across a ramp. The section's last point to
    leave a cylinder need not be a corner of it, so the reach is the
    farthest any point of its hull leaves, found in closed form along the
    hull's sides. This matches the rule for a crossed convex turn: each
    band as it runs out alone, the cuts taken together. A crease running
    on into a boss, and a band meeting another band of the request at that
    vertex without a junction, are capped as before.
  - A straight band between a wall and a drum along its ruling is now set
    aside with the others of its round, as every other band is. Rounded
    first in a crossed junction, it was cut at once, and the band after it
    was built against its cut rather than the solid as given.
  - Measured with r = 0.5 against column-integrated references, the whole
    volume taken. A disc of radius 2.5 up to its chord at 0.5: chord wall
    square 0.262833 → 0.264802 (0.264802), drafted 0.1 refused → 0.208768
    (0.208768), overhanging 0.1 0.330381 → 0.333366 (0.333366). A 10 mm
    cube shaved by a drum of radius 6, its bottom edge and the upright at
    its end in either order: 0.379097 → 0.395766 (0.395766). The cube's
    top edge alone at r = 1: 1.423508 → 1.553030 (1.553030).
- **Two straight rounds meeting across a convex turn take the union of
  each alone.** (#342)
  - This is a convention change. Where the two bands' end sections
    reflected into each other across the plane bisecting the turn (walls
    square, or leaning alike), the two met flush in that plane. With the
    walls leaning alike with their feet out, each band's section reaches
    past that plane under the other's, and the flush meeting left a ridge
    standing that either band alone takes.
  - Two straight creases across a convex turn now always take a crossed
    junction: each band runs out past the vertex as it does alone, and
    the cuts take their union. This is the rule already used where the
    sections differ, and it does not depend on the order the edges are
    named in. Where the walls are square or lean with their feet in, the
    union is the flush meeting, and the figures are unchanged. A curved
    crease at a shallow convex turn still meets flush in the bisecting
    plane.
  - Measured with r = 0.5, both walls leaning by 0.1 with their feet out,
    both orders, against the union integrated a height at a time: 45°
    0.489549 → 0.489918 (0.489918), 60° 0.494839 → 0.495104 (0.495104),
    90° 0.500516 → 0.500669 (0.500669), 120° 0.503794 → 0.503882
    (0.503882).

- **An IGES curve on a surface (142) takes the representation its sender
  preferred.** (#404)
  - The model-space curve was always taken when present. Where the
    preference flag says the parameter-space curve composed with the
    surface is preferred (1), that composition is now lifted and builds
    the trim, as it already did where no model-space curve is given.
    Where the lift is refused or falls short of its target, the
    model-space curve stands in and a warning says so. Any other
    preference reads as before.
  - Both curves now count as read, whichever one builds the trim.
  - Measured on a synthetic deck: a patch trimmed by a parameter-space
    circle whose model-space curve is the inscribed square reads the
    disc's area `6.25π` within its lift and mesh bound when the circle
    is preferred, and the square's 12.5 when the model-space curve is.

- **The kernel's scope and its capability ledger are stated in its own
  terms.** The parity audit against another kernel's public headers is
  gone: `docs/PARITY.md`, `docs/parity/` and the header-profiling tool are
  removed. `docs/capabilities.toml` lists the kernel's 90 capabilities with
  the evidence for each, and `docs/CAPABILITIES.md` is its rendered form,
  held to the code by `tools/capabilities.py check`. Refusal messages that
  pointed at `docs/PARITY.md` point at `docs/CAPABILITIES.md`. Module
  docs no longer map modules to another kernel's packages.
- **`PROVENANCE.md` records where the code and design came from.**
  `CONTRIBUTING.md` now forbids reading another kernel's implementation
  source for work here, by people or tools.
- **A straight-seat fillet's wedge is built without measuring its
  pcurves.** Every pcurve of the wedge a fillet between two planes takes
  away or adds is its edge's closed-form image, lines and arcs on the legs
  and caps and rulings and circles on the blend cylinder, so it stands
  within rounding of its edge and stating the gaps widened nothing. The
  measurement sampled each pcurve some 580 times, about 5% of rounding all
  twelve edges of a 10 mm cube. Wedges that may carry a fitted pcurve
  still measure theirs. The cube rounds to the same 26 faces and volume.
- **A fit at fixed parameters finds each knot span's samples by
  bisection.** Refining a fit sifted every sample for every span, which
  was about half the cost of each fitted pcurve a sphere face takes when a
  boolean recharts it about a clear pole; the samples of a fixed-parameter
  fit come in order, so each span's are one run. The knots chosen, and so
  every fit, are unchanged. Rounding all twelve edges of a 10 mm cube
  fits 24 such images, in 12.4 ms rather than 18.6.
- **Straight bands seated flush on planar faces are cut by trimming those
  faces, not by the general fuse (#384).** `fillet_edges` applies each
  round's wedges through the new `ogeom_bool::cut_seated`. Where each
  wedge is planar seats and one strip of a cylinder, every seat lies on a
  face of the solid facing its way, the wedge's other edges cross those
  faces from vertex to vertex (an edge of the solid is split where a
  contact lands inside it), no other face of the solid reaches into the
  region the seats' planes and the blend's chord bound, and no new edge
  leaves a vertex along another closely enough that the general fuse
  would split the two, the faces are trimmed along the wedge's own edges
  and the blends close what is left. The result is assembled as the
  general fuse assembles its own, and its solid and history are the
  same, except that no face is recorded as generating the shell it is
  sewn into. Everywhere else the cut is the general fuse's. Rounding all
  twelve edges of a 10 mm cube takes 58 ms rather than 103, to the same
  26 faces, 48 edges and volume.
- **Rims seated on a plane and a cylinder about one axis are cut by
  trimming those faces, not by the general fuse (#384).**
  `ogeom_bool::cut_seated` now also takes a rim's wedge: an annulus of a
  plane square to the axis, a band of the cylinder and the quarter of a
  torus's tube between them. Where the annulus lies on a planar face of
  the solid one of whose loops is the rim, the band on the one other face
  walking the rim, a cylinder walking it between the two walks of one
  seam, each facing the seat's way, every other edge of those two faces
  stands clear of the seats, and every other face of the solid stands out
  of the ring the wedge's cross-section sweeps, the two faces are trimmed
  along the contact circles and the tube closes what is left. Closed edges
  holding one vertex are cut at their middle and the contacts where the
  seam and the meridian meet them, as the general fuse cuts them, and the
  solid and its history are the same. Everywhere else the cut is the
  general fuse's. Rounding a cylinder's top rim takes 1.2 ms rather than
  4.1, and all four rims of a tube 5.3 ms rather than 75, to the same
  faces, edges, volume and area.
- **Concave rims are filled by trimming the plane and cylinder they sit
  on, not by the general fuse (#384).** `fillet_edges` applies the wedges
  it fuses on through the new `ogeom_bool::fuse_seated`, which takes a
  concave rim's wedge (a boss's base, a blind hole's floor, a shaft's
  shoulder) under the conditions `cut_seated` takes a convex one's, with
  the two faces facing against the seats. The ring the wedge sweeps then
  lies outside the solid, and the union is the solid with the seats taken
  out and the tube as the wedge faces it. The solid and its history are
  the general fuse's; everywhere else the union is the general fuse's.
  Rounding a boss's base takes 1.35 ms rather than 4.7, a stepped shaft's
  shoulder 1.40 ms rather than 11.7, and a pin's base and its blind
  hole's floor 3.2 ms rather than 21.9, to the same faces, edges, volume
  and area.
- **A rim's wedge applied in one boolean with other wedges is trimmed,
  and only the rest goes to the general fuse (#384).** `fillet_edges` cuts
  wedges whose boxes stand apart in one boolean, and a group mixing a
  rim's wedge with a straight band's or a junction's went to the general
  fuse whole. `ogeom_bool::cut_seated` now splits such a tool where every
  rim's box stands further from every other lump's than the boolean welds
  past the two lumps' tolerances: the rims are trimmed out, keeping the
  faces the general fuse would set aside for the whole tool, and the
  other lumps are cut from what is left, by trimming where they qualify.
  `ogeom_bool::fuse_seated` splits a union of concave rims' wedges with
  other wedges the same way. The solid and its history are the general
  fuse's. Rounding a boss's top rim with the plate's top loop takes
  15.7 ms rather than 18.6, a tube boss's two rims with the plate's
  bottom loop 17.5 ms rather than 27.7, a boss's top rim with one
  straight edge 1.9 ms rather than 6.1, and a boss's base with a rib's
  two base edges 4.8 ms rather than 7.7.
- **The arcs of a tangent chain are rounded by trimming the faces their
  wedges sit on, not by the general fuse (#384).** Rounding a chain of
  straight edges and arcs, such as a block's top loop once its upright
  edges are rounded or a slot's rim, cuts the straight bands' wedges and
  then the arcs', or the other way round. An arc's wedge is a piece of a
  rim's between two planar ends through its axis. Where the straight bands
  went first, each end lies on the face a band's cut left, and the piece
  went through the general fuse in full. `ogeom_bool::cut_seated` now
  trims such a piece as it trims a straight band, at most half a turn,
  with its wall on a cylinder face of the solid and every other face of
  the solid standing out of the sector of the ring the piece lies in.
  `ogeom_bool::fuse_seated` now trims the same wedges for a concave chain
  (a rounded pocket's floor), and a straight band's where it is filled.
  The seat a trim reads may now lie on any loop of its face: a face's other
  loop is read as standing outside the seat where one of its vertices
  stands off the seat's box, where it was declined where the loop's box met
  the seat's. The solid and its history are the general fuse's; everywhere
  else the boolean is the general fuse's. Single-threaded, rounding a
  block's top loop once its upright edges are rounded takes 98 ms rather
  than 495, a stadium slot's rim 10 ms rather than 220, and a rounded
  pocket's floor loop 100 ms rather than 466, to the same faces, edges,
  volume and area.
- **The corners a fillet request rounds are cut in one boolean where
  they stand clear of each other (#384).** `fillet_edges` rounded each
  corner where three or more of its edges meet with its own cut. On a
  solid of planar faces whose corners each take one ball, every corner's
  tool is now built first and the tools are cut together, where each
  tool's block keeps clear of every other corner and every two balls stand
  further apart than the boolean's touching distance for a pole: each
  corner is then read, and its tool built, as it would be after the others
  were cut, and the one cut takes away what the cuts in turn do. Elsewhere,
  and wherever the one cut does not close, the corners are cut in turn as
  before. New `ogeom_bool::POLE_TOUCHING` states that distance. Rounding
  all twelve edges of a 10 mm cube takes 101 ms rather than 105.5, to the
  same 26 faces, 48 edges and volume.
- **A line touching a circle or ellipse in its plane meets it in closed
  form.** The closed form left every tangent line/conic pair to the
  sampling path, at about a millisecond a pair, and a rolling-ball blend
  is tangent to its faces by construction: rounding all twelve edges of a
  10 mm cube at radius 1 ran 992 such pairs and took 796 ms. A line whose
  quadratic along the conic has a double root to within that root's own
  rounding bound is now answered as the one touch the sampling path
  reports, its reach the same valley measured with the exact distance; a
  touch outside a trim answers empty where a lower bound on the curves'
  distance clears the gap. Near-tangent lines that miss or cross, and
  touches the sampling path would read as overlap, still go to it. The
  cube now rounds in 124 ms, to the same 26 faces and volume.
- **A section piece's hug test screens contact edges by their boxes.**
  Whether a piece runs along an edge was measured against every contact
  edge of the face, though own edges were already screened by their
  boxes; a contact edge's box now screens it the same way, with the same
  answer. The four outer edges of a plate with 256 holes fillet in 49 ms
  where they took 72.
- **A rim whose blend's rail passes a sphere's pole closer than a stride
  of the march refuses by name (#365).** The rail's chart turns half
  round within that stride, which no fit through the stations follows;
  the refusal now says so rather than reporting the fit's error. A bore
  whose wall toward the pole is thinner than the ball blends as a wider
  wall does, which a test now holds.
- **PARITY no longer lists the curved corner's straight-band-first order
  as owed (#367).** Rounding the shaved cube's corner with the straight
  band, then the ruling's, then the rim arc's, and the corner tool after,
  leaves no stub: it closes the corner the ruling-first order does, to a
  part in ten million, which the test now holds where it held 1e-4.

- **`unify_same_domain` merges the curved faces it still left split.**
  - Faces that already walk a seam merge. Two cylinder bands stacked end
    to end become one band, its seam running on from one band's into the
    other's. Two half tori that each cover the whole tube become one torus
    face. Bands whose seams start at different angles round the circle
    they share stay split: one face would need a seam running along part
    of that circle. So do bands whose seams meet but run opposite ways
    from the vertex between them.
  - Extrusions on an ellipse, a hyperbola or a parabola are compared, as
    is a basis curve whose parameter runs backwards or is shifted along
    itself. An elliptic wall fused with its mirror image is one face
    again, a tube closed by a seam.
  - Two faces on one surface that close round its axis merge. The edge
    they share may carry one pcurve for both, which one face reads a
    period away from the rest of its boundary.
  - A seam two faces make takes its second side as its first moved a
    whole period across, so the two faces' pcurves may span any ranges.
    Before, both had to be chart lines of one length.

  Each map between charts is exact, and a pair no exact map relates stays
  split. Where a cluster holding such faces does not merge as a whole,
  the clusters its other faces make on their own are merged as before
  (#315).

- **`check` reports an edge whose pcurve stands farther off its curve
  than the edge is long (#314).** Such a pcurve describes another stretch
  of its surface: a fit whose end landed on the far side of a closed
  chart runs the far way round it. A tolerance widened to cover the gap
  hid it from every containment test, and the next boolean welded the
  face's boundary into dust. An exact pcurve is never reported, but one
  only within its tolerance of the edge can be: the test is a heuristic,
  and it is reported as suspect, not broken.

- **The licence changes to LGPL-2.1-only.** It was MIT OR Apache-2.0.
  ogeom was built with Open CASCADE Technology as its reference kernel:
  the public API vocabulary, the data-model conventions and the default
  tolerances follow it, the project's first goal was interoperability
  with OCCT-based software, and for one day OCCT's implementation source
  was read. That makes the kernel a derivative work in the licensing
  sense, whatever the side-by-side audit found (it found no translated
  code). LGPL-2.1, the licence OCCT itself carries, states that
  position, and `PROVENANCE.md` has the full record. Copies of earlier
  releases keep the terms they were distributed under. The rule against
  naming another kernel in committed files is dropped with it: the
  licence and the provenance name OCCT, and code and docs write in the
  field's vocabulary.

### Fixed

- **A profile with any edge in its plane sweeps into a valid solid or
  refuses.** (#466)
  - A prism, and a revolution through less than a turn, keep the profile
    as their end caps. A cap edge with no closed-form pcurve on the
    profile's plane was left bare, and the sweep returned a solid `check`
    rejected ("bounds a face it has no pcurve on"). This happened for a
    hyperbola arc, a reversed trimmed curve and a curve on a surface.
  - A hyperbola in a plane now takes its exact pcurve there, a new
    planar `Hyperbola2d`, at the same parameter. Run backwards, it takes
    one only over a domain symmetric about zero. The native format
    carries the planar hyperbola, and the `.brep` format carries it when
    it runs forward. A STEP file states no pcurve for it, and the reader
    projects the edge again.
  - Any other cap edge with no closed form takes a pcurve fitted at the
    curve's own parameters. The edge and its vertices widen to the
    measured miss of the lifted fit. An edge standing off the profile's
    plane by more than its own tolerance is refused, with the distance
    named. Exact pcurves are unchanged. A full turn keeps no caps and
    fits nothing.
  - Measured: a diamond with sides of every curve type that lies in a
    plane, swept into a prism and turned 2π and 1.5π, is valid every time,
    within 4e-15 relative of the area times the height or Pappus's volume.
    A fitted circle arc widens its edges to 3.0e-6, a fitted hyperbola
    arc to 6.2e-6.

- **A profile with offset or parabola edges sweeps into a valid solid.**
  (#464)
  - A prism, and a revolution through less than a turn, keep the profile
    as their end caps, and each cap edge takes its exact pcurve on the
    profile's plane. An offset curve and a parabola had none, so the caps
    came back with bare edges and `check` rejected the solid ("bounds a
    face it has no pcurve on").
  - An offset whose reference stands square to the plane now takes its
    basis's pcurve offset in the chart. A parabola takes the one quadratic
    Bézier it is in its own parameter. Both are the same point as the
    curve at every parameter. An offset whose reference leans off the
    plane's normal far enough to lift it past the confusion distance
    still has none.
  - Measured: a diamond with sides that are offsets of B-spline lines,
    lines, trimmed B-splines and half circles, trimmed offsets, or
    parabolas, swept into a prism and turned a full turn and three
    quarters, is valid and within 1e-15 relative of the area times the
    height, or of Pappus's volume.

- **A border of one surface lying in the other is a section, whatever
  the pair.** (#462)
  - Where a border of one surface's chart lies in the other surface, the
    two meet along it, and the marcher did not find it: no pair of cells
    crosses there, the border pierces the other nowhere, and a walk along
    it leaves the chart at its first step. #450 answered this for a plane
    against a B-spline patch or an extrusion of a B-spline profile only.
    An extrusion of a trimmed or offset profile met a plane through its
    end ruling nowhere. A plane through the side edges of a prism with
    trimmed sides, a drum through a spline-sided prism's side edges, and
    a plane on the rim between two revolved spline faces were refused
    ("the kept pieces did not close into a shell").
  - The marcher now samples each border of either chart. A stretch lying
    in the other surface within the confusion distance, by projection at
    every point, and leaving it at an angle (the normals' sine past 3e-2,
    the threshold that tells a contact from a crossing) is traced along
    the border as a branch. A border left tangentially is still a
    contact, as before.
  - Measured: the diamond prism with trimmed half-ellipse or trimmed
    B-spline sides, halved through two opposite corners at five turns,
    within 1e-15 relative of half of `2 R² h` and the bulges; the
    B-spline-sided prism and a drum through two of its side edges, at
    five turns, within 2.3e-11 relative of the triangle and the circular
    segment; the spindle halved on its rim, whole and three quarters,
    upright and leaned, within 1.7e-10 relative of `1.125` times its
    turn.

- **Slivers beside a section touching a drum's rim are kept.** (#460)
  - A drum ending in a plane tangent to a disc's wall was still refused
    in 6 of 64 placements ("the kept pieces did not close into a
    shell"), every retry included. Where the drum's seam crossed the
    section within 1e-4 of the rim it touches, the stretch from the seam
    to the touch was read as running along the rim, at the loose width a
    fitted section is measured at, and split the disc's wall but not the
    drum's. Where the drum's wall stood past the disc in a sliver 5e-4
    wide, the touches bounding it welded the whole wall at 6.6e-3, how
    far they run along the rim, and the sliver was dropped.
  - On the last retry a fitted section now runs along an edge only within
    the sum of its tolerance and the edge's, and a junction's run and a
    fitted section's reach weld only the strand ends they concern, not
    the face. Results the earlier attempts build are unchanged.
  - Measured: all 64 placements cut and common either way round, valid,
    within 9e-11 of `∫∫ (√(R² − y²) − x₀) dy dz` over the drum's section.

- **A section touching an edge of a surface ends on the touch.** (#457)
  - Where the section of two surfaces touches an edge of one without
    crossing it (a drum's rim in a plane tangent to a bore), the walk
    either stopped where its step was held on the edge, 1e-5 to 1e-2
    short of the touch, or stepped over the touch and ran on. A chord's
    round capped tangent to a pocket's bore left its section's end that
    far off the rim, and the cut was refused ("the kept pieces did not
    close into a shell").
  - The touch is now solved as one: where the edge, as a curve, touches
    the other surface, from beside it and kept only where Newton's
    linear model holds over the solve. A walk ends there, a branch run
    over a touch is cut at it, and the arc beyond a touch ending one arc
    alone is followed from it.
  - Measured: a disc less a drum capped in a plane tangent to its wall
    cuts and commons to the closed form in elliptic integrals within
    5e-11, either way round. The chord's round capped tangent to the
    bore, at drafts 0.05 and 0.1, cuts and commons either way with cut
    plus common equal to the solid within 1.1e-9.

- **A drum's section with another is found where it lies in a strip of
  the drum's height, both faces part of a turn.** (#456)
  - The boolean marches a pair of surfaces over windows of their faces.
    A window of a drum or a cone less than a turn wide was a trimmed
    view without rims, so the seeding never pierced them. A short arc
    of the section inside one pair of sampling cells was then never
    traced: a half drum ending a hair past a quarter disc's bore, or a
    chord's round ending its cap a few microns past a bore at a pocket's
    corner. The cut was refused ("the kept pieces did not close into a
    shell").
  - The arcs of a windowed drum's or cone's rims that the window spans
    are now pierced as its whole rims were.
  - Measured: the quarter disc less the half drum ending 1e-5, 1e-4 and
    1e-3 past the bore cuts and commons to its closed form within the
    drum's sliver past the bore. The chord's round cut from a drafted
    pocket after the arc's (chord at −0.5, draft 0.05) now closes in
    both orders, to 120.170582 either way, with cut plus common equal to
    the solid to 1e-9.

- **Rounds at a bore meeting a drafted chord through its axis or past it
  take the corner in either order.** (#342)
  - In a disc's pocket up to a chord whose wall is drafted, the corner of
    the arc and the chord is a crossed junction. Before the arc's cap, the
    arc's band was widened by what stands off the bore only over its own
    part, a few degrees round the axis. With the chord through the axis
    or past it, the drafted wall leans inside the bore's cylinder along
    the whole chord, beside the chord's own band. The junction refused
    there, and the request fell back to capped bands. Those bands fell
    short, and with the arc named first, or at a 0.05 draft, the boolean
    refused them ("the kept pieces did not close into a shell").
  - The widening now covers everything off the bore before the arc's cap.
    The corner takes the chord band's run beside it, which lies above
    both bands. A widening that misses what stands past the chord's wall
    and beyond the cap is now read as holding nothing there, where before
    the next boolean refused it.
  - Measured within a 2.4 box round one corner at r = 0.5, against a
    column-integrated reference, in both orders. Chord through the centre:
    draft 0.1 gives 0.127166 (reference 0.127167; main refused, or
    0.117731), and draft 0.05 gives 0.134830 (0.134829; main refused).
    Chord at −0.5: draft 0.1 gives 0.133610 (0.133610; main 0.121966),
    and draft 0.05 gives 0.141499 (0.141499; main refused). A chord at
    0.5, a square chord wall, and an overhanging one, which still falls
    back by convention, keep their figures.

- **A plane through a spline-sided prism's side edges cuts it.** (#450)
  - A plane holding an extrusion's direction and passing through the
    edge where a side face ends met that face nowhere: the extrusion of
    a spline side was marched, and no pair of cells crosses on a line the
    plane holds exactly. A spline patch whose border lay in the plane was
    marched the same way, its border line refused for want of a cell
    beyond it. A prism of a diamond with straight B-spline sides, cut on
    the plane through two opposite corners, was refused ("the kept pieces
    did not close into a shell"); with line sides it cut exactly.
  - A plane meets a spline patch, and now the extrusion of a spline
    profile read as that patch on its own chart, in the whole iso lines
    it holds, exactly, the extrusion's columns stated as its straight
    rulings. A line on an open border of the chart stands where the
    surface leaves the plane there at an angle, the normals' sine past
    the 3e-2 the marcher takes as tangent all along it; a border left
    tangentially is still marched.
  - Measured: the diamond prism, turned five ways, halved on the plane
    through two opposite corners and cut on the plane through one corner
    and the middle of a far side, common and cut, gives the closed-form
    volume to 1e-9 relative, as with line sides. A square twisted
    exactly an eighth turn before a sharp corner, which #432 left refused
    when rounded, now mitres and rounds within 1.5e-4 of the closed form
    (bound 0.19).

- **A mitre where a skinned run meets a straight one, under the
  auxiliary and binormal laws.** (#394)
  - A curved run under those laws is skinned, and ends on the skin's fit
    of each loop of the section. Fits of a hole or a smooth loop wobble
    across the exact loop the law places on the straight run, about
    1.7e-5 apart and crossing it eight times a quarter turn. Run on past
    the corner, both cut the mitre plane, and the fuse refused the
    slivers between them ("the kept pieces did not close into a shell").
  - The straight run's prism now takes, for each loop it shares with the
    skinned run within the skin's tolerance, the skin's own fit mirrored
    across the mitre plane, so both runs cut that plane alike. Where the
    straight run comes first, the skinned run is run back from its own
    start. The straight run stands off the law by at most what the skin
    already does.
  - Measured against a column-integrated reference (Pappus round the arc
    and the corner's wedges): a holed square and a disc round an arc into
    a line, level and leaning binormal and an auxiliary spine, arc or line
    first, within 2.9e-4 of the reference less the arc's own skin error.
    Two skinned runs meeting at a mitre with a loop in common are still
    refused.

- **A prism's or revolution's wall over a stretch of an unbounded line
  ends where the edge does.** (#445)
  - The edge's range was carried from the stored curve's domain to the
    placed one's by an affine map from the domain's far end. On a line's
    default domain of a billion either side, that moved an end by up to
    the 1.19e-7 spacing of `f64` there, past the edge's tolerance of 1e-7.
    The wall's pcurves then left the edge, and `check` rejected the solid.
    The range is now carried by `carried`, which leaves it as it is where
    the domains agree and scales it about zero under a scaling.
  - Measured on a rib's top creases rounded into a boss (r = 0.85) with
    the turn between them swept from -30 to 34 degrees: the prism, and
    each round of it, was invalid at -30, -1, 7, 10, 18, 22 and 32
    degrees (gaps of 1.12e-7 to 1.79e-7) and is valid at all of them now.
    Every removed volume is unchanged to the sixth decimal.

- **`solid_from_mesh` builds no facet whose boundary crosses itself.**
  (#397)
  - A facet all but tangent to a curved face beside it met that face's
    surface along a curve that ran across the facet's other side a
    little short of the corner they share, and the facet's boundary
    closed a figure eight. On NIST FTC-07 a cylinder's section crossed a
    facet's straight side 0.065 mm from their corner, closing a loop
    7.6e-4 mm wide where every tolerance is 1e-7: 31 crossings on 24
    facets. The loops are narrower than the 1e-3 mm chord the fold check
    draws boundaries at, and nothing saw them.
  - The build now asks every face for its chart crossings
    (`chart_crossings`) once the fold check passes. A facet crossing
    itself beside one curved face is built as a fan onto it, as a facet
    folding or turned in already is: the ruled surface from its far
    corner to the seam drawn straight in the curved face's chart, whose
    boundary is its chart's square. A crossing no fan mends facets the
    curved faces at it, as a fold does.
  - Measured: FTC-07 and its six relistings come back with 0 crossings
    (30 to 33 before) and the same faces, 607 with 172 curved and none
    faceted on the plain listing; its volume moves from 1726246.6420 to
    1726246.6577 against the part's 1726246.6525, its area from
    391367.0891 to 391367.0918 against 391365.3670. The rough rounded
    box's 3 crossings go, its 2123 faces and 12 cylinders stand, its
    volume moves by 1.6e-4 mm³ (3838.814744 to 3838.814902 against the
    mesh's 3835.757990) and its area by -1.0e-4 mm².

- **A D-shaped block's arc and chord round alike in either order.**
  (#342)
  - Where a drum's rim meets a straight crease across a convex turn, the
    rim band's curvature keeps the two from meeting in a plane, and no
    junction was taken: the later band ran on through the earlier. With
    the chord through the drum's axis, the arc's cap lies in the chord's
    wall, and the arc rounded after the chord left part of the corner
    standing.
  - Such a turn is now a crossed junction, as one between two straight
    creases whose sections differ: each band is built from the solid as
    given and runs out as it does alone, and the cuts take their union.
  - Measured on a disc of radius 2.5 up to its chord, r = 0.5, within a
    box round one corner against a column-integrated reference: chord
    through the centre, arc then chord 0.119607 and chord then arc
    0.105795 → 0.119607 both (0.119607); chord at 0.5 unchanged at
    0.118374 (0.118375).

- **Booleans whose fitted and exact edges part tangentially, leaving a
  lune between them, build.** (#433)
  - A boolean that refuses is now tried once more, last, reading every
    contact at the radii its curves hold to (their stated tolerances, or
    the measured distance from their faces where wider) instead of the
    1e-4 width a fitted curve was measured within. A stretch two edges
    share ends where they part by those radii, measured once for the pair
    so both sides split at one vertex; a shallow crossing that cannot be
    told from a vertex of either edge is placed at it; a contact near an
    edge is read inside or out at that edge's own curve.
  - A piece thinner than the band the other solid is classified within
    (area at most the band times half its perimeter) is read from the
    other solid's edges that bound it: past a convex edge it lies
    outside, past a concave one inside. Edges that disagree still refuse.
  - The attempt runs only where every earlier one refused, so no boolean
    that built before changes.
  - Measured: a torus elbow and a drum mitred and fused face to face
    built in one direction of seven; all seven now build valid, each
    within 1.9e-9 of its pieces' volume sum. Pads on the unit circle and
    on an ellipse 2e-5 and 3e-5 longer, stacked, fuse to π(2 + w) within
    1e-14, in both orders.

- **A face with no wires is measured in closed form.** Its area was
  read off its mesh: a sphere of radius 2 measured 48.60 at the default
  deflection against 4πr² = 50.27. Where its surface's domain is
  bounded, `surface_properties` now integrates it over that whole
  domain, as the mesh draws it, and reports a deflection of zero.
  Volume still meshes such a face, since no edge of it says its shell
  closes.

- **IGES parents' matrices carry what they own, as IGES 5.3, section
  3.2.3, table 4 says.** (#404)
  - A trimmed or bounded surface (144, 143) ignored its own matrix. It
    now carries its surface and its boundaries.
  - A subfigure definition (308, 320) ignored its own matrix. It now
    carries its members before an instance places them.
  - A drafting entity's matrix placed only its callout's plane. It now
    also carries the leaders, witness lines and symbol geometry it owns.
  - A child owned only logically (subordinate switch 02) was carried by
    its parent's matrix: a ruled surface's rails, a surface of
    revolution's axis and generatrix, a tabulated cylinder's directrix,
    an offset curve's or offset surface's base, and a composite curve's
    pieces. Such a child now stands where its own matrix puts it. A
    composite curve's piece flagged independent (00) is now carried, as
    every other parent's is.
  - An independent bounded plane (108, form 1) was counted as skipped. It
    now reads as the face its bounding curve encloses, the curve carried
    by the plane's matrix.
  - Files whose parents carry no matrix, or whose children are owned
    physically, read as before.
  - Measured on a synthetic deck, each under a lift of 7 or a quarter
    turn: a trimmed square area 16 about (4, 4, 7), a drum on logically
    owned curves area `12π` about (0, 0, 1.5), a turned definition's unit
    patch about (-0.5, 0.5, 50), and a bounded disc area `4π` about
    (0, 0, 12).
- **Where two surfaces cross at a shallow angle, the march places each
  sample within its chord of the section.** (#356)
  - The correction onto both surfaces stopped at a residual of a
    hundredth of the confusion distance. Across a section the surfaces
    part only by the distance times the sine of their angle, so a sample
    could stand that residual over the sine across it: up to 3e-4 off a
    line where a patch crosses a plane at a slope of a few millionths,
    with chords of 1e-6. The samples zigzagged and some doubled back.
  - Where the next Newton step from a landing would move it by more than
    the chord, the correction now settles to the residual that leaves the
    point within the chord, read off Newton's linear model at the landing
    (the chord over the norm of the position gradient times the inverse
    Jacobian), and keeps the further solve only where the Jacobian
    changes over the move by at most half of itself, as seen through its
    inverse. Beside a touch, where two stretches cross within that move,
    the landing stays as it was. Every walk, seed, carry and join takes
    its march's chord; silhouettes and rolling-ball blends walked through
    the same walker take theirs.
  - Over the 36 configurations of the line grid of #356 (s 0.3 to 0.03,
    K 1 to 100, chords 1e-4 to 1e-7), samples clear of the crossings stand
    at most 0.99 chord off their line against up to 110 chords, no section
    doubles back where two did, and the length traced beyond the exact one
    falls from up to 4.6e-4 to at most 1.9e-7, all of it within a crossing's
    tolerance band.
  - A chord whose middle lies within the confusion distance of both
    surfaces is joined straight, as the join's own fallback already took
    it, so a stretch between two touches settled off their exact
    crossings by their own tolerance is still joined.
- **A drum's rim meets a wall drafted out from under the block, and a
  bore's rim meets a drafted chord.** (#342)
  - Past a drum's rim cap the band the block's crease carried on stood off
    the drum where the block's wall leans out from under it, and what
    stands there lies above both bands. Taking it needs the rim's band
    widened past the drum, a widening on the drum's own cylinder, which
    the boolean refused; the request fell back to capped bands. The rim's
    torus band is now widened by what stands off the drum over its part,
    and the open corner is what stands off the drum past the block's wall.
    Within a box round the corner, r = 0.5, wall at -0.1, drum of radius 3:
    0.131918 against 0.138465 mm³ with the block's crease at `y = 1`, and
    0.124301 against 0.128202 at `y = 2`; both orders now take 0.138463
    and 0.128200.
  - A bore, the material outside the cylinder, took no junction. Its open
    corner bends away from the tangent plane into the material, so it is
    cut by the bore's own cylinder. A disc of radius 2.5 up to a chord at
    0.5 drafted by 0.1: 0.117063 against 0.124692 mm³; both orders now
    take 0.124693. A chord wall that overhangs still meets as before.
  - In the boolean, a fitted curve near an edge was read as running along
    it wherever it stayed within the loose width fitted curves are
    measured at, even where it crossed the edge at a shallow angle: a
    band's section across a drum, crossing the drum's edge at 0.04 rad, ran
    along it for a few microns and left the strand after it unattached,
    and the kept pieces did not close. A stretch of a fitted curve that
    ends inside both curves, its distance from the edge growing there in
    proportion to the way from where they meet, is now no overlap; one
    parting tangentially still is.
  - In the boolean's check that an empty cut is right, a face's triangle
    middles lifted from its chart were read as points of its boundary even
    where a trim bending into the face put them off it, and the cut of a
    solid lying wholly in the tool was refused. A middle is now sampled
    only where it lies inside its face's trim.

- **A drum's rim meets a leaning straight crease, and a disc's pocket up
  to a chord meets at both corners.** (#342)
  - Two top creases where a drum fused onto a block meets the block's
    wall leaning in under it had end sections that differ, and the rim's
    crease is a circle, so no junction was taken: the bands met as each was
    built and left the corner standing. The rim's torus band is now taken
    from its cap to where its section stands wholly past the block's wall,
    carried on past its cap as the straight band of its end section, and
    intersected with the block's band in one boolean; the rest is joined
    back on. It is taken where the block's band carried on stands nowhere
    off the drum before the rim's cap, which a closed polygon test in plan
    decides. Within a box round the corner, r = 0.5, wall at 0.1: 0.162018
    against 0.172312 mm³ with the block's crease at `y = 1` on a drum of
    radius 3, and 0.154401 against 0.160586 at `y = 2`, what lies above both
    bands in the material, a column at a time; both orders now take
    0.172311 and 0.160586. A rim with the material without it, or a wall drafted out from
    under the block, meets as before.
  - The two creases of a pocket of a disc up to a chord meet at both their
    ends. The second junction of that loop of two was dropped as though it
    repeated the first, and its corner kept its caps: 0.130616 against
    0.140398 mm³ round one corner, in one order. With the chord through the
    centre the arc's two caps lie in one plane, and the far one was read
    for the near one, so neither corner met: 0.131285 against 0.143268.
    Both orders now take what lies above both bands at each corner.

- **A pocket off square in plan with every wall drafted its own way meets
  at every corner.** (#342) A trapezoid with corners at 63 and 117
  degrees, walls at 0.1, 0.05, -0.05 and 0.08, r = 0.5, fell back to
  capped corners: 0.635925 and 0.640714 mm³ against 0.681475, the sum of
  the creases' sections and each corner's gain above both bands. A pentagon fell back
  the same way, 0.809829 and 0.806727 against 0.856168. Both orders now
  give the ideal within 1e-7. Three defects stood in the way.
  - At the acute corners the ring's wedge holds faces a tenth of a micron
    wide, narrower than the gap the halves of the band are sewn across.
    `sew_within` took both sides of such a face for one edge and the
    short edge across its end for a point, yet kept the face, walking
    each side out and back, and the point as a closed edge. A ring that,
    with the edges shrunk to a point left out, walks each edge once each
    way now bounds nothing and goes with its face, as a ring collapsed
    onto one edge already did, and an edge shrunk to a point is walked by
    no ring.
  - Sewing across a gap widened the edge it kept to reach the one it
    dropped, but not that edge's vertices, so `check` called the ring's
    wedge broken. The vertices now widen with it.
  - The boolean measures whether a contact edge runs along a fitted edge
    of the face it lies in at 48 samples along the contact. A fitted edge
    shorter than three of those spacings, along the pocket's corner, was
    never along it, so the contact split the face again back along its
    own boundary and the face kept no piece. Where the samples find no
    stretch, the edge is now measured from its own side: where all of it
    lies within the width of the contact, it runs along it between the
    feet of its ends. A cube less a prism whose foot ends in a fitted edge
    0.0005 or 0.005 long along the cube's edge was refused and now keeps
    its closed-form volume.
  - Pockets that rounded before give the same figures.

- **IGES applies an instance's matrix and a composite curve's matrix in
  the specification's order.** (#404)
  - A singular subfigure instance (408) applied its matrix to its own
    translation as well. IGES 5.3, section 3.6.2, applies the scale, then
    the matrix, then the translation, which the matrix does not move. A
    cube placed by a quarter turn and `(0, −100, 0)` now centres at
    `(−1, −99, 1)`, not `(99, 1, 1)`.
  - A composite curve's matrix was ignored. Per section 3.2.3, it carries
    its physically dependent pieces on after their own matrices.
  - Files whose instances and composites have no matrix read as before.
- **A marched branch striding past a touch it passes inside the
  confusion band is cut there.** A branch was cut at a touch only where a
  step passed within the touch's radius, sixty chords. Where the surfaces
  stand within the confusion distance `ε` of each other for some way round
  a touch, the walk places its samples across the section only that
  loosely, and a straight stride along one arm can pass the touch further
  off. The branch then ran on past the touch, alongside the arc beyond it.
  A step now passes over a touch it passes within `√(2ε / k)` of it, `k`
  the largest bend there: that is how far a straight chord to the touch
  stays on both surfaces. The reach is held to half the way to the nearest
  other touch. On
  `z = (x³ − s²x)(y³ − s²y)` with `s = 0.1` and a chord of 1e-6, a branch
  down `x = −0.1` passed the touch at (−0.1, 0) 1e-4 off. It ran on 0.003
  beyond it, and one along `x = 0` ran 0.0026 past (0, 0.1). Both now end
  on their touches. The section traced twice goes from 5.6e-3 to 0, and no
  swept case gets worse (#356).

- **A stretch of section marched twice where the surfaces cross shallowly
  comes back once.** Two tracings of one stretch were told apart by the
  chord alone. But a marched sample is held to both surfaces only within
  the confusion distance `ε`, so across the section it is placed only to
  within `ε / sin θ`, `θ` the angle the surfaces cross at. Near a touch,
  where `θ` falls toward zero, a stretch walked twice wandered further off
  the first tracing than the chord, and both were kept. A sample now
  counts as on another tracing within the chord along the section and
  within `2ε / sin θ` across it. An arc whose loose end runs along an arc
  kept after it is also cut back. On a saddle crossing a plane along lines
  4e-3 apart, a 0.0115 stub of `x = 0` no longer comes back beside the
  whole line. On Viviani's curve with a radius of 30 and a chord of 1e-6, a
  0.012 stretch is no longer traced twice. Eight swept cases of Viviani's
  curve lose their overlap, and no swept case gets worse (#356, #327).

- **Two fillets meet at an inside corner whose walls lean alike, their
  feet out under the corner as a moulded pocket's are.** The two bands'
  end sections reflect into each other in the plane bisecting the corner,
  so they meet flush there. But each section reaches past that plane down
  its wall, and a band carried on from its cap cannot slide those points
  back onto the plane. The corner was refused, and both bands were capped
  at the vertex.
  - Along a straight crease between planes, the band is now the straight
    band of its end section on both sides of its cap. It is carried on past
    the plane and cut back to it. The two halves then join across the plane
    as before.
  - Every pocket with all its walls drafted alike fell back this way. A
    4 x 4 pocket, one deep, every wall at 0.1, r = 0.5: 0.676015 mm³ on
    main against 0.708103, the sum of the creases' sections and each
    corner's gain above both bands. A trapezoid pocket with acute and
    obtuse corners: 0.578745 against 0.614621. Both orders now give the
    ideal within 1e-5.
  - Corners that rounded before take the same path. The new path runs only
    where the old one refused; curved creases or curved hosts still
    refuse there, as before.

- **The STEP reader closes a spherical cap across its pole the right way
  round.** An exporter writes a cap as its rim circle plus a meridian arc
  from the rim up to the pole, walked up and back down, with no edge at
  the pole. In the chart, the arc's two walks cross the pole row a whole
  turn apart, and the rim joins them a whole period off. The reader's
  pole closing laid no run where a joint stood a period off, so the
  crossing was read the shortest way. Each such cap then measured short:
  - its area by 3% (13.115 against 2πR²(1 − sin v₀) = 13.525);
  - its volume by 0.0395 and 0.0353 for the two caps of the DoomCube 350
    assembly, each checked by closing the cap with its rim's disc
    (5.8486 against πh²(3R − h)/3 = 5.8881, and 5.2007 against 5π/3);
  - booleans through the cap refused ("a face boundary strand dangles").

  Each face is now solved as it is built, by the new
  `ogeom_algo::pole_crossings`. A joint is at a pole where the whole row
  of the chart through it collapses onto the joint vertex within that
  vertex's tolerance. Each crossing there is the shortest way or a turn
  either way more, and the reader keeps the single choice under which the
  wire's unrolled chart closes and encloses area the right way round. It
  lays a degenerate edge only where that choice is not the shortest way.
  `fix_shape`'s pole step is the same solver, now called from there.

  In the real-world check exactly two solids change, DoomCube 709 and 711,
  by the cap corrections: −0.0395 and +2 × 0.0353. A common and a cut of
  709 with a box round its cap now succeed. The common's volume equals
  its fine mesh to 4e-5. The solid's remaining 0.045 gap to its mesh lies
  in the cut, away from the cap; on main it was hidden by the cap's error
  of the opposite sign. The reader counts each face it closed this way
  (`pole-crossed`).

- **`fix_shape` orients a shell's faces along their shared edges, and no
  longer turns a solid's faces against their neighbours.** Two cases.

  A shell whose faces already shared their edges, with one face flagged
  the wrong way, came back unchanged. Sewing had turned that face, but
  its result was kept only when sewing joined an edge. The face is now
  turned back.

  In a solid, each face that the material probe read as facing in was
  turned on its own. On a boundary that passes through itself, the probe
  misreads a few faces. A Voron 2.4 part has a thin zigzag wall whose two
  sides cross each other. `fix_shape` turned 6 of its 270 faces, which
  left 20 edges walked the same way by both their faces. Now the faces
  are first oriented by the topology: each run of faces joined across
  shared edges keeps the orientation most of its faces have. A run then
  turns as a whole, and only where more than half its faces read facing
  in. The boxes that were turned before are still turned. The Voron part
  comes back as it went in. Its findings name the crossing itself, which
  no orientation mends (#397).
- **A saddle crossing a plane along two lines a few thousandths apart
  comes back with both.** Where the march found one of the two crossings,
  the other was never looked for, and where it was, the search cut the
  patch into squares: across the strip between the lines the patch stays
  within `d²` of the plane, so no square could be set aside until it was
  `d` across, along the whole length of both lines, and the search ran out.
  A single crossing now starts the search too, and a patch is halved
  across the way its normals turn most, a strip along the lines parting
  from its neighbours without being cut along its length. Where halving
  does not complete, as about many crossings close together, the search
  runs again quartering. Of 56 such saddles with the lines 1e-3 to 6e-3
  apart, 55 come back whole, against none (#356).

- **A square pocket's rim rounded in one `fillet_edges` call meets at
  every corner (#387).** The flush junction closing the loop was dropped,
  so two of the four inside corners kept their cells: a 4 by 4 pocket at
  r = 0.5 lost 0.882375 mm³ where its sections and four cells come to
  0.906342. A junction closing a loop of joined bands is now taken, and
  the ring's wedge is sewn to itself across that plane. A crossed corner
  whose band's end arc the boolean returned in two pieces, at its
  circle's seam (which follows the part's placement), is read as the one
  arc it is, so a pocket turned in space rounds as it does built square.
- **A fillet band carried on past its cap holds every vertex at least as
  loosely as the edges it bounds.** Sewing the band carried on to its
  wedge, and stating its pcurves' gaps, could widen an edge past its
  vertices' tolerance, leaving a tool solid `check` calls broken for the
  boolean that crossed junctions cut with. Each vertex is now widened to
  its edges' tolerance.
- **A point every fixed ray reaches only through a tangency or along a
  coplanar face is classified exactly (#390).** The classifier cast six
  rays in fixed world directions and, where each met a tangency, an edge
  or a face it ran along, fell back on a mesh's winding number, which
  can read neither in nor out; the boolean then refused. Which rays
  graze depends on the part's placement, so a fillet that rounded a
  drafted pocket square (0.830051763 mm³) fell back to capped corners
  (0.800339681) once the pocket was turned by a placement. A probe is now
  also decided by a segment to a point inside one face (its own nearest
  point, or the middle of one of its mesh's triangles nearest the
  probe), taken only where it meets that face transversally clear of its
  trim and is shown to miss every other face before it.

- **Crossed junctions of `fillet_edges` stand in a loop or chain of the
  request's edges (#342).** A drafted wall's corners in a pocket rim
  were left capped: a crossed corner was taken only between two edges
  meeting no other. Flush junctions are now taken first and crossed ones
  beside them, and a band crossed at both ends is cut in two halfway
  along, each half joined into its own corner and the two joined again
  in the band's own section. A pocket four square with one wall drafted a
  tenth over its depth now loses 0.856067 mm³ in either order, its
  sections, square cells and drafted corners exactly, where it lost
  0.837197 or 0.835998.
- **A rim whose blend's rail passes a sphere's pole within a stride of the
  march rounds.** On a ball of radius 10 drilled parallel to its axis
  with a bore of radius 1.5 at offset 2.499 or 2.501, rounded at r = 1,
  the ball's rail passes about 1.1e-3 from the pole, where its chart turns
  half round within one stride; the blend refused by name. The leg is now
  built on the same sphere charted with its pole on another axis of its
  frame, and the ball's faces on the sphere are rebuilt on that chart for
  the boolean to melt it with, their old seam and pole edges dropped. Both
  round, through `fillet_edge` and `fillet_edges`, to the volumes of the
  same bore turned onto another axis (4048.121147 and 4048.127658) (#365).
- **A boolean through a face whose ring runs out along an edge and
  straight back builds.** Such an edge is a slit: the face lies on both
  sides of it, and its tip is a vertex no other edge of the ring reaches.
  The arrangement read the slit as boundary and refused, "a face boundary
  strand dangles" (Voron 2.4's solid 977, five cylinder faces whose rings
  hold such slits stated as seams, 20 cuts). A slit is now left out of the
  face's boundary, peeled one turnaround at a time, and the pieces of the
  face are rebuilt without it. A seam walked out to a pole and back is not
  a slit: there the chart's whole row is one point, and the ring is left
  as stated (#344).

- **Two sections leaving a node along one tangent are arranged in the
  order they leave it.** The sections of two faces that touch along a line,
  a thin void between them, part at second order, and the arrangement read
  their order from their first chords, which can lean either way. On a
  plane, where sections are sampled to a known sag, two neighbouring
  sections are read again on a wider circle where they stand apart by
  more than that sag; two curves that meet only at their ends keep their
  order out to there. A box cut through solid 977's void, between a
  cylinder and a spline wall within 0.00065 of it, built faces walking
  an edge the same way (#344).

- **A surface crossing another along close straight lines comes back
  with each line whole and traced once.** On a patch crossing a plane
  along a grid of lines 0.1 apart, the search for the crossings cut the
  plane into quarters as often as the patch and ran out of pairs before it
  reached the crossings, so the march's own cuts were kept and stretches
  went missing (0.69 of the lines at chord 1e-4). A plane is now never
  quartered. A stretch joined between two crossings solved its middle
  from one end's touch, where the two surfaces are tangent, and could land
  on another line: the solve now starts beside the middle. An arc stopped
  short of a crossing it heads at is joined to it, and two arcs between the
  same two crossings along the same arms are kept once. Across lines
  0.3, 0.1 and 0.03 apart and chords 1e-4 to 1e-7 no line is missing; at
  0.03 apart with the patch rising 1e-6 off the plane between crossings a
  stretch can still be traced twice (#356).

- **A boolean through a wall far shorter than its face's outline weld
  builds.** A face's outline in its chart was welded end to end within a
  length meant in space but read in the chart, so a wall a micron tall on
  a chart stretched tenfold, or a step 50 nm deep, read its two long sides
  as one and every point of it as outside: its section was dropped and the
  shell did not close (Voron 2.4's screws, 14 cuts). The weld now measures
  each gap on the surface, never joins the two ends of one edge, and
  reaches no further than half the face's shortest edge (#380).

- **The volume of a solid holding a sliver narrower than its vertex, or
  pcurves stated whole turns off a chart that closes without repeating, is
  measured right.** The chart integration refused a loop whose two rulings
  meet only within their vertex a quarter turn apart in a narrow
  cylinder's chart, and read a seam's far side or a rim stated whole
  periods off its spline chart's window where it was stated; the mesh that
  measured those faces instead read the same loops wrong. Voron 151 measured
  62.29 for 63.08, and the screws holding such slivers 0.5624 for 0.5663.
  The loop now closes through the vertex where the chart between the ends
  lies within it, and is folded onto the window by whole periods (#378).

- **A boolean through a sliver of a narrow cylinder whose rulings meet
  within one vertex builds.** On the Voron 2.4 screws, a quarter turn of
  a cylinder 0.005 mm across has two rulings that end at one vertex whose
  radius holds both ends: a quarter turn apart in the chart, a few microns
  in space. Two ends of one vertex are now also welded where the chart
  between them lies within that vertex on the surface (#345).

- **A face whose ring crosses its periodic chart's join is folded beside
  faces set aside too.** The fuse first sets aside the faces the other
  solid cannot reach. A face sharing edges with them was read as stated, so
  its ring did not close, and the fuse fell back to gathering every face.
  Doom 350 solid 709 has a countersink cone whose ring crosses the join. A
  sphere cap elsewhere on that solid refuses when every face is gathered,
  so a cut through the cone refused. The ring is now folded there as
  everywhere else; along a direction the chart only closes in, the fuse
  still gathers every face (#344).

- **A boolean through a face whose edges stop either side of their vertex
  builds.** A file may close two curves on one vertex with each curve's end
  within the vertex's radius but further apart than either radius, as Doom
  350 solid 709 does. The arrangement welded an end only within its own
  radius, so the joint stayed open and the cut refused with "a face boundary
  strand dangles". Two ends that both lie within one vertex's radius are now
  welded within both their reaches (#373).

- **A wall swept from a whole circle or ellipse goes out of IGES as a
  tabulated cylinder (122) over it.** It went out as a spline patch closed
  where the conic starts, which no face can run across. An elliptic wall
  fused with its mirror image and unified, its seam on the mirror plane,
  read back at 125.332 of 125.664; it now reads back exactly (#346).
- **A boolean through a closed spline whose seam states its far side two
  periods off builds.** Both sides read the same points, but the window
  they bounded was two periods wide and the rims met neither column. On a
  chart that closes without repeating, the far side is now read a period
  from the near one, the one copy the chart covers (#344).

- **A boolean through a face whose ring walks a seam both ways and states
  its rims on another branch of the chart builds.** A drum wall or a
  closed spline with its rim stated a period or more off the window its
  seam bounds left the rim's ends where no seam column began, and every
  boolean through it refused: "a face boundary strand dangles". Such a
  ring is now walked onto the seam's branch: each seam occurrence on the
  side that keeps the face on the ring's left, every other edge moved by
  whole periods to start where the ring has reached, and only where the
  ring then closes. A pole joint between a seam's two sides is read on
  the sides walked. Eight box cuts through such faces in a printer
  assembly now build, each cut and its common adding up to the solid
  (#344).
- **`divide_by_area` cuts a converted cone without leaving a sliver.** The
  disc of a frustum or apex cone converted to splines has its middle line
  read from its sampled trims, a hair off the wall's knot line through
  the same rim point: the rim was cut at both, 1.5e-7 apart, leaving an
  edge whose pcurve stood twice its length off, and the divided solid
  would not convert again ("cannot normalize a vector"). A cut line that
  crosses an edge within an end vertex's radius of it now meets that
  vertex (#350).
- **`apply_draft` drafts a fitted wall a boolean cut out of its patch.**
  The cut wall's rims are split in halves, each rebuilt on the section
  marched round the drafted wall, a curve that closes without repeating.
  A half ending on the section's closure read that end at the start of
  the curve and ran back over the other half, and the drafted half of a
  skinned drum broke its chart ("its boundary does not close in its
  surface's chart"). The old edge's direction now picks which end of the
  section the closure end is read at, as it does on a periodic section;
  the half drafts as the whole drum drafted and then cut does (#359).

- **A marched section finds crossings closer together than its stride.**
  Where two crossings of a section stood within one step of the march, a
  third could hide in a step unseen, and the pair was cut as traced. Both
  charts are now cut into patches, a pair set aside where their boxes do
  not meet or no normal of one can be parallel to a normal of the other,
  and every crossing is solved from the pairs left at the chord's size;
  where that search completes, each crossing's arms are joined to the next
  crossing along the section or followed out, a stretch too shallow to
  march held where it stays within the confusion distance of both
  surfaces. On a saddle crossing a plane along three lines, of 640 sizes
  and shapes, 455 came back cut as traced and wrong, and now 27 do. Where
  the surfaces part so slowly that the search does not complete within
  the points a branch may hold, as on a grid of lines crossing a tenth
  apart, the pair is still cut as traced (#356).

- **Two fillets of one `fillet_edges` call meet at an inside corner off
  square in plan, or between two walls that both lean (#342).** Each
  band carried on is widened past its wall, and the two widenings met
  past both walls' planes, in the open corner, so the corner was refused
  and the bands capped at the vertex. Inside a 2.4 mm box round the
  corner, r = 0.5 and one wall leaning at a slope of 0.1, that left
  0.0136 to 0.0445 mm³ standing at plan angles of 120° to 170°. The open
  corner is now taken out of the two widenings, a wedge's own run-out
  past an oblique far end no longer counts against the corner. Every
  case measured matches the material above both bands to within 1e-6.
- **A boolean folds a ring onto its pole's branch of a closed chart at
  any chart scale (#349).** Beside a pole or apex edge the fold matched a
  ring's ends within a fixed chart distance of 0.001, which a chart
  closing over less than that period reads as no jump at all: a cone
  sector on a patch closing in 0.0005 of `u`, its ruling on the closure
  stated on the other line from the apex's, failed with "a face boundary
  strand dangles". The ends now meet within the edges' radii, the reach
  the arrangement's weld has there, where the chart shrinks along the
  pole's line.

- **A boolean through a cone or ball face whose ring reaches the apex or a
  pole with no edge there builds.** Two rulings or meridians meeting at
  the apex stand apart in the chart, at the columns they arrive along,
  and with no degenerate edge nothing drew the apex row between them: the
  boolean refused with "a face boundary strand dangles". Such a joint now
  closes along the row as a degenerate edge would, where both ends lift to
  the joint's vertex, the surface collapses that row to it within the
  vertex's radius, the row is shorter than a period, and the face lies on
  its left; any other joint still refuses. Eighteen of the twenty box cuts through
  such faces in two printer assemblies now build, each cut and its common
  adding up to the solid (#343).
- **A boolean reads a period jump on a surface smaller than its edges'
  radii as a jump (#345).** The fold took a joint whose ends lie a period
  apart in the chart but close on the surface for one across an apex, and
  a ring on a cylinder smaller round than its edges' radii was left
  unfolded. A joint now stands across an apex only where the surface's
  derivative along the chart between its ends vanishes, and a whole
  period across any other joint is a jump of branch.

- **`apply_draft` gives a drafted planar wall a window that holds it.**
  The leaned plane's chart has its origin on the hinge and its `u` along
  it, but its window was the old plane's numbers grown: on a prism's wall
  the face stood at u in [-8, 0] and the window at [-0.8, 8.8]. The old
  window is now turned with the face and read in the leaned plane's own
  chart. Offsetting the drafted prism was invalid ("a pcurve that leaves
  its curve by 0.517"), its walls meeting the leaned one beside its
  trim; it is now valid and measures its grown planes' volume exactly
  (#332).
- **`to_bezier` leaves a cylinder's or a torus's old seam a plain trim of
  each patch it bounds.** Cut into patches, a closed surface's seam bounds
  the first and the last patch once each, and neither patch closes across
  it, but the edge kept a seam on each with both columns, one a whole
  turn outside the patch. The patch now keeps only the column its
  boundary walks, as a plain pcurve (#331).

- **`check` reports a seam on a surface that does not close across it.**
  A seam's two sides stand a period apart across a direction its surface
  closes in; on a surface that does not, the side a face does not walk
  names another stretch of the surface, which a face reading the edge by
  its orientation can take (#331).

- **A solid fused with its mirror image, or placed by a mirror, reads
  back from IGES, and from STEP with every face facing out.**
  - A ball went out of IGES with its centre and radius alone, and a
    reader takes such a ball's axis as z: one about any other axis came
    back with its poles, and its seam with them, where its edges were
    not. It now goes out with its axis and the direction its angle
    starts from.
  - A wall swept round a closed curve goes out of IGES as a spline
    patch closed where the curve starts; other systems' STEP files hold
    such patches too. An edge on that closing line was read on the wrong
    side of it, and a wall lost a strip. Such an edge now meets the edge
    before it across the line, or, opening its loop, the loop's last.
  - A mirror turns every face's chart over. A face written as a spline
    patch or a swept surface, chart and all, was flagged as facing the
    way its orientation says, against the face's own side; in IGES a
    plane or other analytic surface restated about its axis kept loops
    running about the chart's normal rather than the file's. Both
    writers now flag such a face turned, including a STEP face on a
    trimmed or offset analytic surface, IGES writes its loops about the
    file's normal, and STEP bounds a mirrored face with its interior on
    their left.
  - STEP states an offset of a mirrored analytic surface along the
    normal its basis is read with, which a mirror turned: a mirrored
    drum wall stated as an offset had come back moved twice its
    distance.

  A ball about x fused with its mirror image read back from IGES at
  496.7 of 523.6, a mirrored elliptic prism from STEP with its wall
  facing into the material; both now read back as they were. Eleven
  solids of two printer assemblies read from STEP change volume by 0.1
  to 0.3 %, each a face on a closed spline patch that had lost a strip
  (#336).

- **A boolean brings a hole of one edge onto its outer ring's branch of
  a closed chart (#348).** A hole whose single closed edge states its
  pcurve a period away from the face's outer ring was left there, and the
  boolean read the face without its hole: a quarter wall's sheet cut
  through such a hole kept the hole's area as face. It is now moved the
  whole periods that put it inside the outer ring's extent, as a hole of
  several edges is.

- **A boolean through a face whose ring folds across its closed chart's
  join builds (#333).** An edge on the closure of a closed chart may state
  its pcurve on the other branch from the edges it meets: a quarter-round
  wall's straight edge at `u = 0` while its arcs end at `u = 2 pi`. The
  ring is one closed curve on the surface, and the mesher folds it back,
  but the boolean read each pcurve where it was stated and failed with "a
  face boundary strand dangles". Each edge is now moved by whole periods
  onto the branch that continues its ring, anchored on the ring's pole or
  apex edge where it has one, and each hole onto its outer ring's branch.
  A cone sector, a sphere lune and a drilled wall with any edge restated
  a period off now cut and intersect exactly as when stated on one
  branch, and 62 of 140 box cuts through folded faces of two real-world
  assemblies that failed now build. Rings holding a seam are left as
  stated (#344), as are apex joints with no degenerate edge (#343) and
  tiny cylinders whose period jumps read as apex joints (#345).

- **`to_nurbs` and `to_bezier` convert solids divided off their knots
  exactly.** A drum, ball, ring or frustum cut by `divide_by_angle` into
  pieces that do not end on its spline form's knots converted wrong:
  - Cutting an edge cut every description at the same fraction of its
    own range, though a pcurve paced its own way reaches that fraction
    elsewhere: `to_bezier` of the drum and the frustum left trims 0.012
    and 0.006 off their curves. Each description is now cut where it
    passes the cut's point.
  - The divided pieces' iso trims, stated as linear splines, were refitted
    rather than restated as iso lines, and a fitted end on a pole or a
    closure took any column: the ball read a face both ways, the ring had
    trims 18.4 off. Any trim whose points share one coordinate is now an
    iso line, placed where a polished projection puts it rather than to
    the confusion distance along the surface (a thin tube's volume was
    1.8e-5 off); a fit's pole samples and pole corners take the column
    they arrive along; and a seam on an iso line is that line exactly.
  - A pole row spanned the whole chart whatever part of it the face held.
  Each conversion now checks valid, draws closed and keeps the solid's
  volume; `to_bezier` of the ball and the ring, refused before, builds
  (#339).

- **A lone edge is held to one ball-fit check through either entry
  point.** `fillet_edge` let a ball's contact run off a host face
  anywhere, as long as it fit somewhere along the edge, while
  `fillet_edges` refused any station where it ran off, even for one
  edge. On a wall that carries a boss with its underside 0.2 over the
  floor, `fillet_edge` cut a radius-0.5 band along the whole edge. There
  the ball touched the wall's plane 0.3 inside the material, so the band
  grooved the block under the boss. On a wall that narrows to a tip,
  `fillet_edges` refused a band that `fillet_edge` cut off correctly on
  the faces the edge runs into. A station may now let the contact run
  off its face only where the contact stands outside the solid, through
  either entry point for one edge. The boss is refused by both, and the
  tip is rounded by both to the closed form. A request of two or more
  edges keeps the check at every station (#316).

- **A solid fused with its own mirror image reads back from STEP as the
  solid it was.** A mirror leaves a circle's, ellipse's, hyperbola's or
  parabola's frame left-handed, and STEP states a frame by its axis and
  its `x`, the `y` taken as their cross product: the writer handed a
  reader the mirrored conic's `-y`, so an arc ran the other way round
  between its vertices. A conic now goes out with the axis `x × y`, the
  same curve at every parameter. And a seam walked into a pole and back
  out of it was read with both its sides on one column, the wire walked
  on from the wrong one; a seam's second use is now placed a period from
  its first wherever the walk lands them together, in the STEP and the
  IGES reader alike. A ball, a torus and an elliptic prism each fused
  with its mirror image read back with a volume of 0, and now at their
  own, unified or not (#323).

- **Marched sections through a point where two surfaces touch come back
  as the curves' pieces, not round its corners, where the touches stand
  apart.** Where two surfaces touch at a point and bend apart differently
  across it, the section crosses itself there along two arms. A walk
  coming in along one arm could turn onto the other, or wander there, and
  the section came back running round corners, over itself and short of
  the touch; an arm no walk took was lost. The two arms are now solved at
  each touch, from the difference of the surfaces' second fundamental
  forms, and a branch is cut at a touch unless it comes in along one arm
  and goes on along the same one, however far off the touch a step
  turning the corner runs. An arm no section runs along is followed out
  from the touch, and a stretch traced twice is kept once. A point where
  the surfaces only pass within the confusion distance of each other is
  no touch to cut at. Where two touches stand within a stride of the
  march of each other, a touch between them may go unfound, and the pair
  is cut as before. Equal cylinders crossed square come back as four
  halves at every radius tried, a sphere cut by a cylinder through its
  centre (Viviani's curve) no longer comes back with a whole lobe twice,
  and on a saddle crossing a plane along three lines 38 of the 493 sizes
  and shapes that came back wrong now do not (#327; touches a stride
  apart: #356).

- **A part fused with its mirror image builds where a face's chart closes
  in `v` along its boundary.** The bake a reflecting placement forces
  places an edge on the closure of a closed chart beside its face, but did
  so only for a column on a chart closing in `u`: a row on a chart closing
  in `v` (a torus wall's edge on the circle its chart closes along) took
  whichever side the projection answered, and a column ending on that
  closure was not given the iso line at all, its fallback pcurve standing
  a period away from the rows. A notched profile revolved a quarter turn
  and fused with its mirror failed with "a face boundary strand dangles";
  it now builds, at twice Pappus's volume for each half (#319).
- **`fillet_edge` runs a round out onto a ramp, as `fillet_edges` does.**
  The run-out of #304 was reached only through `fillet_edges`. Rounded
  alone, the same edge still stopped square at its end, leaving a cap
  and a step. On the #304 carriage it removed 57.79 mm³ against 63.62.
  The two entry points now give the same rounding on one edge, except at
  an end against a neighbouring blend: `fillet_edge` stops flush there,
  which is the state `round_vertex` closes. A shape that bounds no
  inside (an open shell) keeps the flush cap rather than refusing.
- **`round_vertex` closes a corner that a band ran out through.**
  Rounded one at a time, an edge's band now runs out past a vertex where
  the face there leans out past its cap. The band of the edge it crosses
  then stops flush against it, short of the corner's block, and the
  corner tool left slivers there. Where one ball touches every face, the
  tool now carries that band back to its block first. It does so only
  for a band that stands against another of the corner's bands; a band
  stopped by a notch is left alone. Across 80 cases (every corner of
  five sheared blocks, with the edges rounded in both orders), three
  fillets one at a time and `round_vertex` now give the solid
  `fillet_edges` gives the three edges in one call. Before, 16 matched
  and 10 refused (#313).

- **A section stalling a long stride short of a touch is carried onto
  it.** Where a marched section runs straight into a point where the two
  surfaces touch, the walk takes long strides, and cutting the stall's
  wander away left its end a whole stride short. The touch was solved
  from that end but taken only within 2400 chords of it, so a longer
  stride left the section short: 2.7 µm on two swept walls meeting at
  their apexes. A touch solved from a stalled end is now that end's at
  any distance, and the end is followed out to it along the section.
  Where another touch lies as near and the end does not run into it, the
  end is carried onto the one solved from it instead: on a saddle
  crossing a plane along three lines, an arm ending 5.2e-2 short of the
  crossing it runs into now reaches it (#321).

- **`offset_shape` builds a closed offset of a face whose seam is two
  edges running opposite ways from a vertex.** A seam edge running
  against its meridian's parameter (from a ball's equator down to its
  south pole) was rebuilt from the angle its new ends stand at, taken in
  the first turn: three quarters of a turn the long way round, through
  rows the sphere's chart does not have. The offset passed `check`,
  measured the exact volume and meshed open. The rebuilt seam now starts
  in the turn the old one's rows stand in and runs the way they ran.
  Every face the offset rebuilds is held to closing in its surface's
  chart, and the offset fails by name where one does not (#325).

- **`check` reports a face whose boundary jumps across rows its chart
  does not have.** Where a pcurve ran past a sphere's pole, the gap to
  the next edge's start was skipped because the surface refused to
  evaluate it, and the face passed while meshing open. The gap is now
  measured from where the surface stops in that direction (#325).

- **Two fillets of one `fillet_edges` call meet at a drafted corner along
  their own intersection (#326).** Where the two walls at a junction
  differ (a drafted wall beside a square one), no plane bounds the
  meeting and each band was capped as without a junction. At an inside
  corner the cell above both bands kept its material, 0.009 to 0.020 mm³
  at r = 0.5. Now, where both creases are straight, each band is carried
  on with what stands past its wall's plane, and the two are intersected
  in one boolean: each wedge is kept to what stands above the other band
  carried on, and the corner is bounded by the bands' own intersection,
  where the material holds all of it. At an outside corner the band
  built second stopped at its vertex where the solid ran on under it,
  so one order came short (0.000493 mm³ at a 0.5 mm draft); each band
  now runs out as it does alone and the cuts take their union, in either
  order. A crossed corner is taken only between two edges that meet no
  other edge of the request. Every case measured matches a semi-analytic
  ideal; a corner the construction does not close rounds with its flush
  junctions alone, as before.

### Performance

- **A whole ball or ring is read without meshing it before a boolean.**
  Whether a solid of one closed face is inside out, a half space, was
  read from the sign of its mesh's volume; and a face's box meshed the
  face to decide which of a sphere's or a torus's closed-form extremes it
  holds. A face covering the whole of its sphere or torus (every edge a
  seam walked once each way, or a pole) is now read from its surface: it
  holds every extreme, and faces in or out as its surface's normal and its
  orientation say, whatever its placement. New `covers_whole_surface`
  answers whether a face is one. Every boolean answer is unchanged; a
  whole ball's or ring's box is now the closed-form one, where the mesh's
  points could stand a few units in the last place past it. Fusing two
  balls takes 1.35 ms instead of 6.0 ms, cutting a ball from a box 0.91 ms
  instead of 3.3 ms, and fusing two rings 49 ms instead of 59 ms.
- **A plane resting beside a torus, and a ball resting against a drum's
  wall, are answered in closed form, not marched (#439).** A plane not
  square to a torus's axis, and a sphere off a cylinder's axis, had no
  closed form, so the general intersector marched them even where they
  never cross. Where they rest against each other at one point it found
  no crossing after seconds of seeding and answered `Apart`. Now
  `surface_surface` measures the plane against the torus's reach along
  its normal (`R` times the sine of the normal's angle to the axis, plus
  `r`) and the sphere's reach from the axis against the cylinder's
  radius. Clear of it, the pair is `Apart`; within the confusion distance
  of it, `Touching` at the one point where it rests. Where they cross,
  the pair is marched as before. Every boolean and fillet answer is
  unchanged: the boolean took the march's `Apart` and takes a touch the
  same way. Rounding the top loop of a 10 mm block whose vertical edges
  are rounded takes 0.54 s instead of 26.0 s, and rounding every edge of a
  D-shaped prism 0.87 s instead of 47.0 s, on one thread.

## [0.11.5] - 2026-10-19

A patch release of fixes, with two changes in behaviour. A part fused
with its own mirror image builds (#302). `fillet_edges` runs a round on
past its edge's end down a ramp (#304), rounds a chain carried on along
its tangent continuations, under an overhanging roof and round a short
rim (#306), meets two fillets of one call where their creases turn
rather than capping each, and stops a marched round where its crease
turns back into the material. `unify_same_domain` merges faces on one
curved carrier as well as one plane (#305), and a boolean passes through
the faces a single-lump tool cannot reach.

### Changed

- **A boolean passes through the faces a single-lump tool cannot reach,
  even inside that tool's box.** A face whose box missed the other
  solid's box was set aside and passed through untouched; one inside a
  lump's box was always gathered, split and rebuilt. Where the other
  solid is a single lump, a face whose box misses every one of that
  lump's faces now lies wholly in or out of it, and one point of it
  classified outside sets it aside too. Two fillet bands joined round a
  plate's corner span the plate in one lump while touching only its
  edges; on a plate drilled with 256 holes they cut in 142 ms where they
  took 197.

- **`unify_same_domain` merges curved faces as well as flat ones.** Faces
  sharing an edge on one cylinder, cone, sphere or torus (one axis, one
  size) or on one extrusion (one basis curve swept along one line) merge
  into one face on the first face's surface, whatever frame each stored
  its surface in. The others' pcurves are carried into that chart by the
  exact map between the two charts, each parameter run either way and
  shifted, by whole periods where it repeats; a pair whose charts no such
  map relates stays split. Halves that close round their axis, such as
  the two halves of a bore, keep one shared edge as the merged face's
  seam. A part fused with its own mirror image no longer keeps every
  curved face the mirror plane crosses split in two: on the carriage of
  the report, 46 cylinder faces become 38 and 20 extrusion faces 16, the
  volume unchanged (#305).

### Fixed

- **A marched round stops where its crease turns back into the material.**
  A crease along a wall swept from a spline, ending where the wall turns
  on into a round boss 1.3° back behind it, was refused: "the marched
  band's end where its crease runs out on another face does not melt with
  that face". The band ran out past the vertex as through an end wall,
  but the crease carried on runs into the material there, and the boss
  curving away let the ball out only 3.5 mm on, along the boss's own rim,
  in a grazing contact the melt could not resolve. A band now stops at
  its vertex where every face there stands behind the crease carried a
  radius on, by more than the face's own tolerance. The carriage's lower
  rib crease rounds alone to its closed form, 2.741 mm³, and the
  carriage's ten-edge fillet comes to 11686.820 mm³ against the source's
  11686.827, `check` clean (#306).
- **Two fillets of one `fillet_edges` call meet where their creases turn.**
  Two convex creases sharing a face and turning at a vertex were rounded
  each to its own cap there. Across an inside corner (an L-shaped top's
  re-entrant corner) both caps stood and the corner cell both rounds hold
  kept its material, short by twice the corner's moment about its wall
  times tan(t/2), t the turn: 0.0959 mm³ at a square corner of radius 1.
  Across a shallow outside turn, a spline wall turning 3.6° into a round
  boss's rim, the marched band ran out along the rim, crossed its crease
  2 mm on and cut the rim's round short there, which kept 0.185 mm³ of the
  rim's corner on the carriage of #306. Carried on past the vertex, each
  band is a straight band of its end section, and two of one radius cross
  in the plane bisecting their centre lines: across an inside turn both now
  run on to that plane, across an outside turn whose bands' sections there
  agree to a marched band's fit both stop in it, and the two are joined
  across it into one wedge. That meeting is taken only where it is exact,
  the reflection in the plane carrying one band's end section onto the
  other's, and wherever it is not (a drafted wall beside a square one, a
  turn folding back past the solid, a wider outside turn of curved bands)
  or does not close, the request rounds as it did before. Each band's
  faces stay credited to its own edge. On the
  carriage's chain the five junctions agree with the source's shape to
  within 0.002 mm³ each, and two pieces of one torus whose edges include a
  section of it meet: the boolean carries an edge between two charts of
  one torus as it does between two of one cylinder (#306).
- **A solid fuses with its own mirror image.** Fusing a part with its
  mirror across a plane it stands on failed with "arrangement left no
  piece of the face" or "a face boundary strand dangles". A mirrored
  circle stands on a left-handed frame, and its closed-form pcurve on a
  cylinder, cone, sphere or torus ran the other way round the axis, so
  the image's arcs split the part's faces on its own side. Baking the
  image, a fitted edge ending on a vertex on a cylinder's seam line was
  pinned to the vertex's image on the other side of the seam, and an
  edge lying along the seam line was put on the chart's edge away from
  its face. Each is now placed where its curve and its face are (#302).
- **A round runs out past its edge's end onto a ramp.** Where a straight
  edge's floor gave way at its end to a face rising obliquely, the wall
  carrying on above it, `fillet_edges` stopped the round square at the
  end vertex: a cap and a step were left, and the corner past the end
  kept its material. The round's cylinder now carries on along the face
  it stays tangent to until the face that follows has cut its whole
  section away, and stops on that face with no cap. A step or a ledge
  past the end, material no ball on the edge touches, still stops the
  round square (#304).
- **A fillet chain rounds on along its tangent continuations.** Two
  things refused a chain carried on along tangent-continuous edges. A
  straight crease running up under a face that overhangs its corner (a
  rib under a slanted roof) was refused, its contacts standing behind
  the roof, and rounded alone its band cut into the roof. The band is now
  trimmed by the roof's plane and stops on it. A short arc of a rim was
  taken for a piece of a rim that closes whenever the rim stood opposite
  its middle, and its round ran the whole turn, through the part where
  the rim stops. A rim now runs the whole turn only where the edges on
  its circle that bound its wall close it (#306).

## [0.11.4] - 2026-10-18

A patch release of fixes. A shape imported with an extrusion face (from
`.brep` or STEP) saves to the native format and reads back (#295).
Rounds and chamfers meeting at an inside corner whose vertex stands a
hair off one edge's wall meet flush rather than refusing (#298), and a
section running into a point where two curved faces touch at their
ends reaches it, so a pad standing on a swept wall's end fuses (#297).

### Fixed

- **The native format reads back an extrusion whose window starts below
  zero.** `native::write` writes an extrusion's window as it stands, and
  `.brep` and STEP imports give extrusions windows reaching both ways
  from their curve, but `native::read` refused any window not starting
  at zero: a shape imported with an extrusion face could not be saved
  and opened again. The reader now takes any finite, increasing window
  (#295).
- **Rounds and chamfers meeting at an inside corner whose vertex stands
  off one edge's wall meet flush.** A straight edge's cap was slid back
  only to its vertex, so where the vertex stood at one edge's line end,
  off the wall that edge runs into by less than its tolerance, the cap
  stood past the wall and cut a sliver off the neighbouring blend; on an
  imported part the boolean's arrangement could not place it and
  `fillet_edges` refused. The cap now also slides back to where the line
  meets the next straight edge of either of its faces at that vertex,
  whichever place lies furthest back, each within the vertex's tolerance
  (#298).
- **A section running into a point where two curved faces touch at
  their ends reaches it.** Where two walls swept from quadratic profiles
  touch at the apexes their profiles end on (a pad's rounded side
  standing on a wall's end), the marched section between them stalls
  short of the touch and is carried onto it once the touch is solved.
  The solve's Jacobian was differenced through the clamp to each face's
  chart, so a step past the bound the touch sits on read as no change
  and the solve gave up: the section stopped 0.2 µm short and the fuse
  failed with "edge 1 ends where edge 2 does not begin". The
  differences now step back into the chart at a bound (#297).

## [0.11.3] - 2026-10-18

A patch release of fixes. `check` names a solid whose boundary passes
through itself, settled on the exact surfaces, and no longer lets one
probe decide which way a face or shell faces (#276). `refine` merges
every split wall (#281), a skinned loft between sections wound opposite
ways or closing to a point no longer folds (#269, #282), surface walkers
close a loop only through a step they solved (#273, #278), and `.brep`
files saved by other applications read as the solids they hold (#283).
Rounds on walls swept from fitted profiles and at vertices their edge
overruns now build (#287). Sweeps, marched rounds and silhouettes state
how far they stray between samples. Mass properties and shape healing
probe each face's orientation once, about 10% faster.

### Fixed

- **A straight edge's cap moves to its vertex only by the line's measured
  overrun.** A round's or a chamfer's cap was slid to the nearest vertex
  in either direction by up to the vertex's tolerance, so a line falling
  short of a tolerant vertex was stretched, and a short edge could pair
  both ends with one vertex. The cap now moves back by exactly how far
  the line runs past its own vertex, and only where that vertex stands
  within its tolerance of the end and the two vertices' tolerances do
  not between them cover the edge.
- **Only a marched round that meets an analytic one at a vertex is built
  ahead of the order asked.** The marched band was put first in every
  set mixing the two kinds, and an edge was taken as marched by its
  curve alone, so a straight edge rounded by the march was not, and a
  spline far from any analytic edge was moved anyway. The order now reads
  the construction each edge is actually rounded by.
- **A spline crease's continuation is never cut short by its refit.** A
  crease continued past its ends is put back on its hosts only where
  both hosts are found to meet along the whole continuation; where a host
  ends first, the continuation stays as itself rather than being refitted
  over the stretch walked, which could end it a fraction of a radius past
  the crease. A wall swept from a fitted profile is continued twice as
  far as the crease, so the crease's continuation meets it to its end.
  An error in that refit is returned rather than taken as no refit.
- **A round on a wall swept from a fitted profile no longer refuses, and
  two rounds meeting at a vertex their edge overruns are capped there
  (#287).** A community printer carriage, rebuilt by printCAD from its
  FreeCAD file, failed to round: a rib's top edges — a straight run, a run
  along a wall swept from a quadratic B-spline profile, and a round end —
  meet at corners of ten and three and a half degrees. The swept wall
  ended where its profile does, so no ball seated past the crease's end
  and the run-out's cap plane cut its host nowhere; and the crease along
  the wall's top, a fitted cubic that does not continue the wall, strayed
  off both hosts by tens of microns where the run-out's cap stands. The
  wall's profile is now continued the same way past its ends and the wall
  swept over it, and a crease continued as itself past its ends is
  refitted onto where its two hosts actually meet, so a run-out's caps
  and crease stand on them. Where a marched band and an analytic one meet
  at a vertex, the marched band is now built first, against the faces as
  they stand, so its run-out melts with untouched hosts rather than
  landing along the other band's rail. On the hex-nut variant, two top
  edges meeting at an inside corner of the top face leave one edge's line
  running 0.24 µm past the vertex it ends at, within the vertex's
  tolerance: a round capped where its edge's curve ends cut the other
  round 0.24 µm deep, a sliver the boolean's arrangement could not place.
  A straight edge's cap is now slid along its line to the vertex it ends
  at, so the two meet flush.
- **`check` names a solid whose boundary passes through itself (#276).**
  Nothing in `check` looked for faces crossing each other, and which way
  a face faces was read at one place, the middle of its largest mesh
  triangle. A wall folded through itself passed as valid wherever that
  place missed the fold: a skinned loft twisted on to a point at most of
  the apex positions tried, and a prism over a profile that loops through
  itself always. Each solid's faces are now meshed, the triangles that
  cross each other found by their boxes, and every crossing settled on
  the surfaces themselves: a point on both where they meet at an angle,
  inside both faces' trims and clear of their boundaries. A face crossing
  itself or another face of the solid is named broken; faces meeting
  along the edges they share, or meshes crossing where the surfaces only
  come close, are not. On the real-world assemblies DoomCube's solid 173
  is now named for twelve crossings of cylinder and plane faces, each
  confirmed by a winding number of -1 beside it, and Voron's solid 907,
  already broken, for its four.

- **One probe no longer decides which way a face or a shell faces
  (#276).** `make_solid` turned a closed shell over on the first face
  probe that read, at the middle of that face's largest mesh triangle; a
  skinned loft folded near the point it closes to had that probe land in
  the fold, and came out inside out. It now turns a shell only where
  three places spread over its faces agree, and otherwise by the volume
  the shell encloses. `check` reads each face at two places, its largest
  triangle and the large triangle farthest from it, and names a face
  that reads out of the material at one and into it at the other rather
  than settling it by whichever half the probe landed on: the top of a
  box half covered by a second box as another shell of one solid now
  reads both ways under every placement. A probe can misread on a sound
  solid, so a face read both ways is named only where the boundary was
  found crossing itself, where the place a single probe asked reads it
  facing in, or where the solid's winding number bears out the inward
  reading.

- **A sampled fit states no less than its miss at a crease** (#270).
  `fit_surface_sampled` read each span at its quarter points and climbed
  at most sixteen peaks, and `fit_curve_sampled`, `fit_curve_2d_sampled`
  and `fit_trace_sampled` read eighth points the same way. Where the
  geometry creases (a polygon read by arc length) the miss peaks sharply
  at the corner, between those points: a triangle twisted up a prism
  stated 9.5e-4 against a tolerance of 1e-3 and stood 1.18e-3 off at a
  corner, and a closed polygon of 25 corners stated 9.98e-4 and stood
  1.11e-3 off. A fit those points pass is now read again at eight points
  of every knot span of its own (a surface along its sample lines and
  midway between its rows), every peak climbed that could reach the
  widest, a surface's widest climbed across its lines too, and refined
  where that misses; the error it states is the one read so. Sweeps and
  lofts no longer re-read their skins' rows themselves: the fit reads
  them at least as closely.

- **A skinned loft between sections wound opposite ways no longer passes
  through itself (#282).** The skin joins its sections' points at equal
  fractions of the way round, so two circles sketched on facing planes
  (the two end faces of a part, each sketch facing out of its face), which
  wind opposite ways, had each point joined to one across the loop: the
  wall folded through the middle and the solid held no volume, and a
  subtractive loft between them failed with "the kept pieces did not
  close into a shell". A section winding against the one before it is now
  walked backwards, from its point that rows up best with that section.
  Lofts to a point and sections moved along a path are unchanged.

- **`unify_same_domain` merges every split wall it meets and no longer
  refuses a valid solid (#281, #282).** Three faults, each met on parts drawn
  elsewhere and padded or pocketed again. Faces on one plane were grouped
  by the normal their plane stored, so a pad whose profile ran the other
  way round, whose walls store their planes facing in and are used
  reversed, kept every wall it shared with the part split (an open box
  came back with fifteen faces, not eleven); faces are now grouped by
  the side they face. The merged boundary was chained into a single
  wire, so a merged face with a hole in it (a pocket's floor, a boss's
  rim across the split) refused the whole refine with "the edges do not
  form one connected path"; the boundary is now walked as the merged
  faces walk it, into an outer ring and one ring per hole, each wound as
  the material lies. And the merged face kept the first face's plane
  window, which a bounded plane (a prism's wall) sized to that face
  alone, so the rest of the merged face read as off it and the solid's
  orientation probe called its neighbours turned in; the window now grows
  to hold every merged face's. A group whose boundary still does not
  chain into rings (two rings through one vertex, an edge a third face
  also bounds) is left split rather than failing the refine.
  `ogeom_algo::edges_turn`, the signed turn of a walk of edges about an
  axis, is now public.

- **A shape saved in the `.brep` text format by another application reads
  as the solid it is (#283).** Project archives keep each part's shape in
  that format as the shape alone, and the reader got them wrong in four
  ways: it required the content-type line those files leave out; it read a
  seam's second chart number and its continuity, written with no space
  between them (`2CN`), as one word; it composed a chain of placements in
  reverse; and it dropped the placement an edge's chart is bound to, so a
  rim used at two placements (the two ends of a padded profile) gave both
  ends the first one's chart and faces came out empty or wrong. A
  10 x 5 mm padded disc measured 520 mm^3 instead of 1570.8. Reading those
  files also needed: placements kept as chains of the file's elementary
  ones, so shared vertices and edges stay shared; wires walked end to end
  (the file lists a wire's edges as a set); an edge's start vertex held
  first; charts derived for edges on planes, which the format leaves out;
  surfaces on left-handed seats rebuilt right-handed with their charts and
  faces turned over; a cone's charts read along the ruling the format
  measures them on; extrusions reaching both ways; and a face whose
  placement cancels its wire's rebuilt where it stands. `read_numbered`
  returns the edges and faces in the file's own numbering, which differs
  from the walked wires' order for edges. The writer now writes placement
  records the format (and the reader) accept, the file's own placement,
  curve placements, each face's use of a chart at its own placement, and
  cone charts along the ruling. Checked on 262 shapes from saved projects
  against an independent reading of the same files: every face, edge and
  solid count agrees, as does the edge and face numbering, and 176 of the
  178 solids measure their volume to a millionth (the other two to 5e-5,
  the share of the measurement taken from a mesh).

- **A skinned loft closing to a point no longer folds its wall through
  itself near the point (#269).** The skin through the sections was
  carried on to the point, and its rows bent to meet it: past a section
  turned against the ones before it (twisted squares and rectangles) they
  shrank through a figure eight and came out wound the other way, so the
  wall near the point faced into the material. Depending on where the
  orientation probe landed the solid came back either visibly turned in
  (`check` rejected it) or quietly self-intersecting. The smooth skin is
  now measured densely between the last section and the point, and kept
  only where it is unfolded there (a graph over the angle round, and the
  height along, the line from the section's middle to the point); where
  it would fold, the last section closes on the point by the straight
  cone over it, which cannot fold, with an edge where the two meet.
  Gentle lofts and pipes to a point stay smooth as before. The patch
  closing on a point is also turned outward by its signed volume over the
  whole patch rather than its normal at one point.

- **A walked curve closes only through a step it solved (#273).** The
  shared walker behind surface intersections, silhouettes and marched
  rounds declared a loop closed whenever a step landed within a step of
  its start, and joined the start on without solving that last chord.
  Across the seam of a closed but non-periodic surface the curve can
  break (a spline sphere's silhouette by about 0.004), and the walk
  joined it anyway: the closing chord missed the curve by 3e-3 to 8e-3 and
  the loop claimed a closed curve that is not there. The start must now
  lie within the step being taken, the condition solved on the plane
  through it must land on it, and the condition solved halfway must keep
  to the chord; otherwise the walk ends open at the seam.

- **A tangential contact closes only through a step it solved (#278).**
  The walker tracing where two surfaces touch had its own closure test:
  a loop closed whenever a step landed within a step of its start, and
  the start was joined on unsolved. Halves cut at the seam of a closed
  but non-periodic patch were joined whenever their ends stood within
  twice the longest step of each other. Where the contact breaks at the
  seam (by 0.004 for a strip resting on a drum along a climbing curve)
  both claimed a closed contact that is not there, the closing chord
  missing it by half the break. The step onto the start is now solved
  and held to the chord as in the shared walker, and halves cut at a
  seam are joined only where their ends, each landed on the seam, meet
  within the chord; otherwise the contact ends open at the seam.

- **A curve whose backward half closes comes back as that loop (#278).**
  Following a curve both ways (the shared walker, the intersection
  tracer and the tangential contact tracer), a forward half that stopped
  short of closing (out of points, or stalled where the other way
  passed) was joined to a backward half that had closed, so the result
  ran round the loop and on past the start a second time, reported
  truncated or as leaving the domain. A closed backward half is now the
  answer, a loop, turned to run the way the forward half set off.

- **A part placed in an assembly is checked as it is on its own (#275).**
  Which way a face faces was probed at the middle of its largest mesh
  triangle, ranked by area from the placed mesh. A placement rounds the
  positions, and triangles of one area (a cylinder strip's, a rectangle's
  two) swapped places; on a face that reads out over part of its extent
  and in over the rest, the rotation placing a part in an assembly turned
  one more face "inside out" than the part read alone, and repair turned
  it back and forth without settling. The face is now meshed in its own
  frame for the probe, so a rigid placement moves the probe points with
  the material and changes no answer: the Voron 2.4r2 body
  `molex_0471510001:1` reads 6 reversed faces placed as alone, and none
  after repair.

- **A marched round's line of contact states how far it strays between
  the march's stations.** Between two faces that share no edge, the line
  along which the ball touches a face was fitted through the march's
  touch points and stated to its target there only. Where the touch point
  runs unevenly along the guide (a floor against a twisted, rippled
  patch) the fit strayed from the ball's true contact between the
  stations by ten times its stated tolerance and more. The ball
  is now seated exactly at the quarters between every two stations, the
  line is read against those sections both ways and fitted again with
  any it misses, and the edge states the farthest it was read.

- **A skinned loft holds its sections round tight corners.** A section is
  read along a polyline, and the polyline was spaced by one sag measured
  over a first pass of 64 points and spread evenly. Round a corner tighter
  than that pass the chords cut across it: a pipe closing a square to a
  point along an arc had its skin 7.2e-3 inside its sections' corners at a
  tolerance of 1e-3, and a loft through rounded squares had its end rings
  1.9 times the tolerance off. Each piece of the polyline is now halved
  until the curve at its quarter points stands within a micron of it.

- **A skin's error is read densely along its rows.** The sampled fit
  reads each span at its quarter points and climbs the peaks it sees
  there, and a peak between them could stand over the target: pipes
  through sections held a blended ring up to 1.06 times the tolerance off
  while stating less. Each skin is now read at eight points a knot span
  along every row, its peaks climbed, and refitted tighter where it
  strays; the error it states is the one read so.

- **A skinned loft, pipe or helical sweep states its edges' pcurve gaps.**
  A border two skins share is fitted on one and read on the other at
  samples, and three rails of a loft through squares and long rectangles
  turning between them stated less than their pcurves stood off. Every
  edge a skinned loft, pipe, pipe shell or helical sweep makes now states
  what its pcurves measure along it; the given sections', profile's and
  spine's own edges are left as they were.

- **A marched silhouette is drawn on its silhouette between the walked
  points.** A torus or spline face's silhouette was fitted through the
  points its walk found and kept whether or not the fit met the walk's
  chord of 1e-5, measured at those points only. Seen at a slant, a bumpy
  spline sheet's outline was drawn up to 1.6 off its silhouette, a fat
  torus's 0.7 and a spline cone's 0.01, in 24 of 60 views measured. The
  fit is now made against the silhouette itself, solved between the
  walked points, and refined until it holds there; one that cannot be
  held is not drawn. All 60 views now keep within 6e-6.

- **`check` does not name a crossing within a face's vertex tolerance.**
  A crossing counted only where it stood clear of both faces' trims by
  the faces' and their edges' tolerances, not their vertices'. Edges meet
  at a vertex only to the vertex's tolerance, and an imported vertex is
  often looser than its edges, so faces overlapping within a loose
  vertex were named as crossing there. The vertices' tolerances now
  count too.

### Performance

- **A face's single orientation probe stops at the first place that
  reads.** The probe that the mass properties, `inside_out_faces`, offset,
  meshing a solid and `fix_shape` ask went on to read up to three places
  and kept only the first; it now reads the first and stops, among the
  same candidate triangles, so every answer is unchanged. The mass
  properties of the largest corpus part take 8 percent less time and
  `fix_shape` on it 12 percent less.

## [0.11.2] - 2026-10-15

A patch release with one fix (#261): a face with no seam edge, such as a
torus band between two circles, now meshes the seam ogeom adds to the
chord asked. Where a rim stood a file's slop off its surface, that seam
came out in eight straight steps whatever the chord, standing 0.09 mm off
a torus of 19 mm section while the mesh reported its chord met. The
README and the exchange guide now list everything the mesh converter
rebuilds.

### Fixed

- **A seamless face meshes its seam to the chord when a rim stands off its
  surface (#261).** A face bounded by two rings that each go once round
  and no seam edge, such as a torus band between two circles, is drawn
  with a seam the mesher puts in between the rings. Where a ring stood
  off the surface by a file's slop the seam came out in eight straight
  steps whatever the chord, each standing 0.09 mm off a torus of 19 mm
  section. The seam is now halved until every step keeps to the chord
  and the angular limit, as the face's grid lines are.
## [0.11.1] - 2026-10-14

A patch release of fixes. Meshes drawn purely in strips, as OpenCascade
and truck write them, convert as ogeom's own do: a narrow band is no
longer read as a ball, a fillet beside a floor stays a cylinder, and a
bore whose rim strays from it takes no open sweep (#253). Drills through
a converted rough corner, grazing the chords where its facets meet the
rounds, cut and fill valid instead of refusing (#249).

### Fixed

- **A narrow band drawn in strips is not read as a ball (#253).**
  nist_ftc_07 drawn in strips has a flat face 2.2 mm by 33 mm between two
  small cones, each drawn as strips from rim to rim. Its eight corners lie
  in two rows of four, each row a trapezoid, so all eight lie on one ball
  of radius 90, and the three faces were recognized as that ball and built
  over the wrong part of it: it and two faces beside it were faceted. A
  ball's region now needs its facets to span about the same angle at its
  centre, as a drum's need about its axis: the middle strip spans more
  than ten times its neighbours', and the band comes back a plane between
  two cones, nothing faceted.

- **A fillet drawn in strips beside a flat floor stays a cylinder (#253).**
  nist_ctc_02's floor meets a fillet's foot through a few slivers, and the
  floor with them is a smooth region that bends only there. With no row
  through the fillet's middle, a surface of revolution about an axis far
  off held the fillet's rims, its straight strips and the floor alike, and
  took the fillet with the floor; its seams with the corner balls beside
  it could not be placed, and the fillet and both balls were faceted. A
  smooth region whose triangles are all but a few those of flat faces
  much larger than its other facets is now not tried as a sweep.

- **A bore whose rim strays from it, drawn in strips, takes no open sweep
  (#253).** nist_ctc_02's two large bores, drawn from rim to rim, have a
  far rim standing up to 0.019 off the bore, so no drum holds them, and an
  extrusion was fitted to the strips: its profile went all the way round
  without closing, the last strip spanned the gap between its ends, and
  no boundary closed on it, so the converter gave up the region as
  unplaced. A sweep whose profile does not close is now refused where a
  triangle of its region spans from one end of the profile to the other;
  the bore's strips stay facets, as nothing holds their vertices.

- **A drill grazing a converted round's chords cuts and fills valid
  (#249).** A converted rough corner meets each recognized round along
  chords that stand off the round by their sag, 0.0096, and lie exactly
  in the facets' planes. An edge's tolerance was read as a radius on both
  its faces, so a drill tangent to the round's rim ran within that of a
  chord for a third of a millimetre on the facet, and the boolean welded
  junctions 0.13 across that swallowed the facets beside it. Each face
  now reads an edge as far as it stands off that face, meeting the round
  along its own image of the chord; the two faces' crossings of a loose
  edge, as far apart along it as its tolerance over the crossing face's
  lean, are one junction; two sections through one face cut each other
  off a loose edge only where both are kept; and a kept piece whose
  boundary is all dust is dropped. The two drills of radius 1.4 beside the
  rough box's corner cut and fill valid again.

- **A converted facet whose seam with a drum folds a micron past its
  other side is built as a fan.** The converter looks for a face whose
  boundary crosses itself in its chart, but drew each face as the
  tessellator does first and looked closer only where that drawing
  crossed, so a fold shallower than the tessellator's chord went unseen.
  The rough rounded box came back with ten facets whose seam with the
  round beside them crossed the facet's other side by up to 1.2 microns
  over 0.07 mm by a corner, and a boolean through one of them found a
  piece of the facet with no inside. Every face is now drawn finely, and
  such a facet's seam is threaded straight.

## [0.11.0] - 2026-10-14

A minor release, since the API changes: `Callout` gains `text` and
`items`, so a callout carries its words and the faces it names through
STEP and the native format. Cylinders and cones now mesh between their
rims, a third fewer triangles on drum-heavy parts (door_handle_guide
2,956 to 2,004 at 0.0125 mm, against OpenCascade's 1,948), and every
triangle of a curved face, a cone's apex included, keeps to the chord.
The mesh converter reads meshes drawn in strips, as OpenCascade and
truck draw them, and a face it builds the wrong way round costs only
itself. A ray meets a spline face wherever it crosses it, and a solid
with thread faces checks in seconds instead of minutes.

### Added

- **`Callout::text` and `Callout::items` (#228, #227).** A callout carries
  the words it says as text of its own, written to STEP as a
  `TEXT_LITERAL` in an annotation text occurrence beside its strokes and
  read back from it, so a note survives a round trip. A note that draws no
  semantic annotation names the faces or edges it is about (a bore's
  thread callout), written as a shape aspect the callout is associated
  with and read back from it. Both travel through the native format.
  `Callout` gains two public fields, so building one with a struct literal
  needs them.

### Changed

- **Drums and cones are drawn between their rims (#213).** A cylinder or
  cone face carried a row of points through its middle that its straight
  rulings never asked for, a point per column. It is now drawn between its
  two ends, its rulings held as constraints. It keeps the row where one
  end closes on an apex, where an edge states it strays from the face by
  more than a hundredth of the chord (a fitted intersection in the file),
  or where it runs smoothly into another curved face; a point on every
  fourth column of it where it runs smoothly into a flat one. At 0.0125,
  OpenCascade in brackets: door_handle_guide 2,956 to 2,004 triangles
  (1,948), rail_guide 8,432 to 5,576 (6,484), sled_base 29,434 to 25,878
  (33,822), every triangle within the chord.

### Fixed

- **A bore running into a drafted wall, drawn in strips, comes back round
  (#253).** Drawn in strips from rim to rim, each strip's corners lie only
  on its two rims, so a first sample took strips of both the bore and the
  wall and fitted neither, and the two were one smooth region no sweep
  held: nist_ftc_06's two bores came back as 64 planes. A free smooth
  region no surface holds whole is now pieced along its strips' long
  sides, and taken where every piece lies on a drum or a cone.

- **A converted face closed the wrong way gives up its recognition alone.**
  When a converted body's volume disagreed with its mesh's, every
  recognized region of the body was faceted. nist_ftc_07 drawn in strips
  and listed from its 1994th triangle came back as 2,792 planes for one
  narrow band built over the wrong part of a ball. A face whose area is
  more than half again, or less than two thirds, of its facets' is now
  given up first, and the body checked again; it comes back with 167
  curved faces and 3 faceted.

- **A small ball's region is a cap only about a true latitude (#253).** A
  converted sphere region whose boundary ring lay within the reach of a
  plane was framed as a cap about that plane's normal, its pole inside.
  On a small ball the reach is loose: nist_ctc_02's corner blend drawn in
  strips came back on a ball nine millimetres across whose ring stood a
  tenth of a millimetre out of its plane, and the face's own mesh folded
  round a pole its boundary did not go round. A ring must now also lie
  within a hundredth of its own size of its plane.

- **A band left between two fillets is taken on its own surface before
  the sweeps.** A narrow band whose seeds all grew onto its still-free
  fillets fitted nothing, and the sweep pass then took it as a surface
  of revolution through its rim vertices whose axis missed the rims, so
  it was faceted. Whether that happened hung on the seed order, which any
  change to the mesh moves. A free region a canonical surface holds
  exactly is now taken on it before any sweep is tried.

- **A ring of flat strips is not read as a cone (#253).** A mesh drawn
  the way OpenCascade and truck draw a drum, each fillet a few flat
  strips from rim to rim, has the innermost strips of a block's four
  bottom fillets with all their corners on one shallow cone about the
  bottom face's middle, each strip a chord of a quarter of it. The
  converter took that ring for a cone, which then kept the fillets beside
  it from being placed: four of a rounded block's twelve fillets came back
  as sixteen planes. A drum or cone whose facet spans more than a sixth of
  a turn about its axis is no longer taken.

- **A band between two fillets is recognized once they are claimed
  (#250).** Recognition grows regions seed by seed, and a seed's sample
  grows across smooth edges: on a narrow band between two tangent fillets
  (a chamfer whose edges are rounded) it runs onto both fillets, fits
  nothing, and its first triangles were retired for good unless its first
  dozen vertices already turned through a tight bend. Whether the band
  was recognized then hung on the order the seeds came in, which any
  change of the triangle count reshuffles: eight more triangles elsewhere
  in grid_point_on_a_diagonal_boundary's mesh left its 45° band a fitted
  sweep whose rims could not be placed, and faceted. Every seed whose
  sample fits nothing is now taken once more after its neighbours are
  claimed, and the band comes back a cone.

- **A cone meshes within its chord at its apex.** Every ruling of a cone
  meets at its apex, and a triangle there joining points on rulings far
  apart cut inside the cone: on cone_point, 22 triangles stood up to 2.06
  chords off at 0.0125. Every curved face now holds the middles of a
  triangle's sides, lifted, to the chord as well as its centre, as #212
  did for surfaces curved both ways alone.

- **Every triangle of a surface curved both ways keeps to the chord
  (#212).** The grid holds its rows and columns to the chord, but a
  surface curved both ways twists between them, and a triangle across a
  cell could stand off it in its middle by twice the chord or more:
  freeform_panel's wall stood 0.22 off at a chord of 0.1, a ball of
  radius 5 meshed alone 0.17 off. A triangle whose middle, lifted onto
  the surface, stands farther than the chord from its plane has its
  middle put in, until none does: the wall keeps to 0.1 with 2% more
  triangles, and to 0.0125 with 0.1% more. A surface straight one way,
  drawn in strips along its rulings, is left as it was.

- **A straight curve meets a spline surface where its patches' control
  hulls say it may (#239).** A ray against a spline face was tried against
  every cell of the face's seeding it passed within the surface's bow of,
  with Newton from every segment beside it: a solid of a user's assembly
  whose rays ran along two thread faces of 49,152 cells took 47 s to check
  on one thread. The surface is cut into its Bézier patches once, each
  patch lies in its control points' box, and a ray is asked only of the
  patches, halved down to small ones, whose boxes it meets within the gap.
  Nothing a ray crosses is passed over: the hull holds every point of its
  patch. The solid checks in 2.4 s.

- **A ray meets a spline face wherever it crosses it (#246).** A cell of a
  spline surface's seeding that a ray passes within the surface's bow of
  was asked of only the ray's segment passing nearest the cell's middle,
  and on a large cell that segment is far from where the ray crosses: a
  ray from beside a thin spline wall of Voron 2.4 solid #2 missed the wall
  and read the wrong side. Every segment passing within the cell's bow of
  the cell's triangle is asked again.
- **A face's inside-out test probes only from points on the face.** It
  probed from the chart middles of the face's largest triangles, and a
  thin ring's triangles span its hole: their middles stood in the hole,
  off the face, and read whatever surrounds it. A middle the face does
  not hold is passed over.

- **A solid with thread faces checks in seconds (#239).** A curve's
  segments find the cells of a spline surface's seeding near them through
  a tree of the cells' boxes, built once per surface, instead of every
  cell being asked of every block of segments. And the finer trimming
  rings a face is drawn with where a ray meets it near its trim are drawn
  at a tenth of the chord only, no longer at a tenth of the angle as well:
  a thread face of 404 edges, each winding a whole turn, was drawn in
  819,000 points. The user's solid that took 6 s checks in 3.3 s.

- **A ray along a large spline face is intersected in time to its cells
  near it (#239).** A curve meeting a spline surface was tried against
  every cell of the surface's seeding along every one of its segments,
  and a cell the curve passes within the surface's bow of was polished by
  Newton from every segment beside it, twice each. A solid of a user's
  assembly with a thread face of 404 edges, its check casting rays along
  that face, took 177 s on one thread. Segments are now taken a block at
  a time against each cell, and a cell the curve does not cross is asked
  only of the segment passing nearest it and only where that segment
  comes within the cell's bow of the cell itself: the same solid checks
  in 6 s.

- **A section crossing a tight edge near its vertex is kept beside a
  loose edge (#235).** A face's weld is its loosest edge's doubt, and a
  stretch of its boundary shorter than that was welded away as dust:
  one edge stating 0.017 on a block's face welded the 0.02 a drill cut
  off the face's tight bottom edge by its corner, which the face across
  that edge kept, and the cut refused. The weld is now local: each strand
  end is welded at the doubt of the edges and vertices near it, onto a
  node of the face's boundary where one stands within reach, and two
  sections hugging an edge are matched by the weld at each point along
  them. A loose edge welds what crosses or ends near it, and a crossing of
  a tight edge elsewhere on the face stays where it is, from a hundredth
  of a micron off the vertex to past the loose edge's doubt.

- **A face curved one way only meshes in strips along its rulings
  (#213).** A drum, a cone or an extrusion took rows across its rulings
  to keep its cells near square, though the surface is straight along
  them: at a chord of 0.0125 step-conformance's cone_point meshed to 5,021
  triangles and sled_base to 60,802, against OpenCascade's 1,563 and
  33,822. Such a face now keeps only the rows its sag asks for, and each
  of its columns is held as a constraint from boundary to boundary
  through the grid's points on it, so no triangle bridges two columns:
  cone_point meshes to 914 triangles and sled_base to 29,434.
- **A mesh drawn in strips reads back on its rounds.** A fillet meshed as
  long strips between its rulings, one strip per facet, comes back as its
  cylinder: a run of facets between two planes is joined across its long
  borders, not across the short ones its ends share with a corner left
  faceted. A curved region in pieces, or a drum one of whose facets spans
  many times the turn of the others, is no surface: on a coarse mesh a flat
  face and the first facets of the fillets either side of it lie on one
  wide circle, and were read back as one cylinder in some triangle orders.

- **A face of a periodic surface splits along a curve entering it across
  the low end of its range (#229).** A cylinder whose side is two half
  faces refused a rectangle facing the junction where one half's angle
  range starts: the path was cut at that junction, and the piece cut off,
  moved a period away, never met the face. The path is now cut half a
  period from the face's middle, off the face, and the face splits across
  either end of its range.

- **An offset of a line or a circle projects onto a plane as a line, an
  arc or an ellipse (#226).** `project_edge_onto_plane` fitted a B-spline
  to an edge on an offset curve, the kind a thickening leaves on a
  cylinder's wall, where a sketch wants the line or arc it is. An offset of
  a line is now projected as the line moved square to it, and an offset
  of a circle in its own plane as the circle with the offset added to its
  radius, over the same angles: a circle onto a parallel plane, an ellipse
  onto a tilted one, every point of the edge on it.

- **A doubly curved face meshes within its chord across its cells
  (#212).** A grid cell whose row and column each sag the chord stood off
  up to twice it across its diagonal: a ball of radius 10 at a chord of
  0.0125 stood 0.016 off. Where any cell of a face sags past the chord,
  every cell of it gets its centre and the grid's lines hold each cell to
  its own triangles, so a ball keeps its symmetry; the ball now stands
  within the chord, and step-conformance's bowl_01 within it at 0.0125,
  0.05 and 0.1. A face of fewer than 128 cells keeps its grid as drawn,
  and a surface straight one way (a drum, a cone) its rows and columns.

- **A mesh drawn finely at a tangent rim reads back on its surfaces.**
  `solid_from_mesh` turned away a rim's facets where two surfaces meet
  tangentially (a dome on a drum, a ball on a pin): at the coplanar
  distance they lie in one plane with the facets across the rim, and a
  region grown before its neighbour refused them with the neighbour's
  still-free facets, leaving the ball in thousands of planes. A patch's
  members another curved region holds no longer stop it, and each region
  is offered its free facets again once all are found. A facet drawn to a
  cell's middle, a hair deeper than an equilateral one, is the surface's
  too. A mesh every coordinate of which single precision holds scatters by
  its quantum, and a coplanar distance asked below it is raised even where
  its facets lie too flat to show the scatter.

- **PMI is anchored to the part it is about (#227).** `write_step` hung
  every shape aspect, datum and datum target from the first part's
  product definition, whatever part owned its faces: a dimension or note
  on a second part's bore named the first part. Each now hangs from the
  part whose topology it names.

- **`split_face` splits along a closed curve across a periodic face's seam
  (#208).** A closed path is walked from where it first meets the face's
  boundary; one that crosses a seam meets the seam's two sides there a
  hair apart, and the hair between them, lying on the seam, was refused
  as a curve running along the boundary. A span whose ends stand within
  the path's own reach of each other is now no stretch, and the face
  splits: the piece inside the curve as two faces meeting at the seam,
  the piece outside wrapping round.

- **A boolean of a plane and a cone keeps the exact conic (#211).** The
  intersection stated the hyperbola or parabola exactly, over a span that
  runs past the cone's window for the face to trim, and the boolean fitted
  its pcurve on the cone by projection, which clamps to the window: the
  points past it stood off the surface, the fit was refused, and the pair
  was marched into a spline. A cone's image is now read in closed form past
  its window, so the section edge is the conic and projects onto its plane
  as `ProjectedCurve::Hyperbola` or `Parabola`.

- **A cone bounded at its apex meshes with no flat triangles (#210).** Since
  0.10.0 a cone's rows near its apex hold points on the same ruling, and
  the chart joined the apex to two of them in a triangle whose corners lie
  on one line. Its edges paired, so the mesh welded closed, but a reader
  that drops triangles with no area, as a 3MF reader may, found a hole. A
  face's flat triangles now split their neighbour across the long side
  instead, leaving every shared edge as it was.

- **Faces meeting at a vertex place it at one point (#217).** A boundary
  point whose chart parameter came out subnormal-tiny was held at zero by
  the triangulation, which then no longer found its exact parameters or
  the vertex it is anchored to, and lifted it off the surface instead: a
  corner 1.3e-14 mm from where its neighbours put it, and a mesh welded on
  exact coordinates open along three edges. The point is now known by
  where the triangulation holds it.

- **A fuzzy boolean at a coarse fuzz keeps both operands (#215).**
  `fuse_fuzzy` loosens every weld and probe the boolean makes, and those
  reach well past the fuzz: at a fuzz of 0.1 mm on ten-millimetre boxes,
  two boxes overlapping by 0.05 mm or touching fused to the first box
  alone, and boxes overlapping by 3 mm at 0.05 came out empty. A fuzzy
  result is now kept only where its volume agrees, within the fuzz over
  the operands' area, with the exact operation's (or with what the two
  can hold where that is refused); otherwise the fuzz is taken a tenth
  smaller, down to the exact operation itself. The same holds for
  `cut_fuzzy` and `common_fuzzy`.

- **A face no closed form speaks for leaves only itself to the mesh
  (#209).** One face whose volume integral has no closed form (a surface
  the charts cannot integrate, or a boundary they cannot walk) sent the
  whole solid to the mesh, so a cylinder's side drawn to the default chord
  read 0.6% short: text embossed on a Ø20 cylinder measured as volume
  lost. That face alone is now measured on its mesh, drawn to a tenth of
  the chord, every other face in closed form, and the strip between each
  of its chords and the edge's arc is added back, so a meshed face meets
  its exact neighbours on their edges. The drill through `nist_ctc_01` whose sliver is meshed now agrees
  with the part less the core to a part in a trillion, from a part in a
  billion.

- **A large face meets a fine chord (#216).** A face's grid stopped
  refining once a parameter direction reached 512 steps, and it refined
  from one end, leaving the rest of the direction at whatever the cap
  left: step-conformance's `bowl_01`, 380 mm across, stood 0.012 mm off at
  chords of 0.0033 and 0.00083 mm. A direction may now take as many steps
  as the deflection's `max_segments` (4096 by default), and refines in
  rounds, every interval that errs split at once, so a direction the cap
  stops is refined evenly. The bowl's rim stands within 1.05 chords, and
  meshes faster and lighter at the finer chord.

## [0.10.0] - 2026-10-12

A minor release, since the API changes: `ProjectedCurve` gains `Parabola`
and `Hyperbola`. An oblique plane now cuts a cone in an exact conic, and
`common_fuzzy` joins `fuse_fuzzy` and `cut_fuzzy`. It closes issue #202:
STEP sphere faces bounded by one parallel read as the half their ring walks
round, B-spline faces closed round with no seam written mesh, an edge the
segment cap stops is sampled all along, and cones mesh close to their chord.
Several silent wrong answers are fixed: a slanted drum's rounded rim that
lost a ring of material, and a NIST part whose balls read inside out. Every
edge of a slanted box rounds, a whole torus and a solid with voids offset,
and booleans near a ball's pole, a closed section's seam and a line's far
domain end no longer refuse. STEP import is about 30% faster on large
assemblies.

### Added

- **An oblique plane cuts a cone in an exact conic.** A plane oblique to a
  cone's axis and clear of its apex meets it in an ellipse where it is
  steeper to the axis than the cone's wall, a parabola where it runs
  parallel to a ruling, and a hyperbola where it is shallower, one branch
  on each nappe; each is now stated in closed form, the branches the
  cone's window reaches only, where it was marched and fitted before.

- **`project_edge_onto_plane` keeps parabolas and hyperbolas.** They come
  back as `ProjectedCurve::Parabola` and `ProjectedCurve::Hyperbola`
  (exact, on their principal axes; a segment where the conic's plane
  stands square to the target) instead of fitted splines, and a spline
  whose knots are not clamped (a periodic one) comes back clamped, still
  exact.
- **`common_fuzzy`** beside `fuse_fuzzy` and `cut_fuzzy`: the common at a
  stated fuzz, for a host retrying a tolerance-sensitive common as it
  retries the other two.

### Changed

- **STEP import is about 30% faster on large assemblies.** A curve the
  readers ask whether it runs through a chart pole (a cone's apex, a
  sphere's pole) is passed by where its guaranteed box stands clear of
  the pole, and `make_solid` settles which way a closed shell faces by
  probing a few of its faces instead of meshing all of them, meshing
  only where no probe settles. The DoomCube 350 assembly reads in 9.1 s
  instead of 12.8, the Voron 2.4r2 assembly in 18 s instead of 22, and
  every solid comes back the same.

### Fixed

- **A cone meshes close to its chord.** A cone's circumference runs from
  nothing at its apex to its widest at the base, and the triangulation,
  Delaunay in a chart scaled by the mean, bridged the wide end's columns
  three at a time between rows up to six cells apart: a pointed cone stood
  0.09 mm off at a 0.0125 mm chord (#202). A cone's rows are now held to a
  cell and a half, and it stands 0.016 mm off.

- **An edge the segment cap stops is sampled all along.** An edge drawn
  to a chord its curve needs more points for than the cap allows was
  split from its start until the cap was reached, and the rest kept only
  its seed points: long chords across a wavy rim, and the face inside it
  meshed past its rim. A freeform panel whose fitted rims a self-crossing
  boundary had drawn finer stood 5.5 mm off itself at a 0.0125 mm chord
  (#202). Such a curve is now split by rounds, every segment that asks at
  once, so the points the cap allows cover it evenly; it still says it
  missed its chord.

- **A ball whose pole a plane passes beside is charted clear of it.** A
  boolean restates a whole ball on a chart whose poles stand clear of the
  other operand's curved faces, and left its planes out: a plane through a
  pole cuts a meridian, which splits there exactly. A plane a hair beside
  the pole cuts a small circle round it instead, whose image the chart
  cannot hold, and the boolean refused ("a branch passing beside a chart's
  pole fits nothing yet"): two fused boxes cut by a ball, at chain stress
  seed 9, twice. Such a plane now counts against the pole too.

- **A B-spline closed round with no seam written meshes.** Some modellers
  write a closed B-spline face (a drum, a dome, a ring of a torus) bounded
  by its rims alone, with no seam edge. Each rim then runs the chart from
  one end column to the other and meets no other, the boundary enclosed
  nothing, and the face meshed to no triangles at all (#202). The seam is
  now put in: the surface's own column at the period's join, exact, closing
  the rims into one loop, or running to the pole where the surface
  collapses a row to a point.

- **A ball face a file bounds by one parallel is the half it walks round.**
  A STEP sphere face bounded by its equator alone, with no vertex at a
  pole, read as whichever half the chart's start fell on: the hemisphere
  above came back as the one below, 20 mm off for a radius of 20 (#202).
  The pole the file leaves unwritten is now put in on the side the ring is
  walked round, its material's side about the face's own normal, and the
  face built as a cap on it, as a cone bounded by one ring already has its
  apex put in. The NIST part `nist_ftc_06` read its two balls inside out
  for the same reason, as pockets where its loops and the concave fillets
  at their feet stand them on the plate: its volume is now 3,323,546 mm3,
  not 3,289,989.

- **A closed section starts at a vertex it passes a hair from.** A circle
  the boolean cuts from a plane and a ball or a drum starts where a world
  axis falls in its plane, and a piece reaching across that start is split
  there. On a solid a hair off the world axes the start fell a few
  ten-thousandths of a radian from a face's vertex on the circle, the
  sliver between was taken for dust, and the arc beside it was left short
  of its neighbours: rounding every edge of such a box refused with "a
  face boundary strand dangles". Such a circle now starts at the vertex.

- **A curve answers an ulp past the end of a domain that ends far out.**
  A curve refused a parameter more than the parametric confusion (1e-9)
  past its domain, but a line's unbounded domain ends at a billion, where
  one ulp is a tenth of a micron, and a projection's or an intersection's
  sample stepped to the end lands an ulp past it. The corner tool's cut
  then failed on a domain error. The slack is now a few ulps of the
  domain's ends where that exceeds the parametric confusion; a micron
  past is refused as before. In the stress round scenario the last
  refusals of every edge of a slanted box went with it: 900 of 900.

- **A slanted drum's rounded rim keeps the material it should.** A rim
  fillet's torus band runs between a ring on the wall, on the tube's outer
  equator, and a ring on the cap. The equator is the seam of the torus
  chart's tube angle, and at a slant rounding read its row at the far end
  of the chart, so the band ran the other three quarters of the tube and
  the drum lost a ring of material, valid by `check` all the same: 8.74
  against 9.37 for a drum of radius 1 and height 3 rounded at 0.2. A row
  on the seam of a chart periodic in its second parameter is now taken
  at the end nearer the other ring. In a stress scenario rounding a
  slanted drum's rims, 65 cases in 100 were wrong; none are now.
- **Every edge of a box at a slant rounds.** A box corner's ball meets
  two of its bands along meridians through the poles of its chart, and the
  boolean splits such a circle at the poles found by projecting them onto
  it. At a slant the projection leaves the poles a nanoradian or so off,
  past the slack the half-meridian's straight pcurve allowed, so the arc
  was marched and fitted instead, and the bands then could not close
  against the fitted rim: three corners of a slanted cube in eight refused.
  The slack is now the angle the ball's radius turns into the confusion
  distance, and the lift back onto the circle still measures the line.
  In the stress round scenario, which rounds every edge of a slanted box
  in three cases in ten, refusals fell from about 17% to about 3%.

- **`offset_shape` offsets a whole torus and a solid with voids.** A
  whole torus is one face bounded by nothing but its seams, and the offset
  assembled faces from their trims: it refused every torus, out or in. It
  is now the moved torus whole. And a solid with a void sewed all its
  faces into one shell and refused when they made two; each shell is now
  sewn on its own and the solid put together from them.
- **An unclamped spline evaluates at its ends and takes a knot there.**
  With its first or last knot repeated, an unclamped spline's span lookup
  landed on a span of no length at that end of its domain and evaluated
  to NaN, and a knot inserted at the end of its domain was worked in the
  wrong span and moved the curve. The lookup passes over empty spans, and
  insertion works in the span the value opens.
- **A face's seam columns meet the rest of its boundary within its radii.**
  A boolean reads whether a face uses both columns of a seam by whether
  each column's ends meet the face's other edges, and it held them to a
  fixed parameter snap. A torus fused with a drum keeps a meridian loop
  whose pcurve stops 7.6e-5 short of the seam column it closes on, two
  hundredths of a micron in space, and cutting the fused solid dropped the
  column and refused. Ends that meet across the chart within both edges'
  radii in space now meet.
- **A ray meeting a curved face near its trim reads the true trim.** The
  exact classifier counts a ray's crossing of a face where the hit lies
  inside the face's rings, polylines straight in the face's chart. Beside
  a sphere's pole the chart stretches, and a ring whose points keep to
  their edges in space strays between them in the chart by many times the
  band a hit is held to. A ray meeting a ball just outside a sliver of it
  two thirds of a millimetre from its pole counted the sliver, and the
  check called a valid common's ball face turned in. A hit within fifty
  bands of a curved face's trim is now read against its rings drawn ten
  times finer.
- **The mesh converter facets a face whose boundary does not close.** A
  printed part of the DoomCube 350 assembly, meshed and converted back,
  came back with a face whose wire's edges stood 4 mm apart. They were
  drawn alike from both sides, so the solid meshed closed and passed the
  converter's checks while it checked broken and a tenth short of the
  part's volume. A built face whose wire does not meet end to end, or
  whose boundary does not close in its chart, is now faceted
  (`FallbackReason::BoundaryOpen`).
- **A thin lump is not taken for a membrane.** A leaning torus whose tube
  pokes a little way past a drum's or a cone's flat end: the cut keeps the
  tip as a lump of its own, bounded by a patch of the torus and a patch of
  the end on one loop. The two face apart and are as thin as a weld near
  the loop, and read at one point there they were dropped as a membrane:
  the cut came back valid and 0.009 cubic millimetres short. Every inner
  point of the one face must now lie on the other.
- **A boolean refuses a face it made that does not close in its chart.**
  A drill whose top line touched `nist_ftc_11`'s bore rim where the bore
  meets the torus above it came back as a cut whose drill wall did not
  close in its chart, a solid that checks broken. Each face the boolean
  makes is now held to the chart test `check` makes of it, and one that
  fails is refused, which lets the boolean's retries try another way.
  `ogeom_algo::chart_breaks` asks the test of one face, passing over the
  wires a caller took over as they stood.
- **A ball whose pole a section passes beside is charted clear of it.** A
  ball and a torus whose tube passed a thousandth and a half of the ball's
  radius from its pole marched their section across the chart's collapse:
  snapped onto the pole on the ball, a run of points on the torus. Every
  result lost a sliver and the fuse's ball face faced into the material.
  A whole ball is now charted about an axis whose poles stand clear
  wherever the other operand comes within a hundredth of its radius of a
  pole, not a thousandth.
- **A face whose integral does not settle is meshed alone.** Where one
  face's closed-form volume integral could not settle (a sliver of a small
  round that a large drill's section meets tangentially), the volume
  meshed the whole solid instead: a drill through `nist_ctc_01` measured
  258 cubic millimetres off the part less the drilled core, and took 70
  seconds at a chord of 2e-4. That face is now measured on its own mesh
  and every other face in closed form; the cut agrees with the part less
  the core to a part in a billion, in a fraction of a second.
- **A short arc between two rims is found.** Two drums across each other,
  or a drum across a cone, cross in one loop that their rims can cut into
  a long arc and a short one tucked into a corner of both walls' charts.
  The short one lay inside one pair of the marcher's sampling cells and
  crossed neither, so it was never seeded, and the boolean refused: the
  kept pieces did not close. The rims of drums and cones are now
  intersected with the other surface, as a spline's borders were, and
  every piercing seeds.
- **`check` does not call a gap broken by rounding.** An edge of a drill
  through `nist_ctc_03` stated its tolerance at the gap its pcurve stood
  off, and `check`, measuring the same gap again at the edge's start, read
  it ten ulps wider and called the cut broken. A gap past its tolerance by
  a billionth of it is within it.
- **A hole named a turn round its wall is drawn as a hole.** A face on a
  periodic surface whose hole's pcurves stand a whole period from the
  face's own bound (a converted bore whose cross bores' loops came out at
  u = 2pi beside a wall drawn from -pi to pi) drew whole, with each hole
  drawn again as a disc over it: more area than the surface has, and the
  solid's mesh, its faces turned to agree, came out open and partly inside
  out. Each ring is now slid by whole periods to stand where the face's
  bound does.
- **A piece read by its corner is read again at its centre.** A face
  leaning a hundred-thousandth of a radian into the other solid, the two
  sharing an edge, was read at its roomiest probe, which stood by that
  edge within rounding of the other's face and read outside: the fuse
  kept it and its pieces did not close. Where kept pieces do not close,
  the boolean now tries again with each piece read at the centre of its
  area where that lies inside it.
- **A blend's wedge states its edges' pcurve gaps before it is applied.**
  A marched band's rails and legs carry fitted pcurves that stand off
  their edges between their samples; one wedge's edge stated 1.3e-4 where
  its pcurve stood 1.8e-4 off. The wedge's edges now state what their
  pcurves measure along them before the boolean reads it.
- **A face split along a projected curve states its new edges' gaps.** A
  curve dropped onto a face along its normals lands through a fitted
  pcurve where no closed form projects it, and the fit strays between its
  samples: a circle dropped onto a drum's wall made edges stating seven
  times less than their pcurves stood off. Every edge the split makes now
  states what its pcurves measure along it.
- **A sheet offset or thickened states its edges' pcurve gaps.** The walls
  between a thickened sheet's layers are ruled between their rims, fitted
  where no closed form joins them, and their pcurves lift to the fit rather
  than to the rim: a drum's wall thickened by one stated 1e-7 on rims its
  end walls' pcurves stood 2.1e-4 off, and a moved fitted face's pcurves
  strayed between their samples. Every edge the offset or the thickening
  makes now states what its pcurves measure along it; the given sheet's
  own edges are left as they were.
- **A drafted face's rebuilt edges state their pcurves' gaps.** An edge
  the draft rebuilds on a fitted support was widened by how far its
  projected pcurve's samples stood off the curve, not by how far the
  pcurve strays between them: on a drum drafted about a neutral tilted
  off its axis an edge stated 1e-4 where its pcurve stood 2.9e-4 off.
  The gap is now measured along the whole edge, and the edge and its
  vertices state it.
- **A drill whose seam grazes a face keeps its seam.** A through drill
  whose top line, its seam, stands a hair above a block's top face cuts
  the face along two lines either side of the seam, closer than the
  weld, and sewing joins them into one line the drill's face meets on
  both sides. The line kept the two plain pcurves the drill's face had
  on either side, and a lookup by surface found the first for both
  walks: the cut came back invalid. Sewing now states such a line as the
  seam it is.
- **A solid's mesh settles its shared edges before it is stitched.** Where
  a face's rings crossed, its edges were drawn finer and the faces across
  them drawn once more, and nothing asked those again: a narrow strip
  whose one rim a neighbour refined could then cross itself, its coarse
  rim sagging over the fine one, and the strip lost a boundary point the
  neighbour kept. Faces drawn again are now asked again, and their edges
  drawn finer in turn, until none comes up short. A face asked which
  chords it wants starts from the chords its edges are already held to.
  One solid of the DoomCube 350 assembly meshed 34 edges open; it now
  meshes 27 open, on a face whose rims touch within a micron.

- **A point on a thread's flank is placed on its own turn.** The exact
  classifier folds a hit into a face's chart by whole periods to find the
  copy inside the face's rings, and tried one period either way. A face
  winding many turns round a drum (a thread's flank, its ring unrolled
  over fourteen turns) holds the hit many periods over, so the hit read
  outside the face, a ray through it was miscounted, and `check` called
  faces of a threaded pin of the Voron 2.4r2 assembly inside out. Every
  period within the rings' extent is now tried, nearest first.
- **`check` no longer calls a sliver face inside out.** The orientation
  probe steps off a face's largest triangles to either side. Off a
  sub-micron step folded back between two nearly coplanar faces, a step
  longer than the sliver is wide passes under one neighbour or over the
  other and reads their sides, and a correctly turned sliver was reported
  reversed. A step now counts only where the segment off the face meets
  no other face; otherwise a shorter step is tried, then another
  triangle. Twenty solids of the DoomCube 350 assembly and one of the
  Voron 2.4r2 assembly checked broken for this since 0.9.1.
- **A face wrapping a closed chart once is trimmed in its chart.** A ring
  that comes out a whole period below a surface closed but not periodic
  (a B-spline band all the way round) is slid back into the chart, by
  whole periods. The count rounded a turn's length over the chart's span,
  a hair over one, up to two, and the ring came to rest a period above
  the chart. Every point of the face then read outside its own trim:
  the exact classifier missed the face, `check` called the faces beside
  it inside out, and booleans asking it which side a piece is on were
  answered wrongly. Three faces of one solid of the Voron 2.4r2 assembly
  were named.
- **A STEP or IGES face on a small chart reads again.** Where a face's
  pcurves missed each other at a corner by more than a millionth of the
  chart's span, the reader took the miss for a pole and laid a pole edge
  in it. On a patch whose chart spans a fortieth of a unit, a miss of
  6e-8 is that, and too short to be an edge: the segment was refused and
  the whole file with it. The Voron 2.4r2 assembly (1563 solids), which
  read before 0.7.0, was refused from then on. A miss too short in the
  chart to be an edge, whose ends stand at one point in space, is now
  the fitting slop it is.
- **A drill tangent to a wall along a round's edge cuts the sliver behind
  it.** Where two faces touch along a curve that runs along a boundary
  edge of one of them, the solid's boundary turns at that edge, and the
  other face can pass from outside the solid into it across the line it
  touches: a drill tangent to a wall exactly where a concave round leaves
  it, the drill flatter than the round, lines the cut from the tangent
  line on. The touch was set aside as splitting nothing, the drill's wall
  stayed one piece, and the cut was refused as not closing. The edge is
  now carried onto the touching face as a contact, which splits it there.
- **A polynomial's leading term is kept where its roots are large.**
  `solve::roots` dropped a leading coefficient small beside the largest
  one, which `(x - 1e6)^4` is, and solved what was left: a line anchored
  some ten thousand units from a cone found no crossing, and the quartic
  came back with a root at a quarter of the true one. The leading term is
  now dropped only where it is nothing at the size of the remaining
  polynomial's roots.
- **A boolean traces every edge and vertex of its arguments.** `fuse`,
  `cut`, `common` and `split_sheet` record each edge and vertex as modified
  into its pieces in the result where it was cut short, split or rebuilt,
  copied where it was rebuilt whole, and deleted where nothing of it is
  left, as they already recorded faces. One the result holds as it stands
  traces to itself. Before, a shortened edge, a removed corner and the
  tool's edges and vertices traced to themselves, which the result does not
  hold.

## [0.9.16] - 2026-10-09

A patch release: no breaking API changes. STEP bands whose rims start on
different columns read seamed, a drill nicking an edge a few microns deep
closes because the two faces of the edge decide a hug once, a thin drum
through a long drill's wall is cut on both sides, an affine rebuild
records every sub-shape in its history, and a round hole thinner than a
micron cut in two halves stays a hole. Mesh conversion no longer reports
a fallback for an empty region.

### Added

- **`ogeom_algo::split_ring_on_column`** splits a closed circle rim of a
  revolution band where it crosses the column of the other rim's vertex,
  so `make_band_of_rings` can seam the band there.

### Fixed

- **A drill nicking an edge a few microns deep closes.** Where a tool
  crosses an edge twice and the two faces of the edge meet it in curves a
  few microns apart, both faces now keep the edge between the crossings
  when the tool's own face cannot hold the sliver between the curves. The
  cut was refused.
- **A shear or an uneven scale traces every sub-shape.**
  `general_transformed_shape` records each wire, vertex and degenerate edge
  of the input as modified into its rebuilt twin, as it already did for
  faces, other edges and containers, so `trace` and `is_affected` answer
  for every sub-shape, alone or composed through `then`. `baked_where_placed` and `restate_geometry`
  rebuild the same way and record the same. The docs of `Model::placed` and
  `Document::occurrences_of` state that they record no history and how a
  sub-shape is found in the placed shape.

- **Mesh conversion reports no fallback for a region it never builds.**
  Where the coplanar distance widens and one region the default distance
  found is split between two regions found wider, putting both back gives
  it back once; the second copy, holding no triangle, used to be counted
  as a curved face left faceted.

- **A STEP band whose rims start on different columns reads seamed.** A
  cylinder or torus face bounded by two whole circles whose vertices stand
  apart has one rim split on the other's column and takes its seam there,
  and the faces across the split rim walk both pieces. Such a face used to
  warn that no seam could be synthesised and mesh through a fallback.
  `make_revolution_band` now says why it refuses such rims.

## [0.9.15] - 2026-10-09

A patch release: no breaking API changes. Shapes far from the origin keep
their facing and features across conversion, sewing, extension, chamfers,
pipes, meshing, booleans and healing, with converted cylinders and spheres
keeping their surfaces there. Mass properties are exact under a scaling
placement, a corner rounds the same whichever edge comes first, and mesh
conversion depends far less on how its triangle list is ordered or which
way a fitted axis points. Timing tests measure thread processor time.

### Fixed

- **Converted cylinders and spheres far from the origin keep their
  surface.** Axes and sphere centres are fitted about the part's own
  points; a cylinder at 1e7 and a sphere at 1e6 came back as facets.

- **Mesh conversion no longer depends on which way a fitted axis points.**
  A band's seam runs through wherever its rims start, a seam through a
  hole's corner counts as crossing it, and a closed piece on a sphere or
  torus is a whole face; such a piece used to leave an empty shell and an
  error.

- **A circle lying in a plane far from the origin keeps its exact image
  there.** At 1e7 its axes were rounded off square, its edges became
  fitted curves, and measuring the solid took seconds.

- **Shapes far from the origin keep their facing and features.** Ring and
  volume sums are taken about a point of the shape instead of the world
  origin in mesh conversion, sewing, face extension, chamfers and wedges,
  fixed-section pipes, slit and crack detection in meshing, and ring areas
  in booleans and healing; at 1e6 and beyond a converted plate no longer
  comes back inside out, a chamfer no longer refuses, and a hole of a few
  microns no longer reads as a slit.

- **Converted faces meet edge to edge on far more listings of one mesh.**
  A drill closed round by facets it took back is built as a band with a
  seam, and a long round fitted leaning is put on the exact cylinder
  between its two largest flanks with its stray rows; on `nist_ftc_07`,
  faces meet in 29 of 30 rotations of the triangle list, up from 9.

- **`solid_from_mesh` keeps rounds whose distance widened elsewhere.**
  When the default coplanar distance widens because a fit pressed against
  it, a curved region found at the wider distance that cannot be built is
  put back as the default distance found it, instead of being faceted; a
  mesh whose triangle list starts elsewhere no longer facets long rounds
  beside small fits the wider distance admits.

- **A corner rounds to the same solid whichever edge is asked first.**
  Where a curved face takes over at the corner, a straight edge's blend no
  longer runs on past its own vertex into material the solid continues
  with; it took off more than the union of the blends, or was refused.

- **Mass properties are exact under a scaling placement.**
  `volume_properties` and `surface_properties` integrate faces under a
  placement that scales on their exact surfaces, rigidly placed and scaled
  back, instead of falling back to the tessellation (a ball scaled by 0.5
  measured about 3e-4 low).

## [0.9.14] - 2026-10-08

A patch release: no breaking API changes. The kernel runs in the browser:
on `wasm32-unknown-unknown` it reads no clock it lacks (`ogeom_core::clock`,
#149), and a host can lend parallel stages its threads, such as web
workers sharing memory (`parallel::set_pool`, #150). A marched section's
stated error bounds its miss against the real crossing, and point-in-solid
classification answers by the winding number where every ray meets a
tangency.

### Added

- **A host can lend parallel stages its threads (#150).**
  `parallel::Pool` and `parallel::set_pool` run `map_ordered`'s workers on
  a pool the host provides (web workers sharing memory on
  `wasm32-unknown-unknown`, through `wasm-bindgen-rayon`); output stays
  bit-identical at any worker count, a panic on a worker resumes on the
  caller, and items a pool leaves unclaimed are run by the caller.

- **The kernel runs on `wasm32-unknown-unknown` (#149).**
  `ogeom_core::clock::Instant` is `std::time::Instant` where there is a
  clock; there, timings read zero, or `performance.now()` with the new
  `web` feature, and parallel stages stay on the calling thread. A CI job
  meshes a box and a cylinder in Node.

### Fixed

- **A marched section's stated error bounds its miss against the real
  crossing.** Its curve and both pcurves are measured against the surface
  crossing found at each station, so a shallow crossing no longer states
  less than it stands off (3.9e-5 where the fit stood 1.9e-3 off), and the
  section edge's tolerance holds.

- **Point-in-solid classification answers where every ray meets a
  tangency.** After all six rays are refused, the point's generalized
  winding number against the faces' meshes answers when the point stands
  clear of each face's mesh band and the number is within a hundredth of 0
  or 1; otherwise it still refuses. Drafted slabs whose pad probes were
  refused at depth 3 now cut and fuse.

## [0.9.13] - 2026-10-08

A patch release: no breaking API changes. A boolean's shell closure no
longer rides on its sections' fit error, so tighter sections no longer flip
a shell from closed to refused, and a merged junction stop splits every
edge it stands on (every band order of a curved corner closes). Marched
sections fit from knots the cubic's error estimate places (crossed
cylinders about a third faster), a boolean on a large solid it has read
before reads it back (`Model::kept_read`), nearest points on B-spline
patches are found by branch and bound (about 11 times faster), and exact
mass properties refine one panel at a time.

### Added

- **`gauss_kronrod_21_rule`** in `ogeom-math`: the ten-point Gauss and
  twenty-one-point Kronrod pair with both weights per node.

- **`Model::kept_read`** keeps what a reader finds of a shape and
  everything below it until any of it changes; `SolidBoundary::of_kept`
  prepares a solid's boundary through it.

### Changed

- **Exact mass properties refine one panel at a time.** Each boundary
  panel of a face integrated in its chart takes a ten-point Gauss and
  twenty-one-point Kronrod pair, and only the panels whose estimates miss
  are halved, instead of doubling every panel; values settle to a part in
  a million million or to rounding. Exact volume of the largest corpus
  part is about 20% faster on one thread, and two converted faces that
  fell back to the mesh are integrated exactly.

- **Marched sections fit from knots placed by the cubic's error
  estimate.** The joint fit of a marched intersection starts its knots
  where a cubic's estimated miss between samples holds, measures every
  sample, and falls back to the free fit only where it misses: a fuse of
  two crossed cylinders is about a third faster, with a tighter true miss.
  `fit_points_joint_from` fits joint samples from given knots.

- **Nearest points on a B-spline patch by branch and bound.**
  `project_on_surface` and `SurfaceSeeds` search the patch's knot spans
  from the nearest bound down instead of a grid of one seed or more per
  span: about 11 times faster on a 100 by 100 span patch, STEP import of
  the largest corpus part about 20% faster, thick spline sheets about 11%.
  On a B-spline patch the projection is the nearest point of the whole
  patch to within a hundredth of the confusion, `samples` is not read
  there, and two equally near feet resolve to the smaller u, then v.

- **A boolean on a large solid it has read before reads it back.** The
  solid's faces and prepared boundary are kept instead of walked again: a
  corner hole into a 582-face plate costs about 1.5 times the same hole
  into a 22-face plate, down from 2.1 to 2.5 times.

### Fixed

- **A boolean's shell closure no longer rides on the sections' fit
  error.** Two crossings either side of a graze count as one junction
  within a fixed reach, a section reading just outside its face along the
  boundary at a node is ordered inside it, a short stub runs along an edge
  only if it stays close to it for its length, and a seam join clamps its
  height to the surface's range. Sections fitted more tightly no longer
  flip a drill's shell from closed to refused.

## [0.9.12] - 2026-10-08

A patch release: no breaking API changes. Fixes: a solid with rounded edges
and corner balls offsets and shells; a spline fitted through every point
stays on the shape between them; projected pcurves state their widest miss
between samples; a refused revolution band leaves its rims untouched. Speed:
exact mass properties about twice as fast, `check` about twice and
`fix_shape` 2.4 times, booleans on drilled parts about 35 to 40%, marched
fillets and booleans about 10%, thick spline sheets about three times,
pipes through circles about thirteen times, IGES reading about 2.7 times,
tessellation and mesh welding faster with identical meshes, and kernel
hash maps with a fixed seed so iteration order is the same in every run.

### Added

- **`ProbeCache`, `check_kept` and `inside_out_faces_kept`.** The
  face-orientation probe keeps its face meshes and ray answers by what they
  read, so a check after a mend recomputes only what the mend changed.
  `newton_system_fixed_lazy` takes the residual and the Jacobian
  separately.

- **Face meshes kept by what a face is made of.** `FaceMeshCache`,
  `triangulate_face_kept`, `triangulate_with_face_meshes`,
  `face_meshes_for` and `ogeom_topo::face_content`: each face's mesh comes
  back beside the welded whole, drawn once, and can be kept across models
  by what each face is made of; a kept mesh is the one drawing it again
  gives.

### Changed

- **A face piece's interior probes past the first are found only when the
  first does not settle it.** `boolean_many_faces` is about 35% faster
  than 0.9.11 in all and `boolean_drill` about 40%, with the same
  results. A new node links to its children without searching their
  other holders.

- **Booleans hold a piece's probes to its exact rings without resampling
  every circle in full:** a circular strand's fine samples are drawn only
  where a probe's ray can cross them (`boolean_many_faces` about 14%
  faster).

- **Exact mass properties are about twice as fast.** A chart face's exact
  samples are summed once rather than in every run, rational spline
  patches take the fast isoline evaluation, and faces are prepared and
  their orientations checked in parallel. Values agree with before to
  5e-13.

- **Self-intersection checks prove most face pairs apart before measuring
  them.** Each piece is held in a guaranteed hull (a spline's control net,
  a plane's corners) and pairs are separated by box or by projection,
  halving the closer one within a budget; only pairs it cannot separate
  are measured exactly. Thickening a spline sheet takes about a third of
  the time. Curve-surface intersection builds its seeding grid with each
  corner evaluated once.

- **Pipes through sections match each section's start without copying
  the loop.** A pipe through two circles is about 13 times faster.

- **Tessellation repeats less work, with identical meshes.** The stored
  tessellation indexes each face mesh once and draws each edge once (13%
  faster on the largest corpus part); grid refinement samples each grid
  value once (whole-part meshing about 10% faster, a torus 12 to 18%); the
  repair passes find border T-junctions through a grid and walk cracks in
  linear time, and welding chains its cells (a 1M-triangle weld 12%
  faster).

- **The orientation probe is about twice as fast, `fix_shape` 2.4 times.**
  A solid's boundary prepares each spline surface for rays once, and the
  curve-surface polish works out derivatives only at points it steps onto.
  `check` on the largest corpus part goes from 176 to 80 ms and
  `fix_shape` from 474 to 196 ms at 20 threads, with the same diagnoses.

- **`check` runs in parallel.** Its per-edge, per-wire, per-face and
  containment checks give the same diagnosis in the same order;
  self-intersection finds candidate face pairs by a box sweep and measures
  them in parallel (11.9 s to 1.3 s on the largest corpus part at 20
  threads).

- **The IGES reader derives edge pcurves in parallel** ahead of the face
  walk, as the STEP reader does, with identical warnings and output.

- **`repair_same_parameter` measures edges in parallel with the dense
  lifted-gap measure,** at matched parameters on the surface where the
  pcurve is placed, so it widens edges the old 24-sample check missed; an
  edge it cannot measure keeps its flag. `fix_shape` builds a face's
  surface only when an edge needs a trim and fits trims in parallel; a
  trim refused at the cap reports its exact offset, and
  `pcurve_fit::fit_projected_pcurve_within` returns the refusal as a value.

- **Faster nearest points on cones, tori and swept surfaces.** Cones,
  tori and trimmed elementary surfaces are projected in closed form, and
  extrusions and revolutions through one scan of their profile instead of
  a seed grid; exact hidden lines are about 15% faster.

- **Marching solves on the stack.** The intersection walker, the marched
  intersections' solves and the rolling-ball fillet march use the
  fixed-size Newton solver and one jet per surface, and the walker reuses
  its last Jacobian for the tangent (`ogeom_intersect::walk::null_tangent`,
  `Condition::tangent_from`): marched fillets about 11% and marched
  booleans about 9% faster.

- **Kernel hash maps are hashbrown with a fixed seed.** `ogeom_core::FastMap`
  and `FastSet` (in `ogeom_core::collections`, with `FastHasher`) replace
  the std maps in every kernel crate, so a map's iteration order is the
  same in every run. Welding a 1M-triangle mesh and reading a
  1M-triangle STL take less than half the time; STEP import is 5 to 30%
  faster, and looks instances up in a table indexed by instance number.

- **Converting a mesh to a solid repeats less work.** Edge fits, face
  meshes, face areas and exact face integrals are kept across the
  conversion's builds. The result is unchanged; the largest corpus part
  converts about 1.5 times faster on one thread, mid-size parts 15 to 40%
  faster.

## [0.9.11] - 2026-10-07

A patch release: no breaking API changes. A prism far from the origin keeps
its rails and its facing, a torus whose hole nearly closes meshes like an
open one, a dropped coincident piece and every sub-shape of a moved shape
trace through history, and `check` reports a seam stored as plain pcurves
while `fix_shape` repairs it (a solid saved from 0.9.9). Building
primitives is about five times faster again, exact hidden lines run in
parallel and ask fewer faces, and hidden line removal splits an edge where
it passes behind a contour. The benches cover mass properties, `check`,
`fix_shape`, mesh conversion, exchange, hidden lines and offsets, and take
`--threads`.

### Added

- **`PreparedSurface` keeps a surface's seeding for many curve
  queries.** `intersect_curve_surface` samples the surface into cells on
  every call; a `PreparedSurface` samples it once and gives every curve
  the same answer.

### Changed

- **Primitives are built without meshing them to check their facing.**
  Boxes, wedges, polyhedra, cylinders, spheres, cones and tori wind their
  faces outward by construction, and skip the probe `make_solid` runs on a
  closed shell to turn one facing inward; building a box is about five
  times faster again. `make_solid` still turns imported and sewn shells.

- **Exact hidden lines run in parallel and ask fewer faces.** Edges,
  silhouettes and their visibility are worked out a curve at a time across
  threads, in the same order at any thread count. Each point asks only the
  faces whose projected box covers it, starting with the face that hid the
  previous point, and each spline face is sampled once rather than once
  per point. The drawing is unchanged; a corpus part of 520 faces draws
  about twice as fast on one thread.

- **Exact hidden lines ask a spline face about the stretch of the ray
  inside the face's box,** not the whole ray. The seeded search no longer
  polishes from every cell of the face, and a 664-face part draws in about
  30 s on one thread instead of 20 minutes. Edges keep their visibility.
  Where a silhouette point's ray grazes its own spline face, or crosses it
  again a fraction of a micron away, a short piece of the silhouette can
  now come out the other way.

### Fixed

- **Projected pcurves state their widest miss between samples.** Quarter
  points and the peaks beside them are read, not only the samples and
  their middles; an edge's measured pcurve gap finds its widest peak where
  the pcurve is paced differently from its curve, so an edge read from a
  file no longer states a tolerance below how far its fitted pcurve stands
  off. A fixed-parameter pcurve fit that reaches interpolation stands on
  averaged knots.

- **Dividing a face along an iso line of a surface with no closed form**
  (an offset of a spline) keeps the new edge within the confusion distance
  of the surface between its samples.

- **A spline fitted through every point stays on the shape between
  them.** Curve and grid fits that refine all the way to interpolation take
  averaged knots, so they no longer swing off between their points (an
  11 x 11 saddle grid from 1.7 off to 1.3e-3). The fit docs state that the
  error of a fit through bare points is measured at those points only.

- **A solid with rounded edges and corner balls offsets and shells.**
  `offset_shape` and `make_thick_solid` on a box with all edges filleted
  refused the result as collapsed: a blend arc whose end fell on its
  circle's start was rebuilt the long way round, and an offset corner
  sphere lost its frame, putting its pole on a face corner. Blend arcs keep
  their quarter turn and offset spheres keep their parameterisation.

- **A refused revolution band leaves its rims untouched.** A band whose
  rims' vertices sit on different columns is refused, but its rims kept a
  pcurve on an unused copy of the surface; `repair_same_parameter` then
  widened them to a chord of the circle (up to 40 mm on the NIST test
  parts). They now stay at the file's own slop.

- **A prism far from the origin keeps its rails and its facing (#144).**
  `make_prism` of a face 1e6 to 2e7 from the origin no longer fails with a
  domain error on its rails or returns a solid wound inward; the facing of
  a closed mesh is measured from a point of the mesh rather than the
  origin.

- **A torus whose hole nearly closes meshes like an open one (#146).** The
  face grid no longer spreads every row to the width of the inner equator,
  so such a torus meshes into about as many triangles as an open one, and
  `check` and `make_revolution` of it take milliseconds instead of
  minutes.

- **A coincident piece dropped for the other argument's traces to it
  (#143).** `common`, `fuse` and `cut` record a face whose piece is dropped
  for the other argument's coincident piece as modified into that piece,
  not deleted, so `trace` finds the face that now stands for its patch.

- **A move traces every sub-shape (#145).** `transformed` records every
  sub-shape of the input as modified into its occurrence in the moved
  shape, so `trace` and `is_affected` answer for faces, edges and
  vertices, alone or composed through `then`.

- **Hidden line removal splits an edge where it passes behind a
  contour.** An edge partly behind a face is drawn part visible and part
  hidden, where a straight one was drawn wholly one or the other; the exact
  path places the change on the outline, not at a sample.

- **`unify_same_domain` adds its merged faces in a fixed order,** so equal
  inputs give equal node ids from run to run.

- **The 2D wire offset measures its surviving pieces to the exact source
  arcs,** so an inward offset beside a wide arc keeps that arc's offset and
  no longer fails with "survivors close no loop".

- **A full revolution of a placed profile stores its seams where they are
  used.** A profile turned about the axis before a full revolve (to keep
  the seam off the sketch plane) got seam pcurves recorded as if it had not
  been turned, and `check` reported each as leaving its curve by the turn's
  chord. The seam records the profile's placement, as a partial turn
  already did.

- **Check reports a seam stored as plain pcurves, and healing repairs
  it.** An edge along a closed surface's seam, used once each way by the
  faces on that surface but holding its sides as plain pcurves, gave both
  uses the same side, and the next boolean failed ("a face boundary strand
  dangles"). `check` reports such an edge, and any face whose boundary does
  not close in its surface's parameter space, as broken. `fix_shape` (and
  the new `ogeom::algo::join_seam_columns`) makes the edge a seam, so a
  solid saved from 0.9.9 is repaired without rebuilding. Fixed pipes of a
  section with a closed edge, and IGES faces whose seam is sewn from two
  edges, come out with a proper seam.

## [0.9.10] - 2026-10-07

A patch release: no breaking API changes. A boolean leaves a seam between
two faces rebuilt on one surface as a seam, so a solid cut that way goes
through the next boolean (a regression in 0.9.9). A model compacts in place
with `Model::retain_reachable`, keeping the handles of what survives. A
placed dome or ball drilled through its pole closes, IGES trimmed surfaces
bounded through a pole read in exact pieces, sub-micron holes stay holes, a
surface fit goes through points its surface cannot evaluate, and nearly
coplanar planar faces that only touch refuse less.

### Added

- **A model compacts in place, keeping the handles of what survives
  (#141).** `Model::retain_reachable(roots)` drops every node, curve,
  pcurve, surface, triangulation and provenance entry the roots do not
  reach; handles into what is kept resolve as before, dropped ones no
  longer resolve, and a clone afterwards costs what is kept (a 582-face
  plate built hole by hole: 7,522 nodes from 25,952, cloned in about 1 ms).
  Adds `Arena::retain`, `Arena::next_index`, `Arena::key_at` and
  `ProvenanceTable::retain`. A native file is written renumbered when the
  model has gaps.

### Fixed

- **A boolean leaves a seam between two faces on one surface as a seam
  (#142).** A face rebuilt on its own surface beside faces set aside took
  one side's position for an edge both its neighbours stood on across the
  surface's seam, so a solid cut that way (a turned, grooved shaft bored
  along its axis) failed the next boolean ("a face boundary strand
  dangles"). Both positions are kept as the seam they are. A solid saved
  from the faulty build still carries the bad edge; building it again
  fixes it.

- **A placed dome or ball drilled through its pole closes.** A sphere
  face whose pole a curved wall passes through is re-charted when its solid
  is moved, rotated, mirrored or scaled, or lies on a left-handed sphere,
  and a mirrored or scaled whole ball drilled that way no longer gives a
  wrong volume.

- **IGES trimmed surfaces bounded through a pole read in exact pieces.** A
  boundary of a trimmed or bounded surface (143, 144) that runs through a
  pole of a surface read with it (a sphere's pole, a cone's apex, a patch's
  collapsed side) is cut at the pole in every face it bounds, so the faces
  still sew and each piece takes an exact pcurve. A half ball bounded by a
  meridian circle reads closed with the exact volume.

- **A boolean between nearly coplanar planar faces that only touch
  refuses less.** A face lying on one side of the other's plane, to
  rounding, is no longer split by the planes' line, which rounding alone
  placed across a converted corner facet meeting a pad's wall at a vertex,
  or along a facet's side edge beside the next wall. A pad into a slab
  drafted 0.00001 cuts and fuses at 45 of 60 row-count sets, against 34.

- **Sub-micron holes stay holes.** An inner ring is dropped as a slit only
  where its edges walk out and back through the same vertices and it is
  thinner than a micron; a ring that goes round (a wall poked through by
  less than a micron) is meshed open however thin.

- **A surface fit goes through points its surface cannot evaluate.** A
  surface fit (`to_nurbs_within`, a patch's degree restriction) takes the
  limit of nearby points where its surface cannot be evaluated, such as an
  offset sphere's poles, and counts how closely they agree in the error it
  states; it still refuses, saying so, where no limit exists, such as an
  offset cone's apex.

## [0.9.9] - 2026-10-07

A patch release: no breaking API changes. A boolean costs what the other
solid reaches: faces it cannot reach pass through as the same nodes, so a
small cut into a 582-face plate costs about what it costs on a 22-face one,
and a fillet's cost no longer depends on faces far from its edges.
Operations that edit in place copy shared nodes first so earlier solids
never change, and a boolean whose faster path fails starts again from the
solids as given. A planar face walks its rings about its plane, a general
draft of a closed wall closes its hinge, and a converted round between two
surfaces on its axis meets them exactly.

### Changed

- **A fillet's cost no longer depends on faces far from its edges
  (#135).** A boolean whose operand is several separate solids sets aside
  every face that misses all of them, so a fillet pass whose wedges sit on
  opposite sides of a part leaves the rest alone. The four outer edges of a
  plate with 256 holes round in under twice the plain plate's time.

- **Operations on a boolean's result leave its operands alone.** Faces a
  boolean passes through stay shared with the operand; heal, sew, fillets
  and chamfers, `remove_faces`, `unify_same_domain` and `tessellate` copy
  those shared nodes before editing and record the copies as exact copies
  in their history (`Model::unshare`, `Model::unshare_each`,
  `ogeom::algo::on_own_nodes`). `Model::widen` still grows a shared
  tolerance in place.

- **A converted round between two surfaces on its axis meets them
  exactly.** Mesh conversion rebuilds a round between cylinders, cones,
  planes square to the axis, or one of them and a sphere centred on it, as
  the torus tangent to both (a tube lip or a round-floored groove
  included), rather than a fit at a near tangency.

- **A boolean costs what the other solid reaches (#134).** Faces whose
  box misses the other solid's pass through as the same nodes
  (`History::copy_of` reports each), a plane's untouched holes are left out
  of its arrangement, and only the region near the tool is split,
  classified and sewn again; a box tree replaces the all-pairs face filter.
  A corner hole into a 582-face plate costs 1.5 to 2 times the same hole
  into a 22-face one, down from 23 times. `SolidBoundary::of_closed` builds
  a classifier without checking closure again.

### Fixed

- **A boolean that falls back to every face starts from the solids as
  given.** Where setting faces aside could not be assembled, what the
  attempt had widened of the solids' own edges and vertices stayed widened
  (sewing a failed attempt grew some by centimetres), and the boolean over
  every face then failed too. `Model::note_widened`, `undo_widened` and
  `keep_widened` journal the tolerances `Model::widen` grows and
  `Model::node_mut` hands out; the boolean puts them back before taking
  every face.

- **A boolean against a solid whose coplanar faces are still split no
  longer refuses.** An edge held by two coplanar faces of one operand was
  laid onto the other operand's face twice, which lost the hole it bounded
  ("kept pieces did not close"); a box minus a frame fused from four bars
  along its top edges is cut.

- **A general draft of a closed wall closes its hinge.** `general_draft`
  dropped the segment that closes a hinge when it came up after the chain's
  ends had grown to meet it, and could start a closed hinge half a turn from
  the seam and double back; a skinned wall drafts about a level plane at
  any height along it.

- **A planar face walks its rings about its plane.** `make_face` on a
  plane keeps the plane's normal and walks the outer ring counter-clockwise
  about it and each hole clockwise, whichever way the wires were walked, so
  a cap sewn onto a sheet faces out. `make_face_with_pcurves` still walks
  its rings as listed.

## [0.9.8] - 2026-10-07

A patch release: no breaking API changes. Booleans take imported parts with
spline fillets over small rings, each face's bounding box is kept in the
model, and a fillet reads only what is near each edge. A native snapshot
carries only the provenance its shapes reach (1.6 MB instead of 39 MB for a
582-face plate) and is written as format version 3 when its ids have gaps.
A mirrored face is shaded the way it faces, two half domes sew into an open
sheet, and a three-arc loft and a revolved dome thicken. Exact volumes
close the slits beside fitted and straying pcurves, so a converted or read
part's volume no longer depends on the face listed first.

### Added

- **`ProvenanceTable::forget` and `Model::forget_provenance`.** An
  entity's provenance entry can be dropped while its id stays issued; the
  id then resolves to nothing, like an id from another document.

- **Each face's bounding box is kept in the model (#136).** It is found
  once (`ogeom::algo::face_bounds`, `tight_bounds`, the booleans' face-pair
  filter), read in constant time with `Model::face_bounds`, and forgotten
  when the face, anything it is made of, or its surface changes in place.
  A second `tight_bounds` of a 582-face part is about 25 times faster.

### Changed

- **A fillet reads only what is near each edge (#135).** Run-out probes ask
  the solid's exact boundary instead of a fine mesh of the whole solid, a
  point on a plane face with many holes is classified against the holes
  near it, and a round's wedges whose boxes meet go in a few booleans
  rather than one each. The outer edges of a plate with 256 holes fillet in
  about 0.22 s (1.85 s before). Results are unchanged.

### Fixed

- **A pcurve that strays off its curve only near its ends gets its
  strip.** The exact volume checked a pcurve against its edge's curve
  inside the span only; one on the curve there but off it towards its ends
  left its slit open. A STEP part with such edges no longer changes volume
  with the face listed first (3.8e-9 relative before, now 2.7e-13).

- **Two half domes sew into an open sheet (#138).** Sewing two half
  domes, about opposite axes or mirrored, joined their half rims as one
  edge and closed the sheet when the arc's centre sat a little off the
  axis. A face holding one edge at two placements keeps its meridians on
  their own columns and its pole as one vertex when sewing rebuilds it, so
  only the meridian seams join.

- **A loft through three arcs and a revolved dome thicken (#140).** A
  smooth loft whose sections are two-piece curves is built without a crease
  along their middle knot, so `make_thick_sheet` and `offset_sheet` fit
  its offset; lofts through such sections get a slightly different surface
  between the sections. A revolved circle in a plane through its axis
  offsets exactly, and offsets no longer go wrong at a pole.

- **A native snapshot carries only the provenance its shapes reach
  (#137).** `native::write` and `compacted` wrote every provenance entry
  the model ever made (39 MB for a 582-face plate cut hole by hole); they
  keep the entries of the written identities and their ancestry: 1.6 MB,
  written in about 17 ms. Entity ids are unchanged, and a file whose ids
  have gaps is written as format version 3, which older readers refuse.

- **A face under a mirroring placement is shaded the way it faces
  (#139).** `triangulate_face` and the whole-shape mesh gave normals
  against their triangles' winding under a reflection; they agree under
  any placement.

- **Exact volumes close the slits beside fitted edges.** Where a pcurve
  stands off its edge's curve by more than a confusion distance, fitted
  edges within their stated tolerance included, the strip between them is
  measured and closed onto the curves' ends at every corner. A converted
  part's volume no longer depends on the point its moments are taken from
  (5e-5 relative before, now 5e-10), and the torus drill commons measure
  within 1e-13 of the exact core (1.5e-9 before).

- **Booleans take imported parts with spline fillets over small rings
  (#133).** A projection onto a surface stopped when its residual times the
  chart's speed fell under the confusion, which on a ring of a sixth of a
  millimetre left a pcurve microns short of the seam, so the face's
  boundary read open ("a face boundary strand dangles") and every boolean
  against the part failed; it converges along the surface itself. Torus
  faces whose trim sits a period off in both directions are cut, a piece of
  an edge whose pcurve lags its curve states the lag as its tolerance, and
  a tool standing in a gap of a part with microns of edge slop leaves it
  whole instead of refusing.

- **A sampled fit states no less than its widest miss.** Surface and curve
  fits made at the true geometry's own parameters (spline conversion,
  degree restriction, sweeps, lofts, fillings, pipe sections) read their
  error at samples and fixed fractions of each span, which could sit under
  the true miss (an offset drum stated 7.83e-5 against 8.45e-5); the
  highest readings are followed to the top of their peak, so the stated
  error holds.

## [0.9.7] - 2026-10-07

A patch release: no breaking API changes. A reversed face's holes are no
longer taken for its outer wire (STEP and IGES export marked a hole as the
outer bound), and STEP and IGES write a reversed circle, ellipse or trimmed
curve at its own points. IGES trimmed-surface boundaries round-trip with
reversed edges, sewing matches closed edges in sense, and `make_solid`
turns a consistent shell facing the wrong way. Fillets over a sphere's pole
build and mesh closed at every turn, a band round crossed bores meshes to
its volume, and points beside an edge near a pole classify in or out.

### Added

- **`CircleCurve::forward` and `EllipseCurve::forward`.** The same curve
  running forward, its frame turned half a turn about its x axis.

- **`Model::outer_wire`.** A face's outer wire, read from the stored order
  (outer first) whichever way the face is used; a reversed face's walk
  lists its holes first.

### Fixed

- **Points beside an edge near a sphere's pole classify in or out.** The
  boundary band was carried into the chart with one scale for both
  directions, the smaller tangent's, so half a millimetre from a pole it
  reached twenty times too far toward it and points beside a rim read as on
  the boundary; the band is measured in space along each chart direction.
  A fillet whose rim's middle is its point nearest the pole seats its ball
  at another point of the rim and builds.

- **Fillets whose ball passes over a sphere's pole.** The sphere leg's seam
  ran along the meridian from the march's first station, through the rim
  it bounds, so the wedge meshed open (5.7 to 64 against 4.44 cubic
  millimetres) and the fillet was refused at most turns of a bore about the
  ball's axis. The rail is split where the seam can run clear of the rim:
  the blend meshes closed and fillets at every turn measured.

- **A reversed circle or ellipse is written and imaged as itself.** STEP
  wrote a circle or ellipse running backwards as its forward conic, and the
  reader took the other arc between its vertices; and its exact pcurve ran
  the forward way on every analytic surface, so a face bounded by one on a
  plane failed `check`. Both follow the curve's own direction.

- **A fillet band round crossed bores meshes to its volume.** The legs of a
  band round a closed loop where two bores cross had their two rings imaged
  a whole turn apart in the wall's chart, so the wedge meshed to -304
  against an exact 8.4 cubic millimetres; both rings are put on one branch.
  `chain_wire_branches` places a wire by points along its images, so a ring
  closed by one edge is moved onto the first wire's branch.

- **`make_solid` turns a consistent shell facing the wrong way.** A closed
  shell whose faces all agree but face in made an inside-out solid; the
  outer shell is turned to face out and a void's shell to face into its
  cavity, so an imported shell with every face reversed becomes a valid
  solid. Sewing keeps such a shell as given, since it may bound a void, and
  builders that want a shell facing in use `Model::add_solid`.

- **IGES writes a reversed trimmed curve at its own points.** A trim running
  against its basis was written over the basis's forward range, so a
  reversed half circle came back as the other half.

- **IGES trimmed-surface boundaries round-trip with reversed edges.** Each
  boundary piece was written in its curve's own direction, so a sheet or
  shell with an edge used reversed (any reversed face) was refused on
  reading; pieces are written the way the wire walks them, a circle or
  ellipse running backwards along its curve at its own points, and the
  reader chains a boundary from either end of its first piece, so decks
  written the old way read too.

- **Sewing matches two closed edges in sense by their quarter points.** A
  cap's rim against a wall's coincides end to end both ways round, so
  sewing could read them as walked the same way: a cylinder read from
  trimmed surfaces came back with its caps facing in.

- **A reversed face's holes are no longer taken for its outer wire.** STEP
  export of a reversed face with holes marked a hole as
  `FACE_OUTER_BOUND`, and IGES export of one as a trimmed surface wrote a
  hole as the outer boundary, so other readers were told the hole bounds
  the face; a fixed-section pipe of such a face faced its side walls into
  the material and stored a cap's hole as its outer wire. Each reads the
  stored order. Point classification uses its lazy hole rings on such
  faces too, with answers unchanged.

## [0.9.6] - 2026-10-06

A patch release: no breaking API changes. It closes issues #127 to #132:
thickening a loft between arcs of different radii, a worm and a bearing
combined on one plane, a bore removed from or moved in an extruded plate,
the middle path of a tube with a sharp corner, and a filleted bracket
measured exactly after unifying its faces. Every producer, the STEP and
IGES exchange and sewing now walk each edge once each way, and `check`
names an edge both its faces walk the same way. A widened edge carries its
vertices, and sewing orients the faces it joins by their walk.

### Fixed

- **A worm and a bearing standing on one plane combine (#128).** An exact
  section on a curved face met a fitted boundary edge only within the
  confusion distance unless the edge was a line or conic, and a near miss
  running a long way along the curves was read as a touch; the section
  meets any boundary edge within that edge's tolerance, and a near miss is
  dropped only where the curves are tangent. A grooved worm against a ring
  whose flat faces cut its lands closes for common, cut and fuse at every
  height tried.

- **Moving a bore sideways, and removing it from an extruded plate (#130,
  #129).** A circle's exact pcurve on a cylinder, cone, sphere or torus
  covered one turn from the circle's origin, so an arc running past the
  origin, as a moved bore's split ring does, failed to evaluate; it covers
  the edge's range. A bore's wall removed from an extruded plate, whose
  base is used reversed, makes the plate whole.

- **A unified face keeps its curved edges' pcurves (#132).**
  `unify_same_domain` gave every boundary edge of a merged plane a chord
  between its ends as its pcurve, which stood 0.9 off a fillet's arc, so
  the face left the exact volume path; an edge takes the pcurve it already
  carries on a merged face, turned into the kept chart. A filleted bracket
  refined by the unify measures exactly, to the same volume as before it.

- **Thickening a loft between arcs of different radii (#127).** Such a
  surface is smooth but parameterised only C0 along its middle knot, and
  the offset's cubic fit stalled at about 5.5e-3 there; the fit keeps a
  triple knot at every knot where the base surface is C1 or less, and
  fits to 1.9e-7.

- **The middle path of a tube with one sharp corner (#131).** A section
  square to one leg reached past the mitre into the other and the march
  read it as turning back; the path is marched from both ends, each leg cut
  back from the corner, and joined at the corner as two edges. A tube with
  two or more sharp corners is still refused.

- **Sewing orients the faces it joins by their walk.** Each edge two faces
  share ends up walked once each way, open sheets included, keeping the
  orientation most faces already have; a closed shell is then made to face
  out. A shell that cannot be oriented (a half-twisted band) keeps its
  faces, and the new `Sewn::edges_walked_one_way` names the edges that
  show it.

- **A widened edge carries its vertices, and a widened face its edges.**
  Marched fillet bands, offset rebuilds, splitting, dividing,
  same-parameter, sweep rails, converted edges, IGES pcurve gaps, and faces
  rejoined at a seam or rebuilt on a ball chart widened an edge (or face)
  in place and left its vertices tighter, which `check` reports broken; they
  widen through the cascade that carries them along.

- **Fillet wedges walk each edge once each way.** Rim arc wedges and
  marched blend wedges (closed, run-out, pinched, face blends) built faces
  whose rings ran the same way as the face beside them; each keeps its
  material on the left of its rings. `make_band_between` walks its lower
  ring forward whichever order the rings are given.

- **Every other producer walks each edge once each way.** Faces built by
  the mesh converter, by baking a reflecting placement or restating a
  surface, by seam healing, ring re-anchoring and surface recognition, by
  face removal, draft, face offset and move, by sheet rounds and by
  straight and rim blend wedges, and faces written to IGES, keep their
  material on the left of their rings. A reversed IGES face was written
  with its loops under the face's sense.

- **Lofts, pipes and sweeps walk each edge once each way.** Lofts between
  polygons, pipes along part of a circle, skinned lofts and pipes, pipes
  through sections or along law-driven frames, helical sweeps and pipe
  shells built caps and walls whose rings ran the same way as the face
  beside them; each face keeps its material on the left of its rings.
  Geometry and volumes are unchanged.

- **Reversed STEP faces read with their loops as the file walks them.** The
  reader stored a reversed face's loops as written and then turned the face
  over, so every such face walked its loops backward; every edge of a
  third-party file now reads walked once each way. Files written by earlier
  ogeom with backward loops on such faces are repaired on reading.

- **Mesh conversion places a tangent seam only through its vertices.** A
  chord seam between two faces meeting almost tangentially was built even
  where its vertices stood far off a surface; a surface standing off its
  vertices now falls back to facets. A curved face that absorbed facets and
  ran over its neighbours gives them back before anything is faceted, so
  every corpus part converts with no faceted region.

- **Sewing turns over a face that faces into the closed shell it joins.**
  A box sewn from sheets with a lid filled the wrong way round made a solid
  `check` called reversed; where two faces walk a shared edge the same way,
  the faces the orientation probe finds facing in are turned, provided every
  shared edge is then walked once each way and the shell encloses a
  positive volume. Faces that already face out, and open shells, keep the
  orientation they were given.

## [0.9.5] - 2026-10-06

A patch release: no breaking API changes. It closes issues #120 to #126:
fillings over crossing semicircles, a four-line saddle and a tube's
tangent cap; a sweep surface through a tangent join; a folded sheet
thickened with a mitre; a sewn box measured exactly; and a STEP file that
0.9.0 read clean reading clean again. STEP carries each edge's pcurves both
ways, sewing keeps a reversed face's rings, swept solids walk each edge
once each way, and the mesh converter converts every NIST solid in the
corpus with no curved face faceted. Fillets the stress run refused now
round, and a boolean on a large solid leaves less work on untouched faces.

### Added

- **`fit_surface_grid_at_with_knots`.** A fixed-parameter grid fit that
  starts from given interior knots and keeps them.

- **`check` names an edge both its faces walk the same way.** In a shell,
  each edge between two faces is walked once each way when every face keeps
  its material on the left of its rings; an edge both faces walk the same
  way is a suspect finding, in closed shells and sheets.

- **A folded sheet thickens into one solid (#124).** `make_thick_sheet`
  joins faces meeting at a crease with a mitre on either side, where the
  crease runs from border to border and those borders stand square to it
  (a profile extruded square to its plane); other creases are refused by
  name. `offset_sheet` still refuses creases.

- **STEP carries each edge's pcurves.** Edges are written as surface curves
  with a pcurve on each face (a seam curve with both sides on a closed
  face), wherever the file can state them exactly in the curve's parameter.

### Changed

- **Less overhead on the faces a boolean leaves alone.** Edge lookups while
  rebuilding kept faces are hashed instead of scanned, widening a node's
  tolerance and querying the vertex bins allocate nothing, and
  `History::followed_by` composes a step that records nothing without
  copying the whole history. A small drill or boss against a large solid
  takes about 8% less CPU, with identical results.

- **One fitter for curves measured between their samples.**
  `ogeom_algo::traced::fit_traced` and `fit_traced_2d` call
  `ogeom_geom::fit::fit_curve_sampled`: the least-squares fit first, the
  spline through every sample where it misses, each interval checked at its
  eighths. Traced fits come out closer and mostly faster (a helix projected
  onto a plane at a fifth of its target instead of all of it, about seventy
  times faster). Added `fit_curve_2d_sampled`, which also takes closed
  curves, and the budgets `SAMPLED_SPANS` and `SAMPLED_LOOP_SPANS`.

### Fixed

- **A STEP file checks clean where a circle's range runs past one turn
  (#126).** The reader measured such a circle's pcurve gap as nothing, so
  the edge kept a tighter tolerance than its pcurve; and the nearest point
  of a closed curve beside either end was missed, so a repair could add a
  broken edge. A file 0.9.0 read clean reads and repairs clean again.

- **Mesh conversion depends less on triangle order.** A fillet torus a few
  facets round is recognized whichever triangle seeds it, a sphere fitted
  with its centre far away no longer swallows neighbouring regions, and the
  facets beside a face reaching past its triangles become fans before the
  face is faceted. Corpus parts come back on fewer faces.

- **A sewn box of planes measures exactly (#125).** `volume_properties`
  sent a shell holding a face whose flag points into the material (a box
  sewn from six sheets with a lid filled facing in) to the mesh, and its
  properties came back approximate. Faces walking a shared edge the same
  way are tied as turned against each other, each tied set is probed
  against the solid, and a face facing in is counted turned back: the box
  measures 500 and 400 exactly.

- **A sweep surface follows a line into a tangent arc (#120).**
  `make_sweep_surface` skinned the whole spine at once across the
  curvature jump at a tangent join, and its error stalled just above the
  target; each spine edge is skinned on its own, so the sheet has one face
  per profile edge and spine edge.

- **Fillings take saddles, crossing semicircles and tube caps (#121, #122,
  #123).** `make_filling_n` read rounding noise along a straight side as
  the loop crossing itself, refusing a four-line saddle; a point within
  rounding of a line counts as on it. A loop whose corners look smooth from
  the plane it spans (two semicircles meeting square in crossing planes),
  and tangent sides whose supports stand square to that plane (a cap on a
  tube), are filled by a patch free in all three coordinates over a chart
  drawn from the loop, without folds; a tube's G1 and G2 cap leaves the rim
  along the wall. A corner where no surface can be tangent to both supports
  is refused naming the angle, and a steep support's refusal states its
  angle and the limit.

- **STEP keeps a file's pcurves, and sewing keeps a reversed face's
  rings.** The reader derived every edge's pcurves itself; it keeps the
  file's pcurve where it follows the curve within twice the projection's
  distance, so a part comes back bounded as it was written (a converted
  part's 1842 pcurves come back as written and its volume unchanged).
  Sewing rebuilt a reversed face, such as a boolean cut's tool faces, with
  its rings walked backward, the same way as the face beside it; it keeps
  them as stored.

- **Swept solids walk each edge once each way.** A prism or revolution
  over a profile edge walked backward, a profile ring wound against its
  role, a tapered prism round a hole, a baked reversed face and a
  revolution band given its upper ring first each built faces that walked
  their rings the same way as the face beside them. Each face keeps its
  material on the left of its rings.

- **Edge tolerances cover a pcurve's peaks between samples.** An edge's
  stated tolerance was measured at 256 samples, so a piece cut from it
  could stand a few microns past it; the widest peaks are searched for
  their tops.

- **Fillets the stress run refused.** A fillet on the loop where two bores
  cross decided each face's side by a step so short it landed inside the
  face on both sides; it steps the whole setback where the short step
  cannot tell. Rounding an edge and then one at its end where a drum takes
  over took the first blend's cap edge for a piece of the second; a piece
  must lie between the same surfaces as the edge. A boolean's material test
  on a flat face with a rotated chart used a ray level with a corner and
  dropped the strip between two cuts; it uses the slanted ray its trims
  use. Thirteen stress fillets over three seeds that were refused now
  match their closed forms.

- **The corpus's last faceted curved faces convert on their surfaces.** A
  fillet torus round a rounded corner, meshed a few facets round, comes
  back as one torus instead of spheres; a disc capping a chamfered end is
  no longer taken by the chamfer's cone; a drill point whose rim is arcs of
  one circle (under a hex socket) closes at its apex as one cone face; a
  flat face meeting a curved one below the crease angle no longer turns it
  into a kinked sweep; a single facet with two sides on curved faces is
  built as a ruled wedge; and a run of facets between two single facets is
  no longer put on a cylinder. Every NIST solid in the corpus converts with
  no curved face faceted, and a socket screw comes back as the 20 faces of
  its source instead of 470.

- **The mesh converter facets a face it cannot build.** `solid_from_mesh`
  failed outright when one face could not be built, meshed or measured (a
  degenerate spline sliver failed with "cannot normalize a vector"); the
  face is faceted and named with the new `FallbackReason::BuildFailed`. The
  converter also probes every face for orientation as `check` does, so a
  planar face turned in among planes is rebuilt from its triangles
  (`FallbackReason::TurnedIn`); `curved_faceted` and `fallbacks` can name
  planar regions too.

- **Fewer converted faces fall back to facets.** A chain of three vertices
  between faces meeting all but tangentially is threaded with a parabola
  instead of faceting both faces; a single facet whose seam with a curved
  face cannot be placed, or that folds it, is built as a fan onto it; a
  drill point whose rim starts off its seam is turned to meet it; a small
  curved region joins a neighbour of another kind it lies on; and coarse
  pieces of a fillet torus fitted as spheres are put back on the torus.
  NIST parts convert with fewer faceted regions (several with none), a
  corpus part that converted invalid is valid, and the truth bench's mean
  volume error falls from 3.5e-6 to 7.1e-7.

## [0.9.4] - 2026-10-06

A patch release: no breaking API changes. Exact mass properties weigh a
void turned inside out as a void and close the sliver a chord leaves
beside its arc, so converted parts bounded by chords measure what their
tessellation does. The mesh converter rebuilds drill points, narrow rounds,
drafted corners, facets folded against curved faces and tangent chords as
exact faces: NIST parts convert to far fewer faces, one from 4104 to 1997
and another with no face left faceted. A pipe shell's wall holds between
rows at a curved corner, and the area of a spline face that folds over
settles at once.

### Changed

- **The area of a spline face that folds over settles at once.** Where a
  trimmed spline face's normal nearly vanishes inside it, the exact area
  graded its panels towards the fold instead of halving them all: a wavy
  wall drafted 0.1 takes 0.07 s instead of about 5 s, and its area is more
  accurate. Other areas and every volume are unchanged.

### Fixed

- **Exact mass properties close the sliver a chord leaves beside its
  arc.** An edge taking a straight chord or another closed-form curve for a
  curved section, stating the sag as its tolerance, left the strip between
  the face's pcurve and the curve open wherever it was wider than a
  thousandth, faces measured as chart rectangles ignored it altogether, and
  areas counted strips by walk direction. Such strips are closed as wide as
  the edge's tolerance allows and count as added surface: a converted part
  bounded by chords measured 1.3e-3 under its volume and its blend tori
  1.8% under their area; it now matches its fine tessellation.

- **Drill points, stand-in tori and tangent chords convert exact.** A cone
  closing at its apex inside one rim (a drill point) had no layout and fell
  to facets; it is built like a sphere's cap, its seam running from the
  apex to the rim. A strip of a cylinder fitted as a torus thousands of
  millimetres in radius comes back as the cylinder the torus stands in for.
  An edge between two curved faces meeting almost tangentially is traced
  through points placed on both surfaces. NIST parts convert to far fewer
  faces (one from 4104 to 1997, its volume error from 1.3e-3 to 3.7e-5;
  another from 345 to 155 with no face faceted).

- **Narrow rounds and drafted corners convert to their own surfaces.** A
  narrow round meshed only a few rows across came back as near-flat
  spheres hundreds of thousands of times its size, which won the chord
  tie-break over the round's exact cylinder; a sphere or torus missing the
  points by more than a closer fit no longer takes part in it. The corners
  of drafted walls, meshed too few facets round to fit on their own, come
  back as cones tangent to the walls, so the fillets at their foot no
  longer fold under them. A NIST part converts to 927 faces instead of
  1460 with every torus built, and most of the others to fewer faces.

- **A facet folded against a curved face is built as a fan.** Mesh
  conversion faceted a curved face where a single-triangle facet beside it
  met it almost tangentially and folded under it; such a facet is a ruled
  fan from its far corner to the seam on the curved surface, so the curved
  face stays curved. NIST parts convert to fewer faces (687 instead of 726,
  4115 instead of 4138).

- **A void turned inside out is weighed as a void.** `volume_properties`
  sent a solid whose void shell had every face reversed to the mesh, which
  added the cavity to the volume (a block of 1000 with a void of 64 measured
  1064, not 936). A whole shell facing in is turned over in the closed
  form, and the mesh path probes every face of a solid with several shells
  as `check` does.

- **A pipe shell's wall at a curved corner holds between rows.** A curved
  leg meeting another leg at a crossing was fitted through its rows only: a
  disc down an arc turning onto a straight leg stood 7.3e-3 off at a 1e-4
  target, and a square section through the same corner was refused. The
  run is checked between its rows (7.7e-5 off) and the square section
  builds.

## [0.9.3] - 2026-10-06

A patch release: no breaking API changes. More silent wrong answers are
closed. `check` probes every face of a solid for orientation, and an exact
volume settles each face's orientation before trusting it, so a face
turned into its solid can no longer pass as valid or give a wrong volume.
Sweeps hold their walls between stations along the path, drafts and loft
end patches hold their geometry between samples, and a helical sweep
through part of a turn ends on an outward cap. A trimmed spline face's
exact volume is up to a hundred times faster, and converted parts keep
their straight-threaded seams.

### Changed

- **A trimmed spline face's exact volume settles on what it returns.** The
  volume integral refined until the face's area settled too, which on a
  nearly folded spline face took millions of samples; it stops when two
  runs agree on the volume and its moments. A drafted wavy wall's volume
  check drops from 10 s to 0.1 s; inertias move within the integrator's
  stated 1e-10 agreement.

### Fixed

- **An exact volume never takes a face's sign on trust.**
  `volume_properties` read a face whose orientation it could not compare
  along its edges (a cone apex, a face placed inside the shape, a separate
  shell) as agreeing, so a face turned into its solid gave a wrong exact
  volume: a drum with a nearly flat cone cap turned in measured 784.09
  where it holds 786.71. Each group of faces the edge walk cannot tie is
  probed against its solid, any of its faces able to settle it (the broad
  face of a very thin sheet passes the question to its rim), and the volume
  falls back to the mesh only where a group faces in or none of its faces
  can be decided. A prism's placed top is compared too.

- **A converted part keeps its straight-threaded seams.** Mesh conversion
  threw away every seam it had threaded straight when any facet beside one
  faced into the material, and the rebuild faceted curved faces for folded
  seams; the threading stays and only the curved faces beside such a facet
  are faceted. A NIST part converts to 726 faces instead of 977, with 144
  curved faces where it kept 127.

- **`check` probes every face of a solid for orientation.** A single face
  turned into its solid passed where the faces' flags could not be compared
  along their edges (a nearly flat conical or domed cap, a face whose
  boundary middle lies outside it): a loft whose end patch faced in kept
  half its volume and read valid. Every face is probed off both its sides,
  the step shortened on solids thinner than it. The mesh converter asks the
  same probe and falls back to facets where a converted face folds under its
  neighbour: three NIST parts come back with a few more faces (ftc_10 from
  704 to 977) and their volumes unchanged within tolerance.

- **Sweeps hold their walls between stations along the path.** Pipe
  shells, law pipes and helical sweeps were checked along the profile but
  not between stations: a disc swept down a spline turning a quarter within
  one span stood 1.1e-2 off at 1e-4, a guided law loft 6.4e-3. They are
  fitted to the swept profile between stations as well. A general draft
  holds its wall to the turned rulings across its whole height (2.0e-3 off,
  now 1.8e-5), and a skinned loft's non-planar end patch is held to its
  cone (9.6e-2 off and invalid, now valid within 6e-5) and faces out of the
  solid.

- **A helical sweep through part of a turn ends on an outward cap.** A
  sweep that was not whole turns ended on a cap facing into the solid. An
  auxiliary-law pipe whose guide starts on the first station's plane no
  longer turns its start section.

- **Drafting an extruded spline wall holds it between stations.** The
  drafted wall was fitted through 65 stations and checked there only: a
  wavy wall drafted 0.1 at 1e-4 stood 2.5e-3 off its turned rulings. It is
  fitted to the rulings between stations as well (3.1e-5), and extends far
  enough past the profile's ends for a neighbour meeting it at a slant.

## [0.9.2] - 2026-10-05

A patch release: no breaking API changes. Approximations that measured
their error only at the points they were fitted through, and so reported a
tolerance met while straying far between them, are measured between their
samples, refined where they miss, and state the worst error found: spline
conversion of surfaces (a sphere stood 0.28 off), sweeps, lofts, pipes and
fillings (up to 0.18), heal boundaries, plane projections, IGES trims and
offsets, canonical recognition, and fillet and boolean images. New
`fit_surface_sampled`, `fit_curve_sampled`, `fit_trace_sampled` and
`ogeom_algo::traced` fit known geometry this way.

### Fixed

- **Sweeps, lofts and fillings hold their tolerance between samples.**
  Skinned pipes, skinned lofts (aligned, closed, apex and law forms),
  fillings, swept strips (pipe shell, cornered loft, helical sweep) and
  closed and multisection pipe rings were fitted through samples and
  reported the tolerance met while straying up to 0.18 between them. They
  are fitted at the true geometry's own parameters, checked at quarter
  points (eighths along curves), refined where they miss, and state the
  worst measured. `normal_projection` measures between its stations and
  widens the edge to the error it states. `fit_surface_sampled`,
  `fit_curve_sampled` and `fit_trace_sampled` fit known geometry this way.

- **Fitted curves are measured between their samples.** A heal boundary
  moved onto its face, a curve projected onto a plane, a band connector, an
  IGES trim given only in parameter space or a varying offset, and a
  boolean section's image were fitted through samples and their error read
  at those samples only: a moved boundary stood 0.1 off its face while
  reporting 3e-12, a projection 2e-2 off while stating under 1e-5, a
  301-corner IGES trim lost its corners by 6e-2. They are fitted at their
  source's own parameters (`ogeom_algo::traced`), checked at every
  interval's quarters, refined where they miss, and their edges state the
  worst measured.

- **Canonical recognition checks the whole face.** `recognize_surface`
  verified a candidate plane, sphere, cylinder or cone on an interior grid
  that never reached the face's border strip, so a patch bumped near its
  edge read as a plane; it verifies over the whole chart, and
  `canonical_simplify` states the certificate on the rebuilt face and its
  edges. A line or circle replacing a curved edge states its deviation, and
  `fix_face_pcurves` widens an edge's vertices with the edge.

- **Surfaces fitted as splines hold their tolerance between samples.**
  `to_nurbs_within` and a patch's degree restriction checked the fitted
  surface only at its grid samples and reported the tolerance met: a
  sphere converted at 1e-4 stood up to 0.28 off between them, an offset
  drum 1.29e-3. The surface is fitted at its own parameters, checked at
  every span's middle as well, and refined where it misses; every case
  measured now holds the tolerance asked. A surface fit takes tens of
  milliseconds where it took about one.

## [0.9.1] - 2026-10-03

A patch release: no public API changes. `check` now compares every edge's
pcurves with its curve, not only where an edge claims `same_parameter`, and
every producer that understated its edges' tolerances states them first:
the mesh converter (a plane's arc pcurve stood thousands of millimetres
off), the STEP and IGES readers (which also cut an edge through a pole into
exact pieces, keeping its colour and PMI on them), `fix_shape`, booleans
under a fuzz, NURBS conversion and projection onto offset surfaces. STEP
export writes edge colours and edge PMI.

### Added

- **STEP export writes edge colours and edge PMI.** A coloured edge's
  `EDGE_CURVE` gets a curve style, and dimensions, tolerances, datums and
  datum targets whose features are edges are written as shape aspects
  naming the written `EDGE_CURVE`s, so both survive a write and read round
  trip; edge features were dropped before.

### Changed

- **`check` compares every edge's pcurves with its curve.** It compared
  them only on edges claiming `same_parameter`, which few producers set, so
  an edge whose pcurve stood thousands of millimetres off its curve read
  valid. Every pcurve, seams included, is now measured along the whole edge
  against the edge's tolerance: at the same parameter where the edge claims
  `same_parameter`, otherwise against the nearest point of the curve. The
  producers that understated their edges' tolerances (the mesh converter,
  the STEP and IGES readers, `fix_shape`, booleans under a fuzz, NURBS
  conversion) state them first.

### Fixed

- **Projection onto an offset surface.** Projecting a point onto an offset
  surface returned the search's starting guess, an offset having no second
  derivative; it now steps on first derivatives there. Faces restated as
  offsets keep their pcurves on their edges (an offset drum's stood 0.11
  off). `to_nurbs_within` keeps the fitted surface's iso line as an edge's
  pcurve only where it follows the edge within the fit target, and every
  re-derived edge states how far its pcurves stand off its curve.

- **Fuzzy booleans keep their pieces' stated tolerances.** Under a fuzz,
  pieces of operand edges, contacts and sections dropped tolerances below
  the fuzz to the minimum, and the sew let twins within the fuzz through
  unstated; pieces keep their source's tolerance, and a contact's piece
  states how far its pcurve on the other surface stands off.

- **Edge attributes follow an edge cut at a pole.** On STEP import, a style
  or a PMI shape aspect (dimension, tolerance, datum) that names an edge the
  reader cut at a pole lands on every piece of that edge, and edge styles
  are read at all; the cut edge is left with nothing referring to it.

- **Edges through a pole read in exact pieces.** The STEP and IGES readers
  cut an edge whose curve runs through a pole of a face it bounds (a
  sphere's pole, a cone's apex, a patch's collapsed side) at that pole, so
  each piece takes an exact image in the chart. A meridian circle through
  both poles of a sphere reads at the confusion tolerance, where its fitted
  pcurve stood millimetres off with a warning; the hemisphere's area is
  exact.

- **Read edges state where their pcurves run.** The STEP and IGES readers
  raise each edge's tolerance to how far its pcurves, exact or fitted, stand
  from its curve, measured densely and searched near poles and ends: a line
  microns off its plane, a sphere pcurve straying near a pole, and spheres
  and converted solids read back stated 1e-7. Beyond a millimetre the
  readers warn instead. `reduce_tolerances`, and so `fix_shape`, no longer
  shrinks an edge below how far its pcurves stand from its curve; it
  measures every pcurve, seams included. Five more drills through a NIST
  part cut in the stress run.

- **A converted plane face's pcurves stay on their edges.** An arc fitted on
  a neighbouring face that only crosses the plane along a few millimetres
  was imaged on the plane as a full-radius circle thousands of millimetres
  off its edge, which `check` did not compare. A plane keeps the closed-form
  image only where it lands within the edge's tolerance and otherwise
  projects the curve; every converted edge's tolerance covers how far its
  pcurves stand off its curve, measured at 256 points.

## [0.9.0] - 2026-10-03

A minor release, since the API changes: `MeshSolidReport` gains
`fallbacks`, `SurfaceKind` and `RegionRefusal` gain variants, and
`FillSide::edge` is the edge bounding the filling. It closes issues #116 to
#119: placed edges and sheets wherever fill, blend, extend and bake read
them, filling sides chained by coincident ends, sections of different edge
counts in lofts and ruled sheets, and offsets of fitted sheets. Several
silent wrong answers are fixed: drills touching a torus or a ball, drills
through a dome's pole, arcs bounded short on a turned frame, a sewn sliver
moving a corner, edges stating a smaller tolerance than their pcurves keep,
and exact volumes that depended on where a loop's integral began. The
stress harness now judges volumes exactly to a millionth. The mesh
converter converts sphere caps with holes, tori whose holes leave no seam
free and patches tangent to their neighbours; booleans on large solids and
exact mass properties run several times faster.

### Added

- **Placed edges and sheets everywhere they are read (#116).**
  `baked_shape`, `to_nurbs`, `to_nurbs_within`, `restate_geometry` and
  `general_transformed_shape` rebuild shells, faces and compounds of them,
  not only solids, history mapping each face and edge occurrence to its
  twin; `baked_where_placed` bakes only where something is placed.
  `make_filling_n`, `make_blend_surface`, `extend_face` and
  `fillet_sheet_edges` take placed edges and faces where they stand, and
  `sew` joins a placed face to faces built from it. A prism's walls, a
  floor and a filling across their top edges sew into a closed shell of
  volume 500.

- **Filling sides chain by coincident ends (#117).** `make_filling_n` joins
  sides whose ends coincide within tolerance, a shared vertex node taking
  precedence; where neighbouring sides carry their own vertices the face
  is bounded by new edges on the given curves (`FillSide::edge` is the
  edge bounding the face), and the history records them.

- **Sphere caps with holes convert to one face.** A sphere bounded by one
  rim and holes (a dome drilled off its pole, a hill on a hemisphere) comes
  back as one sphere face, its pole on the face and its seam clear of the
  holes; a sphere whose hole holds the pole comes back as one face between
  its two rims. A hill on a hemisphere converts to 3 faces instead of 5304,
  and a socket screw's button head in the corpus from 1272 faces to 215.

- **Patches meet a curved neighbour on its surface.** Where a fitted patch
  runs out tangentially into a cylinder, cone, sphere or torus, the seam is
  threaded in that surface's own coordinates, cornered where it turns, and
  mapped back onto it: a hill on a round bar meets the bar within 6e-10
  instead of 5.5e-5. A large canonical face beside a smooth region bounds
  that region instead of being tried with it as one patch.

- **Tori whose holes leave no seam free convert whole.** A torus whose
  holes cross every parallel or every meridian comes back as one exact
  toroidal face, its seam threaded through the holes and the full circle
  the other way as its second seam; a ring grooved all the way round its
  axis or its tube comes back as an exact band. Volumes match the source
  within 1.4e-7.

- **Mesh conversion steps hold or say why not.** `MeshSolidReport::fallbacks`
  names each curved region the build faceted, with a `FallbackReason`.
  Regions changed by `merge`, `split` or a free `fit` go through the
  tangent passes (rounds, blends, corner balls) the way recognized regions
  do, and one- or two-triangle facets a step made are no longer absorbed
  into curved neighbours. `fit` takes `SurfaceKind::Patch` (refused as
  `NotADisk`, `TooNarrowForAPatch` or `PatchDoesNotVerify`), and `merge`
  and `split` fall back to a verified patch where no other surface fits.
  `solid_from_mesh` results are unchanged. The new report field and enum
  variants break code that builds the report or matches these enums
  exhaustively.

- **Free-form regions tangent to a plane convert to a patch.** A smooth
  region that runs out tangentially into a plane comes back as a fitted
  patch bounded by a seam lying on that plane, instead of taking the plane
  with it and falling to facets: a hill tangent to a flat top converts to 7
  faces instead of 1153. A plane of 32 triangles or more bounds a smooth
  region, and a region whose patch does not verify over the mean-value map
  is charted by projection onto its mean plane.

- **Fine meshes of free-form surfaces convert to one patch.** Where
  recognition cut a smooth free-form area of a fine mesh into small
  spheres, cones, cylinders and tori, `solid_from_mesh` fits the area and
  the pieces it encloses as one B-spline patch, kept only if the patch
  verifies. A bumped plate meshed at 80 or 120 cells a side comes back as
  six faces instead of twelve to twenty-four thousand; a canonical face
  tangent to anything outside the smooth area is never taken.

### Changed

- **The stress harness judges volumes exactly, to a millionth.** Volume
  identities (cut and common adding up to the part, fuse and common to both)
  are checked on exact volumes with a slack of 1e-6 relative instead of
  tessellated volumes with 2e-3, which had hidden errors of a few tenths of
  a percent. A result the exact integral refuses is measured at a 2e-4
  chord with the chord's error added to its slack, and listed.

- **Exact mass properties are several times faster on long pcurves.** Short
  boundary panels take as few Gauss points as their error bound allows,
  inner integrals on analytic surfaces are not refined again between runs,
  and a run is summed as it goes instead of held. A coarsely meshed hill on
  a hemisphere and a torus threaded through six holes convert about four
  times faster; volumes are unchanged to within 1e-12.

- **Booleans on large solids cost closer to what the tool touches.**
  Untouched faces pass through the sew as they stand, the classifier draws
  a plane's holes only near where it asks, and a cut face's holes that
  nothing comes near are handled apart from its arrangement. A small drill
  or boss against a 900-hole plate or an 1870-face converted mesh takes
  about 1.8 times less CPU time, with identical results.

### Fixed

- **A sewn edge states where its pcurves run.** Where the sew merged a
  boolean's fitted section with an existing circle it nearly follows, the
  kept edge stated the circle's 1e-7 while the copied pcurves stood up to
  1.2e-4 off it. The kept edge now raises its tolerance, and its vertices',
  to the distance measured, and the exact volume still integrates the strip
  between pcurve and curve on an analytic edge.

- **G1 fillings no longer sag in the middle.** `make_filling_n` accepted
  the first control net whose sides met the tolerance, and weighed bending
  energy too lightly for the interior to follow anything but the boundary's
  fitting error: a dome over four strips tangent to a sphere cratered to
  7.98 where its crown stands near 10. The energy's weight now follows the
  tolerance, and a net whose boundary residual outweighs the energy is
  refined before it is accepted: the same dome crowns at 10.09.

- **Offsetting fitted sheets works (#119).** Where refining a free-form
  offset cannot reach the approximation tolerance, `offset_sheet` and
  `make_thick_sheet` keep the closest fit if it lies within the face's own
  tolerance or a ten-thousandth of the distance, whichever is larger, and
  report the deviation in the moved face's tolerance; they failed at 1e-6.
  A filled saddle offset by 1 lies within 3.4e-6 of that distance.

- **Ruled and lofted sheets between sections of different edge counts
  (#118).** `make_ruled` and `make_loft_surface` match such sections by arc
  length: each section is cut where the others break, the pieces keep their
  exact geometry, and the sheet passes through every section exactly.

- **Exact mass properties take pcurves that run beside their edge.** A face
  whose pcurve runs beside its edge's curve by up to a thousandth of a
  millimetre is integrated with the strip between the two, and a junction
  whose ends lie anywhere within its vertex's tolerance closes; such faces
  went to a tessellation. Drills whose section edge took its curve from a
  nearby circle are measured exactly: cut and common add up to the part
  within 6e-9 over every stress drill on three seeds.

- **Crossings at a section's sharp turn are found.** A marched section that
  wobbles through a sharp turn within its stated error is fitted in pieces
  where its error exceeds the finest step the march took, and the
  curve-curve crossing search continues by alternate projection where its
  Newton solve stalls at a curve that almost stops along its parameter. A
  drill all but touching inside a torus's outer equator, its seam on the
  touching line, cuts at every placement tried instead of refusing.

- **Drills touching a torus or a ball cut it.** Where a drill touches a
  torus inside its outer equator, or a ball at its pole, the section is a
  figure eight, and the marcher traced it as one loop round both lobes or
  stalled short of the touch: depending on the order of the operands the
  cut was refused or came back as the whole torus, which `check` called
  valid. The figure eight is cut into arcs ending exactly on the touch.
  Cut, common, union and cells whose sections leave every piece of both
  solids outside the other while the solids overlap now refuse instead of
  returning a solid. Of 77 drill placements swept round a torus's equator,
  74 cut within 1e-8 (48 were refused or wrong).

- **Exact area and volume close junctions a pcurve misses.** A face's
  chart integral treated each loop as closed, so where a fitted pcurve ended
  a little off its vertex the result depended on where the integral
  started, by up to a strip round the whole face: a frustum drilled along
  its base measured its cut's cone wall 6e-4 square millimetres short, and
  a drill inside a torus's tube measured its common 1.8e-3 short. Such
  junctions are closed in the chart; cut and common add up to the part
  within 1e-10 on drilled cylinders, cones, spheres and tori.

- **A marched blend's rail follows its rail on each face.** A closed blend
  marched between two faces fits its rail's chart image on each face
  through the rail's own points, ending on the rail's vertex, not through
  the march's stations: on a drum leaning on a block it stays within 2e-5
  of the rail instead of 1.3e-4, and the fills at 1 and 2 degrees
  extrapolate to the upright torus ring within 8e-6.

- **Drills whose wall grazes an edge of the part.** A section touching an
  edge very flatly (a drill lying on a frustum's base through its rim, or
  along a round's tangent line) stopped short of the touch point and the
  drill was never split there; the stretch of edge within a micron of the
  drill is carried onto it as a contact. A section ending at a junction on
  the drill's seam goes to the copy of the junction it arrives at, not
  always the far one, and a touching section's junction is held to a
  hundredth of a millimetre instead of swallowing the edge. Eight stress
  drills over three seeds that were refused now cut, each checked against
  an independent integral within 1e-6.

- **Drills lying in a face that runs into a round.** A closed section's
  piece starting a hair short of its curve's seam was read back as one
  vertex and its arc dropped ("wire 0 is open"); a line through a figure
  eight's double point was crossed once instead of twice; and two
  crossings of an edge that the other solid grazes between them, at most a
  tenth of a millimetre apart, are one junction. Two NIST drills the stress
  run refused now cut, cut and common adding up to the part within 1e-6.

- **The sew joins edges by each end's own tolerances.** Two edges were
  taken for one when both their ends lay within the widest tolerance either
  edge held, so a vertex widened where a section was welded could join the
  two sides of a sliver a quarter of a micron wide at their other ends and
  move a corner: a pad fused into a slab drafted 0.00005 came out valid but
  2.1e-3 cubic millimetres off. Each pair of ends is now compared within
  those two ends' own tolerances. A facet leaning off a wall by tens of
  microradians crosses it on their shared edge wherever its far side stands
  off the wall by more than the weld distance: slabs drafted 0.00005 and
  0.00002 under a pad fuse in nearly every case tried.

- **Pads into slabs whose rows lean off the wall by thousandths.** Where a
  converted slab's facet rows lean off a pad's wall by a couple of
  thousandths of a radian, the sections either side of a row edge reach it
  a micron apart; where their gap is within the weld distance they join as
  one junction, its vertex claiming at most two microns. The mesher keeps
  an inner ring about a micron wide as a hole rather than dropping it as a
  slit. Every slab drafted 0.001 to 0.0001 tried fuses, cuts and is taken
  in common.

- **Converted faces no longer fold in space.** A curved face whose own
  edges cross within their tolerance (a sliver of a thin fillet fitted as a
  huge sphere) is faceted instead of built, and an edge between two curved
  faces meeting almost tangentially is threaded through points resting on
  both, so it keeps to them (0.06 mm instead of 1.1 mm off on a NIST part).
  Every NIST part meshed at a thousandth of its diagonal tessellates closed
  without the converter's closure check.

- **The converter's scatter estimate skips a bend's nearly flat edges.**
  The default coplanar distance no longer reads the nearly flat edges of a
  coarse curved or drafted wall as scatter, so such walls convert at the
  distance single precision resolves, each row a face of its own, with no
  volume doubt from vertices standing off merged planes (slabs drafted
  0.0001 to 0.0007: distance 1.7e-5 instead of up to 7.6e-4). Noisy meshes
  beside creases estimate their scatter more closely: a filleted disc
  exported with noise comes back as 5 recognized faces instead of 651
  facets.

- **Arcs on a turned frame are bounded whole.** An edge on a circle or
  ellipse arc whose frame is turned against the world axes was bounded by
  its ends and its frame's quarter points, short of the arc. Classification
  could skip the face such an edge bounds: the common of a half ball seamed
  along its flat face with a drum touching that face was refused as
  overlapping instead of coming back empty. An arc's box now holds its
  extremes along each world axis.

- **A shared arc across the start of a full circle.** A half ball charted
  about a tilted axis, standing on a drum of its radius (the two sharing
  the flat face and its rim), failed fuse, cut and common with
  "arrangement left no piece of the face". Where a shared arc runs across
  the start of a full-circle edge, every stretch the edge covers is read as
  boundary, on both sides of that start, not only the longer one.

- **An upright drum's face blend where its foot runs off the block.**
  `blend_faces` between an upright drum's wall and a block's top, where the
  foot runs off the block's side, builds an exact torus band ending where
  the foot's edge leaves the block, capped in the meridian plane there,
  with the fill the ring's share of the arc to 1e-11; it was refused. A
  closed contact circle starting on a face's seam is no longer cut into a
  sliver and refused when the face is split along it.

- **Drills through a dome's pole.** A trimmed sphere face (a dome, a cap, a
  dome on a drum) drilled by a wall passing through its pole was refused,
  turned inside out or given wrong volumes (a cut more than seven times the
  whole). Such a face is recharted about a clear axis before the boolean,
  as a whole ball already was, its trim edges carried over and its seam
  laid where no trim reaches: cut and common add up to the whole within
  1e-8 and the common matches an independent integral within 1e-9.

- **Face blends end where the crease's edges end.** `blend_faces` between
  faces that share no edge ends a capped round where the solid's edges
  along the crease end, if that comes before the ball leaves a face,
  capped in the ball's section through their end. A drum's foot leaning
  along a block's edge blends at 1, 2 and 15 degrees with the axis on the
  edge or past it, where it was refused; fills tend to the upright closed
  form within 3e-6. A cap whose plane runs almost along a drum's axis meets
  the wall in a chord instead of an ellipse too long to mesh.

- **Pads fuse into slabs drafted under a thousandth.** A planar face
  leaning off another by under a thousandth of a radian, the two sharing
  an edge to their points' rounding, crosses it on that edge rather than
  on the solved line a sliver off it, where the band in which they stand
  within the weld distance is under a tenth of either face. A pad pushed
  into a converted slab drafted 0.0001 to 0.0007 fuses, cuts and is taken
  in common instead of being refused.

- **A round's cap flush with a side face of the solid.** A marched round
  whose cap stands in the plane of a side face (an edge round or face blend
  of a leaning drum's foot running off the block's edge) melts into the
  solid: the band's end row lies in its cap plane and carries its measured
  tolerance. Arcs of one crease split at a seam meet in one cap plane
  instead of overlapping or leaving a sliver, which also raises fills that
  came out short at the seam. Fills tend to the upright closed form within
  3e-6.

## [0.8.0] - 2026-10-03

A minor release, since the API changes: `MeshSolidReport` gains a field,
and the mesh converter's steps are a public API (`MeshRegions`, with
`merge`, `split`, `fit` and `build`). The converter keeps an open mesh's
free boundary on its surfaces, seats a holed torus's seam clear of its
holes, tells a finely drawn curve from scatter and runs 3 to 5 times
faster on large meshes. Three silent wrong answers are fixed: a ball
drilled through its poles, empty boolean results taken as valid, and an
edge round over a drum's foot run off the block. Marched sections state
their measured error, pads fuse into slabs drafted less than their chords
sag, and `intersect_surfaces` answers `Same` for coincident patches.

### Added

- **Mesh conversion's steps as an API.** `MeshRegions::find` returns a
  mesh's regions (triangles, surface, deviation, neighbours); `merge`,
  `split` along a vertex path and `fit` of a `SurfaceKind` with
  `FitConstraints` (a fixed axis or radius) correct them, each verified at
  every vertex, and `build` makes the solid through the same pipeline as
  `solid_from_mesh`, which is now `find` then `build` with unchanged
  results. A step that cannot be taken returns a named `RegionRefusal` and
  leaves the regions as they were.

- **Tori drilled across their equators convert whole.** `solid_from_mesh`
  places a torus's seam round the axis on a parallel no hole crosses, as
  its seam round the tube already was, so a hole across the outer equator
  is an inner wire of one toroidal face instead of leaving the tube in
  facets.

- **Open meshes keep their curved faces' free boundary.** `solid_from_mesh`
  cuts a run of boundary with one face and nothing across it where it
  turns by the crease angle, and places it on the face's surface: an exact
  parallel circle or ruling where one holds, otherwise a curve fitted
  through its vertices' feet, with a pcurve and a tolerance measured
  against the mesh. Open tubes, half cylinder sheets, domes and free-form
  sheets come back as one face each instead of facets.
  `MeshSolidReport::free_edges_fitted` counts the fitted edges.

- **Face blends end where the seat runs off a face.** `blend_faces`
  rounds two curved faces of a solid, marched between them, where the ball
  touches both over only part of its seat: a crease that runs off the
  solid, or a face holding part of the loop. The round ends where the
  ball's line of contact leaves a face, closed by a flat cap in the ball's
  section, and is cut off or fused on as the closed round is. Slanted-drum
  fills tend to the torus ring's share within 2e-5.

- **Mesh conversion fits B-spline patches.** `solid_from_mesh` rebuilds a
  smooth region no plane, canonical surface or sweep fits as one B-spline
  patch, where the region is one disk and the patch verifies: charted by a
  nearly fitting canonical surface or by the mean-value map onto a square,
  fitted with thin-plate fairing and parameter correction, and checked at
  every vertex and inside every triangle, Jacobian and normals included.
  `MeshSolidOptions::patches` (on by default) controls it, and the report
  counts patch faces and the regions refused (`patches_not_disk`,
  `patches_narrow`, `patches_unverified`). `fit::fit_surface_scattered_at`
  fits a scattered surface at the caller's own parameters and knots, with
  thin-plate fairing.

- **Rounds between curved faces with no shared edge or axis.**
  `blend_faces` and `fillet_faces` round curved faces that share no edge
  and whose surfaces share no direction or axis (a B-spline face, a
  cylinder at a slant to a plane). The ball is marched along where the
  surfaces cross, or, between separate faces whose surfaces do not cross,
  along where their offsets by the radius cross; the round is a B-spline
  band fitted through the ball's arcs, tangent within a tenth of a degree.
  On a solid the corner is cut off or filled where the ball touches both
  faces all the way round a closed seat; between separate faces the round
  spans where both faces reach, and with `trim` the faces and the round
  sew into one shell. A seat a solid's faces hold over part of its length
  only, a marched round closing on itself between separate faces, and a
  line of contact crossing its face more than once are refused by name.

- **A coincident patch's edge that weaves across the other's window is
  carried stretch by stretch.** An edge of a B-spline patch that leaves a
  coincident patch's window and comes back onto it was refused. It is
  carried into that patch's chart as one contact per stretch on the
  window, each with its own pcurve; a sheared bar across a sheared drum on
  a shared plane fuses, cuts and intersects to the exact volumes.

- **`blend_faces` blends curved faces.** Faces that meet along edges of
  the solid are rounded along those edges: exact cylinders and tori where
  the seat has a closed form, marched B-spline bands otherwise, B-spline
  faces included. Faces that share no edge are rounded with an exact
  cylinder or torus where their surfaces share a direction or an axis,
  the corner cut off or filled as the solid shows. Faces sharing no edge
  and no direction or axis are refused by name. Closed marched fillet
  bands sample each ball arc more finely, and stay tangent within a tenth
  of a degree.

- **Tangential contact traced as its own curve.** Two surfaces touching
  along a curve no closed form names (a torus resting in a free-form
  cradle) came back as fragments. The trace predicts along the null
  direction of the difference of the two surfaces' second fundamental
  forms and corrects onto the locus where they meet with parallel
  normals, landing on a patch's edges and closing across a seam: one
  curve with pcurves on both surfaces, marked tangential. A seed at an
  isolated touch, or where the surfaces come near without touching, is
  refused by name.

### Changed

- **`intersect_surfaces` answers `Same` for coincident patches.** Two
  patches with no closed form that lie on one another within the requested
  tolerance are reported as coincident for every caller, where fillets,
  offsets and sections used to march noise along them. Surface projection
  (`project_on_surface`, `project_on_surface_from`, `SurfaceSeeds`,
  `SurfaceProjection`) lives in `ogeom-geom`; the `ogeom-algo` paths still
  work.

- **Mesh conversion is 3 to 5 times faster on large meshes.** Each seam
  curve is solved once per conversion and the per-face checks (areas,
  probes, the exact volume) run in parallel, with identical results: a
  23k-triangle part converts in 9 s instead of 48 s.
  `ogeom_mesh::triangulate_with_chords` returns a shape's triangulation
  together with the edge chords its faces agreed on.

### Fixed

- **Marched sections state the error their curves have.** A marched
  intersection's tolerance is measured in space (the curve against its
  trace, each pcurve on its surface against the curve, and the surfaces'
  gap over the sine of their angle where they cross shallowly) instead of
  the fit's residual times the surface's stretch, which on a wide drum
  stated hundreds of times the real error and widened the boolean's
  junctions as far. A branch whose fit misses by a hundred times its
  tolerance is fitted in pieces, and a section reaches an edge within its
  own reach plus the edge's tolerance. Two more drills in the stress run
  succeed.

- **Booleans with a ball whose poles lie on the other solid's wall.** A
  drill whose wall passes through both poles of a ball cut nothing away or
  everything (an empty or whole result `check` called valid), or took the
  common a little short. A whole ball meeting a curved wall at a pole is
  charted about another axis of its frame, clear of the wall, before the
  boolean; cut and common add up to the ball and match an independent
  integral within 1e-6 over drills in six directions and three radii.
  Empty results are checked too: a cut leaves nothing only where no part of
  the solid lies outside the tool, a common only where the solids share no
  volume, and a fuse never; otherwise each is refused by name.

- **Mesh conversion tells a finely drawn curve from scatter.** The
  estimate behind the default coplanar distance weighs each nearly flat
  edge's turn against its neighbours': a steady turn from one edge to the
  next is a curve, not noise. Slabs drafted by fractions of a degree keep
  each row a face of its own at the default distance (pads into them fuse,
  cut and are taken in common at every row count tried), and a free-form
  sheet with no flat stretch comes back one fitted patch instead of
  thousands of facets.

- **Exact area and volume on a face holed by a closed spline.** A face's
  loops are told apart by their areas, integrated span by span; one rule
  over a whole closed spline could read a hole's direction wrong and add
  the hole to the face.

- **Edge rounds stop where a drum's foot runs off the block.**
  `fillet_edges` on an open elliptic crease (a leaning drum whose foot runs
  off the block) rounded the whole ellipse and left material past the
  block's side; the round now covers the arc only, or is refused by name
  where the boolean cannot merge its end with the face the crease runs out
  on. Fillets and chamfers on a circular foot that runs off the block took
  the whole turn, read the wrong side of the wall and removed material;
  they now round only the arc, with exact torus or cone fills.

- **Pads fuse into slabs whose rounded corners lean in by less than
  their chords sag.** Where a pad wall crosses an edge between two facets
  lying almost in the wall, the two facets' sections reach the edge up to
  a few ten-thousandths of a millimetre apart; they now join as one point,
  and a stub of section between two such crossings no longer leaves a face
  piece with no inside. A slab drafted 0.001 fuses, cuts and is taken in
  common at every row count tried, while its rows stay faces of their own.

- **A boolean cuts periodic faces whose outlines sit turns away from
  their chart, and fits sharp-tipped section loops in pieces.** A face on
  a periodic surface whose outline was stated one or two turns outside the
  chart's range read every section piece as outside it, so the section was
  dropped and the cut did not close; its pcurves are moved home by whole
  periods when it is gathered. A closed section loop that turns sharply
  where two drums all but touch, fitted whole, missed its trace by
  millimetres and swelled the junctions it fed; it is fitted in halves.
  Two oblique drills through a NIST part that were refused give valid
  solids whose cut and common add up to the part.

- **A sphere face bounded by a great circle through its poles is
  healed into half meridians.** Its stored pcurve was a fit that strayed
  far from the meridian, so a boolean's section ended away from where the
  boundary was paved and a drill across the hemisphere was refused.
  Healing splits such a circle at the poles into exact half meridians,
  joins pieces meeting on a pole's row by a pole edge, and rebuilds the
  neighbouring face on the pieces; the repair is refused where the rebuilt
  face's area would change. The part's volume then integrates exactly.

## [0.7.0] - 2026-10-02

A minor release, since the API changes: `Canonical` is no longer `Copy`,
`MeshSolid` and `BlendContact` gain fields. Its larger part is surface
modelling on sheets: sewing across a gap and a prism's lids, trimming and
splitting sheets with booleans, splitting a face along curves, growing a
face past an edge, offsetting and thickening sheets, lofts, ruled surfaces
and sweeps as sheets with guide curves, blend surfaces and curves, N-sided
fills on their own edges, fillets on sheets and between separate faces,
curvature continuity in blend analysis, and sheet export to STEP and IGES
(#104 to #115). The mesh converter rebuilds surfaces of revolution all the
way round, fillets and corner balls from the faces they blend, and every
converted part it is measured on now tessellates closed; a boolean keeps a
piece reaching across a periodic face's seam as one face and classifies a
piece touching the other solid only along lines.

### Added

- **Fillets on sheets and between separate faces.** `fillet_sheet_edges`
  rounds edges where two faces of a shell meet, and `fillet_faces` rounds
  the corner between two faces of separate shapes, the ball on the side
  their normals point to, optionally trimming both back into one shell
  with the round. The round is an exact cylinder where the supports share
  a direction, an exact torus where they share an axis (a plane square to
  it, cylinders, cones, spheres, tori), and on a sheet a marched B-spline
  band elsewhere. Each face is rebuilt along its line of contact and
  shares that edge with the round. A B-spline pair given to
  `fillet_faces`, and a closed sheet edge with no closed-form round, are
  refused by name.

- **Sheets offset and thicken.** `offset_sheet` moves a face or a shell
  along its normals, and `make_thick_sheet` thickens a face or an open
  shell into a solid, on one side or both, closing it with side faces
  along its free edges. Analytic faces move to their exact parallel;
  B-spline and other faces are fitted within the approximation tolerance,
  the measured deviation held in the face's tolerance. Trims carry over
  unchanged. Creases between faces, offsets past a radius of curvature,
  edges shared by more than two faces and layers that run into each other
  are refused by name. `fit_surface_grid_at` fits a tensor-product grid at
  given parameters, keeping the caller's chart.

- **N-sided fills on their own edges.** `make_filling_n` fills a hole
  bounded by any number of edges with one face bounded by those same
  edges, with pcurves, so `sew` joins it to its neighbours. Each side
  meets its support face G0, G1 or G2, and the surface passes through
  given interior vertices and edges. The gap, angle and curvature step
  reached are measured per side and returned in `Filled`. A fit that
  cannot meet the tolerance, and a hole that is not a height field over
  its boundary's plane, are refused by name.

- **Blend surfaces and blend curves bridge two edges.**
  `make_blend_surface` bridges two edges of different faces with a
  B-spline face that shares both edges and meets each face C0, G1 or G2:
  exactly against planes, and within a measured budget (1e-5 rad, a
  curvature step of 1e-4 over the gap's width) against curved faces.
  `make_blend_curve` bridges the ends of two edges G1 or G2, `End` naming
  an edge's start or end as traversed. Circle and ellipse edges, closed
  and seam edges, and C1 or C2 asked across unrelated surfaces are refused
  by name.

- **`split_face` cuts a face along curves on it or projected onto it.**
  Curves lie on the face, or are dropped onto it along its normals or
  along a direction (`Projection`); curves meeting end to end are
  followed as one, and a curve crossing a periodic face's seam is followed
  across it. Each stretch from boundary to boundary divides the face, and
  a closed curve inside it cuts out the region it encloses. The pieces
  share the new edges, neighbouring faces keep sharing the boundary edges
  a cut ends on, and history maps the face to its pieces. A curve that
  ends inside the face, misses it, or runs along its boundary is refused
  with a reason.

- **Booleans trim and split sheets.** `common` and `cut` accept a face, a
  shell or a compound of those against a solid or a half space, and return
  what of the sheet lies inside or outside it, sewn into a shell or a
  compound of shells. `split_sheet` splits a sheet where another sheet
  crosses it and keeps every piece. Pieces keep their surfaces, and
  history maps each to the face it came from. A sheet as the tool of
  `cut`, two sheets given to `common`, a sheet given to `fuse`, and a
  piece lying on the tool's boundary are refused by name.

- **Sheets from lofts, ruled surfaces and sweeps.** `make_ruled`,
  `make_loft_surface`, `make_sweep_surface` and `make_sweep_two_rails`
  build sheets from edges or wires, open or closed, planar or not. Lofts
  and ruled surfaces pass through their sections exactly (exact rational
  B-splines, degree one across when ruled, a plane between coplanar
  segments, optionally closed back to the first section). Sweeps follow a
  frame law or two rails, held to ten confusions of the swept profile and
  measured. A sheet has a face per section edge, and bounding section
  edges are the caller's own, so the sheet sews to them. A loft follows
  guide curves as a Gordon surface through single-edge sections and their
  guides: exact on the sections and within ten confusions on each guide,
  measured; a guide missing a section is refused by name.

- **Faces grow past a boundary edge with `extend_face`.** A face is
  rebuilt with the edge moved `length` out along the surface, square to
  it, and the face's other edges stay shared with its neighbours. Planes
  grow across any straight edge; cylinders, cones, spheres and tori grow
  along a parameter line on their own surface; B-spline faces grow past
  the side of their patch, naturally or linearly at G1 or G2
  (`Extension`). History maps the old face to the new one, and the cases
  it does not handle are refused by name.

- **Blend analysis measures curvature continuity, and faces give
  curvature samples.** `BlendContact` carries `curvature_error`, the
  largest step in normal curvature square to the shared edge over the
  stations; an edge that cannot be measured reports it infinite.
  `face_curvature_samples` returns the principal curvatures on a grid
  inside a face's trim, placed where the face stands and signed against
  its own normal, for curvature maps and zebra stripes.
  `SurfaceCurvature::normal_curvature` gives the curvature along any
  tangent.

- **STEP and IGES export write sheet bodies.** `write_step` writes shells
  no solid owns and faces no shell owns as `SHELL_BASED_SURFACE_MODEL`s in
  a `MANIFOLD_SURFACE_SHAPE_REPRESENTATION`, and a part holding solids and
  sheets as both representations in one part. `write_iges` writes sheet
  faces as trimmed surfaces (144) beside the manifold solid B-reps it
  writes for solids, keeping a reversed face's side. Both writers refuse
  wireframe (a free wire, edge or vertex) and empty parts by name instead
  of dropping them.

- **Sewing across a gap.** `sew_within` sews faces whose edges meet up
  to a given distance apart, as surfaces built separately (an imported
  sheet, a fitted fill beside an exact extrusion) do: two edges whose ends
  and middle lie within the gap are one, and the edge kept and its
  vertices widen their tolerances to reach the edge it replaces. An edge
  that runs along two or more shorter edges of other faces is first split
  where their vertices meet it, so each piece has a twin. Two unit squares
  a twentieth of a millimetre apart sew across a tenth into one shell and
  stay two across a hundredth; a 10 mm edge against two 5 mm edges sews
  into one shell.

- **A boolean takes a compound of disjoint solids as an argument.** A
  grid of drums cut from a plate at once is one boolean, not one per drum:
  four hundred bores take 0.2 s where four hundred cuts took over a
  minute.

- **Refine around the faces an operation made.**
  `unify_same_domain_around` merges the given faces with the faces on
  their carriers they reach across shared edges, and leaves every other
  face of the solid as it is, split or not: an application refining after
  a feature unifies what the feature made without undoing splits the part
  already held.

- **Mesh conversion rebuilds extrusions and surfaces of revolution.**
  With `MeshSolidOptions::sweeps`, a smooth region none of the canonical
  surfaces fits is tried as an extrusion of a fitted profile (its normals
  all square to one direction) and as a surface of revolution (its
  normals all meeting one axis), each held to the coplanar distance at
  every vertex. The axis read from the normals is brought close by a
  search on the samples' radii against their heights, then settled by
  Gauss-Newton on the samples themselves: each round fits the profile and
  solves for the axis's offsets and tilts with a correction to the
  profile, which separates a tilt from the profile's shape over a whole
  turn as over part of one. A recognized band the region runs into
  smoothly is offered with it, so a stretch of a wavy profile that passes
  for a cone joins the one sweep. A wall extruded from a wavy spline comes
  back one extrusion and five planes, not nineteen facets; a wavy profile
  turned half round, one surface of revolution, and turned all the way
  round, one surface of revolution closed on itself between two flat
  ends. On by default.

- **Mesh conversion builds fillets and corner balls from their
  neighbours.** A torus between a plane and a cylinder or cone on an axis
  square to it is put on that axis, its tube's centre a radius off both;
  rounds between planes that meet at corner balls take one radius, and
  each ball is centred a radius off the three planes its rounds run
  between, where their axes meet. Only the radius is fitted, and each
  surface is kept only where it holds every vertex of its region. From a
  single-precision mesh of a rounded block the corner balls now meet their
  fillets tangentially to rounding, where the fitted ones stood off their
  axes by a few tenths of a micron.

### Fixed

- **A boolean piece touching the other solid along lines or points is
  classified, not refused.** A face piece whose every probe landed on the
  other solid's boundary (a cylinder inscribed in a box touches its walls
  at the quarter turns the probes use) was refused for want of a
  coincident partner. It is asked again at points off the contact, which
  a contact along lines or points cannot cover. A hole touching its outer
  boundary at isolated points (a disc inscribed in a square) is placed by
  a vote of chord midpoints, not by a point that may be the one it
  touches at.

- **A boolean keeps a piece reaching across a periodic face's seam as
  one face.** The boolean splits each face in its own chart, whose edge
  is the seam, so a piece reaching across the seam came back as two faces
  sewn along it: a tube or drum halved square to its seam gave three
  faces. Pieces of a cylinder, cone, sphere, torus or full revolution
  that meet across its seam are joined into one face, on the surface
  turned so its seam falls where the face does not reach, in solid
  booleans, sheet booleans and `split_sheet` alike. `History` gains
  `without_repeated_images`.

- **A prism's far end edges sew.** `sew` told edges apart by their nodes,
  and a prism's far end edges are its profile's own edge nodes carried
  along it by a location: each was read where the profile stands, so a
  lid at the far end joined nothing. Faces that hold an edge node placed
  twice are baked first, every edge its own node where it stands, and an
  edge's curve is read through its own location as well as the edge's.

- **A converted solid that would mesh open gives up the curved faces that
  open it.** A face whose trim folds within the tolerance its own seams
  claim passes every check on the solid and still meshes over itself.
  Where the solid's tessellation does not close, each face is meshed on
  its own, and the curved faces using a mesh edge other than twice (or,
  where only planes do, with an edge drawn within the reach of it) fall
  back to facets. Every converted part of the corpus and the NIST set now
  tessellates closed: ftc_07 gives up 20 curved faces of 143, Body28 3,
  ctc_02 2, handle-pickup 1. A part that needs it converts slower (ctc_02
  from 27 s to 50 s); one that meshes closed pays one tessellation.

- **A converted torus is no longer built with its seam through a hole.**
  A hole's outline was unwrapped round the torus's axis but not round its
  tube, so a hole across the outer equator, where the tube's angle starts,
  read as clear of the seam there. The whole torus was built with the hole
  as an inner wire and its seam edge running through the hole: a face
  `check` passed and no tessellation could close. Such a ring's tube now
  stays facets; its bore still comes back a cylinder.

- **Faces meet at a vertex their edges end off.** An edge's curve may end
  off its vertex by up to the vertex's tolerance, and the next edge's off
  it on another side. Each face round the vertex closed its own gap
  between two of the ends, and between them the faces left a hole; on a
  curved face the boundary also ran from the curve's end to the vertex and
  back in its chart, a fold its triangles covered twice. Each edge is now
  drawn from its vertices, in space and in every face's chart, so the
  faces round a vertex meet at one point. Converted parts whose faces,
  meshed one by one, left edges unmatched: 77777_1 from 64 to none, Body28
  from 122 to 31, ctc_02 from 25 to 6, ctc_01 and ctc_03 to none. Every
  imported STEP solid in the corpus closes as before.

- **Facets left between curved faces join them.** Where curved regions
  meet, or a fit stops a row short, a mesh leaves planar facets of a
  triangle or two that no surface claims, bounded by seams as wide as
  themselves, whose trims fold. A facet whose corners all lie within the
  seams' reach of a curved neighbour joins it; a region that cannot be
  built with the facets it took gives them back before it is faceted, and
  a conversion that leaves a face turned into the material is made again
  without them. A coarse plate with filleted rims comes back its own 18
  faces instead of 56, the truth bench's faces over the exact count fall
  from a mean log ratio of 0.27 to 0.17, and large real exports lose up
  to a third of their faces with every curved surface kept.

- **A needle on a flat face is part of it.** A triangle no higher than
  the coplanar distance (three corners all but on a line) has no plane of
  its own at that distance, but its rounded normal could lean far enough
  that it became a face of almost no area, turned whichever way its
  rounding pointed, which `check` could then find facing into the
  material. Such a triangle joins a neighbouring plane by its corners'
  distance alone.

- **Mesh conversion raises a coplanar distance the mesh cannot meet.**
  With recognition, a distance asked below two and a half times the
  mesh's measured scatter (or twice its quantum) left regions in
  fragments whose seams opened: a single-precision rounded block asked to
  a tenth of a micron tessellated open while `check` passed. It is raised
  to that floor, `MeshSolidReport::coplanar_distance_raised` says so, and
  the new `MeshSolid::coplanar_distance` is the distance the faces were
  built to.

- **A converted face whose boundary folds over itself is withdrawn.** Each
  built face's boundary is drawn in its surface's chart, and two of its
  chords that cross are a fold. A fold deeper than the tolerance of the
  two edges that cross, outside every corner's tolerance, on edges of the
  face other than a seam, and still there when the boundary is drawn a
  hundred times finer, names its two seams, and the curved faces across
  them fall back to facets. `check` passes such a face; its triangles
  cover part of it twice.
- **A seam threaded along a tangency no longer hooks past its corner.**
  Where a fillet meets a face tangentially the mesh's boundary between
  them wanders, and can step past its corner and back; the curve threaded
  through it ran on past the corner and came back, folding the trims of
  both faces. Only the points that close in on the seam's end are
  threaded, and where a cubic through them still overshoots, the polyline
  through them is.

- **A point on a face of any repeating surface classifies on the right
  side.** Points were folded by a period toward a face's outline only on
  cylinders, cones, spheres and tori; on a whole turn of a revolution or a
  periodic patch, an outline kept at negative angles left the point a
  period away and the solid classified inside out. The surface now says
  whether and how it repeats.

- **Crossings along an edge each known to within its doubt meet where
  both doubts reach.** A pad pushed into a converted slab whose corner's
  rows hold different numbers of facets refused: the facets' diagonal
  edges run all but in the pad's walls, and the walls' sections on the
  facets either side of one reached it tens of microns apart, each within
  its own doubt but not within the larger one alone. Two crossings are now
  one point where they stand no further apart than their doubts together.

- **A coarse mesh comes back the same whatever order its triangles are
  listed in.** A fillet two or three facets wide, first sampled together
  with the corners beside it, fitted nothing and was left faceted, or not,
  by the luck of the order. Its triangles are now sampled again once the
  corners are claimed, and a first sample no longer reaches across a flat
  face into the fillets on its far side. A sphere or a torus fitted to a
  handful of noisy vertices is no longer taken. Several converted parts
  come back with far fewer faces; one of 6268 triangles, 1615 faces before,
  comes back in 13.

- **A converted part tessellates closed where a single facet meets a
  fillet.** A facet left beside a fillet, thinner than the curve its plane
  meets the fillet along bulges, had that curve as its seam: the facet's
  trim folded back across itself, its triangles covered the fold twice,
  and the whole part's tessellation opened there while every check passed.
  Such a seam is now threaded straight through the facet's corners, its
  tolerance how far it stands off the fillet, and the fillet keeps its
  surface; where the straight seam strays too far, the fillet is built
  faceted. A sliver between two curved faces whose seams cross in its
  plane, and any facet whose mesh does not meet its neighbours', is
  threaded the same way. A threaded build that leaves a face turned into
  the material or collapsed onto a line is set aside for the unthreaded
  one. Three converted corpus parts and the rough rounded box now
  tessellate closed.

- **A polynomial's roots are found whatever their pattern.** The real
  roots of a quartic or higher came from the companion matrix's
  eigenvalues, whose iteration had no limit; a quartic whose roots pair
  off by sign (a line through a torus's middle) defeated its shifts and
  the call never returned. The iteration is now capped, and where it does
  not settle the roots are found by Durand-Kerner instead.

- **Halves of a cone sharing their seam ruling mesh closed.** The ruling
  two halves on one surface share has one image in the chart, at the seam,
  and the half running from a full turn back to half a turn needs it a
  turn over. Walked into the apex the two columns are equally near, the
  ring wound onto the other half, and both halves meshed over one side: a
  NIST part read from AP242 tessellated open. The winding is now undone
  from the tie where the walk left it.

- **A converted solid no longer comes back with a facet facing into its
  material.** On a coarse mesh a fitted face could run past the facet it
  should have ended on and enclose a sliver outside it; the solid checked
  sound but for that facet's orientation. Once a build would be accepted,
  the faces round every recognized face are probed as the validity check
  probes orientation, and a recognized face that overlaps another is built
  faceted. The STEP corpus, meshed and converted back, is now a test.

- **A face moved or offset clear of a neighbour it was tangent to is
  refused.** Moving the side of a block rounded all round outward by half
  a millimetre left the face's surface there and its edges where they had
  been, on fillets it no longer met, and the solid came back valid with
  less volume than before. Every corner the edit re-solves is now held to
  lie on each surface it sits on, and the edit is refused by name where it
  does not.

### Changed

- **`recognize::Canonical` is no longer `Copy`, and has a `Swept`
  variant** holding a surface of revolution or an extrusion
  ([`SweptShape`]). Code matching on it needs the new arm; code copying
  it, a `clone`.

- **`tight_bounds` is hundreds of times faster on ruled and round
  parts.** Every face was meshed and searched for an interior extreme, and
  every edge searched numerically. A face on a plane, drum, cone or linear
  extrusion reaches no further than its boundary, a straight edge leads at
  its ends, and a circle at an angle its frame names: those are taken as
  they are. A plate with four hundred bores bounds in 1 ms, not 311 ms, to
  the same box.
- **A small cut into a large solid costs what it touches.** A face the
  tool leaves alone (no section on it, no edge of it split, no coincident
  partner) is neither split nor arranged: its rings are its own wires and
  its state is read once, at the middle of one of its edges. The other
  solid's trims are drawn only where a point or a ray comes near them, a
  crossed face's outline is drawn once, circles are sampled to their sag
  in closed form, outline ends are welded among their neighbours, a
  piece's scanlines come from a heap as the search asks, and an edge
  rebuilt whole is built once for both its faces. Sewing looks up eight
  grid cells per point, not sixty-four, compares a loose edge only with
  edges starting within its reach, and finds a middle's foot on a line or
  circle in closed form. A pocket into a plate with four hundred bores
  takes about 50 ms, from 160.
- **A line meets a cone or a torus in closed form.** Both took the
  general seeded path, a Newton polish per seed; a cone's crossings are
  now a quadratic's roots and a torus's a quartic's, tangencies included,
  each polished once onto the surface, the quartic solved about the
  line's nearest approach to the torus's centre in units of its size, so
  a torus far down the line is not lost to rounding. Every ray a
  classification casts
  through a part with tori is cheaper: converting a part with three tori
  takes 2.5 s, not 27.

- **A chain of edges rounds in a few booleans, not one per edge.** The
  chain is taken in rounds, no two edges of a round sharing a vertex
  (pieces of one circle count as one); each edge's blend is built against
  what the earlier rounds left, and a round's blends are applied together,
  one cut and one fuse. The checks before the blends mesh each face once
  and pass over edge pairs whose boxes stand apart. A hundred bore rims
  round in 1.4 s, from 10, and a box's twelve edges in half the time.

## [0.6.1] - 2026-10-01

A patch release: `cargo semver-checks` against 0.6.0 finds no breaking
change, and every crate packages and verifies in a publish dry run. Most of
it answers an adversarial review of the kernel: mirrors and negative scales
transform curves, surfaces, primitives and booleans the right way round;
spheres come back whole through STEP and IGES; unclamped and periodic
splines reverse, rescale and split unchanged; offsets and fillets refuse what
they cannot build instead of returning it inside out. The boolean settles
near-coincident faces by one rule, welding a sliver thinner than the weld
distance and keeping anything thicker as real geometry, and pads on a part
converted face for facet now fuse back into it. A pipe tee's junction loop
and a slab's full round fillet, and a thread groove cuts the same at any
start angle. A history now says which faces came through a boolean
unchanged (`History::copy_of`), and a boolean on a face with hundreds of
holes is several times faster.

Known limits, each with its diagnosis in `docs/PLAN.md` or
`docs/REVIEW.md`: a converted rounded block's whole-shape tessellation can
stay open along its corner spheres while `check` passes, and a pad into a
slab whose faceted walls lean in by less than their chords sag refuses at
some facet counts.

### Fixed

- **Pads on a part converted face for facet.** A pad pushed down from a
  top face of a converted part, its walls standing on facets a hundredth of
  a radian off them, refused ("the kept pieces did not close into a
  shell"). Two planes sharing an edge only to their points' rounding were
  crossed on a line the solve left a sliver under the edge; they now cross
  on the edge itself. And where a facet row's edge runs all but in a wall's
  plane, the sections on the facets either side reached it microns apart
  and split it twice; they are now one junction.

- **A slab rounded at half its thickness.** Rounding every edge of a box
  at half its thinnest side refused ("an earlier blend in the chain
  consumed this edge"): the balls of the two corners at each short edge's
  ends meet across it and leave nothing of it. An edge of the chain whose
  ends are both corners the chain rounded, and which they consumed, is now
  rounded by them; the result is the full round, its volume Steiner's.

- **A closed crease asked for as a chain of its arcs.** `fillet_edges` on
  every arc of a junction loop refused after the first ("an earlier blend
  in the chain consumed this edge"): the first arc's blend runs the whole
  loop. An arc of the chain lying between the same two faces as one
  already rounded is now done, and the history records it consumed.

- **A pipe tee's junction loop would not round.** Filleting the loop where
  a branch drum meets a main drum refused ("the kept pieces did not close
  into a shell"). The blend's leg on the branch is a band between the loop
  and the ball's contact ring, and the two rings' chart images started on
  one column a period apart, so the band's seam wound a full turn round the
  branch as a helix. The second ring's image is now moved into the first's
  period and the seam runs straight along the column.

- **Drums the weld distance apart in axis.** Two drums whose axes stood
  1e-5 to about 1.07e-5 apart refused every operation. The rim arc across
  the thin crescent one leaves on the other's cap was asked whether it lies
  inside the cap only at its middle, where it lies inside by less than the
  boundary's own doubt; it was dropped and the drum was not split. A piece
  is now asked at points either side of its middle as well.

- **A thread groove cut into a rod depended on its start angle.** A
  helical sweep cut into a body of revolution failed at some angles about
  the axis ("edge N ends where edge N+1 does not begin", or a shell that
  did not close) and passed at others. A marched section stopped a
  fraction of a step short of the edge of the patch it left, so it missed
  the next patch's section by more than either curve's tolerance. The
  marcher now ends such a section on the edge itself; a sphere's pole and
  other bounds that collapse to a point are not edges and are left alone.

- **Blocks tilted by a hair, and slivers at the weld distance.** Two
  blocks stacked and hinged, or standing in one another turned, by 1e-5 to
  1e-4 rad refused. Two planes all but parallel now cross on their true
  line (measured from a point on the first, not the world origin, whose
  offsets cancelled to a third of a micron); an edge of either lying on the
  other's plane splits it as a coincident pair's edge would; and a target
  edge splits where a contact ends on it. The welding retry takes a second,
  wider step, held to the same tolerance limit, and drums whose sliver
  lies at the weld distance itself are sectioned rather than lost between
  the two readings.

- **Slivers under the weld distance are welded.** Drums a few microns apart
  in radius, axis or end, and a box tilted a few microradians on another,
  refused. A boolean that refuses now runs once more at a confusion
  distance ten times wider, where such a sliver is one the boolean welds,
  and keeps that result only if nothing in it states a tolerance past the
  weld distance, a hundred confusion distances, so no thicker sliver is
  lost. Drums 3e-5 apart in axis, whose crossing lines were not sought
  closer than 1e-4, now meet in them.

- **Drums a hundredth of a micron to a tenth of a millimetre apart.** Two
  drums 2e-5 to 1e-4 apart in radius, or 1e-4 apart in axis, refused or
  (offset in axis) came out with a misbuilt face. Where a planar face has a
  curved edge, a section's place inside it and a piece's probe are now read
  against the exact boundary rather than its chords, a point within the
  weld distance of the boundary counting as on it; two circles in one plane
  cross where their closed form says, not along a sampled run millimetres
  long; a closed section no longer has a crossing pulled to its own
  arbitrary start; and an exact section ending at a junction is carried
  onto it rather than widening the face's weld to the junction's span.
  Volumes match the closed forms.

- **A cavity touching its wall is the solid's void.** Cutting a ball from
  a box it touches from inside left the box as one solid and the ball as an
  inside-out solid of its own beside it, where the ball reached a tenth of
  a micron past the wall: shells were nested only when one's box held the
  other's exactly. The container's box now takes the touch tolerance.

- **Faces within tolerance of the other solid's.** A drill whose wall
  pokes a tenth of a micron out of a box side, a ball resting a tenth of a
  micron into a box top, and a box stacked on another and tilted a
  microradian all refused: a piece lay on the other solid's boundary with
  no coincident partner face. Where that refusal stands after the nested
  fallback, the boolean now runs once more, settling each such piece by
  asking just off it on both sides: outside or inside on both, it is
  that; inside only on its own material side, the other's face runs along
  it backed the same way.

- **Small corner rounds and an L-bracket's every edge.** Rounding a
  block's three edges at a corner refused at any radius up to 0.1: the
  ball patch's rim leaves each sharp edge tangentially, and an arc stated
  as a trim of its circle was held to the loose overlap width meant for
  fitted curves, so the boolean read the rim as running along the edge for
  a stretch and split it where the band beside it was whole. A trimmed
  exact curve is now as exact as its basis. Every edge of an L-bracket
  failed in either order: the corner tool took the step's re-entrant
  corner, where a concave edge meets two convex ones, for a convex corner,
  because the planes through it taken whole look like one. It now probes
  beside each corner edge and leaves a re-entrant corner to the blends.

- **A STEP void reads as a cavity whichever way it is written (#102, #103).** A void
  shell whose flags turn it to add material (written inside out and then
  marked reversed as well, as some exporters do) read inside out: the
  cavity counted as extra material and `check` reported its faces broken.
  Each void is now oriented by its geometry as it is read, with a warning
  naming it, and keeps its faces and placement.

- **Drums a thousandth apart.** Two drums 1e-3 apart in radius or axis
  refused: the wider one's cap leaves an annulus or crescent that thin, and
  its outline's chords bowed by more than its width, so the point the
  piece was classified by fell outside it. A curved strand on a plane is
  now sampled until its chords stand within 1e-4 of it. Two drums 1e-4
  apart in axis, whose union came out short by 0.03, now refuse instead.

- **A cut by an identical solid, and a cap that all but meets a wall.**
  Cutting a solid by an identical copy described another way (two spheres
  on turned frames) failed in the arrangement; a solid lying within its
  tool is now cut away entirely. A drum a hair wider than the box it
  stands in refused at most turns of its seam: a planar face's bound was
  its rim's sampled points, and a disc's rim bows out between them by more
  than the drum's cap pokes past the box wall, so the pair was never
  intersected. A plane's bound now carries its rim's sampling slack, as a
  ruled face's does.

- **Three intersector edge cases answer.** A plane tangent to a cone
  along a ruling read `Apart`: a plane through a cone's apex now holds its
  rulings in closed form (one tangent, two crossing, or the apex alone),
  each stated over the cone's window. One circle written twice in the plane
  gave both whole domains as its overlap, whatever the phase and winding;
  `on_b` now names where the second stands at the first's ends, as in
  space. Collinear segments meeting end to end, in the plane or in space,
  returned nothing; they share their end point.

- **A periodic pattern traces its faces and may have gaps.**
  `make_periodic` recorded history only for the whole solid; every face
  and edge of the cell now traces to its own image and generates its
  copies' images through each fuse. A pattern whose copies stand apart
  failed at the second copy, the fuse refusing the compound the first one
  left; the copies are now fused only into the pieces they join.

- **Every outline is drawn once.** HLR `project` read its silhouettes off
  a mesh whose faces were meshed apart, so every face border counted as an
  outline and every outline edge was drawn a second time over the model
  edge. The mesh is welded first, and a silhouette along a drawn edge is
  left to the edge.

- **A bare or placed profile tapers.** `make_prism_tapered` refused a
  face straight from `make_face`, whose edges carry no pcurves yet, and a
  placed profile, whose material side it probed in the wrong frame. The
  bare profile gains its exact pcurves, and a placed one is swept where it
  was built and the prism placed where it stands.

- **A straight run of corners is no obstacle to a face.**
  `make_polyhedron` took a face's normal from its first three corners and
  refused a ring starting on a straight run (and would have turned a
  concave ring starting on a reflex corner inside out); the whole ring's
  winding decides now. `medial_axis` refused a polygon with a corner on a
  straight side; such a corner is dropped, and the axis is the polygon's
  without it.

- **Touching circles have their third tangent line.**
  `lines_tangent_to_two_circles` dropped a pair of tangents that had
  closed into one, so circles touching from outside got two lines instead
  of three and circles touching from inside none instead of one.

- **Polynomial roots keep their double roots.** A tangency's double root
  sits where the discriminant is zero only up to rounding, and the quadratic
  and cubic solvers dropped it when rounding left the discriminant a hair
  below or above; above the cubic, the companion matrix split it into a
  complex pair and it was dropped there too. All three now find it. The doc
  no longer claims a closed-form quartic.

- **Degree elevation keeps the continuity at every knot.**
  `elevate_degree` left each interior knot at full multiplicity, a C0
  joint where the curve was C2, against its doc; the knots it introduced
  now come out exactly, each interior knot one higher than it was.

- **Conics and helices answer every derivative order.** A circle's,
  ellipse's, hyperbola's and helix's derivatives past the second (the
  helix's past the third) came back zero, as did a 2D circle's and
  ellipse's; they are now the closed forms at any order. A parabola's and a
  line's higher orders are zero, as they should be.

- **A near half turn and a small chart keep their digits.**
  `Quaternion::between` refused directions within some 5e-4 rad of a half
  turn (a length tolerance held against a cosine) and lost digits short of
  that; it now reads the angle off the cross product and refuses only within
  the angular tolerance. `SurfacePoint::normal` held the tangents' cross
  product, an area, against a length and refused a healthy chart a tenth of a
  millimetre across.

- **A large box measures at once and exactly.** The exact surface
  integral cut every chart rectangle into quarter-turn panels in both
  parameters, so a plane's millimetre parameters gave a metre cube some
  four hundred thousand panels: seconds per measure, longer at three metres,
  and drifting in the tenth digit. A length parameter now takes one panel
  and an angular one at most sixty-four.

- **`make_wedge` documents its zero top extents.** The doc said a zero
  top extent is refused; it builds a ridge or a point top, five faces either
  way, and the doc now says so.

- **A mesh of positions alone welds.** `Triangulation::welded` indexed
  the normals and parameters of every vertex and panicked on a mesh read
  without them; the optional arrays now follow only where they are full.

- **A converted torus measures exactly and divides into spans.** A torus
  stated as a spline far from the origin had its closing seam sampled back
  onto the chart's start, so its seam sides were fitted with ringing and the
  volume came out off by a few tenths of a percent; a spline seam that is a
  chart line now takes the exact measure, also when it is cut into pieces.
  `divide_by_continuity` walked a torus's ring of seams with its first seam
  on the wrong chart side and never cut it, so `to_bezier` left one
  multi-span face; it now cuts into sixteen, each patch held to its own span.

- **A drum shorter than a unit divides.** `divide_by_angle` and
  `divide_by_area` read a straight cut line's direction off the surface at
  parameters zero and one, past the window of a drum or cone less than a
  unit high, and refused. The line is read inside its own span.

- **A mesh reports the chord it meets.** On a surface curved both ways a
  grid cell's middle can stand off the surface by up to twice the chord,
  while `deflection_met` said the chord was met. Each triangle's middle is
  now measured against the surface, and the flag is false where any
  stands off by more than the chord.

- **A scaled body meshes to the chord where it stands.** A face is drawn
  in its own frame and placed after, so a body scaled by its placement
  was meshed to the caller's chord before the scale stretched it: ten
  times the chord at a tenfold scale, and mass properties two percent off.
  The chord is now taken where the face stands.

- **Tight bounds reach a sphere's and a torus's extremes.** In a general
  frame `tight_bounds` could come out up to 45 microns small on a sphere or
  a torus, its descent stopping short of the extreme. Their extremes are
  closed form and are now taken wherever they lie inside the face.

- **A placed face's medial axis moves with it.** `medial_axis` read a
  moved face's boundary in its surface's unplaced frame and lifted the
  axis back through it, giving an axis neither where the face stands nor
  where it was written. The frame is now placed with the face.

- **Lengths and points along a curve turn its kinks.** `curve_length`
  could fail across a polyline spline's corner, and `parameter_at_length`
  (with `points_by_count` and `points_by_spacing` on it) took the failure
  for a zero length and returned wrong points silently. Lengths are
  integrated span by span between a spline's knots, and a failure is
  reported.

- **A periodic spline is closed.** `is_closed` on a B-spline compared its
  end control points, which are its ends only for clamped knots: a
  periodic spline, and a surface extruded from one, read as open. A
  periodic spline is closed, and any other is closed where its ends meet.

- **A conical helix measures its own length.** `HelixCurve::arc_length`
  ignored the taper and reported a cylindrical helix's length; a conical
  helix's speed grows with its radius, and its length now integrates it
  exactly.

- **Scaled curves name the images of their points.** Under a scaling a
  planar trim of anything but a bare line (a trim of a trim, of an offset)
  kept its old window, a parabola its old domain, and a curve on a surface
  its old chart curve while the surface's length directions rescaled; each
  then named the wrong points. All three now follow the rescaled
  parameter, and the whole-shape conversion carries edge ranges the same
  way.

- **Tangent circles keep their digits far from the origin.**
  `circles_tangent_to_three` and `circles_of_radius_tangent_to_two`
  solved rows holding each target's squared distance from the origin, so
  a small circle ten metres out came back off by a percent of its radius.
  Both now solve about the targets' own middle.

- **A nearly circular ellipse is a circle, not NaN.** An ellipse whose
  minor radius exceeded its major by less than the confusion tolerance
  was accepted and reported a NaN focal distance and eccentricity. Such a
  minor radius is stored as the major.

- **Two circles' bisector bends round the smaller one.** `bisector` of two
  circles of unequal radii put the hyperbola's equidistant branch round
  the second circle whichever was smaller, the wrong branch when the
  first was. Its frame now points toward the smaller circle.

- **A cone narrowing up its axis measures its distance.** A negative half
  angle (a cone narrowing in `+z`, which tapered sweeps and chamfers
  build) is accepted, but `Cone::distance_to` measured such a cone's own
  points as off it. Its nappes now lean out from the axis whichever way
  the cone widens.

- **Cone and torus inversion past the apex and the axis.**
  `cone_parameters` returned a point's own angle and height, the chart
  point only on the cone's own nappe and only for a point already on it;
  `torus_parameters` did not round-trip on a spindle torus's folded half.
  Each now returns the nearest point's parameters, half a turn round where
  that point lies past the apex or the axis.

- **Trimmed planar curves meet as their bases do.** `intersect_curves_2d`
  sampled any trimmed curve, and missed a line tangent to an arc, a line
  crossing one twice just below the tangent, and a circle touching an arc.
  It answers on the bases in closed form and clips to the trims' windows,
  as the space-curve path already did.

- **A curve's tangencies and lying stretches come back once.** Where no
  closed form applies, `intersect_curve_surface` reported a tangency up to
  eighteen times, and a curve lying in the surface (a cone's ruling, a
  circle on a sphere) as hundreds of piercings. Piercings between which
  the curve never leaves the surface are now one contact: a tangency once,
  or a stretch in `lying`.

- **A sphere on a torus's axis meets it in exact parallels.** The pair
  went to the marcher, and a sphere seated against the tube came back as
  two fitted pieces of its contact, a hundred times past their stated
  tolerance. It has a closed form: each place the sphere's meridian meets
  the tube's is a parallel, one tangential where the sphere is seated.

- **A fitted section states the error its images reach.** A section
  fitted beside a cone's apex reported a tolerance two hundred times
  tighter than its pcurve on the cone stood from the curve: the error was
  read at the trace's samples, and the chart turns fast between them
  there. Each pcurve is now lifted through its surface at dense stations,
  and the stated tolerance covers the worst.

- **A cone through its apex is sectioned on both nappes.** The closed
  forms for a cone against a coaxial cylinder, a coaxial cone and a level
  plane answered on the chart's own nappe only: a cone face running past
  its apex lost its far-nappe section, two parallel-sided cones meeting
  past an apex read as apart, and a level plane past the apex went to the
  marcher. Each crossing is now reported, and the cone's pcurve places a
  far-nappe parallel half a turn round; a face stopping short of its apex
  drops it by its window.

- **Arcs of one circle overlap across its seam.** The overlap of two
  trims of one circle was clipped within a single turn, so a stretch
  crossing the seam was lost, and a pair walked opposite ways reported
  parameters outside the second trim. The correspondence is now carried
  round whole turns, and every stretch both trims share comes back inside
  both.

- **A spindle torus is sectioned on its folded half.** A torus whose tube
  swallows its axis has, in every meridian, two tube circles reaching past
  the axis, and the closed forms against a plane square to its axis, a
  coaxial cylinder and a coaxial torus read only one: the parallels on
  the folded half were dropped. Each now reads both, and the pcurve on
  the torus takes such a parallel half a turn round.

- **Endlessly nested glTF and VRML are refused.** Both readers parse by
  recursion, and a document nesting arrays, objects or nodes a hundred
  thousand deep ran the stack out and aborted the process. Nesting past
  256 levels is refused by name.

- **A PLY face naming a missing vertex is refused.** `read_ply` indexed
  the vertices with the file's face indices unchecked and panicked on one
  past the end; it refuses the file by name, as `read_obj` does.

- **Flat caps of a fine mesh far from the origin stay whole.** A finely
  meshed part stored in single precision a long way out broke its flat
  faces into dozens: their small triangles lean past the coplanar angle by
  the rounding of their corners alone. A triangle's allowed lean now
  widens by its corners' rounding (twice the mesh's quantum over its
  smallest altitude).

- **A mesh's void stays a void wherever its ray lands.** The converter
  decided nesting with one ray, and a ray through an edge two of the outer
  piece's triangles share counted the crossing twice, or not at all, so a
  void came back as a separate solid. A ray grazing a triangle's boundary
  is set aside and another direction tried.

- **Unusable converter options are refused.** `solid_from_mesh` read a
  NaN, zero or negative `coplanar_distance`, a negative `quantum`, and a
  NaN or zero `crease` as something else, silently. Each is refused by
  name, as the documentation said.

- **Whole spheres and tori from single-precision meshes check valid.** The
  converter widened a whole face's tolerance to its fit's deviation but
  left its seam, poles and vertices tighter, which `check` refuses. The
  widening now reaches everything the face bounds.

- **A point in a spherical void classifies out of the material.** The
  bound of a sphere face used inside out covered only its seam: its chart
  outline put both seam walks on one column, so `classify_in_solid_exact`
  took points in a void bounded by it as inside the material. The first
  seam walk of an outline takes the image its own record names.

- **A pointed cone converts up to its apex.** `to_nurbs` and
  `general_transformed_shape` widened a cone face's window by a margin
  that ran past the apex into the other nappe, and the converted cone came
  back 15 to 20% too big. The margin stops at the apex.

- **A hole wound like its outline revolves hollow.** `make_revolution` read
  the hand of the outline alone, so a hole wound the same way swept walls
  facing into the material and added the hole's volume. Each wire's turn
  is measured, as the prism does: the outline must turn one way and every
  hole the other.

- **Faces welded within microns are not kept as a membrane.** Two boxes
  stacked a few microns apart, closer than the boolean welds vertices,
  fused into one shell holding both facing walls on the same edges: a
  wall of no thickness that `check` passed. Such a pair, bounded by the
  same edges and lying on each other facing apart, is dropped as the
  exact contact's opposed pair is, and a shell left with nothing else
  goes.

- **A solid poking microns through another is not taken as nested.** Where
  the general boolean refused, the nesting fallback read the boundaries at
  a mesh's samples, a tenth of a millimetre apart, and took a drum two
  microns wider than its box as lying within it. Each inner face now
  climbs its chart from the samples nearest each outer face, and a point
  found outside refuses the nesting.

- **Fillets whose balls overlap are refused.** `fillet_edges` accepted two
  edges of a narrow face (or a short wall's two rims) whose bands overlap
  across it, and built a solid that was neither rounding. A chain whose
  contact on a shared face stands inside another edge's band is refused;
  edges meeting at a corner still meet there.

- **A fillet the faces cannot hold is refused.** `fillet_edge` never
  checked that the ball fits: a radius of 100 on a 2 mm box deleted the
  solid and reported success. It now checks each face reaches the ball's
  contact, as `fillet_edges` did, and the check itself no longer skips a
  radius so large its side test stepped off the face.

- **An offset through the part is refused.** `offset_shape` and
  `make_thick_solid` returned inside-out solids when the offset or the wall
  ran a part's faces through each other; they now refuse, as `offset_faces`
  and `move_faces` did.

- **Located solids offset where they stand.** `offset_shape`,
  `offset_faces`, `move_faces` and `make_thick_solid` rebuilt a solid
  moved by a location in its nodes' own frame: wrong geometry that passed
  `check`, or refusals. Such a solid is baked into world coordinates
  first, as an instanced one already was.

- **Over-repeated end knots are refused.** `KnotVector::new` accepted a
  domain end repeated `degree + 2` times, where evaluation gave NaN; it is
  refused as a malformed vector.

- **Unclamped splines split into Bezier segments and elevate unchanged.**
  `to_bezier_segments` and `elevate_degree` assumed clamped ends and moved
  an unclamped or periodic curve by millimetres. The domain's ends are now
  raised to full multiplicity like every interior knot, and each segment
  takes the control points its span reads.

- **Unclamped and periodic splines reverse and rescale unchanged.**
  `KnotVector::reversed` and `reparameterized` clamped every knot into the
  domain, changing a spline whose outer knots lie beyond it by up to
  millimetres. Only the knots at the domain's ends are snapped now.

- **Sections winding round a periodic surface fit.** The closed sections of
  two crossing cylinders at an angle (an oblique pipe tee) came back
  millimetres off both surfaces with a stated tolerance in the hundreds:
  each loop's image on the branch pipe ends a period from where it
  begins, so the loop fitted as an open trace and stalled. Such a loop is
  now fitted closed in space with its images' join smooth across the seam
  (`fit::fit_points_joint_winding`) wherever the plain fit misses.

- **Spheres come back whole through STEP and IGES.** A sphere face bounded
  by its seam walked both ways, as an exchange file states it, read back
  with no region to mesh: whole spheres came back empty, and a sphere cut
  on its side or fused with another lost most of its face. Where a
  boundary meets a pole at one vertex from two sides of the chart, the
  readers now close it with an edge across the pole row
  (`closed_at_poles`). Both writers also dropped a solid's voids; they now
  write them (`BREP_WITH_VOIDS`, and the void shells of IGES entity 186),
  and the IGES reader honours the shells' orientation flags.

- **glTF read under a mirroring node stays right side out.** A node whose
  transform reflects (a negative scale, or a matrix with a negative
  determinant) left the triangles wound inside out and turned every normal
  inward; the winding is turned back and the normals follow the map.

- **Booleans with mirrored curved solids.** A mirrored cone gave wrong
  results that passed `check`, and mirrored cylinders, spheres and tori
  were refused. Baking a reflected placement restates each analytic surface
  on its frame's right-handed twin, and a seam that runs against its
  chart's direction gets its pcurves inside the chart window.

- **Primitives in a mirrored frame are sound.** A box, cylinder, cone,
  torus or wedge laid out in a left-handed frame came out inside out, and a
  sphere failed to build. A round primitive is built on the frame's
  right-handed twin about the same axis, and a box or a wedge on the frame
  with `x` and `y` exchanged, which spans the same corner; the box's faces
  keep the names the caller's frame gives them.

- **A mirrored cylinder or cone converts to a spline turning its way.** The
  section circles were built on a right-handed frame, so a mirrored
  cylinder's or cone's spline turned the other way and its normals
  pointed in; they are built on the surface's own frame.

- **An offset surface's closed form is the offset.** `OffsetSurface::analytic`
  rebuilt an offset plane on a right-handed frame, turning a mirrored
  plane's normal over, and an offset cone on its basis's frame, naming
  points a little along the axis from the offset's; the plane keeps its
  frame's hand, and the cone's frame rises with the offset, so each names
  the offset's points at the offset's parameters.

- **Mirrored revolutions and offsets are the images of the originals.** A
  surface of revolution carried by a reflection turned the other way round
  its axis, and offset curves and surfaces landed on the other side of
  their basis; each now names, at every parameter, the image of the point
  it named. `Transform2::preserves_handedness` no longer reads a negative
  scale, a half turn in the plane, as a reflection.

- **A negative scale reverses directions.** `Transform::apply_direction`
  (and `Transform2`'s) left a direction unchanged under any uniform scale;
  a negative one is a point mirror with a scale, and reverses it, so
  circles, planes and frames scaled by a negative factor keep the sense
  their points do.

- **`fix_shape` turns a located face in place (#101).** A face turned right
  way out was rebuilt on a new node without its location, so a void whose
  faces are placed by one (an assembly part read from STEP) moved out of its
  solid and still checked broken. The turned face keeps its location.

- **A fuzzy boolean takes its fuzz as the faces' coincidence.** Two faces
  lying within the fuzz of one another over what they share are one
  surface to `fuse_fuzzy` and `cut_fuzzy`, as they are within the
  tolerances the faces state: a prism of a converted part's top face along
  its rounds, left as facets each leaning a few microns, fuses and cuts at
  a fuzz wider than the lean.

- **A pad beside a wall leaning a few microns off it is sound.** Where
  two exact edges (lines, circles) come within a measuring width of one
  another, the boolean now holds them to their own tolerances, and takes
  them for one curve only over a stretch ending where one of them ends:
  a pad's side meeting a wall along its top edge and parting from it by
  a lean was read as running along the wall's edge for a tenth of a
  millimetre, and the pad's bottom edge was left short of the wall.

- **The kernel builds on its declared minimum Rust, 1.90, again.**
  `ogeom-bool` used an `if let` guard on a match arm, which Rust 1.90
  refuses; 0.6.0 built only on newer compilers.

### Added

- **A history says which faces came through a boolean unchanged.**
  `History::copy` records a modification that is an exact copy (a new node
  on the same surface within the same boundary) and `History::copy_of`
  reads it back. A boolean reports so every face the tool never touched,
  and the sew every face it only moves onto shared edges; composition keeps
  the mark only where every step copied the shape or left it alone. An
  application holding names, meshes or bounds for a face carries them
  across without comparing geometry.

### Changed

- **A boolean on a face with hundreds of holes is several times faster.**
  Splitting a face picked its probe points with a scanline per vertex
  height crossing every segment, and placed each hole by asking every
  pair of cycles about every other; on a plate with four hundred bores a
  small pocket took 580 ms, most of it on faces the pocket never touches.
  The scanlines now stop once their choice is settled and meet only the
  segments at their height, each hole's containers are found once, and
  the material test skips polylines the ray cannot cross. The probes and
  pieces are the same; the pocket takes 170 ms.

## [0.6.0] - 2026-09-30

A minor release: `cargo semver-checks` against 0.5.1 finds one breaking
change, a field added to `MeshSolidOptions`, and every crate packages and
verifies in a publish dry run. A mesh now converts in two steps if wanted,
as it is and then refined (`refine_solid`), and recognition is sounder
throughout: noisy flat faces, rounds tangent to their faces, thin rims,
crests and single-precision input. The boolean holds faces to the
tolerances they state, so a pad sketched from a converted part's face
fuses, cuts and commons with it (#99, #100).

A known limit: a pad whose sides run along a solid's walls at a tiny angle
(a converted solid's facets, tilted by the mesh's rounding, under a pad
sketched from them) leaves a sliver the boolean does not yet close; a
refined solid's walls are not facets, and pads along them fuse.

### Changed

- **`MeshSolidOptions` has a new field, `keep_vertices`.** Code building
  the options with a struct literal names it, or takes the rest from
  `MeshSolidOptions::default()`.

### Added

- **A converted solid refines in a second step.** `refine_solid` finds
  the cylinders, cones, spheres and tori among a solid's facets and
  rebuilds them, as `solid_from_mesh` with recognition would: a mesh can
  be converted as it is first (`recognize: false`) and refined on top.
  `MeshSolidOptions::keep_vertices` keeps every mesh vertex a vertex of
  the solid, so the first step keeps the samples along curved stretches
  that merged edges would drop. A recognized region lying on a larger
  neighbour's surface of its own kind now joins it, where two fits of one
  cone parted over their few vertices' slop.

- **STL precision reaches the mesh converter.** `ogeom::io::stl::read_with_quantum`
  returns the mesh with its encoding's rounding: single precision at the
  largest coordinate for binary files, the finest printed digit for ASCII.
  `MeshSolidOptions::quantum` takes it (or `single_precision_quantum` for a
  mesh held in `f32`), and the default coplanar distance is then never
  below twice it. An ASCII file written to six significant digits a
  hundred millimetres from the origin kept its fillets only as facets; it
  now comes back with cylinders.

### Fixed

- **A face met along one line by two faces of the other solid is split
  there once.** Where a section runs along one face's own edge while
  reading inside both faces by a hair, it now counts as that edge's split
  of the other face, as one reading outside does, and a second face
  meeting the other along the same edge (a chamfer beside the face it
  bevels) gives way to it. A prism of a chamfered slot plate's top face
  pushed through its bottom had its slot walls cut twice along the bottom
  rim, and the boolean left an open shell.

- **A pad sketched from a converted part's face fuses, cuts and commons
  with the part.** A face states how far it may stand off its surface (a
  face fitted to a mesh, by its fit's deviation), and the boolean now
  holds two faces to that: they lie on one surface where every point of
  the one over the other stands within what the two state, though their
  surfaces differ (a sketch keeping its points in single precision draws
  an arc through them on another centre and radius, a few microns off the
  wall it was taken from). The edges carried between the two, the
  pieces' partners on the shared surface and the pieces one stands in for
  are all read at that tolerance, and an edge running past the other
  surface's window is carried over the stretch it has on it. A pad from
  such a sketch along a converted part's walls was refused, a piece on the
  other's boundary with no partner face. A box and a sheared copy sharing
  a plane with it, refused as unresolved same-domain contact, now fuse,
  cut and common to their exact volumes.

- **A section leaving a face's boundary all but along it splits the face.**
  The boolean orders the strands round each junction of a face's
  arrangement by where they leave a small circle about it, not by their
  first steps: a section starting at a tangent can set out a hair outside
  the boundary before turning in, and read by its first step it crossed
  the boundary and the face was left with no piece at all.

- **A round's crest meeting a flat face is a circle.** A mesh's rim where
  a torus meets a flat face along its tube's crest, a few slops off the
  crest after rounding, was taken for a free curve: the circle's radius
  turns infinitely fast with its height there. A ring within the reach of
  the crest is now the crest circle, and the torus is built as a band
  between two circles rather than falling to facets.

- **A thin disc's rounded rim comes back a torus from a mesh.** When a
  curved region shed the flat face's facets it had taken in and the refit
  of what was left found nothing, the whole region went back to facets;
  its rim's normals, one-sided once those facets were gone, misled the
  fit. The surface already held now stands where it still holds the rest.

- **A round between two flat faces comes back tangent to them.** The mesh
  converter puts a curved region that meets two non-parallel planes across
  smooth edges on the cylinder tangent to both, its axis along their line
  of meeting and its radius the median of the vertices' own, wherever that
  holds every vertex within the coplanar distance. A narrow round's few
  rows of vertices fitted a cylinder standing a micron off its faces.

- **Noisy flat faces come back whole from a mesh.** The mesh converter's
  default coplanar distance also follows how far the mesh's vertices stand
  off the flat faces they lie on (measured across nearly coplanar edges),
  and a flat face is grown against a plane refitted to what it has
  gathered rather than its first triangle's, whose lean carried a large
  face's far side past the distance. A drilled block exported with noise
  came back with its top split in two; it now comes back its seven faces.

- **Rounds met by flat faces along a tangent come back exact.** A flat
  face's facets drawn across the tangent circle it shares with a round have
  their corners on the round; they are now peeled from the round's region
  (which is fitted again without them) instead of making the whole region
  fail. Two regions left on one surface and sharing an edge are one. A
  plate with rounded corners and a rounded bottom rim comes back on its
  eighteen exact surfaces, four tori among them, where it took a thousand
  faces; a thin disc with a rounded rim on its four.

- **A converted dome on a cylinder is sound.** A sphere centred on the axis
  of the cylinders and cones beside it kept its own angular origin, so the
  cap's seam and the cylinder's started a quarter turn apart on the circle
  they share; the solid failed the checker, and on a coarser mesh the
  sphere was not built at all. The sphere now takes the axis's frame.

- **A curved region no longer takes a flat face's facets.** The mesh is
  first cut into flat patches, and a curved region may take a patch much
  larger than its own facets only whole, with every corner on its surface.
  A round met along a tangent by a long flat wall took the wall's facets
  whose near corners lay on it; a finely meshed plate with rounded rims
  came back in a thousand faces and now in five hundred. A plane beside a
  recognized face takes its boundary from it, so the face must also keep
  its flat neighbours on their own triangles, or it is faceted.

- **A recognized region stays one connected piece.** Dropping the triangles
  of vertices a fit refused could leave a region in pieces joined by
  nothing, each then built as a face on the first fit's surface. The
  largest piece is now the claim, and the others stay free to seed regions
  of their own or fall to the planes. A perforated cup's conversion lost
  three spheres its mesh does not hold, and came to within 0.02% of its
  volume from 0.17%.

- **The mesh converter checks each body it builds, and gives up only what
  fails.**
  - A recognized face pointing against the triangles it replaces, or
    covering a part of its surface they do not (sampled on the face, each
    point within the rise of the triangles near it), is faceted and the
    part rebuilt.
  - A body that does not hold its mesh volume gives up its own recognized
    regions; other bodies keep theirs.
  - The volume is measured at the facets' mean offset rather than a
    millionth of the part: a converted part of four thousand faces took
    55 s and takes 5 s. A thin disc with a rounded rim, meshed coarsely,
    made the conversion fail outright; it now comes back a sound solid.

- **A mesh of bodies glued along an edge converts to closed solids.** An
  edge four triangles use (two blocks touching along it, as exporters
  write glued parts) paired nothing, and both bodies came back open
  shells. The bodies are now what the edges used exactly twice join, and at
  a crowded edge each body's own two triangles pair up.

- **`solid_from_mesh` no longer builds faces over the wrong part of their
  surface.** A recognized region whose fit refused a few vertices lost the
  triangles bringing them, which could leave it in patches joined by
  nothing; its face was built from one patch's boundary and the solid
  passed the checker with the wrong volume (0.85% on a perforated cup).
  - Each patch of a region is now a face of its own.
  - A face whose bounds reach well past its triangles is faceted and the
    part rebuilt.
  - The finished solid must hold the mesh's volume within what its
    recognized surfaces may add over their facets; past that recognition
    is withdrawn and the mesh comes back faceted, flagged by the new
    `MeshSolidReport::recognition_withdrawn`.

- **`solid_from_mesh` unfolds slivers folded under their neighbours.** An
  exporter that moves a vertex across a thin triangle leaves it facing
  against all three of its neighbours while its winding agrees with
  theirs; it became a face pointing into the material, and the solid
  failed the checker. The diagonal it shares with a neighbour is now
  swapped, choosing the best-shaped pair that faces with the rest.

- **A pad through a converted part (#100, in part).** A pad cut from a
  converted part's own face and pushed through its walls and chamfered
  bottom now fuses, cuts and intersects, matching the same booleans on the
  exact part. Six causes, each where the converted part's single-precision
  slop met a rule written for exact geometry:
  - An edge lying on the other solid's surface within the two faces'
    stated tolerances takes its closed-form image at that looseness.
  - One cylinder on two frames, where the fitted image refuses, has the
    owner's pcurve carried across the charts exactly.
  - An exact section meeting a line or conic boundary edge on a curved
    face counts a near miss within that edge's tolerance as a crossing.
  - A contact measured along an edge ends at the edge's vertex, not a
    measuring width past it; the overshoot left one face a sliver its
    neighbour did not have.
  - A contact strand is folded onto the face's branch with the leaning
    ray: its middle sits level with the vertices it paved, where the level
    ray's parity is wrong.
  - A prism's rails take the tolerance of the vertex they sweep.

- **A chamfer along a rounded rim removed too much; a fillet refused it.**
  A chamfer on an arc of a rim (a rounded corner's quarter) cut a wedge
  turned the whole way round the corner's axis, taking material on the
  far side too: the solid passed the checker with the wrong volume. It
  now turns through the arc alone, as the fillet's does. The probe that
  tells a rim that stops from an arc split off a whole rim read the arc
  the wrong way round when the cap faced against the circle's normal, so
  a fillet along the bottom rim of a rounded block was refused; it now
  reads the arc's own middle. Both blends along a chain of lines and
  quarter arcs now match the closed-form volume.

- **`solid_from_mesh` recognizes one-row chamfers and exporter noise.**
  - A chamfer meshed as a single row of triangles has its corners on two
    circles, which a sphere and a torus fit as well as a cone. The cone
    is now taken, and a cone fit reads its lean from radius against
    height.
  - Coaxial surfaces now share the axis of the best determined one: a
    cylinder's before a cone's, and a larger region before a smaller. A
    fit that does not hold on the shared axis is fitted again about it.
    Neighbouring rims then meet on one circle.
  - An axis within 1e-3 of square to a plane of the solid is made square
    to the largest such plane, where the region still fits. Rims on that
    plane are then exact circles, not fitted splines, and a feature
    sketched on it extrudes into exact cylinders.
  - When no `coplanar_distance` is given and the curved fits press
    against the default, the distance widens (up to four times). A mesh
    whose exporter moved its vertices a little further no longer leaves
    thin planar slivers at the foot of each fillet.

- **A pad within a converted solid fuses (#99).** A pad sketched on a
  converted part's own outline and sunk into it runs along the part's
  walls and touches its recognized ones at every vertex, within the
  tolerance the conversion states; the boolean read that contact as
  ambiguous or as edge contact and refused. Where it refuses, `fuse` and
  `common` now ask whether one solid lies within the other at the
  looseness the solids state, and answer with the outer one or the inner
  one. A piece whose probe no ray can read is also asked at its next probe
  before the boolean gives up.

- **A traced section round a whole loop is no longer taken for a wandering
  trace.** Two surfaces meeting all the way round (a near-coaxial cone and
  cylinder), whose faces hold only an arc of the loop, gave a section
  longer than four times the faces' extent, and the guard against traces
  lapping beside a chart's pole refused it. The guard now allows a turn
  round the section's own extent as well.

## [0.5.1] - 2026-09-29

A patch release: `cargo semver-checks` against 0.5.0 finds no change
needing a bump in any of the fifteen crates, and every crate packages and
verifies in a publish dry run. It adds direct editing of a solid's faces
(`offset_faces`, `move_faces`), and a mesh converted with recognition
keeps its recognized faces on the triangles they replace.

### Added

- **Direct editing of a solid's faces (#98).** `offset_faces` offsets the
  named faces along their outward normals and `move_faces` moves them by a
  translation or a rotation; the faces around them follow, each staying on
  its own surface with its edges re-derived where the moved faces now meet
  it. A plane moves parallel, a cylinder, cone, sphere or torus changes
  radius; several faces move together; the history maps every face to the
  face it became. An edit that runs a face past the faces across from it
  or into the rest of the solid is refused by name.

- **`check_self_intersection_near`** tests only the face pairs holding one
  of the faces named: after an edit, the pairs it left alone cannot have
  begun to cross.

### Fixed

- **`offset_shape` of a bored block and of a ball.** A seam reached
  reversed (a bore a boolean cut) was rebuilt the wrong way round, and a
  sphere's poles and two-pole band had no rebuild; both now offset.

- **A surface recognized in a mesh stays on the triangles it replaces
  (#97).** A region was accepted when its vertices lay on the fitted
  surface, so long facets across a flat stretch between two rounded
  corners could be taken for a torus through the corners, whose face then
  stood far outside the part and turned the classifier's rays: the solid
  checked inside out. Every triangle of a region must now stay on the
  surface as closely as its own turn allows, and an edge between two
  recognized faces is not traced where either surface turns through a
  right angle along the mesh chord it replaces. `tight_bounds` no longer
  counts a face's surface beyond the chords its chart triangulation runs
  on, which put a dome a few hundredths below the edge that bounds it.

## [0.5.0] - 2026-09-29

A minor release, since one public enum gains a variant: `cargo
semver-checks` against 0.4.1 finds `PipeLaw::Fixed` added to an exhaustive
enum, and nothing else in the fifteen crates. It adds helical sweeps of a
profile off the axis, pipes and skinned lofts that close to a point, and a
pipe law that carries its section by translation. Unevenly scaled spheres
and revolved solids now mesh, measure and cut as the ellipsoids they are,
through their axis and on their seam's plane included; a profile turned
past a quarter turn revolves outward; and `check` finds a solid inside out
as a whole.

### Changed

- **`PipeLaw` has a new variant, `Fixed` (#90).** A `match` over
  `PipeLaw` needs an arm for it.

### Added

- **Helical sweeps of a profile off the axis (#89).**
  `make_helical_sweep` takes a planar face in any plane clear of the axis,
  not only one through it: a level profile square to the axis climbs as a
  ramp or a coil of a flat section. The caps are the profile's plane and
  that plane turned with it. A taper is taken for a profile through the
  axis or square to it; a profile the axis runs through, or one the screw
  carries partly forward through its plane and partly back, is refused.

- **A multisection pipe closes to a point (#92).** `make_pipe_sections`
  takes a vertex as its first or last section where the spine starts or
  ends: the section shrinks onto the spine and the skin closes on the
  point, uncapped there. Down a straight spine it is the exact pyramid or
  cone over the section.

- **`PipeLaw::Fixed` (#90):** the section keeps the frame it has at the
  spine's start, carried along the spine by translation and never turned.
  Each section edge sweeps the exact B-spline `edge(u) + spine(v)` (the
  control points summed, the weights multiplied), so nothing is fitted. A
  spine that crosses the section's plane one way and back the other is
  refused, as is one running in the plane.

### Fixed

- **A revolved solid scaled on every axis meshes (#96).** Converting a
  face with a pole under a general transform looked for the pole's row on
  the new surface at the pole's old place, so a scaling that moved the
  pole put its degenerate edge on the wrong end of the chart and the face
  enclosed nothing. The row is now found at the pole's mapped place: a
  hemisphere scaled to radii 8, 5, 3 meshes and measures as half the
  ellipsoid. A pcurve refitted by the conversion now holds its target in
  space, through the chart's stretch and between its samples, so a quarter
  turn scaled that way is integrated exactly instead of falling back to a
  mesh 4% light.

- **A profile turned past a quarter turn about the axis revolves
  outward.** A profile placed by a rotation had its missing pcurves read
  from its curves as placed on its surface as stated, so past a quarter
  turn `make_revolution` built its solid inside out, and at exactly a
  quarter turn refused it as lying on the axis. The curves are now read in
  the surface's own frame.

- **Cutting a solid of revolution or an ellipsoid through its axis.** A
  plane through a revolution's axis now meets it in the profile turned
  onto the plane, exactly, where marching stalled at the poles. A plane
  holding whole columns or rows of a spline surface meets it in those iso
  lines, so an ellipsoid from a scaled sphere halves on the plane of its
  seam. A point is classified against a revolution whose profile crosses
  the axis at either of its two places in the chart, so the half of the
  profile a face does not use no longer reads interior points as outside.

- **`check` finds a solid inside out as a whole.** Faces whose flags all
  agree with each other were taken as facing out, so a shell turned
  inside out everywhere read as valid. Where the flags agree, one face of
  each shell is now probed against the material; a shell found inside out
  is reported face by face, and `fix_shape` turns it back.

- **A fillet chain holding closed edges rounds each (#94).** A piece of a
  rim the solid holds in parts is rounded the whole turn round, so the
  rim's other pieces (a second half passed in, or a whole rim a boolean
  split in two on an earlier step) were gone when their turn came and the
  chain refused as though its members interfered. A member lying on a
  circle the chain has already rounded now counts as done. Eleven fillet
  stress cases refused before now round.

- **Booleans of whole revolutions whose seams coincide off the profile
  plane (#95).** A piece of a face was classified at a probe placed by the
  widest scanline interval across it; a scanline in the thin gap just past
  an arc's chords, and inside the arc, crosses no chord of it and reports
  an interval as wide as the piece, with its midpoint in a hole. Room is now
  the lesser of the interval's width and its scanline's gap. Four fillet
  stress cases refused before now round.

- **A skin narrowing to a point along a bend was built inside out.** Its
  wall was turned away from the centroid of every row, which a bent skin
  puts outside itself; it is turned away from its own section instead.

- **A skinned loft through several sections closes to a point (#91).** A
  family fit refining past one control point per datum left spans with no
  data and gave up on its best earlier round, which missed the tolerance;
  it now takes the interpolating spline there. A loft's sections are read
  off chords kept within a micron of their curves instead of 64 per edge,
  so a loft of circles is a solid of revolution to its tolerance.

- **Unevenly scaled spheres and revolved solids are what they look like
  (#93).** A revolution whose profile runs on past its pole (a whole circle
  turned) was converted to a spline over the whole profile, and its pole
  put at the chart's end instead of where the surface reaches it; each
  face is now converted over the rows it uses, a pole is placed where the
  surface reaches it, and under an affine map whether the conversion
  turned the normal is measured before the map. A partial sweep's two
  sides share their edge on the axis through the conversion. A boolean
  with a sphere spelt as a spline no longer tears a closed section at the
  seam, so a scaled sphere cuts as the ellipsoid it is.

## [0.4.1] - 2026-09-28

A patch release with no API change needing a bump; `cargo semver-checks`
finds nothing against 0.4.0. A boolean's history now reaches every face of
its result, and a fillet's or chamfer's history its blend or bevel face.

### Fixed

- **A boolean's history reaches every face of its result (#88).** The
  pieces were recorded, but sewing rebuilds them onto shared edges and
  its record was not carried on, so only faces the sew left alone were
  reached. A fillet or chamfer also records its blend or bevel face as
  generated from the edge it replaces.

## [0.4.0] - 2026-09-28

A minor release, since two public types change shape; `cargo
semver-checks` against 0.3.4 finds those two and nothing else. It follows
an adversarial review of the kernel: hostile files are
refused instead of aborting the process, results of several shells nest
by where they are, faces built inside out are found by `check` and
turned by `fix_shape`, and most of the time went out of booleans,
meshing, mass and fitting (the stress harness runs in 20 s where it took
143 s, every outcome the same or better). It adds revolutions up to a
face, arc joins for thick solids, pipe frame laws, corner styles and
multisection pipes, flat spirals, sections of faces and open shells, and
seams for imported periodic faces, with which booleans cut parts they
refused before.

### Changed

- **`FixReport` has a `faces_turned` field**, so a literal building one
  needs it.
- **`step::parse::Arg::Typed` holds its keyword and arguments boxed
  together,** `Typed(Box<(String, Vec<Arg>)>)` instead of
  `Typed(String, Vec<Arg>)`, which makes every argument 32 bytes instead
  of 56.

### Performance

- **`check` walks the shape once** for every type it examines, instead of
  once per type.

- **A boolean places each edge's curve once,** shared by the faces on
  either side of it, instead of copying and transforming it per face.

- **An exact mass integral holds at most a million samples.** Each run's
  comparison measures are summed as its samples are drawn; a run past the
  bound is drawn a second time into the result instead of held, the same
  samples in the same order.

- **Fewer surface evaluations per Newton step.** A march's correction
  takes the walked point and its gradient from the evaluation its system
  already made (four evaluations a step instead of six, the same values),
  and extrema read each surface's point and derivatives from one jet.

- **A walk's tangent is found without allocating:** the null vector's
  minors are expanded in place, with the same arithmetic, and each accepted
  state is stored once instead of copied.

- **A node's identity is a slot lookup,** not a hash: the table is a
  vector indexed by node, cheaper to copy for undo, and `shape_of` answers
  in node order.

- **Sewing bins its few loose edges apart,** so one edge with a wide
  stated tolerance no longer coarsens the twin search for every other.

- **A planar section skips faces wholly on one side of its plane** by
  their boxes, before intersecting their surfaces.

- **Meshing measures each chart edge's sag once** across the refinement's
  rounds and the triangles sharing it, and a fillet chain finds its next
  edges by projection only among edges whose boxes reach them.

- **A face's arrangement has no quadratic steps left.** Strand ends snap
  through a grid, dangling chains peel through a queue, each dart's angle
  is computed once before sorting, and hole nesting reads each cycle's
  nodes, area and box once and tests boxes before polygons. Every choice
  is the one the old scans made.

- **A march's seeding tests each cell only against cells it could meet,**
  binned by box, in their own order, instead of every cell of the other
  surface.

- **A STEP argument is 32 bytes instead of 56:** the rare typed value's
  keyword and list are boxed rather than sizing every argument of a file.

- **A boolean splits and classifies its faces in parallel,** each face's
  pieces and junctions joined in face order, and intersects a contact
  edge only with target edges whose boxes it can reach.

- **Checks and bounds repeat less.** `shape_bounds` bounds each shared
  edge and vertex once per call; the face-orientation agreement places
  each edge's curve once and keeps one golden-section point per round;
  `check_self_intersection` gathers each face once, skips pairs whose
  boxes stand apart, and keys adjacency by node rather than by a hash
  that could collide.

- **`cells` runs one boolean, not three:** the three cells are three
  choices of the same classified pieces. Pipes and evolved solids fuse
  their pieces in a balanced order, neighbours with neighbours, instead of
  each onto everything before it.

- **The intersectors' small Newton solves allocate nothing.** Curve
  crossings, curve and surface extrema, curve-surface polishing and the
  marcher's contact solve run on `newton_system_fixed`, the same damped
  iteration on stack arrays; a curve crossing reads each point from its
  derivative table.

- **Hidden lines test each sample against the triangles and faces in
  front of it only.** The mesh drawing bins triangles by where they
  project; the exact drawing skips faces whose projected box misses the
  sample or that lie wholly behind it. Mesh-to-solid finds the triangle
  across a sliver from an edge map instead of a scan.

- **The STEP reader borrows its arguments instead of copying them,** and
  the writer finds a coloured product's entities by node instead of
  scanning every written entity per product.

- **Projection onto a plane, cylinder or sphere starts from its closed
  form** where the foot lies in the window, instead of a 33 by 33 grid.
  `SolidMesh` meshes a solid once to classify many points; the fillet's
  run-out probes use it instead of meshing per probe.

- **Smaller savings.** A surface's derivative table is filled only to the
  total order asked for, each control row summed across `v` once; the
  thick-solid edge scan maps edges to faces in one walk; `ancestors_of`
  searches each candidate for the target's type only; `fix_shape` skips
  measuring an edge whose ends already stand further apart than the reach.
  A marched fillet's chord and the offset guard's finer retry follow the
  part's own size, not fixed lengths.

- **Spline basis derivatives run on fixed arrays below degree eight,** bit
  for bit the general recurrence without its nested scratch vectors: a
  third off the thread-groove tests. Curve crossings skip segment pairs
  whose boxes stand further apart than the reach.

- **Exact volumes of spline faces cost a third of what they did.** The
  chart integral reads a polynomial patch along each inner integral's own
  `v` from the net summed across `v` once, and a polynomial patch's first
  and second derivatives skip the rational quotient. `Surface::point_d1_at`
  gives the point with its first derivatives from one evaluation.

- **A boolean intersects its face pairs in parallel,** each pair's
  findings joined in pair order so the result is the same at any thread
  count. A parallel stage inside another now runs on the worker it lands
  on instead of spawning threads per outer item, and the machine's thread
  count is read once.

- **Spline fits solve their normal equations instead of inverting them.**
  An open curve's normal matrix is banded, and a banded Cholesky costs its
  size times the band squared where the inverse cost its size cubed. A
  system so near singular that the data does not decide a span's control
  points now reports itself singular, as an exactly singular one did,
  instead of solving to control points anywhere.

- **Curve crossings seed only where the curves come nearest.** Two curves
  running near each other (a boolean's section beside the edge it was cut
  along) put hundreds of segment pairs in reach, and each was polished by
  Newton to the same few crossings. Seeds are now the pairs nearest along
  either curve, and stretches the curves share are not polished at all.
  The stress harness runs in 16 s instead of 143 s with every outcome
  unchanged; one corpus drill went from 18 s to 1.5 s.

- **Walks, classification and evaluation do less of the same work.**
  `explore` no longer descends below the type it is looking for; a
  placed child is built in one step; a solid's prepared boundary answers a
  point near a face from what it prepared, not by rebuilding the face's
  trim; interior mesh points are tested against a banded ring index;
  de Boor evaluation and mass simplices use stack storage; the identity
  transform returns its input; the walker measures each tangent once; the
  boolean builds each face's chart outline and trim samples once.

### Added

- **`fix_shape` turns faces that face into their solid's material,** and
  reports how many (`faces_turned`); `inside_out_faces` names them.

- **`seam_periodic_faces`** gives imported periodic faces the seam or pole
  edge their chart needs: a band between rings that each go round once
  (the whole circle re-anchored and the other ring's crossing edge split
  where needed), a sphere's cap bounded by one ring, and a sphere's face
  reaching a pole with no edge along the pole's row. Each faces the way
  its loop runs round the chart. `reanchor_periodic_rings` runs it after
  the rings. `make_band_of_rings` builds a band from rings of any number
  of edges.

- **`native::compacted`:** a fresh model holding only what given roots
  reach, for reclaiming what a long session's discarded and failed
  operations left in the model.
- **`Document::set_undo_limit`:** at most 256 checkpoints are kept by
  default, the oldest dropped first.

- **`solve::newton_system_fixed`:** the damped Newton iteration for a
  fixed number of unknowns, allocation-free.

- **`SolidMesh`:** a solid's boundary meshed once, to classify many points
  against.

- **`Surface::point_d1_at`:** the point and first derivatives from one
  evaluation, for callers that need both at every sample.

- **`fit_surface_grid_sections`:** a grid fit with rows placed across by
  chord length, for sections lofted at spacings of their own.

- **`check` finds faces turned inside out.** A flipped face leaves every
  edge used twice, so nothing topological saw it and a mesh volume repaired
  it silently. Where the faces' flags disagree along their shared edges,
  each face is probed a step off both sides; a face whose outside is
  material is reported broken.

- **Revolutions up to a face (#73).** `make_revolution_until` turns each
  point of a profile about its axis until its circle first meets a limit
  face's surface, refused by name where a circle never meets it or meets
  it only past half a turn.
- **Arc joins for thick solids (#75).** `make_thick_solid_with` takes a
  `Join`: the walls extended to meet, or rounded about each edge convex on
  the side they grow toward, a cylinder of the wall's thickness about the
  edge and a ball at a corner.
- **Pipe frame laws (#76).** `make_pipe_shell_law` and
  `make_pipe_shell_with` take a `PipeLaw`: rotation-minimizing, Frenet, an
  auxiliary spine the section's normal points at, or a fixed binormal.
- **Extended and round pipe corners (#77).** `make_pipe_shell_with` takes
  `PipeCorners`: the mitre, legs run on past the corner and fused, or the
  section turned about the corner. Built for a face down an open spine of
  straight legs.
- **Multisection pipes (#78).** `make_pipe_sections` sweeps several planar
  sections down one spine, each standing where the spine crosses its plane
  and blended between them in the spine's moving frame.
- **Flat spirals (#79).** `make_helical_sweep` with a pitch of zero and a
  taper per turn turns the profile in its plane square to the axis while
  moving it out, refused where one turn would meet the last.
- **Pipe shells round straight corners are exact.** A face swept down an
  open spine of straight legs with corners is each leg's prism, run on
  past its corners, trimmed at the mitre planes and fused: planes and
  drums, measured exactly, whatever the profile's seam (#82).
- **Sections of a face or an open shell by a plane (#84).**
  `section_face` gives the edges where a face, shell or compound of faces
  crosses a plane, each trimmed to the face it lies on: exact where the
  surface and the plane have a closed form, fitted otherwise, a face's own
  boundary where it lies in the plane. `section` takes such a sheet
  against a half space bounded by a plane.

### Fixed

- **`canonical_simplify` rebuilds faces whose edges are used reversed.**
  An edge it left unchanged kept its occurrence's direction and had it
  applied again, so a bore spelt as a spline running the other way round
  failed with "edge 0 ends where edge 1 does not begin".

- **Mirror-placed bodies measure exactly.** The exact mass path refused
  every reflecting placement, so a mirrored part was measured from its
  mesh (a small drum 2.5% out at the default chord); a reflection is now
  taken with the chart's normal turned back. The bound-filter audit
  reports an unfiltered fill that refuses instead of failing the boolean.

- **A face's chart outline closes across the openings its edges' ends
  leave**, so a point level with one no longer reads inside the face. The
  boolean's paving kept a section piece off a narrow torus wedge that way;
  the bound-filter audit (`OGEOM_BOOL_AUDIT_BOUNDS`) now passes on the
  curved-corner blends.

- **Booleans on parts with seamless periodic faces.** Every face of a
  part is arranged, so one such face refused any boolean on it; three of
  the six stress drills into `nist_ftc_06` now cut. A spherical pocket
  imported facing out of the sphere comes out of healing facing in.

- **Spline walls whose pcurves stray past the chart's border by rounding
  are measured exactly.** The exact volume path and its orientation check
  bring such points back into the surface's domain instead of falling back
  to the mesh, which read a binormal-law pipe 1.5% high at the default
  deflection (#87).

- **`make_pipe_sections` matches each section's start and sense to the
  section before** in the spine's frame, so sections drawn from different
  seams or running opposite ways blend without a twist (#86).

- **An auxiliary pipe guide ending a hair short of the last station** (a
  sketch's end in single precision) is carried along its end tangent to the
  station's plane when within the sweep's tolerance, instead of refused
  (#85).

- **`check` names each part too tight for its bound once.** A vertex
  reached through both of a face's edges was reported twice.

- **STEP names past ASCII round-trip.** The writer spells them in `\X2\`
  escapes and doubles a backslash; the reader decodes `\X2\`, `\X4\`,
  `\X\`, `\S\` and `\\`, and reads raw bytes as UTF-8 where they are.
  A number that is not finite is refused rather than written as `NaN.0`.
- **Smaller fixes.** `section` shares one vertex where consecutive pieces
  meet, so the section is a wire's worth of edges. A pair of faces marched
  whole because one curve had no chart image no longer drops the pair's
  later curves or repeats its earlier ones. A strand end inside two
  junctions names the nearer. `check_tolerances` holds a face to its edges
  across the wire between them. Area and volume fall back to the mesh
  alike where the closed forms cannot evaluate, and both report a
  cancelled watch. Exact hidden lines draw a closed silhouette to the
  deflection instead of as an octagon. A face edge with only a pcurve is
  meshed at the chord the faces agreed on.

- **A loft through unevenly spaced sections sagged between them.** The
  skin's parameters across the sections ignored their spacing, so a square
  twisting through steps of 0.625 and 0.3125 measured a fifth light or
  would not fit. Sections are now placed across the skin by chord length
  (`fit_surface_grid_sections`). A section of a face whose ends two
  neighbouring faces place a hair apart shares one vertex that reaches
  both.
- **A frame-law pipe asked every edge of its spine for the whole station
  count,** and a law whose frame jumps doubled its stations twelve times.
  Stations are now shared out by length, and a law that outruns 2048
  sections is refused.

- **Trims scaled with their basis only sometimes.** A trimmed surface
  whose basis domain started at zero, and a trim of a trim of a line, kept
  their old parameters under a scaling. Trims are now carried affinely from
  the basis's old domain to its new one.
- **A tangency contact on a full drum read as outside its face.** The
  contact's trim test saw only one seam column; it now reads the same
  welded outline every other trim test does.

- **Faces built inside out, found by the new check.** A fillet that adds
  material (a concave seat) on a marched blend faced its blend band into
  the fill; a draft about an oblique neutral faced its wall inward whenever
  its hinge chained against the old surface's direction; a skinned tube
  round a bend in any plane but XY faced its wall into the bend;
  `to_nurbs` of a mirrored part turned its bore walls; canonical
  recognition of a spline facing its axis kept a flag the cylinder no
  longer agreed with. Each now measures the side it builds on.
- **Newton on an unevaluable point read as a root.** Where a surface or
  curve could not be evaluated mid-solve, the residual fell back to zero
  in the walker and to a made-up point elsewhere; it is now infinite, so
  the damped step backs off.

- **A part standing in another's notch was swallowed as a void.** Results
  of several shells were nested by bounding box, so a disjoint part inside
  another's box became a void of it, and an island in a cavity a second
  void. Shells are now nested by classifying one against the other, by
  depth: a shell inside a void is a solid of its own. `make_volume` shares
  the same nesting.
- **Two planes a microradian apart were called one plane.** Parallel was
  decided on one minus the cosine, a million times looser than the angle it
  was compared with. Planes are now parallel when they part by no more than
  the confusion distance across the region their windows cover.
- **Surface extrema never polished a farthest point, and reported stalls
  as extrema.** Seeds are the local minima and maxima of the near and far
  fields over each surface's lattice, thinned apart, and every polished
  point in curve and surface extrema is checked for stationarity by angle.
- **A line touching a sphere or a cylinder was missed about half the
  time.** A tangent discriminant that rounding leaves a few ulps below zero
  is now the tangency.
- **Smaller fixes.** Periodic curves and surfaces refuse a parameter that
  is not finite instead of answering a point that is not a number. Two
  boolean contact paves cap an imported edge's tolerance as the others do.
  `fuse` decides half spaces at the caller's tolerances. Sewing carries
  pcurves across in a fixed order, so repeated runs build identical tables.

- **Malformed and hostile files are refused, never aborted on.** The STEP
  reader no longer panics on entities with too few arguments, overflows its
  stack on entities that refer in a circle or on deeply nested argument
  lists, or allocates by the largest entity id, a knot multiplicity or a
  degree in the file. IGES reads its fixed columns as bytes (a label past
  ASCII no longer panics), bounds Hollerith lengths and every count by the
  entity's own parameters, and refuses transform and curve cycles. The
  native, BREP and glTF readers reserve no more than the file holds; 3MF
  archive offsets are bounded by the file and a deflated entry is refused
  the moment it inflates past its promised size. A placement's power costs
  its bit length, not its size.

- **Sweeps sampled a placed section where it was stored, not where it
  stood.** A section wire carrying a placement was read at its unplaced
  position when a pipe or loft sampled it.
- **A boolean between a helical sweep and a box whose face holds the helix
  axis did not close (#81).** The sweep's walls were fitted a quarter turn
  at a time with their borders on the profile's own plane and the planes
  square to it, so such a face met the walls only along those borders,
  where no section starts. The borders now stand an eighth of a turn off
  those planes, and the box keeps exactly half the coil either side. A
  junction between fitted sections now closes in the exact integral within
  both sections' tolerances or the vertex's, whichever is wider.
- **A helical sweep trimmed by a boolean measured about 3% out (#80).**
  The fitted sections the boolean leaves on its walls meet within their
  own stated tolerances, and the exact integral refused them for missing
  a millionth of the chart, falling back to a coarse mesh. Chart loops now
  close, and pcurves lie on their edges, within those tolerances, so the
  common and the cut measure their closed forms at the default deflection.
- **A groove cut into a bore drilled by a cylinder primitive failed
  (#74).** The primitive's wall is stored over exactly its own height, so a
  section crossing its rim stopped a hair short of it and the cut refused to
  close. The boolean now stretches a cylinder's window a little past an end
  the other face reaches beyond.
- **A half space clear of a solid sectioned it anyway (#83).** The box
  standing in for a planar half space reached twice the other shape's
  diagonal from the plane, so a plane further off than that put the box's
  far face inside the shape: a section returned edges and a common trimmed
  material. The box now reaches past the shape's far side.
- **Re-anchoring periodic rings turned rebuilt faces inside out.**
  `reanchor_periodic_rings` reversed a reversed face once when rebuilding
  it and again when rebuilding its shell, so the face came back facing
  into the solid. `check` passed it and a mesh volume did not notice, but
  a boolean's exact volume lost the material behind it: a drill through
  an imported plate measured 1.5% light. Found by the new robustness
  harness.

## [0.3.4] - 2026-09-26

A patch release with no API change needing a bump; `cargo semver-checks`
finds nothing against 0.3.3. It adds typed DXF reading, helical sweeps,
exact pipe shells along lines and arcs, the medial axis of any planar
face and the middle path of a pipe-like solid, dividing and remodelling
shapes, NURBS conversion within a tolerance, tight bounds, VRML reading
and more of STEP and IGES. Fillets pinch at tangent poles, round corners
where a curved face meets and close in more orders at apexes and curved
corners, and a fillet or chamfer too wide for its face is refused.

### Added

- **Typed DXF reading (#72).** `read_dxf_entities` reads lines, arcs,
  circles, ellipses, splines and polylines with their bulges and closed
  flags, the drawing's `$INSUNITS`, and which layers and linetypes are
  dashed. `read_dxf` now closes a closed polyline and no longer reads a
  `POLYLINE` header's elevation point as its first vertex.
- **The middle path of a pipe-like solid.** `middle_path` walks a solid
  from one named face to another and returns the curve through its
  cross-section centroids: a line for a straight tube, a fitted spline
  otherwise. Each point is the centroid of a section cut square to the
  path, holes subtracted. The result carries the measured deviation, held
  under the tolerance asked for.
- **Fillets pinched at a tangent pole.** A crease whose two faces turn
  tangent at one of its ends, such as the seam of two equal drums crossing
  where the drums touch, blends: the band narrows to a point at the pole.
  Both the fused drums' notch and their common part's ridge round, one
  seam edge at a time or all in one call.
- **Corners where a curved face meets round.** `round_vertex` rounds a
  convex vertex whose three faces include a drum, cone, ball or spline:
  one ball touches all three and its patch closes the corner. `fillet_edges`
  on the three edges of such a corner closes it with that ball.
- **Exact fillets where a wall meets a drum along a ruling.** A straight
  edge between two planes or drums parallel to it (the vertical edge where
  an extruded profile's line meets its arc) blends with an exact drum band,
  and no longer fails for want of a planar pair.

- **Pipe shells round a skew corner against a curved leg.** Where a
  spine turns a curved leg's end out of its plane, the sweep is built in
  pieces: each side runs on straight past the corner, is trimmed by the
  mitre plane and the pieces are fused. This was refused.

- **The medial axis of any planar face.** `medial_graph` returns the
  axis of a face with holes, reflex corners and curved edges as a graph of
  branches, each the exact bisector of its two boundary elements: lines,
  parabolas, hyperbolas and ellipses for segments, corners and circular
  arcs, and a fitted curve held to the tolerance for any other edge. Each
  branch point carries its clearance, and every branch is checked to keep
  no other boundary element nearer than its own two.

- **More of IGES reads.** A conic arc whose axes turn, parametric
  spline surfaces (114) as exact bicubic B-splines, offset curves whose
  distance varies, and trims given only in a B-spline surface's parameters.
  Model-space notes, leaders, labels, symbols and dimensions read as PMI
  callouts, a dimension whose text states a number also as its value.
  Subfigure instances place their shared definition, and levels and
  groups read as layers.
  Constructive solids read too: the primitives, solids of revolution and
  extrusion, boolean trees, assemblies and instances.
- **Small faces and small solids go.** `fix_small_faces` collapses a face
  that fits in a ball of the given size to a point, and a thin strip to one
  of its long sides, rebuilding the faces around it. `remove_small_solids`
  drops the solids enclosing less than a given volume.
- **More STEP geometry reads.** Curves: Bezier, uniform and quasi-uniform
  splines, hyperbolas, parabolas, trimmed curves, offset curves,
  polylines, composite curves and curve replicas. Surfaces: surfaces of
  revolution, offset surfaces, rectangular trimmed surfaces, curve-bounded
  surfaces, degenerate tori, surface replicas, the Bezier and uniform
  spline forms, and a rectangular grid of patches joined into one spline.
  The writer puts hyperbolas, parabolas, offset curves, surfaces of
  revolution, extrusions and offset surfaces out exactly instead of as
  splines.
- **An offset of an analytic surface is that analytic surface.**
  `OffsetSurface::analytic` returns the plane, drum, cone, ball or torus
  the offset is.
- **A trimmed surface converts to a B-spline** as its basis over the
  trim's window. `BSplineCurve::segment` and `BSplineSurface::segment`
  cut a spline to a piece of itself, parameters kept.
- **Dividing shapes.** `divide_face` cuts a face along a line of its own
  parameters, with the surface's exact iso-curve as the new edge.
  `divide_by_continuity`, `divide_by_angle` and `divide_by_area` cut edges
  and faces at weak knots, past an angle or above an area. `to_bezier`
  leaves every curve and surface a single Bezier span.
- **Every face divides.** A face with no boundary of its own (a whole
  ball or torus) is first bounded by its chart's sides, seams and poles
  included; a face on a surface with no closed-form iso-curve (an offset
  of a spline) is cut along one fitted at its own parameters.
- **IGES trims on analytic surfaces.** A trim given only in the
  parameters of a plane, cylinder, cone, sphere or torus is lifted through
  the format's parameterization, and those surfaces take the reference
  direction the file names for where their angle starts. So is one on a
  tabulated cylinder, and on a surface of revolution of a line or spline.
- **Remodelling.** `restate_geometry` rebuilds a solid with its surfaces
  and curves restated by the caller, every trim re-derived.
  `restrict_degree` refits every spline above a degree at that degree,
  within a tolerance. `swept_to_elementary` names an extruded or revolved
  line or circle as the plane, drum, cone, ball or torus it is.
- **NURBS conversion within a tolerance.** `to_nurbs_within` converts a
  solid exactly where a closed form exists and fits the rest (offset
  surfaces, helices, offset curves) within the tolerance given, where
  `to_nurbs` refuses. `SurfaceGeometry::fitted_bspline` fits any bounded
  surface.
- **A curve continues to a named point.** `BSplineCurve::extended_to`
  joins on a piece that carries the curve's end derivatives to the order
  asked and ends at the point. `BSplineCurve::extended` now meets the
  length asked exactly rather than to first order.
- **Global minima.** `global_minimum` finds the least value of a
  function over a box by branch and bound, certified to a tolerance under
  slope and curvature bounds it estimates; `swarm_minimum` is a particle
  swarm; `minimize_local` is a Nelder-Mead descent in several variables.
- **VRML reads.** `read_vrml` reads VRML 2.0 and 1.0 scenes into placed
  meshes: face sets and the box, ball, drum and cone primitives, with
  `DEF`/`USE`, transforms, switches and material colours honoured.
- **Helical sweeps.** `make_helical_sweep` turns a profile lying in a
  plane through an axis along a helix about it (a thread or a spring):
  every point runs its own helix, cylindrical or tapered, either hand.
- **Pipe shells along lines and arcs are exact.** A spine of lines and
  arcs meeting tangent sweeps as extrusions and revolutions of the
  section, fused, so its walls are planes, drums, cones and tori and the
  solid measures in closed form.
- **A wire swept round a skew corner.** A closed planar wire profile
  sweeps round a corner that turns a curved leg out of its plane as the
  walls of the face it bounds, where it was refused.
- **Skinned lofts keep corners and circles** (#57). Two sections loft
  ruled and exact; more with matching corners loft as one strip per
  edge, corners kept; coaxial circles loft as a solid of revolution. Two
  squares give the prism and two circles the drum.
- **Tight bounds** (#65). `tight_bounds` gives the smallest axis-aligned
  box holding a shape, each side where the shape actually reaches: a
  revolved tube of outer radius 6 bounds at 6, not at its carrier's 8.29.

### Fixed

- **A washer's medial axis did not settle (#71).** A face bounded by two
  full circles has no branch point, and its one branch, the ellipse
  between the circles, was never emitted. It now comes back as one closed
  branch naming a single vertex at both ends.
- **A fillet on a prism's base edge was refused as too wide (#70).** The
  check stepped into the reversed base face on the wrong side of the
  edge. Each face's inner side is now read from the face itself.
- **A rim blend beside a sphere's pole could not be built.** Where the
  ball's rail on the sphere passes over the pole while the rim does not,
  the sphere's leg was asked to be a band between two loops that both
  wind round the pole. It is now the band from the rail to the pole with
  the rim cut from it.
- **A section through a point where two surfaces touch was fitted loosely.**
  The tracer could report one sample's position with its neighbour's
  parameters there, and the joint fit stalled near 1e-5 mm at that
  sample. Such samples are dropped before fitting, and a torus against a
  cylinder now fits to about 3e-7 mm.
- **A straight band after a rim arc's band could not close at a curved
  corner.** Its section through the rim's torus runs tangent to the
  torus's end meridian into the corner, stopped short of the vertex, and
  split the meridian into a sliver. A tangent touch at an edge's end is
  now that end's vertex, and the section's chart image is bent onto it.
- **A flush fillet beside a rounded apex could fail to close.** A straight
  band end tangent to the corner sphere's rim was read as running along
  the rim for a short stretch, which split the rim at a vertex the band
  never had. A line and a circle now share no stretch. The four flush
  fillets close after an oblique apex's corner, and the five after an
  irregular pentagonal one.
- **A holed profile round a closed spine filled its hole.** Where the
  hole's ring wound the other way from the one the sweep assumed, its
  tunnel came out as material and the solid measured the hole's volume
  too many. Each shell is now turned by the volume it encloses.
- **A partial revolve could build its end caps inside out.** Where the
  profile's ring wound against its plane's normal, the walls followed the
  ring and the caps the plane, and the caps faced into the solid; a cap
  in a plane through the axis adds nothing to the volume, so only its
  orientation was wrong. The caps now face the way the prism's do and the
  walls are turned to meet them. Surfaces of revolution of a spline
  curve trimmed by their chart now measure exactly too.
- **Spline faces with any trim measure exactly.** Mass properties
  integrate a spline face, or a spline curve extruded or revolved, round
  its own trims on the exact surface, the panels broken at its knots; a
  converted box drilled through keeps its volume to rounding.
- **A helical sweep measured 0.65% light at the default deflection**
  (#69). Its fitted walls were measured from their mesh. A spline face
  trimmed by its own chart's borders is now integrated on the exact
  surface, knot span by knot span, and the thread measures to Pappus.
- **A chamfer wider than its face cut through it** (#67): a 12 mm bevel
  on a 10 mm face sliced a slab off the box. A setback running past a
  face's far side is now refused, naming the edge and the distance. A
  fillet is held to the same on any edge, curved faces and concave edges
  included: its ball's contacts, set back the radius times the tangent of
  half the faces' turn, must land on the solid's boundary.
- **A refined solid of pads measured from its mesh** (#64). The face
  `unify_same_domain` merged lay on one of the merged planes, whose window
  did not cover the union, and the plane refused to be read past it. A
  plane is now read anywhere; its window only says where its face was
  built.
- **A wall drafted a degree and a half could not be measured** (#66): its
  trim reached past the plane's window. The same change reads it.
- **A pad placed again at a scale measured and fused wrong** (#63): a
  cube padded from a square and copied 1.5 times the size measured 144
  where it is 216. An edge's trim stored under the pad's own placement
  was not found once the whole solid was placed again, and the trim of
  the edge's unplaced twin was read instead.
- **A partial revolve of a profile touching its axis failed** (#62). The
  edge on the axis is its own image; a wire now joins one vertex under two
  placements that land on one point, so the start and end faces close. A
  profile face built without trims gets exact ones before a prism or
  revolution keeps it as an end, and is read at its own middle rather than
  at its plane's origin.
- **A pipe shell failed on a profile ring walked by reversed edges**
  (#58). A reversed edge's ends were swapped twice.
- **A profile square to a curved spine's start was refused as leaning**
  (#59). The test now reads the spine's exact start tangent, to an
  angle's tolerance.
- **A circle piped along a line was 2.6% short of its cylinder** (#60),
  and a square along an arc 0.6% short of Pappus. Both are now exact.
- **A thread groove failed to cut from a bored block at some lengths.**
  Where the sweep had turned a flank just inside the bore, the flank met
  the bore's wall along a short run at a grazing angle, shallower than the
  sampled cells could see, and that section was missed. Sections are now
  also found from where they cross a spline's border, so 5.75 and 6.75
  turns cut.
- **Rebuilding a torus broke its trims.** `to_nurbs` and `baked_shape`
  took every seam for a `u` column, so a torus's `v` seam was given
  pcurves on the wrong edges of the chart and the face could not be
  meshed. A seam on the `v` rows now gets its trims there.
- **A pipe shell failed along an arc that runs on into a straight line.**
  One skin fitted across the join, where the curvature steps, and could not
  reach the tolerance. Each spine edge now skins its own run, the runs
  sharing the section at the join, and a run whose sections share a plane
  is built on that plane.
- **IGES curve-on-surface entities read their model-space curve from the
  wrong field.** The reader took the parameter-space curve as the
  model-space one; it now reads the fields in the order the format gives
  them.
- **A boolean on a prism split its near cap where its far cap was cut.**
  The far cap of an extrusion is its near cap moved, sharing its edges
  under a displacement, and the boolean shared edge splits by edge alone.
  Edges now count as one only at the same placement.
- **`analyse_blend` failed on a corner patch whose corner is its ball's
  pole.** A pcurve ending on the pole stood a rounding past the sphere's
  chart there; the station is now read at the chart's edge.

## [0.3.3] - 2026-09-25

A patch release with no API change; `cargo semver-checks` finds nothing
against 0.3.2. `make_half_space` accepts any face, so a cut or common can
stop on a drum, ball, cone, ring or spline patch. A thread groove cuts
from a plain or bored block. Holes drilled through large converted parts,
beside rounded corners, near a ball's rim and all but parallel to a bore
cut valid, and booleans on parts of tens of thousands of faces spend far
less time on bookkeeping.

### Added

- **Half spaces bounded by curved faces.** `cut`, `common` and `section`
  accept a half space from `make_half_space` whatever its face's surface.
  A drum, cone, ball or ring divides all of space and is exact on its
  surface; a drum and a cone are read as unbounded along their axes. An
  open surface (a spline patch, an extrusion, a revolution, trimmed or
  offset surfaces) divides space across its own extent: the other
  argument must lie across it, and the surface must not fold back along
  the direction into its material, or the operation says which. An
  extrusion can stop exactly on a curved face.

### Fixed

- **A thread groove would not cut from a block.** A section swept along a
  helix skins one wall per profile edge, each fitted at its own pace
  along the sweep, and a wall adopting its neighbour's rail took that
  rail's image on itself to be straight; along a helix it drifted a fifth
  of a millimetre off. The image is now read off the adopting wall. And a
  ray starting a few microns from a wall winding several turns crossed it
  in the gap between the wall and the flat cells seeding its crossings,
  so a point in the groove read as outside; a curve near a cell by its
  surface's bow from it is now seeded there. The groove cuts from a block
  for any number of turns.
- **A thread groove would not cut from a bored block.** Where the groove
  runs in and out of the bore, four things went wrong. A section winding
  several turns round the bore was fitted as one curve that wandered
  hundreds of millimetres between the points it was checked at; fits are
  now checked between their samples too, and a trace no single fit
  follows is fitted in pieces. Such a section crossed the bore's seam
  once a turn and was read as running along it. A section ending along a
  rail, or two meeting a rail that grazes the bore, ended beside the
  rail's split point by more than the weld reached in the face's chart;
  they now meet where the junction there says. And a ray crossing a wall
  near the end of its patch was solved from a seed that stepped off the
  patch, so a point in the groove's end read as outside. The groove cuts
  from a bored block for one to seven turns.
- **`make_half_space` could pick the wrong side of a closed surface.** The
  side was read from the normal at one sample of the face, which on a rod
  or a ball can face the given point across the surface. It is read where
  the surface comes nearest the point.

- **A hole drilled through a large converted part failed, and slowly.**
  One loose edge on a face set how far every end on that face was welded,
  so short edges and short sections collapsed whole. On a small sphere the
  weld was also taken in radians as if they were millimetres. The weld
  now reaches only as far as each piece's own doubt, scaled into the
  face's chart. Where a circle or ellipse crosses another in a different
  plane, the crossing is solved in closed form rather than by sampling.
  On a converted part of 25 000 faces a drill cuts in 6 s where it
  failed after 18 s.
- **A hole drilled beside a rounded corner often failed.** Where a hole
  crosses a fillet's tangent line, the section on the flat and the section
  on the fillet meet tangentially, and their crossing with each other or
  with the line landed up to a couple of microns off. The flat was left
  unsplit or a wire had two vertices where it needed one. Such a crossing
  now stops where the sections cross the shared edge, and a face welds its
  ends as far apart as the stops on its edges were found. A seeded scan of
  a hundred holes beside the corner of a rounded box went from a third
  failing to none.
- **A converted part's volume could not be measured.** Its faces meet on
  curves threaded through the mesh a few hundredths off, and the drawn
  mesh did not weld shut there, so the volume was refused as having no
  closed boundary. Closure is now asked of the topology, and where the
  drawn mesh stays open the volume is summed face by face from each
  face's own mesh. Cracks in a drawn mesh narrower than its chord are
  sealed, and triangles laid twice facing opposite ways are cancelled.
- **Booleans on parts of tens of thousands of faces spent seconds on
  bookkeeping.** Junctions were merged by comparing every pair, deleted
  faces were found by scanning every kept one, sewing projected edge
  middles before checking the ends matched, and a point classified
  against a solid asked every face. Each now looks only at what is near
  it. A line crossing a circle or ellipse is solved in closed form rather
  than by sampling. A drill through a converted part of 25 000 faces
  cuts in 2 s.
- **Sewing could open a face's wire.** Two edges within their own
  tolerance of each other were merged, but when their ends were two
  vertices farther apart than either vertex's tolerance, the faces
  around the dropped edge kept ending on its vertex, and their wires
  broke. Merged edges now end on one vertex, widened to reach both.
- **Holes through a converted part failed or took the wrong volume.**
  Two causes, together half of a scan of random holes through a real
  printer's part. A plane all but parallel to a hole's axis meets it in
  an ellipse kilometres long, whose short stretch inside the hole was
  missed by sampling, so the face was never cut. And a bore the converter
  opens along a slit, its chart starting at the slit rather than at the
  cylinder's zero, had points just short of the slit read as outside it:
  sections were cut short, a ray missed a wall and read material as air,
  and holes were attached a period away from the face they cut. Points
  are now placed on the side of the seam the face's trim is on, and a
  section is kept wherever it may pass through the surfaces' extents.
- **A hole's rim passing close to a loose edge was dropped.** Splitting a
  face, a loop was taken for part of the face's boundary, not a hole in
  it, whenever any point of it came within the face's weld of the
  boundary, and the weld on a face with loose edges is a few hundredths.
  A loop is part of the boundary only when it meets it at a node.
- **Circles and ellipses in all but one plane missed their crossing.**
  The closed form for conics in different planes looked for crossings
  only where one crosses the other's plane. Planes a hair apart, a hole's
  rim fitted to one facet group and a section through another, let the
  curves pass within the gap elsewhere, and the crossing was missed. Such
  pairs are measured by sampling.
- **A hole's rim cut a few microns from its own start stayed open.** A
  rim is cut where it crosses the hole's seam, and when that falls a few
  microns from where the rim's curve starts, the sliver between collapses
  in the chart of a face with loose edges. It did not collapse in space,
  so the rest of the rim ended on two vertices. A strand a face's weld
  collapses now makes its ends one junction.
- **Holes all but parallel to a bore or wall of a converted part failed,
  and slowly.** Two drums whose axes lean a ten-thousandth apart were
  marched, and their fitted crossings cost seconds each to intersect; a
  plane leaning a hundred-thousandth off a drum's axis met it in an
  ellipse tens of metres long, too coarse a ruler for the drum's height.
  Both are now solved in the drum's cross-sections along its height and
  kept as a line, or a cubic through the stations, with their departure
  stated as the section's tolerance. A hole that took 17 s and failed
  cuts in under a second, and a scan of a hundred random holes through a
  converted printer's part fails in one, from four.
- **A mesh with a flat sliver at a T-junction converted open.** A
  triangle with three distinct corners and no area was dropped, and the
  long side it sealed and its two short sides were each left used once.
  The triangle across the long side is now split at the sliver's middle
  corner, and the solid closes.
- **A hole near a ball's rim failed.** A drum passing through a ball
  close to its rim, all but grazing its far side, was marched, and the
  trace wandered along the long thin loops into a curve kilometres long.
  Where every line along the drum meets the ball twice, each loop is a
  function of the angle round the drum; the loops are now sampled exactly
  so, fitted with their images on both surfaces, and checked midway.

## [0.3.2] - 2026-09-24

A patch release with no API change; `cargo semver-checks` finds nothing
against 0.3.1. `project_edge_onto_plane` is new. Meshes convert far more
completely: coarse fillets, curved faces meeting along any curve, balls
and rings with holes, and a real printer's part all come back exact and
valid. Mass properties integrate most faces exactly whatever their trim,
and booleans cut cleanly along a closed surface's seam and through a
countersink whose bore is its narrow end.

### Added

- **`project_edge_onto_plane`** projects an edge orthogonally onto a plane
  as an exact curve in the plane's coordinates, for a sketch to constrain
  against. A line comes back a line or a point; a circle or ellipse a
  circle, an ellipse or a segment, its arc range counter-clockwise; a
  B-spline the spline on its projected control points. Any other curve
  is fitted and the fit's error returned.

### Fixed

- **The rim of a bore cut with a prism would not fillet or chamfer.** A
  circle face pushed down through a block, against its own axis, swept a
  general extrusion rather than the cylinder it is, and the blend refused
  its rim. It now sweeps a cylinder whichever way along its axis it
  travels, and the rim rounds with one torus, or bevels with one cone,
  exactly.

- **A converted printer's part came out invalid, a chamfer and fillet
  corners faceted.** Edge tolerances were recorded equal to the gaps they
  were measured from, and the checker, measuring again, found them a
  rounding over; they are recorded a millionth wider. A band round its
  axis cut open along a slit (a few triangles off its surface, left by
  the exporter) is built as a patch running round to the slit instead of
  falling to facets. And where two faces meet all but tangentially, a
  fillet running into a corner the mesh leaves free-form, the boundary is
  a curve threaded through the chain's own vertices, its tolerance how
  far it strays, instead of faceting the fillet whole.

- **Balls cut by two planes, and bores crossing at a coarse mesh, came
  back wrong or faceted.** A ball's axis could put both poles inside the
  cut-away part, where no loop goes round them, and the face covered the
  wrong side; both poles now stand on the ball's own surface. A flat face
  with two edges was taken for a sliver strip whenever it had two, and is
  now only when both are straight. An edge with an exact curve gets its
  image on a curved face interpolated at the edge's own parameters, close
  enough for its volume to be integrated exactly. Where a coarse mesh
  leaves skewed triangles on a surface next to a hole, or a few slivers
  together round a pole, they are taken into the surface when their
  middles sag no more than its curvature explains.

- **A countersink whose bore was exactly its narrow end refused to cut.**
  The cone's rim lay in the bore's wall, and the rim edge and the section
  along it were split at different points: the section at its middle, the
  rim where another section crossed it. Once every section is paved, a
  piece running along an edge is cut wherever the edge is split, and the
  edge wherever the piece ends.

- **A ball or a ring with holes in it stayed faceted.** Only a ball cut
  by parallel planes, or a ring cut square to its tube, was rebuilt. A
  ball or ring whose boundary only makes holes in it (a boss fused on a
  ball, a ball bored across twice, a ring pierced through its tube) now
  comes back as the whole surface with the holes as inner wires: its axis
  is turned so its poles stand farthest from every hole, and its seam into
  the widest angle the holes leave free.

- **Curved faces meeting along a curve that is no circle stayed
  faceted.** A bar drilled across, or a pipe with a branch, came back as
  hundreds of facets. Where two recognized surfaces meet along any curve,
  the meeting is solved onto both and interpolated once, its image on
  each surface made at the same parameters. A face round its axis between
  two rims of any shape, holed or not, is built with a seam of its own
  placed clear of its holes. Both parts convert exactly at a fine mesh and
  a coarse one.
- **A sliver face left a hole in the mesh of its solid.** A face narrower
  than a millionth of its length was drawn as nothing, and its
  neighbours' meshes did not meet across it. It is drawn as a fan across
  its own boundary points, which its neighbours share.

- **Coarse meshes lost their fillets.** At a printer's usual export
  settings a small fillet is three or four facets across, and it came
  back as flat faces. The first fit now tries a dozen vertices when they
  already turn through a bend, a fit needs half as many samples again as
  its surface has parameters, a failed sample that spanned a bend retires
  only its seed, and a facet spanning two rows of a tight fillet joins it
  when the surface accounts for its lean. A fully filleted box meshed at
  the default deflection comes back as its 26 faces.

- **A cut along a cylinder's or torus's seam refused to close.** A box
  face lying on the seam split the plane along the seam edge in pieces
  the curved face did not share. An open section is no longer mistaken
  for a closed loop and split at its middle, and a plane through a
  torus's axis meets it in its two exact tube circles instead of fitted
  ones.
- **Concave faces sent volumes to the mesh.** The check that faces agree
  about which way is out misread faces such as a three-quarter disc, and
  faces cut along their seam. It now reads each face's side from its
  walked boundary, so these shapes are integrated exactly.

- **Mass properties were a part in a hundred off on elliptic walls at the
  default deflection.** Only faces bounded by a chart rectangle or a full
  circle on an analytic surface were integrated exactly; anything else
  sent the shape to its mesh. A face on an analytic surface, or on a line
  or conic swept or revolved, is now integrated round its own trim, so an
  elliptic pad measures to rounding whatever deflection is asked for.
  Shapes with spline faces are still measured on their mesh, and say so
  in `deflection`. The default deflection no longer claims a part in a
  thousand.

## [0.3.1] - 2026-09-23

A patch release with no API change; `cargo semver-checks` finds nothing
against 0.3.0. Converted meshes rebuild whole spheres and tori, sphere caps
and zones, and bent tubes exactly, and drills through converted parts are
faster and close where they refused. The docs are rewritten shorter and
plainer.

### Added

- **`solid_from_mesh` rebuilds spheres and tori that close on themselves.**
  A whole sphere or torus comes back as one face. A sphere cut by one
  plane comes back as a cap, cut by two parallel planes as a zone, and a
  torus cut across its tube as a bent tube, each bounded by full circles.
  Before, these stayed faceted.
  The thin triangles a mesh often has round a pole, whose normals lean
  far from the surface's, join the region that surrounds them.

### Fixed

- **A drill through a finely faceted part spent minutes pairing
  sections.** Every section on the drill's wall was tried for crossings
  against every other section on it, a curve-curve intersection each, and
  a wall met by a few thousand facets of a converted mesh carries a few
  thousand sections: seven minutes on one part. A crossing between two
  sections that share a face counts only inside the faces they do not
  share, so only sections whose other faces' bounds meet are tried: the
  same part pairs its sections in eighteen seconds.
- **A drill through a face of many holes spent seconds on that one
  face.** A section was tried against every boundary edge of its faces,
  and a flat face bounded by a few thousand edges (a panel with four
  hundred vents) took them all; and choosing a split piece's nine probe
  points compared each of hundreds of thousands of candidates with every
  column already taken. Each boundary edge carries its own bound, and
  one outside the face the section must also lie in is skipped; the probe
  search stops at nine columns. A drill through such a part takes under a
  second where it took twenty.
- **A drill through near-coplanar facets refused to close.** A converted
  mesh merges facets that lie within its coplanar distance of one plane
  into one planar face, whose boundary then stands off that plane by up
  to that distance. A drill's section, lying in the plane, met such a
  boundary only within its own reach and missed it, so a whole face's
  section was dropped; neighbouring sections ended on their shared edge
  microns apart, farther than the drill wall's weld reached; and a
  section's end and the edge's split point were built as two vertices.
  On a plane a section meets an edge within the edge's radius, the weld
  on a face reaches as far as the edges its sections end on, and every
  crossing of a tolerant edge is one junction for both: the part that
  refused cuts to a valid solid.

## [0.3.0] - 2026-09-23

A minor release, for one change a caller can see: `FixReport` gains a
field, and a struct literal that built one no longer compiles. That is the only
break `cargo semver-checks` finds against 0.2.1. Meshes come in: 3MF
packages read, deflated, multi-part and ZIP64 alike, and a triangle mesh
becomes a solid whose planar regions merge into faces and whose curved
regions are rebuilt on the cylinders, cones, spheres and tori they lie on,
a scope addition `docs/SCOPE.md` argues. A boolean on such a solid of
thousands of faces takes a second where it took minutes. Imports hold the
tolerance rules they state, IGES round-trips real parts, and a face drawn
out to its whole plane, or the long way round an ellipse, meshes on
itself.

### Changed

- **`FixReport` gains `tolerances_widened`.** A public field on a struct
  callers can build: a literal that built one names the new field too.

### Added

- **`solid_from_mesh`**, a B-rep from a triangle mesh, with
  `MeshSolidOptions`, `MeshSolid` and `MeshSolidReport`. The topology is
  built from the mesh's own connectivity once its repeated vertices are
  welded, with no sewing search, so a 200 000-triangle mesh converts in
  about a second. Coplanar triangles merge into planar faces bounded by
  their outer loops and holes, and collinear boundary runs into single
  edges: an STL cube becomes six faces and twelve edges. Windings are
  made to agree and face outward, a closed piece inside another becomes a
  void, and a mesh that does not close comes back as open shells, with
  holes, non-manifold edges, dropped triangles and flipped windings
  counted in the report.
- **`solid_from_mesh` recognizes curved regions.** Triangles across which
  the surface turns smoothly grow into regions for as long as their
  vertices lie on one cylinder, cone, sphere or torus, verified at the
  coplanar distance, and each is rebuilt on that surface; its boundary
  with each neighbour is placed on the surface exactly, as a parallel
  circle or a ruling, and a band all the way round its axis gets a seam.
  A meshed bore comes back as a cylinder between two circles, and a
  rounded box as six planes, twelve cylinders and eight spheres. A region
  whose boundary is no such curve, or whose surface is round both ways,
  stays faceted and is counted (`curved_faces`, `curved_faceted`);
  `MeshSolidOptions::recognize` turns it off and `crease` sets the angle
  that counts as an edge. `recognize_points` exposes the recognition
  itself, with its measured deviation as the certificate. The mesh
  conversion is in scope by `docs/SCOPE.md`, which says why.
- **`read_3mf`**, reading a 3MF package into one placed mesh per build
  item, with `ThreeMfImport`, `ThreeMfObject` and `ObjectType`.
  Components are flattened, the production extension's multi-part
  packages followed (every slicer writes one object per part), the
  model's unit scaled to millimetres, a mirroring transform's triangles
  rewound, and a uniform object colour kept; per-triangle colours,
  textures, supports and unknown required extensions come back as
  warnings. Deflated entries are inflated by a decoder in the crate, and
  the archive is read through its central directory, so entries streamed
  with their sizes after the data read too. `read_package` reads deflated
  entries as well, checks every entry's checksum, and reads ZIP64
  archives (which streaming writers emit whatever a package's size, and
  which half the 3MF files downloaded from model sites are), refusing
  encrypted entries and archives spanning several disks by name.
- **`restore_containment`**, the pass that establishes the tolerance
  containment rule the checker enforces: every edge widened to at least
  the faces it bounds, every vertex to at least the edges it bounds,
  only ever growing, returning how many entities grew.

### Fixed

- **A sliver face was drawn as its whole plane.** A planar face narrower
  than a millionth of its length (a strip 163 mm long and 40 nm wide, as
  one CAD export leaves along a mirrored panel) has its two long sides
  merged into one line by the mesher, and no boundary ring survived; the
  mesher then took it for a face without wires and drew the plane's whole
  window, a hundred metres each way, so an imported printer's skirts
  rendered as bars across the scene. A face with wires whose rings all
  collapse encloses nothing, and meshes to nothing.
- **A boolean on a solid of many small faces took minutes.** Rebuilding
  the result compared every strand end with every vertex minted so far,
  and `sew` compared every edge with every other (projecting one's
  middle onto the other's curve for each pair), and every vertex and
  face likewise, so the cost grew with the square of the face count even
  when the tool touched a handful of faces. A thin drill through a
  converted mesh of 6 000 faces took two minutes. Vertices, junctions
  and edge ends are binned by position and compared only with those
  near them, in the order the full scan visited them, so the results are
  the same: the same drill takes under a second, and 12 000 faces two.
- **An edge on a very eccentric ellipse could run the long way round.**
  The STEP and IGES readers placed a vertex on an ellipse by its
  eccentric anomaly, which is exact only for a point on the curve. A
  vertex two microns off an ellipse 6.5 m by 1.8 mm read as nine
  millimetres along it, the end parameter fell before the start, and the
  edge took almost the whole ellipse: `triangulate_face` drew its faces
  out to ten metres. The readers refine the closed form to the nearest
  point on the curve.
- **A STEP read broke the tolerance containment rule, and `fix_shape`
  did not restore it.** The reader widens an edge to how far its pcurves
  sit from its curve, and left the edge's vertices at the confusion
  tolerance: `check` called every such vertex broken, hundreds on an
  ordinary part and thousands on an assembly, and `fix_shape` only ever
  tightened, so it never cleared them. The STEP and IGES readers run the
  containment pass on every body they build, and `fix_shape` runs it
  last, after the reduction.
- **A reshape reversed a wire's walk once too often.** A substitution
  rebuilds everything above the substituted entity, and read each
  occurrence's children with that occurrence's orientation composed in,
  then oriented the new node as the occurrence again. The occurrence being
  rebuilt came out right, but the new node it stored was inverted, and
  every other occurrence of it (a wire shared by a face reached reversed)
  walked its edges tail to head. Children are read through the forward
  occurrence and the result oriented once.
- **Collapsing a degenerate edge opened a gap.** `fix_shape` collapses an
  edge shorter than its vertices' tolerances onto one vertex, and the
  surviving vertex kept its own tolerance, so the neighbouring edges'
  curves stopped the collapsed length short of it. The survivor widens to
  reach every vertex it absorbs.
- **`write_iges` overflowed the eighty-column record, and its own reader
  refused the file.** A near-zero real went out in positional notation:
  a coefficient of 1.5e-51 spelt in sixty-nine characters where a record
  holds sixty-four. Every real takes the shorter of its positional and
  exponent spellings, and a parameter longer than a whole record fills
  every record it crosses, so the reader rejoins it with nothing inserted.
- **`read_iges` refused a solid whose curves missed their vertices by
  nanometres.** IGES states no tolerances, so every vertex was built at the
  confusion tolerance, and a writer's last digit (half a nanometre)
  refused the whole solid; a vertex a hair past a bounded curve's end
  also asked for a parameter outside the curve. A vertex's tolerance now
  grows to state its curves' miss, up to the millimetre past which a
  boundary is not that curve's at all, and the range is held to the
  curve's domain. Every solid in the exchange corpus survives the
  kernel's own IGES write and read, at the volume it went out with.

## [0.2.1] - 2026-09-23

A patch release: no public signature changed in any crate, checked against
0.2.0 by `cargo semver-checks`. It adds, and it fixes output that was wrong.
Fillets reach every host the march can seat on, fitted patches included,
and the corner tool rounds any convex planar vertex. Edges asked together
meet: `chamfer_edges` mitres, and `fillet_edges` closes a shared corner
with the rolling ball's patch rather than leaving the bands' caps. This is a
visible change for any caller that rounds a box corner in one call. A
viewer meshing face by face can now agree with the whole-shape mesh along
every edge, and STEP files from plainer writers and mesh converters read.

### Added

- **The marched blend takes fitted hosts as far as the melt.** A
  B-spline patch, or a swept or revolved surface, is a host the march
  had refused by name. It marches now: the chart inverts by projection
  warm-started from the last station, a patch that meets itself at its
  seam wraps as a periodic one does, the patch is continued past its
  face by a few radii so the ball can run out through a wall, and a
  seat the boolean split into arcs at the hosts' seams is closed back
  into one loop through the neighbours the two hosts share, arc by arc
  where each continues the last tangentially, and marched as the loop it
  is. The wedge is built on the host's own patch, widened in place, its
  rail loop slid into the host's own window, and it melts: a converted
  box's edge, the crease round a converted post, a converted cone's rim
  and a converted drum's rim all round to valid solids. A melt the
  boolean still cannot resolve refuses by name. Under it, `to_nurbs` now carries the
  old vertices' and edges' recorded tolerances into the converted solid,
  which it had dropped, and the boolean inverts a probe on a patch by
  projection and treats two identical patches as one chart.
- **Chamfers bevel a chain of edges as one operation, mitred.**
  `chamfer_edges` and `chamfer_edges_with` (a `Chamfer` per edge, in any
  of the three spellings) build every wedge on the solid as it stands
  before the call and then apply them all, so where two bevels meet at a
  vertex each still reaches the corner and the bevel planes meet along
  their own line. Beveling a box's four top edges one call at a time left
  a small tetrahedron and two extra faces at every corner; as one
  operation it is ten faces and the exact volume. Three bevels at a
  convex vertex meet at one point.
- **The corner tool rounds a vertex no single ball touches.** A
  rectangular pyramid's apex, an irregular pentagonal one, or any convex
  planar vertex whose faces share no tangent ball was refused by name.
  The ball's centre may sit anywhere a radius in from every host plane,
  and at such a vertex that region's tip is a few points joined by short
  ridges; the rounded corner is the exact envelope of the rolling ball, a
  sphere at each tip vertex and a cylinder along each ridge, each cut
  with its own compartment (the sphere by the one-ball tool on its three
  planes, the ridge by the flush fillet of a virtual crease), the
  compartments meeting cap to cap on the planes square to the ridges. A
  ridge seven microns long stays the sliver of cylinder it is. The
  sphere's pole stands along a ridge where there is one, so the ridge
  fillet's cap meets it on a circle the chart images exactly. The flush
  fillets follow at a sharp rectangular apex; where a sphere clears a
  fourth plane by a few hundredths of a millimetre the envelope keeps a
  sliver of that plane, and the flush fillet meeting the sliver still
  dies in the cut.
- **IGES reads more of the 1980s.** A conic arc (104) of any axis-aligned
  kind reads as its own curve: hyperbola and parabola alongside the
  ellipse, the arc's ends read off the file's points whichever way round
  they are listed. A ruled surface (118) reads as the degree-one patch
  between the two curves' exact spline forms, raised to one degree and
  refined to one knot vector, the second curve walked backward where the
  direction flag says. A constant offset curve (130) and an offset
  surface (140) read as this vocabulary's offsets, the file's direction
  conventions carried by the sign. Parametric spline surfaces (114) and
  rotated conics stay refused by name.
- **Gauss–Kronrod quadrature, and an integral over a rectangle of
  parameters.** The adaptive integral estimated its error by halving and
  comparing, thirty evaluations an interval; it now takes the seven-point
  Gauss and fifteen-point Kronrod pair, fifteen evaluations shared, the
  gap between the two the error. `integrate_2d` integrates over a chart
  rectangle the same way, in tensor form, halving a cell along whichever
  direction the pair finds rougher, so a crease running across the chart
  costs a line of cells rather than a field of them.
- **A stated-tolerance fit for what has no exact spline, and a degree
  limit.** `Curve::to_bspline` still refuses a helix, an offset curve or a
  surface curve rather than approximate silently; `Curve::fitted_bspline_over`
  is the explicit ask, fitting the curve at its own parameters to a named
  tolerance, measured between the samples as well as at them, the error
  reported as found. `BSplineCurve::restricted_to_degree` and
  `BSplineSurface::restricted_to_degree` bring a curve or a patch under a
  degree limit the same way, which is what an exchange format with limits
  needs written.
- **B-spline curves and patches extend.** A B-spline used to end where
  its net ended: an analytic surface's window widens over its carrier,
  but a patch had nothing past its last row. `BSplineCurve::extended`
  and `BSplineSurface::extended` continue a curve past either end, or a
  patch past any of its four sides, by about a length in space: the
  polynomial continuation of the curve's own end derivatives to the
  order asked, raised to the degree and joined on. A rational arc
  continued at order two stays on its circle, and a cylinder patch
  continued round its circle or along its axis stays on its cylinder.
  `widened_to_hold` continues a patch, side by side, as far as a point
  stands off the side its projection clamped to, which is what the
  neighbour-extension steps need on spline faces.
- **Marched blends on cone, sphere and torus hosts.** A seat with one
  host a cone, a sphere or a torus (a bore through a cone's wall, a hole
  drilled through a ring's tube, a ball drilled off its centre) was
  refused by name; the chart inversions those surfaces have in closed
  form are carried now, and a band whose rings start on different columns
  gets a fitted connector where only a cylinder had an exact one. A full
  circular rim whose hosts are not a planar cap and its coaxial wall (a
  bore straight down a ball's axis) takes the march too, where the
  revolved blend had refused it by name. Under it, two fixes any seat
  could hit: a looping seat whose join is a
  corner (the two arcs of a boolean's seam joined end to end) steers
  the march by a smooth refit of itself, since a march cannot cross a
  corner in the guide's derivatives; and the band's first station is
  re-solved exactly on the apex column, where a station a fraction of a
  stride off the host's own seam left a sliver the boolean could not
  hold.
- **Draft on any face.** A wall on a raw fitted patch, or a wall of
  revolution about a neutral plane tilted against its axis, used to be
  refused. Both draft now the way a mould-maker drafts: the face becomes
  the ruled surface through its crossing with the neutral plane, each
  ruling the pull direction turned by the draft angle about the
  crossing's tangent. The planar and revolved drafts are the closed
  forms of this construction. The crossing is read off the face's own mesh
  and corrected onto the surface, so it lies inside the face whatever
  the chart does; the support is the ruled surface itself, degree one
  along the ruling, between a cubic fitted through the crossing and one
  fitted through the rulings' tips at the same parameters. The rebuild
  underneath learned what a fitted support needs: a band on one is
  assembled wire by wire, its rings marched against the caps and
  re-seamed at the vertex the seam starts from, its seam the support's
  own iso-curve (`BSplineSurface::iso_u_curve`, new), and a seam's end
  vertex the crossing of that column with the other seat rather than
  the nearest point two seats agree on.
- **A skinned loft caps a non-planar end section.** A wavy rim, on which
  no plane can stand, is capped by a patch skinned from the rim down to a
  point inside it, sharing the wall's ring edge; the loft used to refuse
  it. A planar end keeps its plane.

### Fixed

- **A prism swept from a clockwise profile was inside out.** A closed
  wire on a plane bounds one region however it is walked, but the prism
  read each wall's side off its edge's direction, so a square walked
  clockwise about the travel swept its walls facing into the material
  while the caps faced out, and the volume was refused as wound inward.
  Each ring's walls now follow the ring's winding about the travel, the
  outer ring turning positively and every hole the other way, whichever
  way the caller walked them.
- **Fillets asked together at a corner close it with the ball's patch.**
  Three fillets at a box corner through `fillet_edges` left the bands'
  flat caps standing where a ball rolls round the corner. Where three or
  more edges of one call meet at a vertex, the corner tool now rounds
  the vertex first and the bands stop flush against its patch: an
  octant of a sphere at a box corner, the envelope of spheres and
  cylinders at an apex no single ball touches. A corner the tool does
  not speak keeps the caps. Under it, the corner tool read a vertex's
  point where its node was built rather than where it stands, and at a
  prism's far end it rounded the near corner instead.
- **A band meeting two bands at an apex left its triple point as two vertices.** Three
  descriptions of one point arrived a tenth of a micron apart in turn,
  each within the last's reach and none within the first's, and the
  junctions the arrangement adds after the first merge were never merged
  with each other; they are, an end welds into a vertex within the
  vertex's own tolerance and the weld, and a crossing at an edge's end is
  the end on the edge as well as on the section. What is still owed is
  named in the ledger: a band's section touching an edge describes the
  touch point twice, a tenth of a micron apart, and the third flush fillet
  at a sharp apex still dies on it.
- **STEP faces on a surface of linear extrusion were skipped, and the
  body drew open.** A curve swept along a vector is how some writers
  spell a drum's wall, and the reader skipped every such face by name,
  so a part with a slot's rounded ends on them showed its inside through
  the gaps. A circle swept along its own axis now reads as the cylinder
  it is, a line swept as the plane it is, and anything else as the swept
  surface itself over a window the face's own edges size, its trims
  fitted by projection.
- **An open edge whose vertices stand on its curve short of the ends.**
  A file that writes the whole spline and lets the vertices say where
  the edge stops, millimetres in, had the edge held to the whole curve,
  overshooting its neighbours, and the faces it bounded drew as nothing.
  The window between the vertices' own feet on the curve is taken, and
  the report tallies it as `vertex-window`.
- **A marched section between two patches wandered, and a loop cut at
  a seam stayed open.** The boolean fits a marched section jointly with
  its two chart images, and parameterised the samples centripetally: a
  walk that starts on a window's rim halves its step to a micron and
  grows it back only twofold a point, so most of its samples crowd one
  end, and a cubic through them could not be the straight line they lay
  on: a section came back twelve hundred millimetres off it, and the
  boolean refused the wedge for running beyond a turn. A joint fit that
  misses its target centripetally is fitted again by chord length, and
  the closer of the two stands. A loop
  walked round a converted drum reached the seam from both sides and was
  flagged as having left the domain; its ends coincide and it is closed,
  its chart image unwrapped across the seam as across a period.
- **A fused post converted to patches drew open, and the conversion
  failed the checker.** The second half of a rim split by a fuse was
  written as the chart segment between its projected endpoints, the seam
  vertex projecting to either column, and it drew the front half
  mirrored; the wall's ring lost its far side. A seam endpoint's column
  is decided by the edge's own interior. The fit's slop widened each
  converted edge past the vertices bounding it; the vertices widen with
  the edge.
- **A fitted pcurve hooked off its curve between the fitter's samples.**
  The reader fits an edge's pcurve by projecting samples spaced evenly in
  the curve's parameter and measuring the fit at those samples alone. A
  curve is free to run far faster at one end than the other (a fan
  blade's root meets its hub through most of its bend inside the first
  sample interval), and there the cubic hooked four tenths of a
  millimetre past the curve while every sample sat within a hundredth:
  the ring crossed the face's neighbouring ring in the chart and the hub
  drew with overlapping triangles and open edges at every blade. The fit
  is now measured between its samples as well, through the surface, and
  where it leaves the curve by more than a micron it is asked again at
  twice the samples, a few rounds at most, the closest round kept; the
  reported error and `fit-short` count what was found between the
  samples too.
- **A mesh assembled face by face cracked along every narrow face's
  edges.** A face narrower than a few chords draws its edges finer than
  asked, as does one whose boundary crosses itself at the chord, and
  `triangulate` tells the faces across those edges to draw them the
  same. `triangulate_face` could not be told, so a viewer keeping one
  mesh per face (and the kernel's own stored tessellation, built the
  same way) drew the two sides of such an edge to different points: a
  seam of cracks round every fillet and along every thin plate.
  `edge_chords_for` answers the agreement once per shape and
  `triangulate_face_with` draws a face to it; `tessellate` stores its
  faces and its edge polylines by the same answer.
- **A bore crossed by a hole at its wall drew with a ring missing.** The
  hole's two rims wind the bore's chart between them, and pairing them
  along the first rim's last column could run into the hole's own
  opening, or leave a rim's tail on the far side of the seam from the
  ring it was folded into. The pairing column is chosen in the widest
  gap free of the other rings, either rim is cut whichever way it runs,
  and the rings the band passes are slid into the same stretch of the
  period; a cross hole through a bore wall draws closed.
- **STEP files spelt in older or plainer vocabulary.** A face written as
  `FACE_SURFACE` rather than `ADVANCED_FACE`, a shape carried by a
  `SHAPE_REPRESENTATION` or a manifold surface or faceted representation
  rather than an advanced B-rep one, and a product whose formation is
  written with its source read as they were meant; the source-spelt
  formation's product is taken from its third slot whether it stands
  alone or inside a complex instance, and a product whose name slot is
  blank (a mesh converter fills the id and leaves the name) reads by
  its id rather than its entity number. An edge whose
  vertices stand at descending parameters on its curve whatever its
  sense flag says (every edge written forward and the line left to run
  the other way) is reversed to run with the curve, tallied as
  `edge-against-curve`.
- **Removing a tangent chain of blends.** A stadium's top rim rounded
  in one call, its four bands named for removal together, failed by
  name: at a tangent junction neither band's crease pierces the other's
  wall (the straight crease grazes the round wall and the round crease
  grazes the flat one), so no corner was placed and the straight crease
  was taken for a wrapping band. The junction is now where two creases
  touch, placed exactly by the cross-section edge the two bands share:
  the foot of that edge on either crease. And a band standing across a
  circular crease's seam is read as one run about its own centre rather
  than its complement, so an end wall grows back its outer half and not
  its inner.
- **A torus band between two parallels drew as two whole tori laid over
  each other.** Each wire of such a face is a single edge winding the
  chart once. Walked one wire at a time, each rim closed on its own
  translate a tube-period over, which is the whole torus cut along that
  rim, and two of those cancel where they overlap: the band drew to
  more area than its torus has. Wound rims are paired into one ring
  before any rim is closed alone, on either periodic axis, so a ball's
  belt between two latitude wires holds the same way. The parity check
  that guards a boundary had been flagging exactly these faces and
  drawing them again, eight times finer, every time; the benchmark part
  that carries two such bands tessellates in a twentieth of the time.
- **A face narrower than the chord shaded as a saw of fins.** Drawn at
  the caller's chord its boundary sags between points by more than the
  face is wide (a thread flank a tenth of a millimetre wide at a chord
  three times that), and every triangle across the width stands off the
  surface by the sag. A face narrower than four chords has its edges
  drawn to a quarter of its width, and the faces across those edges draw
  them the same. The whole-shape pass reads every face's width off its
  rings before drawing anything, so the big faces round a thousand small
  fillets are drawn once; on a real assembly the area standing more than
  45° off the surface falls from 540 to 80 mm², for 7% more triangles
  and a quarter more time.
- **A vertex at a cone's apex, a sphere's pole or a patch's collapsed
  corner had no normal.** The surface has none there, and the mesh carried
  a zero, which shaded the tip black and dragged every normal it was
  welded with towards nothing. The vertex now carries the limit normal
  along its own column, found a step inside the domain: the apex of a
  cone seen up a ruling has that ruling's normal. On a real assembly,
  1,199 such vertices; none now.

## [0.2.0] - 2026-09-22

A minor bump before 1.0, which is to say it breaks one thing: the default
angular deflection is half a radian. Everything else is additions and
fixes, most of them to meshing, measured on a community printer assembly
of 1,563 solids: every one of them meshes watertight at the default
deflection, with no triangle standing off its surface, in a third of the
time 0.1.0 took.

### Added

- **`heal::fix_shape`**, one pass over a shape nobody promised was
  well-formed: diagnose; put a wire's edges end to end where an order
  exists; collapse an edge shorter than its own vertices' tolerances;
  give a pcurve-less edge the trim projection can honestly fit; sew a
  compound of loose faces or an open shell; tighten tolerances; diagnose
  again. `Fixed` carries the result, its history, and a `FixReport` of
  what was done and what the checker still sees. Small faces and small
  solids are not removed (that is defeaturing), and the report says so
  by leaving them in `after`.

- **`Surface::curvature_at`**: the principal curvatures and directions at
  a point, signed against the surface's own normal, with the mean and
  Gaussian curvatures they combine to and whether the point is an
  umbilic: what curvature display and zebra analysis ask, from the two
  fundamental forms of the jet the trait already carried.
- **Ruled lofts take skew walls.** A wall between two segments that are
  not coplanar (a square lofted to the same square turned an eighth of a
  turn) is the bilinear patch through its four corners, the ruled
  surface between them and exact, where it was refused by name and
  routed to the skinned loft.

### Fixed

- **Meshing a real assembly closed.** 1,499 of 1,563 solids meshed
  watertight before; 1,562 do now, and none refuse. Six causes, each with a
  fixture cut from the file that showed it:
  - `Triangulation::is_closed` counted an edge's uses and wanted exactly
    two. It now requires every edge to be crossed as often each way, which
    is what the divergence theorem needs: it catches a face wound inside
    out (which counting passed) and accepts four triangles round one edge
    (which counting refused).
  - A face narrower than the chord error its boundary is drawn with had a
    boundary that crossed itself; the triangulator returned fragments. The
    boundary is redrawn finer, keyed to the edge so both faces sharing it
    agree.
  - A run along a cone's apex row was dropped when its rulings stood
    exactly a quarter turn apart. Whether a gap is an apex run is decided
    by its width and by whether its chart midpoint lifts to the shared
    vertex, not by a fraction of the period.
  - A chart direction was called degenerate by ratio to the other
    direction, so a densely parameterised `v` made every `u` look weak and
    two half-millimetre edges fitted as a point. It is judged by what
    crossing the whole span moves, in length.
  - Ring folding across a chart's join engaged on periodicity only; a
    B-spline tube that closes without repeating never folded. Closure is
    the test now, folded rings are slid back into the domain, and
    `Surface` evaluation wraps a parameter past a closed join instead of
    refusing it.
  - An inner loop thinner than a micron (arcs out, fitted splines back)
    is a slit, not a hole, and is no longer handed to the triangulator.
- **A long bore drew square between its cross holes.** A cylinder never
  sags along its axis, so its grid got one interior row, and the Delaunay
  triangulation bridged from each rim to it with triangles a quarter turn
  wide, under the three chords the repair pass fires at, so they stayed.
  Grid cells are held to a bounded aspect: rows close enough, measured
  through the surface, that no triangle reaches across more than a few
  columns. Fewer triangles on the bore than before, not more.
- **Faces shaded as quilts of creases.** Delaunay in the chart is not
  Delaunay on the surface when the chart's units differ by axis (a
  cylinder's `u` in radians against its `v` in millimetres, a fitted
  strip's `u` over a fiftieth of a unit against a `v` over one), and the
  slivers it makes across the narrow way lift folded, flat across a bend
  the surface takes in between, their normals pointing where none of
  their vertices' do. Three things, measured on a real assembly by the
  area of triangles standing more than 45° off their vertices' normals,
  61,000 mm² before and 540 after:
  - the triangulation runs in the chart scaled to the surface's own
    metric, the mean tangent length each way, so Delaunay sees distances
    as space does;
  - an interior grid point keeps a third of a cell clear of the boundary,
    where a point hugging a boundary chord makes a sliver that stands off
    the surface as a fin;
  - the angular deflection is not asked of a segment both shorter than
    the chord tolerance and a sixteenth of its edge, such as a fitted edge's
    end hook a few microns long, which bisection chased down to the
    resolution of the parameter and handed the face a fan of hairs.
  The same assembly meshes in a third of the time with a fifth fewer
  triangles, the slivers the repair pass used to chase now never made.
- **A bore's inside was coarser than its rims.** Edges are discretized
  to the chord *and* the angular deflection; the interior grid was held
  to the chord alone. A viewer scaling its chord to a body's size gave a
  long extrusion's bore thirteen-sided rims and a seven-sided inside. Grid
  cells are held to the normal's turn as well now, so a surface's inside
  is as round as its boundary.
- **A face whose chart sat far from its origin drew as two triangles.**
  The scale a degenerate chart triangle was measured against was taken
  from the ring's coordinates rather than its span, so a cylinder whose
  axis point the file placed half a metre away had every honest cell
  called a hair. Fixture cut from the file that showed it.
- **At a coarser angular deflection, three more ways a face came apart**,
  each with a fixture cut from the file that showed it:
  - A boundary drawn coarsely enough to cross itself and enclose *nothing*
    (an annulus narrower than its rims' sag) was refused, where one that
    enclosed fragments was drawn again with finer edges. Empty is short
    too; a face is refused only when finer edges still enclose nothing.
  - The sag repair inserted a sliver's centre that was its own apex to the
    last bits, round after round, each a hair on the last; the degenerate
    filter dropped the hairs and left a hole. A repair point that lands on
    a vertex already there is not inserted.
  - A grid point that fell exactly on a boundary segment running
    diagonally across the chart was inserted and split that constraint on
    one face alone: a T-junction against the face across the edge. Grid
    points on the boundary stay out.
- **A crossing buried under interior points, a spike, and slop at the
  closing vertex.** Whether a face's boundary crossed itself was read
  from the count of triangles against boundary points, which over a face
  with hundreds of interior points misses a crossing that costs a
  handful; a slotted cylinder wall drew with two holes more than it has.
  The boundary is triangulated on its own first, where the count is
  exact and told by constraint parity rather than by where a hair's
  centre rounds to. What that exposed: a ring that runs out along an edge
  and straight back is a spike that bounds nothing and comes off; and a
  ring whose last point is its first a file's slop away (the two edges'
  own ends of the vertex they share) closed with a fold over its first
  segment, a crossing a fraction of a micron deep. Merged, every solid in
  the real assembly meshes watertight at the default deflection, 1,563 of
  1,563; at half a radian, 1,562.
- **A closed edge's seam is moved to its vertex.** A fitted loop written
  with its start wherever the fit began, the edge's one vertex millimetres
  along it, was held to the curve's own seam, and the vertex's tolerance
  widened to the miss: 2.18 mm in a real assembly, a reach a solid's
  border weld then used. The reader moves the seam to the vertex: the same
  curve, begun where the edge does. `BSplineCurve::reseamed_at` and
  `bspline::join` are new.
- **`triangulate_face` drew every face twice** since the sliver refinement
  landed. It draws once, and again only where the first came up short.
  Face by face, the assembly above meshes in 9.4 s where it took 18.6.

### Changed

- **The default angular deflection is half a radian**, twenty-eight
  degrees, a circle in thirteen segments (what B-rep kernels have long
  defaulted to), where it was 0.2, eleven degrees and thirty-two. With
  the interior of a face now held to the angular deflection as its
  edges are, the old default cost four times the triangles on every
  cylinder; a real assembly meshes to 4.0M triangles at the new default
  where the old gave 19.5M. Ask for `angular: 0.2` to have what the old
  default drew.
- `BSplineSurface` settles whether its net closes at construction, so
  evaluation past a closed join costs what evaluation inside costs.
- CI verifies the declared `rust-version` on every push, reading it from
  the manifest so the two cannot drift.

## [0.1.0] - 2026-09-18

First public release.

### Added

- **Geometry.** Parametric curves and surfaces in two and three dimensions
  (lines, conics, B-splines rational and not, extrusions, revolutions, offsets
  and trims) behind adaptor traits, on a B-spline substrate with knot
  insertion, splitting, degree elevation and Bézier decomposition.
- **Topology.** One shared B-rep model: geometry and topology in arenas, a
  shape a cheap handle into them, the same node placed, mirrored or instanced
  many times without being copied. Per-entity tolerances, location chains,
  orientation composed through the tree, and edges carrying a list of
  representations rather than one curve.
- **Construction and measurement.** Primitives, polyhedra, sewing, shape
  validity, mass properties exactly where a closed form exists and from the
  mesh otherwise, bounds, projection and classification.
- **Booleans.** A general fuse with the filters over it (union, difference,
  intersection, section) plus defeaturing over the same machinery.
- **Blends.** Constant and variable radius fillets and chamfers on single
  edges, tangent chains and full rims, closed forms where they exist and a
  marched rolling ball where they do not, and the corner where three blends
  meet.
- **Offsets and sweeps.** Offsetting, shelling, sweeping along spines open and
  closed, lofting through sections, and draft.
- **Tessellation.** Edge discretization to chord and angular tolerances, and
  constrained Delaunay triangulation per face in its own chart.
- **Healing.** Validity diagnosis, sewing, reanchoring periodic rings, and
  instructed fixes for trims a reader refused and boundaries sitting off the
  surface they bound.
- **Drawings.** Hidden line removal, sections and hatching.
- **Documents.** Assemblies and product structure, appearance, PMI, saved
  views, transactions and persistence.
- **Exchange.** STEP and IGES in both directions (STEP carrying assemblies,
  colours, semantic PMI and saved views) plus the native format, `.brep`,
  STL, DXF, glTF, OBJ, PLY, VRML and 3MF.

### Known restrictions

This release is measured against 97 named capabilities in
`docs/PARITY.md`, and the ledger is part of the build: 66 are covered, 16
carry a stated restriction, 9 diverge by design, 5 do not apply and 1 is
unreviewed. `tools/check.sh` fails if the audit and the code drift apart.
The `partial` rows say exactly what each one does not do; the largest are
blend hosts beyond planes and cylinders, the general N-edged setback vertex,
shape healing's wire reordering and small-feature removal, the medial axis
beyond convex polygons, and the IGES entities listed as refused by name.

Nothing here is a silent gap. A capability that is not implemented refuses
by name rather than returning an answer it cannot stand behind.

[Unreleased]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.11.5...HEAD
[0.11.5]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.11.4...v0.11.5
[0.11.4]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.11.3...v0.11.4
[0.11.3]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.11.2...v0.11.3
[0.11.2]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.11.1...v0.11.2
[0.11.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.11.0...v0.11.1
[0.11.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.16...v0.10.0
[0.9.16]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.15...v0.9.16
[0.9.15]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.14...v0.9.15
[0.9.14]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.13...v0.9.14
[0.9.13]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.12...v0.9.13
[0.9.12]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.11...v0.9.12
[0.9.11]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.10...v0.9.11
[0.9.10]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.9...v0.9.10
[0.9.9]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.8...v0.9.9
[0.9.8]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.7...v0.9.8
[0.9.7]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.6...v0.9.7
[0.9.6]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.5...v0.9.6
[0.9.5]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.4...v0.9.5
[0.9.4]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.3...v0.9.4
[0.9.3]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.2...v0.9.3
[0.9.2]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.6.1...v0.7.0
[0.6.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.4...v0.4.0
[0.3.4]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.3...v0.3.4
[0.3.3]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/gilbertorconde/ogeom-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/gilbertorconde/ogeom-rs/releases/tag/v0.1.0
