# Plan

What is left to build, and the rules the work follows.

## Where the remaining work is tracked

- **GitHub issues** hold bugs and requests, one per issue.
- **The capability ledger** (`docs/CAPABILITIES.md`, generated from
  `docs/capabilities.toml`) records the kernel's named capabilities. The restriction on each `partial` row is the
  list of what that capability does not do yet. `tools/check.sh` holds the
  ledger to its evidence.
- **This file** holds the open items that are neither: known limits of the
  code, refusals by name, design items and speed work.

Finished work is recorded in `CHANGELOG.md` and in the tests that pin it, not
here.

## Open

Within each group the items are ordered by risk: a possible wrong answer
first, then refusals, then cost.

### Booleans

**A section passing a pole within its tolerance, beyond its chord, is
refused.** A marched branch through a pole of either surface's chart, the
pole lying on the other surface within the trace's chord, is cut there
(`cut_at_poles` in `ogeom-intersect/src/march.rs`): each piece ends on the
pole, its image there on the line of the chart it leaves along, and the
arrangement cuts the face's pole row where the pieces end on it. A branch
passing beside a pole, farther than the chord, is fitted whole; where the
joint fit of curve and images misses its tolerance a hundredfold and
either chart has a pole, the curve is fitted in space alone and each image
read off it by inversion, measured lifted against the curve (`projected`
in `ogeom-intersect/src/approx.rs`), the closer fit standing. A face left
on its own chart (where `ball_chart.rs` cannot restate it: a face with one
new pole inside whose loop round it the seam meridian crosses at no
vertex, a face placed inside its solid, a face with an edge that has no
curve in space, a face whose trim holds no pole of its chart) therefore
takes a section through its pole or beside it. A branch fitted whole that
passes a pole within its stated tolerance (the chord and the fit's error)
is traced again, cut at the pole within that tolerance
(`branches_cut_at_poles`, `passes_a_pole`): fitted whole there it states a
tolerance about its own distance from the pole, since no image follows the
chart's half turn, and the pole is a point of the section as far as the fit
can say. Left on their own charts, a ball and a ball slab drilled with the
drill's wall through the poles or a ten-millionth of the radius beside
them, and a ball drilled with it a hundred-thousandth to a hundredth of the
radius off them either side, cut at the closed forms. Still refused: a
section whose image runs past the pole row by more than its tolerance ("a
section through a sphere chart's pole is not resolved"). A loop crossing
the sphere's seam within the window that reads a crossing as the seam's
end vertex, the pole, ends in the chart off the pole row by that window
over the radius and anywhere along it, the chart's free parameter
shrinking to nothing there; the end is carried onto the pole's line at its
foot, held to the line, and the pole is cut for it. A march round a pole
that starts beside it closes on the step that passes over its start (the
step onto the start, solved there, lands anywhere along the pole's line),
and a closed branch passing over its own start mid-way is held to its
first lap. All six placements at 2e-5 off cut at the closed forms. "A branch
passing beside a chart's pole fits nothing yet" stands where both fits
miss. Only spheres' poles are measured; a patch's collapsed side takes the
same path untested.

**Nearly coplanar faces cross beside a curve only where the rest resolves.**
Two planar faces leaning a hair apart cross at a small angle, and a strip
of one between the crossing and its edge, narrower than the confusion
distance's reach divided by the angle, stands on the other at every point.
On the attempt that settles such pieces from both sides (`settled` in
`ogeom-bool/src/lib.rs`), the edges bounding the strip are carried onto
the other face from where a plane-plane section line crosses them
(`graze_contacts`), so the other face's matching strip is split off and
the two are kept as one contact. A line edge is crossed first; where that
refuses, any edge in the face's plane is (`section_crossings`: a circle in
closed form, any other curve where its distance across the line changes
sign at 64 steps), and the stretch carried between two crossings is not
welded shut as a graze. A cube, a drum's disc and an arched block's arc
against a box whose wall leans 1e-5 to 1e-3 radians and crosses a hundredth of
a micron to a micron inside the edge fuse to the closed form within the
weld distance over the band (`boolean_contact.rs`). An arched block whose
arc is the same arc as a rational B-spline meets the box's wall in a
section that crosses the arch's top rim all but tangentially, within the
confusion distance of it over tenths of a millimetre, and the marcher
stops short of the rim there. Where every retry refuses, the boolean runs
again with that section split where the arc's height to the wall's plane
reaches the rim, in closed form (`plane_extrusion_at_rims` in
`ogeom-intersect/src/extrusion.rs`), an edge crossing a section beside a
plane held to its own distance off the plane over its lean out of it,
and an exact section's end touching an edge held to its own distance
off the edge. Over leans of ±1e-5, ±1e-4, 1e-6, 3e-4 and 1e-3 and
crossings 1e-6 to 1e-2 inside the crown, 183 of 189 operations answer
within 7.3e-7 of the closed form; 6 refuse, at a lean of 1e-5 with the
crossing 1e-5 or 3e-5 inside ("a piece lies on the other solid's
boundary" for the cut, fuse and common at 3e-5, "the kept pieces did not
close" at 1e-5). Traced on the cut at 3e-5: the end face's cap above the
crossing, 3e-5 high and 0.035 wide, lies inside the band where the two
planes stand within the confusion distance (1e-2 wide at this lean); its
interior probe falls on the crossing line, and the first attempt finds it
on the box's boundary with no partner. Settled from both sides it reads as
the aligned contact, and the kept pieces then do not close. The attempts
that do build a solid (at a wider confusion, and split at the rims) weld
the crown's junctions across the band, three vertices widened to 1.3e-2,
2.2e-2 and 3.9e-2, and are refused by the weld limit. Resolving it needs
the arrangement to hold the cap's own short edges and vertices apart
inside the band. At a lean of 1e-5 the faces the wall leaves on the
arch's wall, a sliver 5e-5 thick at most, are measured on their own
meshes, within 5.2e-9 of the closed form at a chord of 1e-3 as at 1e-4
(#683). The cut
at a lean of 1e-5 keeps the sliver above the wall, a retry at a wider
confusion dropping a shell only within the caller's weld distance (#674). Where
the strip is wider the ordinary crossing splits both faces exactly, and
pairing the band as same-domain would move the result by up to the weld
distance over it, so the band is not paired there.

**A section whose fit misses widens its junctions as far.** A marched
section keeps its fit's measured error as its tolerance, and the paving's
reach and the junctions along an edge grow with it, unbounded. Over the
workspace suite (heavy tier included) and the stress suite at seeds 1 to
4, a boolean's sections no longer state more than 2.1e-3, with a sphere
chart's pole met only where the chart is left unturned, which refuses.
Two still miss by more than a few chords. A skin crossing a mitre
plane ever more shallowly toward the corner
(`an_off_centre_section_down_an_l_mitres_between_skinned_runs` in
`crates/ogeom/tests/pipe_corners.rs`, the surfaces' normals parallel to
1e-14 at the corner) states 2.2e-4: within the confusion distance of both
surfaces the walk's correction can land a sample off the crossing, one by
2.9e-5 in the skin's chart, and the fit follows it. Holding each span
between stations to the curve at its middle (the condition solved halfway
within the chord of the cubic the two stations' positions and tangents
make, else the step halved) brings that section to 1.8e-5, but breaks the
walk beside branch points: a band crossing a cone through its node
(`a_band_crossing_a_cone_in_a_window_cut_by_its_node_has_no_section_along_the_cut`)
and two equal cylinders crossed square
(`large_equal_cylinders_crossed_square_come_back_as_four_halves`) lose their
ends at the node and the touch, the halfway solve there landing on the other
branch, the step halved to the floor and the walk creeping through. Judging
only halfway solves within a quarter step of the chord keeps the branch
points and loses the mitre (1.5e-4): the stations the check refuses there
are the ones whose halfway solve lands as far off. Telling a station landed
off a shallow crossing from the other branch at a touch needs the halfway
solve read against the condition's own conditioning there, not a distance.
A drum meeting a
face blend's patch at a few ten-thousandths of a radian states up to
5.2e-4 (the stress face-blend scenarios), where the section solved
between two samples lands tens of chords off their chord and is not
filled in.

### Fillets and blends

**Curved face blends refused by name.** `blend_faces` and `fillet_faces`
round curved faces that share no edge where their surfaces share a direction
or an axis, and otherwise march the ball round where the surfaces cross;
between separate faces a march closed round its seat is a band fitted closed
round the loop. Where both faces let the ball go round a closed seat, never
at the same station, each face's stretch is read with the loop opened in the
middle of that face's own gap (its line read over the loop less the middle
half of the gap, so the cut is an open line), the stretches are intersected
round the loop, and the round is a run over each overlap: capped and applied
in turn by `blend_faces`, and by `fillet_faces` a round over each run, with
`trim` or without
(`a_face_blend_rounds_two_runs_where_the_faces_let_the_ball_go_apart`,
`a_plane_and_a_leaning_wall_piece_round_in_two_runs`). A marched line of
contact crossing a face more than once (a face letting the ball go at two
places round its seat, or an open seat whose line runs out of a face and
back) holds the ball over a stretch each time, and every overlap of the two
faces' stretches is a run; with `trim` each face is cut along every
stretch and keeps every piece beyond one
(`a_face_blend_rounds_two_runs_where_a_face_lets_the_ball_go_at_two_places`,
`a_floor_letting_the_ball_go_at_two_places_rounds_in_two_runs`,
`a_floor_its_line_of_contact_runs_out_of_and_back_rounds_in_two_runs`).
Edges a line is cut into that meet end to end, across a face's seam, are
one stretch (`a_round_whose_line_crosses_the_drum_s_seam_stands_untrimmed`).
They refuse (`curved_face_blends_refuse_by_name` in `tests/blends.rs`,
`face_rounds_refuse_by_name` and `curved_rounds_refuse_by_name` in
`tests/sheet_fillets.rs`):

- with `trim`, a marched line of contact crossing a face's seam (the rail
  would be two edges), or crossing a face over a stretch the other face
  does not reach (the face would keep a piece no round meets); and a line
  of contact crossing a face more than once where the seat has a closed
  form (two planes, or surfaces sharing a direction or an axis);
- corners where the surfaces cross in more than one place the ball rolls
  in (two drums fused, a drum's flat and its wall), rounded corner by
  corner by `blend_faces` and `fillet_faces`
  (`a_face_blend_fills_both_creases_over_older_smaller_rounds`,
  `two_drums_round_both_creases_their_walls_meet_in`), are told apart by
  the sum of their reaches, each the larger of its crease's distances to
  the contacts, or else by a line parting the convex polygons that hold
  each corner's section: its crease, its contacts and the meeting points
  of each arc's end tangents
  (`a_face_blend_takes_off_two_corners_whose_balls_overlap_but_whose_corners_stand_apart`).
  Owed: corners whose polygons overlap while the corners do not, told
  apart by their own boundaries; and with `trim` over several corners,
  each later cut keeps the earlier rails, which a cut rebuilding them is
  refused by name. Both are refusals, not wrong answers. A sweep of the
  flat and drum blend (4 flat positions by 12 radii) builds valid and to
  its Green's-theorem volume within 4e-14 wherever it builds. A sweep of
  two fused drums rounded with and without `trim` (6 separations by 5
  radii) gives every face's closed-form area exactly, and the trimmed
  shells are usable;
- a round crossing a face that stands between the two (a ball crossing a
  chamfer between the faces it rounds), where the solid beside the round is
  open where the seat says material, or the reverse. No ball rolls in that
  corner touching only the two faces: on the test's drum, chamfered 1, the
  ball of radius 1.5 touching the wall and the lid has its centre 1.41 from
  the chamfer's plane. The round would leave the material through the
  chamfer, giving round, chamfer, round meeting at creases. It needs a ball
  touching the chamfer as well (a three-face blend). On every path, two
  planes, swept and marched, `blend_faces` refuses by name where another
  face of the solid, on any surface, holds a point inside the ball at a
  station of the run (`between.rs`;
  `a_face_blend_refuses_a_ball_that_crosses_a_bevel_between_the_faces`,
  `a_face_blend_refuses_a_spline_bevel_or_a_bevel_over_part_of_the_edge`,
  `a_face_blend_refuses_a_wedge_its_ball_crosses_over_part_of_the_foot`,
  `a_face_blend_refuses_a_wedge_on_the_foot_s_edges_its_ball_crosses`):
  a chamfer the ball crosses, uneven or not, plane or B-spline, a round
  smaller than the one it would replace, a wedge fused into a drum's foot
  over a short stretch, the foot parted from the top or its edges shared
  either side of the wedge. The run is read whole, halved between stations
  until the centre's travel fits under the face's clearance. Between two
  planes and on the swept path a face inside the ball counts only inside
  its lune, the part between its centre and its arc
  (`a_face_blend_rounds_a_part_thinner_than_its_ball_where_the_lune_stays_inside`),
  and a face on a plane along a straight run is cleared of every lune at
  once from the run's two end sections
  (`a_face_blend_rounds_a_bar_whose_far_faces_pass_just_behind_the_ball_s_centre`).
  The marched path, where the section turns along the run, reads the
  lunes too (`between.rs`, `TurningLunes`;
  `a_marched_face_blend_rounds_a_part_thinner_than_its_ball_where_the_lune_stays_inside`,
  `a_marched_face_blend_fills_over_an_older_round_just_smaller`): each
  lune, from the centre to the contacts in the section through them, is
  cut into cells of radius, of the way between the contacts' directions
  and of the run between two stations. A cell's points lie within half
  its radial stretch and its outer radius times its half-angle of its
  middle point in the middle section, and move over its stretch of the
  run by at most the centre's travel and the outer radius times the turn
  of the direction, each contact's and the centre's travel bounded by
  their arcs' chord and twice the height (the bound the centre's travel
  is read by between stations). A cell whose middle stands farther from
  the face is clear; else it is halved until the face's nearest point is
  placed on some station's section plane by the secant and found inside
  its sector. Owed: the travel bound assumes each arc between two
  stations is convex, as the centre's already does; a proof needs the
  march's step control to bound the arcs' curvature. On every face-blend
  run of `blends.rs`, `sheet_fillets.rs` and `consumed_marched.rs` (77
  runs) the centre's true travel from the nearer station, at 64 samples
  an interval, stays within the bound (worst 0.9999 of it on a curved
  path, equal on a straight one), so no case breaks it yet. Owed with it:
  a second-order clearance near the contacts. Eight faceblend stress
  refusals at seeds 1 to 8 where a ball rolls (tee/76, tee/102, tee/115,
  leaning/23, /75 and /114) are a crease first rounded just smaller than
  the ball (0.1597 under 0.1625 at tee/76): the older round stands about
  6e-6 off the lune, which the cells' first-order bounds cannot show in
  65536 cells, and the blend is refused by name. An older marched round just smaller than the ball (within a few
  hundredths of it) is fitted to 2e-4 and reaches into the ball by a few
  micrometres beside a contact, where its border runs within its own
  tolerance of the face; the blend is refused there, as the fitted round
  does stand inside the ball. Faces the
  crease or a line of contact meets, and faces holding a capped run's
  ends, are where the round runs off (a crease running off a block's
  side) and are not asked. Owed: a face that ends the run without
  meeting either line or holding an end (a narrow face past the one a
  line of contact leaves by) is asked as a face between. The sound rule
  would leave out a face lying wholly beyond a cap plane of the run, as
  no part of the capped round reaches it. That holds only where every
  ball's lune of the run stays on the run's side of that plane, which a
  marched run curving back past its own end does not. The proof needs
  each station's section read against both caps, and the case wants a
  reproducer first. A face's distance
  is the kernel's point-to-face distance: the foot on its surface where it
  lands inside the trim, and where it lands outside on a curved face, the
  local feet inside it sought from samples of the face, as well as its
  edges.

**An edge round whose ball is wider than a face it rounds.** `fillet_edge`
refuses where the ball's contact on one of the edge's faces stands past
that face's far side at every station: the face is narrower than the
setback all along the edge (a sliver face where two oblique boxes are
fused, a triangle face whose apex stands under the setback). Where the
face is a planar strip whose far edge turns the same way as the edge
(convex both: a chamfer, a narrow top), the ball rolls between the faces
either side of the strip and, where the strip lies wholly in the corner
it takes off, the strip goes with the corner (`consumed.rs`). The concave
mirror (a flat in an inside corner) rounds where the strips are planes
whose edges run parallel to the seam over its length and end on planar
faces, square to it or slanted: the ball touches the faces either side,
and the fill is the prism of one section over the strips, cut by the end
planes (`stepped.rs`). The convex seat's ends may be slanted to the edge
too: the corner it takes off is the section's prism between the end
planes, each face's stretch and each end face's trace asked there, and
the cut is cut back to the planes (`consumed_slanted` in `consumed.rs`). Where the far edge turns the other way (a concave seam where two
solids are fused at a slant), the ball touches the other face and rests
on the far edge, its centre on the ellipse where the plane a radius off
that face meets the cylinder a radius about the edge, and where it rests
on that edge alone over the whole crease the round is the band through
its arcs (`rolled.rs`); an end face meeting the far edge's end, touched
there by the end ball, lies behind the plane through that contact square
to the ball's radius and is clear of every ball. Where the ball reaches
the far edge of a strip beyond the narrow face (a wall on a staircase with
a deep drop past its lower step, #720), and every strip between is a plane
whose far edge runs parallel to the seam over its length, ending on two
planar faces square to it or slanted, the round is the cylinder about its
centre's line and the fill the prism of one section, cut by the end
planes and fused on over the strips (`stepped.rs`). Every ball over the
span is cleared of each planar face by its least distance from the
centre's segment, in closed form, so the end faces the end balls touch are
cleared exactly.
A ball that comes down on the plane face the strips run into (the floor
under a shallow staircase) touches it and the face it touches, the steps
going into the fill. Where the strips' far edges are slanted to the seam, the ball may rest
on a convex far edge of a strip beyond (`rolled.rs`, `stepped_rolled`):
its centre runs on the ellipse about that edge's line, the round is the
band through its arcs, and the fill is the band's wedge with a leg across
several hosts, one face on each strip between (`StripPiece` in
`marched.rs`). The walk takes each strip's far edge as the first edge it
crosses square from the middle of the edge before, so a strip that runs
past the seam's ends is walked too. A face is cleared only of its part in
front of the touched face's plane, since every point on or behind that
plane stands a radius off every centre; an end face the touch line meets
is then cleared. Over each interval halved, the part of a face behind the
planes through the far edge's line square to the ball's radius at both
its ends is behind every station's plane between (the interval under a
half turn), so only the part in front of the touched plane and of one of
those is read; a held far edge ending where the end ball touches an end
face is cleared so, each plane moved back by the confusion. Past a
concave far edge over slanted strips the ball may come down on the face
the strips run into (a floor), its centre on the line a radius off that
face and the touched one: the floor goes in as the last strip, its far
line the contact line, and a face on or behind the floor's plane stands a
radius off every centre. Owed there: a ball that changes what it holds
partway along the crease over strips beyond the narrow face (a chain of
`carried.rs` pieces held on edges and faces the seam does not bound): no
run proves its balls over the whole crease, so it is not built. Where the ball changes what it holds along the
crease, the round is a chain of bands meeting tangent where it holds
three things at once (`carried.rs`): a face on each side (the cylinder
between the planes), a face and the other face's edge (the ellipse), or a
face and the other face's corner (the circle over the face about the
corner). Of the 11 stress blend refusals of that kind at seeds 1 to 4,
five round alone: seed 1 cases 12 and 48 and seed 4 case 31 on the edge
all along, seed 3 case 33 from one strip's far edge down onto both faces
and onto the other's far edge, and seed 2 case 54 from the far edge round
the triangular strip's corner. Seed 8 case 46 rounds in its request, edge
to corner to edge. The rest are refused by name. In seed 3 cases 17, 32
and 54 the ball comes onto the face beyond the far edge, where it touches
neither of the edge's faces on that side and rounds the corner between
the face it touches and the face beyond: that is the round of the seam
between those two, which the request does not name, and needs that seam
in the request with a corner closing the vertex where the seams meet
(seed 3 case 32 also starts where its strip's end edge turns the ball's
sections back). In seed 1 case 36 and seed 3 case 4 the ball rests on
the far edges of both faces at once, its centre on the quartic where the
cylinders a radius about the two edges meet, which has no closed-form
piece here. In seed 4 case 48 the ball's contact on the face it touches
runs over a concave seam of that face onto the face beyond. Seed 1 case
12 rounds alone, but another edge of its request shares a face with it:
a band the ball rests on a far edge for, all along or in part, is taken
in a request only for an edge sharing no face and no vertex with another,
and owes its contacts on the face it touches and on the far edge to the
check of each band against the others'. Within a request, a ball whose
contact runs off a face along part of an edge, out of the solid, is cut
off there as alone, its band checked against the others' on the stretch
the face holds, up to the face's side, and on the face across that side
(`carried_contacts` in `support.rs`). Of the 4 stress blend refusals of
that kind at seeds 1 to 4, seed 1 case 6 and seed 3 cases 13 and 14
round, each to its edges rounded alone with the joint at a shared corner
as at a radius that fits (and six more at seeds 5 to 8). Seed 4 case 32 is refused by name: the side
its band is cut off at stands inside another seam's band on the face
across. Owed, alone and in a request alike: past a face's side,
a concave band's fill stands as a flat on the face's plane, where the
ball rests on the side instead (`carried.rs`), a run read only from a
ball holding that side at the crease's middle.

**A pole leg whose edge winds round the pole and whose rail does not.** A
marched band whose rail passes over a sphere's pole has a leg from the
winding loop to the pole with the other loop cut out as a hole. When the
winding loop is the rail, `pole_leg` in `marched.rs` splits it clear of the
hole and starts the seam there. When it is the edge, the seam leaves from the
edge's own vertex and can cross the hole. No measured part builds one.

**Variable radius on chains.** `fillet_edges_variable` rounds a straight
edge between planes under any law, a closed loop of curved edges marched
with the law's radius per station, and an open curved edge whose ball runs
out past its ends: the law is held at the end's radius past each end, so
the run-out is the constant round's at that radius. An open chain of arcs of
one curve (a crease split at a vertex) is that one seat over the window its
arcs cover. A chain of different seats continuing each other end to end
(a rounded-rectangle pad's top loop) gives each edge the law over its own
stretch, stops each band in its end sections and joins the wedges across
each tangent junction's plane. A chain of straight edges between planes
turning on one shared planar host carries each band past a turning vertex
at the law's radius there (the constant band of that radius, which the
junction is read for): run out across a convex turn, carried to the
bisecting plane across a re-entrant one. The law is kept as given over
each edge, its value and slope at the edge's ends its own: a turn's mitre
is G0 whatever the slope, so pinning the slope level there would change
the law and buy nothing. A loop joined round one convex turn is opened at
a straight band built in two halves. A marched band whose crease ends on
a planar face holding its end section stops at the vertex, capped in that
face's own plane. A chain turning where marched bands meet is met in
`fillet_edges`' own rounds, each piece under its window of the law. A straight edge between planes rounded alone whose
contact runs off a face beside it (a strip beside a groove) is read
station by station at the law's radius as the constant round reads its
own, cut off in the open and refused in the material, and owed the same
beside a narrow face in a chain (the constant round carries such a band's
contacts past the face's side for the other bands to be checked against,
`carried_contacts`, which under a law change along the edge) and on a
marched seat. Owed:
a turn where a marched band meets a straight band between planes (the
constant round's band there is an exact torus band; under a law it is
fitted, and its end meets the straight band in what the boolean does not
melt, built first or second, capped in its wall or run out); a convex turn whose sections differ where the constant
band beside it is trimmed by a face it runs into (its run-out is no
straight band ended square, and the law's band is joined to it at the
vertex as before; untried). Across a convex turn whose sections differ
otherwise (a drafted wall's foot leaning out under the far face), a band
running out past the vertex is built on past it as one surface, the law
held level, as far as the constant band runs out: a seam at the vertex
would stand inside the material where the other band crosses it, which
the cut does not close. A re-entrant turn whose sections
differ keeps each band under the law to the vertex and carries its end
section on level past it, the widening past its wall a prism from behind
where the other band reaches, so no face of it lies in the other band's
wall.
A straight edge along a drum's ruling
under a law is marched, the exact sweep of `ruled.rs` holding one section
only; a seat pinched where its hosts turn tangent solves each station and
its pole at the law's radius there.

**Setback corners beyond the one-ball planar vertex.** `round_vertex_setback`
fits one trimmed patch G1 to the fit target at a convex planar corner one
ball touches. Owed: G2 (needs band sections curvature continuous with their
faces, so the end arcs stop bending the patch where the faces do not), the
patch as N quads from the sides' midpoints to a centre fitted by
`fit_network` so tangency holds exactly, vertices where a curved face
meets or the planes share no tangent ball, and setbacks asked through
`fillet_edges`' corners. `fillet_edges_with` closes corners of three
different radii with the same patch, a box corner at 1, 3 and 2 included (its
narrow band's end turns into a face's trim, and the patch's tangency there
closes about twofold a net; it meets the target at 135 controls). The
quads round a centre with every side held exactly (`fit_loop_quads`) meet a
pyramid's apex but not yet a box corner's six-sided opening: the centre's
corner holding cannot follow the larger opening (branch
`feat/corner-quads`, #718).
Where a convex and a concave band of one request meet at a vertex
(`setback.rs`, `signed_corners`; #712), the request is rounded in one
order whatever its own, and a vertex of three planar faces and three
straight edges, all rounded, is closed as one corner: each band stopped
its own width past the least setback at which its section clears the
others, the box square to the three edges from each section back to
twice the widest band behind the vertex cut away, and the opening of
three arcs and three conics filled by one patch G1 to every band end and
face. On a face reflex at the vertex (an L-bracket's front) the face's
part kept in the box is the side of its conic away from the vertex, out
to the box's walls; coplanar faces a fusion left at the vertex count as
one. Where an edge of the vertex is left sharp, the request is refused by
name in either order: the opening between the band ends crosses that
edge, and a patch tangent to both its faces there needs two tangent
planes at one point. A patch with a collapsed side at the crossing closes
such an opening on a convex corner (`sharp.rs`,
`round_vertex_setback_sharp`): the sharp edge cut at its own setback at
`Q`, a three-sided patch there with its fourth side shrunk to `Q`, fitted
by `fit_network` with its halves of the trims from `Q` held and held across
the two faces, its normal turning from one face's to the other's along the
shrunk side, and one `make_filling_n` patch G1 to it, the bands and the
faces over the rest. `blend_edges` closes a convex corner of one round,
one bevel and one edge left sharp the same way (`close_sharp_corner`): the
sharp edge cut as far out as the farther band, and the rest meeting the
bevel's end and the trims from it in position only, as #803's corners do.
Owed: the same at a signed corner (the L-bracket of #712, where the cell is
not convex) and through `fillet_edges`' own corners; signed corners of
more than three edges or with a curved face, refused by name.

### Construction and exchange

**STEP pcurves the file cannot state exactly.** The reader keeps a file's
pcurve where its lifted gap is within twice the projection's, restated in
this kernel's units (the file's length and plane angle units, applied to
each surface parameter as ISO 10303-42 parameterizes the surface) and,
where the file states it at another pace or phase than its curve's, read
at the affine map of the curve's parameter that puts it on the edge's
ends: a conic either way round and turned, a line at any pace, a spline
over its own domain either way. Still derived: a trig or offset trace and
a pcurve of a curve written as a spline conversion, none of which STEP
states at a linear map of the curve's parameter; an ellipse turned other
than by half turns; a conic on a chart stretched unevenly or on a line swept
into a plane (whose chart here is the file's sheared where the sweep leans,
taking a circle to an ellipse at a phase no conic here states; lines and
splines on it are carried through the map and kept). A surface replica's
chart is its parent's with each length parameter scaled by the replica's
scale, and its pcurves are restated so. What blocks each still derived: a
trace and a spline conversion have no linear map of the curve's parameter to
state; the conics need a conic stated at a phase other than half turns from
its major axis, which no conic here carries. The writer still leaves out
a reversed conic pcurve: written forward, it would ask other readers to
take a pcurve at another pace than its curve.

**Edges through a pole on reading.** Both readers cut an edge through a pole
of a face it bounds (a sphere's pole, a cone's apex, a patch's collapsed side)
before any face is built. IGES trimmed surfaces (143, 144) cut every boundary
segment through a pole of any surface read with them, so neighbours sew piece
to piece. Not cut: a trim placed by one subfigure instance against a pole of
a surface placed by another, and pcurves other producers build through a
pole. A subfigure definition is built once and shared by every placement,
so the first would cut the definition for poles only some placements meet,
and needs such definitions built per placement; the second has no file at
hand that writes one. The cut edge stays in the model unused, and a STEP style or shape
aspect naming its file id lands on every piece.

**A sketch curve through a pole of an integrated chart (#729).**
`SurfaceSheet::chart_pieces` cuts a curve where it lands on a sphere's pole
or a cone's apex, both found in closed form; through the collapsed side of
a revolved or spline surface it is refused by name ("crosses a pole of
the surface's chart"): a geodesic integrated in the chart that would run
on across the collapsed side is stopped there, the chart's metric having
degenerated, and no landing is integrated past the pole. Those landings
are integrated, so the pole's preimage on the sheet is not known in closed
form. Building it needs that preimage located by the integrated map's own
distance to the pole, with the geodesic near it followed in a chart
regular there.

**A network's mixed point, held tangent where that fits better (#653).**
`fit_network` leaves a point where a three-sided patch meets four-sided
ones near tangent, its angle in `Network::worst_tangent`: the gluing
there misses the fit's cross derivatives by 2 to 4.5 per cent however fine
the spans. Holding it anyway, with the points' parameters corrected under
the hold (branch `wip/network-mixed-gluing`), is exactly tangent and fits
better on a noisy ball's cap (9.3e-3 off against 1.5e-2 near tangent) but
stalls on the clean one (1.6e-3 against 2.3e-4). Owed: run the held fit
beside the near-tangent one and keep it only where it is exactly tangent
and no further off the points than the near-tangent fit, so no result
gets worse.

### Mesh conversion

`solid_from_mesh` rebuilds planes, the four canonical surfaces, extrusions,
surfaces of revolution and fitted B-spline patches, and leaves every other
region faceted; `MeshSolidReport::fallbacks` names each curved region it
faceted and why. It reads meshes drawn in strips along rulings, as
other kernels commonly draw drums and cones, as it reads ogeom's own. Each
change is held to the converter's standard (verified against every sample,
edges placed on both surfaces, facets as the fallback) and measured on
`mesh_corpus` and the truth bench before and after.

**Seams between fitted regions.** Each region is fitted on its own and each
seam solved point by point along the mesh boundary, or threaded as a chord
through the boundary's vertices. Folds and faces that mesh open are caught
after the build and their curved faces withdrawn to facets, so every
measured part tessellates closed at the cost of those faces. A chord where no
point carries onto both surfaces (a torus and a cylinder fitted to pieces of
one spline blend) can stand 0.2 off, stated as its tolerance. A single facet
whose seam with a curved face would be a chord wider than the reach is built
as a fan onto that face instead (a wedge split at its corner first), its
seam lifted onto the face. A plane of two facets gives up the one along such
a seam into a plane of its own, which is then fanned the same way, the other
then a facet alone. A facet whose seam the plan took within the reach and the
build placed past it is fanned too, and a facet with a curved face across
every side whose seam stands past the reach is split at its middle into
three, each a facet with one curved side, and the conversion made again. On
`mesh_corpus` every edge between a curved face and a plane of one or two
facets then holds within the reach except:

- on the sliver part, where a plane of two facets beside a torus keeps
  seams of 0.019 against a reach of 0.004.

Design items, in order:

- **The sliver part's plane of two facets**, not yet diagnosed.
- **Seams between two curved faces fitted to one spline.** On `mesh_corpus`
  one such seam stands past its reach: on `nist_ctc_02`, a cone (half-angle
  3°) beside a torus (major radius 112, tube 10), at 0.097 against a reach
  of 0.093. They meet nearly tangent with no point on both, and the torus's
  axis stands 2.0° off the cone's. Holding them tangent along a circle
  needs them coaxial, and turning the torus onto the cone's axis moves its
  samples millimetres off, past its fit; the fix lies in fits that share
  their neighbours' tangency or take the spline whole. Every other seam
  between two curved faces on the corpus holds within its reach.
- **Junction rules** (Benière et al. 2012: a junction stays only if its
  faces are pairwise adjacent, dangling edges dropped). No part of
  `mesh_corpus` shows a junction they would change; they wait for one.

Tracing a seam as a branch of the two fitted surfaces' intersection has
nothing to trace on the measured parts. Of the 427 chords on `mesh_corpus`,
258 stand wider than the reach (twenty coplanar distances), and 254 of those
because one surface stands farther than the reach off the other at a point
of the mesh path. A signed distance moves no faster than the point it is
taken at, so no point on both surfaces lies within the reach of that point,
and Newton's step from the path lands 0.5 to 15 off it. The widest are all
of this kind: `nist_ctc_02`'s planes and tori fitted to pieces of B-spline
faces (0.63 over a span of 28, the torus's equator bowing off the plane by
its sag), its cylinders beside planes (0.46), and its tori beside cylinders
(0.11, meeting at a quarter of a degree to seven). The other four are cubics
through three vertices standing a little wider than the polyline would. Of
the 169 within the reach, 3 have a crossing within it all along. The fix is
upstream, in fits that share their neighbours' tangency or take the spline
face whole; the chord stays as the fallback. The planes beside those tori
and cylinders were single facets, now built as fans (above): their seams lie
on the curved faces, and of the main corpus with the step-conformance parts
the widest seam on a curved face is 0.35 (p99 0.046), against 0.63 (0.16).

**Fillets fitted apart from their supports.** Derived from their supports:
rounds between two planes, a torus between any two supports on its axis
(cylinders, cones and planes square to it, or one of them and a sphere
centred on it), corner balls, and a sphere put on the torus a round turns on
beside a plane. Not derived: a torus between two spheres, a piece of a torus
between curved supports met as a sphere, variable radius, and fillets between
curved faces off a common axis, which keep their fitted surfaces and meet
their supports at a near tangency. A torus corner a few facets round between
two fillet cylinders and a floor fits nothing (every sample reaches into the
cylinders) and stays facets.

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
  edge, under 1 inside;
- a ring (a spline blend round a pocket) is no disk, and cutting it into
  disks does not help on a coarse mesh. On `nist_ctc_02` its blends are
  rings four facets across (radius 10, chord 1.17). Halved, then halved
  again, the pieces fail: the mean-value map flattens triangles with all
  three corners on one side of the square, a half ring folds over its mean
  plane, and the quarter rings' patches stand 3 to 5 percent past a
  triangle's sag at its middle, or hold too few vertices inside. One piece
  of eight verified, and the facets left beside it widened the worst seam
  from 0.35 to 0.73. A ring needs a chart that goes round it (periodic in
  one direction) and enough rows across to verify against.

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

**Converted parts of the real assemblies.** Fifty solids of the DoomCube 350
assembly and fifty-three of the Voron 2.4, each meshed at a thousandth of
its diagonal and converted back, come back valid within 3e-3 of their
volume. Without recognition (`MeshSolidOptions::recognize` off) the build
asks nothing of the solid it returns, so a mesh passing through itself comes
back a solid crossing itself there; with recognition the conversion is
refused, naming where.

**A face whose ring runs back along itself is meshed, not healed.** Where a
face's boundary runs out along a path and back along it within the widest
tolerance of its edges (a spline edge rising to a circle and running on
along it to the vertex they share, a hair past it in the chart), the ring
crosses itself at any chord. On its last pass the mesher takes the tail out
of the ring and draws it as a strip at its edges' points
(`ogeom-mesh/src/tails.rs`), and the rest triangulates whole. The model
keeps the tail: the face still bounds a stretch of no width between two
edges its neighbours do not share with each other, and `fix_shape` strips
such a pair only where no other face walks either edge. A tail whose way
back leaves the way out across a segment wider than that tolerance is left
crossing, as is a ring that crosses itself anywhere else.

### Speed

**A per-face state cache.** The boolean's build phase could read a piece's
classification instead of probing for it. No case needs it for correctness.

**The rest of the local boolean.** Faces whose box misses the other solid's
are set aside (`ogeom-bool/src/aside.rs`): not gathered, split, classified or
rebuilt, and passed through as the same nodes; a plane's holes clear of the
other solid are left out of its arrangement and put back into the piece
holding them. What a boolean reads of a solid of 64 faces or more (the
setting aside's walk of every face, the classifier's prepared faces and
the rings drawn of them) is kept in the model (`Model::kept_read`) and
read back by the next boolean on the solid unchanged; any edit below it
forgets it. A corner hole into a plate of 582 faces costs about 1.5 times
the same hole into one of 22 (heavy test `local_boolean`, bound 2.0).
Still linear in the solid: `side`, the history copies and the shells of
the faces passed through, about 0.4 ms at one thread on the 582-face
plate, and the first boolean on a solid reads all of it. A solid read for
the first time takes each face it shares with a solid whose read is kept
from that read (`Model::kept_reads`): a boolean's result takes the faces
it passed through from its operand's setting aside and prepared boundary,
rings drawn included, and reads only the faces made; the nodes below it
are walked down to the nodes a kept read holds. The setting aside asks
whether a face of the tool lies outside the solid of the boundary the
general fuse classifies against, prepared once. Along a chain of 200
holes and bosses into a plate (heavy test `local_boolean`) a step at 306
faces costs 1.31 ms against 1.48, the chain 8 percent fewer cycles and 12
percent at 361 steps, to the same result; at 66 faces nothing changes.
Taking the faces' reads is under 2 percent of that at 200 steps: the
edges' holders are rebuilt in the solid's own order from each face's walk
as a fresh read builds them. Still linear along the chain: those holders,
`side`, the history copies, the shells' closed check and each new face's
rim walk (`explore_unique` per face; a plain walk saves 2 percent).
`edge_pcurve_gap` answers a line, a circle or an ellipse against a line,
a circle or an ellipse on a plane, or an iso line of a cylinder, a cone,
a sphere or a torus, in closed form where that stands within the edge's
tolerance (`pcurve_gap/closed.rs`), within rounding of what its sampling
measures; the chain then costs 22 percent fewer cycles at 200 steps and
17 at 361, and 80 percent of the gaps the stress suite measures are
answered so. Still sampled: B-spline pcurves, the chart traces of
oblique sections, helices and any gap past the tolerance. A STEP read
measures the pcurves it derives with the same forms where it acts on the
gap only past the confusion distance
(`pcurve_fit::lifted_gap_past_confusion`): the Voron 2.4 assembly reads
in 23 percent fewer cycles at one thread. Still sampled there: the
derived pcurve's gap where the file states pcurves of its own, which
bounds the ones kept.
The sampled curve crossing in space passes over a segment whose box
stands further from a sample than its nearest so far
(`curves::Segments`): `fillet_box_all` 26 percent fewer cycles,
`boolean_marched` 9. Measured and not landed: the classifier draws a
cylinder wall's rings at the full chord though its pcurves are straight
in the chart, about 15 percent of `boolean_drill` and of a step of the
chain and 40 percent of `offset_part`; drawn at the fewest points the
walk's period choices (rings placed by the mean of their points) put a
NIST part's wall a period apart in places, so the fewer points are not
the same rings. Boxing the section against each edge before crossing
them rejects almost no pair on the benches. Where the assembly around the faces set aside does not close (the
sew welds a shared vertex, a seam join reaches a shared edge, a strand would
widen a shared vertex), the boolean runs again over every face; converted
slabs with tolerant vertices take that path.

**Operations that edit in place.** A boolean notes the nodes it passes
through as held by its operands too (`Model::note_held`). Heal, sew, the
fillets and chamfers, `remove_faces`, `unify_same_domain` and `tessellate`
copy what is held below the shape they are given before editing
(`Model::unshare`, `ogeom_algo::on_own_nodes`). Still edited where they
stand: a node held whole and handed in itself (a face set aside, given to
`fix_face_pcurves`); a profile's edges, which a sweep (`make_prism`,
`make_revolution`, the pipes and lofts) and `split_face` describe on the
surfaces they build; a boolean's own vertex widening and re-charting of a
ball. Each grows a tolerance, or adds a pcurve on a surface none of the
operand's faces lie on (clearing the edge's same-parameter flag). Copying
there instead parts the profile from the solid it came from: a decision
for those operations, not taken here. Sharing made by
anything but the boolean (`Reshape`, compounds, placed instances) is not
noted, and the notes are not saved with a document.

**Fillets in fewer passes.** A chain is taken in rounds, no two edges of a
round sharing a vertex, each round in one boolean each way (wedges whose
boxes meet in as few booleans as keep each one's boxes apart). Everything
else in a call reads only what is near each edge, so a call costs its
booleans plus a few milliseconds per edge: the four outer edges of a plate
with 256 holes take about 2.5 times the plain plate's time, every rim of
it about 25 ms a hole. What is left is the boolean's own cost over the
whole solid, and a polygon's outline still takes two rounds. Building the
later round against the solid as given, as though the earlier cuts were
made, and cutting every wedge at once (fused first, or as disjoint pieces)
gives a tool the boolean cannot close against even a plain box ("the kept
pieces did not close into a shell"): its legs lie on the box's faces all
round the outline.

### Quality and performance

From a read-only audit of every crate (2026-10-07). Nothing below is
measured yet: each item starts with a bench, and lands only where paired
runs beat the threshold it declares (the performance-optimizer skill).
Kernel threading already runs through `ogeom_core::parallel::map_ordered`
(scoped threads, results in item order, bit-identical at any thread count).

**Benches.** `tools/ogeom-bench` (see its README) reports min, median
and MAD per bench after a warm-up, loops short benches to 50 ms of
samples, takes `--threads N` and `--filter`, and `threads()` honours
`OGEOM_THREADS`. The baseline is recorded at one thread so its ratios
compare across core counts. Benches now cover mass properties, `check`
and `fix_shape` on the largest corpus part, sew of a 1536-face shell,
`solid_from_mesh` on a tessellated part, STEP write, IGES read (of the
kernel's own IGES of a corpus part: the corpus has no IGES file), weld and
STL read at 1M triangles, `project_exact` and the polygonal HLR, a
thickened spline face, a pipe through circular sections, a guided pipe, a
marched fillet and a box with all 12 edges filleted. Still missing: a STEP
file of 10 MB or more (the largest corpus file is 1.2 MB). On the
largest part `solid_from_mesh` does not finish a warm-up and three samples
in two minutes, so its bench runs on the small part. `project_exact` takes
about 30 s at one thread there (7 s at 20), nearly all of it Newton on its
spline faces, so its heavier bench runs on a 520-face part.

An instruction-count gate in CI would be deterministic on a loaded runner
and looks feasible: `iai-callgrind` needs valgrind on the runner (an apt
install on ubuntu-latest) and a bench target with a few of the benches
above at reduced size, since callgrind runs them 20 to 50 times slower;
`perf stat -e instructions` needs `perf_event_paranoid` at 1 or below,
which hosted runners do not promise. The work: an `iai` bench crate, a
recorded count per bench, and a gate at a few percent.

**Cheap and certain.**
- `map_ordered` starts fresh threads for every stage, two items included.
  Measured, and left alone: the items are heavy (faces, edges, face pairs)
  and few (a drill hands over 2, 2 and 5), so the default count already
  beats `--threads 1` on the small benches (boolean_drill 1.0 against
  1.45 ms, fillet_block 20 against 37 ms). A minimum batch per worker
  serialises those stages: 2 items per worker cost fillet_box_all 3%,
  4 cost it 10%, 16 nearly doubled fillet_block. A scoped spawn costs
  8 to 16 us per thread here, about 8% of boolean_drill and under 1% of
  the fillets. Letting the caller take items as one of the workers gained
  3% on boolean_drill and 7% on boolean_local, under the 10% bar; a
  persistent pool could win at most that 8% on the smallest booleans.
- Closed-form feet on curves: `project_on_curve` in
  `ogeom-algo/src/measure.rs` samples a line or a circle like any curve
  (`line_parameter` and `circle_parameter` exist in
  `ogeom-math/src/elementary.rs`), and refines with Brent's minimizer
  where Newton on `(C - P) . C'` after the bracket would do. Surfaces
  (cone, torus, extrusions and revolutions through their profile) already
  have closed-form feet. Reading the DoomCube and Voron 2.4 assemblies
  projects only onto splines this way (33,000 calls, 0.2 s of 27 s), so
  the gain is for other callers.
- The intersection walker: six or more unknowns still go to the heap
  solver, `Condition::system_at` allocates on every evaluation, and each
  step corrects from the last point rather than a tangent prediction.
  Measured, and left alone: the correction evaluating the surface pair
  into fixed storage, with no vectors per Newton step, takes 3.2 percent
  fewer user instructions on `boolean_marched` and 0.8 percent on
  `fillet_marched` (pinned, deterministic), to the same result; under a
  5 percent bar. The blend's condition evaluates its guide through
  `derivatives_at`, which allocates as well.
- One jet call where separate point, first and second derivative calls
  are made: `CurveOnSurface::derivatives_at`, the mesh lift (point then
  normal), and the side probes of mass integration; allocation-free curve
  derivatives (a curve jet; `derivatives_at` returns a `Vec`).
- Box rejection: `general_3d` in the curve crossings samples both curves
  before any whole-curve box test, and `section_face.rs` does not filter
  edges by box; section-against-edge paving tests the other face's box,
  not the section's own; a spline's sub-range is bounded by its whole
  control polygon.
- The guided sweep law samples each guide at 257 points at every station
  and bisects each sign change 60 times.
- The mesh still draws an edge again for each face's rings and chords and
  in each settling round; its sealing passes build the edge-use map three
  times. Measured, and left alone: on `nist_ctc_02`, `nist_ctc_04` and
  `nist_ftc_10` at a ten-thousandth of their diagonals, drawing each edge's
  polyline once per call instead of once per face saves 2.4 to 2.8 percent
  of the triangulation's instructions, under a 5 percent bar, and
  everything after the faces (weld, folds, seal) is a tenth of it, the
  edge-use maps about a hundredth.
- Primitive refinement in the converter builds a new matrix every
  iteration. Piece nesting answers a piece outside another's box without a
  ray (`mesh_to_solid_pieces`).

**Algorithms worth the work.**
- Classify a boolean's pieces once per connected region, not per piece
  (Requicha and Voelcker 1985): pieces meeting across an unpaved sub-edge
  share a state. Measured: classification is about 5% of
  `boolean_many_faces` (36 probes over four drills, faces the tool misses
  already set aside) and 0.1% of the fillets, so it waits for a case where
  pieces are many; it needs the split and classify of `general_fuse_as`
  apart first. The boolean's time is the arrangement (about 2.9 ms of the
  four drills) and assembly (about 2.1 ms).
- Surface feet on a trimmed or offset spline still seed from a grid; the
  span search on a plain patch could clip its spans to a trim window, and
  bound an offset by its basis's boxes grown by the offset distance.
- Curve-surface intersection prepares its sample grid once per face and
  clips the ray to the face's box; point-in-solid rays pay it on every
  call now.
- March sections at the cubic's own interpolation error and check the fit
  between stations (Bajaj et al. 1988; Barnhill and Kersey 1990), instead
  of a straight-chord step about ten times finer than the fit needs. The
  fillet march takes stations eight times denser than its fit for the
  same reason. Measured on two crossed cylinders: the joint fit is half
  the fuse (5 to 7 ms a branch from one span through 1556 points); a fit
  started from knots placed by the cubic's error estimate is a third
  faster with a tighter true miss, but moves sections by sub-micron
  amounts and the shells of a few stress drills and two tests close or
  refuse the other way (patch kept as `hermite-stations-fit.patch` in the
  scratch crate). Shell closure that flips on sub-micron section changes
  is the item to settle first.
- March general pairs over the faces' parameter boxes, and cache a pair of
  surfaces marched once within one fill.
- Marched seeds search the cell pairs that come close without crossing
  for a crossing or for the stationary approach a small loop runs round
  (Sederberg and Meyers 1988). Where the surfaces lie close along a whole
  curve (a blend against its support) the squares around a search that
  found nothing are passed over, to keep the cost within the grid's, so a
  small loop beside such a band can still be missed. Normal cones over the
  spline patches' hulls would prove those squares empty instead.
- The winding number fallback of the exact classifier (where every ray
  meets a tangency) sums the faces' meshes triangle by triangle; a
  hierarchical far field (Barill et al. 2018) once a solid of many faces
  asks it often. Exact solid angles of trimmed planes would spare meshing
  them. A probe within a curved face's mesh band at the finest chord still
  refuses.
- Arrangement darts ordered by exact orientation, ties walked apart.
- Primitive fits with analytic gradients and Levenberg-Marquardt (Lukacs,
  Marshall and Martin 1998).
- Fitting a band or a pcurve: factor the normal matrix once per round for
  all rows, a cyclic system as band plus border; warm-start knots from
  the last round.
- Exact offsets of extrusions and revolutions through their profile's
  offset; spline offsets through the local-refinement fit.
- Mesh `simplify` by quadric error metrics with a heap and the link
  condition (Garland and Heckbert 1997); the triangulation keeps triangles
  by the exact flood rather than a centroid parity.
- Fillet corners cut together where their tools stand apart, not one at a
  time against the whole solid; the ruled fillet through the round's
  batching; the curved corner probes the exact boundary instead of meshing
  the whole solid.

**Structure.** `fill` (2,300 lines) and `general_fuse_as` (1,500) in the
boolean, `mesh_solid.rs` (12,900 lines), and the largest builders in the
fillet and offset crates split into named phases with typed outputs, which
the region classification and the parallel read-only phases need. Mesh
conversion moved to its own crate would stop it rebuilding the whole
graph.

**Consumers' builds.** Measured and documented in the README ("In your own
crate"): on a dependent crate, `x86-64-v3` takes 23 percent fewer
instructions and CPU time converting and refining a 29,003-face solid, thin
LTO with one codegen unit 3 percent fewer instructions and 6 percent less
time, both 26 and 25 percent, to the same result on that solid and on
`mesh_corpus`. The workspace keeps the baseline target: a consumer chooses
what its machines run.

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
- **Scope is what `docs/SCOPE.md` admits.** It is normative and says how
  to decide a case. Scope is about capability, not structure.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against,
  bundles or imports another kernel, no one working here reads another
  kernel's implementation source, and the design is worked out here, not
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
- **Each face's box is kept in the model.** The model keeps a face's tight
  box in the face's own frame, found the first time it is asked for
  (`ogeom_algo::face_bounds`, `tight_bounds`, the boolean's pair filter) and
  read in constant time after (`Model::face_bounds`). It is forgotten when a
  node of the face is handed out by `node_mut` or a surface by
  `surface_mut`. `tight_bounds` of a 582-face plate a second time costs a
  twentieth of the first. Faces a boolean rebuilds are new nodes, so a
  result's boxes are found again (see "Sharing unchanged faces' nodes").
- **The converter uses no learned model.** Published learned reconstruction
  (2024 to 2026) gives a valid solid for 70 to 76% of parts under a hundred
  faces, every face a B-spline, and places its own bottleneck in rebuilding
  the b-rep, which is the converter's verified build. A learned proposal is
  also a result nothing here measures. Global selection by an integer program
  is out for the same reason: its time has no bound.
- **SAT, X\_T and JT are not read or written.** Their specifications are
  unpublished, and implementing them would mean reverse engineering files
  instead of reading a standard. If a specification is published, this lifts.
