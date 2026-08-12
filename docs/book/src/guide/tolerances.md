# Tolerances

Real geometry is inexact: imported models carry gaps, intersections are
computed numerically, and two points that should coincide rarely do to the
last bit. A kernel's tolerance model is how it stays honest about that,
and ogeom's has three rules.

**There is no global epsilon.** Every operation takes a `Tolerances`
argument — `Tolerances::millimetres()` in every example in this guide —
which sets the scale-appropriate thresholds: `confusion()` is the distance
below which two points are the same point, and the rest derive from it.
Passing tolerances explicitly is deliberate friction: it makes the unit
system and the precision expectations of a call site visible at the call
site.

**Tolerances are per entity, and they only grow.** Beyond the baseline,
each vertex, edge and face carries its own tolerance — the radius within
which its stated geometry is trusted. A healthy primitive's entities sit
at the baseline; an imported or heavily-operated-on model carries wider
ones where the geometry genuinely is less certain. Operations may *widen*
an entity's tolerance to record honest uncertainty (a sew that merges two
vertices a micron apart widens the survivor to cover both); nothing ever
narrows one silently, because narrowing is a claim of precision that was
not measured. `ogeom::heal::reduce_tolerances` narrows them the only
honest way — by re-measuring.

**The containment rule ties it together.** An edge must lie within its
faces' tolerance regions, a vertex within its edges'. The validity checker
(`ogeom::algo::check`) verifies exactly what the builders promise — it
accepts a gap that is within stated tolerance and rejects one that is not,
and it is deliberately no stricter than the builders, because a checker
that rejects what the kernel legitimately builds teaches people to ignore
it.

The full semantics — what the tolerance of each entity kind means, and why
the model is not negotiable — are §5 and §9 of
[the data model](data-model.md).
