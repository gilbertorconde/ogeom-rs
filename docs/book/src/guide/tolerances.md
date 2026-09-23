# Tolerances

Imported models have gaps, intersections are computed numerically, and
points that should coincide rarely match exactly. ogeom handles this with
three rules.

## 1. No global epsilon

Every operation takes a `Tolerances` argument. The examples in this guide
use `Tolerances::millimetres()`.

- `confusion()` is the distance below which two points count as the same
  point.
- The other thresholds derive from it.

This keeps the units and expected precision visible at each call site.

## 2. Tolerances are per entity and only grow

Each vertex, edge and face carries its own tolerance: the radius within
which its stated geometry is trusted.

- Entities of a clean primitive sit at the baseline.
- Imported or heavily modified models have wider tolerances where the
  geometry is less certain.
- Operations may widen a tolerance to record real uncertainty. Example:
  when `sew` merges two vertices a micron apart, the surviving vertex
  widens to cover both.
- Nothing narrows a tolerance silently, because that would claim precision
  nobody measured. `ogeom::heal::reduce_tolerances` narrows tolerances by
  re-measuring the geometry.

## 3. The containment rule

- An edge must lie within the tolerance regions of its faces.
- A vertex must lie within the tolerance regions of its edges.

The validity checker (`ogeom::algo::check`) accepts a gap within the
stated tolerance and rejects one outside it. It is no stricter than the
builders, so it never rejects what the kernel legitimately builds.

The full semantics (what each entity kind's tolerance means, and why) are
in §5 and §9 of [the data model](data-model.md).
