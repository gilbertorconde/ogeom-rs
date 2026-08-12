# Booleans

The four classics live in `ogeom::boolean`, all with the same shape:

```rust,ignore
let out = ogeom::boolean::fuse(&mut model, &a, &b, tol)?.shape;    // union
let out = ogeom::boolean::common(&mut model, &a, &b, tol)?.shape;  // intersection
let out = ogeom::boolean::cut(&mut model, &a, &b, tol)?.shape;     // difference
let out = ogeom::boolean::section(&mut model, &a, &b, tol)?.shape; // the curves where they meet
```

The [getting-started example](getting-started.md#a-first-solid) is a
`cut`, measured against its closed form.

Beyond the classics:

- **`cells`** computes the full cellular decomposition of two solids —
  every region classified against both inputs — which is what the four
  classics select from, exposed for callers that need a different
  selection.
- **`fuse_fuzzy` / `cut_fuzzy`** take an explicit fuzz distance for
  inputs whose faces almost coincide — imported geometry, mostly — where
  the exact operation would produce sliver faces along the near-contact.
- **`make_volume`** builds the solids enclosed by an arbitrary set of
  faces.
- **`remove_faces`** deletes a feature's faces from a solid and heals the
  wound — the defeaturing operation: neighbours extend to fill, or the
  band is re-intersected where extension cannot close it.
- **`make_periodic`** prepares shapes for pattern-repetition along an
  axis.

## What the boolean promises

Tangencies are handled, not wished away: a tool tangent to a face — even
at a vertex of its own surface's parametrisation — produces the section
curve it should, and the suite holds those cases to closed forms (sphere
octants, blend corners). Where two inputs genuinely interfere in a way
the algorithm cannot resolve honestly, the operation refuses by name
rather than returning a shape that looks right until it is measured.
