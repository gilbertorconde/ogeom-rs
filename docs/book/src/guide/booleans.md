# Booleans

The four standard booleans are in `ogeom::boolean` and share one
signature:

```rust,ignore
let out = ogeom::boolean::fuse(&mut model, &a, &b, tol)?.shape;    // union
let out = ogeom::boolean::common(&mut model, &a, &b, tol)?.shape;  // intersection
let out = ogeom::boolean::cut(&mut model, &a, &b, tol)?.shape;     // difference
let out = ogeom::boolean::section(&mut model, &a, &b, tol)?.shape; // the curves where they meet
```

The [getting-started example](getting-started.md#a-first-solid) is a
`cut`, checked against its exact volume.

## Other operations

| Function | What it does |
|---|---|
| `cells` | Full cellular decomposition of two solids: every region classified against both inputs. The four booleans select from this. Use it when you need a different selection. |
| `fuse_fuzzy`, `cut_fuzzy` | Take an explicit fuzz distance, for inputs whose faces almost coincide (mostly imported geometry). Avoids the sliver faces an exact operation would create along the near-contact. |
| `make_volume` | Builds the solids enclosed by an arbitrary set of faces. |
| `remove_faces` | Defeaturing: deletes a feature's faces from a solid and closes the gap. Neighbouring faces extend to fill it, or the band is re-intersected where extension cannot close it. |
| `make_periodic` | Prepares shapes for repetition along an axis. |

## Guarantees

- **Tangent cases work.** A tool tangent to a face, even at a
  parametrisation singularity of that face's surface, produces the correct
  section curve. The test suite checks these cases against closed forms
  (sphere octants, blend corners).
- **Unresolvable cases are refused.** If two inputs interfere in a way the
  algorithm cannot resolve correctly, the operation returns a named error
  instead of a shape that is wrong.
