# outside

Code that was written here but is not part of the kernel.

`docs/SCOPE.md` sets the kernel's bounds: parity with the reference kernel's
modelling modules (FoundationClasses, ModelingData, ModelingAlgorithms,
DataExchange) and nothing else. Everything in this directory is a real
capability outside those bounds. It works and it is tested. It is kept rather
than deleted, because deleting working code to make a point is a bad trade.

It is **not** part of the kernel's workspace:

- The root `Cargo.toml` has `exclude = ["outside"]`. Nothing here can be pulled
  back in by a path dependency unless that exclusion is deliberately removed.
- This directory has its own workspace and its own lockfile. Nothing here
  affects the kernel's dependency graph or what `cargo deny` checks.

```
cd outside && cargo test --workspace
```

## What is here, and why it is not kernel

| Crate | Why it is outside |
|---|---|
| `ogeom-sketch` | 2D geometric constraint solving. No conventional CAD kernel ships this; each application built on one supplies its own. Its module doc said so before it moved. |
| `ogeom-recognize` | Feature recognition (reading design intent back out of raw topology) and machining process planning. Recognition has no counterpart in the reference kernel; process planning is not geometry at all. |
| `ogeom-select` | BVH picking, marquee selection, and the draft and thickness analyses built on them. Selection belongs to the reference kernel's Visualization module, which is out of scope. |

## Maintenance

`tools/check.sh` does not compile this directory, and the kernel will change
underneath it.

- A non-gating CI job builds it, so breakage is visible without blocking the
  kernel on code the kernel has dropped.
- If that job has been red for a while, treat it as information, not an
  emergency. This is a snapshot of working code, not a maintained product.

To bring something back into the kernel, go through `docs/SCOPE.md` first:
argue the scope change in writing, then move the code. Mesh-to-b-rep with
canonical surface recognition came back this way: `docs/SCOPE.md` admits it,
and the kernel's `solid_from_mesh` replaced the crate that held it here.
