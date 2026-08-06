# outside

Code that was written here and is not part of the kernel.

`docs/SCOPE.md` sets the kernel's bounds: parity with the reference kernel's
modelling modules — FoundationClasses, ModelingData, ModelingAlgorithms,
DataExchange — and nothing else. Everything in this directory is a real
capability that falls outside those bounds. It works, it is tested, and it is
kept rather than deleted, because deleting working code to make a point is a
poor trade.

It is **not** part of the kernel's workspace. The root `Cargo.toml` carries
`exclude = ["outside"]`, so nothing here can be pulled back in by a path
dependency without that exclusion being removed deliberately. This directory has
its own workspace and its own lockfile; the kernel's dependency graph, and what
`cargo deny` reasons about, are unaffected by anything here.

```
cd outside && cargo test --workspace
```

## What is here, and why it is not kernel

| Crate | Why it is outside |
|---|---|
| `ogeom-sketch` | 2D geometric constraint solving. No conventional CAD kernel ships this; applications built on them each supply their own. Its own module doc said so before it moved. |
| `ogeom-recognize` | Feature recognition — reading design intent back out of raw topology — plus machining process planning. Recognition has no counterpart in the reference; process planning is not geometry at all. |
| `ogeom-select` | BVH picking, marquee selection, and the draft and thickness analyses that ride on them. Selection lives in the reference's Visualization module, which is out of scope. |
| `ogeom-reverse` | Mesh → b-rep with canonical surface recognition. Reverse engineering; the reference does not do it. |

## Rot

Nothing in `tools/check.sh` compiles this directory, and the kernel will move
underneath it. A non-gating CI job builds it so breakage is visible without
holding the kernel hostage to code the kernel disowned. If that job has been red
for a while, that is information rather than an emergency: this is a snapshot of
working code, not a maintained product.

If something here is ever wanted back inside, the way in is `docs/SCOPE.md` —
argue the scope change first, in writing, then move the code.
