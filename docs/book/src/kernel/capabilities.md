# Capabilities

The kernel's capability ledger, included verbatim from
[`docs/CAPABILITIES.md`](https://github.com/gilbertorconde/ogeom-rs/blob/main/docs/CAPABILITIES.md).
It is generated from the hand-written `docs/capabilities.toml`, and
`tools/check.sh` fails if it goes stale, so this chapter is as current as the
commit the book was built from.

Every capability has one of six verdicts, and the gate enforces the evidence
each requires:

| Verdict | Means | Must cite |
|---|---|---|
| `covered` | built and tested | symbols that resolve in the built rustdoc, and tests that exist in the tree |
| `partial` | built, with a stated restriction | the restriction in words, plus symbols and tests |
| `divergent` | deliberately done another way | the reasoning |
| `absent` | not built | a plan entry |
| `n/a` | not a kernel capability for us | the reasoning |
| `unreviewed` | not yet reviewed | counted against a ratchet that can only go down |

Rename a cited symbol or delete a cited test, and the build fails.

---

{{#include ../../../CAPABILITIES.md}}
