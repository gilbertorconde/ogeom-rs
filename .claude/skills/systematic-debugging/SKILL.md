---
name: systematic-debugging
description: "Use when encountering any bug, test failure, wrong shape, refusal or unexpected behaviour, before proposing a fix. Four-phase root-cause process for a geometry kernel: reproduce, trace to the source, one hypothesis at a time, a pinning test before the fix. For a red PR check loop use ci-iterate; for a measured slowness use performance-optimizer."
---

# Systematic debugging

Find the root cause before attempting a fix. A symptom fix in a kernel is
worse than no fix: the wrong answer moves to a shape nobody is looking at.

## The rule

No fix without a root-cause investigation first. If phase 1 is not done,
there is no fix to propose.

Use this for every technical issue: a failing test, a wrong volume, a
refusal that should not happen, a `check` failure, a panic, a stress case
that went `invalid` or `wrong`, a build failure. Use it especially when the
fix seems obvious, when a previous fix did not work, or when the issue is
not fully understood.

## Phase 1: root cause

Before any fix:

1. **Read the error completely.** The refusal text names the thing and the
   reason; the `check` diagnosis names the entity; a panic names the line.
   Note every number in it.
2. **Reproduce deterministically.** A test, a replay with
   `ogeom-stress --case <scenario>/<part>/<n>`, a corpus file through
   `ogeom-cli`, or a repro program in a scratch crate outside the repository
   that depends on the workspace by path. The repro must give the same
   answer every run; property tests record their failing seed in
   `*.proptest-regressions`, so read it.
3. **Check recent changes.** `git log` on the touched crates, a baseline
   drift in `tools/ogeom-stress/baseline.json`, a tolerance that moved.
4. **Gather evidence at each boundary.** The kernel is a pipeline, and the
   wrong value usually enters it several stages before it is noticed. Name
   the suspected boundary, then print what crosses it, with its tolerance:
   - the pcurve against its 3D curve, and the edge's representations
     against one another;
   - the section the marcher produced against both surfaces;
   - the 2D stage of the boolean (the face's trim) against the 3D pieces;
   - the pieces kept against the pieces classified;
   - the location chain at the point a result is built;
   - the file record read against the entity built from it, and back.
   Run once to see where it breaks, then investigate that stage.
5. **Trace the data flow backwards.** Where does the bad value originate?
   What called this with it? Keep going up until the source. Fix at the
   source, not where the symptom appears. See
   [root-cause-tracing.md](root-cause-tracing.md).

## Phase 2: pattern

1. Find a working example in the same codebase: the sibling case that does
   not fail (the cylinder where the cone fails, the planar where the spline
   fails, the right-handed where the mirrored fails).
2. Compare against the reference: the paper or specification the algorithm
   follows, read completely, not skimmed.
3. List every difference between working and broken, however small. A sign,
   a domain endpoint, a tolerance scale, a seam.
4. Name the dependencies: which tolerance, which parameterization, which
   invariant the code assumes.

## Phase 3: hypothesis

1. State one hypothesis: "X is the root cause because Y." Write it down.
2. Test it with the smallest change: one variable at a time.
3. Did it hold? Then phase 4. Did it not? Form a new hypothesis. Do not
   stack fixes.
4. When you do not know, say so, and gather more evidence instead of
   guessing.

## Phase 4: implementation

1. **Write the pinning test first.** The simplest reproduction, measured
   against a closed form, a published value or an independent computation.
   It must fail on the pre-fix code for the intended reason. Follow the
   authoring gate in the test-audit skill. A test that never failed proves
   nothing.
2. **One fix.** The root cause only. No "while I am here", no bundled
   refactor.
3. **Verify.** The pinning test passes, nothing else broke, the issue is
   gone in the original reproduction too. Run `cargo test -p <crate>` for
   the owner, then `tools/check.sh` before claiming success.
4. **If the fix does not work:** stop and count. Fewer than three attempts:
   return to phase 1 with the new information. Three or more: the pattern
   is wrong, not the hypothesis. Stop and discuss the design with the user
   before any further attempt.

## What a fix may not do

- Loosen a tolerance to make a test pass without the reason in the same
  commit.
- Turn a wrong answer into a silent one: an empty shape, a flag, a result
  that skips the hard case.
- Add a refusal without naming the thing and the reason, and without a
  test that pins the text. A refusal is a correct answer; a stub is not.
- Describe the change in a comment. Comments describe present behaviour.
  The investigation goes in the commit message.
- Name the file the user was testing with, in a comment or a test name. The
  test describes the geometric condition.

## Red flags

Stop and return to phase 1 when you catch yourself thinking: "quick fix for
now", "just try changing X", "add a few changes and run the tests", "skip the
test, I checked by hand", "it is probably X", "I do not fully understand but
this might work", "one more attempt" after two failures, or when each fix
reveals a new problem somewhere else.

The user's signals that the process has slipped: "is that not happening?"
(assumed, not verified), "will it show us...?" (no evidence gathered), "stop
guessing", "we're stuck?".

## When there is no root cause

If the investigation shows the issue is truly in the environment (a
toolchain, a loaded machine, an upstream crate), then document what was
investigated, handle it with a clear error or a named refusal, and never
with a retry or a silenced result. Most "no root cause" cases are an
incomplete investigation.

## Quick reference

| Phase | Activities | Done when |
|---|---|---|
| 1. Root cause | read, reproduce, check changes, probe the boundaries, trace back | what and why are understood |
| 2. Pattern | working sibling, reference text, list differences | the difference is named |
| 3. Hypothesis | one theory, smallest test | confirmed, or a new theory |
| 4. Implementation | pinning test, one fix, verify | the test passes and the suite is green |
