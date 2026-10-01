---
name: comment-cleanup
description: "Diff-scoped cleanup before review: removes narration, history, plan references and war stories from comments, keeps the ones that state what the code cannot, and strips slop (quieted warnings, swallowed results, one-use helpers, dual paths) from the same hunks. Use before committing, when the comment-rot lint fails, or when asked to deslop, trim or unify comments. A crate-wide pass is a campaign (below). Not a correctness review (use review) and not a test-portfolio audit (use test-audit)."
---

# Comment cleanup

Clean only the diff before review. Preserve behaviour absolutely.

## The rule

`CONTRIBUTING.md`: comments describe current behaviour, never the change
that produced it. No "used to", "formerly", "now returns", "since X
landed", and no issue, PR or commit as the reason for a behaviour; those
belong in the commit message. `lint-comment-rot: ignore` on a line opts it
out when a reference is genuinely needed.

This codebase explains itself in prose: module docs (`//!`) that state the
design, and comments that say why a thing is shaped as it is, are its
style. There is no percentage budget. The test is what a comment states,
not how many there are. What has gone wrong is not the amount but the
voice: some comments state a condition, others tell the story of the case
that found it, cite a plan that no longer exists, or pitch the design.

## Measure

```bash
node tools/lint-comment-rot.mjs                 # lines added against origin/main
node tools/lint-comment-rot.mjs --pedantic      # plus the advisory tier
node tools/lint-comment-rot.mjs --staged        # what the pre-commit hook sees
node tools/lint-comment-rot.mjs --all           # the whole tree, what CI runs
```

The gating tier names constructions that cannot describe present
behaviour. The advisory tier fires on correct comments often enough that
it is a review aid, not a rule: read each hit, do not delete on sight. The
patterns below that the lint does not catch (plan references, war stories,
pitch, kernel names) are found by reading.

## What a comment is for

A comment earns its place when it states what the code cannot:

- **the geometric condition** a branch handles, and what goes wrong
  without it: "a seam walked both ways bounds no region";
- **a tolerance and its reason**, which `CONTRIBUTING.md` requires next to
  every tolerance in a test and every factor in the code: "a hundred
  confusions: below any feature the pipeline resolves, above the slop a
  fitted section carries";
- **an approximation and its stated error**;
- **an invariant** of `docs/DATA_MODEL.md`, cited by section (`§5`), and
  what relies on it here;
- **a refusal's reason**: why this case is refused rather than attempted;
- **a side effect, an ordering, or a failure behaviour** a caller cannot
  see from the signature.

## Checklist

1. Scope to `git diff origin/main...HEAD` (or the range the user names).
   Never clean the whole tree in a feature commit.
2. For each added or changed comment, decide one action: KEEP, DELETE,
   UPDATE or MERGE. KEEP only when it states one of the things above.
3. Delete on sight:
   - narration of what the next line does, and syntax explanation;
   - history: "changed from", "used to", "now", "we", issue or commit
     numbers, dates, decision logs, paths not taken, "this crate found";
   - **plan references**: `§N` or a letter-number code that points at
     anything but `docs/DATA_MODEL.md` (the plan has no numbered sections),
     "milestone", "phase", "the deferred entry", "owed", "comes first
     because"; the plan is `docs/PLAN.md` and the ledger is
     `docs/parity/parity.toml`, and a comment cites neither as a reason;
   - **war stories**: the case that found the bug told as a story, with
     its numbers ("three descriptions arrived a tenth of a micron apart in
     turn, and the wire had two vertices where it needed one"). Rewrite as
     the condition and the rule it implies, present tense. The story goes
     in the commit message, and the case in a test named for the
     geometric condition;
   - **the user's file**: no name of a model or file someone tested with,
     in a comment or a test name; the geometric condition, instead;
   - **pitch**: claims about other software, "the insight worth
     preserving", "critically", "there is no X today", why the design is
     clever. A module doc states what the module does and the design it
     follows, not its case for existing;
   - **another kernel's name**, product or project, anywhere. The one
     sanctioned place for a reference *package* name is a module doc's
     `*Elsewhere:*` line, which maps the module to the reference packages
     it answers for, the way the parity ledger does; it never names a
     product, and nothing below the module doc names a reference class;
   - banners and dividers that restate the next symbol;
   - doc comments on small private functions that restate the name;
   - promises the code does not keep, and names or types restated as
     prose.
4. Keep each surviving comment timeless and direct, a why and not a what,
   in the voice of the file around it:
   - present tense, the condition first, then the consequence;
   - an inline comment is one paragraph at most; a longer argument moves
     to the function's doc comment or the module doc;
   - a number appears only when it is the contract (a factor, a bound), as
     the rule, never as the anecdote;
   - no em-dash character, no semicolons joining clauses.
5. In the same hunks, remove slop that is abnormal for the file:
   - an `#[allow(...)]` added to quiet a warning instead of fixing it, or
     without a `reason`;
   - `.ok()`, `unwrap_or_default()`, `let _ =` or `unwrap_or(0.0)` on a
     result that was a refusal: a failure is a value and is returned;
   - defensive checks for states the types already exclude;
   - a tolerance literal where `Tolerances` has the value;
   - `clone()`s that dodge a borrow a reference would satisfy;
   - one-use helpers and intermediate variables that add no domain
     meaning;
   - a `pub` widened only for a test (see test-audit);
   - aliases, re-exports, flags and dual paths kept for a caller that no
     longer exists;
   - style that conflicts with the surrounding file.
6. Make no functional edit. If a cleanup could change behaviour, leave it
   and report it.
7. Rerun the lint, `cargo fmt --all -- --check`, and the owning crate's
   clippy and tests.

## Campaign

A campaign cleans the comments of one crate in one commit. Read every
comment in the crate, give each an action and a one-line reason, apply the
ledger, comments only (no behaviour change, no renames, no reformatting),
and run `node tools/lint-comment-rot.mjs --all` and the crate's tests
unchanged. Have one independent reviewer read the deleted comments for a
lost contract (a tolerance reason, an invariant, a refusal's reason), and
restore only those with source evidence. Start with the crates whose
inline comments run longest: the boolean, the fillet corner, the mesher.

## Report

One to three sentences: what changed, the lint result before and after,
and any non-trivial item left for the author. Run this skill before
`review`, never in place of it.
