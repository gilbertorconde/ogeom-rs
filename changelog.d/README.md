# Changelog fragments

Each pull request that changes behaviour adds its `CHANGELOG.md` entry here
as a file of its own, instead of editing `CHANGELOG.md`. Two pull requests
then never touch the same lines, and neither has to be rebased for the
other's entry.

- **One file per change**, named after its issue: `1234-short-slug.md`
  (a change with no issue: `0-short-slug.md`).
- **The file is the entry**, exactly as it will read under the version's
  heading: one top-level `- **Added: …**` / `- **Fixed: …**` /
  `- **Changed: …**` bullet with its sub-bullets, the issue number in the
  bold line as usual.
- A follow-up to a change not yet released edits that change's fragment.

At a release, `tools/changelog.py release <version>` moves every fragment
into `CHANGELOG.md` under `## [<version>] - <date>`, the highest issue number
first, and deletes the fragments. `tools/changelog.py preview` prints the
section it would write. `tools/changelog.py check` fails on a fragment that
is not a single top-level bullet.
