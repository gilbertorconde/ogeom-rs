#!/usr/bin/env python3
"""Changelog fragments: changelog.d/<issue>-<slug>.md, one entry each.

  changelog.py check              every fragment is one top-level bullet
  changelog.py preview            print the section a release would write
  changelog.py release <version> [--date YYYY-MM-DD]
                                  move the fragments into CHANGELOG.md
"""
import datetime
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FRAGMENTS = ROOT / "changelog.d"
CHANGELOG = ROOT / "CHANGELOG.md"
NAME = re.compile(r"^(\d+)-[a-z0-9][a-z0-9-]*\.md$")


def fragments():
    """The fragments, the highest issue number first, then by name."""
    found = []
    for path in FRAGMENTS.glob("*.md"):
        if path.name == "README.md":
            continue
        m = NAME.match(path.name)
        if not m:
            sys.exit(f"changelog.d/{path.name}: name it <issue>-<slug>.md")
        found.append((-int(m.group(1)), path.name, path))
    return [p for _, _, p in sorted(found)]


def entry(path):
    text = path.read_text().strip("\n")
    lines = text.splitlines()
    if not lines or not lines[0].startswith("- "):
        sys.exit(f"changelog.d/{path.name}: an entry starts with '- '")
    if any(l.startswith("- ") for l in lines[1:]):
        sys.exit(f"changelog.d/{path.name}: one top-level bullet per fragment")
    return text


def section():
    return "\n\n".join(entry(p) for p in fragments())


def release(version, date):
    body = section()
    old = CHANGELOG.read_text()
    marker = "## [Unreleased]\n"
    if marker not in old:
        sys.exit("CHANGELOG.md has no '## [Unreleased]' heading")
    head, rest = old.split(marker, 1)
    nxt = re.search(r"^## \[", rest, re.M)
    pending = (rest[: nxt.start()] if nxt else rest).strip("\n")
    tail = rest[nxt.start():] if nxt else ""
    parts = [p for p in (body, pending) if p]
    new = f"{head}{marker}\n## [{version}] - {date}\n\n" + "\n\n".join(parts) + "\n\n" + tail
    CHANGELOG.write_text(new)
    for p in fragments():
        p.unlink()


def main(argv):
    if len(argv) < 2 or argv[1] not in ("check", "preview", "release"):
        sys.exit(__doc__)
    if argv[1] == "check":
        section()
        print(f"changelog.d: {len(fragments())} fragment(s) ok")
    elif argv[1] == "preview":
        print(section())
    else:
        if len(argv) < 3:
            sys.exit(__doc__)
        date = datetime.date.today().isoformat()
        if "--date" in argv:
            date = argv[argv.index("--date") + 1]
        release(argv[2], date)


if __name__ == "__main__":
    main(sys.argv)
