#!/usr/bin/env python3
"""The facade gate: `crates/ogeom` names every public item of every crate root.

  facade.py check    fail when a crate root exports a name the facade omits
  facade.py list     print each crate's root names, one `crate name` per line

`crates/ogeom/src/lib.rs` re-exports items by name, so the `ogeom::<module>`
paths stay put when an item moves between crates. A name the facade never
lists is unreachable for users, so a new public item has to be added there,
or listed in INTERNAL below as deliberately kept out.

Read from the source, not from rustdoc, so the gate needs no nightly
toolchain. It reads only what rustfmt'd crate roots contain: `pub` items and
`pub use` trees at the top level of `src/lib.rs`, and `#[macro_export]`
macros anywhere in the crate. Items under `cfg(test)` or the `test-hooks`
feature are not part of the API and are skipped.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FACADE = ROOT / "crates" / "ogeom" / "src" / "lib.rs"

# Root names the facade leaves out on purpose, as `crate::name`.
INTERNAL: set[str] = {
    # Reached as `ogeom::boolean::POLE_TOUCHING`, and the booleans' own step.
    "ogeom_algo::POLE_TOUCHING",
    "ogeom_algo::operands_clear_of_poles",
    # The blend engine's entry points into the N-sided filling.
    "ogeom_algo::FillRegion",
    "ogeom_algo::make_filling_n_fixed",
    "ogeom_algo::make_filling_n_region",
    "ogeom_algo::make_filling_n_within",
    # The upgrade's steps the healing and the many-operand union call.
    "ogeom_algo::merge_edges_at",
    "ogeom_algo::measure_edge",
}

LEXEME = re.compile(
    r"""//[^\n]*"""
    r"""|/\*.*?\*/"""
    r"""|r(?P<hashes>\#*)".*?"(?P=hashes)"""
    r"""|b?"(?:\\.|[^"\\])*\""""
    r"""|b?'(?:\\.|[^\\'\n])'""",
    re.S,
)
ITEM = re.compile(
    r"\bpub\s+(?:(?:unsafe|async|const|extern)\s+)*fn\s+(\w+)"
    r"|\bpub\s+(?:struct|enum|trait|type|const|static|mod|union)\s+(\w+)"
    r"|\bpub\s+use\s+([^;]*);"
)
SKIP_CFG = re.compile(r'#\[cfg\((?:test|feature\s*=\s*"test-hooks")\)\]')
MACRO = re.compile(r"#\[macro_export\][^!]*?macro_rules!\s*(\w+)", re.S)


def strip(text: str) -> str:
    """Comments and literals blanked in place, so braces and keywords inside
    them vanish and every offset still points at the same source character."""

    def blank(m: re.Match) -> str:
        return re.sub(r"[^\n]", " ", m.group(0))

    return LEXEME.sub(blank, text)


def depths(text: str) -> list[int]:
    out, d = [], 0
    for ch in text:
        if ch == "}":
            d -= 1
        out.append(d)
        if ch == "{":
            d += 1
    return out


def use_leaves(tree: str, prefix: list[str]) -> list[tuple[list[str], str]]:
    """(path, bound name) for each leaf of a use tree."""
    tree = tree.strip()
    leaves = []
    m = re.match(r"^(.*?)::\s*\{(.*)\}$", tree, re.S)
    if m:
        base = prefix + [s for s in m.group(1).split("::") if s.strip()]
        base = [s.strip() for s in base]
        depth, start, body = 0, 0, m.group(2)
        parts = []
        for i, ch in enumerate(body):
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
            elif ch == "," and depth == 0:
                parts.append(body[start:i])
                start = i + 1
        parts.append(body[start:])
        for p in parts:
            if p.strip():
                leaves += use_leaves(p, base)
        return leaves
    if tree.startswith("{") and tree.endswith("}"):
        return use_leaves("self::" + tree, prefix[:])
    alias = None
    am = re.match(r"^(.*?)\s+as\s+(\w+)$", tree, re.S)
    if am:
        tree, alias = am.group(1).strip(), am.group(2)
    path = prefix + [s.strip() for s in tree.split("::")]
    if path[-1] == "*":
        sys.exit(f"glob re-export `{'::'.join(path)}` at a crate root: name the items")
    if path[-1] == "self":
        path = path[:-1]
    if alias == "_":
        return []
    return [(path, alias or path[-1])]


def top_level(source: str):
    """(kind, payload) for each `pub` item at the top level of a file."""
    text = strip(source)
    d = depths(text)
    for m in ITEM.finditer(text):
        if d[m.start()] != 0:
            continue
        head = text[: m.start()]
        cut = max(head.rfind(";"), head.rfind("}"))
        if SKIP_CFG.search(source[cut + 1 : m.start()]):
            continue
        if m.group(3) is not None:
            yield "use", m.group(3)
        else:
            yield "item", m.group(1) or m.group(2)


def crates() -> dict[str, pathlib.Path]:
    found = {}
    for toml in sorted((ROOT / "crates").glob("*/Cargo.toml")):
        name = re.search(r'^name\s*=\s*"([^"]+)"', toml.read_text(), re.M).group(1)
        if name != "ogeom":
            found[name.replace("-", "_")] = toml.parent / "src"
    return found


def root_names(src: pathlib.Path) -> set[str]:
    names = set()
    for kind, payload in top_level((src / "lib.rs").read_text()):
        if kind == "item":
            names.add(payload)
        else:
            names.update(name for _, name in use_leaves(payload, []))
    for path in src.rglob("*.rs"):
        names.update(MACRO.findall(strip(path.read_text())))
    return names


def facade_names() -> dict[str, set[str]]:
    """Names the facade re-exports straight from a crate root, per crate."""
    named: dict[str, set[str]] = {}
    for m in re.finditer(r"\bpub\s+use\s+([^;]*);", strip(FACADE.read_text())):
        for path, _ in use_leaves(m.group(1), []):
            if len(path) == 2 and path[0].startswith("ogeom_"):
                named.setdefault(path[0], set()).add(path[1])
    return named


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    roots = {c: root_names(src) for c, src in crates().items()}
    if cmd == "list":
        for c, names in roots.items():
            for n in sorted(names):
                print(c, n)
        return
    if cmd != "check":
        sys.exit(__doc__)
    named = facade_names()
    missing = sorted(
        f"{c}::{n}"
        for c, names in roots.items()
        for n in names - named.get(c, set())
        if f"{c}::{n}" not in INTERNAL
    )
    if missing:
        print(f"{FACADE.relative_to(ROOT)} does not re-export:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        print("Add each to its facade module, or to INTERNAL in tools/facade.py.", file=sys.stderr)
        sys.exit(1)
    total = sum(len(n) for n in roots.values())
    print(f"facade: {total} root names across {len(roots)} crates, all re-exported")


if __name__ == "__main__":
    main()
