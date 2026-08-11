#!/usr/bin/env python3
"""apisurf — read a reference CAD kernel's shape, and how applications use it.

**Analysis tool. Not part of the build.** `scope` derives the parity index that
`docs/SCOPE.md` defines; `profile` measures usage, which sequences work and
never sets its bounds.

Given a reference kernel's headers and an application that consumes them, this
reports which types and members the application actually touches, how often, and
from which of its modules. The output is a *usage profile*.

What it is good for: **sequencing**. It shows which parts of a kernel sit on the
hot paths of real work, and therefore which deserve the earliest design
attention, the earliest benchmarks, and the most careful API shape.

What it is emphatically not good for: deciding what to build. Reference counts
understate anything reached through a narrow facade — boolean operations show a
few dozen references and represent perhaps a fifth of a kernel's implementation
effort — and they say nothing at all about capabilities the sampled application
implements for itself, such as constraint solving. See `docs/SCOPE.md`.

Neither input is a dependency of this repository. Point it at whatever you have
locally; nothing is vendored and nothing is committed.

Two passes:

1. **Index the reference.** Collect every header the consumer includes, emit one
   translation unit including all of them, and parse it once with libclang. One
   parse rather than several hundred keeps this to seconds. Records classes,
   methods with signatures, enums and type aliases.

2. **Scan the consumer.** Regex over its sources for includes, type references
   and member calls. Without a compile database this is textual rather than
   semantic, so member names owned by more than one class are flagged ambiguous
   and attributed to each candidate; the unambiguous counts are the ones to
   trust.

Two subcommands:

``profile`` (the original behaviour) — needs libclang and both checkouts::

    python3 tools/apisurf/apisurf.py profile \\
        --reference /path/to/kernel/src \\
        --consumer  /path/to/application \\
        --out docs/api_surface.json

``scope`` — needs only the reference checkout, no libclang. Walks the
reference's own module and toolkit manifests, applies the triage rules from
`docs/SCOPE.md`, and emits the parity index every audit row is checked
against::

    python3 tools/apisurf/apisurf.py scope \\
        --reference vendor/OCCT \\
        --out docs/parity/reference-index.tsv
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import subprocess
import sys
import tempfile

# libclang is imported lazily inside cmd_profile: the scope subcommand is
# deliberately runnable on a machine with nothing but Python and a checkout.
ci = None


def _need_libclang():
    global ci
    try:
        import clang.cindex as _ci
    except ImportError:
        sys.exit("apisurf profile needs python-clang: pacman -S python-clang / pip install clang")
    ci = _ci


# Third-party headers vendored inside the consumer's own tree. Excluding them
# keeps the profile honest: code the application bundles is not code the kernel
# is being asked to provide.
VENDORED_PREFIXES = (
    "SMESH", "SMDS", "SMESHDS", "StdMeshers", "NETGENPlugin", "MED", "Driver",
    "Utils", "Basics", "SALOME", "Jt", "OCC_", "netgen",
)

# Directories under the consumer's src/ that hold vendored third-party code.
SKIP_DIRS = ("3rdParty", "zipios++", "CXX")

INCLUDE_RE = re.compile(r'#\s*include\s*<([A-Za-z][A-Za-z0-9_]*)\.hxx>')
SOURCE_EXTS = (".cpp", ".cxx", ".cc", ".h", ".hpp", ".hxx", ".pyi", ".inl")


def clang_resource_dir() -> str | None:
    for exe in ("clang++", "clang"):
        try:
            out = subprocess.run(
                [exe, "-print-resource-dir"], capture_output=True, text=True, check=True
            )
            return out.stdout.strip()
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue
    return None


def find_libclang() -> str | None:
    import glob
    for pat in ("/usr/lib/libclang.so*", "/usr/lib64/libclang.so*",
                "/usr/lib/llvm*/lib/libclang.so*"):
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[-1]
    return None


# ==========================================================================
# scope — derive the parity index from the reference tree's own manifests
# ==========================================================================
#
# The primary key of the audit is OUR capabilities, not these headers; this
# index is the completeness check those capabilities are balanced against.
# Every rule here has a stable id, materialised on every row it touches, so a
# rule change shows up in review as a diff naming exactly which headers moved.

IN_SCOPE_MODULES = ("FoundationClasses", "ModelingData", "ModelingAlgorithms",
                    "DataExchange")

# --- package-level rules --------------------------------------------------

# C-STEP / C-IGES: one class per EXPRESS entity (plus its RW twin). Auditing
# them one by one is absurd and dropping them loses the coverage number the
# exchange layer wants, so they collapse to per-package rows whose coverage
# figure is regenerated from the entity tables in ogeom-io. The *translator*
# packages (StepToGeom, GeomToStep, IGESToBRep, ...) and the file models
# (StepData, StepFile, IGESData, IGESFile) are capabilities and stay in the
# keep pool.
COLLAPSE_STEP = re.compile(
    r"^(RWStep[A-Z]\w*|Step(Basic|Geom|Shape|Repr|Visual|DimTol|Element|FEA|"
    r"Kinematics|AP203|AP209|AP214|AP242))$")
COLLAPSE_IGES = re.compile(
    r"^IGES(Basic|Geom|Dimen|Draw|Appli|Defs|Graph|Solid)$")

# R2: generic containers. Rust's generics provide every one of these for free;
# there is no capability to audit, only monomorphization spelled out by hand.
R2_CONTAINER_PKGS = frozenset({
    "NCollection", "TColStd", "TColgp", "TColGeom", "TColGeom2d",
    "TColQuantity", "TShort", "TCollection",
})

# R3: operating-system and infrastructure surface. Memory, threads, streams,
# strings, resources, persistence plumbing — the standard library's job, not a
# geometry kernel's. The carve-out is cancellable progress, which is a kernel
# API shape (ogeom-core's Watch) and not an OS abstraction.
R3_INFRA_PKGS = frozenset({
    "Standard", "OSD", "Plugin", "Resource", "FSD", "Storage", "StdFail",
    "Message", "FlexLexer",
})
R3_CARVE_OUT = re.compile(r"^Message_Progress")

# R4: the TopOpeBRep* boolean family, superseded *inside the reference* by the
# BOPAlgo pipeline. Parity with a kernel does not include parity with the
# parts it keeps only for compatibility.
R4_SUPERSEDED = re.compile(r"^TopOpeBRep")

# R5: OCAF persistence drivers that ride along inside DataExchange. The
# exchange *document* (XCAFDoc) is in scope; serializing it through the
# application framework's Bin/Xml driver stack is the framework's concern.
R5_OCAF_DRIVERS = frozenset({
    "BinXCAFDrivers", "BinMXCAFDoc", "XmlXCAFDrivers", "XmlMXCAFDoc",
})

# --- header-level rule ----------------------------------------------------

# R1: CDL generic instantiations. The reference predates C++ templates in its
# own dialect; every instantiation got a header whose name chains the
# arguments together with "Of". Three shapes give it away: a container word
# before the Of, a The-/My- prefixed segment (the dialect's naming for an
# algorithm's own internals), or two Ofs chained. What survives is protected
# by SEMANTIC_OF: names where "...Of..." is English — a surface of revolution,
# an arc of a circle, a degree of freedom — are capabilities, not
# instantiations.
R1_OF = re.compile(r"Of[A-Z]")
R1_CONTAINER_WORD = re.compile(
    r"(Array[12]?|List|Map|Sequence|Set|Queue|Stack|Iterator|Tree|Buffer|Row"
    r"|Column|H?Seq|TSeq|Vector|GlobalNode|Node|Actor|Binder|Hasher)Of[A-Z]")
R1_THE_MY = re.compile(r"(?:^|_)(?:The|My)[A-Z]|[a-z0-9](?:The|My)[A-Z]")
SEMANTIC_OF = re.compile(
    r"(Surface|Arc|Type|Name|Out|Degree|Pair|Couple|Distribution|Energy"
    r"|Elements|Representation|Centre|Center|Week|Coefficient|Structure"
    r"|Selector|Axis|Moment)Of[A-Z]")

KIND_RES = (
    ("enum", re.compile(r"\benum\s+(?:class\s+)?{}\b")),
    ("class", re.compile(r"\b(?:class|struct)\s+{}\b")),
    ("alias", re.compile(r"(?:\btypedef\b[^;]*\b{}\s*;|\busing\s+{}\s*=)")),
)


def scope_rows(reference: str) -> list[dict]:
    """One row per public header of the four in-scope modules."""
    modules = {}
    with open(os.path.join(reference, "adm", "MODULES")) as fh:
        for line in fh:
            parts = line.split()
            if parts:
                modules[parts[0]] = parts[1:]
    for m in IN_SCOPE_MODULES:
        if m not in modules:
            sys.exit(f"module {m} not in {reference}/adm/MODULES")

    rows = []
    for module in IN_SCOPE_MODULES:
        for tk in modules[module]:
            pkg_file = os.path.join(reference, "src", tk, "PACKAGES")
            if not os.path.exists(pkg_file):
                sys.exit(f"toolkit {tk} has no PACKAGES file")
            for pkg in open(pkg_file).read().split():
                pkgdir = os.path.join(reference, "src", pkg)
                if not os.path.isdir(pkgdir):
                    continue
                for f in sorted(os.listdir(pkgdir)):
                    if not f.endswith(".hxx"):
                        continue
                    stem = f[:-4]
                    text = open(os.path.join(pkgdir, f), errors="replace").read()
                    kind = "other"
                    for k, pat in KIND_RES:
                        if re.search(pat.pattern.format(stem, stem), text):
                            kind = k
                            break
                    rows.append({
                        "module": module, "toolkit": tk, "package": pkg,
                        "header": stem, "kind": kind,
                        "deprecated": int("Standard_DEPRECATED" in text),
                    })
    return rows


def triage(rows: list[dict]) -> None:
    """Assign bucket and rule in place. Order matters: package rules first, so
    the exchange collapse removes its semantic-Of entity names (SurfaceOf
    Revolution and friends) before R1 ever sees them."""
    for r in rows:
        pkg, stem = r["package"], r["header"]
        classpart = stem.split("_", 1)[-1]
        if COLLAPSE_STEP.match(pkg):
            r["bucket"], r["rule"] = "collapse", "C-STEP"
        elif COLLAPSE_IGES.match(pkg):
            r["bucket"], r["rule"] = "collapse", "C-IGES"
        elif pkg in R2_CONTAINER_PKGS:
            r["bucket"], r["rule"] = "drop", "R2"
        elif pkg in R3_INFRA_PKGS:
            if R3_CARVE_OUT.match(stem):
                r["bucket"], r["rule"] = "keep", "V-PROGRESS"
            else:
                r["bucket"], r["rule"] = "drop", "R3"
        elif R4_SUPERSEDED.match(pkg):
            r["bucket"], r["rule"] = "drop", "R4"
        elif pkg in R5_OCAF_DRIVERS:
            r["bucket"], r["rule"] = "drop", "R5"
        elif R1_THE_MY.search(classpart) or (
                R1_OF.search(classpart) and (
                    R1_CONTAINER_WORD.search(stem)
                    or len(R1_OF.findall(classpart)) > 1
                    or not SEMANTIC_OF.search(stem))):
            r["bucket"], r["rule"] = "drop", "R1"
        else:
            r["bucket"], r["rule"] = "keep", ""

    # Safety valve: a package the rules emptied still gets one keep row — a
    # package that is nothing but instantiations of one algorithm (AppDef,
    # BRepApprox) is still a capability to audit, and silence here would hide
    # it. The namesake header is re-kept where there is one, else the first.
    by_pkg = collections.defaultdict(list)
    for r in rows:
        by_pkg[r["package"]].append(r)
    for pkg, prs in by_pkg.items():
        if any(x["bucket"] != "drop" for x in prs):
            continue
        if all(x["rule"] == "R1" for x in prs):
            namesake = next((x for x in prs if x["header"] == pkg), prs[0])
            namesake["bucket"], namesake["rule"] = "keep", "V-EMPTIED"


def cmd_scope(args) -> int:
    rows = scope_rows(args.reference)
    triage(rows)

    usage = {}
    if os.path.exists(args.profile):
        with open(args.profile) as fh:
            usage = json.load(fh).get("include_counts", {})
    elif args.profile != SCOPE_DEFAULT_PROFILE:
        sys.exit(f"profile {args.profile} not found")
    for r in rows:
        r["usage"] = usage.get(r["header"], 0)

    rows.sort(key=lambda r: (r["module"], r["package"], r["header"]))
    cols = ("module", "toolkit", "package", "header", "kind", "deprecated",
            "bucket", "rule", "usage")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    buckets = collections.Counter(r["bucket"] for r in rows)
    rules = collections.Counter(r["rule"] for r in rows if r["rule"])
    packages = {r["package"] for r in rows}
    print(f"wrote {args.out}")
    print(f"  {len(packages)} packages, {len(rows)} headers")
    print(f"  keep {buckets['keep']}  collapse {buckets['collapse']}  "
          f"drop {buckets['drop']}")
    for rule, n in sorted(rules.items()):
        print(f"    {rule:<11} {n}")

    if args.sample:
        # Deterministic stratified sample of the drops, for the hand read that
        # the rules must survive before they are trusted.
        import random
        rng = random.Random(0)
        dropped = collections.defaultdict(list)
        for r in rows:
            if r["bucket"] == "drop":
                dropped[r["rule"]].append(r)
        total = sum(len(v) for v in dropped.values())
        print(f"\nsample of {args.sample} dropped headers "
              f"(deterministic, stratified over {total}):")
        for rule in sorted(dropped):
            share = max(1, round(args.sample * len(dropped[rule]) / total))
            for r in rng.sample(dropped[rule], min(share, len(dropped[rule]))):
                print(f"  {rule:<4} {r['package']:<20} {r['header']}")
    return 0


SCOPE_DEFAULT_PROFILE = "docs/api_surface.json"


# ==========================================================================
# profile — how one application exercises the reference (original behaviour)
# ==========================================================================

# --------------------------------------------------------------------------
# Pass 1 — what the consumer includes
# --------------------------------------------------------------------------

def consumer_sources(consumer: str):
    src = os.path.join(consumer, "src")
    for root, dirs, files in os.walk(src):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for f in files:
            if f.endswith(SOURCE_EXTS):
                yield os.path.join(root, f)


def reference_header_index(reference: str) -> dict[str, tuple[str, str]]:
    """header stem -> (package, absolute path)"""
    index = {}
    src = os.path.join(reference, "src")
    for pkg in sorted(os.listdir(src)):
        pkgdir = os.path.join(src, pkg)
        if not os.path.isdir(pkgdir):
            continue
        for f in os.listdir(pkgdir):
            if f.endswith(".hxx"):
                index[f[:-4]] = (pkg, os.path.join(pkgdir, f))
    return index


def collect_includes(consumer: str, ref_headers: dict) -> tuple[dict, dict, dict]:
    """Returns (header -> count, header -> {module: count}, file -> [headers])."""
    counts = collections.Counter()
    per_module = collections.defaultdict(collections.Counter)
    files = {}
    src_root = os.path.join(consumer, "src")

    for path in consumer_sources(consumer):
        try:
            text = open(path, "r", errors="replace").read()
        except OSError:
            continue
        hits = [h for h in INCLUDE_RE.findall(text)
                if h in ref_headers and not h.startswith(VENDORED_PREFIXES)]
        if not hits:
            continue
        rel = os.path.relpath(path, src_root)
        module = "/".join(rel.split(os.sep)[:3]) if rel.startswith("Mod") else rel.split(os.sep)[0]
        files[rel] = sorted(set(hits))
        for h in hits:
            counts[h] += 1
            per_module[h][module] += 1
    return counts, per_module, files


# --------------------------------------------------------------------------
# Pass 2 — parse those headers, once
# --------------------------------------------------------------------------

def parse_reference(reference: str, headers: list[str], ref_headers: dict) -> dict:
    import glob

    incs = ["-I" + d for d in sorted(glob.glob(os.path.join(reference, "src", "*")))]
    args = ["-x", "c++", "-std=c++17", "-fsyntax-only", "-DHAVE_CONFIG_H=0"] + incs
    res = clang_resource_dir()
    if res:
        args += ["-isystem", os.path.join(res, "include")]

    with tempfile.NamedTemporaryFile("w", suffix=".cpp", delete=False) as tmp:
        for h in sorted(headers):
            tmp.write(f"#include <{h}.hxx>\n")
        unit = tmp.name

    try:
        idx = ci.Index.create()
        tu = idx.parse(unit, args=args, options=ci.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES)
        fatal = [d for d in tu.diagnostics if d.severity >= 4]
        errors = [d for d in tu.diagnostics if d.severity == 3]
        if fatal:
            print(f"  fatal: {fatal[0].spelling}", file=sys.stderr)
        if errors:
            print(f"  {len(errors)} parse errors (continuing; first: {errors[0].spelling})",
                  file=sys.stderr)

        classes: dict[str, dict] = {}
        enums: dict[str, dict] = {}
        aliases: list[tuple[str, str, str, object]] = []

        def is_kernel_name(name: str, pkg: str) -> bool:
            """Many kernels use a strict ``Pkg_Name`` convention, with the
            package's own static-utility class named plainly ``Pkg``.

            Enforcing it filters out nested helpers with generic names —
            ``Hasher``, ``Mesh``, ``Iterator`` — that would otherwise collide with
            unrelated identifiers in the consumer and inflate counts badly.
            """
            return bool(name) and (name == pkg or name.startswith(pkg + "_"))

        def members_of(cur) -> tuple[dict, dict, list]:
            members, statics, ctors = {}, {}, []
            for ch in cur.get_children():
                if ch.access_specifier != ci.AccessSpecifier.PUBLIC:
                    continue
                if ch.kind == ci.CursorKind.CXX_METHOD:
                    target = statics if ch.is_static_method() else members
                    target.setdefault(ch.spelling, []).append(ch.displayname)
                elif ch.kind == ci.CursorKind.CONSTRUCTOR:
                    ctors.append(ch.displayname)
            return members, statics, ctors

        for cur in tu.cursor.walk_preorder():
            if not cur.location.file:
                continue
            fname = os.path.basename(cur.location.file.name)
            if not fname.endswith(".hxx"):
                continue
            stem = fname[:-4]
            pkg = ref_headers.get(stem, ("?", ""))[0]

            if cur.kind in (ci.CursorKind.CLASS_DECL, ci.CursorKind.STRUCT_DECL) and cur.is_definition():
                name = cur.spelling
                if not is_kernel_name(name, pkg) or name in classes:
                    continue
                members, statics, ctors = members_of(cur)
                classes[name] = {
                    "package": pkg,
                    "header": stem,
                    "kind": "class",
                    "bases": [b.type.spelling for b in cur.get_children()
                              if b.kind == ci.CursorKind.CXX_BASE_SPECIFIER],
                    "constructors": ctors,
                    "methods": {k: sorted(set(v)) for k, v in sorted(members.items())},
                    "static_methods": {k: sorted(set(v)) for k, v in sorted(statics.items())},
                }

            elif cur.kind in (ci.CursorKind.TYPEDEF_DECL, ci.CursorKind.TYPE_ALIAS_DECL):
                # Shape maps, lists and coordinate arrays are frequently typedefs
                # of container templates rather than classes. Missing them would
                # drop some of the most heavily used types in the profile.
                name = cur.spelling
                if is_kernel_name(name, pkg) and name not in classes:
                    aliases.append((name, pkg, stem, cur))

        for name, pkg, stem, cur in aliases:
            if name in classes:
                continue
            under = cur.underlying_typedef_type
            decl = under.get_declaration()
            members, statics, ctors = members_of(decl) if decl.is_definition() else ({}, {}, [])
            classes[name] = {
                "package": pkg,
                "header": stem,
                "kind": "alias",
                "aliases": under.spelling,
                "bases": [b.type.spelling for b in decl.get_children()
                          if b.kind == ci.CursorKind.CXX_BASE_SPECIFIER] if decl else [],
                "constructors": ctors,
                "methods": {k: sorted(set(v)) for k, v in sorted(members.items())},
                "static_methods": {k: sorted(set(v)) for k, v in sorted(statics.items())},
            }

        for cur in tu.cursor.walk_preorder():
            if cur.kind == ci.CursorKind.ENUM_DECL and cur.is_definition() and cur.location.file:
                fname = os.path.basename(cur.location.file.name)
                if not fname.endswith(".hxx"):
                    continue
                stem = fname[:-4]
                pkg = ref_headers.get(stem, ("?", ""))[0]
                if is_kernel_name(cur.spelling, pkg):
                    enums.setdefault(cur.spelling, {
                        "package": pkg,
                        "header": stem,
                        "values": [c.spelling for c in cur.get_children()],
                    })

        return {"classes": classes, "enums": enums}
    finally:
        os.unlink(unit)


# --------------------------------------------------------------------------
# Pass 3 — which members does the consumer call?
# --------------------------------------------------------------------------

def scan_usage(consumer: str, index: dict) -> tuple[dict, dict, dict]:
    classes = index["classes"]
    enums = index["enums"]

    # member name -> set of owning classes, for ambiguity accounting
    owners = collections.defaultdict(set)
    for cname, c in classes.items():
        for m in list(c["methods"]) + list(c["static_methods"]):
            owners[m].add(cname)

    all_types = set(classes) | set(enums)
    enum_values = {v: e for e, d in enums.items() for v in d["values"]}

    type_refs = collections.Counter()
    member_calls = collections.Counter()
    enum_refs = collections.Counter()
    qualified = collections.Counter()   # Class::member — unambiguous by construction

    token_re = re.compile(r'\b([A-Za-z_][A-Za-z0-9_]*)\b')
    qual_re = re.compile(r'\b([A-Za-z_][A-Za-z0-9_]*)\s*::\s*([A-Za-z_][A-Za-z0-9_]*)')
    call_re = re.compile(r'(?:\.|->)\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(')

    for path in consumer_sources(consumer):
        try:
            text = open(path, "r", errors="replace").read()
        except OSError:
            continue
        # cheap prefilter: skip files that include nothing from the reference
        if ".hxx>" not in text:
            continue

        for tok in token_re.findall(text):
            if tok in all_types:
                type_refs[tok] += 1
            elif tok in enum_values:
                enum_refs[tok] += 1

        for cls, mem in qual_re.findall(text):
            if cls in classes and (mem in classes[cls]["methods"]
                                   or mem in classes[cls]["static_methods"]):
                qualified[(cls, mem)] += 1

        for mem in call_re.findall(text):
            if mem in owners:
                member_calls[mem] += 1

    return type_refs, member_calls, {"enum_refs": enum_refs, "qualified": qualified,
                                     "owners": owners}


def cmd_profile(args) -> int:
    _need_libclang()
    lib = args.libclang or find_libclang()
    if lib:
        ci.Config.set_library_file(lib)

    print("indexing reference headers ...")
    ref_headers = reference_header_index(args.reference)
    print(f"  {len(ref_headers)} headers in {args.reference}/src")

    print("scanning consumer includes ...")
    inc_counts, per_module, per_file = collect_includes(args.consumer, ref_headers)
    print(f"  {len(inc_counts)} distinct reference headers included, "
          f"{sum(inc_counts.values())} include lines, {len(per_file)} files")

    print("parsing reference (one translation unit) ...")
    index = parse_reference(args.reference, list(inc_counts), ref_headers)
    print(f"  {len(index['classes'])} classes, {len(index['enums'])} enums reachable")

    print("scanning consumer usage ...")
    type_refs, member_calls, extra = scan_usage(args.consumer, index)
    owners = extra["owners"]
    qualified = extra["qualified"]
    print(f"  {len(type_refs)} distinct types referenced")

    # ---- assemble ----
    all_classes = index["classes"]

    def ancestry(cname: str, seen: set | None = None) -> list[str]:
        """Transitive base classes. Algorithm classes typically declare almost
        nothing of their own, inheriting their whole useful surface from a base
        several levels up, so a per-class member list alone badly understates
        what a type actually offers."""
        seen = seen if seen is not None else set()
        out = []
        for b in all_classes.get(cname, {}).get("bases", []):
            b = b.split("<")[0].strip().removeprefix("class ").removeprefix("struct ")
            if b in all_classes and b not in seen:
                seen.add(b)
                out.append(b)
                out.extend(ancestry(b, seen))
        return out

    packages = collections.defaultdict(lambda: {"headers": 0, "includes": 0,
                                                "classes": [], "type_refs": 0})
    out_classes = {}
    for cname, c in sorted(all_classes.items()):
        refs = type_refs.get(cname, 0)
        if refs == 0 and inc_counts.get(c["header"], 0) == 0:
            continue
        used, unused = {}, []
        for m, sigs in list(c["methods"].items()) + list(c["static_methods"].items()):
            n = member_calls.get(m, 0) + qualified.get((cname, m), 0)
            if n:
                used[m] = {
                    "calls": n,
                    "ambiguous": len(owners[m]) > 1,
                    "owners": len(owners[m]),
                    "signatures": sigs,
                }
            else:
                unused.append(m)

        bases = ancestry(cname)
        inherited = {}
        for b in bases:
            bc = all_classes[b]
            for m in list(bc["methods"]) + list(bc["static_methods"]):
                n = member_calls.get(m, 0) + qualified.get((cname, m), 0)
                if n and m not in used:
                    inherited.setdefault(m, {"calls": n, "from": b})

        entry = {
            "package": c["package"],
            "header": c["header"],
            "kind": c.get("kind", "class"),
            "bases": c["bases"],
            "ancestry": bases,
            "includes": inc_counts.get(c["header"], 0),
            "type_refs": refs,
            "modules": dict(per_module.get(c["header"], {})),
            "constructors": c["constructors"],
            "used_members": dict(sorted(used.items(), key=lambda kv: -kv[1]["calls"])),
            "inherited_used_members": dict(sorted(inherited.items(),
                                                  key=lambda kv: -kv[1]["calls"])),
            "unused_members": sorted(unused),
            "coverage": f"{len(used)}/{len(used) + len(unused)}",
        }
        if "aliases" in c:
            entry["aliases"] = c["aliases"]
        out_classes[cname] = entry
        p = packages[c["package"]]
        p["classes"].append(cname)
        p["type_refs"] += refs
        p["includes"] += inc_counts.get(c["header"], 0)

    for pkg in packages:
        packages[pkg]["headers"] = len(packages[pkg]["classes"])
        packages[pkg]["classes"].sort()

    doc = {
        "_comment": "Usage profile from tools/apisurf/apisurf.py. Sequencing input "
                    "only, NOT a scope specification. See docs/SCOPE.md.",
        "reference_ref": subprocess.run(["git", "-C", args.reference, "describe", "--tags"],
                                        capture_output=True, text=True).stdout.strip(),
        "consumer_ref": subprocess.run(["git", "-C", args.consumer, "rev-parse", "--short", "HEAD"],
                                       capture_output=True, text=True).stdout.strip(),
        "totals": {
            "headers_included": len(inc_counts),
            "include_lines": sum(inc_counts.values()),
            "files_using_reference": len(per_file),
            "classes": len(out_classes),
            "enums": len(index["enums"]),
        },
        "packages": dict(sorted(packages.items(),
                                key=lambda kv: -kv[1]["type_refs"])),
        "classes": out_classes,
        "enums": {k: v for k, v in sorted(index["enums"].items())},
        "include_counts": dict(inc_counts.most_common()),
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(doc, fh, indent=1, sort_keys=False)
    print(f"wrote {args.out}")

    print("\ntop packages by type references:")
    for pkg, d in list(doc["packages"].items())[:20]:
        print(f"  {pkg:<20} refs={d['type_refs']:<6} classes={d['headers']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command")

    prof = sub.add_parser("profile", help="profile how an application exercises the reference")
    prof.add_argument("--reference", required=True,
                      help="root of a reference kernel's source tree (not vendored, not a build dep)")
    prof.add_argument("--consumer", required=True,
                      help="root of an application that consumes it")
    prof.add_argument("--out", default="docs/api_surface.json")
    prof.add_argument("--libclang", default=None)
    prof.set_defaults(func=cmd_profile)

    scope = sub.add_parser("scope", help="derive the parity index from the reference tree")
    scope.add_argument("--reference", required=True,
                       help="root of a reference kernel's checkout (adm/MODULES must exist)")
    scope.add_argument("--out", default="docs/parity/reference-index.tsv")
    scope.add_argument("--profile", default=SCOPE_DEFAULT_PROFILE,
                       help="api_surface.json for the usage column (0s if absent)")
    scope.add_argument("--sample", type=int, default=0, metavar="N",
                       help="print a deterministic stratified sample of N dropped headers")
    scope.set_defaults(func=cmd_scope)

    # Back-compat: the tool was a flat CLI before it had subcommands, and the
    # docstring of every previously generated api_surface.json points at that
    # invocation.
    argv = sys.argv[1:]
    if argv and argv[0].startswith("--") and argv[0] not in ("-h", "--help"):
        argv = ["profile"] + argv
    args = ap.parse_args(argv)
    if not args.command:
        ap.print_help()
        return 2
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
