"""build_clean.py — turn raw SUGRA catalog archives into a clean,
per-T / per-combo tree (one folder per tensor number, one file per combo).

It reads ONLY the "canonical" catalog directories and ignores the redundant
auxiliary ones, so no deduplication is needed:

    KEEP : cat_phase1_nhc_merged_*          (single externals + NHC clusters)
           cat_nhc_ext_nhc_ext_unified_*sext  (multi-external combos; covers
                                                both ...sext and ...nsext)
    SKIP : cat_1ext_nosu2_*   (mirror subset of phase1)
           cat_nhc_ext_t*     (non-unified; mirror subset)
           *_nonsugra         (non-SUGRA variants)

Output layout (same as the old "catalogs_clean"):

    <out>/T<n>_<n>/<combo>.cat

USAGE
-----
    # source = a folder of per-T archives  (e.g. "sugra blocks/T8_8.zip ...")
    python build_clean.py "sugra blocks" catalog_clean

    # source = a single per-T zip, or an already-extracted tree, also work:
    python build_clean.py T8_8.zip catalog_clean
    python build_clean.py extracted_tree/ catalog_clean

A body-level safety dedup is applied and the number removed is reported; for
canonical-only input it is 0 (proof there were no real duplicates).
"""
from __future__ import annotations

import os
import re
import sys
import zipfile

ENTRY_RE = re.compile(r"ENTRY\b.*?\nEND", re.S)
_T_IN_PATH = re.compile(r"(?:^|/)T(\d+)_\d+/")


def is_canonical(catdir: str) -> bool:
    """True for the directories that make up the canonical catalog."""
    catdir = catdir.rstrip("/").rsplit("/", 1)[-1]
    if catdir.endswith("_nonsugra"):
        return False
    if catdir.startswith("cat_phase1_nhc_merged_"):
        return True
    return bool(re.match(r"cat_nhc_ext_nhc_ext_unified_.*sext$", catdir))


def _combo_of(member_path: str) -> str:
    return os.path.basename(member_path)[:-4]  # strip ".cat"


def _entries(text: str):
    return ENTRY_RE.findall(text)


def _add(groups: dict, combo: str, text: str) -> int:
    """Add entries from one .cat text to the per-combo group; return # added."""
    g = groups.setdefault(combo, {"seen": set(), "order": []})
    added = 0
    for blk in _entries(text):
        body = "\n".join(blk.split("\n")[1:])  # drop the "ENTRY <id>" line
        if body not in g["seen"]:
            g["seen"].add(body)
            g["order"].append(body)
            added += 1
    return added


def _write(out_t_dir: str, groups: dict) -> int:
    os.makedirs(out_t_dir, exist_ok=True)
    n = 0
    for combo, g in groups.items():
        if not g["order"]:        # skip combos with zero bases (empty .cat)
            continue
        with open(os.path.join(out_t_dir, combo + ".cat"), "w") as f:
            for i, body in enumerate(g["order"]):
                f.write(f"ENTRY {i}\n{body}\n")
                n += 1
    written = sum(1 for g in groups.values() if g["order"])
    return n, written


# ── per-T processing from a zip or an extracted dir ───────────────────────
def _process_zip(zip_path: str, out_root: str, stats: dict):
    with zipfile.ZipFile(zip_path) as zf:
        # group canonical .cat members by T_H
        byT: dict[int, list[str]] = {}
        for name in zf.namelist():
            if not name.endswith(".cat"):
                continue
            parent = name.rsplit("/", 1)[0] if "/" in name else ""
            if not is_canonical(parent):
                continue
            m = _T_IN_PATH.search(name) or re.search(r"_t(\d+?)\1?_su2", parent)
            t = int(m.group(1)) if m else _infer_t(name)
            if t is None:
                continue
            byT.setdefault(t, []).append(name)
        for t, members in byT.items():
            groups: dict = {}
            raw = 0
            for mem in members:
                txt = zf.read(mem).decode("utf-8", errors="replace")
                blks = _entries(txt)
                raw += len(blks)
                _add(groups, _combo_of(mem), txt)
            kept, ncombo = _write(os.path.join(out_root, f"T{t}_{t}"), groups)
            stats["raw"] += raw
            stats["kept"] += kept
            stats["combos"] += ncombo
            print(f"  T{t}: combos={ncombo} entries={kept} "
                  f"(dups removed={raw - kept})", flush=True)


def _process_dir(tree_root: str, out_root: str, stats: dict):
    byT: dict[int, list[str]] = {}
    for dirpath, _dirs, files in os.walk(tree_root):
        if not is_canonical(dirpath):
            continue
        rel = "/" + os.path.relpath(dirpath, tree_root).replace(os.sep, "/") + "/"
        m = _T_IN_PATH.search(rel) or re.search(r"_t(\d+?)\1?_su2", dirpath)
        t = int(m.group(1)) if m else None
        for fn in files:
            if fn.endswith(".cat") and t is not None:
                byT.setdefault(t, []).append(os.path.join(dirpath, fn))
    for t, paths in byT.items():
        groups: dict = {}
        raw = 0
        for p in paths:
            txt = open(p, encoding="utf-8", errors="replace").read()
            raw += len(_entries(txt))
            _add(groups, _combo_of(p), txt)
        kept, ncombo = _write(os.path.join(out_root, f"T{t}_{t}"), groups)
        stats["raw"] += raw
        stats["kept"] += kept
        stats["combos"] += ncombo
        print(f"  T{t}: combos={ncombo} entries={kept} "
              f"(dups removed={raw - kept})", flush=True)


def _infer_t(name: str):
    m = re.search(r"t(\d+)_su2", name)
    if m:
        s = m.group(1)
        return int(s[: len(s) // 2]) if len(s) % 2 == 0 else None
    return None


def build(source: str, out_root: str):
    stats = {"raw": 0, "kept": 0, "combos": 0}
    if os.path.isdir(source):
        zips = sorted(f for f in os.listdir(source)
                      if re.match(r"T\d+_\d+\.zip$", f))
        if zips:                                  # folder of per-T zips
            for z in sorted(zips, key=lambda s: int(re.match(r"T(\d+)_", s).group(1))):
                print(f"[{z}]")
                _process_zip(os.path.join(source, z), out_root, stats)
        else:                                     # already-extracted tree
            _process_dir(source, out_root, stats)
    elif source.endswith(".zip"):
        _process_zip(source, out_root, stats)
    else:
        raise SystemExit(f"don't know how to read source: {source!r}")
    print(f"\nDONE -> {out_root}/")
    print(f"  combos(files)={stats['combos']:,}  entries={stats['kept']:,}  "
          f"dups removed={stats['raw'] - stats['kept']:,}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        raise SystemExit("usage: python build_clean.py <source> <out_dir>")
    build(sys.argv[1], sys.argv[2])
