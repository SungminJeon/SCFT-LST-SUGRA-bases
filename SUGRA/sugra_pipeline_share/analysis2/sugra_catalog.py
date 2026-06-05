"""Python port of the SUGRACatalog query interface.

Mirrors the Mathematica example (example_for_hamada.m): load entries from a
catalog tree (directory **or zip archive**), then filter by T_H, T_min, and
external types.

Key design choices:
- Reads .cat files directly from zip archives, no extraction step.
- Uses numpy for the T_min linear-algebra (much faster than Mathematica).
- Parallel parsing via ProcessPoolExecutor for large catalogs.

CAT file grammar (as emitted by the pipeline):
    ENTRY <id>
    BASE <catalogId> <catalogType> <baseT> <ifSize>     (optional)
    SPEC <specId> <tag> <extSI> <targetSI> <intNum> <isHat1>   (one per ext, before EXTERNALS)
    COMBO <combo-string>                                 (optional)
    EXTERNALS <n>
        <curveIdx> <specId> <tag> <extSI> <targetSI> <intNum> <isHat1>
        ...
    ATTACH <attachIdx> <ifSize>                          (optional)
    PHYSICS <T> <Hc> <V> <Hn> <det> <sigPos> <sigNeg> <sigZero>
    IF <n>
        <n integers per row>  x n
    REMAINING <n>                                        (optional)
        <curveIdx> <selfInt> <availCC>  x n
    GLUE <a> <Ta> <b> <Tb>                               (optional)
    BLOWDOWN <n>  ... BDSIG ...                          (optional)
    END
"""
from __future__ import annotations

import os
import re
import zipfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Iterable, Iterator

import numpy as np


# ── External SI → gauge name ────────────────────────────────────────────
EXT_SI_TO_GAUGE = {
    -2: "su2", -3: "su3", -4: "so8", -5: "f4",
    -6: "e6", -7: "e7'", -8: "e7", -12: "e8",
}


# ── Entry data class ────────────────────────────────────────────────────
@dataclass
class Entry:
    id: int = 0
    combo: str = ""
    externals: list[dict] = field(default_factory=list)
    catalog_id: int = 0
    catalog_type: str = ""
    base_t: int = 0
    T: int = 0
    Hc: int = 0
    V: int = 0
    Hn: int = 0
    det: int = 0
    sig_pos: int = 0
    sig_neg: int = 0
    sig_zero: int = 0
    IF: np.ndarray = field(default_factory=lambda: np.zeros((0, 0), dtype=np.int64))
    remaining: list[dict] = field(default_factory=list)
    source: str = ""  # which file/zip the entry came from

    @property
    def ext_tags(self) -> list[str]:
        return [e["tag"] for e in self.externals]

    @property
    def ext_tag_set(self) -> list[str]:
        return sorted(self.ext_tags)

    def ext_indices(self) -> list[int]:
        """0-indexed curve indices that are externals."""
        return [e["curveIdx"] for e in self.externals if e.get("curveIdx", -1) >= 0]


# ── Parser ──────────────────────────────────────────────────────────────
def _parse_cat_text(text: str, source: str = "") -> list[Entry]:
    entries: list[Entry] = []
    lines = text.split("\n")
    n_lines = len(lines)
    i = 0
    while i < n_lines:
        line = lines[i].strip()
        if not line.startswith("ENTRY"):
            i += 1
            continue
        try:
            ent_id = int(line[6:].strip())
        except ValueError:
            ent_id = 0
        e = Entry(id=ent_id, source=source)
        i += 1
        while i < n_lines:
            line = lines[i].strip()
            if line == "END":
                entries.append(e)
                i += 1
                break
            if line.startswith("SPEC"):
                p = line.split()
                if len(p) >= 7:
                    e.externals.append({
                        "curveIdx": -1, "specId": int(p[1]), "tag": p[2],
                        "extSI": int(p[3]), "targetSI": int(p[4]),
                        "intNum": int(p[5]), "isHat1": p[6] == "1",
                    })
                    if not e.combo:
                        e.combo = p[2]
            elif line.startswith("COMBO"):
                e.combo = line[6:].strip()
            elif line.startswith("EXTERNALS"):
                n_ext = int(line.split()[1])
                e.externals = []  # full list overrides any SPEC stubs
                for _ in range(n_ext):
                    i += 1
                    if i >= n_lines:
                        break
                    p = lines[i].split()
                    if len(p) >= 7:
                        e.externals.append({
                            "curveIdx": int(p[0]), "specId": int(p[1]), "tag": p[2],
                            "extSI": int(p[3]), "targetSI": int(p[4]),
                            "intNum": int(p[5]), "isHat1": p[6] == "1",
                        })
            elif line.startswith("BASE"):
                p = line.split()
                if len(p) >= 5:
                    e.catalog_id = int(p[1])
                    e.catalog_type = p[2]
                    e.base_t = int(p[3])
            elif line.startswith("ATTACH"):
                # ATTACH <attachIdx> <ifSize> -- updates the last SPEC's curveIdx
                p = line.split()
                if len(p) >= 3 and e.externals and e.externals[-1]["curveIdx"] == -1:
                    e.externals[-1]["curveIdx"] = int(p[2])
            elif line.startswith("PHYSICS"):
                p = line.split()
                if len(p) >= 9:
                    e.T = int(p[1])
                    e.Hc = int(p[2])
                    e.V = int(p[3])
                    e.Hn = int(p[4])
                    e.det = int(p[5])
                    e.sig_pos = int(p[6])
                    e.sig_neg = int(p[7])
                    e.sig_zero = int(p[8])
            elif line.startswith("IF"):
                n = int(line[3:].strip())
                mat = np.zeros((n, n), dtype=np.int64)
                for r in range(n):
                    i += 1
                    if i >= n_lines:
                        break
                    row = lines[i].split()
                    if len(row) == n:
                        mat[r] = [int(x) for x in row]
                e.IF = mat
            elif line.startswith("REMAINING"):
                n_rem = int(line.split()[1])
                for _ in range(n_rem):
                    i += 1
                    if i >= n_lines:
                        break
                    p = lines[i].split()
                    if len(p) >= 3:
                        e.remaining.append({
                            "curveIdx": int(p[0]),
                            "selfInt": int(p[1]),
                            "availCC": float(p[2]),
                        })
            # BLOWDOWN/BDSIG/GLUE blocks: skip; not used by the query helpers
            i += 1
    return entries


def parse_cat_file(path: str) -> list[Entry]:
    with open(path, "rb") as f:
        text = f.read().decode("utf-8", errors="replace")
    return _parse_cat_text(text, source=os.path.basename(path))


# ── T_min ───────────────────────────────────────────────────────────────
def t_min(IF: np.ndarray, h1: Iterable[int] = ()) -> int | None:
    """Return T_min (T_crit). None if IF is unimodular (b0 column is degenerate)."""
    n = len(IF)
    if n == 0:
        return None
    h1_set = set(h1)
    M = np.zeros((n + 1, n + 1), dtype=np.float64)
    M[:n, :n] = IF
    for i in range(n):
        b = -1.0 if i in h1_set else float(IF[i, i] + 2)
        M[i, n] = b
        M[n, i] = b
    M[n, n] = 0.0
    B = float(np.linalg.det(M))
    M[n, n] = 1.0
    A = float(np.linalg.det(M)) - B
    if abs(A) < 1e-10:
        return None
    tc = 9.0 - (-B / A)
    return int(np.ceil(tc - 1e-9))


def entry_t_min(e: Entry) -> int | None:
    h1 = [x["curveIdx"] for x in e.externals if x.get("isHat1")]
    return t_min(e.IF, h1)


# ── Discovery helpers ───────────────────────────────────────────────────
_CAT_DIR_KEEP = re.compile(
    r"(cat_phase1_nhc_merged_[^/]+|cat_nhc_ext_nhc_ext_unified_[^/]+sext)$"
)


def _is_kept_cat_dir(name: str) -> bool:
    base = name.rstrip("/").rsplit("/", 1)[-1]
    return bool(_CAT_DIR_KEEP.match(base))


def _t_h_from_path(path: str) -> int | None:
    """Extract the T<n>_<n> tensor number from a path; None if not found."""
    m = re.search(r"(?:^|/)T(\d+)_\1(?=/)", path)
    return int(m.group(1)) if m else None


# ── Loading ─────────────────────────────────────────────────────────────
def _load_zip_batch(args):
    """Open the zip once per worker and parse a batch of members."""
    zip_path, members = args
    base = os.path.basename(zip_path)
    out: list[Entry] = []
    with zipfile.ZipFile(zip_path) as zf:
        for member in members:
            with zf.open(member) as f:
                text = f.read().decode("utf-8", errors="replace")
            out.extend(_parse_cat_text(text, source=f"{base}!{member}"))
    return out


def _load_file_batch(paths: list[str]) -> list[Entry]:
    out: list[Entry] = []
    for p in paths:
        out.extend(parse_cat_file(p))
    return out


def load_catalog(
    path: str,
    t_h_range: tuple[int, int] | None = None,
    workers: int | None = None,
    verbose: bool = True,
) -> list[Entry]:
    """Load SUGRA entries from a directory or a zip archive.

    t_h_range = (lo, hi) inclusive restricts which T<n>_<n> subdirs are read,
    so you can start with a small slice (fast) and enlarge once happy.

    workers = number of parallel processes (default: os.cpu_count()).
              Pass workers=1 to disable parallelism.
    """
    if path.endswith(".zip") or zipfile.is_zipfile(path):
        return _load_from_zip(path, t_h_range, workers, verbose)
    return _load_from_dir(path, t_h_range, workers, verbose)


def _load_from_zip(zip_path: str, t_h_range, workers, verbose) -> list[Entry]:
    with zipfile.ZipFile(zip_path) as zf:
        members = []
        for info in zf.infolist():
            name = info.filename
            if not name.endswith(".cat") or info.is_dir():
                continue
            parent = name.rsplit("/", 1)[0] if "/" in name else ""
            if not _is_kept_cat_dir(parent):
                continue
            if t_h_range is not None:
                t = _t_h_from_path(name)
                if t is None or not (t_h_range[0] <= t <= t_h_range[1]):
                    continue
            members.append(name)
    if verbose:
        print(f"[zip] {os.path.basename(zip_path)}: {len(members)} .cat files")
    if not members:
        return []
    return _run_batched(members, zip_path, workers, verbose, kind="zip")


def _load_from_dir(root: str, t_h_range, workers, verbose) -> list[Entry]:
    cat_files: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        # keep walking into per-tensor T<n>_<n> dirs honoring the range
        if t_h_range is not None:
            rel = os.path.relpath(dirpath, root)
            t = _t_h_from_path("/" + rel.replace(os.sep, "/") + "/")
            if t is not None and not (t_h_range[0] <= t <= t_h_range[1]):
                dirnames[:] = []
                continue
        if not _is_kept_cat_dir(dirpath):
            continue
        for fn in filenames:
            if fn.endswith(".cat"):
                cat_files.append(os.path.join(dirpath, fn))
    if verbose:
        print(f"[dir] {root}: {len(cat_files)} .cat files")
    if not cat_files:
        return []
    return _run_batched(cat_files, None, workers, verbose, kind="dir")


def _run_batched(items, zip_path, workers, verbose, kind) -> list[Entry]:
    if workers is None:
        workers = min(os.cpu_count() or 4, 8)
    if workers == 1 or len(items) < 16:
        if kind == "zip":
            return _load_zip_batch((zip_path, items))
        return _load_file_batch(items)

    # split items into ~workers*2 batches so faster workers can grab more.
    n_batches = max(1, workers * 2)
    chunk = max(1, (len(items) + n_batches - 1) // n_batches)
    batches = [items[i:i + chunk] for i in range(0, len(items), chunk)]
    if kind == "zip":
        args = [(zip_path, b) for b in batches]
        fn = _load_zip_batch
    else:
        args = batches
        fn = _load_file_batch
    out: list[Entry] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for chunk_out in ex.map(fn, args):
            out.extend(chunk_out)
    if verbose:
        print(f"loaded {len(out)} entries (workers={workers}, batches={len(batches)})")
    return out


# ── Queries ─────────────────────────────────────────────────────────────
def pick_bases(
    entries: list[Entry],
    t_h: int | None = None,
    t_min_val: int | None = None,
    externals: Iterable[str] | None = None,
) -> list[Entry]:
    """Filter entries by any combination of T_H, T_min, and (sorted) externals."""
    ext_target = sorted(externals) if externals is not None else None
    out: list[Entry] = []
    for e in entries:
        if t_h is not None and e.base_t != t_h:
            continue
        if t_min_val is not None:
            tm = entry_t_min(e)
            if tm != t_min_val:
                continue
        if ext_target is not None and e.ext_tag_set != ext_target:
            continue
        out.append(e)
    return out


def find_by_combo(entries: list[Entry], combo: str) -> list[Entry]:
    return [e for e in entries if e.combo == combo]


def find_containing_tag(entries: list[Entry], tag: str) -> list[Entry]:
    return [e for e in entries if tag in e.ext_tags]


# ── Pretty printing ─────────────────────────────────────────────────────
def show_entry(e: Entry) -> str:
    lines = [
        f"ENTRY {e.id}  [{e.catalog_type}]  T_H={e.base_t}  T={e.T}",
        f"  combo: {e.combo}",
        f"  externals: {', '.join(e.ext_tags) or '(none)'}",
        f"  H={e.Hc}+{e.Hn}n  V={e.V}  det={e.det}  sig={e.sig_pos},{e.sig_neg},{e.sig_zero}",
    ]
    tm = entry_t_min(e)
    lines.append(f"  T_min: {tm}")
    return "\n".join(lines)


def show_if(e: Entry) -> str:
    if e.IF.size == 0:
        return "(empty IF)"
    return "\n".join(" ".join(f"{v:3d}" for v in row) for row in e.IF)


def table_form(
    entries: list[Entry],
    columns: tuple[str, ...] = ("id", "catalog_type", "base_t", "t_min", "combo", "T"),
) -> str:
    """Aligned text table of the given attribute columns."""
    def cell(e: Entry, col: str):
        if col == "t_min":
            return entry_t_min(e)
        return getattr(e, col, "")
    rows = [[str(cell(e, c)) for c in columns] for e in entries]
    widths = [max(len(c), *(len(r[i]) for r in rows)) if rows else len(c)
              for i, c in enumerate(columns)]
    fmt = " | ".join(f"{{:<{w}}}" for w in widths)
    out = [fmt.format(*columns), "-+-".join("-" * w for w in widths)]
    out.extend(fmt.format(*r) for r in rows)
    return "\n".join(out)
