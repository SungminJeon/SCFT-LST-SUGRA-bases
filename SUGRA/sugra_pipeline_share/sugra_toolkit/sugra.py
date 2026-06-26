"""sugra.py — friendly loader for the 6D SUGRA base catalog.

The catalog is laid out one folder per LST tensor number, one file per
external combo:

    <root>/T<n>_<n>/<combo>.cat        e.g.  catalog/T8_8/e6+su3.cat

Each .cat file holds many "bases" (one ENTRY..END block each). You almost
never need to know the on-disk format — just use the API below.

QUICK START
-----------
    import sugra
    cat = sugra.Catalog("catalog")        # point at the data root folder

    cat.tensors()                         # -> [2, 3, 4, ...] T_H values present
    cat.combos(8)[:10]                    # -> combo names available at T_H=8

    bases = cat.load(8, "e6+su3")         # -> list[Base] for that T_H + combo
    b = bases[0]
    b.T_H, b.T, b.combo                   # tensor count, extended T, combo
    b.tags                                # ['e6', 'su3']
    b.physics                             # dict: Hc, V, Hn, det, sig_*
    b.IF                                  # numpy intersection-form matrix
    b.t_min()                             # T_min (None if undefined)

    # Don't know the combo? Filter across everything:
    hits = cat.find(t_h=8, externals=["e6", "su3"])
    hits = cat.find(externals=["e8"])           # any base containing an e8
    hits = cat.find(t_h=5, t_min=8)

No special setup, no multiprocessing — safe to `import sugra` and use in any
script, notebook, or REPL.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

__all__ = ["Catalog", "Base"]


# ── one SUGRA base ────────────────────────────────────────────────────────
@dataclass
class Base:
    """A single 6D SUGRA base (one ENTRY in a .cat file)."""

    id: int
    combo: str                       # "+"-joined external tags, e.g. "e6+su3"
    externals: list[dict]            # each: curveIdx, specId, tag, extSI,
                                     #       targetSI, intNum, isHat1
    catalog_id: int                  # underlying LST id
    catalog_type: str                # LST family label, e.g. "NN", "DM:A(3)"
    T_H: int                         # tensor count of the underlying LST
    T: int                           # extended-base T (= sig_neg)
    Hc: int                          # charged hypers
    V: int                           # vectors
    Hn: int                          # neutral hypers
    det: int
    sig_pos: int
    sig_neg: int
    sig_zero: int
    IF: np.ndarray                   # intersection-form matrix (n x n)
    remaining: list[dict] = field(default_factory=list)
    source: str = ""                 # file the base came from

    # convenience views ----------------------------------------------------
    @property
    def tags(self) -> list[str]:
        """External tags as they appear, e.g. ['e6', 'su3']."""
        return [e["tag"] for e in self.externals]

    @property
    def tag_set(self) -> list[str]:
        """Sorted external tags (order-independent comparison)."""
        return sorted(self.tags)

    @property
    def physics(self) -> dict:
        return {"T": self.T, "Hc": self.Hc, "V": self.V, "Hn": self.Hn,
                "det": self.det, "sig": (self.sig_pos, self.sig_neg, self.sig_zero)}

    def t_min(self) -> int | None:
        """Minimum T admitting a SUGRA b0 lift; None if undefined."""
        return _t_min(self.IF, [e["curveIdx"] for e in self.externals
                                if e.get("isHat1") and e.get("curveIdx", -1) >= 0])

    def __repr__(self) -> str:
        return (f"Base(id={self.id}, T_H={self.T_H}, T={self.T}, "
                f"combo='{self.combo}', IF={self.IF.shape[0]}x{self.IF.shape[0]})")


# ── the catalog ───────────────────────────────────────────────────────────
class Catalog:
    """Read-only view over a per-T / per-combo .cat tree."""

    def __init__(self, root: str):
        root = os.path.expanduser(root)
        if not os.path.isdir(root):
            raise FileNotFoundError(f"catalog root not found: {root!r}")
        self.root = root
        # map T_H -> folder path (handles T8_8, T8, etc.)
        self._tdir: dict[int, str] = {}
        for name in os.listdir(root):
            p = os.path.join(root, name)
            if not os.path.isdir(p):
                continue
            m = re.match(r"T(\d+)(?:_\d+)?$", name)
            if m:
                self._tdir[int(m.group(1))] = p
        if not self._tdir:
            raise FileNotFoundError(
                f"no T<n> folders under {root!r}. Expected e.g. {root}/T8_8/")

    # discovery ------------------------------------------------------------
    def tensors(self) -> list[int]:
        """Sorted list of T_H values present."""
        return sorted(self._tdir)

    def combos(self, t_h: int) -> list[str]:
        """Sorted combo names available at this T_H."""
        d = self._tdir.get(t_h)
        if d is None:
            raise KeyError(f"T_H={t_h} not found. Available: {self.tensors()}")
        return sorted(f[:-4] for f in os.listdir(d) if f.endswith(".cat"))

    # loading --------------------------------------------------------------
    def load(self, t_h: int, combo: str | None = None) -> list[Base]:
        """Load bases for a T_H, optionally a single combo.

        cat.load(8)            -> every base at T_H=8
        cat.load(8, "e6+su3")  -> bases for that exact combo only
        """
        d = self._tdir.get(t_h)
        if d is None:
            raise KeyError(f"T_H={t_h} not found. Available: {self.tensors()}")
        if combo is not None:
            path = os.path.join(d, combo + ".cat")
            if not os.path.isfile(path):
                raise FileNotFoundError(
                    f"combo {combo!r} not found at T_H={t_h}. "
                    f"Try cat.combos({t_h}) to list them.")
            return _parse_file(path)
        out: list[Base] = []
        for f in sorted(os.listdir(d)):
            if f.endswith(".cat"):
                out.extend(_parse_file(os.path.join(d, f)))
        return out

    # filtering ------------------------------------------------------------
    def find(self, t_h: int | None = None, t_min: int | None = None,
             externals=None, combo: str | None = None,
             contains: str | None = None) -> list[Base]:
        """Filter bases by any combination of criteria.

        t_h        : exact tensor count (or None = all T_H — can be slow)
        combo      : exact combo string  (fast: reads one file per T_H)
        externals  : multiset match, order-independent, e.g. ["e6", "su3"]
        contains   : at least one external with this tag, e.g. "e8"
        t_min      : exact T_min value
        """
        t_list = [t_h] if t_h is not None else self.tensors()
        ext_target = sorted(externals) if externals is not None else None
        out: list[Base] = []
        for t in t_list:
            bases = self.load(t, combo) if combo else self.load(t)
            for b in bases:
                if ext_target is not None and b.tag_set != ext_target:
                    continue
                if contains is not None and contains not in b.tags:
                    continue
                if t_min is not None and b.t_min() != t_min:
                    continue
                out.append(b)
        return out


# ── .cat parsing (you don't need to call these directly) ──────────────────
def _parse_file(path: str) -> list[Base]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    return _parse_text(text, os.path.basename(path))


def _parse_text(text: str, source: str = "") -> list[Base]:
    lines = text.split("\n")
    n = len(lines)
    out: list[Base] = []
    i = 0
    while i < n:
        if not lines[i].startswith("ENTRY"):
            i += 1
            continue
        try:
            ent_id = int(lines[i].split()[1])
        except (IndexError, ValueError):
            ent_id = 0
        combo = ""
        externals: list[dict] = []
        cat_id = base_t = T = Hc = V = Hn = det = sp = sn = sz = 0
        cat_type = ""
        IF = np.zeros((0, 0), dtype=np.int64)
        remaining: list[dict] = []
        i += 1
        while i < n:
            ln = lines[i]
            if ln.startswith("END"):
                i += 1
                break
            if ln.startswith("COMBO"):
                combo = ln[6:].strip()
            elif ln.startswith("EXTERNALS"):
                k = int(ln.split()[1])
                for _ in range(k):
                    i += 1
                    p = lines[i].split()
                    if len(p) >= 7:
                        externals.append({
                            "curveIdx": int(p[0]), "specId": int(p[1]),
                            "tag": p[2], "extSI": int(p[3]),
                            "targetSI": int(p[4]), "intNum": int(p[5]),
                            "isHat1": p[6] == "1"})
            elif ln.startswith("SPEC"):  # phase1-style single external
                p = ln.split()
                if len(p) >= 7:
                    externals.append({
                        "curveIdx": -1, "specId": int(p[1]), "tag": p[2],
                        "extSI": int(p[3]), "targetSI": int(p[4]),
                        "intNum": int(p[5]), "isHat1": p[6] == "1"})
                    if not combo:
                        combo = p[2]
            elif ln.startswith("ATTACH"):
                p = ln.split()
                if len(p) >= 3 and externals and externals[-1]["curveIdx"] == -1:
                    externals[-1]["curveIdx"] = int(p[2])
            elif ln.startswith("BASE"):
                p = ln.split()
                if len(p) >= 5:
                    cat_id, cat_type, base_t = int(p[1]), p[2], int(p[3])
            elif ln.startswith("PHYSICS"):
                p = ln.split()
                if len(p) >= 9:
                    T, Hc, V, Hn, det, sp, sn, sz = (int(x) for x in p[1:9])
            elif ln.startswith("IF"):
                k = int(ln[3:].strip())
                mat = np.zeros((k, k), dtype=np.int64)
                for r in range(k):
                    i += 1
                    row = lines[i].split()
                    if len(row) == k:
                        mat[r] = [int(x) for x in row]
                IF = mat
            elif ln.startswith("REMAINING"):
                k = int(ln.split()[1])
                for _ in range(k):
                    i += 1
                    p = lines[i].split()
                    if len(p) >= 3:
                        remaining.append({"curveIdx": int(p[0]),
                                          "selfInt": int(p[1]),
                                          "availCC": float(p[2])})
            i += 1
        if not combo and externals:
            combo = "+".join(sorted(e["tag"] for e in externals))
        out.append(Base(ent_id, combo, externals, cat_id, cat_type, base_t,
                        T, Hc, V, Hn, det, sp, sn, sz, IF, remaining, source))
    return out


def _t_min(IF: np.ndarray, hat1=()) -> int | None:
    """T_min (T_crit). None if the IF is unimodular (degenerate b0 column)."""
    n = len(IF)
    if n == 0:
        return None
    h1 = set(hat1)
    M = np.zeros((n + 1, n + 1), dtype=np.float64)
    M[:n, :n] = IF
    for i in range(n):
        b = -1.0 if i in h1 else float(IF[i, i] + 2)
        M[i, n] = M[n, i] = b
    M[n, n] = 0.0
    B = float(np.linalg.det(M))
    M[n, n] = 1.0
    A = float(np.linalg.det(M)) - B
    if abs(A) < 1e-10:
        return None
    return int(np.ceil(9.0 - (-B / A) - 1e-9))
