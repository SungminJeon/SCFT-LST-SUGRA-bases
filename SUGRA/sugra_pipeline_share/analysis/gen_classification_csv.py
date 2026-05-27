#!/usr/bin/env python3
"""Generate classification_t010.csv from reg chain + nhc chain + 1-ext results.
Follows the same logic as classification_gb_v2.nb:
- LSTGauge: NHC decomposition of non-(-1), non-external curves
- SUGRAGauge: NHC decomposition of ALL non-(-1) curves + hat1/su2 enhancement
- ExtGauge: external tags (su2, su3, e6, etc.)
"""

import os, sys, csv
import numpy as np
from collections import Counter

# ── NHC Table (matches C++ and Mathematica) ──
NHC_TABLE = []

def _build_nhc(name, self_ints, H, V, ints=None):
    n = len(self_ints)
    IF = np.zeros((n, n))
    for i in range(n):
        IF[i, i] = self_ints[i]
        if i + 1 < n:
            k = ints[i] if (ints and i < len(ints)) else 1
            IF[i, i+1] = k
            IF[i+1, i] = k
    eigs = sorted(np.linalg.eigvalsh(IF))
    NHC_TABLE.append({"name": name, "selfInts": self_ints, "H": H, "V": V, "eigs": eigs})

_build_nhc("-3", [-3], 0, 8)
_build_nhc("-4", [-4], 0, 28)
_build_nhc("-5", [-5], 0, 52)
_build_nhc("-6", [-6], 0, 78)
_build_nhc("-7", [-7], 28, 133)
_build_nhc("-8", [-8], 0, 133)
_build_nhc("-12", [-12], 0, 248)
_build_nhc("-2-3", [-2, -3], 8, 17)
_build_nhc("-2-4", [-2, -4], 256, 183, [2])
_build_nhc("-2-3-2", [-2, -3, -2], 16, 27)
_build_nhc("-2-2-3", [-2, -2, -3], 8, 17)

# NHC name -> gauge string (from Mathematica nhcToGauge)
NHC_TO_GAUGE = {
    "-3": "su3",
    "-4": "so8",
    "-5": "f4",
    "-6": "e6",
    "-7": "e7'",
    "-8": "e7",
    "-12": "e8",
    "-2-3": "su2+g2",
    "-2-4": "su8+so16",
    "-2-3-2": "su2+so7+su2",
    "-2-2-3": "sp1+g2",
    "pure-2": "",
    "???": "???",
}

# External gauge mapping
EXT_SI_TO_GAUGE = {
    -2: "su2", -3: "su3", -4: "so8", -5: "f4",
    -6: "e6", -7: "e7'", -8: "e7", -12: "e8",
}


def nhc_lookup(sub_if):
    """Match submatrix to NHC table by eigenvalue comparison. Returns name string."""
    eigs = sorted(np.linalg.eigvalsh(sub_if))
    for nhc in NHC_TABLE:
        if len(eigs) == len(nhc["eigs"]) and np.max(np.abs(np.array(eigs) - np.array(nhc["eigs"]))) < 1e-6:
            return nhc["name"]
    # Pure -2 check
    if all(sub_if[i, i] == -2 for i in range(len(sub_if))):
        return "pure-2"
    return "???"


def find_connected_components(matrix, indices):
    """BFS connected components among given indices."""
    if not indices:
        return []
    idx_set = set(indices)
    visited = set()
    components = []
    for start in sorted(indices):
        if start in visited:
            continue
        comp = []
        queue = [start]
        visited.add(start)
        while queue:
            cur = queue.pop(0)
            comp.append(cur)
            for nb in sorted(idx_set):
                if nb not in visited and matrix[cur][nb] != 0:
                    visited.add(nb)
                    queue.append(nb)
        components.append(sorted(comp))
    return components


def nhc_gauge_for_component(matrix, comp):
    """Get gauge string for a connected component via NHC lookup."""
    n = len(comp)
    sub = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            sub[a, b] = matrix[comp[a]][comp[b]]
    name = nhc_lookup(sub)
    return NHC_TO_GAUGE.get(name, name)


def format_gauge_list(gauges):
    """Format gauge list: sorted, with counts. e.g. '2 su3 + e6'"""
    if not gauges:
        return "(none)"
    counts = Counter(gauges)
    parts = []
    for g in sorted(counts.keys()):
        c = counts[g]
        if "+" in g:
            wrap = f"({g})"
        else:
            wrap = g
        if c > 1:
            parts.append(f"{c} {wrap}")
        else:
            parts.append(g)
    return " + ".join(parts) if parts else "(none)"


# ── Entry parsing ──

def parse_cat_file(filepath):
    """Parse a .cat file, yield entries as dicts."""
    with open(filepath) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("ENTRY"):
            entry = {"externals": [], "if_matrix": []}
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if line == "END":
                    yield entry
                    i += 1
                    break
                if line.startswith("SPEC"):
                    parts = line.split()
                    entry["externals"].append({
                        "curve_idx": None,  # filled later from ATTACH
                        "tag": parts[2],
                        "ext_si": int(parts[3]),
                        "target_si": int(parts[4]),
                        "int_num": int(parts[5]),
                        "is_hat1": int(parts[6]) == 1,
                    })
                elif line.startswith("COMBO"):
                    entry["combo"] = line.split(None, 1)[1]
                elif line.startswith("EXTERNALS"):
                    n_ext = int(line.split()[1])
                    for j in range(n_ext):
                        i += 1
                        parts = lines[i].strip().split()
                        entry["externals"].append({
                            "curve_idx": int(parts[0]),
                            "spec_id": int(parts[1]),
                            "tag": parts[2],
                            "ext_si": int(parts[3]),
                            "target_si": int(parts[4]),
                            "int_num": int(parts[5]),
                            "is_hat1": int(parts[6]) == 1,
                        })
                elif line.startswith("BASE"):
                    parts = line.split()
                    entry["base_si"] = int(parts[1])
                    entry["catalog_type"] = parts[2]
                    entry["base_T"] = int(parts[3])
                    entry["if_size"] = int(parts[4])
                elif line.startswith("ATTACH"):
                    parts = line.split()
                    entry["attach_idx"] = int(parts[1])
                    entry["attach_if_size"] = int(parts[2])
                elif line.startswith("PHYSICS"):
                    parts = line.split()
                    entry["T"] = int(parts[1])
                    entry["H_charged"] = int(parts[2])
                    entry["V"] = int(parts[3])
                    entry["H_neutral"] = int(parts[4])
                elif line.startswith("IF"):
                    n = int(line.split()[1])
                    matrix = []
                    for j in range(n):
                        i += 1
                        row = list(map(int, lines[i].strip().split()))
                        matrix.append(row)
                    entry["if_matrix"] = matrix
                i += 1
        else:
            i += 1


def get_ext_indices(entry):
    """Get set of curve indices that are externals."""
    indices = set()
    for ext in entry["externals"]:
        if ext.get("curve_idx") is not None:
            indices.add(ext["curve_idx"])
    # Phase1: SPEC format, external is last curve in IF
    if not indices and "attach_if_size" in entry:
        indices.add(len(entry["if_matrix"]) - 1)
    return indices


# ── Classification functions (matching Mathematica) ──

# NHC tag -> cluster gauge name
NHC_TAG_TO_GAUGE = {
    "nhc_2_3": "su2+g2",
    "nhc_2_3_2": "su2+so7+su2",
    "nhc_2_2_3": "sp1+g2",
    "nhc_2_4": "su8+so16",
}

def get_ext_gauge(entry):
    """ExtGauge: gauge names from external specs, NHC clusters as units."""
    gauges = []
    # Count NHC cluster curves to determine how many clusters
    nhc_curve_counts = {}  # tag -> count of curves seen
    for ext in entry["externals"]:
        tag = ext.get("tag", "")
        if tag in NHC_TAG_TO_GAUGE:
            nhc_curve_counts[tag] = nhc_curve_counts.get(tag, 0) + 1
    # Convert curve counts to cluster counts
    nhc_cluster_added = {}  # tag -> clusters already added
    for tag, count in nhc_curve_counts.items():
        curves_per_cluster = {"nhc_2_3": 2, "nhc_2_3_2": 3, "nhc_2_2_3": 3, "nhc_2_4": 2}.get(tag, 2)
        n_clusters = count // curves_per_cluster
        for _ in range(n_clusters):
            gauges.append(NHC_TAG_TO_GAUGE[tag])
    nhc_tags_done = set(nhc_curve_counts.keys())
    for ext in entry["externals"]:
        tag = ext.get("tag", "")
        if tag in nhc_tags_done:
            continue
        # Mixed externals: use actual gauge name
        TAG_TO_GAUGE = {"su2n3": "su2", "su3n2": "g2", "su8": "su8", "so16n2": "so16",
                        "su2n3mix": "su2", "su3n2mix": "g2", "su8mix": "su8", "so16n2mix": "so16",
                        "so7": "so7", "so7mix": "so7"}
        if tag in TAG_TO_GAUGE:
            gauges.append(TAG_TO_GAUGE[tag])
            continue
        si = ext["ext_si"]
        if ext["is_hat1"]:
            gauges.append("hat1")
        elif si in EXT_SI_TO_GAUGE:
            gauges.append(EXT_SI_TO_GAUGE[si])
        else:
            gauges.append(f"si{si}")
    return format_gauge_list(gauges)


def get_lst_gauge(entry):
    """LSTGauge: NHC decomposition of non-(-1), non-external curves."""
    matrix = entry["if_matrix"]
    n = len(matrix)
    ext_indices = get_ext_indices(entry)

    non_m1 = [i for i in range(n) if matrix[i][i] != -1 and i not in ext_indices]
    if not non_m1:
        return "(none)"

    components = find_connected_components(matrix, non_m1)
    gauges = []
    for comp in components:
        g = nhc_gauge_for_component(matrix, comp)
        if g and g != "":
            gauges.append(g)
    return format_gauge_list(gauges)


def get_sugra_gauge(entry):
    """SUGRAGauge: NHC decomposition of ALL non-(-1) curves + hat1 + standalone su2 enhancement."""
    matrix = entry["if_matrix"]
    n = len(matrix)
    ext_indices = get_ext_indices(entry)

    # 1. NHC decomposition on all non-(-1) curves
    non_m1 = [i for i in range(n) if matrix[i][i] != -1]
    components = find_connected_components(matrix, non_m1)
    gauges = []
    for comp in components:
        g = nhc_gauge_for_component(matrix, comp)
        if g and g != "":
            gauges.append(g)

    # 2. Hat1 enhancement
    for ext in entry["externals"]:
        if ext["is_hat1"]:
            if ext["target_si"] == -1:
                gauges.append("su8")
            else:
                gauges.append("su16")

    # 3. Standalone su2 on (-1): external -2 targeting -1, only connects to -1 curves
    for ext in entry["externals"]:
        if not ext["is_hat1"] and ext["ext_si"] == -2 and ext["target_si"] == -1:
            ci = ext.get("curve_idx")
            if ci is None:
                # Phase1: last curve
                ci = len(matrix) - 1
            standalone = True
            for j in range(n):
                if j != ci and matrix[ci][j] != 0 and matrix[j][j] != -1:
                    standalone = False
                    break
            if standalone:
                gauges.append("su2")

    return format_gauge_list(gauges)


# ── Directory scanning ──

def scan_directories():
    dirs = []
    base = "/Users/seongminjeon/Desktop/sugra_pipeline"

    # 1-ext regular (no-su2)
    d = os.path.join(base, "cat_1ext_nosu2_t010")
    if os.path.isdir(d):
        dirs.append(d)

    # 1-ext NHC
    d = os.path.join(base, "cat_nhc_ext")
    if os.path.isdir(d):
        dirs.append(d)

    # unified chain (SUGRA only) — handle both naming conventions
    for r in range(2, 11):
        for prefix in ["cat_nhc_ext_unified", "cat_nhc_ext_nhc_ext_unified"]:
            for suffix in [f"{prefix}_r{r}sext", f"{prefix}_r{r}nsext"]:
                d = os.path.join(base, suffix)
                if os.path.isdir(d):
                    dirs.append(d)

    return dirs


def main():
    dirs = scan_directories()
    print(f"Scanning {len(dirs)} directories...")

    classification = Counter()
    total_entries = 0

    for d in dirs:
        cat_files = [f for f in os.listdir(d) if f.endswith(".cat")]
        for cat_file in cat_files:
            filepath = os.path.join(d, cat_file)
            try:
                for entry in parse_cat_file(filepath):
                    # Skip entries with su2 externals (tag="su2", not NHC)
                    has_su2 = any(ext["tag"] == "su2"
                                 for ext in entry["externals"])
                    if has_su2:
                        continue

                    # Skip entries with total ext cc > 16
                    total_cc = 0.0
                    for ext in entry["externals"]:
                        si, k = ext["ext_si"], ext["int_num"]
                        if ext["is_hat1"]:
                            if ext["target_si"] == -1:
                                total_cc += 63/8  # su8: 63*1/(1+8)=7.0
                            else:
                                total_cc += 255/17  # su16: 255*1/(1+16)=15.0
                        else:
                            gauge = {-2:(3,2), -3:(8,3), -4:(28,6), -5:(52,9),
                                     -6:(78,12), -7:(133,18), -8:(133,18), -12:(248,30)}
                            if si in gauge:
                                dim, h = gauge[si]
                                total_cc += k*dim/(k+h)
                    if total_cc > 16.0 + 1e-9:
                        continue

                    T_H = entry.get("base_T", entry.get("T", 0))
                    ext_gauge = get_ext_gauge(entry)
                    lst_gauge = get_lst_gauge(entry)
                    sugra_gauge = get_sugra_gauge(entry)

                    classification[(T_H, ext_gauge, lst_gauge, sugra_gauge)] += 1
                    total_entries += 1

                    if total_entries % 500000 == 0:
                        print(f"  Processed {total_entries} entries...")
            except Exception as e:
                print(f"  Error in {filepath}: {e}", file=sys.stderr)

    print(f"Total entries: {total_entries}")
    print(f"Unique classifications: {len(classification)}")

    rows = sorted(classification.items(), key=lambda x: (x[0][0], -x[1], x[0][1]))

    outfile = os.path.join("/Users/seongminjeon/Desktop/sugra_pipeline", "classification_t010.csv")
    with open(outfile, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["T_H", "External", "LST_Gauge", "SUGRA_Gauge", "Count"])
        for (t, ext, lst, sugra), count in rows:
            w.writerow([t, ext, lst, sugra, count])

    print(f"Written to {outfile}")


if __name__ == "__main__":
    main()
