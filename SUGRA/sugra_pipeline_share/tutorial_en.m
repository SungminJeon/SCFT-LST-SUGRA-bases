(* ::Package:: *)

(*  SUGRA Catalog Tutorial (English)
    ──────────────────────────────────────────────────────────────
    Step 1: run the C++ pipeline first, e.g.
              ./run_pipeline.sh 0 10
            which produces:
              cat_1ext_nosu2_t010/                (Phase 1 output)
              cat_nhc_ext_t010/                   (NHC ext output)
              cat_nhc_ext_nhc_ext_unified_t010_r<N>{s,ns}ext/  (Phase 2 chain)
    Step 2: load any of those directories into Mathematica below.
    ─────────────────────────────────────────────────────────────── *)


(* ── 0. Setup ─────────────────────────────────────────────────── *)
SetDirectory[NotebookDirectory[]]    (* or: SetDirectory["<your path>"] *)
<< "SUGRACatalog.m"


(* ── 1. Load a catalog ────────────────────────────────────────── *)

(* (a) Phase 1 catalog: single-curve externals (e6, e7, e8, su3, ...) *)
entries = LoadCatalog["cat_1ext_nosu2_t010"];
CatalogSummary[entries]

(* (b) NHC ext catalog: full NHC chain externals (nhc_2_3, nhc_2_3_2, ...) *)
entriesNHC = LoadCatalog["cat_nhc_ext_t010"];
CatalogSummary[entriesNHC]

(* (c) Phase 2 chain rounds: e.g. round 4 SUGRA output *)
(* entriesR4 = LoadCatalog["cat_nhc_ext_nhc_ext_unified_t010_r4sext"]; *)


(* ── 2. Inspect entries ───────────────────────────────────────── *)

(* List of entries matching a spec/combo tag *)
e8List = FindEntries[entries, "e8"];     (* e8 external *)
Length[e8List]

su2List = FindEntries[entries, "su2"];   (* su2 external *)
nhc23 = FindEntries[entriesNHC, "nhc_2_3"];

(* Display a single entry: quiver, IF, anomaly numbers, fiber class *)
ShowEntry[e8List[[1]]]

(* Accessors *)
IF      = EntryIF[e8List[[1]]];        (* intersection-form matrix *)
extIdx  = ExtIndices[e8List[[1]]];     (* external curve indices (0-indexed) *)
e8List[[1]]["catalogId"]                (* LST base catalog_id from unified.cat *)
e8List[[1]]["catalogType"]              (* "NN", "N1", "DM:A(n)", ... *)


(* ── 3. Geometric invariants ──────────────────────────────────── *)

(* Fiber class: positive integer null vector of IF (if unimodular: returns {}) *)
FiberClass[IF]

(* T_min / T_max from fiber-class / anomaly constraints *)
TMin[e8List[[3]]]      (* $Failed if IF is unimodular *)
TMax[e8List[[3]]]      (* anomaly bound on T *)

(* Determinant + signature already pre-computed *)
e8List[[1]]["det"]
{e8List[[1]]["sigPos"], e8List[[1]]["sigNeg"], e8List[[1]]["sigZero"]}

(* Non-(-1) connected components (= candidate NHC clusters) *)
Clusters[IF]


(* ── 4. Anomaly cancellation check ────────────────────────────── *)

(* H_charged, V, H_neutral come pre-computed from the C++ pipeline *)
e8List[[1]]["Hc"]
e8List[[1]]["V"]
e8List[[1]]["Hn"]   (* should be 273 - 29 T + V - H_c, must be >= 0 *)

(* Δ = H_charged - V + 29 T : matches the Hamada comparison key *)
delta = e8List[[1]]["Hc"] - e8List[[1]]["V"] + 29 e8List[[1]]["T"]


(* ── 5. Filtering ─────────────────────────────────────────────── *)

(* Pick entries whose |det| is a perfect square (geometric constraint) *)
SelectDetSquare[e8List]

(* Manual filter: e.g. non-unimodular su2-ext entries *)
Select[FindEntries[entries, "su2"], Abs[#["det"]] > 1 &]


(* ── 6. Blowdown (for displaying the unblown base) ────────────── *)

BlowdownEntry[e8List[[1]]]               (* uni: T fixed *)
BlowdownEntry[e8List[[3]]]               (* non-uni: T_min auto *)
BlowdownEntry[e8List[[3]], 12]           (* override T *)


(* ── 7. Gluing two phase1 entries via a shared external ───────── *)

a = e8List[[1]]; b = e8List[[3]];
glued = GlueEntries[a, b];               (* combined IF + metadata *)


(* ── 8. Manual IF input ───────────────────────────────────────── *)

(* You can copy-paste an IF matrix that ShowEntry printed and play with it: *)
myIF = {{-1, 1, 1}, {1, -1, 0}, {1, 0, -12}};
FiberClass[myIF]
TMin[myIF]
Clusters[myIF]


(* ── 9. Full function reference ───────────────────────────────── *)

? SUGRACatalog`*

(* Or per-function: *)
? LoadCatalog
? ShowEntry
? FiberClass
? TMin
? Clusters
? GlueIFs
