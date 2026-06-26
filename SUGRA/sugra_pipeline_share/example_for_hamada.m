(* ============================================================ *)
(*  example_for_hamada.m                                        *)
(*  6D SUGRA base selection by T_H, T_min, and external types.  *)
(*  Run cell-by-cell in a Mathematica front-end.                *)
(* ============================================================ *)

(* ─── Setup ─────────────────────────────────────────────────── *)
(* Run from the folder that contains SUGRACatalog.m + catalog_analysis.m.
   Works for both notebook (.nb) and script (wolframscript) use.       *)
SetDirectory[
  If[StringQ[$InputFileName] && $InputFileName =!= "",
     DirectoryName[$InputFileName],
     NotebookDirectory[]]
];
<< "SUGRACatalog.m";
<< "catalog_analysis.m";


(* ─── 1. Load all SUGRA bases for T = 0..10 ────────────────── *)
(* A run of the pipeline produces one top-level directory
   T_<MIN>_<MAX>/ containing several sub-catalogs:
     cat_phase1_nhc_merged_<sfx>      <- 1-external bases (incl. NHC ext)
     cat_nhc_ext_nhc_ext_unified_<sfx>_r2sext   <- 2 externals
     ...                                          ...
     cat_nhc_ext_nhc_ext_unified_<sfx>_r9sext   <- 9 externals
   The helper below loads them all and merges into one list.    *)

LoadAllSUGRA[root_String] := Flatten[
  LoadCatalog /@ Join[
    FileNames["cat_phase1_nhc_merged_*"      , root],
    FileNames["cat_nhc_ext_nhc_ext_unified_*sext", root]
  ]];

entries = LoadAllSUGRA["T_0_10"];
Length[entries]
(* Each entry is an Association with the following keys:        *)
(*    id, combo, externals, catalogId, catalogType,             *)
(*    baseT, T, Hc, V, Hn, det, sigPos, sigNeg, sigZero,        *)
(*    IF, remaining                                              *)


(* ─── 2. Field meanings ────────────────────────────────────── *)
(*   baseT   = T_H   ( tensor count of the underlying LST )      *)
(*   T       = full sig_neg of the extended base                 *)
(*   TMin[e] = T_min ( minimum T that admits a SUGRA b0 lift )    *)
(*   externals = list of <|"tag" -> ..., "extSI" -> ...,         *)
(*                          "targetSI" -> ..., "intNum" -> ...,  *)
(*                          "isHat1" -> ..., "curveIdx" -> ...|> *)


(* ─── 3. Convenience selectors ─────────────────────────────── *)
(*  Helper: list of external tags attached to an entry.          *)
ExtTags[e_]  := e["externals"][[All, "tag"]];

(*  Helper: multiset (sorted) of external tags.                  *)
ExtTagSet[e_] := Sort[ExtTags[e]];


(* ─── 4. Pick by T_H ───────────────────────────────────────── *)
(*  All bases whose underlying LST has T_H = 3.                  *)
sel = Select[entries, #["baseT"] == 3 &];
Length[sel]


(* ─── 5. Pick by T_min ─────────────────────────────────────── *)
(*  Bases with T_min == 7. (TMin returns $Failed for unimodular  *)
(*  bases, hence the IntegerQ guard.)                            *)
sel = Select[entries, IntegerQ[TMin[#]] && TMin[#] == 7 &];
Length[sel]


(* ─── 6. Pick by external types ────────────────────────────── *)
(*  6a) Combo string match (fastest, requires exact ordering).   *)
sel = FindEntries[entries, "e6+su3"];
Length[sel]

(*  6b) Multiset match (order-independent).                      *)
target = Sort[{"e6", "su3"}];
sel = Select[entries, ExtTagSet[#] == target &];
Length[sel]

(*  6c) Contains at least one e8 external.                       *)
sel = Select[entries, MemberQ[ExtTags[#], "e8"] &];
Length[sel]

(*  6d) Bases that contain exactly two externals, both e6.       *)
sel = Select[entries, ExtTagSet[#] == {"e6", "e6"} &];
Length[sel]


(* ─── 7. Combined filter (T_H + T_min + externals) ─────────── *)
(*  Example: T_H = 5, T_min = 8, externals = { e7, e7, su3 }.    *)
target = Sort[{"e7", "e7", "su3"}];

sel = Select[entries,
        #["baseT"] == 5 &&
        IntegerQ[TMin[#]] && TMin[#] == 8 &&
        ExtTagSet[#] == target &
      ];

Length[sel]
ShowEntry /@ sel       (* pretty-print each hit *)


(* ─── 8. Reusable function form ────────────────────────────── *)
(*  Wrap the combined filter so you can call it ad hoc.          *)
PickBases[entries_, TH_:Automatic, Tmin_:Automatic, exts_:Automatic] :=
  Select[entries,
    And @@ {
      TH    === Automatic || #["baseT"] == TH,
      Tmin  === Automatic || (IntegerQ[TMin[#]] && TMin[#] == Tmin),
      exts  === Automatic || ExtTagSet[#] == Sort[exts]
    } &
  ];

(*  Usage examples:                                              *)
hits = PickBases[entries, 3];                       (* T_H = 3 *)
hits = PickBases[entries, 3, 6];                    (* T_H = 3, T_min = 6 *)
hits = PickBases[entries, 3, 6, {"e6", "su3"}];     (* + externals = {e6,su3} *)
hits = PickBases[entries, Automatic, 8, {"e8"}];    (* T_min = 8, single e8 *)


(* ─── 9. Tabulate results ──────────────────────────────────── *)
TableForm[
  {#["id"], #["catalogType"], #["baseT"], TMin[#], #["combo"], #["T"]} & /@
    PickBases[entries, 3, Automatic, {"e6", "su3"}],
  TableHeadings -> {None,
    {"id", "LST", "T_H", "T_min", "combo", "T"}}
]


(* ─── 10. Show full entry (IF + physics + externals) ──────── *)
PrettyEntry[sel[[1]]]
ShowIF[sel[[1]]]
