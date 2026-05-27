(* ::Package:: *)

(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  SUGRA Base Construction Guide                                 *)
(*  Step 1: C++ pipeline\:c73c\:b85c cat_1ext/ \:c0dd\:c131                       *)
(*  Step 2: \:c774 \:d30c\:c77c\:c744 Mathematica\:c5d0\:c11c \:c140 \:b2e8\:c704\:b85c \:c2e4\:d589               *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)


(* \[HorizontalLine]\[HorizontalLine] 0. Setup \[HorizontalLine]\[HorizontalLine] *)
SetDirectory["/Users/seongminjeon/Desktop/sugra_pipeline"]     
<< "SUGRACatalog.m"
entries = LoadCatalog["cat_1ext"];
CatalogSummary[entries]


(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  PART 1: Entry \:d0d0\:c0c9                                            *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* \[HorizontalLine]\[HorizontalLine] 1.1 External type\:bcc4 \:bcf4\:ae30 \[HorizontalLine]\[HorizontalLine] *)
e8s = FindEntries[entries, "e8"];
su2s = FindEntries[entries, "su2"];
su2n3s = FindEntries[entries, "su2n3"];

(* \[HorizontalLine]\[HorizontalLine] 1.2 Entry \:c0c1\:c138 \[HorizontalLine]\[HorizontalLine] *)
ShowEntry[e8s[[1]]]
(* \:cd9c\:b825:
   e8  |  DM:A(2)  |  T=2
   det=11  sig={1,2,0}  \[CapitalDelta]=-215  T_max=17
   Externals: {e8} @ curves {2}
   Fiber class: {1,1,0}  (ext @ {{2, \[RightArrow] f\[CenterDot]ext=0}})
   [Graph]
   IF (copy-pasteable):
   {{-1, 1, 1}, {1, -1, 0}, {1, 0, -12}}
*)

(* \[HorizontalLine]\[HorizontalLine] 1.3 IF, fiber class, T range \[HorizontalLine]\[HorizontalLine] *)
IF = EntryIF[e8s[[1]]];
FiberClass[IF]                 (* \:c815\:c218 null vector *)
TMin[e8s[[3]]]                 (* non-uni\:c758 T_min *)
TMax[e8s[[3]]]                 (* anomaly bound *)

(* \[HorizontalLine]\[HorizontalLine] 1.4 Non-unimodular entry \:baa8\:c544\:bcf4\:ae30 \[HorizontalLine]\[HorizontalLine] *)
su2n3nu = Select[FindEntries[entries, "su2n3"], Abs[#["det"]] > 1 &];
Length[su2n3nu]

(* quiver \:cb49 \:bcf4\:ae30 *)
Do[
  Print["\[HorizontalLine]\[HorizontalLine] #", i, " T=", su2n3nu[[i]]["T"],
    " det=", su2n3nu[[i]]["det"], " \[HorizontalLine]\[HorizontalLine]"];
  Print[EntryQuiver[su2n3nu[[i]]]],
  {i, Min[10, Length[su2n3nu]]}
]

(* \[HorizontalLine]\[HorizontalLine] 1.5 Blowdown \[HorizontalLine]\[HorizontalLine] *)
BlowdownEntry[e8s[[1]]]              (* uni: T \:c790\:b3d9 *)
BlowdownEntry[su2n3nu[[1]]]          (* non-uni: T_min \:c790\:b3d9 *)
BlowdownEntry[su2n3nu[[1]], 12]      (* T=12 \:c9c0\:c815 *)


(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  PART 2: \:b450 LST\:b97c \:bd99\:c5ec\:c11c SUGRA base \:b9cc\:b4e4\:ae30                    *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* \[HorizontalLine]\[HorizontalLine] 2.1 \:ac19\:c740 spec\:c5d0\:c11c \:b450 entry \:c120\:d0dd \[HorizontalLine]\[HorizontalLine] *)

(* \:c608: e8 external *)
a = e8s[[1]];
b = e8s[[3]];

Print["Entry A:"]; ShowEntry[a];
Print["Entry B:"]; ShowEntry[b];

(* \[HorizontalLine]\[HorizontalLine] 2.2 Glue \[HorizontalLine]\[HorizontalLine] *)
glued = GlueEntries[a, b];
(* \:cd9c\:b825:
   Glued: 3+5-1 = 7 curves
   A: T=2 DM:A(2)  B: T=4 DM:A(4)
   Shared ext: e8 @ curve 2
   det=0  sig={1,5,1}
   Fiber class: {1,1,0,-1,-1,-1,-1}
   ...
*)

(* \[HorizontalLine]\[HorizontalLine] 2.3 Anomaly \:acc4\:c0b0 \[HorizontalLine]\[HorizontalLine] *)
(* external tag\:c744 \:b118\:aca8\:c57c enhancement\:ac00 \:bc18\:c601\:b428 *)
ComputeAnomaly[glued, {"e8"}]
(* \:cd9c\:b825:
   NHC Clusters:
     {2} {-12} \[RightArrow] -12  H=0 V=248
   External enhancements:
     e8 \[RightArrow] in NHC  dH=0 dV=0
   Result: T=5  H_c=0  V=248  H_n=76
   Anomaly OK
*)


(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  PART 3: T \:c870\:ac74 \:d655\:c778                                           *)
(*                                                                *)
(*  \:d575\:c2ec \:c870\:ac74:                                                    *)
(*    glued IF\:c758 negative eigenvalue \:ac1c\:c218 (= T_glued)\:ac00           *)
(*    T_min(A) + T_min(B) - 1 \:c774\:c0c1\:c774\:c5b4\:c57c \:d568.                      *)
(*    (-1 \:c740 shared external\:b85c \:d569\:ccd0\:c9c4 fiber class \:bc29\:d5a5)            *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* \[HorizontalLine]\[HorizontalLine] 3.1 T_glued \:d655\:c778 \[HorizontalLine]\[HorizontalLine] *)
eigsGlued = Eigenvalues[N@glued];
Tglued = Count[eigsGlued, _?(# < -10^-8 &)];
Print["T_glued = ", Tglued];

(* \[HorizontalLine]\[HorizontalLine] 3.2 \:ac01 entry\:c758 T_min \:d655\:c778 \[HorizontalLine]\[HorizontalLine] *)
tminA = TMin[a];    (* uni\:ba74 $Failed \[RightArrow] T = sig_neg \:adf8 \:c790\:ccb4 *)
tminB = TMin[b];

(* uni entry\:b294 T_min = T (\:ace0\:c815) *)
If[tminA === $Failed, tminA = a["T"]];
If[tminB === $Failed, tminB = b["T"]];

Print["T_min(A) = ", tminA];
Print["T_min(B) = ", tminB];
Print["T_min(A) + T_min(B) - 1 = ", tminA + tminB - 1];
Print["T_glued = ", Tglued];

If[Tglued >= tminA + tminB - 1,
  Print[Style["\[Checkmark] T condition satisfied", Darker[Green], Bold]],
  Print[Style["\:2717 T too small", Red, Bold]]
]

(* \[HorizontalLine]\[HorizontalLine] 3.3 Extended IF + blowdown \[HorizontalLine]\[HorizontalLine] *)
(* T\:b97c \:c9c1\:c811 \:c9c0\:c815\:d574\:c11c extended IF\:b97c \:ad6c\:c131 *)
Ttarget = tminA + tminB - 1;     (* \:cd5c\:c18c T: shared external\:b85c -1 *)
ext = ExtendedIF[glued, Ttarget];
{bd, nBlown} = Blowdown[ext];
Print["Extended IF @ T=", Ttarget, ":"];
Print["  blown down: ", nBlown, " curves"];
Print["  sig = ", Module[{ev = Eigenvalues[N@bd]},
  {Count[ev,_?(#>10^-8&)], Count[ev,_?(#<-10^-8&)], Count[ev,_?(-10^-8<=#<=10^-8&)]}]];
Print["  det = ", Det[bd]];
MatrixForm[bd]


(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  PART 4: Systematic scan \[LongDash] \:ac19\:c740 spec\:c758 \:baa8\:b4e0 pair \:c2dc\:b3c4          *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* \:c608: e8 entry\:b4e4\:c758 \:baa8\:b4e0 pair\:c5d0\:c11c anomaly OK\:c778 \:ac83 \:cc3e\:ae30 *)
Print["=== Scanning e8 pairs ==="];
Do[
  Module[{ea = e8s[[i]], eb = e8s[[j]], g, eigs, Tg, ta, tb},
    g = GlueIFs[EntryIF[ea], ExtIndices[ea][[1]],
                EntryIF[eb], ExtIndices[eb][[1]]];
    eigs = Eigenvalues[N@g];
    Tg = Count[eigs, _?(# < -10^-8 &)];
    ta = If[TMin[ea] === $Failed, ea["T"], TMin[ea]];
    tb = If[TMin[eb] === $Failed, eb["T"], TMin[eb]];
    If[Tg >= ta + tb - 1,
      (* Anomaly check *)
      Module[{cls = Clusters[g], totalH = 0, totalV = 0, Hn},
        Do[Module[{nhc = NHCLookup[cl["subIF"]]},
          If[nhc["H"] >= 0, totalH += nhc["H"]; totalV += nhc["V"]]],
          {cl, cls}];
        Hn = 273 + totalV - totalH - 29*Tg;
        If[Hn >= 0,
          Print["  (", i, ",", j, ") T_A=", ea["T"], " T_B=", eb["T"],
            " T_glued=", Tg, " H_n=", Hn, Style[" \[Checkmark]", Darker[Green]]]
        ]
      ]
    ]
  ],
  {i, Length[e8s]}, {j, i, Length[e8s]}
]


(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  PART 5: \:b2e4\:b978 spec\:b4e4\:b3c4 \:c2dc\:b3c4                                    *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)

(* su2n3 non-uni pair *)
Print["=== su2n3 non-uni example ==="];
a2 = su2n3nu[[1]];
b2 = su2n3nu[[2]];
ShowEntry[a2];
ShowEntry[b2];

glued2 = GlueEntries[a2, b2];
ComputeAnomaly[glued2, {"su2n3"}]

(* T condition *)
Tg2 = Count[Eigenvalues[N@glued2], _?(# < -10^-8 &)];
Print["T_glued = ", Tg2];
Print["T_min(A) = ", TMin[a2], "  T_min(B) = ", TMin[b2]];


(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)
(*  Quick Reference                                               *)
(*                                                                *)
(*  LoadCatalog["cat_1ext"]       \:ce74\:d0c8\:b85c\:adf8 \:b85c\:b4dc                    *)
(*  CatalogSummary[entries]       \:d1b5\:acc4                            *)
(*  FindEntries[entries, "e8"]    external type\:c73c\:b85c \:d544\:d130           *)
(*  ShowEntry[e]                  \:c0c1\:c138 (copy-pasteable IF)        *)
(*  EntryIF[e]                    IF matrix                       *)
(*  ExtIndices[e]                 external curve index (0-idx)    *)
(*  FiberClass[IF]                fiber class (\:c815\:c218 null vector)  *)
(*  TMin[e], TMax[e]              T range                         *)
(*  GlueEntries[a, b]             \:b450 entry\:c758 external \:d569\:ce68       *)
(*  ComputeAnomaly[IF, {"tag"}]   NHC + enhancement \[RightArrow] anomaly    *)
(*  ExtendedIF[IF, T]             IF + b\:2080                        *)
(*  Blowdown[extIF]               exceptional curve blowdown     *)
(*  BlowdownEntry[e]              entry\:c5d0\:c11c \:d55c\:bc29\:c5d0                *)
(* \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550 *)


f4s = FindEntries[entries, "f4"];
Select[Union[entries[[All, "combo"]]], StringContainsQ[#, "f4"] &]

Length[f4s]


f4s = FindEntries[entries, "f4"];
fa = f4s[[1]]; fb = f4s[[1]];
g = GlueIFs[EntryIF[fa], ExtIndices[fa][[1]], EntryIF[fb], ExtIndices[fb][[1]]];
Tg = Count[Eigenvalues[N@g], _?(# < -10^-8 &)];
ta = If[TMin[fa] === $Failed, fa["T"], TMin[fa]];
tb = If[TMin[fb] === $Failed, fb["T"], TMin[fb]];
Print["Tg=", Tg, " ta=", ta, " tb=", tb];


f4s = FindEntries[entries, "f4"];
Print["=== Scanning f4 pairs ==="];
Do[
  Module[{fa = f4s[[i]], fb = f4s[[j]], g, Tg, ta, tb},
    g = GlueIFs[EntryIF[fa], ExtIndices[fa][[1]],
                EntryIF[fb], ExtIndices[fb][[1]]];
    Tg = Count[Eigenvalues[N@g], _?(# < -10^-8 &)];
    ta = If[TMin[fa] === $Failed, fa["T"], TMin[fa]];
    tb = If[TMin[fb] === $Failed, fb["T"], TMin[fb]];
    If[Tg >= ta + tb - 1,
      Module[{cls = Clusters[g], totalH = 0, totalV = 0, Hn},
        Do[Module[{nhc = NHCLookup[cl["subIF"]]},
          If[nhc["H"] >= 0, totalH += nhc["H"]; totalV += nhc["V"]]],
          {cl, cls}];
        Hn = 273 + totalV - totalH - 29*Tg;
        If[Hn >= 0,
          Print["(", i, ",", j, ") T_A=", fa["T"], " T_B=", fb["T"],
            " Tg=", Tg, " Hn=", Hn]
        ]
      ]
    ]
  ],
  {i, Length[f4s]}, {j, i, Length[f4s]}
]


f4s = FindEntries[entries, "f4"];
Print["=== Scanning f4 pairs ==="];
Do[
  Module[{fa = f4s[[i]], fb = f4s[[j]], g, Tg, ta, tb},
    g = GlueIFs[EntryIF[fa], ExtIndices[fa][[1]],
                EntryIF[fb], ExtIndices[fb][[1]]];
    Tg = Count[Eigenvalues[N@g], _?(# < -10^-8 &)];
    ta = If[TMin[fa] === $Failed, fa["T"], TMin[fa]];
    tb = If[TMin[fb] === $Failed, fb["T"], TMin[fb]];
    If[Tg >= ta + tb - 1,
      Module[{cls = Clusters[g], totalH = 0, totalV = 0, Hn,
              ext, bd, n, sigExt, ccOK = True, nn = Length[g]},
        (* NHC *)
        Do[Module[{nhc = NHCLookup[cl["subIF"]]},
          If[nhc["H"] >= 0, totalH += nhc["H"]; totalV += nhc["V"]]],
          {cl, cls}];
        Hn = 273 + totalV - totalH - 29*Tg;
        If[Hn < 0, Continue[]];
        
        (* CC budget: each (-1) curve, sum cc of non-(-1) neighbors \[LessEqual] 8 *)
        Do[If[g[[k,k]] == -1,
          Module[{cc = 0.},
            Do[If[m != k && g[[k,m]] != 0 && g[[m,m]] != -1,
              Module[{si = g[[m,m]], intN = Abs[g[[k,m]]], dim, hdual},
                {dim, hdual} = Switch[si,
                  -2, {3, 2}, -3, {8, 3}, -4, {28, 6}, -5, {52, 9},
                  -6, {78, 12}, -7, {133, 18}, -8, {133, 18}, -12, {248, 30},
                  _, {0, 1}];
                If[dim > 0, cc += intN * dim / (intN + hdual)]]],
              {m, nn}];
            If[cc > 8.01, ccOK = False]]],
          {k, nn}];
        If[!ccOK, Continue[]];
        
        (* Extended IF sig check *)
        ext = ExtendedIF[g, Tg];
        {bd, n} = Blowdown[ext];
        sigExt = Module[{ev = Eigenvalues[N@bd]},
          {Count[ev, _?(# > 10^-8 &)], 
           Count[ev, _?(# < -10^-8 &)], 
           Count[ev, _?(-10^-8 <= # <= 10^-8 &)]}];
        If[sigExt[[1]] != 1, Continue[]];
        
        Print["(", i, ",", j, ") T_A=", fa["T"], " T_B=", fb["T"],
          " Tg=", Tg, " Hn=", Hn, " sig=", sigExt, " bd: ", MatrixForm[bd]]
      ]
    ]
  ],
  {i, Length[f4s]}, {j, i, Length[f4s]}
]



