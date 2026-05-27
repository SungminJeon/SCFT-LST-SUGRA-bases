(* ::Package:: *)

(* SUGRA Catalog Tutorial \[LongDash] \:c140 \:b2e8\:c704\:b85c \:c2e4\:d589 *)

SetDirectory["/Users/seongminjeon/Desktop/sugra_pipeline"]  (* \:bcf8\:c778 \:acbd\:b85c *)
<< "SUGRACatalog.m"

(* \[HorizontalLine]\[HorizontalLine] \:b85c\:b4dc \[HorizontalLine]\[HorizontalLine] *)
entries = LoadCatalog["cat_1ext"];
CatalogSummary[entries]

(* \[HorizontalLine]\[HorizontalLine] \:cc3e\:ae30 \[HorizontalLine]\[HorizontalLine] *)
e8s = FindEntries[entries, "e8"];
ShowEntry[e8s[[1]]]
(* \[RightArrow] Fiber class: {1,1,0}  (ext @ {{2, \[RightArrow] f\[CenterDot]ext=0}})  *)
(* \[RightArrow] \:c774 \:acbd\:c6b0 external(curve 2)\:c740 fiber\:c640 intersection 0 *)

(* \[HorizontalLine]\[HorizontalLine] Fiber class \:c9c1\:c811 \:b2e4\:b8e8\:ae30 \[HorizontalLine]\[HorizontalLine] *)
IF = EntryIF[e8s[[1]]];
fc = FiberClass[IF]             (* \:c815\:c218 null vector *)
fc = FiberClass[e8s[[1]]]      (* entry \:b118\:aca8\:b3c4 \:b428 *)

(* external\:c774 \:bd99\:b294 target curve\:c758 fiber multiplicity *)
extIdx = ExtIndices[e8s[[1]]][[1]];   (* 0-indexed *)
fc[[extIdx + 1]]                       (* fiber\[CenterDot]external *)

(* \[HorizontalLine]\[HorizontalLine] non-uni\:b3c4 \:b3d9\:c77c \[HorizontalLine]\[HorizontalLine] *)
ShowEntry[e8s[[3]]]
FiberClass[e8s[[3]]]
TMin[e8s[[3]]]
TMax[e8s[[3]]]

(* \[HorizontalLine]\[HorizontalLine] IF \:bcf5\:bd99 \[HorizontalLine]\[HorizontalLine] *)
(* ShowEntry\:ac00 \:cd9c\:b825\:d55c IF\:b97c \:bcf5\:c0ac\:d574\:c11c: *)
(* myIF = {{-1,1,1},{1,-1,0},{1,0,-12}}; *)
(* FiberClass[myIF] *)

(* \[HorizontalLine]\[HorizontalLine] Blowdown (\:d544\:c694\:d558\:ba74) \[HorizontalLine]\[HorizontalLine] *)
BlowdownEntry[e8s[[1]]]
BlowdownEntry[e8s[[3]]]        (* non-uni: T_min\:c5d0\:c11c \:c790\:b3d9 *)
BlowdownEntry[e8s[[1]], 5]     (* \:d2b9\:c815 T *)


su2s = FindEntries[entries,"su2"];
ShowEntry[su2s[[1054]]]
ShowEntry[su2s[[1053]]]
hat1m1 = FindEntries[entries,"hat1m1"];
ShowEntry[hat1m1[[1]]]

ExtIndices[su2s[[1]]]
(* \:c608: {2} \:ac00 \:b098\:c640\:c57c \:d568 *)

(* \:c9c1\:c811 \:d655\:c778 *)
su2s[[1]]["externals"]



Select[FindEntries[entries, "su2n3"], Abs[#["det"]] > 1 &]


su2n3nu = Select[FindEntries[entries, "su2n3"], Abs[#["det"]] > 1 &];
Length[su2n3nu]

(* \:c804\:bd80 \:bcf4\:ae30 *)
Do[
  Print["\[HorizontalLine]\[HorizontalLine] #", i, " T=", su2n3nu[[i]]["T"], " det=", su2n3nu[[i]]["det"], " \[HorizontalLine]\[HorizontalLine]"];
  Print[EntryQuiver[su2n3nu[[i]]]],
  {i, Length[su2n3nu]}
]

su2snu = Select[FindEntries[entries,"su2"], Abs[#["det"]] > 1 &];

Do[
  Print["\[HorizontalLine]\[HorizontalLine] #", i, " T=", su2snu[[i]]["T"], " det=", su2snu[[i]]["det"], " \[HorizontalLine]\[HorizontalLine]"];
  Print[EntryQuiver[su2snu[[i]]]],
  {i, Length[su2snu]}
]






(* 1. IF\:b07c\:b9ac glue *)(*5,5 \[Rule] ok *) 
glued = GlueEntries[su2n3nu[[5]], su2snu[[5]]];

(* 2. glued IF\:b85c extended \:b9cc\:b4e4\:ae30 *)
ext = ExtendedIF[glued, 8];    (* \:c6d0\:d558\:b294 T *)

(* 3. blowdown *)
{bd, n} = Blowdown[ext];
MatrixForm[bd]
Det[bd]

ShowEntry[su2n3nu[[5]]] 
ShowEntry[su2snu[[5]]]
ShowEntry[so8n2[[1]]]


(* ::Print:: *)
(*Graph[{1, 2, 3, 4, 5, 6}, {Null, {{1, 2}, {1, 4}, {1, 5}, {2, 3}, {3, 6}}}, {FormatType -> TraditionalForm, ImageSize -> 350, EdgeStyle -> {GrayLevel[0.5]}, GraphLayout -> "SpringElectricalEmbedding", ImageSize -> 350, VertexLabels -> {6 -> Placed[Style[2, GrayLevel[1], Bold, 16], Center], 3 -> Placed[Style[1, GrayLevel[0], Bold, 16], Center], 5 -> Placed[Style[2, GrayLevel[0], Bold, 16], Center], 4 -> Placed[Style[2, GrayLevel[0], Bold, 16], Center], 2 -> Placed[Style[2, GrayLevel[0], Bold, 16], Center], 1 -> Placed[Style[2, GrayLevel[0], Bold, 16], Center]}, VertexShapeFunction -> {5 -> "Circle", 6 -> "Diamond", 1 -> "Circle", 4 -> "Circle", 2 -> "Circle", 3 -> "Circle"}, VertexSize -> {0.6}, VertexStyle -> {6 -> RGBColor[1, 0, 0], 1 -> RGBColor[0.87, 0.94, 1], 2 -> RGBColor[0.87, 0.94, 1], 3 -> RGBColor[0.87, 0.94, 1], 4 -> RGBColor[0.87, 0.94, 1], 5 -> RGBColor[0.87, 0.94, 1]}}]c*)


(* 1. IF\:b07c\:b9ac glue *)
glued = GlueEntries[su2n3nu[[39]], su2n3nu[[19]]];

(* 2. glued IF\:b85c extended \:b9cc\:b4e4\:ae30 *)
ext = ExtendedIF[glued, 8];    (* \:c6d0\:d558\:b294 T *)

(* 3. blowdown *)
{bd, n} = Blowdown[ext];
MatrixForm[bd]
Det[bd]

ShowEntry[su2n3nu[[39]]] 
ShowEntry[su2n3nu[[19]]]



