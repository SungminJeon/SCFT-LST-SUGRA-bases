(* SUGRACatalog.m — Minimal SUGRA catalog tools *)

BeginPackage["SUGRACatalog`"];

LoadCatalog::usage = "LoadCatalog[dir] loads all .cat files.";
EntryIF::usage = "EntryIF[entry] returns IF matrix.";
EntryQuiver::usage = "EntryQuiver[entry] returns Graph.";
ExtIndices::usage = "ExtIndices[entry] returns external curve indices (0-indexed).";
ShowEntry::usage = "ShowEntry[entry] displays entry with copy-pasteable IF.";
FindEntries::usage = "FindEntries[entries, spec] or FindEntries[entries, spec, T].";
DetSquareQ::usage = "DetSquareQ[entry] returns True if |det| is a perfect square.";
SelectDetSquare::usage = "SelectDetSquare[entries] filters entries with |det|=n².";
CatalogSummary::usage = "CatalogSummary[entries] prints stats.";
TMin::usage = "TMin[entry] or TMin[IF] returns T_min (T_crit). $Failed if unimodular.";
TMax::usage = "TMax[entry] returns T_max from anomaly bound.";
BlowdownEntry::usage = "BlowdownEntry[entry] or BlowdownEntry[entry, T] blows down and shows labeled result.";
FiberClass::usage = "FiberClass[IF] returns integer null vector (fiber class). {} if full rank.";
ExtendedIF::usage = "ExtendedIF[IF] or ExtendedIF[IF, T] builds IF+b0. Auto T from eigenvalues if omitted.";
Blowdown::usage = "Blowdown[IF] blows down all (-1) curves. Returns {bdIF, nBlown}.";
GlueEntries::usage = "GlueEntries[a, b] glues two entries by identifying their external curves.";
GlueIFs::usage = "GlueIFs[IF1, ep1, IF2, ep2] glues two IFs at endpoints (0-indexed). Returns matrix.";
Clusters::usage = "Clusters[IF] returns non-(-1) connected components.";
ComputeAnomaly::usage = "ComputeAnomaly[IF, Hc, V] or ComputeAnomaly[IF] or ComputeAnomaly[IF, {tags}].";
NHCLookup::usage = "NHCLookup[subIF] matches submatrix to NHC table. Returns <|name, H, V|>.";

Begin["`Private`"];

(* ── Parser ── *)
parseCatFile[text_String] := Module[
  {lines, entries={}, cur, inEntry=False, ifLeft=0, ifRow=0, ifN=0,
   remLeft=0, extLeft=0, mat, remaining, externals,
   bdLeft=0, bdRow=0, bdN=0, bdMat},
  lines = StringSplit[text, "\n"];
  Do[Which[
    StringMatchQ[line,"ENTRY *"],
      cur=<|"id"->ToExpression@StringDrop[line,6],"combo"->"","externals"->{},"catalogId"->0,
        "catalogType"->"","baseT"->0,"T"->0,"Hc"->0,"V"->0,"Hn"->0,"det"->0,
        "sigPos"->0,"sigNeg"->0,"sigZero"->0,"IF"->{},"remaining"->{},
        "glueA"->"","glueB"->"","glueTA"->0,"glueTB"->0,
        "blowdown"->{},"bdSig"->{}|>;
      externals={};remaining={};inEntry=True;ifLeft=0;remLeft=0;extLeft=0,
    inEntry&&StringMatchQ[line,"SPEC *"],
      Module[{p=StringSplit[line]},If[Length[p]>=7,cur["combo"]=p[[3]];
        AppendTo[externals,<|"curveIdx"->-1,"specId"->ToExpression@p[[2]],"tag"->p[[3]],
          "extSI"->ToExpression@p[[4]],"targetSI"->ToExpression@p[[5]],
          "intNum"->ToExpression@p[[6]],"isHat1"->(p[[7]]=="1")|>]]],
    inEntry&&StringMatchQ[line,"COMBO *"],
      cur["combo"]=StringTrim@StringDrop[line,6],
    inEntry&&StringMatchQ[line,"EXTERNALS *"],
      extLeft=ToExpression@StringDrop[line,10];externals={},
    inEntry&&extLeft>0,
      Module[{p=StringSplit[line]},If[Length[p]>=7,
        AppendTo[externals,<|"curveIdx"->ToExpression@p[[1]],"specId"->ToExpression@p[[2]],
          "tag"->p[[3]],"extSI"->ToExpression@p[[4]],"targetSI"->ToExpression@p[[5]],
          "intNum"->ToExpression@p[[6]],"isHat1"->(p[[7]]=="1")|>]];extLeft--],
    inEntry&&StringMatchQ[line,"BASE *"],
      Module[{p=StringSplit[line]},If[Length[p]>=5,
        cur["catalogId"]=ToExpression@p[[2]];cur["catalogType"]=p[[3]];cur["baseT"]=ToExpression@p[[4]]]],
    inEntry&&StringMatchQ[line,"ATTACH *"],
      Module[{p=StringSplit[line]},If[Length[p]>=3&&Length[externals]>0&&externals[[-1]]["curveIdx"]==-1,
        externals[[-1]] = Append[KeyDrop[externals[[-1]], "curveIdx"], "curveIdx" -> ToExpression@p[[3]]]]],
    inEntry&&StringMatchQ[line,"PHYSICS *"],
      Module[{p=StringSplit[line]},If[Length[p]>=9,
        cur["T"]=ToExpression@p[[2]];cur["Hc"]=ToExpression@p[[3]];cur["V"]=ToExpression@p[[4]];
        cur["Hn"]=ToExpression@p[[5]];cur["det"]=ToExpression@p[[6]];cur["sigPos"]=ToExpression@p[[7]];
        cur["sigNeg"]=ToExpression@p[[8]];cur["sigZero"]=ToExpression@p[[9]]]],
    inEntry&&StringMatchQ[line,"IF *"],
      ifN=ToExpression@StringDrop[line,3];mat=ConstantArray[0,{ifN,ifN}];ifLeft=ifN;ifRow=1,
    inEntry&&ifLeft>0,
      Module[{v=ToExpression/@StringSplit[line]},If[Length[v]==ifN,mat[[ifRow]]=v;ifRow++;ifLeft--]],
    inEntry&&StringMatchQ[line,"GLUE *"],
      Module[{p=StringSplit[line]},If[Length[p]>=5,
        cur["glueA"]=p[[2]];cur["glueTA"]=ToExpression@p[[3]];
        cur["glueB"]=p[[4]];cur["glueTB"]=ToExpression@p[[5]]]],
    inEntry&&StringMatchQ[line,"BLOWDOWN *"],
      bdN=ToExpression@StringDrop[line,9];bdMat=ConstantArray[0,{bdN,bdN}];bdLeft=bdN;bdRow=1,
    inEntry&&bdLeft>0,
      Module[{v=ToExpression/@StringSplit[line]},If[Length[v]==bdN,bdMat[[bdRow]]=v;bdRow++;bdLeft--]],
    inEntry&&StringMatchQ[line,"BDSIG *"],
      Module[{p=StringSplit[line]},If[Length[p]>=5,
        cur["bdSig"]={ToExpression@p[[2]],ToExpression@p[[3]],ToExpression@p[[4]],ToExpression@p[[5]]}]],
    inEntry&&StringMatchQ[line,"REMAINING *"],remLeft=ToExpression@StringDrop[line,10];remaining={},
    inEntry&&remLeft>0,
      Module[{p=StringSplit[line]},If[Length[p]>=3,
        AppendTo[remaining,<|"curveIdx"->ToExpression@p[[1]],"selfInt"->ToExpression@p[[2]],
          "availCC"->ToExpression@p[[3]]|>]];remLeft--],
    inEntry&&line=="END",
      cur["IF"]=mat;cur["remaining"]=remaining;cur["externals"]=externals;
      If[bdN>0, cur["blowdown"]=bdMat]; bdN=0;
      AppendTo[entries,cur];inEntry=False
  ],{line,lines}];entries];

LoadCatalog[dir_String] := Module[{files, all={}},
  files=FileNames["*.cat",dir];
  Do[all=Join[all,parseCatFile[Import[f,"Text"]]],{f,files}];
  Print["Loaded ",Length[all]," entries from ",dir];all];

(* ── Accessors ── *)
EntryIF[e_Association] := e["IF"];
ExtIndices[e_Association] := e["externals"][[All,"curveIdx"]];

(* ── Quiver ── *)
EntryQuiver[e_Association] := Module[{IF=e["IF"],n,edges={},labels,extSet},
  n=Length[IF];labels=Table[-IF[[i,i]],{i,n}];extSet=ExtIndices[e]+1;
  Do[If[IF[[i,j]]!=0,AppendTo[edges,i<->j]],{i,n},{j,i+1,n}];
  Graph[Range[n],edges,
    VertexLabels->Table[i->Placed[Style[labels[[i]],
      If[MemberQ[extSet,i],White,Black],Bold,16],Center],{i,n}],
    VertexSize->0.6,
    VertexStyle->Table[i->If[MemberQ[extSet,i],Red,LightBlue],{i,n}],
    VertexShapeFunction->Table[i->If[MemberQ[extSet,i],"Diamond","Circle"],{i,n}],
    EdgeStyle->{Thick,Gray},
    GraphLayout->"SpringElectricalEmbedding",
    ImageSize->350]];

(* ── IF to copy-pasteable string ── *)
ifString[IF_?MatrixQ] := Module[{n=Length[IF]},
  "{"<>StringRiffle[
    Table["{"<>StringRiffle[ToString/@IF[[i]],","]<>"}",{i,n}],
    ","
  ]<>"}"
];

(* ── T_min ── *)
TMin[IF_?MatrixQ, h1_List:{}] := Module[{n=Length[IF], M, B, A, Tc},
  M = ConstantArray[0., {n+1,n+1}]; M[[1;;n,1;;n]] = N@IF;
  Do[Module[{b=If[MemberQ[h1,i-1],-1.,IF[[i,i]]+2.]},M[[i,n+1]]=b;M[[n+1,i]]=b],{i,n}];
  M[[n+1,n+1]]=0.; B=Det[M]; M[[n+1,n+1]]=1.; A=Det[M]-B;
  If[Abs[A]<10^-10, Return[$Failed]];
  Tc = 9. - (-B/A);
  Ceiling[Tc - 10^-9]];

TMin[e_Association] := Module[{h1},
  h1 = Cases[e["externals"], _?(#["isHat1"]&)][[All,"curveIdx"]];
  TMin[e["IF"], h1]];

(* ── T_max ── *)
TMax[e_Association] := Min[Floor[(273 + e["V"] - e["Hc"])/29], 193];

(* ── ExtendedIF: add b0 row/col ── *)
ExtendedIF[IF_?MatrixQ] := Module[{n=Length[IF], T, ext},
  T = Count[Eigenvalues[N@IF], _?(# < -10^-8 &)];
  ext = ConstantArray[0, {n+1, n+1}]; ext[[1;;n, 1;;n]] = IF;
  Do[ext[[i,n+1]] = IF[[i,i]]+2; ext[[n+1,i]] = IF[[i,i]]+2, {i,n}];
  ext[[n+1,n+1]] = 9-T; ext];

ExtendedIF[IF_?MatrixQ, T_Integer] := Module[{n=Length[IF], ext},
  ext = ConstantArray[0, {n+1, n+1}]; ext[[1;;n, 1;;n]] = IF;
  Do[ext[[i,n+1]] = IF[[i,i]]+2; ext[[n+1,i]] = IF[[i,i]]+2, {i,n}];
  ext[[n+1,n+1]] = 9-T; ext];

(* with hat1: b0·C = -1 for hat1, b0·C = C²+2 for others *)
ExtendedIF[IF_?MatrixQ, T_Integer, h1_List] := Module[{n=Length[IF], ext},
  ext = ConstantArray[0, {n+1, n+1}]; ext[[1;;n, 1;;n]] = IF;
  Do[Module[{b = If[MemberQ[h1, i-1], -1, IF[[i,i]]+2]},
    ext[[i,n+1]] = b; ext[[n+1,i]] = b], {i,n}];
  ext[[n+1,n+1]] = 9-T; ext];

(* ── Blowdown: blow down exceptional curves only ── *)
(* Exceptional = C²=-1 AND b₀·C=1 (last col = 1) *)
Blowdown[IF_?MatrixQ] := Module[{cur=IF, n, m1, bd=0, lastCol},
  While[True,
    n = Length[cur]; m1 = -1; lastCol = n;
    Do[If[cur[[i,i]] == -1 && cur[[i, lastCol]] == 1,
      m1 = i; Break[]], {i, n}];
    If[m1 == -1, Break[]];
    Module[{keep = Delete[Range[n], m1]},
      cur = Table[cur[[i,j]] + cur[[i,m1]]*cur[[m1,j]], {i,keep}, {j,keep}]];
    bd++];
  {cur, bd}];


(* ── GlueEntries: identify external curves ── *)
GlueEntries[a_Association, b_Association] := Module[
  {IF1 = a["IF"], IF2 = b["IF"], ep1, ep2, n1, n2, nt, G, mb, sharedTag},
  ep1 = ExtIndices[a][[1]];
  ep2 = ExtIndices[b][[1]];
  sharedTag = a["externals"][[1]]["tag"];
  n1 = Length[IF1]; n2 = Length[IF2];
  nt = n1 + n2 - 1;
  G = ConstantArray[0, {nt, nt}];
  G[[1 ;; n1, 1 ;; n1]] = IF1;
  mb[j_] := If[j == ep2, ep1 + 1, n1 + If[j < ep2, j + 1, j]];
  Do[If[j != ep2, G[[mb[j], mb[j]]] = IF2[[j + 1, j + 1]]], {j, 0, n2 - 1}];
  Do[If[IF2[[j + 1, k + 1]] != 0,
    G[[mb[j], mb[k]]] += IF2[[j + 1, k + 1]];
    G[[mb[k], mb[j]]] += IF2[[k + 1, j + 1]]],
    {j, 0, n2 - 1}, {k, j + 1, n2 - 1}];
  Print["Glued: ", n1, "+", n2, "-1 = ", nt, " curves"];
  Print["A: T=", a["T"], " ", a["catalogType"], "  B: T=", b["T"], " ", b["catalogType"]];
  Print["Shared ext: ", sharedTag, " @ curve ", ep1];
  Print["det=", Det[G], "  sig=",
    Module[{ev = Eigenvalues[N@G]},
      {Count[ev,_?(#>10^-8&)], Count[ev,_?(#<-10^-8&)], Count[ev,_?(-10^-8<=#<=10^-8&)]}]];
  Module[{fc = FiberClass[G]},
    If[Length[fc] > 0, Print["Fiber class: ", fc]]];
  Print[""];
  Print[Style["Anomaly with shared external:", Bold]];
  Print[Style["  ComputeAnomaly[glued, {\"" <> sharedTag <> "\"}]", "Input"]];
  Print[InputForm[G]];
  G];

(* ── GlueIFs: quiet matrix-level glue ── *)
GlueIFs[IF1_?MatrixQ, ep1_Integer, IF2_?MatrixQ, ep2_Integer] := Module[
  {n1=Length[IF1], n2=Length[IF2], nt, G, mb},
  nt = n1+n2-1; G = ConstantArray[0, {nt,nt}]; G[[1;;n1, 1;;n1]] = IF1;
  mb[j_] := If[j==ep2, ep1+1, n1+If[j<ep2, j+1, j]];
  Do[If[j!=ep2, G[[mb[j], mb[j]]] = IF2[[j+1, j+1]]], {j, 0, n2-1}];
  Do[If[IF2[[j+1,k+1]]!=0, G[[mb[j],mb[k]]] += IF2[[j+1,k+1]]; G[[mb[k],mb[j]]] += IF2[[k+1,j+1]]], {j,0,n2-1},{k,j+1,n2-1}];
  G];
(* ── Clusters: find non-(-1) connected components ── *)
Clusters[IF_?MatrixQ] := Module[
  {n = Length[IF], nonM1, adj, visited, comps = {}, queue, comp},
  nonM1 = Select[Range[n], IF[[#, #]] != -1 &];
  adj = Table[Select[nonM1, # != i && IF[[i, #]] != 0 &], {i, nonM1}];
  adj = Association[Thread[nonM1 -> adj]];
  visited = Association[Thread[nonM1 -> False]];
  Do[
    If[!visited[s],
      comp = {}; queue = {s}; visited[s] = True;
      While[Length[queue] > 0,
        Module[{cur = First[queue]},
          queue = Rest[queue];
          AppendTo[comp, cur];
          Do[If[!visited[nb], visited[nb] = True; AppendTo[queue, nb]], {nb, adj[cur]}]]];
      comp = Sort[comp];
      AppendTo[comps, <|
        "indices" -> comp - 1,  (* 0-indexed *)
        "selfInts" -> (IF[[#, #]] & /@ comp),
        "size" -> Length[comp],
        "subIF" -> IF[[comp, comp]]|>]],
    {s, nonM1}];
  comps];

(* ── NHC Table: build reference IFs and match by eigenvalue spectrum ── *)

buildNHCEntry[name_String, selfInts_List, H_Integer, V_Integer, ints_List:{}] := Module[
  {n = Length[selfInts], IF},
  IF = DiagonalMatrix[selfInts];
  Do[Module[{k = If[i <= Length[ints], ints[[i]], 1]},
    IF[[i, i+1]] = k; IF[[i+1, i]] = k], {i, n-1}];
  <|"name" -> name, "selfInts" -> selfInts, "H" -> H, "V" -> V,
    "IF" -> IF, "eigs" -> Sort[Eigenvalues[N@IF]]|>];

nhcTable = {
  (* Single curves *)
  buildNHCEntry["-3", {-3}, 0, 8],
  buildNHCEntry["-4", {-4}, 0, 28],
  buildNHCEntry["-5", {-5}, 0, 52],
  buildNHCEntry["-6", {-6}, 0, 78],
  buildNHCEntry["-7", {-7}, 28, 133],
  buildNHCEntry["-8", {-8}, 0, 133],
  buildNHCEntry["-12", {-12}, 0, 248],
  (* Two-curve chains *)
  buildNHCEntry["-2-3", {-2, -3}, 8, 17],
  buildNHCEntry["-2-4", {-2, -4}, 256, 183, {2}],  (* int=2 *)
  (* Three-curve chains *)
  buildNHCEntry["-2-3-2", {-2, -3, -2}, 16, 27],
  buildNHCEntry["-2-2-3", {-2, -2, -3}, 8, 17]
};

eigsMatch[a_List, b_List] := Length[a] == Length[b] && Max[Abs[a - b]] < 10^-6;

nhcLookup[subIF_?MatrixQ] := Module[{eigs = Sort[Eigenvalues[N@subIF]], result = None},
  Do[If[eigsMatch[eigs, nhc["eigs"]], result = nhc; Break[]], {nhc, nhcTable}];
  If[result =!= None, Return[result]];
  If[AllTrue[Diagonal[subIF], # == -2 &],
    Return[<|"name" -> "pure-2", "H" -> 0, "V" -> 0|>]];
  <|"name" -> "???", "H" -> -1, "V" -> -1|>];

NHCLookup[subIF_?MatrixQ] := nhcLookup[subIF];

(* ── ComputeAnomaly ── *)

(* Enhancement rules for externals *)
externalHV[tag_String] := Which[
  (* su(2)_k → -1: standalone -2 gets su(2) enhancement *)
  StringMatchQ[tag, "su2" ~~ ___] && !StringContainsQ[tag, "n"],
    <|"dH" -> 4, "dV" -> 3, "desc" -> "su(2) on ext -2"|>,
  (* hat1 → -1: su(8) *)
  tag == "hat1m1",
    <|"dH" -> 36, "dV" -> 63, "desc" -> "su(8) on hat1"|>,
  (* hat1 → -2: su(16) *)
  tag == "hat1m2",
    <|"dH" -> 264, "dV" -> 255, "desc" -> "su(16) on hat1"|>,
  (* NHC-forming externals: enhancement is in the NHC cluster itself *)
  True, <|"dH" -> 0, "dV" -> 0, "desc" -> "in NHC"|>
];

(* Manual H,V *)
ComputeAnomaly[IF_?MatrixQ, Hc_Integer, V_Integer] := Module[{T, Hn},
  T = Count[Eigenvalues[N@IF], _?(# < -10^-8 &)];
  Hn = 273 + V - Hc - 29*T;
  Print["T=", T, "  H_c=", Hc, "  V=", V, "  H_n=", Hn];
  If[Hn >= 0, Print[Style["Anomaly OK", Darker[Green]]], Print[Style["H_n < 0!", Red]]];
  <|"T" -> T, "Hc" -> Hc, "V" -> V, "Hn" -> Hn|>];

(* Auto from single entry *)
ComputeAnomaly[e_Association] := Module[{},
  Print["T=", e["T"], "  H_c=", e["Hc"], "  V=", e["V"], "  H_n=", e["Hn"]];
  <|"T" -> e["T"], "Hc" -> e["Hc"], "V" -> e["V"], "Hn" -> e["Hn"]|>];

(* Auto from IF, with external tags for enhancement *)
ComputeAnomaly[IF_?MatrixQ, extTags_List] := Module[
  {cls, totalH = 0, totalV = 0, T, Hn, unknown = False},
  cls = Clusters[IF];
  Print[Style["NHC Clusters:", Bold]];
  Do[
    Module[{nhc = nhcLookup[cl["subIF"]]},
      Print["  ", cl["indices"], " ", cl["selfInts"], " \[Rule] ", nhc["name"],
        "  H=", nhc["H"], " V=", nhc["V"]];
      If[nhc["H"] == -1,
        unknown = True;
        Print["    subIF: ", MatrixForm[cl["subIF"]]],
        totalH += nhc["H"]; totalV += nhc["V"]]],
    {cl, cls}];
  (* External enhancements *)
  If[Length[extTags] > 0,
    Print[Style["External enhancements:", Bold]];
    Do[
      Module[{enh = externalHV[tag]},
        Print["  ", tag, " \[Rule] ", enh["desc"], "  dH=", enh["dH"], " dV=", enh["dV"]];
        totalH += enh["dH"]; totalV += enh["dV"]],
      {tag, extTags}]];
  If[unknown, Print[Style["WARNING: Unknown cluster! H,V incomplete.", Red]]];
  T = Count[Eigenvalues[N@IF], _?(# < -10^-8 &)];
  Hn = 273 + totalV - totalH - 29*T;
  Print[Style["Result: ", Bold], "T=", T, "  H_c=", totalH, "  V=", totalV, "  H_n=", Hn];
  If[Hn >= 0, Print[Style["Anomaly OK", Darker[Green]]], Print[Style["H_n < 0!", Red]]];
  <|"T" -> T, "Hc" -> totalH, "V" -> totalV, "Hn" -> Hn, "clusters" -> cls|>];

(* Plain IF, no external info *)
ComputeAnomaly[IF_?MatrixQ] := ComputeAnomaly[IF, {}];
FiberClass[IF_?MatrixQ] := Module[{ns, v, g},
  ns = NullSpace[IF];
  If[Length[ns] == 0, Return[{}]];
  v = ns[[1]];
  (* Make integer: multiply by LCM of denominators *)
  v = v * LCM @@ Denominator /@ v;
  (* Normalize: divide by GCD, ensure first nonzero is positive *)
  g = GCD @@ Abs /@ v;
  If[g > 0, v = v/g];
  If[First[Select[v, # != 0 &]] < 0, v = -v];
  v];

FiberClass[e_Association] := FiberClass[e["IF"]];

(* ── BlowdownEntry: blowdown with index tracking ── *)
BlowdownEntry[e_Association] := Module[{IF=e["IF"], h1, tmin},
  h1 = Cases[e["externals"], _?(#["isHat1"]&)][[All,"curveIdx"]];
  If[Abs[e["det"]]==1,
    blowdownTracked[IF, ExtIndices[e], h1, Automatic],
    tmin = TMin[e];
    If[tmin === $Failed, Return[$Failed]];
    blowdownTracked[IF, ExtIndices[e], h1, tmin]
  ]];

BlowdownEntry[e_Association, T_Integer] := Module[{h1},
  h1 = Cases[e["externals"], _?(#["isHat1"]&)][[All,"curveIdx"]];
  blowdownTracked[e["IF"], ExtIndices[e], h1, T]];

blowdownTracked[IF_?MatrixQ, extIdxs0_List, h1_List, Tval_] := Module[
  {n=Length[IF], ext, cur, labels, curH1=h1, m1, newLabels, idx, bd=0, T},
  
  (* Build extended IF *)
  T = If[Tval === Automatic,
    Count[Eigenvalues[N@IF], _?(# < -10^-8 &)],
    Tval];
  ext = ConstantArray[0, {n+1,n+1}]; ext[[1;;n,1;;n]] = IF;
  Do[Module[{b=If[MemberQ[h1,i-1],-1,IF[[i,i]]+2]},ext[[i,n+1]]=b;ext[[n+1,i]]=b],{i,n}];
  ext[[n+1,n+1]] = 9-T;
  
  (* Labels: 0-indexed original curve id, last = b0 *)
  labels = Append[Range[0, n-1], "b0"];
  cur = ext;
  
  (* Blowdown exceptional curves only: C²=-1 and b₀·C=1 *)
  While[True,
    n = Length[cur]; m1 = -1;
    Do[If[cur[[i,i]] == -1 && cur[[i, n]] == 1,
      m1 = i; Break[]], {i, n}];
    If[m1 == -1, Break[]];
    Module[{keep = Delete[Range[n], m1], newCur},
      newCur = Table[cur[[i,j]] + cur[[i,m1]]*cur[[m1,j]], {i,keep}, {j,keep}];
      labels = Delete[labels, m1];
      cur = newCur];
    bd++];
  
  (* Identify which rows are externals, b0 *)
  Module[{extRows, b0Row, extTags, result},
    extRows = Flatten[Position[labels, #] & /@ extIdxs0];
    b0Row = FirstPosition[labels, "b0"][[1]];
    extTags = Table[
      Module[{origIdx = extIdxs0[[i]], row = extRows[[i]]},
        <|"origIdx" -> origIdx,
          "row" -> row,
          "selfInt" -> cur[[row, row]],
          "b0dot" -> cur[[row, b0Row]],
          "intersections" -> cur[[row]]|>],
      {i, Length[extIdxs0]}];
    
    result = <|
      "matrix" -> cur,
      "labels" -> labels,
      "b0row" -> b0Row,
      "b0sq" -> cur[[b0Row, b0Row]],
      "externals" -> extTags,
      "sig" -> Module[{ev=Eigenvalues[N@cur]},
        {Count[ev,_?(#>10^-8&)], Count[ev,_?(#<-10^-8&)], Count[ev,_?(-10^-8<=#<=10^-8&)]}],
      "det" -> Det[cur],
      "T" -> T,
      "blowdowns" -> bd|>;
    
    (* Print summary *)
    Print[Style["Blowdown result (T=", Bold], Style[T, Bold], Style[")", Bold]];
    Print["b0\[Squared] = ", 9-T, ",  ", bd, " curves blown down,  sig = ", result["sig"], ",  det = ", result["det"]];
    Print["Remaining labels (0-idx): ", labels];
    Print[MatrixForm[cur]];
    Do[
      Module[{et = extTags[[i]]},
        Print[Style["External ", Bold], extIdxs0[[i]], 
          Style[" (C\[Squared]=", Gray], et["selfInt"], Style[")", Gray],
          ":  b0\[CenterDot]C = ", et["b0dot"],
          ",  row = ", cur[[et["row"]]]]],
      {i, Length[extTags]}];
    
    result
  ]];
ShowEntry[e_Association] := PrettyEntry[e];

(* IF만 보기: 행렬 + copy-pasteable 형태 *)
ShowIF[e_Association] := Module[{},
  Print[Style["IF (" <> ToString[Length[e["IF"]]] <> "\[Times]" <> ToString[Length[e["IF"]]] <> "):", Bold]];
  Print[MatrixForm[e["IF"]]];
  Print[""];
  Print[Style["Copy-pasteable:", Gray]];
  Print[InputForm[e["IF"]]];
];

(* ── Filter ── *)
FindEntries[entries_List,spec_String] := Select[entries,#["combo"]==spec&];
FindEntries[entries_List,spec_String,T_Integer] := Select[entries,#["combo"]==spec&&#["T"]==T&];

DetSquareQ[e_Association] := IntegerQ[Sqrt[Abs[e["det"]]]];
SelectDetSquare[entries_List] := Select[entries, DetSquareQ];

(* ── Summary ── *)
CatalogSummary[entries_List] := Module[{combos,Tdist,nUni,comboData},
  nUni=Count[entries,_?(Abs[#["det"]]==1&)];
  Tdist=SortBy[Tally[entries[[All,"T"]]],First];
  comboData = Association @@ (# -> <|"uni"->0,"nonuni"->0|> & /@ Union[entries[[All,"combo"]]]);
  Do[If[Abs[e["det"]]==1, comboData[e["combo"]]["uni"]++, comboData[e["combo"]]["nonuni"]++], {e,entries}];
  combos = SortBy[{#, comboData[#]["uni"], comboData[#]["nonuni"], comboData[#]["uni"]+comboData[#]["nonuni"]} & /@ Keys[comboData], -Last[#]&];
  Print[Style["Catalog: ",Bold],Length[entries]," entries (",nUni," uni, ",Length[entries]-nUni," non-uni)"];
  Print[Grid[Prepend[Tdist,{"T","Count"}],Frame->All]];
  Print[Grid[Prepend[combos,{"Combo","Uni","Non-uni","Total"}],Frame->All]]];

(* ── Chain Notation (cat2tex style) ── *)

ChainNotation[e_Association] := Module[
  {IF = e["IF"], n, extIdxs, adj, visited, start, chain, extSet, result},
  n = Length[IF];
  extIdxs = ExtIndices[e] + 1; (* 1-indexed *)
  extSet = Association[Thread[extIdxs -> Range[Length[extIdxs]]]];

  (* Build adjacency excluding externals *)
  adj = Table[
    Select[Range[n], # != i && IF[[i, #]] != 0 && !KeyExistsQ[extSet, #] &],
    {i, n}];

  (* Find start: leaf that is NOT a target of external *)
  Module[{leaves, targets},
    targets = {};
    Do[Do[If[IF[[ext, j]] != 0 && !KeyExistsQ[extSet, j], AppendTo[targets, j]], {j, n}], {ext, extIdxs}];
    leaves = Select[Range[n], !KeyExistsQ[extSet, #] && Length[adj[[#]]] <= 1 &];
    leaves = Select[leaves, !MemberQ[targets, #] &];
    start = If[Length[leaves] > 0, First[leaves],
      First[Select[Range[n], !KeyExistsQ[extSet, #] &], 1]]
  ];

  (* Walk chain *)
  visited = ConstantArray[False, n];
  chain = {};
  Module[{walk},
    walk[cur_] := Module[{si, extMark = "", children},
      visited[[cur]] = True;
      si = -IF[[cur, cur]];
      (* Check if external connects here *)
      Do[If[IF[[ext, cur]] != 0,
        extMark = extMark <> Superscript[Style[ToString[Abs[IF[[ext, cur]]]], Red], ""]],
        {ext, extIdxs}];
      AppendTo[chain, Row[{si, extMark}]];
      children = Select[adj[[cur]], !visited[[#]] &];
      Do[walk[c], {c, children}]];
    walk[start]];

  Row[chain, ""]
];

(* ── Gauge Formatting ── *)

FormatGaugeName[s_String] := Module[{parts, formatted},
  If[s == "--" || s == "", return[Style["--", Gray]]];
  parts = StringSplit[s, "+"];
  formatted = Table[
    Switch[p,
      "su2", Subscript[Style["\[GothicS]\[GothicU]", Italic], 2],
      "su3", Subscript[Style["\[GothicS]\[GothicU]", Italic], 3],
      "su8", Subscript[Style["\[GothicS]\[GothicU]", Italic], 8],
      "sp1", Subscript[Style["\[GothicS]\[GothicP]", Italic], 1],
      "g2", Subscript[Style["\[GothicG]", Italic], 2],
      "so7", Subscript[Style["\[GothicS]\[GothicO]", Italic], 7],
      "so8", Subscript[Style["\[GothicS]\[GothicO]", Italic], 8],
      "so16", Subscript[Style["\[GothicS]\[GothicO]", Italic], 16],
      "f4", Subscript[Style["\[GothicF]", Italic], 4],
      "e6", Subscript[Style["\[GothicE]", Italic], 6],
      "e7", Subscript[Style["\[GothicE]", Italic], 7],
      "e7'", Subsuperscript[Style["\[GothicE]", Italic], 7, "\[Prime]"],
      "e8", Subscript[Style["\[GothicE]", Italic], 8],
      "hat1", OverHat["1"],
      _, p],
    {p, parts}];
  Row[Riffle[formatted, "\[DirectSum]"]]
];

(* ── Pretty Display ── *)

PrettyEntry[e_Association] := Module[
  {sig, delta, Tmax, tmin, extGauge, lstGauge, chainStr},
  sig = {e["sigPos"], e["sigNeg"], e["sigZero"]};
  delta = e["Hc"] - e["V"] + 29*e["T"] - 273;
  Tmax = TMax[e];

  (* External gauge *)
  extGauge = StringRiffle[Sort[Table[
    Switch[ext["extSI"], -2,"su2", -3,"su3", -4,"so8", -5,"f4",
           -6,"e6", -7,"e7'", -8,"e7", -12,"e8", -1,"hat1", _, ToString[ext["extSI"]]],
    {ext, e["externals"]}]], "+"];

  (* LST gauge using NHCLookup *)
  lstGauge = LSTGauge[e];

  Print[Style["━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", Gray]];
  Print[Style["Chain: ", Bold], ChainNotation[e]];
  Print[Style["T = ", Bold], e["T"], "   ",
        Style["External: ", Bold], FormatGaugeName[extGauge], "   ",
        Style["LST Gauge: ", Bold], FormatGaugeName[lstGauge]];
  Print["det = ", e["det"], "  sig = ", sig, "  \[CapitalDelta] = ", delta,
        "  T_max = ", Tmax];
  Print["H_c = ", e["Hc"], "  V = ", e["V"], "  H_n = ", e["Hn"]];
  If[Abs[e["det"]] > 1,
    tmin = TMin[e];
    If[tmin =!= $Failed, Print["T_min = ", tmin, "  T \[Element] [", tmin, ", ", Tmax, "]"]]];
  Print[EntryQuiver[e]];
];

(* ── Batch Pretty Display ── *)
PrettyList[entries_List, max_:20] := Do[
  PrettyEntry[entries[[i]]];
  Print[""],
  {i, Min[Length[entries], max]}];

End[];
EndPackage[];

Print["SUGRACatalog loaded."];
Print["  entries = LoadCatalog[\"cat_1ext\"]"];
Print["  CatalogSummary[entries]"];
Print["  ShowEntry[entries[[1]]]  (* copy-pasteable IF *)"];
Print["  FindEntries[entries, \"e8\"]"];
Print["  TMin[entry], TMax[entry]"];
