(* catalog_analysis.m — SUGRA Catalog Analysis Notebook *)
(* 카탈로그에서 LST gauge, external gauge, T로 검색 *)

SetDirectory["/Users/seongminjeon/Desktop/sugra_pipeline"]
<< "SUGRACatalog.m"

(* ═══════════════════════════════════════════════════ *)
(* 1. 카탈로그 로드 *)
(* ═══════════════════════════════════════════════════ *)

(* 원하는 T 범위 선택 *)
entries = LoadCatalog["cat_1ext_t010"];
(* entries = LoadCatalog["cat_2ext_t010"]; *)
(* entries = LoadCatalog["cat_1ext_t1120"]; *)

(* 여러 디렉토리 합치기 *)
(*
allEntries = Join[
    LoadCatalog["cat_1ext_t010"],
    LoadCatalog["cat_1ext_t1120"],
    LoadCatalog["cat_1ext_t2130"]
];
*)

(* ═══════════════════════════════════════════════════ *)
(* 2. NHC Gauge Algebra 계산 *)
(* ═══════════════════════════════════════════════════ *)

(* IF에서 external 제외한 base의 gauge algebra 계산 (NHCLookup 사용) *)
LSTGauge[entry_Association] := Module[
  {IF = entry["IF"], n, extIdxs, nonM1, visited, comps, gauges = {}},
  n = Length[IF];
  extIdxs = ExtIndices[entry] + 1; (* 1-indexed *)

  (* non-(-1), non-external curves *)
  nonM1 = Select[Range[n], IF[[#, #]] != -1 && !MemberQ[extIdxs, #] &];

  (* BFS connected components *)
  visited = Association[Thread[nonM1 -> False]];
  Do[
    If[!visited[s],
      Module[{comp = {}, queue = {s}},
        visited[s] = True;
        While[Length[queue] > 0,
          Module[{cur = First[queue]},
            queue = Rest[queue];
            AppendTo[comp, cur];
            Do[If[!visited[nb] && IF[[cur, nb]] != 0,
              visited[nb] = True; AppendTo[queue, nb]], {nb, nonM1}]]];
        comp = Sort[comp];

        (* Extract submatrix and use NHCLookup for eigenvalue-based matching *)
        Module[{subIF, nhcResult},
          subIF = IF[[comp, comp]];
          (* Skip pure -2 chains *)
          If[AllTrue[Diagonal[subIF], # == -2 &], Null,
            nhcResult = NHCLookup[subIF];
            AppendTo[gauges, nhcResult["name"]]
          ]]]],
    {s, nonM1}];

  StringRiffle[Sort[DeleteCases[gauges, Null]], " + "]
];

(* External gauge 이름 *)
ExtGauge[entry_Association] := Module[
  {exts = entry["externals"], gauges},
  gauges = Table[
    Switch[ext["extSI"],
      -2, "su2", -3, "su3", -4, "so8", -5, "f4",
      -6, "e6", -7, "e7'", -8, "e7", -12, "e8", -1, "hat1",
      _, ToString[ext["extSI"]]],
    {ext, exts}];
  StringRiffle[Sort[gauges], "+"]
];

(* ═══════════════════════════════════════════════════ *)
(* 3. 검색 함수들 *)
(* ═══════════════════════════════════════════════════ *)

(* T값으로 검색 *)
FindByT[entries_List, T_Integer] := Select[entries, #["T"] == T &];

(* External gauge로 검색 *)
FindByExtGauge[entries_List, gauge_String] :=
  Select[entries, ExtGauge[#] == gauge &];

(* LST gauge로 검색 *)
FindByLSTGauge[entries_List, gauge_String] :=
  Select[entries, LSTGauge[#] == gauge &];

(* 조합 검색: T + External + LST gauge *)
FindBases[entries_List, T_Integer, extGauge_String, lstGauge_String] :=
  Select[entries,
    #["T"] == T && ExtGauge[#] == extGauge && LSTGauge[#] == lstGauge &];

(* 부분 매칭: LST gauge에 특정 algebra가 포함된 것 *)
FindByLSTGaugeContains[entries_List, substr_String] :=
  Select[entries, StringContainsQ[LSTGauge[#], substr] &];

(* ═══════════════════════════════════════════════════ *)
(* 4. 통계 함수들 *)
(* ═══════════════════════════════════════════════════ *)

(* 전체 분류 테이블 생성 *)
ClassificationTable[entries_List] := Module[
  {groups, table},
  groups = GatherBy[entries, {#["T"], ExtGauge[#], LSTGauge[#]} &];
  table = Table[
    Module[{e = First[g]},
      <|"T" -> e["T"], "ExtGauge" -> ExtGauge[e],
        "LSTGauge" -> LSTGauge[e], "Count" -> Length[g]|>],
    {g, groups}];
  SortBy[table, {#["T"], -#["Count"]} &]
];

(* T별 요약 *)
SummaryByT[entries_List] := Module[{groups},
  groups = GroupBy[entries, #["T"] &];
  Table[{T, Length[groups[T]]}, {T, Sort[Keys[groups]]}] //
    TableForm[#, TableHeadings -> {None, {"T", "Count"}}] &
];

(* External gauge별 요약 *)
SummaryByExtGauge[entries_List] := Module[{groups},
  groups = GroupBy[entries, ExtGauge];
  Table[{g, Length[groups[g]]}, {g, Sort[Keys[groups]]}] //
    TableForm[#, TableHeadings -> {None, {"ExtGauge", "Count"}}] &
];

(* LST gauge별 요약 *)
SummaryByLSTGauge[entries_List] := Module[{groups},
  groups = GroupBy[entries, LSTGauge];
  SortBy[
    Table[{g, Length[groups[g]]}, {g, Keys[groups]}],
    -#[[2]] &] //
    TableForm[#, TableHeadings -> {None, {"LSTGauge", "Count"}}] &
];

(* ═══════════════════════════════════════════════════ *)
(* 5. 사용 예시 *)
(* ═══════════════════════════════════════════════════ *)

(* --- 기본 통계 --- *)
Print["Total entries: ", Length[entries]];
SummaryByT[entries]
SummaryByExtGauge[entries]
SummaryByLSTGauge[entries]

(* --- 특정 조건 검색 --- *)

(* T=10, su2 external, e6+f4+su3 gauge를 가진 base들 *)
hits = FindBases[entries, 10, "su2", "e6+f4+su3"];
Length[hits]
ShowEntry[hits[[1]]]

(* T=7, e6 LST gauge를 가진 모든 base *)
e6bases = Select[entries, #["T"] == 7 && LSTGauge[#] == "e6" &];
Length[e6bases]

(* LST gauge에 "e8"이 포함된 모든 base *)
e8bases = FindByLSTGaugeContains[entries, "e8"];
Length[e8bases]

(* --- 특정 base의 상세 정보 --- *)
(*
e = hits[[1]];
Print["T = ", e["T"]];
Print["External: ", ExtGauge[e]];
Print["LST Gauge: ", LSTGauge[e]];
Print["det = ", e["det"], "  sig = ", {e["sigPos"], e["sigNeg"], e["sigZero"]}];
Print["H = ", e["Hc"], "  V = ", e["V"], "  Hn = ", e["Hn"]];
ShowEntry[e]
BlowdownEntry[e]
*)

(* --- 분류 테이블 (전체) --- *)
(* ct = ClassificationTable[entries]; *)
(* TableForm[ct[[1;;20]], TableHeadings -> {None, {"T","Ext","LST Gauge","Count"}}] *)
