payloadRowsToStateSummaryRows[payload_Association] := 
    Module[{schema, rows, onMinus1Dict, onIslandDict}, 
     schema = Lookup[payload, "schema", {}]; 
      rows = Lookup[payload, "rows", {}]; onMinus1Dict = 
       Lookup[payload, "onMinus1Dict", {}]; onIslandDict = 
       Lookup[payload, "onIslandDict", {}]; If[schema === {} || rows === {}, 
       Return[{}]]; Function[row, Module[{a, i1, i2}, 
         a = AssociationThread[schema, row]; i1 = Lookup[a, "OnMinus1Index", 
            Missing["NoIndex"]]; i2 = Lookup[a, "OnIslandIndex", 
            Missing["NoIndex"]]; Join[KeyDrop[a, {"OnMinus1Index", 
             "OnIslandIndex"}], Association["OnMinus1" -> 
             If[IntegerQ[i1] && 1 <= i1 <= Length[onMinus1Dict], 
              onMinus1Dict[[i1]], ""], "OnIsland" -> If[IntegerQ[i2] && 1 <= 
                i2 <= Length[onIslandDict], onIslandDict[[i2]], ""]]]]] /@ 
       rows]
 
loadCompactStateSummary[file_] := Import[file, "WXF"]
 
LoadCompactStateSummaryDirectory[dir_, Tsugra_:All] := 
    Module[{pattern, files}, pattern = If[Tsugra === All, "*.wxf", 
        StringJoin["T", ToString[Tsugra], "_*.wxf"]]; 
      files = FileNames[pattern, dir]; AssociationThread[
       FileBaseName /@ files, loadCompactStateSummary /@ files]]
 
StateSummaryFromCompactDirectory[dir_, Tsugra_, which_:All] := 
    Module[{allData, selected, rows, keys}, 
     allData = LoadCompactStateSummaryDirectory[dir, Tsugra]; 
      keys = Switch[which, All, Keys[allData], _Integer, 
        {StringJoin["T", ToString[Tsugra], "_i", ToString[which]]}, 
        {_Integer..}, (StringJoin["T", ToString[Tsugra], "_i", 
           ToString[#1]] & ) /@ which, _String, {which}, _, Keys[allData]]; 
      selected = KeyTake[allData, Intersection[keys, Keys[allData]]]; 
      rows = Flatten[KeyValueMap[Function[{fname, payload}, 
          payloadRowsToStateSummaryRows[payload]], selected], 1]; 
      Dataset[SortBy[rows, {Replace[Lookup[#1, "Tmin", Missing["NoTmin"]], 
           _Missing -> Infinity] & , Replace[Lookup[#1, "Delta", 
            Missing["NoDelta"]], _Missing -> Infinity] & , 
         Lookup[#1, "n2", 0] & , Lookup[#1, "n3", 0] & , 
         Lookup[#1, "n4", 0] & , Lookup[#1, "n5", 0] & , 
         Lookup[#1, "n6", 0] & , Lookup[#1, "n7", 0] & , 
         Lookup[#1, "n8", 0] & , Lookup[#1, "n12", 0] & , 
         Lookup[#1, "OnMinus1", ""] & , Lookup[#1, "OnIsland", ""] & }]]]
 
groupStateSummaryByTotalAdded[ds_Dataset] := 
    Module[{rows, totalAdded, grouped}, rows = Normal[ds]; 
      totalAdded[row_] := Total[Lookup[row, {"n2", "n3", "n4", "n5", "n6", 
          "n7", "n8", "n12"}, 0]]; grouped = GroupBy[rows, totalAdded]; 
      KeySort[Association[KeyValueMap[Function[{k, vals}, 
          k -> Dataset[SortBy[vals, {Replace[Lookup[#1, "Tmin", Missing[
                  "NoTmin"]], _Missing -> Infinity] & , Replace[Lookup[#1, 
                 "Delta", Missing["NoDelta"]], _Missing -> Infinity] & , 
              Lookup[#1, "n2", 0] & , Lookup[#1, "n3", 0] & , 
              Lookup[#1, "n4", 0] & , Lookup[#1, "n5", 0] & , 
              Lookup[#1, "n6", 0] & , Lookup[#1, "n7", 0] & , 
              Lookup[#1, "n8", 0] & , Lookup[#1, "n12", 0] & , 
              Lookup[#1, "OnMinus1", ""] & , Lookup[#1, "OnIsland", 
                ""] & }]]], grouped]]]]
