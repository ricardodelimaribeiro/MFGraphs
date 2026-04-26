Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

scList = {{1,2,3,5},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}};

Print["ListQ: ", ListQ[scList]];
Print["All Lists: ", AllTrue[scList, ListQ]];
Print["All Length 4: ", AllTrue[scList, Length[#] === 4 &]];
Print["All Take 3 IntegerQ: ", AllTrue[scList, AllTrue[Take[#, 3], IntegerQ] &]];
Print["All Last NumericQ: ", AllTrue[Last /@ scList, NumericQ]];

(* Check what makeScenario sees *)
model = <|
    "Vertices" -> {1, 2, 3, 4},
    "Adjacency" -> {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}},
    "Entries" -> {{1, 100}},
    "Exits" -> {{3, 0}, {4, 0}},
    "Switching" -> scList
|>;

Print["SwitchingCostsNumericQ result: ", MFGraphs`Private`SwitchingCostsNumericQ[model]];
