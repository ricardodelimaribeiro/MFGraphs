Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

Print["MFGraphs loaded"];

Print["Testing makeScenario..."];
raw = <|
    "Model" -> <|
        "Vertices List" -> {1, 2, 3},
        "Adjacency Matrix" -> SparseArray[{{1, 2} -> 1, {2, 3} -> 1}, {3, 3}],
        "Entrance Vertices and Flows" -> {{1, 100}},
        "Exit Vertices and Terminal Costs" -> {{2, 0}, {3, 20}},
        "Switching Costs" -> <|{1, 2, 3} -> 10|>
    |>
|>;

sc = makeScenario[raw];
Print["Scenario head: ", Head[sc]];
If[FailureQ[sc], Print["Failure: ", sc]];

Print["Testing makeSystem..."];
sys = makeSystem[sc];
Print["System head: ", Head[sys]];

Print["Testing SystemDataFlatten..."];
d2e = SystemDataFlatten[sys];
Print["D2E keys: ", Keys[d2e]];
Print["CONSISTENT_COSTS=", d2e["consistentCosts"]];

Print["Testing solver..."];
result = ReduceSystem[sys];
Print["Result keys: ", If[AssociationQ[result], Keys[result], Head[result]]];
If[AssociationQ[result],
  Print["Rules: ", result["Rules"]];
  Print["Equations: ", result["Equations"]];
];
