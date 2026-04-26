Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

Print["MFGraphs loaded"];

Print["Testing makeScenario..."];
raw = <|
    "Model" -> <|
        "Vertices" -> {1, 2, 3},
        "Adjacency" -> SparseArray[{{1, 2} -> 1, {2, 3} -> 1}, {3, 3}],
        "Entries" -> {{1, 100}},
        "Exits" -> {{2, 0}, {3, 20}},
        "Switching" -> <|{1, 2, 3} -> 10|>
    |>
|>;

sc = makeScenario[raw];
Print["Scenario head: ", Head[sc]];
If[FailureQ[sc], Print["Failure: ", sc]];

Print["Testing makeSystem..."];
sys = makeSystem[sc];
Print["System head: ", Head[sys]];

Print["Testing systemDataFlatten..."];
d2e = systemDataFlatten[sys];
Print["D2E keys: ", Keys[d2e]];
Print["CONSISTENT_COSTS=", d2e["ConsistentCosts"]];

Print["Testing solver..."];
result = reduceSystem[sys];
Print["Result keys: ", If[AssociationQ[result], Keys[result], Head[result]]];
If[AssociationQ[result],
  Print["Rules: ", result["Rules"]];
  Print["Equations: ", result["Equations"]];
];
