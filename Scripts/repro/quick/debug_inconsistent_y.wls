Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

Print["Testing Inconsistent Y shortcut..."];
(* Use the 3-arg version which resolves Automatic *)
entries = {{1, 100}};
exits = {{3, 0}, {4, 0}};

s = getExampleScenario["Inconsistent Y shortcut", entries, exits, Automatic];
Print["Scenario valid: ", ScenarioQ[s]];
If[FailureQ[s], Print["Failure: ", s]];

sys = makeSystem[s];
d2e = systemDataFlatten[sys];
Print["CONSISTENT_COSTS=", d2e["ConsistentCosts"]];

Print["Testing solver on inconsistent case..."];
result = reduceSystem[sys];
Print["Result keys: ", If[AssociationQ[result], Keys[result], Head[result]]];
If[AssociationQ[result],
  Print["Feasibility: ", If[ListQ[result["Rules"]] && result["Rules"] =!= {}, "Feasible", "Infeasible"]]
];
