
$RepoRoot = DirectoryName[$InputFileName];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

$testModel = <|
    "Vertices List"                    -> {1, 2},
    "Adjacency Matrix"                 -> {{0, 1}, {0, 0}}, (* 1 -> 2 only *)
    "Entrance Vertices and Flows"      -> {{1, 10}},
    "Exit Vertices and Terminal Costs" -> {{2, 0}},
    "Switching Costs"                  -> {}
|>;

Print["Building scenario..."];
s = makeScenario[<|"Model" -> $testModel|>];
If[FailureQ[s], Print["Failed to make scenario: ", s]; Exit[1]];

Print["Building system..."];
sys = makeSystem[s];
If[FailureQ[sys], Print["Failed to make system: ", sys]; Exit[1]];

data = SystemDataFlatten[sys];

Print["\n--- SignedFlows ---"];
Print[Keys[data["SignedFlows"]]];

Print["\n--- BalanceSplittingFlows ---"];
Print[data["BalanceSplittingFlows"]];

Exit[0];
