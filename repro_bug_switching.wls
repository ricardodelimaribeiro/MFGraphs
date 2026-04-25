
$RepoRoot = DirectoryName[$InputFileName];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

$testModel = <|
    "Vertices List"                    -> {1, 2, 3, 4},
    "Adjacency Matrix"                 -> {{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
    "Entrance Vertices and Flows"      -> {{1, 10}},
    "Exit Vertices and Terminal Costs" -> {{3, 0}, {4, 0}},
    "Switching Costs"                  -> <|{1, 2, 3} -> 5, {1, 2, 4} -> 10|>
|>;

Print["Building scenario..."];
s = makeScenario[<|"Model" -> $testModel|>];
If[FailureQ[s], Print["Failed to make scenario: ", s]; Exit[1]];

(* Manually fix the functions to test the fix *)
Print["Redefining IneqSwitch..."];
Unprotect[IneqSwitch];
IneqSwitch[u_, Switching_Association][r_, i_, w_] :=
    u[i, r] <= u[i, w] + Switching[{r, i, w}];

Unprotect[AltSwitch];
AltSwitch[j_, u_, Switching_][r_, i_, w_] :=
    (j[r, i, w] == 0) || (u[i, r] == u[i, w] + Switching[{r, i, w}]);

Print["Information for IneqSwitch:"];
Print[Information[IneqSwitch]];

Print["Direct test of IneqSwitch:"];
Print[IneqSwitch[u, <|{1, 2, 3} -> 5|>][1, 2, 3]];

Print["Building system (with manual fix)..."];
sys = makeSystem[s];
If[FailureQ[sys], Print["Failed to make system: ", sys]; Exit[1]];

data = SystemDataFlatten[sys];

Print["\n--- Keys in data ---"];
Print[Keys[data]];

Print["\n--- Switching Costs ---"];
Print[data["SwitchingCosts"]];

Print["\n--- consistentCosts ---"];
Print[data["consistentCosts"]];

Print["\n--- IneqSwitchingByVertex ---"];
Print[data["IneqSwitchingByVertex"]];

Print["\n--- AltOptCond ---"];
Print[data["AltOptCond"]];

Print["\n--- activeAuxTriples ---"];
Print[data["activeAuxTriples"]];

Exit[0];
