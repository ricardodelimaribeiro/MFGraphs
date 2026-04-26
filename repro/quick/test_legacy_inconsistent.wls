
(* Load archives to support legacy tests *)
mfgDir = "/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs";
Get[FileNameJoin[{mfgDir, "MFGraphs.wl"}]];
Get[FileNameJoin[{mfgDir, "archive", "DataToEquations.wl"}]];
Get[FileNameJoin[{mfgDir, "archive", "Solvers.wl"}]];
Get[FileNameJoin[{mfgDir, "Examples", "archive", "ExamplesData.wl"}]];

Print["Loaded legacy environment"];

data = Append[
    MFGraphs`GetExampleData["chain with two exits"] /. {
        MFGraphs`I1 -> 100,
        MFGraphs`U1 -> 0,
        MFGraphs`U2 -> 0
    },
    "Switching" -> {{1, 2, 3, 10}}
];

Print["DATA_SC=" , data["Switching"]];

d2e = Quiet[
    MFGraphs`DataToEquations[data],
    {MFGraphs`DataToEquations::switchingcosts}
];

Print["D2E_SC_KEYS=", Keys[d2e["SwitchingCosts"]]];
Print["D2E_SC_VALUES=", Values[d2e["SwitchingCosts"]]];

consistency = MFGraphs`IsSwitchingCostConsistent[Normal @ Lookup[d2e, "SwitchingCosts", <||>]];
Print["CONSISTENCY=", consistency];

result = Quiet[
    MFGraphs`CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 60],
    {MFGraphs`DataToEquations::switchingcosts}
];

Print["SOLVER=", result["Solver"]];
Print["FEASIBILITY=", result["Feasibility"]];
