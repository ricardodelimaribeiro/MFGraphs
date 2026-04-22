#!/usr/bin/env wolframscript
(* Extract and validate paper tier benchmark results *)

(* Load latest benchmark CSV *)
$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
csvFile = FileNameJoin[{$RepoRoot, "Results", "benchmark_latest.csv"}];
If[!FileExistsQ[csvFile],
    Print["ERROR: ", csvFile, " not found"];
    Exit[1]
];

(* Import CSV as Association *)
data = Import[csvFile, {"CSV", "Dataset"}];

(* Filter for paper tier *)
paperTierData = data[Select[#"Tier" == "paper" &]];

If[Length[paperTierData] == 0,
    Print["WARNING: No paper tier data found yet"];
    Exit[0]
];

Print["=== Paper Tier (HRF Scenario 1) Validation Results ==="];
Print[""];

(* Extract by solver *)
Do[
    solverData = paperTierData[Select[#"Solver" == solver &]];
    If[Length[solverData] > 0,
        row = solverData[1];
        Print["Solver: ", solver];
        Print["  Status: ", row["Status"]];
        Print["  ResultKind: ", row["ResultKind"]];
        Print["  Feasibility: ", row["Feasibility"]];
        Print["  WallTime: ", row["WallTime"], "s"];
        Print["  Memory: ", row["ChildMaxMemoryUsed"] / 1024^2, " MB"];
        If[row["Status"] === "OK",
            Print["  ✅ PASS"],
            Print["  ❌ FAIL: ", row["Message"]]
        ];
        Print[""]
    ],
    {solver, {"CriticalCongestion"}}
]

(* Validate critical solver produced feasible solution *)
criticalData = paperTierData[Select[#"Solver" == "CriticalCongestion" &]];
If[Length[criticalData] > 0,
    row = criticalData[1];
    Print["=== Validation Summary ==="];
    Print[""];
    If[row["Status"] === "OK" && row["ResultKind"] === "Success" && row["Feasibility"] === "Feasible",
        Print["✅ HRF Scenario 1 VALIDATED"];
        Print["   CriticalCongestion found feasible equilibrium"];
        Print["   Wall-clock time: ", row["WallTime"], "s"];
        Print["   Safe for paper publication"],
        Print["⚠️  HRF Scenario 1 validation INCOMPLETE"];
        Print["   Status: ", row["Status"], " | ResultKind: ", row["ResultKind"], " | Feasibility: ", row["Feasibility"]]
    ]
]

Print[""];
Print["Results file: ", csvFile];
Print["Results timestamp: ", data[1]["Timestamp"]];
