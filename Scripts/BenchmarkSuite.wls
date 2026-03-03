#!/usr/bin/env wolframscript
(* BenchmarkSuite.wls -- Main benchmark runner for MFGraphs *)
(* Usage: wolframscript -file Scripts/BenchmarkSuite.wls [tier] *)
(* tier: "small", "medium", "large", "vlarge", or "all" (default: "all") *)

Print["MFGraphs Benchmark Suite"];
Print["========================"];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print["OS: ", $SystemID];
Print[];

(* --- Setup package path --- *)
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

(* --- Load helpers --- *)
Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

(* --- Parse command-line arguments --- *)
$SelectedTier = If[Length[$ScriptCommandLine] >= 2,
  $ScriptCommandLine[[2]],
  "all"
];

Print["Selected tier: ", $SelectedTier];
Print[];

(* --- Determine which cases to run --- *)
$CasesToRun = If[$SelectedTier === "all",
  Flatten[Values[$BenchmarkTiers]],
  Lookup[$BenchmarkTiers, $SelectedTier, {}]
];

GetTier[key_] := SelectFirst[Keys[$BenchmarkTiers], MemberQ[$BenchmarkTiers[#], key] &, "unknown"];

(* --- Benchmark a single case with a single solver --- *)

BenchmarkOneSolver[key_, tier_String, solverName_String, d2e_Association, timeout_Integer] :=
  Module[{result, meta, residual = Missing["N/A"]},
    meta = GraphMetadata[d2e];
    Print["  Solver: ", solverName, " (timeout: ", timeout, "s)"];

    result = Switch[solverName,
      "CriticalCongestion",
        SafeExecute[CriticalCongestionSolver[d2e], timeout],
      "NonLinear",
        Module[{critResult},
          (* NonLinear needs critical solution first *)
          critResult = SafeExecute[CriticalCongestionSolver[d2e], timeout];
          If[critResult["Status"] === "OK",
            Block[{V = Function[{x, edge}, 0]},
              SafeExecute[NonLinear[critResult["Result"]], timeout]
            ],
            <|"Time" -> 0, "Memory" -> 0, "Status" -> "SKIPPED",
              "Result" -> "Critical solver failed"|>
          ]
        ],
      "Monotone",
        SafeExecute[MonotoneSolver[d2e], timeout],
      _,
        <|"Time" -> 0, "Memory" -> 0, "Status" -> "SKIPPED", "Result" -> "Unknown solver"|>
    ];

    (* Compute residual for successful solves *)
    If[result["Status"] === "OK" && AssociationQ[result["Result"]],
      residual = Switch[solverName,
        "CriticalCongestion", ComputeResidual[result["Result"], "AssoCritical"],
        "NonLinear", ComputeResidual[result["Result"], "AssoNonCritical"],
        "Monotone", Missing["N/A"],
        _, Missing["N/A"]
      ]
    ];

    Print["    Status: ", result["Status"],
          " | Time: ", FormatTime[result["Time"]], "s",
          " | Memory: ", FormatMemoryMB[result["Memory"]], " MB"];

    <|
      "Key" -> ToString[key],
      "Tier" -> tier,
      "Solver" -> solverName,
      "NumVertices" -> meta["NumVertices"],
      "NumEdges" -> meta["NumEdges"],
      "NumEntries" -> meta["NumEntries"],
      "NumExits" -> meta["NumExits"],
      "HasSwitching" -> meta["HasSwitching"],
      "D2ETime" -> Missing["Included"],
      "SolveTime" -> result["Time"],
      "SolveMemory" -> result["Memory"],
      "Status" -> result["Status"],
      "Residual" -> residual,
      "Timestamp" -> DateString[]
    |>
  ];

(* --- Benchmark a single case across all applicable solvers --- *)

BenchmarkCase[key_] :=
  Module[{tier, timeout, data, d2eResult, d2e, results = {}, params},
    tier = GetTier[key];
    timeout = Lookup[$TierTimeouts, tier, 300];
    params = GetParams[key];

    Print["\n--- Case: ", key, " (tier: ", tier, ") ---"];

    (* Load data *)
    data = Quiet @ Check[DataG[key], $Failed];
    If[data === $Failed,
      Print["  FAILED to load DataG[", key, "]"];
      Return[{<|"Key" -> ToString[key], "Tier" -> tier, "Solver" -> "ALL",
               "Status" -> "FAILED", "Timestamp" -> DateString[]|>}]
    ];

    (* Apply parameter substitutions *)
    If[params =!= {},
      data = data /. params
    ];

    (* Convert to equations (timed separately) *)
    Print["  Running Data2Equations..."];
    d2eResult = SafeExecute[Data2Equations[data], timeout];
    Print["    D2E Status: ", d2eResult["Status"],
          " | Time: ", FormatTime[d2eResult["Time"]], "s"];

    If[d2eResult["Status"] =!= "OK",
      Return[{<|"Key" -> ToString[key], "Tier" -> tier, "Solver" -> "Data2Equations",
               "D2ETime" -> d2eResult["Time"], "Status" -> d2eResult["Status"],
               "Timestamp" -> DateString[]|>}]
    ];

    d2e = d2eResult["Result"];

    (* Run each solver *)
    AppendTo[results,
      BenchmarkOneSolver[key, tier, "CriticalCongestion", d2e, timeout]];

    AppendTo[results,
      BenchmarkOneSolver[key, tier, "NonLinear", d2e, timeout]];

    AppendTo[results,
      BenchmarkOneSolver[key, tier, "Monotone", d2e, timeout]];

    (* Store D2E time in all results *)
    results = Map[Append[#, "D2ETime" -> d2eResult["Time"]] &, results];

    results
  ];

(* --- Main execution --- *)

Print["Running benchmarks for ", Length[$CasesToRun], " cases...\n"];

$AllResults = Flatten[BenchmarkCase /@ $CasesToRun];

(* --- Export results --- *)
$ResultsDir = FileNameJoin[{Directory[], "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];

$Timestamp = DateString[{"Year", "Month", "Day", "-", "Hour", "Minute"}];
$CSVPath = FileNameJoin[{$ResultsDir, "benchmark_" <> $Timestamp <> ".csv"}];
$JSONPath = FileNameJoin[{$ResultsDir, "benchmark_" <> $Timestamp <> ".json"}];
$LatestCSV = FileNameJoin[{$ResultsDir, "benchmark_latest.csv"}];
$LatestJSON = FileNameJoin[{$ResultsDir, "benchmark_latest.json"}];

ExportResultsCSV[$AllResults, $CSVPath];
ExportResultsJSON[$AllResults, $JSONPath];
ExportResultsCSV[$AllResults, $LatestCSV];
ExportResultsJSON[$AllResults, $LatestJSON];

(* --- Print summary --- *)
PrintResultsSummary[$AllResults];

Print["\nBenchmark suite completed at ", DateString[]];
