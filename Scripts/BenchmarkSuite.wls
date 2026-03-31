(* Usage: wolframscript -file Scripts/BenchmarkSuite.wls [tier] *)
(* tier: "small", "medium", "large", "vlarge", or "all" (default: "all") *)

(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol::Print *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
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
(* Support: wolframscript -file Scripts/BenchmarkSuite.wls [tier] [--tag "description"] *)
Module[{args = Rest[$ScriptCommandLine], tagPos},
  tagPos = FirstPosition[args, "--tag", Missing[], {1}];
  If[MissingQ[tagPos],
    $HistoryTag = None;
    $SelectedTier = SelectFirst[args, !StringStartsQ[#, "--"]&, "all"],
    (* has --tag *)
    $HistoryTag = If[Length[args] >= First[tagPos] + 1,
      args[[First[tagPos] + 1]], "untagged"];
    $SelectedTier = SelectFirst[args[[;; First[tagPos] - 1]],
      !StringStartsQ[#, "--"]&, "all"]
  ]
];

Print["Selected tier: ", $SelectedTier];
If[$HistoryTag =!= None, Print["History tag:   ", $HistoryTag]];
Print["Subkernels:    ", $KernelCount, " (of ", $ProcessorCount, " processors)"];
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
      "NonLinearSolver",
        Module[{critResult},
          (* NonLinearSolver needs critical solution first *)
          critResult = SafeExecute[CriticalCongestionSolver[d2e], timeout];
          If[critResult["Status"] === "OK",
            Block[{V = Function[{x, edge}, 0]},
              SafeExecute[NonLinearSolver[critResult["Result"]], timeout]
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
        "NonLinearSolver", ComputeResidual[result["Result"], "AssoNonCritical"],
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
    data = Quiet @ Check[GetExampleData[key], $Failed];
    If[data === $Failed,
      Print["  FAILED to load GetExampleData[", key, "]"];
      Return[{<|"Key" -> ToString[key], "Tier" -> tier, "Solver" -> "ALL",
               "Status" -> "FAILED", "Timestamp" -> DateString[]|>}]
    ];

    (* Apply parameter substitutions *)
    If[params =!= {},
      data = data /. params
    ];

    (* Convert to equations (timed separately) *)
    Print["  Running DataToEquations..."];
    d2eResult = SafeExecute[DataToEquations[data], timeout];
    Print["    D2E Status: ", d2eResult["Status"],
          " | Time: ", FormatTime[d2eResult["Time"]], "s"];

    If[d2eResult["Status"] =!= "OK",
      Return[{<|"Key" -> ToString[key], "Tier" -> tier, "Solver" -> "DataToEquations",
               "D2ETime" -> d2eResult["Time"], "Status" -> d2eResult["Status"],
               "Timestamp" -> DateString[]|>}]
    ];

    d2e = d2eResult["Result"];

    (* Run each solver *)
    AppendTo[results,
      BenchmarkOneSolver[key, tier, "CriticalCongestion", d2e, timeout]];

    AppendTo[results,
      BenchmarkOneSolver[key, tier, "NonLinearSolver", d2e, timeout]];

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

(* --- Append to PARALLEL_PERFORMANCE_HISTORY.md if --tag was given --- *)
If[$HistoryTag =!= None,
  Module[{histPath, commitResult, commitHash, date, tier, rows, tierRows,
          entry, lines, existingContent, insertPos},

    histPath = FileNameJoin[{Directory[], "PARALLEL_PERFORMANCE_HISTORY.md"}];
    commitResult = RunProcess[{"git", "rev-parse", "--short", "HEAD"}];
    commitHash = StringTrim[commitResult["StandardOutput"]];
    date = DateString[{"Year", "-", "Month", "-", "Day"}];

    (* Build per-tier tables *)
    tier = $SelectedTier;
    rows = Select[$AllResults, #["Tier"] === tier || tier === "all"&];

    (* Group by case key *)
    BuildTierTable[tierName_] := Module[{tierRows2, caseKeys, tableLines},
      tierRows2 = Select[$AllResults, #["Tier"] === tierName &];
      If[Length[tierRows2] === 0, Return[""]];
      caseKeys = DeleteDuplicates[#["Key"] & /@ tierRows2];
      tableLines = {"#### Tier: " <> tierName, "",
        "| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Total (s) | Status |",
        "|------|---------|---------------|-------------|-----------|--------|"};
      Do[
        Module[{caseRows, d2eT, critT, nlT, totalT, status},
          caseRows = Select[tierRows2, #["Key"] === k &];
          d2eT = SelectFirst[caseRows, True&, <||>]["D2ETime"];
          critT = SelectFirst[Select[caseRows, #["Solver"] === "CriticalCongestion"&], True&, <||>];
          nlT  = SelectFirst[Select[caseRows, #["Solver"] === "NonLinearSolver"&],   True&, <||>];
          fmtT[t_] := If[NumericQ[t], ToString[CForm[N[t]]], "—"];
          d2eT   = fmtT[d2eT];
          critST = fmtT[If[KeyExistsQ[critT, "SolveTime"], critT["SolveTime"], Missing[]]];
          nlST   = fmtT[If[KeyExistsQ[nlT,   "SolveTime"], nlT["SolveTime"],   Missing[]]];
          totalN = Plus @@ Select[{
            If[KeyExistsQ[critT, "SolveTime"] && NumericQ[critT["SolveTime"]], critT["SolveTime"], 0],
            If[KeyExistsQ[nlT,   "SolveTime"] && NumericQ[nlT["SolveTime"]],   nlT["SolveTime"],   0]
          }, NumericQ];
          totalT = If[totalN > 0, fmtT[totalN], "—"];
          status = SelectFirst[caseRows, True&, <||>]["Status"];
          AppendTo[tableLines,
            "| " <> k <> " | " <> d2eT <> " | " <> critST <> " | " <>
            nlST <> " | " <> totalT <> " | " <> ToString[status] <> " |"]
        ],
        {k, caseKeys}
      ];
      StringRiffle[tableLines, "\n"]
    ];

    tiers = If[$SelectedTier === "all",
      Keys[$BenchmarkTiers],
      {$SelectedTier}
    ];
    tierTables = StringRiffle[
      Select[BuildTierTable /@ tiers, StringLength[#] > 0&],
      "\n\n"
    ];

    entry = StringRiffle[{
      "### " <> date <> " \[LongDash] " <> $HistoryTag,
      "",
      "**Commit:** `" <> commitHash <> "`",
      "**Date:** " <> DateString[],
      "**Kernels:** " <> ToString[$KernelCount],
      "",
      tierTables,
      "",
      "#### Rationale",
      "_" <> $HistoryTag <> "_",
      "",
      "#### Changes",
      "_TODO_",
      "",
      "#### Interpretation",
      "_TODO_",
      ""
    }, "\n"];

    (* Insert before the "## Future entries" section *)
    existingContent = Import[histPath, "Text"];
    insertPos = StringPosition[existingContent, "## Future entries"];
    If[Length[insertPos] > 0,
      existingContent = StringInsert[existingContent, entry <> "\n---\n\n", insertPos[[1,1]]],
      existingContent = existingContent <> "\n---\n\n" <> entry
    ];
    Export[histPath, existingContent, "Text"];
    Print["History entry appended to PARALLEL_PERFORMANCE_HISTORY.md (tag: ", $HistoryTag, ")"];
  ]
(* :!CodeAnalysis::EndBlock:: *)
];
