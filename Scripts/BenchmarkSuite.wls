(* Usage: wolframscript -file Scripts/BenchmarkSuite.wls [tier] [case=key] [timeout=seconds] [tag=label] *)
(* tier: "small", "core", "stress", "large", "vlarge", legacy "medium", or "all" (default: "all") *)

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
(* Support: wolframscript -file Scripts/BenchmarkSuite.wls [tier] [case=key] [timeout=seconds] [tag=description] *)
NormalizeCaseArg[arg_String] := If[StringMatchQ[arg, NumberString], ToExpression[arg], arg];
NormalizeCaseArg[arg_] := arg;

$HistoryTag = None;
$SelectedTier = "all";
$SelectedCase = All;
$TimeoutOverride = Automatic;

Module[{args = Rest[$ScriptCommandLine], i = 1, arg, parsedTimeout},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--tag",
        i++;
        $HistoryTag = If[i <= Length[args], args[[i]], "untagged"],
      "--case",
        i++;
        If[i <= Length[args], $SelectedCase = NormalizeCaseArg[args[[i]]]],
      "--timeout",
        i++;
        If[i <= Length[args],
          parsedTimeout = Quiet @ Check[ToExpression[args[[i]]], $Failed];
          If[NumericQ[parsedTimeout] && parsedTimeout > 0,
            $TimeoutOverride = parsedTimeout
          ]
        ],
      _ /; StringStartsQ[arg, "case="],
        $SelectedCase = NormalizeCaseArg[StringDrop[arg, StringLength["case="]]],
      _ /; StringStartsQ[arg, "timeout="],
        parsedTimeout = Quiet @ Check[ToExpression[StringDrop[arg, StringLength["timeout="]]], $Failed];
        If[NumericQ[parsedTimeout] && parsedTimeout > 0,
          $TimeoutOverride = parsedTimeout
        ],
      _ /; StringStartsQ[arg, "tag="],
        $HistoryTag = StringDrop[arg, StringLength["tag="]],
      _,
        If[!StringStartsQ[arg, "--"], $SelectedTier = arg]
    ];
    i++;
  ];
];

$RequestedTier = $SelectedTier;
$SelectedTier = Lookup[$BenchmarkTierAliases, $SelectedTier, $SelectedTier];

Print["Selected tier: ", $SelectedTier];
If[$RequestedTier =!= $SelectedTier, Print["Requested alias: ", $RequestedTier]];
If[$SelectedCase =!= All, Print["Selected case: ", $SelectedCase]];
If[NumericQ[$TimeoutOverride], Print["Timeout override: ", $TimeoutOverride, "s"]];
If[$HistoryTag =!= None, Print["History tag:   ", $HistoryTag]];
Print["Subkernels:    ", $KernelCount, " (of ", $ProcessorCount, " processors)"];
Print[];

(* --- Determine which cases to run --- *)
$CasesToRun = If[$SelectedCase =!= All,
  {$SelectedCase},
  If[$SelectedTier === "all",
    Flatten[Lookup[$BenchmarkTiers, $DefaultBenchmarkTiers]],
    Lookup[$BenchmarkTiers, $SelectedTier, {}]
  ]
];

GetTier[key_] := SelectFirst[Keys[$BenchmarkTiers], MemberQ[$BenchmarkTiers[#], key] &, "unknown"];

(* --- Benchmark a single case with a single solver --- *)

BenchmarkOneSolver[key_, tier_String, solverName_String, d2e_Association, timeout_?NumericQ] :=
  Module[{result, meta, residual = Missing["N/A"], equationResidual = Missing["N/A"]},
    meta = GraphMetadata[d2e];
    Print["  Solver: ", solverName, " (timeout: ", timeout, "s)"];

    result = Switch[solverName,
      "CriticalCongestion",
        SafeExecute[CriticalCongestionSolver[d2e, "ReturnShape" -> "Standard"], timeout],
      "NonLinearSolver",
        Module[{critResult},
          (* NonLinearSolver needs critical solution first *)
          critResult = SafeExecute[CriticalCongestionSolver[d2e, "ReturnShape" -> "Standard"], timeout];
          If[critResult["Status"] === "OK",
            Block[{V = Function[{x, edge}, 0]},
              SafeExecute[NonLinearSolver[critResult["Result"], "ReturnShape" -> "Standard"], timeout]
            ],
            <|"Time" -> 0, "Memory" -> 0, "Status" -> "SKIPPED",
              "Result" -> "Critical solver failed"|>
          ]
        ],
      "Monotone",
        SafeExecute[MonotoneSolver[d2e, "ReturnShape" -> "Standard"], timeout],
      _,
        <|"Time" -> 0, "Memory" -> 0, "Status" -> "SKIPPED", "Result" -> "Unknown solver"|>
    ];

    (* Compute residual for successful solves *)
    If[result["Status"] === "OK" && AssociationQ[result["Result"]],
      residual = Lookup[result["Result"], "KirchhoffResidual", Missing["N/A"]];
      equationResidual = Switch[solverName,
        "CriticalCongestion", ComputeResidual[result["Result"], "AssoCritical"],
        "NonLinearSolver", ComputeResidual[result["Result"], "AssoNonCritical"],
        "Monotone", Missing["N/A"],
        _, Missing["N/A"]
      ]
    ];

    (* Verify against reference hashes for CriticalCongestion *)
    Module[{verification = "N/A"},
      If[result["Status"] === "OK" && solverName === "CriticalCongestion",
        verification = VerifySolution[key,
          Lookup[result["Result"], "AssoCritical", Null]];
      ];
      Print["    Status: ", result["Status"],
            " | Time: ", FormatTime[result["Time"]], "s",
            " | Memory: ", FormatMemoryMB[result["Memory"]], " MB",
            If[verification =!= "N/A", " | Ref: " <> verification, ""]];
    ];

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
      "EquationResidual" -> equationResidual,
      "Timestamp" -> DateString[]
    |>
  ];

CaseStatusSummary[caseRows_List] :=
  Module[{badRows},
    badRows = Select[caseRows, Lookup[#, "Status", "UNKNOWN"] =!= "OK" &];
    If[badRows === {},
      "OK",
      StringRiffle[
        (Lookup[#, "Solver", "Unknown"] <> ":" <> ToString[Lookup[#, "Status", "UNKNOWN"]]) & /@ badRows,
        ", "
      ]
    ]
  ];

(* --- Benchmark a single case across all applicable solvers --- *)

BenchmarkCase[key_] :=
  Module[{tier, timeout, data, d2eResult, d2e, results = {}, params},
    tier = GetTier[key];
    timeout = If[NumericQ[$TimeoutOverride], $TimeoutOverride, Lookup[$TierTimeouts, tier, 300]];
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

(* --- Export results --- *)
$ResultsDir = FileNameJoin[{Directory[], "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];

$Timestamp = DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];
$CSVPath = FileNameJoin[{$ResultsDir, "benchmark_" <> $Timestamp <> ".csv"}];
$JSONPath = FileNameJoin[{$ResultsDir, "benchmark_" <> $Timestamp <> ".json"}];
$LatestCSV = FileNameJoin[{$ResultsDir, "benchmark_latest.csv"}];
$LatestJSON = FileNameJoin[{$ResultsDir, "benchmark_latest.json"}];

ExportBenchmarkSnapshot[results_List] := Module[{},
  ExportResultsCSV[results, $CSVPath];
  ExportResultsJSON[results, $JSONPath];
  ExportResultsCSV[results, $LatestCSV];
  ExportResultsJSON[results, $LatestJSON];
];

$AllResults = {};
Do[
  $AllResults = Join[$AllResults, BenchmarkCase[key]];
  ExportBenchmarkSnapshot[$AllResults],
  {key, $CasesToRun}
];

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
        "| Case | D2E (s) | CritSolver (s) | NLSolver (s) | Monotone (s) | Total (s) | Case Status |",
        "|------|---------|---------------|-------------|--------------|-----------|-------------|"};
      Do[
        Module[{caseRows, d2eT, critT, nlT, monoT, totalT, caseStatus},
          caseRows = Select[tierRows2, #["Key"] === k &];
          d2eT = SelectFirst[caseRows, True&, <||>]["D2ETime"];
          critT = SelectFirst[Select[caseRows, #["Solver"] === "CriticalCongestion"&], True&, <||>];
          nlT  = SelectFirst[Select[caseRows, #["Solver"] === "NonLinearSolver"&], True&, <||>];
          monoT = SelectFirst[Select[caseRows, #["Solver"] === "Monotone"&], True&, <||>];
          fmtT[t_] := If[NumericQ[t], ToString[CForm[N[t]]], "—"];
          d2eT   = fmtT[d2eT];
          critST = fmtT[If[KeyExistsQ[critT, "SolveTime"], critT["SolveTime"], Missing[]]];
          nlST   = fmtT[If[KeyExistsQ[nlT, "SolveTime"], nlT["SolveTime"], Missing[]]];
          monoST = fmtT[If[KeyExistsQ[monoT, "SolveTime"], monoT["SolveTime"], Missing[]]];
          totalN = Plus @@ Select[{
            If[KeyExistsQ[critT, "SolveTime"] && NumericQ[critT["SolveTime"]], critT["SolveTime"], 0],
            If[KeyExistsQ[nlT, "SolveTime"] && NumericQ[nlT["SolveTime"]], nlT["SolveTime"], 0],
            If[KeyExistsQ[monoT, "SolveTime"] && NumericQ[monoT["SolveTime"]], monoT["SolveTime"], 0]
          }, NumericQ];
          totalT = If[totalN > 0, fmtT[totalN], "—"];
          caseStatus = CaseStatusSummary[caseRows];
          AppendTo[tableLines,
            "| " <> k <> " | " <> d2eT <> " | " <> critST <> " | " <>
            nlST <> " | " <> monoST <> " | " <> totalT <> " | " <> caseStatus <> " |"]
        ],
        {k, caseKeys}
      ];
      StringRiffle[tableLines, "\n"]
    ];

    tiers = If[$SelectedTier === "all",
      $DefaultBenchmarkTiers,
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
      "_Describe the code or benchmark configuration change being measured in this run._",
      "",
      "#### Interpretation",
      "_Summarize the main performance signal, regressions, or follow-up actions from this run._",
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
