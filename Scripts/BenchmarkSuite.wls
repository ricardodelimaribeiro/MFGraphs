(* Usage: wolframscript -file Scripts/BenchmarkSuite.wls [tier] [case=key] [timeout=seconds] [isolation=mode] [tag=label] [rationale=...] [changes=...] [interpretation=...] *)
(* tier: "small", "core", "stress", "paper", "inconsistent-switching", "large", "vlarge", legacy "medium", or "all" (default: "all") *)

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
(* Support: wolframscript -file Scripts/BenchmarkSuite.wls [tier] [case=key] [timeout=seconds] [isolation=process|inprocess|auto] [tag=description] *)
NormalizeCaseArg[arg_String] := If[StringMatchQ[arg, NumberString], ToExpression[arg], arg];
NormalizeCaseArg[arg_] := arg;
NormalizeIsolationMode[arg_String] := Switch[ToLowerCase[arg],
  "process" | "isolated", "Process",
  "inprocess" | "in-process" | "none", "InProcess",
  "auto", "Auto",
  _, $Failed
];
NormalizeIsolationMode[_] := $Failed;

$HistoryTag = None;
$HistoryRationale = None;
$HistoryChanges = None;
$HistoryInterpretation = None;
$SelectedTier = "all";
$SelectedCase = All;
$TimeoutOverride = Automatic;
$IsolationModeArg = "inprocess";

Module[{args = Rest[$ScriptCommandLine], i = 1, arg, parsedTimeout, parsedIsolation},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--tag",
        i++;
        $HistoryTag = If[i <= Length[args], args[[i]], "untagged"],
      "--rationale",
        i++;
        If[i <= Length[args], $HistoryRationale = args[[i]]],
      "--changes",
        i++;
        If[i <= Length[args], $HistoryChanges = args[[i]]],
      "--interpretation",
        i++;
        If[i <= Length[args], $HistoryInterpretation = args[[i]]],
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
      "--isolation",
        i++;
        If[i <= Length[args], $IsolationModeArg = args[[i]]],
      _ /; StringStartsQ[arg, "case="],
        $SelectedCase = NormalizeCaseArg[StringDrop[arg, StringLength["case="]]],
      _ /; StringStartsQ[arg, "timeout="],
        parsedTimeout = Quiet @ Check[ToExpression[StringDrop[arg, StringLength["timeout="]]], $Failed];
        If[NumericQ[parsedTimeout] && parsedTimeout > 0,
          $TimeoutOverride = parsedTimeout
        ],
      _ /; StringStartsQ[arg, "isolation="],
        $IsolationModeArg = StringDrop[arg, StringLength["isolation="]],
      _ /; StringStartsQ[arg, "tag="],
        $HistoryTag = StringDrop[arg, StringLength["tag="]],
      _ /; StringStartsQ[arg, "rationale="],
        $HistoryRationale = StringDrop[arg, StringLength["rationale="]],
      _ /; StringStartsQ[arg, "changes="],
        $HistoryChanges = StringDrop[arg, StringLength["changes="]],
      _ /; StringStartsQ[arg, "interpretation="],
        $HistoryInterpretation = StringDrop[arg, StringLength["interpretation="]],
      _,
        If[!StringStartsQ[arg, "--"], $SelectedTier = arg]
    ];
    i++;
  ];

  parsedIsolation = NormalizeIsolationMode[$IsolationModeArg];
  If[parsedIsolation === $Failed,
    Print["Invalid isolation mode: ", $IsolationModeArg];
    Print["Valid values: inprocess, process, auto"];
    Exit[2]
  ];
  $RequestedIsolationMode = parsedIsolation;
  $SafeExecuteMode = If[$RequestedIsolationMode === "Auto", "InProcess", $RequestedIsolationMode];
];

$RequestedTier = $SelectedTier;
$SelectedTier = Lookup[$BenchmarkTierAliases, $SelectedTier, $SelectedTier];

Print["Selected tier: ", $SelectedTier];
If[$RequestedTier =!= $SelectedTier, Print["Requested alias: ", $RequestedTier]];
If[$SelectedCase =!= All, Print["Selected case: ", $SelectedCase]];
If[NumericQ[$TimeoutOverride], Print["Timeout override: ", $TimeoutOverride, "s"]];
Print["Isolation mode: ", If[$RequestedIsolationMode === "Auto", "Auto (tier-routed)", $RequestedIsolationMode]];
If[$HistoryTag =!= None, Print["History tag:   ", $HistoryTag]];
If[$HistoryTag =!= None &&
   MemberQ[{$HistoryRationale, $HistoryChanges, $HistoryInterpretation}, None],
  Print["History entries require explicit rationale, changes, and interpretation text."];
  Exit[2]
];
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
ResolveCaseIsolationMode[tier_String] :=
  If[
    $RequestedIsolationMode === "Auto",
    If[MemberQ[{"large", "vlarge"}, tier], "Process", "InProcess"],
    $RequestedIsolationMode
  ];

(* --- Benchmark a single case with a single solver --- *)

BenchmarkOneSolver[key_, tier_String, solverName_String, d2e_Association, timeout_?NumericQ] :=
  Module[{result, meta, residual = Missing["N/A"], equationResidual = Missing["N/A"],
      convergence, stopReason = Missing["N/A"], convergenceResidual = Missing["N/A"],
      iterations = Missing["N/A"], solveTime = Missing["N/A"],
      solveMemory = Missing["N/A"], executionMode = Missing["N/A"],
      rawResult = Missing["NotAvailable"], resultKind = Missing["NotAvailable"],
      feasibility = Missing["NotAvailable"], solverMessage = Missing["NotAvailable"]},
    meta = GraphMetadata[d2e];
    Print["  Solver: ", solverName, " (timeout: ", timeout, "s)"];

    result = Switch[solverName,
      "CriticalCongestion",
        SafeExecute[CriticalCongestionSolver[d2e], timeout],
      "NonLinearSolver",
        SafeExecute[NonLinearSolver[d2e], timeout],
      "Monotone",
        SafeExecute[MonotoneSolver[d2e], timeout],
      _,
        <|"Time" -> 0, "Memory" -> 0, "Status" -> "SKIPPED", "Result" -> "Unknown solver"|>
    ];
    solveTime =
      If[
        NumericQ[Lookup[result, "InnerSolveTime", Missing["N/A"]]],
        Lookup[result, "InnerSolveTime", Missing["N/A"]],
        Lookup[result, "Time", Missing["N/A"]]
      ];
    solveMemory =
      If[
        NumericQ[Lookup[result, "ChildMaxMemoryUsed", Missing["N/A"]]],
        Lookup[result, "ChildMaxMemoryUsed", Missing["N/A"]],
        Lookup[result, "Memory", Missing["N/A"]]
      ];
    executionMode = Lookup[result, "ExecutionMode", If[$SafeExecuteMode === "Process", "Process", "InProcess"]];
    rawResult = Lookup[result, "Result", Missing["NotAvailable"]];
    If[AssociationQ[rawResult],
      resultKind = Lookup[rawResult, "ResultKind", Missing["NotAvailable"]];
      feasibility = Lookup[rawResult, "Feasibility", Missing["NotAvailable"]];
      solverMessage = Lookup[rawResult, "Message", Missing["NotAvailable"]];
    ];

    (* Compute residual for successful solves *)
    If[result["Status"] === "OK" && AssociationQ[rawResult],
      residual = Lookup[rawResult, "KirchhoffResidual", Missing["N/A"]];
      convergence = Lookup[rawResult, "Convergence", Missing["N/A"]];
      If[AssociationQ[convergence],
        stopReason = Lookup[convergence, "StopReason", Missing["N/A"]];
        convergenceResidual = Lookup[convergence, "FinalResidual", Missing["N/A"]];
        iterations = Lookup[convergence, "Iterations", Missing["N/A"]];
      ];
      equationResidual = Switch[solverName,
        "CriticalCongestion", ComputeResidual[Join[d2e, rawResult], "AssoCritical"],
        "NonLinearSolver", ComputeResidual[Join[d2e, rawResult], "AssoNonCritical"],
        "Monotone", ComputeResidual[Join[d2e, rawResult], "AssoMonotone"],
        _, Missing["N/A"]
      ]
    ];

    (* Verify against reference hashes for CriticalCongestion *)
    Module[{verification = "N/A"},
      If[result["Status"] === "OK" && solverName === "CriticalCongestion",
        verification = VerifySolution[key,
          Lookup[rawResult, "AssoCritical", Null]];
      ];
      Print["    Status: ", result["Status"],
            " | Time: ", FormatTime[solveTime], "s",
            " | Memory: ", FormatMemoryMB[solveMemory], " MB",
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
      "ExecutionMode" -> executionMode,
      "WallTime" -> Lookup[result, "WallTime", Lookup[result, "Time", Missing["N/A"]]],
      "InnerSolveTime" -> Lookup[result, "InnerSolveTime", Missing["N/A"]],
      "ChildMaxMemoryUsed" -> Lookup[result, "ChildMaxMemoryUsed", Missing["N/A"]],
      "D2ETime" -> Missing["Included"],
      "SolveTime" -> solveTime,
      "SolveMemory" -> solveMemory,
      "Status" -> result["Status"],
      "ResultKind" -> resultKind,
      "Feasibility" -> feasibility,
      "Message" -> solverMessage,
      "Residual" -> residual,
      "EquationResidual" -> equationResidual,
      "StopReason" -> stopReason,
      "ConvergenceResidual" -> convergenceResidual,
      "Iterations" -> iterations,
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
  Module[{tier, timeout, data, d2eResult, d2e, results = {}, params, d2eTime, caseIsolationMode},
    tier = GetTier[key];
    timeout = If[NumericQ[$TimeoutOverride], $TimeoutOverride, Lookup[$TierTimeouts, tier, 300]];
    params = GetParams[key];
    caseIsolationMode = ResolveCaseIsolationMode[tier];
    $SafeExecuteMode = caseIsolationMode;

    Print["\n--- Case: ", key, " (tier: ", tier, ") ---"];
    Print["  Isolation: ", caseIsolationMode];

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
    d2eTime =
      If[
        NumericQ[Lookup[d2eResult, "InnerSolveTime", Missing["N/A"]]],
        Lookup[d2eResult, "InnerSolveTime", Missing["N/A"]],
        Lookup[d2eResult, "Time", Missing["N/A"]]
      ];
    Print["    D2E Status: ", d2eResult["Status"],
          " | Time: ", FormatTime[d2eTime], "s"];

    If[d2eResult["Status"] =!= "OK",
      Return[{<|"Key" -> ToString[key], "Tier" -> tier, "Solver" -> "DataToEquations",
               "D2ETime" -> d2eTime, "Status" -> d2eResult["Status"],
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
    results = Map[Append[#, "D2ETime" -> d2eTime] &, results];

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
      $HistoryRationale,
      "",
      "#### Changes",
      $HistoryChanges,
      "",
      "#### Interpretation",
      $HistoryInterpretation,
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
