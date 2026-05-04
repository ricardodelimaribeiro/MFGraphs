#!/usr/bin/env wolframscript
(* Usage: wolframscript -file Scripts/BenchmarkSystemSolver.wls [--solver dnf|optimizeddnf|activeset|reduce|boolean|findinstance] [--tag label] [--timeout s] [--case name] *)
(* Benchmarks the active scenario-kernel solver path and stores solutions to Results/. *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

Print["MFGraphs System Solver Benchmark"];
Print["================================"];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print["OS: ", $SystemID];
Print[];

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$HistoryTag = None;
$TimeoutSecs = 60;
$SelectedCase = All;
$SolverName = "dnf";

Module[{args = Rest[$ScriptCommandLine], i = 1, arg},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--tag",
        i++; If[i <= Length[args], $HistoryTag = args[[i]]],
      "--timeout",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[NumericQ[v] && v > 0, $TimeoutSecs = v]]],
      "--case",
        i++; If[i <= Length[args], $SelectedCase = args[[i]]],
      "--solver",
        i++; If[i <= Length[args], $SolverName = ToLowerCase[args[[i]]]],
      _, Null
    ];
    i++;
  ]
];

$SolverNames = {"dnf", "optimizeddnf", "activeset", "reduce", "boolean", "findinstance"};

If[!MemberQ[$SolverNames, $SolverName],
  Print["Unknown solver: ", $SolverName];
  Print["Available solvers: ", StringRiffle[$SolverNames, ", "]];
  Exit[2]
];

RunSelectedSolver[sys_] :=
  Switch[$SolverName,
    "dnf", solversTools`dnfReduceSystem[sys],
    "optimizeddnf", solversTools`optimizedDNFReduceSystem[sys],
    "activeset", solversTools`activeSetReduceSystem[sys],
    "reduce", solversTools`reduceSystem[sys],
    "boolean", solversTools`booleanReduceSystem[sys],
    "findinstance", solversTools`findInstanceSystem[sys]
  ];

Print["Solver: ", $SolverName];
Print["Timeout per case: ", $TimeoutSecs, "s"];
If[$HistoryTag =!= None, Print["Tag: ", $HistoryTag]];
If[$SelectedCase =!= All, Print["Case filter: ", $SelectedCase]];
Print[];

$BenchCases = {
  {"chain-2v", Function[gridScenario[{2}, {{1, 10}}, {{2, 0}}]]},
  {"chain-3v-1exit", Function[gridScenario[{3}, {{1, 10}}, {{3, 0}}]]},
  {"chain-3v-2exit", Function[gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}]]},
  {"example-7", Function[getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}]]},
  {"chain-5v-1exit", Function[gridScenario[{5}, {{1, 10}}, {{5, 0}}]]},
  {"grid-2x3", Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]]},
  {"grid-3x2", Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]]},
  {"example-12", Function[getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}]]}
};
$CaseNames = $BenchCases[[All, 1]];

If[$SelectedCase =!= All,
  $BenchCases = Select[$BenchCases, #[[1]] === $SelectedCase &];
  If[$BenchCases === {},
    Print["No case named '", $SelectedCase, "'. Available: ",
      StringRiffle[$CaseNames, ", "]];
    Exit[2]
  ]
];

$ResultsDir = FileNameJoin[{$RepoRoot, "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];

$Timestamp = DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];
$ResultStem = "system_solver_" <> $SolverName;
$SolFile = FileNameJoin[{$ResultsDir, $ResultStem <> "_solutions_" <> $Timestamp <> ".wl"}];
$LatestSolFile = FileNameJoin[{$ResultsDir, $ResultStem <> "_solutions_latest.wl"}];
$CSVFile = FileNameJoin[{$ResultsDir, $ResultStem <> "_" <> $Timestamp <> ".csv"}];
$LatestCSVFile = FileNameJoin[{$ResultsDir, $ResultStem <> "_latest.csv"}];

FmtMs[t_?NumericQ] := ToString[Round[N[t * 1000], 0.01]] <> " ms";
FmtMs[_] := "-";
FmtMsCell[t_?NumericQ] := ToString[Round[N[t * 1000], 0.01]];
FmtMsCell[_] := "-";

ResultKind[sol_] := solversTools`Private`solutionResultKind[sol] /. "Unknown" -> "Other";

RunCase[{name_String, buildFn_}] :=
  Module[{sys, warmupTime, repTime, solWarmup, repSol, sol, buildMs,
          warmupMs, repMs, kind, valid, timeout, varDiag, transitionStatus,
          transitionResiduals},
    timeout = $TimeoutSecs;
    Print["Case: ", name];

    {buildMs, sys} = AbsoluteTiming[makeSystem[buildFn[]]];
    Print["  Build: ", FmtMs[buildMs]];

    {warmupTime, solWarmup} = AbsoluteTiming[
      TimeConstrained[RunSelectedSolver[sys], timeout, $TimedOut]
    ];
    If[solWarmup === $TimedOut,
      Print["  Warmup: TIMED OUT after ", FmtMs[warmupTime]];
      Return[<|"Name" -> name, "BuildMs" -> buildMs * 1000,
        "WarmupMs" -> warmupTime * 1000, "RepMs" -> Missing["TimedOut"],
        "Status" -> "TIMEOUT", "Kind" -> "-", "Valid" -> Missing["TimedOut"],
        "Solution" -> $TimedOut|>]
    ];
    warmupMs = warmupTime * 1000;

    {repTime, repSol} = AbsoluteTiming[
      TimeConstrained[RunSelectedSolver[sys], timeout, $TimedOut]
    ];
    sol = solWarmup;
    kind = ResultKind[sol];
    varDiag = solversTools`Private`solutionVariableDiagnostics[sys, sol];
    transitionStatus = Lookup[varDiag, "TransitionFlowStatus", "-"];
    transitionResiduals = Lookup[varDiag, "TransitionFlowResidualCount", "-"];
    repMs = If[repSol === $TimedOut, Missing["TimedOut"], repTime * 1000];
    valid = If[sol === $TimedOut, Missing["TimedOut"], isValidSystemSolution[sys, sol]];

    Print["  Warmup: ", FmtMs[warmupTime],
      " | Repeat: ", If[repSol === $TimedOut, "TIMED OUT", FmtMs[repTime]],
      " | Kind: ", kind,
      " | TransitionFlows: ", transitionStatus, " (residual ", transitionResiduals, ")",
      " | Valid: ", valid];

    <|"Name" -> name, "BuildMs" -> buildMs * 1000,
      "WarmupMs" -> warmupMs, "RepMs" -> repMs,
      "Status" -> If[sol === $TimedOut, "TIMEOUT", "OK"],
      "Kind" -> kind, "TransitionFlowStatus" -> transitionStatus,
      "TransitionFlowResidualCount" -> transitionResiduals,
      "Valid" -> valid, "Solution" -> sol|>
  ];

Print["Running ", Length[$BenchCases], " cases...\n"];
$AllResults = RunCase /@ $BenchCases;

$SolutionStore = Association[
  Flatten[{
    "Timestamp" -> DateString[],
    "Commit" -> StringTrim[RunProcess[{"git", "rev-parse", "--short", "HEAD"}]["StandardOutput"]],
    "Solver" -> $SolverName,
    "Cases" -> Association[
      Table[
        r["Name"] -> <|"Status" -> r["Status"], "Kind" -> r["Kind"],
          "TransitionFlowStatus" -> r["TransitionFlowStatus"],
          "TransitionFlowResidualCount" -> r["TransitionFlowResidualCount"],
          "Valid" -> r["Valid"], "Solution" -> r["Solution"]|>,
        {r, $AllResults}
      ]
    ]
  }]
];

Export[$SolFile, $SolutionStore, "WL"];
Export[$LatestSolFile, $SolutionStore, "WL"];
Print["\nSolutions written to: ", $SolFile];

$CSVRows = Join[
  {{"Case", "Solver", "BuildMs", "WarmupMs", "RepMs", "Status", "Kind", "TransitionFlowStatus", "TransitionFlowResidualCount", "Valid"}},
  Table[
    {r["Name"], $SolverName, r["BuildMs"], r["WarmupMs"],
      Replace[r["RepMs"], Missing[_] -> ""],
      r["Status"], r["Kind"], r["TransitionFlowStatus"], r["TransitionFlowResidualCount"],
      Replace[r["Valid"], Missing[_] -> ""]},
    {r, $AllResults}
  ]
];
Export[$CSVFile, $CSVRows, "CSV"];
Export[$LatestCSVFile, $CSVRows, "CSV"];
Print["CSV written to:       ", $CSVFile];

Print["\n| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |"];
Print["|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|"];
Do[
  Print["| ", r["Name"], " | ", $SolverName, " | ",
    FmtMsCell[r["BuildMs"]/1000], " | ",
    FmtMsCell[r["WarmupMs"]/1000], " | ",
    Replace[r["RepMs"], {Missing[_] -> "-", v_?NumericQ :> FmtMsCell[v/1000]}], " | ",
    r["Status"], " | ", r["Kind"], " | ", r["TransitionFlowStatus"], " | ",
    r["TransitionFlowResidualCount"], " | ",
    Replace[r["Valid"], {Missing[_] -> "-", v_ :> ToString[v]}], " |"],
  {r, $AllResults}
];

If[$HistoryTag =!= None,
  Module[{histPath, commitHash, date, tableLines, entry, existing,
          insertPos},
    histPath = FileNameJoin[{$RepoRoot, "BENCHMARKS.md"}];
    commitHash = StringTrim[RunProcess[{"git", "rev-parse", "--short", "HEAD"}]["StandardOutput"]];
    date = DateString[{"Year", "-", "Month", "-", "Day"}];

    tableLines = {"| Case | Solver | Build (ms) | Warmup (ms) | Repeat (ms) | Status | Kind | Transition flows | Residual Jts | Valid |",
      "|------|--------|-----------|------------|--------------|--------|------|------------------|--------------|-------|"};
    Do[
      AppendTo[tableLines,
        "| " <> r["Name"] <> " | " <> $SolverName <> " | " <>
        FmtMsCell[r["BuildMs"]/1000] <> " | " <>
        FmtMsCell[r["WarmupMs"]/1000] <> " | " <>
        Replace[r["RepMs"], {Missing[_] -> "-", v_?NumericQ :> FmtMsCell[v/1000]}] <> " | " <>
        r["Status"] <> " | " <> r["Kind"] <> " | " <>
        ToString[r["TransitionFlowStatus"]] <> " | " <>
        ToString[r["TransitionFlowResidualCount"]] <> " | " <>
        Replace[r["Valid"], {Missing[_] -> "-", v_ :> ToString[v]}] <> " |"],
      {r, $AllResults}
    ];

    entry = StringRiffle[{
      "### " <> date <> " - " <> $HistoryTag,
      "",
      "**Commit:** `" <> commitHash <> "`  ",
      "**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  ",
      "**Solver:** `" <> $SolverName <> "`  ",
      "**Timeout per case:** " <> ToString[$TimeoutSecs] <> "s  ",
      "**Solutions file:** `Results/" <> $ResultStem <> "_solutions_" <> $Timestamp <> ".wl`",
      "",
      StringRiffle[tableLines, "\n"],
      ""
    }, "\n"];

    existing = Import[histPath, "Text"];
    insertPos = StringPosition[existing, "## system solver benchmark history"];
    If[Length[insertPos] > 0,
      Module[{afterHeader, nextEntry},
        afterHeader = StringDrop[existing, insertPos[[1, 2]]];
        nextEntry = StringPosition[afterHeader, "\n### "];
        If[Length[nextEntry] > 0,
          existing = StringInsert[existing, "\n" <> entry <> "---\n",
            insertPos[[1, 2]] + nextEntry[[1, 1]]],
          existing = existing <> "\n---\n\n" <> entry
        ]
      ],
      existing = existing <> "\n\n## system solver benchmark history\n\n---\n\n" <> entry
    ];
    Export[histPath, existing, "Text"];
    Print["\nHistory entry appended to BENCHMARKS.md (tag: ", $HistoryTag, ")"];
  ]
];

Print["\nDone at ", DateString[]];

(* :!CodeAnalysis::EndBlock:: *)
