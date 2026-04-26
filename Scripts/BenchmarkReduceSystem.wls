#!/usr/bin/env wolframscript
(* Usage: wolframscript -file Scripts/BenchmarkReduceSystem.wls [--tag label] [--timeout s] [--case name] *)
(* Benchmarks reduceSystem across representative cases and stores solutions to a file. *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

Print["MFGraphs reduceSystem Benchmark"];
Print["================================"];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print["OS: ", $SystemID];
Print[];

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

(* --- Parse arguments --- *)
$HistoryTag  = None;
$TimeoutSecs = 60;
$SelectedCase = All;

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
      _, Null
    ];
    i++;
  ]
];

Print["Timeout per case: ", $TimeoutSecs, "s"];
If[$HistoryTag =!= None, Print["Tag: ", $HistoryTag]];
If[$SelectedCase =!= All, Print["Case filter: ", $SelectedCase]];
Print[];

(* --- Case definitions --- *)
(* Each entry: {name, buildFn}  where buildFn[] returns an mfgSystem. *)

$BenchCases = {
  {"chain-2v",
    Function[gridScenario[{2}, {{1, 10}}, {{2, 0}}]] },
  {"chain-3v-1exit",
    Function[gridScenario[{3}, {{1, 10}}, {{3, 0}}]] },
  {"chain-3v-2exit",
    Function[gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}]] },
  {"example-7",
    Function[getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}]] },
  {"chain-5v-1exit",
    Function[gridScenario[{5}, {{1, 10}}, {{5, 0}}]] },
  {"example-12",
    Function[getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}]] }
};

If[$SelectedCase =!= All,
  $BenchCases = Select[$BenchCases, #[[1]] === $SelectedCase &];
  If[$BenchCases === {},
    Print["No case named '", $SelectedCase, "'. Available: ",
          StringRiffle[$BenchCases[[All, 1]], ", "]];
    Exit[2]
  ]
];

(* --- Results storage --- *)
$ResultsDir = FileNameJoin[{$RepoRoot, "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];

$Timestamp    = DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];
$SolFile      = FileNameJoin[{$ResultsDir, "reduce_system_solutions_" <> $Timestamp <> ".wl"}];
$LatestSolFile = FileNameJoin[{$ResultsDir, "reduce_system_solutions_latest.wl"}];
$CSVFile      = FileNameJoin[{$ResultsDir, "reduce_system_" <> $Timestamp <> ".csv"}];
$LatestCSVFile = FileNameJoin[{$ResultsDir, "reduce_system_latest.csv"}];

(* --- Benchmark helpers --- *)

FmtMs[t_?NumericQ] := ToString[Round[N[t * 1000], 0.01]] <> " ms";
FmtMs[_] := "—";
FmtMsCell[t_?NumericQ] := ToString[Round[N[t * 1000], 0.01]];
FmtMsCell[_] := "—";

ResultKind[sol_] := Which[
  ListQ[sol] && AllTrue[sol, MatchQ[#, _Rule] &], "Rules",
  AssociationQ[sol] && KeyExistsQ[sol, "Rules"], "Underdetermined",
  sol === {}, "NoSolution",
  True, "Other"
];

(* --- Run one case --- *)

RunCase[{name_String, buildFn_}] :=
  Module[{sys, warmupTime, repTime, solWarmup, sol,
          buildMs, warmupMs, repMs, kind, valid, timeout, result},
    timeout = $TimeoutSecs;
    Print["Case: ", name];

    (* Build system (not timed in benchmark — one-time cost) *)
    {buildMs, sys} = AbsoluteTiming[makeSystem[buildFn[]]];
    Print["  Build: ", FmtMs[buildMs]];

    (* Warmup run with timeout *)
    {warmupTime, solWarmup} = AbsoluteTiming[
      TimeConstrained[reduceSystem[sys], timeout, $TimedOut]
    ];
    If[solWarmup === $TimedOut,
      Print["  Warmup: TIMED OUT after ", FmtMs[warmupTime]];
      Return[<|"Name" -> name, "BuildMs" -> buildMs * 1000,
               "WarmupMs" -> warmupTime * 1000, "RepMs" -> Missing["TimedOut"],
               "Status" -> "TIMEOUT", "Kind" -> "—", "Valid" -> Missing["TimedOut"],
               "Solution" -> $TimedOut|>]
    ];
    warmupMs = warmupTime * 1000;

    (* Repeated timing — 3 reps, short timeout per rep *)
    repTime = RepeatedTiming[
      TimeConstrained[reduceSystem[sys], timeout, $TimedOut],
      3
    ][[1]];
    sol  = reduceSystem[sys];  (* final solution to store *)
    kind = ResultKind[sol];
    repMs = repTime * 1000;

    valid = isReduceSystemSolution[sys, sol];

    Print["  Warmup: ", FmtMs[warmupTime],
          " | Repeated: ", FmtMs[repTime],
          " | Kind: ", kind,
          " | Valid: ", valid];

    <|"Name" -> name, "BuildMs" -> buildMs * 1000,
      "WarmupMs" -> warmupMs, "RepMs" -> repMs,
      "Status" -> "OK", "Kind" -> kind, "Valid" -> valid, "Solution" -> sol|>
  ];

(* --- Main loop --- *)

Print["Running ", Length[$BenchCases], " cases...\n"];

$AllResults = RunCase /@ $BenchCases;

(* --- Store solutions --- *)

$SolutionStore = Association[
  Flatten[{
    "Timestamp" -> DateString[],
    "Commit"    -> StringTrim[RunProcess[{"git", "rev-parse", "--short", "HEAD"}]["StandardOutput"]],
    "Cases"     -> Association[
      Table[
        r["Name"] -> <|"Status" -> r["Status"], "Kind" -> r["Kind"],
                       "Valid" -> r["Valid"], "Solution" -> r["Solution"]|>,
        {r, $AllResults}
      ]
    ]
  }]
];

Export[$SolFile, $SolutionStore, "WL"];
Export[$LatestSolFile, $SolutionStore, "WL"];
Print["\nSolutions written to: ", $SolFile];

(* --- Store CSV --- *)

$CSVRows = Join[
  {{"Case", "BuildMs", "WarmupMs", "RepMs", "Status", "Kind", "Valid"}},
  Table[
    {r["Name"], r["BuildMs"], r["WarmupMs"],
     Replace[r["RepMs"], Missing[_] -> ""],
     r["Status"], r["Kind"],
     Replace[r["Valid"], Missing[_] -> ""]},
    {r, $AllResults}
  ]
];
Export[$CSVFile, $CSVRows, "CSV"];
Export[$LatestCSVFile, $CSVRows, "CSV"];
Print["CSV written to:       ", $CSVFile];

(* --- Print summary table --- *)

Print["\n| Case | Build (ms) | Warmup (ms) | Repeated (ms) | Status | Kind | Valid |"];
Print["|------|-----------|------------|--------------|--------|------|-------|"];
Do[
  Print["| ", r["Name"], " | ",
        FmtMsCell[r["BuildMs"]/1000], " | ",
        FmtMsCell[r["WarmupMs"]/1000], " | ",
        Replace[r["RepMs"], {Missing[_] -> "—", v_?NumericQ :> FmtMsCell[v/1000]}], " | ",
        r["Status"], " | ", r["Kind"], " | ",
        Replace[r["Valid"], {Missing[_] -> "—", v_ :> ToString[v]}], " |"],
  {r, $AllResults}
];

(* --- Append to BENCHMARKS.md if --tag given --- *)

If[$HistoryTag =!= None,
  Module[{histPath, commitHash, date, tableLines, entry,
          existing, insertPos},
    histPath   = FileNameJoin[{$RepoRoot, "BENCHMARKS.md"}];
    commitHash = StringTrim[RunProcess[{"git", "rev-parse", "--short", "HEAD"}]["StandardOutput"]];
    date       = DateString[{"Year", "-", "Month", "-", "Day"}];

    tableLines = {"| Case | Build (ms) | Warmup (ms) | Repeated (ms) | Status | Kind | Valid |",
                  "|------|-----------|------------|--------------|--------|------|-------|"};
    Do[
      AppendTo[tableLines,
        "| " <> r["Name"] <> " | " <>
        FmtMsCell[r["BuildMs"]/1000] <> " | " <>
        FmtMsCell[r["WarmupMs"]/1000] <> " | " <>
        Replace[r["RepMs"], {Missing[_] -> "—",
                             v_?NumericQ :> FmtMsCell[v/1000]}] <> " | " <>
        r["Status"] <> " | " <> r["Kind"] <> " | " <>
        Replace[r["Valid"], {Missing[_] -> "—", v_ :> ToString[v]}] <> " |"],
      {r, $AllResults}
    ];

    entry = StringRiffle[{
      "### " <> date <> " \[LongDash] " <> $HistoryTag,
      "",
      "**Commit:** `" <> commitHash <> "`  ",
      "**Environment:** local `wolframscript`, `$MFGraphsVerbose=False`  ",
      "**Timeout per case:** " <> ToString[$TimeoutSecs] <> "s  ",
      "**Solutions file:** `Results/reduce_system_solutions_" <> $Timestamp <> ".wl`",
      "",
      StringRiffle[tableLines, "\n"],
      ""
    }, "\n"];

    existing = Import[histPath, "Text"];
    insertPos = StringPosition[existing, "## reduceSystem benchmark history"];
    If[Length[insertPos] > 0,
      (* Insert after the section header and the next blank line *)
      Module[{afterHeader, nextEntry},
        afterHeader = StringDrop[existing, insertPos[[1, 2]]];
        nextEntry   = StringPosition[afterHeader, "\n### "];
        If[Length[nextEntry] > 0,
          existing = StringInsert[existing,
            "\n" <> entry <> "---\n",
            insertPos[[1, 2]] + nextEntry[[1, 1]]],
          existing = existing <> "\n---\n\n" <> entry
        ]
      ],
      existing = existing <> "\n\n## reduceSystem benchmark history (manual)\n\n---\n\n" <> entry
    ];
    Export[histPath, existing, "Text"];
    Print["\nHistory entry appended to BENCHMARKS.md (tag: ", $HistoryTag, ")"];
  ]
];

Print["\nDone at ", DateString[]];
(* :!CodeAnalysis::EndBlock:: *)
