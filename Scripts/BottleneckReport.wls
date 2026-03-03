#!/usr/bin/env wolframscript
(* BottleneckReport.wls -- Run profiled benchmarks and generate a bottleneck report *)
(* Usage: wolframscript -file Scripts/BottleneckReport.wls *)

Print["MFGraphs Bottleneck Analysis"];
Print["=============================="];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print[];

(* --- Setup --- *)
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];
Get[FileNameJoin[{Directory[], "Scripts", "ProfileInstrument.wls"}]];

(* --- Select representative cases for profiling --- *)
(* One from each complexity class *)
$ProfileCases = <|
  "small_simple" -> 3,                   (* 2-edge linear chain *)
  "medium_nosw" -> 7,                    (* Y-network, no switching *)
  "medium_sw" -> 8,                      (* Y-network, with switching *)
  "medium_attraction" -> 12,             (* Attraction, no switching *)
  "medium_triangle" -> 14,               (* Triangle with switching *)
  "large_multipath" -> 20,               (* 9 vertices, 2 in, 3 out *)
  "large_braess" -> "Braess congest"     (* Braess paradox *)
|>;

(* --- Run profiled benchmarks --- *)

$AllProfiles = <||>;
$AllTimings = {};

Do[
  Module[{key, data, d2e, params, critResult, t0, totalTime, report},
    key = $ProfileCases[label];
    params = GetParams[key];
    Print["\n=== Profiling: ", label, " (case ", key, ") ==="];

    (* Load and prepare data *)
    data = DataG[key];
    If[params =!= {}, data = data /. params];

    (* Initialize profiling *)
    InitProfile[];
    InstallProfileWrappers[];

    (* Profile CriticalCongestionSolver *)
    Print["  Running CriticalCongestionSolver..."];
    t0 = AbsoluteTime[];
    d2e = Data2Equations[data];
    critResult = Quiet @ Check[CriticalCongestionSolver[d2e], $Failed];
    totalTime = AbsoluteTime[] - t0;
    Print["  Total time: ", NumberForm[totalTime, {6, 3}], "s"];

    report = Normal @ ProfileReport[];
    $AllProfiles[label] = report;
    AppendTo[$AllTimings, <|
      "Label" -> label,
      "Key" -> ToString[key],
      "TotalTime" -> totalTime,
      "ZAndTime" -> $ProfileTimings["ZAnd"],
      "ZAndCalls" -> $ProfileCounts["ZAnd"],
      "ZAndMaxDepth" -> $ProfileExtras["ZAnd_MaxDepth"],
      "TripleStepTime" -> $ProfileTimings["TripleStep"],
      "TripleStepCalls" -> $ProfileCounts["TripleStep"],
      "GradProjTime" -> $ProfileTimings["GradientProjection"],
      "GradProjCalls" -> $ProfileCounts["GradientProjection"],
      "Status" -> If[critResult === $Failed, "FAILED", "OK"]
    |>];

    UninstallProfileWrappers[];
  ],
  {label, Keys[$ProfileCases]}
];

(* --- Generate Markdown report --- *)

$ReportLines = {
  "# MFGraphs Bottleneck Analysis Report",
  "",
  "Generated: " <> DateString[],
  "Mathematica: " <> $Version,
  "",
  "## Summary",
  "",
  "| Label | Case | Total (s) | ZAnd (s) | ZAnd Calls | ZAnd Depth | TripleStep (s) | TripleStep Calls | GradProj (s) | Status |",
  "|---|---|---|---|---|---|---|---|---|---|"
};

Do[
  AppendTo[$ReportLines,
    "| " <> r["Label"] <> " | " <> r["Key"] <> " | " <>
    ToString[NumberForm[r["TotalTime"], {6, 3}]] <> " | " <>
    ToString[NumberForm[r["ZAndTime"], {6, 3}]] <> " | " <>
    ToString[r["ZAndCalls"]] <> " | " <>
    ToString[r["ZAndMaxDepth"]] <> " | " <>
    ToString[NumberForm[r["TripleStepTime"], {6, 3}]] <> " | " <>
    ToString[r["TripleStepCalls"]] <> " | " <>
    ToString[NumberForm[r["GradProjTime"], {6, 3}]] <> " | " <>
    r["Status"] <> " |"
  ],
  {r, $AllTimings}
];

AppendTo[$ReportLines, ""];
AppendTo[$ReportLines, "## Detailed Profiles"];
AppendTo[$ReportLines, ""];

Do[
  Module[{profile},
    AppendTo[$ReportLines, "### " <> label <> " (case " <> ToString[$ProfileCases[label]] <> ")"];
    AppendTo[$ReportLines, ""];
    AppendTo[$ReportLines, "| Function | Total Time (s) | Call Count | Avg Time (s) |"];
    AppendTo[$ReportLines, "|---|---|---|---|"];
    profile = $AllProfiles[label];
    Do[
      AppendTo[$ReportLines,
        "| " <> ToString[r["Function"]] <> " | " <>
        ToString[NumberForm[r["TotalTime"], {6, 4}]] <> " | " <>
        ToString[r["CallCount"]] <> " | " <>
        ToString[NumberForm[r["AvgTime"], {6, 6}]] <> " |"
      ],
      {r, profile}
    ];
    AppendTo[$ReportLines, ""];
  ],
  {label, Keys[$ProfileCases]}
];

(* Key observations *)
AppendTo[$ReportLines, "## Key Observations"];
AppendTo[$ReportLines, ""];

(* Find the case with the most ZAnd calls *)
Module[{maxZAnd, maxTriple},
  maxZAnd = MaximalBy[$AllTimings, #["ZAndCalls"] &];
  maxTriple = MaximalBy[$AllTimings, #["TripleStepCalls"] &];
  If[Length[maxZAnd] > 0,
    AppendTo[$ReportLines,
      "- **Highest ZAnd call count**: " <> First[maxZAnd]["Label"] <>
      " with " <> ToString[First[maxZAnd]["ZAndCalls"]] <> " calls, " <>
      "taking " <> ToString[NumberForm[First[maxZAnd]["ZAndTime"], {6, 3}]] <> "s"
    ];
  ];
  If[Length[maxTriple] > 0,
    AppendTo[$ReportLines,
      "- **Highest TripleStep call count**: " <> First[maxTriple]["Label"] <>
      " with " <> ToString[First[maxTriple]["TripleStepCalls"]] <> " calls, " <>
      "taking " <> ToString[NumberForm[First[maxTriple]["TripleStepTime"], {6, 3}]] <> "s"
    ];
  ];
];

AppendTo[$ReportLines, ""];
AppendTo[$ReportLines, "## Optimization Targets"];
AppendTo[$ReportLines, ""];
AppendTo[$ReportLines, "1. **ZAnd Solve Memoization**: Cache results of `Solve` within `ReZAnd` to avoid redundant symbolic solves across Or-branches."];
AppendTo[$ReportLines, "2. **ZAnd Branch Pruning**: Early-exit when `xp === False` to avoid processing inconsistent branches."];
AppendTo[$ReportLines, "3. **PseudoInverse Caching**: Cache the PseudoInverse result in `GradientProjection` to avoid O(n^3) recomputation per ODE step."];
AppendTo[$ReportLines, "4. **M Interpolation**: Pre-compute `M[j, x, edge]` on a grid and interpolate instead of calling `FindRoot` per point."];

(* --- Export report --- *)
$ResultsDir = FileNameJoin[{Directory[], "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];

$ReportPath = FileNameJoin[{$ResultsDir, "bottleneck_report.md"}];
Export[$ReportPath, StringRiffle[$ReportLines, "\n"], "Text"];
Print["\nBottleneck report exported to ", $ReportPath];

(* Also export raw profiling data as JSON *)
$DataPath = FileNameJoin[{$ResultsDir, "bottleneck_data.json"}];
Export[$DataPath, $AllTimings, "JSON"];
Print["Raw data exported to ", $DataPath];

Print["\nBottleneck analysis completed at ", DateString[]];
