#!/usr/bin/env wolframscript
(* CompareDNF.wls -- Compare BooleanConvert["DNF"] vs DNFReduce *)
(*
   Usage:
     wolframscript -file Scripts/CompareDNF.wls
     wolframscript -file Scripts/CompareDNF.wls <case_key>
     wolframscript -file Scripts/CompareDNF.wls --tag "Description of change" \
       --rationale "Why this change exists" \
       --changes "What code changed" \
       --interpretation "How to read the results"
     wolframscript -file Scripts/CompareDNF.wls <case_key> --tag "Description" \
       rationale="..." changes="..." interpretation="..."

   --tag "..."   Run all cases and append a new entry to docs/history/DNF_PERFORMANCE_HISTORY.md.
                 The entry includes the current git commit hash, date, results table,
                 and the supplied rationale / changes / interpretation prose.
*)
Print["=== BooleanConvert vs DNFReduce Comparison ==="];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print[];

(* --- Setup --- *)
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

(* --- Comparison runner --- *)

CompareMethods[key_, params_] :=
  Module[{data, d2e, xp, sorted, initRules, tag,
          bcTime, bcResult, bcMem,
          dnfTime, dnfResult, dnfMem,
          pipeTime, pipeResult, pipeMem,
          bcDisjuncts, dnfDisjuncts, pipeDisjuncts,
          bcLeafCount, dnfLeafCount, pipeLeafCount,
          bcConjuncts, dnfConjuncts, pipeConjuncts,
          equivTime, equiv},

    Print["--- Case: ", key, " ---"];

    (* Load data and build equations *)
    data = GetExampleData[key];
    If[params =!= {}, data = data /. params];
    d2e = DataToEquations[data];
    ClearSolveCache[];

    (* Extract the input that would go to DNFReduce *)
    {xp, sorted, initRules, tag} = Quiet @ Check[
      ExtractDNFInput[d2e],
      {$Failed, $Failed, <||>, "failed"}
    ];

    If[tag === "failed",
      Print["  FAILED to extract DNF input."];
      Return[<|"Key" -> ToString[key], "Status" -> "FAILED"|>]
    ];

    If[tag === "trivial",
      Print["  System is trivial (no disjunctions). Skipping comparison."];
      Return[<|"Key" -> ToString[key], "Status" -> "TRIVIAL"|>]
    ];

    Print["  Input: ", LeafCount[xp], " leaves in xp, ",
          CountDisjuncts[sorted], " disjuncts in sorted"];

    (* --- Method 1: BooleanConvert[expr, "DNF"] --- *)
    Print["  Running BooleanConvert[..., \"DNF\"]..."];
    bcMem = MemoryInUse[];
    {bcTime, bcResult} = AbsoluteTiming[
      TimeConstrained[
        Quiet @ Check[BooleanConvert[xp && sorted, "DNF"], $Failed],
        300, $TimedOut
      ]
    ];
    bcMem = MaxMemoryUsed[] - bcMem;

    If[bcResult === $TimedOut,
      Print["  BooleanConvert TIMED OUT (300s)"];
      bcDisjuncts = -1; bcLeafCount = -1; bcConjuncts = {};,
      If[bcResult === $Failed,
        Print["  BooleanConvert FAILED"];
        bcDisjuncts = -1; bcLeafCount = -1; bcConjuncts = {};,
        bcDisjuncts = CountDisjuncts[bcResult];
        bcLeafCount = LeafCount[bcResult];
        bcConjuncts = ConjunctsPerDisjunct[bcResult];
        Print["  BooleanConvert: ", ToString@N[bcTime, 6], "s, ",
              bcDisjuncts, " disjuncts, ", bcLeafCount, " leaves"];
      ]
    ];

    (* --- Method 2: DNFReduce --- *)
    Print["  Running DNFReduce..."];
    ClearSolveCache[];
    dnfMem = MemoryInUse[];
    {dnfTime, dnfResult} = AbsoluteTiming[
      TimeConstrained[
        Quiet @ Check[DNFReduce[xp, sorted], $Failed],
        300, $TimedOut
      ]
    ];
    dnfMem = MaxMemoryUsed[] - dnfMem;

    If[dnfResult === $TimedOut,
      Print["  DNFReduce TIMED OUT (300s)"];
      dnfDisjuncts = -1; dnfLeafCount = -1; dnfConjuncts = {};,
      If[dnfResult === $Failed,
        Print["  DNFReduce FAILED"];
        dnfDisjuncts = -1; dnfLeafCount = -1; dnfConjuncts = {};,
        dnfDisjuncts = CountDisjuncts[dnfResult];
        dnfLeafCount = LeafCount[dnfResult];
        dnfConjuncts = ConjunctsPerDisjunct[dnfResult];
        Print["  DNFReduce:      ", ToString@N[dnfTime, 6], "s, ",
              dnfDisjuncts, " disjuncts, ", dnfLeafCount, " leaves"];
      ]
    ];

    (* --- Method 3: BooleanConvert + Reduce (the pipeline approach) --- *)
    Print["  Running Reduce on BooleanConvert output..."];
    If[bcResult === $TimedOut || bcResult === $Failed,
      Print["  Skipping Method 3 because BooleanConvert failed."];
      pipeTime = 0; pipeResult = $Failed; pipeMem = 0;
      pipeDisjuncts = -1; pipeLeafCount = -1; pipeConjuncts = {};,

      ClearSystemCache[];
      pipeMem = MemoryInUse[];
      {pipeTime, pipeResult} = AbsoluteTiming[
        TimeConstrained[
          Quiet @ Check[Reduce[bcResult, Reals], $Failed],
          300, $TimedOut
        ]
      ];
      pipeMem = MaxMemoryUsed[] - pipeMem;

      If[pipeResult === $TimedOut,
        Print["  Method 3 TIMED OUT (300s)"];
        pipeDisjuncts = -1; pipeLeafCount = -1; pipeConjuncts = {};,
        If[pipeResult === $Failed,
          Print["  Method 3 FAILED"];
          pipeDisjuncts = -1; pipeLeafCount = -1; pipeConjuncts = {};,
          pipeDisjuncts = CountDisjuncts[pipeResult];
          pipeLeafCount = LeafCount[pipeResult];
          pipeConjuncts = ConjunctsPerDisjunct[pipeResult];
          Print["  Method 3:       ", ToString@N[pipeTime, 6], "s, ",
                pipeDisjuncts, " disjuncts, ", pipeLeafCount, " leaves"];
        ]
      ];
    ];

    (* --- Equivalence check --- *)
    equiv = "SKIPPED";
    equivTime = 0;
    If[bcResult =!= $TimedOut && bcResult =!= $Failed &&
       dnfResult =!= $TimedOut && dnfResult =!= $Failed,
      Print["  Checking equivalence..."];
      {equivTime, equiv} = AbsoluteTiming[
        TimeConstrained[
          Quiet @ Check[
            Reduce[bcResult \[Equivalent] dnfResult, Reals],
            "ERROR"
          ],
          60, "EQUIV_TIMEOUT"
        ]
      ];
      equiv = Which[
        equiv === True, "EQUIVALENT",
        equiv === "EQUIV_TIMEOUT", "EQUIV_TIMEOUT",
        equiv === "ERROR", "EQUIV_ERROR",
        True, "NOT_EQUIVALENT: " <> ToString[InputForm[equiv]]
      ];
      Print["  Equivalence: ", equiv, " (", ToString@N[equivTime, 6], "s)"];
    ];

    (* --- Summary --- *)
    Print["  Speedup (BC/DNF):   ",
      If[dnfTime > 0 && bcTime > 0 && bcResult =!= $TimedOut && dnfResult =!= $TimedOut,
        ToString@N[bcTime / dnfTime, 2] <> "x",
        "N/A"
      ]
    ];
    Print["  Speedup (Pipe/DNF): ",
      If[dnfTime > 0 && pipeTime > 0 && pipeResult =!= $TimedOut && dnfResult =!= $TimedOut,
        ToString@N[(bcTime + pipeTime) / dnfTime, 2] <> "x",
        "N/A"
      ]
    ];
    Print[];

    <|
      "Key" -> ToString[key],
      "Status" -> "OK",
      "InputLeafCount" -> LeafCount[xp && sorted],
      "InputDisjuncts" -> CountDisjuncts[sorted],
      "BC_Time" -> bcTime,
      "BC_Disjuncts" -> bcDisjuncts,
      "BC_LeafCount" -> bcLeafCount,
      "BC_AvgConjuncts" -> If[Length[bcConjuncts] > 0, N@Mean[bcConjuncts], -1],
      "BC_Memory" -> bcMem,
      "DNF_Time" -> dnfTime,
      "DNF_Disjuncts" -> dnfDisjuncts,
      "DNF_LeafCount" -> dnfLeafCount,
      "DNF_AvgConjuncts" -> If[Length[dnfConjuncts] > 0, N@Mean[dnfConjuncts], -1],
      "DNF_Memory" -> dnfMem,
      "Pipe_Time" -> pipeTime,
      "Pipe_Disjuncts" -> pipeDisjuncts,
      "Pipe_LeafCount" -> pipeLeafCount,
      "Pipe_AvgConjuncts" -> If[Length[pipeConjuncts] > 0, N@Mean[pipeConjuncts], -1],
      "Pipe_Memory" -> pipeMem,
      "Equivalence" -> equiv,
      "EquivCheckTime" -> equivTime,
      "Speedup_BC_over_DNF" -> If[dnfTime > 0 && bcTime > 0, bcTime / dnfTime, -1],
      "Speedup_Pipe_over_DNF" -> If[dnfTime > 0 && bcTime > 0 && pipeTime >= 0, (bcTime + pipeTime) / dnfTime, -1]
    |>
  ];

(* --- CLI argument parsing --- *)

Module[{args = Rest[$ScriptCommandLine], tagPos},
  $historyTag = None;
  $historyRationale = None;
  $historyChanges = None;
  $historyInterpretation = None;
  $requestedCase = "all";
  While[Length[args] > 0,
    Switch[First[args],
      "--tag",
        If[Length[args] >= 2, $historyTag = args[[2]]];
        args = Drop[args, 2],
      "--rationale",
        If[Length[args] >= 2, $historyRationale = args[[2]]];
        args = Drop[args, 2],
      "--changes",
        If[Length[args] >= 2, $historyChanges = args[[2]]];
        args = Drop[args, 2],
      "--interpretation",
        If[Length[args] >= 2, $historyInterpretation = args[[2]]];
        args = Drop[args, 2],
      _ /; StringStartsQ[First[args], "tag="],
        $historyTag = StringDrop[First[args], StringLength["tag="]];
        args = Rest[args],
      _ /; StringStartsQ[First[args], "rationale="],
        $historyRationale = StringDrop[First[args], StringLength["rationale="]];
        args = Rest[args],
      _ /; StringStartsQ[First[args], "changes="],
        $historyChanges = StringDrop[First[args], StringLength["changes="]];
        args = Rest[args],
      _ /; StringStartsQ[First[args], "interpretation="],
        $historyInterpretation = StringDrop[First[args], StringLength["interpretation="]];
        args = Rest[args],
      _,
        If[$requestedCase === "all", $requestedCase = First[args]];
        args = Rest[args]
    ]
  ];
];

If[$historyTag =!= None,
  Print["History tag: ", $historyTag];
  If[MemberQ[{$historyRationale, $historyChanges, $historyInterpretation}, None],
    Print["History entries require explicit rationale, changes, and interpretation text."];
    Exit[2]
  ];
  Print["Will append entry to docs/history/DNF_PERFORMANCE_HISTORY.md"]
];

$results = {};

Do[
  Module[{key, params, result},
    {key, params} = $NamedBenchmarkCases[label];
    If[$requestedCase =!= "all" && label =!= $requestedCase, Continue[]];
    result = CompareMethods[key, params];
    AppendTo[$results, result];
  ],
  {label, Keys[$NamedBenchmarkCases]}
];

(* --- Generate Markdown report --- *)

$lines = {
  "# BooleanConvert vs DNFReduce Comparison",
  "",
  "Generated: " <> DateString[],
  "Mathematica: " <> $Version,
  "",
  "## Summary",
  "",
  "| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Pipe Time (s) | **Speedup (Pipe)** | BC Disjuncts | DNF Disjuncts | Pipe Disjuncts | BC Leaves | DNF Leaves | Pipe Leaves | Equiv |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|---|"
};

Do[
  If[r["Status"] === "OK",
    AppendTo[$lines,
      "| " <> r["Key"] <> " | " <>
      ToString[r["InputDisjuncts"]] <> " | " <>
      If[r["BC_Time"] >= 0, ToString@N[r["BC_Time"], 3], "TIMEOUT"] <> " | " <>
      If[r["DNF_Time"] >= 0, ToString@N[r["DNF_Time"], 3], "TIMEOUT"] <> " | " <>
      If[r["Pipe_Time"] >= 0, ToString@N[r["Pipe_Time"], 3], "TIMEOUT"] <> " | " <>
      If[r["Speedup_Pipe_over_DNF"] > 0, ToString@N[r["Speedup_Pipe_over_DNF"], 2] <> "x", "N/A"] <> " | " <>

      ToString[r["BC_Disjuncts"]] <> " | " <>
      ToString[r["DNF_Disjuncts"]] <> " | " <>
      ToString[r["Pipe_Disjuncts"]] <> " | " <>
      ToString[r["BC_LeafCount"]] <> " | " <>
      ToString[r["DNF_LeafCount"]] <> " | " <>
      ToString[r["Pipe_LeafCount"]] <> " | " <>
      Which[r["Equivalence"] === "EQUIVALENT", "\[Checkmark]",
            r["Equivalence"] === "EQUIV_TIMEOUT", "\[Dagger]",
            True, r["Equivalence"]] <> " |"
    ],
    AppendTo[$lines, "| " <> r["Key"] <> " | - | - | - | - | - | - | - | - | - | - | - | " <> r["Status"] <> " |"];
  ],
  {r, $results}
];

AppendTo[$lines, ""];
AppendTo[$lines, "## Interpretation"];
AppendTo[$lines, ""];
AppendTo[$lines, "- **Speedup (Pipe)**: How much faster `DNFReduce` is compared to the `BooleanConvert` + `Reduce` pipeline."];
AppendTo[$lines, "- **BC Time**: Time for pure boolean conversion."];
AppendTo[$lines, "- **Pipe Time**: Time for boolean conversion followed by algebraic simplification (`Reduce`). This is the fair baseline for `DNFReduce`."];
AppendTo[$lines, "- **Disjuncts**: Fewer disjuncts = simpler output (fewer solution branches)"];
AppendTo[$lines, "- **Leaves**: Smaller LeafCount = more compact expression"];
AppendTo[$lines, ""];
AppendTo[$lines, "## Method details"];
AppendTo[$lines, ""];
AppendTo[$lines, "- **BooleanConvert**: Mathematica's built-in `BooleanConvert[expr, \"DNF\"]`. Treats the expression as a pure boolean formula."];
AppendTo[$lines, "- **DNFReduce**: Our custom solver that interleaves algebraic solving (substitution and pruning) with DNF conversion."];
AppendTo[$lines, "- **Pipe (Method 3)**: `BooleanConvert` followed by `Reduce[..., Reals]`. This isolates whether the benefit of `DNFReduce` comes from the simplification itself or the interleaved strategy."];
AppendTo[$lines, ""];
AppendTo[$lines, "## Key differences"];
AppendTo[$lines, ""];
AppendTo[$lines, "DNFReduce is *not* a pure boolean conversion. It interleaves algebraic solving with DNF conversion:"];
AppendTo[$lines, "- Equalities are solved and substituted, reducing the number of free variables"];
AppendTo[$lines, "- Branches that become False after substitution are pruned immediately"];
AppendTo[$lines, "- The output is a simplified DNF where each disjunct has had equalities resolved"];
AppendTo[$lines, ""];
AppendTo[$lines, "BooleanConvert treats the entire expression as a boolean formula and does not perform algebraic simplification. Method 3 (Pipe) adds that simplification as a post-processing step."];

(* --- Export --- *)
$ResultsDir = FileNameJoin[{Directory[], "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];

$reportPath = FileNameJoin[{$ResultsDir, "compare_dnf.md"}];
Export[$reportPath, StringRiffle[$lines, "\n"], "Text"];
Print["Report exported to ", $reportPath];

$dataPath = FileNameJoin[{$ResultsDir, "compare_dnf.json"}];
Export[$dataPath, $results, "JSON"];
Print["Data exported to ", $dataPath];

(* --- Append to performance history if --tag was given --- *)

If[$historyTag =!= None,
  Module[{commitHash, historyPath, entry, tableRows, okResults},

    (* Resolve current git commit hash via RunProcess *)
    commitHash = Quiet @ Check[
      Module[{proc = RunProcess[{"git", "rev-parse", "--short", "HEAD"}]},
        If[proc["ExitCode"] === 0, StringTrim[proc["StandardOutput"]], "unknown"]
      ],
      "unknown"
    ];

    historyPath = FileNameJoin[{Directory[], "docs", "history", "DNF_PERFORMANCE_HISTORY.md"}];
    (* Fallback to root if not found in docs/history *)
    If[!FileExistsQ[historyPath] && FileExistsQ[FileNameJoin[{Directory[], "DNF_PERFORMANCE_HISTORY.md"}]],
      historyPath = FileNameJoin[{Directory[], "DNF_PERFORMANCE_HISTORY.md"}]
    ];

    (* Build the table rows *)
    okResults = Select[$results, #["Status"] === "OK" &];

    tableRows = "| Case | In Lv | BC (s) | DNF (s) | Pipe (s) | **Speedup** | BC Dj | DNF Dj | Pipe Dj | BC Lv | DNF Lv | Pipe Lv | Eq |\n" <>
                "|------|-------|--------|---------|----------|-------------|-------|--------|---------|-------|--------|---------|----|\n" <>
                StringRiffle[
                  Table[
                    "| " <> r["Key"] <> " | " <>
                    ToString[r["InputLeafCount"]] <> " | " <>
                    If[r["BC_Time"] >= 0, ToString@N[r["BC_Time"], 3], "T/O"] <> " | " <>
                    If[r["DNF_Time"] >= 0, ToString@N[r["DNF_Time"], 3], "T/O"] <> " | " <>
                    If[r["Pipe_Time"] >= 0, ToString@N[r["Pipe_Time"], 3], "T/O"] <> " | " <>
                    If[r["Speedup_Pipe_over_DNF"] > 0,
                       ToString@N[r["Speedup_Pipe_over_DNF"], 2] <> "x", "N/A"] <> " | " <>
                    ToString[r["BC_Disjuncts"]] <> " | " <>
                    ToString[r["DNF_Disjuncts"]] <> " | " <>
                    ToString[r["Pipe_Disjuncts"]] <> " | " <>
                    ToString[r["BC_LeafCount"]] <> " | " <>
                    ToString[r["DNF_LeafCount"]] <> " | " <>
                    ToString[r["Pipe_LeafCount"]] <> " | " <>
                    Which[r["Equivalence"] === "EQUIVALENT", "\[Checkmark]",
                          r["Equivalence"] === "EQUIV_TIMEOUT", "~\[Dagger]",
                          r["Equivalence"] === "EQUIV_ERROR", "ERR",
                          True, r["Equivalence"]] <> " |",
                    {r, okResults}
                  ],
                  "\n"
                ];

    entry = StringJoin[
      "\n---\n\n",
      "### ", DateString[{"Year", "-", "Month", "-", "Day"}], " \[LongDash] ", $historyTag, "\n\n",
      "**Commit:** `", commitHash, "`  \n",
      "**Date:** ", DateString[], "  \n",
      "**Generated by:** `wolframscript -file Scripts/CompareDNF.wls --tag \"", $historyTag, "\"`\n\n",
      "#### Measured performance\n\n",
      tableRows, "\n\n",
      "\[Dagger] *Equivalence check timed out (60 s); outputs are expected to be equivalent.*\n\n",
      "#### Rationale\n\n",
      $historyRationale, "\n\n",
      "#### Changes\n\n",
      $historyChanges, "\n\n",
      "#### Interpretation\n\n",
      $historyInterpretation, "\n\n"
    ];

    (* Append to history file *)
    If[FileExistsQ[historyPath],
      (* Append before the template section if present, otherwise at end *)
      Module[{existing, templateIdx, newContent},
        existing = Import[historyPath, "Text"];
        templateIdx = StringPosition[existing, "## Future entries (template)"];
        If[Length[templateIdx] > 0,
          newContent = StringInsert[existing, entry, First[First[templateIdx]]],
          newContent = existing <> entry
        ];
        Export[historyPath, newContent, "Text"]
      ],
      Export[historyPath, entry, "Text"]
    ];
    Print["History entry appended to ", historyPath];
  ]
];

Print["\nComparison completed at ", DateString[]];
