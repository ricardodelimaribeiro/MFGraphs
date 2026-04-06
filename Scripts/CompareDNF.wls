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

   --tag "..."   Run all cases and append a new entry to DNF_PERFORMANCE_HISTORY.md.
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

(* --- Utility functions --- *)

(* Count the number of disjuncts in a DNF expression *)
CountDisjuncts[expr_Or] := Length[expr];
CountDisjuncts[False] := 0;
CountDisjuncts[True] := 1;
CountDisjuncts[_] := 1;

(* Count conjuncts per disjunct *)
ConjunctsPerDisjunct[expr_Or] := Length /@ (List @@ expr) /. {n_Integer :> n, x_ :> 1};
ConjunctsPerDisjunct[expr_And] := {Length[expr]};
ConjunctsPerDisjunct[_] := {1};

(* Extract the system that DNFSolveStep would pass to DNFReduce *)
ExtractDNFInput[eqs_Association] :=
  Module[{preEqs, newSystem, initRules, sorted, xp},
    (* Run preprocessing to get the system *)
    preEqs = If[KeyExistsQ[eqs, "InitRules"],
      eqs,
      MFGPreprocessing[eqs]
    ];
    initRules = Lookup[preEqs, "InitRules", <||>];
    newSystem = Lookup[preEqs, "NewSystem", {True, True, True}];

    (* Simulate what MFGSystemSolver does before DNFSolveStep: *)
    (* Apply zero-flow (critical congestion) substitution *)
    Module[{js, approxJs, costpluscurrents, Ncpc},
      js = Lookup[preEqs, "js", {}];
      approxJs = AssociationThread[js, 0 js];
      costpluscurrents = Lookup[preEqs, "costpluscurrents", {}];
      Ncpc = Expand /@ (costpluscurrents /. approxJs);
      initRules = Expand /@ (initRules /. Ncpc);
      newSystem = newSystem /. Ncpc;
    ];

    (* TripleClean *)
    {newSystem, initRules} = TripleClean[{newSystem, initRules}];

    (* Build the arguments that DNFSolveStep passes to DNFReduce *)
    If[newSystem[[3]] === True,
      (* No disjunctions — nothing interesting to compare *)
      {True, True, initRules, "trivial"},
      sorted = DeduplicateByComplexity[newSystem[[3]]];
      xp = And @@ Most[newSystem];
      {xp, sorted, initRules, "nontrivial"}
    ]
  ];

(* --- Comparison runner --- *)

CompareMethods[key_, params_] :=
  Module[{data, d2e, xp, sorted, initRules, tag,
          bcTime, bcResult, bcMem,
          dnfTime, dnfResult, dnfMem,
          bcDisjuncts, dnfDisjuncts,
          bcLeafCount, dnfLeafCount,
          bcConjuncts, dnfConjuncts,
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

    (* --- Method 3: Reduce then BooleanConvert (the pipeline approach) --- *)
    (* This is what the code currently does after DNFReduce *)
    (* We skip this for now if either method timed out *)

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
    Print["  Speedup (BC/DNF): ",
      If[dnfTime > 0 && bcTime > 0 && bcResult =!= $TimedOut && dnfResult =!= $TimedOut,
        ToString@N[bcTime / dnfTime, 2] <> "x",
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
      "Equivalence" -> equiv,
      "EquivCheckTime" -> equivTime,
      "Speedup_BC_over_DNF" -> If[dnfTime > 0 && bcTime > 0, bcTime / dnfTime, -1]
    |>
  ];

(* --- Test cases --- *)
(* Select cases that have nontrivial Or-branches (switching costs) *)
$ComparisonCases = <|
  "Y_noswitch" -> {7, {I1 -> 50, U1 -> 0, U2 -> 0}},
  "Y_switch" -> {8, {I1 -> 100, U1 -> 0, U2 -> 0, S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1}},
  "Attraction_noswitch" -> {12, {I1 -> 100, U1 -> 0}},
  "Attraction_switch" -> {11, {I1 -> 100, U1 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1}},
  "Triangle_switch" -> {14, {I1 -> 100, U1 -> 0, U2 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1, S5 -> 1, S6 -> 1}},
  "Braess" -> {"Braess congest", {}},
  "Paper" -> {"Paper example", {I1 -> 100, I2 -> 50, U1 -> 0, U2 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1, S5 -> 1, S6 -> 1, S7 -> 1, S8 -> 1, S9 -> 1, S10 -> 1, S11 -> 1, S12 -> 1, S13 -> 1, S14 -> 1, S15 -> 1, S16 -> 1}},
  "Multi_20" -> {20, {I1 -> 80, I2 -> 60, U1 -> 0, U2 -> 0, U3 -> 0}},
  "Grid0303" -> {"Grid0303", {}},
  "Grid0505" -> {"Grid0505", {}},
  "Grid0710" -> {"Grid0710", {}},
  "Grid1010" -> {"Grid1010", {}},
  "Grid1020" -> {"Grid1020", {}}
|>;

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
  Print["Will append entry to DNF_PERFORMANCE_HISTORY.md"]
];

$results = {};

Do[
  Module[{key, params, result},
    {key, params} = $ComparisonCases[label];
    If[$requestedCase =!= "all" && label =!= $requestedCase, Continue[]];
    result = CompareMethods[key, params];
    AppendTo[$results, result];
  ],
  {label, Keys[$ComparisonCases]}
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
  "| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Speedup | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equivalence |",
  "|---|---|---|---|---|---|---|---|---|---|"
};

Do[
  If[r["Status"] === "OK",
    AppendTo[$lines,
      "| " <> r["Key"] <> " | " <>
      ToString[r["InputDisjuncts"]] <> " | " <>
      If[r["BC_Time"] >= 0, ToString@N[r["BC_Time"], 3], "TIMEOUT"] <> " | " <>
      If[r["DNF_Time"] >= 0, ToString@N[r["DNF_Time"], 3], "TIMEOUT"] <> " | " <>
      If[r["Speedup_BC_over_DNF"] > 0, ToString@N[r["Speedup_BC_over_DNF"], 2] <> "x", "N/A"] <> " | " <>
      ToString[r["BC_Disjuncts"]] <> " | " <>
      ToString[r["DNF_Disjuncts"]] <> " | " <>
      ToString[r["BC_LeafCount"]] <> " | " <>
      ToString[r["DNF_LeafCount"]] <> " | " <>
      r["Equivalence"] <> " |"
    ],
    AppendTo[$lines, "| " <> r["Key"] <> " | - | - | - | - | - | - | - | - | " <> r["Status"] <> " |"];
  ],
  {r, $results}
];

AppendTo[$lines, ""];
AppendTo[$lines, "## Interpretation"];
AppendTo[$lines, ""];
AppendTo[$lines, "- **Speedup > 1**: BooleanConvert is slower than DNFReduce (our method is faster)"];
AppendTo[$lines, "- **Speedup < 1**: BooleanConvert is faster than DNFReduce"];
AppendTo[$lines, "- **Disjuncts**: Fewer disjuncts = simpler output (fewer solution branches)"];
AppendTo[$lines, "- **Leaves**: Smaller LeafCount = more compact expression"];
AppendTo[$lines, ""];
AppendTo[$lines, "## Method details"];
AppendTo[$lines, ""];
AppendTo[$lines, "- **BooleanConvert**: Mathematica's built-in `BooleanConvert[expr, \"DNF\"]`. Treats the expression as a pure boolean formula and converts to disjunctive normal form without exploiting equation structure."];
AppendTo[$lines, "- **DNFReduce**: Our custom solver that (1) solves equalities by substitution, (2) caches Solve/Reduce results, (3) prunes inconsistent branches early, (4) short-circuits when a branch already covers the expression."];
AppendTo[$lines, ""];
AppendTo[$lines, "## Key differences"];
AppendTo[$lines, ""];
AppendTo[$lines, "DNFReduce is *not* a pure boolean conversion. It interleaves algebraic solving with DNF conversion:"];
AppendTo[$lines, "- Equalities are solved and substituted, reducing the number of free variables"];
AppendTo[$lines, "- Branches that become False after substitution are pruned immediately"];
AppendTo[$lines, "- The output is a simplified DNF where each disjunct has had equalities resolved"];
AppendTo[$lines, ""];
AppendTo[$lines, "BooleanConvert treats the entire expression as a boolean formula and does not perform algebraic simplification. This can lead to:"];
AppendTo[$lines, "- More disjuncts (redundant branches not eliminated)"];
AppendTo[$lines, "- Larger expressions (equalities not resolved)"];
AppendTo[$lines, "- Potentially faster for small problems where the boolean structure is simple"];

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

    historyPath = FileNameJoin[{Directory[], "DNF_PERFORMANCE_HISTORY.md"}];

    (* Build the table rows *)
    okResults = Select[$results, #["Status"] === "OK" &];

    tableRows = "| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |\n" <>
                "|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|\n" <>
                StringRiffle[
                  Table[
                    "| " <> r["Key"] <> " | " <>
                    ToString[r["InputLeafCount"]] <> " | " <>
                    If[r["BC_Time"] >= 0, ToString@N[r["BC_Time"], 4], "TIMEOUT"] <> " | " <>
                    If[r["DNF_Time"] >= 0, ToString@N[r["DNF_Time"], 4], "TIMEOUT"] <> " | " <>
                    If[r["Speedup_BC_over_DNF"] > 0,
                       ToString@N[r["Speedup_BC_over_DNF"], 3] <> "x", "N/A"] <> " | " <>
                    ToString[r["BC_Disjuncts"]] <> " | " <>
                    ToString[r["DNF_Disjuncts"]] <> " | " <>
                    ToString[r["BC_LeafCount"]] <> " | " <>
                    ToString[r["DNF_LeafCount"]] <> " | " <>
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
