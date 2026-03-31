#!/usr/bin/env wolframscript
(* ProfileReductionStrategies.wls -- Compare disjunct reduction strategies *)
(* Usage: wolframscript -file Scripts/ProfileReductionStrategies.wls <case_key> *)
(* Keys: Attraction_noswitch, Attraction_switch, Braess, Multi_20 *)

$caseKey = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "Multi_20"];

Print["=== Disjunct Reduction Strategy Profiler ==="];
Print["Case: ", $caseKey];
Print["Date: ", DateString[]];
Print[];

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

$ComparisonCases = <|
  "Attraction_noswitch" -> {12, {I1 -> 100, U1 -> 0}},
  "Attraction_switch" -> {11, {I1 -> 100, U1 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1}},
  "Braess" -> {"Braess congest", {}},
  "Multi_20" -> {20, {I1 -> 80, I2 -> 60, U1 -> 0, U2 -> 0, U3 -> 0}}
|>;

(* --- Extract DNF input (same as CompareDNF.wls) --- *)
ExtractDNFInput[eqs_Association] :=
  Module[{preEqs, newSystem, initRules, sorted, xp},
    preEqs = If[KeyExistsQ[eqs, "InitRules"], eqs, MFGPreprocessing[eqs]];
    initRules = Lookup[preEqs, "InitRules", <||>];
    newSystem = Lookup[preEqs, "NewSystem", {True, True, True}];
    Module[{js, approxJs, costpluscurrents, Ncpc},
      js = Lookup[preEqs, "js", {}];
      approxJs = AssociationThread[js, 0 js];
      costpluscurrents = Lookup[preEqs, "costpluscurrents", {}];
      Ncpc = Expand /@ (costpluscurrents /. approxJs);
      initRules = Expand /@ (initRules /. Ncpc);
      newSystem = newSystem /. Ncpc;
    ];
    {newSystem, initRules} = TripleClean[{newSystem, initRules}];
    If[newSystem[[3]] === True,
      {True, True, initRules, "trivial"},
      sorted = DeduplicateByComplexity[newSystem[[3]]];
      xp = And @@ Most[newSystem];
      {xp, sorted, initRules, "nontrivial"}
    ]
  ];

(* ================================================================ *)
(* Strategy 1: Subsumption pruning                                   *)
(* If branch A => branch B (A implies B), then B is redundant.       *)
(* Keep the more general branches.                                   *)
(* ================================================================ *)
SubsumptionPrune[branches_List, perCheckTimeout_:2] :=
  Module[{n = Length[branches], keep, subsumed, t, result, checks = 0, pruned = 0},
    keep = ConstantArray[True, n];
    Print["    Subsumption: ", n, " branches, up to ", n(n-1)/2, " pairwise checks (", perCheckTimeout, "s each)"];
    Do[
      If[!keep[[i]], Continue[]];
      Do[
        If[!keep[[j]], Continue[]];
        checks++;
        (* Check if branch i implies branch j (making j redundant) *)
        result = Quiet@TimeConstrained[
          Reduce[Implies[branches[[i]], branches[[j]]], Reals],
          perCheckTimeout, $TimedOut
        ];
        If[result === True,
          keep[[j]] = False; pruned++;
          Print["    Pruned branch ", j, " (subsumed by ", i, ") [", pruned, " pruned so far, ", checks, " checks]"];
          Continue[]
        ];
        (* Check if branch j implies branch i (making i redundant) *)
        result = Quiet@TimeConstrained[
          Reduce[Implies[branches[[j]], branches[[i]]], Reals],
          perCheckTimeout, $TimedOut
        ];
        If[result === True,
          keep[[i]] = False; pruned++;
          Print["    Pruned branch ", i, " (subsumed by ", j, ") [", pruned, " pruned so far, ", checks, " checks]"];
          Break[]  (* i is gone, skip to next i *)
        ],
        {j, i + 1, n}
      ],
      {i, 1, n}
    ];
    Print["    Subsumption done: ", checks, " checks, ", pruned, " pruned, ", Count[keep, True], " remaining"];
    Pick[branches, keep]
  ];

(* ================================================================ *)
(* Strategy 2: Full Reduce on the entire Or expression               *)
(* Let Mathematica's Reduce merge/simplify all branches at once.     *)
(* ================================================================ *)
FullReduceStrategy[branches_List, timeout_:300] :=
  Module[{expr, result, t},
    expr = Or @@ branches;
    Print["    FullReduce: ", Length[branches], " branches, ", LeafCount[expr], " leaves, timeout ", timeout, "s"];
    {t, result} = AbsoluteTiming[
      TimeConstrained[Quiet@Reduce[expr, Reals], timeout, $TimedOut]
    ];
    If[result === $TimedOut,
      Print["    FullReduce: TIMED OUT at ", timeout, "s"];
      {$TimedOut, t, -1},
      Module[{nDisjuncts},
        nDisjuncts = If[Head[result] === Or, Length[result], 1];
        Print["    FullReduce: ", t, "s, ", nDisjuncts, " disjuncts, ", LeafCount[result], " leaves"];
        {result, t, nDisjuncts}
      ]
    ]
  ];

(* ================================================================ *)
(* Strategy 3: Group by equality skeleton, merge inequalities        *)
(* Extract equalities from each branch, group branches sharing the   *)
(* same equalities, then try to simplify within each group.          *)
(* ================================================================ *)
ExtractEqualities[branch_And] :=
  Module[{parts = List @@ branch, eqs, ineqs},
    eqs = SortBy[Select[parts, Head[#] === Equal &], LeafCount];
    ineqs = Select[parts, Head[#] =!= Equal &];
    {eqs, ineqs}
  ];
ExtractEqualities[branch_Equal] := {{branch}, {}};
ExtractEqualities[branch_] := {{}, {branch}};

GroupAndMergeStrategy[branches_List, perGroupTimeout_:30] :=
  Module[{grouped, keys, result = {}, t, groupResult},
    Print["    GroupAndMerge: ", Length[branches], " branches"];
    grouped = GroupBy[branches, First@ExtractEqualities[#] &];
    keys = Keys[grouped];
    Print["    Found ", Length[keys], " equality groups"];
    Do[
      Module[{group, eqs, ineqSets, mergedIneqs},
        group = grouped[k];
        eqs = k;
        ineqSets = Last@ExtractEqualities[#] & /@ group;
        Print["    Group (", Length[eqs], " eqs): ", Length[group], " branches"];
        If[Length[group] === 1,
          AppendTo[result, group[[1]]],
          (* Try to simplify the Or of inequality conjunctions *)
          mergedIneqs = Or @@ (And @@@ ineqSets);
          {t, groupResult} = AbsoluteTiming[
            TimeConstrained[Quiet@Reduce[mergedIneqs, Reals], perGroupTimeout, $TimedOut]
          ];
          If[groupResult === $TimedOut,
            Print["      Merge timed out (", perGroupTimeout, "s), keeping ", Length[group], " branches"];
            result = Join[result, group],
            Module[{mergedBranches},
              mergedBranches = If[Head[groupResult] === Or,
                (And @@ Join[eqs, If[Head[#] === And, List @@ #, {#}]]) & /@ (List @@ groupResult),
                {And @@ Join[eqs, If[Head[groupResult] === And, List @@ groupResult, {groupResult}]]}
              ];
              Print["      Merged ", Length[group], " -> ", Length[mergedBranches], " in ", t, "s"];
              result = Join[result, mergedBranches]
            ]
          ]
        ]
      ],
      {k, keys}
    ];
    Print["    GroupAndMerge done: ", Length[branches], " -> ", Length[result], " branches"];
    result
  ];

(* ================================================================ *)
(* Main: load case, run DNFReduce, then try each strategy            *)
(* ================================================================ *)

$caseSpec = $ComparisonCases[$caseKey];
If[!ListQ[$caseSpec],
  Print["ERROR: Unknown case key '", $caseKey, "'"];
  Quit[1]
];

$exampleKey = $caseSpec[[1]];
$params = $caseSpec[[2]];

Print["Loading case ", $exampleKey, "..."];
$data = GetExampleData[$exampleKey];
If[$params =!= {}, $data = $data /. $params];
$d2e = DataToEquations[$data];
ClearSolveCache[];

{$xp, $sorted, $initRules, $tag} = Quiet@Check[ExtractDNFInput[$d2e], {$Failed, $Failed, <||>, "failed"}];
If[$tag === "trivial" || $tag === "failed",
  Print["  Case is ", $tag, ". Nothing to profile."];
  Quit[]
];

Print["  Input: ", LeafCount[$xp], " leaves, ", If[Head[$sorted] === Or, Length[$sorted], 1], " disjuncts"];

(* Run DNFReduce *)
Print["\n--- Running DNFReduce ---"];
ClearSolveCache[];
{$dnfTime, $dnfResult} = AbsoluteTiming[
  TimeConstrained[Quiet@Check[DNFReduce[$xp, $sorted], $Failed], 300, $TimedOut]
];
If[$dnfResult === $TimedOut || $dnfResult === $Failed,
  Print["  DNFReduce ", If[$dnfResult === $TimedOut, "TIMED OUT", "FAILED"]];
  Quit[]
];

$branches = If[Head[$dnfResult] === Or, List @@ $dnfResult, {$dnfResult}];
Print["  DNFReduce: ", $dnfTime, "s, ", Length[$branches], " disjuncts, ", LeafCount[$dnfResult], " leaves"];

If[Length[$branches] <= 1,
  Print["\n  Only 1 disjunct — nothing to reduce further."];
  Quit[]
];

(* --- Strategy 1: Subsumption --- *)
Print["\n=== Strategy 1: Subsumption Pruning ==="];
If[Length[$branches] > 500,
  Print["  Skipping: too many branches (", Length[$branches], ") for pairwise checks"];
  Print["  Estimated checks: ", Length[$branches] (Length[$branches] - 1) / 2];
  (* Try on a random sample *)
  Print["  Running on first 50 branches as sample..."];
  {$s1Time, $s1Result} = AbsoluteTiming[SubsumptionPrune[$branches[[;;50]], 1]];
  Print["  Sample: ", 50, " -> ", Length[$s1Result], " in ", $s1Time, "s"],

  {$s1Time, $s1Result} = AbsoluteTiming[SubsumptionPrune[$branches, 2]];
  Print["  Result: ", Length[$branches], " -> ", Length[$s1Result], " in ", $s1Time, "s"]
];

(* --- Strategy 2: Full Reduce --- *)
Print["\n=== Strategy 2: Full Reduce ==="];
{$s2Time, $s2Data} = AbsoluteTiming[FullReduceStrategy[$branches, 120]];

(* --- Strategy 3: Group and Merge --- *)
Print["\n=== Strategy 3: Group by Equality Skeleton ==="];
{$s3Time, $s3Result} = AbsoluteTiming[GroupAndMergeStrategy[$branches, 30]];
Print["  Result: ", Length[$branches], " -> ", Length[$s3Result], " in ", $s3Time, "s"];

(* --- Summary --- *)
Print["\n=== SUMMARY ==="];
Print["Case: ", $caseKey, " (", $exampleKey, ")"];
Print["DNFReduce output: ", Length[$branches], " disjuncts"];
Print["Strategy 1 (Subsumption):     ", If[ValueQ[$s1Result], ToString[Length[$s1Result]] <> " disjuncts, " <> ToString@Round[N[$s1Time], 0.01] <> "s", "skipped/sampled"]];
Print["Strategy 2 (Full Reduce):     ", If[$s2Data[[1]] === $TimedOut, "TIMED OUT", ToString[$s2Data[[3]]] <> " disjuncts, " <> ToString@Round[N[$s2Data[[2]]], 0.01] <> "s"]];
Print["Strategy 3 (Group & Merge):   ", ToString[Length[$s3Result]], " disjuncts, ", ToString@Round[N[$s3Time], 0.01], "s"];
