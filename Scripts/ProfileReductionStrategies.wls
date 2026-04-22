#!/usr/bin/env wolframscript
(* ProfileReductionStrategies.wls -- Compare disjunct reduction strategies *)
(* Usage: wolframscript -file Scripts/ProfileReductionStrategies.wls <case_key> *)
(* Keys: Attraction_noswitch, Attraction_switch, Braess, Multi_20 *)

$caseKey = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "Multi_20"];

Print["=== Disjunct Reduction Strategy Profiler ==="];
Print["Case: ", $caseKey];
Print["Date: ", DateString[]];
Print[];

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

Get[FileNameJoin[{$RepoRoot, "Scripts", "BenchmarkHelpers.wls"}]];

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
      Print["    FAILED: FullReduce timed out"];
      Return[{$TimedOut, t}],
      Print["    FullReduce done: ", LeafCount[result], " leaves in ", t, "s"];
      Return[{result, t}]
    ]
  ];

(* ================================================================ *)
(* Main profiling logic                                             *)
(* ================================================================ *)

RunProfile[label_String] := Module[
  {caseSpec, caseID, params, data, d2e, xp, sorted, initRules, tag,
   branches, result1, result2, time1, time2, final1, final2},

  caseSpec = Lookup[$NamedBenchmarkCases, label, $Failed];
  If[caseSpec === $Failed,
    Print["Unknown case: ", label];
    Return[]
  ];
  {caseID, params} = caseSpec;

  Print["Loading data for ", label, " (ID: ", caseID, ")..."];
  data = GetExampleData[caseID] /. params;
  d2e = DataToEquations[data];

  Print["Extracting DNF input..."];
  {xp, sorted, initRules, tag} = ExtractDNFInput[d2e];
  If[tag === "trivial", Print["Trivial case, skipping."]; Return[]];

  branches = List @@ sorted;
  Print["Input branches: ", Length[branches]];

  Print["\n--- Running Strategy 1: Subsumption Pruning ---"];
  {time1, result1} = AbsoluteTiming[SubsumptionPrune[branches]];
  final1 = Or @@ result1;
  Print["Strategy 1 Result: ", Length[result1], " branches, ", LeafCount[final1], " leaves"];

  Print["\n--- Running Strategy 2: Full Reduce ---"];
  {result2, time2} = FullReduceStrategy[branches];
  If[result2 =!= $TimedOut,
    final2 = result2;
    Print["Strategy 2 Result: ", LeafCount[final2], " leaves"]
  ];

  Print["\n--- Summary ---"];
  Print["Strategy 1 Time: ", time1, "s"];
  Print["Strategy 2 Time: ", time2, "s"];
  If[result2 =!= $TimedOut,
    Print["LeafCount Ratio (S2/S1): ", N[LeafCount[final2]/LeafCount[final1]]];
  ];
];

RunProfile[$caseKey];
