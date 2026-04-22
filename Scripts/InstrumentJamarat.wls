#!/usr/bin/env wolframscript
(* LEGACY_GETEXAMPLEDATA: This script still depends on GetExampleData and is excluded from the active scenario-first workflow. *)
(*
  InstrumentJamarat.wls -- Instrumentation wrapper for Jamaratv9 solver run

  This script wraps the MFGraphs solver functions to collect:
  - Branch pruning statistics (immediate, elimination, survival)
  - Per-pass timing breakdown (TripleClean iterations)
  - Final statistics and LeafCount reduction

  Usage:
    wolframscript -file Scripts/InstrumentJamarat.wls

  Output:
    - Console: detailed timing and counter reports
    - Results/jamaratv9_pruning_stats.json: machine-readable statistics
*)

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Get[FileNameJoin[{$RepoRoot, "MFGraphs", "MFGraphs.wl"}]];
Get[FileNameJoin[{$RepoRoot, "Scripts", "BenchmarkHelpers.wls"}]];

(* ========== Global Statistics Accumulator ========== *)

prunStats = <|
  "timestamp" -> Now,
  "network" -> "Jamaratv9",
  "q" -> Null,
  "branchesEntered" -> 0,
  "branchesPrunedImmediate" -> 0,
  "branchesPrunedElimination" -> 0,
  "leavesSurviving" -> 0,
  "totalWallTime" -> 0.0,
  "passTimings" -> {},
  "dnfReduceStats" -> <|"branchesIn" -> 0, "branchesOut" -> 0|>,
  "reduceDisjunctsStats" -> <|"input" -> 0, "output" -> 0|>
|>;

(* ========== Wrapper: CriticalCongestionSolver ========== *)

Module[{origCriticalCongestionSolver = CriticalCongestionSolver},
  CriticalCongestionSolver[Eqs_, opts___] := Module[{t, r},
    {t, r} = AbsoluteTiming[origCriticalCongestionSolver[Eqs, opts]];
    prunStats["totalWallTime"] = t;
    (* Extract q from preprocessing if available *)
    If[KeyExistsQ[Eqs, "NewSystem"],
      Module[{newSys = Eqs["NewSystem"]},
        If[Head[newSys] === List && Length[newSys] > 2,
          Module[{orBranches = newSys[[3]]},
            If[Head[orBranches] === Or,
              prunStats["q"] = Length[List @@ orBranches],
              If[orBranches =!= True,
                prunStats["q"] = 1,
                prunStats["q"] = 0
              ]
            ]
          ]
        ]
      ]
    ];
    Print["\n>>> CriticalCongestionSolver completed"];
    Print["    Total wall time: ", N[t, 5], " seconds"];
    Print["    q (uncertain triples): ", prunStats["q"]];
    If[prunStats["q"] > 0,
      Print["    Total sign alternatives: 2^", prunStats["q"]]
    ];
    r
  ]
];

(* ========== Wrapper: TripleClean (per-pass timing) ========== *)

Module[{origTripleClean = TripleClean},
  TripleClean[args___] := Module[{t, r, passRecord},
    {t, r} = AbsoluteTiming[origTripleClean[args]];
    passRecord = <|
      "elapsed" -> N[t, 6],
      "timestamp" -> Now
    |>;
    AppendTo[prunStats["passTimings"], passRecord];
    If[Mod[Length[prunStats["passTimings"]], 5] == 0,
      Print["  TripleClean pass ", Length[prunStats["passTimings"]],
            " completed in ", N[t, 5], "s"]
    ];
    r
  ]
];

(* ========== Wrapper: DNFReduce (branch tracking) ========== *)

Module[{origDNFReduce = DNFReduce},
  DNFReduce[xp_, sys___] := Module[{r},
    (* Increment branch counter on entry *)
    prunStats["branchesEntered"]++;
    (* Call original *)
    r = origDNFReduce[xp, sys];
    (* Check for False return (immediate pruning) *)
    If[r === False,
      prunStats["branchesPrunedImmediate"]++
    ];
    r
  ]
];

(* ========== Wrapper: ReduceDisjuncts (post-reduction tracking) ========== *)

Module[{origReduceDisjuncts = ReduceDisjuncts},
  ReduceDisjuncts[expr_] := Module[{r, inputCount, outputCount},
    (* Count input branches *)
    inputCount = If[Head[expr] === Or, Length[List @@ expr], 1];
    prunStats["reduceDisjunctsStats"]["input"] = inputCount;
    (* Call original *)
    r = origReduceDisjuncts[expr];
    (* Count output branches *)
    outputCount = If[Head[r] === Or, Length[List @@ r], If[r === False, 0, 1]];
    prunStats["reduceDisjunctsStats"]["output"] = outputCount;
    (* Track elimination *)
    prunStats["branchesPrunedElimination"] = inputCount - outputCount;
    prunStats["leavesSurviving"] = outputCount;
    Print["\n>>> ReduceDisjuncts: ", inputCount, " -> ", outputCount,
          " branches (pruned: ", inputCount - outputCount, ")"];
    r
  ]
];

(* ========== Main Execution ========== *)

Print["\n", Style["========== JAMARATV9 PRUNING STATISTICS RUN ==========", Bold]];
Print["Timestamp: ", prunStats["timestamp"]];
Print["Network: ", prunStats["network"]];
Print[];

(* Load data with canonical parameters *)
data = GetExampleData["Jamaratv9"];
params = GetParams["Jamaratv9"];
Print["Parameters: I1 -> 100, I2 -> 50, U1 -> 0, U2 -> 0, U3 -> 0"];

data = data /. {I1 -> 100, I2 -> 50, U1 -> 0, U2 -> 0, U3 -> 0};
Print["Entrances: ", data["Entrance Vertices and Flows"]];
Print["Exits: ", data["Exit Vertices and Terminal Costs"]];
Print["Switching costs: (none)"];
Print[];

(* Enable verbose mode *)
$MFGraphsVerbose = True;

(* Run the solver *)
Print[Style["Running solver...", Italic]];
Print[];
d2e = DataToEquations[data];
result = CriticalCongestionSolver[d2e];

(* Print final statistics *)
Print[];
Print[Style["========== FINAL STATISTICS ==========", Bold]];
Print["Total wall time: ", N[prunStats["totalWallTime"], 5], " seconds"];
Print["q (uncertain triples): ", prunStats["q"]];
If[prunStats["q"] > 0,
  Print["Total sign alternatives: 2^", prunStats["q"]]
];
Print[];
Print["Branch tracking:"];
Print["  Branches entered: ", prunStats["branchesEntered"]];
Print["  Pruned (immediate): ", prunStats["branchesPrunedImmediate"]];
Print["  Pruned (elimination): ", prunStats["branchesPrunedElimination"]];
Print["  Surviving leaves: ", prunStats["leavesSurviving"]];
Print[];
Print["TripleClean passes: ", Length[prunStats["passTimings"]]];
If[Length[prunStats["passTimings"]] > 0,
  Print["  Total pass time: ",
    N[Total[Map[#["elapsed"] &, prunStats["passTimings"]]], 5], " seconds"]
];
Print[];

(* Check feasibility *)
If[IsFeasible[result],
  Print[Style["✓ FEASIBLE SOLUTION FOUND", Green]],
  Print[Style["✗ INFEASIBLE (no solution)", Red]]
];
Print[];

(* Export results *)
exportPath = FileNameJoin[{$RepoRoot, "Results", "jamaratv9_pruning_stats.json"}];
If[!DirectoryQ[FileNameJoin[{$RepoRoot, "Results"}]],
  CreateDirectory[FileNameJoin[{$RepoRoot, "Results"}]]
];
Export[exportPath, prunStats, "JSON"];
Print["Results exported to: ", exportPath];
Print[];

(* Return the statistics *)
prunStats
