#!/usr/bin/env wolframscript
(* LEGACY_GETEXAMPLEDATA: This script still depends on GetExampleData and is excluded from the active scenario-first workflow. *)
(* ProfilePreprocessing.wls -- Profile each preprocessing step and test optimization strategies *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
SetDirectory[$RepoRoot];

Needs["MFGraphs`"];

Get[FileNameJoin[{$RepoRoot, "Scripts", "BenchmarkHelpers.wls"}]];

$MFGraphsVerbose = False;

(* ============================================================ *)
(* Utility: run expr N times and return {mean, stddev}          *)
(* ============================================================ *)

TimeMean[expr_, n_:3] := Module[{times},
  times = Table[First@AbsoluteTiming[expr], {n}];
  {Mean[times], StandardDeviation[times]}
];

(* ============================================================ *)
(* Test cases: 8 (6 switching costs), 11 (16 switching costs),  *)
(* "Paper example" (4 switching costs)                          *)
(* ============================================================ *)

$TestCases = {8, 19, "Paper example"};

(* ============================================================ *)
(* PHASE 1: Baseline Profiling                                  *)
(* Manually replicate MFGPreprocessing steps with timing        *)
(* ============================================================ *)

ProfilePreprocessingBaseline[caseKey_] := Module[
  {data, d2e, Eqs, InitRules, RuleBalanceGatheringFlows, EqEntryIn,
   RuleEntryOut, RuleExitFlowsIn, RuleExitValues, EqValueAuxiliaryEdges,
   IneqSwitchingByVertex, AltOptCond, EqBalanceSplittingFlows, AltFlows,
   AltTransitionFlows, IneqJs, IneqJts, NewSystem, Rules, EqGeneral,
   RuleEntryValues, timings, items, ineqTriple,
   time1, time2, time3, time4, time5, time6, time7},

  timings = <||>;

  (* Step 0: DataToEquations *)
  data = GetExampleData[caseKey] /. GetParams[caseKey];
  {time1, d2e} = AbsoluteTiming[DataToEquations[data]];
  timings["DataToEquations"] = time1;
  If[d2e === $Failed, Print["  DataToEquations returned $Failed — skipping case."]; Return[timings, Module]];
  Eqs = d2e;

  (* Extract fields *)
  ClearSolveCache[];
  RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows", $Failed];
  EqEntryIn = Lookup[Eqs, "EqEntryIn", $Failed];
  RuleEntryOut = Lookup[Eqs, "RuleEntryOut", $Failed];
  RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", $Failed];
  RuleExitValues = Lookup[Eqs, "RuleExitValues", $Failed];
  RuleEntryValues = Lookup[Eqs, "RuleEntryValues", $Failed];
  AltOptCond = Lookup[Eqs, "AltOptCond", $Failed];
  AltFlows = Lookup[Eqs, "AltFlows", $Failed];
  AltTransitionFlows = Lookup[Eqs, "AltTransitionFlows", $Failed];
  IneqJs = Lookup[Eqs, "IneqJs", $Failed];
  IneqJts = Lookup[Eqs, "IneqJts", $Failed];
  EqGeneral = Lookup[Eqs, "EqGeneral", $Failed];
  IneqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex", $Failed];
  EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", $Failed];
  EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];

  (* Step 1: Build InitRules *)
  {time1, InitRules} = AbsoluteTiming[
    Module[{ir},
      ir = Association[Flatten[ToRules /@ EqEntryIn]];
      AssociateTo[ir, RuleEntryOut];
      AssociateTo[ir, RuleExitFlowsIn];
      AssociateTo[ir, RuleEntryValues];
      AssociateTo[ir, RuleBalanceGatheringFlows];
      AssociateTo[ir, RuleExitValues];
      ir
    ]
  ];
  timings["BuildInitRules"] = time1;

  (* Step 2: Simplify switching cost inequalities (BASELINE: per-element) *)
  {time2, IneqSwitchingByVertex} = AbsoluteTiming[
    With[{it = IneqSwitchingByVertex /. InitRules},
      And @@ MFGParallelMap[Simplify, it]
    ]
  ];
  timings["SimplifySwitchingCosts_PerElement"] = time2;

  (* Step 3: Substitute + Solve balance equations *)
  EqBalanceSplittingFlows = EqBalanceSplittingFlows /. InitRules;
  EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
  {time3, Rules} = AbsoluteTiming[First @ Solve[EqBalanceSplittingFlows && EqValueAuxiliaryEdges] // Quiet];
  InitRules = Join[InitRules /. Rules, Association @ Rules];
  timings["SolveBalanceEquations"] = time3;

  (* Step 4: SystemToTriple *)
  {time4, ineqTriple} = AbsoluteTiming[SystemToTriple[IneqSwitchingByVertex]];
  NewSystem = {ineqTriple[[1]], IneqJts && IneqJs && ineqTriple[[2]],
               AltFlows && AltTransitionFlows && AltOptCond && ineqTriple[[3]]};
  timings["SystemToTriple"] = time4;

  (* Step 5: TripleClean #1 *)
  {time5, {NewSystem, InitRules}} = AbsoluteTiming[TripleClean[{NewSystem, InitRules}]];
  timings["TripleClean_1"] = time5;

  (* Step 6: Add EqGeneral + TripleClean #2 *)
  NewSystem[[1]] = EqGeneral;
  {time6, {NewSystem, InitRules}} = AbsoluteTiming[TripleClean[{NewSystem, InitRules}]];
  timings["TripleClean_2_WithEqGeneral"] = time6;

  (* Step 7: MFGSystemSolver (critical congestion: all flows = 0) *)
  Module[{PreEqs, js, AssoCritical},
    PreEqs = Join[Eqs, <|"InitRules" -> InitRules, "NewSystem" -> NewSystem|>];
    js = Lookup[PreEqs, "js", $Failed];
    ClearSolveCache[];
    {time7, AssoCritical} = AbsoluteTiming[MFGSystemSolver[PreEqs][AssociationThread[js, 0 js]]];
    timings["MFGSystemSolver"] = time7;
  ];

  timings["TOTAL_Preprocessing"] = timings["BuildInitRules"] + timings["SimplifySwitchingCosts_PerElement"] +
    timings["SolveBalanceEquations"] + timings["SystemToTriple"] +
    timings["TripleClean_1"] + timings["TripleClean_2_WithEqGeneral"];

  timings
];

(* ============================================================ *)
(* PHASE 2: Test optimization strategies                        *)
(* ============================================================ *)

(* Strategy A: Simplify whole system instead of per-element *)
ProfileStrategyA[caseKey_] := Module[
  {data, Eqs, InitRules, IneqSwitchingByVertex, timeBaseline, timeOptimized,
   resultBaseline, resultOptimized, items},

  data = GetExampleData[caseKey] /. GetParams[caseKey];
  Eqs = DataToEquations[data];

  (* Build InitRules *)
  InitRules = Association[Flatten[ToRules /@ Lookup[Eqs, "EqEntryIn"]]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleEntryOut"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleExitFlowsIn"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleEntryValues"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleBalanceGatheringFlows"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleExitValues"]];

  IneqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex"];
  items = IneqSwitchingByVertex /. InitRules;

  (* Baseline: per-element Simplify *)
  {timeBaseline, resultBaseline} = AbsoluteTiming[And @@ MFGParallelMap[Simplify, items]];

  (* Strategy A: whole-system Simplify *)
  {timeOptimized, resultOptimized} = AbsoluteTiming[Simplify[And @@ items]];

  <|
    "Baseline_Time" -> timeBaseline,
    "StrategyA_Time" -> timeOptimized,
    "Speedup" -> timeBaseline / Max[timeOptimized, 0.0001],
    "Equivalent" -> TimeConstrained[Reduce[Equivalent[resultBaseline, resultOptimized], Reals] === True, 30, "Timeout"]
  |>
];

(* Strategy B: Single TripleClean with EqGeneral from the start *)
ProfileStrategyB[caseKey_] := Module[
  {data, Eqs, InitRules, IneqSwitchingByVertex, EqBalanceSplittingFlows,
   EqValueAuxiliaryEdges, AltOptCond, AltFlows, AltTransitionFlows,
   IneqJs, IneqJts, EqGeneral, Rules, ineqTriple,
   NewSystem1, NewSystem2, InitRules1, InitRules2,
   timeBaseline, timeOptimized, items},

  data = GetExampleData[caseKey] /. GetParams[caseKey];
  Eqs = DataToEquations[data];

  (* Build InitRules *)
  InitRules = Association[Flatten[ToRules /@ Lookup[Eqs, "EqEntryIn"]]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleEntryOut"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleExitFlowsIn"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleEntryValues"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleBalanceGatheringFlows"]];
  AssociateTo[InitRules, Lookup[Eqs, "RuleExitValues"]];

  IneqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex"];
  items = IneqSwitchingByVertex /. InitRules;
  IneqSwitchingByVertex = And @@ (Simplify /@ items);

  EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows"] /. InitRules;
  EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges"] /. InitRules;
  Rules = First @ Solve[EqBalanceSplittingFlows && EqValueAuxiliaryEdges] // Quiet;
  InitRules = Join[InitRules /. Rules, Association @ Rules];

  AltOptCond = Lookup[Eqs, "AltOptCond"];
  AltFlows = Lookup[Eqs, "AltFlows"];
  AltTransitionFlows = Lookup[Eqs, "AltTransitionFlows"];
  IneqJs = Lookup[Eqs, "IneqJs"];
  IneqJts = Lookup[Eqs, "IneqJts"];
  EqGeneral = Lookup[Eqs, "EqGeneral"];

  ineqTriple = SystemToTriple[IneqSwitchingByVertex];

  (* BASELINE: Two TripleClean passes *)
  ClearSolveCache[];
  InitRules1 = InitRules;
  NewSystem1 = {ineqTriple[[1]], IneqJts && IneqJs && ineqTriple[[2]],
                AltFlows && AltTransitionFlows && AltOptCond && ineqTriple[[3]]};
  {timeBaseline, {NewSystem1, InitRules1}} = AbsoluteTiming[
    Module[{ns, ir},
      {ns, ir} = TripleClean[{NewSystem1, InitRules1}];
      ns[[1]] = EqGeneral;
      TripleClean[{ns, ir}]
    ]
  ];

  (* Strategy B: Single TripleClean with EqGeneral included from start *)
  ClearSolveCache[];
  InitRules2 = InitRules;
  NewSystem2 = {ineqTriple[[1]] && EqGeneral, IneqJts && IneqJs && ineqTriple[[2]],
                AltFlows && AltTransitionFlows && AltOptCond && ineqTriple[[3]]};
  {timeOptimized, {NewSystem2, InitRules2}} = AbsoluteTiming[
    TripleClean[{NewSystem2, InitRules2}]
  ];

  <|
    "Baseline_Time" -> timeBaseline,
    "StrategyB_Time" -> timeOptimized,
    "Speedup" -> timeBaseline / Max[timeOptimized, 0.0001],
    "BaselineRuleCount" -> Length[InitRules1],
    "OptimizedRuleCount" -> Length[InitRules2],
    "RulesEquivalent" -> (Sort[Keys[InitRules1]] === Sort[Keys[InitRules2]])
  |>
];

(* Strategy C: Index-based inequality grouping in MFGSystemSolver *)
(* We test grouping inequalities using GroupBy + FreeQ vs the current Select approach *)
ProfileStrategyC[caseKey_] := Module[
  {data, Eqs, PreEqs, js, jts, NewSystem, InitRules, costpluscurrents,
   Ncpc, approxJs, ineqs, timeBaseline, timeOptimized,
   resultBaseline, resultOptimized},

  data = GetExampleData[caseKey] /. GetParams[caseKey];
  Eqs = DataToEquations[data];
  ClearSolveCache[];
  PreEqs = MFGPreprocessing[Eqs];
  js = Lookup[PreEqs, "js"];
  jts = Lookup[PreEqs, "jts"];
  NewSystem = Lookup[PreEqs, "NewSystem"];
  InitRules = Lookup[PreEqs, "InitRules"];
  costpluscurrents = Lookup[PreEqs, "costpluscurrents"];

  (* Simulate critical congestion substitution *)
  approxJs = AssociationThread[js, 0 js];
  Ncpc = RoundValues @ (Expand /@ (costpluscurrents /. approxJs));
  InitRules = Expand /@ (InitRules /. Ncpc);
  NewSystem = NewSystem /. Ncpc;

  If[NewSystem[[2]] === True,
    Return[<|"Baseline_Time" -> 0., "StrategyC_Time" -> 0., "Speedup" -> 1.,
             "Note" -> "No inequalities to group"|>]
  ];

  ineqs = NewSystem[[2]];

  (* Baseline: Select approach (current code) *)
  {timeBaseline, resultBaseline} = AbsoluteTiming[
    MFGParallelMap[Function[jt, Select[ineqs, !FreeQ[jt][#]&]], jts]
  ];

  (* Strategy C: Build index once, then lookup *)
  {timeOptimized, resultOptimized} = AbsoluteTiming[
    Module[{ineqList, index},
      ineqList = If[Head[ineqs] === And, List @@ ineqs, {ineqs}];
      (* Build reverse index: for each inequality, which jts appear in it *)
      index = Association[# -> {} & /@ jts];
      Do[
        Do[
          If[!FreeQ[ineq, jt],
            index[jt] = Append[index[jt], ineq]
          ],
          {jt, jts}
        ],
        {ineq, ineqList}
      ];
      (* Convert back to And expressions *)
      If[Length[#] == 0, True, And @@ #] & /@ Values[index]
    ]
  ];

  <|
    "Baseline_Time" -> timeBaseline,
    "StrategyC_Time" -> timeOptimized,
    "Speedup" -> timeBaseline / Max[timeOptimized, 0.0001],
    "NumInequalities" -> If[Head[ineqs] === And, Length[ineqs], 1],
    "NumTransitionFlows" -> Length[jts]
  |>
];

(* ============================================================ *)
(* PHASE 3: Full pipeline comparison (all strategies combined)  *)
(* ============================================================ *)

ProfileFullPipeline[caseKey_] := Module[
  {data, timeBaseline, timeOptimized, resultBaseline, resultOptimized},

  data = GetExampleData[caseKey] /. GetParams[caseKey];

  (* Baseline: standard CriticalCongestionSolver *)
  ClearSolveCache[];
  {timeBaseline, resultBaseline} = AbsoluteTiming[
    CriticalCongestionSolver[DataToEquations[data]]
  ];

  <|
    "CriticalCongestionSolver_Total" -> timeBaseline,
    "Status" -> Lookup[resultBaseline, "Status", "Unknown"]
  |>
];

(* ============================================================ *)
(* MAIN: Run all profiles                                       *)
(* ============================================================ *)

Print["╔══════════════════════════════════════════════════════╗"];
Print["║    MFGraphs Preprocessing Performance Profiler      ║"];
Print["╚══════════════════════════════════════════════════════╝"];
Print[];

(* Phase 1: Baseline *)
Print["━━━ PHASE 1: Baseline Step-by-Step Profiling ━━━"];
Print[];

Do[
  Print["──── Case: ", caseKey, " ────"];
  Module[{timings},
    timings = ProfilePreprocessingBaseline[caseKey];
    Do[
      Print["  ", StringPadRight[key, 40], " ", ToString[CForm[N[timings[key]]]], " s"],
      {key, Keys[timings]}
    ];
    Print[];
  ],
  {caseKey, $TestCases}
];

(* Helper for formatting *)
FmtTime[t_?NumericQ] := ToString[CForm[N[t]]];
FmtTime[x_] := ToString[x];

(* Phase 2: Strategy comparisons *)
Print["━━━ PHASE 2: Strategy A — Whole-System Simplify ━━━"];
Print[];

Do[
  Print["──── Case: ", caseKey, " ────"];
  Module[{result},
    result = ProfileStrategyA[caseKey];
    Print["  Baseline (per-element):  ", FmtTime[result["Baseline_Time"]], " s"];
    Print["  Strategy A (whole):      ", FmtTime[result["StrategyA_Time"]], " s"];
    Print["  Speedup:                 ", FmtTime[result["Speedup"]], "x"];
    Print["  Equivalent:              ", result["Equivalent"]];
    Print[];
  ],
  {caseKey, $TestCases}
];

Print["━━━ PHASE 2: Strategy B — Single TripleClean ━━━"];
Print[];

Do[
  Print["──── Case: ", caseKey, " ────"];
  Module[{result},
    result = ProfileStrategyB[caseKey];
    Print["  Baseline (two passes):   ", FmtTime[result["Baseline_Time"]], " s"];
    Print["  Strategy B (single):     ", FmtTime[result["StrategyB_Time"]], " s"];
    Print["  Speedup:                 ", FmtTime[result["Speedup"]], "x"];
    Print["  Rules equivalent:        ", result["RulesEquivalent"]];
    Print["  Baseline rules:          ", result["BaselineRuleCount"]];
    Print["  Optimized rules:         ", result["OptimizedRuleCount"]];
    Print[];
  ],
  {caseKey, $TestCases}
];

Print["━━━ PHASE 2: Strategy C — Index-Based Inequality Grouping ━━━"];
Print[];

Do[
  Print["──── Case: ", caseKey, " ────"];
  Module[{result},
    result = ProfileStrategyC[caseKey];
    Print["  Baseline (Select):       ", FmtTime[result["Baseline_Time"]], " s"];
    Print["  Strategy C (Index):      ", FmtTime[result["StrategyC_Time"]], " s"];
    Print["  Speedup:                 ", FmtTime[result["Speedup"]], "x"];
    If[KeyExistsQ[result, "NumInequalities"],
      Print["  Inequalities:            ", result["NumInequalities"]];
      Print["  Transition flows:        ", result["NumTransitionFlows"]];
    ];
    Print[];
  ],
  {caseKey, $TestCases}
];

(* Phase 3: Full pipeline baseline *)
Print["━━━ PHASE 3: Full Pipeline Baseline ━━━"];
Print[];

Do[
  Print["──── Case: ", caseKey, " ────"];
  Module[{result},
    result = ProfileFullPipeline[caseKey];
    Print["  CriticalCongestionSolver: ", FmtTime[result["CriticalCongestionSolver_Total"]], " s"];
    Print["  Status:                   ", result["Status"]];
    Print[];
  ],
  {caseKey, $TestCases}
];

Print["Done."];

(* :!CodeAnalysis::EndBlock:: *)
