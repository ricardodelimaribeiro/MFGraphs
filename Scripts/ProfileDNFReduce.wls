#!/usr/bin/env wolframscript
(* Usage:
   wolframscript -file Scripts/ProfileDNFReduce.wls [--case name|all] [--timeout s]
     [--order original|nonor-stable|or-small|or-large|j-before-u|u-before-j|
              transition-before-edge|edge-before-transition|few-vars|path-aware|all]
     [--recursion-limit n]
*)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$SelectedCase = "grid-3x2";
$TimeoutSecs = 300;
$SelectedOrder = "original";
$RecursionLimitOpt = Automatic;

Module[{args = Rest[$ScriptCommandLine], i = 1, arg},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--case",
        i++; If[i <= Length[args], $SelectedCase = ToLowerCase[args[[i]]]],
      "--timeout",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[NumericQ[v] && v >= 0, $TimeoutSecs = v]]],
      "--order",
        i++; If[i <= Length[args], $SelectedOrder = ToLowerCase[args[[i]]]],
      "--recursion-limit",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[IntegerQ[v] && v > 0, $RecursionLimitOpt = v]]],
      _, Null
    ];
    i++;
  ]
];

$Orders = {
  "original",
  "nonor-stable",
  "or-small",
  "or-large",
  "j-before-u",
  "u-before-j",
  "transition-before-edge",
  "edge-before-transition",
  "few-vars",
  "path-aware"
};

If[$SelectedOrder =!= "all" && !MemberQ[$Orders, $SelectedOrder],
  Print["Unknown order: ", $SelectedOrder];
  Print["Available orders: ", StringRiffle[Append[$Orders, "all"], ", "]];
  Exit[2]
];

$RunOrders = If[$SelectedOrder === "all", $Orders, {$SelectedOrder}];

$Cases = {
  {"grid-2x2", Function[gridScenario[{2, 2}, {{1, 100}}, {{4, 0}}]]},
  {"grid-3x2", Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]]},
  {"grid-2x3", Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]]},
  {"grid-2x4", Function[gridScenario[{2, 4}, {{1, 100}}, {{8, 0}}]]},
  {"grid-3x3", Function[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]]},
  {"case-12",  Function[getExampleScenario[12, {{1, 100}}, {{4, 0}}]]},
  {"case-21",  Function[getExampleScenario[21, {{1, 50}, {2, 50}}, {{10, 0}, {11, 0}, {12, 0}}]]}
};
$CaseNames = $Cases[[All, 1]];

If[$SelectedCase =!= "all",
  $Cases = Select[$Cases, #[[1]] === $SelectedCase &];
  If[$Cases === {},
    Print["No case named '", $SelectedCase, "'. Available: ", StringRiffle[Append[$CaseNames, "all"], ", "]];
    Exit[2]
  ]
];

Fmt[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 3}, NumberPadding -> {"", "0"}]];
Fmt[x_] := ToString[x, InputForm];
Ms[t_?NumericQ] := N[1000 t];
AndLen[expr_] := If[Head[expr] === And, Length[List @@ expr], If[TrueQ[expr], 0, 1]];

ModelEdgeCount[sys_] :=
  With[{graph = systemData[sys, "Graph"]},
    If[MatchQ[graph, _Graph], EdgeCount[graph], Missing["Unavailable"]]
  ];

Print["MFGraphs DNF reducer diagnostic profile"];
Print["Date\t", DateString[]];
Print["TimeoutSeconds\t", $TimeoutSecs];
Print["Orders\t", StringRiffle[$RunOrders, ","]];
If[$RecursionLimitOpt =!= Automatic, Print["RecursionLimit\t", $RecursionLimitOpt]];
If[$SelectedCase =!= "all", Print["CaseFilter\t", $SelectedCase]];
Print[
  "Header\tCase\tOrder\tBuildMs\tStatus\tKind\tTransitionFlowStatus\tTransitionFlowRules\tTransitionFlowResiduals\tValid\tVertices\tEdges\tVars\tPreRules\tPreConjuncts\tOrs\tOrLeafCounts\t" <>
  "SolverInputMs\tDNFMs\tProcessed\tBranchesStarted\tBranchesKept\tBranchesDropped\tMaxDepth\tEqualitySolves\tSubstitutions\t" <>
  "ExpressionLeaves\tResultLeaves\tTimeoutPhase\tLastHead\tLastLeaves\tLastConjunct"
];

RunOne[name_String, buildFn_, order_String] :=
  Module[{tBuild, sys, report, summary, ordering, status, result, valid, timeoutLocation,
          lastConjunct, lastHead, lastLeaves, kind},
    {tBuild, sys} = AbsoluteTiming[makeSystem[buildFn[]]];
    report = Block[{$RecursionLimit = If[$RecursionLimitOpt === Automatic, $RecursionLimit, $RecursionLimitOpt]},
      solversTools`Private`dnfReduceDiagnosticReport[
        sys,
        "Order" -> order,
        "Timeout" -> $TimeoutSecs,
        "TraceLength" -> 12
      ]
    ];
    summary = Lookup[report, "Summary", <||>];
    ordering = Lookup[report, "OrderingSummary", <||>];
    status = Lookup[report, "Status", "Unknown"];
    result = Lookup[report, "Result", Missing["NoResult"]];
    kind = Lookup[summary, "ResultKind", If[result === $TimedOut, "TIMEOUT", "Unknown"]];
    valid = If[result === $TimedOut || status =!= "OK", "TimedOut", ToString[isValidSystemSolution[sys, result]]];
    timeoutLocation = Lookup[report, "TimeoutLocation", <||>];
    lastConjunct = Lookup[report, "LastConjunct", Missing["None"]];
    lastHead = Lookup[timeoutLocation, "LastConjunctHead", ToString[Head[lastConjunct]]];
    lastLeaves = Lookup[timeoutLocation, "LastConjunctLeafCount", Missing["Unavailable"]];
    Print[StringRiffle[
      {
        name,
        order,
        Fmt[Ms[tBuild]],
        status,
        kind,
        ToString[Lookup[summary, "TransitionFlowStatus", Missing["Unavailable"]]],
        ToString[Lookup[summary, "TransitionFlowRuleCount", Missing["Unavailable"]]],
        ToString[Lookup[summary, "TransitionFlowResidualCount", Missing["Unavailable"]]],
        valid,
        ToString[Length[systemData[sys, "Vertices"]]],
        ToString[ModelEdgeCount[sys]],
        ToString[Lookup[summary, "VariableCount", Missing["Unavailable"]]],
        ToString[Lookup[summary, "PreRuleCount", Missing["Unavailable"]]],
        ToString[Lookup[summary, "PreprocessingConjunctCount", Missing["Unavailable"]]],
        ToString[Lookup[ordering, "OrCount", Missing["Unavailable"]]],
        ToString[Lookup[ordering, "OrLeafCounts", Missing["Unavailable"]], InputForm],
        Fmt[Ms[Lookup[summary, "SolverInputSeconds", Missing["Unavailable"]]]],
        Fmt[Ms[Lookup[summary, "DNFSeconds", Missing["Unavailable"]]]],
        ToString[Lookup[summary, "ConjunctsProcessed", Missing["Unavailable"]]],
        ToString[Lookup[summary, "BranchesStarted", Missing["Unavailable"]]],
        ToString[Lookup[summary, "BranchesKept", Missing["Unavailable"]]],
        ToString[Lookup[summary, "BranchesDropped", Missing["Unavailable"]]],
        ToString[Lookup[summary, "MaxRecursionDepth", Missing["Unavailable"]]],
        ToString[Lookup[summary, "EqualitySolves", Missing["Unavailable"]]],
        ToString[Lookup[summary, "Substitutions", Missing["Unavailable"]]],
        ToString[Lookup[summary, "ExpressionLeafCount", Missing["Unavailable"]]],
        ToString[Lookup[summary, "ResultLeafCount", Missing["Unavailable"]]],
        ToString[Lookup[timeoutLocation, "Phase", Missing["Unavailable"]]],
        ToString[lastHead],
        ToString[lastLeaves],
        StringTake[ToString[lastConjunct, InputForm], UpTo[240]]
      },
      "\t"
    ]];
  ];

Scan[
  Function[case,
    Scan[
      Function[order, RunOne[case[[1]], case[[2]], order]],
      $RunOrders
    ]
  ],
  $Cases
];

(* :!CodeAnalysis::EndBlock:: *)
