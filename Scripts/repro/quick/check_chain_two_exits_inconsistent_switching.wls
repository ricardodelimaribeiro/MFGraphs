Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

ClearAll[data, consistency, d2e, result, assocCritical, printable];

printable[expr_] := StringReplace[
  ToString[InputForm[expr]],
  "MFGraphs`Private`" -> ""
];

data = GetExampleData["chain with two exits"] /. {I1 -> 100, U1 -> 0, U2 -> 20};
data = AssociateTo[data, "Switching" -> {{1, 2, 3, 10}}];

Print["DATA=" <> printable[data]];

consistency = IsSwitchingCostConsistent[Normal @ Association[Most /@ data["Switching"], Last /@ data["Switching"]]];
Print["CONSISTENCY=" <> printable[consistency]];

Quiet[
  d2e = DataToEquations[data];
  , {DataToEquations::switchingcosts}
];

Print["D2E_KEYS=" <> printable[Sort @ Keys[d2e]]];
Print["D2E_SWITCHING=" <> printable[d2e["SwitchingCosts"]]];
If[KeyExistsQ[d2e, "AllTransitionFlowsConstraint"],
  Print["HAS_ALL_TRANSITION_FLOWS_CONSTRAINT=True"],
  Print["HAS_ALL_TRANSITION_FLOWS_CONSTRAINT=False"]
];

Quiet[
  result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 60];
  , {DataToEquations::switchingcosts}
];

Print["RESULT_KIND=" <> printable[Lookup[result, "ResultKind", Missing[]]]];
Print["FEASIBILITY=" <> printable[Lookup[result, "Feasibility", Missing[]]]];
Print["MESSAGE=" <> printable[Lookup[result, "Message", Missing[]]]];
Print["SOLVER=" <> printable[Lookup[result, "Solver", Missing[]]]];

assocCritical = Lookup[result, "AssoCritical", Missing[]];
Print["ASSO_CRITICAL_HEAD=" <> printable[Head[assocCritical]]];
Print["ASSO_CRITICAL=" <> printable[assocCritical]];
