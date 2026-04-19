Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

ClearAll[printable, runCase];
printable[expr_] := StringReplace[
  ToString[InputForm[expr]],
  "MFGraphs`Private`" -> ""
];

runCase[label_, switchingCosts_] := Module[{data, d2e, consistency, result},
  data = GetExampleData["chain with two exits"] /. {I1 -> 100, U1 -> 0, U2 -> 0};
  data = AssociateTo[data, "Switching Costs" -> switchingCosts];
  Print["CASE=" <> label];
  Print["DATA_SWITCHING=" <> printable[data["Switching Costs"]]];
  d2e = Quiet[DataToEquations[data], {DataToEquations::switchingcosts}];
  consistency = IsSwitchingCostConsistent[Normal @ d2e["SwitchingCosts"]];
  Print["CONSISTENCY_D2E=" <> printable[consistency]];
  Print["ALT_TRANSITION_FLOWS_HEAD=" <> printable[Head[d2e["AltTransitionFlows"]]]];
  result = Quiet[CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 60], {DataToEquations::switchingcosts}];
  Print["FEASIBILITY=" <> printable[Lookup[result, "Feasibility", Missing[]]]];
  Print["MESSAGE=" <> printable[Lookup[result, "Message", Missing[]]]];
  Print["ASSO_CRITICAL=" <> printable[Lookup[result, "AssoCritical", Missing[]]]];
  Print["---"];
];

runCase["baseline_no_switching", {}];
runCase["requested_single_switch_123_cost10", {{1, 2, 3, 10}}];
runCase["symmetric_switch_123_and_321_cost10", {{1, 2, 3, 10}, {3, 2, 1, 10}}];
