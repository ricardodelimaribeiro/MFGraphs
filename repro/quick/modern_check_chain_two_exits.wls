Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

ClearAll[scenario, sys, d2e, result, assocCritical, printable];

printable[expr_] := StringReplace[
  ToString[InputForm[expr]],
  "MFGraphs`Private`" -> ""
];

(* Use the new Scenario and System API *)
(* Chain 1->2->3, entry at 1 flow 100, exits at 2 (cost 0) and 3 (cost 20) *)
(* Add inconsistent switching cost: 1->2->3 cost 10 *)
scenario = GridScenario[{3}, {{1, 100}}, {{2, 0}, {3, 20}}, {{1, 2, 3, 10}}];

Print["SCENARIO_VALID=" <> printable[ScenarioQ[scenario]]];

sys = makeSystem[scenario];
d2e = SystemDataFlatten[sys];

Print["D2E_KEYS=" <> printable[Sort @ Keys[d2e]]];
Print["D2E_SWITCHING=" <> printable[d2e["SwitchingCosts"]]];

result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 60];

Print["RESULT_KIND=" <> printable[Lookup[result, "ResultKind", Missing[]]]];
Print["FEASIBILITY=" <> printable[Lookup[result, "Feasibility", Missing[]]]];
Print["MESSAGE=" <> printable[Lookup[result, "Message", Missing[]]]];

assocCritical = Lookup[result, "AssoCritical", Missing[]];
Print["ASSO_CRITICAL=" <> printable[assocCritical]];
