SetDirectory[DirectoryName[$InputFileName]];
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];

ClearAll[FormatExpr, ClauseTriples, StandardTransitionClauses, ActualTransitionClauses, AnalyzeExample];

FormatExpr[expr_] := ToString[Unevaluated[expr], InputForm];

ClauseTriples[expr_] :=
  SortBy[
    Cases[Unevaluated[expr], Equal[j[a_, b_, c_], 0] :> {a, b, c}, Infinity],
    FormatExpr
  ];

StandardTransitionClauses[auxTriples_List] :=
  DeleteDuplicates @ Map[
    ClauseTriples,
    AltFlowOp[j] /@ Select[auxTriples, OrderedQ[{First[#], Last[#]}] &]
  ];

ActualTransitionClauses[altTransitionFlows_] :=
  DeleteDuplicates @ Map[
    ClauseTriples,
    If[Head[altTransitionFlows] === And, List @@ altTransitionFlows, {altTransitionFlows}]
  ];

AnalyzeExample[name_String] :=
  Module[{data, d2e, switchingCosts, consistentQ, altTransitionFlows, auxTriples,
      standardClauses, actualClauses, extraClauses, critical, criticalResult,
      uAssignments, criticalResultKind, criticalFeasibility, criticalMessage,
      renderedExtraClauses},
    data = GetExampleData[name];
    d2e = DataToEquations[data];
    switchingCosts = Normal @ Lookup[d2e, "SwitchingCosts", <||>];
    consistentQ = IsSwitchingCostConsistent[switchingCosts];
    altTransitionFlows = Lookup[d2e, "AltTransitionFlows", Missing["NotAvailable"]];
    auxTriples = Lookup[d2e, "auxTriples", {}];
    standardClauses = StandardTransitionClauses[auxTriples];
    actualClauses = ActualTransitionClauses[altTransitionFlows];
    extraClauses = Select[actualClauses, FreeQ[standardClauses, #] &];
    renderedExtraClauses =
      FormatExpr /@ Map[
        HoldForm[(j @@ #[[1]] == 0) || (j @@ #[[2]] == 0)] &,
        Select[extraClauses, Length[#] == 2 &]
      ];
    critical = Quiet @ Check[CriticalCongestionSolver[d2e], $Failed];
    criticalResultKind = If[AssociationQ[critical], Lookup[critical, "ResultKind", Missing["NotAvailable"]], Missing["NotAvailable"]];
    criticalFeasibility = If[AssociationQ[critical], Lookup[critical, "Feasibility", Missing["NotAvailable"]], Missing["NotAvailable"]];
    criticalMessage = If[AssociationQ[critical], Lookup[critical, "Message", Missing["NotAvailable"]], "No symbolic solution returned"];
    criticalResult = If[AssociationQ[critical], Lookup[critical, "AssoCritical", <||>], <||>];
    uAssignments = If[AssociationQ[criticalResult], Select[Normal @ criticalResult, MatchQ[First[#], _u] &], {}];
    <|
      "Name" -> name,
      "ConsistentQ" -> consistentQ,
      "SwitchingCosts" -> (FormatExpr /@ switchingCosts),
      "StandardClauseCount" -> Length[standardClauses],
      "ActualClauseCount" -> Length[actualClauses],
      "ExtraCrossPairClauseCount" -> Length[extraClauses],
      "ExtraCrossPairClauses" -> renderedExtraClauses,
      "AltTransitionFlows" -> FormatExpr[altTransitionFlows],
      "CriticalResultKind" -> FormatExpr[criticalResultKind],
      "CriticalFeasibility" -> FormatExpr[criticalFeasibility],
      "CriticalMessage" -> FormatExpr[criticalMessage],
      "CriticalUAssignments" -> (FormatExpr /@ (Rule @@@ uAssignments))
    |>
  ];

result = AnalyzeExample /@ {
  "Inconsistent Y shortcut",
  "Inconsistent attraction shortcut"
};

Export[FileNameJoin[{Directory[], "Results", "inconsistent_examples_verification.json"}], result, "RawJSON"];
Print[ExportString[result, "RawJSON", "Compact" -> False]];
