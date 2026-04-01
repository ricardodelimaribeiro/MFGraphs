#!/usr/bin/env wolframscript
(* Compare stationary solver outputs on standardized result contracts. *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
MFGraphs`$MFGraphsVerbose = False;

$HamiltonianOptions = {
  "PotentialFunction" -> Function[{x, edge}, 0],
  "CongestionExponentFunction" -> Function[edge, 1],
  "InteractionFunction" -> Function[{m, edge}, -1/m^2]
};

$ComparisonCases = <|
  "case3" -> <|
    "Key" -> 3,
    "Params" -> {I1 -> 80, U1 -> 0},
    "CriticalVsNonLinearTolerance" -> 10^-5,
    "CriticalVsMonotoneTolerance" -> 10^-4,
    "MonotoneOptions" -> {"ResidualTolerance" -> 10^-6, "MaxTime" -> 10, "MaxSteps" -> 2000}
  |>,
  "case7" -> <|
    "Key" -> 7,
    "Params" -> {I1 -> 80, U1 -> 0, U2 -> 10},
    "CriticalVsNonLinearTolerance" -> 10^-5,
    "CriticalVsMonotoneTolerance" -> 10^-4,
    "MonotoneOptions" -> {"ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000}
  |>,
  "case8" -> <|
    "Key" -> 8,
    "Params" -> {I1 -> 100, U1 -> 0, U2 -> 10, S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1},
    "CriticalVsNonLinearTolerance" -> 10^-5,
    "CriticalVsMonotoneTolerance" -> 10^-4,
    "MonotoneOptions" -> {"ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000}
  |>,
  "case18" -> <|
    "Key" -> 18,
    "Params" -> {I1 -> 80, U1 -> 0, S1 -> 2, S2 -> 1},
    "CriticalVsNonLinearTolerance" -> 10^-5,
    "CriticalVsMonotoneTolerance" -> 10^-4,
    "MonotoneOptions" -> {"ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000}
  |>,
  "hrf-scenario1" -> <|
    "Key" -> "HRF Scenario 1",
    "Params" -> {I1 -> 100, I2 -> 100, U1 -> 0, U2 -> 0},
    "CriticalVsNonLinearTolerance" -> 10^-5,
    "CriticalVsMonotoneTolerance" -> 10^-4,
    "MonotoneOptions" -> {"ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000}
  |>
|>;

$SmokeComparisonLabels = {"case3", "case7", "case8", "case18"};
$PaperComparisonLabels = {"hrf-scenario1"};
$RoutineComparisonLabels = Complement[Keys[$ComparisonCases], $PaperComparisonLabels];

NormalizeLabel[arg_String] := ToLowerCase[arg];
NormalizeLabel[arg_] := arg;

RequestedLabel = NormalizeLabel @ If[Length[Rest[$ScriptCommandLine]] > 0,
  First[Rest[$ScriptCommandLine]],
  "all"
];

SuccessResultQ[result_Association] :=
  Lookup[result, "ResultKind", "Failure"] === "Success" &&
  AssociationQ[Lookup[result, "Solution", Missing["NotAvailable"]]];

VectorDifference[a_, b_] /; VectorQ[a, NumericQ] && VectorQ[b, NumericQ] :=
  Max[Abs[a - b]];
VectorDifference[_, _] := Infinity;

CompareCase[label_String, spec_Association] :=
  Module[{data, d2e, critical, nonlinear, monotone, criticalVsNonLinear,
      criticalVsMonotone, resultKinds, kirchhoffResiduals, edgesAgree, pass},
    Print["--- ", label, " (", spec["Key"], ") ---"];
    data = GetExampleData[spec["Key"]] /. spec["Params"];
    d2e = DataToEquations[data];

    critical = Quiet[CriticalCongestionSolver[d2e]];
    nonlinear = Quiet[NonLinearSolver[d2e, Sequence @@ $HamiltonianOptions, "Tolerance" -> 10^-8, "MaxIterations" -> 20]];
    monotone = Quiet[MonotoneSolverFromData[data, Sequence @@ $HamiltonianOptions, Sequence @@ spec["MonotoneOptions"]]];

    resultKinds = Lookup[#, "ResultKind", "Missing"] & /@ {critical, nonlinear, monotone};
    Print["  Result kinds: ", resultKinds];
    If[!And @@ (SuccessResultQ /@ {critical, nonlinear, monotone}),
      Print["  FAIL: at least one solver did not converge successfully."];
      Return[False]
    ];

    edgesAgree =
      critical["ComparableEdges"] === nonlinear["ComparableEdges"] &&
      nonlinear["ComparableEdges"] === monotone["ComparableEdges"];
    criticalVsNonLinear =
      VectorDifference[critical["ComparableFlowVector"], nonlinear["ComparableFlowVector"]];
    criticalVsMonotone =
      VectorDifference[critical["ComparableFlowVector"], monotone["ComparableFlowVector"]];
    kirchhoffResiduals = Lookup[#, "KirchhoffResidual", Infinity] & /@ {critical, nonlinear, monotone};

    Print["  Critical vs NonLinear: ", criticalVsNonLinear];
    Print["  Critical vs Monotone: ", criticalVsMonotone];
    Print["  Kirchhoff residuals: ", kirchhoffResiduals];

    pass =
      edgesAgree &&
      criticalVsNonLinear < spec["CriticalVsNonLinearTolerance"] &&
      criticalVsMonotone < spec["CriticalVsMonotoneTolerance"] &&
      Max[kirchhoffResiduals] < 10^-6;

    If[pass,
      Print["  PASS"],
      Print["  FAIL: standardized outputs are not comparable within tolerance."]
    ];
    pass
  ];

SelectedLabels =
  Which[
    RequestedLabel === "all",
      $RoutineComparisonLabels,
    RequestedLabel === "smoke",
      $SmokeComparisonLabels,
    RequestedLabel === "paper",
      $PaperComparisonLabels,
    RequestedLabel === "exhaustive",
      Keys[$ComparisonCases],
    True,
      Select[Keys[$ComparisonCases], # === RequestedLabel &]
  ];

If[SelectedLabels === {},
  Print["Unknown case label: ", RequestedLabel];
  Print["Available labels: ", StringRiffle[Keys[$ComparisonCases], ", "], ", smoke, paper, all, exhaustive"];
  Exit[1];
];

AllPass = And @@ (CompareCase[#, $ComparisonCases[#]] & /@ SelectedLabels);

If[AllPass,
  Print["SUCCESS: standardized stationary solver outputs agree across all requested cases."];
  Exit[0],
  Print["FAILURE: at least one standardized solver comparison failed."];
  Exit[1]
];

(* :!CodeAnalysis::EndBlock:: *)
