#!/usr/bin/env wolframscript
(* Compare stationary solver outputs on standardized result contracts. *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
MFGraphs`$MFGraphsVerbose = False;

Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

$HamiltonianOptions = {
  "PotentialFunction" -> Function[{x, edge}, 0],
  "CongestionExponentFunction" -> Function[edge, 1],
  "InteractionFunction" -> Function[{m, edge}, -1/m^2]
};

(* Use centralized case definitions where possible *)
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
    "Params" -> (Last @ Lookup[$NamedBenchmarkCases, "Y_noswitch"]),
    "CriticalVsNonLinearTolerance" -> 10^-5,
    "CriticalVsMonotoneTolerance" -> 10^-4,
    "MonotoneOptions" -> {"ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000}
  |>,
  "case8" -> <|
    "Key" -> 8,
    "Params" -> (Last @ Lookup[$NamedBenchmarkCases, "Y_switch"]),
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
    "Params" -> (Lookup[$CaseParams, "HRF Scenario 1", {}]),
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

    Print["  Kirchhoff residuals: ", kirchhoffResiduals];
    Print["  Max diff (Crit vs NL): ", criticalVsNonLinear];
    Print["  Max diff (Crit vs Mono): ", criticalVsMonotone];

    pass = True;
    If[criticalVsNonLinear > spec["CriticalVsNonLinearTolerance"],
      Print["  [FAIL] Critical vs NonLinear difference exceeds tolerance!"];
      pass = False;
    ];
    If[criticalVsMonotone > spec["CriticalVsMonotoneTolerance"],
      Print["  [FAIL] Critical vs Monotone difference exceeds tolerance!"];
      pass = False;
    ];
    If[!edgesAgree,
      Print["  [FAIL] Comparable edges do not match across solvers!"];
      pass = False;
    ];

    pass
  ];

results = Reap[
  Do[
    If[RequestedLabel === "all" || RequestedLabel === label,
      Sow[CompareCase[label, $ComparisonCases[label]]]
    ],
    {label, Keys[$ComparisonCases]}
  ]
][[2, 1]];

If[And @@ results,
  Print["\nPASS: All compared cases match within tolerance."];
  Exit[0],
  Print["\nFAIL: At least one case failed comparison."];
  Exit[1]
];
