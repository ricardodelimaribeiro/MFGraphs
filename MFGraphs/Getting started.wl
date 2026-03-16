(* ::Package:: *)

(* MFGraphs getting-started script.
   Evaluate sections in order from a fresh kernel. Legacy solver outputs remain
   the default; add "ReturnShape" -> "Standard" to receive a normalized result
   envelope with keys such as "Solver", "ResultKind", "Feasibility", and
   "Solution". *)

<< MFGraphs`

$MFGraphsVerbose = False;

(* --- Critical congestion example --- *)

criticalData = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
criticalD2E = DataToEquations[criticalData];

criticalLegacy = CriticalCongestionSolver[criticalD2E];
criticalLegacy["AssoCritical"];

criticalStandard = CriticalCongestionSolver[
    criticalD2E,
    "ReturnShape" -> "Standard"
];
criticalStandard["Solution"];

(* --- Non-linear example --- *)

Block[{V = Function[{x, edge}, 0],
       alpha = Function[edge, 1],
       g = Function[{m, edge}, -1/m^2}},
    nonlinearData = GetExampleData[12] /. {I1 -> 2, U1 -> 0};
    nonlinearD2E = DataToEquations[nonlinearData];
    nonlinearLegacy = NonLinearSolver[nonlinearD2E, "MaxIterations" -> 5];
    nonlinearLegacy["AssoNonCritical"];

    nonlinearStandard = NonLinearSolver[
        nonlinearD2E,
        "MaxIterations" -> 5,
        "ReturnShape" -> "Standard"
    ];
    nonlinearStandard["Solution"];
];

(* --- Monotone example --- *)

Block[{V = Function[{x, edge}, 0],
       alpha = Function[edge, 1],
       g = Function[{m, edge}, -1/m^2}},
    monotoneData = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
    monotoneLegacy = MonotoneSolverFromData[monotoneData, "TimeSteps" -> 20];
    monotoneLegacy;

    monotoneStandard = MonotoneSolverFromData[
        monotoneData,
        "TimeSteps" -> 20,
        "ReturnShape" -> "Standard"
    ];
    monotoneStandard["Solution"];
];
