(* Wolfram Language Test file *)
(* Thin-routing tests for SolveMFG (phase 1 unified API). *)

Test[
    Module[{data, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        result = Quiet[SolveMFG[data]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None}
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: Automatic defaults to critical congestion"
]

Test[
    Module[{data, d2e, direct, routed, drift},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        direct = Quiet[CriticalCongestionSolver[d2e]];
        routed = Quiet[SolveMFG[d2e, Method -> "CriticalCongestion"]];
        drift = Max[Abs[direct["ComparableFlowVector"] - routed["ComparableFlowVector"]]];
        Lookup[direct, {"ResultKind", "Feasibility", "Message"}] ===
            Lookup[routed, {"ResultKind", "Feasibility", "Message"}] &&
        NumericQ[drift] && drift <= 10^-8
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: CriticalCongestion method matches direct solver on compiled input"
]

Test[
    Module[{data, result},
        data = GetExampleData[3] /. {I1 -> 2, U1 -> 0};
        result = Quiet[
            SolveMFG[
                data,
                Method -> "NonLinear",
                "MaxIterations" -> 20,
                "Tolerance" -> 10^-8,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"NonLinear", "Success", "Feasible", None}
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: NonLinear method dispatches from raw data"
]

Test[
    Module[{data, result},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        result = Quiet[
            SolveMFG[
                data,
                Method -> "Monotone",
                "ResidualTolerance" -> 10^-6,
                "MaxTime" -> 10,
                "MaxSteps" -> 2000,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, "Solver", None] === "Monotone" &&
        MemberQ[{"Success", "NonConverged"}, Lookup[result, "ResultKind", None]]
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: Monotone method dispatches from raw data"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            SolveMFG[
                d2e,
                Method -> "Monotone",
                "ResidualTolerance" -> 10^-6,
                "MaxTime" -> 10,
                "MaxSteps" -> 2000,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, "Solver", None] === "Monotone" &&
        AssociationQ[Lookup[result, "Convergence", Missing["NotAvailable"]]]
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: Monotone method accepts compiled equation associations"
]

Test[
    Module[{result},
        result = SolveMFG[<||>, Method -> "NotAMethod"];
        Lookup[result, {"Solver", "ResultKind", "Message", "Method"}] ===
            {"SolveMFG", "Failure", "UnknownMethod", "NotAMethod"}
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: unknown method returns failure envelope"
]

(* Automatic cascade tests (issue #23): trace + strict fallback/abort behavior. *)

Test[
    Module[{compiled, result, trace},
        compiled = <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}|>;
        result = Block[{CriticalCongestionSolver, MonotoneSolver, NonLinearSolver},
            CriticalCongestionSolver[___] =
                <|
                    "Solver" -> "CriticalCongestion",
                    "ResultKind" -> "Failure",
                    "Feasibility" -> "Infeasible",
                    "Message" -> "NoSolution",
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> Missing["NotApplicable"]
                |>;
            MonotoneSolver[___] =
                <|
                    "Solver" -> "Monotone",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>
                |>;
            NonLinearSolver[___] =
                <|
                    "Solver" -> "NonLinear",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>
                |>;
            SolveMFG[compiled, Method -> "Automatic"]
        ];
        trace = Lookup[result, "MethodTrace", {}];
        Lookup[result, {"Solver", "MethodUsed", "ResultKind", "Feasibility", "Message"}] ===
            {"Monotone", "Monotone", "Success", "Feasible", None} &&
        Length[trace] === 2 &&
        Lookup[trace[[1]], "Decision", None] === "Fallback" &&
        Lookup[trace[[2]], "Decision", None] === "Done"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: recoverable critical failure cascades to monotone"
]

Test[
    Module[{compiled, result, trace},
        compiled = <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}|>;
        result = Block[{CriticalCongestionSolver, MonotoneSolver, NonLinearSolver},
            CriticalCongestionSolver[___] =
                <|
                    "Solver" -> "CriticalCongestion",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Infeasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>
                |>;
            MonotoneSolver[___] =
                <|
                    "Solver" -> "Monotone",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>
                |>;
            NonLinearSolver[___] =
                <|
                    "Solver" -> "NonLinear",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>
                |>;
            SolveMFG[compiled, Method -> "Automatic"]
        ];
        trace = Lookup[result, "MethodTrace", {}];
        Lookup[result, {"Solver", "MethodUsed", "ResultKind", "Feasibility"}] ===
            {"Monotone", "Monotone", "Success", "Feasible"} &&
        Length[trace] === 2 &&
        Lookup[trace[[1]], "ResultKind", None] === "Success" &&
        Lookup[trace[[1]], "Feasibility", None] === "Infeasible" &&
        Lookup[trace[[1]], "Decision", None] === "Fallback"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: infeasible false-success critical result falls back"
]

Test[
    Module[{compiled, result, trace},
        compiled = <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}|>;
        result = Block[{CriticalCongestionSolver, MonotoneSolver, NonLinearSolver},
            CriticalCongestionSolver[___] =
                <|
                    "Solver" -> "CriticalCongestion",
                    "ResultKind" -> "Failure",
                    "Feasibility" -> "Infeasible",
                    "Message" -> "NoSolution",
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> Missing["NotApplicable"]
                |>;
            MonotoneSolver[___] =
                <|
                    "Solver" -> "Monotone",
                    "ResultKind" -> "NonConverged",
                    "Feasibility" -> "Infeasible",
                    "Message" -> "ResidualExceedsTolerance",
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> <|"StopReason" -> "MaxTimeReached"|>
                |>;
            NonLinearSolver[___] =
                <|
                    "Solver" -> "NonLinear",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>
                |>;
            SolveMFG[compiled, Method -> "Automatic"]
        ];
        trace = Lookup[result, "MethodTrace", {}];
        Lookup[result, {"Solver", "MethodUsed", "ResultKind", "Feasibility", "Message"}] ===
            {"NonLinear", "NonLinear", "Success", "Feasible", None} &&
        Length[trace] === 3 &&
        Lookup[trace[[1]], "Decision", None] === "Fallback" &&
        Lookup[trace[[2]], "Decision", None] === "Fallback" &&
        Lookup[trace[[3]], "Decision", None] === "Done"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: cascades through nonlinear on dual recoverable failures"
]

Test[
    Module[{compiled, result, trace},
        compiled = <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}|>;
        result = Block[{CriticalCongestionSolver},
            CriticalCongestionSolver[___] = $Failed;
            SolveMFG[compiled, Method -> "Automatic"]
        ];
        trace = Lookup[result, "MethodTrace", {}];
        Lookup[result, {"Solver", "ResultKind", "Message", "MethodUsed"}] ===
            {"SolveMFG", "Failure", "InvalidSolverReturnShape", "CriticalCongestion"} &&
        Length[trace] === 1 &&
        Lookup[trace[[1]], "Decision", None] === "Abort"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: malformed solver return triggers strict abort"
]

Test[
    Module[{compiled, result, trace, nonLinearCalls = 0},
        compiled = <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}|>;
        result = Block[{CriticalCongestionSolver, MonotoneSolver, NonLinearSolver},
            CriticalCongestionSolver[___] =
                <|
                    "Solver" -> "CriticalCongestion",
                    "ResultKind" -> "Failure",
                    "Feasibility" -> "Infeasible",
                    "Message" -> "NoSolution",
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> Missing["NotApplicable"]
                |>;
            MonotoneSolver[___] =
                <|
                    "Solver" -> "Monotone",
                    "ResultKind" -> "Degenerate",
                    "Feasibility" -> Missing["NotApplicable"],
                    "Message" -> "DegenerateCase",
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> <|"StopReason" -> "DegenerateCase"|>
                |>;
            NonLinearSolver[___] :=
                (
                    nonLinearCalls++;
                    <|
                        "Solver" -> "NonLinear",
                        "ResultKind" -> "Success",
                        "Feasibility" -> "Feasible",
                        "Message" -> None,
                        "Solution" -> <||>,
                        "Convergence" -> <||>
                    |>
                );
            SolveMFG[compiled, Method -> "Automatic"]
        ];
        trace = Lookup[result, "MethodTrace", {}];
        Lookup[result, {"Solver", "MethodUsed", "ResultKind", "Message"}] ===
            {"Monotone", "Monotone", "Degenerate", "DegenerateCase"} &&
        nonLinearCalls === 0 &&
        Length[trace] === 2 &&
        Lookup[trace[[2]], "Decision", None] === "Done"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: monotone degenerate case is terminal"
]

Test[
    Module[{data, result, trace, methodUsed},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        result = Quiet[SolveMFG[data, Method -> "Automatic"]];
        trace = Lookup[result, "MethodTrace", {}];
        methodUsed = Lookup[result, "MethodUsed", Missing["NotAvailable"]];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        ListQ[trace] && Length[trace] >= 1 &&
        StringQ[methodUsed]
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: real-data run returns method trace and method used"
]
