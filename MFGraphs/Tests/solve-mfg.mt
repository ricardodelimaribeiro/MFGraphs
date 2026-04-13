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

(* ================================================================ *)
(* Alpha-aware dispatch tests (issue #23)                           *)
(* ================================================================ *)

(* Helper: run SolveMFG with mocked solvers, returning {result, criticalCalls}. *)
solveMFGWithMockedSolvers[alphaSpec_] :=
    Module[{compiled, result, criticalCalls = 0},
        compiled = <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}|>;
        result = Block[{CriticalCongestionSolver, MonotoneSolver, NonLinearSolver},
            CriticalCongestionSolver[___] :=
                (
                    criticalCalls++;
                    <|
                        "Solver" -> "CriticalCongestion",
                        "ResultKind" -> "Success",
                        "Feasibility" -> "Feasible",
                        "Message" -> None,
                        "Solution" -> <||>,
                        "Convergence" -> Missing["NotApplicable"]
                    |>
                );
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
            SolveMFG[compiled, Method -> "Automatic",
                "CongestionExponentFunction" -> alphaSpec]
        ];
        {result, criticalCalls}
    ];

solveMFGExplicitCriticalWithMockCompile[alphaSpec_] :=
    Module[{data, result, criticalCalls = 0, compileCallCount = 0},
        data = <|"Mock" -> True|>;
        result = Block[{SolveMFGCompileInput, CriticalCongestionSolver},
            SolveMFGCompileInput[input_, opts_List:{}] :=
                (
                    compileCallCount++;
                    <|"EqGeneral" -> True, "js" -> {}, "jts" -> {}, "CompiledInput" -> input, "CompileOpts" -> opts|>
                );
            CriticalCongestionSolver[eqs_, ___] :=
                (
                    criticalCalls++;
                    <|
                        "Solver" -> "CriticalCongestion",
                        "ResultKind" -> "Success",
                        "Feasibility" -> "Feasible",
                        "Message" -> None,
                        "Solution" -> eqs,
                        "Convergence" -> Missing["NotApplicable"]
                    |>
                );
            SolveMFG[
                data,
                Method -> "CriticalCongestion",
                "CongestionExponentFunction" -> alphaSpec
            ]
        ];
        {result, compileCallCount, criticalCalls}
    ];

Test[
    Module[{result, criticalCalls, trace},
        {result, criticalCalls} = solveMFGWithMockedSolvers[2&];
        trace = Lookup[result, "MethodTrace", {}];
        criticalCalls === 0 &&
        Lookup[result, "MethodUsed", None] === "Monotone" &&
        Length[trace] === 1 &&
        Lookup[trace[[1]], "Method", None] === "Monotone"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: alpha!=1 skips CriticalCongestion"
]

Test[
    Module[{result, criticalCalls},
        {result, criticalCalls} = solveMFGWithMockedSolvers[<|e[1, 2] -> 1.5|>];
        criticalCalls === 0 &&
        Lookup[result, "MethodUsed", None] === "Monotone"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: Association alpha!=1 skips CriticalCongestion"
]

Test[
    Module[{result, criticalCalls},
        {result, criticalCalls} = solveMFGWithMockedSolvers[<|e[1, 2] -> 1, e[2, 3] -> 1|>];
        criticalCalls === 1 &&
        Lookup[result, "MethodUsed", None] === "CriticalCongestion"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: uniform alpha=1 Association still routes to Critical"
]

Test[
    Module[{result, criticalCalls},
        {result, criticalCalls} = solveMFGWithMockedSolvers[1&];
        criticalCalls === 1 &&
        Lookup[result, "MethodUsed", None] === "CriticalCongestion"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: pure-function alpha=1 still routes to Critical"
]

Test[
    Module[{result, criticalCalls},
        {result, criticalCalls} = solveMFGWithMockedSolvers[Function[edge, 1]];
        criticalCalls === 1 &&
        Lookup[result, "MethodUsed", None] === "CriticalCongestion"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: named Function alpha=1 still routes to Critical"
]

Test[
    Module[{result, compileCallCount, criticalCalls},
        {result, compileCallCount, criticalCalls} = solveMFGExplicitCriticalWithMockCompile[1.5];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message", "Solution"}] === {
            "CriticalCongestion",
            "Failure",
            Missing["NotAvailable"],
            "Method -> 'CriticalCongestion' only supports CongestionExponentFunction -> 1.",
            Missing["NotAvailable"]
        } &&
        compileCallCount === 0 &&
        criticalCalls === 0
    ]
    ,
    True
    ,
    TestID -> "SolveMFG critical routing: scalar alpha!=1 fails before compile"
]

Test[
    Module[{result, compileCallCount, criticalCalls},
        {result, compileCallCount, criticalCalls} =
            solveMFGExplicitCriticalWithMockCompile[<|e[1, 2] -> 2|>];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message", "Solution"}] === {
            "CriticalCongestion",
            "Failure",
            Missing["NotAvailable"],
            "Method -> 'CriticalCongestion' only supports CongestionExponentFunction -> 1.",
            Missing["NotAvailable"]
        } &&
        compileCallCount === 0 &&
        criticalCalls === 0
    ]
    ,
    True
    ,
    TestID -> "SolveMFG critical routing: association alpha!=1 fails before compile"
]

Test[
    Module[{result, compileCallCount, criticalCalls},
        {result, compileCallCount, criticalCalls} = solveMFGExplicitCriticalWithMockCompile[1&];
        compileCallCount > 0 &&
        criticalCalls === 1 &&
        Lookup[result, "Solver", None] === "CriticalCongestion"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG critical routing: pure-function alpha=1 compiles and runs critical"
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
                "CongestionExponentFunction" -> 1.5
            ]
        ];
        Lookup[result, "Solver", None] === "NonLinear" &&
        MemberQ[{"Success", "NonConverged"}, Lookup[result, "ResultKind", None]]
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: Hamiltonian options forwarded to NonLinear with scalar alpha"
]

Test[
    Module[{data, result, criticalCalls = 0},
        data = <|"Mock" -> True|>;
        result = Block[{DataToEquations, CriticalCongestionSolver, MonotoneSolver, NonLinearSolver},
            DataToEquations[_] :=
                <|
                    "EqGeneral" -> alpha[{1, 2}],
                    "js" -> {},
                    "jts" -> {},
                    "CompiledAlpha" -> alpha[{1, 2}]
                |>;
            CriticalCongestionSolver[___] := (criticalCalls++; $Failed);
            MonotoneSolver[___] =
                <|
                    "Solver" -> "Monotone",
                    "ResultKind" -> "Failure",
                    "Feasibility" -> "Infeasible",
                    "Message" -> "MockFailure",
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> Missing["NotApplicable"]
                |>;
            NonLinearSolver[d2e_, ___] :=
                <|
                    "Solver" -> "NonLinear",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>,
                    "CompiledAlpha" -> Lookup[d2e, "CompiledAlpha", Missing["NotAvailable"]]
                |>;
            SolveMFG[
                data,
                Method -> "Automatic",
                "CongestionExponentFunction" -> 1.5
            ]
        ];
        criticalCalls === 0 &&
        Lookup[result, "MethodUsed", None] === "NonLinear" &&
        Lookup[result, "CompiledAlpha", Missing["NotAvailable"]] === 1.5
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: raw alpha override recompiles equations before NonLinear fallback"
]

Test[
    Module[{data, result},
        data = <|"Mock" -> True|>;
        result = Block[{DataToEquations, NonLinearSolver},
            DataToEquations[_] :=
                <|
                    "EqGeneral" -> alpha[e[1, 2]],
                    "js" -> {},
                    "jts" -> {},
                    "CompiledAlpha" -> alpha[e[1, 2]]
                |>;
            NonLinearSolver[d2e_, ___] :=
                <|
                    "Solver" -> "NonLinear",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>,
                    "CompiledAlpha" -> Lookup[d2e, "CompiledAlpha", Missing["NotAvailable"]]
                |>;
            SolveMFG[
                data,
                Method -> "NonLinear",
                "CongestionExponentFunction" -> <|e[1, 2] -> 1.5|>
            ]
        ];
        Lookup[result, {"Solver", "CompiledAlpha"}] === {"NonLinear", 1.5}
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: raw NonLinear compile honors Association alpha override"
]

Test[
    Module[{data, result},
        data = <|"Mock" -> True|>;
        result = Block[{DataToEquations, MonotoneSolver},
            DataToEquations[_] :=
                <|
                    "EqGeneral" -> alpha[e[1, 2]],
                    "js" -> {},
                    "jts" -> {},
                    "CompiledAlpha" -> alpha[e[1, 2]]
                |>;
            MonotoneSolver[d2e_, ___] :=
                <|
                    "Solver" -> "Monotone",
                    "ResultKind" -> "Success",
                    "Feasibility" -> "Feasible",
                    "Message" -> None,
                    "Solution" -> <||>,
                    "Convergence" -> <||>,
                    "CompiledAlpha" -> Lookup[d2e, "CompiledAlpha", Missing["NotAvailable"]]
                |>;
            MonotoneSolverFromData[
                data,
                "CongestionExponentFunction" -> <|e[1, 2] -> 1.5|>
            ]
        ];
        Lookup[result, {"Solver", "CompiledAlpha"}] === {"Monotone", 1.5}
    ]
    ,
    True
    ,
    TestID -> "MonotoneSolverFromData: raw compile honors Association alpha override"
]
(* NormalizeEdgeFunction unit tests *)

Test[
    Module[{f},
        f = NormalizeEdgeFunction[2.5];
        f[{1, 2}]
    ]
    ,
    2.5
    ,
    TestID -> "NormalizeEdgeFunction: scalar input returns constant function"
]

Test[
    Module[{f},
        f = NormalizeEdgeFunction[<|e[1, 2] -> 1.5, e[2, 3] -> 2.0|>];
        {f[e[1, 2]], f[e[2, 3]], f[e[3, 4]]}
    ]
    ,
    {1.5, 2.0, 1}
    ,
    TestID -> "NormalizeEdgeFunction: Association with default fallback"
]

Test[
    Module[{f},
        f = NormalizeEdgeFunction[Automatic];
        f[{1, 2}]
    ]
    ,
    1
    ,
    TestID -> "NormalizeEdgeFunction: Automatic returns default-valued function"
]

Test[
    Module[{f, custom},
        custom = Function[edge, Length[edge]];
        f = NormalizeEdgeFunction[custom];
        f[{1, 2}]
    ]
    ,
    2
    ,
    TestID -> "NormalizeEdgeFunction: Function passthrough"
]
