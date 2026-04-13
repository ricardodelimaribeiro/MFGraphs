(* Wolfram Language Test file *)
(* Contract tests for the standardized stationary-solver return shape. *)

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None} &&
        AssociationQ[result["AssoCritical"]] &&
        result["Solution"] === result["AssoCritical"] &&
        AssociationQ[result["FlowAssociation"]] &&
        AssociationQ[result["SignedEdgeFlows"]] &&
        VectorQ[result["ComparableFlowVector"], NumericQ] &&
        NumericQ[result["KirchhoffResidual"]] &&
        NumericQ[result["BoundaryMassResidual"]] &&
        KeyExistsQ[result, "UtilityClassResidual"] &&
        IsFeasible[result]
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: critical congestion solver"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData["chain with two exits"] /. {I1 -> 80, U1 -> 0, U2 -> 10};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None} &&
        Length[Lookup[data, "Exit Vertices and Terminal Costs", {}]] === 2 &&
        NumericQ[Lookup[result, "BoundaryMassResidual", Missing["NotAvailable"]]] &&
        Lookup[result, "BoundaryMassResidual", Infinity] <= 10^-10
    ]
    ,
    True
    ,
    TestID -> "Critical solver: chain with two exits example is feasible and mass-conservative"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        IsCriticalSolution[result]
    ]
    ,
    True
    ,
    TestID -> "Critical solution checker: feasible case satisfies full MFG constraints"
]

Test[
    Module[{data, d2e, result, report},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        report = IsCriticalSolution[result, "ReturnReport" -> True];
        TrueQ[report["Valid"]] &&
        AssociationQ[report["BlockChecks"]] &&
        report["EquationResidualPass"]
    ]
    ,
    True
    ,
    TestID -> "Critical solution checker: report mode returns block checks and residual pass"
]

Test[
    Module[{report},
        report = IsCriticalSolution[
            <|"AssoCritical" -> <||>|>,
            "ReturnReport" -> True
        ];
        !TrueQ[report["Valid"]] &&
        report["Reason"] === "MissingRequiredFields" &&
        ListQ[report["MissingKeys"]] &&
        Length[report["MissingKeys"]] > 0
    ]
    ,
    True
    ,
    TestID -> "Critical solution checker: strict mode rejects incomplete inputs"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[d2e],
            {MFGSystemSolver::nosolution}
        ];
        !IsCriticalSolution[result]
    ]
    ,
    True
    ,
    TestID -> "Critical solution checker: infeasible case is rejected"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[3] /. {I1 -> 2, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            NonLinearSolver[
                d2e,
                "MaxIterations" -> 20,
                "Tolerance" -> 10^-8,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"NonLinear", "Success", "Feasible", None} &&
        AssociationQ[result["AssoNonCritical"]] &&
        result["Solution"] === result["AssoNonCritical"] &&
        AssociationQ[result["FlowAssociation"]] &&
        AssociationQ[result["SignedEdgeFlows"]] &&
        VectorQ[result["ComparableFlowVector"], NumericQ] &&
        NumericQ[result["KirchhoffResidual"]] &&
        NumericQ[result["BoundaryMassResidual"]] &&
        KeyExistsQ[result, "UtilityClassResidual"] &&
        AssociationQ[result["Convergence"]] &&
        NumericQ[result["Convergence"]["SolveTime"]] &&
        IsFeasible[result]
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: nonlinear solver"
]

Test[
    Module[{data, result},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        result = Quiet[
            MonotoneSolverFromData[
                data,
                "ResidualTolerance" -> 10^-6,
                "MaxTime" -> 10,
                "MaxSteps" -> 2000,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, "Solver", None] === "Monotone" &&
        MemberQ[{"Success", "NonConverged"}, Lookup[result, "ResultKind", None]] &&
        If[
            Lookup[result, "ResultKind", None] === "NonConverged",
            Lookup[result, "Message", None] === "ResidualExceedsTolerance",
            Lookup[result, "Message", None] === None
        ] &&
        AssociationQ[result["AssoMonotone"]] &&
        result["Solution"] === result["AssoMonotone"] &&
        AssociationQ[result["FlowAssociation"]] &&
        AssociationQ[result["SignedEdgeFlows"]] &&
        VectorQ[result["ComparableFlowVector"], NumericQ] &&
        NumericQ[result["KirchhoffResidual"]] &&
        KeyExistsQ[result, "BoundaryMassResidual"] &&
        KeyExistsQ[result, "UtilityClassResidual"] &&
        AssociationQ[result["Convergence"]] &&
        NumericQ[result["Convergence"]["FinalResidual"]] &&
        IntegerQ[result["ReducedStateDimension"]] &&
        IntegerQ[result["FullStateDimension"]]
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: monotone solver success"
]

Test[
    Module[{data, result},
        data = <|
            "Vertices List" -> {1},
            "Adjacency Matrix" -> {{0}},
            "Entrance Vertices and Flows" -> {},
            "Exit Vertices and Terminal Costs" -> {},
            "Switching Costs" -> {}
        |>;
        result = Quiet[MonotoneSolverFromData[data]];
        Lookup[result, {"Solver", "ResultKind", "Message"}] ===
            {"Monotone", "Degenerate", "DegenerateCase"} &&
        result["Solution"] === Missing["NotAvailable"] &&
        result["AssoMonotone"] === Missing["NotAvailable"] &&
        result["Convergence"]["StopReason"] === "DegenerateCase"
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: monotone degenerate case"
]

Test[
    Module[{data, result},
        data = GetExampleData[3] /. {I1 -> 0, U1 -> 0};
        result = Quiet[
            MonotoneSolverFromData[
                data,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, {"Solver", "ResultKind", "Message"}] ===
            {"Monotone", "Success", None} &&
        result["Convergence"]["StopReason"] === "QuadraticCriticalSolve" &&
        VectorQ[result["ComparableFlowVector"], NumericQ] &&
        Max[Abs[result["ComparableFlowVector"]]] < 10^-9
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: monotone zero-flow critical solve"
]

Test[
    Module[{data, d2e, b1, k1, vars1, b2, k2, costPlaceholder, vars2},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        {b1, k1, vars1} = GetKirchhoffLinearSystem[d2e];
        {b2, k2, costPlaceholder, vars2} = GetKirchhoffMatrix[d2e];
        b1 === b2 && Normal[k1] === Normal[k2] && vars1 === vars2
    ]
    ,
    True
    ,
    TestID -> "GetKirchhoffLinearSystem matches GetKirchhoffMatrix core outputs"
]

Test[
    Module[{data, d2e, critical, nonlinear, monotone, cVec, nVec, mVec},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        d2e = DataToEquations[data];
        critical = Quiet[CriticalCongestionSolver[d2e]];
        nonlinear = Quiet[
            NonLinearSolver[
                d2e,
                "MaxIterations" -> 20,
                "Tolerance" -> 10^-8,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        monotone = Quiet[
            MonotoneSolverFromData[
                data,
                "ResidualTolerance" -> 10^-6,
                "MaxTime" -> 10,
                "MaxSteps" -> 2000,
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        cVec = critical["ComparableFlowVector"];
        nVec = nonlinear["ComparableFlowVector"];
        mVec = monotone["ComparableFlowVector"];
        critical["ComparableEdges"] === nonlinear["ComparableEdges"] &&
        nonlinear["ComparableEdges"] === monotone["ComparableEdges"] &&
        MemberQ[{"Success", "NonConverged"}, Lookup[monotone, "ResultKind", None]] &&
        If[
            Lookup[monotone, "ResultKind", None] === "NonConverged",
            Lookup[monotone, "Message", None] === "ResidualExceedsTolerance",
            True
        ] &&
        Max[Abs[cVec - nVec]] < 10^-5 &&
        Max[Abs[cVec - mVec]] < 10^-4 &&
        monotone["KirchhoffResidual"] <= 10^-6 + 10^-12
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: comparable signed edge flows across solvers"
]

Test[
    Module[{data, d2e, b, k, jj, seed, reduced},
        data = GetExampleData[7] /. {I1 -> 80, U1 -> 0, U2 -> 10};
        d2e = DataToEquations[data];
        {b, k, jj} = GetKirchhoffLinearSystem[d2e];
        seed = First @ FindInstance[k . jj == b && And @@ Thread[jj >= 10^-8], jj, Reals];
        reduced = MFGraphs`Private`BuildReducedKirchhoffCoordinates[d2e, N[jj /. seed]];
        0 <= reduced["StateDimension"] <= Length[jj] - MatrixRank[Normal[k]] &&
        reduced["FullDimension"] === Length[jj] &&
        Max[Abs[k . reduced["BasePoint"] - b]] < 10^-7 &&
        Dimensions[reduced["BasisMatrix"]] === {Length[jj], reduced["StateDimension"]} &&
        reduced["CostInvisibleDimension"] >= 0
    ]
    ,
    True
    ,
    TestID -> "Monotone reduced state: observable basis stays feasible and bounded by the Kirchhoff nullspace"
]

Test[
    Module[{data, d2e, result},
        (* Case 8 has switching costs, so the direct solver is skipped *)
        data = GetExampleData[8] /. {I1 -> 100, U1 -> 0, U2 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1, S5 -> 1, S6 -> 1};
        d2e = DataToEquations[data];
        (* Disable numeric backend to force symbolic path *)
        d2e = Append[d2e, "CriticalNumericBackendMode" -> False];
        result = Quiet[CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 0.0001]];
        AssociationQ[result] &&
        Lookup[result, "ResultKind"] === "Failure" &&
        Lookup[result, "Feasibility"] === Missing["NotAvailable"] &&
        Lookup[result, "Message"] === "SymbolicTimeout" &&
        TrueQ[Lookup[result, "SymbolicSolverTimedOut", False]]
    ]
    ,
    True
    ,
    TestID -> "CriticalCongestionSolver: symbolic timeout returns structured envelope, not $Failed"
]
