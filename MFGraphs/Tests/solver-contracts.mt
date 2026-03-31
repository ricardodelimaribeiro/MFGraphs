(* Wolfram Language Test file *)
(* Contract tests for the opt-in standardized solver return shape. *)

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e, "ReturnShape" -> "Standard"]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None} &&
        AssociationQ[result["AssoCritical"]] &&
        result["Solution"] === result["AssoCritical"] &&
        AssociationQ[result["FlowAssociation"]] &&
        AssociationQ[result["SignedEdgeFlows"]] &&
        VectorQ[result["ComparableFlowVector"], NumericQ] &&
        NumericQ[result["KirchhoffResidual"]] &&
        IsFeasible[result]
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: critical congestion solver"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[3] /. {I1 -> 2, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            NonLinearSolver[
                d2e,
                "MaxIterations" -> 2,
                "ReturnShape" -> "Standard",
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
                "TimeSteps" -> 20,
                "ReturnShape" -> "Standard",
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, {"Solver", "ResultKind", "Message"}] ===
            {"Monotone", "Success", None} &&
        AssociationQ[result["AssoMonotone"]] &&
        result["Solution"] === result["AssoMonotone"] &&
        AssociationQ[result["FlowAssociation"]] &&
        AssociationQ[result["SignedEdgeFlows"]] &&
        VectorQ[result["ComparableFlowVector"], NumericQ] &&
        NumericQ[result["KirchhoffResidual"]] &&
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
        result = Quiet[MonotoneSolverFromData[data, "ReturnShape" -> "Standard"]];
        Lookup[result, {"Solver", "ResultKind", "Message"}] ===
            {"Monotone", "Degenerate", "DegenerateCase"} &&
        result["Solution"] === Missing["NotAvailable"] &&
        result["AssoMonotone"] === Missing["NotAvailable"]
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
                "ReturnShape" -> "Standard",
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        Lookup[result, {"Solver", "ResultKind", "Message"}] ===
            {"Monotone", "Failure", "SeedFindInstanceFailed"} &&
        result["Solution"] === Missing["NotAvailable"] &&
        result["AssoMonotone"] === Missing["NotAvailable"]
    ]
    ,
    True
    ,
    TestID -> "Standard return shape: monotone seed failure"
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
        critical = Quiet[CriticalCongestionSolver[d2e, "ReturnShape" -> "Standard"]];
        nonlinear = Quiet[
            NonLinearSolver[
                d2e,
                "MaxIterations" -> 2,
                "ReturnShape" -> "Standard",
                "PotentialFunction" -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction" -> Function[{m, edge}, -1/m^2]
            ]
        ];
        monotone = Quiet[
            MonotoneSolverFromData[
                data,
                "TimeSteps" -> 20,
                "ReturnShape" -> "Standard",
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
        Max[Abs[cVec - nVec]] < 10^-5 &&
        Max[Abs[cVec - mVec]] < 10^-4 &&
        monotone["KirchhoffResidual"] < 10^-6
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
        reduced["StateDimension"] === Length[jj] - MatrixRank[Normal[k]] &&
        reduced["FullDimension"] === Length[jj] &&
        Max[Abs[k . reduced["BasePoint"] - b]] < 10^-7 &&
        Dimensions[reduced["BasisMatrix"]] === {Length[jj], reduced["StateDimension"]}
    ]
    ,
    True
    ,
    TestID -> "Monotone reduced state: affine basis matches Kirchhoff dimension"
]
