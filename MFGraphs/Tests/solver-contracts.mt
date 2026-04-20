(* Wolfram Language Test file *)
(* Contract tests for the standardized critical-solver return shape. *)

If[!MemberQ[$Packages, "MFGraphs`"],
    Get[
        FileNameJoin[
            {
                DirectoryName[$InputFileName],
                "..",
                "MFGraphs.wl"
            }
        ]
    ]
];

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
        KeyExistsQ[result, "UnresolvedEquations"] &&
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
        data = MFGraphs`GetExampleData["Jamaratv9"] /. {
            MFGraphs`I1 -> 130,
            MFGraphs`I2 -> 128,
            MFGraphs`U1 -> 20,
            MFGraphs`U2 -> 100,
            MFGraphs`U3 -> 0
        };
        d2e = MFGraphs`DataToEquations[data];
        d2e = Append[d2e, "CriticalNumericBackendMode" -> False];
        result = Quiet[MFGraphs`CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 0.01]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Failure", Missing["NotAvailable"], "SymbolicTimeout"} &&
        TrueQ[Lookup[result, "SymbolicSolverTimedOut", False]] &&
        Lookup[result, "AssoCritical", Missing["NotAvailable"]] === Null &&
        result["Solution"] === Missing["NotAvailable"]
    ]
    ,
    True
    ,
    TestID -> "Critical solver: timed out symbolic path returns timeout envelope instead of infeasibility"
]

Test[
    Module[{data, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0, S1 -> 3, S2 -> 1};
        result = Quiet[CriticalCongestionSolver[data]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None} &&
        AssociationQ[Lookup[result, "UnresolvedEquations", Missing["NotAvailable"]]]
    ]
    ,
    True
    ,
    TestID -> "Critical solver: underdetermined inconsistent-switching is Success/Feasible with symbolic UnresolvedEquations"
]
