Test[
    Module[{data, d2e, result},
        data = Append[
            MFGraphs`GetExampleData["chain with two exits"] /. {
                MFGraphs`I1 -> 100,
                MFGraphs`U1 -> 0,
                MFGraphs`U2 -> 0
            },
            "Switching Costs" -> {{1, 2, 3, 10}}
        ];
        d2e = Quiet[
            MFGraphs`DataToEquations[data],
            {MFGraphs`DataToEquations::switchingcosts}
        ];
        result = Quiet[
            MFGraphs`CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 60],
            {MFGraphs`DataToEquations::switchingcosts}
        ];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None} &&
        AssociationQ[Lookup[result, "AssoCritical", Missing["NotAvailable"]]] &&
        AssociationQ[Lookup[result, "UnresolvedEquations", Missing["NotAvailable"]]] &&
        IsFeasible[result]
    ]
    ,
    True
    ,
    TestID -> "Critical solver: inconsistent switching on chain-with-two-exits remains feasible"
]

Test[
    Module[{data, d2e},
        data = Append[
            MFGraphs`GetExampleData["chain with two exits"] /. {
                MFGraphs`I1 -> 100,
                MFGraphs`U1 -> 0,
                MFGraphs`U2 -> 0
            },
            "Switching Costs" -> {{1, 2, 3, 10}}
        ];
        d2e = Quiet[
            MFGraphs`DataToEquations[data],
            {MFGraphs`DataToEquations::switchingcosts}
        ];
        Lookup[d2e, "AltTransitionFlows", Missing["NotAvailable"]] =!= Missing["NotAvailable"] &&
        MFGraphs`IsSwitchingCostConsistent[Normal @ Lookup[d2e, "SwitchingCosts", <||>]] === False
    ]
    ,
    True
    ,
    TestID -> "Critical preprocessing: inconsistent switching is recorded while preserving transition constraints"
]
