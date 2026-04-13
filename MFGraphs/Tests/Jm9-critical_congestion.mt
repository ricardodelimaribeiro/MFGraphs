(* Wolfram Language Test file *)

Test[
    Quiet[
        d2e = DataToEquations[
            GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}
        ];
        result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 120.];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        AssociationQ[Lookup[result, "AssoCritical", Missing["NotAvailable"]]]
    ]
    ,
    True
    ,
    TestID -> "Jamaratv9: baseline case remains feasible with a 120-second symbolic time budget"
]

Test[
    Quiet[
        d2e = DataToEquations[
            GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}
        ];
        result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 0.01];
        Lookup[result, {"ResultKind", "Feasibility", "Message"}] ===
            {"Failure", Missing["NotAvailable"], "SymbolicTimeout"} &&
        TrueQ[Lookup[result, "SymbolicSolverTimedOut", False]]
    ]
    ,
    True
    ,
    TestID -> "Jamaratv9: symbolic timeout is reported as timeout, not infeasibility"
]
