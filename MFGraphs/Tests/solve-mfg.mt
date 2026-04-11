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
