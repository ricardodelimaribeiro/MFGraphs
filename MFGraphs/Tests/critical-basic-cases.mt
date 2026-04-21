(* Wolfram Language Test file *)
(* Consolidated baseline critical-feasibility checks for representative scenarios. *)

Test[
    Module[{cases, runCase},
        cases = {
            <|
                "Data" -> (GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0}),
                "ExpectedFlowKey" -> j[en1, 1],
                "ExpectedFlowKey2" -> j[1, 2]
            |>,
            <|
                "Data" -> (GetExampleData[7] /. {I1 -> 11, U1 -> 0, U2 -> 10}),
                "ExpectedFlowKey" -> j[en1, 1],
                "ExpectedFlowKey2" -> j[1, 2]
            |>,
            <|
                "Data" -> (GetExampleData[12] /. {I1 -> 2, U1 -> 0}),
                "ExpectedFlowKey" -> j[en1, 1],
                "ExpectedFlowKey2" -> j[4, ex4]
            |>
        };
        runCase = Function[{case},
            Module[{eqs, result, asso},
                eqs = DataToEquations[case["Data"]];
                result = CriticalCongestionSolver[
                    Join[eqs, <|"CriticalNumericBackendMode" -> False|>]
                ];
                asso = result["AssoCritical"];
                And[
                    AssociationQ[asso],
                    IsFeasible[result],
                    Lookup[result, "Feasibility"] === "Feasible",
                    AnyTrue[Keys[asso], MatchQ[#, case["ExpectedFlowKey"]] &],
                    AnyTrue[Keys[asso], MatchQ[#, case["ExpectedFlowKey2"]] &]
                ]
            ]
        ];
        And @@ (runCase /@ cases)
    ]
    ,
    True
    ,
    TestID -> "Critical baseline cases: representative symbolic feasibility checks"
]
