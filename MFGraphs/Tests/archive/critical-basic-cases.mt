(* Wolfram Language Test file *)
(* Baseline critical-feasibility checks for representative scenarios. *)

Test[
    Module[{data, eqs, result, asso},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        eqs = DataToEquations[data];
        result = CriticalCongestionSolver[
            Join[eqs, <|"CriticalNumericBackendMode" -> False|>]
        ];
        asso = result["AssoCritical"];
        And[
            AssociationQ[asso],
            IsFeasible[result],
            Lookup[result, "Feasibility"] === "Feasible",
            AnyTrue[Keys[asso], MatchQ[#, j[en1, 1]] &],
            AnyTrue[Keys[asso], MatchQ[#, j[1, 2]] &]
        ]
    ]
    ,
    True
    ,
    TestID -> "Critical baseline: case 7 symmetric exits is feasible"
]

Test[
    Module[{data, eqs, result, asso},
        data = GetExampleData[7] /. {I1 -> 11, U1 -> 0, U2 -> 10};
        eqs = DataToEquations[data];
        result = CriticalCongestionSolver[
            Join[eqs, <|"CriticalNumericBackendMode" -> False|>]
        ];
        asso = result["AssoCritical"];
        And[
            AssociationQ[asso],
            IsFeasible[result],
            Lookup[result, "Feasibility"] === "Feasible",
            AnyTrue[Keys[asso], MatchQ[#, j[en1, 1]] &],
            AnyTrue[Keys[asso], MatchQ[#, j[1, 2]] &]
        ]
    ]
    ,
    True
    ,
    TestID -> "Critical baseline: case 7 asymmetric exits is feasible"
]

Test[
    Module[{data, eqs, result, asso},
        data = GetExampleData[12] /. {I1 -> 2, U1 -> 0};
        eqs = DataToEquations[data];
        result = CriticalCongestionSolver[
            Join[eqs, <|"CriticalNumericBackendMode" -> False|>]
        ];
        asso = result["AssoCritical"];
        And[
            AssociationQ[asso],
            IsFeasible[result],
            Lookup[result, "Feasibility"] === "Feasible",
            AnyTrue[Keys[asso], MatchQ[#, j[en1, 1]] &],
            AnyTrue[Keys[asso], MatchQ[#, j[4, ex4]] &]
        ]
    ]
    ,
    True
    ,
    TestID -> "Critical baseline: case 12 attraction is feasible"
]
