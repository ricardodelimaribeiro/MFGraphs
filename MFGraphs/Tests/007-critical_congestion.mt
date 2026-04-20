(* Wolfram Language Test file *)
Test[
	MFGEquations = DataToEquations[GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0}];
    result = CriticalCongestionSolver[
            Join[MFGEquations, <|"CriticalNumericBackendMode" -> False|>]
        ];
    asso = result["AssoCritical"];
    And[
        AssociationQ[asso],
        IsFeasible[result],
        Lookup[result, "Feasibility"] === "Feasible",
        (* Verify that flow variables are present *)
        AnyTrue[Keys[asso], MatchQ[#, j[en1, 1]] &],
        AnyTrue[Keys[asso], MatchQ[#, j[1, 2]] &]
    ]
	,
	True
    ,
    TestID -> "Case 7: Y-network symmetric exits (symbolic region)"
]
