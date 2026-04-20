(* Wolfram Language Test file *)
Test[
	MFGEquations = DataToEquations[GetExampleData[12] /. {I1 -> 2, U1 -> 0}];
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
        AnyTrue[Keys[asso], MatchQ[#, j[4, ex4]] &]
    ]
	,
	True
    ,
    TestID -> "Case 12: attraction problem (symbolic region)"
]
