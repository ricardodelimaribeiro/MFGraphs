(* Wolfram Language Test file *)
Test[
	MFGEquations = DataToEquations[GetExampleData[7] /. {I1 -> 11, U1 -> 0, U2 -> 10}];
    Association[
        Normal[CriticalCongestionSolver[
            Join[MFGEquations, <|"CriticalNumericBackendMode" -> False|>]
        ]["AssoCritical"]] /.
            {MFGraphs`Private`u -> u, MFGraphs`Private`j -> j}
    ]
	,
	(*currents result*)
	(*<|u1 -> 43/2, u2 -> 21/2, u3 -> 21/2, u4 -> 0, u5 -> 21/2, u6 -> 10, u7 -> 0, u8 -> 0, u9 -> 10, u10 -> 10, u11 -> 43/2, u12 -> 43/2, j1 -> 0, j2 -> 11, j3 -> 0, j4 -> 21/2, j5 -> 0, j6 -> 1/2, j7 -> 0, j8 -> 21/2, j9 -> 0, j10 -> 1/2, j11 -> 0, j12 -> 11|>*)
	(*vertices indices notation*)
	<|u[en1, 1] -> 43/2, u[1, en1] -> 43/2, u[ex3, 3] -> 0, u[ex4, 4] -> 10, u[3, ex3] -> 0, u[4, ex4] -> 10, u[1, 2] -> 21/2, u[2, 3] -> 0, u[2, 4] -> 10, u[2, 1] -> 43/2, u[3, 2] -> 21/2, u[4, 2] -> 21/2, j[en1, 1] -> 11, j[1, en1] -> 0, j[ex3, 3] -> 0, j[ex4, 4] -> 0, j[3, ex3] -> 21/2, j[4, ex4] -> 1/2, j[1, 2] -> 11, j[2, 3] -> 21/2, j[2, 4] -> 1/2, j[2, 1] -> 0, j[3, 2] -> 0, j[4, 2] -> 0, j[2, 1, en1] -> 0, j[en1, 1, 2] -> 11, j[1, 2, 3] -> 21/2, j[1, 2, 4] -> 1/2, j[3, 2, 1] -> 0, j[3, 2, 4] -> 0, j[4, 2, 1] -> 0, j[4, 2, 3] -> 0, j[2, 3, ex3] -> 21/2, j[ex3, 3, 2] -> 0, j[2, 4, ex4] -> 1/2, j[ex4, 4, 2] -> 0|>
    ,
    TestID -> "Case 7: Y-network asymmetric exits"
]
