(* Wolfram Language Test file *)
Test[
	MFGEquations = Data2Equations[DataG[7] /. {I1 -> 100, U1 -> 0, U2 -> 0}];
    CriticalCongestionSolver[MFGEquations]["AssoCritical"]
	,
	(*Seems fine!*)
	<|u[en1, 1] -> 150, u[1, en1] -> 150, u[ex3, 3] -> 0, u[ex4, 4] -> 0, u[3, ex3] -> 0, u[4, ex4] -> 0, u[1, 2] -> 50, u[2, 3] -> 0, u[2, 4] -> 0, u[2, 1] -> 150, u[3, 2] -> 50, u[4, 2] -> 50, j[en1, 1] -> 100, j[1, en1] -> 0, j[ex3, 3] -> 0, j[ex4, 4] -> 0, j[3, ex3] -> 50, j[4, ex4] -> 50, j[1, 2] -> 100, j[2, 3] -> 50, j[2, 4] -> 50, j[2, 1] -> 0, j[3, 2] -> 0, j[4, 2] -> 0, j[2, 1, en1] -> 0, j[en1, 1, 2] -> 100, j[1, 2, 3] -> 50, j[1, 2, 4] -> 50, j[3, 2, 1] -> 0, j[3, 2, 4] -> 0, j[4, 2, 1] -> 0, j[4, 2, 3] -> 0, j[2, 3, ex3] -> 50, j[ex3, 3, 2] -> 0, j[2, 4, ex4] -> 50, j[ex4, 4, 2] -> 0|>
    ,
    TestID -> "Case 7: Y-network symmetric exits"
]
