(* Wolfram Language Test file *)

Test[
	MFGEquations = D2E[DataG[7] /. {I1 -> 1, U1 -> 0, U2 -> 10}];
    CriticalCongestionSolver2[MFGEquations]
	,
	Null
    ,
    TestID -> "split no solution"
]
