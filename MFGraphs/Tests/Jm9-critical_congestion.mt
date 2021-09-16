(* Wolfram Language Test file *)

Test[
	MFGEquations = D2E[DataG["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
    CriticalCongestionSolver2[MFGEquations]
	,
	Null
    ,
    TestID -> "Jamarat not feasible"
]
