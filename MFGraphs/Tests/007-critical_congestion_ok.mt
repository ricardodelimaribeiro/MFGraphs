(* Wolfram Language Test file *)
Get["ExamplesData.m"];
Test[
	MFGEquations = Data2Equations[DataG[7] /. {I1 -> 11, U1 -> 0, U2 -> 10}];
    CriticalCongestionSolver2[MFGEquations]
	,
	<|u1 -> 43/2, u2 -> 21/2, u3 -> 21/2, u4 -> 0, u5 -> 21/2, u6 -> 10, u7 -> 0, u8 -> 0, u9 -> 10, u10 -> 10, u11 -> 43/2, u12 -> 43/2, j1 -> 0, j2 -> 11, j3 -> 0, j4 -> 21/2, j5 -> 0, j6 -> 1/2, j7 -> 0, j8 -> 21/2, j9 -> 0, j10 -> 1/2, j11 -> 0, j12 -> 11|>
    ,
    TestID -> "split with solution"
]
