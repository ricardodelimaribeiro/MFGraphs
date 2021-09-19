(* Wolfram Language Test file *)

Test[
	MFGEquations = D2E[DataG[12] /. {I1 -> 2, U1 -> 0}];
    CriticalCongestionSolver2[MFGEquations]
	,
	<|u1 -> 2, u2 -> 1, u3 -> 2, u4 -> 1, u5 -> 1, u6 -> 1, u7 -> 1, u8 -> 0, u9 -> 1, u10 -> 0, u11 -> 0, u12 -> 0, u13 -> 2, u14 -> 2, j1 -> 0, j2 -> 1, j3 -> 0, j4 -> 1, j5 -> 0, j6 -> 0, j7 -> 0, j8 -> 1, j9 -> 0, j10 -> 1, j11 -> 0, j12 -> 2, j13 -> 0, j14 -> 2|>
    ,
    TestID -> "attraction problem"
]
