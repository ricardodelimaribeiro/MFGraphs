(* Wolfram Language Test file *)

Test[
	MFGEquations = D2E[DataG[12] /. {I1 -> 1, U1 -> 0}];
    CriticalCongestionSolver2[MFGEquations]
	,
	<|u1 -> 1, u2 -> 1/2, u3 -> 1, u4 -> 1/2, u5 -> 1/2, u6 -> 1/2, u7 -> 1/2, u8 -> 0, u9 -> 1/2, u10 -> 0, u11 -> 0, u12 -> 0, u13 -> 1, u14 -> 1, j1 -> 0, j2 -> 1/2, j3 -> 0, j4 -> 1/2, j5 -> 0, j6 -> 0, j7 -> 0, j8 -> 1/2, j9 -> 0, j10 -> 1/2, j11 -> 0, j12 -> 1, j13 -> 0, j14 -> 1|>
    ,
    TestID -> "attraction problem"
]
