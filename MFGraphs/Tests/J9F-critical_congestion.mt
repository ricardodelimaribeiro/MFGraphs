(* Wolfram Language Test file *)

Test[
	MFGEquations = D2E[DataG["Jamaratv9"] /. {I1 -> 2, I2 -> 2000, U1 -> 20, U2 -> 100, U3 -> 0}];
    CriticalCongestionSolver2[MFGEquations]
	,
	<|u1 -> 110558/41, u2 -> 140608/41, u3 -> 110558/41, u4 -> 80426/41, u5 -> 140608/41, u6 -> 88658/41, u7 -> 80426/41, u8 -> 88658/41, u9 -> 80426/41, u10 -> 42062/41, u11 -> 88658/41, u12 -> 44940/41, u13 -> 42062/41, u14 -> 44940/41, u15 -> 42062/41, u16 -> 20, u17 -> 44940/41, u18 -> 100, u19 -> 20, u20 -> 0, u21 -> 20, u22 -> 20, u23 -> 100, u24 -> 0, u25 -> 100, u26 -> 100, u27 -> 0, u28 -> 0, u29 -> 110558/41, u30 -> 110558/41, u31 -> 140608/41, u32 -> 140608/41, j1 -> 30050/41, j2 -> 0, j3 -> 0, j4 -> 30132/41, j5 -> 0, j6 -> 51950/41, j7 -> 8232/41, j8 -> 0, j9 -> 0, j10 -> 38364/41, j11 -> 0, j12 -> 43718/41, j13 -> 2878/41, j14 -> 0, j15 -> 0, j16 -> 41242/41, j17 -> 0, j18 -> 40840/41, j19 -> 0, j20 -> 20, j21 -> 0, j22 -> 40422/41, j23 -> 0, j24 -> 100, j25 -> 0, j26 -> 36740/41, j27 -> 0, j28 -> 120, j29 -> 0, j30 -> 2, j31 -> 0, j32 -> 2000|>
    ,
    TestID -> "Jamarat feasible!"
]
