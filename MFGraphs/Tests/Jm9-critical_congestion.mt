(* Wolfram Language Test file *)
(* Infeasible case: the solver returns a solution with negative flows *)

Test[
	Quiet[
		d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
		result = CriticalCongestionSolver[d2e]["AssoCritical"];
		Min[Select[Values[result], NumericQ]] < 0
	]
	,
	True
    ,
    TestID -> "Jamaratv9: infeasible (negative flows)"
]
