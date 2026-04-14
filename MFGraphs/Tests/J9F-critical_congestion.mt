(* Wolfram Language Test file *)

Test[
	d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 2, I2 -> 2000, U1 -> 20, U2 -> 100, U3 -> 0}];
    result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 300];
    AssociationQ[result["AssoCritical"]]
	,
	True
	,
    TestID -> "Jamaratv9: feasible"
]
