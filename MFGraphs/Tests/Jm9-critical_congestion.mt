(* Wolfram Language Test file *)
(* Infeasible case: the solver now short-circuits instead of returning a
   negative-flow association. *)

Test[
	Quiet[
		d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
		result = CriticalCongestionSolver[d2e];
		Lookup[result, "ResultKind"] === "Failure" &&
		Lookup[result, "AssoCritical", Missing["NotAvailable"]] === Null
	]
	,
	True
    ,
    TestID -> "Jamaratv9: infeasible cases short-circuit without a critical-flow association"
]

Test[
	Quiet[
		d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
		result = CriticalCongestionSolver[d2e];
		result["Status"]
	]
	,
	"Infeasible"
    ,
    TestID -> "Jamaratv9: Status is Infeasible"
]
