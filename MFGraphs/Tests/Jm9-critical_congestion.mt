(* Wolfram Language Test file *)
(* Jamaratv9: large network solver tests including timeout envelope. *)

Test[
	Quiet[
		d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
		result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 300];
		Lookup[result, "Feasibility"]
	]
	,
	"Feasible"
    ,
    TestID -> "Jamaratv9-CriticalCongestion-Feasible"
]

Test[
	Quiet[
		d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
		(* Disable numeric backend to force symbolic path, then timeout immediately *)
		d2e = Append[d2e, "CriticalNumericBackendMode" -> False];
		result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 0.0001];
		{Lookup[result, "ResultKind"], Lookup[result, "Feasibility"]}
	]
	,
	{"Failure", Missing["NotAvailable"]}
    ,
    TestID -> "Jamaratv9-CriticalCongestion-TimeoutEnvelope"
]
