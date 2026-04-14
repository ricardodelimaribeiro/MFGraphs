(* Wolfram Language Test file *)

Test[
	d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 2, I2 -> 2000, U1 -> 20, U2 -> 100, U3 -> 0}];
    result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 300];
    Module[{asso, flowVars},
        asso = result["AssoCritical"];
        flowVars = Select[Normal[asso], StringMatchQ[ToString[First[#]], "j[*"] &];
        And[
            AssociationQ[asso],                    (* Result structure is valid *)
            IsFeasible[result],                     (* Solution is feasible *)
            Length[flowVars] > 0,                   (* At least some flows computed *)
            Lookup[result, "Feasibility"] === "Feasible"  (* Explicit feasibility flag *)
        ]
    ]
	,
	True
	,
    TestID -> "Jamaratv9: feasible with structural and feasibility validation"
]
