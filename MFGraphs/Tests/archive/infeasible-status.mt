(* Wolfram Language Test file *)
(* Tests for the Status key in solver output *)

(* Test: CriticalCongestionSolver returns "Feasible" for a valid case *)
Test[
    Module[{d2e, result},
        Quiet[
            d2e = DataToEquations[GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0}];
            result = CriticalCongestionSolver[d2e];
            result["Status"]
        ]
    ]
    ,
    "Feasible"
    ,
    TestID -> "Status: case 7 is Feasible"
]

(* Test: Jamaratv9 baseline remains Infeasible (no feasible solution exists for this variant) *)
Test[
    Module[{d2e, result},
        Quiet[
            d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
            result = CriticalCongestionSolver[d2e];
            result["Status"]
        ]
    ]
    ,
    "Infeasible"
    ,
    TestID -> "Status: Jamaratv9 baseline case remains Infeasible"
]

(* Test: forced symbolic timeout reports missing status instead of infeasibility *)
Test[
    Module[{d2e, result},
        Quiet[
            d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
            d2e = Append[d2e, "CriticalNumericBackendMode" -> False];
            result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 0.0001];
            result["Status"]
        ]
    ]
    ,
    Missing["NotAvailable"]
    ,
    TestID -> "Status: Jamaratv9 symbolic timeout reports missing status"
]

(* Test: IsFeasible helper *)
Test[
    Module[{d2e, result},
        Quiet[
            d2e = DataToEquations[GetExampleData[12] /. {I1 -> 2, U1 -> 0}];
            result = CriticalCongestionSolver[d2e];
            IsFeasible[result]
        ]
    ]
    ,
    True
    ,
    TestID -> "IsFeasible: case 12 is feasible"
]

(* Test: MFGSystemSolver short-circuits as soon as the reduced system is infeasible *)
Test[
    Module[{fakeSystem},
        Quiet[
            fakeSystem = <|
                "us" -> {},
                "js" -> {},
                "jts" -> {fakeJT},
                "InitRules" -> <||>,
                "NewSystem" -> {True, False, True},
                "costpluscurrents" -> <||>
            |>;
            MFGSystemSolver[fakeSystem][<||>]["Solution"] === Null
        ]
    ]
    ,
    True
    ,
    TestID -> "MFGSystemSolver: early exit when inequality block is False"
]

(* Test: DataToEquations handles inconsistent switching costs gracefully
   (violates triangle inequality) instead of returning $Failed.
   The solver should issue a warning but attempt to find a solution
   using disjunctive constraint strategy. *)
Test[
    Module[{data, d2e, result},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {{0, 1, 1}, {0, 0, 1}, {0, 0, 0}},
            "Entrance Vertices and Flows" -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            (* Inconsistent switching costs: violate triangle inequality
               Cost(1→2→3) = 5+5=10 < Cost(1→3) = 3, violating triangle inequality *)
            "Switching Costs" -> {{1, 2, 3, 5}, {1, 3, 2, 5}}
        |>;
        Quiet[
            d2e = DataToEquations[data];
            (* Check that DataToEquations didn't fail *)
            result = CriticalCongestionSolver[d2e, "SymbolicTimeLimit" -> 60];
            (* Result should be either Feasible or Infeasible, not $Failed *)
            And[
                result =!= $Failed,
                MemberQ[{"Feasible", "Infeasible"},
                    Lookup[result, "Feasibility", Missing[]]]
            ]
        ]
    ]
    ,
    True
    ,
    TestID -> "Inconsistent switching costs: solver handles gracefully without $Failed"
]
