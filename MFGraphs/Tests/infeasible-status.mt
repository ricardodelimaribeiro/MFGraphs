(* Wolfram Language Test file *)
(* Tests for the Status key in solver output *)

(* Test: CriticalCongestionSolver returns "Feasible" for a valid case *)
Test[
    Quiet[
        d2e = DataToEquations[GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0}];
        result = CriticalCongestionSolver[d2e];
        result["Status"]
    ]
    ,
    "Feasible"
    ,
    TestID -> "Status: case 7 is Feasible"
]

(* Test: Infeasible case returns "Infeasible" status *)
Test[
    Quiet[
        d2e = DataToEquations[GetExampleData["Jamaratv9"] /. {I1 -> 130, I2 -> 128, U1 -> 20, U2 -> 100, U3 -> 0}];
        result = CriticalCongestionSolver[d2e];
        result["Status"]
    ]
    ,
    "Infeasible"
    ,
    TestID -> "Status: Jamaratv9 infeasible case returns Infeasible"
]

(* Test: IsFeasible helper *)
Test[
    Quiet[
        d2e = DataToEquations[GetExampleData[12] /. {I1 -> 2, U1 -> 0}];
        result = CriticalCongestionSolver[d2e];
        IsFeasible[result]
    ]
    ,
    True
    ,
    TestID -> "IsFeasible: case 12 is feasible"
]

(* Test: MFGSystemSolver short-circuits as soon as the reduced system is infeasible *)
Test[
    Quiet[
        fakeSystem = <|
            "us" -> {},
            "js" -> {},
            "jts" -> {fakeJT},
            "InitRules" -> <||>,
            "NewSystem" -> {True, False, True},
            "costpluscurrents" -> <||>
        |>;
        MFGSystemSolver[fakeSystem][<||>] === Null
    ]
    ,
    True
    ,
    TestID -> "MFGSystemSolver: early exit when inequality block is False"
]
