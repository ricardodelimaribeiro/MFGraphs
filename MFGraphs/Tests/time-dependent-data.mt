(* Wolfram Language Test file *)
(* Tests for time-dependent data detection and validation. *)

(* --- IsTimeDependentQ --- *)

Test[
    IsTimeDependentQ[GetExampleData[2] /. {I1 -> 100, U1 -> 0}]
    ,
    False
    ,
    TestID -> "IsTimeDependentQ: static example returns False"
]

Test[
    IsTimeDependentQ[GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0}]
    ,
    False
    ,
    TestID -> "IsTimeDependentQ: Y-network static returns False"
]

Test[
    IsTimeDependentQ[GetTimeDependentExampleData["TD-1"]]
    ,
    True
    ,
    TestID -> "IsTimeDependentQ: TD-1 returns True"
]

Test[
    IsTimeDependentQ[GetTimeDependentExampleData["TD-2"]]
    ,
    True
    ,
    TestID -> "IsTimeDependentQ: TD-2 returns True"
]

Test[
    IsTimeDependentQ[<|"Vertices List" -> {1, 2}|>]
    ,
    False
    ,
    TestID -> "IsTimeDependentQ: no Time Horizon returns False"
]

Test[
    IsTimeDependentQ["not an association"]
    ,
    False
    ,
    TestID -> "IsTimeDependentQ: non-Association returns False"
]

(* --- ValidateTimeDependentData --- *)

Test[
    ValidateTimeDependentData[GetTimeDependentExampleData["TD-1"]]
    ,
    True
    ,
    TestID -> "ValidateTimeDependentData: TD-1 is valid"
]

Test[
    ValidateTimeDependentData[GetTimeDependentExampleData["TD-5"]]
    ,
    True
    ,
    TestID -> "ValidateTimeDependentData: TD-5 is valid"
]

Test[
    StringQ[ValidateTimeDependentData[<|"Vertices List" -> {1, 2}|>]]
    ,
    True
    ,
    TestID -> "ValidateTimeDependentData: missing Time Horizon fails"
]

Test[
    StringQ[ValidateTimeDependentData[<|"Time Horizon" -> -1.0|>]]
    ,
    True
    ,
    TestID -> "ValidateTimeDependentData: negative Time Horizon fails"
]

Test[
    StringQ[ValidateTimeDependentData[<|
        "Time Horizon" -> 1.0, "Time Steps" -> 0|>]]
    ,
    True
    ,
    TestID -> "ValidateTimeDependentData: zero Time Steps fails"
]

(* --- GetTimeDependentExampleData --- *)

Test[
    AssociationQ[GetTimeDependentExampleData["TD-1"]]
    ,
    True
    ,
    TestID -> "GetTimeDependentExampleData: TD-1 returns Association"
]

Test[
    KeyExistsQ[GetTimeDependentExampleData["TD-3"], "Time Dependent Entrance Flows"]
    ,
    True
    ,
    TestID -> "GetTimeDependentExampleData: TD-3 has time-dependent flows"
]

Test[
    Quiet[GetTimeDependentExampleData["NONEXISTENT"]]
    ,
    $Failed
    ,
    TestID -> "GetTimeDependentExampleData: bad key returns $Failed"
]

(* --- Backward compatibility: static solvers ignore time-dependent keys --- *)

Test[
    Module[{data, d2e, result},
        data = GetTimeDependentExampleData["TD-1"];
        d2e = Quiet[DataToEquations[data]];
        result = Quiet[CriticalCongestionSolver[d2e]];
        result["Status"] === "Feasible"
    ]
    ,
    True
    ,
    TestID -> "Backward compat: CriticalCongestionSolver works on TD data"
]
