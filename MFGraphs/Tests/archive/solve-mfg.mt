(* Wolfram Language Test file *)
(* Critical-only routing tests for SolveMFG. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "MFGraphs.wl"}]];

Test[
    Module[{data, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        result = Quiet[SolveMFG[data]];
        Lookup[result, {"Solver", "ResultKind", "Feasibility", "Message"}] ===
            {"CriticalCongestion", "Success", "Feasible", None}
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: Automatic defaults to critical congestion"
]

Test[
    Module[{data, d2e, direct, routed, drift},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        direct = Quiet[CriticalCongestionSolver[d2e]];
        routed = Quiet[SolveMFG[d2e, Method -> "CriticalCongestion"]];
        drift = Max[Abs[direct["ComparableFlowVector"] - routed["ComparableFlowVector"]]];
        Lookup[direct, {"ResultKind", "Feasibility", "Message"}] ===
            Lookup[routed, {"ResultKind", "Feasibility", "Message"}] &&
        NumericQ[drift] && drift <= 10^-8
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: CriticalCongestion method matches direct solver on compiled input"
]

Test[
    Module[{result},
        result = SolveMFG[<||>, Method -> "NotAMethod"];
        Lookup[result, {"Solver", "ResultKind", "Message", "Method"}] ===
            {"SolveMFG", "Failure", "UnknownMethod", "NotAMethod"}
    ]
    ,
    True
    ,
    TestID -> "SolveMFG routing: unknown method returns failure envelope"
]

Test[
    Module[{data, result, trace},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        result = Quiet[SolveMFG[data, Method -> "Automatic"]];
        trace = Lookup[result, "MethodTrace", {}];
        Lookup[result, {"Solver", "MethodUsed", "ResultKind", "Feasibility"}] ===
            {"CriticalCongestion", "CriticalCongestion", "Success", "Feasible"} &&
        Length[trace] === 1 &&
        Lookup[First[trace], "Decision", None] === "Done"
    ]
    ,
    True
    ,
    TestID -> "SolveMFG automatic: trace records single critical stage"
]
