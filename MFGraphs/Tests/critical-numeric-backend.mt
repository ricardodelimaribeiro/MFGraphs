(* Wolfram Language Test file *)
(* Critical numeric backend integration and fallback behavior. *)

Test[
    Module[{data, d2e, symbolic, numeric, drift},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        symbolic = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> False|>]
            ]
        ];
        numeric = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>]
            ]
        ];
        drift = Max[Abs[symbolic["ComparableFlowVector"] - numeric["ComparableFlowVector"]]];
        Lookup[symbolic, {"ResultKind", "Feasibility", "Message"}] ===
            Lookup[numeric, {"ResultKind", "Feasibility", "Message"}] &&
        TrueQ[Lookup[numeric, "NumericBackendUsed", False]] &&
        Lookup[numeric, "NumericBackendFallbackReason", Missing["NotAvailable"]] === None &&
        NumericQ[Lookup[numeric, "NumericBackendSolveTime", Missing["NotAvailable"]]] &&
        NumericQ[drift] && drift <= 10^-6 &&
        IsCriticalSolution[numeric]
    ],
    True,
    TestID -> "Critical numeric backend parity: forced backend matches symbolic result"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[KeyDrop[d2e, "NumericState"], <|"CriticalNumericBackendMode" -> "Force"|>]
            ]
        ];
        !TrueQ[Lookup[result, "NumericBackendUsed", False]] &&
        Lookup[result, "NumericBackendFallbackReason", Missing["NotAvailable"]] === "IneligibleInput" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "Critical numeric backend fallback: ineligible input stays symbolic"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[
                    d2e,
                    <|
                        "CriticalNumericBackendMode" -> "Force",
                        "ForceNumericBackendValidationFailure" -> True
                    |>
                ]
            ]
        ];
        !TrueQ[Lookup[result, "NumericBackendUsed", False]] &&
        Lookup[result, "NumericBackendFallbackReason", Missing["NotAvailable"]] === "ForcedValidationFailure" &&
        Lookup[result, "ResultKind", Missing["NotAvailable"]] === "Success" &&
        Lookup[result, "Feasibility", Missing["NotAvailable"]] === "Feasible" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "Critical numeric backend fallback: forced validation reject demotes numeric candidate"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[
                    d2e,
                    <|
                        "CriticalNumericBackendMode" -> "Force",
                        "CriticalNumericBackendTimeLimit" -> 10^-9
                    |>
                ]
            ]
        ];
        !TrueQ[Lookup[result, "NumericBackendUsed", False]] &&
        Lookup[result, "NumericBackendFallbackReason", Missing["NotAvailable"]] === "NumericBackendTimeout" &&
        Lookup[result, "ResultKind", Missing["NotAvailable"]] === "Success" &&
        Lookup[result, "Feasibility", Missing["NotAvailable"]] === "Feasible" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "Critical numeric backend fallback: timeout fast-fails to symbolic solver"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        And @@ (KeyExistsQ[result, #] & /@ {
            "NumericBackendUsed",
            "NumericBackendFallbackReason",
            "NumericBackendSolveTime"
        })
    ],
    True,
    TestID -> "Critical numeric backend telemetry: keys always present"
]
