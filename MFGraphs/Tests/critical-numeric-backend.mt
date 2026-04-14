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

Test[
    Module[{data, d2e, symbolic, numeric, drift},
        data = GetExampleData[1] /. {I1 -> 100, U1 -> 0};
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
        drift = Quiet @ Check[
            Module[{deltaVec},
                deltaVec = symbolic["ComparableFlowVector"] - numeric["ComparableFlowVector"];
                If[ListQ[deltaVec] && deltaVec === {},
                    0.,
                    Norm[deltaVec, Infinity]
                ]
            ],
            Infinity
        ];
        Lookup[numeric, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        Lookup[numeric, "NumericBackendStrategy", Missing["NotAvailable"]] === "JFirst" &&
        TrueQ[Lookup[numeric, "JFirstBackendUsed", False]] &&
        Lookup[numeric, "JFirstBackendFallbackReason", Missing["NotAvailable"]] === None &&
        NumericQ[drift] && drift <= 10^-6 &&
        IsCriticalSolution[numeric]
    ],
    True,
    TestID -> "Critical j-first backend: DAG baseline uses j-first strategy with parity"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[27] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>]
            ]
        ];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        Lookup[result, "NumericBackendStrategy", Missing["NotAvailable"]] =!= "JFirst" &&
        !TrueQ[Lookup[result, "JFirstBackendUsed", False]] &&
        Lookup[result, "JFirstBackendFallbackReason", Missing["NotAvailable"]] === "NonDAG" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "Critical j-first backend: non-DAG case skips j-first and falls back cleanly"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[9] /. {I1 -> 100, I2 -> 40, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            Block[
                {
                    MFGraphs`Private`SolveCriticalJFirstUtilities =
                        (Function[{eqs, flowAssoc, uVars, tol},
                            Failure["CriticalJFirstBackend", <|"Reason" -> "InjectedURecoveryFailure"|>]
                        ])
                },
                CriticalCongestionSolver[
                    Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>]
                ]
            ]
        ];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        Lookup[result, "NumericBackendStrategy", Missing["NotAvailable"]] === "Coupled" &&
        !TrueQ[Lookup[result, "JFirstBackendUsed", False]] &&
        Lookup[result, "JFirstBackendFallbackReason", Missing["NotAvailable"]] === "InjectedURecoveryFailure" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "Critical j-first backend: U-recovery failure falls back to coupled numeric"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        And @@ (KeyExistsQ[result, #] & /@ {
            "NumericBackendUsed",
            "NumericBackendFallbackReason",
            "NumericBackendSolveTime",
            "NumericBackendStrategy",
            "JFirstBackendUsed",
            "JFirstBackendFallbackReason"
        })
    ],
    True,
    TestID -> "Critical j-first backend telemetry: additive keys are always present"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[1] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            Block[
                {
                    QuadraticOptimization = (Function[{args___}, $Failed]),
                    ConvexOptimization = (Function[{args___}, $Failed])
                },
                CriticalCongestionSolver[
                    Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>]
                ]
            ]
        ];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        Lookup[result, "NumericBackendStrategy", Missing["NotAvailable"]] === "JFirst" &&
        Lookup[result, "JFirstOptimizerUsed", Missing["NotAvailable"]] === "LinearOptimization" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "LinearOptimization path: J-First fallback after QP optimizers fail"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>]
            ]
        ];
        (* Test that LP-related telemetry keys are boolean and properly initialized *)
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        BooleanQ[Lookup[result, "NumericBackendLPUsed", False]] &&
        BooleanQ[Lookup[result, "NumericBackendLPSimplexRetryUsed", False]] &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "LP telemetry validation: LP-related keys are properly initialized as booleans"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[1] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            Block[
                {
                    QuadraticOptimization = (Function[{___}, $Failed]),
                    ConvexOptimization    = (Function[{___}, $Failed]),
                    LinearOptimization    = (Function[{___}, $Failed])
                },
                CriticalCongestionSolver[
                    Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>]
                ]
            ]
        ];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        Lookup[result, "JFirstBackendFallbackReason", Missing["NotAvailable"]] === "FlowInfeasible" &&
        Lookup[result, "NumericBackendStrategy", Missing["NotAvailable"]] === "Coupled" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "J-First FlowInfeasible: all optimizers fail, fallback to Coupled succeeds"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        KeyExistsQ[result, "JFirstOptimizerUsed"]
    ],
    True,
    TestID -> "Telemetry: JFirstOptimizerUsed key always present"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[CriticalCongestionSolver[d2e]];
        KeyExistsQ[result, "NumericBackendLPUsed"] && KeyExistsQ[result, "NumericBackendLPSimplexRetryUsed"]
    ],
    True,
    TestID -> "Telemetry: LP and Simplex retry telemetry keys always present"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[27] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>],
                "ForceOptimizer" -> "LinearOptimization"
            ]
        ];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        TrueQ[Lookup[result, "NumericBackendLPUsed", False]] &&
        !TrueQ[Lookup[result, "NumericBackendLPSimplexRetryUsed", False]] &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "ForceOptimizer: LinearOptimization path gives correct result"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[27] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>],
                "ForceOptimizer" -> "Simplex"
            ]
        ];
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        TrueQ[Lookup[result, "NumericBackendLPSimplexRetryUsed", False]] &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "ForceOptimizer: Simplex path gives correct result"
]

Test[
    Module[{data, d2e, result},
        data = GetExampleData[1] /. {I1 -> 100, U1 -> 0};
        d2e = DataToEquations[data];
        result = Quiet[
            CriticalCongestionSolver[
                Join[d2e, <|"CriticalNumericBackendMode" -> "Force"|>],
                "ForceOptimizer" -> None
            ]
        ];
        (* Numeric backend is blocked; solver falls through to symbolic engine *)
        Lookup[result, {"ResultKind", "Feasibility"}] === {"Success", "Feasible"} &&
        !TrueQ[Lookup[result, "NumericBackendUsed", True]] &&
        Lookup[result, "NumericBackendStrategy", Missing["NotAvailable"]] === "Symbolic" &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "ForceOptimizer: None bypasses numeric tier, symbolic fallback succeeds"
]
