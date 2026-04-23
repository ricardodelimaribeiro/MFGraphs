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
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
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
    Module[{eqs, reduction},
        eqs = <|
            "us" -> {u[1, 2], u[3, 2]},
            "SwitchingCosts" -> <|
                {1, 2, 3} -> 0,
                {3, 2, 1} -> 0
            |>,
            "RuleExitValues" -> <||>,
            "RuleEntryValues" -> <||>
        |>;
        reduction = MFGraphs`Private`BuildUEquivalenceReduction[eqs];
        Lookup[reduction, "EliminatedCount", 0] >= 1
    ],
    True,
    TestID -> "Critical u-equivalence reduction: list-key switching costs eliminate equivalent u vars"
]

Test[
    Module[{data, d2e, patched, result},
        data = GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        patched = Join[d2e, <|"SwitchingCosts" -> AssociationMap[N, d2e["SwitchingCosts"]]|>];
        result = Quiet[
            CriticalCongestionSolver[
                Join[patched, <|"CriticalNumericBackendMode" -> "Force"|>]
            ]
        ];
        Lookup[result, "ResultKind", Missing["NotAvailable"]] === "Success" &&
        TrueQ[Lookup[result, "NumericBackendUsed", False]] &&
        IsCriticalSolution[result]
    ],
    True,
    TestID -> "Critical numeric backend eligibility: machine-zero switching costs keep numeric path enabled"
]

Test[
    Module[{eqs},
        eqs = <|
            "NumericState" -> <|
                "KirchhoffKM" -> SparseArray[{{1, 1} -> 1., {2, 2} -> 1.}, {2, 2}],
                "KirchhoffB" -> {0., 0.},
                "KirchhoffVariables" -> {j[a, b], j[b, a]}
            |>,
            "SwitchingCosts" -> <|{a, b, c} -> 0.|>,
            "auxiliaryGraph" -> Graph[{a, b}, {}],
            "Entrance Vertices and Flows" -> {{"entryA", 1.}},
            "Exit Vertices and Terminal Costs" -> {{"exitA", 0.}}
        |>;
        TrueQ[MFGraphs`Private`CriticalNumericBackendEligibleQ[eqs]]
    ],
    True,
    TestID -> "Critical numeric backend eligibility: symbolic vertex labels are allowed when values are numeric"
]
