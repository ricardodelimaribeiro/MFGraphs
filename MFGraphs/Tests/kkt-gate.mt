(* Wolfram Language Test file *)
(* Contract tests for ClassifyKKT — the three-way KKT gate for NonLinearSolver results. *)
(*
   Test strategy: synthetic result envelopes are built directly so these tests run
   in milliseconds without calling any solver.  One live integration test (the last)
   exercises the end-to-end path on a small known-feasible case.
*)

If[!MemberQ[$Packages, "MFGraphs`"],
    Get[
        FileNameJoin[
            {
                DirectoryName[$InputFileName],
                "..",
                "MFGraphs.wl"
            }
        ]
    ]
];

(* ------------------------------------------------------------------ *)
(* Helper: build a minimal synthetic NonLinearSolver result envelope   *)
(* ------------------------------------------------------------------ *)
syntheticResult[minFlow_, kirchhoff_, convResidual_] :=
    Module[{sol},
        sol = If[NumericQ[minFlow],
            Association[j[1, 2] -> minFlow, j[2, 3] -> 1.0],
            Missing["NotAvailable"]
        ];
        <|
            "Solver"             -> "NonLinear",
            "ResultKind"         -> "Success",
            "Feasibility"        -> "Feasible",
            "Message"            -> None,
            "Solution"           -> sol,
            "AssoNonCritical"    -> sol,
            "KirchhoffResidual"  -> kirchhoff,
            "Convergence"        -> <|"FinalResidual" -> convResidual|>
        |>
    ];

(* ------------------------------------------------------------------ *)
(* Return shape contract                                               *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{r, gate},
        r    = syntheticResult[1.0, 1.*^-8, 1.*^-10];
        gate = ClassifyKKT[r];
        AssociationQ[gate] &&
        KeyExistsQ[gate, "KKTClass"] &&
        KeyExistsQ[gate, "KKTReason"] &&
        KeyExistsQ[gate, "KKTMetrics"] &&
        KeyExistsQ[gate, "KKTViolations"] &&
        AssociationQ[gate["KKTMetrics"]] &&
        ListQ[gate["KKTViolations"]]
    ]
    ,
    True
    ,
    TestID -> "ClassifyKKT: return shape contract"
]

(* ------------------------------------------------------------------ *)
(* Feasible path                                                       *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{r},
        r = syntheticResult[5.0, 1.*^-8, 1.*^-10];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Feasible"
    ,
    TestID -> "ClassifyKKT: all conditions within Feasible thresholds -> Feasible"
]

Test[
    Module[{r, gate},
        r    = syntheticResult[5.0, 1.*^-8, 1.*^-10];
        gate = ClassifyKKT[r];
        gate["KKTReason"] === None && gate["KKTViolations"] === {}
    ]
    ,
    True
    ,
    TestID -> "ClassifyKKT: Feasible result has None reason and empty violations"
]

(* Mixed numeric/non-numeric flow values must not be partially aggregated.
   PrimalMinFlow is marked Missing -> Borderline/metrics-unavailable path. *)
Test[
    Module[{r, gate},
        r = syntheticResult[5.0, 1.*^-8, 1.*^-10];
        r["AssoNonCritical"] = <|j[1, 2] -> 1.0, j[2, 3] -> x|>;
        r["Solution"] = r["AssoNonCritical"];
        gate = ClassifyKKT[r];
        MissingQ[gate["KKTMetrics"]["PrimalMinFlow"]] &&
        gate["KKTClass"] === "Borderline" &&
        gate["KKTReason"] === "MetricsUnavailable"
    ]
    ,
    True
    ,
    TestID -> "ClassifyKKT: mixed numeric/non-numeric flow values do not yield partial primal min"
]

(* Exact fixed-point: ConvergenceResidual is Missing — should count as Feasible *)
Test[
    Module[{r},
        r = syntheticResult[5.0, 1.*^-8, Missing["NotAvailable"]];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Feasible"
    ,
    TestID -> "ClassifyKKT: missing ConvergenceResidual (exact fixed-point) counts as Feasible"
]

(* ------------------------------------------------------------------ *)
(* Borderline paths                                                    *)
(* ------------------------------------------------------------------ *)

(* Small negative flow: within PrimalBorderlineThreshold *)
Test[
    Module[{r},
        r = syntheticResult[-5.*^-8, 1.*^-8, 1.*^-10];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Borderline"
    ,
    TestID -> "ClassifyKKT: small negative flow within borderline threshold -> Borderline"
]

Test[
    Module[{r},
        r = syntheticResult[-5.*^-8, 1.*^-8, 1.*^-10];
        ClassifyKKT[r]["KKTReason"]
    ]
    ,
    "SmallNegativeFlow"
    ,
    TestID -> "ClassifyKKT: small negative flow reason is SmallNegativeFlow"
]

(* Kirchhoff slightly above Feasible threshold but below Borderline *)
Test[
    Module[{r},
        r = syntheticResult[1.0, 5.*^-5, 1.*^-10];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Borderline"
    ,
    TestID -> "ClassifyKKT: Kirchhoff between Feasible and Borderline thresholds -> Borderline"
]

(* Slow convergence: FinalResidual between Feasible and Borderline thresholds *)
Test[
    Module[{r},
        r = syntheticResult[1.0, 1.*^-8, 5.*^-5];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Borderline"
    ,
    TestID -> "ClassifyKKT: ConvergenceResidual between thresholds -> Borderline"
]

(* ------------------------------------------------------------------ *)
(* Infeasible paths                                                    *)
(* ------------------------------------------------------------------ *)

(* Clearly negative flows (beyond PrimalBorderlineThreshold = -1e-6) *)
Test[
    Module[{r},
        r = syntheticResult[-0.5, 1.*^-8, 1.*^-10];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Infeasible"
    ,
    TestID -> "ClassifyKKT: clearly negative flow -> Infeasible"
]

Test[
    Module[{r},
        r = syntheticResult[-0.5, 1.*^-8, 1.*^-10];
        ClassifyKKT[r]["KKTReason"]
    ]
    ,
    "NegativeFlow"
    ,
    TestID -> "ClassifyKKT: negative flow reason is NegativeFlow"
]

(* Kirchhoff residual above Borderline threshold *)
Test[
    Module[{r},
        r = syntheticResult[1.0, 0.5, 1.*^-10];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Infeasible"
    ,
    TestID -> "ClassifyKKT: KirchhoffResidual above borderline threshold -> Infeasible"
]

Test[
    Module[{r},
        r = syntheticResult[1.0, 0.5, 1.*^-10];
        ClassifyKKT[r]["KKTReason"]
    ]
    ,
    "KirchhoffViolation"
    ,
    TestID -> "ClassifyKKT: Kirchhoff violation reason is KirchhoffViolation"
]

(* Convergence residual above Borderline threshold: not converged *)
Test[
    Module[{r},
        r = syntheticResult[1.0, 1.*^-8, 0.5];
        ClassifyKKT[r]["KKTClass"]
    ]
    ,
    "Infeasible"
    ,
    TestID -> "ClassifyKKT: ConvergenceResidual above borderline threshold -> Infeasible"
]

Test[
    Module[{r},
        r = syntheticResult[1.0, 1.*^-8, 0.5];
        ClassifyKKT[r]["KKTReason"]
    ]
    ,
    "NotConverged"
    ,
    TestID -> "ClassifyKKT: not-converged reason is NotConverged"
]

(* ------------------------------------------------------------------ *)
(* Violations list contract                                            *)
(* ------------------------------------------------------------------ *)

(* Two simultaneous violations are both reported *)
Test[
    Module[{r},
        r = syntheticResult[-5.*^-8, 5.*^-5, 1.*^-10];
        Sort[ClassifyKKT[r]["KKTViolations"]]
    ]
    ,
    Sort[{"Primal", "Kirchhoff"}]
    ,
    TestID -> "ClassifyKKT: simultaneous Primal+Kirchhoff violations both appear in list"
]

(* ------------------------------------------------------------------ *)
(* Custom threshold options                                            *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{r},
        (* Raise the Kirchhoff feasible threshold so 1e-4 passes *)
        r = syntheticResult[1.0, 1.*^-4, 1.*^-10];
        ClassifyKKT[r, "KirchhoffFeasibleThreshold" -> 10^-3]["KKTClass"]
    ]
    ,
    "Feasible"
    ,
    TestID -> "ClassifyKKT: custom KirchhoffFeasibleThreshold widens Feasible band"
]

(* ------------------------------------------------------------------ *)
(* KKTMetrics content                                                  *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{r, m},
        r = syntheticResult[3.7, 2.1*^-9, 4.5*^-11];
        m = ClassifyKKT[r]["KKTMetrics"];
        NumericQ[m["PrimalMinFlow"]] &&
        NumericQ[m["KirchhoffResidual"]] &&
        NumericQ[m["ConvergenceResidual"]] &&
        (* Primal: min flow among {3.7, 1.0} is 1.0 *)
        Abs[m["PrimalMinFlow"] - 1.0] < 10^-10 &&
        Abs[m["KirchhoffResidual"]  - 2.1*^-9]  < 10^-15 &&
        Abs[m["ConvergenceResidual"] - 4.5*^-11] < 10^-20
    ]
    ,
    True
    ,
    TestID -> "ClassifyKKT: KKTMetrics values match what was injected into the result envelope"
]

(* ------------------------------------------------------------------ *)
(* Live integration: small known-feasible NonLinear case              *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{data, result, gate},
        data   = GetExampleData[3] /. {I1 -> 2, U1 -> 0};
        result = Quiet[
            NonLinearSolver[
                data,
                "MaxIterations" -> 20,
                "Tolerance"     -> 10^-8,
                "PotentialFunction"          -> Function[{x, edge}, 0],
                "CongestionExponentFunction" -> Function[edge, 1],
                "InteractionFunction"        -> Function[{m, edge}, -1/m^2]
            ]
        ];
        gate = ClassifyKKT[result];
        gate["KKTClass"] === "Feasible" &&
        gate["KKTReason"] === None &&
        gate["KKTViolations"] === {} &&
        NumericQ[gate["KKTMetrics"]["PrimalMinFlow"]] &&
        gate["KKTMetrics"]["PrimalMinFlow"] >= 0
    ]
    ,
    True
    ,
    TestID -> "ClassifyKKT: live NonLinear run on case 3 classifies as Feasible"
]
