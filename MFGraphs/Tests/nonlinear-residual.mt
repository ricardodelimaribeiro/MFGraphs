(* Wolfram Language Test file *)
(* Manual diagnostics for non-linear residual gating in NonLinearSolver. *)

Test[
    Module[{data, d2e, result, residual},
        data = GetExampleData[8] /. {
            I1 -> 30, U1 -> 0, U2 -> 0,
            S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1
        };
        d2e = DataToEquations[data];
        result = Quiet @ NonLinearSolver[
            d2e,
            "CongestionExponentFunction" -> Function[edge, 1.1],
            "MaxIterations" -> 15,
            "Tolerance" -> 10^-6
        ];
        residual = Lookup[result, "NonLinearResidual", Missing["NotAvailable"]];
        !IsFeasible[result] &&
        NumericQ[residual] &&
        residual > 10^-3
    ],
    True,
    TestID -> "NonLinear residual trap case is not feasible"
]

Test[
    Module[{data, d2e, result, residual},
        data = GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        result = Quiet @ NonLinearSolver[
            d2e,
            "CongestionExponentFunction" -> Function[edge, 1],
            "MaxIterations" -> 15,
            "Tolerance" -> 10^-6
        ];
        residual = Lookup[result, "NonLinearResidual", Missing["NotAvailable"]];
        IsFeasible[result] &&
        NumericQ[residual] &&
        residual <= 10^-4
    ],
    True,
    TestID -> "Linear baseline remains feasible with tiny residual"
]
