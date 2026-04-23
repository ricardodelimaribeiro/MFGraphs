(* Wolfram Language Package *)
(*
   SolveMFGDispatch: Unified SolveMFG entrypoint and automatic solver cascade.

   Holds Options[SolveMFG], SolveMFGRunStage (the per-stage dispatcher) and
   SolveMFGAutomaticDispatch (the CriticalCongestion -> Monotone -> NonLinear
   cascade), plus the SolveMFG[input, ...] top-level definitions and small
   envelope helpers.

   Load order dependency (enforced by MFGraphs.wl): all solver submodules
   (Solvers, NonLinearSolver, Monotone) must load first; this file references
   their public symbols and Options[].
*)

Options[SolveMFG] = DeleteDuplicatesBy[
    Join[
        {Method -> "Automatic"},
        Options[CriticalCongestionSolver]
    ],
    First
];

SolveMFGUnknownMethodResult[method_] :=
    <|
        "Solver" -> "SolveMFG",
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotApplicable"],
        "Message" -> "UnknownMethod",
        "Method" -> method,
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGCompiledInputQ[input_Association] :=
    KeyExistsQ[input, "EqGeneral"] &&
    KeyExistsQ[input, "js"] &&
    KeyExistsQ[input, "jts"];

$SolveMFGRequiredEnvelopeKeys = {"Solver", "ResultKind", "Feasibility", "Message"};

SolveMFGClassifyOutcome[result_] :=
    Module[{resultKind, feasibility, message},
        If[!AssociationQ[result], Return["Abort"]];
        If[!And @@ (KeyExistsQ[result, #] & /@ $SolveMFGRequiredEnvelopeKeys), Return["Abort"]];
        resultKind = Lookup[result, "ResultKind", Missing["NotAvailable"]];
        feasibility = Lookup[result, "Feasibility", Missing["NotAvailable"]];
        message = Lookup[result, "Message", Missing["NotAvailable"]];
        Which[
            resultKind === "Success" && feasibility === "Feasible", "Done",
            resultKind === "Degenerate" && message === "DegenerateCase", "Done",
            feasibility === "Infeasible" || MemberQ[{"Failure", "NonConverged"}, resultKind], "Fallback",
            True, "Abort"
        ]
    ];

SolveMFGBuildTraceEntry[method_, result_, decision_] :=
    Module[{solver, resultKind, feasibility, message},
        solver = If[AssociationQ[result], Lookup[result, "Solver", "InvalidSolver"], "InvalidSolver"];
        resultKind = If[AssociationQ[result], Lookup[result, "ResultKind", Missing["NotAvailable"]], Missing["NotAvailable"]];
        feasibility = If[AssociationQ[result], Lookup[result, "Feasibility", Missing["NotAvailable"]], Missing["NotAvailable"]];
        message = If[AssociationQ[result], Lookup[result, "Message", Missing["NotAvailable"]], "InvalidSolverReturnShape"];
        <|
            "Method" -> method,
            "Solver" -> solver,
            "ResultKind" -> resultKind,
            "Feasibility" -> feasibility,
            "Message" -> message,
            "Decision" -> decision
        |>
    ];

SolveMFGFinalizeWithTrace[result_Association, methodUsed_, trace_List] :=
    Join[result, <|"MethodUsed" -> methodUsed, "MethodTrace" -> trace|>];

SolveMFGAbortResult[message_, methodUsed_, trace_List] :=
    <|
        "Solver" -> "SolveMFG",
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotApplicable"],
        "Message" -> message,
        "Method" -> "Automatic",
        "MethodUsed" -> methodUsed,
        "MethodTrace" -> trace,
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGHeldOptionValue[heldOpts_HoldComplete, key_] :=
    Module[{matches},
        matches = Cases[
            heldOpts,
            HoldPattern[Rule[key, value_]] :> HoldComplete[value],
            Infinity
        ];
        If[matches === {}, HoldComplete[Automatic], Last[matches]]
    ];

SolveMFGUnsupportedCriticalAlphaResult[] :=
    <|
        "Solver" -> "CriticalCongestion",
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotAvailable"],
        "Message" -> "Method -> 'CriticalCongestion' only supports CongestionExponentFunction -> 1.",
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGMissingOptionalSolverResult[solverName_String] :=
    <|
        "Solver" -> solverName,
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotAvailable"],
        "Message" -> "SolverNotAvailable",
        "UnavailableSolver" -> solverName,
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGSymbolAvailableQ[sym_Symbol] :=
    DownValues[sym] =!= {} || SubValues[sym] =!= {} || OwnValues[sym] =!= {};

SolveMFGIsCriticalQ[heldOpts_HoldComplete] :=
    Module[{heldAlphaOpt, alphaOpt},
        heldAlphaOpt = SolveMFGHeldOptionValue[heldOpts, "CongestionExponentFunction"];
        If[
            MatchQ[
                heldAlphaOpt,
                HoldComplete[
                    HoldPattern[1 &] | HoldPattern[Function[_, 1]]
                ]
            ],
            Return[True, Module]
        ];
        alphaOpt = ReleaseHold[heldAlphaOpt];
        Which[
            alphaOpt === Automatic, True,
            NumericQ[alphaOpt], alphaOpt == 1,
            AssociationQ[alphaOpt], AllTrue[Values[alphaOpt], # == 1 &],
            True, False
        ]
    ];

SolveMFGCompileInput[input_Association, opts_List:{}] :=
    Module[{potentialFunction, congestionExponentFunction, interactionFunction},
        If[SolveMFGCompiledInputQ[input],
            Return[input, Module]
        ];
        potentialFunction = Replace["PotentialFunction" /. opts, "PotentialFunction" -> Automatic];
        congestionExponentFunction = Replace[
            "CongestionExponentFunction" /. opts,
            "CongestionExponentFunction" -> Automatic
        ];
        interactionFunction = Replace["InteractionFunction" /. opts, "InteractionFunction" -> Automatic];
        WithHamiltonianFunctions[
            potentialFunction,
            congestionExponentFunction,
            interactionFunction,
            DataToEquations[input]
        ]
    ];

(* SolveMFGRunStage — central dispatcher for one solver stage.

   dispatchMode semantics
   ----------------------
   "Direct"    — called from SolveMFG with user-supplied input that may be raw or
                 pre-compiled. The CriticalCongestion branch enforces the alpha guard here.

   "Automatic" — called from SolveMFGAutomaticDispatch, which always pre-compiles input
                 to a d2e before invoking this function. Two invariants hold:
                   • SolveMFGCompiledInputQ[input] is True for every branch, so the
                     SolveMFGCompileInput fallbacks below are dead code in this path.
                   • "CriticalCongestion" only appears in the attempt list when
                     isCritical=True, so the alpha guard is intentionally skipped.
                 The SolveMFGCompileInput branches exist so Direct callers work correctly
                 without a separate pre-compilation step. Do not remove them. *)
SolveMFGRunStage[
    method_String,
    input_Association,
    opts_List,
    heldOpts_HoldComplete,
    dispatchMode_String : "Direct"
] :=
    Module[{d2e, isCritical},
        Switch[method,
            "CriticalCongestion",
                If[
                    SolveMFGCompiledInputQ[input],
                    d2e = input,
                    (* Direct only: guard against alpha != 1 before compiling.
                       Automatic callers guarantee isCritical=True before placing this
                       stage in the attempt list, so the guard is skipped there. *)
                    If[dispatchMode === "Direct",
                        isCritical = SolveMFGIsCriticalQ[heldOpts];
                        If[!TrueQ[isCritical],
                            Return[SolveMFGUnsupportedCriticalAlphaResult[], Module]
                        ]
                    ];
                    d2e = SolveMFGCompileInput[input, opts]
                ];
                CriticalCongestionSolver[
                    d2e,
                    Sequence @@ FilterRules[opts, Options[CriticalCongestionSolver]]
                ]
            ,
            "Monotone",
                (* Direct mode: raw input goes to MonotoneSolverFromData (which calls
                   DataToEquations internally); pre-compiled input goes to MonotoneSolver.
                   These are distinct entry points — MonotoneSolver requires a compiled d2e.
                   Automatic mode: input is always pre-compiled, so MonotoneSolver is used
                   directly and the SolveMFGCompileInput branch is dead code. *)
                If[
                    dispatchMode === "Automatic",
                    If[!SolveMFGSymbolAvailableQ[MonotoneSolver],
                        Return[SolveMFGMissingOptionalSolverResult["MonotoneSolver"], Module]
                    ];
                    d2e = If[SolveMFGCompiledInputQ[input], input, SolveMFGCompileInput[input, opts]];
                    MonotoneSolver[
                        d2e,
                        Sequence @@ FilterRules[opts, Options[MonotoneSolver]]
                    ],
                    If[
                        SolveMFGCompiledInputQ[input],
                        If[!SolveMFGSymbolAvailableQ[MonotoneSolver],
                            Return[SolveMFGMissingOptionalSolverResult["MonotoneSolver"], Module]
                        ];
                        MonotoneSolver[
                            input,
                            Sequence @@ FilterRules[opts, Options[MonotoneSolver]]
                        ],
                        If[!SolveMFGSymbolAvailableQ[MonotoneSolverFromData],
                            Return[SolveMFGMissingOptionalSolverResult["MonotoneSolverFromData"], Module]
                        ];
                        MonotoneSolverFromData[
                            input,
                            Sequence @@ FilterRules[opts, Options[MonotoneSolverFromData]]
                        ]
                    ]
                ]
            ,
            "NonLinear",
                (* SolveMFGCompileInput is a no-op when input is already compiled. *)
                If[!SolveMFGSymbolAvailableQ[NonLinearSolver],
                    Return[SolveMFGMissingOptionalSolverResult["NonLinearSolver"], Module]
                ];
                d2e = If[SolveMFGCompiledInputQ[input], input, SolveMFGCompileInput[input, opts]];
                NonLinearSolver[
                    d2e,
                    Sequence @@ FilterRules[opts, Options[NonLinearSolver]]
                ]
            ,
            _,
                SolveMFGUnknownMethodResult[method]
        ]
    ];

SolveMFGAutomaticDispatch[input_Association, opts_List] :=
    Module[{d2e, result = Missing["NotAvailable"], decision, trace = {}, methodUsed = "Automatic",
            attempts, isCritical, heldOpts},
        heldOpts = HoldComplete[{opts}];
        d2e = If[SolveMFGCompiledInputQ[input], input, Quiet @ Check[SolveMFGCompileInput[input, opts], $Failed]];
        If[!AssociationQ[d2e],
            Return[SolveMFGAbortResult["DataToEquationsFailed", methodUsed, trace], Module]
        ];
        isCritical = SolveMFGIsCriticalQ[heldOpts];
        attempts = Join[
            (* alpha=1: try CriticalCongestion first; alpha!=1: skip it *)
            If[isCritical, {"CriticalCongestion"}, {}],
            {"Monotone", "NonLinear"}
        ];
        Do[
            result = SolveMFGRunStage[attempt, d2e, opts, heldOpts, "Automatic"];
            decision = SolveMFGClassifyOutcome[result];
            methodUsed = attempt;
            trace = Append[trace, SolveMFGBuildTraceEntry[methodUsed, result, decision]];
            Switch[decision,
                "Done",
                    Return[SolveMFGFinalizeWithTrace[result, methodUsed, trace], Module],
                "Fallback",
                    Null,
                _,
                    Return[SolveMFGAbortResult["InvalidSolverReturnShape", methodUsed, trace], Module]
            ],
            {attempt, attempts}
        ];
        If[AssociationQ[result],
            SolveMFGFinalizeWithTrace[result, methodUsed, trace],
            SolveMFGAbortResult["InvalidSolverReturnShape", methodUsed, trace]
        ]
    ];

SolveMFG[input_Association, opts:OptionsPattern[]] :=
    Module[{method, normalizedMethod, heldOpts},
        heldOpts = HoldComplete[{opts}];
        method = OptionValue[Method];
        normalizedMethod = ToString[method];

        Switch[normalizedMethod,
            "Monotone" | "MonotoneSolver",
                SolveMFGRunStage["Monotone", input, {opts}, heldOpts, "Direct"]
            ,
            "CriticalCongestion" | "CriticalCongestionSolver",
                SolveMFGRunStage["CriticalCongestion", input, {opts}, heldOpts, "Direct"]
            ,
            "Automatic",
                SolveMFGAutomaticDispatch[input, {opts}]
            ,
            "NonLinear" | "NonLinearSolver",
                SolveMFGRunStage["NonLinear", input, {opts}, heldOpts, "Direct"]
            ,
            _,
                SolveMFGUnknownMethodResult[method]
        ]
    ];

SolveMFG[_, opts:OptionsPattern[]] :=
    SolveMFGUnknownMethodResult[OptionValue[Method]];
