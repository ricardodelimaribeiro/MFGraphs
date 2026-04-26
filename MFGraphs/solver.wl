(* Wolfram Language package *)
(*
   Solver.wl: symbolic solver for mfgSystem objects.

   Current scope: ReduceSystem collects the structural equations,
   flow-balance constraints, non-negativity inequalities, and
   complementarity alternatives from an mfgSystem, substitutes the
   exit-value boundary conditions, then calls Reduce over the Reals.
   Switching-cost inequalities are not included.
*)

BeginPackage["MFGraphs`"];

ReduceSystem::usage =
"ReduceSystem[sys] reduces the structural equations, flow balance, \
non-negativity constraints, and complementarity conditions of the \
mfgSystem sys using Reduce over the Reals. Includes AltOptCond \
(switching-cost optimality complementarity). IneqSwitchingByVertex \
is not included in this pass. \
Returns a list of rules when the system is fully determined, or \
<|\"Rules\" -> rules, \"Equations\" -> residual|> when underdetermined.";

IsReduceSystemSolution::usage =
"IsReduceSystemSolution[sys, sol] checks whether sol (the output of \
ReduceSystem[sys]) satisfies the constraint blocks of sys. Returns True or \
False. With option \"ReturnReport\" -> True, returns a detailed association \
with per-block results. Tolerance for numeric checks is set via \
\"Tolerance\" (default 10^-6). For underdetermined solutions the partial \
rules are checked; blocks that remain symbolic after substitution are \
reported as Indeterminate, not False.";

Begin["`Private`"];

ReduceSystem[sys_?mfgSystemQ] :=
    Module[{eqEntryIn, eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts, altFlows, altTrans, altOptCond, eqExitJunctions,
            ruleExitVals, constraints, allVars},

        eqEntryIn    = SystemData[sys, "EqEntryIn"];
        eqSplit      = SystemData[sys, "EqBalanceSplittingFlows"];
        eqGather     = SystemData[sys, "EqBalanceGatheringFlows"];
        eqHJ         = SystemData[sys, "EqGeneral"];
        ineqJs       = SystemData[sys, "IneqJs"];
        ineqJts      = SystemData[sys, "IneqJts"];
        altFlows     = SystemData[sys, "AltFlows"];
        altTrans     = SystemData[sys, "AltTransitionFlows"];
        altOptCond      = SystemData[sys, "AltOptCond"];
        eqExitJunctions = SystemData[sys, "EqExitJunctions"];
        ruleExitVals    = Normal @ SystemData[sys, "RuleExitValues"];

        constraints = And[
            And @@ eqEntryIn,
            eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts,
            altFlows, altTrans,
            altOptCond,
            eqExitJunctions
        ] /. ruleExitVals;

        allVars = Select[
            Variables[constraints /. {Equal -> List, Or -> List, And -> List}],
            MatchQ[#, _j | _u] &
        ];

        With[{reduced = Reduce[constraints, allVars, Reals]},
            iParseReduceResult[reduced, allVars]
        ]
    ];

iBranchToRules[branch_, allVars_] :=
    With[{conjuncts = Switch[Head[branch], And, List @@ branch, _, {branch}]},
        Cases[conjuncts,
            HoldPattern[v_ == val_] /;
                MemberQ[allVars, v] && FreeQ[val, Alternatives @@ allVars] :>
                (v -> val)
        ]
    ];

iParseReduceResult[reduced_, allVars_] :=
    Module[{conjuncts, rules, residual},
        If[Head[reduced] === Or,
            Return[
                With[{common = Fold[Intersection,
                                    iBranchToRules[#, allVars] & /@ List @@ reduced]},
                    <|"Rules" -> common, "Equations" -> reduced|>
                ]
            ]
        ];
        conjuncts = Switch[Head[reduced],
            And,   List @@ reduced,
            True,  {},
            False, Return[{}],
            _,     {reduced}
        ];
        rules = Cases[conjuncts,
            HoldPattern[v_ == val_] /;
                MemberQ[allVars, v] && FreeQ[val, Alternatives @@ allVars] :>
                (v -> val)
        ];
        residual = DeleteCases[conjuncts,
            HoldPattern[v_ == val_] /;
                MemberQ[allVars, v] && FreeQ[val, Alternatives @@ allVars]
        ];
        If[residual === {},
            rules,
            <|"Rules" -> rules, "Equations" -> And @@ residual|>
        ]
    ];

(* --- Solution validation --- *)

Options[IsReduceSystemSolution] = {
    "Tolerance"    -> 10^-6,
    "ReturnReport" -> False
};

IsReduceSystemSolution[sys_?mfgSystemQ, sol_, OptionsPattern[]] :=
    Module[{tol, returnReportQ, kind, rules, ruleExitVals,
            blocks, blockResults, concretelyFailed, overall, report},
        tol          = OptionValue["Tolerance"];
        returnReportQ = TrueQ[OptionValue["ReturnReport"]];

        kind = Which[
            ListQ[sol] && (sol === {} || MatchQ[sol, {__Rule}]), "Determined",
            AssociationQ[sol] && KeyExistsQ[sol, "Rules"],       "Underdetermined",
            True,                                                  "Unknown"
        ];

        If[kind === "Unknown",
            Return[If[returnReportQ,
                <|"Valid" -> False, "Kind" -> kind,
                  "Reason" -> "UnrecognizedSolutionFormat", "BlockChecks" -> <||>|>,
                False], Module]
        ];

        rules        = If[kind === "Determined", sol, Lookup[sol, "Rules", {}]];
        ruleExitVals = Normal @ SystemData[sys, "RuleExitValues"];

        (* For underdetermined: fail fast if residual equations are False *)
        If[kind === "Underdetermined",
            With[{eqs = Lookup[sol, "Equations", True]},
                If[eqs === False,
                    Return[If[returnReportQ,
                        <|"Valid" -> False, "Kind" -> kind,
                          "Reason" -> "InconsistentResidual", "BlockChecks" -> <||>|>,
                        False], Module]
                ]
            ]
        ];

        blocks = <|
            "EqEntryIn"               -> (And @@ SystemData[sys, "EqEntryIn"]),
            "EqBalanceSplittingFlows" -> SystemData[sys, "EqBalanceSplittingFlows"],
            "EqBalanceGatheringFlows" -> SystemData[sys, "EqBalanceGatheringFlows"],
            "EqGeneral"               -> SystemData[sys, "EqGeneral"],
            "IneqJs"                  -> SystemData[sys, "IneqJs"],
            "IneqJts"                 -> SystemData[sys, "IneqJts"],
            "AltFlows"                -> SystemData[sys, "AltFlows"],
            "AltTransitionFlows"      -> SystemData[sys, "AltTransitionFlows"],
            "AltOptCond"              -> SystemData[sys, "AltOptCond"],
            "EqExitJunctions"         -> SystemData[sys, "EqExitJunctions"]
        |>;

        blockResults = Association @ KeyValueMap[
            #1 -> iCheckBlock[#2 /. ruleExitVals /. rules, tol] &,
            blocks
        ];

        concretelyFailed = Keys @ Select[blockResults, (# === False) &];
        overall = concretelyFailed === {};

        report = <|
            "Valid"                  -> overall,
            "Kind"                   -> kind,
            "Reason"                 -> If[overall, None, "ConstraintViolation"],
            "ConcretelyFailedBlocks" -> concretelyFailed,
            "BlockChecks"            -> blockResults,
            "Tolerance"              -> tol
        |>;

        If[returnReportQ, report, overall]
    ];

(* Recursively evaluate a constraint expression after solution substitution.
   Returns True, False, or Indeterminate. *)
iCheckBlock[True,  _] := True;
iCheckBlock[False, _] := False;

iCheckBlock[expr_And, tol_] :=
    With[{results = iCheckBlock[#, tol] & /@ List @@ expr},
        Which[
            MemberQ[results, False],       False,
            MemberQ[results, Indeterminate], Indeterminate,
            True,                            True
        ]
    ];

iCheckBlock[expr_Or, tol_] :=
    With[{results = iCheckBlock[#, tol] & /@ List @@ expr},
        Which[
            MemberQ[results, True],          True,
            MemberQ[results, Indeterminate], Indeterminate,
            True,                            False
        ]
    ];

iCheckBlock[lhs_ == rhs_, tol_] :=
    With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
        Which[
            d === $Failed,     Indeterminate,
            NumericQ[d],       Abs[d] <= tol,
            True,              With[{s = Quiet @ Check[Simplify[lhs == rhs], $Failed]},
                                   Which[s === True, True, s === False, False, True, Indeterminate]]
        ]
    ];

iCheckBlock[lhs_ >= rhs_, tol_] :=
    With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d >= -tol, True, Indeterminate]
    ];

iCheckBlock[lhs_ <= rhs_, tol_] :=
    With[{d = Quiet @ Check[N[rhs - lhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d >= -tol, True, Indeterminate]
    ];

iCheckBlock[lhs_ > rhs_, tol_] :=
    With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d > tol, True, Indeterminate]
    ];

iCheckBlock[lhs_ < rhs_, tol_] :=
    With[{d = Quiet @ Check[N[rhs - lhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d > tol, True, Indeterminate]
    ];

iCheckBlock[Inequality[a_, Less, b_, Less, c_], tol_] :=
    With[{d1 = Quiet @ Check[N[b - a], $Failed], d2 = Quiet @ Check[N[c - b], $Failed]},
        Which[
            d1 === $Failed || d2 === $Failed, Indeterminate,
            NumericQ[d1] && NumericQ[d2],     d1 > tol && d2 > tol,
            True,                             Indeterminate
        ]
    ];

iCheckBlock[expr_, tol_] :=
    With[{s = Quiet @ Check[Simplify[expr], $Failed]},
        Which[s === True, True, s === False, False, True, Indeterminate]
    ];

End[];

EndPackage[];
