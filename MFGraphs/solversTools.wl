(* Wolfram Language package *)
(*
   solversTools.wl: symbolic solver for mfgSystem objects.

   Current scope: reduceSystem collects the structural equations,
   flow-balance constraints, non-negativity inequalities, and
   complementarity alternatives from an mfgSystem, seeds the rule set
   with exit-value boundary conditions, iteratively eliminates linear
   equations, then calls Reduce over the Reals on the residual system.
*)

BeginPackage["solversTools`", {"primitives`", "systemTools`"}];

reduceSystem::usage =
"reduceSystem[sys] reduces the structural equations, flow balance, \
non-negativity constraints, and complementarity conditions of the \
mfgSystem sys using Reduce over the Reals. Includes AltOptCond \
(switching-cost optimality complementarity) and \
IneqSwitchingByVertex (switching-cost optimality inequalities). \
Returns a list of rules when the system is fully determined, or \
<|\"Rules\" -> rules, \"Equations\" -> residual|> when underdetermined.";

isReduceSystemSolution::usage =
"isReduceSystemSolution[sys, sol] checks whether sol (the output of \
reduceSystem[sys]) satisfies the constraint blocks of sys. Returns True or \
False. With option \"ReturnReport\" -> True, returns a detailed association \
with per-block results. Tolerance for numeric checks is set via \
\"Tolerance\" (default 10^-6). For underdetermined solutions the partial \
rules are checked; blocks that remain symbolic after substitution are \
reported as Indeterminate, not False.";

Begin["`Private`"];

reduceSystem[sys_?mfgSystemQ] :=
    Module[{eqEntryIn, eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts, ineqSwitchingByVertex, altFlows, altTrans, altOptCond,
            ruleExitVals, baseConstraints, constraints, allVars, rulesAcc},

        eqEntryIn    = systemData[sys, "EqEntryIn"];
        eqSplit      = systemData[sys, "EqBalanceSplittingFlows"];
        eqGather     = systemData[sys, "EqBalanceGatheringFlows"];
        eqHJ         = systemData[sys, "EqGeneral"];
        ineqJs       = systemData[sys, "IneqJs"];
        ineqJts      = systemData[sys, "IneqJts"];
        ineqSwitchingByVertex = systemData[sys, "IneqSwitchingByVertex"];
        altFlows     = systemData[sys, "AltFlows"];
        altTrans     = systemData[sys, "AltTransitionFlows"];
        altOptCond      = systemData[sys, "AltOptCond"];
        ruleExitVals    = Normal @ systemData[sys, "RuleExitValues"];

        baseConstraints = And[
            And @@ eqEntryIn,
            eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts,
            ineqSwitchingByVertex,
            altFlows, altTrans,
            altOptCond
        ];

        allVars = Select[
            Variables[(baseConstraints /. ruleExitVals) /. {Equal -> List, Or -> List, And -> List}],
            trackedVarQ
        ];
        {constraints, rulesAcc} = accumulateLinearRules[baseConstraints, allVars, ruleExitVals];
        allVars = Select[
            Variables[constraints /. {Equal -> List, Or -> List, And -> List}],
            trackedVarQ
        ];

        With[{reduced = Reduce[constraints, allVars, Reals]},
            attachAccumulatedRules[parseReduceResult[reduced, allVars], rulesAcc]
        ]
    ];

trackedVarQ[var_] :=
    MatchQ[var, _[__]] && MemberQ[{"j", "u"}, Quiet @ Check[SymbolName[Head[var]], ""]];

mergeRules[oldRules_List, newRules_List] :=
    Reverse @ DeleteDuplicatesBy[Reverse @ Join[oldRules, newRules], First];

normalizeRules[rules_List] :=
    Module[{normalized = rules, next, iter = 0, maxIter},
        maxIter = Max[10, 2 Length[rules] + 5];
        While[iter < maxIter,
            next = Map[
                Rule[First[#], (Last[#] /. normalized)] &,
                normalized
            ];
            If[next === normalized, Break[]];
            normalized = next;
            iter++
        ];
        normalized
    ];

linearEquationQ[eq_Equal, vars_List] :=
    Module[{expr = Subtract @@ List @@ eq},
        PolynomialQ[expr, vars] && Exponent[expr, vars] <= 1
    ];

topLevelEquations[constraints_] :=
    Select[
        Switch[Head[constraints], And, List @@ constraints, _, {constraints}],
        MatchQ[#, _Equal] &
    ];

extractLinearRules[equations_List, vars_List, existingRules_List] :=
    Module[{sol, solRules},
        If[equations === {}, Return[{}, Module]];
        sol = Check[
            Quiet[Solve[equations, vars, Reals], Solve::svars],
            $Failed
        ];
        If[sol === $Failed || !ListQ[sol] || Length[sol] =!= 1 || !ListQ[First[sol]],
            Return[{}, Module]
        ];
        solRules = Cases[First[sol],
            r : Rule[v_, rhs_] /; MemberQ[vars, v] && FreeQ[rhs, v] :> r
        ];
        solRules = Select[solRules, FreeQ[Last[#], C] &];
        Select[solRules, FreeQ[existingRules[[All, 1]], First[#]] &]
    ];

accumulateLinearRules[baseConstraints_, vars_List, initRules_List] :=
    Module[{rulesAcc = normalizeRules[initRules], constraints, eqs, newRules, iter = 0, maxIter},
        maxIter = Max[10, 2 Length[vars] + 5];
        constraints = baseConstraints /. rulesAcc;
        While[iter < maxIter,
            eqs = topLevelEquations[constraints];
            If[eqs === {}, Break[]];
            newRules = extractLinearRules[eqs, vars, rulesAcc];
            If[newRules === {}, Break[]];
            rulesAcc = normalizeRules[mergeRules[rulesAcc, newRules]];
            constraints = baseConstraints /. rulesAcc;
            iter++
        ];
        {constraints, rulesAcc}
    ];

attachAccumulatedRules[result_, rulesAcc_List] /; ListQ[result] :=
    normalizeRules[mergeRules[rulesAcc, result]];

attachAccumulatedRules[result_, rulesAcc_List] /; AssociationQ[result] :=
    Association[result, "Rules" -> normalizeRules[mergeRules[rulesAcc, Lookup[result, "Rules", {}]]]];

attachAccumulatedRules[result_, _] := result;

branchToRules[branch_, allVars_] :=
    With[{conjuncts = Switch[Head[branch], And, List @@ branch, _, {branch}]},
        Cases[conjuncts,
            HoldPattern[v_ == val_] /;
                MemberQ[allVars, v] && FreeQ[val, Alternatives @@ allVars] :>
                (v -> val)
        ]
    ];

parseReduceResult[reduced_, allVars_] :=
    Module[{conjuncts, rules, residual},
        If[Head[reduced] === Or,
            Return[
                With[{common = Fold[Intersection,
                                    branchToRules[#, allVars] & /@ List @@ reduced]},
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

Options[isReduceSystemSolution] = {
    "Tolerance"    -> 10^-6,
    "ReturnReport" -> False
};

isReduceSystemSolution[sys_?mfgSystemQ, sol_, OptionsPattern[]] :=
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
        ruleExitVals = Normal @ systemData[sys, "RuleExitValues"];

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
            "EqEntryIn"               -> (And @@ systemData[sys, "EqEntryIn"]),
            "EqBalanceSplittingFlows" -> systemData[sys, "EqBalanceSplittingFlows"],
            "EqBalanceGatheringFlows" -> systemData[sys, "EqBalanceGatheringFlows"],
            "EqGeneral"               -> systemData[sys, "EqGeneral"],
            "IneqJs"                  -> systemData[sys, "IneqJs"],
            "IneqJts"                 -> systemData[sys, "IneqJts"],
            "IneqSwitchingByVertex"   -> systemData[sys, "IneqSwitchingByVertex"],
            "AltFlows"                -> systemData[sys, "AltFlows"],
            "AltTransitionFlows"      -> systemData[sys, "AltTransitionFlows"],
            "AltOptCond"              -> systemData[sys, "AltOptCond"]
        |>;

        blockResults = Association @ KeyValueMap[
            #1 -> checkBlock[#2 /. ruleExitVals /. rules, tol] &,
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
checkBlock[True,  _] := True;
checkBlock[False, _] := False;

checkBlock[expr_And, tol_] :=
    With[{results = checkBlock[#, tol] & /@ List @@ expr},
        Which[
            MemberQ[results, False],       False,
            MemberQ[results, Indeterminate], Indeterminate,
            True,                            True
        ]
    ];

checkBlock[expr_Or, tol_] :=
    With[{results = checkBlock[#, tol] & /@ List @@ expr},
        Which[
            MemberQ[results, True],          True,
            MemberQ[results, Indeterminate], Indeterminate,
            True,                            False
        ]
    ];

checkBlock[lhs_ == rhs_, tol_] :=
    With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
        Which[
            d === $Failed,     Indeterminate,
            NumericQ[d],       Abs[d] <= tol,
            True,              With[{s = Quiet @ Check[Simplify[lhs == rhs], $Failed]},
                                   Which[s === True, True, s === False, False, True, Indeterminate]]
        ]
    ];

checkBlock[lhs_ >= rhs_, tol_] :=
    With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d >= -tol, True, Indeterminate]
    ];

checkBlock[lhs_ <= rhs_, tol_] :=
    With[{d = Quiet @ Check[N[rhs - lhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d >= -tol, True, Indeterminate]
    ];

checkBlock[lhs_ > rhs_, tol_] :=
    With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d > tol, True, Indeterminate]
    ];

checkBlock[lhs_ < rhs_, tol_] :=
    With[{d = Quiet @ Check[N[rhs - lhs], $Failed]},
        Which[d === $Failed, Indeterminate, NumericQ[d], d > tol, True, Indeterminate]
    ];

checkBlock[Inequality[a_, Less, b_, Less, c_], tol_] :=
    With[{d1 = Quiet @ Check[N[b - a], $Failed], d2 = Quiet @ Check[N[c - b], $Failed]},
        Which[
            d1 === $Failed || d2 === $Failed, Indeterminate,
            NumericQ[d1] && NumericQ[d2],     d1 > tol && d2 > tol,
            True,                             Indeterminate
        ]
    ];

checkBlock[expr_, tol_] :=
    With[{s = Quiet @ Check[Simplify[expr], $Failed]},
        Which[s === True, True, s === False, False, True, Indeterminate]
    ];

End[];

EndPackage[];
