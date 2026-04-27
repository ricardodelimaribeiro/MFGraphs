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
Fails for non-critical congestion systems where Alpha != 1 on any edge. \
Returns a list of rules when the system is fully determined, or \
<|\"Rules\" -> rules, \"Equations\" -> residual|> when underdetermined.";

isValidSystemSolution::usage =
"isValidSystemSolution[sys, sol] checks whether sol (the output of \
reduceSystem[sys]) satisfies the constraint blocks of sys. Returns True or \
False. With option \"ReturnReport\" -> True, returns a detailed association \
with per-block results. Tolerance for numeric checks is set via \
\"Tolerance\" (default 10^-6). For underdetermined solutions the partial \
rules are checked; blocks that remain symbolic after substitution are \
reported as Indeterminate, not False.";

dnfReduce::usage =
"dnfReduce[xp, sys] simplifies xp && sys by solving equalities, substituting \
their solutions throughout the system, and distributing over disjunctions. \
Returns a DNF expression with all equalities eliminated where possible. \
dnfReduce[xp, sys, elem] is the 3-argument form used internally to process \
one conjunct elem from sys.";

dnfReduceSystem::usage =
"dnfReduceSystem[sys] solves the mfgSystem sys using linear preprocessing \
followed by dnfReduce instead of Reduce. Handles cases where Reduce times \
out by using equality-substitution and disjunction-distribution. Returns a \
list of rules when fully determined, or \
<|\"Rules\" -> rules, \"Equations\" -> residual|> when underdetermined. \
Fails for non-critical congestion systems where Alpha != 1 on any edge.";

reduceSystem::noncritical =
"reduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

dnfReduceSystem::noncritical =
"dnfReduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

Begin["`Private`"];

trackedVarQ::usage =
"trackedVarQ[var] returns True for solver variables tracked during Reduce post-processing.";

mergeRules::usage =
"mergeRules[oldRules, newRules] joins rule lists while keeping the latest rule for each left-hand side.";

normalizeRules::usage =
"normalizeRules[rules] rewrites right-hand sides through the full rule set.";

topLevelEquations::usage =
"topLevelEquations[constraints] extracts top-level Equal expressions from an And expression or single constraint.";

extractLinearRules::usage =
"extractLinearRules[equations, vars, existingRules] solves linear equations for new eliminable rules.";

accumulateLinearRules::usage =
"accumulateLinearRules[baseConstraints, vars, initRules] repeatedly extracts linear rules from the constraint system.";

attachAccumulatedRules::usage =
"attachAccumulatedRules[result, rulesAcc] merges accumulated rules into a reduceSystem result.";

branchToRules::usage =
"branchToRules[branch, allVars] extracts explicit variable rules from one Reduce branch.";

parseReduceResult::usage =
"parseReduceResult[reduced, allVars] converts a Reduce result to rules or a rules-plus-residual association.";

checkBlock::usage =
"checkBlock[expr, tol] evaluates a substituted constraint block using tolerance tol for numeric comparisons.";

reduceSystem[sys_?mfgSystemQ] :=
    Module[{eqEntryIn, eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts, ineqSwitchingByVertex, altFlows, altTrans, altOptCond,
            ruleExitVals, baseConstraints, constraints, allVars, rulesAcc},

        If[!criticalCongestionSystemQ[sys],
            Message[reduceSystem::noncritical];
            Return[
                Failure["reduceSystem", <|
                    "Message" -> "reduceSystem supports only critical congestion systems with Alpha == 1 on every edge."
                |>],
                Module
            ]
        ];

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

dnfReduceSystem[sys_?mfgSystemQ] :=
    Module[{eqEntryIn, eqSplit, eqGather, eqHJ,
            ineqJs, ineqJts, ineqSwitchingByVertex, altFlows, altTrans, altOptCond,
            ruleExitVals, baseConstraints, constraints, allVars, rulesAcc},

        If[!criticalCongestionSystemQ[sys],
            Message[dnfReduceSystem::noncritical];
            Return[
                Failure["dnfReduceSystem", <|
                    "Message" -> "dnfReduceSystem supports only critical congestion systems with Alpha == 1 on every edge."
                |>],
                Module
            ]
        ];

        eqEntryIn    = systemData[sys, "EqEntryIn"];
        eqSplit      = systemData[sys, "EqBalanceSplittingFlows"];
        eqGather     = systemData[sys, "EqBalanceGatheringFlows"];
        eqHJ         = systemData[sys, "EqGeneral"];
        ineqJs       = systemData[sys, "IneqJs"];
        ineqJts      = systemData[sys, "IneqJts"];
        ineqSwitchingByVertex = systemData[sys, "IneqSwitchingByVertex"];
        altFlows     = systemData[sys, "AltFlows"];
        altTrans     = systemData[sys, "AltTransitionFlows"];
        altOptCond   = systemData[sys, "AltOptCond"];
        ruleExitVals = Normal @ systemData[sys, "RuleExitValues"];

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

        With[{reduced = dnfReduce[True, constraints]},
            attachAccumulatedRules[parseReduceResult[reduced, allVars], rulesAcc]
        ]
    ];

trackedVarQ[var_] :=
    MatchQ[var, j[__] | u[__]];

flattenConjuncts::usage = "flattenConjuncts[expr] unpacks expr into a conjunct list for iterative processing.";
substituteSolution::usage = "substituteSolution[rst, sol] substitutes solution rules sol into rst.";

flattenConjuncts[True]  := {}
flattenConjuncts[False] := {False}
flattenConjuncts[x_And] := List @@ x
flattenConjuncts[x_]    := {x}

substituteSolution[rst_?BooleanQ, _] := rst
substituteSolution[rst_, sol_]       := rst /. sol

dnfReduce[_,    False] := False
dnfReduce[False, _]    := False
dnfReduce[xp_,  True]  := xp

dnfReduce[xp_, eq_Equal] := dnfReduce[xp, True, eq]

dnfReduce[xp_, sys_And] :=
    Module[{conjuncts = List @@ sys, xpAcc = xp, i = 1, elem,
            sol, fsol, newxp, rest, newRest},
        While[i <= Length[conjuncts],
            If[xpAcc === False, Return[False, Module]];
            elem = conjuncts[[i]];
            Which[
                elem === True,
                    i++,
                elem === False,
                    Return[False, Module],
                Head[elem] === Or,
                    rest = If[i >= Length[conjuncts], True, And @@ conjuncts[[i + 1 ;;]]];
                    Return[
                        Module[{results = {}, r},
                            Do[
                                r = dnfReduce[xpAcc, rest, b];
                                If[r =!= False, AppendTo[results, r]],
                                {b, List @@ elem}
                            ];
                            Which[
                                results === {}, False,
                                Length[results] === 1, First[results],
                                True, Or @@ results
                            ]
                        ],
                        Module
                    ],
                Head[elem] === Equal,
                    sol  = Quiet[Solve[elem], Solve::svars];
                    fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
                    newxp = substituteSolution[xpAcc, fsol];
                    If[newxp === False, Return[False, Module]];
                    xpAcc = newxp && elem;
                    If[i < Length[conjuncts],
                        newRest = substituteSolution[And @@ conjuncts[[i + 1 ;;]], fsol];
                        conjuncts = Join[conjuncts[[;; i]], flattenConjuncts[newRest]]
                    ];
                    i++,
                True,
                    xpAcc = xpAcc && elem;
                    i++
            ]
        ];
        xpAcc
    ]

dnfReduce[xp_, sys_Or] :=
    Module[{results = {}, r},
        Do[
            r = dnfReduce[xp, b];
            If[r =!= False, AppendTo[results, r]],
            {b, List @@ sys}
        ];
        Which[
            results === {}, False,
            Length[results] === 1, First[results],
            True, Or @@ results
        ]
    ]

dnfReduce[xp_, leq_] := xp && leq

dnfReduce[xp_, rst_, fst_Equal] :=
    Module[{sol, fsol, newxp, newrst},
        sol  = Quiet[Solve[fst], Solve::svars];
        fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
        newxp = substituteSolution[xp, fsol];
        If[newxp === False,
            False,
            newrst = substituteSolution[rst, fsol];
            dnfReduce[newxp && fst, newrst]
        ]
    ]

dnfReduce[xp_, rst_, fst_] := dnfReduce[xp && fst, rst]

(* reduceSystem only handles the critical-congestion structural system. The
   system object carries the scenario Hamiltonian so this check can reject any
   non-1 default Alpha or per-edge EdgeAlpha before Reduce sees nonlinear
   placeholder equations. *)
criticalCongestionSystemQ[sys_?mfgSystemQ] :=
    Module[{halfPairs, hamiltonian, alphaDefault, edgeAlpha, alphaAtEdge},
        halfPairs   = systemData[sys, "HalfPairs"];
        hamiltonian = systemData[sys, "Hamiltonian"];
        If[MissingQ[hamiltonian] || !AssociationQ[hamiltonian], hamiltonian = <||>];
        alphaDefault = Lookup[hamiltonian, "Alpha", 1];
        edgeAlpha    = Lookup[hamiltonian, "EdgeAlpha", <||>];
        If[!ListQ[halfPairs] || !AssociationQ[edgeAlpha], Return[False, Module]];
        alphaAtEdge[edge_List] :=
            Lookup[
                edgeAlpha,
                Key[edge],
                Lookup[edgeAlpha, Key[Reverse[edge]], alphaDefault]
            ];
        alphaDefault === 1 && AllTrue[halfPairs, alphaAtEdge[#] === 1 &]
    ];

mergeRules[oldRules_List, newRules_List] :=
    Reverse @ DeleteDuplicatesBy[Reverse @ Join[oldRules, newRules], First];

normalizeRules[rules_List] :=
    Map[Rule[First[#], ReplaceRepeated[Last[#], rules]] &, rules];

topLevelEquations[constraints_] :=
    Select[
        Switch[Head[constraints], And, List @@ constraints, _, {constraints}],
        MatchQ[#, _Equal] &
    ];

extractLinearRules[equations_List, vars_List, existingRules_List] :=
    Module[{sol, solRules},
        If[equations === {}, Return[{}, Module]];
        sol = Quiet[Solve[equations, vars], Solve::svars];
        If[!ListQ[sol] || Length[sol] === 0, Return[{}, Module]];
        
        solRules = Cases[First[sol],
            r : Rule[v_, rhs_] /; MemberQ[vars, v] && FreeQ[rhs, v] && FreeQ[rhs, C] :> r
        ];
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
            False, Return[<|"Rules" -> {}, "Equations" -> False|>],
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

Options[isValidSystemSolution] = {
    "Tolerance"    -> 10^-6,
    "ReturnReport" -> False
};

isValidSystemSolution[sys_?mfgSystemQ, sol_, OptionsPattern[]] :=
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
