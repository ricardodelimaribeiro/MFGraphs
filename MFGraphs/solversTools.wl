(* Wolfram Language package *)
(*
   solversTools.wl: symbolic solvers for mfgSystem objects.

   Three solvers share a common preprocessing pipeline (buildSolverInputs):
   collect structural equations, apply exit-value rules, eliminate linear
   equations via accumulateLinearRules. They differ only in the final step:
     reduceSystem         -- Reduce[constraints, allVars, Reals]
     dnfReduceSystem      -- dnfReduce[True, constraints]
     optimizedDNFReduceSystem
                          -- branch-state exact DNF reduction
     activeSetReduceSystem
                          -- branch-state exact active-set reduction
     booleanReduceSystem  -- BooleanConvert to DNF, Reduce per disjunct
     findInstanceSystem   -- FindInstance[constraints, allVars, Reals]
*)

BeginPackage["solversTools`", {"primitives`", "utilities`", "systemTools`"}];

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

computeKirchhoffResidual::usage =
"computeKirchhoffResidual[sys, sol] calculates the Kirchhoff Residual (the maximum absolute divergence \
in the transition flux conservation equations) for the solution sol. \
Returns the maximum numerical residual, or Missing[\"NotAvailable\"] if the \
solution is symbolic or incomplete.";

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

optimizedDNFReduceSystem::usage =
"optimizedDNFReduceSystem[sys] is an opt-in exact solver for critical \
congestion systems. It follows the same preprocessing and output contract as \
dnfReduceSystem. It carries small DNF branch families directly as rules plus \
residual constraints, and falls back to the proven exact DNF reducer for \
larger residual variable sets. \
Returns a list of rules when fully determined, or \
<|\"Rules\" -> rules, \"Equations\" -> residual|> when underdetermined. \
Fails for non-critical congestion systems where Alpha != 1 on any edge.";

activeSetReduceSystem::usage =
"activeSetReduceSystem[sys] is an opt-in exact active-set solver for the \
critical-congestion linear complementarity structure. It enumerates small \
complementarity alternatives incrementally with exact linear substitution and \
falls back to the proven exact DNF reducer for larger residual variable sets. \
Returns the same rule/residual shape as dnfReduceSystem. Fails for \
non-critical congestion systems where Alpha != 1 on any edge.";

reduceSystem::noncritical =
"reduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

dnfReduceSystem::noncritical =
"dnfReduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

optimizedDNFReduceSystem::noncritical =
"optimizedDNFReduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

activeSetReduceSystem::noncritical =
"activeSetReduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

booleanReduceSystem::usage =
"booleanReduceSystem[sys] solves the mfgSystem sys by converting the \
preprocessed constraint system to DNF via BooleanConvert, then calling \
Reduce independently on each disjunct. Each disjunct is a pure conjunction \
(no Or), so Reduce avoids case-splitting. Non-False results are collected; \
if the system has a unique equilibrium all non-False results are equivalent. \
Returns a list of rules when fully determined, or \
<|\"Rules\" -> rules, \"Equations\" -> residual|> when underdetermined. \
Fails for non-critical congestion systems where Alpha != 1 on any edge. \
Options: \"DisjunctTimeout\" (default 30s per Reduce call), \
\"ReturnAll\" (default False; True returns all non-False parsed results).";

booleanReduceSystem::noncritical =
"booleanReduceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

booleanReduceSystem::multisol =
"booleanReduceSystem found `1` non-False disjuncts with differing rules. \
Returning first; use \"ReturnAll\" -> True to inspect all results.";

findInstanceSystem::usage =
"findInstanceSystem[sys] solves the mfgSystem sys by collecting and \
linearly preprocessing constraints, then calling FindInstance over the \
remaining variables. Returns one feasible list of rules. If no instance is \
found or the final solve times out, returns \
<|\"Rules\" -> accumulatedRules, \"Equations\" -> False|>. \
Fails for non-critical congestion systems where Alpha != 1 on any edge. \
Options: \"Timeout\" (default Infinity).";

findInstanceSystem::noncritical =
"findInstanceSystem supports only critical congestion systems with Alpha == 1 on every edge.";

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

parseDNFReduceResult::usage =
"parseDNFReduceResult[reduced, allVars] converts dnfReduce output to rules \
or a rules-plus-residual association, harvesting equality solutions from \
surviving DNF branches.";

branchStateReduceResult::usage =
"branchStateReduceResult[constraints, allVars] reduces constraints by carrying \
surviving branches directly as exact rules plus residual constraints.";

harvestDNFBranch::usage =
"harvestDNFBranch[branch, allVars] solves explicit branch equalities, \
substitutes them into residual constraints, and returns a branch association \
or False when the branch is contradictory.";

buildSolverInputs::usage =
"buildSolverInputs[sys] collects all constraint blocks from sys, applies \
exit-value rules, and runs accumulateLinearRules. Returns \
{constraints, allVars, rulesAcc} ready for the final solve step.";

checkBlock::usage =
"checkBlock[expr, tol] evaluates a substituted constraint block using tolerance tol for numeric comparisons.";

solutionResultKind::usage =
"solutionResultKind[sol] classifies solver output by its residual content.";

dnfResidualDiagnostics::usage =
"dnfResidualDiagnostics[sol] returns lightweight branch and tracked-variable diagnostics for solver residuals.";

dnfReduceDiagnosticReport::usage =
"dnfReduceDiagnosticReport[sys, opts] runs the private DNF reducer with \
instrumentation and returns a diagnostic association. This helper is intended \
for scripts and tests; dnfReduceSystem behavior is unchanged.";

reduceSystem[sys_?mfgSystemQ] :=
    withCriticalCongestionSolver[sys, "reduceSystem",
        Function[{constraints, allVars},
            parseReduceResult[Reduce[constraints, allVars, Reals], allVars]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

dnfReduceSystem[sys_?mfgSystemQ] :=
    withCriticalCongestionSolver[sys, "dnfReduceSystem",
        Function[{constraints, allVars},
            parseDNFReduceResult[dnfReduce[True, constraints], allVars]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

optimizedDNFReduceSystem[sys_?mfgSystemQ] :=
    withCriticalCongestionSolver[sys, "optimizedDNFReduceSystem",
        Function[{constraints, allVars},
            If[Length[allVars] <= 8,
                branchStateReduceResult[constraints, allVars],
                parseDNFReduceResult[dnfReduce[True, constraints], allVars]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

activeSetReduceSystem[sys_?mfgSystemQ] :=
    Module[{ruleExitVals, deterministicConstraints, activeConstraints, allVars,
            detReduced, rulesAcc, activeReduced, result, branches},
        If[!criticalCongestionSystemQ[sys],
            Message[activeSetReduceSystem::noncritical];
            Return[Failure["activeSetReduceSystem", <|"Message" -> "activeSetReduceSystem supports only critical congestion systems with Alpha == 1 on every edge."|>], Module]
        ];
        ruleExitVals = Normal @ systemData[sys, "RuleExitValues"];
        deterministicConstraints = And[
            And @@ systemData[sys, "EqEntryIn"],
            systemData[sys, "EqBalanceSplittingFlows"],
            systemData[sys, "EqBalanceGatheringFlows"],
            systemData[sys, "EqGeneral"],
            systemData[sys, "IneqJs"],
            systemData[sys, "IneqJts"],
            systemData[sys, "IneqSwitchingByVertex"]
        ];
        activeConstraints = And[
            systemData[sys, "AltFlows"],
            systemData[sys, "AltTransitionFlows"],
            systemData[sys, "AltOptCond"]
        ];
        allVars = Select[
            Variables[(And[deterministicConstraints, activeConstraints] /. ruleExitVals) /. {Equal -> List, Or -> List, And -> List}],
            trackedVarQ
        ];
        {detReduced, rulesAcc} = accumulateLinearRules[deterministicConstraints, allVars, ruleExitVals];
        activeReduced = activeConstraints /. rulesAcc;
        allVars = Select[
            Variables[And[detReduced, activeReduced] /. {Equal -> List, Or -> List, And -> List}],
            trackedVarQ
        ];
        result = If[Length[allVars] <= 8,
            branches = branchStateReduce[activeReduced, allVars];
            branches = branchStateReduceFromBranches[branches, detReduced, allVars];
            branchStateFinalizeResult[branches, allVars],
            parseDNFReduceResult[dnfReduce[True, And[detReduced, activeReduced]], allVars]
        ];
        attachAccumulatedRules[result, rulesAcc]
    ];

Options[booleanReduceSystem] = {
    "DisjunctTimeout" -> 30,
    "ReturnAll"       -> False
};

booleanReduceSystem[sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    withCriticalCongestionSolver[sys, "booleanReduceSystem",
        Function[{constraints, allVars},
            Module[{dnf, disjuncts, reducedList, nonFalse, parsed, timeout, returnAll},
                timeout   = OptionValue[booleanReduceSystem, {opts}, "DisjunctTimeout"];
                returnAll = TrueQ[OptionValue[booleanReduceSystem, {opts}, "ReturnAll"]];
                dnf       = BooleanConvert[constraints, "DNF"];
                disjuncts = If[Head[dnf] === Or, List @@ dnf, {dnf}];
                reducedList = TimeConstrained[Reduce[#, allVars, Reals], timeout, $Aborted] & /@ disjuncts;
                nonFalse = Select[reducedList, # =!= False && # =!= $Aborted &];
                If[nonFalse === {},
                    Return[<|"Rules" -> {}, "Equations" -> False|>, Module]
                ];
                parsed = parseReduceResult[#, allVars] & /@ nonFalse;
                If[Length[parsed] > 1,
                    With[{ruleSets = (If[ListQ[#], #, Lookup[#, "Rules", {}]] & /@ parsed)},
                        If[!SameQ @@ (Sort /@ ruleSets), Message[booleanReduceSystem::multisol, Length[parsed]]]
                    ]
                ];
                If[returnAll, parsed, First[parsed]]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

Options[findInstanceSystem] = {
    "Timeout" -> Infinity
};

findInstanceSystem[sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    withCriticalCongestionSolver[sys, "findInstanceSystem",
        Function[{constraints, allVars},
            Module[{timeout, instance},
                timeout = OptionValue[findInstanceSystem, {opts}, "Timeout"];
                If[allVars === {},
                    Return[
                        If[TrueQ[Simplify[constraints]],
                            {},
                            <|"Rules" -> {}, "Equations" -> False|>
                        ],
                        Module
                    ]
                ];
                instance = TimeConstrained[
                    FindInstance[constraints, allVars, Reals],
                    timeout,
                    $Aborted
                ];
                If[MatchQ[instance, {{___Rule}, ___}],
                    First[instance],
                    <|"Rules" -> {}, "Equations" -> False|>
                ]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

trackedVarQ[var_] :=
    MatchQ[var, j[__] | u[__]];

residualTrackedVariables[residual_] :=
    DeleteDuplicates @ Cases[Unevaluated[residual], j[__] | u[__], {0, Infinity}];

dnfResidualDiagnostics[sol_] :=
    Module[{residual, topLevelBranches, nestedOrQ, trackedVars},
        residual = If[AssociationQ[sol], Lookup[sol, "Equations", True], True];
        topLevelBranches = If[Head[residual] === Or, Length[List @@ residual], 0];
        nestedOrQ = If[Head[residual] === Or,
            !FreeQ[List @@ residual, Or],
            !FreeQ[residual, Or]
        ];
        trackedVars = residualTrackedVariables[residual];
        <|
            "TopLevelBranchCount" -> topLevelBranches,
            "NestedOrQ" -> nestedOrQ,
            "TrackedVariablesQ" -> (trackedVars =!= {}),
            "TrackedVariables" -> trackedVars
        |>
    ];

solutionResultKind[sol_] :=
    Module[{residual, diagnostics},
        Which[
            sol === $TimedOut,
                "TIMEOUT",
            sol === {},
                "NoSolution",
            ListQ[sol] && AllTrue[sol, MatchQ[#, _Rule] &],
                "Rules",
            AssociationQ[sol] && KeyExistsQ[sol, "Rules"],
                residual = Lookup[sol, "Equations", True];
                Which[
                    residual === True, "Rules",
                    residual === False, "NoSolution",
                    !FreeQ[residual, Or], "Branched",
                    True,
                        diagnostics = dnfResidualDiagnostics[sol];
                        If[TrueQ[diagnostics["TrackedVariablesQ"]], "Parametric", "Residual"]
                ],
            True,
                "Unknown"
        ]
    ];

buildSolverInputs[sys_?mfgSystemQ] :=
    Module[{data, ruleExitVals, baseConstraints, allVars, constraints, rulesAcc},
        data = systemDataFlatten[sys];
        ruleExitVals = Normal @ Lookup[data, "RuleExitValues"];
        (* Zero-cost symmetric triple pairs force u-equalities by antisymmetry of >=.
           Injecting them here applies the substitution to all constraints uniformly
           before the DNF stage, rather than only partially at construction time. *)
        ruleExitVals = Join[ruleExitVals, Lookup[data, "ZeroSwitchUEqualities", {}]];
        baseConstraints = And[
            And @@ Lookup[data, "EqEntryIn"],
            Lookup[data, "EqBalanceSplittingFlows"],
            Lookup[data, "EqBalanceGatheringFlows"],
            Lookup[data, "EqGeneral"],
            Lookup[data, "IneqJs"],
            Lookup[data, "IneqJts"],
            Lookup[data, "IneqSwitchingByVertex"],
            Lookup[data, "AltFlows"],
            Lookup[data, "AltTransitionFlows"],
            Lookup[data, "AltOptCond"]
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
        {constraints, allVars, rulesAcc}
    ];

flattenConjuncts::usage = "flattenConjuncts[expr] unpacks expr into a conjunct list for iterative processing.";
substituteSolution::usage = "substituteSolution[rst, sol] substitutes solution rules sol into rst.";

flattenConjuncts[True]  := {}
flattenConjuncts[False] := {False}
flattenConjuncts[x_And] := List @@ x
flattenConjuncts[x_]    := {x}

substituteSolution[rst_?BooleanQ, _] := rst
substituteSolution[rst_, sol_]       := rst /. sol

constraintPriority[expr_] :=
    Which[
        MatchQ[expr, _Equal], 0,
        FreeQ[expr, Or],      1,
        True,                 2
    ];

orderedConjuncts[expr_] :=
    SortBy[flattenConjuncts[expr], {constraintPriority, ToString[#, InputForm] &}];

orderedAlternatives[expr_Or] := SortBy[List @@ expr, ToString[#, InputForm] &];
orderedAlternatives[expr_]   := {expr};

dnfTopLevelConjuncts[expr_] := flattenConjuncts[expr];

dnfVarVertices[j[a_, b_]] := {a, b};
dnfVarVertices[j[a_, b_, c_]] := {a, b, c};
dnfVarVertices[u[a_, b_]] := {a, b};
dnfVarVertices[_] := {};

dnfExpressionTrackedVariables[expr_] :=
    DeleteDuplicates @ Cases[Unevaluated[expr], j[__] | u[__], {0, Infinity}];

dnfExpressionVertices[expr_] :=
    DeleteDuplicates @ Cases[Flatten[dnfVarVertices /@ dnfExpressionTrackedVariables[expr]], _Integer];

dnfTrackedVariableKind[var_] :=
    Which[
        MatchQ[var, j[_, _, _]], "TransitionFlow",
        MatchQ[var, j[_, _]],    "EdgeFlow",
        MatchQ[var, u[_, _]],    "Value",
        True,                    "Other"
    ];

dnfPathInfo[sys_?mfgSystemQ] :=
    Module[{graph, entries, exits, entryVertices, exitVertices, paths, pathVertices,
            distances, exitDistances, entryDistances, vertices},
        graph = systemData[sys, "Graph"];
        entries = systemData[sys, "Entries"];
        exits = systemData[sys, "Exits"];
        If[!MatchQ[graph, _Graph] || !ListQ[entries] || !ListQ[exits],
            Return[<|"PathVertices" -> {}, "Distances" -> <||>, "ExitDistances" -> <||>, "EntryDistances" -> <||>|>, Module]
        ];
        entryVertices = Cases[entries, {v_, _} :> v];
        exitVertices = Cases[exits, {v_, _} :> v];
        paths = DeleteCases[
            Flatten[Table[FindShortestPath[graph, en, ex], {en, entryVertices}, {ex, exitVertices}], 1],
            {} | _FindShortestPath
        ];
        pathVertices = DeleteDuplicates @ Flatten[paths];
        vertices = VertexList[graph];
        distances = Association @ Table[
            v -> If[pathVertices === {}, Infinity, Min[GraphDistance[graph, v, #] & /@ pathVertices]],
            {v, vertices}
        ];
        exitDistances = Association @ Table[
            v -> If[exitVertices === {}, Infinity, Min[GraphDistance[graph, v, #] & /@ exitVertices]],
            {v, vertices}
        ];
        entryDistances = Association @ Table[
            v -> If[entryVertices === {}, Infinity, Min[GraphDistance[graph, v, #] & /@ entryVertices]],
            {v, vertices}
        ];
        <|
            "PathVertices" -> pathVertices,
            "Distances" -> distances,
            "ExitDistances" -> exitDistances,
            "EntryDistances" -> entryDistances
        |>
    ];

dnfPathRank[expr_, pathInfo_Association] :=
    Module[{verts, pathVerts, distances, exitDistances, entryDistances, onPath,
            minDistance, exitDistance, entryDistance},
        verts = dnfExpressionVertices[expr];
        pathVerts = Lookup[pathInfo, "PathVertices", {}];
        distances = Lookup[pathInfo, "Distances", <||>];
        exitDistances = Lookup[pathInfo, "ExitDistances", <||>];
        entryDistances = Lookup[pathInfo, "EntryDistances", <||>];
        onPath = If[Intersection[verts, pathVerts] === {}, 1, 0];
        minDistance = If[verts === {}, Infinity, Min[Lookup[distances, #, Infinity] & /@ verts]];
        exitDistance = If[verts === {}, Infinity, Min[Lookup[exitDistances, #, Infinity] & /@ verts]];
        entryDistance = If[verts === {}, Infinity, Min[Lookup[entryDistances, #, Infinity] & /@ verts]];
        {onPath, minDistance, exitDistance, entryDistance, LeafCount[Unevaluated[expr]], ToString[Unevaluated[expr], InputForm]}
    ];

dnfOrderingSummary[conjuncts_List] :=
    Module[{ors, vars},
        ors = Select[conjuncts, !FreeQ[#, Or] &];
        vars = dnfExpressionTrackedVariables /@ ors;
        <|
            "ConjunctCount" -> Length[conjuncts],
            "EqualityCount" -> Count[conjuncts, _Equal],
            "NonOrCount" -> Count[conjuncts, x_ /; FreeQ[x, Or] && !MatchQ[x, _Equal]],
            "OrCount" -> Length[ors],
            "OrLeafCounts" -> (LeafCount /@ ors),
            "OrTrackedVariableKinds" -> Counts[Flatten[dnfTrackedVariableKind /@ Flatten[vars]]]
        |>
    ];

dnfOrderConjuncts[expr_, "original", ___] := dnfTopLevelConjuncts[expr];

dnfOrderConjuncts[expr_, order_String, sys_: Missing["NoSystem"]] :=
    Module[{conjuncts, indexed, pathInfo},
        conjuncts = dnfTopLevelConjuncts[expr];
        indexed = MapIndexed[{#1, First[#2]} &, conjuncts];
        pathInfo = If[StringStartsQ[order, "path"], dnfPathInfo[sys], <||>];
        First /@ SortBy[indexed,
            Switch[order,
                "nonor-stable",
                    {constraintPriority[#[[1]]], #[[2]]} &,
                "or-small",
                    {If[FreeQ[#[[1]], Or], 0, 1], If[FreeQ[#[[1]], Or], 0, LeafCount[#[[1]]]], #[[2]]} &,
                "or-large",
                    {If[FreeQ[#[[1]], Or], 0, 1], If[FreeQ[#[[1]], Or], 0, -LeafCount[#[[1]]]], #[[2]]} &,
                "j-before-u",
                    {If[FreeQ[#[[1]], Or], 0, 1],
                     If[!FreeQ[#[[1]], j[__]], 0, If[!FreeQ[#[[1]], u[__]], 1, 2]], #[[2]]} &,
                "u-before-j",
                    {If[FreeQ[#[[1]], Or], 0, 1],
                     If[!FreeQ[#[[1]], u[__]], 0, If[!FreeQ[#[[1]], j[__]], 1, 2]], #[[2]]} &,
                "transition-before-edge",
                    {If[FreeQ[#[[1]], Or], 0, 1],
                     If[!FreeQ[#[[1]], j[_, _, _]], 0, If[!FreeQ[#[[1]], j[_, _]], 1, 2]], #[[2]]} &,
                "edge-before-transition",
                    {If[FreeQ[#[[1]], Or], 0, 1],
                     If[!FreeQ[#[[1]], j[_, _]] && FreeQ[#[[1]], j[_, _, _]], 0, If[!FreeQ[#[[1]], j[_, _, _]], 1, 2]], #[[2]]} &,
                "few-vars",
                    {If[FreeQ[#[[1]], Or], 0, 1], Length[dnfExpressionTrackedVariables[#[[1]]]], #[[2]]} &,
                "path-aware",
                    {If[FreeQ[#[[1]], Or], 0, 1], Sequence @@ dnfPathRank[#[[1]], pathInfo], #[[2]]} &,
                _,
                    {#[[2]]} &
            ]
        ]
    ];

dnfOrderExpression[expr_, order_String, sys_: Missing["NoSystem"]] :=
    With[{conjuncts = dnfOrderConjuncts[expr, order, sys]},
        Which[
            conjuncts === {}, True,
            Length[conjuncts] === 1, First[conjuncts],
            True, And @@ conjuncts
        ]
    ];

branchKey[branch_Association] :=
    ToString[
        InputForm[
            <|
                "Rules" -> SortBy[Lookup[branch, "Rules", {}], ToString[First[#], InputForm] &],
                "Residuals" -> SortBy[Lookup[branch, "Residuals", {}], ToString[#, InputForm] &]
            |>
        ]
    ];

deduplicateBranches[branches_List] := DeleteDuplicatesBy[branches, branchKey];

trackedFreeQ[expr_] := FreeQ[expr, j[__] | u[__]];

cheapConstraintSimplify[expr_] :=
    Module[{evaluated = expr},
        Which[
            evaluated === True || evaluated === False,
                evaluated,
            trackedFreeQ[evaluated],
                Quiet @ Check[Simplify[evaluated], evaluated],
            MatchQ[evaluated, _Equal] && TrueQ[Subtract @@ (List @@ evaluated) === 0],
                True,
            True,
                evaluated
        ]
    ];

linearRuleFromEquality[eq_Equal, vars_List] :=
    Module[{expr, presentVars, varAlt, coeff, rest, rule = None},
        expr = Expand[Subtract @@ (List @@ eq)];
        presentVars = Select[vars, !FreeQ[expr, #] &];
        If[presentVars === {}, Return[None, Module]];
        varAlt = Alternatives @@ vars;
        Do[
            coeff = Simplify[Coefficient[expr, v]];
            rest  = Simplify[expr - coeff v];
            If[coeff =!= 0 &&
               FreeQ[coeff, varAlt] &&
               FreeQ[rest, v] &&
               Simplify[expr - coeff v - rest] === 0,
                rule = v -> Simplify[-rest/coeff];
                Break[]
            ],
            {v, presentVars}
        ];
        rule
    ];

simplifyResiduals[residuals_List, rules_List] :=
    DeleteCases[cheapConstraintSimplify /@ (residuals /. rules), True];

addBranchConstraint[branch_Association, elem_, vars_List] :=
    Module[{rules, residuals, simplified, rule, newRules, newResiduals},
        rules = Lookup[branch, "Rules", {}];
        residuals = Lookup[branch, "Residuals", {}];
        simplified = cheapConstraintSimplify[elem /. rules];
        Which[
            simplified === True,
                {branch},
            simplified === False,
                {},
            Head[simplified] === And,
                Fold[
                    Function[{branches, conjunct},
                        Join @@ (addBranchConstraint[#, conjunct, vars] & /@ branches)
                    ],
                    {branch},
                    orderedConjuncts[simplified]
                ],
            Head[simplified] === Or,
                Join @@ (addBranchConstraint[branch, #, vars] & /@ orderedAlternatives[simplified]),
            Head[simplified] === Equal,
                rule = linearRuleFromEquality[simplified, vars];
                If[MatchQ[rule, _Rule],
                    newRules = normalizeRules[mergeRules[rules, {rule}]];
                    newResiduals = simplifyResiduals[residuals, newRules];
                    If[MemberQ[newResiduals, False],
                        {},
                        {<|"Rules" -> newRules, "Residuals" -> newResiduals|>}
                    ],
                    newResiduals = simplifyResiduals[Append[residuals, simplified], rules];
                    If[MemberQ[newResiduals, False],
                        {},
                        {<|"Rules" -> rules, "Residuals" -> newResiduals|>}
                    ]
                ],
            True,
                newResiduals = simplifyResiduals[Append[residuals, simplified], rules];
                If[MemberQ[newResiduals, False],
                    {},
                    {<|"Rules" -> rules, "Residuals" -> newResiduals|>}
                ]
        ]
    ];

branchStateReduce[constraints_, allVars_List] :=
    Fold[
        Function[{branches, elem},
            deduplicateBranches[Join @@ (addBranchConstraint[#, elem, allVars] & /@ branches)]
        ],
        {<|"Rules" -> {}, "Residuals" -> {}|>},
        orderedConjuncts[constraints]
    ];

branchStateReduceFromBranches[branches_List, constraints_, allVars_List] :=
    Fold[
        Function[{stateBranches, elem},
            deduplicateBranches[Join @@ (addBranchConstraint[#, elem, allVars] & /@ stateBranches)]
        ],
        branches,
        orderedConjuncts[constraints]
    ];

groundRuleQ[rule_Rule, allVars_List] :=
    allVars === {} || FreeQ[Last[rule], Alternatives @@ allVars];

splitGroundRules[rules_List, allVars_List] :=
    Module[{normalized, ground, symbolic},
        normalized = normalizeRules[rules];
        ground = Select[normalized, groundRuleQ[#, allVars] &];
        symbolic = Complement[normalized, ground];
        {ground, Equal @@@ symbolic}
    ];

branchStateFinalizeResult[branches_List, allVars_List] :=
    Module[{ruleSets, common, commonGround, commonResiduals, branchExprs,
            branchRules, branchGround, branchSymbolic, residuals, residualExpr},
        If[branches === {},
            Return[<|"Rules" -> {}, "Equations" -> False|>, Module]
        ];
        ruleSets = Lookup[#, "Rules", {}] & /@ branches;
        common = If[ruleSets === {}, {}, Fold[Intersection, First[ruleSets], Rest[ruleSets]]];
        {commonGround, commonResiduals} = splitGroundRules[common, allVars];
        branchExprs = Simplify /@ Map[
            Function[branch,
                branchRules = Complement[Lookup[branch, "Rules", {}], common];
                {branchGround, branchSymbolic} = splitGroundRules[branchRules, allVars];
                residuals = Lookup[branch, "Residuals", {}];
                And[
                    And @@ (Equal @@@ branchGround),
                    And @@ branchSymbolic,
                    And @@ residuals
                ] /. commonGround
            ],
            branches
        ];
        branchExprs = DeleteDuplicates @ DeleteCases[DeleteCases[branchExprs, False], True];
        residualExpr = Simplify @ And[
            And @@ (commonResiduals /. commonGround),
            If[branchExprs === {}, True, If[Length[branchExprs] === 1, First[branchExprs], Or @@ branchExprs]]
        ];
        If[residualExpr === True,
            commonGround,
            <|"Rules" -> commonGround, "Equations" -> residualExpr|>
        ]
    ];

branchStateReduceResult[constraints_, allVars_List] :=
    branchStateFinalizeResult[branchStateReduce[constraints, allVars], allVars];

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

ClearAttributes[dnfTraceEvent, HoldFirst];
ClearAttributes[dnfCounterIncrement, HoldFirst];
ClearAttributes[dnfReduceInstrumented, HoldAll];
ClearAttributes[dnfReduceInstrumentedAnd, HoldAll];

dnfTraceEvent[_, event_Association] :=
    Module[{trace, limit},
        $dnfCurrentMonitor["CurrentPhase"] = Lookup[event, "Phase", Lookup[$dnfCurrentMonitor, "CurrentPhase", "DNF"]];
        If[KeyExistsQ[event, "Expression"],
            $dnfCurrentMonitor["LastConjunct"] = Lookup[event, "Expression"];
            $dnfCurrentMonitor["LastConjunctLeafCount"] = LeafCount[Lookup[event, "Expression"]];
            $dnfCurrentMonitor["LastConjunctHead"] = ToString[Head[Lookup[event, "Expression"]]]
        ];
        $dnfCurrentMonitor["TimeoutLocation"] = <|
            "Phase" -> $dnfCurrentMonitor["CurrentPhase"],
            "LastConjunct" -> Lookup[$dnfCurrentMonitor, "LastConjunct", None],
            "LastConjunctHead" -> Lookup[$dnfCurrentMonitor, "LastConjunctHead", None],
            "LastConjunctLeafCount" -> Lookup[$dnfCurrentMonitor, "LastConjunctLeafCount", Missing["Unavailable"]]
        |>;
        limit = Lookup[$dnfCurrentMonitor, "TraceLength", 20];
        trace = Append[Lookup[$dnfCurrentMonitor, "LastEvents", {}], KeyDrop[event, {"Expression"}]];
        $dnfCurrentMonitor["LastEvents"] = If[Length[trace] > limit, Take[trace, -limit], trace];
    ];

dnfCounterIncrement[_, key_String, n_: 1] :=
    ($dnfCurrentMonitor[key] = Lookup[$dnfCurrentMonitor, key, 0] + n);

dnfReduceInstrumented[_, False, monitor_, _String, _] :=
    (dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "FalseSystem"|>]; False);

dnfReduceInstrumented[False, _, monitor_, _String, _] :=
    (dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "FalsePrefix"|>]; False);

dnfReduceInstrumented[xp_, True, monitor_, _String, _] :=
    (dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "TrueSystem"|>]; xp);

dnfReduceInstrumented[xp_, eq_Equal, monitor_, order_String, sys_] :=
    dnfReduceInstrumented[xp, True, eq, monitor, order, sys, 1];

dnfReduceInstrumented[xp_, sys_And, monitor_, order_String, sysObj_] :=
    With[{conjunctsValue = dnfOrderConjuncts[sys, order, sysObj]},
        dnfReduceInstrumentedAnd[xp, conjunctsValue, monitor, order, sysObj, 1]
    ];

dnfReduceInstrumented[xp_, sys_Or, monitor_, order_String, sysObj_] :=
    Module[{results = {}, r, alternatives = List @@ sys},
        dnfCounterIncrement[monitor, "OrsProcessed"];
        dnfCounterIncrement[monitor, "BranchesStarted", Length[alternatives]];
        dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "TopLevelOr", "BranchCount" -> Length[alternatives], "Expression" -> sys|>];
        Scan[
            Function[branch,
            With[{branchValue = branch},
                r = dnfReduceInstrumented[xp, branchValue, monitor, order, sysObj]
            ];
            If[r =!= False,
                AppendTo[results, r]; dnfCounterIncrement[monitor, "BranchesKept"],
                dnfCounterIncrement[monitor, "BranchesDropped"]
            ]
            ],
            alternatives
        ];
        Which[
            results === {}, False,
            Length[results] === 1, First[results],
            True, Or @@ results
        ]
    ];

dnfReduceInstrumented[xp_, leq_, monitor_, _String, _] :=
    (dnfCounterIncrement[monitor, "ConjunctsProcessed"];
     dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "AppendConstraint", "Expression" -> leq|>];
     xp && leq);

dnfReduceInstrumentedAnd[xp_, conjuncts_List, monitor_, order_String, sysObj_, depth_Integer] :=
    Module[{xpAcc = xp, conjunctList = conjuncts, i = 1, elem, sol, fsol, newxp, rest, newRest},
        $dnfCurrentMonitor["MaxRecursionDepth"] = Max[Lookup[$dnfCurrentMonitor, "MaxRecursionDepth", 0], depth];
        While[i <= Length[conjunctList],
            If[xpAcc === False, Return[False, Module]];
            elem = conjunctList[[i]];
            dnfCounterIncrement[monitor, "ConjunctsProcessed"];
            dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "ProcessConjunct", "Index" -> i, "Depth" -> depth, "Expression" -> elem|>];
            Which[
                elem === True,
                    i++,
                elem === False,
                    dnfCounterIncrement[monitor, "BranchesDropped"];
                    Return[False, Module],
                Head[elem] === Or,
                    dnfCounterIncrement[monitor, "OrsProcessed"];
                    dnfCounterIncrement[monitor, "BranchesStarted", Length[List @@ elem]];
                    rest = If[i >= Length[conjunctList], True, And @@ conjunctList[[i + 1 ;;]]];
                    Return[
                        Module[{results = {}, r},
                            Scan[
                                Function[branch,
                                dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "EnterBranch", "Depth" -> depth + 1, "Expression" -> branch|>];
                                With[{prefixValue = xpAcc, restValue = rest, branchValue = branch, depthValue = depth + 1},
                                    r = dnfReduceInstrumented[prefixValue, restValue, branchValue, monitor, order, sysObj, depthValue]
                                ];
                                If[r =!= False,
                                    AppendTo[results, r]; dnfCounterIncrement[monitor, "BranchesKept"],
                                    dnfCounterIncrement[monitor, "BranchesDropped"]
                                ]
                                ],
                                List @@ elem
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
                    dnfCounterIncrement[monitor, "EqualityAttempts"];
                    sol  = Quiet[Solve[elem], Solve::svars];
                    fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
                    If[fsol =!= {},
                        dnfCounterIncrement[monitor, "EqualitySolves"];
                        dnfCounterIncrement[monitor, "Substitutions", Length[fsol]]
                    ];
                    newxp = substituteSolution[xpAcc, fsol];
                    If[newxp === False, Return[False, Module]];
                    xpAcc = newxp && elem;
                    If[i < Length[conjunctList],
                        newRest = substituteSolution[And @@ conjunctList[[i + 1 ;;]], fsol];
                        conjunctList = Join[conjunctList[[;; i]], dnfOrderConjuncts[newRest, order, sysObj]]
                    ];
                    i++,
                True,
                    xpAcc = xpAcc && elem;
                    i++
            ]
        ];
        xpAcc
    ];

dnfReduceInstrumented[xp_, rst_, fst_Equal, monitor_, order_String, sysObj_, depth_Integer] :=
    Module[{sol, fsol, newxp, newrst},
        $dnfCurrentMonitor["MaxRecursionDepth"] = Max[Lookup[$dnfCurrentMonitor, "MaxRecursionDepth", 0], depth];
        dnfCounterIncrement[monitor, "ConjunctsProcessed"];
        dnfCounterIncrement[monitor, "EqualityAttempts"];
        dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "ProcessBranchEquality", "Depth" -> depth, "Expression" -> fst|>];
        sol  = Quiet[Solve[fst], Solve::svars];
        fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
        If[fsol =!= {},
            dnfCounterIncrement[monitor, "EqualitySolves"];
            dnfCounterIncrement[monitor, "Substitutions", Length[fsol]]
        ];
        newxp = substituteSolution[xp, fsol];
        If[newxp === False,
            False,
            newrst = substituteSolution[rst, fsol];
            With[{prefixValue = newxp && fst, restValue = dnfOrderExpression[newrst, order, sysObj]},
                dnfReduceInstrumented[prefixValue, restValue, monitor, order, sysObj]
            ]
        ]
    ];

dnfReduceInstrumented[xp_, rst_, fst_, monitor_, order_String, sysObj_, depth_Integer] :=
    ($dnfCurrentMonitor["MaxRecursionDepth"] = Max[Lookup[$dnfCurrentMonitor, "MaxRecursionDepth", 0], depth];
     dnfCounterIncrement[monitor, "ConjunctsProcessed"];
     dnfTraceEvent[monitor, <|"Phase" -> "DNF", "Action" -> "ProcessBranchConstraint", "Depth" -> depth, "Expression" -> fst|>];
     With[{prefixValue = xp && fst, restValue = dnfOrderExpression[rst, order, sysObj]},
         dnfReduceInstrumented[prefixValue, restValue, monitor, order, sysObj]
     ]);

Options[dnfReduceDiagnosticReport] = {
    "Order" -> "original",
    "Timeout" -> Infinity,
    "TraceLength" -> 20
};

dnfReduceDiagnosticReport[sys_?mfgSystemQ, OptionsPattern[]] :=
    Module[{order, timeout, traceLength, tInputs, inputs, constraints, allVars,
            rulesAcc, conjuncts, orderedConstraints, orderingSummary, monitor,
            tDNF, reduced, status, parsed, result, summary},
        order = ToLowerCase[ToString[OptionValue["Order"]]];
        timeout = OptionValue["Timeout"];
        traceLength = OptionValue["TraceLength"];
        If[!criticalCongestionSystemQ[sys],
            Return[<|"Status" -> "Failure", "Reason" -> "NonCriticalCongestionSystem"|>, Module]
        ];

        {tInputs, inputs} = AbsoluteTiming[buildSolverInputs[sys]];
        {constraints, allVars, rulesAcc} = inputs;
        conjuncts = dnfTopLevelConjuncts[constraints];
        orderedConstraints = dnfOrderExpression[constraints, order, sys];
        orderingSummary = Association[
            dnfOrderingSummary[conjuncts],
            "Order" -> order,
            "OrderedSameMultisetQ" -> (Sort[ToString[#, InputForm] & /@ conjuncts] ===
                Sort[ToString[#, InputForm] & /@ dnfTopLevelConjuncts[orderedConstraints]])
        ];
        monitor = <|
            "TraceLength" -> traceLength,
            "LastEvents" -> {},
            "CurrentPhase" -> "DNF",
            "ConjunctsProcessed" -> 0,
            "OrsProcessed" -> 0,
            "BranchesStarted" -> 0,
            "BranchesKept" -> 0,
            "BranchesDropped" -> 0,
            "EqualityAttempts" -> 0,
            "EqualitySolves" -> 0,
            "Substitutions" -> 0,
            "MaxRecursionDepth" -> 0,
            "ExpressionLeafCount" -> LeafCount[constraints],
            "OrderedExpressionLeafCount" -> LeafCount[orderedConstraints]
        |>;
        $dnfCurrentMonitor = monitor;

        {tDNF, reduced} = AbsoluteTiming[
            TimeConstrained[
                With[{orderedConstraintsValue = orderedConstraints, orderValue = order, sysValue = sys},
                    Quiet[
                        dnfReduceInstrumented[True, orderedConstraintsValue, monitor, orderValue, sysValue],
                        Set::write
                    ]
                ],
                timeout,
                $TimedOut
            ]
        ];
        monitor = $dnfCurrentMonitor;
        status = If[reduced === $TimedOut, "Timeout", "OK"];
        parsed = If[status === "OK",
            $dnfCurrentMonitor["CurrentPhase"] = "Parse";
            parseDNFReduceResult[reduced, allVars],
            Missing["TimedOut"]
        ];
        result = If[status === "OK", attachAccumulatedRules[parsed, rulesAcc], $TimedOut];
        summary = <|
            "SolverInputSeconds" -> tInputs,
            "DNFSeconds" -> tDNF,
            "PreRuleCount" -> Length[rulesAcc],
            "VariableCount" -> Length[allVars],
            "PreprocessingConjunctCount" -> Length[conjuncts],
            "ConjunctsProcessed" -> Lookup[monitor, "ConjunctsProcessed", 0],
            "OrsProcessed" -> Lookup[monitor, "OrsProcessed", 0],
            "BranchesStarted" -> Lookup[monitor, "BranchesStarted", 0],
            "BranchesKept" -> Lookup[monitor, "BranchesKept", 0],
            "BranchesDropped" -> Lookup[monitor, "BranchesDropped", 0],
            "EqualityAttempts" -> Lookup[monitor, "EqualityAttempts", 0],
            "EqualitySolves" -> Lookup[monitor, "EqualitySolves", 0],
            "Substitutions" -> Lookup[monitor, "Substitutions", 0],
            "MaxRecursionDepth" -> Lookup[monitor, "MaxRecursionDepth", 0],
            "ExpressionLeafCount" -> Lookup[monitor, "ExpressionLeafCount", Missing["Unavailable"]],
            "OrderedExpressionLeafCount" -> Lookup[monitor, "OrderedExpressionLeafCount", Missing["Unavailable"]],
            "ResultLeafCount" -> If[status === "OK", LeafCount[reduced], Missing["TimedOut"]],
            "ResultKind" -> If[status === "OK", solutionResultKind[result], "TIMEOUT"]
        |>;
        <|
            "Status" -> status,
            "Order" -> order,
            "Result" -> result,
            "ReducedExpression" -> If[status === "OK", reduced, Missing["TimedOut"]],
            "Summary" -> summary,
            "OrderingSummary" -> orderingSummary,
            "LastConjunct" -> Lookup[monitor, "LastConjunct", Missing["None"]],
            "LastEvents" -> Lookup[monitor, "LastEvents", {}],
            "TimeoutLocation" -> Lookup[monitor, "TimeoutLocation", <||>]
        |>
    ];

(* reduceSystem only handles the critical-congestion structural system. The
   system object carries the scenario Hamiltonian so this check can reject any
   non-1 default Alpha or per-edge EdgeAlpha before Reduce sees nonlinear
   placeholder equations. *)

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

harvestDNFBranch[False, _] := False;

harvestDNFBranch[branch_, allVars_List] :=
    Module[{conjuncts, eqs, sol, rules, reducedConjuncts, residual,
            residualVars},
        conjuncts = flattenConjuncts[branch];
        eqs = Select[conjuncts, MatchQ[#, _Equal] &];
        sol = If[eqs === {} || allVars === {},
            {{}},
            Quiet[Solve[eqs, allVars], Solve::svars]
        ];
        If[!ListQ[sol] || sol === {}, Return[False, Module]];
        rules = Cases[First[sol],
            r : Rule[v_, rhs_] /;
                MemberQ[allVars, v] &&
                FreeQ[rhs, Alternatives @@ allVars] :> r
        ];
        reducedConjuncts = Simplify /@ (conjuncts /. rules);
        reducedConjuncts = DeleteCases[reducedConjuncts, True];
        residual = If[reducedConjuncts === {}, True, And @@ reducedConjuncts];
        residual = Simplify[residual];
        If[residual === False, Return[False, Module]];
        residualVars = Select[
            Variables[residual /. {Equal -> List, Or -> List, And -> List}],
            MemberQ[allVars, #] &
        ];
        If[residualVars === {},
            If[TrueQ[Simplify[residual]],
                <|"Rules" -> rules, "Equations" -> True|>,
                False
            ],
            <|"Rules" -> rules, "Equations" -> residual|>
        ]
    ];

parseDNFReduceResult[reduced_, allVars_List] :=
    Module[{branches, harvested, ruleSets, common, branchExprs},
        If[reduced === False, Return[<|"Rules" -> {}, "Equations" -> False|>, Module]];
        branches = If[Head[reduced] === Or, List @@ reduced, {reduced}];
        harvested = DeleteCases[harvestDNFBranch[#, allVars] & /@ branches, False];
        If[harvested === {},
            Return[<|"Rules" -> {}, "Equations" -> False|>, Module]
        ];
        If[Length[harvested] === 1,
            With[{rules = Lookup[First[harvested], "Rules", {}],
                  residual = Lookup[First[harvested], "Equations", True]},
                Return[
                    If[residual === True,
                        rules,
                        <|"Rules" -> rules, "Equations" -> residual|>
                    ],
                    Module
                ]
            ]
        ];
        ruleSets = Lookup[#, "Rules", {}] & /@ harvested;
        common = If[ruleSets === {}, {}, Fold[Intersection, First[ruleSets], Rest[ruleSets]]];
        branchExprs = Simplify /@ Map[
            With[{branchRules = Complement[Lookup[#, "Rules", {}], common],
                  residual = Lookup[#, "Equations", True]},
                And[
                    And @@ (Equal @@@ branchRules),
                    residual
                ] /. common
            ] &,
            harvested
        ];
        branchExprs = DeleteDuplicates @ DeleteCases[DeleteCases[branchExprs, False], True];
        If[branchExprs === {},
            common,
            <|"Rules" -> common,
              "Equations" -> Simplify @ If[Length[branchExprs] === 1, First[branchExprs], Or @@ branchExprs]|>
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

        kind = solutionResultKind[sol];

        If[!MemberQ[{"Rules", "Branched", "Parametric", "Residual", "NoSolution"}, kind],
            Return[If[returnReportQ,
                <|"Valid" -> False, "Kind" -> kind,
                  "Reason" -> "UnrecognizedSolutionFormat", "BlockChecks" -> <||>|>,
                False], Module]
        ];

        rules        = If[ListQ[sol], sol, Lookup[sol, "Rules", {}]];
        ruleExitVals = Normal @ systemData[sys, "RuleExitValues"];

        (* For residual-bearing results: fail fast if residual equations are inconsistent. *)
        If[AssociationQ[sol],
            With[{eqs = Lookup[sol, "Equations", True]},
                If[eqs === False || TrueQ[Quiet @ Check[Simplify[eqs] === False, False]],
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

computeKirchhoffResidual[sys_?mfgSystemQ, sol_] :=
    Module[{rules, balanceSplitting, balanceGathering, diffs, numDiffs},
        rules = If[ListQ[sol], sol, Lookup[sol, "Rules", {}]];
        If[!ListQ[rules], Return[Missing["NotAvailable"], Module]];

        balanceSplitting = systemData[sys, "BalanceSplittingFlows"];
        balanceGathering = systemData[sys, "BalanceGatheringFlows"];

        If[MissingQ[balanceSplitting], balanceSplitting = {}];
        If[MissingQ[balanceGathering], balanceGathering = {}];

        diffs = Join[balanceSplitting, balanceGathering] /. rules;
        numDiffs = Quiet @ Check[N[diffs], $Failed];

        If[numDiffs === $Failed || !VectorQ[numDiffs, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];

        If[Length[numDiffs] === 0, 0., Max[Abs[numDiffs]]]
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
    With[{s = Quiet @ Check[Simplify[lhs == rhs], $Failed]},
        Which[
            s === True,  True,
            s === False, False,
            True,        With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
                             Which[
                                 d === $Failed, Indeterminate,
                                 NumericQ[d],   Abs[d] <= tol,
                                 True,          Indeterminate
                             ]]
        ]
    ];

checkBlock[lhs_ >= rhs_, tol_] :=
    With[{s = Quiet @ Check[Simplify[lhs >= rhs], $Failed]},
        Which[
            s === True,  True,
            s === False, False,
            True,        With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
                             Which[d === $Failed, Indeterminate, NumericQ[d], d >= -tol, True, Indeterminate]]
        ]
    ];

checkBlock[lhs_ <= rhs_, tol_] :=
    With[{s = Quiet @ Check[Simplify[lhs <= rhs], $Failed]},
        Which[
            s === True,  True,
            s === False, False,
            True,        With[{d = Quiet @ Check[N[rhs - lhs], $Failed]},
                             Which[d === $Failed, Indeterminate, NumericQ[d], d >= -tol, True, Indeterminate]]
        ]
    ];

checkBlock[lhs_ > rhs_, tol_] :=
    With[{s = Quiet @ Check[Simplify[lhs > rhs], $Failed]},
        Which[
            s === True,  True,
            s === False, False,
            True,        With[{d = Quiet @ Check[N[lhs - rhs], $Failed]},
                             Which[d === $Failed, Indeterminate, NumericQ[d], d > tol, True, Indeterminate]]
        ]
    ];

checkBlock[lhs_ < rhs_, tol_] :=
    With[{s = Quiet @ Check[Simplify[lhs < rhs], $Failed]},
        Which[
            s === True,  True,
            s === False, False,
            True,        With[{d = Quiet @ Check[N[rhs - lhs], $Failed]},
                             Which[d === $Failed, Indeterminate, NumericQ[d], d > tol, True, Indeterminate]]
        ]
    ];

checkBlock[Inequality[a_, Less, b_, Less, c_], tol_] :=
    With[{s = Quiet @ Check[Simplify[a < b < c], $Failed]},
        Which[
            s === True,  True,
            s === False, False,
            True,        With[{d1 = Quiet @ Check[N[b - a], $Failed], d2 = Quiet @ Check[N[c - b], $Failed]},
                             Which[
                                 d1 === $Failed || d2 === $Failed, Indeterminate,
                                 NumericQ[d1] && NumericQ[d2],     d1 > tol && d2 > tol,
                                 True,                             Indeterminate
                             ]]
        ]
    ];

checkBlock[expr_, tol_] :=
    With[{s = Quiet @ Check[Simplify[expr], $Failed]},
        Which[s === True, True, s === False, False, True, Indeterminate]
    ];

End[];

EndPackage[];
