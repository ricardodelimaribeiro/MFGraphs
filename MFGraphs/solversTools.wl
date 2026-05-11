(* Wolfram Language package *)
(*
   solversTools.wl: symbolic solvers for mfgSystem objects.

   Three solvers share a common preprocessing pipeline (buildSolverInputs):
   collect structural equations, apply exit-value rules, eliminate equalities
   via accumulateEqualityRules. They differ only in the final step:
     reduceSystem         -- Wolfram Reduce over Reals (CAD/real
                             quantifier-elimination backend for polynomial
                             systems)
     dnfReduceSystem      -- package recursive DNF reducer:
                             equality substitution + disjunction expansion
     optimizedDNFReduceSystem
                          -- package branch-state exact DNF reduction,
                             falling back to dnfReduce
     activeSetReduceSystem
                          -- package active-set/complementarity branch
                             enumeration, falling back to dnfReduce
     booleanReduceSystem  -- Wolfram BooleanConvert to DNF, then Reduce per
                             disjunct
     booleanMinimizeSystem
                          -- Wolfram BooleanMinimize to minimal DNF, then
                             Reduce per disjunct
     booleanMinimizeReduceSystem
                          -- package arm-prune + component decomposition,
                             BooleanMinimize per component, Reduce per
                             disjunct
     findInstanceSystem   -- Wolfram FindInstance over Reals
*)

BeginPackage["solversTools`", {"primitives`", "utilities`", "systemTools`"}];

reduceSystem::usage =
"reduceSystem[sys] reduces the structural equations, flow balance, \
non-negativity constraints, and complementarity conditions of the \
mfgSystem sys using Wolfram Reduce over the Reals; for polynomial real \
systems this is the CAD/real quantifier-elimination backend. Includes AltOptCond \
(switching-cost optimality complementarity) and \
IneqSwitchingByVertex (switching-cost optimality inequalities). \
Fails for non-critical congestion systems where Alpha != 1 on any edge. \
Returns a list of rules when the system is fully determined, or \
<|\"Rules\" -> rules, \"Residual\" -> residual|> when underdetermined.";

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
one conjunct elem from sys. Implemented with direct recursion: the Or case \
spawns one recursive call per branch and the Equal case recurses on the \
substituted remainder. Deeply nested Or-chains may approach $RecursionLimit; \
see dnfReduceProcedural for a stack-based iterative alternative.";

dnfReduceSystem::usage =
"dnfReduceSystem[sys] solves the mfgSystem sys using linear preprocessing \
followed by dnfReduce instead of Reduce. Handles cases where Reduce times \
out by using equality-substitution and disjunction-distribution. Returns a \
list of rules when fully determined, or \
<|\"Rules\" -> rules, \"Residual\" -> residual|> when underdetermined. \
Fails for non-critical congestion systems where Alpha != 1 on any edge.";

optimizedDNFReduceSystem::usage =
"optimizedDNFReduceSystem[sys] is an opt-in exact solver for critical \
congestion systems. It follows the same preprocessing and output contract as \
dnfReduceSystem. It carries small DNF branch families directly as rules plus \
residual constraints, and falls back to the proven exact DNF reducer for \
larger residual variable sets. \
Returns a list of rules when fully determined, or \
<|\"Rules\" -> rules, \"Residual\" -> residual|> when underdetermined. \
Fails for non-critical congestion systems where Alpha != 1 on any edge.";

activeSetReduceSystem::usage =
"activeSetReduceSystem[sys] is an opt-in exact active-set solver for the \
critical-congestion linear complementarity structure. It enumerates small \
complementarity alternatives incrementally with exact linear substitution and \
falls back to the proven exact DNF reducer for larger residual variable sets. \
Returns the same rule/residual shape as dnfReduceSystem. Fails for \
non-critical congestion systems where Alpha != 1 on any edge.";

booleanReduceSystem::usage =
"booleanReduceSystem[sys] solves the mfgSystem sys by converting the \
preprocessed constraint system to DNF via BooleanConvert, then calling \
Reduce independently on each disjunct. This is DNF conversion followed by \
real quantifier elimination / CAD-style solving per pure conjunction. Each \
disjunct has no Or, so Reduce avoids case-splitting. Non-False results are collected; \
if the system has a unique equilibrium all non-False results are equivalent. \
Returns a list of rules when fully determined, or \
<|\"Rules\" -> rules, \"Residual\" -> residual|> when underdetermined. \
Fails for non-critical congestion systems where Alpha != 1 on any edge. \
Options: \"DisjunctTimeout\" (default 30s per Reduce call), \
\"ReturnAll\" (default False; True returns all non-False parsed results).";

booleanReduceSystem::multisol =
"booleanReduceSystem found `1` non-False disjuncts with differing rules. \
Returning first; use \"ReturnAll\" -> True to inspect all results.";

booleanMinimizeSystem::usage =
"booleanMinimizeSystem[sys] is a head-to-head variant of booleanReduceSystem \
that calls BooleanMinimize[constraints, \"DNF\"] in place of \
BooleanConvert[constraints, \"DNF\"]. This is exact minimal-DNF Boolean \
minimization, analogous to classical two-level SOP minimization such as \
Quine-McCluskey/Petrick-style methods, followed by Reduce per disjunct. \
Wolfram does not document BooleanMinimize as a specific QMC/Petrick \
implementation. Same preprocessing and return shape as booleanReduceSystem. \
Use to compare the two Boolean-stage operations on identical input. Fails for \
non-critical congestion systems where Alpha != 1 on any edge. Options: \
\"DisjunctTimeout\" (default 30s per Reduce call), \
\"ReturnAll\" (default False).";

booleanMinimizeSystem::multisol =
"booleanMinimizeSystem found `1` non-False disjuncts with differing rules. \
Returning first; use \"ReturnAll\" -> True to inspect all results.";

booleanMinimizeReduceSystem::usage =
"booleanMinimizeReduceSystem[sys] solves the mfgSystem sys by attacking the \
disjunctive structure of the preprocessed constraint system before DNF \
expansion. It (1) prunes individual complementarity arms that are infeasible \
against the linear part via FindInstance; (2) decomposes the surviving \
disjunctive atoms into connected components by shared variables; \
(3) BooleanMinimizes each component to minimal DNF and Reduces per disjunct. \
Returns the same rule/residual shape as booleanReduceSystem. Fails for \
non-critical congestion systems where Alpha != 1 on any edge. \
Options: \"ArmTimeout\" (default 2s per FindInstance arm check), \
\"DisjunctTimeout\" (default 30s per Reduce call), \
\"ReturnAll\" (default False).";

booleanMinimizeReduceSystem::multisol =
"booleanMinimizeReduceSystem found `1` non-False disjuncts in a component with \
differing rules. Returning first; use \"ReturnAll\" -> True to inspect all results.";

findInstanceSystem::usage =
"findInstanceSystem[sys] solves the mfgSystem sys by collecting and \
linearly preprocessing constraints, then calling FindInstance over the \
remaining real variables. This is real satisfiability / instance finding using \
Wolfram's real-system solver backend. Returns one feasible list of rules. If no instance is \
found or the final solve times out, returns \
<|\"Rules\" -> accumulatedRules, \"Residual\" -> False|>. \
Fails for non-critical congestion systems where Alpha != 1 on any edge. \
Options: \"Timeout\" (default Infinity).";

solutionReport::usage =
"solutionReport[sys, sol] returns a read-only diagnostic association for an \
existing solver result. The report includes result kind, validation report, \
Kirchhoff residual, residual branch diagnostics, primary residual variables, \
and transition-flow determinacy diagnostics. It does not rewrite sol.";

solutionBranchCostReport::usage =
"solutionBranchCostReport[sys, sol] ranks top-level residual branches by a \
diagnostic critical-congestion objective. It reports terminal exit cost, \
switching cost, edge-flow quadratic cost, total objective, determined and \
residual variables, and validation status for each branch.";

directCriticalSystem::usage =
"directCriticalSystem[sys] is an explicit opt-in solver for critical \
congestion systems with all numeric zero switching costs and equal numeric \
exit costs. It uses graph distance to exits to forbid flow directions that \
move farther from every exit, then solves the resulting feasibility system. \
Returns the standard raw solver shape or Failure[\"directCriticalSystem\", ...].";

flowFirstCriticalSystem::usage =
"flowFirstCriticalSystem[sys] is an explicit opt-in solver for critical \
congestion systems. It minimizes the critical quadratic flow objective over \
Kirchhoff/nonnegative flow constraints, then recovers value variables from \
the remaining system. Returns the standard raw solver shape or \
Failure[\"flowFirstCriticalSystem\", ...].";

Begin["`Private`"];

(* The following functions are declared and implemented in utilities.wl:
   - mergeRules
   - normalizeRules
*)

(* Variable-count cutoff below which branch-state enumeration beats full DNF
   reduction. Re-tuned 2026-05-10 to 0 (branch-state only when there are zero
   tracked vars to enumerate): once branchStateFinalizeResult was taught to
   promote residual equalities into rules, branch-state stopped winning on any
   non-trivial benchmark case while still costing 25-35% on chain-3v-2exit and
   chain-3-midentry. Kept as a knob rather than deleting the path so a future
   case profile with high branching factor can re-enable it via the sweep
   script. Shared by optimizedDNFReduceSystem and activeSetReduceSystem. *)
$optimizedDNFVarThreshold = 0;

trackedVarQ::usage =
"trackedVarQ[var] returns True for solver variables tracked during Reduce post-processing.";

collectTrackedVars::usage =
"collectTrackedVars[sys, expr] returns the tracked solver variables (j[__] | u[__]) appearing in expr, in u-first construction order. Pulls the pre-grouped Us, Js, Jts lists straight from sys instead of re-deriving via Variables, so the column-pivot priority required by rulesFromEqualities is established once at the system level rather than re-sorted per call.";

topLevelEquations::usage =
"topLevelEquations[constraints] extracts top-level Equal expressions from an And expression or single constraint.";

rulesFromEqualities::usage =
"rulesFromEqualities[equations, vars, existingRules] solves the equalities for new eliminable rules. Assumes linear inputs; see makeSystem contract.";

accumulateEqualityRules::usage =
"accumulateEqualityRules[baseConstraints, vars, initRules] repeatedly extracts substitution rules from equalities in the constraint system. Assumes linear inputs; see makeSystem contract.";

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
exit-value rules, and runs accumulateEqualityRules. Returns \
{constraints, allVars, rulesAcc} ready for the final solve step.";

booleanDnfReducePipeline::usage =
"booleanDnfReducePipeline[dnfFn, multisolEmit][constraints, allVars, timeout, returnAll] \
is the shared body of booleanReduceSystem and booleanMinimizeSystem. dnfFn is the Boolean \
DNF operation (BooleanConvert or BooleanMinimize); multisolEmit is a 1-arg callback that \
emits the solver-specific ::multisol message when distinct disjuncts produce differing \
rule sets.";

checkBlock::usage =
"checkBlock[expr, tol] evaluates a substituted constraint block using tolerance tol for numeric comparisons.";

solutionResultKind::usage =
"solutionResultKind[sol] classifies solver output by its residual content.";

primarySolutionVarQ::usage =
"primarySolutionVarQ[var] returns True for primary solver variables (j[a,b] or u[a,b]).";

transitionFlowVarQ::usage =
"transitionFlowVarQ[var] returns True for transition flow variables (j[a,b,c]).";

residualPrimaryVariables::usage =
"residualPrimaryVariables[residual] returns a list of primary variables (j[a,b], u[a,b]) present in the residual expression.";

residualTransitionFlows::usage =
"residualTransitionFlows[residual] returns a list of transition flow variables (j[a,b,c]) present in the residual expression.";

residualTrackedVariables::usage =
"residualTrackedVariables[residual] returns a list of all tracked variables present in the residual expression.";

orInvolvingPrimaryVariableQ::usage =
"orInvolvingPrimaryVariableQ[residual] returns True if the residual contains an Or expression involving primary variables.";

solutionRules::usage =
"solutionRules[sol] extracts the rule list from a solver result (list or association).";

solutionVariableDiagnostics::usage =
"solutionVariableDiagnostics[sys, sol] reports primary solution classification \
and transition-flow determinacy for a solver result.";

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
    Module[{result},
        resetSolveCacheCounters[];
        result = withCriticalCongestionSolver[sys, "dnfReduceSystem",
            Function[{constraints, allVars},
                parseDNFReduceResult[dnfReduce[True, constraints], allVars]
            ],
            buildSolverInputs,
            attachAccumulatedRules
        ];
        If[TrueQ[$MFGraphsVerbose], printSolveCacheSummary["dnfReduceSystem"]];
        result
    ];

optimizedDNFReduceSystem[sys_?mfgSystemQ] :=
    withCriticalCongestionSolver[sys, "optimizedDNFReduceSystem",
        Function[{constraints, allVars},
            If[Length[allVars] <= $optimizedDNFVarThreshold,
                branchStateReduceResult[constraints, allVars],
                parseDNFReduceResult[dnfReduce[True, constraints], allVars]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

activeSetReduceSystem[sys_?mfgSystemQ] :=
    withCriticalCongestionGuard[sys, "activeSetReduceSystem",
    Module[{initRules, deterministicConstraints, activeConstraints, allVars,
            detReduced, rulesAcc, activeReduced, result, branches},
        initRules = Join[
            Normal @ With[{v = systemData[sys, "RuleBalanceGatheringFlows"]}, If[AssociationQ[v], v, <||>]],
            With[{v = systemData[sys, "ZeroSwitchUEqualities"]}, If[ListQ[v], v, {}]]
        ];
        deterministicConstraints = And[
            And @@ systemData[sys, "EqEntryIn"],
            systemData[sys, "EqBalanceSplittingFlows"],
            systemData[sys, "EqBalanceGatheringFlows"],
            And @@ systemData[sys, "EqGeneral"],
            And @@ systemData[sys, "IneqJs"],
            And @@ systemData[sys, "IneqJts"],
            And @@ systemData[sys, "IneqSwitchingByVertex"],
            And @@ systemData[sys, "IneqExitValues"]
        ];
        activeConstraints = And[
            And @@ systemData[sys, "AltFlows"],
            And @@ systemData[sys, "AltTransitionFlows"],
            And @@ systemData[sys, "AltOptCond"],
            And @@ systemData[sys, "AltExitCond"]
        ];
        allVars = collectTrackedVars[sys, And[deterministicConstraints, activeConstraints] /. initRules];
        {detReduced, rulesAcc} = accumulateEqualityRules[deterministicConstraints, allVars, initRules];
        activeReduced = activeConstraints /. rulesAcc;
        (* Eliminated vars are exactly the LHS of rulesAcc that were not in initRules.
           Set-subtract instead of re-walking the reduced expression tree. *)
        allVars = DeleteCases[allVars, Alternatives @@ (First /@ rulesAcc)];
        result = If[Length[allVars] <= $optimizedDNFVarThreshold,
            branches = branchStateReduce[activeReduced, allVars];
            branches = branchStateReduceFromBranches[branches, detReduced, allVars];
            branchStateFinalizeResult[branches, allVars],
            parseDNFReduceResult[dnfReduce[True, And[detReduced, activeReduced]], allVars]
        ];
        attachAccumulatedRules[result, rulesAcc]
    ]
    ];

(* Shared body for booleanReduceSystem / booleanMinimizeSystem.
   dnfFn is the Boolean DNF operation (BooleanConvert or BooleanMinimize);
   multisolEmit is a 1-arg callback that emits the solver-specific
   ::multisol message when distinct disjuncts produce differing rule sets. *)
booleanDnfReducePipeline[dnfFn_, multisolEmit_][constraints_, allVars_, timeout_, returnAll_] :=
    Module[{dnf, disjuncts, reducedList, nonFalse, parsed},
        dnf         = dnfFn[constraints, "DNF"];
        disjuncts   = If[Head[dnf] === Or, List @@ dnf, {dnf}];
        reducedList = TimeConstrained[Reduce[#, allVars, Reals], timeout, $Aborted] & /@ disjuncts;
        nonFalse    = Select[reducedList, # =!= False && # =!= $Aborted &];
        If[nonFalse === {},
            Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
        ];
        parsed = parseReduceResult[#, allVars] & /@ nonFalse;
        If[Length[parsed] > 1,
            With[{ruleSets = (If[ListQ[#], #, Lookup[#, "Rules", {}]] & /@ parsed)},
                If[!SameQ @@ (Sort /@ ruleSets), multisolEmit[Length[parsed]]]
            ]
        ];
        If[returnAll, parsed, First[parsed]]
    ];

Options[booleanReduceSystem] = {
    "DisjunctTimeout" -> 30,
    "ReturnAll"       -> False
};

booleanReduceSystem[sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    withCriticalCongestionSolver[sys, "booleanReduceSystem",
        Function[{constraints, allVars},
            booleanDnfReducePipeline[
                BooleanConvert,
                Function[{n}, Message[booleanReduceSystem::multisol, n]]
            ][
                constraints, allVars,
                OptionValue[booleanReduceSystem, {opts}, "DisjunctTimeout"],
                TrueQ[OptionValue[booleanReduceSystem, {opts}, "ReturnAll"]]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

(* ---------- booleanMinimizeSystem (plain) ---------- *)

Options[booleanMinimizeSystem] = {
    "DisjunctTimeout" -> 30,
    "ReturnAll"       -> False
};

booleanMinimizeSystem[sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    withCriticalCongestionSolver[sys, "booleanMinimizeSystem",
        Function[{constraints, allVars},
            booleanDnfReducePipeline[
                BooleanMinimize,
                Function[{n}, Message[booleanMinimizeSystem::multisol, n]]
            ][
                constraints, allVars,
                OptionValue[booleanMinimizeSystem, {opts}, "DisjunctTimeout"],
                TrueQ[OptionValue[booleanMinimizeSystem, {opts}, "ReturnAll"]]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

(* ---------- booleanMinimizeReduceSystem (levered) ---------- *)

(* Lever 1: prune complementarity arms that are infeasible against the linear part. *)
pruneDisjunctiveArms[linearAtoms_List, disjAtoms_List, allVars_List, armTimeout_] :=
    Module[{lin = linearAtoms, disj = disjAtoms, changed = True, status = "Unchanged",
            newDisj, atom, arms, kept, fi, infeasibleAll = False},
        While[changed && !infeasibleAll,
            changed = False;
            newDisj = {};
            Do[
                atom = disj[[k]];
                arms = If[Head[atom] === Or, List @@ atom, {atom}];
                kept = Select[arms,
                    Function[arm,
                        fi = TimeConstrained[
                            FindInstance[And @@ Append[lin, arm], allVars, Reals, 1],
                            armTimeout, $Aborted
                        ];
                        Which[
                            fi === $Aborted,             True,    (* keep on timeout *)
                            MatchQ[fi, {{___Rule}, ___}], True,
                            True,                         False
                        ]
                    ]
                ];
                Which[
                    kept === {},
                        infeasibleAll = True;
                        Break[],
                    Length[kept] === 1,
                        lin = Append[lin, First[kept]];
                        changed = True;
                        status = "Reduced",
                    Length[kept] < Length[arms],
                        AppendTo[newDisj, Or @@ kept];
                        changed = True;
                        status = "Reduced",
                    True,
                        AppendTo[newDisj, atom]
                ],
                {k, Length[disj]}
            ];
            If[!infeasibleAll, disj = newDisj]
        ];
        If[infeasibleAll,
            {lin, disj, "Infeasible"},
            {lin, disj, status}
        ]
    ];

(* Lever 2: decompose into connected components of the variable-cooccurrence graph
   built from BOTH linear and disjunctive atoms. Returns one record per component:
   {compDisjAtoms, compLinearAtoms, compVars}. Variables are coupled if they appear
   together in any atom; this preserves linkages through linear equalities/inequalities
   that would otherwise split coupled disjuncts into different components. *)
varsOf[expr_, varSet_] :=
    DeleteDuplicates@Cases[expr, v_ /; KeyExistsQ[varSet, v], {0, Infinity}, Heads -> False];

decomposeByVariableSharing[linearAtoms_List, disjAtoms_List, allVars_List] :=
    Module[{varSet, linVars, disjVars, edges, varComponents, varToComp,
            compIndex, n, results},
        n = Length[allVars];
        If[n === 0, Return[{{disjAtoms, linearAtoms, {}}}]];
        varSet = Association[# -> True & /@ allVars];
        linVars  = varsOf[#, varSet] & /@ linearAtoms;
        disjVars = varsOf[#, varSet] & /@ disjAtoms;
        edges = DeleteDuplicates@Flatten[
            Function[vs,
                If[Length[vs] >= 2,
                    Table[vs[[1]] <-> vs[[k]], {k, 2, Length[vs]}],
                    {}
                ]
            ] /@ Join[linVars, disjVars]
        ];
        varComponents = ConnectedComponents[Graph[allVars, edges]];
        varToComp = Association[];
        Do[
            (varToComp[#] = k) & /@ varComponents[[k]],
            {k, Length[varComponents]}
        ];
        compIndex[atomVars_] := If[atomVars === {}, Missing[], varToComp[First[atomVars]]];
        results = Table[
            With[{compVars = varComponents[[k]]},
                {
                    Pick[disjAtoms,   (compIndex[#] === k) & /@ disjVars],
                    Pick[linearAtoms, (compIndex[#] === k) & /@ linVars],
                    compVars
                }
            ],
            {k, Length[varComponents]}
        ];
        results
    ];

Options[booleanMinimizeReduceSystem] = {
    "ArmTimeout"      -> 2,
    "DisjunctTimeout" -> 30,
    "ReturnAll"       -> False
};

booleanMinimizeReduceSystem[sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    withCriticalCongestionSolver[sys, "booleanMinimizeReduceSystem",
        Function[{constraints, allVars},
            Module[{atoms, linearAtoms, disjAtoms, status, components,
                    armTimeout, disjTimeout, returnAll, mergedRules,
                    componentResults, anyFalse = False, infeasible = False},
                armTimeout  = OptionValue[booleanMinimizeReduceSystem, {opts}, "ArmTimeout"];
                disjTimeout = OptionValue[booleanMinimizeReduceSystem, {opts}, "DisjunctTimeout"];
                returnAll   = TrueQ[OptionValue[booleanMinimizeReduceSystem, {opts}, "ReturnAll"]];

                atoms = If[Head[constraints] === And, List @@ constraints, {constraints}];
                {linearAtoms, disjAtoms} = {
                    Select[atoms, Head[#] =!= Or &],
                    Select[atoms, Head[#] === Or &]
                };

                {linearAtoms, disjAtoms, status} =
                    pruneDisjunctiveArms[linearAtoms, disjAtoms, allVars, armTimeout];
                If[status === "Infeasible",
                    Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
                ];

                components = decomposeByVariableSharing[linearAtoms, disjAtoms, allVars];

                componentResults = Map[
                    Function[comp,
                        Module[{compDisjAtoms, compLinear, compVars,
                                dnf, disjuncts, reducedList, nonFalse, parsed,
                                fullSystem, redOne},
                            {compDisjAtoms, compLinear, compVars} = comp;
                            If[compDisjAtoms === {},
                                (* Pure-linear component: one Reduce call. *)
                                redOne = TimeConstrained[
                                    Reduce[
                                        If[compLinear === {}, True, And @@ compLinear],
                                        compVars, Reals
                                    ],
                                    disjTimeout, $Aborted
                                ];
                                If[redOne === False || redOne === $Aborted,
                                    anyFalse = True;
                                    Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
                                ];
                                Return[parseReduceResult[redOne, compVars], Module]
                            ];
                            dnf = BooleanMinimize[And @@ compDisjAtoms, "DNF"];
                            disjuncts = If[Head[dnf] === Or, List @@ dnf, {dnf}];
                            reducedList = TimeConstrained[
                                Reduce[
                                    And @@ Join[compLinear, {#}],
                                    compVars, Reals
                                ],
                                disjTimeout, $Aborted
                            ] & /@ disjuncts;
                            nonFalse = Select[reducedList, # =!= False && # =!= $Aborted &];
                            If[nonFalse === {},
                                anyFalse = True;
                                Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
                            ];
                            parsed = parseReduceResult[#, compVars] & /@ nonFalse;
                            If[Length[parsed] > 1,
                                With[{ruleSets = (If[ListQ[#], #, Lookup[#, "Rules", {}]] & /@ parsed)},
                                    If[!SameQ @@ (Sort /@ ruleSets),
                                        Message[booleanMinimizeReduceSystem::multisol, Length[parsed]]
                                    ]
                                ]
                            ];
                            If[returnAll, parsed, First[parsed]]
                        ]
                    ],
                    components
                ];

                If[anyFalse,
                    Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
                ];

                (* Merge per-component rule sets: components are variable-disjoint. *)
                If[returnAll,
                    componentResults,
                    Module[{rulesList, residuals},
                        rulesList = If[ListQ[#], #, Lookup[#, "Rules", {}]] & /@ componentResults;
                        residuals = Cases[componentResults,
                            a_Association :> Lookup[a, "Residual", True]];
                        mergedRules = Flatten[rulesList];
                        If[residuals === {} || And @@ (TrueQ /@ residuals),
                            mergedRules,
                            <|"Rules" -> mergedRules,
                              "Residual" -> Simplify[And @@ residuals]|>
                        ]
                    ]
                ]
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
                            <|"Rules" -> {}, "Residual" -> False|>
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
                    <|"Rules" -> {}, "Residual" -> False|>
                ]
            ]
        ],
        buildSolverInputs,
        attachAccumulatedRules
    ];

directCriticalSystem[sys_?mfgSystemQ] :=
    Module[{eqs, exits, switching, exitCosts, sol},
        If[!criticalCongestionSystemQ[sys],
            Message[MFGraphs::noncritical, "directCriticalSystem"];
            Return[Failure["directCriticalSystem", <|"Reason" -> "NonCritical"|>]]
        ];
        
        exits = systemData[sys, "Exits"];
        exitCosts = If[ListQ[exits], exits[[All, 2]], {}];
        If[!AllTrue[exitCosts, NumericQ] || (Length[DeleteDuplicates[exitCosts]] > 1),
            Return[Failure["directCriticalSystem", <|"Tag" -> "directCriticalSystem", "Reason" -> "ExitCostsNotEqualNumeric"|>]]
        ];
        
        switching = Values[systemData[sys, "SwitchingCosts"]];
        If[MissingQ[switching] || !AllTrue[switching, NumericQ] || AnyTrue[switching, # != 0 &],
            Return[Failure["directCriticalSystem", <|"Tag" -> "directCriticalSystem", "Reason" -> "NonZeroOrNonNumericSwitchingCosts"|>]]
        ];

        eqs = sysToEqsInternal[sys];
        sol = directCriticalSolverInternal[eqs];
        If[AssociationQ[sol],
            sol,
            Failure["directCriticalSystem", <|"Reason" -> "SolveFailed"|>]
        ]
    ];

flowFirstCriticalSystem[sys_?mfgSystemQ] :=
    Module[{eqs, numericState, sol},
        If[!criticalCongestionSystemQ[sys],
            Message[MFGraphs::noncritical, "flowFirstCriticalSystem"];
            Return[Failure["flowFirstCriticalSystem", <|"Reason" -> "NonCritical"|>]]
        ];
        eqs = sysToEqsInternal[sys];
        numericState = Lookup[eqs, "NumericState", <||>];
        sol = SolveCriticalJFirstBackend[eqs, numericState];
        If[AssociationQ[sol],
            sol,
            Failure["flowFirstCriticalSystem", <|"Reason" -> "SolveFailed"|>]
        ]
    ];

trackedVarQ[var_] :=
    MatchQ[var, j[__] | u[__]];

collectTrackedVars[sys_?mfgSystemQ, expr_] :=
    Select[
        Join[
            systemData[sys, "Us"],
            systemData[sys, "Js"],
            systemData[sys, "Jts"]
        ],
        !FreeQ[expr, #] &
    ];

primarySolutionVarQ[var_] :=
    MatchQ[var, j[_, _] | u[__]];

transitionFlowVarQ[var_] :=
    MatchQ[var, j[_, _, _]];

residualPrimaryVariables[residual_] :=
    DeleteDuplicates @ Cases[residual, (j[_, _] | u[__]), {0, Infinity}];

residualTransitionFlows[residual_] :=
    DeleteDuplicates @ Cases[residual, j[_, _, _], {0, Infinity}];

residualTrackedVariables[residual_] :=
    DeleteDuplicates @ Join[residualPrimaryVariables[residual], residualTransitionFlows[residual]];

orInvolvingPrimaryVariableQ[residual_] :=
    !FreeQ[
        residual,
        expr_Or /; residualPrimaryVariables[expr] =!= {}
    ];

solutionRules[sol_] :=
    If[ListQ[sol], sol, If[AssociationQ[sol], Lookup[sol, "Rules", {}], {}]];

dnfResidualDiagnostics[sol_] :=
    Module[{residual, topLevelBranches, nestedOrQ, primaryVars, transitionVars, trackedVars},
        residual = If[AssociationQ[sol], Lookup[sol, "Residual", True], True];
        topLevelBranches = If[Head[residual] === Or, Length[List @@ residual], 0];
        nestedOrQ = If[Head[residual] === Or,
            !FreeQ[List @@ residual, Or],
            !FreeQ[residual, Or]
        ];
        primaryVars = residualPrimaryVariables[residual];
        transitionVars = residualTransitionFlows[residual];
        trackedVars = residualTrackedVariables[residual];
        <|
            "TopLevelBranchCount" -> topLevelBranches,
            "NestedOrQ" -> nestedOrQ,
            "TrackedVariablesQ" -> (trackedVars =!= {}),
            "TrackedVariables" -> trackedVars,
            "PrimaryVariablesQ" -> (primaryVars =!= {}),
            "PrimaryVariables" -> primaryVars,
            "TransitionFlowVariablesQ" -> (transitionVars =!= {}),
            "TransitionFlowVariables" -> transitionVars
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
                residual = Lookup[sol, "Residual", True];
                Which[
                    residual === True, "Rules",
                    residual === False, "NoSolution",
                    orInvolvingPrimaryVariableQ[residual], "Branched",
                    True,
                        diagnostics = dnfResidualDiagnostics[sol];
                        Which[
                            TrueQ[diagnostics["PrimaryVariablesQ"]], "Underdetermined",
                            TrueQ[diagnostics["TransitionFlowVariablesQ"]], "Rules",
                            True, "ResidualLogic"
                        ]
                ],
            True,
                "Unknown"
        ]
    ];

solutionVariableDiagnostics[sys_?mfgSystemQ, sol_] :=
    Module[{rules, residual, transitionFlows, transitionRuleFlows,
            residualTransitions, residualTransitionSet, determinedTransitions},
        rules = solutionRules[sol];
        residual = If[AssociationQ[sol], Lookup[sol, "Residual", True], True];
        residual = residual /. rules;
        transitionFlows = Replace[systemData[sys, "Jts"], Except[_List] -> {}];
        transitionRuleFlows = DeleteDuplicates @ Cases[rules, Rule[lhs_?transitionFlowVarQ, _] :> lhs];
        residualTransitions = Select[residualTransitionFlows[residual], MemberQ[transitionFlows, #] &];
        residualTransitionSet = AssociationThread[residualTransitions, True];
        determinedTransitions = Select[transitionFlows, !KeyExistsQ[residualTransitionSet, #] &];
        <|
            "PrimaryResultKind" -> solutionResultKind[sol],
            "PrimaryResidualVariables" -> residualPrimaryVariables[residual],
            "TransitionFlowStatus" -> If[residualTransitions === {}, "Unique", "Underdetermined"],
            "TransitionFlowCount" -> Length[transitionFlows],
            "TransitionFlowRuleCount" -> Length[Select[transitionRuleFlows, MemberQ[transitionFlows, #] &]],
            "TransitionFlowResidualCount" -> Length[residualTransitions],
            "DeterminedTransitionFlows" -> determinedTransitions,
            "ResidualTransitionFlows" -> residualTransitions
        |>
    ];

solutionReport[sys_?mfgSystemQ, sol_] :=
    Module[{kind, validation, residual, varDiagnostics, dnfDiagnostics},
        kind = solutionResultKind[sol];
        validation = isValidSystemSolution[sys, sol, "ReturnReport" -> True];
        residual = computeKirchhoffResidual[sys, sol];
        varDiagnostics = solutionVariableDiagnostics[sys, sol];
        dnfDiagnostics = dnfResidualDiagnostics[sol];
        <|
            "ResultKind" -> kind,
            "ValidationReport" -> validation,
            "KirchhoffResidual" -> residual,
            "TransitionFlowDiagnostics" -> varDiagnostics,
            "ResidualTransitionFlows" -> Lookup[varDiagnostics, "ResidualTransitionFlows", {}],
            "DNFDiagnostics" -> dnfDiagnostics
        |>
    ];

solutionBranchCostReport[sys_?mfgSystemQ, sol_] :=
    Module[{rules, residual, branches, branchData, sortedBranches},
        rules = solutionRules[sol];
        residual = If[AssociationQ[sol], Lookup[sol, "Residual", True], True];
        branches = If[Head[residual] === Or, List @@ residual, {residual}];
        
        branchData = MapIndexed[
            Function[{branch, idx},
                Module[{branchRules, branchSol, exitCost, switchCost, edgeFlowCost, totalObjective, validation},
                    branchRules = branchToRules[branch, residualTrackedVariables[branch]];
                    branchSol = normalizeRules[mergeRules[rules, branchRules]];
                    
                    (* Simplified diagnostic objective: Exit costs + Switching costs + 0.5 * Sum[j^2] *)
                    exitCost = Total[Cases[branchSol, Rule[u[_, v_], val_?NumericQ] /; MemberQ[systemData[sys, "Exits"][[All, 1]], v] :> val]];
                    switchCost = Total[Cases[branchSol, Rule[j[r_, i_, w_], val_?NumericQ] :> val * switchingCostLookup[systemData[sys, "SwitchingCosts"]][r, i, w]]];
                    edgeFlowCost = 0.5 * Total[Cases[branchSol, Rule[j[_, _], val_?NumericQ] :> val^2]];
                    totalObjective = N[exitCost + switchCost + edgeFlowCost];
                    
                    validation = isValidSystemSolution[sys, <|"Rules" -> branchSol, "Residual" -> True|>];
                    
                    <|
                        "Index" -> idx[[1]],
                        "TotalObjective" -> If[NumericQ[totalObjective], totalObjective, 10^10],
                        "ExitCost" -> exitCost,
                        "SwitchingCost" -> switchCost,
                        "EdgeFlowQuadraticCost" -> edgeFlowCost,
                        "DeterminedVariables" -> Keys[branchRules],
                        "ResidualVariables" -> residualTrackedVariables[branch],
                        "Valid" -> validation,
                        "Expression" -> branch
                    |>
                ]
            ],
            branches
        ];
        
        sortedBranches = SortBy[branchData, #TotalObjective &];
        
        <|
            "BranchCount" -> Length[branches],
            "BestBranches" -> Take[sortedBranches, UpTo[3]],
            "Branches" -> branchData
        |>
    ];

buildSolverInputs[sys_?mfgSystemQ] :=
    Module[{data, initRules, baseConstraints, allVars, constraints, rulesAcc},
        data = systemDataFlatten[sys];
        (* Seed with two rule families before the DNF stage:
           1. RuleBalanceGatheringFlows: j[a,b] -> Sum j[a,b,c] eliminates all
              non-transition edge flows upfront so residuals are j[a,b,c]-only.
           2. ZeroSwitchUEqualities: zero-cost symmetric triple pairs force
              u-equalities by antisymmetry of >=; applied uniformly here. *)
        initRules = Join[
            Normal @ Lookup[data, "RuleBalanceGatheringFlows", <||>],
            Lookup[data, "ZeroSwitchUEqualities", {}]
        ];
        baseConstraints = And[
            And @@ Lookup[data, "EqEntryIn"],
            Lookup[data, "EqBalanceSplittingFlows"],
            Lookup[data, "EqBalanceGatheringFlows"],
            And @@ Lookup[data, "EqGeneral"],
            And @@ Lookup[data, "IneqJs"],
            And @@ Lookup[data, "IneqJts"],
            And @@ Lookup[data, "IneqSwitchingByVertex"],
            And @@ Lookup[data, "IneqExitValues"],
            And @@ Lookup[data, "AltFlows"],
            And @@ Lookup[data, "AltTransitionFlows"],
            And @@ Lookup[data, "AltOptCond"],
            And @@ Lookup[data, "AltExitCond"]
        ];
        allVars = collectTrackedVars[sys, baseConstraints /. initRules];
        {constraints, rulesAcc} = accumulateEqualityRules[baseConstraints, allVars, initRules];
        (* Eliminated vars are LHS of rulesAcc; set-subtract instead of re-walking. *)
        allVars = DeleteCases[allVars, Alternatives @@ (First /@ rulesAcc)];
        {constraints, allVars, rulesAcc}
    ];

flattenConjuncts::usage = "flattenConjuncts[expr] unpacks expr into a conjunct list for iterative processing.";
substituteSolution::usage = "substituteSolution[rst, sol] substitutes solution rules sol into rst.";

constraintPriority::usage =
"constraintPriority[expr] returns an integer priority (0 for Equal, 1 for non-Or, 2 for Or) used for ordering conjuncts.";

orderedConjuncts::usage =
"orderedConjuncts[expr] returns a list of conjuncts from expr sorted by constraintPriority and string form.";

orderedAlternatives::usage =
"orderedAlternatives[expr] returns the list of disjuncts from an Or expression sorted by string form.";

dnfTopLevelConjuncts::usage =
"dnfTopLevelConjuncts[expr] returns the top-level conjuncts of a boolean expression.";

dnfVarVertices::usage =
"dnfVarVertices[var] returns the vertex indices associated with a solver variable.";

dnfExpressionTrackedVariables::usage =
"dnfExpressionTrackedVariables[expr] returns the list of all tracked variables (j[__], u[__]) in an expression.";

dnfExpressionVertices::usage =
"dnfExpressionVertices[expr] returns the list of all vertex indices mentioned in an expression's variables.";

dnfTrackedVariableKind::usage =
"dnfTrackedVariableKind[var] returns a string describing the variable type (TransitionFlow, EdgeFlow, Value, or Other).";

dnfPathInfo::usage =
"dnfPathInfo[sys] returns an association with path and distance metadata used for path-aware DNF ordering.";

dnfPathRank::usage =
"dnfPathRank[expr, pathInfo] returns a ranking tuple for an expression based on its distance to optimal paths.";

dnfOrderingSummary::usage =
"dnfOrderingSummary[conjuncts] returns diagnostic counts of conjunct types (Equal, Or, etc.) for a list of conjuncts.";

dnfOrderConjuncts::usage =
"dnfOrderConjuncts[expr, order, sys] reorders the conjuncts of expr according to the specified strategy.";

dnfOrderExpression::usage =
"dnfOrderExpression[expr, order, sys] reorders and reassembles a boolean expression according to the specified strategy.";

branchKey::usage =
"branchKey[branch] returns a canonical string key for a branch association for deduplication.";

deduplicateBranches::usage =
"deduplicateBranches[branches] removes duplicate branch associations based on their canonical keys.";

trackedFreeQ::usage =
"trackedFreeQ[expr] returns True if the expression contains no tracked solver variables.";

cheapConstraintSimplify::usage =
"cheapConstraintSimplify[expr] performs lightweight simplification of boolean constraints.";

(* linearRuleFromEquality removed: subsumed by rulesFromEqualities[{eq}, vars, existing]
   under makeSystem's linearity guarantee. *)

simplifyResiduals::usage =
"simplifyResiduals[residuals, rules] substitutes rules into a list of residuals and simplifies.";

addBranchConstraint::usage =
"addBranchConstraint[branch, elem, vars] incorporates a new constraint into an existing branch, potentially splitting it.";

branchStateReduce::usage =
"branchStateReduce[constraints, allVars] reduces a constraint system into a list of branch associations.";

branchStateReduceFromBranches::usage =
"branchStateReduceFromBranches[branches, constraints, allVars] continues branch-state reduction from an existing set of branches.";

groundRuleQ::usage =
"groundRuleQ[rule, allVars] returns True if the RHS of a rule is free of any tracked variables.";

splitGroundRules::usage =
"splitGroundRules[rules, allVars] splits a rule list into ground (numeric/constant) rules and symbolic residual equalities.";

branchStateFinalizeResult::usage =
"branchStateFinalizeResult[branches, allVars] converts a list of branch associations into the final solver result format.";

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
    (* Canonical-order tiebreaker via Sort instead of ToString[..., InputForm] —
       both inputs evaluate through the same kernel so structural order suffices. *)
    SortBy[flattenConjuncts[expr], {constraintPriority, Identity}];

orderedAlternatives[expr_Or] := Sort[List @@ expr];
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
        (* Hash is a stable totally-ordered tiebreaker; avoids InputForm pretty-print cost. *)
        {onPath, minDistance, exitDistance, entryDistance, LeafCount[Unevaluated[expr]], Hash[Unevaluated[expr]]}
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
    (* Hash of canonically-sorted Rules + Residuals replaces a full InputForm
       pretty-print. DeleteDuplicatesBy only requires equal keys for equal
       branches, which Hash satisfies for structurally-equal inputs. *)
    Hash[{
        Sort @ Lookup[branch, "Rules", {}],
        Sort @ Lookup[branch, "Residuals", {}]
    }];

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

simplifyResiduals[residuals_List, rules_List] :=
    DeleteCases[cheapConstraintSimplify /@ (residuals /. rules), True];

addBranchConstraint[branch_Association, elem_, vars_List] :=
    Module[{rules, residuals, simplified, solvedRules, newRules, newResiduals},
        rules = Lookup[branch, "Rules", {}];
        residuals = Lookup[branch, "Residuals", {}];
        simplified = cheapConstraintSimplify[elem /. rules];
        Which[
            simplified === True,
                {branch},
            simplified === False,
                {},
            Head[simplified] === And,
                (* Batch all top-level Equal conjuncts through one rulesFromEqualities
                   call instead of recursing per conjunct. RowReduce is complete on
                   linear systems, so one CoefficientArrays + RowReduce subsumes N
                   per-equation calls. The non-equation conjuncts (and any equations
                   that survive batch substitution as non-trivial) fall through the
                   per-conjunct fold for contradiction detection. *)
                Module[{conjuncts, eqs, rest, batchedRules, mergedRules,
                        seedResiduals, seedBranch, residualEqs},
                    conjuncts = orderedConjuncts[simplified];
                    eqs  = Cases[conjuncts, _Equal];
                    rest = DeleteCases[conjuncts, _Equal];
                    batchedRules = If[eqs === {}, {}, rulesFromEqualities[eqs, vars, rules]];
                    If[batchedRules === {},
                        Fold[
                            Function[{branches, conjunct},
                                Join @@ (addBranchConstraint[#, conjunct, vars] & /@ branches)
                            ],
                            {branch},
                            conjuncts
                        ],
                        mergedRules   = normalizeRules[mergeRules[rules, batchedRules]];
                        seedResiduals = simplifyResiduals[residuals, batchedRules];
                        If[MemberQ[seedResiduals, False],
                            {},
                            seedBranch  = <|"Rules" -> mergedRules, "Residuals" -> seedResiduals|>;
                            (* Drop equations whose batched substitution makes them True
                               (pivoted or redundant); keep the rest for the fold to
                               either consume into residuals or detect as contradictions. *)
                            residualEqs = DeleteCases[
                                cheapConstraintSimplify /@ (eqs /. batchedRules),
                                True
                            ];
                            If[MemberQ[residualEqs, False],
                                {},
                                Fold[
                                    Function[{branches, conjunct},
                                        Join @@ (addBranchConstraint[#, conjunct, vars] & /@ branches)
                                    ],
                                    {seedBranch},
                                    Join[residualEqs, rest]
                                ]
                            ]
                        ]
                    ]
                ],
            Head[simplified] === Or,
                Join @@ (addBranchConstraint[branch, #, vars] & /@ orderedAlternatives[simplified]),
            Head[simplified] === Equal,
                solvedRules = rulesFromEqualities[{simplified}, vars, rules];
                If[solvedRules =!= {},
                    newRules = normalizeRules[mergeRules[rules, solvedRules]];
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
            branchRules, branchGround, branchSymbolic, residuals, residualExpr,
            finalResidual, finalRules},
        If[branches === {},
            Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
        ];
        ruleSets = Lookup[#, "Rules", {}] & /@ branches;
        common = commonRulesFromBranches[ruleSets];
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
        (* Promote top-level equalities surviving in residualExpr into rules
           (ground or parametric — the package treats both the same; see
           RuleBalanceGatheringFlows in systemTools.wl which seeds initRules
           with parametric rules). The per-branch reducer leaves equalities in
           Residuals rather than Rules; without this pass a fully-solved result
           misreports as Underdetermined. *)
        {finalResidual, finalRules} = accumulateEqualityRules[residualExpr, allVars, commonGround];
        If[finalResidual === True,
            finalRules,
            <|"Rules" -> finalRules, "Residual" -> finalResidual|>
        ]
    ];

branchStateReduceResult[constraints_, allVars_List] :=
    branchStateFinalizeResult[branchStateReduce[constraints, allVars], allVars];

(* --- cachedSolve: memoized Quiet[Solve[expr], Solve::svars] ------------- *)
(* Used inside dnfReduce / dnfReduceProcedural / dnfReduceInstrumented.
   On a baseline measurement (2026-05-11) grid-4x4 made 71 051 non-trivial
   Solve calls on only 114 distinct keys (99.8% hit rate). The cache is
   keyed by Hash[expr] (32-bit, fast); collisions are negligible at our
   scale. Reset before each dnfReduceSystem invocation. *)

$dnfSolveCache    = <||>;
$dnfSolveHits     = 0;
$dnfSolveMiss     = 0;
$dnfSolveMissTag  = Unique["dnfSolveMissTag$"];  (* sentinel — guaranteed not to collide with any Solve result *)

resetSolveCacheCounters[] := (
    $dnfSolveHits = 0;
    $dnfSolveMiss = 0;
);

clearSolveCache[] := ($dnfSolveCache = <||>);

cachedSolve[expr_] :=
    Module[{h = Hash[expr], cached},
        cached = Lookup[$dnfSolveCache, h, $dnfSolveMissTag];
        If[cached === $dnfSolveMissTag,
            $dnfSolveMiss++;
            cached = Quiet[Solve[expr], Solve::svars];
            $dnfSolveCache[h] = cached;
            ,
            $dnfSolveHits++;
        ];
        cached
    ];

printSolveCacheSummary[label_String] :=
    Module[{total, hitRate},
        total   = $dnfSolveHits + $dnfSolveMiss;
        hitRate = If[total === 0, 0., N[$dnfSolveHits / total]];
        Print["SOLVE_CACHE|", label,
              "|total=", total,
              "|hits=",  $dnfSolveHits,
              "|miss=",  $dnfSolveMiss,
              "|hit_rate=", hitRate,
              "|cache_size=", Length[$dnfSolveCache]];
    ];

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
                    sol  = cachedSolve[elem];
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
        sol  = cachedSolve[fst];
        fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
        newxp = substituteSolution[xp, fsol];
        If[newxp === False,
            False,
            newrst = substituteSolution[rst, fsol];
            dnfReduce[newxp && fst, newrst]
        ]
    ]

dnfReduce[xp_, rst_, fst_] := dnfReduce[xp && fst, rst]

dnfReduceProcedural::usage =
"dnfReduceProcedural[xp, sys, allVars] is a stack-based iterative equivalent \
of dnfReduce. Uses an explicit worklist instead of recursive calls, avoiding \
$RecursionLimit on deeply nested Or-chains. Produces the same output as \
dnfReduce for well-formed inputs. Or-branches are pushed in original order \
(Prepend + Reverse) so branch priority matches the recursive version. \
Not wired into dnfReduceSystem by default; call directly for deep Or-chains.";

dnfReduceProcedural[xp0_, sys0_, allVars_List] :=
    Module[{stack = {<|"xp" -> xp0, "rst" -> True, "fst" -> sys0|>},
            results = {}, item, xp, rst, fst,
            conjuncts, xpAcc, i, elem, sol, fsol, newxp, rest, newRest},
        While[stack =!= {},
            item  = First[stack];
            stack = Rest[stack];
            xp = item["xp"]; rst = item["rst"]; fst = item["fst"];
            Which[
                xp  === False || fst === False,
                    Null,
                fst === True,
                    AppendTo[results, xp && rst],
                Head[fst] === Or,
                    Do[stack = Prepend[stack, <|"xp" -> xp, "rst" -> rst, "fst" -> b|>],
                       {b, Reverse[List @@ fst]}],
                Head[fst] === And,
                    conjuncts = List @@ fst;
                    xpAcc     = xp;
                    i         = 1;
                    Module[{outcome = "normal"},
                        While[i <= Length[conjuncts] && outcome === "normal",
                            If[xpAcc === False, outcome = "false"; Break[]];
                            elem = conjuncts[[i]];
                            Which[
                                elem === True,
                                    i++,
                                elem === False,
                                    outcome = "false",
                                Head[elem] === Or,
                                    rest = If[i >= Length[conjuncts], rst,
                                              And[rst, And @@ conjuncts[[i + 1 ;;]]]];
                                    Do[stack = Prepend[stack, <|"xp" -> xpAcc, "rst" -> rest, "fst" -> b|>],
                                       {b, Reverse[List @@ elem]}];
                                    outcome = "branched",
                                Head[elem] === Equal,
                                    sol  = cachedSolve[elem];
                                    fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
                                    newxp = substituteSolution[xpAcc, fsol];
                                    If[newxp === False,
                                        outcome = "false",
                                        xpAcc = newxp && elem;
                                        If[i < Length[conjuncts],
                                            newRest = substituteSolution[And @@ conjuncts[[i + 1 ;;]], fsol];
                                            conjuncts = Join[conjuncts[[;; i]], flattenConjuncts[newRest]]
                                        ]
                                    ];
                                    i++,
                                True,
                                    xpAcc = xpAcc && elem;
                                    i++
                            ]
                        ];
                        If[outcome === "normal",
                            AppendTo[results, xpAcc && rst]]
                    ],
                Head[fst] === Equal,
                    sol  = cachedSolve[fst];
                    fsol = If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}];
                    newxp = substituteSolution[xp, fsol];
                    If[newxp =!= False,
                        stack = Prepend[stack,
                            <|"xp" -> newxp && fst, "rst" -> True,
                              "fst" -> substituteSolution[rst, fsol]|>]],
                True,
                    AppendTo[results, xp && rst && fst]
            ]
        ];
        Which[
            results === {}, False,
            Length[results] === 1, First[results],
            True, Or @@ results
        ]
    ]

dnfTraceEvent::usage =
"dnfTraceEvent[monitor, event] logs an event to the DNF instrumentation monitor.";

dnfCounterIncrement::usage =
"dnfCounterIncrement[monitor, key, n] increments a diagnostic counter in the DNF monitor.";

dnfReduceInstrumented::usage =
"dnfReduceInstrumented[xp, sys, monitor, order, sysObj] is the instrumented \
recursive engine used by dnfReduceDiagnosticReport. Directly recursive with \
an explicit depth counter: each Or-branch and Equal-continuation increments \
depth by 1. Mirrors dnfReduce with branch-ordering and monitoring side effects \
added. Non-reentrant: writes into $dnfCurrentMonitor via mutable assignments.";

dnfReduceInstrumentedAnd::usage =
"dnfReduceInstrumentedAnd[xp, conjuncts, monitor, order, sysObj, depth] implements the And-case for the instrumented DNF reducer.";

sysToEqsInternal::usage =
"sysToEqsInternal[sys] extracts a flat association of all relevant system data for internal solver backends.";

directCriticalSolverInternal::usage =
"directCriticalSolverInternal[Eqs] implements the core logic for directCriticalSystem using graph distance heuristics.";

SolveCriticalJFirstBackend::usage =
"SolveCriticalJFirstBackend[Eqs, numericState] implements the core logic for flowFirstCriticalSystem.";

BuildLinearSystemFromEqualities::usage =
"BuildLinearSystemFromEqualities[equalities, vars] extracts a coefficient matrix and RHS vector from a list of linear equalities.";

LinearSolveCandidate::usage =
"LinearSolveCandidate[aMat, bVec] computes a least-squares candidate solution for a linear system.";

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
                    sol  = cachedSolve[elem];
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
        sol  = cachedSolve[fst];
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
            tDNF, reduced, status, parsed, result, variableDiagnostics, summary},
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
            "OrderedSameMultisetQ" -> (Sort[conjuncts] === Sort[dnfTopLevelConjuncts[orderedConstraints]])
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
        variableDiagnostics = If[status === "OK",
            solutionVariableDiagnostics[sys, result],
            <|
                "PrimaryResultKind" -> "TIMEOUT",
                "PrimaryResidualVariables" -> {},
                "TransitionFlowStatus" -> Missing["TimedOut"],
                "TransitionFlowCount" -> Length[Replace[systemData[sys, "Jts"], Except[_List] -> {}]],
                "TransitionFlowRuleCount" -> Missing["TimedOut"],
                "TransitionFlowResidualCount" -> Missing["TimedOut"],
                "DeterminedTransitionFlows" -> Missing["TimedOut"],
                "ResidualTransitionFlows" -> Missing["TimedOut"]
            |>
        ];
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
            "ResultKind" -> If[status === "OK", solutionResultKind[result], "TIMEOUT"],
            "PrimaryResultKind" -> Lookup[variableDiagnostics, "PrimaryResultKind", Missing["Unavailable"]],
            "PrimaryResidualVariables" -> Lookup[variableDiagnostics, "PrimaryResidualVariables", Missing["Unavailable"]],
            "TransitionFlowStatus" -> Lookup[variableDiagnostics, "TransitionFlowStatus", Missing["Unavailable"]],
            "TransitionFlowCount" -> Lookup[variableDiagnostics, "TransitionFlowCount", Missing["Unavailable"]],
            "TransitionFlowRuleCount" -> Lookup[variableDiagnostics, "TransitionFlowRuleCount", Missing["Unavailable"]],
            "TransitionFlowResidualCount" -> Lookup[variableDiagnostics, "TransitionFlowResidualCount", Missing["Unavailable"]]
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

rulesFromEqualities[equations_List, vars_List, existingRules_List] :=
    Module[{n = Length[vars], arrays, c0, c1, augmented, rref, candidateRules, existingLHS},
        If[equations === {} || n === 0, Return[{}, Module]];
        (* Equations represent c1.vars + c0 == 0; RowReduce on the augmented
           matrix [c1 | -c0] reads off solved rules in one pass. Safe because
           (a) makeSystem guarantees linearity in {j, u, z} and (b) eqGeneral
           no longer carries Cost[__] coefficients, so coefficients here are
           pure rationals. Pivot priority is column order; vars arrive u-first
           via collectTrackedVars so RowReduce eliminates value variables ahead
           of flow variables, keeping Alt* disjunctions collapsible downstream. *)
        arrays = Quiet @ Check[
            CoefficientArrays[equations /. Equal[a_, b_] :> a - b, vars],
            $Failed
        ];
        If[!ListQ[arrays] || Length[arrays] < 2, Return[{}, Module]];
        {c0, c1} = Normal /@ Take[arrays, 2];
        augmented = ArrayFlatten[{{c1, Transpose[{-c0}]}}];
        rref = RowReduce[augmented];
        candidateRules = Table[
            Module[{row = rref[[i]], leadIdx},
                leadIdx = SelectFirst[Range[n], row[[#]] =!= 0 &, 0];
                If[leadIdx === 0,
                    Nothing,
                    vars[[leadIdx]] ->
                        -Sum[row[[k]] vars[[k]], {k, leadIdx + 1, n}] + row[[n + 1]]
                ]
            ],
            {i, Length[rref]}
        ];
        (* O(1) LHS membership via Association instead of O(n) FreeQ scan per candidate. *)
        existingLHS = If[existingRules === {}, <||>,
            AssociationThread[existingRules[[All, 1]] -> True]];
        Select[candidateRules, !KeyExistsQ[existingLHS, First[#]] &]
    ];

(* RowReduce is complete on a linear system, so one rulesFromEqualities call
   extracts every implied rule from the current top-level equation set in one
   shot. Recursion is needed only when substitution exposes a brand-new
   top-level equality (Or-branch collapse, paired-inequality fusion such as
   a<=b && a>=b -> a==b). Termination: each step's rulesFromEqualities
   LHS-dedup filter forbids re-solving the same variable, so the residual
   equation set strictly shrinks. The depth guard catches a broken linearity
   invariant or non-shrinking constraint that would otherwise loop silently. *)
$accumulateEqualityRulesDepthCap = 8;

accumulateEqualityRules[baseConstraints_, vars_List, initRules_List, depth_Integer:0] :=
    Module[{constraints, eqs, newRules, nextConstraints, nextEqs, mergedRules},
        constraints = baseConstraints /. initRules;
        eqs = topLevelEquations[constraints];
        newRules = If[eqs === {}, {}, rulesFromEqualities[eqs, vars, initRules]];
        If[newRules === {},
            Return[{constraints, normalizeRules[initRules]}, Module]
        ];
        nextConstraints = constraints /. newRules;
        nextEqs = topLevelEquations[nextConstraints];
        mergedRules = mergeRules[initRules, newRules];
        Which[
            Complement[nextEqs, eqs] === {},
                {nextConstraints, normalizeRules[mergedRules]},
            depth >= $accumulateEqualityRulesDepthCap,
                mfgPrint["accumulateEqualityRules: depth cap hit; check linearity invariant."];
                {nextConstraints, normalizeRules[mergedRules]},
            True,
                accumulateEqualityRules[nextConstraints, vars, mergedRules, depth + 1]
        ]
    ];

mergeAccumulatedRules[rulesAcc_, rules_] := normalizeRules[mergeRules[rulesAcc, rules]];

(* Rules common to every branch's rule set. Empty input -> {}. *)
commonRulesFromBranches[ruleSets_List] :=
    If[ruleSets === {}, {}, Fold[Intersection, First[ruleSets], Rest[ruleSets]]];

attachAccumulatedRules[result_, rulesAcc_List] /; ListQ[result] :=
    mergeAccumulatedRules[rulesAcc, result];

attachAccumulatedRules[result_, rulesAcc_List] /; AssociationQ[result] :=
    Association[result, "Rules" -> mergeAccumulatedRules[rulesAcc, Lookup[result, "Rules", {}]]];

attachAccumulatedRules[result_, _] := result;

(* O(1) KeyExistsQ replaces per-match O(n) MemberQ scan. branchToRules
   accepts an optional pre-built (allVarsSet, varsAlt) pair so the Or-branch
   fan-out in parseReduceResult builds them once instead of per branch. *)
branchToRules[branch_, allVars_List] :=
    branchToRules[branch, AssociationThread[allVars -> True], Alternatives @@ allVars];

branchToRules[branch_, allVarsSet_Association, varsAlt_] :=
    With[{conjuncts = Switch[Head[branch], And, List @@ branch, _, {branch}]},
        Cases[conjuncts,
            HoldPattern[v_ == val_] /;
                KeyExistsQ[allVarsSet, v] && FreeQ[val, varsAlt] :>
                (v -> val)
        ]
    ];

parseReduceResult[reduced_, allVars_] :=
    Module[{conjuncts, rules, residual, allVarsSet, varsAlt},
        allVarsSet = AssociationThread[allVars -> True];
        varsAlt = Alternatives @@ allVars;
        If[Head[reduced] === Or,
            Return[
                With[{common = commonRulesFromBranches[
                                    branchToRules[#, allVarsSet, varsAlt] & /@ List @@ reduced]},
                    <|"Rules" -> common, "Residual" -> reduced|>
                ]
            ]
        ];
        conjuncts = Switch[Head[reduced],
            And,   List @@ reduced,
            True,  {},
            False, Return[<|"Rules" -> {}, "Residual" -> False|>],
            _,     {reduced}
        ];
        rules = Cases[conjuncts,
            HoldPattern[v_ == val_] /;
                KeyExistsQ[allVarsSet, v] && FreeQ[val, varsAlt] :>
                (v -> val)
        ];
        residual = DeleteCases[conjuncts,
            HoldPattern[v_ == val_] /;
                KeyExistsQ[allVarsSet, v] && FreeQ[val, varsAlt]
        ];
        If[residual === {},
            rules,
            <|"Rules" -> rules, "Residual" -> And @@ residual|>
        ]
    ];

harvestDNFBranch[False, _] := False;

harvestDNFBranch[branch_, allVars_List] :=
    Module[{conjuncts, eqs, sol, rules, reducedConjuncts, residual, allVarsSet, varsAlt},
        conjuncts = flattenConjuncts[branch];
        eqs = Select[conjuncts, MatchQ[#, _Equal] &];
        sol = If[eqs === {} || allVars === {},
            {{}},
            Quiet[Solve[eqs, allVars], Solve::svars]
        ];
        If[!ListQ[sol] || sol === {}, Return[False, Module]];
        allVarsSet = AssociationThread[allVars -> True];
        varsAlt = Alternatives @@ allVars;
        rules = Cases[First[sol],
            r : Rule[v_, rhs_] /;
                KeyExistsQ[allVarsSet, v] && FreeQ[rhs, varsAlt] :> r
        ];
        reducedConjuncts = Simplify /@ (conjuncts /. rules);
        reducedConjuncts = DeleteCases[reducedConjuncts, True];
        residual = If[reducedConjuncts === {}, True, And @@ reducedConjuncts];
        residual = Simplify[residual];
        (* A residual that is neither True nor False is live by definition.
           Variables[...] cannot see inequality atoms, so do not use it here. *)
        Which[
            residual === False,
                False,
            residual === True,
                <|"Rules" -> rules, "Residual" -> True|>,
            True,
                <|"Rules" -> rules, "Residual" -> residual|>
        ]
    ];

parseDNFReduceResult[reduced_, allVars_List] :=
    Module[{branches, harvested, ruleSets, common, branchExprs},
        If[reduced === False, Return[<|"Rules" -> {}, "Residual" -> False|>, Module]];
        branches = If[Head[reduced] === Or, List @@ reduced, {reduced}];
        harvested = DeleteCases[harvestDNFBranch[#, allVars] & /@ branches, False];
        If[harvested === {},
            Return[<|"Rules" -> {}, "Residual" -> False|>, Module]
        ];
        If[Length[harvested] === 1,
            With[{rules = Lookup[First[harvested], "Rules", {}],
                  residual = Lookup[First[harvested], "Residual", True]},
                Return[
                    If[residual === True,
                        rules,
                        <|"Rules" -> rules, "Residual" -> residual|>
                    ],
                    Module
                ]
            ]
        ];
        ruleSets = Lookup[#, "Rules", {}] & /@ harvested;
        common = commonRulesFromBranches[ruleSets];
        branchExprs = Simplify /@ Map[
            With[{branchRules = Complement[Lookup[#, "Rules", {}], common],
                  residual = Lookup[#, "Residual", True]},
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
              "Residual" -> Simplify @ If[Length[branchExprs] === 1, First[branchExprs], Or @@ branchExprs]|>
        ]
    ];

(* --- Solution validation --- *)

Options[isValidSystemSolution] = {
    "Tolerance"    -> 10^-6,
    "ReturnReport" -> False
};

isValidSystemSolution[sys_?mfgSystemQ, sol_, OptionsPattern[]] :=
    Module[{tol, returnReportQ, kind, rules,
            blocks, blockResults, concretelyFailed, overall, report},
        tol          = OptionValue["Tolerance"];
        returnReportQ = TrueQ[OptionValue["ReturnReport"]];

        kind = solutionResultKind[sol];

        If[!MemberQ[{"Rules", "Branched", "Underdetermined", "ResidualLogic", "NoSolution"}, kind],
            Return[If[returnReportQ,
                <|"Valid" -> False, "Kind" -> kind,
                  "Reason" -> "UnrecognizedSolutionFormat", "BlockChecks" -> <||>|>,
                False], Module]
        ];

        rules = If[ListQ[sol], sol, Lookup[sol, "Rules", {}]];

        (* For residual-bearing results: fail fast if residual equations are inconsistent. *)
        If[AssociationQ[sol],
            With[{eqs = Lookup[sol, "Residual", True]},
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
            "EqGeneral"               -> And @@ systemData[sys, "EqGeneral"],
            "IneqJs"                  -> And @@ systemData[sys, "IneqJs"],
            "IneqJts"                 -> And @@ systemData[sys, "IneqJts"],
            "IneqSwitchingByVertex"   -> And @@ systemData[sys, "IneqSwitchingByVertex"],
            "IneqExitValues"          -> With[{v = systemData[sys, "IneqExitValues"]},
                                            If[ListQ[v] && v =!= {}, And @@ v, True]],
            "AltFlows"                -> And @@ systemData[sys, "AltFlows"],
            "AltTransitionFlows"      -> And @@ systemData[sys, "AltTransitionFlows"],
            "AltOptCond"              -> And @@ systemData[sys, "AltOptCond"],
            "AltExitCond"             -> With[{v = systemData[sys, "AltExitCond"]},
                                            If[ListQ[v] && v =!= {}, And @@ v, True]]
        |>;

        blockResults = Association @ KeyValueMap[
            #1 -> checkBlock[#2 /. rules, tol] &,
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

(* --- Internal helpers for direct/flow-first solvers --- *)

sysToEqsInternal[sys_] :=
    Module[{data, Eqs},
        data = systemDataFlatten[sys];
        Eqs = <|
            "Entrance Vertices and Flows" -> systemData[sys, "Entries"],
            "Exit Vertices and Terminal Costs" -> systemData[sys, "Exits"],
            "SwitchingCosts" -> systemData[sys, "SwitchingCosts"],
            "RuleEntryValues" -> Lookup[data, "RuleEntryValues", <||>],
            "IneqExitValues" -> Lookup[data, "IneqExitValues", {}],
            "AltExitCond"    -> Lookup[data, "AltExitCond", {}],
            "RuleBalanceGatheringFlows" -> Lookup[data, "RuleBalanceGatheringFlows", <||>],
            "EqGeneral" -> Lookup[data, "EqGeneral", True],
            "EqEntryIn" -> Lookup[data, "EqEntryIn", True],
            "EqBalanceSplittingFlows" -> Lookup[data, "EqBalanceSplittingFlows", True],
            "EqBalanceGatheringFlows" -> Lookup[data, "EqBalanceGatheringFlows", True],
            "EqValueAuxiliaryEdges" -> Lookup[data, "EqValueAuxiliaryEdges", True],
            "IneqJs" -> Lookup[data, "IneqJs", True],
            "IneqJts" -> Lookup[data, "IneqJts", True],
            "IneqSwitchingByVertex" -> Lookup[data, "IneqSwitchingByVertex", True],
            "AltFlows" -> Lookup[data, "AltFlows", True],
            "AltTransitionFlows" -> Lookup[data, "AltTransitionFlows", True],
            "AltOptCond" -> Lookup[data, "AltOptCond", True],
            "costpluscurrents" -> Lookup[data, "costpluscurrents", <||>],
            "us" -> systemData[sys, "Us"],
            "js" -> systemData[sys, "Js"],
            "jts" -> systemData[sys, "Jts"],
            "auxiliaryGraph" -> Graph[systemData[sys, "AuxVertices"], systemData[sys, "AuxEdges"]],
            "auxExitVertices" -> systemData[sys, "AuxExitVertices"],
            "auxEntryVertices" -> systemData[sys, "AuxEntryVertices"],
            "NumericState" -> systemData[sys, "NumericState"]
        |>;
        Eqs
    ];

directCriticalSolverInternal[Eqs_] :=
    Module[{graph, exitVerts, distToExit, us, js, jts, auxPairs, zeroRules, eqSystem, allVars, result, flowVars, ineqs, inst},
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        If[graph === None || exitVerts === {}, Return[$Failed]];
        
        distToExit = Association @ Table[v -> Min[GraphDistance[graph, v, #] & /@ exitVerts], {v, VertexList[graph]}];
        us = Lookup[Eqs, "us", {}];
        js = Lookup[Eqs, "js", {}];
        jts = Lookup[Eqs, "jts", {}];
        auxPairs = Cases[us, u[a_, b_] :> {a, b}];
        
        zeroRules = Association[];
        Do[
            With[{a = p[[1]], b = p[[2]]},
                If[KeyExistsQ[distToExit, a] && KeyExistsQ[distToExit, b] && distToExit[a] < distToExit[b],
                    If[MemberQ[js, j[a, b]], zeroRules[j[a, b]] = 0]
                ]
            ],
            {p, Join[Cases[js, j[a_, b_] :> {a, b}], Cases[jts, j[a_, b_, c_] :> {a, c}]]}
        ];
            
        eqSystem = And @@ Flatten[{
            (And @@ Flatten[{Lookup[Eqs, "EqGeneral", True]}]) /. Thread[Keys[Lookup[Eqs, "costpluscurrents", <||>]] -> 0],
            Lookup[Eqs, "EqEntryIn", True],
            Lookup[Eqs, "EqBalanceSplittingFlows", True],
            Lookup[Eqs, "EqBalanceGatheringFlows", True],
            Lookup[Eqs, "EqValueAuxiliaryEdges", True],
            KeyValueMap[Equal, Lookup[Eqs, "RuleEntryValues", <||>]],
            Lookup[Eqs, "IneqExitValues", {}],
            Lookup[Eqs, "AltExitCond", {}],
            Equal @@@ Normal[Lookup[Eqs, "RuleBalanceGatheringFlows", <||>]],
            Equal[#, 0] & /@ Keys[zeroRules],
            Module[{byVertex = GroupBy[auxPairs, Last]},
                Flatten @ KeyValueMap[
                    Function[{v, pairs},
                        If[Length[pairs] > 1,
                            (u @@ # == u @@ pairs[[1]]) & /@ Rest[pairs],
                            {}
                        ]
                    ],
                    byVertex
                ]
            ]
        }];
        
        allVars = Join[us, js, jts];
        flowVars = Join[js, jts];
        ineqs = And @@ (# >= 0 & /@ flowVars);
        
        inst = Quiet @ FindInstance[eqSystem && ineqs, allVars, Reals];
        If[MatchQ[inst, {{___Rule}, ___}],
            <|"Rules" -> First[inst], "Residual" -> True|>,
            $Failed
        ]
    ];

(* --- SolveCriticalJFirstBackend and utilities --- *)

SolveCriticalJFirstBackend[Eqs_, numericState_] :=
    Module[{equalities, ineqs, flowVars, allVars, linearSystem, aMat, bVec, xCandidate, inst, allConstraints},
        (* Collect all constraints for FindInstance *)
        allConstraints = And @@ Flatten[{
            Lookup[Eqs, "EqEntryIn", True], 
            Lookup[Eqs, "EqBalanceSplittingFlows", True], 
            Lookup[Eqs, "EqBalanceGatheringFlows", True], 
            Lookup[Eqs, "EqGeneral", True],
            Lookup[Eqs, "EqValueAuxiliaryEdges", True],
            Lookup[Eqs, "IneqJs", True],
            Lookup[Eqs, "IneqJts", True],
            Lookup[Eqs, "AltFlows", True],
            Lookup[Eqs, "AltTransitionFlows", True],
            Lookup[Eqs, "AltOptCond", True],
            Lookup[Eqs, "IneqExitValues", {}],
            Lookup[Eqs, "AltExitCond", {}],
            KeyValueMap[Equal, Lookup[Eqs, "RuleEntryValues", <||>]]
        }];
        
        equalities = Select[flattenConjuncts[allConstraints], MatchQ[#, _Equal] &];
        flowVars = Join[Lookup[Eqs, "js", {}], Lookup[Eqs, "jts", {}]];
        allVars = Join[Lookup[Eqs, "us", {}], flowVars];
        
        linearSystem = BuildLinearSystemFromEqualities[equalities, allVars];
        If[AssociationQ[linearSystem],
            xCandidate = LinearSolveCandidate[linearSystem["A"], linearSystem["b"]];
            If[VectorQ[xCandidate, NumericQ] && AllTrue[xCandidate[[Lookup[AssociationThread[allVars, Range[Length[allVars]]], flowVars]]], # >= -10^-6 &],
                With[{candidateRules = Normal[AssociationThread[allVars, xCandidate]]},
                    If[checkBlock[allConstraints /. candidateRules, 10^-5] === True,
                        Return[<|"Rules" -> candidateRules, "Residual" -> True|>]
                    ]
                ]
            ];
        ];
        
        (* Fallback to FindInstance with full system *)
        inst = Quiet @ FindInstance[allConstraints, allVars, Reals];
        If[MatchQ[inst, {{___Rule}, ___}],
            <|"Rules" -> First[inst], "Residual" -> True|>,
            $Failed
        ]
    ];

BuildLinearSystemFromEqualities[equalities_List, vars_List] :=
    Module[{exprs, placeholders, subExprs, coeffData, bVec, aMat},
        If[equalities === {} || vars === {}, Return[$Failed]];
        exprs = equalities /. Equal[l_, r_] :> (l - r);
        placeholders = Table[Unique["lv"], {Length[vars]}];
        subExprs = exprs /. Thread[vars -> placeholders];
        coeffData = Quiet @ Check[CoefficientArrays[subExprs, placeholders], $Failed];
        If[coeffData === $Failed || Length[coeffData] < 2, Return[$Failed]];
        bVec = Developer`ToPackedArray @ N @ (-Normal[coeffData[[1]]]);
        aMat = SparseArray[coeffData[[2]]];
        <|"A" -> aMat, "b" -> bVec|>
    ];

LinearSolveCandidate[aMat_, bVec_] :=
    Module[{x},
        x = Quiet @ Check[LeastSquares[N[aMat], N[bVec]], $Failed];
        If[VectorQ[x, NumericQ], Developer`ToPackedArray[x], $Failed]
    ];

End[];

EndPackage[];
