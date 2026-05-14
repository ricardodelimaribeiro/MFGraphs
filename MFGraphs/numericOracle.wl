(* ::Package:: *)

(* numericOracle.wl
   ----------------
   OPT-IN subpackage. NOT loaded by Needs["MFGraphs`"].
   Load explicitly: Needs["numericOracle`"] (after Needs["MFGraphs`"]).

   This is the ONLY file in MFGraphs that introduces floating-point work.
   It exposes a single function `numericOracleClassify[sys]` that takes an
   exact mfgSystem, runs a single NMinimize over the complementarity
   objective, and returns an Association classifying each flow variable
   as "Active", "Inactive", or "Ambiguous" based on a two-threshold safety
   band. The result contains only symbolic variable names plus a few
   diagnostic scalars; no numerical rules cross the package boundary.

   Pair with `addOracleEqualities` (in systemTools, symbolic-only) which
   consumes the classification and produces a new mfgSystem with j[e] == 0
   appended to EqGeneral for each Inactive variable. The exact symbolic
   solver then runs on the pruned system and produces a rational answer.

   Architecture mirrors the archived FictitiousPlayBackend.wl + the Phase 5/6
   oracle bridge, but with the iteration loop replaced by a single global
   NMinimize call. The known vulnerabilities of the archived iteration
   (cost-scale fragility, oscillation, DAG cycles) are bypassed because
   there is no iteration. The float-to-exact chasm (vulnerability #3) is
   addressed by the two-threshold safety band.
*)

BeginPackage["numericOracle`", {"primitives`", "utilities`", "scenarioTools`", "unknownsTools`", "systemTools`", "solversTools`"}];

numericOracleClassify::usage = "numericOracleClassify[sys] runs FindInstance over the LP relaxation (drops complementarity Or-disjunctions), classifies each flow variable as Active / Inactive / Ambiguous from one feasible vertex, and returns an Association: <|\"Inactive\" -> {...}, \"Active\" -> {...}, \"Ambiguous\" -> {...}, \"Converged\" -> True|False|>. Options: \"TightThreshold\" (default 10^-8), \"SafeThreshold\" (default 10^-3), \"Timeout\" (default 30s for FindInstance). All floats live inside this subpackage; the returned Association contains only symbolic variable names.

Caveat: a one-shot LP vertex's zero coordinates aren't necessarily zero at equilibrium. Pruning by these classifications can produce an infeasible system (Residual -> False) on scenarios where the LP vertex doesn't coincide with the equilibrium support. Use solveScenarioWithOracle for an automatic fallback.";

numericFictitiousPlayClassify::usage = "numericFictitiousPlayClassify[sys, opts] runs an iterative LP-refinement classifier inspired by the archived FictitiousPlayBackend. Each iteration solves the LP relaxation, classifies edges by current flow, and tightens constraints by enforcing complementarity at confidently-classified edges; the loop terminates when the support classification is stable for StableIterationsLimit consecutive iterations or MaxIterations is exhausted. Returns the same Association shape as numericOracleClassify, so the same addOracleEqualities bridge can consume it.

Options: MaxIterations (30), TightThreshold (10^-8), SafeThreshold (10^-3), StableIterationsLimit (3), Timeout (30s per LP), TieBreakThreshold (10^-5).

Designed to be more reliable than the one-shot numericOracleClassify because it converges on a single self-consistent classification rather than reading off one arbitrary LP vertex. Same caveat applies: the safety band keeps an edge as Ambiguous when its flow falls between TightThreshold and SafeThreshold, so addOracleEqualities never pins a 'maybe-zero' variable.";

solveScenarioWithOracle::usage = "solveScenarioWithOracle[s] composes the full opt-in oracle pipeline on scenario s: makeSymbolicUnknowns -> makeSystem -> addSymmetryEqualities -> numericOracleClassify -> addOracleEqualities -> activeSetReduceSystem. If the pruned solve returns Residual -> False (the oracle over-pruned past the safety probe in addOracleEqualities), the wrapper transparently falls back to activeSetReduceSystem on the symmetry-augmented but unpruned system and returns that result. solveScenarioWithOracle[s, solver] uses the supplied solver (must accept an mfgSystem) in place of activeSetReduceSystem. Options TightThreshold / SafeThreshold / Timeout are forwarded to numericOracleClassify.";

(* The wrapper below codifies this pipeline; the comment is retained as a
   reference for callers who want to invoke the stages by hand:
     sys     = makeSystem[scenario];
     sym     = addSymmetryEqualities[sys, scenario];
     oracle  = numericOracleClassify[sym];
     pruned  = addOracleEqualities[sym, oracle];
     result  = activeSetReduceSystem[pruned];
     (* If Lookup[result, "Residual", True] === False the oracle over-pruned;
        rerun activeSetReduceSystem[sym] to fall back. *)
*)

Begin["`Private`"];

(* Build the LP-relaxation constraints (drop complementarity Or-disjunctions).
   Returns the symbolic constraints + the variable list. NO floats yet. *)
buildLPRelaxation[sys_] :=
    Module[{deterministic, allVars},
        deterministic = And[
            And @@ systemData[sys, "EqEntryIn"],
            systemData[sys, "EqBalanceSplittingFlows"],
            systemData[sys, "EqBalanceGatheringFlows"],
            And @@ systemData[sys, "EqGeneral"],
            And @@ systemData[sys, "IneqJs"],
            And @@ systemData[sys, "IneqJts"],
            And @@ systemData[sys, "IneqSwitchingByVertex"],
            And @@ systemData[sys, "IneqExitValues"]
        ];
        allVars = Join[
            systemData[sys, "Js"],
            systemData[sys, "Jts"],
            systemData[sys, "Us"]
        ];
        <|"Constraints" -> deterministic, "Vars" -> allVars|>
    ];

(* ---- Matrix-form LP path (opt-in via "UseMatrixSolver" -> True) ----
   Symbolic FindInstance / Minimize on the LP relaxation does not scale
   past ~300 vars (FindInstance times out at 90s on Grid1010's 1696 vars).
   These helpers convert the same constraints to numeric matrix form via
   CoefficientArrays and call LinearProgramming directly, bypassing the
   symbolic CAD overhead. Same float-isolation invariant: returned values
   are still flow-variable point assignments, classified into Inactive /
   Active / Ambiguous lists by the same downstream code. *)

$matrixLPTimeout = 30;

buildLPMatrices[sys_] :=
    Module[{vars, jVars, eqList, ineqList, eqExprs, ineqExprs,
            eqArrays, ineqArrays, AEq, bEq, AIneq, bIneq},
        vars = Join[
            systemData[sys, "Js"],
            systemData[sys, "Jts"],
            systemData[sys, "Us"]
        ];
        jVars = Cases[vars, _j];
        (* Equalities: collect from EqEntryIn, EqBalance*, EqGeneral.
           EqBalance* are stored as a single And expression; flatten. *)
        eqList = DeleteCases[
            Flatten[{
                systemData[sys, "EqEntryIn"],
                If[Head[systemData[sys, "EqBalanceSplittingFlows"]] === And,
                    List @@ systemData[sys, "EqBalanceSplittingFlows"],
                    {systemData[sys, "EqBalanceSplittingFlows"]}],
                If[Head[systemData[sys, "EqBalanceGatheringFlows"]] === And,
                    List @@ systemData[sys, "EqBalanceGatheringFlows"],
                    {systemData[sys, "EqBalanceGatheringFlows"]}],
                systemData[sys, "EqGeneral"]
            }],
            True
        ];
        ineqList = DeleteCases[
            Flatten[{
                systemData[sys, "IneqJs"],
                systemData[sys, "IneqJts"],
                systemData[sys, "IneqSwitchingByVertex"],
                systemData[sys, "IneqExitValues"]
            }],
            True
        ];
        (* Convert l == r to (l - r), then CoefficientArrays. *)
        eqExprs   = eqList   /. (l_ == r_) :> l - r;
        (* l >= r becomes l - r (nonneg quantity). *)
        ineqExprs = ineqList /. {(l_ >= r_) :> l - r, (l_ <= r_) :> r - l};
        eqArrays   = If[eqExprs   === {}, Null, CoefficientArrays[eqExprs, vars]];
        ineqArrays = If[ineqExprs === {}, Null, CoefficientArrays[ineqExprs, vars]];
        (* For sparse equality, A_eq.x + c_eq = 0 ⇒ A_eq.x = -c_eq. *)
        bEq   = If[eqArrays   === Null, {}, -Normal[eqArrays[[1]]]];
        AEq   = If[eqArrays   === Null, {}, Normal[eqArrays[[2]]]];
        bIneq = If[ineqArrays === Null, {}, -Normal[ineqArrays[[1]]]];
        AIneq = If[ineqArrays === Null, {}, Normal[ineqArrays[[2]]]];
        (* Sanity: all eqArrays must be degree-1 (linear). *)
        If[(eqArrays =!= Null && Length[eqArrays] != 2) ||
           (ineqArrays =!= Null && Length[ineqArrays] != 2),
           Return[$Failed, Module]];
        <|
            "Vars"     -> vars,
            "JVars"    -> jVars,
            "JIndices" -> Flatten[Position[vars, _j, {1}, Heads -> False]],
            "AEq"      -> AEq,
            "bEq"      -> bEq,
            "AIneq"    -> AIneq,
            "bIneq"    -> bIneq
        |>
    ];

(* Solve LinearProgramming with given objective vector c (over numeric vars).
   Returns Association[var -> value] on success, $Aborted on timeout, $Failed
   on infeasibility. *)
solveLPMatrix[matrices_, c_] :=
    Module[{n, A, bSign, bounds, jIdx, sol},
        n     = Length[matrices["Vars"]];
        jIdx  = matrices["JIndices"];
        A     = If[matrices["AEq"]   === {}, {}, matrices["AEq"]];
        A     = If[matrices["AIneq"] === {}, A, Join[A, matrices["AIneq"]]];
        bSign = Join[
            Transpose[{matrices["bEq"],   ConstantArray[0, Length[matrices["bEq"]]]}],
            Transpose[{matrices["bIneq"], ConstantArray[1, Length[matrices["bIneq"]]]}]
        ];
        bounds = Table[{-Infinity, Infinity}, {n}];
        Do[bounds[[i]] = {0., Infinity}, {i, jIdx}];
        sol = TimeConstrained[
            Quiet @ Check[
                LinearProgramming[N[c], N[A], N[bSign], N[bounds]],
                $Failed],
            $matrixLPTimeout,
            $Aborted
        ];
        Which[
            sol === $Aborted,                $Aborted,
            sol === $Failed,                 $Failed,
            ! VectorQ[sol, NumericQ],        $Failed,
            True,                            AssociationThread[matrices["Vars"], sol]
        ]
    ];

(* Random-objective LP solve for the K-fold consensus path.
   sign = +1 ⇒ minimize ∑ wᵢ jᵢ (small-flow vertex)
   sign = -1 ⇒ maximize ∑ wᵢ jᵢ (large-flow vertex; LinearProgramming
   minimizes c.x, so we negate c). *)
solveLPMatrixRandom[matrices_, sign_:1] :=
    Module[{n, c, jIdx, coeffs},
        n      = Length[matrices["Vars"]];
        jIdx   = matrices["JIndices"];
        c      = ConstantArray[0., n];
        coeffs = sign * RandomReal[{0.1, 1.0}, Length[jIdx]];
        Do[c[[ jIdx[[k]] ]] = coeffs[[k]], {k, Length[jIdx]}];
        solveLPMatrix[matrices, c]
    ];

Options[numericOracleClassify] = {
    "TightThreshold"  -> 10^-8,
    "SafeThreshold"   -> 10^-3,
    "Timeout"         -> 30,
    "ConfirmTimeout"  -> 0.5,
    "UseMatrixSolver" -> False,
    "MatrixLPTimeout" -> 30
};

numericOracleClassify[sys_?mfgSystemQ, OptionsPattern[]] :=
    Module[{lp, constraints, vars, jVars, fi, point, classify,
            inactive, active, ambiguous, eps0, eps1, timeout, t0, elapsed,
            useMatrix, matrices},
        eps0      = OptionValue["TightThreshold"];
        eps1      = OptionValue["SafeThreshold"];
        timeout   = OptionValue["Timeout"];
        useMatrix = TrueQ[OptionValue["UseMatrixSolver"]];
        lp        = buildLPRelaxation[sys];
        constraints = lp["Constraints"];
        vars        = lp["Vars"];
        jVars       = Cases[vars, _j];
        t0 = AbsoluteTime[];
        If[useMatrix,
            (* Matrix path: build matrices once, solve with LinearProgramming. *)
            matrices = buildLPMatrices[sys];
            If[matrices === $Failed,
                Return[<|"Inactive" -> {}, "Active" -> {}, "Ambiguous" -> {},
                         "Converged" -> False, "Reason" -> "MatrixBuildFailed",
                         "Elapsed" -> AbsoluteTime[] - t0|>, Module]
            ];
            Block[{$matrixLPTimeout = OptionValue["MatrixLPTimeout"]},
                fi = solveLPMatrix[matrices,
                    (* zero objective ≡ FindInstance: any feasible point *)
                    ConstantArray[0., Length[matrices["Vars"]]]]
            ];
            elapsed = AbsoluteTime[] - t0;
            If[fi === $Aborted || fi === $Failed,
                Return[<|"Inactive" -> {}, "Active" -> {}, "Ambiguous" -> {},
                         "Converged" -> False,
                         "Reason" -> If[fi === $Aborted, "MatrixLPTimeout", "MatrixLPInfeasible"],
                         "Elapsed" -> elapsed|>, Module]
            ];
            point = Normal[fi],
            (* Symbolic path (default): FindInstance on the And expression. *)
            fi = TimeConstrained[
                Quiet @ FindInstance[constraints, vars, Reals, 1],
                timeout,
                $Aborted
            ];
            elapsed = AbsoluteTime[] - t0;
            If[fi === $Aborted || ! MatchQ[fi, {{___Rule}, ___}] || fi === {},
                Return[<|"Inactive" -> {}, "Active" -> {}, "Ambiguous" -> {},
                         "Converged" -> False,
                         "Reason" -> If[fi === $Aborted, "FindInstanceTimeout", "FindInstanceEmpty"],
                         "Elapsed" -> elapsed|>, Module]
            ];
            point = First[fi]
        ];
        classify[v_] :=
            With[{x = Abs[v /. point]},
                Which[
                    ! NumericQ[x],   "Ambiguous",
                    x < eps0,        "Inactive",
                    x > eps1,        "Active",
                    True,            "Ambiguous"
                ]
            ];
        inactive  = Select[jVars, classify[#] === "Inactive" &];
        active    = Select[jVars, classify[#] === "Active"   &];
        ambiguous = Select[jVars, classify[#] === "Ambiguous" &];
        (* WARNING: a one-shot LP FindInstance returns ONE feasible vertex.
           Variables that are zero at that vertex are not necessarily zero
           at the equilibrium. The downstream addOracleEqualities runs a
           safety check to detect over-pruning; if the pruned system turns
           out to be infeasible the caller sees Residual -> False. *)
        <|"Inactive"  -> inactive,
          "Active"    -> active,
          "Ambiguous" -> ambiguous,
          "Converged" -> True,
          "Method"    -> "LPRelaxation",
          "Elapsed"   -> elapsed|>
    ];

(* ---- K-fold consensus oracle ("FictitiousPlay" by analogy) ----
   Solves the LP relaxation K times with K different random linear
   objectives, then classifies a flow variable as Inactive only if it
   was zero in EVERY LP solve. This is the consensus / voting analogue
   of Fictitious Play for our setting: rather than iterating a softmax
   best-response, we sample K different polytope vertices and trust the
   intersection of their zero-sets.

   Why this instead of the archived softmax FP: the archive's softmax
   needs Bellman-potential extraction at every iteration, which entangles
   with helpers (RecoverCriticalFlowAssociation, BuildMonotoneValueSystem)
   that don't exist in the current package. K-fold consensus achieves
   the same goal -- a more reliable Inactive classification than one-shot
   LP -- using only the LP solver primitives we already have.

   Mitigations from the Phase 5/6 technical-debt memory:
   - Cost-scale fragility: N/A (we solve LP, not softmax).
   - Symmetric oscillation: handled by the safety band -- tied edges
     stay Ambiguous, and consensus across K random objectives breaks
     symmetric ties statistically.
   - Float-to-exact chasm: same two-threshold safety band as
     numericOracleClassify. *)

Options[numericFictitiousPlayClassify] = {
    "MaxIterations"          -> 5,        (* K, the number of LP solves *)
    "TightThreshold"         -> 10^-8,
    "SafeThreshold"          -> 10^-3,
    "Timeout"                -> 30,       (* per-LP timeout *)
    "ObjectiveSeed"          -> 42,       (* deterministic random seed *)
    "UseMatrixSolver"        -> False,
    "MatrixLPTimeout"        -> 30
};

numericFictitiousPlayClassify[sys_?mfgSystemQ, OptionsPattern[]] :=
    Module[{vars, jVars, baseConstraints, eps0, eps1, K, perLPTimeout, seed,
            points, classifications, allInactive, alwaysInactive, alwaysActive,
            ambiguous, totalElapsed, t0, classify, kSucceeded, useMatrix, matrices},
        eps0          = OptionValue["TightThreshold"];
        eps1          = OptionValue["SafeThreshold"];
        K             = OptionValue["MaxIterations"];
        perLPTimeout  = OptionValue["Timeout"];
        seed          = OptionValue["ObjectiveSeed"];
        useMatrix     = TrueQ[OptionValue["UseMatrixSolver"]];
        vars   = Join[
            systemData[sys, "Js"], systemData[sys, "Jts"], systemData[sys, "Us"]];
        jVars  = Cases[vars, _j];
        matrices = If[useMatrix, buildLPMatrices[sys], Null];
        If[useMatrix && matrices === $Failed,
            Return[<|"Inactive" -> {}, "Active" -> {}, "Ambiguous" -> jVars,
                     "Converged" -> False, "Method" -> "FictitiousPlayConsensus",
                     "Iterations" -> 0, "StopReason" -> "MatrixBuildFailed",
                     "Elapsed" -> 0|>, Module]
        ];
        baseConstraints = And[
            And @@ systemData[sys, "EqEntryIn"],
            systemData[sys, "EqBalanceSplittingFlows"],
            systemData[sys, "EqBalanceGatheringFlows"],
            And @@ systemData[sys, "EqGeneral"],
            And @@ systemData[sys, "IneqJs"],
            And @@ systemData[sys, "IneqJts"],
            And @@ systemData[sys, "IneqSwitchingByVertex"],
            And @@ systemData[sys, "IneqExitValues"]
        ];
        classify[v_, pt_] :=
            With[{x = Abs[v /. pt]},
                Which[
                    ! NumericQ[x],   "Ambiguous",
                    x < eps0,        "Inactive",
                    x > eps1,        "Active",
                    True,            "Ambiguous"
                ]
            ];
        SeedRandom[seed];
        t0     = AbsoluteTime[];
        points = {};
        Do[
            Module[{coeffs, obj, res},
                (* Diverse vertex sampling: k=1 plain FindInstance (one
                   feasible vertex), k=2..K alternate Minimize/Maximize on
                   random nonnegative-coefficient objectives to push toward
                   different polytope corners without going unbounded.
                   Matrix path uses LinearProgramming directly. *)
                If[useMatrix,
                    Block[{$matrixLPTimeout = OptionValue["MatrixLPTimeout"]},
                        (* Always positive random coefficients: minimizing a
                           positive-weighted j-sum is bounded (j-vars >= 0),
                           and different random weights pick different vertices
                           of the polytope. Negative coefficients would make
                           LinearProgramming unbounded because j-vars have no
                           upper bound. *)
                        res = solveLPMatrixRandom[matrices, 1]
                    ];
                    If[AssociationQ[res], AppendTo[points, Normal[res]]],
                    Which[
                        k == 1,
                            res = TimeConstrained[
                                Quiet @ FindInstance[baseConstraints, vars, Reals, 1],
                                perLPTimeout, $Aborted];
                            If[MatchQ[res, {{___Rule}, ___}] && res =!= {},
                                AppendTo[points, First[res]]
                            ],
                        EvenQ[k],
                            coeffs = RandomReal[{0.1, 1.0}, Length[jVars]];
                            obj    = coeffs . jVars;
                            res = TimeConstrained[
                                Quiet @ Minimize[{obj, baseConstraints}, vars],
                                perLPTimeout, $Aborted];
                            If[MatchQ[res, {_?NumericQ, {___Rule}}] && NumericQ[First[res]] && Abs[First[res]] < Infinity,
                                AppendTo[points, Last[res]]
                            ],
                        True,
                            coeffs = RandomReal[{0.1, 1.0}, Length[jVars]];
                            obj    = coeffs . jVars;
                            res = TimeConstrained[
                                Quiet @ Maximize[{obj, baseConstraints}, vars],
                                perLPTimeout, $Aborted];
                            If[MatchQ[res, {_?NumericQ, {___Rule}}] && NumericQ[First[res]] && Abs[First[res]] < Infinity,
                                AppendTo[points, Last[res]]
                            ]
                    ]
                ]
            ],
            {k, K}
        ];
        totalElapsed = AbsoluteTime[] - t0;
        kSucceeded   = Length[points];
        If[kSucceeded == 0,
            Return[<|"Inactive" -> {}, "Active" -> {}, "Ambiguous" -> jVars,
                     "Converged" -> False, "Method" -> "FictitiousPlayConsensus",
                     "Iterations" -> 0, "StopReason" -> "AllLPsFailed",
                     "Elapsed" -> totalElapsed|>, Module]
        ];
        (* Per-LP classifications. *)
        classifications = Function[pt,
            AssociationThread[jVars, classify[#, pt] & /@ jVars]] /@ points;
        (* Consensus rules: Inactive iff Inactive in ALL points; Active iff
           Active in at least one point AND not Inactive in any point. *)
        alwaysInactive = Select[jVars,
            Function[v, AllTrue[classifications, #[v] === "Inactive" &]]];
        alwaysActive   = Select[jVars,
            Function[v,
                AnyTrue[classifications, #[v] === "Active" &] &&
                AllTrue[classifications, #[v] =!= "Inactive" &]
            ]];
        ambiguous      = Complement[jVars, alwaysInactive, alwaysActive];
        <|"Inactive"   -> alwaysInactive,
          "Active"     -> alwaysActive,
          "Ambiguous"  -> ambiguous,
          "Converged"  -> True,
          "Method"     -> "FictitiousPlayConsensus",
          "Iterations" -> kSucceeded,
          "StopReason" -> If[kSucceeded == K, "AllLPsCompleted", "PartialConsensus"],
          "Elapsed"    -> totalElapsed|>
    ];

(* ---- solveScenarioWithOracle: high-level wrapper ----
   Runs the full oracle pipeline and silently falls back to the unpruned
   active-set solve when the pruned solve ends with Residual -> False.
   addOracleEqualities (in systemTools) already runs a FindInstance probe
   to detect over-pruning before pruning; the fallback here handles the
   secondary case where the active-set solver itself bottoms out on a
   pruned-but-feasible system.

   Forwarded options match numericOracleClassify so callers can tune the
   classifier from one call site. *)
Options[solveScenarioWithOracle] = {
    "TightThreshold" -> 10^-8,
    "SafeThreshold"  -> 10^-3,
    "Timeout"        -> 30
};

solveScenarioWithOracle[s_?scenarioQ, opts:OptionsPattern[]] :=
    solveScenarioWithOracle[s, activeSetReduceSystem, opts];

solveScenarioWithOracle[s_?scenarioQ, solver_, OptionsPattern[]] :=
    Module[{unk, sys, sym, oracle, pruned, result},
        unk    = makeSymbolicUnknowns[s];
        sys    = makeSystem[s, unk];
        sym    = addSymmetryEqualities[sys, s];
        oracle = numericOracleClassify[sym,
            "TightThreshold" -> OptionValue["TightThreshold"],
            "SafeThreshold"  -> OptionValue["SafeThreshold"],
            "Timeout"        -> OptionValue["Timeout"]];
        pruned = addOracleEqualities[sym, oracle];
        result = solver[pruned];
        If[AssociationQ[result] && Lookup[result, "Residual", True] === False,
            solver[sym],
            result
        ]
    ];

End[];
EndPackage[];
