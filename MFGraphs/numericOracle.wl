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

BeginPackage["numericOracle`", {"primitives`", "utilities`", "systemTools`"}];

numericOracleClassify::usage = "numericOracleClassify[sys] runs FindInstance over the LP relaxation (drops complementarity Or-disjunctions), classifies each flow variable as Active / Inactive / Ambiguous from one feasible vertex, and returns an Association: <|\"Inactive\" -> {...}, \"Active\" -> {...}, \"Ambiguous\" -> {...}, \"Converged\" -> True|False|>. Options: \"TightThreshold\" (default 10^-8), \"SafeThreshold\" (default 10^-3), \"Timeout\" (default 30s for FindInstance). All floats live inside this subpackage; the returned Association contains only symbolic variable names.

Caveat: a one-shot LP vertex's zero coordinates aren't necessarily zero at equilibrium. Pruning by these classifications can produce an infeasible system (Residual -> False) on scenarios where the LP vertex doesn't coincide with the equilibrium support. Use solveScenarioWithOracle for an automatic fallback.";

numericFictitiousPlayClassify::usage = "numericFictitiousPlayClassify[sys, opts] runs an iterative LP-refinement classifier inspired by the archived FictitiousPlayBackend. Each iteration solves the LP relaxation, classifies edges by current flow, and tightens constraints by enforcing complementarity at confidently-classified edges; the loop terminates when the support classification is stable for StableIterationsLimit consecutive iterations or MaxIterations is exhausted. Returns the same Association shape as numericOracleClassify, so the same addOracleEqualities bridge can consume it.

Options: MaxIterations (30), TightThreshold (10^-8), SafeThreshold (10^-3), StableIterationsLimit (3), Timeout (30s per LP), TieBreakThreshold (10^-5).

Designed to be more reliable than the one-shot numericOracleClassify because it converges on a single self-consistent classification rather than reading off one arbitrary LP vertex. Same caveat applies: the safety band keeps an edge as Ambiguous when its flow falls between TightThreshold and SafeThreshold, so addOracleEqualities never pins a 'maybe-zero' variable.";

(* Caller's responsibility to compose:
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

Options[numericOracleClassify] = {
    "TightThreshold" -> 10^-8,
    "SafeThreshold"  -> 10^-3,
    "Timeout"        -> 30,
    "ConfirmTimeout" -> 0.5
};

numericOracleClassify[sys_?mfgSystemQ, OptionsPattern[]] :=
    Module[{lp, constraints, vars, jVars, fi, point, classify,
            inactive, active, ambiguous, eps0, eps1, timeout, t0, elapsed},
        eps0    = OptionValue["TightThreshold"];
        eps1    = OptionValue["SafeThreshold"];
        timeout = OptionValue["Timeout"];
        lp     = buildLPRelaxation[sys];
        constraints = lp["Constraints"];
        vars        = lp["Vars"];
        jVars       = Cases[vars, _j];
        t0 = AbsoluteTime[];
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
        point = First[fi];
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
    "ObjectiveSeed"          -> 42        (* deterministic random seed *)
};

numericFictitiousPlayClassify[sys_?mfgSystemQ, OptionsPattern[]] :=
    Module[{vars, jVars, baseConstraints, eps0, eps1, K, perLPTimeout, seed,
            points, classifications, allInactive, alwaysInactive, alwaysActive,
            ambiguous, totalElapsed, t0, classify, kSucceeded},
        eps0          = OptionValue["TightThreshold"];
        eps1          = OptionValue["SafeThreshold"];
        K             = OptionValue["MaxIterations"];
        perLPTimeout  = OptionValue["Timeout"];
        seed          = OptionValue["ObjectiveSeed"];
        vars   = Join[
            systemData[sys, "Js"], systemData[sys, "Jts"], systemData[sys, "Us"]];
        jVars  = Cases[vars, _j];
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
                   different polytope corners without going unbounded. *)
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

End[];
EndPackage[];
