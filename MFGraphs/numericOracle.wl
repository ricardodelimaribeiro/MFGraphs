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

End[];
EndPackage[];
