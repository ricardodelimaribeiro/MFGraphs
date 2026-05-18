#!/usr/bin/env wolframscript
(* PR 3 / Phase A.3 of the disjunct-influence study: for each scenario,
   per-disjunct sensitivity. Bound to scenarios where:
     - total disjunct count <= $maxDisjuncts (default 30)
     - baseline solve completes within $baselineTimeout (default 10s)

   For each Or-disjunct in the system's Alt* buckets:
     - removeResult: solve the system with the disjunct replaced by True
     - forceResults: per branch of the Or, solve the system with the
                     branch's equation forced in place of the Or
   Compare against the baseline (L-infinity over flow-variable RHS values
   present in both rule sets). Emit JSON one file per scenario.

   This is the most expensive Research script; cap total runtime via the
   per-call $perCallTimeout. Skipped scenarios produce a marker JSON with
   reason = "skipped:<why>" so the workbook can still enumerate them. *)

$RepoRoot = ParentDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
PrependTo[$Path, FileNameJoin[{$RepoRoot, "MFGraphs"}]];
Needs["MFGraphs`"];
Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];

$resultsDir = FileNameJoin[{$RepoRoot, "Scripts", "Research", "results", "sensitivity"}];
If[!DirectoryQ[$resultsDir],
    CreateDirectory[$resultsDir, CreateIntermediateDirectories -> True]];

$maxDisjuncts    = 30;
$baselineTimeout = 10;
$perCallTimeout  = 30;

sanitizeKey[k_Integer] := "case_" <> ToString[k];
sanitizeKey[k_String]  := StringReplace[k, " " -> "_"];

(* Replace target Or in the four Alt* buckets with `replacement`. *)
substituteAlt[sys_, target_, replacement_] :=
    Module[{assoc = systemData[sys], buckets = {"AltFlows", "AltTransitionFlows",
                                                 "AltOptCond", "AltExitCond"}, swap},
        swap = Function[bucket,
            Module[{lst = Lookup[assoc, bucket, {}]},
                If[ListQ[lst],
                    bucket -> Replace[lst, target -> replacement, {1}],
                    bucket -> lst]
            ]];
        mfgSystem[Join[assoc, Association[swap /@ buckets]]]
    ];

(* activeSetReduceSystem returns a flat List of Rule[var, value] (not an
   Association). Treat both shapes (list-of-rules, or association with
   "Rules" key) defensively for forward-compat. *)
ruleAssoc[res_] :=
    Which[
        ListQ[res] && AllTrue[res, MatchQ[#, _Rule | _RuleDelayed] &],
            Association[res],
        AssociationQ[res] && KeyExistsQ[res, "Rules"],
            Association @ res["Rules"],
        AssociationQ[res],
            res,
        True,
            <||>
    ];

(* Result is "good" when it is either a list of rules or an Association
   with rules. Either shape is acceptable for the sensitivity comparison. *)
goodResultQ[res_] :=
    (ListQ[res] && AllTrue[res, MatchQ[#, _Rule | _RuleDelayed] &]) ||
    AssociationQ[res];

solveBounded[sys_] :=
    Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
        CheckAbort[
            TimeConstrained[activeSetReduceSystem[sys], $perCallTimeout, $TimedOut],
            $Aborted]];

(* L-infinity over keys common to both rule sets, after numerical evaluation.
   Returns a real or one of: "Infeasible", "NonComparable", "Timeout". *)
deltaVsBaseline[modifiedRes_, baseAssoc_Association] :=
    Module[{modAssoc, common, diffs},
        Which[
            modifiedRes === $TimedOut,                 Return["Timeout"],
            modifiedRes === $Aborted,                  Return["NonComparable"],
            !goodResultQ[modifiedRes],                 Return["NonComparable"],
            AssociationQ[modifiedRes] &&
              Lookup[modifiedRes, "Residual", True] === False, Return["Infeasible"]
        ];
        modAssoc = ruleAssoc[modifiedRes];
        common = Intersection[Keys[baseAssoc], Keys[modAssoc]];
        If[common === {}, Return["NonComparable"]];
        diffs = Function[k,
            Module[{a = baseAssoc[k], b = modAssoc[k]},
                If[NumericQ[a] && NumericQ[b],
                    Abs[N[a - b]],
                    Missing[]]]] /@ common;
        diffs = Cases[diffs, _?NumericQ];
        If[diffs === {}, "NonComparable", N[Max[diffs]]]
    ];

measureScenario[label_, scenario_] :=
    Module[{sys, ors, baseRes, baseAssoc, baseTime, t0, rows, json, outFile},
        sys = makeSystem[scenario];
        ors = Cases[
            Flatten[{
                systemData[sys, "AltFlows"],
                systemData[sys, "AltTransitionFlows"],
                systemData[sys, "AltOptCond"],
                systemData[sys, "AltExitCond"]}],
            _Or];
        If[Length[ors] > $maxDisjuncts,
            Export[FileNameJoin[{$resultsDir, label <> ".json"}],
                <|"scenario" -> label, "skipped" -> True,
                  "reason" -> "disjunctCount>" <> ToString[$maxDisjuncts],
                  "disjunctCount" -> Length[ors]|>, "JSON"];
            Print["[skip] ", label, " disjuncts=", Length[ors]];
            Return[]];
        t0 = AbsoluteTime[];
        baseRes = Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
            CheckAbort[
                TimeConstrained[activeSetReduceSystem[sys], $baselineTimeout, $TimedOut],
                $Aborted]];
        baseTime = AbsoluteTime[] - t0;
        If[baseRes === $TimedOut || baseRes === $Aborted || !goodResultQ[baseRes],
            Export[FileNameJoin[{$resultsDir, label <> ".json"}],
                <|"scenario" -> label, "skipped" -> True,
                  "reason" -> "baselineTimeout",
                  "baselineSec" -> N[baseTime]|>, "JSON"];
            Print["[skip] ", label, " baseline=", baseTime, "s"];
            Return[]];
        baseAssoc = ruleAssoc[baseRes];
        rows = Function[orExpr,
            Module[{branches = List @@ orExpr, removeRes, removeDelta,
                    forceResults, dId},
                dId = ToString[Sort[branches], InputForm];
                removeRes = solveBounded[substituteAlt[sys, orExpr, True]];
                removeDelta = deltaVsBaseline[removeRes, baseAssoc];
                forceResults = Function[branch,
                    Module[{res = solveBounded[substituteAlt[sys, orExpr, branch]],
                            d},
                        d = deltaVsBaseline[res, baseAssoc];
                        <|"branch" -> ToString[branch, InputForm],
                          "delta"  -> d,
                          "tag"    -> Which[
                              d === "Infeasible", "Infeasible",
                              d === "Timeout",    "Timeout",
                              d === "NonComparable", "NonComparable",
                              NumericQ[d] && d < 10^-8, "Matches",
                              True, "Differs"]|>
                    ]] /@ branches;
                <|"disjunctId" -> dId,
                  "branchCount" -> Length[branches],
                  "removeDelta" -> removeDelta,
                  "removeTag"   -> Which[
                      removeDelta === "Infeasible",   "Infeasible",
                      removeDelta === "Timeout",      "Timeout",
                      removeDelta === "NonComparable","NonComparable",
                      NumericQ[removeDelta] && removeDelta < 10^-8, "Matches",
                      True, "Differs"],
                  "forceResults" -> forceResults|>
            ]] /@ ors;
        json = <|"scenario" -> label,
                 "skipped"  -> False,
                 "baselineSec" -> N[baseTime],
                 "disjunctCount" -> Length[ors],
                 "rows" -> rows|>;
        outFile = FileNameJoin[{$resultsDir, label <> ".json"}];
        Export[outFile, json, "JSON"];
        Print["[done] ", label, " disjuncts=", Length[ors], " base=", N[baseTime], "s"]
    ];

Print["Sweep: bounds disjunctCount<=", $maxDisjuncts,
      " baseline<=", $baselineTimeout, "s"];
Scan[
    Function[k,
        Module[{ee = $ExampleBoundaryDefaults[k]},
            Quiet @ Check[
                measureScenario[sanitizeKey[k], getExampleScenario[k, ee[[1]], ee[[2]]]],
                Print["[FAIL] ", k]]]],
    Keys[$ExampleBoundaryDefaults]
];
Print["Done."];
