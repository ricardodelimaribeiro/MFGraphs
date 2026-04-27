#!/usr/bin/env wolframscript
(* Benchmark: accumulateLinearRules vs TripleClean preprocessing
   Measures time, rules extracted, variables remaining, and disjunct count
   for each approach on the same raw constraint systems. *)

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

(* Load archive helpers: CachedSolve, SystemToTriple, TripleClean *)
Get[FileNameJoin[{$RepoRoot, "MFGraphs", "archive", "DNFReduce.wl"}]];
Get[FileNameJoin[{$RepoRoot, "MFGraphs", "archive", "DataToEquations.wl"}]];

SetAttributes[timed, HoldFirst];
timed[expr_] := Module[{t, r}, {t, r} = AbsoluteTiming[expr]; <|"Time" -> t, "Result" -> r|>];

(* Count Or-branches across the entire expression tree *)
countDisjuncts[expr_] :=
    Total @ Cases[expr, _Or :> Length[Unevaluated @@ {expr}], {0, Infinity}] +
    Total @ Cases[expr, or_Or :> Length[or], {0, Infinity}];

(* Better: count all Or subexpressions and sum their branch counts *)
countDisjuncts[expr_] :=
    Total[Length[#] & /@ Cases[expr, _Or, {0, Infinity}]];

(* Extract raw constraints from an mfgSystem, same logic as reduceSystem *)
extractConstraints[sys_] :=
    Module[{eqEntryIn, eqSplit, eqGather, eqHJ, ineqJs, ineqJts,
            ineqSBV, altFlows, altTrans, altOpt, ruleExit, base, allVars},
        eqEntryIn = systemData[sys, "EqEntryIn"];
        eqSplit    = systemData[sys, "EqBalanceSplittingFlows"];
        eqGather   = systemData[sys, "EqBalanceGatheringFlows"];
        eqHJ       = systemData[sys, "EqGeneral"];
        ineqJs     = systemData[sys, "IneqJs"];
        ineqJts    = systemData[sys, "IneqJts"];
        ineqSBV    = systemData[sys, "IneqSwitchingByVertex"];
        altFlows   = systemData[sys, "AltFlows"];
        altTrans   = systemData[sys, "AltTransitionFlows"];
        altOpt     = systemData[sys, "AltOptCond"];
        ruleExit   = Normal @ systemData[sys, "RuleExitValues"];
        base = And[And @@ eqEntryIn, eqSplit, eqGather, eqHJ,
                   ineqJs, ineqJts, ineqSBV, altFlows, altTrans, altOpt];
        allVars = Select[
            Variables[(base /. ruleExit) /. {Equal -> List, Or -> List, And -> List}],
            MatchQ[#, j[__] | u[__]] &
        ];
        <|"Base" -> base, "ExitRules" -> ruleExit, "AllVars" -> allVars,
          "IneqsBefore"     -> countIneqs[base /. ruleExit],
          "DisjunctsBefore" -> countDisjuncts[base /. ruleExit],
          "OrSummaryBefore" -> orBranchSummary[base /. ruleExit]|>
    ];

residualVarCount[constraints_, vars_List] :=
    Length @ Select[
        Variables[constraints /. {Equal -> List, Or -> List, And -> List}],
        MemberQ[vars, #] &
    ];

countIneqs[constraints_] :=
    Length @ Cases[
        Switch[Head[constraints], And, List @@ constraints, _, {constraints}],
        _LessEqual | _GreaterEqual | _Less | _Greater | _Inequality
    ];

(* Or-branch analysis: count top-level Or conjuncts and their branch sizes *)
orBranchSummary[constraints_] :=
    Module[{orExprs, branchCounts},
        orExprs     = Cases[Switch[Head[constraints], And, List @@ constraints, _, {constraints}], _Or];
        branchCounts = Length /@ orExprs;
        <|
            "OrConjunctCount" -> Length[orExprs],
            "BranchCountTally" -> Sort @ Tally[branchCounts],
            "DNFExpansionSize" -> If[orExprs === {}, 1, Times @@ branchCounts]
        |>
    ];

runAccumulateLinear[data_Association] :=
    Module[{m, constraints, rulesAcc, varsAfter},
        m = timed[
            solversTools`Private`accumulateLinearRules[
                data["Base"], data["AllVars"], data["ExitRules"]]
        ];
        {constraints, rulesAcc} = m["Result"];
        varsAfter = residualVarCount[constraints, data["AllVars"]];
        <|
            "Method"          -> "accumulateLinearRules",
            "TimeSeconds"     -> m["Time"],
            "RulesExtracted"  -> Length[rulesAcc] - Length[data["ExitRules"]],
            "VarsBefore"      -> Length[data["AllVars"]],
            "VarsAfter"       -> varsAfter,
            "IneqsAfter"      -> countIneqs[constraints],
            "DisjunctsAfter"  -> countDisjuncts[constraints],
            "OrSummaryAfter"  -> orBranchSummary[constraints]
        |>
    ];

(* TripleClean adapter: decompose flat And into {EE, NN, OR} then run TripleClean *)
runTripleClean[data_Association] :=
    Module[{after, triple, initRules, m, eeFinal, nnFinal, orFinal, rulesFinal,
            residual, varsAfter},
        after = data["Base"] /. data["ExitRules"];
        triple = SystemToTriple[after];
        initRules = Association @ data["ExitRules"];
        m = timed[TripleClean[{triple, initRules}]];
        {{eeFinal, nnFinal, orFinal}, rulesFinal} = m["Result"];
        residual = And[eeFinal, nnFinal, orFinal];
        varsAfter = residualVarCount[residual, data["AllVars"]];
        <|
            "Method"          -> "TripleClean",
            "TimeSeconds"     -> m["Time"],
            "RulesExtracted"  -> Length[rulesFinal] - Length[data["ExitRules"]],
            "VarsBefore"      -> Length[data["AllVars"]],
            "VarsAfter"       -> varsAfter,
            "IneqsAfter"      -> countIneqs[residual],
            "DisjunctsAfter"  -> countDisjuncts[residual],
            "OrSummaryAfter"  -> orBranchSummary[residual]
        |>
    ];

(* Scenarios *)
scenarioSpecs = {
    <|"Name" -> "chain_3v",
      "Build" -> Function[{}, gridScenario[{3}, {{1, 120}}, {{2, 0}, {3, 10}}]]|>,
    <|"Name" -> "y_case_8",
      "Build" -> Function[{}, getExampleScenario[8, {{1, 80}}, {{3, 0}, {4, 10}}]]|>,
    <|"Name" -> "attraction_11",
      "Build" -> Function[{}, getExampleScenario[11, {{1, 100}}, {{4, 0}}]]|>,
    <|"Name" -> "grid_2x2",
      "Build" -> Function[{}, gridScenario[{2, 2}, {{1, 100}}, {{4, 0}}]]|>
};

runCase[spec_Association] :=
    Module[{name, s, sys, data, resAL, resTC},
        name = spec["Name"];
        Print["--- ", name, " ---"];
        s   = spec["Build"][];
        sys = makeSystem[s];
        data = extractConstraints[sys];
        resAL = runAccumulateLinear[data];
        resTC = runTripleClean[data];
        <|"Case" -> name, "Data" -> data, "accumulateLinearRules" -> resAL, "TripleClean" -> resTC|>
    ];

results = runCase /@ scenarioSpecs;

col[s_, w_] := StringPadRight[ToString[s], w];

Print[];
Print["=== Summary ==="];
Print[col["Case", 15], col["Method", 24],
      col["Time(s)", 9], col["Rules+", 7],
      col["Vars", 14], col["Ineqs", 14], col["Or-conjuncts", 14], "DNF expansion"];
Print[col["", 15], col["", 24],
      col["", 9], col["", 7],
      col["bef->aft", 14], col["bef->aft", 14], col["bef->aft", 14], "bef->aft"];
Print[StringRepeat["-", 110]];
Do[
    With[{al = r["accumulateLinearRules"], tc = r["TripleClean"], d = r["Data"]},
        Print[col[r["Case"], 15], col["accumulateLinearRules", 24],
              col[NumberForm[al["TimeSeconds"], {4,3}], 9],
              col[al["RulesExtracted"], 7],
              col[ToString[al["VarsBefore"]] <> "->" <> ToString[al["VarsAfter"]], 14],
              col[ToString[d["IneqsBefore"]] <> "->" <> ToString[al["IneqsAfter"]], 14],
              col[ToString[d["OrSummaryBefore"]["OrConjunctCount"]] <> "->" <> ToString[al["OrSummaryAfter"]["OrConjunctCount"]], 14],
              ToString[d["OrSummaryBefore"]["DNFExpansionSize"]] <> "->" <> ToString[al["OrSummaryAfter"]["DNFExpansionSize"]]];
        Print[col["", 15], col["TripleClean", 24],
              col[NumberForm[tc["TimeSeconds"], {4,3}], 9],
              col[tc["RulesExtracted"], 7],
              col[ToString[tc["VarsBefore"]] <> "->" <> ToString[tc["VarsAfter"]], 14],
              col[ToString[d["IneqsBefore"]] <> "->" <> ToString[tc["IneqsAfter"]], 14],
              col[ToString[d["OrSummaryBefore"]["OrConjunctCount"]] <> "->" <> ToString[tc["OrSummaryAfter"]["OrConjunctCount"]], 14],
              ToString[d["OrSummaryBefore"]["DNFExpansionSize"]] <> "->" <> ToString[tc["OrSummaryAfter"]["DNFExpansionSize"]]];
        Print[]
    ],
    {r, results}
];

exportResults = Map[KeyDrop[#, "Data"] &, results];
Export[FileNameJoin[{$RepoRoot, "Scripts", "benchmark_preprocessing_results.json"}],
       exportResults, "JSON"];
Print["Results saved to Scripts/benchmark_preprocessing_results.json"];
