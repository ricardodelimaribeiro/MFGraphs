#!/usr/bin/env wolframscript
(* Analyze residual solution branches by entry costs.

   Usage:
     wolframscript -file Scripts/AnalyzeSolutionBranches.wls chain-one-exit
     wolframscript -file Scripts/AnalyzeSolutionBranches.wls chain-two-exits
     wolframscript -file Scripts/AnalyzeSolutionBranches.wls example-7
     wolframscript -file Scripts/AnalyzeSolutionBranches.wls jamarat-captured

   The jamarat-captured mode reads the commented jamaratSol assignment from
   MFGraphs/Jamarat.wl instead of rerunning the expensive solver.
*)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

ClearAll[
    optionValue,
    makeCase,
    extractCapturedJamaratSolution,
    solutionRules,
    solutionResidual,
    dnfBranches,
    entryCostData,
    branchVariables,
    branchEntryCostSummary,
    analyzeSolutionBranches,
    printBranchSummary
];

optionValue[args_List, name_String, default_] :=
    Module[{prefix = "--" <> name <> "="},
        FirstCase[
            args,
            s_String /; StringStartsQ[s, prefix] :> StringDrop[s, StringLength[prefix]],
            default
        ]
    ];

makeCase["chain-one-exit"] :=
    Module[{s, sys, sol},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        <|"Scenario" -> s, "System" -> sys, "Solution" -> sol|>
    ];

makeCase["chain-two-exits"] :=
    Module[{s, sys, sol},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        <|"Scenario" -> s, "System" -> sys, "Solution" -> sol|>
    ];

makeCase["example-7"] :=
    Module[{s, sys, sol},
        s = gridScenario[{3}, {{2, 100}}, {{1, 0}, {3, 10}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        <|"Scenario" -> s, "System" -> sys, "Solution" -> sol|>
    ];

makeCase[name_String] /; StringMatchQ[name, "example-" ~~ DigitCharacter ..] :=
    Module[{key, s, sys, sol},
        key = ToExpression @ StringDrop[name, StringLength["example-"]];
        s = getExampleScenario[key, {{1, 100}}, {{3, 0}, {4, 10}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        <|"Scenario" -> s, "System" -> sys, "Solution" -> sol|>
    ];

extractCapturedJamaratSolution[] :=
    Module[{path, text, start, stop, snippet},
        path = FileNameJoin[{$RepoRoot, "MFGraphs", "Jamarat.wl"}];
        text = Import[path, "Text"];
        start = StringPosition[text, "(*jamaratSol="];
        If[start === {}, Return[$Failed, Module]];
        stop = StringPosition[text, "|>;*)", IgnoreCase -> False];
        stop = Select[stop, #[[1]] > start[[1, 1]] &];
        If[stop === {}, Return[$Failed, Module]];
        snippet = StringTake[text, {start[[1, 1]] + StringLength["(*jamaratSol="], stop[[1, 1]] + 2}];
        snippet = StringTrim[snippet];
        If[StringEndsQ[snippet, ";"], snippet = StringDrop[snippet, -1]];
        ToExpression[snippet]
    ];

makeCase["jamarat-captured"] :=
    Module[{s, sys, sol},
        s = getExampleScenario[
            "Jamaratv9",
            {{1, 10}, {2, 50}},
            {{7, 0}, {8, 0}, {9, 55}}
        ];
        sys = makeSystem[s];
        sol = extractCapturedJamaratSolution[];
        If[sol === $Failed,
            Print["Could not find commented jamaratSol assignment in MFGraphs/Jamarat.wl"];
            Exit[2]
        ];
        <|"Scenario" -> s, "System" -> sys, "Solution" -> sol|>
    ];

makeCase[name_String] := (
    Print["Unknown case: ", name];
    Print["Known cases: chain-one-exit, chain-two-exits, example-7, jamarat-captured"];
    Exit[2]
);

solutionRules[sol_] :=
    If[ListQ[sol], sol, If[AssociationQ[sol], Lookup[sol, "Rules", {}], {}]];

solutionResidual[sol_] :=
    If[AssociationQ[sol], Lookup[sol, "Equations", True], True];

dnfBranches[expr_, timeout_] :=
    Module[{dnf},
        dnf = TimeConstrained[
            Quiet @ Check[BooleanConvert[expr, "DNF"], $Failed],
            timeout,
            $TimedOut
        ];
        Which[
            dnf === $TimedOut || dnf === $Failed, dnf,
            dnf === False, {},
            Head[dnf] === Or, List @@ dnf,
            True, {dnf}
        ]
    ];

entryCostData[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{pairs, vars, flows},
        pairs = systemData[sys, "InAuxEntryPairs"];
        vars = u @@@ pairs;
        flows = Last /@ scenarioData[s, "Model"]["Entries"];
        <|"Pairs" -> pairs, "Variables" -> vars, "Flows" -> flows|>
    ];

branchVariables[expr_] :=
    Select[
        Variables[expr /. {Equal -> List, Unequal -> List, LessEqual -> List,
            GreaterEqual -> List, Less -> List, Greater -> List, And -> List, Or -> List}],
        MatchQ[#, j[__] | u[__] | z[__]] &
    ];

branchEntryCostSummary[branch_, rules_List, entry_Association, minimizeTimeout_] :=
    Module[{entryVars, entryFlows, entryValues, objective, constraints, vars, min},
        entryVars = entry["Variables"];
        entryFlows = entry["Flows"];
        constraints = branch /. rules;
        entryValues = Simplify[entryVars /. rules];
        objective = Simplify[Total[entryFlows * entryValues]];
        vars = DeleteDuplicates @ branchVariables[{constraints, objective}];
        min = If[vars === {},
            If[TrueQ[constraints], {objective, {}}, $Failed],
            TimeConstrained[
                Quiet @ Check[
                    Minimize[{objective, constraints}, vars, Reals],
                    $Failed
                ],
                minimizeTimeout,
                $TimedOut
            ]
        ];
        <|
            "EntryValues" -> entryValues,
            "WeightedEntryCost" -> objective,
            "MinimizeResult" -> min,
            "ResidualVariableCount" -> Length[vars],
            "ResidualVariables" -> vars
        |>
    ];

analyzeSolutionBranches[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts_Association] :=
    Module[{rules, residual, entry, branches, summaries, dnfTimeout, minimizeTimeout, maxBranches},
        dnfTimeout = Lookup[opts, "DNFTimeout", 30];
        minimizeTimeout = Lookup[opts, "MinimizeTimeout", 10];
        maxBranches = Lookup[opts, "MaxBranches", 25];
        rules = solutionRules[sol];
        residual = Simplify[solutionResidual[sol] /. rules];
        entry = entryCostData[s, sys];
        branches = dnfBranches[residual, dnfTimeout];
        If[branches === $TimedOut || branches === $Failed,
            Return[<|
                "Status" -> branches,
                "ResultKind" -> solversTools`Private`solutionResultKind[sol],
                "EntryVariables" -> entry["Variables"],
                "RulesEntryValues" -> Simplify[entry["Variables"] /. rules],
                "RulesWeightedEntryCost" -> Simplify[Total[entry["Flows"] * (entry["Variables"] /. rules)]],
                "ResidualLeafCount" -> LeafCount[residual]
            |>]
        ];
        summaries = branchEntryCostSummary[#, rules, entry, minimizeTimeout] & /@ Take[branches, UpTo[maxBranches]];
        <|
            "Status" -> "OK",
            "ResultKind" -> solversTools`Private`solutionResultKind[sol],
            "EntryVariables" -> entry["Variables"],
            "EntryFlows" -> entry["Flows"],
            "BranchCount" -> Length[branches],
            "AnalyzedBranchCount" -> Length[summaries],
            "ResidualLeafCount" -> LeafCount[residual],
            "Branches" -> summaries
        |>
    ];

printBranchSummary[analysis_Association, printLimit_] :=
    Module[{branches, rows, ranked},
        Print["Result kind: ", Lookup[analysis, "ResultKind", Missing["NA"]]];
        Print["Entry variables: ", Lookup[analysis, "EntryVariables", {}]];
        If[Lookup[analysis, "Status", "OK"] =!= "OK",
            Print["Branch extraction status: ", analysis["Status"]];
            Print["Residual leaf count: ", Lookup[analysis, "ResidualLeafCount", Missing["NA"]]];
            Print["Entry values from rules only: ", Lookup[analysis, "RulesEntryValues", Missing["NA"]]];
            Print["Weighted entry cost from rules only: ", Lookup[analysis, "RulesWeightedEntryCost", Missing["NA"]]];
            Return[]
        ];
        Print["Branch count: ", analysis["BranchCount"], " (analyzed ", analysis["AnalyzedBranchCount"], ")"];
        Print["Residual leaf count: ", analysis["ResidualLeafCount"]];
        branches = analysis["Branches"];
        rows = MapIndexed[
            With[{b = #1, ix = First[#2]},
                <|
                    "Branch" -> ix,
                    "EntryValues" -> b["EntryValues"],
                    "WeightedEntryCost" -> b["WeightedEntryCost"],
                    "MinimizeResult" -> b["MinimizeResult"],
                    "ResidualVariableCount" -> b["ResidualVariableCount"]
                |>
            ] &,
            branches
        ];
        ranked = SortBy[
            rows,
            Replace[#["MinimizeResult"], {
                {val_, ___} :> {0, N[val]},
                _ :> {1, Infinity}
            }] &
        ];
        Print["Ranked branch summaries:"];
        Scan[Print, Take[ranked, UpTo[printLimit]]];
    ];

Module[{args, caseName, maxBranches, dnfTimeout, minimizeTimeout, printLimit, data, analysis},
    args = Rest[$ScriptCommandLine];
    caseName = If[Length[args] >= 1 && !StringStartsQ[args[[1]], "--"], args[[1]], "chain-two-exits"];
    maxBranches = ToExpression @ optionValue[args, "max-branches", "25"];
    dnfTimeout = ToExpression @ optionValue[args, "dnf-timeout", "30"];
    minimizeTimeout = ToExpression @ optionValue[args, "minimize-timeout", "10"];
    printLimit = ToExpression @ optionValue[args, "print-limit", ToString[maxBranches]];
    data = makeCase[caseName];
    analysis = analyzeSolutionBranches[
        data["Scenario"],
        data["System"],
        data["Solution"],
        <|
            "MaxBranches" -> maxBranches,
            "DNFTimeout" -> dnfTimeout,
            "MinimizeTimeout" -> minimizeTimeout
        |>
    ];
    printBranchSummary[analysis, printLimit];
];
