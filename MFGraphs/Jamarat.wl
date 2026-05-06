(* ::Package:: *)

Quit[]


(* ::Subsection:: *)
(*Initialization*)


(* Notebook-friendly MFGraphs workbook for the Jamaratv9 scenario only. *)
(* Evaluate cells one at a time or section by section \[LongDash] do not evaluate the entire file at once. *)
(* Force clean reload \[LongDash] safe to re-evaluate without restarting the kernel. *)
(*Quiet[
    With[
        {
            reloadContexts = {
                "MFGraphs`",
                "primitives`",
                "scenarioTools`",
                "examples`",
                "unknownsTools`",
                "systemTools`",
                "solversTools`",
                "graphicsTools`"
            }
        },
        Remove["Global`*"];
        Scan[
            Remove[# <> "*", # <> "Private`*"]&,
            reloadContexts
        ];
        $Packages = DeleteCases[
            $Packages,
            Alternatives @@ Join[reloadContexts, (# <> "Private`") /@ reloadContexts]
        ];
        $ContextPath = DeleteCases[
            $ContextPath,
            Alternatives @@ Join[reloadContexts, (# <> "Private`") /@ reloadContexts]
        ];
    ];
];*)

(* This file lives alongside MFGraphs.wl in the same directory. *)
mfgDir = If[$InputFileName === "", 
    NotebookDirectory[], 
    DirectoryName[$InputFileName]
];

(* Robustness: ensure mfgDir is a string for ParentDirectory *)
If[!StringQ[mfgDir] || mfgDir === "", 
    mfgDir = ExpandFileName["."];
];

mfgParentDir = ParentDirectory[mfgDir];
If[!MemberQ[$Path, mfgParentDir],
    PrependTo[$Path, mfgParentDir]
];
Needs["MFGraphs`"];

ClearAll[DescribeOutput];
DescribeOutput[title_String, description_String, expr_] :=
    Column[
        {
            Style[title, 15, Bold, RGBColor[0.15, 0.26, 0.45]],
            Style[description, 11, GrayLevel[0.35]],
            expr
        },
        Alignment -> Left,
        Spacings -> {0.2, 0.7}
    ];

ClearAll[
    WorkbookFilePath,
    CapturedSolutionFromWorkbook,
    DNFBranches,
    EntryCostData,
    BranchVariables,
    BranchEntryCostSummary,
    RankSolutionBranches,
    JamaratScenarioSummary
];

WorkbookFilePath[] :=
    FileNameJoin[{mfgDir, "Jamarat.wl"}];

CapturedSolutionFromWorkbook[symbolName_String] :=
    Module[{path, text, start, stop, snippet},
        path = WorkbookFilePath[];
        text = Import[path, "Text"];
        start = StringPosition[text, "(*" <> symbolName <> "="];
        If[start === {}, Return[$Failed, Module]];
        stop = StringPosition[text, "|>;*)"];
        stop = Select[stop, #[[1]] > start[[1, 1]] &];
        If[stop === {}, Return[$Failed, Module]];
        snippet = StringTake[
            text,
            {start[[1, 1]] + StringLength["(*" <> symbolName <> "="], stop[[1, 1]] + 2}
        ];
        snippet = StringTrim[snippet];
        If[StringEndsQ[snippet, ";"], snippet = StringDrop[snippet, -1]];
        ToExpression[snippet]
    ];

DNFBranches[expr_, timeout_:60] :=
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

EntryCostData[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{pairs, vars, flows},
        pairs = systemData[sys, "InAuxEntryPairs"];
        vars = u @@@ pairs;
        flows = Last /@ scenarioData[s, "Model"]["Entries"];
        <|"Pairs" -> pairs, "Variables" -> vars, "Flows" -> flows|>
    ];

BranchVariables[expr_] :=
    Select[
        Variables[
            expr /. {Equal -> List, Unequal -> List, LessEqual -> List,
                GreaterEqual -> List, Less -> List, Greater -> List, And -> List, Or -> List}
        ],
        MatchQ[#, j[__] | u[__] | z[__]] &
    ];

BranchEntryCostSummary[branch_, rules_List, entry_Association, minimizeTimeout_:10] :=
    Module[{entryVars, entryFlows, entryValues, objective, constraints, vars, min},
        entryVars = entry["Variables"];
        entryFlows = entry["Flows"];
        constraints = branch /. rules;
        entryValues = Simplify[entryVars /. rules];
        objective = Simplify[Total[entryFlows * entryValues]];
        vars = DeleteDuplicates @ BranchVariables[{constraints, objective}];
        min = If[vars === {},
            If[TrueQ[constraints], {objective, {}}, $Failed],
            TimeConstrained[
                Quiet @ Check[Minimize[{objective, constraints}, vars, Reals], $Failed],
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

RankSolutionBranches[s_?scenarioQ, sys_?mfgSystemQ, sol_Association,
        maxBranches_:Infinity, dnfTimeout_:60, minimizeTimeout_:10] :=
    Module[{rules, residual, entry, branches, summaries, rows},
        rules = Lookup[sol, "Rules", {}];
        residual = Simplify[Lookup[sol, "Residual", True] /. rules];
        entry = EntryCostData[s, sys];
        branches = DNFBranches[residual, dnfTimeout];
        If[branches === $TimedOut || branches === $Failed,
            Return[
                <|
                    "Status" -> branches,
                    "EntryVariables" -> entry["Variables"],
                    "EntryFlows" -> entry["Flows"],
                    "ResidualLeafCount" -> LeafCount[residual]
                |>,
                Module
            ]
        ];
        summaries = BranchEntryCostSummary[#, rules, entry, minimizeTimeout] & /@
            Take[branches, UpTo[maxBranches]];
        rows = MapIndexed[
            Join[<|"Branch" -> First[#2]|>, #1] &,
            summaries
        ];
        <|
            "Status" -> "OK",
            "EntryVariables" -> entry["Variables"],
            "EntryFlows" -> entry["Flows"],
            "BranchCount" -> Length[branches],
            "AnalyzedBranchCount" -> Length[rows],
            "ResidualLeafCount" -> LeafCount[residual],
            "RankedBranches" -> SortBy[
                rows,
                Replace[#["MinimizeResult"], {
                    {val_, ___} :> {0, N[val]},
                    _ :> {1, Infinity}
                }] &
            ]
        |>
    ];

JamaratScenarioSummary[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{model, entries, exits, augmented},
        model = scenarioData[s, "Model"];
        entries = model["Entries"];
        exits = model["Exits"];
        augmented = augmentAuxiliaryGraph[sys];
        <|
            "Entries" -> entries,
            "EntryFlowTotal" -> Total[Last /@ entries],
            "Exits" -> exits,
            "ExitCostTotal" -> Total[Last /@ exits],
            "EntryFlowTotalGreaterThanExitCostTotal" -> (Total[Last /@ entries] > Total[Last /@ exits]),
            "NetworkEdges" -> Length[systemData[sys, "Edges"]],
            "FlowVariables" -> Length[systemData[sys, "Js"]],
            "TransitionFlowVariables" -> Length[systemData[sys, "Jts"]],
            "ValueVariables" -> Length[systemData[sys, "Us"]],
            "AugmentedFlowEdges" -> Length[augmented["FlowEdges"]],
            "AugmentedTransitionEdges" -> Length[augmented["TransitionEdges"]]
        |>
    ];

(* Use the global flag from primitives context *)
$MFGraphsVerbose = False;


(* ::Subsection:: *)
(*Simplified Jamarat cycle*)


(* ::Subsubsection:: *)
(*Simplified Jamarat end*)


(* ::Text:: *)
(*This is a shorter version of the Jamarat example for when the flows are at the last pillar. *)


jamaratEnd=cycleScenario[5, {{1,100},{5,100}},{{3,0},{4,40},{2,30}}];
jamaratEndSystem=makeSystem[jamaratEnd];
AbsoluteTiming[jamaratEndSol=solveScenario[jamaratEnd];]
Column[{
    DescribeOutput[
        "Simplified Jamarat augmented infrastructure",
        "Augmented graph before solving \[LongDash] structure only.",
        mfgAugmentedPlot[jamaratEnd, jamaratEndSystem, <||>,
            PlotLabel -> "Simplified Jamarat infrastructure",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Simplified Jamarat augmented solution",
        "Augmented graph with solved flow, transition, and u values.",
        mfgAugmentedPlot[jamaratEnd, jamaratEndSystem, jamaratEndSol,
            PlotLabel -> "Simplified Jamarat solution",
            ImageSize -> Large]
    ]
}]


jamaratEndSol


(* ::Subsubsection::Closed:: *)
(*Simplified Jamarat end*)


(* ::Text:: *)
(*This is a shorter version of the Jamarat example for when the flows are at the last pillar. *)


jamaratEnd=cycleScenario[5, {{1,200},{2,100}},{{3,10},{4,0},{5,10}}];


jamaratEndSystem = makeSystem[jamaratEnd];


AbsoluteTiming[jamaratEndSol=solveScenario[jamaratEnd];]


Column[{
    DescribeOutput[
        "Simplified Jamarat augmented infrastructure",
        "Augmented graph before solving \[LongDash] structure only.",
        mfgAugmentedPlot[jamaratEnd, jamaratEndSystem, <||>,
            PlotLabel -> "Simplified Jamarat infrastructure",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Simplified Jamarat augmented solution",
        "Augmented graph with solved flow, transition, and u values.",
        mfgAugmentedPlot[jamaratEnd, jamaratEndSystem, jamaratEndSol,
            PlotLabel -> "Simplified Jamarat solution",
            ImageSize -> Large]
    ]
}]


(* ::Subsection:: *)
(*Jamarat 9 vertices*)


(* --- 1. Captured Jamaratv9 run from the named scenario registry --- *)

jamaratScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 40}, {2, 50}},
    {{7, 3}, {8, 4}, {9, 2}}
];

jamaratSystem = makeSystem[jamaratScenario];
jamaratAugmented = augmentAuxiliaryGraph[jamaratSystem];

Column[{
    DescribeOutput[
        "Jamaratv9 scenario",
        "Named registry example with two entrances and three exits.",
        <|
            "scenarioQ" -> scenarioQ[jamaratScenario],
            "Entries" -> scenarioData[jamaratScenario, "Model"]["Entries"],
            "Exits" -> scenarioData[jamaratScenario, "Model"]["Exits"],
            "Network edges" -> systemData[jamaratSystem, "Edges"]
        |>
    ],
    DescribeOutput[
        "Jamaratv9 augmented road-traffic graph",
        "Augmented graph showing flow and transition edges before solving.",
        mfgAugmentedPlot[jamaratScenario, jamaratSystem, <||>,
            PlotLabel -> "Jamaratv9 augmented infrastructure",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Jamaratv9 augmented graph metadata",
        "augmentAuxiliaryGraph exposes the graph structure used by the plot.",
        <|
            "Vertex count" -> Length[jamaratAugmented["Vertices"]],
            "Flow edge count" -> Length[jamaratAugmented["FlowEdges"]],
            "Transition edge count" -> Length[jamaratAugmented["TransitionEdges"]],
            "First flow edges" -> Take[jamaratAugmented["FlowEdges"], UpTo[8]],
            "First transition edges" -> Take[jamaratAugmented["TransitionEdges"], UpTo[8]]
        |>
    ]
}]


jamaratScenarioSummary = JamaratScenarioSummary[jamaratScenario, jamaratSystem];

DescribeOutput[
    "Jamaratv9 scalar totals",
    "The captured run uses entry-flow total 60 and exit-cost total 55.",
    jamaratScenarioSummary
]


(* Optional rerun of the captured scenario. This may be slow; evaluate explicitly. *)
AbsoluteTiming[jamaratSol=solveScenario[jamaratScenario]]


If[Head[jamaratSol] =!= Association, jamaratSol = <|"Rules" -> jamaratSol, "Residual" -> True|>]


jamaratSolutionSummary = <|
    "CapturedSolutionQ" -> AssociationQ[jamaratSol],
    "ResultKind" -> solversTools`Private`solutionResultKind[jamaratSol],
    "RuleCount" -> Length[Lookup[jamaratSol, "Rules", {}]],
    "ResidualLeafCount" -> LeafCount[jamaratEquations],
    "EntryVariables" -> EntryCostData[jamaratScenario, jamaratSystem]["Variables"],
    "EntryValuesAfterRules" -> Simplify[
        EntryCostData[jamaratScenario, jamaratSystem]["Variables"] /. Lookup[jamaratSol, "Rules", {}]
    ]
|>;

DescribeOutput[
    "Captured Jamarat solution",
    "The large commented jamaratSol block above is loaded without rerunning solveScenario.",
    jamaratSolutionSummary
]


jamaratBranchRanking = RankSolutionBranches[
    jamaratScenario,
    jamaratSystem,
    jamaratSol,
    Infinity,
    60,
    10
];

DescribeOutput[
    "Jamarat residual branches ranked by entry cost",
    "The score is the fixed-flow weighted entry value: 10 u[auxEntry1,1] + 50 u[auxEntry2,2].",
    <|
        "Status" -> jamaratBranchRanking["Status"],
        "BranchCount" -> jamaratBranchRanking["BranchCount"],
        "AnalyzedBranchCount" -> jamaratBranchRanking["AnalyzedBranchCount"],
        "BestBranches" -> Take[jamaratBranchRanking["RankedBranches"], UpTo[10]]
    |>
]


jamaratBestBranchRules = Replace[
    First[Lookup[First[jamaratBranchRanking["RankedBranches"]], "MinimizeResult", {$Failed, {}}]],
    {
        _?NumericQ :> Last[First[jamaratBranchRanking["RankedBranches"]]["MinimizeResult"]],
        _ :> {}
    }
];

jamaratBestBranchSol = <|
    "Rules" -> Join[Lookup[jamaratSol, "Rules", {}], jamaratBestBranchRules],
    "Residual" -> True
|>;

Column[{
    DescribeOutput[
        "Jamaratv9 augmented infrastructure",
        "Augmented graph showing flow and transition edges before solving.",
        mfgAugmentedPlot[
            jamaratScenario,
            jamaratSystem,
            <||>,
            PlotLabel -> "Jamaratv9 augmented infrastructure",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, augmented",
        "Augmented graph with solved flow, transition, and u values for the best ranked branch.",
        mfgAugmentedPlot[
            jamaratScenario,
            jamaratSystem,
            jamaratBestBranchSol,
            PlotLabel -> "Jamaratv9 augmented solution, best ranked branch",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution validation",
        "The original captured branched solution remains the solver output being validated.",
        isValidSystemSolution[jamaratSystem, jamaratSol]
    ]
}]


(* --- 2. Jamaratv9 higher-entry-flow run --- *)

jamaratHighEntryScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 20}, {2, 50}},
    {{7, 0}, {8, 0}, {9, 55}}
];

jamaratHighEntrySystem = makeSystem[jamaratHighEntryScenario];
jamaratHighEntryAugmented = augmentAuxiliaryGraph[jamaratHighEntrySystem];
jamaratHighEntryScenarioSummary = JamaratScenarioSummary[
    jamaratHighEntryScenario,
    jamaratHighEntrySystem
];

Column[{
    DescribeOutput[
        "Jamaratv9 higher-entry-flow scenario",
        "This second run has entry-flow total 70, which is larger than exit-cost total 55.",
        jamaratHighEntryScenarioSummary
    ],
    DescribeOutput[
        "Jamaratv9 higher-entry-flow augmented graph",
        "Augmented graph \[LongDash] same shape as before, only boundary data changes.",
        mfgAugmentedPlot[
            jamaratHighEntryScenario,
            jamaratHighEntrySystem,
            <||>,
            PlotLabel -> "Jamaratv9 augmented infrastructure, entries {20,50}",
            ImageSize -> Large
        ]
    ]
}]


(* Optional expensive solve for the higher-entry-flow scenario. Evaluate this cell explicitly. *)
jamaratHighEntryRun[timeout_:Infinity] :=
    AbsoluteTiming[
        If[timeout === Infinity,
            solveScenario[jamaratHighEntryScenario],
            TimeConstrained[solveScenario[jamaratHighEntryScenario], timeout, $TimedOut]
        ]
    ];

(* Example:
jamaratHighEntryTimedSol = jamaratHighEntryRun[3600];
jamaratHighEntrySol = Last[jamaratHighEntryTimedSol];
jamaratHighEntryBranchRanking = If[
    AssociationQ[jamaratHighEntrySol],
    RankSolutionBranches[jamaratHighEntryScenario, jamaratHighEntrySystem, jamaratHighEntrySol, Infinity, 60, 10],
    Missing["NoSolution", jamaratHighEntrySol]
];
*)


(* --- 1. Captured Jamaratv9 run from the named scenario registry --- *)

jamaratScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 100}, {2, 100}},
    {{7, 0}, {8, 0}, {9, 0}}
];

jamaratSystem = makeSystem[jamaratScenario];
jamaratAugmented = augmentAuxiliaryGraph[jamaratSystem];

Column[{
    DescribeOutput[
        "Jamaratv9 scenario",
        "Named registry example with two entrances and three exits.",
        <|
            "scenarioQ" -> scenarioQ[jamaratScenario],
            "Entries" -> scenarioData[jamaratScenario, "Model"]["Entries"],
            "Exits" -> scenarioData[jamaratScenario, "Model"]["Exits"],
            "Network edges" -> systemData[jamaratSystem, "Edges"]
        |>
    ],
    DescribeOutput[
        "Jamaratv9 augmented road-traffic graph",
        "Augmented graph showing flow and transition edges before solving.",
        mfgAugmentedPlot[jamaratScenario, jamaratSystem, <||>,
            PlotLabel -> "Jamaratv9 augmented infrastructure",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Jamaratv9 augmented graph metadata",
        "augmentAuxiliaryGraph exposes the graph structure used by the plot.",
        <|
            "Vertex count" -> Length[jamaratAugmented["Vertices"]],
            "Flow edge count" -> Length[jamaratAugmented["FlowEdges"]],
            "Transition edge count" -> Length[jamaratAugmented["TransitionEdges"]],
            "First flow edges" -> Take[jamaratAugmented["FlowEdges"], UpTo[8]],
            "First transition edges" -> Take[jamaratAugmented["TransitionEdges"], UpTo[8]]
        |>
    ]
}]


jamaratRun[]
