(* ::Package:: *)

(* Notebook-friendly MFGraphs workbook for the Jamaratv9 scenario only. *)

(* Evaluate cells one at a time or section by section \[LongDash] do not evaluate the entire file at once. *)

(* Force clean reload \[LongDash] safe to re-evaluate without restarting the kernel. *)
Quiet[
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
];

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
        residual = Simplify[Lookup[sol, "Equations", True] /. rules];
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


(* --- 1. Captured Jamaratv9 run from the named scenario registry --- *)

jamaratScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 10}, {2, 50}},
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
        "Jamaratv9 topology plot",
        "scenarioTopologyPlot shows the original network with entry/exit coloring.",
        scenarioTopologyPlot[jamaratScenario, jamaratSystem,
            PlotLabel -> "Jamaratv9 topology",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Jamaratv9 augmented road-traffic graph",
        "mfgAugmentedPlot shows flow edges and transition edges without requiring a solved system.",
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
jamaratRun[] := AbsoluteTiming[solveScenario[jamaratScenario]]


jamaratRun[]


(* ::Input:: *)
(*jamaratSol=<|"Rules"->{u["auxExit7",7]->0,u["auxExit8",8]->0,u["auxExit9",9]->55,u[3,1]->u[2,1],u[4,2]->u[1,2],u[4,3]->u[1,3],u[5,3]->u[1,3],u[6,4]->u[2,4],u[3,4]->u[2,4],u[6,5]->u[3,5],u[7,5]->u[3,5],u[8,6]->u[4,6],u[5,6]->u[4,6],u[9,7]->u[5,7],u[9,8]->u[6,8],u[8,9]->u[7,9],j[1,2]->-50+j[2,1]+j[4,3]+j[2,4,6],j[2,4]->j[4,3]+j[2,4,6],j[3,1]->-60+j[1,3]+j[4,3]+j[2,4,6],j[3,4]->j[6,5]-j[2,4,6]+j[4,6,8],j[4,6]->j[6,5]+j[4,6,8],j[5,3]->-60+j[3,5]+j[6,5]+j[4,6,8],j[5,6]->j[6,8]-j[4,6,8],j[7,5]->-60+j[5,7]+j[6,8],j[7,"auxExit7"]->60-j[8,"auxExit8"],j[8,6]->0,j[8,9]->j[6,8]-j[8,"auxExit8"],j[9,7]->j[6,8]-j[8,"auxExit8"],j["auxEntry1",1]->10,j["auxEntry2",2]->50,j[1,2,4]->-50+j[2,1]+j[4,3]+j[2,4,6],j[1,3,4]->j[1,3]-j[1,3,5],j[2,1,3]->j[2,1],j[2,4,3]->j[4,3],j[3,1,2]->-60+j[1,3]+j[4,3]+j[2,4,6],j[3,4,2]->0,j[3,4,6]->j[6,5]-j[2,4,6]+j[4,6,8],j[3,5,6]->j[3,5]-j[3,5,7],j[4,2,1]->0,j[4,3,1]->-j[3,5]+j[4,3]+j[1,3,5],j[4,3,5]->j[3,5]-j[1,3,5],j[4,6,5]->j[6,5],j[5,3,1]->-60+j[1,3]+j[3,5]-j[1,3,5]+j[2,4,6],j[5,3,4]->-j[1,3]+j[6,5]+j[1,3,5]-j[2,4,6]+j[4,6,8],j[5,6,4]->0,j[5,6,8]->j[6,8]-j[4,6,8],j[5,7,9]->0,j[5,7,"auxExit7"]->j[5,7],j[6,4,3]->0,j[6,5,3]->-j[5,7]+j[6,5]+j[3,5,7],j[6,5,7]->j[5,7]-j[3,5,7],j[6,8,9]->j[6,8]-j[8,"auxExit8"],j[6,8,"auxExit8"]->j[8,"auxExit8"],j[7,5,3]->-60+j[3,5]+j[5,7]-j[3,5,7]+j[4,6,8],j[7,5,6]->-j[3,5]+j[6,8]+j[3,5,7]-j[4,6,8],j[7,9,8]->0,j[7,9,"auxExit9"]->0,j[8,6,5]->0,j[8,9,7]->j[6,8]-j[8,"auxExit8"],j[8,9,"auxExit9"]->0,j[9,7,5]->-60+j[5,7]+j[6,8],j[9,7,"auxExit7"]->60-j[5,7]-j[8,"auxExit8"],j[9,8,"auxExit8"]->0,j["auxEntry1",1,2]->10-j[1,3]+j[2,1],j["auxEntry1",1,3]->j[1,3]-j[2,1],j["auxEntry2",2,1]->j[2,1],j["auxEntry2",2,4]->50-j[2,1],j[4,2]->0,j[6,4]->0,j[7,9]->0,j[9,8]->0,j[9,"auxExit9"]->0,j[6,4,2]->0,j[8,6,4]->0,j[9,8,6]->0},"Equations"->(j[4,3]==0&&((j[6,5]==0&&((u[6,8]==0&&((j[5,7]==0&&j[6,8]==60&&j[8,"auxExit8"]==60&&j[3,5,7]==0&&u[4,6]==60&&((j[3,5]==0&&j[1,3,5]==0&&j[4,6,8]==60&&u[2,4]==120&&((j[2,1]==0&&((j[1,3]==0&&j[2,4,6]==60&&u[1,2]==180&&u[2,1]==190&&u["auxEntry1",1]==190&&u["auxEntry2",2]==180&&((u[1,3]==u[3,5]&&((u[1,3]==u[5,7]&&u[1,3]<=0&&(u[1,3]==u[7,9]||u[7,9]<=55))||(u[5,7]<=0&&u[7,9]<=55)))||(u[3,5]==u[5,7]&&u[3,5]<=0&&(u[3,5]==u[7,9]||u[7,9]<=55))||(u[5,7]==0&&u[7,9]==0)||(u[5,7]==u[7,9]&&u[5,7]<=0)))||(j[1,3]==10&&j[2,4,6]==50&&u[1,2]==170&&u[1,3]==130&&u[2,1]==140&&u["auxEntry1",1]==140&&u["auxEntry2",2]==170&&((u[5,7]<=0&&u[5,7]==u[7,9])||(u[5,7]==0&&u[7,9]==0)||(u[3,5]<=0&&u[3,5]==u[5,7]&&(u[3,5]==u[7,9]||u[7,9]<=55))))))||(((j[1,3]==60&&j[2,1]==50&&j[2,4,6]==0&&u[1,2]==290&&u[1,3]==180&&u[2,1]==240&&u["auxEntry1",1]==240&&u["auxEntry2",2]==290)||(2 j[1,3]==35&&2 j[2,1]==15&&2 j[2,4,6]==85&&2 u[1,2]==325&&2 u[1,3]==275&&u[2,1]==155&&u["auxEntry1",1]==155&&2 u["auxEntry2",2]==325))&&((u[5,7]<=0&&u[5,7]==u[7,9])||(u[5,7]==0&&u[7,9]==0)||(u[3,5]==u[5,7]&&u[3,5]<=0&&(u[3,5]==u[7,9]||u[7,9]<=55))))))||(((u[5,7]<=0&&u[5,7]==u[7,9])||(u[5,7]==0&&u[7,9]==0))&&((j[1,3]==60&&j[2,1]==50&&j[2,4,6]==0&&((j[3,5]==60&&j[1,3,5]==60&&j[4,6,8]==0&&u[1,2]==290&&u[1,3]==180&&u[2,1]==240&&u[3,5]==120&&u["auxEntry1",1]==240&&u["auxEntry2",2]==290)||(j[3,5]==30&&j[1,3,5]==30&&j[4,6,8]==30&&u[1,2]==230&&u[1,3]==120&&u[2,1]==180&&u[2,4]==90&&u[3,5]==90&&u["auxEntry1",1]==180&&u["auxEntry2",2]==230)))||(j[1,3]==10&&j[2,1]==0&&j[3,5]==10&&j[1,3,5]==10&&j[2,4,6]==50&&j[4,6,8]==50&&u[1,2]==160&&u[1,3]==80&&u[2,1]==90&&u[2,4]==110&&u[3,5]==70&&u["auxEntry1",1]==90&&u["auxEntry2",2]==160)||(3 u[1,2]==410&&3 u["auxEntry2",2]==410&&((3 j[1,3]==65&&3 j[2,1]==35&&3 j[3,5]==65&&3 j[1,3,5]==65&&3 j[2,4,6]==115&&3 j[4,6,8]==115&&3 u[1,3]==310&&u[2,1]==125&&3 u[2,4]==295&&3 u[3,5]==245&&u["auxEntry1",1]==125)||(3 j[1,3]==68&&3 j[2,1]==38&&3 j[3,5]==62&&3 j[1,3,5]==62&&3 j[2,4,6]==112&&3 j[4,6,8]==118&&3 u[1,3]==304&&u[2,1]==124&&3 u[2,4]==298&&3 u[3,5]==242&&u["auxEntry1",1]==124)))))))||(u[5,7]==0&&u[7,9]==0&&((j[1,3]==60&&j[2,1]==50&&j[2,4,6]==0&&((j[3,5]==60&&j[1,3,5]==60&&j[4,6,8]==0&&((j[5,7]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[3,5,7]==60&&u[1,2]==230&&u[1,3]==120&&u[2,1]==180&&u[3,5]==60&&u["auxEntry1",1]==180&&u["auxEntry2",2]==230)||(j[5,7]==40&&j[6,8]==20&&j[8,"auxExit8"]==20&&j[3,5,7]==40&&u[1,2]==210&&u[1,3]==100&&u[2,1]==160&&u[3,5]==40&&u[4,6]==20&&u["auxEntry1",1]==160&&u["auxEntry2",2]==210)))||(j[3,5]==36&&j[5,7]==36&&j[6,8]==24&&j[8,"auxExit8"]==24&&j[1,3,5]==36&&j[3,5,7]==36&&j[4,6,8]==24&&u[1,2]==182&&u[1,3]==72&&u[2,1]==132&&u[2,4]==48&&u[3,5]==36&&u[4,6]==24&&u["auxEntry1",1]==132&&u["auxEntry2",2]==182)||(11 j[3,5]==420&&11 j[5,7]==360&&11 j[6,8]==300&&11 j[8,"auxExit8"]==300&&11 j[1,3,5]==420&&11 j[3,5,7]==360&&11 j[4,6,8]==240&&11 u[1,2]==1990&&11 u[1,3]==780&&11 u[2,1]==1440&&11 u[2,4]==540&&11 u[3,5]==360&&11 u[4,6]==300&&11 u["auxEntry1",1]==1440&&11 u["auxEntry2",2]==1990)))||(j[1,3]==10&&j[2,1]==0&&j[3,5]==10&&j[5,7]==10&&j[6,8]==50&&j[8,"auxExit8"]==50&&j[1,3,5]==10&&j[2,4,6]==50&&j[3,5,7]==10&&j[4,6,8]==50&&u[1,2]==150&&u[1,3]==20&&u[2,1]==30&&u[2,4]==100&&u[3,5]==10&&u[4,6]==50&&u["auxEntry1",1]==30&&u["auxEntry2",2]==150)||(7 j[1,3]==190&&7 j[2,1]==120&&7 j[3,5]==190&&7 j[5,7]==190&&7 j[6,8]==230&&7 j[8,"auxExit8"]==230&&7 j[1,3,5]==190&&7 j[2,4,6]==230&&7 j[3,5,7]==190&&7 j[4,6,8]==230&&7 u[1,2]==690&&7 u[1,3]==380&&7 u[2,1]==570&&7 u[2,4]==460&&7 u[3,5]==190&&7 u[4,6]==230&&7 u["auxEntry1",1]==570&&7 u["auxEntry2",2]==690)))))||(j[1,3]==60&&j[2,1]==50&&j[3,5]==60&&j[5,7]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[1,3,5]==60&&j[2,4,6]==0&&j[3,5,7]==60&&j[4,6,8]==0&&u[1,2]==230&&u[1,3]==120&&u[2,1]==180&&u[3,5]==60&&u[5,7]==0&&u["auxEntry1",1]==180&&u["auxEntry2",2]==230&&((u[2,4]==u[4,6]&&((u[2,4]==u[6,8]&&u[2,4]<=0&&(u[2,4]==u[7,9]||u[7,9]<=55))||(u[6,8]<=0&&u[7,9]<=55)))||(u[4,6]==u[6,8]&&u[4,6]<=0&&(u[4,6]==u[7,9]||u[7,9]<=55))||(u[6,8]==u[7,9]&&u[6,8]<=0)))))||(u[5,7]==0&&((j[3,5]==0&&j[1,3,5]==0&&j[3,5,7]==0&&((j[2,1]==0&&((j[1,3]==0&&j[2,4,6]==60&&((u[6,8]==0&&u[7,9]==0&&((j[5,7]==60&&j[6,5]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[1,2]==240&&u[2,1]==250&&u[2,4]==180&&u[3,5]==60&&u[4,6]==120&&u["auxEntry1",1]==250&&u["auxEntry2",2]==240)||(j[5,7]==20&&j[6,5]==20&&j[6,8]==40&&j[8,"auxExit8"]==40&&j[4,6,8]==40&&u[1,2]==160&&u[2,1]==170&&u[2,4]==100&&u[3,5]==20&&u[4,6]==40&&u["auxEntry1",1]==170&&u["auxEntry2",2]==160)))||(j[5,7]==60&&j[6,5]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[1,2]==240&&u[2,1]==250&&u[2,4]==180&&u[3,5]==60&&u[4,6]==120&&u[6,8]<=0&&u[6,8]==u[7,9]&&u["auxEntry1",1]==250&&u["auxEntry2",2]==240)))||(j[1,3]==10&&j[2,4,6]==50&&((u[6,8]==0&&u[7,9]==0&&((j[5,7]==60&&j[6,5]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[1,2]==230&&u[1,3]==190&&u[2,1]==200&&u[2,4]==180&&u[3,5]==60&&u[4,6]==120&&u["auxEntry1",1]==200&&u["auxEntry2",2]==230)||(j[5,7]==20&&j[6,5]==20&&j[6,8]==40&&j[8,"auxExit8"]==40&&j[4,6,8]==40&&u[1,2]==150&&u[1,3]==110&&u[2,1]==120&&u[2,4]==100&&u[3,5]==20&&u[4,6]==40&&u["auxEntry1",1]==120&&u["auxEntry2",2]==150)))||(j[5,7]==60&&j[6,5]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[1,2]==230&&u[1,3]==190&&u[2,1]==200&&u[2,4]==180&&u[3,5]==60&&u[4,6]==120&&u[6,8]==u[7,9]&&u["auxEntry1",1]==200&&u["auxEntry2",2]==230&&u[6,8]<=0)))))||(u[6,8]==0&&u[7,9]==0&&((j[5,7]==60&&j[6,5]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[2,4]==180&&u[3,5]==60&&u[4,6]==120&&((j[1,3]==60&&j[2,1]==50&&j[2,4,6]==0&&u[1,2]==350&&u[1,3]==240&&u[2,1]==300&&u["auxEntry1",1]==300&&u["auxEntry2",2]==350)||(2 j[1,3]==35&&2 j[2,1]==15&&2 j[2,4,6]==85&&2 u[1,2]==445&&2 u[1,3]==395&&u[2,1]==215&&u["auxEntry1",1]==215&&2 u["auxEntry2",2]==445)))||(j[5,7]==20&&j[6,5]==20&&j[6,8]==40&&j[8,"auxExit8"]==40&&j[4,6,8]==40&&u[2,4]==100&&u[3,5]==20&&u[4,6]==40&&((j[1,3]==60&&j[2,1]==50&&j[2,4,6]==0&&u[1,2]==270&&u[1,3]==160&&u[2,1]==220&&u["auxEntry1",1]==220&&u["auxEntry2",2]==270)||(2 j[1,3]==35&&2 j[2,1]==15&&2 j[2,4,6]==85&&2 u[1,2]==285&&2 u[1,3]==235&&u[2,1]==135&&u["auxEntry1",1]==135&&2 u["auxEntry2",2]==285)))))||(j[5,7]==60&&j[6,5]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[2,4]==180&&u[3,5]==60&&u[4,6]==120&&u[6,8]<=0&&u[6,8]==u[7,9]&&((j[1,3]==60&&j[2,1]==50&&j[2,4,6]==0&&u[1,2]==350&&u[1,3]==240&&u[2,1]==300&&u["auxEntry1",1]==300&&u["auxEntry2",2]==350)||(2 j[1,3]==35&&2 j[2,1]==15&&2 j[2,4,6]==85&&2 u[1,2]==445&&2 u[1,3]==395&&u[2,1]==215&&u["auxEntry1",1]==215&&2 u["auxEntry2",2]==445)))))||(u[6,8]==0&&u[7,9]==0&&((j[1,3]==10&&j[2,1]==0&&j[3,5]==10&&j[1,3,5]==10&&j[2,4,6]==50&&j[3,5,7]==10&&((j[5,7]==60&&j[6,5]==50&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[1,2]==210&&u[1,3]==70&&u[2,1]==80&&u[2,4]==160&&u[3,5]==60&&u[4,6]==110&&u["auxEntry1",1]==80&&u["auxEntry2",2]==210)||(3 j[5,7]==70&&3 j[6,5]==40&&3 j[6,8]==110&&3 j[8,"auxExit8"]==110&&3 j[4,6,8]==110&&3 u[1,2]==410&&3 u[1,3]==100&&3 u[2,1]==130&&3 u[2,4]==260&&3 u[3,5]==70&&3 u[4,6]==110&&3 u["auxEntry1",1]==130&&3 u["auxEntry2",2]==410)))||(j[5,7]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[3,5]==60&&((j[1,3]==60&&j[2,1]==50&&j[3,5]==45&&j[6,5]==15&&j[1,3,5]==45&&j[2,4,6]==0&&j[3,5,7]==45&&u[1,2]==215&&u[1,3]==105&&u[2,1]==165&&u[2,4]==90&&u[4,6]==75&&u["auxEntry1",1]==165&&u["auxEntry2",2]==215)||(3 j[1,3]==95&&3 j[2,1]==65&&3 j[3,5]==95&&3 j[6,5]==85&&3 j[1,3,5]==95&&3 j[2,4,6]==85&&3 j[3,5,7]==95&&u[1,2]==145&&3 u[1,3]==275&&3 u[2,1]==370&&3 u[2,4]==350&&3 u[4,6]==265&&3 u["auxEntry1",1]==370&&u["auxEntry2",2]==145)))||(17 j[1,3]==450&&17 j[2,1]==280&&17 j[3,5]==450&&17 j[5,7]==490&&17 j[6,5]==40&&17 j[6,8]==530&&17 j[8,"auxExit8"]==530&&17 j[1,3,5]==450&&17 j[2,4,6]==570&&17 j[3,5,7]==450&&17 j[4,6,8]==530&&17 u[1,2]==1670&&17 u[1,3]==940&&17 u[2,1]==1390&&17 u[2,4]==1100&&17 u[3,5]==490&&17 u[4,6]==530&&17 u["auxEntry1",1]==1390&&17 u["auxEntry2",2]==1670)))||(j[5,7]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[3,5]==60&&u[6,8]<=0&&u[6,8]==u[7,9]&&((j[1,3]==10&&j[2,1]==0&&j[3,5]==10&&j[6,5]==50&&j[1,3,5]==10&&j[2,4,6]==50&&j[3,5,7]==10&&u[1,2]==210&&u[1,3]==70&&u[2,1]==80&&u[2,4]==160&&u[4,6]==110&&u["auxEntry1",1]==80&&u["auxEntry2",2]==210)||(j[1,3]==60&&j[2,1]==50&&j[3,5]==45&&j[6,5]==15&&j[1,3,5]==45&&j[2,4,6]==0&&j[3,5,7]==45&&u[1,2]==215&&u[1,3]==105&&u[2,1]==165&&u[2,4]==90&&u[4,6]==75&&u["auxEntry1",1]==165&&u["auxEntry2",2]==215)||(3 j[1,3]==95&&3 j[2,1]==65&&3 j[3,5]==95&&3 j[6,5]==85&&3 j[1,3,5]==95&&3 j[2,4,6]==85&&3 j[3,5,7]==95&&u[1,2]==145&&3 u[1,3]==275&&3 u[2,1]==370&&3 u[2,4]==350&&3 u[4,6]==265&&3 u["auxEntry1",1]==370&&u["auxEntry2",2]==145)))))))||(j[6,5]==0&&((j[3,5]==60&&j[2,4,6]==0&&j[4,6,8]==0&&((u[5,7]==0&&((j[5,7]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[3,5,7]==60&&u[1,3]==120&&u[3,5]==60&&((j[2,1]==0&&((j[1,3]==0&&j[4,3]==60&&j[1,3,5]==0&&u[1,2]==240&&u[2,1]==250&&u[2,4]==180&&u["auxEntry1",1]==250&&u["auxEntry2",2]==240)||(j[1,3]==10&&j[4,3]==50&&j[1,3,5]==10&&u[1,2]==220&&u[2,1]==130&&u[2,4]==170&&u["auxEntry1",1]==130&&u["auxEntry2",2]==220)))||(2 j[1,3]==65&&2 j[2,1]==45&&2 j[4,3]==55&&2 j[1,3,5]==65&&u[1,2]==175&&2 u[2,1]==305&&2 u[2,4]==295&&2 u["auxEntry1",1]==305&&u["auxEntry2",2]==175))&&((u[6,8]==0&&u[7,9]==0)||(u[6,8]<=0&&u[6,8]==u[7,9])||(u[4,6]==u[6,8]&&u[4,6]<=0&&(u[4,6]==u[7,9]||u[7,9]<=55))))||(u[6,8]==0&&u[7,9]==0&&((j[5,7]==0&&j[6,8]==60&&j[8,"auxExit8"]==60&&j[3,5,7]==0&&u[1,3]==180&&u[3,5]==120&&u[4,6]==60&&((j[2,1]==0&&((j[1,3]==0&&j[4,3]==60&&j[1,3,5]==0&&u[1,2]==300&&u[2,1]==310&&u[2,4]==240&&u["auxEntry1",1]==310&&u["auxEntry2",2]==300)||(j[1,3]==10&&j[4,3]==50&&j[1,3,5]==10&&u[1,2]==280&&u[2,1]==190&&u[2,4]==230&&u["auxEntry1",1]==190&&u["auxEntry2",2]==280)))||(2 j[1,3]==65&&2 j[2,1]==45&&2 j[4,3]==55&&2 j[1,3,5]==65&&u[1,2]==235&&2 u[2,1]==425&&2 u[2,4]==415&&2 u["auxEntry1",1]==425&&u["auxEntry2",2]==235)))||(j[5,7]==40&&j[6,8]==20&&j[8,"auxExit8"]==20&&j[3,5,7]==40&&u[1,3]==100&&u[3,5]==40&&u[4,6]==20&&((j[2,1]==0&&((j[1,3]==0&&j[4,3]==60&&j[1,3,5]==0&&u[1,2]==220&&u[2,1]==230&&u[2,4]==160&&u["auxEntry1",1]==230&&u["auxEntry2",2]==220)||(j[1,3]==10&&j[4,3]==50&&j[1,3,5]==10&&u[1,2]==200&&u[2,1]==110&&u[2,4]==150&&u["auxEntry1",1]==110&&u["auxEntry2",2]==200)))||(2 j[1,3]==65&&2 j[2,1]==45&&2 j[4,3]==55&&2 j[1,3,5]==65&&u[1,2]==155&&2 u[2,1]==265&&2 u[2,4]==255&&2 u["auxEntry1",1]==265&&u["auxEntry2",2]==155)))))))||(j[5,7]==0&&j[6,8]==60&&j[8,"auxExit8"]==60&&j[3,5,7]==0&&u[1,3]==180&&u[3,5]==120&&u[4,6]==60&&u[5,7]<=0&&u[6,8]==0&&u[5,7]==u[7,9]&&((j[2,1]==0&&((j[1,3]==0&&j[4,3]==60&&j[1,3,5]==0&&u[1,2]==300&&u[2,1]==310&&u[2,4]==240&&u["auxEntry1",1]==310&&u["auxEntry2",2]==300)||(j[1,3]==10&&j[4,3]==50&&j[1,3,5]==10&&u[1,2]==280&&u[2,1]==190&&u[2,4]==230&&u["auxEntry1",1]==190&&u["auxEntry2",2]==280)))||(2 j[1,3]==65&&2 j[2,1]==45&&2 j[4,3]==55&&2 j[1,3,5]==65&&u[1,2]==235&&2 u[2,1]==425&&2 u[2,4]==415&&2 u["auxEntry1",1]==425&&u["auxEntry2",2]==235)))))||(u[6,8]==0&&((j[2,1]==0&&((j[5,7]==0&&j[6,8]==60&&j[8,"auxExit8"]==60&&j[3,5,7]==0&&u[4,6]==60&&((u[5,7]<=0&&u[5,7]==u[7,9])||(u[5,7]==0&&u[7,9]==0))&&((j[1,3]==0&&j[3,5]==15&&j[4,3]==15&&j[1,3,5]==0&&j[2,4,6]==45&&j[4,6,8]==45&&u[1,2]==165&&u[1,3]==90&&u[2,1]==175&&u[2,4]==105&&u[3,5]==75&&u["auxEntry1",1]==175&&u["auxEntry2",2]==165)||(j[1,3]==10&&2 j[3,5]==35&&2 j[4,3]==15&&j[1,3,5]==10&&2 j[2,4,6]==85&&2 j[4,6,8]==85&&2 u[1,2]==305&&u[1,3]==95&&u[2,1]==105&&2 u[2,4]==205&&2 u[3,5]==155&&u["auxEntry1",1]==105&&2 u["auxEntry2",2]==305)))||(u[5,7]==0&&u[7,9]==0&&((j[1,3]==0&&j[3,5]==24&&j[4,3]==24&&j[5,7]==24&&j[6,8]==36&&j[8,"auxExit8"]==36&&j[1,3,5]==0&&j[2,4,6]==36&&j[3,5,7]==24&&j[4,6,8]==36&&u[1,2]==132&&u[1,3]==48&&u[2,1]==142&&u[2,4]==72&&u[3,5]==24&&u[4,6]==36&&u["auxEntry1",1]==142&&u["auxEntry2",2]==132)||(j[1,3]==10&&j[3,5]==26&&j[4,3]==16&&j[5,7]==26&&j[6,8]==34&&j[8,"auxExit8"]==34&&j[1,3,5]==10&&j[2,4,6]==34&&j[3,5,7]==26&&j[4,6,8]==34&&u[1,2]==118&&u[1,3]==52&&u[2,1]==62&&u[2,4]==68&&u[3,5]==26&&u[4,6]==34&&u["auxEntry1",1]==62&&u["auxEntry2",2]==118)))))||(19 j[1,3]==470&&19 j[2,1]==280&&19 j[3,5]==550&&19 j[4,3]==80&&19 j[5,7]==550&&19 j[6,8]==590&&19 j[8,"auxExit8"]==590&&19 j[1,3,5]==470&&19 j[2,4,6]==590&&19 j[3,5,7]==550&&19 j[4,6,8]==590&&19 u[1,2]==1850&&19 u[1,3]==1100&&19 u[2,1]==1570&&19 u[2,4]==1180&&19 u[3,5]==550&&19 u[4,6]==590&&u[5,7]==0&&u[7,9]==0&&19 u["auxEntry1",1]==1570&&19 u["auxEntry2",2]==1850)))))||(u[5,7]==0&&((u[6,8]==0&&u[7,9]==0&&((j[2,1]==0&&((j[1,3]==0&&j[1,3,5]==0&&((j[3,5]==30&&j[4,3]==30&&j[5,7]==60&&j[6,5]==30&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[2,4,6]==30&&j[3,5,7]==30&&j[4,6,8]==0&&u[1,2]==180&&u[1,3]==90&&u[2,1]==190&&u[2,4]==120&&u[3,5]==60&&u[4,6]==90&&u["auxEntry1",1]==190&&u["auxEntry2",2]==180)||(11 j[3,5]==240&&11 j[4,3]==240&&11 j[5,7]==300&&11 j[6,5]==60&&11 j[6,8]==360&&11 j[8,"auxExit8"]==360&&11 j[2,4,6]==420&&11 j[3,5,7]==240&&11 j[4,6,8]==360&&11 u[1,2]==1440&&11 u[1,3]==540&&11 u[2,1]==1550&&11 u[2,4]==780&&11 u[3,5]==300&&11 u[4,6]==360&&11 u["auxEntry1",1]==1550&&11 u["auxEntry2",2]==1440)))||(j[1,3]==10&&j[1,3,5]==10&&((2 j[3,5]==65&&2 j[4,3]==45&&j[5,7]==60&&2 j[6,5]==55&&j[6,8]==0&&j[8,"auxExit8"]==0&&2 j[2,4,6]==55&&2 j[3,5,7]==65&&j[4,6,8]==0&&u[1,2]==165&&2 u[1,3]==185&&2 u[2,1]==205&&u[2,4]==115&&u[3,5]==60&&2 u[4,6]==175&&2 u["auxEntry1",1]==205&&u["auxEntry2",2]==165)||(11 j[3,5]==270&&11 j[4,3]==160&&11 j[5,7]==310&&11 j[6,5]==40&&11 j[6,8]==350&&11 j[8,"auxExit8"]==350&&11 j[2,4,6]==390&&11 j[3,5,7]==270&&11 j[4,6,8]==350&&11 u[1,2]==1290&&11 u[1,3]==580&&11 u[2,1]==690&&11 u[2,4]==740&&11 u[3,5]==310&&11 u[4,6]==350&&11 u["auxEntry1",1]==690&&11 u["auxEntry2",2]==1290)))))||(3 j[1,3]==80&&3 j[2,1]==50&&3 j[3,5]==110&&j[4,3]==10&&j[5,7]==60&&3 j[6,5]==70&&j[6,8]==0&&j[8,"auxExit8"]==0&&3 j[1,3,5]==80&&3 j[2,4,6]==70&&3 j[3,5,7]==110&&j[4,6,8]==0&&u[1,2]==140&&3 u[1,3]==290&&3 u[2,1]==370&&3 u[2,4]==320&&u[3,5]==60&&3 u[4,6]==250&&3 u["auxEntry1",1]==370&&u["auxEntry2",2]==140)||(41 j[1,3]==1010&&41 j[2,1]==600&&41 j[3,5]==1170&&41 j[4,3]==160&&41 j[5,7]==1210&&41 j[6,5]==40&&41 j[6,8]==1250&&41 j[8,"auxExit8"]==1250&&41 j[1,3,5]==1010&&41 j[2,4,6]==1290&&41 j[3,5,7]==1170&&41 j[4,6,8]==1250&&41 u[1,2]==3990&&41 u[1,3]==2380&&41 u[2,1]==3390&&41 u[2,4]==2540&&41 u[3,5]==1210&&41 u[4,6]==1250&&41 u["auxEntry1",1]==3390&&41 u["auxEntry2",2]==3990)))||(j[5,7]==60&&j[6,8]==0&&j[8,"auxExit8"]==0&&j[4,6,8]==0&&u[3,5]==60&&u[6,8]<=0&&u[6,8]==u[7,9]&&((j[2,1]==0&&((j[1,3]==0&&j[3,5]==30&&j[4,3]==30&&j[6,5]==30&&j[1,3,5]==0&&j[2,4,6]==30&&j[3,5,7]==30&&u[1,2]==180&&u[1,3]==90&&u[2,1]==190&&u[2,4]==120&&u[4,6]==90&&u["auxEntry1",1]==190&&u["auxEntry2",2]==180)||(j[1,3]==10&&2 j[3,5]==65&&2 j[4,3]==45&&2 j[6,5]==55&&j[1,3,5]==10&&2 j[2,4,6]==55&&2 j[3,5,7]==65&&u[1,2]==165&&2 u[1,3]==185&&2 u[2,1]==205&&u[2,4]==115&&2 u[4,6]==175&&2 u["auxEntry1",1]==205&&u["auxEntry2",2]==165)))||(3 j[1,3]==80&&3 j[2,1]==50&&3 j[3,5]==110&&j[4,3]==10&&3 j[6,5]==70&&3 j[1,3,5]==80&&3 j[2,4,6]==70&&3 j[3,5,7]==110&&u[1,2]==140&&3 u[1,3]==290&&3 u[2,1]==370&&3 u[2,4]==320&&3 u[4,6]==250&&3 u["auxEntry1",1]==370&&u["auxEntry2",2]==140)))))|>;*)


jamaratSol = CapturedSolutionFromWorkbook["jamaratSol"];

jamaratEquations = Lookup[jamaratSol, "Equations", True];

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
    "Equations" -> True
|>;

Column[{
    DescribeOutput[
        "Jamaratv9 problem topology",
        "scenarioTopologyPlot shows the real network with entry and exit vertices highlighted.",
        scenarioTopologyPlot[
            jamaratScenario,
            jamaratSystem,
            PlotLabel -> "Jamaratv9 topology",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 problem augmented infrastructure",
        "mfgAugmentedPlot shows edge-flow and transition-flow variables on the augmented road-traffic graph.",
        mfgAugmentedPlot[
            jamaratScenario,
            jamaratSystem,
            <||>,
            PlotLabel -> "Jamaratv9 augmented infrastructure",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, combined",
        "mfgSolutionPlot combines flow, value, and density views for the best ranked branch of the captured solution.",
        mfgSolutionPlot[
            jamaratScenario,
            jamaratSystem,
            jamaratBestBranchSol,
            PlotLabel -> "Jamaratv9 captured solution, best ranked branch",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, flow",
        "mfgFlowPlot shows directed edge and auxiliary flows for the best ranked branch.",
        mfgFlowPlot[
            jamaratScenario,
            jamaratSystem,
            jamaratBestBranchSol,
            PlotLabel -> "Jamaratv9 flow, best ranked branch",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, values",
        "mfgValuePlot shows endpoint and interpolated value samples for the best ranked branch.",
        mfgValuePlot[
            jamaratScenario,
            jamaratSystem,
            jamaratBestBranchSol,
            PlotLabel -> "Jamaratv9 values, best ranked branch",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, densities",
        "mfgDensityPlot shows inferred edge densities for the best ranked branch.",
        mfgDensityPlot[
            jamaratScenario,
            jamaratSystem,
            jamaratBestBranchSol,
            PlotLabel -> "Jamaratv9 densities, best ranked branch",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, transitions",
        "mfgTransitionPlot shows transition-flow movement between directed edge states.",
        mfgTransitionPlot[
            jamaratScenario,
            jamaratSystem,
            jamaratBestBranchSol,
            PlotLabel -> "Jamaratv9 transitions, best ranked branch",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 captured solution, augmented",
        "mfgAugmentedPlot overlays solved flow and transition values on the augmented infrastructure graph.",
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
        "Jamaratv9 higher-entry-flow topology",
        "Same Jamaratv9 network with a larger first entrance flow.",
        scenarioTopologyPlot[
            jamaratHighEntryScenario,
            jamaratHighEntrySystem,
            PlotLabel -> "Jamaratv9 topology, entries {20,50}",
            ImageSize -> Large
        ]
    ],
    DescribeOutput[
        "Jamaratv9 higher-entry-flow augmented graph",
        "The augmented graph shape is unchanged; only boundary data changes.",
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
        "Jamaratv9 topology plot",
        "scenarioTopologyPlot shows the original network with entry/exit coloring.",
        scenarioTopologyPlot[jamaratScenario, jamaratSystem,
            PlotLabel -> "Jamaratv9 topology",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Jamaratv9 augmented road-traffic graph",
        "mfgAugmentedPlot shows flow edges and transition edges without requiring a solved system.",
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


symJamarat=%[[2]]
