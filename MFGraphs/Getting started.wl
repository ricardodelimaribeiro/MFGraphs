(* ::Package:: *)

(* Notebook-friendly MFGraphs workbook covering the typed scenario kernels
   and the solversTools solver. *)

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
                "graphics`"
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

(* Use the global flag from primitives context *)
$MFGraphsVerbose = False;

(* --- 1. Build scenarios with ExampleScenarios constructors --- *)

sGrid = gridScenario[
    {3},
    {{1, 120.0}},
    {{3, 0.0}}
];

sCycle = cycleScenario[
    3,
    {{1, 50.0}},
    {{2, 0.0}, {3, 10.0}},
    {
        {1, 2, 3, 2.0},
        {3, 2, 1, 1.0}
    }
];

sAm = amScenario[
    {1, 2, 3, 4},
    {
        {0, 1, 1, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 1},
        {0, 0, 0, 0}
    },
    {{1, 80.0}},
    {{4, 0.0}}
];

Column[{
    DescribeOutput[
        "gridScenario output",
        "A typed scenario object with completed defaults and derived topology.",
        sGrid
    ],
    DescribeOutput[
        "cycleScenario output",
        "Cycle topology with explicit switching costs.",
        sCycle
    ],
    DescribeOutput[
        "amScenario output",
        "Adjacency-matrix constructor for explicit benchmark topologies.",
        sAm
    ]
}]


(* --- 2. Use named factory examples from ExampleScenarios.wl --- *)

exY = getExampleScenario[
    7,
    {{1, 100.0}},
    {{3, 0.0}, {4, 10.0}}
];

exGrid = getExampleScenario[
    "Grid0303",
    {{1, 30.0}},
    {{9, 0.0}}
];

Column[{
    DescribeOutput[
        "Named example (7)",
        "Factory-backed scenario with canonical topology and defaults.",
        exY
    ],
    DescribeOutput[
        "Named example (Grid0303)",
        "Grid example from the scenario registry.",
        exGrid
    ]
}]


(* --- 3. Inspect scenario structure from Scenario.wl --- *)

scenarioChecks = <|
    "scenarioQ[exY]" -> scenarioQ[exY],
    "Identity" -> scenarioData[exY, "Identity"],
    "Benchmark" -> scenarioData[exY, "Benchmark"],
    "Hamiltonian" -> scenarioData[exY, "Hamiltonian"],
    "Model keys" -> Keys @ scenarioData[exY, "Model"],
    "Topology keys" -> Keys @ scenarioData[exY, "Topology"]
|>;

DescribeOutput[
    "scenarioData overview",
    "Typed scenarios expose canonical blocks used by the downstream kernels.",
    scenarioChecks
]


(* --- 4. Build unknown bundles from Unknowns.wl --- *)

exampleUnknowns = makeUnknowns[exY];

unknownSummary = <|
    "unknownsQ" -> unknownsQ[exampleUnknowns],
    "js count" -> Length @ unknownsData[exampleUnknowns, "Js"],
    "jts count" -> Length @ unknownsData[exampleUnknowns, "Jts"],
    "us count" -> Length @ unknownsData[exampleUnknowns, "Us"],
    "first js" -> Take[unknownsData[exampleUnknowns, "Js"], UpTo[6]],
    "first jts" -> Take[unknownsData[exampleUnknowns, "Jts"], UpTo[6]],
    "first us" -> Take[unknownsData[exampleUnknowns, "Us"], UpTo[6]]
|>;

DescribeOutput[
    "Unknown bundle",
    "makeUnknowns derives flow, transition-flow, and value-function symbolic variables from scenario topology.",
    unknownSummary
]


(* --- 5. Build structural systems from System.wl --- *)

exampleSystem = makeSystem[exY, exampleUnknowns];

systemSummary = <|
    "mfgSystemQ" -> mfgSystemQ[exampleSystem],
    "System keys" -> Keys @ systemData[exampleSystem],
    "# EqEntryIn" -> Length @ systemData[exampleSystem, "EqEntryIn"],
    "# EqGeneral" -> Length @ systemData[exampleSystem, "EqGeneral"],
    "# AltTransitionFlows" -> Length @ systemData[exampleSystem, "AltTransitionFlows"],
    "# IneqSwitchingByVertex" -> Length @ systemData[exampleSystem, "IneqSwitchingByVertex"],
    "# js" -> Length @ systemData[exampleSystem, "Js"],
    "# jts" -> Length @ systemData[exampleSystem, "Jts"],
    "# us" -> Length @ systemData[exampleSystem, "Us"],
    "Switching-cost consistency" -> isSwitchingCostConsistent[
        Normal @ systemData[exampleSystem, "SwitchingCosts"]
    ]
|>;

DescribeOutput[
    "mfgSystem overview",
    "makeSystem builds structural equations and inequalities without running any solver.",
    systemSummary
]


(* --- 6. Chain with two exits: equations without and with switching costs --- *)

chain2ExNoSC = gridScenario[
    {3},
    {{1, 120}},
    {{2, 10}, {3, 0}}
];
chain2ExWithSC = gridScenario[
    {3},
    {{1, 120}},
    {{2, 10}, {3, 0}},
    {{1, 2, 3, 2}}
];

unkNoSC   = makeUnknowns[chain2ExNoSC];
sysNoSC   = makeSystem[chain2ExNoSC, unkNoSC];
unkWithSC = makeUnknowns[chain2ExWithSC];
sysWithSC = makeSystem[chain2ExWithSC, unkWithSC];

Column[{
    DescribeOutput[
        "Chain 1\[Rule]2\[Rule]3, exits at {2,3} \[LongDash] no switching costs: HJ equations",
        "One Hamilton-Jacobi equation per directed edge.",
        systemData[sysNoSC, "EqGeneral"]
    ],
    DescribeOutput[
        "Chain 1\[Rule]2\[Rule]3, exits at {2,3} \[LongDash] with SC {1,2,3}=2.0: HJ equations",
        "Same HJ equations; switching costs enter via complementarity, not HJ.",
        systemData[sysWithSC, "EqGeneral"]
    ],
    DescribeOutput[
        "Entry/exit boundary conditions (no SC)",
        "Entry flow pinned via EqEntryIn; exit value substitution rules from RuleExitValues.",
        Column[{systemData[sysNoSC, "EqEntryIn"], Normal @ systemData[sysNoSC, "RuleExitValues"]}]
    ],
    DescribeOutput[
        "Flow complementarity \[LongDash] no switching costs",
        "AltFlows: j[i,k]*j[k,i]=0 per edge pair. AltTransitionFlows: trivial (no jts).",
        Column[{systemData[sysNoSC, "AltFlows"], systemData[sysNoSC, "AltTransitionFlows"]}]
    ],
    DescribeOutput[
        "Flow complementarity \[LongDash] with switching costs",
        "IneqSwitchingByVertex: optimality at vertex 2. AltOptCond: combined condition.",
        Column[{systemData[sysWithSC, "IneqSwitchingByVertex"], systemData[sysWithSC, "AltOptCond"]}]
    ]
}]


(* --- 7. Topology visualization (from MFGraphs`Graphics`) --- *)

Column[{
    DescribeOutput[
        "Chain topology \[LongDash] no switching costs",
        "Green = entry, red = exits, gray = internal.",
        scenarioTopologyPlot[chain2ExNoSC, sysNoSC,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3 (no SC)"]
    ],
    DescribeOutput[
        "Chain topology \[LongDash] with switching cost {1,2,3}=2.0",
        "Same topology; SC at vertex 2 penalises continuing from edge 1\[Rule]2 to 2\[Rule]3.",
        scenarioTopologyPlot[chain2ExWithSC, sysWithSC,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3 (SC at 2)"]
    ]
}]


(* --- 8. Solve with reduceSystem --- *)

(* Chain 1->2->3, single exit at 3 (cost=0), entry flow=10.
   All variables are uniquely determined. *)
chain1Ex = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
sys1Ex   = makeSystem[chain1Ex];

DescribeOutput[
    "reduceSystem \[LongDash] chain with one exit",
    "Unique solution: all j and u values pinned by flow balance + HJ + complementarity.",
    reduceSystem[sys1Ex]
]


(* Chain 1->2->3, two exits at {2,0} and {3,10}, entry flow=120.
   Flow split between exits is under-determined; Reduce returns a parametric solution. *)
DescribeOutput[
    "reduceSystem \[LongDash] chain with two exits (no switching costs)",
    "Parametric solution: one free variable governs how flow splits between exits.",
    reduceSystem[sysNoSC]
]


(* --- 9. Visualize reduceSystem solutions (from MFGraphs`Graphics`) --- *)

(* --- Apply to chain with one exit --- *)

sol1Ex = reduceSystem[sys1Ex];

Column[{
    DescribeOutput[
        "Solution validation",
        "isValidSystemSolution confirms that the solved rules satisfy all system constraints.",
        isValidSystemSolution[sys1Ex, sol1Ex]
    ],
    DescribeOutput[
        "Combined solution plot \[LongDash] chain 1\[Rule]2\[Rule]3, single exit",
        "Directed edges show j-flow direction/magnitude; labels show both j and u. Auxiliary edges are included.",
        mfgSolutionPlot[chain1Ex, sys1Ex, sol1Ex,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: combined solution (exit at 3, cost=0)"]
    ],
    DescribeOutput[
        "Flow-only plot \[LongDash] chain 1\[Rule]2\[Rule]3, single exit",
        "Original vertices use one color. Edges have fixed width, and nearby labels show only j-flow values.",
        mfgFlowPlot[chain1Ex, sys1Ex, sol1Ex,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: flow values only"]
    ]
}]


(* --- Apply to chain with two exits (partial rules from underdetermined solution) --- *)

solNoSC = reduceSystem[sysNoSC];

Column[{
    DescribeOutput[
        "Combined solution plot \[LongDash] chain with two exits (no SC)",
        "Edge {1,2} is directed 1\[Rule]2 (j[1,2]>0). Labels show j and u on real + auxiliary edges.",
        mfgSolutionPlot[chain2ExNoSC, sysNoSC, solNoSC,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: combined solution (exits at 2 and 3)"]
    ],
    DescribeOutput[
        "Flow-only plot \[LongDash] chain with two exits (no SC)",
        "Auxiliary entry/exit edges are included. Edges have fixed width, and nearby labels show only j-flow values.",
        mfgFlowPlot[chain2ExNoSC, sysNoSC, solNoSC,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: flow values only (exits at 2 and 3)"]
    ]
}]


(* --- 10. Advanced Solution Visualization (Paper Scheme) --- *)

augChain1 = augmentAuxiliaryGraph[sys1Ex];

Column[{
    DescribeOutput[
        "Transition graph plot",
        "Nodes represent network edges (pairs); edges represent transition flows (j[a,b,c]).",
        mfgTransitionPlot[chain1Ex, sys1Ex, sol1Ex,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: transition graph"]
    ],
    DescribeOutput[
        "Augmented infrastructure helper",
        "augmentAuxiliaryGraph builds the road-traffic graph from AuxPairs and AuxTriples, with explicit j variables for every edge.",
        <|
            "Vertices" -> augChain1["Vertices"],
            "FlowEdges" -> augChain1["FlowEdges"],
            "TransitionEdges" -> augChain1["TransitionEdges"],
            "EdgeVariables" -> augChain1["EdgeVariables"]
        |>
    ],
    DescribeOutput[
        "Augmented infrastructure graph (Paper scheme)",
        "Nodes are road-traffic edge-vertex states. Blue edges are j[a,b] flows; red edges are j[r,i,w] transitions.",
        mfgAugmentedPlot[chain1Ex, sys1Ex, sol1Ex,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: augmented infrastructure"]
    ]
}]


(* --- 11. Paper-style network: Figure 3 \[Rule] Figure 4 transformation --- *)

paperFig3Scenario = graphScenario[
    Graph[{1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3, 3 \[UndirectedEdge] 4}],
    {{1, 10}, {3, 5}},
    {{2, 0}, {4, 0}}
];

paperFig3System = makeSystem[paperFig3Scenario];
paperFig3Augmented = augmentAuxiliaryGraph[paperFig3System];

paperFig3AuxGraph = Graph[
    VertexList @ systemData[paperFig3System, "AuxiliaryGraph"],
    EdgeList @ systemData[paperFig3System, "AuxiliaryGraph"],
    VertexLabels -> Placed["Name", Center],
    VertexStyle -> Normal @ Join[
        AssociationThread[scenarioData[paperFig3Scenario, "Model"]["Vertices"], GrayLevel[0.72]],
        AssociationThread[systemData[paperFig3System, "AuxEntryVertices"], RGBColor[0.38, 0.74, 0.9]],
        AssociationThread[systemData[paperFig3System, "AuxExitVertices"], RGBColor[0.95, 0.7, 0.4]]
    ],
    EdgeStyle -> Directive[GrayLevel[0.45], AbsoluteThickness[2]],
    (*GraphLayout -> "LayeredDigraphEmbedding",*)
    PlotLabel -> Style["Figure 3-style MFG network with entrance and exit edges", 14, Bold],
    ImageSize -> Large
];

Column[{
    DescribeOutput[
        "Figure 3-style auxiliary MFG network",
        "Original vertices are gray. Auxiliary entry vertices are blue, and auxiliary exit vertices are orange.",
        paperFig3AuxGraph
    ],
    DescribeOutput[
        "Figure 4-style augmented road-traffic graph",
        "Blue edges are flow variables j[a,b]. Red edges are transition variables j[r,i,w].",
        mfgAugmentedPlot[paperFig3Scenario, paperFig3System, <||>,
            PlotLabel -> "Paper-style augmented road-traffic graph"]
    ],
    DescribeOutput[
        "Augmented graph metadata",
        "The helper exposes the graph and the exact j variable attached to each road-traffic edge.",
        <|
            "FlowEdges" -> paperFig3Augmented["FlowEdges"],
            "TransitionEdges" -> paperFig3Augmented["TransitionEdges"],
            "EdgeVariables" -> paperFig3Augmented["EdgeVariables"]
        |>
    ]
}]
