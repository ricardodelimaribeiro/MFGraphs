(* ::Package:: *)

(* Notebook-friendly MFGraphs workbook covering the typed scenario kernels
   and the ReduceSystem solver. *)

(* Evaluate cells one at a time or section by section \[LongDash] do not evaluate the entire file at once. *)

(* This file lives alongside MFGraphs.wl in the same directory. *)
mfgDir = If[$InputFileName === "", NotebookDirectory[], DirectoryName[$InputFileName]];

(* Force clean reload \[LongDash] safe to re-evaluate without restarting the kernel. *)
Quiet[
    Unprotect /@ Names["MFGraphs`*"];
    Remove /@ Names["MFGraphs`*"];
    Unprotect[$Packages];
    $Packages = DeleteCases[$Packages, "MFGraphs`"];
    Protect[$Packages]
];
PrependTo[$Path, ParentDirectory[mfgDir]];
Get[FileNameJoin[{mfgDir, "MFGraphs.wl"}]];



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

MFGraphs`$MFGraphsVerbose = False;

(* --- 1. Build scenarios with ExampleScenarios constructors --- *)

gridScenario = GridScenario[
    {3},
    {{1, 120.0}},
    {{3, 0.0}}
];

cycleScenario = CycleScenario[
    3,
    {{1, 50.0}},
    {{2, 0.0}, {3, 10.0}},
    {
        {1, 2, 3, 2.0},
        {3, 2, 1, 1.0}
    }
];

amScenario = AMScenario[
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
        "GridScenario output",
        "A typed scenario object with completed defaults and derived topology.",
        gridScenario
    ],
    DescribeOutput[
        "CycleScenario output",
        "Cycle topology with explicit switching costs.",
        cycleScenario
    ],
    DescribeOutput[
        "AMScenario output",
        "Adjacency-matrix constructor for explicit benchmark topologies.",
        amScenario
    ]
}]


(* --- 2. Use named factory examples from ExampleScenarios.wl --- *)

exampleY = GetExampleScenario[
    7,
    {{1, 100.0}},
    {{3, 0.0}, {4, 10.0}}
];

exampleGrid = GetExampleScenario[
    "Grid0303",
    {{1, 30.0}},
    {{9, 0.0}}
];

Column[{
    DescribeOutput[
        "Named example (7)",
        "Factory-backed scenario with canonical topology and defaults.",
        exampleY
    ],
    DescribeOutput[
        "Named example (Grid0303)",
        "Grid example from the scenario registry.",
        exampleGrid
    ]
}]


(* --- 3. Inspect scenario structure from Scenario.wl --- *)

scenarioChecks = <|
    "scenarioQ[exampleY]" -> scenarioQ[exampleY],
    "Identity" -> ScenarioData[exampleY, "Identity"],
    "Benchmark" -> ScenarioData[exampleY, "Benchmark"],
    "Hamiltonian" -> ScenarioData[exampleY, "Hamiltonian"],
    "Model keys" -> Keys @ ScenarioData[exampleY, "Model"],
    "Topology keys" -> Keys @ ScenarioData[exampleY, "Topology"]
|>;

DescribeOutput[
    "ScenarioData overview",
    "Typed scenarios expose canonical blocks used by the downstream kernels.",
    scenarioChecks
]


(* --- 4. Build unknown bundles from Unknowns.wl --- *)

exampleUnknowns = makeUnknowns[exampleY];

unknownSummary = <|
    "unknownsQ" -> unknownsQ[exampleUnknowns],
    "js count" -> Length @ UnknownsData[exampleUnknowns, "js"],
    "jts count" -> Length @ UnknownsData[exampleUnknowns, "jts"],
    "us count" -> Length @ UnknownsData[exampleUnknowns, "us"],
    "first js" -> Take[UnknownsData[exampleUnknowns, "js"], UpTo[6]],
    "first jts" -> Take[UnknownsData[exampleUnknowns, "jts"], UpTo[6]],
    "first us" -> Take[UnknownsData[exampleUnknowns, "us"], UpTo[6]]
|>;

DescribeOutput[
    "Unknown bundle",
    "makeUnknowns derives flow, transition-flow, and value-function symbolic variables from scenario topology.",
    unknownSummary
]


(* --- 5. Build structural systems from System.wl --- *)

exampleSystem = makeSystem[exampleY, exampleUnknowns];

systemSummary = <|
    "mfgSystemQ" -> mfgSystemQ[exampleSystem],
    "System keys" -> Keys @ SystemData[exampleSystem],
    "# EqEntryIn" -> Length @ SystemData[exampleSystem, "EqEntryIn"],
    "# EqGeneral" -> Length @ SystemData[exampleSystem, "EqGeneral"],
    "# AltTransitionFlows" -> Length @ SystemData[exampleSystem, "AltTransitionFlows"],
    "# IneqSwitchingByVertex" -> Length @ SystemData[exampleSystem, "IneqSwitchingByVertex"],
    "# js" -> Length @ SystemData[exampleSystem, "js"],
    "# jts" -> Length @ SystemData[exampleSystem, "jts"],
    "# us" -> Length @ SystemData[exampleSystem, "us"],
    "Switching-cost consistency" -> IsSwitchingCostConsistent[
        Normal @ SystemData[exampleSystem, "SwitchingCosts"]
    ]
|>;

DescribeOutput[
    "mfgSystem overview",
    "makeSystem builds structural equations and inequalities without running any solver.",
    systemSummary
]


(* --- 6. Chain with two exits: equations without and with switching costs --- *)

chain2ExNoSC = GridScenario[
    {3},
    {{1, 120.0}},
    {{2, 0.0}, {3, 10.0}}
];
chain2ExWithSC = GridScenario[
    {3},
    {{1, 120.0}},
    {{2, 0.0}, {3, 10.0}},
    {{1, 2, 3, 2.0}}
];

unkNoSC   = makeUnknowns[chain2ExNoSC];
sysNoSC   = makeSystem[chain2ExNoSC, unkNoSC];
unkWithSC = makeUnknowns[chain2ExWithSC];
sysWithSC = makeSystem[chain2ExWithSC, unkWithSC];

Column[{
    DescribeOutput[
        "Chain 1\[Rule]2\[Rule]3, exits at {2,3} \[LongDash] no switching costs: HJ equations",
        "One Hamilton-Jacobi equation per directed edge.",
        SystemData[sysNoSC, "EqGeneral"]
    ],
    DescribeOutput[
        "Chain 1\[Rule]2\[Rule]3, exits at {2,3} \[LongDash] with SC {1,2,3}=2.0: HJ equations",
        "Same HJ equations; switching costs enter via complementarity, not HJ.",
        SystemData[sysWithSC, "EqGeneral"]
    ],
    DescribeOutput[
        "Entry/exit boundary conditions (no SC)",
        "Entry flow pinned via EqEntryIn; exit value substitution rules from RuleExitValues.",
        Column[{SystemData[sysNoSC, "EqEntryIn"], Normal @ SystemData[sysNoSC, "RuleExitValues"]}]
    ],
    DescribeOutput[
        "Flow complementarity \[LongDash] no switching costs",
        "AltFlows: j[i,k]*j[k,i]=0 per edge pair. AltTransitionFlows: trivial (no jts).",
        Column[{SystemData[sysNoSC, "AltFlows"], SystemData[sysNoSC, "AltTransitionFlows"]}]
    ],
    DescribeOutput[
        "Flow complementarity \[LongDash] with switching costs",
        "IneqSwitchingByVertex: optimality at vertex 2. AltOptCond: combined condition.",
        Column[{SystemData[sysWithSC, "IneqSwitchingByVertex"], SystemData[sysWithSC, "AltOptCond"]}]
    ]
}]


(* --- 7. Topology visualization --- *)

ClearAll[ScenarioTopologyPlot];

ScenarioTopologyPlot[s_?scenarioQ, sys_?mfgSystemQ, title_: Automatic] :=
    Module[{model, entryV, exitV, internalV, allV, edges, plotTitle},
        model    = ScenarioData[s, "Model"];
        entryV   = First /@ model["Entrance Vertices and Flows"];
        exitV    = First /@ model["Exit Vertices and Terminal Costs"];
        allV     = model["Vertices List"];
        internalV = Complement[allV, entryV, exitV];
        edges    = SystemData[sys, "edgeList"];
        plotTitle = Replace[title, Automatic -> "Network topology"];
        Graph[allV, edges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> Normal @ Join[
                AssociationThread[entryV,    RGBColor[0.22, 0.6, 0.3]],
                AssociationThread[exitV,     RGBColor[0.82, 0.27, 0.2]],
                AssociationThread[internalV, GrayLevel[0.75]]
            ],
            VertexSize -> 0.3,
            EdgeStyle  -> Directive[GrayLevel[0.5], AbsoluteThickness[2]],
            GraphLayout -> "LayeredDigraphEmbedding",
            PlotLabel  -> Style[plotTitle, 14, Bold],
            ImageSize  -> Medium
        ]
    ];

Column[{
    DescribeOutput[
        "Chain topology \[LongDash] no switching costs",
        "Green = entry, red = exits, gray = internal.",
        ScenarioTopologyPlot[chain2ExNoSC, sysNoSC, "Chain 1\[Rule]2\[Rule]3 (no SC)"]
    ],
    DescribeOutput[
        "Chain topology \[LongDash] with switching cost {1,2,3}=2.0",
        "Same topology; SC at vertex 2 penalises continuing from edge 1\[Rule]2 to 2\[Rule]3.",
        ScenarioTopologyPlot[chain2ExWithSC, sysWithSC, "Chain 1\[Rule]2\[Rule]3 (SC at 2)"]
    ]
}]


(* --- 8. Solve with ReduceSystem --- *)

(* Chain 1->2->3, single exit at 3 (cost=0), entry flow=10.
   All variables are uniquely determined. *)
chain1Ex = GridScenario[{3}, {{1, 10}}, {{3, 0}}];
sys1Ex   = makeSystem[chain1Ex];

DescribeOutput[
    "ReduceSystem \[LongDash] chain with one exit",
    "Unique solution: all j and u values pinned by flow balance + HJ + complementarity.",
    ReduceSystem[sys1Ex]
]


(* Chain 1->2->3, two exits at {2,0} and {3,10}, entry flow=120.
   Flow split between exits is under-determined; Reduce returns a parametric solution. *)
DescribeOutput[
    "ReduceSystem \[LongDash] chain with two exits (no switching costs)",
    "Parametric solution: one free variable governs how flow splits between exits.",
    ReduceSystem[sysNoSC]
]


