(* ::Package:: *)

(* Notebook-friendly MFGraphs workbook covering the typed scenario kernels
   and the reduceSystem solver. *)

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

gridScenario = gridScenario[
    {3},
    {{1, 120.0}},
    {{3, 0.0}}
];

cycleScenario = cycleScenario[
    3,
    {{1, 50.0}},
    {{2, 0.0}, {3, 10.0}},
    {
        {1, 2, 3, 2.0},
        {3, 2, 1, 1.0}
    }
];

amScenario = amScenario[
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
        gridScenario
    ],
    DescribeOutput[
        "cycleScenario output",
        "Cycle topology with explicit switching costs.",
        cycleScenario
    ],
    DescribeOutput[
        "amScenario output",
        "Adjacency-matrix constructor for explicit benchmark topologies.",
        amScenario
    ]
}]


(* --- 2. Use named factory examples from ExampleScenarios.wl --- *)

exampleY = getExampleScenario[
    7,
    {{1, 100.0}},
    {{3, 0.0}, {4, 10.0}}
];

exampleGrid = getExampleScenario[
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
    "Identity" -> scenarioData[exampleY, "Identity"],
    "Benchmark" -> scenarioData[exampleY, "Benchmark"],
    "Hamiltonian" -> scenarioData[exampleY, "Hamiltonian"],
    "Model keys" -> Keys @ scenarioData[exampleY, "Model"],
    "Topology keys" -> Keys @ scenarioData[exampleY, "Topology"]
|>;

DescribeOutput[
    "scenarioData overview",
    "Typed scenarios expose canonical blocks used by the downstream kernels.",
    scenarioChecks
]


(* --- 4. Build unknown bundles from Unknowns.wl --- *)

exampleUnknowns = makeUnknowns[exampleY];

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

exampleSystem = makeSystem[exampleY, exampleUnknowns];

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
    "Switching-cost consistency" -> IsSwitchingCostConsistent[
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
        scenarioTopologyPlot[chain2ExNoSC, sysNoSC, "Chain 1\[Rule]2\[Rule]3 (no SC)"]
    ],
    DescribeOutput[
        "Chain topology \[LongDash] with switching cost {1,2,3}=2.0",
        "Same topology; SC at vertex 2 penalises continuing from edge 1\[Rule]2 to 2\[Rule]3.",
        scenarioTopologyPlot[chain2ExWithSC, sysWithSC, "Chain 1\[Rule]2\[Rule]3 (SC at 2)"]
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
        "Combined solution plot \[LongDash] chain 1\[Rule]2\[Rule]3, single exit",
        "Directed edges show j-flow direction/magnitude; labels show both j and u. Auxiliary edges are included.",
        mfgSolutionPlot[chain1Ex, sys1Ex, sol1Ex,
            "Chain 1\[Rule]2\[Rule]3: combined solution (exit at 3, cost=0)"]
    ]
}]


(* --- Apply to chain with two exits (partial rules from underdetermined solution) --- *)

solNoSC = reduceSystem[sysNoSC];

Column[{
    DescribeOutput[
        "Combined solution plot \[LongDash] chain with two exits (no SC)",
        "Edge {1,2} is directed 1\[Rule]2 (j[1,2]>0). Labels show j and u on real + auxiliary edges.",
        mfgSolutionPlot[chain2ExNoSC, sysNoSC, solNoSC,
            "Chain 1\[Rule]2\[Rule]3: combined solution (exits at 2 and 3)"]
    ]
}]
