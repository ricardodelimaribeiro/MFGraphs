(* ::Package:: *)

(* Notebook-friendly MFGraphs workbook focused on the typed scenario kernels.
   This version intentionally excludes solver examples. *)

(* Evaluate cells one at a time or section by section — do not evaluate the entire file at once. *)

(* This file lives alongside MFGraphs.wl in the same directory. *)
mfgDir = If[$InputFileName === "", NotebookDirectory[], DirectoryName[$InputFileName]];

(* Force clean reload — safe to re-evaluate without restarting the kernel. *)
Quiet[
    Unprotect /@ Names["MFGraphs`*"];
    Remove /@ Names["MFGraphs`*"];
    Unprotect[$Packages];
    $Packages = DeleteCases[$Packages, "MFGraphs`"];
    Protect[$Packages]
];
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
    "# Equ" -> Length @ SystemData[exampleSystem, "Equ"],
    "# EqGeneral" -> Length @ SystemData[exampleSystem, "EqGeneral"],
    "# EqTransitionFlow" -> Length @ SystemData[exampleSystem, "EqTransitionFlow"],
    "# EqJunctionPotentials" -> Length @ SystemData[exampleSystem, "EqJunctionPotentials"],
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

(* Optional: inspect raw equation blocks directly if needed. *)
Take[SystemData[exampleSystem, "EqGeneral"], UpTo[10]]
Take[SystemData[exampleSystem, "EqTransitionFlow"], UpTo[10]]
Take[SystemData[exampleSystem, "EqJunctionPotentials"], UpTo[10]]
