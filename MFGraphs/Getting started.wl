(* ::Package:: *)

(*Quit[]*)


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
    EdgeEndpoints,
    EdgeHamiltonianValue,
    FormatHamiltonianTerm,
    EdgeModelSummary,
    EdgeModelPlot
];

EdgeEndpoints[UndirectedEdge[a_, b_]] := {a, b};
EdgeEndpoints[DirectedEdge[a_, b_]] := {a, b};
EdgeEndpoints[{a_, b_}] := {a, b};

EdgeHamiltonianValue[edgeValues_Association, edge_List, default_] :=
    Lookup[
        edgeValues,
        Key[edge],
        Lookup[edgeValues, Key[Reverse[edge]], default]
    ];

FormatHamiltonianTerm[term_] :=
    Which[
        MatchQ[term, _Function], ToString[Unevaluated[term], InputForm],
        True, term
    ];

EdgeModelSummary[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{hamiltonian, edges, edgeModel},
        hamiltonian = scenarioData[s, "Hamiltonian"];
        edges = EdgeEndpoints /@ systemData[sys, "Edges"];
        edgeModel[edge_] := <|
            "Alpha" -> EdgeHamiltonianValue[
                Lookup[hamiltonian, "EdgeAlpha", <||>],
                edge,
                Lookup[hamiltonian, "Alpha", 1]
            ],
            "V" -> EdgeHamiltonianValue[
                Lookup[hamiltonian, "EdgeV", <||>],
                edge,
                Lookup[hamiltonian, "V", 0]
            ],
            "G" -> FormatHamiltonianTerm @ EdgeHamiltonianValue[
                Lookup[hamiltonian, "EdgeG", <||>],
                edge,
                Lookup[hamiltonian, "G", Function[z, -1/z]]
            ]
        |>;
        AssociationThread[edges, edgeModel /@ edges]
    ];

Options[EdgeModelPlot] = {
    GraphLayout -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

EdgeModelPlot[s_?scenarioQ, sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    Module[{summary, modelEdges, edgeLabels},
        summary = EdgeModelSummary[s, sys];
        modelEdges = systemData[sys, "Edges"];
        edgeLabels = Association @ Map[
            With[
                {
                    edge = #,
                    model = Lookup[summary, Key[EdgeEndpoints[#]], <||>]
                },
                edge -> Placed[
                    Style[
                        Row[{
                            "\[Alpha]=", Lookup[model, "Alpha", "?"],
                            ", V=", Lookup[model, "V", "?"],
                            ", G=", Lookup[model, "G", "?"]
                        }],
                        9,
                        Black,
                        Background -> White
                    ],
                    Center
                ]
            ] &,
            modelEdges
        ];
        Graph[
            scenarioTopologyPlot[s, sys, PlotLabel -> None,
                GraphLayout -> OptionValue[GraphLayout],
                ImageSize -> OptionValue[ImageSize]],
            EdgeLabels -> Normal[edgeLabels],
            PlotLabel -> Replace[
                OptionValue[PlotLabel],
                {
                    Automatic -> Style["Per-edge Hamiltonian model", 14, Bold],
                    None -> None,
                    other_ :> Style[other, 14, Bold]
                }
            ],
            ImageSize -> OptionValue[ImageSize]
        ]
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

exY = gridScenario[{3}, {{2, 100.0}}, {{1, 0.0}, {3, 10.0}}];

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


(* --- 4. Show the Hamiltonian model stored on each edge --- *)

edgeModelScenario = makeScenario[
    <|
        "Model" -> <|
            "Vertices" -> {1, 2, 3, 4},
            "Adjacency" -> {
                {0, 1, 1, 0},
                {0, 0, 0, 1},
                {0, 0, 0, 1},
                {0, 0, 0, 0}
            },
            "Entries" -> {{1, 80.0}},
            "Exits" -> {{4, 0.0}},
            "Switching" -> {}
        |>,
        "Hamiltonian" -> <|
            "Alpha" -> 1,
            "V" -> 0,
            "G" -> Function[z, -1/z],
            "EdgeAlpha" -> <|{1, 2} -> 1.0, {1, 3} -> 2.0, {2, 4} -> 1.5|>,
            "EdgeV" -> <|{1, 3} -> 0.25, {3, 4} -> -0.5|>,
            "EdgeG" -> <|{2, 4} -> Function[z, -2/z], {3, 4} -> Function[z, -1/(2 z)]|>
        |>
    |>
];

edgeModelSystem = makeSystem[edgeModelScenario];

Column[{
    DescribeOutput[
        "Topology plot from graphicsTools",
        "scenarioTopologyPlot shows the network while preserving entry/exit coloring.",
        scenarioTopologyPlot[edgeModelScenario, edgeModelSystem,
            PlotLabel -> "Example network for per-edge model parameters",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Per-edge Hamiltonian model",
        "Each edge shows the effective \[Alpha], V, and G after applying per-edge overrides. Current system construction applies \[Alpha]/EdgeAlpha; V/G are stored for visualization and future density work.",
        EdgeModelPlot[edgeModelScenario, edgeModelSystem,
            PlotLabel -> "Hamiltonian model on each edge",
            ImageSize -> Large]
    ],
    DescribeOutput[
        "Per-edge model table",
        "Association keyed by network edge {u,v}; values are the effective model terms used for display.",
        EdgeModelSummary[edgeModelScenario, edgeModelSystem]
    ],
    DescribeOutput[
        "Hamilton-Jacobi equations with EdgeAlpha",
        "The structural equations reflect Alpha/EdgeAlpha on each edge.",
        systemData[edgeModelSystem, "EqGeneral"]
    ]
}]


(* --- 5. Build exact symbolic unknown bundles from unknownsTools.wl --- *)

exampleUnknowns = makeSymbolicUnknowns[exY];

unknownSummary = <|
    "symbolicUnknownsQ" -> symbolicUnknownsQ[exampleUnknowns],
    "js count" -> Length @ symbolicUnknownsData[exampleUnknowns, "Js"],
    "jts count" -> Length @ symbolicUnknownsData[exampleUnknowns, "Jts"],
    "us count" -> Length @ symbolicUnknownsData[exampleUnknowns, "Us"],
    "first js" -> Take[symbolicUnknownsData[exampleUnknowns, "Js"], UpTo[6]],
    "first jts" -> Take[symbolicUnknownsData[exampleUnknowns, "Jts"], UpTo[6]],
    "first us" -> Take[symbolicUnknownsData[exampleUnknowns, "Us"], UpTo[6]]
|>;

DescribeOutput[
    "Symbolic unknown bundle",
    "makeSymbolicUnknowns derives flow, transition-flow, and value-function symbolic variables from scenario topology.",
    unknownSummary
]


(* --- 6. Build structural systems from System.wl --- *)

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
    "Switching-cost consistency" -> systemData[exampleSystem, "ConsistentCosts"]
|>;

DescribeOutput[
    "mfgSystem overview",
    "makeSystem builds structural equations and inequalities without running any solver.",
    systemSummary
]


(* --- 7. Chain with two exits: equations without and with switching costs --- *)

chain2ExNoSC = gridScenario[
    {3},
    {{1, 6}},
    {{2, 10}, {3, 0}}
];
chain2ExWithSC = gridScenario[
    {3},
    {{1, 120}},
    {{2, 10}, {3, 0}},
    {{1, 2, 3, 2}}
];

unkNoSC   = makeSymbolicUnknowns[chain2ExNoSC];
sysNoSC   = makeSystem[chain2ExNoSC, unkNoSC];
unkWithSC = makeSymbolicUnknowns[chain2ExWithSC];
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


(* --- 8. Topology visualization (from MFGraphs`graphicsTools`) --- *)

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


(* --- 9. Solve with solveScenario (DNF-first default) --- *)

(* Chain 1->2->3, single exit at 3 (cost=0), entry flow=10.
   All variables are uniquely determined. *)
chain1Ex = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
sys1Ex   = makeSystem[chain1Ex];

DescribeOutput[
    "solveScenario \[LongDash] chain with one exit",
    "Unique solution: all j and u values pinned by flow balance + HJ + complementarity.",
    solveScenario[chain1Ex]
]


(* Chain 1->2->3, two exits at {2,0} and {3,10}, entry flow=120.
   Flow split between exits is under-determined; the solver returns a parametric solution. *)
DescribeOutput[
    "solveScenario \[LongDash] chain with two exits (no switching costs)",
    "Parametric solution: one free variable governs how flow splits between exits.",
    solveScenario[chain2ExNoSC]
]


(* --- 10. Visualize solveScenario solutions (from MFGraphs`graphicsTools`) --- *)

(* --- Apply to chain with one exit --- *)

sol1Ex = solveScenario[chain1Ex];

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

solNoSC = solveScenario[chain2ExNoSC];

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


(* --- 11. Advanced Solution Visualization (Paper Scheme) --- *)

augChain1 = augmentAuxiliaryGraph[sys1Ex];

Column[{
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
        "Transition graph \[LongDash] gradient node coloring",
        "mfgTransitionPlot now colors nodes on a Red\[Rule]Blue u-value gradient and draws edges as quadratic Bezier arcs. Anti-parallel arcs curve to opposite sides. A color bar legend is shown by default.",
        mfgTransitionPlot[chain1Ex, sys1Ex, sol1Ex,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: transition graph"]
    ],
    DescribeOutput[
        "Transition graph \[LongDash] BendFactor comparison",
        "BendFactor controls arc curvature as a fraction of edge length. BendFactor\[Rule]0 gives straight edges; higher values increase the arc.",
        Row[{
            mfgTransitionPlot[chain1Ex, sys1Ex, sol1Ex,
                PlotLabel -> "BendFactor \[Rule] 0 (straight)", ShowLegend -> False,
                BendFactor -> 0, ImageSize -> 300],
            mfgTransitionPlot[chain1Ex, sys1Ex, sol1Ex,
                PlotLabel -> "BendFactor \[Rule] 0.15 (default)", ShowLegend -> False,
                BendFactor -> 0.15, ImageSize -> 300],
            mfgTransitionPlot[chain1Ex, sys1Ex, sol1Ex,
                PlotLabel -> "BendFactor \[Rule] 0.35", ShowLegend -> False,
                BendFactor -> 0.35, ImageSize -> 300]
        }, Spacer[12]]
    ],
    DescribeOutput[
        "Augmented infrastructure graph \[LongDash] with solution",
        "Nodes are road-traffic edge-vertex states. Blue arcs are j[a,b] flows; red arcs are j[r,i,w] transitions. Anti-parallel flow arcs curve to opposite sides. Color gradient on nodes shows u-values.",
        mfgAugmentedPlot[chain1Ex, sys1Ex, sol1Ex,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: augmented infrastructure (solved)"]
    ],
    DescribeOutput[
        "Augmented infrastructure graph \[LongDash] pre-solve with ShowBoundaryValues\[Rule]False",
        "When ShowBoundaryValues\[Rule]False, boundary values are hidden but entry (green) and exit (red) nodes keep their category color. Internal nodes are gray until a solution is provided.",
        mfgAugmentedPlot[chain1Ex, sys1Ex, <||>,
            PlotLabel -> "Chain 1\[Rule]2\[Rule]3: structure only (ShowBoundaryValues\[Rule]False)",
            ShowBoundaryValues -> False]
    ]
}]


(* --- 12. Paper-style network: Figure 3 \[Rule] Figure 4 transformation --- *)

paperFig3Scenario = graphScenario[
    Graph[{1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3, 3 \[UndirectedEdge] 4}],
    {{1, 10}, {3, 5}},
    {{2, 0}, {4, 0}}
];

paperFig3System = makeSystem[paperFig3Scenario];
paperFig3Augmented = augmentAuxiliaryGraph[paperFig3System];

paperFig3AuxGraph = Graph[
    systemData[paperFig3System, "AuxVertices"],
    systemData[paperFig3System, "AuxEdges"],
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
            PlotLabel -> "Paper-style augmented road-traffic graph",
            ShowBoundaryValues -> False]
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


(* --- 13. Jamaratv9 example from the named scenario registry --- *)

jamaratScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 100}, {2, 50}},
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
        "Jamaratv9 augmented road-traffic graph (structure only)",
        "ShowBoundaryValues\[Rule]False hides u-values but entry (green) and exit (red) nodes keep their category color so boundary topology remains readable.",
        mfgAugmentedPlot[jamaratScenario, jamaratSystem, <||>,
            PlotLabel -> "Jamaratv9 augmented infrastructure (structure only)",
            ShowBoundaryValues -> False,
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
