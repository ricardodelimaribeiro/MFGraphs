(* ::Package:: *)

(* Notebook-friendly MFGraphs workbook.
   This file is meant to be evaluated section-by-section inside a notebook.
   It does not assume that the notebook directory is already on $Path. *)

ClearAll[
    FindMFGraphsRoot,
    ClearMFGraphsNotebookShadows,
    EnsureMFGraphsLoaded,
    DescribeOutput,
    AssociationValue,
    NetEdgeFlows,
    NetworkVisualData,
    FlowStyleDirective,
    BaseNetworkGraphic,
    SolutionNetworkGraphic,
    DensityNetworkGraphic,
    MassDensityCurve,
    ValueFunctionCurve,
    ExitFlowChart,
    JamaratScenario
];

FindMFGraphsRoot[] :=
    Module[{origins, candidates},
        origins = DeleteDuplicates @ Select[
            {
                Quiet @ Check[NotebookDirectory[], Nothing],
                If[StringQ[$InputFileName] && $InputFileName =!= "",
                    DirectoryName[$InputFileName],
                    Nothing
                ],
                Directory[]
            },
            StringQ
        ];
        candidates = DeleteDuplicates @ Flatten[{
            origins,
            FileNameJoin[{#, "MFGraphs"}] & /@ origins,
            FileNameJoin[{#, "..", "MFGraphs"}] & /@ origins
        }];
        SelectFirst[
            candidates,
            FileExistsQ[FileNameJoin[{#, "MFGraphs.wl"}]] &,
            $Failed
        ]
    ];

ClearMFGraphsNotebookShadows[] :=
    Module[{shadowNames, shadowed},
        shadowNames = {
            "GetExampleData",
            "DataToEquations",
            "CriticalCongestionSolver",
            "NonLinearSolver",
            "MonotoneSolverFromData",
            "IsSwitchingCostConsistent",
            "PlotMassDensity",
            "PlotValueFunction",
            "I1", "I2", "I3",
            "U1", "U2", "U3",
            "S1", "S2", "S3", "S4", "S5", "S6",
            "V", "alpha", "g"
        };
        shadowed = Select[shadowNames, NameQ["Global`" <> #] &];
        Scan[ToExpression["Global`" <> #, InputForm, Remove] &, shadowed];
        shadowed
    ];

EnsureMFGraphsLoaded[] :=
    Module[{packageDir = FindMFGraphsRoot[], searchRoot, cleared},
        If[packageDir === $Failed,
            Print["Could not locate MFGraphs.wl relative to the notebook or current directory."];
            Abort[]
        ];
        searchRoot = DirectoryName[packageDir];
        cleared = ClearMFGraphsNotebookShadows[];
        If[FreeQ[$Path, searchRoot], PrependTo[$Path, searchRoot]];
        Get[FileNameJoin[{packageDir, "MFGraphs.wl"}]];
        If[cleared =!= {},
            Print["Removed shadowing Global` symbols before loading MFGraphs: ", cleared];
        ];
        packageDir
    ];

mfGraphsRoot = EnsureMFGraphsLoaded[];
MFGraphs`$MFGraphsVerbose = False;

(* Public MFGraphs symbols are now on $ContextPath, so the examples below can
   use GetExampleData, DataToEquations, I1, U1, ... without a context prefix.
   V, alpha, and g are now public MFGraphs symbols as well, so you can
   override them directly with Block[...] or with WithHamiltonianFunctions[...]. *)

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

AssociationValue[assoc_Association, key_, default_: Missing["NotAvailable"]] :=
    If[KeyExistsQ[assoc, key], assoc[key], default];

(* Net flow on each displayed edge, after eliminating dependent current variables. *)
NetEdgeFlows[d2e_Association, solution_Association, pairs_: Automatic] :=
    Module[{selectedPairs, rules},
        selectedPairs = Replace[
            pairs,
            Automatic -> (List @@@ Lookup[d2e, "edgeList", {}])
        ];
        rules = Join[
            Lookup[d2e, "RuleBalanceGatheringFlows", <||>],
            Lookup[d2e, "RuleExitFlowsIn", <||>],
            Lookup[d2e, "RuleEntryOut", <||>]
        ];
        AssociationThread[
            selectedPairs,
            N[((Lookup[d2e, "SignedFlows", <||>] /@ selectedPairs) /. rules /. solution)]
        ]
    ];

vertexPropertyMap[entryVertices_, exitVertices_, internalVertices_, entryVal_, exitVal_, internalVal_] :=
    Join[
        AssociationThread[entryVertices, ConstantArray[entryVal, Length[entryVertices]]],
        AssociationThread[exitVertices, ConstantArray[exitVal, Length[exitVertices]]],
        AssociationThread[internalVertices, ConstantArray[internalVal, Length[internalVertices]]]
    ];

NetworkVisualData[d2e_Association] :=
    Module[{graph, layoutGraph, coords, coordArray, extent,
      entryVertices, exitVertices, internalVertices},
        graph = Lookup[d2e, "graph", Graph[{}]];
        layoutGraph = Graph[graph, GraphLayout -> "SpringElectricalEmbedding"];
        coords = AssociationThread[VertexList[layoutGraph], N @ GraphEmbedding[layoutGraph]];
        coordArray = Values[coords];
        extent = If[coordArray === {},
            1.,
            Max[(Max /@ Transpose[coordArray]) - (Min /@ Transpose[coordArray])]
        ];
        extent = Max[extent, 1.];
        entryVertices = Lookup[d2e, "entryVertices", {}];
        exitVertices = Lookup[d2e, "exitVertices", {}];
        internalVertices = Complement[VertexList[graph], entryVertices, exitVertices];
        <|
            "graph" -> graph,
            "coords" -> coords,
            "vertexColors" -> vertexPropertyMap[entryVertices, exitVertices, internalVertices,
                RGBColor[0.22, 0.6, 0.3], RGBColor[0.82, 0.27, 0.2], GrayLevel[0.75]],
            "vertexSizes" -> vertexPropertyMap[entryVertices, exitVertices, internalVertices,
                0.34, 0.34, 0.28],
            "vertexRadii" -> vertexPropertyMap[entryVertices, exitVertices, internalVertices,
                0.045 extent, 0.045 extent, 0.035 extent]
        |>
    ];

FlowStyleDirective[flow_?NumericQ, maxFlow_?NumericQ] :=
    Directive[
        If[flow >= 0, RGBColor[0.12, 0.45, 0.78], RGBColor[0.82, 0.39, 0.2]],
        AbsoluteThickness[
            Rescale[Abs[flow], {0, Max[maxFlow, 10^-9]}, {1.5, 8}]
        ],
        Opacity[0.9]
    ];

BaseNetworkGraphic[d2e_Association, title_: Automatic] :=
    Module[{visual, plotTitle},
        visual = NetworkVisualData[d2e];
        plotTitle = Replace[title, Automatic -> "Network structure"];
        Graph[
            visual["graph"],
            GraphLayout -> "SpringElectricalEmbedding",
            VertexLabels -> Placed["Name", Center],
            VertexStyle -> Normal[visual["vertexColors"]],
            VertexSize -> Normal[visual["vertexSizes"]],
            EdgeStyle -> Directive[GrayLevel[0.6], AbsoluteThickness[2]],
            PlotLabel -> Style[plotTitle, 14, Bold],
            ImageSize -> Large
        ]
    ];

SolutionNetworkGraphic[d2e_Association, solution_Association, title_: Automatic] :=
    Module[{visual, pairs, graphEdges, flows, flowValues, maxFlow, edgeStyles,
      edgeLabels, plotTitle},
        visual = NetworkVisualData[d2e];
        graphEdges = Lookup[d2e, "edgeList", {}];
        pairs = List @@@ graphEdges;
        flows = NetEdgeFlows[d2e, solution, pairs];
        flowValues = N[Lookup[flows, pairs]];
        maxFlow = Max[Append[Abs[flowValues], 0]];
        edgeStyles = AssociationThread[
            graphEdges,
            FlowStyleDirective[#, maxFlow] & /@ flowValues
        ];
        edgeLabels = AssociationThread[
            graphEdges,
            Placed[
                Style[
                    NumberForm[Chop[#, 10^-8], {Infinity, 1}],
                    11,
                    Black,
                    Background -> White
                ],
                Center
            ] & /@ flowValues
        ];
        plotTitle = Replace[title, Automatic -> "Net edge flows"];
        Graph[
            visual["graph"],
            GraphLayout -> "SpringElectricalEmbedding",
            VertexLabels -> Placed["Name", Center],
            VertexStyle -> Normal[visual["vertexColors"]],
            VertexSize -> Normal[visual["vertexSizes"]],
            EdgeStyle -> Normal[edgeStyles],
            EdgeLabels -> Normal[edgeLabels],
            PlotLabel -> Style[plotTitle, 14, Bold],
            ImageSize -> Large
        ]
    ];

DensityNetworkGraphic[d2e_Association, solution_Association, title_: Automatic, samples_Integer: 24] :=
    Module[{visual, coords, pairs, flows, sampleXs, densitiesByPair, allDensities,
      maxDensity, colorScale, edgeSegments, edgeLabels, vertexDisks, vertexLabels,
      plotTitle, legendRange},
        visual = NetworkVisualData[d2e];
        coords = visual["coords"];
        pairs = List @@@ Lookup[d2e, "edgeList", {}];
        flows = NetEdgeFlows[d2e, solution, pairs];
        sampleXs = N @ Subdivide[0., 1., samples];
        densitiesByPair = AssociationThread[
            pairs,
            Table[
                Quiet @ Check[
                    MFGraphs`Private`M[AssociationValue[flows, pair, 0.], x, UndirectedEdge @@ pair],
                    0.
                ],
                {pair, pairs}, {x, sampleXs}
            ]
        ];
        allDensities = Select[Flatten[Values[densitiesByPair]], NumericQ];
        maxDensity = Max[Append[allDensities, 0.]];
        legendRange = {0., Max[maxDensity, 1.]};
        colorScale[value_] :=
            ColorData["SolarColors"][Rescale[value, legendRange]];
        edgeSegments = Flatten @ Map[
            Function[pair,
                Module[{edgePoints, densities},
                    edgePoints = ((1 - #) coords[pair[[1]]] + # coords[pair[[2]]]) & /@ sampleXs;
                    densities = AssociationValue[densitiesByPair, pair, ConstantArray[0., Length[sampleXs]]];
                    Table[
                        {
                            colorScale[Mean[densities[[k ;; k + 1]]]],
                            AbsoluteThickness[8],
                            Line[edgePoints[[k ;; k + 1]]]
                        },
                        {k, 1, Length[edgePoints] - 1}
                    ]
                ]
            ],
            pairs
        ];
        edgeLabels = Map[
            Function[pair,
                Inset[
                    Framed[
                        Style[
                            Row[{"flow = ", NumberForm[AssociationValue[flows, pair, 0.], {Infinity, 1}]}],
                            11,
                            Black
                        ],
                        Background -> White,
                        FrameStyle -> GrayLevel[0.8],
                        RoundingRadius -> 3,
                        FrameMargins -> Tiny
                    ],
                    Mean[{coords[pair[[1]]], coords[pair[[2]]]}]
                ]
            ],
            pairs
        ];
        vertexDisks = Map[
            Function[vertex,
                {
                    Lookup[visual["vertexColors"], vertex, GrayLevel[0.75]],
                    EdgeForm[Directive[White, AbsoluteThickness[1.5]]],
                    Disk[
                        coords[vertex],
                        Lookup[visual["vertexRadii"], vertex, 0.04]
                    ]
                }
            ],
            VertexList[visual["graph"]]
        ];
        vertexLabels = Map[
            Function[vertex,
                Inset[
                    Style[vertex, 11, Bold, Black],
                    coords[vertex]
                ]
            ],
            VertexList[visual["graph"]]
        ];
        plotTitle = Replace[title, Automatic -> "Edge density profile"];
        Legended[
            Graphics[
                {
                    edgeSegments,
                    edgeLabels,
                    vertexDisks,
                    vertexLabels
                },
                PlotRange -> All,
                PlotRangePadding -> Scaled[0.15],
                ImageSize -> Large,
                PlotLabel -> Style[plotTitle, 14, Bold]
            ],
            Placed[
                BarLegend[
                    {ColorData["SolarColors"], legendRange},
                    LegendLabel -> Style["Local density", 11]
                ],
                Right
            ]
        ]
    ];

MassDensityCurve[result_Association, pair_List] :=
    Show[
        MFGraphs`Private`PlotMassDensity[result, "AssoNonCritical", pair],
        PlotLabel -> Style[
            Row[{"Mass density along edge ", pair}],
            14,
            Bold
        ],
        AxesLabel -> {"Position on edge", "Density"}
    ];

ValueFunctionCurve[result_Association, pair_List] :=
    Show[
        MFGraphs`Private`PlotValueFunction[result, "AssoNonCritical", pair],
        PlotLabel -> Style[
            Row[{"Value function along edge ", pair}],
            14,
            Bold
        ],
        AxesLabel -> {"Position on edge", "Value"}
    ];

ExitFlowChart[exitFlows_Association, title_: Automatic] :=
    Module[{vertices, values, plotTitle},
        vertices = Keys[exitFlows];
        values = N[Values[exitFlows]];
        plotTitle = Replace[title, Automatic -> "Exit flow totals"];
        BarChart[
            values,
            ChartLabels -> Placed[vertices, Below],
            ChartStyle -> Table[ColorData[97][k], {k, Length[values]}],
            AxesLabel -> {"Exit", "Flow"},
            PlotLabel -> Style[plotTitle, 14, Bold],
            GridLines -> Automatic,
            ImageSize -> Large
        ]
    ];

(* Jamarat helper: solve one release/cost scenario and summarize exit usage. *)
JamaratScenario[params_List] :=
    Module[{data, d2e, result, exitPairs, exitFlows},
        data = GetExampleData["Jamaratv9"] /. params;
        d2e = DataToEquations[data];
        result = Quiet @ CriticalCongestionSolver[d2e];
        exitPairs = List @@@ Lookup[d2e, "exitEdges", {}];
        exitFlows = If[AssociationQ[result["Solution"]],
            AssociationThread[
                Lookup[d2e, "exitVertices", {}],
                N[((Lookup[d2e, "SignedFlows", <||>] /@ exitPairs) /. result["Solution"])]
            ],
            Missing["NotAvailable"]
        ];
        <|
            "Parameters" -> params,
            "Model" -> d2e,
            "Solution" -> result["Solution"],
            "Feasibility" -> result["Feasibility"],
            "ExitFlows" -> exitFlows,
            "InteriorFlows" ->
                If[AssociationQ[result["Solution"]],
                    NetEdgeFlows[d2e, result["Solution"]],
                    Missing["NotAvailable"]
                ]
        |>
    ];


(* --- 1. Explore built-in data --- *)

examplePreview = <|
    "Linear path" -> GetExampleData[3],
    "Y network" -> GetExampleData[7],
    "Attraction problem" -> GetExampleData[12],
    "Jamarat" -> GetExampleData["Jamaratv9"]
|>;

examplePreview["Jamarat"]


(* --- 2. Define your own network --- *)

customData = <|
    "Vertices List" -> {1, 2, 3, 4},
    "Adjacency Matrix" -> {
        {0, 1, 1, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 1},
        {0, 0, 0, 0}
    },
    "Entrance Vertices and Flows" -> {{1, 120}},
    "Exit Vertices and Terminal Costs" -> {{4, 0}},
    "Switching Costs" -> {}
|>;

customD2E = DataToEquations[customData];
Column[{
    DescribeOutput[
        "Custom network edge list",
        "The undirected edges below are the interior corridors of the model before any solve is run.",
        customD2E["edgeList"]
    ],
    DescribeOutput[
        "Custom network topology",
        "Green vertices are entrances, red vertices are exits, and gray vertices are interior junctions.",
        BaseNetworkGraphic[customD2E, "Custom network topology"]
    ]
}]


(* --- 3. Critical congestion solver --- *)

criticalData = GetExampleData[7] /. {
    I1 -> 100,
    U1 -> 0,
    U2 -> 0
};
criticalD2E = DataToEquations[criticalData];

criticalLegacy = CriticalCongestionSolver[criticalD2E];
criticalLegacy["AssoCritical"];

criticalStandard = CriticalCongestionSolver[criticalD2E];
Column[{
    DescribeOutput[
        "Critical solver solution",
        "This association stores the equilibrium directional currents, switching currents, and edge values for the critical-congestion model.",
        criticalStandard["Solution"]
    ],
    DescribeOutput[
        "Net edge flows",
        "Edge thickness and labels show the signed net flow on each interior edge. Positive values follow the stored edge orientation.",
        SolutionNetworkGraphic[
            criticalD2E,
            criticalStandard["Solution"],
            "Critical-congestion net flows"
        ]
    ],
    DescribeOutput[
        "Density along each edge",
        "The color scale shows local agent density along each edge under a simple baseline Hamiltonian with zero potential and the default congestion law.",
        WithHamiltonianFunctions[
            Function[{x, edge}, 0],
            Function[edge, 1],
            Function[{m, edge}, -1 / m^2],
            DensityNetworkGraphic[
                criticalD2E,
                criticalStandard["Solution"],
                "Critical-congestion edge densities"
            ]
        ]
    ],
    DescribeOutput[
        "Net edge flow values",
        "This association is often the easiest summary to inspect numerically when checking route splits.",
        NetEdgeFlows[criticalD2E, criticalStandard["Solution"]]
    ]
}]


(* --- 4. Switching-cost validation --- *)

switchingData = GetExampleData[8] /. {
    I1 -> 100, U1 -> 0, U2 -> 0,
    S1 -> 2, S2 -> 3, S3 -> 2,
    S4 -> 1, S5 -> 3, S6 -> 1
};
switchingD2E = DataToEquations[switchingData];
DescribeOutput[
    "Switching-cost consistency check",
    "A value of True means the local switching penalties are compatible with a globally consistent potential at each junction.",
    IsSwitchingCostConsistent[Normal @ switchingD2E["SwitchingCosts"]]
]


(* --- 5. Non-linear congestion solver and plots --- *)

Module[{potentialFunction, congestionExponentFunction, interactionFunction},
    potentialFunction = Function[{x, edge},
        Which[
            edge === UndirectedEdge[1, 2], 0.25 x,
            edge === UndirectedEdge[2, 4], 0.5,
            True, 0
        ]
    ];
    congestionExponentFunction = Function[edge,
        If[edge === UndirectedEdge[2, 4], 1.3, 1.0]
    ];
    interactionFunction = Function[{m, edge},
        -If[edge === UndirectedEdge[2, 4], 1.5, 1.0] / m^2
    ];
    nonlinearData = GetExampleData[12] /. {
        I1 -> 2,
        U1 -> 0
    };
    nonlinearD2E = DataToEquations[nonlinearData];
    nonlinearResult = NonLinearSolver[
        nonlinearD2E,
        "MaxIterations" -> 5,
        "Tolerance" -> 10^-8,
        "PotentialFunction" -> potentialFunction,
        "CongestionExponentFunction" -> congestionExponentFunction,
        "InteractionFunction" -> interactionFunction
    ];
    Column[{
        DescribeOutput[
            "Non-linear solver solution",
            "This fixed-point solve includes the edge Hamiltonian and congestion interaction, so densities and values respond to the chosen V, alpha, and g.",
            nonlinearResult["Solution"]
        ],
        DescribeOutput[
            "Non-linear equilibrium net flows",
            "Thickness and labels summarize the net current on each edge after the non-linear congestion update converges.",
            SolutionNetworkGraphic[
                nonlinearD2E,
                nonlinearResult["Solution"],
                "Non-linear equilibrium net flows"
            ]
        ],
        DescribeOutput[
            "Non-linear edge densities",
            "Warm colors mark locally denser parts of the network. The legend is pointwise density, not total edge flow.",
            WithHamiltonianFunctions[
                potentialFunction,
                congestionExponentFunction,
                interactionFunction,
                DensityNetworkGraphic[
                    nonlinearD2E,
                    nonlinearResult["Solution"],
                    "Non-linear density map"
                ]
            ]
        ],
        DescribeOutput[
            "Mass density on edge {1, 2}",
            "This one-dimensional slice shows how density changes with normalized position along a selected edge.",
            WithHamiltonianFunctions[
                potentialFunction,
                congestionExponentFunction,
                interactionFunction,
                MassDensityCurve[nonlinearResult, {1, 2}]
            ]
        ],
        DescribeOutput[
            "Value function on edge {1, 2}",
            "This plot shows the value perceived by agents as they move along the same edge under the non-linear equilibrium.",
            WithHamiltonianFunctions[
                potentialFunction,
                congestionExponentFunction,
                interactionFunction,
                ValueFunctionCurve[nonlinearResult, {1, 2}]
            ]
        ]
    }]
]


(* --- 6. Monotone solver --- *)

Module[{potentialFunction, congestionExponentFunction, interactionFunction},
    potentialFunction = Function[{x, edge}, 0];
    congestionExponentFunction = Function[edge, 1];
    interactionFunction = Function[{m, edge}, -1 / m^2];
    monotoneData = GetExampleData[3] /. {
        I1 -> 80,
        U1 -> 0
    };
    monotoneResult = MonotoneSolverFromData[
        monotoneData,
        "ResidualTolerance" -> 10^-6,
        "MaxTime" -> 10,
        "MaxSteps" -> 2000,
        "PotentialFunction" -> potentialFunction,
        "CongestionExponentFunction" -> congestionExponentFunction,
        "InteractionFunction" -> interactionFunction
    ];
    Column[{
        DescribeOutput[
            "Monotone solver solution",
            "The Hessian-Riemannian flow evolves the Kirchhoff variables until it reaches a monotone equilibrium.",
            monotoneResult["Solution"]
        ],
        DescribeOutput[
            "Monotone solver net flows",
            "This graph shows the final net currents produced by the monotone solver on the simple path example.",
            SolutionNetworkGraphic[
                DataToEquations[monotoneData],
                monotoneResult["Solution"],
                "Monotone solver net flows"
            ]
        ],
        DescribeOutput[
            "Monotone edge densities",
            "With the same baseline Hamiltonian, the edge color gradient shows how density varies along the path under the monotone solution.",
            WithHamiltonianFunctions[
                potentialFunction,
                congestionExponentFunction,
                interactionFunction,
                DensityNetworkGraphic[
                    DataToEquations[monotoneData],
                    monotoneResult["Solution"],
                    "Monotone density map"
                ]
            ]
        ]
    }]
]


(* --- 7. Jamarat / pilgrimage crowd routing --- *)

(* The built-in Jamarat network has two entrances and three exits. *)
jamaratTemplate = GetExampleData["Jamaratv9"];
jamaratTemplate["Vertices List"];
jamaratTemplate["Exit Vertices and Terminal Costs"];

(* Feasible release plan from the existing regression tests. *)
jamaratBase = JamaratScenario[
    {
        I1 -> 2,
        I2 -> 2000,
        U1 -> 20,
        U2 -> 100,
        U3 -> 0
    }
];
Column[{
    DescribeOutput[
        "Jamarat feasibility status",
        "This baseline release plan is the first check: it tells you whether the symbolic critical-congestion model finds a feasible equilibrium.",
        jamaratBase["Feasibility"]
    ],
    DescribeOutput[
        "Jamarat baseline net flows",
        "The network view shows where the pilgrims are routed through the interior corridors under the baseline release and exit-cost assumptions.",
        SolutionNetworkGraphic[
            jamaratBase["Model"],
            jamaratBase["Solution"],
            "Jamarat baseline net flows"
        ]
    ],
    DescribeOutput[
        "Jamarat baseline densities",
        "The edge color gradient highlights where the baseline routing would concentrate agents most strongly under the default Hamiltonian model.",
        WithHamiltonianFunctions[
            Function[{x, edge}, 0],
            Function[edge, 1],
            Function[{m, edge}, -1 / m^2],
            DensityNetworkGraphic[
                jamaratBase["Model"],
                jamaratBase["Solution"],
                "Jamarat baseline density map"
            ]
        ]
    ],
    DescribeOutput[
        "Jamarat exit split",
        "Bar heights summarize the total predicted outflow through each exit in the feasible baseline scenario.",
        ExitFlowChart[jamaratBase["ExitFlows"], "Jamarat baseline exit split"]
    ],
    DescribeOutput[
        "Jamarat exit totals",
        "This association is the numeric form of the exit split shown in the chart above.",
        jamaratBase["ExitFlows"]
    ]
}]


(* Overloaded plan: the symbolic critical-congestion solver detects infeasibility. *)
jamaratOverload = JamaratScenario[
    {
        I1 -> 130,
        I2 -> 128,
        U1 -> 20,
        U2 -> 100,
        U3 -> 0
    }
];
Column[{
    DescribeOutput[
        "Overloaded Jamarat scenario",
        "This release plan is intentionally infeasible. Use it as a comparison point when stress-testing gate releases or exit penalties.",
        jamaratOverload["Feasibility"]
    ],
    DescribeOutput[
        "Jamarat network structure",
        "When a scenario is infeasible, it is still useful to inspect the topology alone before changing releases or costs.",
        BaseNetworkGraphic[jamaratOverload["Model"], "Jamarat network structure"]
    ]
}]

(* Optional extension for later:
   To compare guidance policies, create a small list of cost scenarios and map
   JamaratScenario over them. Start from the tested baseline above, then change
   U1, U2, and U3 one at a time and inspect ExitFlows. Keep the first sweep
   conservative so you can tell whether a surprising split comes from your
   policy choice or from a numerically delicate symbolic branch.

   Example template:
   jamaratScenarios = Association[
       "Baseline" -> JamaratScenario[{I1 -> 2, I2 -> 2000, U1 -> 20, U2 -> 100, U3 -> 0}],
       "Discourage Exit 8" -> JamaratScenario[{I1 -> 2, I2 -> 2000, U1 -> 20, U2 -> 130, U3 -> 0}]
   ];

   KeyValueMap[
       <|
           "Scenario" -> #1,
           "Feasibility" -> #2["Feasibility"],
           "Exit7" -> Lookup[#2["ExitFlows"], 7, Missing["NotAvailable"]],
           "Exit8" -> Lookup[#2["ExitFlows"], 8, Missing["NotAvailable"]],
           "Exit9" -> Lookup[#2["ExitFlows"], 9, Missing["NotAvailable"]]
       |>&,
       jamaratScenarios
   ];

   Once you have calibrated edge classes for wide corridors, ramps, and merges,
   reuse the same pattern from the non-linear section above on the Jamarat graph.
   In practice, start with the critical-congestion symbolic solve to screen release
   plans and destination costs, then move to NonLinearSolver only after the graph
   geometry and demand windows are calibrated. *)
