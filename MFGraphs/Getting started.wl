(* ::Package:: *)

Quit[]


$Path


 Directory[]


AppendTo[$Path, FileNameJoin[Directory[],"/"]]


(* Notebook-friendly MFGraphs workbook.
   This file is meant to be evaluated section-by-section inside a notebook.
   It does not assume that the notebook directory is already on $Path. *)

ClearAll[
    FindMFGraphsRoot,
    PackageUsageSymbols,
    ClearMFGraphsNotebookShadows,
    EnsureMFGraphsLoaded,
    DescribeOutput,
    AssociationValue,
    NetEdgeFlows,
    NetworkVisualData,
    FlowStyleDirective,
    NetworkGraphPlot,
    SolutionFlowPlot,
    DensityNetworkGraphic,
    MassDensityCurve,
    ValueFunctionCurve,
    ExitFlowPlot,
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

PackageUsageSymbols[packageDir_String] :=
    Module[{sourceFiles, usageNames},
        sourceFiles = FileNames["*.wl", packageDir, Infinity];
        usageNames = DeleteDuplicates @ Flatten[
            StringCases[
                Import[#, "Text"],
                RegularExpression["(?m)^\\s*([A-Za-z$][A-Za-z0-9$]*)::usage\\s*="] :> "$1"
            ] & /@ sourceFiles
        ];
        usageNames
    ];

ClearMFGraphsNotebookShadows[] :=
    Module[{packageDir = FindMFGraphsRoot[], exported, shadowed},
        exported = If[packageDir === $Failed, {}, PackageUsageSymbols[packageDir]];
        shadowed = Select[exported, NameQ["Global`" <> #] &];
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

(* Reuse canonical visualization helpers from Graphics.wl to avoid workbook drift. *)
AssociationValue[assoc_Association, key_, default_: Missing["NotAvailable"]] :=
    MFGraphs`AssociationValue[assoc, key, default];

(* Net flow on each displayed edge, after eliminating dependent current variables. *)
NetEdgeFlows[d2e_Association, solution_Association, pairs_: Automatic] :=
    MFGraphs`NetEdgeFlows[d2e, solution, pairs];

NetworkVisualData[d2e_Association] :=
    MFGraphs`NetworkVisualData[d2e];

FlowStyleDirective[flow_?NumericQ, maxFlow_?NumericQ] :=
    MFGraphs`FlowStyleDirective[flow, maxFlow];

(* NetworkGraphPlot and SolutionFlowPlot now live in Graphics.wl. *)

DensityNetworkGraphic[d2e_Association, solution_Association, title_: Automatic, samples_Integer: 24] :=
    Module[{visual, coords, pairs, flows, sampleXs, densitiesByPair, allDensities,
      maxDensity, colorScale, edgeSegments, edgeLabels, vertexDisks, vertexLabels,
      plotTitle, legendRange, hasDensityModelQ, usesProxyDensityQ = False},
        visual = NetworkVisualData[d2e];
        coords = visual["coords"];
        pairs = List @@@ Lookup[d2e, "edgeList", {}];
        flows = NetEdgeFlows[d2e, solution, pairs];
        sampleXs = N @ Subdivide[0., 1., samples];
        hasDensityModelQ = DownValues[MFGraphs`Private`M] =!= {};
        densitiesByPair = AssociationThread[
            pairs,
            Table[
                Module[{flow, sampled},
                    flow = N @ AssociationValue[flows, pair, 0.];
                    sampled = If[
                        hasDensityModelQ,
                        Quiet @ Check[
                            N @ (MFGraphs`Private`M[flow, #, UndirectedEdge @@ pair] & /@ sampleXs),
                            {}
                        ],
                        {}
                    ];
                    If[VectorQ[sampled, NumericQ],
                        sampled,
                        usesProxyDensityQ = True;
                        ConstantArray[N @ Max[0., Abs[flow]], Length[sampleXs]]
                    ]
                ],
                {pair, pairs}
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
        If[usesProxyDensityQ,
            plotTitle = plotTitle <> " (flow-magnitude proxy)"
        ];
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
                    LegendLabel -> Style[
                        If[usesProxyDensityQ, "Flow-magnitude proxy", "Local density"],
                        11
                    ]
                ],
                Right
            ]
        ]
    ];

(* ExitFlowPlot now lives in Graphics.wl. *)

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


 Remove["Global`makeUnknowns", "Global`scenario", "Global`unknowns", "Global`UnknownsData"];
  Needs["MFGraphs`"];



Quit[]


(* Notebook-friendly MFGraphs workbook.
   This file is meant to be evaluated section-by-section inside a notebook.
   It does not assume that the notebook directory is already on $Path. *)

ClearAll[
    FindMFGraphsRoot,
    PackageUsageSymbols,
    ClearMFGraphsNotebookShadows,
    EnsureMFGraphsLoaded,
    DescribeOutput,
    AssociationValue,
    NetEdgeFlows,
    NetworkVisualData,
    FlowStyleDirective,
    NetworkGraphPlot,
    SolutionFlowPlot,
    DensityNetworkGraphic,
    MassDensityCurve,
    ValueFunctionCurve,
    ExitFlowPlot,
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

PackageUsageSymbols[packageDir_String] :=
    Module[{sourceFiles, usageNames},
        sourceFiles = FileNames["*.wl", packageDir, Infinity];
        usageNames = DeleteDuplicates @ Flatten[
            StringCases[
                Import[#, "Text"],
                RegularExpression["(?m)^\\s*([A-Za-z$][A-Za-z0-9$]*)::usage\\s*="] :> "$1"
            ] & /@ sourceFiles
        ];
        usageNames
    ];

ClearMFGraphsNotebookShadows[] :=
    Module[{packageDir = FindMFGraphsRoot[], exported, shadowed},
        exported = If[packageDir === $Failed, {}, PackageUsageSymbols[packageDir]];
        shadowed = Select[exported, NameQ["Global`" <> #] &];
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

(* Reuse canonical visualization helpers from Graphics.wl to avoid workbook drift. *)
AssociationValue[assoc_Association, key_, default_: Missing["NotAvailable"]] :=
    MFGraphs`AssociationValue[assoc, key, default];

(* Net flow on each displayed edge, after eliminating dependent current variables. *)
NetEdgeFlows[d2e_Association, solution_Association, pairs_: Automatic] :=
    MFGraphs`NetEdgeFlows[d2e, solution, pairs];

NetworkVisualData[d2e_Association] :=
    MFGraphs`NetworkVisualData[d2e];

FlowStyleDirective[flow_?NumericQ, maxFlow_?NumericQ] :=
    MFGraphs`FlowStyleDirective[flow, maxFlow];

(* NetworkGraphPlot and SolutionFlowPlot now live in Graphics.wl. *)

DensityNetworkGraphic[d2e_Association, solution_Association, title_: Automatic, samples_Integer: 24] :=
    Module[{visual, coords, pairs, flows, sampleXs, densitiesByPair, allDensities,
      maxDensity, colorScale, edgeSegments, edgeLabels, vertexDisks, vertexLabels,
      plotTitle, legendRange, hasDensityModelQ, usesProxyDensityQ = False},
        visual = NetworkVisualData[d2e];
        coords = visual["coords"];
        pairs = List @@@ Lookup[d2e, "edgeList", {}];
        flows = NetEdgeFlows[d2e, solution, pairs];
        sampleXs = N @ Subdivide[0., 1., samples];
        hasDensityModelQ = DownValues[MFGraphs`Private`M] =!= {};
        densitiesByPair = AssociationThread[
            pairs,
            Table[
                Module[{flow, sampled},
                    flow = N @ AssociationValue[flows, pair, 0.];
                    sampled = If[
                        hasDensityModelQ,
                        Quiet @ Check[
                            N @ (MFGraphs`Private`M[flow, #, UndirectedEdge @@ pair] & /@ sampleXs),
                            {}
                        ],
                        {}
                    ];
                    If[VectorQ[sampled, NumericQ],
                        sampled,
                        usesProxyDensityQ = True;
                        ConstantArray[N @ Max[0., Abs[flow]], Length[sampleXs]]
                    ]
                ],
                {pair, pairs}
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
        If[usesProxyDensityQ,
            plotTitle = plotTitle <> " (flow-magnitude proxy)"
        ];
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
                    LegendLabel -> Style[
                        If[usesProxyDensityQ, "Flow-magnitude proxy", "Local density"],
                        11
                    ]
                ],
                Right
            ]
        ]
    ];

(* ExitFlowPlot now lives in Graphics.wl. *)

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
(* --- 2. Typed Scenarios and Unknowns --- *)

(* Scenarios provide a typed wrapper for model data and parameters. *)
typedScenario = makeScenario[<|
    "Identity" -> <|"Name" -> "Y-junction benchmark"|>,
    "Model" -> (GetExampleData[7] /. {I1 -> 100, U1 -> 0, U2 -> 0})
|>]
Head@%


(* makeUnknowns derives the symbolic variables (js, jts, us) from the scenario. *)
unknowns = makeUnknowns[typedScenario]


Column[{
    DescribeOutput[
        "Typed scenario identity",
        "Identity metadata is carried through scenario construction.",
        ScenarioData[typedScenario, "Identity"]
    ],
    DescribeOutput[
        "Symbolic unknowns",
        "js are edge flows, jts are transition flows at junctions, and us are value-function variables.",
        unknowns
    ]
}]


(* --- 3. Define your own network --- *)

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
        NetworkGraphPlot[customD2E, "Custom network topology"]
    ]
}]


(* --- 4. Critical congestion solver --- *)

criticalData = GetExampleData[7] /. {
    I1 -> 100,
    U1 -> 0,
    U2 -> 0
};
criticalD2E = DataToEquations[criticalData];

criticalResult = CriticalCongestionSolver[criticalD2E];
criticalResult["AssoCritical"];
Column[{
    DescribeOutput[
        "Critical solver solution",
        "This association stores the equilibrium directional currents, switching currents, and edge values for the critical-congestion model.",
        criticalResult["Solution"]
    ],
    DescribeOutput[
        "Net edge flows",
        "Edge thickness and labels show the signed net flow on each interior edge. Positive values follow the stored edge orientation.",
        SolutionFlowPlot[
            criticalD2E,
            criticalResult["Solution"],
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
                criticalResult["Solution"],
                "Critical-congestion edge densities"
            ]
        ]
    ],
    DescribeOutput[
        "Net edge flow values",
        "This association is often the easiest summary to inspect numerically when checking route splits.",
        NetEdgeFlows[criticalD2E, criticalResult["Solution"]]
    ]
}]


(* --- 5. Switching-cost validation --- *)

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


(* --- 6. Jamarat / pilgrimage crowd routing --- *)

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
        SolutionFlowPlot[
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
        ExitFlowPlot[jamaratBase["ExitFlows"], "Jamarat baseline exit split"]
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
        NetworkGraphPlot[jamaratOverload["Model"], "Jamarat network structure"]
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
   reuse the same visualization and comparison pattern on the Jamarat graph.
   In practice, start with the critical-congestion symbolic solve to screen release
   plans and destination costs, then iterate on graph geometry and demand windows. *)
