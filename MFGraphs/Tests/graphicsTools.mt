(* Tests for Graphics helpers *)

Test[
    NameQ["graphicsTools`scenarioTopologyPlot"] && NameQ["graphicsTools`mfgSolutionPlot"] &&
    NameQ["graphicsTools`mfgFlowPlot"] && NameQ["graphicsTools`mfgValuePlot"] &&
    NameQ["graphicsTools`mfgDensityPlot"] && NameQ["graphicsTools`mfgTransitionPlot"] &&
    NameQ["graphicsTools`mfgAugmentedPlot"] && NameQ["graphicsTools`augmentAuxiliaryGraph"],
    True,
    TestID -> "GraphicsTools: public symbols exist"
]

Test[
    Module[{names, globalNames},
        names = {"scenarioTopologyPlot", "mfgSolutionPlot", "mfgFlowPlot",
            "mfgValuePlot", "mfgDensityPlot", "mfgTransitionPlot",
            "mfgAugmentedPlot", "augmentAuxiliaryGraph"};
        globalNames = ("Global`" <> #) & /@ names;
        Scan[If[NameQ[#], Remove[#]] &, globalNames];
        Needs["MFGraphs`"];
        And @@ (NameQ /@ names)
    ],
    True,
    TestID -> "GraphicsTools: unqualified symbols are available after Needs"
]

Test[
    Module[{s, sys, sol, topoPlot, solPlot, flowPlot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        topoPlot = scenarioTopologyPlot[s, sys];
        solPlot = mfgSolutionPlot[s, sys, sol];
        flowPlot = mfgFlowPlot[s, sys, sol];
        MatchQ[topoPlot, _Graph] && MatchQ[solPlot, _Graphics | _GraphicsGrid] &&
        MatchQ[flowPlot, _Graph]
    ],
    True,
    TestID -> "GraphicsTools: topology, solution grid, and flow plots return renderable expressions"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        edges = EdgeList[mfgFlowPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[1, 2]]
    ],
    True,
    TestID -> "GraphicsTools: mfgFlowPlot positive j[1,2] displays as DirectedEdge[1,2]"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{3, 120.0}}, {{1, 0.0}, {2, 10.0}}];
        sys = makeSystem[s];
        sol = {j[1, 2] -> 0, j[2, 1] -> 5};
        edges = EdgeList[mfgFlowPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[2, 1]]
    ],
    True,
    TestID -> "GraphicsTools: mfgFlowPlot negative net flow flips display direction to DirectedEdge[2,1]"
]

Test[
    Module[{s, sys, sol, graph, edges},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        graph = mfgFlowPlot[s, sys, sol];
        edges = EdgeList[graph];
        AllTrue[edges, MatchQ[#, _DirectedEdge] &] &&
        Length[edges] == Length[systemData[sys, "AuxEdges"]]
    ],
    True,
    TestID -> "GraphicsTools: mfgFlowPlot displays every auxiliary edge as directed"
]

Test[
    Module[{s, sys, sol, graph, labelText},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        graph = mfgFlowPlot[s, sys, sol];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        StringContainsQ[labelText, "EdgeLabels"] && !StringContainsQ[labelText, "u="]
    ],
    True,
    TestID -> "GraphicsTools: mfgFlowPlot edge labels omit value functions"
]

Test[
    Module[{s, sys, sol, graph, edgeStyleText},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        graph = mfgFlowPlot[s, sys, sol];
        edgeStyleText = ToString[InputForm[AbsoluteOptions[graph, EdgeStyle]]];
        StringContainsQ[edgeStyleText, "AbsoluteThickness[2."]
    ],
    True,
    TestID -> "GraphicsTools: mfgFlowPlot uses fixed edge thickness"
]

Test[
    Module[{s, sys, topoPlot, solPlot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        topoPlot = scenarioTopologyPlot[s, sys, "Topology custom title"];
        solPlot = mfgSolutionPlot[s, sys, <||>, "Solution custom title"];
        !FreeQ[AbsoluteOptions[topoPlot, PlotLabel], "Topology custom title"] &&
        !FreeQ[AbsoluteOptions[solPlot, PlotLabel], "Solution custom title"]
    ],
    True,
    TestID -> "GraphicsTools: custom titles propagate to PlotLabel"
]

Test[
    Module[{s, sys, sol, plots},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        plots = {
            scenarioTopologyPlot[s, sys, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Topology option", ImageSize -> Medium],
            mfgSolutionPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Solution option", ImageSize -> Medium],
            mfgFlowPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Flow option", ImageSize -> Medium],
            mfgValuePlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Value option", ImageSize -> Medium],
            mfgDensityPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Density option", ImageSize -> Medium],
            mfgTransitionPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Transition option", ImageSize -> Medium,
                ShowLegend -> False],
            mfgAugmentedPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Augmented option", ImageSize -> Medium,
                ShowLegend -> False]
        };
        MatchQ[plots[[2]], _Graphics | _GraphicsGrid] &&
        AllTrue[Delete[plots, 2], MatchQ[#, _Graph] &] &&
        And @@ MapThread[
            !FreeQ[AbsoluteOptions[#1, PlotLabel], #2] &,
            {plots, {"Topology option", "Solution option", "Flow option",
                "Value option", "Density option", "Transition option", "Augmented option"}}
        ] &&
        AllTrue[
            Delete[plots, 2],
            !FreeQ[AbsoluteOptions[#, GraphLayout], "LayeredDigraphEmbedding"] &
        ]
    ],
    True,
    TestID -> "GraphicsTools: graph helpers accept PlotLabel GraphLayout and ImageSize options"
]

Test[
    Module[{s, sys, sol, plots},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        plots = {
            scenarioTopologyPlot[s, sys, PlotLabel -> None],
            mfgSolutionPlot[s, sys, sol, PlotLabel -> None],
            mfgFlowPlot[s, sys, sol, PlotLabel -> None],
            mfgValuePlot[s, sys, sol, PlotLabel -> None],
            mfgDensityPlot[s, sys, sol, PlotLabel -> None],
            mfgTransitionPlot[s, sys, sol, PlotLabel -> None, ShowLegend -> False],
            mfgAugmentedPlot[s, sys, sol, PlotLabel -> None, ShowLegend -> False]
        };
        MatchQ[plots[[2]], _Graphics | _GraphicsGrid] &&
        AllTrue[Delete[plots, 2], MatchQ[#, _Graph] &] &&
        AllTrue[plots, FreeQ[AbsoluteOptions[#, PlotLabel], _Style] &]
    ],
    True,
    TestID -> "GraphicsTools: PlotLabel -> None suppresses styled label for all plot helpers"
]

Test[
    Module[{s, sys, sol, transitionPlot, augmentedPlot, transitionFile,
            augmentedFile, transitionDims, augmentedDims},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        transitionPlot = mfgTransitionPlot[s, sys, sol,
            GraphLayout -> "LayeredDigraphEmbedding"];
        augmentedPlot = mfgAugmentedPlot[s, sys, sol,
            GraphLayout -> "LayeredDigraphEmbedding"];
        transitionFile = FileNameJoin[{$TemporaryDirectory, "mfg-transition-render-test.png"}];
        augmentedFile = FileNameJoin[{$TemporaryDirectory, "mfg-augmented-render-test.png"}];
        Export[transitionFile, transitionPlot];
        Export[augmentedFile, augmentedPlot];
        transitionDims = ImageDimensions @ Import[transitionFile];
        augmentedDims = ImageDimensions @ Import[augmentedFile];
        MatchQ[transitionPlot, _Graph | _Legended] &&
        MatchQ[augmentedPlot, _Graph | _Legended] &&
        Min[transitionDims] > 100 &&
        Min[augmentedDims] > 100
    ],
    True,
    TestID -> "GraphicsTools: transition and augmented plots render to nontrivial images"
]

Test[
    Module[{s, sys, graph, labelText},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        graph = mfgValuePlot[s, sys, {u[1, 2] -> 0., u[2, 1] -> 4.}];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        MatchQ[graph, _Graph] &&
        StringContainsQ[labelText, "1/2"] &&
        StringContainsQ[labelText, "2."]
    ],
    True,
    TestID -> "GraphicsTools: mfgValuePlot labels include interpolated midpoint values"
]

Test[
    Module[{s, sys, graph, labelText},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        graph = mfgDensityPlot[s, sys, {j[1, 2] -> 0, j[2, 1] -> 0}];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        MatchQ[graph, _Graph] && StringContainsQ[labelText, "m="] &&
        StringContainsQ[labelText, "0."]
    ],
    True,
    TestID -> "GraphicsTools: mfgDensityPlot displays zero density for zero flow"
]

Test[
    Module[{s, sys, graph, labelText},
        s = makeScenario[<|"Model" -> <|
            "Vertices" -> {1, 2},
            "Adjacency" -> {{0, 1}, {1, 0}},
            "Entries" -> {{1, 1}},
            "Exits" -> {{2, 0}},
            "Switching" -> {}
        |>|>];
        sys = makeSystem[s];
        graph = mfgDensityPlot[s, sys, {j[1, 2] -> 2, j[2, 1] -> 0}];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        MatchQ[graph, _Graph] && StringContainsQ[labelText, "3."]
    ],
    True,
    TestID -> "GraphicsTools: mfgDensityPlot solves positive density with default Hamiltonian"
]

Test[
    Module[{s, sys, graph, labelText},
        s = makeScenario[<|
            "Model" -> <|
                "Vertices" -> {1, 2},
                "Adjacency" -> {{0, 1}, {1, 0}},
                "Entries" -> {{1, 1}},
                "Exits" -> {{2, 0}},
                "Switching" -> {}
            |>,
            "Hamiltonian" -> <|
                "Alpha" -> 1,
                "V" -> -2,
                "G" -> Function[z, -2/z]
            |>
        |>];
        sys = makeSystem[s];
        graph = mfgDensityPlot[s, sys, {j[1, 2] -> 2, j[2, 1] -> 0}];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        MatchQ[graph, _Graph] && StringContainsQ[labelText, "2."]
    ],
    True,
    TestID -> "GraphicsTools: mfgDensityPlot uses scenario Hamiltonian defaults"
]

Test[
    Module[{s, sBroken, sys, graph, labelText},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        sBroken = scenario[Join[scenarioData[s], <|"Hamiltonian" -> <||>|>]];
        graph = mfgDensityPlot[sBroken, sys, {j[1, 2] -> 2, j[2, 1] -> 0}];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        MatchQ[graph, _Graph] && StringContainsQ[labelText, "?"]
    ],
    True,
    TestID -> "GraphicsTools: mfgDensityPlot reports malformed Hamiltonian as unknown density"
]

Test[
    Module[{s, sys, graph, labelText},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        graph = mfgDensityPlot[s, sys, <||>];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        MatchQ[graph, _Graph] && StringContainsQ[labelText, "?"]
    ],
    True,
    TestID -> "GraphicsTools: mfgDensityPlot renders missing solution rules as unknown density"
]

Test[
    Module[{s, sys, sol, rules},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        rules = Replace[sol, a_Association :> Lookup[a, "Rules", {}]];
        MatchQ[mfgSolutionPlot[s, sys, rules], _Graphics | _GraphicsGrid] &&
        MatchQ[mfgFlowPlot[s, sys, rules], _Graph] &&
        MatchQ[mfgValuePlot[s, sys, rules], _Graph] &&
        MatchQ[mfgDensityPlot[s, sys, rules], _Graph]
    ],
    True,
    TestID -> "GraphicsTools: solution panels accept raw rule list solution"
]

Test[
    Module[{s, sys, graph, flowGraph, valueGraph, densityGraph},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        graph = mfgSolutionPlot[s, sys, <||>];
        flowGraph = mfgFlowPlot[s, sys, <||>];
        valueGraph = mfgValuePlot[s, sys, <||>];
        densityGraph = mfgDensityPlot[s, sys, <||>];
        MatchQ[graph, _Graphics | _GraphicsGrid] &&
        MatchQ[flowGraph, _Graph] && Length[EdgeList[flowGraph]] > 0 &&
        MatchQ[valueGraph, _Graph] && Length[EdgeList[valueGraph]] > 0 &&
        MatchQ[densityGraph, _Graph] && Length[EdgeList[densityGraph]] > 0
    ],
    True,
    TestID -> "GraphicsTools: solution panels handle missing Rules association without failure"
]

Test[
    Module[{s, sys, augmented},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        augmented = augmentAuxiliaryGraph[sys];
        AssociationQ[augmented] &&
        MatchQ[augmented["Graph"], _Graph] &&
        KeyExistsQ[augmented, "Vertices"] &&
        KeyExistsQ[augmented, "FlowEdges"] &&
        KeyExistsQ[augmented, "TransitionEdges"] &&
        KeyExistsQ[augmented, "EdgeVariables"] &&
        KeyExistsQ[augmented, "EdgeKinds"]
    ],
    True,
    TestID -> "GraphicsTools: augmentAuxiliaryGraph returns graph metadata"
]

Test[
    Module[{s, sys, augmented, flowEdges, transitionEdges},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        augmented = augmentAuxiliaryGraph[sys];
        flowEdges = augmented["FlowEdges"];
        transitionEdges = augmented["TransitionEdges"];
        MemberQ[flowEdges, DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}]] &&
        MemberQ[flowEdges, DirectedEdge[{"auxExit3", 3}, {3, "auxExit3"}]] &&
        MemberQ[transitionEdges, DirectedEdge[{"auxEntry1", 1}, {2, 1}]] &&
        MemberQ[transitionEdges, DirectedEdge[{2, 3}, {"auxExit3", 3}]] &&
        !MemberQ[transitionEdges, DirectedEdge[{1, "auxEntry1"}, {1, 2}]]
    ],
    True,
    TestID -> "GraphicsTools: augmentAuxiliaryGraph uses road-traffic edge directions"
]

Test[
    Module[{s, sys, augmented, edgeVariables, entryFlowEdge, exitFlowEdge, transitionEdge},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        augmented = augmentAuxiliaryGraph[sys];
        edgeVariables = augmented["EdgeVariables"];
        entryFlowEdge = DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}];
        exitFlowEdge = DirectedEdge[{"auxExit3", 3}, {3, "auxExit3"}];
        transitionEdge = DirectedEdge[{"auxEntry1", 1}, {2, 1}];
        edgeVariables[entryFlowEdge] === j["auxEntry1", 1] &&
        edgeVariables[exitFlowEdge] === j[3, "auxExit3"] &&
        edgeVariables[transitionEdge] === j["auxEntry1", 1, 2]
    ],
    True,
    TestID -> "GraphicsTools: augmentAuxiliaryGraph maps edges to canonical j variables"
]

Test[
    Module[{s, sys, augmented, edgeKinds, entryFlowEdge, transitionEdge},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        augmented = augmentAuxiliaryGraph[sys];
        edgeKinds = augmented["EdgeKinds"];
        entryFlowEdge = DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}];
        transitionEdge = DirectedEdge[{"auxEntry1", 1}, {2, 1}];
        edgeKinds[entryFlowEdge] === "Flow" &&
        edgeKinds[transitionEdge] === "Transition"
    ],
    True,
    TestID -> "GraphicsTools: augmentAuxiliaryGraph records flow and transition edge kinds"
]

Test[
    Module[{s, sys, graph, edges},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        graph = mfgAugmentedPlot[s, sys, <||>, ShowBoundaryValues -> False];
        edges = EdgeList[graph];
        MatchQ[graph, _Graph] &&
        MemberQ[edges, DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}]] &&
        MemberQ[edges, DirectedEdge[{"auxEntry1", 1}, {2, 1}]] &&
        !MemberQ[edges, DirectedEdge[{1, "auxEntry1"}, {1, 2}]]
    ],
    True,
    TestID -> "GraphicsTools: mfgAugmentedPlot uses augmented road-traffic graph"
]
