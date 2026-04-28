(* Tests for Graphics helpers *)

Test[
    NameQ["graphics`scenarioTopologyPlot"] && NameQ["graphics`mfgSolutionPlot"] &&
    NameQ["graphics`mfgFlowPlot"] && NameQ["graphics`mfgTransitionPlot"] &&
    NameQ["graphics`mfgAugmentedPlot"] && NameQ["graphics`augmentAuxiliaryGraph"],
    True,
    TestID -> "Graphics: public symbols exist"
]

Test[
    Module[{names, globalNames},
        names = {"scenarioTopologyPlot", "mfgSolutionPlot", "mfgFlowPlot",
            "mfgTransitionPlot", "mfgAugmentedPlot", "augmentAuxiliaryGraph"};
        globalNames = ("Global`" <> #) & /@ names;
        Scan[If[NameQ[#], Remove[#]] &, globalNames];
        Needs["MFGraphs`"];
        And @@ (NameQ /@ names)
    ],
    True,
    TestID -> "Graphics: unqualified symbols are available after Needs"
]

Test[
    Module[{s, sys, sol, topoPlot, solPlot, flowPlot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        topoPlot = scenarioTopologyPlot[s, sys];
        solPlot = mfgSolutionPlot[s, sys, sol];
        flowPlot = mfgFlowPlot[s, sys, sol];
        MatchQ[topoPlot, _Graph] && MatchQ[solPlot, _Graph] && MatchQ[flowPlot, _Graph]
    ],
    True,
    TestID -> "Graphics: topology, solution, and flow plots return Graph"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        edges = EdgeList[mfgSolutionPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[1, 2]]
    ],
    True,
    TestID -> "Graphics: positive j[1,2] displays as DirectedEdge[1,2]"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        edges = EdgeList[mfgFlowPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[1, 2]]
    ],
    True,
    TestID -> "Graphics: mfgFlowPlot positive j[1,2] displays as DirectedEdge[1,2]"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{3, 120.0}}, {{1, 0.0}, {2, 10.0}}];
        sys = makeSystem[s];
        sol = {j[1, 2] -> 0, j[2, 1] -> 5};
        edges = EdgeList[mfgSolutionPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[2, 1]]
    ],
    True,
    TestID -> "Graphics: negative net flow flips display direction to DirectedEdge[2,1]"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{3, 120.0}}, {{1, 0.0}, {2, 10.0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        edges = EdgeList[mfgFlowPlot[s, sys, sol]];
        ListQ[sol] && MemberQ[edges, DirectedEdge[2, 1]]
    ],
    True,
    TestID -> "Graphics: mfgFlowPlot negative net flow flips display direction to DirectedEdge[2,1]"
]

Test[
    Module[{s, sys, sol, graph, edges},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        graph = mfgFlowPlot[s, sys, sol];
        edges = EdgeList[graph];
        AllTrue[edges, MatchQ[#, _DirectedEdge] &] &&
        Length[edges] == Length[systemData[sys, "AuxEdges"]]
    ],
    True,
    TestID -> "Graphics: mfgFlowPlot displays every auxiliary edge as directed"
]

Test[
    Module[{s, sys, sol, graph, labelText},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        graph = mfgFlowPlot[s, sys, sol];
        labelText = ToString[InputForm[AbsoluteOptions[graph, EdgeLabels]]];
        StringContainsQ[labelText, "EdgeLabels"] && !StringContainsQ[labelText, "u="]
    ],
    True,
    TestID -> "Graphics: mfgFlowPlot edge labels omit value functions"
]

Test[
    Module[{s, sys, sol, graph, edgeStyleText},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        graph = mfgFlowPlot[s, sys, sol];
        edgeStyleText = ToString[InputForm[AbsoluteOptions[graph, EdgeStyle]]];
        StringContainsQ[edgeStyleText, "AbsoluteThickness[2."]
    ],
    True,
    TestID -> "Graphics: mfgFlowPlot uses fixed edge thickness"
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
    TestID -> "Graphics: custom titles propagate to PlotLabel"
]

Test[
    Module[{s, sys, sol, plots},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        plots = {
            scenarioTopologyPlot[s, sys, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Topology option", ImageSize -> Medium],
            mfgSolutionPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Solution option", ImageSize -> Medium],
            mfgFlowPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Flow option", ImageSize -> Medium],
            mfgTransitionPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Transition option", ImageSize -> Medium],
            mfgAugmentedPlot[s, sys, sol, GraphLayout -> "LayeredDigraphEmbedding",
                PlotLabel -> "Augmented option", ImageSize -> Medium]
        };
        AllTrue[plots, MatchQ[#, _Graph] &] &&
        And @@ MapThread[
            !FreeQ[AbsoluteOptions[#1, PlotLabel], #2] &,
            {plots, {"Topology option", "Solution option", "Flow option",
                "Transition option", "Augmented option"}}
        ] &&
        AllTrue[
            plots,
            !FreeQ[AbsoluteOptions[#, GraphLayout], "LayeredDigraphEmbedding"] &
        ]
    ],
    True,
    TestID -> "Graphics: graph helpers accept PlotLabel GraphLayout and ImageSize options"
]

Test[
    Module[{s, sys, sol, transitionPlot, augmentedPlot, transitionFile,
            augmentedFile, transitionDims, augmentedDims},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
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
        MatchQ[transitionPlot, _Graph] &&
        MatchQ[augmentedPlot, _Graph] &&
        Min[transitionDims] > 100 &&
        Min[augmentedDims] > 100
    ],
    True,
    TestID -> "Graphics: transition and augmented plots render to nontrivial images"
]

Test[
    Module[{s, sys, sol, rules},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        rules = Replace[sol, a_Association :> Lookup[a, "Rules", {}]];
        MatchQ[mfgSolutionPlot[s, sys, rules], _Graph] &&
        MatchQ[mfgFlowPlot[s, sys, rules], _Graph]
    ],
    True,
    TestID -> "Graphics: solution and flow plots accept raw rule list solution"
]

Test[
    Module[{s, sys, graph, flowGraph},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        graph = mfgSolutionPlot[s, sys, <||>];
        flowGraph = mfgFlowPlot[s, sys, <||>];
        MatchQ[graph, _Graph] && Length[EdgeList[graph]] > 0 &&
        MatchQ[flowGraph, _Graph] && Length[EdgeList[flowGraph]] > 0
    ],
    True,
    TestID -> "Graphics: solution and flow plots handle missing Rules association without failure"
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
    TestID -> "Graphics: augmentAuxiliaryGraph returns graph metadata"
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
    TestID -> "Graphics: augmentAuxiliaryGraph uses road-traffic edge directions"
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
    TestID -> "Graphics: augmentAuxiliaryGraph maps edges to canonical j variables"
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
    TestID -> "Graphics: augmentAuxiliaryGraph records flow and transition edge kinds"
]

Test[
    Module[{s, sys, graph, edges},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        graph = mfgAugmentedPlot[s, sys, <||>];
        edges = EdgeList[graph];
        MatchQ[graph, _Graph] &&
        MemberQ[edges, DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}]] &&
        MemberQ[edges, DirectedEdge[{"auxEntry1", 1}, {2, 1}]] &&
        !MemberQ[edges, DirectedEdge[{1, "auxEntry1"}, {1, 2}]]
    ],
    True,
    TestID -> "Graphics: mfgAugmentedPlot uses augmented road-traffic graph"
]
