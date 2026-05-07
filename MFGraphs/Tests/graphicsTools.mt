(* Tests for the consolidated graphics surface: rawNetworkPlot, richNetworkPlot,
   and augmentAuxiliaryGraph. *)

Test[
    NameQ["graphicsTools`rawNetworkPlot"] &&
    NameQ["graphicsTools`richNetworkPlot"] &&
    NameQ["graphicsTools`augmentAuxiliaryGraph"],
    True,
    TestID -> "GraphicsTools: public symbols exist"
]

Test[
    Module[{names, globalNames},
        names = {"rawNetworkPlot", "richNetworkPlot", "augmentAuxiliaryGraph"};
        globalNames = ("Global`" <> #) & /@ names;
        Scan[If[NameQ[#], Remove[#]] &, globalNames];
        Needs["MFGraphs`"];
        And @@ (NameQ /@ names)
    ],
    True,
    TestID -> "GraphicsTools: unqualified symbols are available after Needs"
]

Test[
    Module[{names},
        names = {"scenarioTopologyPlot", "mfgSolutionPlot", "mfgFlowPlot",
                 "mfgValuePlot", "mfgDensityPlot", "mfgTransitionPlot",
                 "mfgAugmentedPlot", "mfgAugmentedBoundaryPlot"};
        AllTrue[names, !NameQ["graphicsTools`" <> #] &]
    ],
    True,
    TestID -> "GraphicsTools: deprecated public symbols are removed"
]

Test[
    Module[{s, sys, sol, plot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        plot = rawNetworkPlot[s, sys, sol];
        MatchQ[plot, _Graph | _Legended]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot returns renderable expression"
]

Test[
    Module[{s, sys, sol, plot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        plot = richNetworkPlot[s, sys, sol];
        MatchQ[plot, _Graph | _Legended]
    ],
    True,
    TestID -> "GraphicsTools: richNetworkPlot returns renderable expression"
]

Test[
    Module[{s, sys, plot, vlist},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys];
        vlist = VertexList[plot];
        AllTrue[{1, 2, 3}, MemberQ[vlist, #] &] &&
        AllTrue[vlist, NumericQ]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot without sol shows real vertices only"
]

Test[
    Module[{s, sys, plot, vlist},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys, ShowAuxiliaryVertices -> True];
        vlist = VertexList[plot];
        MemberQ[vlist, "auxEntry1"] && MemberQ[vlist, "auxExit2"] &&
        MemberQ[vlist, "auxExit3"]
    ],
    True,
    TestID -> "GraphicsTools: ShowAuxiliaryVertices includes auxiliary vertices"
]

Test[
    Module[{s, sys, plot, vlist},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        (* ShowBoundaryData should force aux visibility *)
        plot = rawNetworkPlot[s, sys, ShowBoundaryData -> True,
                              ShowAuxiliaryVertices -> False];
        vlist = VertexList[plot];
        MemberQ[vlist, "auxEntry1"]
    ],
    True,
    TestID -> "GraphicsTools: ShowBoundaryData forces aux visibility"
]

Test[
    Module[{s, sys, plot, styles, vstyle, entryColor, exitColor, internalColor},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys];
        vstyle = AbsoluteOptions[plot, VertexStyle][[1, 2]];
        styles = Association[vstyle];
        entryColor = Lookup[styles, 1, Missing[]];
        exitColor = Lookup[styles, 2, Missing[]];
        internalColor = Lookup[styles, 3, Missing[]];
        (* All real vertices should be gray (GrayLevel) *)
        MatchQ[entryColor, _GrayLevel] &&
        MatchQ[exitColor, _GrayLevel] &&
        MatchQ[internalColor, _GrayLevel]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot real vertices are gray (no solution)"
]

Test[
    Module[{s, sys, sol, plot, vstyle, styles, entryColor},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        plot = rawNetworkPlot[s, sys, sol];
        vstyle = AbsoluteOptions[
            Replace[plot, Legended[g_, _] :> g], VertexStyle][[1, 2]];
        styles = Association[vstyle];
        entryColor = Lookup[styles, 1, Missing[]];
        (* Real vertex 1 must remain gray even with a solution *)
        MatchQ[entryColor, _GrayLevel]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot real vertices stay gray with solution"
]

Test[
    Module[{s, sys, plot, vstyle, styles, entryColor, exitColor},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys, ShowAuxiliaryVertices -> True];
        vstyle = AbsoluteOptions[plot, VertexStyle][[1, 2]];
        styles = Association[vstyle];
        entryColor = Lookup[styles, "auxEntry1", Missing[]];
        exitColor = Lookup[styles, "auxExit3", Missing[]];
        MatchQ[entryColor, _RGBColor] && MatchQ[exitColor, _RGBColor]
    ],
    True,
    TestID -> "GraphicsTools: aux vertices receive category colors when shown"
]

Test[
    Module[{s, sys, sol, ruleList, assoc, plotA, plotB},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        ruleList = Replace[sol, a_Association :> Lookup[a, "Rules", {}]];
        assoc = <|"Rules" -> ruleList, "Residual" -> True|>;
        plotA = rawNetworkPlot[s, sys, ruleList];
        plotB = rawNetworkPlot[s, sys, assoc];
        MatchQ[plotA, _Graph | _Legended] && MatchQ[plotB, _Graph | _Legended]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot accepts rule list and association forms"
]

Test[
    Module[{s, sys, plot, edges},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = richNetworkPlot[s, sys];
        edges = EdgeList[Replace[plot, Legended[g_, _] :> g]];
        MemberQ[edges, DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}]] &&
        MemberQ[edges, DirectedEdge[{"auxEntry1", 1}, {2, 1}]]
    ],
    True,
    TestID -> "GraphicsTools: richNetworkPlot includes flow and transition edges by default"
]

Test[
    Module[{s, sys, plot, edges, kinds},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = richNetworkPlot[s, sys, ShowFlowEdges -> False];
        edges = EdgeList[Replace[plot, Legended[g_, _] :> g]];
        (* Without flow edges, the entry-flow self-loop arc should not appear *)
        !MemberQ[edges, DirectedEdge[{1, "auxEntry1"}, {"auxEntry1", 1}]] &&
        Length[edges] > 0
    ],
    True,
    TestID -> "GraphicsTools: richNetworkPlot ShowFlowEdges->False yields transition-only graph"
]

Test[
    Module[{s, sys, plot, vstyle, styles, inNode, otherNode},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = richNetworkPlot[s, sys];
        vstyle = AbsoluteOptions[
            Replace[plot, Legended[g_, _] :> g], VertexStyle][[1, 2]];
        styles = Association[vstyle];
        (* {realV, auxEntry} is the "in" node — should be green RGBColor *)
        inNode = Lookup[styles, Key[{1, "auxEntry1"}], Missing[]];
        (* {auxEntry, realV} has aux at position 1 — should be gray *)
        otherNode = Lookup[styles, Key[{"auxEntry1", 1}], Missing[]];
        MatchQ[inNode, _RGBColor] && MatchQ[otherNode, _GrayLevel]
    ],
    True,
    TestID -> "GraphicsTools: richNetworkPlot only colors label-position aux nodes"
]

Test[
    Module[{s, sys, sol, plot, ruleList},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        ruleList = Replace[sol, a_Association :> Lookup[a, "Rules", {}]];
        plot = richNetworkPlot[s, sys, ruleList];
        MatchQ[plot, _Graph | _Legended]
    ],
    True,
    TestID -> "GraphicsTools: richNetworkPlot accepts rule list solution"
]

Test[
    Module[{s, sys, plotA, plotB},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plotA = rawNetworkPlot[s, sys, "Custom raw title"];
        plotB = richNetworkPlot[s, sys, "Custom rich title"];
        !FreeQ[AbsoluteOptions[Replace[plotA, Legended[g_, _] :> g], PlotLabel],
               "Custom raw title"] &&
        !FreeQ[AbsoluteOptions[Replace[plotB, Legended[g_, _] :> g], PlotLabel],
               "Custom rich title"]
    ],
    True,
    TestID -> "GraphicsTools: positional title argument propagates to PlotLabel"
]

Test[
    Module[{s, sys, sol, plot},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        sol = solveScenario[s];
        plot = rawNetworkPlot[s, sys, sol,
            GraphLayout -> "LayeredDigraphEmbedding",
            PlotLabel -> "Layout test",
            ImageSize -> Medium,
            ShowLegend -> False];
        !FreeQ[AbsoluteOptions[Replace[plot, Legended[g_, _] :> g], GraphLayout],
               "LayeredDigraphEmbedding"]
    ],
    True,
    TestID -> "GraphicsTools: GraphLayout / PlotLabel / ImageSize options propagate"
]

Test[
    Module[{s, sys, plot},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys, PlotLabel -> None];
        FreeQ[AbsoluteOptions[plot, PlotLabel], _Style]
    ],
    True,
    TestID -> "GraphicsTools: PlotLabel -> None suppresses styled label"
]

Test[
    Module[{s, sys, plot, file, dims},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = richNetworkPlot[s, sys, GraphLayout -> "LayeredDigraphEmbedding"];
        file = FileNameJoin[{$TemporaryDirectory, "rich-render-test.png"}];
        Export[file, plot];
        dims = ImageDimensions @ Import[file];
        Min[dims] > 100
    ],
    True,
    TestID -> "GraphicsTools: richNetworkPlot renders to nontrivial image"
]

Test[
    Module[{s, sys, plot, labelText},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys, {u[1, 2] -> 0., u[2, 1] -> 4.},
                              ShowFlowLabels -> False, ShowValueLabels -> True];
        labelText = ToString[InputForm[
            AbsoluteOptions[Replace[plot, Legended[g_, _] :> g], EdgeLabels]]];
        StringContainsQ[labelText, "1/2"]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot ShowValueLabels shows interpolated u midpoints"
]

Test[
    Module[{s, sys, plot, labelText},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        plot = rawNetworkPlot[s, sys, {j[1, 2] -> 0, j[2, 1] -> 0},
                              ShowFlowLabels -> False, ShowDensityLabels -> True];
        labelText = ToString[InputForm[
            AbsoluteOptions[Replace[plot, Legended[g_, _] :> g], EdgeLabels]]];
        StringContainsQ[labelText, "m="]
    ],
    True,
    TestID -> "GraphicsTools: rawNetworkPlot ShowDensityLabels shows m= prefix"
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
