(* Wolfram Language package *)
(* Graphics: public visualization helpers for MFGraphs systems and solutions. *)

BeginPackage["graphicsTools`", {"primitives`", "utilities`", "scenarioTools`", "systemTools`"}];

scenarioTopologyPlot::usage =
"scenarioTopologyPlot[s, sys, opts] plots the scenario topology using vertex coloring for entry, exit, and internal vertices. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

mfgSolutionPlot::usage =
"mfgSolutionPlot[s, sys, sol, opts] plots coordinated solution views: flows, value samples, and agent density. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

mfgFlowPlot::usage =
"mfgFlowPlot[s, sys, sol, opts] plots a flow-only solution graph with real and auxiliary edges. Edges are displayed as directed, and edge labels show only j flow values. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

mfgValuePlot::usage =
"mfgValuePlot[s, sys, sol, opts] plots real network edges with endpoint value variables u[a,b], u[b,a] and linearly interpolated interior samples at 1/4, 1/2, and 3/4. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

mfgDensityPlot::usage =
"mfgDensityPlot[s, sys, sol, opts] plots real network edges with agent density inferred from solved signed spatial flow and the scenario Hamiltonian density equation. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

mfgTransitionPlot::usage =
"mfgTransitionPlot[s, sys, sol, opts] plots the transition graph of the solution. Nodes are AuxPair states {r,i}->{i,w}; directed edges represent transition flows j[r,i,w]. Nodes are labeled with internal values u where available. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

mfgAugmentedPlot::usage =
"mfgAugmentedPlot[s, sys, sol, opts] plots the augmented road-traffic infrastructure graph built by augmentAuxiliaryGraph. Blue edges represent flow variables j[a,b]; red edges represent transition variables j[r,i,w]. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

augmentAuxiliaryGraph::usage =
"augmentAuxiliaryGraph[sys] constructs the road-traffic augmented infrastructure graph from a system's AuxPairs and AuxTriples. Returns an Association containing the Graph, Vertices, FlowEdges, TransitionEdges, EdgeVariables, and EdgeKinds.";

Begin["`Private`"];

Options[scenarioTopologyPlot] = {
    GraphLayout -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

Options[mfgSolutionPlot] = {
    GraphLayout -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Large
};

Options[mfgFlowPlot] = Options[mfgSolutionPlot];

Options[mfgValuePlot] = Options[mfgSolutionPlot];

Options[mfgDensityPlot] = Options[mfgSolutionPlot];

Options[mfgTransitionPlot] = {
    GraphLayout -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Large
};

Options[mfgAugmentedPlot] = Options[mfgTransitionPlot];

netEdgeFlow::usage =
"netEdgeFlow[a, b, rules] returns the solved net flow j[a,b]-j[b,a], defaulting absent flow variables to zero.";

edgeJValue::usage =
"edgeJValue[a, b, rules] returns the solved flow j[a,b], defaulting absent flow variables to zero.";

edgeUValue::usage =
"edgeUValue[a, b, rules] returns the solved value variable for edge endpoints, or Missing[] when absent.";

directedDisplayEdges::usage =
"directedDisplayEdges[edges, rules] orients display edges using solved net flow values.";

stateLabel::usage =
"stateLabel[pair] returns a compact display label for an augmented graph state pair.";

plotLabelValue::usage =
"plotLabelValue[label, default] returns a styled plot label unless label is None.";

realNetworkVertexStyle::usage =
"realNetworkVertexStyle[s] returns vertex styles for real entry, exit, and internal network vertices.";

realNetworkVertexSize::usage =
"realNetworkVertexSize[s] returns vertex sizes for real network vertices.";

edgeHamiltonianValue::usage =
"edgeHamiltonianValue[ham, key, edge] returns scenario Hamiltonian data for an edge with reverse-edge fallback.";

edgeDensityValue::usage =
"edgeDensityValue[s, sys, edge, rules] computes inferred density for a real edge.";

formatPlotNumber::usage =
"formatPlotNumber[value] formats numeric plot labels and returns ? for unavailable values.";

netEdgeFlow[a_, b_, rules_List] :=
    (j[a, b] - j[b, a]) /. rules /. {_j -> 0};

edgeJValue[a_, b_, rules_List] :=
    (j[a, b] /. rules /. {_j -> 0});

edgeUValue[a_, b_, rules_List] :=
    Module[{v1, v2},
        v1 = u[a, b] /. rules;
        v2 = u[b, a] /. rules;
        Which[
            v1 =!= u[a, b], v1,
            v2 =!= u[b, a], v2,
            True, Missing[]
        ]
    ];

directedDisplayEdges[edges_List, rules_List] :=
    edges /. {
        UndirectedEdge[a_, b_] :> Module[{f},
            f = netEdgeFlow[a, b, rules];
            If[NumericQ[f] && f < 0, DirectedEdge[b, a], DirectedEdge[a, b]]
        ],
        DirectedEdge[a_, b_] :> DirectedEdge[a, b]
    };

stateAtomLabel[x_] :=
    StringReplace[ToString[x], {"auxEntry" -> "in", "auxExit" -> "out"}];

stateLabel[{a_, b_}] :=
    Row[{"{", stateAtomLabel[a], ", ", stateAtomLabel[b], "}"}];

plotLabelValue[label_, default_String] :=
    Replace[label, {
        Automatic -> Style[default, 14, Bold],
        None -> None,
        other_ :> Style[other, 14, Bold]
    }];

realNetworkVertexStyle[s_?scenarioQ] :=
    Module[{model, entryV, exitV, realV, internalV},
        model     = scenarioData[s, "Model"];
        entryV    = First /@ model["Entries"];
        exitV     = First /@ model["Exits"];
        realV     = model["Vertices"];
        internalV = Complement[realV, entryV, exitV];
        Normal @ Join[
            AssociationThread[entryV,    RGBColor[0.22, 0.6, 0.3]],
            AssociationThread[exitV,     RGBColor[0.82, 0.27, 0.2]],
            AssociationThread[internalV, GrayLevel[0.72]]
        ]
    ];

realNetworkVertexSize[s_?scenarioQ] :=
    AssociationThread[scenarioData[s, "Model"]["Vertices"], 0.38];

formatPlotNumber[value_?NumericQ] := NumberForm[N[value], {5, 2}];
formatPlotNumber[_] := "?";

edgeHamiltonianValue[ham_Association, key_String, edge_List] :=
    Module[{edgeAssoc, globalValue},
        edgeAssoc = Lookup[ham, "Edge" <> key, <||>];
        globalValue = Lookup[ham, key, Missing["MissingHamiltonianKey", key]];
        If[!AssociationQ[edgeAssoc] || MissingQ[globalValue],
            Missing["InvalidHamiltonian"],
            Lookup[
                edgeAssoc,
                Key[edge],
                Lookup[edgeAssoc, Key[Reverse[edge]], globalValue]
            ]
        ]
    ];

edgeGExpression[g_?NumericQ, m_] := g;
edgeGExpression[g_Function, m_] := g[m];
edgeGExpression[_, _] := Missing["InvalidG"];

edgeSignedFlowValue[sys_?mfgSystemQ, edge_List, rules_List] :=
    Module[{signedFlows, expr},
        signedFlows = systemData[sys, "SignedFlows"];
        If[!AssociationQ[signedFlows], Return[Missing["MissingSignedFlows"], Module]];
        expr = Lookup[signedFlows, Key[edge], Missing["MissingSignedFlow"]];
        If[MissingQ[expr],
            expr = Lookup[signedFlows, Key[Reverse[edge]], Missing["MissingSignedFlow"]];
            If[MissingQ[expr], Return[expr, Module], expr = -expr]
        ];
        expr /. rules
    ];

edgeDensityValue[s_?scenarioQ, sys_?mfgSystemQ, edge_List, rules_List] :=
    Module[{ham, alpha, v, g, jval, m, gexpr, equation, roots},
        ham = scenarioData[s, "Hamiltonian"];
        If[!AssociationQ[ham], Return[Missing["InvalidHamiltonian"], Module]];
        alpha = edgeHamiltonianValue[ham, "Alpha", edge];
        v     = edgeHamiltonianValue[ham, "V",     edge];
        g     = edgeHamiltonianValue[ham, "G",     edge];
        If[MissingQ[alpha] || MissingQ[v] || MissingQ[g],
            Return[Missing["InvalidHamiltonian"], Module]
        ];
        jval  = edgeSignedFlowValue[sys, edge, rules];

        If[!NumericQ[jval], Return[Missing["MissingFlow"], Module]];
        If[Chop[N[jval]] == 0, Return[0, Module]];
        If[!NumericQ[alpha] || !NumericQ[v], Return[Missing["InvalidHamiltonian"], Module]];

        gexpr = edgeGExpression[g, m];
        If[MissingQ[gexpr], Return[gexpr, Module]];
        equation = N[jval]^2/(2 m^(2 - N[alpha])) - gexpr == -N[v];
        roots = Quiet @ Check[
            m /. NSolve[{equation, m > 0}, m, Reals],
            $Failed
        ];
        If[roots === $Failed || roots === {} || !ListQ[roots],
            Missing["DensityNotSolved"],
            FirstCase[N @ roots, value_?NumericQ /; value > 0,
                Missing["DensityNotSolved"]]
        ]
    ];

augmentAuxiliaryGraph[sys_?mfgSystemQ] :=
    Module[{auxPairs, triples, flowEdges, transitionEdges, allEdges, edgeVariables,
            edgeKinds, vertices},
        auxPairs = systemData[sys, "AuxPairs"];
        triples  = systemData[sys, "AuxTriples"];

        flowEdges = DirectedEdge[{#[[2]], #[[1]]}, {#[[1]], #[[2]]}] & /@ auxPairs;
        transitionEdges = DirectedEdge[{#[[1]], #[[2]]}, {#[[3]], #[[2]]}] & /@ triples;
        allEdges = Join[flowEdges, transitionEdges];
        vertices = DeleteDuplicates @ Flatten[List @@@ allEdges, 1];

        edgeVariables = Association @ Join[
            MapThread[#1 -> (j @@ #2) &, {flowEdges, auxPairs}],
            MapThread[#1 -> (j @@ #2) &, {transitionEdges, triples}]
        ];
        edgeKinds = Association @ Join[
            Thread[flowEdges -> "Flow"],
            Thread[transitionEdges -> "Transition"]
        ];

        <|
            "Graph" -> Graph[vertices, allEdges],
            "Vertices" -> vertices,
            "FlowEdges" -> flowEdges,
            "TransitionEdges" -> transitionEdges,
            "EdgeVariables" -> edgeVariables,
            "EdgeKinds" -> edgeKinds
        |>
    ];

scenarioTopologyPlot[s_?scenarioQ, sys_?mfgSystemQ, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    scenarioTopologyPlot[s, sys, PlotLabel -> title, opts];

scenarioTopologyPlot[s_?scenarioQ, sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    Module[{model, entryV, exitV, internalV, allV, edges},
        model    = scenarioData[s, "Model"];
        entryV   = First /@ model["Entries"];
        exitV    = First /@ model["Exits"];
        allV     = model["Vertices"];
        internalV = Complement[allV, entryV, exitV];
        edges    = systemData[sys, "Edges"];
        Graph[allV, edges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> Normal @ Join[
                AssociationThread[entryV,    RGBColor[0.22, 0.6, 0.3]],
                AssociationThread[exitV,     RGBColor[0.82, 0.27, 0.2]],
                AssociationThread[internalV, GrayLevel[0.75]]
            ],
            VertexSize -> 0.3,
            EdgeStyle  -> Directive[GrayLevel[0.5], AbsoluteThickness[2]],
            GraphLayout -> OptionValue[GraphLayout],
            PlotLabel  -> plotLabelValue[OptionValue[PlotLabel], "Network topology"],
            ImageSize  -> OptionValue[ImageSize]
        ]
    ];

mfgSolutionPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgSolutionPlot[s, sys, sol, PlotLabel -> title, opts];

mfgSolutionPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Show[
        GraphicsGrid[
            {{
                mfgFlowPlot[s, sys, sol,
                    GraphLayout -> OptionValue[GraphLayout],
                    PlotLabel -> "Flows",
                    ImageSize -> Medium],
                mfgValuePlot[s, sys, sol,
                    GraphLayout -> OptionValue[GraphLayout],
                    PlotLabel -> "Values",
                    ImageSize -> Medium],
                mfgDensityPlot[s, sys, sol,
                    GraphLayout -> OptionValue[GraphLayout],
                    PlotLabel -> "Agent density",
                    ImageSize -> Medium]
            }},
            ImageSize -> OptionValue[ImageSize]
        ],
        PlotLabel -> plotLabelValue[OptionValue[PlotLabel], "Solution views"],
        ImageSize -> OptionValue[ImageSize]
    ];

mfgValuePlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgValuePlot[s, sys, sol, PlotLabel -> title, opts];

mfgValuePlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{model, realV, edges, rules, edgeLabels, valueAt},
        model = scenarioData[s, "Model"];
        realV = model["Vertices"];
        edges = systemData[sys, "Edges"];
        rules = extractRules[sol];

        valueAt[a_, b_, t_] :=
            Module[{ua, ub},
                ua = u[a, b] /. rules;
                ub = u[b, a] /. rules;
                If[NumericQ[ua] && NumericQ[ub],
                    (1 - t) N[ua] + t N[ub],
                    Missing["MissingValue"]
                ]
            ];

        edgeLabels = Association @ Map[
            Module[{a, b, ua, ub, q1, q2, q3},
                {a, b} = List @@ #;
                ua = u[a, b] /. rules;
                ub = u[b, a] /. rules;
                q1 = valueAt[a, b, 1/4];
                q2 = valueAt[a, b, 1/2];
                q3 = valueAt[a, b, 3/4];
                # -> Placed[
                    Style[
                        Column[{
                            Row[{"u[", a, ",", b, "]=", formatPlotNumber[ua]}],
                            Row[{"1/4=", formatPlotNumber[q1],
                                 "  1/2=", formatPlotNumber[q2],
                                 "  3/4=", formatPlotNumber[q3]}],
                            Row[{"u[", b, ",", a, "]=", formatPlotNumber[ub]}]
                        }, Center, Spacings -> 0.15],
                        8, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            edges
        ];

        Graph[realV, edges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> realNetworkVertexStyle[s],
            VertexSize   -> Normal @ realNetworkVertexSize[s],
            EdgeStyle    -> Directive[RGBColor[0.36, 0.38, 0.42], AbsoluteThickness[2.2], Opacity[0.9]],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Values"],
            ImageSize    -> OptionValue[ImageSize]
        ]
    ];

mfgDensityPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgDensityPlot[s, sys, sol, PlotLabel -> title, opts];

mfgDensityPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{model, realV, edges, rules, densities, numericDensities, maxDensity,
            edgeStyles, edgeLabels},
        model = scenarioData[s, "Model"];
        realV = model["Vertices"];
        edges = systemData[sys, "Edges"];
        rules = extractRules[sol];

        densities = Association @ Map[
            With[{edge = List @@ #},
                # -> edgeDensityValue[s, sys, edge, rules]
            ] &,
            edges
        ];
        numericDensities = Cases[Values[densities], _?NumericQ, Infinity];
        maxDensity = Max[Append[Abs @ numericDensities, 1]];

        edgeStyles = Association @ Map[
            With[{density = densities[#]},
                # -> Directive[
                    If[NumericQ[density] && density > 0,
                        RGBColor[0.1, 0.5, 0.42],
                        GrayLevel[0.48]
                    ],
                    AbsoluteThickness[
                        If[NumericQ[density],
                            Rescale[Abs[density], {0, maxDensity}, {1.5, 7}],
                            2.0
                        ]
                    ],
                    Opacity[0.9]
                ]
            ] &,
            edges
        ];

        edgeLabels = Association @ Map[
            With[{density = densities[#]},
                # -> Placed[
                    Style[
                        Row[{"m=", formatPlotNumber[density]}],
                        9, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            edges
        ];

        Graph[realV, edges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> realNetworkVertexStyle[s],
            VertexSize   -> Normal @ realNetworkVertexSize[s],
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Agent density"],
            ImageSize    -> OptionValue[ImageSize]
        ]
    ];

mfgFlowPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgFlowPlot[s, sys, sol, PlotLabel -> title, opts];

mfgFlowPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{model, realV, auxV, edges, dispEdges, rules, edgeStyles, edgeLabels,
            auxEntryV, auxExitV},
        model     = scenarioData[s, "Model"];
        realV     = model["Vertices"];
        auxV      = Complement[systemData[sys, "AuxVertices"], realV];
        edges     = systemData[sys, "AuxEdges"];
        rules     = extractRules[sol];
        dispEdges = directedDisplayEdges[edges, rules];

        edgeStyles = Association @ Map[
            Module[{a, b, auxQ},
                {a, b} = List @@ #;
                auxQ = MemberQ[auxV, a] || MemberQ[auxV, b];
                # -> Directive[
                    If[auxQ, GrayLevel[0.35], RGBColor[0.12, 0.45, 0.78]],
                    AbsoluteThickness[2.0],
                    Opacity[0.9]
                ]
            ] &,
            dispEdges
        ];

        edgeLabels = Association @ Map[
            Module[{a, b, jv},
                {a, b} = List @@ #;
                jv = edgeJValue[a, b, rules];
                # -> Placed[
                    Style[
                        If[NumericQ[jv], NumberForm[N[jv], {5, 1}], "?"],
                        9, Black, Background -> White
                    ],
                    0.6
                ]
            ] &,
            dispEdges
        ];

        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];

        Graph[Join[realV, auxV], dispEdges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> Normal @ Join[
                Association @ realNetworkVertexStyle[s],
                AssociationThread[Intersection[auxV, auxEntryV], RGBColor[0.38, 0.74, 0.9]],
                AssociationThread[Intersection[auxV, auxExitV],  RGBColor[0.95, 0.7, 0.4]]
            ],
            VertexSize   -> Normal @ Join[
                AssociationThread[realV, 0.38],
                AssociationThread[auxV, 0.28]
            ],
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Flow graph"],
            ImageSize    -> OptionValue[ImageSize]
        ]
    ];

mfgTransitionPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgTransitionPlot[s, sys, sol, PlotLabel -> title, opts];

mfgTransitionPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{triples, rules, transitionEdges, nodeLabels, edgeStyles, edgeLabels,
            numericJ, maxJ, auxPairs},
        triples  = systemData[sys, "AuxTriples"];
        auxPairs = systemData[sys, "AuxPairs"];
        rules    = extractRules[sol];

        transitionEdges = DirectedEdge[#[[1 ;; 2]], #[[2 ;; 3]]] & /@ triples;

        numericJ = Cases[(j @@ # /. rules) & /@ triples, _?NumericQ, Infinity];
        maxJ = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            With[{trip = {#[[1, 1]], #[[1, 2]], #[[2, 2]]}, jv = (j @@ {#[[1, 1]], #[[1, 2]], #[[2, 2]]}) /. rules},
                # -> Directive[
                    Which[
                        NumericQ[jv] && jv > 0, RGBColor[0.12, 0.45, 0.78],
                        NumericQ[jv] && jv < 0, RGBColor[0.82, 0.39, 0.2],
                        True, GrayLevel[0.45]
                    ],
                    AbsoluteThickness[
                        If[NumericQ[jv], Rescale[Abs[jv], {0, maxJ}, {1.5, 7}], 2.0]
                    ],
                    Opacity[0.9]
                ]
            ] &,
            transitionEdges
        ];

        edgeLabels = Association @ Map[
            With[{trip = {#[[1, 1]], #[[1, 2]], #[[2, 2]]}, jv = (j @@ {#[[1, 1]], #[[1, 2]], #[[2, 2]]}) /. rules},
                # -> Placed[
                    Style[
                        If[NumericQ[jv], NumberForm[N[jv], {5, 1}], "?"],
                        9, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            transitionEdges
        ];

        nodeLabels = Association @ Map[
            With[{uv = edgeUValue[#[[1]], #[[2]], rules]},
                # -> Placed[
                    Style[
                        If[NumericQ[uv],
                            Row[{stateLabel[#], "\n", "u=", NumberForm[N[uv], {5, 1}]}],
                            stateLabel[#]
                        ],
                        9, Black, Background -> RGBColor[0.96, 0.96, 0.96]
                    ],
                    Center
                ]
            ] &,
            auxPairs
        ];

        Graph[auxPairs, transitionEdges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.55,
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Transition graph: flows and values"],
            ImageSize    -> OptionValue[ImageSize]
        ]
    ];

mfgAugmentedPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgAugmentedPlot[s, sys, sol, PlotLabel -> title, opts];

mfgAugmentedPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{augmented, vertices, flowEdges, transitionEdges, allAugEdges,
            edgeStyles, edgeLabels, nodeLabels, numericJ, maxJ,
            edgeVariables, edgeKinds, rules, getU, auxEntryV, auxExitV},

        augmented = augmentAuxiliaryGraph[sys];
        vertices = augmented["Vertices"];
        flowEdges = augmented["FlowEdges"];
        transitionEdges = augmented["TransitionEdges"];
        allAugEdges = Join[flowEdges, transitionEdges];
        edgeVariables = augmented["EdgeVariables"];
        edgeKinds = augmented["EdgeKinds"];
        rules    = extractRules[sol];
        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];

        getU[p_] := (u @@ p) /. rules;

        numericJ = Cases[Lookup[edgeVariables, allAugEdges] /. rules, _?NumericQ, Infinity];
        maxJ = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            With[{jv = edgeVariables[#] /. rules, kind = edgeKinds[#]},
                # -> Directive[
                    Which[
                        kind === "Flow", RGBColor[0.12, 0.45, 0.78],
                        kind === "Transition", RGBColor[0.86, 0.25, 0.22],
                        True, GrayLevel[0.45]
                    ],
                    AbsoluteThickness[
                        If[NumericQ[jv], Rescale[Abs[jv], {0, maxJ}, {1.5, 7}], 2.0]
                    ],
                    Opacity[0.8]
                ]
            ] &,
            allAugEdges
        ];

        edgeLabels = Association @ Map[
            With[{jv = edgeVariables[#] /. rules},
                # -> Placed[
                    Style[
                        If[NumericQ[jv] && Abs[jv] > 10^-5, NumberForm[N[jv], {5, 1}], ""],
                        8, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            allAugEdges
        ];

        nodeLabels = Association @ Map[
            With[{uv = getU[#]},
                # -> Placed[
                    Style[
                        If[NumericQ[uv],
                            Row[{stateLabel[#], "\n", NumberForm[N[uv], {5, 2}]}],
                            stateLabel[#]
                        ],
                        9, Black
                    ],
                    Center
                ]
            ] &,
            vertices
        ];

        Graph[vertices, allAugEdges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.58,
            VertexStyle  -> Normal @ Join[
                AssociationThread[Select[vertices, Intersection[#, auxEntryV] =!= {} &], RGBColor[0.38, 0.74, 0.9]],
                AssociationThread[Select[vertices, Intersection[#, auxExitV]  =!= {} &], RGBColor[0.95, 0.7, 0.4]]
            ],
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Augmented infrastructure graph (Paper scheme)"],
            ImageSize    -> OptionValue[ImageSize]
        ]
    ];

End[];

EndPackage[];
