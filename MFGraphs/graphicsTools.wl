(* Wolfram Language package *)
(* Graphics: public visualization helpers for MFGraphs systems and solutions. *)

BeginPackage["graphicsTools`", {"primitives`", "utilities`", "scenarioTools`", "systemTools`"}];

scenarioTopologyPlot::usage =
"scenarioTopologyPlot[s, sys, opts] plots the scenario topology using vertex coloring for entry (green), exit (red), and internal (gray) vertices. Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Medium).";

mfgSolutionPlot::usage =
"mfgSolutionPlot[s, sys, sol, opts] plots coordinated solution views: flows, value samples, and agent density. Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Large).";

mfgFlowPlot::usage =
"mfgFlowPlot[s, sys, sol, opts] plots a flow-only solution graph with real and auxiliary directed edges; edge labels show j flow values. Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Large).";

mfgValuePlot::usage =
"mfgValuePlot[s, sys, sol, opts] plots real network edges with endpoint value variables u[a,b], u[b,a] and linearly interpolated interior samples at 1/4, 1/2, and 3/4. Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Large).";

mfgDensityPlot::usage =
"mfgDensityPlot[s, sys, sol, opts] plots real network edges with agent density inferred from solved signed spatial flow and the scenario Hamiltonian density equation. Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Large).";

mfgTransitionPlot::usage =
"mfgTransitionPlot[s, sys, sol, opts] plots the transition graph. Nodes are AuxPair states {r,i}; directed edges represent transition flows j[r,i,w] labeled with their solved value. Nodes are labeled with internal values u where available. Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Large).";

mfgAugmentedPlot::usage =
"mfgAugmentedPlot[s, sys, sol, opts] plots the augmented infrastructure graph. Blue arcs are flow variables j[a,b]; red arcs are transition variables j[r,i,w]. Anti-parallel flow arcs curve to opposite sides so both directions are visible. Node colors show u-values on a gradient when a solution is provided; zero-flow arcs are invisible (Opacity 0). Options: PlotLabel (default Automatic), GraphLayout (default Automatic), ImageSize (default Large), ColorFunction (default Automatic, a Red\[Rule]Blue blend applied to u-values), ShowBoundaryValues (default True, shows entry/exit u-values on boundary nodes), ShowLegend (default True, shows color bar when u-values are numeric), BendFactor (default 0.15, controls arc curvature as a fraction of edge length; 0 gives straight edges).";

mfgAugmentedBoundaryPlot::usage =
"mfgAugmentedBoundaryPlot[s, sys, sol, opts] plots the augmented infrastructure graph with boundary data overlaid. In the augmented representation, edges correspond to flow and transition flow variables; vertices correspond to value variables. Entry flow edges are labeled with their fixed supply from boundary data; exit source vertices {exN,N} are annotated with their exit cost bound. Use PlotLabel, GraphLayout, and ImageSize options to control display.";

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

Options[mfgAugmentedPlot] = Join[Options[mfgTransitionPlot], {ColorFunction -> Automatic, ShowBoundaryValues -> True, ShowLegend -> True, BendFactor -> 0.15}];

Options[mfgAugmentedBoundaryPlot] = Options[mfgTransitionPlot];

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

stateAtomLabel::usage =
"stateAtomLabel[x] returns a compact string representation of an auxiliary or real vertex index.";

edgeGExpression::usage =
"edgeGExpression[g, m] evaluates the Hamiltonian congestion cost term G at mass m.";

edgeSignedFlowValue::usage =
"edgeSignedFlowValue[sys, edge, rules] returns the solved signed spatial flow for a real network edge.";

stateVertexStyle::usage =
"stateVertexStyle[vertices, sys] returns plot styles for augmented transition graph vertices.";

formatStateNodeLabel::usage =
"formatStateNodeLabel[state, val] formats a node label for the transition graph showing state and value.";

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
    stateAtomLabel[b];

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

stateVertexStyle[vertices_, sys_?mfgSystemQ] :=
    Module[{auxEntryV, auxExitV},
        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];
        Normal @ Join[
            AssociationThread[vertices, GrayLevel[0.72]],
            AssociationThread[Select[vertices, Intersection[#, auxEntryV] =!= {} &], RGBColor[0.22, 0.6, 0.3]],
            AssociationThread[Select[vertices, Intersection[#, auxExitV]  =!= {} &], RGBColor[0.82, 0.27, 0.2]]
        ]
    ];

formatStateNodeLabel[state_, val_] :=
    Placed[
        Style[
            If[NumericQ[val],
                Column[{stateLabel[state], formatPlotNumber[val]}, Center, Spacings -> 0],
                stateLabel[state]
            ],
            9, Black
        ],
        Center
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
    Module[{triples, rules, transitionEdges, vertices, nodeLabels, edgeStyles, edgeLabels,
            numericJ, maxJ},
        triples = systemData[sys, "AuxTriples"];
        rules   = extractRules[sol];

        transitionEdges = DirectedEdge[{#[[1]], #[[2]]}, {#[[3]], #[[2]]}] & /@ triples;
        vertices = DeleteDuplicates @ Flatten[List @@@ transitionEdges, 1];

        numericJ = Cases[(j @@ # /. rules) & /@ triples, _?NumericQ, Infinity];
        maxJ = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            With[{jv = (j @@ {#[[1, 1]], #[[1, 2]], #[[2, 1]]}) /. rules},
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
            With[{jv = (j @@ {#[[1, 1]], #[[1, 2]], #[[2, 1]]}) /. rules},
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
            # -> formatStateNodeLabel[#, (u @@ #) /. rules] &,
            vertices
        ];

        Graph[vertices, transitionEdges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.55,
            VertexStyle  -> stateVertexStyle[vertices, sys],
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
            edgeLabels, nodeLabels, numericJ, maxJ,
            edgeVariables, edgeKinds, rules, getU, getEffectiveU, auxEntryV, auxExitV,
            exitCosts, effectiveFlows, uValues, numericU, uMin, uMax, vertexStyle,
            colorFn, nSeg = 20, legend, graph},

        augmented = augmentAuxiliaryGraph[sys];
        vertices = augmented["Vertices"];
        flowEdges = augmented["FlowEdges"];
        transitionEdges = augmented["TransitionEdges"];
        allAugEdges = Join[flowEdges, transitionEdges];
        edgeVariables = augmented["EdgeVariables"];
        edgeKinds = augmented["EdgeKinds"];
        rules     = extractRules[sol];
        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];
        exitCosts = systemData[sys, "ExitCosts"];

        getU[p_] := (u @@ p) /. rules;
        (* For aux exit vertices the solution may not include u — fall back to the exit cost.
           Exit vertex label can sit at either position in the pair. *)
        getEffectiveU[p_] := With[{uv = getU[p]},
            If[NumericQ[uv], uv,
                If[TrueQ[OptionValue[ShowBoundaryValues]],
                    With[{cost = Lookup[exitCosts, p[[1]], Lookup[exitCosts, p[[2]], Missing[]]]},
                        If[!MissingQ[cost], cost,
                            If[MemberQ[auxEntryV, p[[1]]] || MemberQ[auxEntryV, p[[2]]],
                                With[{adjU = getU[Reverse[p]]},
                                    If[NumericQ[adjU], adjU, Missing[]]
                                ],
                                Missing[]
                            ]
                        ]
                    ],
                    Missing[]
                ]
            ]
        ];

        effectiveFlows = Association @ Map[
            With[{v1 = If[rules =!= {},
                edgeVariables[#] /. rules /. _j -> 0,
                edgeVariables[#]
            ]},
                # -> v1
            ] &,
            allAugEdges
        ];

        numericJ = Cases[Values[effectiveFlows], _?NumericQ, Infinity];
        maxJ = Max[Append[numericJ, 1]];

        edgeLabels = Association @ Map[
            With[{jv = edgeVariables[#] /. rules},
                # -> Placed[
                    Style[
                        If[NumericQ[jv] && jv != 0, NumberForm[N[jv], {5, 1}], ""],
                        8, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            allAugEdges
        ];

        nodeLabels = Association @ Map[
            # -> formatStateNodeLabel[#, getEffectiveU[#]] &,
            vertices
        ];

        colorFn = Replace[OptionValue[ColorFunction],
            Automatic -> (Blend[{Red, Blue}, #] &)];

        uValues  = AssociationMap[getEffectiveU, vertices];
        numericU = Select[uValues, NumericQ];
        If[Length[numericU] > 0,
            uMin = Min[Values[numericU]];
            uMax = Max[Values[numericU]];
            vertexStyle = Normal @ Map[
                If[NumericQ[#],
                    colorFn[If[uMin == uMax, 0.5, Rescale[#, {uMin, uMax}]]],
                    GrayLevel[0.72]
                ] &,
                uValues
            ];
            legend = BarLegend[{colorFn[Rescale[#, {uMin, uMax}]] &, {uMin, uMax}},
                LegendLabel -> Placed["u", Right],
                LegendMarkerSize -> 200],
            vertexStyle = stateVertexStyle[vertices, sys];
            legend = None
        ];

        graph = Graph[vertices, allAugEdges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.58,
            VertexStyle  -> vertexStyle,
            EdgeShapeFunction -> Function[{pts, e},
                With[{
                    uS   = uValues[e[[1]]],
                    uT   = uValues[e[[2]]],
                    effJ = effectiveFlows[e],
                    kind = edgeKinds[e],
                    p1   = N[pts[[1]]],
                    p2   = N[pts[[-1]]]
                },
                    Module[{
                        op    = If[NumericQ[effJ] && effJ == 0, 0, 0.9],
                        thick = AbsoluteThickness[If[NumericQ[effJ], Rescale[effJ, {0, maxJ}, {1.5, 7}], 2.0]],
                        edgeColorAt, d, perp, ctrl, pAt
                    },
                        d    = p2 - p1;
                        perp = With[{len = Norm[d]}, If[len > 0, {-d[[2]], d[[1]]} / len, {0, 1}]];
                        ctrl = 0.5*(p1 + p2) + OptionValue[BendFactor]*Norm[d]*perp;
                        (* Quadratic Bezier curves each edge to the left of its direction;
                           anti-parallel pairs curve to opposite sides and appear as distinct arcs. *)
                        pAt  = Function[t, (1-t)^2*p1 + 2*(1-t)*t*ctrl + t^2*p2];
                        edgeColorAt = If[NumericQ[uS] && NumericQ[uT] && Length[numericU] > 0,
                            colorFn[Rescale[(1 - #)*uS + #*uT, {uMin, uMax}]] &,
                            (Which[kind === "Flow", RGBColor[0.12, 0.45, 0.78],
                                   kind === "Transition", RGBColor[0.86, 0.25, 0.22],
                                   True, GrayLevel[0.45]]) &
                        ];
                        Join[
                            Table[
                                With[{t0 = k/nSeg, t1 = (k + 1)/nSeg},
                                    {edgeColorAt[t0 + 0.5/nSeg], thick, Opacity[op],
                                     Line[{pAt[t0], pAt[t1]}]}],
                                {k, 0, nSeg - 1}],
                            {{GrayLevel[0.3], Opacity[op], Arrowheads[{{0.02, 1}}],
                              Arrow[{pAt[0.45], pAt[0.55]}]}}
                        ]
                    ]
                ]
            ],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Augmented infrastructure graph (Paper scheme)"],
            ImageSize    -> OptionValue[ImageSize]
        ];

        If[legend === None || !TrueQ[OptionValue[ShowLegend]], graph, Legended[graph, legend]]
    ];

entryFlowEdge[pair_List] :=
    DirectedEdge[{pair[[2]], pair[[1]]}, {pair[[1]], pair[[2]]}];

formatBoundaryExitLabel[state_, val_, cost_] :=
    Placed[
        Style[
            If[NumericQ[val] && Abs[val] > 10^-5,
                Column[{
                    Row[{stateLabel[state], "  u=", formatPlotNumber[val]}],
                    Style[Row[{"c\[LessEqual]", formatPlotNumber[cost]}], Italic]
                }, Center, Spacings -> 0.1],
                Column[{
                    stateLabel[state],
                    Style[Row[{"c\[LessEqual]", formatPlotNumber[cost]}], Italic]
                }, Center, Spacings -> 0.1]
            ],
            9, Black
        ],
        Center
    ];

mfgAugmentedBoundaryPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    mfgAugmentedBoundaryPlot[s, sys, sol, PlotLabel -> title, opts];

mfgAugmentedBoundaryPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{augmented, vertices, flowEdges, transitionEdges, allAugEdges,
            edgeStyles, edgeLabels, nodeLabels, numericJ, maxJ,
            edgeVariables, edgeKinds, rules, getU, auxEntryV, auxExitV,
            exitCosts, entryDataAssoc, inAuxEntryPairs, entryEdgeSet},

        augmented = augmentAuxiliaryGraph[sys];
        vertices        = augmented["Vertices"];
        flowEdges       = augmented["FlowEdges"];
        transitionEdges = augmented["TransitionEdges"];
        allAugEdges     = Join[flowEdges, transitionEdges];
        edgeVariables   = augmented["EdgeVariables"];
        edgeKinds       = augmented["EdgeKinds"];
        rules           = extractRules[sol];
        auxEntryV       = systemData[sys, "AuxEntryVertices"];
        auxExitV        = systemData[sys, "AuxExitVertices"];
        exitCosts       = systemData[sys, "ExitCosts"];
        entryDataAssoc  = systemData[sys, "EntryDataAssociation"];
        inAuxEntryPairs = systemData[sys, "InAuxEntryPairs"];
        entryEdgeSet    = Association @ Thread[(entryFlowEdge /@ inAuxEntryPairs) -> inAuxEntryPairs];

        getU[p_] := (u @@ p) /. rules;

        numericJ = Cases[Lookup[edgeVariables, allAugEdges] /. rules, _?NumericQ, Infinity];
        maxJ = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            With[{jv = EchoLabel["jv: "][edgeVariables[#] /. rules /. _j -> 0], kind = edgeKinds[#]},
                Module[{op = If[NumericQ[jv] && jv == 0, 0, 0.8]},
                    EchoLabel["op: "][op];
                    # -> Directive[
                        Which[
                            kind === "Flow", RGBColor[0.12, 0.45, 0.78],
                            kind === "Transition", RGBColor[0.86, 0.25, 0.22],
                            True, GrayLevel[0.45]
                        ],
                        AbsoluteThickness[
                            If[NumericQ[jv], Rescale[Abs[jv], {0, maxJ}, {1.5, 7}], 2.0]
                        ],
                        Opacity[op]
                    ]
                ]
            ] &,
            allAugEdges
        ];

        edgeLabels = Association @ Map[
            With[{jv = (edgeVariables[#] /. rules /. _j -> 0),
                  entryPair = Lookup[entryEdgeSet, Key[#], Missing[]]},
                # -> Placed[
                    Style[
                        Which[
                            !MissingQ[entryPair],
                                With[{supply = Lookup[entryDataAssoc, Key[entryPair], Missing[]]},
                                    If[!MissingQ[supply],
                                        Row[{"f=", formatPlotNumber[supply]}],
                                        If[NumericQ[jv] && jv != 0, NumberForm[N[jv], {5, 1}], ""]
                                    ]
                                ],
                            NumericQ[jv] && jv != 0,
                                NumberForm[N[jv], {5, 1}],
                            True, ""
                        ],
                        8, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            allAugEdges
        ];

        nodeLabels = Association @ Map[
            With[{uval = getU[#], head = First[#]},
                # -> If[MemberQ[auxExitV, head],
                    formatBoundaryExitLabel[#, uval, Lookup[exitCosts, head, Missing[]]],
                    formatStateNodeLabel[#, uval]
                ]
            ] &,
            vertices
        ];

        Graph[vertices, allAugEdges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.58,
            VertexStyle  -> stateVertexStyle[vertices, sys],
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Augmented graph with boundary data"],
            ImageSize    -> OptionValue[ImageSize]
        ]
    ];

End[];

EndPackage[];
