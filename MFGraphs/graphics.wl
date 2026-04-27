(* Wolfram Language package *)
(* Graphics: public visualization helpers for MFGraphs systems and solutions. *)

BeginPackage["graphics`", {"primitives`", "scenarioTools`", "systemTools`"}];

scenarioTopologyPlot::usage =
"scenarioTopologyPlot[s, sys] plots the scenario topology using vertex coloring for entry, exit, and internal vertices. An optional third argument sets the title.";

mfgSolutionPlot::usage =
"mfgSolutionPlot[s, sys, sol] plots a combined solution graph with real and auxiliary edges. Edge labels include both j and u values, and real-edge directions follow solved net flow orientation. An optional fourth argument sets the title.";

mfgTransitionPlot::usage =
"mfgTransitionPlot[s, sys, sol] plots the transition graph of the solution. Nodes represent network edges (or entry/exit), and directed edges represent transition flows (j[a, b, c]). Nodes are labeled with their internal values (u).";

mfgAugmentedPlot::usage =
"mfgAugmentedPlot[s, sys, sol] plots the augmented infrastructure graph (Paper scheme). Nodes represent edge-vertex pairs (e, v), and edges represent both flows within edges and transitions at vertices. Value functions (u) are vertex labels, and transition flows (j) are edge labels.";

Begin["`Private`"];

extractRules::usage =
"extractRules[sol] extracts replacement rules from a solution list or solution association.";

netEdgeFlow::usage =
"netEdgeFlow[a, b, rules] returns the solved net flow j[a,b]-j[b,a], defaulting absent flow variables to zero.";

edgeJValue::usage =
"edgeJValue[a, b, rules] returns the solved flow j[a,b], defaulting absent flow variables to zero.";

edgeUValue::usage =
"edgeUValue[a, b, rules] returns the solved value variable for edge endpoints, or Missing[] when absent.";

directedDisplayEdges::usage =
"directedDisplayEdges[edges, rules] orients display edges using solved net flow values.";

extractRules[sol_List] := Select[sol, MatchQ[#, _Rule | _RuleDelayed] &];
extractRules[sol_Association] := extractRules[Lookup[sol, "Rules", {}]];
extractRules[_] := {};

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

scenarioTopologyPlot[s_?scenarioQ, sys_?mfgSystemQ, title_: Automatic] :=
    Module[{model, entryV, exitV, internalV, allV, edges, plotTitle},
        model    = scenarioData[s, "Model"];
        entryV   = First /@ model["Entries"];
        exitV    = First /@ model["Exits"];
        allV     = model["Vertices"];
        internalV = Complement[allV, entryV, exitV];
        edges    = systemData[sys, "Edges"];
        plotTitle = Replace[title, Automatic -> "Network topology"];
        Graph[allV, edges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> Normal @ Join[
                AssociationThread[entryV,    RGBColor[0.22, 0.6, 0.3]],
                AssociationThread[exitV,     RGBColor[0.82, 0.27, 0.2]],
                AssociationThread[internalV, GrayLevel[0.75]]
            ],
            VertexSize -> 0.3,
            EdgeStyle  -> Directive[GrayLevel[0.5], AbsoluteThickness[2]],
            GraphLayout -> "LayeredDigraphEmbedding",
            PlotLabel  -> Style[plotTitle, 14, Bold],
            ImageSize  -> Medium
        ]
    ];

mfgSolutionPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title_: Automatic] :=
    Module[{model, entryV, exitV, realV, auxV, internalV, edges, dispEdges, rules,
            edgeStyles, edgeLabels, plotTitle, numericJ, maxJ},
        model     = scenarioData[s, "Model"];
        entryV    = First /@ model["Entries"];
        exitV     = First /@ model["Exits"];
        realV     = model["Vertices"];
        auxV      = Complement[systemData[sys, "AuxVertices"], realV];
        internalV = Complement[realV, entryV, exitV];
        edges     = systemData[sys, "AuxEdges"];
        rules     = extractRules[sol];
        dispEdges = directedDisplayEdges[edges, rules];
        numericJ  = Cases[edgeJValue @@@ (List @@@ dispEdges), _?NumericQ, Infinity];
        maxJ      = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            Module[{a, b, jv, auxQ},
                {a, b} = List @@ #;
                jv = edgeJValue[a, b, rules];
                auxQ = MemberQ[auxV, a] || MemberQ[auxV, b];
                # -> Directive[
                    Which[
                        auxQ, RGBColor[0.35, 0.35, 0.35],
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
            dispEdges
        ];

        edgeLabels = Association @ Map[
            Module[{a, b, jv, uv},
                {a, b} = List @@ #;
                jv = edgeJValue[a, b, rules];
                uv = edgeUValue[a, b, rules];
                # -> Placed[
                    Style[
                        Row[{
                            "j=", If[NumericQ[jv], NumberForm[N[jv], {5, 1}], "?"],
                            "  |  u=", If[NumericQ[uv], NumberForm[N[uv], {5, 1}], "?"]
                        }],
                        9, Black, Background -> White
                    ],
                    Center
                ]
            ] &,
            dispEdges
        ];

        plotTitle = Replace[title, Automatic -> "Solution graph: flows and values"];
        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];

        Graph[Join[realV, auxV], dispEdges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> Normal @ Join[
                AssociationThread[entryV,    RGBColor[0.22, 0.6, 0.3]],
                AssociationThread[exitV,     RGBColor[0.82, 0.27, 0.2]],
                AssociationThread[internalV, GrayLevel[0.72]],
                AssociationThread[Intersection[auxV, auxEntryV], RGBColor[0.38, 0.74, 0.9]],
                AssociationThread[Intersection[auxV, auxExitV],  RGBColor[0.95, 0.7, 0.4]]
            ],
            VertexSize   -> Normal @ Join[
                AssociationThread[realV, 0.38],
                AssociationThread[auxV, 0.28]
            ],
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> "LayeredDigraphEmbedding",
            PlotLabel    -> Style[plotTitle, 14, Bold],
            ImageSize    -> Large
        ]
    ];

mfgTransitionPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title_: Automatic] :=
    Module[{triples, rules, transitionEdges, nodeLabels, edgeStyles, edgeLabels,
            numericJ, maxJ, plotTitle, auxPairs},
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
                        Row[{
                            ToString[#],
                            " (u=", If[NumericQ[uv], NumberForm[N[uv], {5, 1}], "?"], ")"
                        }],
                        10, Black, Background -> RGBColor[0.95, 0.95, 0.95]
                    ],
                    Center
                ]
            ] &,
            auxPairs
        ];

        plotTitle = Replace[title, Automatic -> "Transition graph: flows and values"];
        Graph[auxPairs, transitionEdges,
            VertexShapeFunction -> "RoundRectangle",
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> {"Scaled", 0.15},
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> "LayeredDigraphEmbedding",
            PlotLabel    -> Style[plotTitle, 14, Bold],
            ImageSize    -> Large
        ]
    ];

mfgAugmentedPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title_: Automatic] :=
    Module[{auxPairs, triples, rules, flowEdges, transitionEdges, allAugEdges,
            edgeStyles, edgeLabels, nodeLabels, plotTitle, numericJ, maxJ,
            getJ, getU, auxEntryV, auxExitV},

        auxPairs = systemData[sys, "AuxPairs"];
        triples  = systemData[sys, "AuxTriples"];
        rules    = extractRules[sol];
        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];

        (* 1. Flow edges: {u, v} -> {v, u} (magnitude j[v, u]) *)
        flowEdges = Map[DirectedEdge[#, Reverse[#]] &, auxPairs];

        (* 2. Transition edges: {i, r} -> {i, w} (magnitude j[i, r, w]) *)
        transitionEdges = DirectedEdge[{#[[2]], #[[1]]}, {#[[2]], #[[3]]}] & /@ triples;

        allAugEdges = Join[flowEdges, transitionEdges];

        getJ[DirectedEdge[p1_, p2_]] :=
            If[p1 === Reverse[p2],
                j @@ p2, (* Flow within edge: magnitude j[target, source] *)
                j[p1[[1]], p1[[2]], p2[[2]]] (* Transition at vertex: magnitude j[v, e1, e2] *)
            ];

        getU[p_] := (u @@ p) /. rules;

        numericJ = Cases[getJ /@ allAugEdges /. rules, _?NumericQ, Infinity];
        maxJ = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            With[{jv = getJ[#] /. rules},
                # -> Directive[
                    Which[
                        NumericQ[jv] && jv > 0, RGBColor[0.12, 0.45, 0.78],
                        NumericQ[jv] && jv < 0, RGBColor[0.82, 0.39, 0.2],
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
            With[{jv = getJ[#] /. rules},
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
                        Row[{
                            ToString[#],
                            "\n",
                            If[NumericQ[uv], NumberForm[N[uv], {5, 2}], "?"]
                        }],
                        9, Black
                    ],
                    Center
                ]
            ] &,
            auxPairs
        ];

        plotTitle = Replace[title, Automatic -> "Augmented infrastructure graph (Paper scheme)"];
        Graph[auxPairs, allAugEdges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.4,
            VertexStyle  -> Normal @ Join[
                AssociationThread[Select[auxPairs, MemberQ[auxEntryV, #[[1]]] &], RGBColor[0.38, 0.74, 0.9]],
                AssociationThread[Select[auxPairs, MemberQ[auxExitV,  #[[1]]] &], RGBColor[0.95, 0.7, 0.4]]
            ],
            EdgeStyle    -> Normal[edgeStyles],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> "LayeredDigraphEmbedding",
            PlotLabel    -> Style[plotTitle, 14, Bold],
            ImageSize    -> Large
        ]
    ];

End[];

EndPackage[];
