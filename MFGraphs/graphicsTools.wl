(* Wolfram Language package *)
(* Graphics: public visualization helpers for MFGraphs systems and solutions. *)

BeginPackage["graphicsTools`", {"primitives`", "utilities`", "scenarioTools`", "systemTools`"}];

rawNetworkPlot::usage =
"rawNetworkPlot[s, sys, opts] and rawNetworkPlot[s, sys, sol, opts] render the physical network. Real network vertices are gray. Auxiliary entry/exit vertices and boundary edges are hidden by default. Options: ShowAuxiliaryVertices (default False), ShowBoundaryData (default False; when True forces ShowAuxiliaryVertices), ShowBoundaryValues (default True; shows u-values on auxiliary exit vertices when boundary data is shown), ShowFlowLabels (default Automatic; True when sol provided), ShowValueLabels (default Automatic; True when sol provided; shows u-values with 1/4, 1/2, 3/4 interpolations), ShowDensityLabels (default False; shows inferred density m), ColorFunction (default Automatic; Blue\[Rule]Red blend over u-values), ShowLegend (default True), GraphLayout (default Automatic), PlotLabel (default Automatic), ImageSize (default Large). With a solution provided, auxiliary vertices (when shown) and shown vertex states are colored by u-value gradient; physical entry/exit vertices remain gray.";

richNetworkPlot::usage =
"richNetworkPlot[s, sys, opts] and richNetworkPlot[s, sys, sol, opts] render the augmented state-space graph. Nodes are pairs {a,b} representing oriented edge states; edges are flow arcs j[a,b] (blue) and transition arcs j[r,i,w] (red), drawn as quadratic Bezier curves so anti-parallel pairs separate. Only nodes whose label position (b) is an auxiliary entry/exit vertex receive boundary colors. Options: ShowFlowEdges (default True; False yields the transition-only graph), ShowBoundaryData (default False; overlays entry-flow and exit-cost labels), ShowBoundaryValues (default True; shows u/cost on boundary nodes), ShowFlowLabels (default Automatic; True when sol provided), ShowValueLabels (default Automatic; True when sol provided), UseColorFunction (default False; when True colors nodes and eligible edges by u-values using ColorFunction), ColorFunction (default Automatic), ShowLegend (default True), BendFactor (default 0.15; arc curvature as fraction of edge length), GraphLayout (default Automatic), PlotLabel (default Automatic), ImageSize (default Large).";

augmentAuxiliaryGraph::usage =
"augmentAuxiliaryGraph[sys] constructs the road-traffic augmented infrastructure graph from a system's AuxPairs and AuxTriples. Returns an Association containing the Graph, Vertices, FlowEdges, TransitionEdges, EdgeVariables, and EdgeKinds.";

ShowAuxiliaryVertices::usage = "ShowAuxiliaryVertices is an option for rawNetworkPlot.";
ShowBoundaryData::usage      = "ShowBoundaryData is an option for rawNetworkPlot and richNetworkPlot.";
ShowBoundaryValues::usage    = "ShowBoundaryValues is an option for rawNetworkPlot and richNetworkPlot.";
ShowFlowEdges::usage         = "ShowFlowEdges is an option for richNetworkPlot.";
ShowFlowLabels::usage        = "ShowFlowLabels is an option for rawNetworkPlot and richNetworkPlot.";
ShowValueLabels::usage       = "ShowValueLabels is an option for rawNetworkPlot and richNetworkPlot.";
ShowDensityLabels::usage     = "ShowDensityLabels is an option for rawNetworkPlot.";
ShowLegend::usage            = "ShowLegend is an option for rawNetworkPlot and richNetworkPlot.";
UseColorFunction::usage      = "UseColorFunction is an option for richNetworkPlot.";
BendFactor::usage            = "BendFactor is an option for richNetworkPlot.";

Begin["`Private`"];

(* ===========================================================
   Options
   =========================================================== *)

Options[rawNetworkPlot] = {
    ShowAuxiliaryVertices -> False,
    ShowBoundaryData      -> False,
    ShowBoundaryValues    -> True,
    ShowFlowLabels        -> Automatic,
    ShowValueLabels       -> Automatic,
    ShowDensityLabels     -> False,
    ColorFunction         -> Automatic,
    ShowLegend            -> True,
    GraphLayout           -> Automatic,
    PlotLabel             -> Automatic,
    ImageSize             -> Large
};

Options[richNetworkPlot] = {
    ShowFlowEdges         -> True,
    ShowBoundaryData      -> False,
    ShowBoundaryValues    -> True,
    ShowFlowLabels        -> Automatic,
    ShowValueLabels       -> Automatic,
    UseColorFunction      -> False,
    ColorFunction         -> Automatic,
    ShowLegend            -> True,
    BendFactor            -> 0.15,
    GraphLayout           -> Automatic,
    PlotLabel             -> Automatic,
    ImageSize             -> Large
};

(* ===========================================================
   Small formatting helpers
   =========================================================== *)

formatPlotNumber[value_?NumericQ] := NumberForm[N[value], {5, 2}];
formatPlotNumber[_] := "?";

stateAtomLabel[x_] :=
    StringReplace[ToString[x], {"auxEntry" -> "in", "auxExit" -> "out"}];

stateLabel[{a_, b_}] := stateAtomLabel[b];

formatRawVertexLabel[v_] :=
    Placed[Style[stateAtomLabel[v], 9, Black], Center];

formatRawBoundaryExitVertexLabel[v_, val_] :=
    Placed[
        Style[
            If[NumericQ[val],
                Column[{stateAtomLabel[v], Row[{"u=", formatPlotNumber[val]}]},
                       Center, Spacings -> 0],
                stateAtomLabel[v]
            ],
            9, Black
        ],
        Center
    ];

plotLabelValue[label_, default_String] :=
    Replace[label, {
        Automatic -> Style[default, 14, Bold],
        None -> None,
        other_ :> Style[other, 14, Bold]
    }];

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

(* Solution rules — accepts {rules}, <|"Rules" -> rules, ...|>, or <||> *)
extractSolutionRules[sol_] := extractRules[sol];

resolveAutomaticOption[val_, autoCondition_] :=
    Replace[val, {Automatic :> autoCondition, x_ :> TrueQ[x]}];

(* ===========================================================
   Edge-data helpers
   =========================================================== *)

netEdgeFlow[a_, b_, rules_List] :=
    (j[a, b] - j[b, a]) /. rules /. {_j -> 0};

edgeJValue[a_, b_, rules_List] :=
    j[a, b] /. rules /. {_j -> 0};

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

(* Force aux entry edges to point auxEntry -> N and aux exit edges to point
   N -> auxExit, regardless of how the edge currently appears (UndirectedEdge
   or either DirectedEdge orientation). Internal edges are left untouched. *)
forceAuxEdgeDirections[edges_List, inAuxEntryPairs_List, outAuxExitPairs_List] :=
    Module[{entryKeys, exitKeys, canonical},
        entryKeys = Association @ Map[(Sort[#] -> #) &, inAuxEntryPairs];
        exitKeys  = Association @ Map[(Sort[#] -> #) &, outAuxExitPairs];
        canonical[a_, b_] := Module[{key = Sort[{a, b}], pair},
            pair = Lookup[entryKeys, Key[key], Lookup[exitKeys, Key[key], Missing[]]];
            If[MissingQ[pair], Missing[], DirectedEdge @@ pair]
        ];
        Map[
            With[{forced = canonical[#[[1]], #[[2]]]},
                If[MissingQ[forced], #, forced]
            ] &,
            edges
        ]
    ];

edgeHamiltonianValue[ham_Association, key_String, edge_List] :=
    Module[{edgeAssoc, globalValue},
        edgeAssoc = Lookup[ham, "Edge" <> key, <||>];
        globalValue = Lookup[ham, key, Missing["MissingHamiltonianKey", key]];
        If[!AssociationQ[edgeAssoc] || MissingQ[globalValue],
            Missing["InvalidHamiltonian"],
            Lookup[edgeAssoc, Key[edge],
                Lookup[edgeAssoc, Key[Reverse[edge]], globalValue]]
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
            Return[Missing["InvalidHamiltonian"], Module]];
        jval = edgeSignedFlowValue[sys, edge, rules];
        If[!NumericQ[jval], Return[Missing["MissingFlow"], Module]];
        If[Chop[N[jval]] == 0, Return[0, Module]];
        If[!NumericQ[alpha] || !NumericQ[v], Return[Missing["InvalidHamiltonian"], Module]];
        gexpr = edgeGExpression[g, m];
        If[MissingQ[gexpr], Return[gexpr, Module]];
        equation = N[jval]^2/(2 m^(2 - N[alpha])) - gexpr == -N[v];
        roots = Quiet @ Check[m /. NSolve[{equation, m > 0}, m, Reals], $Failed];
        If[roots === $Failed || roots === {} || !ListQ[roots],
            Missing["DensityNotSolved"],
            FirstCase[N @ roots, value_?NumericQ /; value > 0, Missing["DensityNotSolved"]]
        ]
    ];

(* ===========================================================
   Augmented graph construction (public)
   =========================================================== *)

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

(* ===========================================================
   Vertex styling helpers (uniform coloring rules)
   =========================================================== *)

(* Real network vertices: always gray, regardless of entry/exit role. *)
realVertexGrayMap[vertices_List] :=
    AssociationThread[vertices, GrayLevel[0.72]];

(* Auxiliary vertices: shared entry/exit category colors. *)
auxiliaryVertexColorMap[auxV_List, sys_?mfgSystemQ] :=
    Module[{auxEntryV, auxExitV},
        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];
        Join[
            AssociationThread[Intersection[auxV, auxEntryV], RGBColor[0.38, 0.74, 0.9]],
            AssociationThread[Intersection[auxV, auxExitV],  RGBColor[0.95, 0.7, 0.4]]
        ]
    ];

(* Augmented state-pair colors: only nodes where p[[2]] (label position) is aux entry/exit. *)
stateNodeColorMap[vertices_List, sys_?mfgSystemQ] :=
    Module[{base, auxColors, colored},
        base = AssociationThread[vertices, GrayLevel[0.72]];
        auxColors = auxiliaryVertexColorMap[vertices[[All, 2]], sys];
        colored = Association @ Cases[
            Map[
                With[{vtx = #, color = Lookup[auxColors, Key[#[[2]]], Missing[]]},
                    If[MissingQ[color], Nothing, vtx -> color]
                ] &,
                vertices
            ],
            _Rule
        ];
        Join[base, colored]
    ];

(* Gradient style from u-values, with a fallback association for vertices without u. *)
gradientVertexStyle[vertices_List, uValues_Association, fallback_Association,
                     colorFn_, uMin_, uMax_] :=
    Association @ Map[
        With[{vtx = #, uval = uValues[#]},
            vtx -> If[NumericQ[uval],
                colorFn[If[uMin == uMax, 0.5, Rescale[uval, {uMin, uMax}]]],
                Lookup[fallback, Key[vtx], GrayLevel[0.72]]
            ]
        ] &,
        vertices
    ];

rawVertexCoordinatesFromStateLayout[sys_?mfgSystemQ, vertices_List, layout_] :=
    Module[{augmented, stateVertices, stateEdges, stateGraph, embedding,
            coordinateMap, coordinatesFor, coordinateRules},
        augmented = augmentAuxiliaryGraph[sys];
        stateVertices = augmented["Vertices"];
        stateEdges = Join[augmented["FlowEdges"], augmented["TransitionEdges"]];
        stateGraph = Graph[stateVertices, stateEdges, GraphLayout -> layout];
        embedding = Quiet @ Check[GraphEmbedding[stateGraph], $Failed];
        If[embedding === $Failed || Length[embedding] != Length[VertexList[stateGraph]],
            Return[Automatic, Module]
        ];

        coordinateMap = AssociationThread[VertexList[stateGraph], embedding];
        coordinatesFor[v_] := DeleteMissing @ Map[
            If[MatchQ[#, {_, _}] && #[[2]] === v,
                Lookup[coordinateMap, Key[#], Missing[]],
                Nothing
            ] &,
            stateVertices
        ];

        coordinateRules = Normal @ Association @ Cases[
            Map[
                With[{v = #, coords = coordinatesFor[#]},
                    If[coords === {}, Nothing, v -> Mean[coords]]
                ] &,
                vertices
            ],
            _Rule
        ];
        If[coordinateRules === {}, Automatic, coordinateRules]
    ];

barLegendForU[uMin_, uMax_, colorFn_] :=
    BarLegend[{colorFn[Rescale[#, {uMin, uMax}]] &, {uMin, uMax}},
        LegendLabel -> Placed["u", Right],
        LegendMarkerSize -> 200];

(* ===========================================================
   Edge-label builders (raw plot)
   =========================================================== *)

buildFlowLabels[edges_List, rules_List] :=
    Association @ Map[
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
        edges
    ];

buildValueLabels[edges_List, rules_List] :=
    Module[{valueAt},
        valueAt[a_, b_, t_] :=
            Module[{ua, ub},
                ua = u[a, b] /. rules;
                ub = u[b, a] /. rules;
                If[NumericQ[ua] && NumericQ[ub],
                    (1 - t) N[ua] + t N[ub],
                    Missing["MissingValue"]
                ]
            ];
        Association @ Map[
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
        ]
    ];

buildDensityLabels[s_?scenarioQ, sys_?mfgSystemQ, edges_List, rules_List] :=
    Association @ Map[
        With[{edge = List @@ #, density = edgeDensityValue[s, sys, List @@ #, rules]},
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

(* Combine labels from multiple sources (each per-edge). *)
mergeEdgeLabels[parts___] :=
    Module[{nonempty},
        nonempty = Cases[{parts}, _Association?(Length[#] > 0 &)];
        If[nonempty === {}, <||>,
            Association @ Map[
                Function[edge,
                    edge -> Placed[
                        Column[
                            Cases[
                                Lookup[#, edge, Placed["", Center]] & /@ nonempty,
                                Placed[content_, _] :> content
                            ],
                            Center, Spacings -> 0.2
                        ],
                        Center
                    ]
                ],
                Union @@ (Keys /@ nonempty)
            ]
        ]
    ];

(* ===========================================================
   rawNetworkPlot
   =========================================================== *)

rawNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    rawNetworkPlot[s, sys, <||>, opts];

rawNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, title_String, opts : OptionsPattern[]] :=
    rawNetworkPlot[s, sys, <||>, PlotLabel -> title, opts];

rawNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title_String, opts : OptionsPattern[]] :=
    rawNetworkPlot[s, sys, sol, PlotLabel -> title, opts];

rawNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{rules, model, realV, auxEntryV, auxExitV, auxV, showAux, showBoundary,
            showBoundaryValues, showFlow, showValues, showDensity, showLegend, colorFn,
            edges, vertices, edgeLabels, dispEdges, vertexCoords,
            entryDataAssoc, inAuxEntryPairs, exitCosts, entryEdgeMap,
            uValues, numericU, uMin, uMax,
            baseStyleMap, gradientStyle, vertexStyle, vertexLabels,
            graph, legend, flowLabels, valueLabels, densityLabels, boundaryLabels,
            edgeStyle},

        rules = extractSolutionRules[sol];
        model = scenarioData[s, "Model"];
        realV = model["Vertices"];

        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];

        showBoundary = TrueQ[OptionValue[ShowBoundaryData]];
        showBoundaryValues = TrueQ[OptionValue[ShowBoundaryValues]];
        showAux = TrueQ[OptionValue[ShowAuxiliaryVertices]] || showBoundary;
        showFlow = resolveAutomaticOption[OptionValue[ShowFlowLabels], rules =!= {}];
        showValues = resolveAutomaticOption[OptionValue[ShowValueLabels], rules =!= {}];
        showDensity = TrueQ[OptionValue[ShowDensityLabels]];
        showLegend = TrueQ[OptionValue[ShowLegend]];
        colorFn = Replace[OptionValue[ColorFunction], Automatic -> (Blend[{Blue, Red}, #] &)];

        (* Choose vertex / edge sets based on whether aux is shown. *)
        If[showAux,
            edges = systemData[sys, "AuxEdges"];
            auxV = Complement[systemData[sys, "AuxVertices"], realV];
            vertices = Join[realV, auxV];
            ,
            edges = systemData[sys, "Edges"];
            auxV = {};
            vertices = realV
        ];

        (* Project raw vertex coordinates from the same state-label positions
           used by richNetworkPlot, keeping the raw and rich views aligned. *)
        vertexCoords = rawVertexCoordinatesFromStateLayout[
            sys, vertices, OptionValue[GraphLayout]];

        (* Direct edges by net flow when flow info is being shown; else keep undirected. *)
        dispEdges = If[TrueQ[showFlow], directedDisplayEdges[edges, rules], edges];

        (* Aux entry/exit edges always render with their canonical orientation
           (auxEntry -> N, N -> auxExit), matching richNetworkPlot's augmented
           graph. Independent of flow sign and applies even without a solution,
           as long as aux vertices are shown. *)
        If[showAux,
            dispEdges = forceAuxEdgeDirections[
                dispEdges,
                systemData[sys, "InAuxEntryPairs"],
                systemData[sys, "OutAuxExitPairs"]
            ]
        ];

        (* Vertex base styles: real vertices gray; aux vertices get category color. *)
        baseStyleMap = Join[
            realVertexGrayMap[realV],
            auxiliaryVertexColorMap[auxV, sys]
        ];

        (* For aux vertices only, compute a u-value if available. Real vertices have no
           single canonical u (u is per-edge orientation), so they stay gray. *)
        uValues = Association @ Map[
            Function[v,
                v -> Module[{candidates},
                    candidates = Cases[
                        Join[(u[v, #] /. rules) & /@ realV,
                             (u[#, v] /. rules) & /@ realV],
                        _?NumericQ
                    ];
                    If[candidates === {}, Missing[], First[candidates]]
                ]
            ],
            auxV
        ];
        numericU = Select[uValues, NumericQ];
        If[Length[numericU] > 0,
            uMin = Min[Values[numericU]];
            uMax = Max[Values[numericU]];
            gradientStyle = gradientVertexStyle[
                auxV, uValues, baseStyleMap, colorFn, uMin, uMax];
            (* Real vertices stay gray; gradient overlays aux vertices only. *)
            vertexStyle = Join[
                realVertexGrayMap[realV],
                gradientStyle
            ];
            legend = barLegendForU[uMin, uMax, colorFn];
            ,
            vertexStyle = baseStyleMap;
            legend = None
        ];

        exitCosts = systemData[sys, "ExitCosts"];
        vertexLabels = Association @ Map[
            With[{v = #},
                v -> If[
                    showBoundary && showBoundaryValues && MemberQ[auxExitV, v],
                    formatRawBoundaryExitVertexLabel[
                        v,
                        With[{uv = Lookup[uValues, Key[v], Missing[]]},
                            If[NumericQ[uv], uv, Lookup[exitCosts, Key[v], Missing[]]]
                        ]
                    ],
                    formatRawVertexLabel[v]
                ]
            ] &,
            vertices
        ];

        (* Edge labels per option. *)
        flowLabels = If[TrueQ[showFlow], buildFlowLabels[dispEdges, rules], <||>];
        valueLabels = If[showValues, buildValueLabels[dispEdges, rules], <||>];
        densityLabels = If[showDensity, buildDensityLabels[s, sys, dispEdges, rules], <||>];

        (* Boundary labels for entry-flow / exit-cost (raw plot only labels real edges; we
           overlay onto entry aux edges when ShowBoundaryData and aux are shown). *)
        boundaryLabels = If[showBoundary && showAux,
            entryDataAssoc = systemData[sys, "EntryDataAssociation"];
            inAuxEntryPairs = systemData[sys, "InAuxEntryPairs"];
            entryEdgeMap = Association @ Map[
                With[{p = #}, DirectedEdge[p[[1]], p[[2]]] -> p] &,
                inAuxEntryPairs
            ];
            Association @ Map[
                Module[{e = #, a, b, pair, supply},
                    {a, b} = List @@ e;
                    pair = Lookup[entryEdgeMap, Key[DirectedEdge[a, b]], Missing[]];
                    If[!MissingQ[pair],
                        supply = Lookup[entryDataAssoc, Key[pair], Missing[]];
                        If[!MissingQ[supply],
                            e -> Placed[Style[Row[{"j=", formatPlotNumber[supply]}],
                                              8, Black, Background -> White], Center],
                            Nothing
                        ],
                        Nothing
                    ]
                ] &,
                dispEdges
            ]
            ,
            <||>
        ];

        edgeLabels = mergeEdgeLabels[flowLabels, valueLabels, densityLabels, boundaryLabels];

        edgeStyle = Directive[GrayLevel[0.45], AbsoluteThickness[2], Opacity[0.9]];

        graph = Graph[vertices, dispEdges,
            VertexLabels -> Normal[vertexLabels],
            VertexStyle  -> Normal[vertexStyle],
            VertexSize   -> Normal @ Join[
                AssociationThread[realV, 0.38],
                AssociationThread[auxV, 0.28]
            ],
            VertexCoordinates -> vertexCoords,
            EdgeStyle    -> edgeStyle,
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel], "Network"],
            ImageSize    -> OptionValue[ImageSize]
        ];

        If[legend === None || !showLegend, graph, Legended[graph, legend]]
    ];

(* ===========================================================
   richNetworkPlot
   =========================================================== *)

entryFlowEdge[pair_List] :=
    DirectedEdge[{pair[[2]], pair[[1]]}, {pair[[1]], pair[[2]]}];

richNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, opts : OptionsPattern[]] :=
    richNetworkPlot[s, sys, <||>, opts];

richNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, title_String, opts : OptionsPattern[]] :=
    richNetworkPlot[s, sys, <||>, PlotLabel -> title, opts];

richNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, title_String, opts : OptionsPattern[]] :=
    richNetworkPlot[s, sys, sol, PlotLabel -> title, opts];

richNetworkPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_, opts : OptionsPattern[]] :=
    Module[{rules, augmented, allVertices, flowEdges, transitionEdges, edgeVariables,
            edgeKinds, showFlowEdges, showBoundary, showBoundaryValues,
            showFlowLabels, showValueLabels, useColorFunction, showLegend,
            colorFn, bendFactor,
            auxEntryV, auxExitV, exitCosts, entryDataAssoc, inAuxEntryPairs,
            entryEdgeSet, edges, vertices, effectiveFlows, numericJ, maxJ,
            getU, getEffectiveU, uValues, numericU, uMin, uMax, baseStyleMap,
            gradientStyle, vertexStyle, edgeLabels, nodeLabels,
            graph, legend, nSeg = 20},

        rules = extractSolutionRules[sol];
        augmented = augmentAuxiliaryGraph[sys];
        allVertices = augmented["Vertices"];
        flowEdges = augmented["FlowEdges"];
        transitionEdges = augmented["TransitionEdges"];
        edgeVariables = augmented["EdgeVariables"];
        edgeKinds = augmented["EdgeKinds"];

        showFlowEdges = TrueQ[OptionValue[ShowFlowEdges]];
        showBoundary = TrueQ[OptionValue[ShowBoundaryData]];
        showBoundaryValues = TrueQ[OptionValue[ShowBoundaryValues]];
        showFlowLabels = resolveAutomaticOption[OptionValue[ShowFlowLabels], rules =!= {}];
        showValueLabels = resolveAutomaticOption[OptionValue[ShowValueLabels], rules =!= {}];
        useColorFunction = TrueQ[OptionValue[UseColorFunction]];
        showLegend = TrueQ[OptionValue[ShowLegend]];
        colorFn = Replace[OptionValue[ColorFunction], Automatic -> (Blend[{Blue, Red}, #] &)];
        bendFactor = OptionValue[BendFactor];

        auxEntryV = systemData[sys, "AuxEntryVertices"];
        auxExitV  = systemData[sys, "AuxExitVertices"];
        exitCosts = systemData[sys, "ExitCosts"];

        (* Choose edge set: with/without flow arcs. Recompute vertex set from edges shown. *)
        edges = If[showFlowEdges, Join[flowEdges, transitionEdges], transitionEdges];
        vertices = DeleteDuplicates @ Flatten[List @@@ edges, 1];

        getU[p_] := (u @@ p) /. rules;
        getEffectiveU[p_] := With[{uv = getU[p]},
            If[NumericQ[uv], uv,
                If[showBoundaryValues,
                    With[{cost = Lookup[exitCosts, p[[1]],
                                        Lookup[exitCosts, p[[2]], Missing[]]]},
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
            With[{var = edgeVariables[#],
                  jv = If[rules =!= {}, edgeVariables[#] /. rules /. _j -> 0,
                                        edgeVariables[#]]},
                # -> jv
            ] &,
            edges
        ];
        numericJ = Cases[Values[effectiveFlows], _?NumericQ, Infinity];
        maxJ = Max[Append[numericJ, 1]];

        (* Vertex coloring: state pairs with label-position aux reuse raw aux
           category colors. ColorFunction is opt-in for u-value gradients. *)
        baseStyleMap = stateNodeColorMap[vertices, sys];
        uValues = AssociationMap[getEffectiveU, vertices];
        numericU = Select[uValues, NumericQ];
        If[useColorFunction && Length[numericU] > 0,
            uMin = Min[Values[numericU]];
            uMax = Max[Values[numericU]];
            gradientStyle = gradientVertexStyle[
                vertices, uValues, baseStyleMap, colorFn, uMin, uMax];
            vertexStyle = gradientStyle;
            legend = barLegendForU[uMin, uMax, colorFn];
            ,
            vertexStyle = baseStyleMap;
            legend = None
        ];

        (* Node labels: u-value and (optionally) boundary cost overlays. *)
        nodeLabels = Association @ Map[
            Module[{vtx = #, head, uval, cost},
                head = First[vtx];
                uval = If[showValueLabels, getEffectiveU[vtx], Missing[]];
                If[showBoundary && MemberQ[auxExitV, head],
                    cost = Lookup[exitCosts, head, Missing[]];
                    vtx -> formatBoundaryExitLabel[vtx, uval, cost],
                    vtx -> formatStateNodeLabel[vtx, uval]
                ]
            ] &,
            vertices
        ];

        (* Edge labels: j-values plus optional entry-flow boundary labels. *)
        entryDataAssoc = If[showBoundary, systemData[sys, "EntryDataAssociation"], <||>];
        inAuxEntryPairs = If[showBoundary, systemData[sys, "InAuxEntryPairs"], {}];
        entryEdgeSet = If[showBoundary,
            Association @ Thread[(entryFlowEdge /@ inAuxEntryPairs) -> inAuxEntryPairs],
            <||>
        ];

        edgeLabels = Association @ Map[
            With[{e = #, jv = effectiveFlows[#]},
                Module[{entryPair, supply, body},
                    entryPair = Lookup[entryEdgeSet, Key[e], Missing[]];
                    body = Which[
                        showBoundary && !MissingQ[entryPair],
                            supply = Lookup[entryDataAssoc, Key[entryPair], Missing[]];
                            If[!MissingQ[supply],
                                Row[{"j=", formatPlotNumber[supply]}],
                                If[showFlowLabels && NumericQ[jv] && jv != 0,
                                    NumberForm[N[jv], {5, 1}], ""]],
                        showFlowLabels && NumericQ[jv] && jv != 0,
                            NumberForm[N[jv], {5, 1}],
                        True, ""
                    ];
                    e -> Placed[
                        Style[body, 8, Black, Background -> White],
                        Center
                    ]
                ]
            ] &,
            edges
        ];

        graph = Graph[vertices, edges,
            VertexLabels -> Normal[nodeLabels],
            VertexSize   -> 0.58,
            VertexStyle  -> Normal[vertexStyle],
            EdgeShapeFunction -> Function[{pts, e},
                With[{
                    uS   = uValues[e[[1]]],
                    uT   = uValues[e[[2]]],
                    effJ = effectiveFlows[e],
                    kind = edgeKinds[e],
                    p1   = N[pts[[1]]],
                    p2   = N[pts[[-1]]]
                },
                    If[NumericQ[effJ] && effJ == 0, {},
                    Module[{op, thick, edgeColorAt, d, perp, ctrl, pAt},
                        op    = 0.9;
                        thick = AbsoluteThickness[
                            If[NumericQ[effJ], Rescale[effJ, {0, maxJ}, {1.5, 7}], 2.0]];
                        d    = p2 - p1;
                        perp = With[{len = Norm[d]},
                            If[len > 0, {-d[[2]], d[[1]]}/len, {0, 1}]];
                        ctrl = 0.5*(p1 + p2) + bendFactor*Norm[d]*perp;
                        pAt  = Function[t, (1-t)^2*p1 + 2*(1-t)*t*ctrl + t^2*p2];
                        edgeColorAt = If[
                            useColorFunction && NumericQ[uS] && NumericQ[uT] &&
                                Length[numericU] > 0,
                            colorFn[Rescale[(1 - #)*uS + #*uT, {uMin, uMax}]] &,
                            (Which[
                                kind === "Flow", RGBColor[0.12, 0.45, 0.78],
                                kind === "Transition", RGBColor[0.86, 0.25, 0.22],
                                True, GrayLevel[0.45]
                            ]) &
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
                    ]]
                ]
            ],
            EdgeLabels   -> Normal[edgeLabels],
            GraphLayout  -> OptionValue[GraphLayout],
            PlotLabel    -> plotLabelValue[OptionValue[PlotLabel],
                                            "Augmented infrastructure graph"],
            ImageSize    -> OptionValue[ImageSize]
        ];

        If[legend === None || !showLegend, graph, Legended[graph, legend]]
    ];

End[];

EndPackage[];
