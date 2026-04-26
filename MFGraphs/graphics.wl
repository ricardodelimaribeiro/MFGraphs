(* Wolfram Language package *)
(* Graphics: public visualization helpers for MFGraphs systems and solutions. *)

BeginPackage["graphics`", {"primitives`", "scenarioTools`", "systemTools`"}];

scenarioTopologyPlot::usage =
"scenarioTopologyPlot[s, sys] plots the scenario topology using vertex coloring for entry, exit, and internal vertices. An optional third argument sets the title.";

mfgSolutionPlot::usage =
"mfgSolutionPlot[s, sys, sol] plots a combined solution graph with real and auxiliary edges. Edge labels include both j and u values, and real-edge directions follow solved net flow orientation. An optional fourth argument sets the title.";

Begin["`Private`"];

extractRules[sol_List] := sol;
extractRules[sol_Association] := Lookup[sol, "Rules", {}];
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
            True, Missing["NotAvailable"]
        ]
    ];

directedDisplayEdges[edges_List, rules_List] :=
    edges /. {
        UndirectedEdge[a_, b_] :> Module[{f = netEdgeFlow[a, b, rules]},
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
        numericJ  = Cases[edgeJValue @@ (List @@ #) & /@ dispEdges, _?NumericQ, Infinity];
        maxJ      = Max[Append[Abs @ numericJ, 1]];

        edgeStyles = Association @ Map[
            With[{a = #[[1]], b = #[[2]], jv = edgeJValue[#[[1]], #[[2]], rules],
                  auxQ = StringStartsQ[ToString[#[[1]]], "aux"] || StringStartsQ[ToString[#[[2]]], "aux"]},
                DirectedEdge[a, b] -> Directive[
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
            List @@@ dispEdges
        ];

        edgeLabels = Association @ Map[
            With[{a = #[[1]], b = #[[2]], jv = edgeJValue[#[[1]], #[[2]], rules], uv = edgeUValue[#[[1]], #[[2]], rules]},
                DirectedEdge[a, b] -> Placed[
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
            List @@@ dispEdges
        ];

        plotTitle = Replace[title, Automatic -> "Solution graph: flows and values"];
        Graph[Join[realV, auxV], dispEdges,
            VertexLabels -> Placed["Name", Center],
            VertexStyle  -> Normal @ Join[
                AssociationThread[entryV,    RGBColor[0.22, 0.6, 0.3]],
                AssociationThread[exitV,     RGBColor[0.82, 0.27, 0.2]],
                AssociationThread[internalV, GrayLevel[0.72]],
                AssociationThread[Select[auxV, StringStartsQ[ToString[#], "auxEntry"] &], RGBColor[0.38, 0.74, 0.9]],
                AssociationThread[Select[auxV, StringStartsQ[ToString[#], "auxExit"] &],  RGBColor[0.95, 0.7, 0.4]]
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

End[];

EndPackage[];
