(*
   Graphics: public visualization helpers for MFGraphs.

   Extracted from Getting started.wl as part of Issue #25.
   This first wave exposes the core network and flow visualizations while
   keeping workbook-specific scenario helpers and density/value plots in the
   workbook for later extraction.

   Load order dependency: plotting helpers expect a compiled equation/system-style
   association with keys such as edgeList, SignedFlows, and boundary-flow rules.
*)

BeginPackage["MFGraphs`"];

(* --- Public API declarations --- *)
NetworkGraphPlot::usage =
"NetworkGraphPlot[d2e] plots the network structure contained in a compiled equation/system-style association d2e. An optional second argument sets the plot title. The function is intended for workbook and package-level visualization of the directed network topology.";

SolutionFlowPlot::usage =
"SolutionFlowPlot[d2e, solution] plots the network with edge labels and edge styling determined by the net edge flows implied by solution. An optional third argument sets the plot title. The input d2e should be a compiled equation/system-style association and solution should be an association of replacement rules or solved values compatible with its symbolic flow expressions.";

ExitFlowPlot::usage =
"ExitFlowPlot[exitFlows] produces a bar chart of exit-flow totals from an association mapping exit vertices to numeric flow values. An optional second argument sets the plot title.";

AssociationValue::usage =
"AssociationValue[assoc, key, default] returns assoc[key] when key exists, otherwise default (Missing[\"NotAvailable\"] by default).";

NetEdgeFlows::usage =
"NetEdgeFlows[d2e, solution, pairs] returns the net signed flow on each requested edge pair after applying balance and boundary flow rules. pairs defaults to all model edges.";

NetworkVisualData::usage =
"NetworkVisualData[d2e] builds reusable graph layout and vertex styling metadata used by plotting helpers.";

FlowStyleDirective::usage =
"FlowStyleDirective[flow, maxFlow] returns an edge style directive (color/thickness/opacity) scaled by flow magnitude and sign.";

Begin["`Private`"];

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

Options[NetworkGraphPlot] = {};

NetworkGraphPlot[d2e_Association, title_: Automatic, OptionsPattern[]] :=
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

Options[SolutionFlowPlot] = {};

SolutionFlowPlot[d2e_Association, solution_Association, title_: Automatic, OptionsPattern[]] :=
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

Options[ExitFlowPlot] = {};

ExitFlowPlot[exitFlows_Association, title_: Automatic, OptionsPattern[]] :=
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

End[];

EndPackage[];
