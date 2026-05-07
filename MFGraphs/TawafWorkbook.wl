(* ::Package:: *)

Quit[]


(* ::Subsection:: *)
(*Initialization*)


(* Notebook-friendly MFGraphs workbook for Tawaf scenarios. *)
(* Evaluate cells one at a time or section by section \[LongDash] do not evaluate the entire file at once. *)

(* This file lives alongside MFGraphs.wl in the same directory. *)
mfgDir = If[$InputFileName === "",
    NotebookDirectory[],
    DirectoryName[$InputFileName]
];

If[!StringQ[mfgDir] || mfgDir === "",
    mfgDir = ExpandFileName["."]
];

mfgParentDir = ParentDirectory[mfgDir];
If[!MemberQ[$Path, mfgParentDir], PrependTo[$Path, mfgParentDir]];
Needs["MFGraphs`"];

ClearAll[DescribeOutput];
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

ClearAll[TawafScenarioSummary, TawafCouplingPreview];

(* Construction parameters and topology counts at a glance. *)
TawafScenarioSummary[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{meta, model, edges, js},
        meta  = scenarioData[s, "Tawaf"];
        model = scenarioData[s, "Model"];
        edges = systemData[sys, "Edges"];
        js    = systemData[sys, "Js"];
        <|
            "Rounds"          -> meta["Rounds"],
            "NodesPerRound"   -> meta["NodesPerRound"],
            "Layers"          -> meta["Layers"],
            "TotalVertices"   -> VertexCount[model["Graph"]],
            "DirectedEdges"   -> EdgeCount[model["Graph"]],
            "Entries"         -> model["Entries"],
            "Exits"           -> model["Exits"],
            "FlowVariables"   -> Length[js],
            "TransitionFlows" -> Length[systemData[sys, "Jts"]],
            "ValueVariables"  -> Length[systemData[sys, "Us"]]
        |>
    ];

(* Show one HJ equation per logical edge so the coupling is easy to spot.
   Equations that share a physical segment will mention the same Plus[...] term. *)
TawafCouplingPreview[sys_?mfgSystemQ, maxEqs_:8] :=
    With[{eqs = systemData[sys, "EqGeneral"]},
        Take[eqs, UpTo[maxEqs]]
    ];


(* ::Subsection:: *)
(*Smallest coupled case: 2 rounds, 3 nodes per round, 1 layer*)


(* ::Text:: *)
(*Six logical nodes: round 1 = {1, 2, 3}, round 2 = {4, 5, 6}.*)
(*Tangential forward edges 1\[Rule]2 and 4\[Rule]5 share the physical segment "position 1\[Rule]2"; coupling rewrites either flow as their sum.*)


tawaf2x3 = makeTawafScenario[2, 3, 1];
tawaf2x3System = makeTawafSystem[tawaf2x3];

DescribeOutput[
    "2\[Times]3\[Times]1 scenario summary",
    "Smallest coupled Tawaf system. Two rounds of 3 positions; 5 directed tangential edges (last position of round 2 is the exit).",
    TawafScenarioSummary[tawaf2x3, tawaf2x3System]
]


DescribeOutput[
    "2\[Times]3\[Times]1 physical network",
    "Six logical nodes arranged as two rounds. Entry at vertex 1, exit at vertex 6.",
    scenarioTopologyPlot[tawaf2x3, tawaf2x3System,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] physical topology",
        ImageSize -> Medium]
]


DescribeOutput[
    "2\[Times]3\[Times]1 augmented infrastructure (structure only)",
    "Augmented road-traffic graph before solving. Anti-parallel arcs separate so both directions are visible.",
    mfgAugmentedPlot[tawaf2x3, tawaf2x3System, <||>,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


DescribeOutput[
    "2\[Times]3\[Times]1 EqGeneral preview (after coupling)",
    "Equations that share the same physical position should reference the coupled flow sum. For 2\[Times]3, j[1,2] and j[4,5] are paired (both are forward 1\[Rule]2 segments).",
    TawafCouplingPreview[tawaf2x3System, 6]
]


AbsoluteTiming[tawaf2x3Sol = solveScenario[tawaf2x3];]


DescribeOutput[
    "2\[Times]3\[Times]1 solution validity",
    "isValidSystemSolution checks the solver output against the structural constraints.",
    isValidSystemSolution[tawaf2x3System, tawaf2x3Sol]
]


DescribeOutput[
    "2\[Times]3\[Times]1 augmented infrastructure with solved flows",
    "Blue arcs are flow variables j[a,b]; red arcs are transition flows j[r,i,w]. Node colors show u-values on a Red\[Rule]Blue gradient.",
    mfgAugmentedPlot[tawaf2x3, tawaf2x3System, tawaf2x3Sol,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] solved",
        ImageSize -> Large]
]


DescribeOutput[
    "2\[Times]3\[Times]1 flow plot",
    "Real network with auxiliary entry/exit; edge labels show solved j-values.",
    mfgFlowPlot[tawaf2x3, tawaf2x3System, tawaf2x3Sol,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] flow values",
        ImageSize -> Large]
]


(* ::Subsection:: *)
(*Mid case: 3 rounds, 4 nodes per round, 1 layer*)


(* ::Text:: *)
(*Twelve logical nodes; coupling groups have up to 3 logical flows per physical segment.*)


tawaf3x4 = makeTawafScenario[3, 4, 1];
tawaf3x4System = makeTawafSystem[tawaf3x4];

DescribeOutput[
    "3\[Times]4\[Times]1 scenario summary",
    "Three rounds of 4 positions. Each forward tangential segment (1\[Rule]2, 2\[Rule]3, 3\[Rule]4, 4\[Rule]1) is shared across all three rounds.",
    TawafScenarioSummary[tawaf3x4, tawaf3x4System]
]


DescribeOutput[
    "3\[Times]4\[Times]1 augmented infrastructure (structure only)",
    "Pre-solve augmented graph.",
    mfgAugmentedPlot[tawaf3x4, tawaf3x4System, <||>,
        PlotLabel -> "Tawaf 3\[Times]4\[Times]1 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


(* This solve may take longer than the 2\[Times]3 case; evaluate explicitly. *)
AbsoluteTiming[tawaf3x4Sol = solveScenario[tawaf3x4];]


DescribeOutput[
    "3\[Times]4\[Times]1 solution validity",
    "Validation against the structural constraints.",
    isValidSystemSolution[tawaf3x4System, tawaf3x4Sol]
]


DescribeOutput[
    "3\[Times]4\[Times]1 augmented infrastructure with solved flows",
    "u-value gradient and flow-magnitude edge thickness reveal where congestion concentrates.",
    mfgAugmentedPlot[tawaf3x4, tawaf3x4System, tawaf3x4Sol,
        PlotLabel -> "Tawaf 3\[Times]4\[Times]1 \[LongDash] solved",
        ImageSize -> Large]
]


(* ::Subsection:: *)
(*Multi-layer case: 2 rounds, 3 nodes per round, 2 layers*)


(* ::Text:: *)
(*Twelve nodes; radial edges connect adjacent layers in both directions, so radial flows have outward/inward grouping in addition to tangential.*)


tawaf2x3x2 = makeTawafScenario[2, 3, 2];
tawaf2x3x2System = makeTawafSystem[tawaf2x3x2];

DescribeOutput[
    "2\[Times]3\[Times]2 scenario summary",
    "Two layers; radial edges connect each (round, position) across layers in both directions.",
    TawafScenarioSummary[tawaf2x3x2, tawaf2x3x2System]
]


DescribeOutput[
    "2\[Times]3\[Times]2 augmented infrastructure (structure only)",
    "Pre-solve augmented graph for the multi-layer case.",
    mfgAugmentedPlot[tawaf2x3x2, tawaf2x3x2System, <||>,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]2 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


AbsoluteTiming[tawaf2x3x2Sol = TimeConstrained[solveScenario[tawaf2x3x2], 120, $TimedOut];]


(* If the solve completes in time, render the solved augmented plot. *)
If[Head[tawaf2x3x2Sol] =!= Symbol || tawaf2x3x2Sol =!= $TimedOut,
    DescribeOutput[
        "2\[Times]3\[Times]2 augmented infrastructure with solved flows",
        "Multi-layer solved system. Radial coupling distinguishes outward vs inward.",
        mfgAugmentedPlot[tawaf2x3x2, tawaf2x3x2System, tawaf2x3x2Sol,
            PlotLabel -> "Tawaf 2\[Times]3\[Times]2 \[LongDash] solved",
            ImageSize -> Large]
    ],
    Print["2\[Times]3\[Times]2 solve timed out; skipping solved plot."]
]


(* ::Subsection:: *)
(*Canonical Tawaf: 7 rounds, 8 nodes per round (expensive)*)


(* ::Text:: *)
(*Mirrors the canonical Tawaf reference (seven counter-clockwise circumambulations of an 8-station circuit).*)
(*Solve is expensive; the cell wraps the call in TimeConstrained so it can be aborted safely.*)


tawaf7x8 = makeTawafScenario[7, 8, 1];
tawaf7x8System = makeTawafSystem[tawaf7x8];

DescribeOutput[
    "7\[Times]8\[Times]1 scenario summary",
    "Canonical Tawaf parameters. 56 logical nodes; each forward segment is shared by 7 rounds.",
    TawafScenarioSummary[tawaf7x8, tawaf7x8System]
]


DescribeOutput[
    "7\[Times]8\[Times]1 augmented infrastructure (structure only)",
    "Pre-solve augmented graph for the canonical case.",
    mfgAugmentedPlot[tawaf7x8, tawaf7x8System, <||>,
        PlotLabel -> "Tawaf 7\[Times]8\[Times]1 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


(* Optional expensive solve. Generous timeout; abort manually if needed. *)
(* AbsoluteTiming[tawaf7x8Sol = TimeConstrained[solveScenario[tawaf7x8], 600, $TimedOut];] *)


(* If you computed tawaf7x8Sol above and it is not $TimedOut, render with: *)
(*
mfgAugmentedPlot[tawaf7x8, tawaf7x8System, tawaf7x8Sol,
    PlotLabel -> "Tawaf 7\[Times]8\[Times]1 \[LongDash] solved",
    ImageSize -> Large]
*)
