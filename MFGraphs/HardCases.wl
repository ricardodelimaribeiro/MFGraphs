(* ::Package:: *)

(*Quit[]*)


(* ::Title:: *)
(*MFGraphs hard cases workbook*)


(* ::Subsection:: *)
(*Overview*)


(* ::Text:: *)
(*The seven registered example scenarios that the current symbolic solver does not finish on within a 10-minute budget. Each section shows the typed scenario + the augmented road-traffic graph with boundary data, so the topology is visible without running solveScenario.*)


(* ::Text:: *)
(*Cases included:*)
(*	1. case 21 \[LongDash] 12-vertex multi-entrance / multi-exit; crashes the kernel mid-Reduce*)
(*	2. HRF Scenario 1 \[LongDash] 10-vertex multi-route benchmark; Reduce never returns*)
(*	3. Grid0505 \[LongDash] 5\[Times]5 grid; clean 10-min TimeConstrained*)
(*	4. Grid0707 \[LongDash] 7\[Times]7 grid; clean 10-min TimeConstrained*)
(*	5. Grid0710 \[LongDash] 7\[Times]10 grid; clean 10-min TimeConstrained*)
(*	6. Grid1010 \[LongDash] 10\[Times]10 grid; clean 10-min TimeConstrained*)
(*	7. Grid1020 \[LongDash] 10\[Times]20 grid; clean 10-min TimeConstrained (also tripped a license error in one run)*)


(* ::Text:: *)
(*Companion to solutions/README.md and MFGraphs/Tests/example-coverage.mt: the cache fails to cover exactly these seven scenarios.*)


(* ::Subsection:: *)
(*Initialization*)


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


(* ::Subsection:: *)
(*Presentation helpers*)


ClearAll[DescribeOutput, HardCaseSummary];

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

HardCaseSummary[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{model, augmented},
        model = scenarioData[s, "Model"];
        augmented = augmentAuxiliaryGraph[sys];
        <|
            "Vertices" -> VertexCount[model["Graph"]],
            "DirectedEdges" -> EdgeCount[model["Graph"]],
            "Entries" -> model["Entries"],
            "Exits" -> model["Exits"],
            "EntryFlowTotal" -> Total[Last /@ model["Entries"]],
            "FlowVariables" -> Length[systemData[sys, "Js"]],
            "TransitionFlowVariables" -> Length[systemData[sys, "Jts"]],
            "ValueVariables" -> Length[systemData[sys, "Us"]],
            "AugmentedFlowEdges" -> Length[augmented["FlowEdges"]],
            "AugmentedTransitionEdges" -> Length[augmented["TransitionEdges"]]
        |>
    ];

$MFGraphsVerbose = False;


(* ::Section:: *)
(*Section 1 \[LongDash] case 21 (12-vertex multi-entrance / multi-exit)*)


(* ::Text:: *)
(*Reason for inclusion: solveScenario crashes the Wolfram kernel mid-Reduce on this scenario. Symbolic LCP reduction explodes on a topology that has multiple sources, multiple sinks, and intermediate routing choices.*)


case21 = getExampleScenario[21, {{1, 50}, {2, 50}}, {{10, 0}, {11, 0}, {12, 0}}];
case21System = makeSystem[case21];

DescribeOutput[
    "case 21 summary",
    "12 vertices, two entries, three exits.",
    HardCaseSummary[case21, case21System]
]


DescribeOutput[
    "case 21 augmented road-traffic graph (with boundary data)",
    "ShowBoundaryData -> True overlays entry-flow and exit-cost labels on the auxiliary boundary nodes.",
    richNetworkPlot[case21, case21System,
        PlotLabel -> "case 21 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Section:: *)
(*Section 2 \[LongDash] HRF Scenario 1*)


(* ::Text:: *)
(*Reason for inclusion: Reduce never returns. TimeConstrained does not interrupt because the bottleneck is a single C-level Reduce call.*)


hrfScenario = getExampleScenario["HRF Scenario 1", {{1, 100}}, {{8, 0}, {10, 0}}];
hrfSystem = makeSystem[hrfScenario];

DescribeOutput[
    "HRF Scenario 1 summary",
    "10 vertices, one entry, two exits.",
    HardCaseSummary[hrfScenario, hrfSystem]
]


DescribeOutput[
    "HRF Scenario 1 augmented road-traffic graph (with boundary data)",
    "Multi-route benchmark with several parallel paths between the entry and the two exits.",
    richNetworkPlot[hrfScenario, hrfSystem,
        PlotLabel -> "HRF Scenario 1 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Section:: *)
(*Section 3 \[LongDash] Grid scenarios*)


(* ::Text:: *)
(*All five grids share the same boundary spec (single entry at vertex 1, single exit at the last vertex, entry flow 100) so the timeout pattern can be read off vertex count alone.*)


(* ::Subsubsection:: *)
(*3.1 Grid0505 (25 vertices)*)


grid0505 = getExampleScenario["Grid0505", {{1, 100}}, {{25, 0}}];
grid0505System = makeSystem[grid0505];

DescribeOutput[
    "Grid0505 summary",
    "5\[Times]5 grid, vertices 1..25 in row-major order.",
    HardCaseSummary[grid0505, grid0505System]
]


DescribeOutput[
    "Grid0505 augmented road-traffic graph (with boundary data)",
    "Single entry at vertex 1 (blue), single exit at vertex 25 (orange).",
    richNetworkPlot[grid0505, grid0505System,
        PlotLabel -> "Grid0505 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Subsubsection:: *)
(*3.2 Grid0707 (49 vertices)*)


grid0707 = getExampleScenario["Grid0707", {{1, 100}}, {{49, 0}}];
grid0707System = makeSystem[grid0707];

DescribeOutput[
    "Grid0707 summary",
    "7\[Times]7 grid, vertices 1..49 in row-major order.",
    HardCaseSummary[grid0707, grid0707System]
]


DescribeOutput[
    "Grid0707 augmented road-traffic graph (with boundary data)",
    "Single entry at vertex 1, single exit at vertex 49.",
    richNetworkPlot[grid0707, grid0707System,
        PlotLabel -> "Grid0707 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Subsubsection:: *)
(*3.3 Grid0710 (70 vertices)*)


grid0710 = getExampleScenario["Grid0710", {{1, 100}}, {{70, 0}}];
grid0710System = makeSystem[grid0710];

DescribeOutput[
    "Grid0710 summary",
    "7\[Times]10 grid, vertices 1..70 in row-major order. First non-square hard case.",
    HardCaseSummary[grid0710, grid0710System]
]


DescribeOutput[
    "Grid0710 augmented road-traffic graph (with boundary data)",
    "Single entry at vertex 1, single exit at vertex 70.",
    richNetworkPlot[grid0710, grid0710System,
        PlotLabel -> "Grid0710 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Subsubsection:: *)
(*3.4 Grid1010 (100 vertices)*)


grid1010 = getExampleScenario["Grid1010", {{1, 100}}, {{100, 0}}];
grid1010System = makeSystem[grid1010];

DescribeOutput[
    "Grid1010 summary",
    "10\[Times]10 grid, vertices 1..100 in row-major order.",
    HardCaseSummary[grid1010, grid1010System]
]


DescribeOutput[
    "Grid1010 augmented road-traffic graph (with boundary data)",
    "Single entry at vertex 1, single exit at vertex 100. Plot is dense; consider zooming in interactively.",
    richNetworkPlot[grid1010, grid1010System,
        PlotLabel -> "Grid1010 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Subsubsection:: *)
(*3.5 Grid1020 (200 vertices)*)


grid1020 = getExampleScenario["Grid1020", {{1, 100}}, {{200, 0}}];
grid1020System = makeSystem[grid1020];

DescribeOutput[
    "Grid1020 summary",
    "10\[Times]20 grid, vertices 1..200 in row-major order. Largest registered grid; render is intentionally heavy.",
    HardCaseSummary[grid1020, grid1020System]
]


DescribeOutput[
    "Grid1020 augmented road-traffic graph (with boundary data)",
    "Single entry at vertex 1, single exit at vertex 200. Plot may take several seconds to render in the front end.",
    richNetworkPlot[grid1020, grid1020System,
        PlotLabel -> "Grid1020 \[LongDash] augmented (boundary data)",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


(* ::Subsection:: *)
(*Summary*)


(* ::Text:: *)
(*All seven scenarios construct (scenarioQ + mfgSystemQ both pass) and have a fully-populated boundary data layer; the failure mode is purely solver-side. The augmented graphs render quickly, so even the 10\[Times]20 case is interactively explorable in the front end.*)


(* ::Text:: *)
(*Next steps if you want to solve these: try the active-set reducer (activeSetReduceSystem) instead of the DNF-first default; consider a numeric backend; or split the LCP by route to keep per-branch Reduce calls tractable.*)
