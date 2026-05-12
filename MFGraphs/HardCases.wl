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
(*	3. Grid0505 \[LongDash] 5*5 grid; clean 10-min TimeConstrained*)
(*	4. Grid0707 \[LongDash] 7*7 grid; clean 10-min TimeConstrained*)
(*	5. Grid0710 \[LongDash] 7*10 grid; clean 10-min TimeConstrained*)
(*	6. Grid1010 \[LongDash] 10*10 grid; clean 10-min TimeConstrained*)
(*	7. Grid1020 \[LongDash] 10*20 grid; clean 10-min TimeConstrained (also tripped a license error in one run)*)


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


ClearAll[DescribeOutput, HardCaseSummary, HardCaseDNFDiagnostic,
    hardCaseVertexCount, hardCaseEdgeCount];

(* Per-case DNF reducer timeout in seconds. dnfReduce is package-level recursive
   code so TimeConstrained inside dnfReduceDiagnosticReport actually interrupts
   it (unlike the C-level Reduce inside reduceSystem). Bump if you want a deeper
   trace; budget ~ $hardCaseDNFTimeout * 7 per evaluate-all run. *)
$hardCaseDNFTimeout = 30;

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

(* model["Graph"] is only present for scenarios built from a Graph object
   (graphScenario, gridScenario, cycleScenario). amScenario-based scenarios
   store Vertices + Adjacency instead. Read whichever representation exists. *)
hardCaseVertexCount[model_Association] :=
    Which[
        KeyExistsQ[model, "Vertices"], Length[model["Vertices"]],
        KeyExistsQ[model, "Graph"],    VertexCount[model["Graph"]],
        True,                          Missing["Unknown"]
    ];

hardCaseEdgeCount[model_Association] :=
    Which[
        KeyExistsQ[model, "Adjacency"], Total[Flatten[model["Adjacency"]]],
        KeyExistsQ[model, "Graph"],     EdgeCount[model["Graph"]],
        True,                           Missing["Unknown"]
    ];

(* Run the recursive DNF reducer with an internal timeout and return the
   diagnostic association: Status, Summary (branch/conjunct counts + timings),
   OrderingSummary, LastConjunct, LastEvents (recent trace), TimeoutLocation. *)
HardCaseDNFDiagnostic[sys_?mfgSystemQ, timeout_:Automatic] :=
    dnfReduceDiagnosticReport[sys,
        "Timeout" -> If[timeout === Automatic, $hardCaseDNFTimeout, timeout]];

HardCaseSummary[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{model, augmented},
        model = scenarioData[s, "Model"];
        augmented = augmentAuxiliaryGraph[sys];
        <|
            "Vertices" -> hardCaseVertexCount[model],
            "DirectedEdges" -> hardCaseEdgeCount[model],
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

$MFGraphsVerbose = True;


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


DescribeOutput[
    "case 21 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout. Status=OK means the DNF phase finished; Status=Timeout still reports how many branches/conjuncts were processed before the budget ran out.",
    HardCaseDNFDiagnostic[case21System]
]


(* ::Text:: *)
(*Observed zero-flow edges (by inspection of the augmented topology and boundary data):*)
(*	\[Bullet] j[3, 1] == 0  \[LongDash] no return flow from vertex 3 to entry 1*)
(*	\[Bullet] j[4, 2] == 0  \[LongDash] no return flow from vertex 4 to entry 2*)
(*	\[Bullet] j[4, 3, 1] == 0  \[LongDash] transition flow 4\[Rule]3\[Rule]1 forbidden (U-turn back toward entry 1)*)
(*	\[Bullet] j[3, 4, 2] == 0  \[LongDash] transition flow 3\[Rule]4\[Rule]2 forbidden (U-turn back toward entry 2)*)
(*	\[Bullet] j[10, 7] == 0  \[LongDash] no flow from exit 10 back upstream*)
(*	\[Bullet] j[11, 8] == 0  \[LongDash] no flow from exit 11 back upstream*)
(*	\[Bullet] j[12, 9] == 0  \[LongDash] no flow from exit 12 back upstream*)


(* ::Text:: *)
(*Sanity-check these symbolically with: Select[systemData[case21System, "Js"], MemberQ[{j[3,1], j[4,2], j[10,7], j[11,8], j[12,9]}, #] &] (likewise for Jts).*)


 (*solveScenario[case21]*)


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


DescribeOutput[
    "HRF Scenario 1 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout.",
    HardCaseDNFDiagnostic[hrfSystem]
]


(* solveScenario[hrfScenario] *)


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


DescribeOutput[
    "Grid0505 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout.",
    HardCaseDNFDiagnostic[grid0505System]
]


(* solveScenario[grid0505] *)


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


DescribeOutput[
    "Grid0707 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout.",
    HardCaseDNFDiagnostic[grid0707System]
]


(* solveScenario[grid0707] *)


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


DescribeOutput[
    "Grid0710 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout.",
    HardCaseDNFDiagnostic[grid0710System]
]


(* solveScenario[grid0710] *)


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


DescribeOutput[
    "Grid1010 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout.",
    HardCaseDNFDiagnostic[grid1010System]
]


(* solveScenario[grid1010] *)


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


DescribeOutput[
    "Grid1020 dnfReduce diagnostic",
    "Recursive DNF reducer with " <> ToString[$hardCaseDNFTimeout] <> "s timeout. Largest grid \[LongDash] often times out before BranchesStarted moves off zero.",
    HardCaseDNFDiagnostic[grid1020System]
]


(* solveScenario[grid1020] *)


(* ::Subsection:: *)
(*Summary*)


(* ::Text:: *)
(*All seven scenarios construct (scenarioQ + mfgSystemQ both pass) and have a fully-populated boundary data layer; the failure mode is purely solver-side. The augmented graphs render quickly, so even the 10*20 case is interactively explorable in the front end.*)


(* ::Text:: *)
(*Next steps if you want to solve these: try the active-set reducer (activeSetReduceSystem) instead of the DNF-first default; consider a numeric backend; or split the LCP by route to keep per-branch Reduce calls tractable.*)
