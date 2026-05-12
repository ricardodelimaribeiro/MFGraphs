(* ::Package:: *)

(*Quit[]*)


(* ::Title:: *)
(*MFGraphs Jamarat workbook*)


(* ::Subsection:: *)
(*Overview*)


(* ::Text:: *)
(*A guided tour of MFGraphs on the Jamarat scenario family \[LongDash] multi-entrance / multi-exit infrastructure problems. Each section names the package capability it demonstrates. Evaluate cells one at a time or section by section; the Quit[] at the top is a safety device.*)


(* ::Text:: *)
(*Capabilities demonstrated, in order:*)
(*	1. Direct cycle scenario construction (cycleScenario) and a simple multi-exit case*)
(*	2. Named registry retrieval (getExampleScenario["Jamaratv9", ...])*)
(*	3. Augmented state-space graph (richNetworkPlot, augmentAuxiliaryGraph)*)
(*	4. Symbolic solving (solveScenario, dnf-first default)*)
(*	5. Two-mode solution visualisation \[LongDash] flow-coloured edges and value-function-coloured nodes (UseColorFunction -> True)*)
(*	6. Congestion behaviour when entry flow exceeds exit budget (Section 3, optional, expensive)*)


(* ::Subsection:: *)
(*Initialization*)


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


(* ::Subsection:: *)
(*Presentation helpers*)


ClearAll[DescribeOutput, JamaratScenarioSummary];

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

JamaratScenarioSummary[s_?scenarioQ, sys_?mfgSystemQ] :=
    Module[{model, entries, exits, augmented},
        model = scenarioData[s, "Model"];
        entries = model["Entries"];
        exits = model["Exits"];
        augmented = augmentAuxiliaryGraph[sys];
        <|
            "Entries" -> entries,
            "EntryFlowTotal" -> Total[Last /@ entries],
            "Exits" -> exits,
            "ExitCostTotal" -> Total[Last /@ exits],
            "EntryFlowExceedsExitBudget" -> (Total[Last /@ entries] > Total[Last /@ exits]),
            "NetworkEdges" -> Length[systemData[sys, "Edges"]],
            "FlowVariables" -> Length[systemData[sys, "Js"]],
            "TransitionFlowVariables" -> Length[systemData[sys, "Jts"]],
            "ValueVariables" -> Length[systemData[sys, "Us"]],
            "AugmentedFlowEdges" -> Length[augmented["FlowEdges"]],
            "AugmentedTransitionEdges" -> Length[augmented["TransitionEdges"]]
        |>
    ];

$MFGraphsVerbose = False;


(* ::Section:: *)
(*Section 1 \[LongDash] Simplified Jamarat (5-cycle)*)


(* ::Text:: *)
(*A 5-vertex cycle with two entrances (vertices 1 and 5) and three exits (vertex 3 cost 0, vertex 4 cost 40, vertex 2 cost 30). The smallest case where multiple entries compete for shared downstream capacity. Capability: cycleScenario direct constructor + the standard makeSystem / solveScenario pipeline.*)


jamaratEnd = cycleScenario[5, {{1, 100}, {5, 100}}, {{3, 0}, {4, 40}, {2, 30}}];
jamaratEndSystem = makeSystem[jamaratEnd];

DescribeOutput[
    "Simplified Jamarat scenario summary",
    "Two entries, three exits on a 5-cycle. Entry-flow total 200; exit-cost total 70.",
    JamaratScenarioSummary[jamaratEnd, jamaratEndSystem]
]


DescribeOutput[
    "Simplified Jamarat physical topology (rawNetworkPlot)",
    "ShowBoundaryData -> True overlays entry-flow and exit-cost annotations on the auxiliary nodes.",
    rawNetworkPlot[jamaratEnd, jamaratEndSystem,
        PlotLabel -> "Simplified Jamarat \[LongDash] physical topology",
        ShowBoundaryData -> True,
        ImageSize -> Large]
]


DescribeOutput[
    "Simplified Jamarat augmented state-space graph (richNetworkPlot, structure only)",
    "Augmented graph before solving. Anti-parallel arcs separate so both directions are visible.",
    richNetworkPlot[jamaratEnd, jamaratEndSystem,
        PlotLabel -> "Simplified Jamarat \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


AbsoluteTiming[jamaratEndSol = solveScenario[jamaratEnd];]


DescribeOutput[
    "Simplified Jamarat augmented solution (default colouring)",
    "Augmented graph with solved flow, transition, and u values; flow-coloured edges.",
    richNetworkPlot[jamaratEnd, jamaratEndSystem, jamaratEndSol,
        PlotLabel -> "Simplified Jamarat \[LongDash] solved (flow colours)",
        ImageSize -> Large]
]


DescribeOutput[
    "Simplified Jamarat augmented solution (value-coloured)",
    "Same solution rendered with UseColorFunction -> True: nodes coloured by their value-function value on a Blue\[Rule]Red gradient.",
    richNetworkPlot[jamaratEnd, jamaratEndSystem, jamaratEndSol,
        PlotLabel -> "Simplified Jamarat \[LongDash] solved (value colours)",
        UseColorFunction -> True,
        ImageSize -> Large]
]


(* ::Section:: *)
(*Section 2 \[LongDash] Jamaratv9 from the named registry*)


(* ::Text:: *)
(*The canonical Jamaratv9 topology: 9 vertices, two upstream entrances (vertices 1 and 2), three downstream exits (vertices 7, 8, 9). Capability: getExampleScenario from the named registry, end-to-end solve, two-mode solution visualisation.*)


jamaratScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 100}, {2, 100}},
    {{7, 0}, {8, 0}, {9, 0}}
];
jamaratSystem = makeSystem[jamaratScenario];

DescribeOutput[
    "Jamaratv9 scenario summary",
    "Equal-cost exits and equal entries: a balanced baseline. The solver should split flows symmetrically.",
    JamaratScenarioSummary[jamaratScenario, jamaratSystem]
]


DescribeOutput[
    "Jamaratv9 physical topology (rawNetworkPlot)",
    "Real network with auxiliary entry/exit nodes shown. Entries (blue) at 1 and 2; exits (orange) at 7, 8, 9.",
    rawNetworkPlot[jamaratScenario, jamaratSystem,
        PlotLabel -> "Jamaratv9 \[LongDash] physical topology",
        ShowAuxiliaryVertices -> True,
        ImageSize -> Large]
]


DescribeOutput[
    "Jamaratv9 augmented state-space graph (structure only)",
    "Augmented graph before solving; node colours indicate boundary status only.",
    richNetworkPlot[jamaratScenario, jamaratSystem,
        PlotLabel -> "Jamaratv9 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


AbsoluteTiming[jamaratSol = solveScenario[jamaratScenario];]
If[Head[jamaratSol] =!= Association,
    jamaratSol = <|"Rules" -> jamaratSol, "Residual" -> True|>];


DescribeOutput[
    "Jamaratv9 solution (default colouring)",
    "Edges show solved j-values; node colours indicate boundary category.",
    richNetworkPlot[jamaratScenario, jamaratSystem, jamaratSol,
        PlotLabel -> "Jamaratv9 \[LongDash] solved (flow colours)",
        ImageSize -> Large]
]


DescribeOutput[
    "Jamaratv9 solution (value-coloured)",
    "UseColorFunction -> True: node colours show the equilibrium value function on a Blue\[Rule]Red gradient. Vertices closer to the cheapest exits are on the cool end.",
    richNetworkPlot[jamaratScenario, jamaratSystem, jamaratSol,
        PlotLabel -> "Jamaratv9 \[LongDash] solved (value colours)",
        UseColorFunction -> True,
        ImageSize -> Large]
]


(* ::Section:: *)
(*Section 3 \[LongDash] High cost exit avoided variant (optional, 36 seconds)*)


(* ::Text:: *)
(*Same Jamaratv9 topology, but exit costs are unbalanced \[LongDash] vertex 9 has cost 55 while 7 and 8 are free. Entries are smaller (20 + 50). The interesting feature is that the cheapest exit is downstream of a longer route, while the costly exit is closer; congestion vs path length forces a non-trivial split. The solve is wrapped in TimeConstrained because the search space grows with the asymmetry.*)


jamaratHighEntryScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 20}, {2, 50}},
    {{7, 0}, {8, 0}, {9, 55}}
];
jamaratHighEntrySystem = makeSystem[jamaratHighEntryScenario];

DescribeOutput[
    "High-entry-flow Jamaratv9 scenario summary",
    "Asymmetric entries (20 + 50) and asymmetric exit costs (free, free, 55).",
    JamaratScenarioSummary[jamaratHighEntryScenario, jamaratHighEntrySystem]
]


(* Generous timeout. Abort manually if needed. *)
AbsoluteTiming[
    jamaratHighEntrySol =
        TimeConstrained[solveScenario[jamaratHighEntryScenario], 600, $TimedOut];
]


If[AssociationQ[jamaratHighEntrySol] || ListQ[jamaratHighEntrySol],
    DescribeOutput[
        "High-entry-flow Jamaratv9 solution (value-coloured)",
        "Solved augmented graph with value-coloured nodes. Compare to Section 2: the value-function gradient now reflects the asymmetric exit costs.",
        richNetworkPlot[jamaratHighEntryScenario, jamaratHighEntrySystem, jamaratHighEntrySol,
            PlotLabel -> "Jamaratv9 high-entry \[LongDash] solved (value colours)",
            UseColorFunction -> True,
            ImageSize -> Large]
    ],
    Print["High-entry Jamaratv9 solve did not return a solution (",
        jamaratHighEntrySol, "); skipping solved plot."]
]


(* ::Subsection:: *)
(*Summary \[LongDash] capabilities demonstrated*)


(* ::Text:: *)
(*Direct constructors. cycleScenario built the simplified 5-cycle case; getExampleScenario pulled the named "Jamaratv9" topology from the registry. Both expose the same downstream pipeline.*)


(* ::Text:: *)
(*Multi-entrance / multi-exit handling. makeSystem builds entry/exit balance equations transparently when entries and exits are lists of pairs; no Jamarat-specific code path.*)


(* ::Text:: *)
(*Two-mode solution visualisation. richNetworkPlot defaults to flow-coloured edges; UseColorFunction -> True switches to value-function-coloured nodes on a Blue->Red gradient. The two views answer different questions about the same solution.*)


(* ::Text:: *)
(*Pointers: package-side API in CLAUDE.md and the auto-generated API_REFERENCE.md. The named registry behind getExampleScenario is documented at MFGraphs/examples.wl; use listExampleScenarios[] to discover all keys.*)
