(* ::Package:: *)

(*Quit[]*)


(* ::Title:: *)
(*MFGraphs Capabilities Tour \[LongDash] the Tawaf scenario*)


(* ::Subsection:: *)
(*Overview*)


(* ::Text:: *)
(*This document presents the MFGraphs package by following one scenario family \[LongDash] Tawaf circumambulation \[LongDash] through the full pipeline. Each section names the package capability it demonstrates so the reader can see what MFGraphs offers without leaving the workbook.*)


(* ::Text:: *)
(*Capabilities demonstrated, in order:*)
(*	1. Typed scenario construction with metadata storage (makeTawafScenario, scenarioQ, scenarioData)*)
(*	2. Symbolic unknown bundle and structural system (makeSystem behind makeTawafSystem; symbolicUnknowns; mfgSystem)*)
(*	3. Tawaf-specific physical-edge coupling \[LongDash] the package's specialised builder pattern (makeTawafSystem rewrites EqGeneral and AltOptCond on top of an ordinary mfgSystem)*)
(*	4. Symbolic solving via the DNF-first default (dnfReduceSystem applied to the Tawaf-coupled system)*)
(*	5. Solution validation (isValidSystemSolution)*)
(*	6. Visualisation surface (rawNetworkPlot, richNetworkPlot, plus a custom 3D helix view defined locally)*)
(*	7. Scaling from a 6-node smoke test through 12, 12, and 56 nodes*)
(*	8. Density-dependent cost extension (\"Density\" -> True, m family, EqDensityFlow, tawafDensities)*)


(* ::Text:: *)
(*Ritual context, the "unroll then couple" modelling idea, the scope of the coupling rewrite, and the parameter provenance (Black Stone potential, layer interpretation, entry flow) live in docs/research/notes/tawaf-model.md. This document focuses on the package-side surface; that note focuses on what the model actually represents.*)


(* ::Subsubsection:: *)
(*How to read this document*)


(* ::Text:: *)
(*Evaluate cells one at a time or section by section \[LongDash] do not evaluate the entire file at once. The Quit[] at the top is a safety device; if it fires, evaluation stops cleanly. Larger scenarios (3*4 and 2*3*2) are wrapped or noted as expensive; the canonical 7*8 case ships with its solve commented out.*)


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


(* ::Text:: *)
(*DescribeOutput wraps any expression with a title and a one-line description so each cell reads as a slide. TawafScenarioSummary, TawafCouplingPreview, and TawafCouplingComparison surface the package's structural data; tawafHelixPlot is a custom 3D layout that takes advantage of the fact that the unrolled rounds stack vertically into a helix. None of these are part of the package surface \[LongDash] they are presentation helpers defined locally.*)


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

ClearAll[TawafScenarioSummary, TawafCouplingPreview, TawafCouplingComparison];

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

(* Side-by-side comparison: the same EqGeneral equations before and after
   the Tawaf coupling rewrite. The "before" column comes from the bare
   makeSystem; the "after" column from makeTawafSystem. The contrast is
   the central package novelty for this scenario family. *)
TawafCouplingComparison[s_?scenarioQ, maxEqs_:4] :=
    Module[{rawSys, coupledSys, rawEqs, coupledEqs, n, rows},
        rawSys     = makeSystem[s];
        coupledSys = makeTawafSystem[s];
        rawEqs     = systemData[rawSys,     "EqGeneral"];
        coupledEqs = systemData[coupledSys, "EqGeneral"];
        n = Min[maxEqs, Length[rawEqs], Length[coupledEqs]];
        rows = Table[
            {i, rawEqs[[i]], coupledEqs[[i]]},
            {i, n}
        ];
        Grid[
            Prepend[rows,
                Style[#, Bold, GrayLevel[0.3]] & /@
                    {"index", "before coupling (makeSystem)", "after coupling (makeTawafSystem)"}],
            Frame    -> All,
            Alignment -> {{Center, Left, Left}, Top},
            Background -> {None, {RGBColor[0.95, 0.95, 0.97], None}},
            Spacings -> {1.2, 0.8}
        ]
    ];


(* ::Subsubsection:: *)
(*Helix plot helper*)


ClearAll[tawafHelixPlot];

Options[tawafHelixPlot] = {
    "ShowEquivalenceLinks" -> False,
    "ShowVertexLabels"     -> Automatic,
    "Pitch"                -> 0.6,
    "Radius"               -> 3.0,
    "LayerGap"             -> 1.2,
    "ColorByRound"         -> True,
    "EdgeThicknessByFlow"  -> True,
    PlotLabel              -> Automatic,
    ImageSize              -> Large
};

(* 3D helix visualization for a Tawaf scenario.
   Rounds stack vertically (one full turn of angular advance per round, height
   advancing one tick per position-step), so equivalent (round-shifted) nodes
   sit directly above each other. Layers become concentric helices.

   Defined locally in this presentation workbook; could be promoted into
   graphicsTools` if the helix view becomes a recurring need. *)
tawafHelixPlot[s_?scenarioQ, sys_?mfgSystemQ, sol_:<||>, opts:OptionsPattern[]] :=
    Module[{meta, rounds, npr, layers, dz, r0, dr, showLinks, showLabels,
            colorByR, thickByJ, encode, coord, modelEdges, integerEdges,
            equivEdges, vertexCoords, ruleList, flowVals, jMax, thicknessFor,
            edgeStyleRules, vertexStyles, vertices, plotLabel, imgSize,
            roundColor, totalNodes},
        meta   = scenarioData[s, "Tawaf"];
        If[!AssociationQ[meta],
            Return[Failure["tawafHelixPlot",
                <|"MessageTemplate" -> "Scenario lacks Tawaf metadata."|>]]];
        rounds = meta["Rounds"]; npr = meta["NodesPerRound"]; layers = meta["Layers"];
        totalNodes = rounds * npr * layers;

        dz        = OptionValue["Pitch"];
        r0        = OptionValue["Radius"];
        dr        = OptionValue["LayerGap"];
        showLinks = TrueQ[OptionValue["ShowEquivalenceLinks"]];
        showLabels = Replace[OptionValue["ShowVertexLabels"],
            Automatic :> (totalNodes <= 24)];
        colorByR  = TrueQ[OptionValue["ColorByRound"]];
        thickByJ  = TrueQ[OptionValue["EdgeThicknessByFlow"]];
        plotLabel = Replace[OptionValue[PlotLabel], Automatic ->
            StringTemplate["Tawaf helix ``\[Times]``\[Times]``"][rounds, npr, layers]];
        imgSize   = OptionValue[ImageSize];

        encode[r_, p_, l_] := tawafEncode[r, p, l, rounds, npr];
        coord[r_, p_, l_]  := {
            (r0 + (l - 1) dr) Cos[2 Pi (p - 1)/npr],
            (r0 + (l - 1) dr) Sin[2 Pi (p - 1)/npr],
            ((r - 1) npr + (p - 1)) dz
        };

        vertices = Flatten @ Table[encode[r, p, l],
            {l, layers}, {r, rounds}, {p, npr}];

        vertexCoords = Flatten[
            Table[encode[r, p, l] -> coord[r, p, l],
                {l, layers}, {r, rounds}, {p, npr}],
            2];

        modelEdges = EdgeList[scenarioData[s, "Model"]["Graph"]];
        integerEdges = Select[modelEdges,
            IntegerQ[#[[1]]] && IntegerQ[#[[2]]] &];

        equivEdges = If[showLinks && rounds >= 2,
            Flatten @ Table[
                UndirectedEdge[encode[r, p, l], encode[r + 1, p, l]],
                {l, layers}, {p, npr}, {r, rounds - 1}],
            {}];

        ruleList = Which[
            sol === <||> || sol === {}, {},
            AssociationQ[sol], Lookup[sol, "Rules", {}],
            ListQ[sol], sol,
            True, {}
        ];

        flowVals = If[ruleList === {}, <||>,
            Association @ Cases[ruleList,
                HoldPattern[j[a_Integer, b_Integer] -> v_?NumericQ] :>
                    ({a, b} -> Abs[N[v]])]
        ];
        jMax = If[Length[flowVals] === 0, 1, Max[Values[flowVals], 1]];

        thicknessFor[edge_] := With[{key = {edge[[1]], edge[[2]]}},
            If[thickByJ && KeyExistsQ[flowVals, key],
                AbsoluteThickness[1 + 5 flowVals[key]/jMax],
                AbsoluteThickness[1.2]]];

        edgeStyleRules = Join[
            Map[Function[e,
                e -> Directive[RGBColor[0.20, 0.40, 0.75], thicknessFor[e]]],
                integerEdges],
            Map[Function[e,
                e -> Directive[GrayLevel[0.65], Dashed, AbsoluteThickness[0.6]]],
                equivEdges]
        ];

        roundColor[r_] := Blend[
            {RGBColor[0.85, 0.25, 0.25], RGBColor[0.25, 0.35, 0.85]},
            If[rounds === 1, 0.5, (r - 1)/(rounds - 1)]];

        vertexStyles = Flatten @ Table[
            encode[r, p, l] -> If[colorByR, roundColor[r], GrayLevel[0.5]],
            {l, layers}, {r, rounds}, {p, npr}];

        Graph3D[
            vertices,
            Join[integerEdges, equivEdges],
            VertexCoordinates -> vertexCoords,
            VertexStyle       -> vertexStyles,
            VertexSize        -> 0.18,
            VertexLabels      -> If[showLabels,
                                    Placed["Name", Center],
                                    None],
            EdgeStyle         -> edgeStyleRules,
            PlotLabel         -> plotLabel,
            ImageSize         -> imgSize,
            Boxed             -> False
        ]
    ];


(* ::Subsection:: *)
(*Section 1 \[LongDash] Smallest coupled case (2*3*1): a guided tour*)


(* ::Text:: *)
(*Two rounds of three positions, one layer. The smallest scenario where coupling has anything to do: tangential forward edges 1->2 and 4->5 share the physical segment "position 1->2", so the package rewrites either flow as their sum.*)


(* ::Subsubsection:: *)
(*1.1 Building the typed scenario*)


(* ::Text:: *)
(*Capability: makeTawafScenario constructs a typed scenario object that downstream kernels accept without further conversion. Construction parameters are stored in scenarioData[s, "Tawaf"] for later builders to read.*)


tawaf2x3 = makeTawafScenario[2, 3, 1];
tawaf2x3System = makeTawafSystem[tawaf2x3];

DescribeOutput[
    "Scenario object and metadata",
    "scenarioQ confirms the typed wrapper. scenarioData[\"Tawaf\"] retrieves the construction parameters that makeTawafSystem reads later.",
    <|
        "scenarioQ"            -> scenarioQ[tawaf2x3],
        "Tawaf metadata"       -> scenarioData[tawaf2x3, "Tawaf"],
        "Hamiltonian metadata" -> scenarioData[tawaf2x3, "Hamiltonian"]
    |>
]


(* ::Subsubsection:: *)
(*1.2 Inspecting structural data*)


(* ::Text:: *)
(*Capability: typed accessors (scenarioData, systemData) and a flattened topology summary. The 2*3*1 case has 6 logical vertices, 5 directed tangential edges (last position of round 2 has no outgoing tangential edge \[LongDash] that is the exit), and an entry-flow total of 100.*)


DescribeOutput[
    "2\[Times]3\[Times]1 scenario summary",
    "TawafScenarioSummary reads scenarioData and systemData and assembles a one-glance view of the construction.",
    TawafScenarioSummary[tawaf2x3, tawaf2x3System]
]


(* ::Subsubsection:: *)
(*1.3 The Black Stone potential*)


(* ::Text:: *)
(*Capability: per-edge Hamiltonian metadata via Hamiltonian["EdgeV"]. The Tawaf builder sets V = -5.0 on every edge originating from a position-1 node (the Black Stone), 0.0 elsewhere. The negative sign is attractive in the Hamilton-Jacobi formulation and biases trajectories toward saluting the Black Stone on every pass.*)


DescribeOutput[
    "EdgeV distribution on the 2\[Times]3\[Times]1 scenario",
    "Two of the five tangential edges originate at a position-1 node (vertex 1 in round 1; vertex 4 in round 2) and carry V = -5.0; the remaining three carry V = 0.0.",
    Module[{ham, edgeV},
        ham   = scenarioData[tawaf2x3, "Hamiltonian"];
        edgeV = ham["EdgeV"];
        <|
            "Total edges with V" -> Length[edgeV],
            "Edges with V == -5" -> Count[Values[edgeV], -5.0],
            "Edges with V == 0"  -> Count[Values[edgeV], 0.0],
            "All EdgeV pairs"    -> edgeV
        |>
    ]
]


(* ::Subsubsection:: *)
(*1.4 Visualising the unrolled topology*)


(* ::Text:: *)
(*Capability: rawNetworkPlot for the physical/logical network, richNetworkPlot for the augmented state-space graph. Both are pre-solve here; identical calls with a solution argument later produce the solved views.*)


DescribeOutput[
    "2\[Times]3\[Times]1 physical network (rawNetworkPlot)",
    "Six logical vertices arranged as two rounds. Entry at vertex 1, exit at vertex 6.",
    rawNetworkPlot[tawaf2x3, tawaf2x3System,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] physical topology",
        ImageSize -> Medium]
]


DescribeOutput[
    "2\[Times]3\[Times]1 augmented infrastructure (richNetworkPlot, structure only)",
    "Augmented road-traffic graph before solving. Anti-parallel arcs separate so both directions are visible.",
    richNetworkPlot[tawaf2x3, tawaf2x3System,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


(* ::Subsubsection:: *)
(*1.5 The coupling rewrite (the unique Tawaf capability)*)


(* ::Text:: *)
(*Capability: makeTawafSystem builds an ordinary mfgSystem via makeSystem, then rewrites EqGeneral and AltOptCond so logical flows on the same physical segment share congestion. The contrast below is the central novelty for the Tawaf scenario family. The "before" column is the raw makeSystem output; the "after" column is the makeTawafSystem rewrite \[LongDash] notice that occurrences of j[1,2] in any equation become j[1,2] + j[4,5], because both flows traverse the same physical "position 1->2" segment.*)


DescribeOutput[
    "EqGeneral before vs after the coupling rewrite",
    "First four Hamilton-Jacobi equations from each side. Where the same physical segment appears, the rewritten equation references the coupled flow sum.",
    TawafCouplingComparison[tawaf2x3, 4]
]


DescribeOutput[
    "All EqGeneral entries from the coupled system",
    "The rewrite is local: only EqGeneral and AltOptCond change. Other system blocks (balance, complementarity, switching margins) keep the logical-flow form deliberately \[LongDash] each logical node has its own balance and complementarity. See the modelling note's 'Scope of coupling' table.",
    TawafCouplingPreview[tawaf2x3System, Length[systemData[tawaf2x3System, "EqGeneral"]]]
]


(* ::Subsubsection:: *)
(*1.6 Solving symbolically*)


(* ::Text:: *)
(*Capability: dnfReduceSystem (the DNF-first default solver) applied to the Tawaf-coupled system. On the 2*3*1 case it returns a fully-determined rule list in well under a second. Note: solveScenario[s] cannot be used here because it dispatches to makeSystem, not makeTawafSystem \[LongDash] solving the uncoupled system would not satisfy the coupled validator below.*)


AbsoluteTiming[tawaf2x3Sol = dnfReduceSystem[tawaf2x3System];]


(* ::Subsubsection:: *)
(*1.7 Validating the solution*)


(* ::Text:: *)
(*Capability: isValidSystemSolution substitutes the solved rules into every block of the structural system and checks consistency within a numerical tolerance. On a Tawaf-rewritten system it validates the post-rewrite EqGeneral and AltOptCond automatically \[LongDash] no Tawaf-specific code path needed.*)


DescribeOutput[
    "2\[Times]3\[Times]1 solution validity",
    "isValidSystemSolution returns True iff every structural constraint holds.",
    isValidSystemSolution[tawaf2x3System, tawaf2x3Sol]
]


(* ::Subsubsection:: *)
(*1.8 Visualising the solved flow field*)


(* ::Text:: *)
(*Capability: the same rawNetworkPlot and richNetworkPlot used for structure also accept a solution argument and render edge thickness, flow values, value-function gradients, and density labels. The custom 3D helix view (tawafHelixPlot, defined above) takes advantage of the unrolled topology: rounds stack vertically so equivalent positions sit on top of each other.*)


DescribeOutput[
    "2\[Times]3\[Times]1 augmented infrastructure with solved flows (richNetworkPlot)",
    "Blue arcs are flow variables j[a,b]; red arcs are transition flows j[r,i,w]. Node colours show u-values on a Blue\[Rule]Red gradient.",
    richNetworkPlot[tawaf2x3, tawaf2x3System, tawaf2x3Sol,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] solved",
        ImageSize -> Large]
]


DescribeOutput[
    "2\[Times]3\[Times]1 flow plot (rawNetworkPlot)",
    "Real network with auxiliary entry/exit; edge labels show solved j-values.",
    rawNetworkPlot[tawaf2x3, tawaf2x3System, tawaf2x3Sol,
        ShowAuxiliaryVertices -> True,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] flow values",
        ImageSize -> Large]
]


DescribeOutput[
    "2\[Times]3\[Times]1 helix view (structure)",
    "3D layout: angle = position, height = (round, position) in flattened sequence. Equivalent (same position, different round) nodes sit directly above each other on the helix, making the unroll visually interpretable.",
    tawafHelixPlot[tawaf2x3, tawaf2x3System, <||>,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] helix (structure)"]
]


DescribeOutput[
    "2\[Times]3\[Times]1 helix view with solved flows",
    "Edge thickness encodes |j[a,b]|; round is encoded by the red\[Rule]blue node gradient.",
    tawafHelixPlot[tawaf2x3, tawaf2x3System, tawaf2x3Sol,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]1 \[LongDash] helix (solved)"]
]


(* ::Subsubsection:: *)
(*1.9 Density extension (opt-in)*)


(* ::Text:: *)
(*Capability: opt into the density-dependent cost model with "Density" -> True. The scenario builder sets a strictly-negative baseline V on every non-boundary edge ("BaselinePotential" option, default -1.0); the Black-Stone -5 discount is added on top. makeTawafSystem then augments the system with a per-physical-edge density family m[{a,b}] (one symbol per physical undirected edge, keyed by the lex-smallest logical undirected edge in the cohort), an inequality block IneqMs (m > 0), and a consistency block EqDensityFlow whose entries are jPhys^2 + 2 V m + 1 == 0 \[LongDash] derived from the critical-congestion HJB H = 0 with g(m) = -1/(2m). After solving the j/u part of the system, tawafDensities[sys, sol] derives m for each physical edge by inverting the linear consistency equation.*)


tawaf2x3D = makeTawafScenario[2, 3, 1, "Density" -> True];
tawaf2x3DSystem = makeTawafSystem[tawaf2x3D];

DescribeOutput[
    "2\[Times]3\[Times]1 density scenario \[LongDash] EdgeV with strictly-negative baseline",
    "Tangential edges leaving position 1 carry V = baseline + (-5) = -6.0; all other non-boundary edges carry V = baseline = -1.0 (required for the consistency relation to be well-defined).",
    Module[{ham, edgeV},
        ham   = scenarioData[tawaf2x3D, "Hamiltonian"];
        edgeV = ham["EdgeV"];
        <|
            "Tawaf metadata"     -> scenarioData[tawaf2x3D, "Tawaf"],
            "EdgeV"              -> edgeV,
            "Edges with V == -6" -> Count[Values[edgeV], -6.0],
            "Edges with V == -1" -> Count[Values[edgeV], -1.0]
        |>
    ]
]


DescribeOutput[
    "2\[Times]3\[Times]1 density system \[LongDash] new symbolic blocks",
    "Three m's (one per physical undirected edge: cohort {1,2}\[LeftRightArrow]{4,5} canonical = {1,2}; cohort {2,3}\[LeftRightArrow]{5,6} canonical = {2,3}; round-boundary singleton {3,4}). Three IneqMs and three EqDensityFlow equations \[LongDash] one per physical edge by construction.",
    <|
        "Ms"            -> systemData[tawaf2x3DSystem, "Ms"],
        "IneqMs"        -> systemData[tawaf2x3DSystem, "IneqMs"],
        "EqDensityFlow" -> systemData[tawaf2x3DSystem, "EqDensityFlow"]
    |>
]


AbsoluteTiming[tawaf2x3DSol = dnfReduceSystem[tawaf2x3DSystem];]


DescribeOutput[
    "2\[Times]3\[Times]1 density solution \[LongDash] derived per-edge densities",
    "tawafDensities post-processes the j/u solution and returns m[{a,b}] -> value rules by solving the linear consistency equation jPhys^2 + 2 V m + 1 == 0 for each physical edge. All m values must be strictly positive (consistent with IneqMs).",
    Module[{mRules},
        mRules = tawafDensities[tawaf2x3DSystem, tawaf2x3DSol];
        <|
            "j/u solution head" -> Head[tawaf2x3DSol],
            "j/u solution length" -> Length[tawaf2x3DSol],
            "Derived densities" -> mRules,
            "All m > 0?" -> AllTrue[Values[mRules], NumericQ[#] && # > 0 &]
        |>
    ]
]


(* ::Subsection:: *)
(*Section 2 \[LongDash] Mid case 3*4*1: scaling the coupling*)


(* ::Text:: *)
(*Twelve logical nodes; coupling groups have up to 3 logical flows per physical segment (one per round). Each forward tangential segment is shared across all three rounds. Same capabilities as Section 1, exercised at a slightly larger scale.*)


tawaf3x4 = makeTawafScenario[3, 4, 1];
tawaf3x4System = makeTawafSystem[tawaf3x4];

DescribeOutput[
    "3\[Times]4\[Times]1 scenario summary",
    "Three rounds of four positions. Each forward tangential segment (1\[Rule]2, 2\[Rule]3, 3\[Rule]4, 4\[Rule]1) is shared across all three rounds.",
    TawafScenarioSummary[tawaf3x4, tawaf3x4System]
]


DescribeOutput[
    "3\[Times]4\[Times]1 augmented infrastructure (structure only)",
    "Pre-solve augmented graph. Compare visual density to the 2\[Times]3 case to see the cost-of-scaling.",
    richNetworkPlot[tawaf3x4, tawaf3x4System,
        PlotLabel -> "Tawaf 3\[Times]4\[Times]1 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


(* This solve may take longer than the 2\[Times]3 case; evaluate explicitly. *)
AbsoluteTiming[tawaf3x4Sol = dnfReduceSystem[tawaf3x4System];]


DescribeOutput[
    "3\[Times]4\[Times]1 solution validity",
    "Same isValidSystemSolution call; the Tawaf rewrite is transparent to the validator.",
    isValidSystemSolution[tawaf3x4System, tawaf3x4Sol]
]


DescribeOutput[
    "3\[Times]4\[Times]1 augmented infrastructure with solved flows",
    "u-value gradient and flow-magnitude edge thickness reveal where congestion concentrates.",
    richNetworkPlot[tawaf3x4, tawaf3x4System, tawaf3x4Sol,
        PlotLabel -> "Tawaf 3\[Times]4\[Times]1 \[LongDash] solved",
        ImageSize -> Large]
]


DescribeOutput[
    "3\[Times]4\[Times]1 helix view with solved flows",
    "Three rounds, four positions per round. Same-position nodes from successive rounds stack directly above each other along the helix.",
    tawafHelixPlot[tawaf3x4, tawaf3x4System, tawaf3x4Sol,
        PlotLabel -> "Tawaf 3\[Times]4\[Times]1 \[LongDash] helix (solved)"]
]


(* ::Subsection:: *)
(*Section 3 \[LongDash] Multi-layer 2*3*2: radial coupling*)


(* ::Text:: *)
(*Twelve nodes split into two concentric rings (one per layer); radial edges connect adjacent layers in both directions. Demonstrates that the same builder pattern extends to multi-lane topologies. The solve is wrapped in TimeConstrained because the multi-layer search space is larger.*)


tawaf2x3x2 = makeTawafScenario[2, 3, 2];
tawaf2x3x2System = makeTawafSystem[tawaf2x3x2];

DescribeOutput[
    "2\[Times]3\[Times]2 scenario summary",
    "Two layers; radial edges connect each (round, position) across layers in both directions. Entry flow is split uniformly across layers.",
    TawafScenarioSummary[tawaf2x3x2, tawaf2x3x2System]
]


DescribeOutput[
    "2\[Times]3\[Times]2 augmented infrastructure (structure only)",
    "Pre-solve augmented graph for the multi-layer case. The radial edges are visible as the additional connections between concentric clusters.",
    richNetworkPlot[tawaf2x3x2, tawaf2x3x2System,
        PlotLabel -> "Tawaf 2\[Times]3\[Times]2 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


AbsoluteTiming[tawaf2x3x2Sol = TimeConstrained[dnfReduceSystem[tawaf2x3x2System], 120, $TimedOut];]


(* If the solve completes in time, render the solved augmented plot. *)
If[AssociationQ[tawaf2x3x2Sol] || ListQ[tawaf2x3x2Sol],
    DescribeOutput[
        "2\[Times]3\[Times]2 augmented infrastructure with solved flows",
        "Multi-layer solved system. Radial coupling distinguishes outward vs inward.",
        richNetworkPlot[tawaf2x3x2, tawaf2x3x2System, tawaf2x3x2Sol,
            PlotLabel -> "Tawaf 2\[Times]3\[Times]2 \[LongDash] solved",
            ImageSize -> Large]
    ],
    Print["2\[Times]3\[Times]2 solve did not return a solution (",
        tawaf2x3x2Sol, "); skipping solved plot."]
]


DescribeOutput[
    "2\[Times]3\[Times]2 helix view (concentric helices for layers)",
    "Two concentric helices (one per layer). Radial edges connect adjacent layers at matching (round, position).",
    tawafHelixPlot[tawaf2x3x2, tawaf2x3x2System,
        If[AssociationQ[tawaf2x3x2Sol] || ListQ[tawaf2x3x2Sol],
            tawaf2x3x2Sol, <||>],
        PlotLabel -> "Tawaf 2\[Times]3\[Times]2 \[LongDash] helix"]
]


(* ::Subsection:: *)
(*Section 4 \[LongDash] Canonical Tawaf 7*8*1: structure at scale*)


(* ::Text:: *)
(*Mirrors the canonical Tawaf reference: seven counter-clockwise circumambulations of an 8-station circuit. Fifty-six logical nodes; each forward segment is shared by 7 rounds. The full symbolic solve is expensive; the cell that performs it is commented out by default. Structure-only views complete in seconds and showcase how the package's visualisation surface scales.*)


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
    richNetworkPlot[tawaf7x8, tawaf7x8System,
        PlotLabel -> "Tawaf 7\[Times]8\[Times]1 \[LongDash] augmented (structure)",
        ShowBoundaryValues -> False,
        ImageSize -> Large]
]


DescribeOutput[
    "7\[Times]8\[Times]1 helix view (structure only)",
    "Canonical 7-turn helix.",
    tawafHelixPlot[tawaf7x8, tawaf7x8System, <||>,
        PlotLabel -> "Tawaf 7\[Times]8\[Times]1 \[LongDash] helix (structure)"]
]


(* Optional expensive solve. Generous timeout; abort manually if needed. *)
 AbsoluteTiming[tawaf7x8Sol = TimeConstrained[dnfReduceSystem[tawaf7x8System], 600, $TimedOut];] 


DescribeOutput[
    "7\[Times]8\[Times]1 helix view (solved)",
    "Canonical 7-turn helix with solved flows.",
    tawafHelixPlot[tawaf7x8, tawaf7x8System, tawaf7x8Sol,
        PlotLabel -> "Tawaf 7\[Times]8\[Times]1 \[LongDash] helix (solved)",
        "ShowEquivalenceLinks" -> False]
]


(* ::Subsection:: *)
(*Section 5 \[LongDash] Canonical Tawaf 7*8*9 : structure at scale*)


tawaf7x8x9 = makeTawafScenario[7, 8, 9];
tawaf7x8x9System = makeTawafSystem[tawaf7x8x9];

DescribeOutput[
    "7\[Times]8\[Times]9 scenario summary",
    "Tawaf parameters. ",
    TawafScenarioSummary[tawaf7x8x9, tawaf7x8x9System]
]


DescribeOutput[
    "7\[Times]8\[Times]9 helix view (structure only)",
    "7-turn 9-lane helix.",
    tawafHelixPlot[tawaf7x8x9, tawaf7x8x9System, <||>,
        PlotLabel -> "Tawaf 7\[Times]8\[Times]9 \[LongDash] helix (structure)",
        "ShowEquivalenceLinks" -> False]
]


(* ::Subsection:: *)
(*Summary \[LongDash] capabilities demonstrated*)


(* ::Text:: *)
(*By following one scenario family through this document the reader has now seen, in concrete form, every layer of the MFGraphs package surface that the core scenario-kernel phase exposes:*)


(* ::Text:: *)
(*Typed scenario kernel. makeTawafScenario, scenarioQ, scenarioData were demonstrated for one-, two-, and three-dimensional scenarios. Construction parameters survive intact in scenarioData[s, "Tawaf"] for downstream builders to consume.*)


(* ::Text:: *)
(*Specialised builder pattern. makeTawafSystem layers a scenario-specific symbolic rewrite on top of the generic makeSystem output, modifying only the blocks where physical-edge sharing matters (EqGeneral, AltOptCond) and leaving the rest in logical form. The before/after comparison in section 1.5 makes the rewrite visible.*)


(* ::Text:: *)
(*Symbolic system construction and solving. makeSystem, dnfReduceSystem, and isValidSystemSolution apply unchanged to the Tawaf-rewritten system. No code path knows about Tawaf. (This workbook calls dnfReduceSystem directly on the Tawaf-coupled system because solveScenario currently builds the system via makeSystem rather than makeTawafSystem.)*)


(* ::Text:: *)
(*Visualisation surface. rawNetworkPlot and richNetworkPlot accept any mfgSystem and (optionally) a solution; the Tawaf workbook combines them with a scenario-specific 3D helix view to show how a custom plot can sit naturally on top of the shared graphics primitives.*)


(* ::Text:: *)
(*Scaling. The pipeline ran end-to-end on 6, 12, and 12 logical nodes and rendered structure for the canonical 56-node case. The expensive solve at 56 nodes is left as an optional cell.*)


(* ::Text:: *)
(*Pointers: package-side API in CLAUDE.md \[Section] "Tawaf scenario builder" and the auto-generated API_REFERENCE.md. Modelling context (what the scenario represents, why the construction is shaped this way) in docs/research/notes/tawaf-model.md.*)
