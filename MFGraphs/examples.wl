(* Wolfram Language package *)
(* Scenario constructors and named-example registry for MFGraphs.

   Direct constructors — topology clear from arguments:
     gridScenario[{n}, entries, exits]          chain, vertices 1..n
     gridScenario[{r,c}, entries, exits]        grid, vertices 1..r*c row-major
     cycleScenario[n, entries, exits]           cycle connections 1-2-...-n-1
     graphScenario[graph, entries, exits]       any WL Graph object; integer vertex labels required
     amScenario[vl, am, entries, exits]         explicit vertices list + adjacency matrix; integer vertex labels required

   Named examples (benchmark registry):
     getExampleScenario[3, {{1,80}}, {{3,0}}]             chain-3 scenario
     getExampleScenario[12, {{1,100}}, {{4,0}}]           attraction-4 scenario

   All constructors accept optional trailing args: sc, alpha, V, g.
   Hamiltonian defaults are supplied by makeScenario.
   Edges are network connections; topology is symmetrized before system construction,
   and restrictions on movement should be encoded with switching costs such as Infinity.
   Numeric benchmark defaults are in Scripts/BenchmarkHelpers.wls. *)

BeginPackage["examples`", {"primitives`", "scenarioTools`"}];

(* --- Public API declarations --- *)

gridScenario::usage =
"gridScenario[dims, entries, exits] creates a scenario on GridGraph[dims] connections. \
{n} gives a chain with vertices 1..n; {r,c} gives an r\[Times]c grid with vertices 1..r*c (row-major). \
Topology is symmetrized into undirected network edges before system construction. \
Optional: sc (switching costs, default {}), alpha, V, g (Hamiltonian defaults are supplied by makeScenario; V/G are preserved but not applied by current system construction).";

cycleScenario::usage =
"cycleScenario[n, entries, exits] creates a scenario on n-cycle connections, vertices 1..n. \
Topology is symmetrized into undirected network edges before system construction. Optional: sc, alpha, V, g.";

graphScenario::usage =
"graphScenario[graph, entries, exits] creates a scenario from any WL Graph object. \
Graph edges are treated as network connections and symmetrized before system construction; \
movement restrictions should be encoded with switching costs such as Infinity. \
Vertex labels must be positive integers; non-integer labels are not accepted (makeScenario \
returns Failure). Optional: sc, alpha, V, g.";

amScenario::usage =
"amScenario[vl, am, entries, exits] creates a scenario from an explicit vertices list vl \
and adjacency matrix am. Vertex labels in vl must be positive integers; non-integer labels \
are not accepted. The adjacency matrix is treated as network connections and symmetrized before \
system construction. Optional: sc, alpha, V, g.";

getExampleScenario::usage =
"getExampleScenario[n] returns a 6-arg factory Function[{entries,exits,sc,alpha,V,g}, scenario[...]] \
for built-in example n (see listExampleScenarios[] for valid keys). Topology is baked in; all parameters are caller-supplied. \
getExampleScenario[n, entries, exits] calls the factory with canonical defaults: \
  sc=Automatic resolves via $CaseDefaultSC (falls back to {}), \
  alpha=1 (critical congestion), V=0, g=0 (Hamiltonian parameters passed through but not yet applied in system construction). \
getExampleScenario[n, entries, exits, sc] overrides switching costs (pass {} for none). \
getExampleScenario[n, entries, exits, sc, alpha] also overrides alpha. \
getExampleScenario[n, entries, exits, sc, alpha, V] also overrides V. \
getExampleScenario[n, entries, exits, sc, alpha, V, g] overrides all parameters. \
entries={{vertex,flow},...}; exits={{vertex,cost},...}; sc={{i,k,j,cost},...} or Association with 3-tuple keys. \
Returns $Failed for unknown keys. \
Numeric keys (3, 11-18, 20-23, 27, 104, 105) are retained for backward compatibility with \
older scripts; new code should prefer named keys (\"Jamaratv9\", \"Grid0303\", \"Camilli 2015 simple\", etc.). \
Use listExampleScenarios[] to discover all registered keys.";

getExampleScenarioMetadata::usage =
"getExampleScenarioMetadata[key] returns metadata for a built-in example scenario, \
including source/provenance details when available. Returns $Failed for unknown keys.";

listExampleScenarios::usage =
"listExampleScenarios[] returns the sorted list of all keys recognised by getExampleScenario \
(both numeric legacy keys and named keys). Use this for discovery and for iterating over \
the registry in tests.";

Begin["`Private`"];

scenarioFromGraph::usage =
"scenarioFromGraph[graph, entries, exits, sc, alpha, V, g] builds a scenario from a Graph and Hamiltonian parameters.";

makeGridFactory::usage =
"makeGridFactory[dims] returns a built-in example factory backed by gridScenario[dims, ...].";

makeCycleFactory::usage =
"makeCycleFactory[n] returns a built-in example factory backed by cycleScenario[n, ...].";

makeAmFactory::usage =
"makeAmFactory[vl, am] returns a built-in example factory backed by amScenario[vl, am, ...].";

edgePairsToAdjacency::usage =
"edgePairsToAdjacency[vertices, pairs] builds an adjacency matrix from directed edge pairs.";

(* --- Shared topology constants --- *)

$Y1In2OutAM    = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
$Attraction4AM = {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}};

$Camilli2015SourcePath =
    "docs/research/papers/Camilli et al. - 2015 - A model problem for Mean Field Games on networks.pdf";

$Achdou2023SourcePath =
    "docs/research/papers/Achdou et al. - 2023 - First order Mean Field Games on networks.pdf";

$AchdouJunctionVertices = Range[5];
$AchdouJunctionCoordinates = <|
    1 -> {0, 0},
    2 -> {0, 1},
    3 -> {1, 0},
    4 -> {0, -1},
    5 -> {-1, 0}
|>;
$AchdouJunctionEdges = {
    {1, 2},
    {1, 3},
    {1, 4},
    {1, 5}
};

$CamilliSimpleVertices = Range[4];
$CamilliSimpleCoordinates = <|
    1 -> {0, 0},
    2 -> {0, 1},
    3 -> {-1, 2},
    4 -> {1, 2}
|>;
$CamilliSimpleEdges = {
    {1, 2},
    {2, 3},
    {2, 4},
    {3, 4}
};

$CamilliGeneralVertices = Range[17];
$CamilliGeneralCoordinates = <|
    1 -> {0, 0},
    2 -> {0, 1/2},
    3 -> {-1, 1},
    4 -> {-1/2, 3/2},
    5 -> {-1, 2},
    6 -> {-1/2, 5/2},
    7 -> {-1, 3},
    8 -> {-3/2, 7/2},
    9 -> {-2, 3},
    10 -> {-5/2, 5/2},
    11 -> {-3/2, 5/2},
    12 -> {-2, 2},
    13 -> {-3/2, 3/2},
    14 -> {1/2, 3/2},
    15 -> {1, 2},
    16 -> {3/2, 3/2},
    17 -> {1, 1}
|>;
$CamilliGeneralEdges = {
    {1, 2},
    {2, 3},
    {2, 14},
    {2, 17},
    {3, 4},
    {3, 13},
    {4, 5},
    {4, 6},
    {5, 6},
    {5, 11},
    {5, 13},
    {6, 7},
    {7, 8},
    {7, 11},
    {8, 9},
    {9, 10},
    {9, 12},
    {10, 11},
    {11, 12},
    {12, 13},
    {14, 15},
    {15, 16},
    {16, 17},
    {14, 17}
};

$CamilliGeneralImportedEdges = DeleteCases[$CamilliGeneralEdges, {2, 14} | {5, 11}];

(* Canonical switching costs looked up by getExampleScenario when sc=Automatic. *)
$CaseDefaultSC = <|
    11 -> {{1,2,4,1},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
           {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
           {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}},
    13 -> {{1,2,4,1},{1,3,4,1},{4,2,1,1},{4,3,1,1}},
    14 -> {{1,2,3,2},{3,2,1,1}},
    15 -> {{2,1,3,2},{3,1,2,1}},
    16 -> {{1,3,2,1}},
    17 -> {{1,2,3,2},{3,2,1,1}},
    18 -> {{1,2,3,2},{3,2,1,1}},
    "Braess split"       -> {{1,2,4,1},{5,7,8,1}},
    "Braess congest"     -> {{1,2,4,1},{4,6,7,1}},
    "Big Braess split"   -> {{1,3,6,1},{5,7,10,1}},
    "Big Braess congest" -> {{1,3,5,1},{5,7,9,1}},
    "Paper example"      -> {{1,2,3,2},{3,2,1,1},{2,3,4,3},{4,3,2,1}},
    (* SC violates triangle inequality — infeasible by design *)
    "Inconsistent Y shortcut" ->
        {{1,2,3,5},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}},
    "Inconsistent attraction shortcut" ->
        {{1,2,4,5},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
         {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
         {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}}
|>;

$ExamplesDefaultHamiltonianValue = Automatic;

exampleHamiltonianSpec[alpha_, V_, g_] :=
    Association @ Select[
        {"Alpha" -> alpha, "V" -> V, "G" -> g},
        Last[#] =!= $ExamplesDefaultHamiltonianValue &
    ];

scenarioFromGraph[graph_, entries_, exits_, sc_, alpha_, V_, g_] :=
    makeScenario[<|
        "Model" -> <|
            "Graph"                           -> graph,
            "Entries"     -> entries,
            "Exits" -> exits,
            "Switching"                 -> sc
        |>,
        "Hamiltonian" -> exampleHamiltonianSpec[alpha, V, g]
    |>];

gridScenario[dims_List, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultHamiltonianValue,
        V_     : $ExamplesDefaultHamiltonianValue,
        g_     : $ExamplesDefaultHamiltonianValue] :=
    scenarioFromGraph[GridGraph[dims, DirectedEdges -> True], entries, exits, sc, alpha, V, g];

cycleScenario[n_Integer, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultHamiltonianValue,
        V_     : $ExamplesDefaultHamiltonianValue,
        g_     : $ExamplesDefaultHamiltonianValue] :=
    scenarioFromGraph[CycleGraph[n, DirectedEdges -> True], entries, exits, sc, alpha, V, g];

graphScenario[graph_, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultHamiltonianValue,
        V_     : $ExamplesDefaultHamiltonianValue,
        g_     : $ExamplesDefaultHamiltonianValue] :=
    scenarioFromGraph[graph, entries, exits, sc, alpha, V, g];

amScenario[vl_, am_, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultHamiltonianValue,
        V_     : $ExamplesDefaultHamiltonianValue,
        g_     : $ExamplesDefaultHamiltonianValue] :=
    makeScenario[<|
        "Model" -> <|
            "Vertices"                   -> vl,
            "Adjacency"                -> am,
            "Entries"     -> entries,
            "Exits" -> exits,
            "Switching"                 -> sc
        |>,
        "Hamiltonian" -> exampleHamiltonianSpec[alpha, V, g]
    |>];

(* --- Private factory helpers — thin wrappers used by $ExampleScenarios --- *)

makeGridFactory[dims_List]  := gridScenario[dims, ##]&;
makeCycleFactory[n_Integer] := cycleScenario[n, ##]&;
makeAmFactory[vl_, am_]     := amScenario[vl, am, ##]&;

edgePairsToAdjacency[vertices_List, pairs_List] :=
    Module[{index = AssociationThread[vertices -> Range[Length[vertices]]]},
        Normal @ SparseArray[
            ({index[#[[1]]], index[#[[2]]]} -> 1) & /@ pairs,
            {Length[vertices], Length[vertices]}
        ]
    ];

$ExampleScenarioMetadata = <|
    "Camilli 2015 simple" -> <|
        "Title" -> "A model problem for Mean Field Games on networks",
        "SourcePaperPath" -> $Camilli2015SourcePath,
        "Source" -> "Camilli, Carlini, and Marchi (2015), Figure 1 and Section 5.1",
        "PackageInterpretation" -> "Stationary analog only: the meeting point v0 is represented as the unique exit vertex with cost 0; the original time-dependent stochastic MFG and heat-equation fixed-point solver are not implemented here.",
        "VertexCoordinates" -> $CamilliSimpleCoordinates,
        "MeetingVertex" -> 1,
        "ExitDefault" -> {{1, 0}},
        "EdgeList" -> $CamilliSimpleEdges,
        "DeclaredVertexCount" -> 4,
        "DeclaredEdgeCount" -> 4,
        "TimeDependentData" -> <|
            "t0" -> 0.5,
            "Tmax" -> 10,
            "Theta" -> 0.5,
            "CostFormula" -> "c_T(s) = 0.1 max(s - t0, 0) + 0.1 max(T - s, 0)",
            "InitialDensityFormula" -> "m0(x) = g(x) / Integral_Gamma g(y) dy, with g(x) = |x| restricted to the graph",
            "ReportedNumericalResult" -> <|"T2" -> 5.62|>
        |>,
        "RecommendedStationaryEntries" -> <|
            "Description" -> "Split total inflow across the two non-exit upper vertices to reflect g(x)=|x| concentrating mass away from v0.",
            "Vertices" -> {3, 4},
            "ExampleForTotalInflow100" -> {{3, 50}, {4, 50}}
        |>,
        "ImportNotes" -> {
            "Edges are network connections; MFGraphs symmetrizes the topology before system construction.",
            "Switching costs default to zero, matching the paper graph as connections rather than directed movement restrictions."
        }
    |>,
    "Camilli 2015 general" -> <|
        "Title" -> "A model problem for Mean Field Games on networks",
        "SourcePaperPath" -> $Camilli2015SourcePath,
        "Source" -> "Camilli, Carlini, and Marchi (2015), Figure 3 and Section 5.2",
        "PackageInterpretation" -> "Stationary analog only: the meeting point v0 is represented as the unique exit vertex with cost 0; the original time-dependent stochastic MFG and heat-equation fixed-point solver are not implemented here.",
        "VertexCoordinates" -> $CamilliGeneralCoordinates,
        "MeetingVertex" -> 1,
        "ExitDefault" -> {{1, 0}},
        "EdgeList" -> $CamilliGeneralImportedEdges,
        "DeclaredVertexCount" -> 17,
        "DeclaredEdgeCount" -> 22,
        "TimeDependentData" -> <|
            "t0" -> 0.5,
            "Tmax" -> 25,
            "Theta" -> 0.7,
            "CostFormula" -> "c_T(s) = 0.1 max(s - t0, 0) + 0.1 max(T - s, 0)",
            "InitialDensityFormula" -> "m0(x) = g(x) / Integral_Gamma g(y) dy, with g(x)=max(1/4 - |x - p1|^2, 0) + max(1/4 - |x - p2|^2, 0), p1=(1,3/2), p2=(-3/2,3)",
            "ReportedNumericalResult" -> <|"T" -> 23.99, "ErrorEhT" -> 2.35*10^-2|>
        |>,
        "RecommendedStationaryEntries" -> <|
            "Description" -> "Split inflow around neighborhoods of p1=(1,3/2) and p2=(-3/2,3), reflecting the two-bump initial density.",
            "VerticesNearP1" -> {14, 15, 16, 17},
            "VerticesNearP2" -> {7, 8, 9, 11},
            "ExampleForTotalInflow100" -> {{14, 25}, {15, 25}, {7, 25}, {9, 25}}
        |>,
        "ImportNotes" -> {
            "The paper states 17 vertices and 22 edges for Figure 3; this registry preserves that declared topology.",
            "Two visually ambiguous drawn segments, {2,14} and {5,11}, are omitted from the imported edge set to match the declared edge count.",
            "Edges are network connections; MFGraphs symmetrizes the topology before system construction.",
            "Switching costs default to zero, matching the paper graph as connections rather than directed movement restrictions."
        }
    |>,
    3 -> <|
        "Title" -> "Linear chain, 3 vertices",
        "Description" -> "Three-vertex chain 1-2-3 built via gridScenario[{3}, ...]. The smallest non-trivial scenario; useful as a smoke test."
    |>,
    11 -> <|
        "Title" -> "Attraction 4-vertex diamond (case 11)",
        "Description" -> "Diamond topology 1-{2,3}-4 with shared interior vertices and the canonical attraction switching costs."
    |>,
    12 -> <|
        "Title" -> "Attraction 4-vertex diamond (case 12)",
        "Description" -> "Same diamond topology as case 11; reserved for variant boundary data."
    |>,
    13 -> <|
        "Title" -> "Attraction 4-vertex diamond, case 13 variant",
        "Description" -> "1->{2,3}->4 with reduced edges; lighter switching cost set than case 11."
    |>,
    14 -> <|
        "Title" -> "Triangle 3-cycle 1->2->3->1 (case 14)",
        "Description" -> "Three-vertex directed cycle with the canonical Y-shortcut switching costs."
    |>,
    15 -> <|
        "Title" -> "3-vertex misc (case 15)",
        "Description" -> "Two-edge in-tree 2->1, 3->1 used to probe the EqEntryIn/EqExit boundary blocks."
    |>,
    16 -> <|
        "Title" -> "3-vertex chain (case 16)",
        "Description" -> "1->2->3 chain probing single-shortcut switching at vertex 2."
    |>,
    17 -> <|
        "Title" -> "3-vertex chain (case 17)",
        "Description" -> "Same topology as case 16 with a different SC default."
    |>,
    18 -> <|
        "Title" -> "3-vertex chain (case 18)",
        "Description" -> "Same topology as cases 16/17, alternative SC profile."
    |>,
    20 -> <|
        "Title" -> "9-vertex multi-entrance/multi-exit (case 20)",
        "Description" -> "Two-entrance, three-exit topology used as an early Jamarat-like benchmark."
    |>,
    21 -> <|
        "Title" -> "12-vertex multi-entrance/multi-exit (case 21)",
        "Description" -> "Larger benchmark variant of the multi-entrance/multi-exit family."
    |>,
    22 -> <|
        "Title" -> "7-vertex two-entrance/two-exit (case 22)",
        "Description" -> "Compact multi-entrance benchmark."
    |>,
    23 -> <|
        "Title" -> "6-vertex two-entrance/two-exit (case 23)",
        "Description" -> "Smallest multi-entrance benchmark in the registry."
    |>,
    27 -> <|
        "Title" -> "2-vertex undirected edge (case 27)",
        "Description" -> "Trivial scenario: a single undirected edge. Useful for unit-test scaffolding."
    |>,
    104 -> <|
        "Title" -> "Triangle 3-cycle alias (case 104)",
        "Description" -> "Alias for the triangle 1->2->3->1 scenario; same factory as case 14."
    |>,
    105 -> <|
        "Title" -> "3-vertex chain with two exits (case 105)",
        "Description" -> "Alias for the named scenario \"chain with two exits\"."
    |>,
    "triangle with two exits" -> <|
        "Title" -> "Triangle with two exits",
        "Description" -> "Three-vertex directed cycle exposed under a descriptive name; same factory as case 14."
    |>,
    "chain with two exits" -> <|
        "Title" -> "3-vertex chain with two exits",
        "Description" -> "Three-vertex chain 1->2->3 with exits at vertices 2 and 3. Common parametric-solution test case."
    |>,
    "Jamaratv9" -> <|
        "Title" -> "Jamarat v9: 9-vertex two-entrance / three-exit",
        "Description" -> "Multi-entrance, multi-exit infrastructure scenario used in the Jamarat workbook (MFGraphs/Jamarat.wl). Two entries (vertices 1 and 2) feed three downstream exits (vertices 7, 8, 9)."
    |>,
    "Braess split" -> <|
        "Title" -> "Braess paradox: split variant",
        "Description" -> "8-vertex Braess-style graph with parallel routes; switching cost {1,2,4,1} biases the split."
    |>,
    "Braess congest" -> <|
        "Title" -> "Braess paradox: congested variant",
        "Description" -> "7-vertex Braess-style graph that exhibits the congestion paradox under the canonical switching cost."
    |>,
    "New Braess" -> <|
        "Title" -> "Braess with custom congestion field",
        "Description" -> "6-vertex Braess variant carrying an extra \"a\" congestion field on edges 1->4 and 3->6 (j/100). Demonstrates per-edge Hamiltonian extensions outside the standard Alpha/V/G triple."
    |>,
    "Big Braess split" -> <|
        "Title" -> "Larger Braess: split variant",
        "Description" -> "10-vertex extension of the Braess split topology."
    |>,
    "Big Braess congest" -> <|
        "Title" -> "Larger Braess: congested variant",
        "Description" -> "9-vertex extension of the Braess congested topology."
    |>,
    "HRF Scenario 1" -> <|
        "Title" -> "HRF Scenario 1",
        "Description" -> "10-vertex multi-route benchmark from the HRF case study."
    |>,
    "Paper example" -> <|
        "Title" -> "Paper example: 4-vertex chain",
        "Description" -> "Four-vertex chain used as a worked example in early write-ups; shipped with a non-trivial SC profile."
    |>,
    "Inconsistent Y shortcut" -> <|
        "Title" -> "Y-topology with infeasible switching costs",
        "Description" -> "4-vertex Y graph whose canonical SC violates the triangle inequality. Used to validate infeasibility detection."
    |>,
    "Inconsistent attraction shortcut" -> <|
        "Title" -> "Attraction diamond with infeasible switching costs",
        "Description" -> "Attraction-4 topology whose canonical SC violates the triangle inequality. Used to validate infeasibility detection."
    |>,
    "Grid0303" -> <|
        "Title" -> "3x3 grid",
        "Description" -> "GridGraph[{3,3}] connections; vertices 1..9 in row-major order."
    |>,
    "Grid0404" -> <|
        "Title" -> "4x4 grid",
        "Description" -> "GridGraph[{4,4}] connections; vertices 1..16 in row-major order."
    |>,
    "Grid0505" -> <|
        "Title" -> "5x5 grid",
        "Description" -> "GridGraph[{5,5}] connections; vertices 1..25 in row-major order."
    |>,
    "Grid0707" -> <|
        "Title" -> "7x7 grid",
        "Description" -> "GridGraph[{7,7}] connections; vertices 1..49 in row-major order."
    |>,
    "Grid0710" -> <|
        "Title" -> "7x10 grid",
        "Description" -> "GridGraph[{7,10}] connections; vertices 1..70 in row-major order."
    |>,
    "Grid1010" -> <|
        "Title" -> "10x10 grid",
        "Description" -> "GridGraph[{10,10}] connections; vertices 1..100 in row-major order."
    |>,
    "Grid1020" -> <|
        "Title" -> "10x20 grid",
        "Description" -> "GridGraph[{10,20}] connections; vertices 1..200 in row-major order. Largest grid in the registry."
    |>,
    "Achdou 2023 junction" -> <|
        "Title" -> "First order Mean Field Games on networks",
        "Authors" -> {"Yves Achdou", "Paola Mannucci", "Claudio Marchi", "Nicoletta Tchou"},
        "HALId" -> "hal-03729443v3",
        "SourcePaperPath" -> $Achdou2023SourcePath,
        "Source" -> "Achdou, Mannucci, Marchi, and Tchou (2023), deterministic first-order MFGs on junctions",
        "Keywords" -> {"mean field games", "networks", "junctions", "relaxed equilibria", "viscosity solutions", "continuity equation"},
        "ModelClass" -> "Stationary finite graph analog of a first-order deterministic junction model",
        "FiniteAnalogNote" -> "This registry entry is a stationary finite graph analog: the paper studies finite-horizon deterministic relaxed equilibria on junctions, with probability measures over admissible trajectories, mild solutions (u,m), viscosity Hamilton-Jacobi theory, and weak continuity equations. Those formulations are not implemented here.",
        "VertexCoordinates" -> $AchdouJunctionCoordinates,
        "CentralVertex" -> 1,
        "LeafVertices" -> {2, 3, 4, 5},
        "EdgeList" -> $AchdouJunctionEdges,
        "DeclaredVertexCount" -> 5,
        "DeclaredEdgeCount" -> 4,
        "RecommendedStationaryBoundaryExamples" -> <|
            "Description" -> "Metadata-only suggestions for stationary MFGraphs experiments: inject mass at one or more leaves and use the central vertex as the exit.",
            "SingleLeafToCenter" -> <|"Entries" -> {{2, 100}}, "Exits" -> {{1, 0}}|>,
            "BalancedLeavesToCenter" -> <|"Entries" -> {{2, 25}, {3, 25}, {4, 25}, {5, 25}}, "Exits" -> {{1, 0}}|>
        |>,
        "ImportNotes" -> {
            "Edges are finite network connections replacing the paper's semi-infinite junction half-lines.",
            "Entries and exits remain caller-supplied through getExampleScenario[key, entries, exits].",
            "Switching costs default to zero via the existing example registry behavior.",
            "This entry does not implement relaxed equilibria, trajectory measures, viscosity solutions, or finite-horizon first-order MFG solvers."
        }
    |>
|>;

(* --- Scenario registry --- *)

$ExampleScenarios = Association[

    (* ------------------------------------------------------------------ *)
    (* Linear chain, 3 vertices                                           *)
    (* ------------------------------------------------------------------ *)

    3 -> makeGridFactory[{3}],

    (* ------------------------------------------------------------------ *)
    (* Attraction 4-vertex diamond (cases 11, 12, 13)                     *)
    (* ------------------------------------------------------------------ *)

    11 -> makeAmFactory[{1,2,3,4}, $Attraction4AM],
    12 -> makeAmFactory[{1,2,3,4}, $Attraction4AM],
    13 -> makeAmFactory[{1,2,3,4}, {{0,1,1,0},{0,0,0,1},{0,0,0,1},{0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Triangle 3-vertex directed cycle: 1->2->3->1                       *)
    (* ------------------------------------------------------------------ *)

    14                        -> makeCycleFactory[3],
    104                       -> makeCycleFactory[3],
    "triangle with two exits" -> makeCycleFactory[3],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex directed chain with two exits (cases 105, alias)          *)
    (* ------------------------------------------------------------------ *)

    "chain with two exits" -> makeGridFactory[{3}],
    105                    -> makeGridFactory[{3}],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex misc (cases 15–18)                                        *)
    (* ------------------------------------------------------------------ *)

    15 -> makeAmFactory[{1,2,3}, {{0,0,0},{1,0,0},{1,0,0}}],
    16 -> makeAmFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],
    17 -> makeAmFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],
    18 -> makeAmFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Multi-entrance/exit (cases 20–23, "Jamaratv9")                     *)
    (* ------------------------------------------------------------------ *)

    20 -> makeAmFactory[Range[9], {
            {0,0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},
            {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,1,1,1},
            {0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0}}],

    21 -> makeAmFactory[Range[12], {
            {0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,0},
            {0,0,0,1,1,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,0},
            {0,0,0,0,0,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0},
            {0,0,0,0,0,0,0,1,1,1,0,0},{0,0,0,0,0,0,0,0,1,0,1,0},
            {0,0,0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}}],

    22 -> makeAmFactory[Range[7], {
            {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,1,0,0},
            {0,0,0,0,0,1,0},{0,0,0,0,0,1,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}}],

    "Jamaratv9" -> makeAmFactory[Range[9], {
            {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},
            {0,0,0,0,0,1,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
            {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}}],

    23 -> makeAmFactory[Range[6], {
            {0,1,1,0,0,0},{0,0,1,0,0,0},{0,0,0,1,1,0},
            {0,0,0,0,1,1},{0,0,0,0,0,1},{0,0,0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* 2-vertex undirected edge (case 27)                                  *)
    (* ------------------------------------------------------------------ *)

    27 -> makeAmFactory[{1,2}, {{0,1},{1,0}}],

    (* ------------------------------------------------------------------ *)
    (* Braess variants                                                     *)
    (* ------------------------------------------------------------------ *)

    "Braess split" -> makeAmFactory[Range[8], {
            {0,1,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},
            {0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1},
            {0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0}}],

    "Braess congest" -> makeAmFactory[Range[7], {
            {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,0,0,0},
            {0,0,0,0,1,1,0},{0,0,0,0,0,0,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}}],

    (* "New Braess" carries an extra "a" congestion field — not expressible
       via the standard constructors, so defined inline. *)
    "New Braess" -> With[{newBraessAM = {
                {0,1,0,1,0,0},{0,0,1,0,0,0},{0,0,0,0,0,1},
                {0,0,0,0,1,0},{0,0,0,0,0,1},{0,0,0,0,0,0}}},
        Function[{entries, exits, sc, alpha, V, g},
            makeScenario[<|
                "Model" -> <|
                    "Vertices"                   -> Range[6],
                    "Adjacency"                -> newBraessAM,
                    "Entries"     -> entries,
                    "Exits" -> exits,
                    "Switching"                 -> sc,
                    "a" -> Function[{j, edge},
                            Which[
                                edge === DirectedEdge[3,6] || edge === DirectedEdge[1,4], j/100,
                                True, 0]]
                |>,
                "Hamiltonian" -> exampleHamiltonianSpec[alpha, V, g]
            |>]]],

    "Big Braess split" -> makeAmFactory[Range[10], {
            {0,1,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},
            {0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},
            {0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},
            {0,0,0,0,0,0,0,0,0,0}}],

    "Big Braess congest" -> makeAmFactory[Range[9], {
            {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},
            {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
            {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Paper / benchmark cases                                             *)
    (* ------------------------------------------------------------------ *)

    "HRF Scenario 1" -> makeAmFactory[Range[10], {
            {0,1,0,0,0,0,0,0,0,0},{0,0,1,1,0,0,0,0,0,0},{0,0,0,1,1,0,1,0,0,0},
            {0,0,0,0,1,1,1,0,0,0},{0,0,0,0,0,1,1,0,0,0},{0,0,0,0,0,0,1,0,0,0},
            {0,0,0,0,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0}}],

    "Paper example" -> makeGridFactory[{4}],

    "Camilli 2015 simple" -> makeAmFactory[
        $CamilliSimpleVertices,
        edgePairsToAdjacency[$CamilliSimpleVertices, $CamilliSimpleEdges]
    ],

    "Camilli 2015 general" -> makeAmFactory[
        $CamilliGeneralVertices,
        edgePairsToAdjacency[$CamilliGeneralVertices, $CamilliGeneralImportedEdges]
    ],

    "Achdou 2023 junction" -> makeAmFactory[
        $AchdouJunctionVertices,
        edgePairsToAdjacency[$AchdouJunctionVertices, $AchdouJunctionEdges]
    ],

    (* ------------------------------------------------------------------ *)
    (* Inconsistent switching (feature validation — infeasible by design)  *)
    (* ------------------------------------------------------------------ *)

    "Inconsistent Y shortcut"          -> makeAmFactory[{1,2,3,4}, $Y1In2OutAM],
    "Inconsistent attraction shortcut" -> makeAmFactory[{1,2,3,4}, $Attraction4AM],

    (* ------------------------------------------------------------------ *)
    (* Grid cases: GridGraph[{r,c}] connections                           *)
    (* ------------------------------------------------------------------ *)

    "Grid0303" -> makeGridFactory[{3,3}],
    "Grid0404" -> makeGridFactory[{4,4}],
    "Grid0505" -> makeGridFactory[{5,5}],
    "Grid0707" -> makeGridFactory[{7,7}],
    "Grid0710" -> makeGridFactory[{7,10}],
    "Grid1010" -> makeGridFactory[{10,10}],
    "Grid1020" -> makeGridFactory[{10,20}]
];

(* --- Accessors --- *)

getExampleScenario[n_] := Lookup[$ExampleScenarios, n, $Failed];

getExampleScenarioMetadata[n_] := Lookup[$ExampleScenarioMetadata, n, $Failed];

listExampleScenarios[] := SortBy[Keys[$ExampleScenarios], {Head[#] === String, #} &];

(* sc=Automatic resolves to the canonical SC via $CaseDefaultSC, or {} if undefined. *)
getExampleScenario[n_, entries_, exits_,
        sc_    : Automatic,
        alpha_ : $ExamplesDefaultHamiltonianValue,
        V_     : $ExamplesDefaultHamiltonianValue,
        g_     : $ExamplesDefaultHamiltonianValue] :=
    Module[{f = Lookup[$ExampleScenarios, n, $Failed]},
        If[f === $Failed, Return[$Failed]];
        f[entries, exits,
            If[sc === Automatic, Lookup[$CaseDefaultSC, n, {}], sc],
            alpha, V, g]
    ];

End[];

EndPackage[];
