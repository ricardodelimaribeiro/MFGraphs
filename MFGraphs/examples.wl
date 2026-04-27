(* Wolfram Language package *)
(* Scenario constructors and named-example registry for MFGraphs.

   Direct constructors — topology clear from arguments:
     gridScenario[{n}, entries, exits]          chain, vertices 1..n
     gridScenario[{r,c}, entries, exits]        grid, vertices 1..r*c row-major
     cycleScenario[n, entries, exits]           directed n-cycle 1->2->...->n->1
     graphScenario[graph, entries, exits]       any WL directed Graph object; integer vertex labels required
     amScenario[vl, am, entries, exits]         explicit vertices list + adjacency matrix; integer vertex labels required

   Named examples (benchmark registry):
     getExampleScenario[7, {{1,80}}, {{3,0},{4,10}}]      canonical SC for case 7
     getExampleScenario[8, {{1,80}}, {{3,0},{4,10}}, {}]  override: no SC

   All constructors accept optional trailing args: sc, alpha, V, g.
   Defaults: sc={}, alpha=1, V=0, g=Function[z,-1/z] (from $DefaultHamiltonian).
   Numeric benchmark defaults are in Scripts/BenchmarkHelpers.wls. *)

BeginPackage["examples`", {"primitives`", "scenarioTools`"}];

(* --- Public API declarations --- *)

gridScenario::usage =
"gridScenario[dims, entries, exits] creates a scenario on a directed GridGraph[dims]. \
{n} gives a chain with vertices 1..n; {r,c} gives an r\[Times]c grid with vertices 1..r*c (row-major). \
Optional: sc (switching costs, default {}), alpha, V, g (Hamiltonian defaults from $DefaultHamiltonian).";

cycleScenario::usage =
"cycleScenario[n, entries, exits] creates a scenario on a directed n-cycle (1->2->...->n->1), \
vertices 1..n. Optional: sc, alpha, V, g.";

graphScenario::usage =
"graphScenario[graph, entries, exits] creates a scenario from any WL directed Graph object. \
Vertex labels must be positive integers; non-integer labels are not accepted (makeScenario \
returns Failure). Optional: sc, alpha, V, g.";

amScenario::usage =
"amScenario[vl, am, entries, exits] creates a scenario from an explicit vertices list vl \
and adjacency matrix am. Vertex labels in vl must be positive integers; non-integer labels \
are not accepted. Optional: sc, alpha, V, g.";

getExampleScenario::usage =
"getExampleScenario[n] returns a 6-arg factory Function[{entries,exits,sc,alpha,V,g}, scenario[...]] \
for built-in example n. Topology is baked in; all parameters are caller-supplied. \
getExampleScenario[n, entries, exits] calls the factory using the canonical switching costs \
for that case (sc=Automatic resolves via $CaseDefaultSC, defaulting to {} if none defined) \
and standard Hamiltonian defaults (alpha=1, V=0, g=Function[z,-1/z]). \
Additional optional arguments override each default in order: sc, alpha, V, g. \
Pass sc={} explicitly to force no switching costs. \
entries={{vertex,flow},...}, exits={{vertex,cost},...}, sc={{i,k,j,cost},...}. \
Returns $Failed for unknown keys.";

Begin["`Private`"];

scenarioFromGraph::usage =
"scenarioFromGraph[graph, entries, exits, sc, alpha, V, g] builds a scenario from a Graph and Hamiltonian parameters.";

makeGridFactory::usage =
"makeGridFactory[dims] returns a built-in example factory backed by gridScenario[dims, ...].";

makeCycleFactory::usage =
"makeCycleFactory[n] returns a built-in example factory backed by cycleScenario[n, ...].";

makeAmFactory::usage =
"makeAmFactory[vl, am] returns a built-in example factory backed by amScenario[vl, am, ...].";

(* --- Shared topology constants --- *)

$Y1In2OutAM    = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
$Y2In1OutAM    = {{0,1,0,0},{0,0,0,1},{0,1,0,0},{0,0,0,0}};
$Attraction4AM = {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}};

(* Canonical switching costs looked up by getExampleScenario when sc=Automatic. *)
$CaseDefaultSC = <|
    8  -> {{1,2,3,2},{1,2,4,3},{3,2,1,2},{3,2,4,1},{4,2,1,3},{4,2,3,1}},
    10 -> {{1,2,4,2},{1,2,3,3},{3,2,1,2},{3,2,4,1},{4,2,1,3},{4,2,3,1}},
    11 -> {{1,2,4,1},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
           {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
           {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}},
    13 -> {{1,2,4,1},{1,3,4,1},{4,2,1,1},{4,3,1,1}},
    14 -> {{1,2,3,2},{3,2,1,1}},
    15 -> {{2,1,3,2},{3,1,2,1}},
    16 -> {{1,3,2,1}},
    17 -> {{1,2,3,2},{3,2,1,1}},
    18 -> {{1,2,3,2},{3,2,1,1}},
    19 -> {{1,2,3,1},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}},
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

$ExamplesDefaultAlpha = 1;
$ExamplesDefaultV = 0;
$ExamplesDefaultG = Function[z, -1/z];

scenarioFromGraph[graph_, entries_, exits_, sc_, alpha_, V_, g_] :=
    makeScenario[<|
        "Model" -> <|
            "Graph"                           -> graph,
            "Entries"     -> entries,
            "Exits" -> exits,
            "Switching"                 -> sc
        |>,
        "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
    |>];

gridScenario[dims_List, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultAlpha,
        V_     : $ExamplesDefaultV,
        g_     : $ExamplesDefaultG] :=
    scenarioFromGraph[GridGraph[dims, DirectedEdges -> True], entries, exits, sc, alpha, V, g];

cycleScenario[n_Integer, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultAlpha,
        V_     : $ExamplesDefaultV,
        g_     : $ExamplesDefaultG] :=
    scenarioFromGraph[CycleGraph[n, DirectedEdges -> True], entries, exits, sc, alpha, V, g];

graphScenario[graph_, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultAlpha,
        V_     : $ExamplesDefaultV,
        g_     : $ExamplesDefaultG] :=
    scenarioFromGraph[graph, entries, exits, sc, alpha, V, g];

amScenario[vl_, am_, entries_, exits_,
        sc_    : {},
        alpha_ : $ExamplesDefaultAlpha,
        V_     : $ExamplesDefaultV,
        g_     : $ExamplesDefaultG] :=
    makeScenario[<|
        "Model" -> <|
            "Vertices"                   -> vl,
            "Adjacency"                -> am,
            "Entries"     -> entries,
            "Exits" -> exits,
            "Switching"                 -> sc
        |>,
        "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
    |>];

(* --- Private factory helpers — thin wrappers used by $ExampleScenarios --- *)

makeGridFactory[dims_List]  := gridScenario[dims, ##]&;
makeCycleFactory[n_Integer] := cycleScenario[n, ##]&;
makeAmFactory[vl_, am_]     := amScenario[vl, am, ##]&;

(* --- Scenario registry --- *)

$ExampleScenarios = Association[

    (* ------------------------------------------------------------------ *)
    (* Linear chains (cases 1–6): directed GridGraph[{n}]                 *)
    (* ------------------------------------------------------------------ *)

    1 -> makeGridFactory[{1}],
    2 -> makeGridFactory[{2}],
    3 -> makeGridFactory[{3}],
    4 -> makeGridFactory[{4}],
    5 -> makeGridFactory[{5}],
    6 -> makeGridFactory[{10}],

    (* ------------------------------------------------------------------ *)
    (* Y 1-in 2-out, 4 vertices (cases 7, 8, 19)                          *)
    (* ------------------------------------------------------------------ *)

    7  -> makeAmFactory[{1,2,3,4}, $Y1In2OutAM],
    8  -> makeAmFactory[{1,2,3,4}, $Y1In2OutAM],
    19 -> makeAmFactory[{1,2,3,4}, $Y1In2OutAM],

    (* ------------------------------------------------------------------ *)
    (* Y 2-in 1-out, 4 vertices (cases 9, 10)                             *)
    (* ------------------------------------------------------------------ *)

    9  -> makeAmFactory[{1,2,3,4}, $Y2In1OutAM],
    10 -> makeAmFactory[{1,2,3,4}, $Y2In1OutAM],

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
                "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
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

    (* ------------------------------------------------------------------ *)
    (* Inconsistent switching (feature validation — infeasible by design)  *)
    (* ------------------------------------------------------------------ *)

    "Inconsistent Y shortcut"          -> makeAmFactory[{1,2,3,4}, $Y1In2OutAM],
    "Inconsistent attraction shortcut" -> makeAmFactory[{1,2,3,4}, $Attraction4AM],

    (* ------------------------------------------------------------------ *)
    (* Grid cases: directed GridGraph[{r,c}]                              *)
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

(* sc=Automatic resolves to the canonical SC via $CaseDefaultSC, or {} if undefined. *)
getExampleScenario[n_, entries_, exits_,
        sc_    : Automatic,
        alpha_ : $ExamplesDefaultAlpha,
        V_     : $ExamplesDefaultV,
        g_     : $ExamplesDefaultG] :=
    Module[{f = Lookup[$ExampleScenarios, n, $Failed]},
        If[f === $Failed, Return[$Failed]];
        f[entries, exits,
            If[sc === Automatic, Lookup[$CaseDefaultSC, n, {}], sc],
            alpha, V, g]
    ];

End[];

EndPackage[];
