(* Wolfram Language package *)
(* Scenario constructors and named-example registry for MFGraphs.

   Direct constructors — topology clear from arguments:
     GridScenario[{n}, entries, exits]          chain, vertices 1..n
     GridScenario[{r,c}, entries, exits]        grid, vertices 1..r*c row-major
     CycleScenario[n, entries, exits]           directed n-cycle 1->2->...->n->1
     GraphScenario[graph, entries, exits]       any WL directed Graph object
     AMScenario[vl, am, entries, exits]         explicit vertices list + adjacency matrix

   Named examples (benchmark registry):
     GetExampleScenario[7, {{1,80}}, {{3,0},{4,10}}]      canonical SC for case 7
     GetExampleScenario[8, {{1,80}}, {{3,0},{4,10}}, {}]  override: no SC

   All constructors accept optional trailing args: sc, alpha, V, g.
   Defaults: sc={}, alpha=1, V=0, g=Function[z,-1/z] (from $DefaultHamiltonian).
   Numeric benchmark defaults are in Scripts/BenchmarkHelpers.wls. *)

Begin["`Private`"];

(* --- Shared topology constants --- *)

$Y1In2OutAM    = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
$Y2In1OutAM    = {{0,1,0,0},{0,0,0,1},{0,1,0,0},{0,0,0,0}};
$Attraction4AM = {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}};

(* Canonical switching costs looked up by GetExampleScenario when sc=Automatic. *)
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

GridScenario[dims_List, entries_, exits_,
        sc_    : {},
        alpha_ : $DefaultHamiltonian["Alpha"],
        V_     : $DefaultHamiltonian["V"],
        g_     : $DefaultHamiltonian["G"]] :=
    makeScenario[<|
        "Model" -> <|
            "Graph"                           -> GridGraph[dims, DirectedEdges -> True],
            "Entrance Vertices and Flows"     -> entries,
            "Exit Vertices and Terminal Costs" -> exits,
            "Switching Costs"                 -> sc
        |>,
        "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
    |>];

CycleScenario[n_Integer, entries_, exits_,
        sc_    : {},
        alpha_ : $DefaultHamiltonian["Alpha"],
        V_     : $DefaultHamiltonian["V"],
        g_     : $DefaultHamiltonian["G"]] :=
    makeScenario[<|
        "Model" -> <|
            "Graph"                           -> CycleGraph[n, DirectedEdges -> True],
            "Entrance Vertices and Flows"     -> entries,
            "Exit Vertices and Terminal Costs" -> exits,
            "Switching Costs"                 -> sc
        |>,
        "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
    |>];

GraphScenario[graph_, entries_, exits_,
        sc_    : {},
        alpha_ : $DefaultHamiltonian["Alpha"],
        V_     : $DefaultHamiltonian["V"],
        g_     : $DefaultHamiltonian["G"]] :=
    makeScenario[<|
        "Model" -> <|
            "Graph"                           -> graph,
            "Entrance Vertices and Flows"     -> entries,
            "Exit Vertices and Terminal Costs" -> exits,
            "Switching Costs"                 -> sc
        |>,
        "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
    |>];

AMScenario[vl_, am_, entries_, exits_,
        sc_    : {},
        alpha_ : $DefaultHamiltonian["Alpha"],
        V_     : $DefaultHamiltonian["V"],
        g_     : $DefaultHamiltonian["G"]] :=
    makeScenario[<|
        "Model" -> <|
            "Vertices List"                   -> vl,
            "Adjacency Matrix"                -> am,
            "Entrance Vertices and Flows"     -> entries,
            "Exit Vertices and Terminal Costs" -> exits,
            "Switching Costs"                 -> sc
        |>,
        "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
    |>];

(* --- Private factory helpers — thin wrappers used by $ExampleScenarios --- *)

MakeGridFactory[dims_List]  := Function[{entries, exits, sc, alpha, V, g}, GridScenario[dims,  entries, exits, sc, alpha, V, g]];
MakeCycleFactory[n_Integer] := Function[{entries, exits, sc, alpha, V, g}, CycleScenario[n,    entries, exits, sc, alpha, V, g]];
MakeAMFactory[vl_, am_]     := Function[{entries, exits, sc, alpha, V, g}, AMScenario[vl, am, entries, exits, sc, alpha, V, g]];

(* --- Scenario registry --- *)

$ExampleScenarios = Association[

    (* ------------------------------------------------------------------ *)
    (* Linear chains (cases 1–6): directed GridGraph[{n}]                 *)
    (* ------------------------------------------------------------------ *)

    1 -> MakeGridFactory[{1}],
    2 -> MakeGridFactory[{2}],
    3 -> MakeGridFactory[{3}],
    4 -> MakeGridFactory[{4}],
    5 -> MakeGridFactory[{5}],
    6 -> MakeGridFactory[{10}],

    (* ------------------------------------------------------------------ *)
    (* Y 1-in 2-out, 4 vertices (cases 7, 8, 19)                          *)
    (* ------------------------------------------------------------------ *)

    7  -> MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],
    8  -> MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],
    19 -> MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],

    (* ------------------------------------------------------------------ *)
    (* Y 2-in 1-out, 4 vertices (cases 9, 10)                             *)
    (* ------------------------------------------------------------------ *)

    9  -> MakeAMFactory[{1,2,3,4}, $Y2In1OutAM],
    10 -> MakeAMFactory[{1,2,3,4}, $Y2In1OutAM],

    (* ------------------------------------------------------------------ *)
    (* Attraction 4-vertex diamond (cases 11, 12, 13)                     *)
    (* ------------------------------------------------------------------ *)

    11 -> MakeAMFactory[{1,2,3,4}, $Attraction4AM],
    12 -> MakeAMFactory[{1,2,3,4}, $Attraction4AM],
    13 -> MakeAMFactory[{1,2,3,4}, {{0,1,1,0},{0,0,0,1},{0,0,0,1},{0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Triangle 3-vertex directed cycle: 1->2->3->1                       *)
    (* ------------------------------------------------------------------ *)

    14                        -> MakeCycleFactory[3],
    104                       -> MakeCycleFactory[3],
    "triangle with two exits" -> MakeCycleFactory[3],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex directed chain with two exits (cases 105, alias)          *)
    (* ------------------------------------------------------------------ *)

    "chain with two exits" -> MakeGridFactory[{3}],
    105                    -> MakeGridFactory[{3}],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex misc (cases 15–18)                                        *)
    (* ------------------------------------------------------------------ *)

    15 -> MakeAMFactory[{1,2,3}, {{0,0,0},{1,0,0},{1,0,0}}],
    16 -> MakeAMFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],
    17 -> MakeAMFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],
    18 -> MakeAMFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Multi-entrance/exit (cases 20–23, "Jamaratv9")                     *)
    (* ------------------------------------------------------------------ *)

    20 -> MakeAMFactory[Range[9], {
            {0,0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},
            {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,1,1,1},
            {0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0}}],

    21 -> MakeAMFactory[Range[12], {
            {0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,0},
            {0,0,0,1,1,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,0},
            {0,0,0,0,0,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0},
            {0,0,0,0,0,0,0,1,1,1,0,0},{0,0,0,0,0,0,0,0,1,0,1,0},
            {0,0,0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}}],

    22 -> MakeAMFactory[Range[7], {
            {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,1,0,0},
            {0,0,0,0,0,1,0},{0,0,0,0,0,1,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}}],

    "Jamaratv9" -> MakeAMFactory[Range[9], {
            {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},
            {0,0,0,0,0,1,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
            {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}}],

    23 -> MakeAMFactory[Range[6], {
            {0,1,1,0,0,0},{0,0,1,0,0,0},{0,0,0,1,1,0},
            {0,0,0,0,1,1},{0,0,0,0,0,1},{0,0,0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* 2-vertex undirected edge (case 27)                                  *)
    (* ------------------------------------------------------------------ *)

    27 -> MakeAMFactory[{1,2}, {{0,1},{1,0}}],

    (* ------------------------------------------------------------------ *)
    (* Braess variants                                                     *)
    (* ------------------------------------------------------------------ *)

    "Braess split" -> MakeAMFactory[Range[8], {
            {0,1,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},
            {0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1},
            {0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0}}],

    "Braess congest" -> MakeAMFactory[Range[7], {
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
                    "Vertices List"                   -> Range[6],
                    "Adjacency Matrix"                -> newBraessAM,
                    "Entrance Vertices and Flows"     -> entries,
                    "Exit Vertices and Terminal Costs" -> exits,
                    "Switching Costs"                 -> sc,
                    "a" -> Function[{j, edge},
                            Which[
                                edge === DirectedEdge[3,6] || edge === DirectedEdge[1,4], j/100,
                                True, 0]]
                |>,
                "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> g|>
            |>]]],

    "Big Braess split" -> MakeAMFactory[Range[10], {
            {0,1,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},
            {0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},
            {0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},
            {0,0,0,0,0,0,0,0,0,0}}],

    "Big Braess congest" -> MakeAMFactory[Range[9], {
            {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},
            {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
            {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Paper / benchmark cases                                             *)
    (* ------------------------------------------------------------------ *)

    "HRF Scenario 1" -> MakeAMFactory[Range[10], {
            {0,1,0,0,0,0,0,0,0,0},{0,0,1,1,0,0,0,0,0,0},{0,0,0,1,1,0,1,0,0,0},
            {0,0,0,0,1,1,1,0,0,0},{0,0,0,0,0,1,1,0,0,0},{0,0,0,0,0,0,1,0,0,0},
            {0,0,0,0,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0}}],

    "Paper example" -> MakeGridFactory[{4}],

    (* ------------------------------------------------------------------ *)
    (* Inconsistent switching (feature validation — infeasible by design)  *)
    (* ------------------------------------------------------------------ *)

    "Inconsistent Y shortcut"          -> MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],
    "Inconsistent attraction shortcut" -> MakeAMFactory[{1,2,3,4}, $Attraction4AM],

    (* ------------------------------------------------------------------ *)
    (* Grid cases: directed GridGraph[{r,c}]                              *)
    (* ------------------------------------------------------------------ *)

    "Grid0303" -> MakeGridFactory[{3,3}],
    "Grid0404" -> MakeGridFactory[{4,4}],
    "Grid0505" -> MakeGridFactory[{5,5}],
    "Grid0707" -> MakeGridFactory[{7,7}],
    "Grid0710" -> MakeGridFactory[{7,10}],
    "Grid1010" -> MakeGridFactory[{10,10}],
    "Grid1020" -> MakeGridFactory[{10,20}]
];

(* --- Accessors --- *)

GetExampleScenario[n_] := Lookup[$ExampleScenarios, n, $Failed];

(* sc=Automatic resolves to the canonical SC via $CaseDefaultSC, or {} if undefined. *)
GetExampleScenario[n_, entries_, exits_,
        sc_    : Automatic,
        alpha_ : $DefaultHamiltonian["Alpha"],
        V_     : $DefaultHamiltonian["V"],
        g_     : $DefaultHamiltonian["G"]] :=
    Module[{f = Lookup[$ExampleScenarios, n, $Failed]},
        If[f === $Failed, Return[$Failed]];
        f[entries, exits,
            If[sc === Automatic, Lookup[$CaseDefaultSC, n, {}], sc],
            alpha, V, g]
    ];

End[];
