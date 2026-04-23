(* Wolfram Language package *)
(* Self-contained typed scenario factory library for all built-in MFGraphs examples.
   Each entry in $ExampleScenarios is a 6-arg Function:
     Function[{entries, exits, sc, alpha, V, gFunc}, makeScenario[...]]
   Topology is baked in at definition time; all boundary and Hamiltonian parameters
   are caller-supplied. Convenience accessor GetExampleScenario applies defaults.

   Usage:
     f = GetExampleScenario[7]
     f[{{1,80}}, {{3,0},{4,10}}, {}, 1, 0, Function[z,-1/z]]

     GetExampleScenario[7, {{1,80}}, {{3,0},{4,10}}]            (* uses defaults *)
     GetExampleScenario[8, {{1,80}}, {{3,0},{4,10}}]            (* uses canonical SC for case 8 *)
     GetExampleScenario[8, {{1,80}}, {{3,0},{4,10}}, {}]        (* override: no SC *)

   Numeric benchmark defaults are in Scripts/BenchmarkHelpers.wls ($DefaultParams/$CaseParams).
   NormalizeScenarioModel derives "Vertices List" + "Adjacency Matrix" from "Graph" key.
   Default Hamiltonian parameters are taken from $DefaultHamiltonian in Scenario.wl. *)

Begin["`Private`"];

(* --- Shared adjacency-matrix constants for custom topologies --- *)

$Y1In2OutAM    = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
$Y2In1OutAM    = {{0,1,0,0},{0,0,0,1},{0,1,0,0},{0,0,0,0}};
$Attraction4AM = {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}};

(* Canonical switching costs for cases that have a well-defined default SC.
   Looked up by GetExampleScenario when sc argument is Automatic. *)
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
    "Braess split"     -> {{1,2,4,1},{5,7,8,1}},
    "Braess congest"   -> {{1,2,4,1},{4,6,7,1}},
    "Big Braess split" -> {{1,3,6,1},{5,7,10,1}},
    "Big Braess congest" -> {{1,3,5,1},{5,7,9,1}},
    "Paper example"    -> {{1,2,3,2},{3,2,1,1},{2,3,4,3},{4,3,2,1}},
    (* SC that violates the triangle inequality — infeasible by design *)
    "Inconsistent Y shortcut" ->
        {{1,2,3,5},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}},
    "Inconsistent attraction shortcut" ->
        {{1,2,4,5},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
         {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
         {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}}
|>;

$MakeGraphFactory[graph_] :=
    With[{g = graph},
        Function[{entries, exits, sc, alpha, V, gFunc},
            makeScenario[<|
                "Model" -> <|
                    "Graph"                           -> g,
                    "Entrance Vertices and Flows"     -> entries,
                    "Exit Vertices and Terminal Costs" -> exits,
                    "Switching Costs"                 -> sc
                |>,
                "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> gFunc|>
            |>]]];

$MakeAMFactory[vl_, am_] :=
    With[{vertices = vl, matrix = am},
        Function[{entries, exits, sc, alpha, V, gFunc},
            makeScenario[<|
                "Model" -> <|
                    "Vertices List"                   -> vertices,
                    "Adjacency Matrix"                -> matrix,
                    "Entrance Vertices and Flows"     -> entries,
                    "Exit Vertices and Terminal Costs" -> exits,
                    "Switching Costs"                 -> sc
                |>,
                "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> gFunc|>
            |>]]];

(* --- Scenario registry --- *)

$ExampleScenarios = Association[

    (* ------------------------------------------------------------------ *)
    (* Linear chains (cases 1–6): directed GridGraph[{n}]                 *)
    (* ------------------------------------------------------------------ *)

    1 -> $MakeGraphFactory[GridGraph[{1},  DirectedEdges -> True]],
    2 -> $MakeGraphFactory[GridGraph[{2},  DirectedEdges -> True]],
    3 -> $MakeGraphFactory[GridGraph[{3},  DirectedEdges -> True]],
    4 -> $MakeGraphFactory[GridGraph[{4},  DirectedEdges -> True]],
    5 -> $MakeGraphFactory[GridGraph[{5},  DirectedEdges -> True]],
    6 -> $MakeGraphFactory[GridGraph[{10}, DirectedEdges -> True]],

    (* ------------------------------------------------------------------ *)
    (* Y 1-in 2-out, 4 vertices (cases 7, 8, 19)                          *)
    (* ------------------------------------------------------------------ *)

    7  -> $MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],
    8  -> $MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],
    19 -> $MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],

    (* ------------------------------------------------------------------ *)
    (* Y 2-in 1-out, 4 vertices (cases 9, 10)                             *)
    (* ------------------------------------------------------------------ *)

    9  -> $MakeAMFactory[{1,2,3,4}, $Y2In1OutAM],
    10 -> $MakeAMFactory[{1,2,3,4}, $Y2In1OutAM],

    (* ------------------------------------------------------------------ *)
    (* Attraction 4-vertex diamond (cases 11, 12, 13)                     *)
    (* ------------------------------------------------------------------ *)

    11 -> $MakeAMFactory[{1,2,3,4}, $Attraction4AM],
    12 -> $MakeAMFactory[{1,2,3,4}, $Attraction4AM],
    13 -> $MakeAMFactory[{1,2,3,4}, {{0,1,1,0},{0,0,0,1},{0,0,0,1},{0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Triangle 3-vertex directed cycle: 1->2->3->1                       *)
    (* CycleGraph[3, DirectedEdges->True]                                  *)
    (* ------------------------------------------------------------------ *)

    14                        -> $MakeGraphFactory[CycleGraph[3, DirectedEdges -> True]],
    104                       -> $MakeGraphFactory[CycleGraph[3, DirectedEdges -> True]],
    "triangle with two exits" -> $MakeGraphFactory[CycleGraph[3, DirectedEdges -> True]],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex directed chain with two exits (cases 105, alias)          *)
    (* ------------------------------------------------------------------ *)

    "chain with two exits" -> $MakeGraphFactory[GridGraph[{3}, DirectedEdges -> True]],
    105                    -> $MakeGraphFactory[GridGraph[{3}, DirectedEdges -> True]],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex misc (cases 15–18)                                        *)
    (* ------------------------------------------------------------------ *)

    15 -> $MakeAMFactory[{1,2,3}, {{0,0,0},{1,0,0},{1,0,0}}],
    16 -> $MakeAMFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],
    17 -> $MakeAMFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],
    18 -> $MakeAMFactory[{1,2,3}, {{0,1,0},{0,0,1},{0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Multi-entrance/exit (cases 20–23, "Jamaratv9")                     *)
    (* ------------------------------------------------------------------ *)

    20 -> $MakeAMFactory[Range[9], {
            {0,0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},
            {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,1,1,1},
            {0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0}}],

    21 -> $MakeAMFactory[Range[12], {
            {0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,0},
            {0,0,0,1,1,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,0},
            {0,0,0,0,0,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0},
            {0,0,0,0,0,0,0,1,1,1,0,0},{0,0,0,0,0,0,0,0,1,0,1,0},
            {0,0,0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}}],

    22 -> $MakeAMFactory[Range[7], {
            {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,1,0,0},
            {0,0,0,0,0,1,0},{0,0,0,0,0,1,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}}],

    "Jamaratv9" -> $MakeAMFactory[Range[9], {
            {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},
            {0,0,0,0,0,1,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
            {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}}],

    23 -> $MakeAMFactory[Range[6], {
            {0,1,1,0,0,0},{0,0,1,0,0,0},{0,0,0,1,1,0},
            {0,0,0,0,1,1},{0,0,0,0,0,1},{0,0,0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* 2-vertex undirected edge (case 27)                                  *)
    (* ------------------------------------------------------------------ *)

    27 -> $MakeAMFactory[{1,2}, {{0,1},{1,0}}],

    (* ------------------------------------------------------------------ *)
    (* Braess variants                                                     *)
    (* ------------------------------------------------------------------ *)

    "Braess split" -> $MakeAMFactory[Range[8], {
            {0,1,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},
            {0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1},
            {0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0}}],

    "Braess congest" -> $MakeAMFactory[Range[7], {
            {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,0,0,0},
            {0,0,0,0,1,1,0},{0,0,0,0,0,0,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}}],

    "New Braess" -> With[{newBraessAM = {
                {0,1,0,1,0,0},{0,0,1,0,0,0},{0,0,0,0,0,1},
                {0,0,0,0,1,0},{0,0,0,0,0,1},{0,0,0,0,0,0}}},
        Function[{entries, exits, sc, alpha, V, gFunc},
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
                "Hamiltonian" -> <|"Alpha" -> alpha, "V" -> V, "G" -> gFunc|>
            |>]]],

    "Big Braess split" -> $MakeAMFactory[Range[10], {
            {0,1,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},
            {0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},
            {0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},
            {0,0,0,0,0,0,0,0,0,0}}],

    "Big Braess congest" -> $MakeAMFactory[Range[9], {
            {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},
            {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
            {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}}],

    (* ------------------------------------------------------------------ *)
    (* Paper / benchmark cases                                             *)
    (* ------------------------------------------------------------------ *)

    "HRF Scenario 1" -> $MakeAMFactory[Range[10], {
            {0,1,0,0,0,0,0,0,0,0},{0,0,1,1,0,0,0,0,0,0},{0,0,0,1,1,0,1,0,0,0},
            {0,0,0,0,1,1,1,0,0,0},{0,0,0,0,0,1,1,0,0,0},{0,0,0,0,0,0,1,0,0,0},
            {0,0,0,0,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0}}],

    "Paper example" -> $MakeGraphFactory[GridGraph[{4}, DirectedEdges -> True]],

    (* ------------------------------------------------------------------ *)
    (* Inconsistent switching (feature validation — infeasible by design)  *)
    (* ------------------------------------------------------------------ *)

    "Inconsistent Y shortcut" ->
        $MakeAMFactory[{1,2,3,4}, $Y1In2OutAM],

    "Inconsistent attraction shortcut" ->
        $MakeAMFactory[{1,2,3,4}, $Attraction4AM],

    (* ------------------------------------------------------------------ *)
    (* Grid cases: directed GridGraph[{r,c}]                              *)
    (* ------------------------------------------------------------------ *)

    "Grid0303" -> $MakeGraphFactory[GridGraph[{3,3},   DirectedEdges -> True]],
    "Grid0404" -> $MakeGraphFactory[GridGraph[{4,4},   DirectedEdges -> True]],
    "Grid0505" -> $MakeGraphFactory[GridGraph[{5,5},   DirectedEdges -> True]],
    "Grid0707" -> $MakeGraphFactory[GridGraph[{7,7},   DirectedEdges -> True]],
    "Grid0710" -> $MakeGraphFactory[GridGraph[{7,10},  DirectedEdges -> True]],
    "Grid1010" -> $MakeGraphFactory[GridGraph[{10,10}, DirectedEdges -> True]],
    "Grid1020" -> $MakeGraphFactory[GridGraph[{10,20}, DirectedEdges -> True]]
];

(* --- Accessors --- *)

GetExampleScenario[n_] := Lookup[$ExampleScenarios, n, $Failed];

(* sc=Automatic resolves to the canonical SC for the case (from $CaseDefaultSC),
   or {} if the case has no canonical SC. Hamiltonian defaults from $DefaultHamiltonian. *)
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
