(* Wolfram Language package *)
(* Self-contained typed scenario factory library for all built-in MFGraphs examples.
   Each entry is a Function[{entries, exits}, makeScenario[...]] that binds topology
   and switching costs but leaves boundary conditions to the caller.
   Numeric benchmark defaults are in Scripts/BenchmarkHelpers.wls ($DefaultParams/$CaseParams).
   Uses WL built-in graph constructors where possible; NormalizeScenarioModel
   derives "Vertices List" and "Adjacency Matrix" from the "Graph" key.

   Usage:
     GetExampleScenario[7]                          returns the factory Function
     GetExampleScenario[7][{{1,80}}, {{3,0},{4,10}}] returns a scenario[...] *)

Begin["`Private`"];

(* --- Shared adjacency-matrix constants for custom topologies --- *)

$Y1In2OutAM    = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
$Y2In1OutAM    = {{0,1,0,0},{0,0,0,1},{0,1,0,0},{0,0,0,0}};
$Attraction4AM = {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}};

(* --- Scenario factories --- *)
(* Each value is Function[{entries, exits}, makeScenario[...]] where
   entries = {{vertex, flow}, ...} and exits = {{vertex, cost}, ...}. *)

$ExampleScenarios = Association[

    (* ------------------------------------------------------------------ *)
    (* Linear chains (cases 1–6): GridGraph[{n}] directed chain            *)
    (* ------------------------------------------------------------------ *)

    1 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{1}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    2 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{2}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    3 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{3}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    4 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{4}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    5 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{5}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    6 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{10}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Y 1-in 2-out, 4 vertices (cases 7, 8, 19)                          *)
    (* ------------------------------------------------------------------ *)

    7 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Y1In2OutAM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    8 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Y1In2OutAM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,3,2},{1,2,4,3},{3,2,1,2},{3,2,4,1},{4,2,1,3},{4,2,3,1}}|>|>]],

    19 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Y1In2OutAM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,3,1},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Y 2-in 1-out, 4 vertices (cases 9, 10)                             *)
    (* ------------------------------------------------------------------ *)

    9 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Y2In1OutAM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    10 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Y2In1OutAM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,4,2},{1,2,3,3},{3,2,1,2},{3,2,4,1},{4,2,1,3},{4,2,3,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Attraction 4-vertex diamond (cases 11, 12, 13)                     *)
    (* ------------------------------------------------------------------ *)

    11 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Attraction4AM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,4,1},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
                    {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
                    {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}}|>|>]],

    12 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Attraction4AM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    13 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> {{0,1,1,0},{0,0,0,1},{0,0,0,1},{0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,4,1},{1,3,4,1},{4,2,1,1},{4,3,1,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Triangle 3-vertex directed cycle (cases 14, 104, "triangle with two exits") *)
    (* CycleGraph[3, DirectedEdges->True]: 1->2->3->1                     *)
    (* ------------------------------------------------------------------ *)

    14 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> CycleGraph[3, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,3,2},{3,2,1,1}}|>|>]],

    104 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> CycleGraph[3, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "triangle with two exits" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> CycleGraph[3, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex linear chain with two exits (cases 105, "chain with two exits") *)
    (* ------------------------------------------------------------------ *)

    "chain with two exits" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{3}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    105 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{3}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex misc (cases 15–18)                                        *)
    (* ------------------------------------------------------------------ *)

    15 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3},
                "Adjacency Matrix"               -> {{0,0,0},{1,0,0},{1,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{2,1,3,2},{3,1,2,1}}|>|>]],

    16 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3},
                "Adjacency Matrix"               -> {{0,1,0},{0,0,1},{0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,3,2,1}}|>|>]],

    17 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3},
                "Adjacency Matrix"               -> {{0,1,0},{0,0,1},{0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,3,2},{3,2,1,1}}|>|>]],

    18 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3},
                "Adjacency Matrix"               -> {{0,1,0},{0,0,1},{0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,3,2},{3,2,1,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Multi-entrance/exit (cases 20–23, "Jamaratv9")                     *)
    (* ------------------------------------------------------------------ *)

    20 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[9],
                "Adjacency Matrix"               -> {
                    {0,0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},
                    {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,1,1,1},
                    {0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    21 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[12],
                "Adjacency Matrix"               -> {
                    {0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,0},
                    {0,0,0,1,1,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,0},
                    {0,0,0,0,0,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0},
                    {0,0,0,0,0,0,0,1,1,1,0,0},{0,0,0,0,0,0,0,0,1,0,1,0},
                    {0,0,0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,0},
                    {0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    22 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[7],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,1,0,0},
                    {0,0,0,0,0,1,0},{0,0,0,0,0,1,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Jamaratv9" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[9],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},
                    {0,0,0,0,0,1,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
                    {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    23 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[6],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0},{0,0,1,0,0,0},{0,0,0,1,1,0},
                    {0,0,0,0,1,1},{0,0,0,0,0,1},{0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* 2-vertex undirected edge (case 27)                                  *)
    (* ------------------------------------------------------------------ *)

    27 -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2},
                "Adjacency Matrix"               -> {{0,1},{1,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Braess variants                                                     *)
    (* ------------------------------------------------------------------ *)

    "Braess split" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[8],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},
                    {0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1},
                    {0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,4,1},{5,7,8,1}}|>|>]],

    "Braess congest" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[7],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,0,0,0},
                    {0,0,0,0,1,1,0},{0,0,0,0,0,0,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,4,1},{4,6,7,1}}|>|>]],

    "New Braess" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[6],
                "Adjacency Matrix"               -> {
                    {0,1,0,1,0,0},{0,0,1,0,0,0},{0,0,0,0,0,1},
                    {0,0,0,0,1,0},{0,0,0,0,0,1},{0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,3,1},{3,2,1,1},{4,5,6,1},{6,5,4,1}},
                "a" -> Function[{j, edge},
                        Which[
                            edge === DirectedEdge[3,6] || edge === DirectedEdge[1,4], j/100,
                            True, 0
                        ]]|>|>]],

    "Big Braess split" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[10],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},
                    {0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},
                    {0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},
                    {0,0,0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,3,6,1},{5,7,10,1}}|>|>]],

    "Big Braess congest" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[9],
                "Adjacency Matrix"               -> {
                    {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},
                    {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
                    {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,3,5,1},{5,7,9,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Paper / benchmark cases                                             *)
    (* ------------------------------------------------------------------ *)

    "HRF Scenario 1" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> Range[10],
                "Adjacency Matrix"               -> {
                    {0,1,0,0,0,0,0,0,0,0},{0,0,1,1,0,0,0,0,0,0},{0,0,0,1,1,0,1,0,0,0},
                    {0,0,0,0,1,1,1,0,0,0},{0,0,0,0,0,1,1,0,0,0},{0,0,0,0,0,0,1,0,0,0},
                    {0,0,0,0,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},
                    {0,0,0,0,0,0,0,0,0,0}},
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Paper example" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{4}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {{1,2,3,2},{3,2,1,1},{2,3,4,3},{4,3,2,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Inconsistent switching (feature validation — infeasible by design)  *)
    (* ------------------------------------------------------------------ *)

    "Inconsistent Y shortcut" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Y1In2OutAM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,3,5},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}}|>|>]],

    "Inconsistent attraction shortcut" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Vertices List"                  -> {1, 2, 3, 4},
                "Adjacency Matrix"               -> $Attraction4AM,
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {
                    {1,2,4,5},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
                    {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
                    {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}}|>|>]],

    (* ------------------------------------------------------------------ *)
    (* Grid cases: GridGraph[{r,c}, DirectedEdges->True]                  *)
    (* ------------------------------------------------------------------ *)

    "Grid0303" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{3,3}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Grid0404" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{4,4}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Grid0505" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{5,5}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Grid0707" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{7,7}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Grid0710" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{7,10}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Grid1010" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{10,10}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]],

    "Grid1020" -> Function[{entries, exits},
            makeScenario[<|"Model" -> <|
                "Graph"                          -> GridGraph[{10,20}, DirectedEdges->True],
                "Entrance Vertices and Flows"    -> entries,
                "Exit Vertices and Terminal Costs" -> exits,
                "Switching Costs"                -> {}|>|>]]
];

GetExampleScenario[n_] := Lookup[$ExampleScenarios, n, $Failed];

End[];
