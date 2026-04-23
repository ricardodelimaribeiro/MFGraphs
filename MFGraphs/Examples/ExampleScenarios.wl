(* Wolfram Language package *)
(* Self-contained typed scenario library for all built-in MFGraphs examples.
   All boundary values are concrete numbers — no symbolic parameters.
   Numeric values match the benchmark defaults in Scripts/BenchmarkHelpers.wls.
   Uses WL built-in graph constructors where possible; NormalizeScenarioModel
   derives "Vertices List" and "Adjacency Matrix" from the "Graph" key. *)

Begin["`Private`"];

(* --- Shared adjacency-matrix constants for custom topologies --- *)

$Y1In2OutAM    = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
$Y2In1OutAM    = {{0,1,0,0},{0,0,0,1},{0,1,0,0},{0,0,0,0}};
$Attraction4AM = {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}};

(* --- Scenario definitions --- *)

$ExampleScenarios = Association[

    (* ------------------------------------------------------------------ *)
    (* Linear chains (cases 1–6): GridGraph[{n}] directed chain            *)
    (* ------------------------------------------------------------------ *)

    1 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{1}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{1, 0}},
            "Switching Costs"                -> {}|>|>],

    2 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{2}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}},
            "Switching Costs"                -> {}|>|>],

    3 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{3}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs"                -> {}|>|>],

    4 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{4}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {}|>|>],

    5 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{5}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{5, 0}},
            "Switching Costs"                -> {}|>|>],

    6 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{10}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{10, 0}},
            "Switching Costs"                -> {}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Y 1-in 2-out, 4 vertices (cases 7, 8, 19)                          *)
    (* ------------------------------------------------------------------ *)

    7 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Y1In2OutAM,
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}, {4, 10}},
            "Switching Costs"                -> {}|>|>],

    8 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Y1In2OutAM,
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}, {4, 10}},
            "Switching Costs"                -> {
                {1,2,3,2},{1,2,4,3},{3,2,1,2},{3,2,4,1},{4,2,1,3},{4,2,3,1}}|>|>],

    19 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Y1In2OutAM,
            "Entrance Vertices and Flows"    -> {{1, 100}, {4, 50}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs"                -> {
                {1,2,3,1},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Y 2-in 1-out, 4 vertices (cases 9, 10)                             *)
    (* ------------------------------------------------------------------ *)

    9 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Y2In1OutAM,
            "Entrance Vertices and Flows"    -> {{1, 100}, {3, 40}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {}|>|>],

    10 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Y2In1OutAM,
            "Entrance Vertices and Flows"    -> {{1, 100}, {3, 40}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {
                {1,2,4,2},{1,2,3,3},{3,2,1,2},{3,2,4,1},{4,2,1,3},{4,2,3,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Attraction 4-vertex diamond (cases 11, 12, 13)                     *)
    (* ------------------------------------------------------------------ *)

    11 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Attraction4AM,
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {
                {1,2,4,1},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
                {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
                {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}}|>|>],

    12 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Attraction4AM,
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {}|>|>],

    13 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> {{0,1,1,0},{0,0,0,1},{0,0,0,1},{0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {
                {1,2,4,1},{1,3,4,1},{4,2,1,1},{4,3,1,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Triangle 3-vertex directed cycle (cases 14, 104, "triangle with two exits") *)
    (* CycleGraph[3, DirectedEdges->True]: 1->2->3->1                     *)
    (* ------------------------------------------------------------------ *)

    14 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> CycleGraph[3, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs"                -> {{1,2,3,2},{3,2,1,1}}|>|>],

    104 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> CycleGraph[3, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 80}, {2, 40}},
            "Exit Vertices and Terminal Costs" -> {{1, 0}, {2, 0}, {3, 0}},
            "Switching Costs"                -> {}|>|>],

    "triangle with two exits" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> CycleGraph[3, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}, {3, 10}},
            "Switching Costs"                -> {}|>|>],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex linear chain with two exits (cases 105, "chain with two exits") *)
    (* ------------------------------------------------------------------ *)

    "chain with two exits" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{3}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}, {3, 10}},
            "Switching Costs"                -> {}|>|>],

    105 -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{3}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}, {3, 10}},
            "Switching Costs"                -> {}|>|>],

    (* ------------------------------------------------------------------ *)
    (* 3-vertex misc (cases 15–18)                                        *)
    (* ------------------------------------------------------------------ *)

    15 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3},
            "Adjacency Matrix"               -> {{0,0,0},{1,0,0},{1,0,0}},
            "Entrance Vertices and Flows"    -> {{2, 80}, {3, 40}},
            "Exit Vertices and Terminal Costs" -> {{1, 0}},
            "Switching Costs"                -> {{2,1,3,2},{3,1,2,1}}|>|>],

    16 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3},
            "Adjacency Matrix"               -> {{0,1,0},{0,0,1},{0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs"                -> {{1,3,2,1}}|>|>],

    17 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3},
            "Adjacency Matrix"               -> {{0,1,0},{0,0,1},{0,0,0}},
            "Entrance Vertices and Flows"    -> {{2, 80}},
            "Exit Vertices and Terminal Costs" -> {{1, 0}, {3, 10}},
            "Switching Costs"                -> {{1,2,3,2},{3,2,1,1}}|>|>],

    18 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3},
            "Adjacency Matrix"               -> {{0,1,0},{0,0,1},{0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 80}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs"                -> {{1,2,3,2},{3,2,1,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Multi-entrance/exit (cases 20–23, "Jamaratv9")                     *)
    (* ------------------------------------------------------------------ *)

    20 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[9],
            "Adjacency Matrix"               -> {
                {0,0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},
                {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,1,1,1},
                {0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}, {2, 50}},
            "Exit Vertices and Terminal Costs" -> {{7, 0}, {8, 0}, {9, 0}},
            "Switching Costs"                -> {}|>|>],

    21 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[12],
            "Adjacency Matrix"               -> {
                {0,0,1,0,0,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0,0,0},
                {0,0,0,1,1,0,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0,0,0},
                {0,0,0,0,0,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,0},
                {0,0,0,0,0,0,0,1,1,1,0,0},{0,0,0,0,0,0,0,0,1,0,1,0},
                {0,0,0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,0},
                {0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}, {2, 50}},
            "Exit Vertices and Terminal Costs" -> {{10, 0}, {11, 0}, {12, 0}},
            "Switching Costs"                -> {}|>|>],

    22 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[7],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,1,0,0},
                {0,0,0,0,0,1,0},{0,0,0,0,0,1,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}, {2, 50}},
            "Exit Vertices and Terminal Costs" -> {{5, 0}, {6, 0}, {7, 0}},
            "Switching Costs"                -> {}|>|>],

    "Jamaratv9" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[9],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},
                {0,0,0,0,0,1,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
                {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}, {2, 50}},
            "Exit Vertices and Terminal Costs" -> {{7, 0}, {8, 0}, {9, 0}},
            "Switching Costs"                -> {}|>|>],

    23 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[6],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0},{0,0,1,0,0,0},{0,0,0,1,1,0},
                {0,0,0,0,1,1},{0,0,0,0,0,1},{0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}, {2, 50}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}, {5, 0}, {6, 0}},
            "Switching Costs"                -> {}|>|>],

    (* ------------------------------------------------------------------ *)
    (* 2-vertex undirected edge (case 27)                                  *)
    (* ------------------------------------------------------------------ *)

    27 -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2},
            "Adjacency Matrix"               -> {{0,1},{1,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}},
            "Switching Costs"                -> {}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Braess variants                                                     *)
    (* ------------------------------------------------------------------ *)

    "Braess split" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[8],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},
                {0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1},
                {0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{8, 0}},
            "Switching Costs"                -> {{1,2,4,1},{5,7,8,1}}|>|>],

    "Braess congest" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[7],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,0,0,0},
                {0,0,0,0,1,1,0},{0,0,0,0,0,0,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{7, 0}},
            "Switching Costs"                -> {{1,2,4,1},{4,6,7,1}}|>|>],

    "New Braess" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[6],
            "Adjacency Matrix"               -> {
                {0,1,0,1,0,0},{0,0,1,0,0,0},{0,0,0,0,0,1},
                {0,0,0,0,1,0},{0,0,0,0,0,1},{0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{6, 0}},
            "Switching Costs"                -> {{1,2,3,1},{3,2,1,1},{4,5,6,1},{6,5,4,1}},
            "a" -> Function[{j, edge},
                    Which[
                        edge === DirectedEdge[3,6] || edge === DirectedEdge[1,4], j/100,
                        True, 0
                    ]]|>|>],

    "Big Braess split" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[10],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},
                {0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},
                {0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},
                {0,0,0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{10, 0}},
            "Switching Costs"                -> {{1,3,6,1},{5,7,10,1}}|>|>],

    "Big Braess congest" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[9],
            "Adjacency Matrix"               -> {
                {0,1,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},
                {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,1,0},
                {0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{9, 0}},
            "Switching Costs"                -> {{1,3,5,1},{5,7,9,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Paper / benchmark cases                                             *)
    (* ------------------------------------------------------------------ *)

    "HRF Scenario 1" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> Range[10],
            "Adjacency Matrix"               -> {
                {0,1,0,0,0,0,0,0,0,0},{0,0,1,1,0,0,0,0,0,0},{0,0,0,1,1,0,1,0,0,0},
                {0,0,0,0,1,1,1,0,0,0},{0,0,0,0,0,1,1,0,0,0},{0,0,0,0,0,0,1,0,0,0},
                {0,0,0,0,0,0,0,1,0,1},{0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},
                {0,0,0,0,0,0,0,0,0,0}},
            "Entrance Vertices and Flows"    -> {{1, 100}, {9, 100}},
            "Exit Vertices and Terminal Costs" -> {{8, 0}, {10, 0}},
            "Switching Costs"                -> {}|>|>],

    "Paper example" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{4}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 80}, {3, 40}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}, {4, 0}},
            "Switching Costs"                -> {{1,2,3,2},{3,2,1,1},{2,3,4,3},{4,3,2,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Inconsistent switching (feature validation — not feasible)         *)
    (* ------------------------------------------------------------------ *)

    "Inconsistent Y shortcut" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Y1In2OutAM,
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}, {4, 0}},
            "Switching Costs"                -> {
                {1,2,3,5},{1,2,4,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1}}|>|>],

    "Inconsistent attraction shortcut" -> makeScenario[<|"Model" -> <|
            "Vertices List"                  -> {1, 2, 3, 4},
            "Adjacency Matrix"               -> $Attraction4AM,
            "Entrance Vertices and Flows"    -> {{1, 20}},
            "Exit Vertices and Terminal Costs" -> {{4, 0}},
            "Switching Costs"                -> {
                {1,2,4,5},{1,2,3,1},{3,2,1,1},{3,2,4,1},{4,2,1,1},{4,2,3,1},
                {1,3,4,1},{4,3,1,1},{1,3,2,1},{2,3,1,1},{3,4,2,1},{2,4,3,1},
                {2,3,4,1},{4,3,2,1},{3,1,2,1},{2,1,3,1}}|>|>],

    (* ------------------------------------------------------------------ *)
    (* Grid cases: GridGraph[{r,c}, DirectedEdges->True]                  *)
    (* ------------------------------------------------------------------ *)

    "Grid0303" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{3,3}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{9, 0}},
            "Switching Costs"                -> {}|>|>],

    "Grid0404" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{4,4}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{16, 0}},
            "Switching Costs"                -> {}|>|>],

    "Grid0505" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{5,5}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{25, 0}},
            "Switching Costs"                -> {}|>|>],

    "Grid0707" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{7,7}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{49, 0}},
            "Switching Costs"                -> {}|>|>],

    "Grid0710" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{7,10}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{70, 0}},
            "Switching Costs"                -> {}|>|>],

    "Grid1010" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{10,10}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{100, 0}},
            "Switching Costs"                -> {}|>|>],

    "Grid1020" -> makeScenario[<|"Model" -> <|
            "Graph"                          -> GridGraph[{10,20}, DirectedEdges->True],
            "Entrance Vertices and Flows"    -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{200, 0}},
            "Switching Costs"                -> {}|>|>]
];

GetExampleScenario[n_] := Lookup[$ExampleScenarios, n, $Failed];

End[];
