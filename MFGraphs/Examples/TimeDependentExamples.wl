(* ::Package:: *)

(* Wolfram Language package *)
(* Time-dependent test cases for MFGraphs *)

GetTimeDependentExampleData::usage =
"GetTimeDependentExampleData[key] returns an Association with all standard network keys
plus time-dependent keys (\"Time Horizon\", \"Time Steps\", etc.).
Available keys: \"TD-1\" through \"TD-5\".";

GetTimeDependentExampleData::badkey =
"Unknown time-dependent example key: `1`.";

Begin["`Private`"];

(* ============================================================ *)
(* TD-1: Single edge, constant entrance, zero terminal cost      *)
(* Simplest case: should converge to the stationary solution     *)
(* at each time step since boundary conditions are constant.     *)
(* ============================================================ *)

tdTest["TD-1"] = <|
    "Vertices List" -> {1, 2},
    "Adjacency Matrix" -> {{0, 1}, {0, 0}},
    "Entrance Vertices and Flows" -> {{1, 100}},
    "Exit Vertices and Terminal Costs" -> {{2, 0}},
    "Switching Costs" -> {},
    "Time Horizon" -> 1.0,
    "Time Steps" -> 5,
    "Initial Mass Distribution" -> Function[{x, edge}, 0],
    "Terminal Cost Function" -> Function[{x, edge}, 0]
|>;

(* ============================================================ *)
(* TD-2: Linear 3-edge chain with time-varying entrance flow     *)
(* Entrance flow oscillates: I(t) = 100 (1 + 0.5 Sin[2 Pi t])  *)
(* Tests that different time steps produce different solutions.  *)
(* ============================================================ *)

tdTest["TD-2"] = <|
    "Vertices List" -> {1, 2, 3},
    "Adjacency Matrix" -> {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
    "Entrance Vertices and Flows" -> {{1, 100}},
    "Exit Vertices and Terminal Costs" -> {{3, 0}},
    "Switching Costs" -> {},
    "Time Horizon" -> 1.0,
    "Time Steps" -> 10,
    "Initial Mass Distribution" -> Function[{x, edge}, 0],
    "Terminal Cost Function" -> Function[{x, edge}, 0],
    "Time Dependent Entrance Flows" ->
        Function[t, {{1, 100 (1 + 0.5 Sin[2 Pi t])}}]
|>;

(* ============================================================ *)
(* TD-3: Y-network (case 7 topology) with pulsed entrance       *)
(* Entry at vertex 1, exits at vertices 3 and 4.                *)
(* Flow pulse: I(t) = 200 for t < 0.5, then 50.                *)
(* Tests flow splitting dynamics over time.                      *)
(* ============================================================ *)

tdTest["TD-3"] = <|
    "Vertices List" -> {1, 2, 3, 4},
    "Adjacency Matrix" -> {{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
    "Entrance Vertices and Flows" -> {{1, 100}},
    "Exit Vertices and Terminal Costs" -> {{3, 0}, {4, 0}},
    "Switching Costs" -> {},
    "Time Horizon" -> 1.0,
    "Time Steps" -> 10,
    "Initial Mass Distribution" -> Function[{x, edge}, 0],
    "Terminal Cost Function" -> Function[{x, edge}, 0],
    "Time Dependent Entrance Flows" ->
        Function[t, {{1, If[t < 0.5, 200, 50]}}]
|>;

(* ============================================================ *)
(* TD-4: Attraction network (case 12 topology) with ramp        *)
(* 1 entrance at vertex 1, 1 exit at vertex 4.                  *)
(* Two paths: 1->2->4 and 1->3->4 (plus middle edge 2->3).     *)
(* Entry flow ramps up linearly: I(t) = 50 + 150 t.            *)
(* Tests complementarity + time interaction.                     *)
(* ============================================================ *)

tdTest["TD-4"] = <|
    "Vertices List" -> {1, 2, 3, 4},
    "Adjacency Matrix" -> {{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
    "Entrance Vertices and Flows" -> {{1, 100}},
    "Exit Vertices and Terminal Costs" -> {{4, 0}},
    "Switching Costs" -> {},
    "Time Horizon" -> 1.0,
    "Time Steps" -> 10,
    "Initial Mass Distribution" -> Function[{x, edge}, 0],
    "Terminal Cost Function" -> Function[{x, edge}, 0],
    "Time Dependent Entrance Flows" ->
        Function[t, {{1, 50 + 150 t}}]
|>;

(* ============================================================ *)
(* TD-5: Single edge, analytical validation case                 *)
(* V=0, alpha=1, constant entrance I=100, terminal cost=0.      *)
(* With zero potential and constant data, the time-dependent     *)
(* solution should match the stationary solution at every step.  *)
(* This serves as a regression/validation test.                  *)
(* ============================================================ *)

tdTest["TD-5"] = <|
    "Vertices List" -> {1, 2},
    "Adjacency Matrix" -> {{0, 1}, {0, 0}},
    "Entrance Vertices and Flows" -> {{1, 100}},
    "Exit Vertices and Terminal Costs" -> {{2, 0}},
    "Switching Costs" -> {},
    "Time Horizon" -> 2.0,
    "Time Steps" -> 20,
    "Initial Mass Distribution" -> Function[{x, edge}, 0],
    "Terminal Cost Function" -> Function[{x, edge}, 0]
|>;

(* ============================================================ *)
(* Accessor function                                             *)
(* ============================================================ *)

$TimeDependentExampleKeys = {"TD-1", "TD-2", "TD-3", "TD-4", "TD-5"};

GetTimeDependentExampleData[key_String] :=
    If[MemberQ[$TimeDependentExampleKeys, key],
        tdTest[key],
        Message[GetTimeDependentExampleData::badkey, key]; $Failed
    ]

GetTimeDependentExampleData[key_] :=
    (Message[GetTimeDependentExampleData::badkey, key]; $Failed)

End[];
