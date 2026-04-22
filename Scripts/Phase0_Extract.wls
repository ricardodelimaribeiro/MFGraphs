#!/usr/bin/env wolframscript
(* LEGACY_GETEXAMPLEDATA: This script still depends on GetExampleData and is excluded from the active scenario-first workflow. *)
(* Phase 0: Extract Jamaratv9 basic data for verification *)

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Get[FileNameJoin[{$RepoRoot, "MFGraphs", "MFGraphs.wl"}]];
Get[FileNameJoin[{$RepoRoot, "Scripts", "BenchmarkHelpers.wls"}]];

(* Load Jamaratv9 data *)
data = GetExampleData["Jamaratv9"];
Print["=== JAMARATV9 BASIC DATA ==="];
Print[];

(* 0.1: Verify edge count *)
Print["0.1: EDGE COUNT VERIFICATION"];
Print["Vertices: ", data["Vertices List"]];
adj = data["Adjacency Matrix"];
edgePositions = Position[adj, 1];
Print["Directed edge count: ", Length[edgePositions]];
Print["Edges (as row,col):"];
Do[Print["  ", edgePositions[[i]][[1]], " -> ", edgePositions[[i]][[2]]], {i, 1, Length[edgePositions]}];
Print[];

(* 0.2: Extract data with default parameters *)
Print["0.2: PARAMETER SUBSTITUTION"];
params = GetParams["Jamaratv9"];
Print["Parameters to use: ", params];
dataSubstituted = data /. params;
Print["Entrances: ", dataSubstituted["Entrance Vertices and Flows"]];
Print["Exits: ", dataSubstituted["Exit Vertices and Terminal Costs"]];
Print["Switching Costs: ", dataSubstituted["Switching Costs"]];
Print[];

(* 0.3: Run preprocessing to find q *)
Print["0.3: PREPROCESSING - EXTRACT q VALUE"];
$MFGraphsVerbose = True;
d2e = DataToEquations[dataSubstituted];
initRules = d2e["InitRules"];
newSystem = d2e["NewSystem"];
Print["InitRules keys: ", Keys[initRules]];
Print["InitRules length: ", Length[initRules]];

(* Count OR branches in NewSystem *)
newSystemList = newSystem;
If[Head[newSystem] === List, newSystemList = newSystem, newSystemList = {newSystem}];
Print["NewSystem structure: ", Map[Head, newSystemList]];

(* Try to estimate q from NewSystem *)
If[Length[newSystemList] > 2,
  orBranches = newSystemList[[3]];
  If[Head[orBranches] === Or,
    Print["Or-branches count: ", Length[List @@ orBranches]],
    If[orBranches === True,
      Print["Or-branches: None (True)"],
      Print["Or-branches: Single constraint"]
    ]
  ]
];
Print[];

Print["=== END PHASE 0 ==="];
