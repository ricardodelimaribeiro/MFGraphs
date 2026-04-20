#!/usr/bin/env wolframscript

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];

failures = {};

assert[cond_, label_] := If[TrueQ[cond], Null, AppendTo[failures, label]];

(* README quick start smoke *)
data = MFGraphs`GetExampleData[12] /. {MFGraphs`I1 -> 100, MFGraphs`U1 -> 0};
d2e = MFGraphs`DataToEquations[data];
result = MFGraphs`CriticalCongestionSolver[d2e];

assert[AssociationQ[d2e], "DataToEquations returns an association"];
assert[AssociationQ[result], "CriticalCongestionSolver returns an association"];
assert[KeyExistsQ[result, "Feasibility"], "Result contains Feasibility key"];
assert[KeyExistsQ[result, "Solution"], "Result contains Solution key"];

(* README plotting snippet smoke *)
networkPlot = MFGraphs`NetworkGraphPlot[d2e];
flowPlot = MFGraphs`SolutionFlowPlot[d2e, result["Solution"]];

assert[Head[networkPlot] === Graph, "NetworkGraphPlot returns Graph"];
assert[Head[flowPlot] === Graph, "SolutionFlowPlot returns Graph"];

If[failures === {},
  Print["README smoke check: PASS"];
  Exit[0],
  Print["README smoke check: FAIL"];
  Scan[Print["  - ", #] &, failures];
  Exit[1]
];
