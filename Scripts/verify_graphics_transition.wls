#!/usr/bin/env wolframscript
(* Usage: wolframscript -file Scripts/verify_graphics_transition.wls *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

fail[label_String, value_] := (
    Print["FAILURE: ", label, " returned ", Head[value]];
    Exit[1]
);

Print["--- Manual Verification: graphics transition plots ---"];
s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
sys = makeSystem[s];

Print["Solving known-fast scenario..."];
sol = TimeConstrained[solveScenario[s], 20, $TimedOut];
If[sol === $TimedOut, (
    Print["FAILURE: solveScenario timed out."];
    Exit[1]
)];

Print["Calling mfgTransitionPlot..."];
plot = mfgTransitionPlot[s, sys, sol];
Print["Head of plot: ", Head[plot]];
If[!MatchQ[plot, _Graph], fail["mfgTransitionPlot", plot]];
Print["SUCCESS: mfgTransitionPlot returned a Graph."];

Print["\n--- Manual Verification: mfgSolutionPlot styling ---"];
plotSol = mfgSolutionPlot[s, sys, sol];
Print["Head of solution plot: ", Head[plotSol]];
If[!MatchQ[plotSol, _Graph], fail["mfgSolutionPlot", plotSol]];
Print["SUCCESS: mfgSolutionPlot returned a Graph."];
