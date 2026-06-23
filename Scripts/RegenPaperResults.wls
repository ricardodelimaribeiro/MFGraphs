#!/usr/bin/env wolframscript
(* RegenPaperResults.wls -- Regenerate the §5 result plots for the paper
   "An Exact Solver for Stationary Mean-Field Games on Networks" by replaying
   each EXE-block scenario through the current solver and exporting the
   density/value-profile grids to PDF.

   Companion to RegenPaperFigures.wls (which regenerates the network diagrams).
   Each plot reproduces its EXE block's exact data (entries/exits/switching
   costs/g/edges/PlotRange) via amScenario -> makeSystem -> solveScenario ->
   edgeDensityProfile/edgeValueProfile.

   Outputs (to the paper's executables/ directory):
     merge1result.pdf       -- road merge, g(m)=m
     merge3result.pdf       -- road merge, g(m)=log m
     fork1result.pdf        -- road fork, cheap exit at v3
     fork2result.pdf        -- road fork, cheap exit at v4
     fork3result.pdf        -- road fork, equal exit costs
     corridorresult.pdf     -- corridor, direct route vs bypass
     corridorpsiresult.pdf  -- corridor with a turn penalty psi=40 at v2

   Usage: wolframscript -file Scripts/RegenPaperResults.wls
*)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

exeDir = "/Users/ribeirrd/Dropbox/Apps/Overleaf/Stationary solutions of a congested MFG on Networks/executables";
If[!DirectoryQ[exeDir],
  Print["Executables directory not found: ", exeDir];
  Exit[1]
];

V = Function[x, (1/2) Sin[2 Pi (x + 1/4)]^2];
Id = (# &); Lg = (Log[#] &);
mfAdj  = {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};   (* merge & fork edge set *)
corAdj = {{0,1,0,1},{0,0,1,0},{0,0,0,1},{0,0,0,0}};   (* corridor: direct + bypass *)

res[label_, entries_, exits_, psi_, gfun_, adj_, edges_, densR_, valR_] :=
  Module[{s, sys, sol, rules, grid},
    s = amScenario[{1,2,3,4}, adj, entries, exits, psi, 1, V, gfun];
    sys = makeSystem[s, makeSymbolicUnknowns[s]];
    sol = solveScenario[s];
    rules = If[AssociationQ[sol], sol["Rules"], sol];
    grid = Grid[{
      Plot[edgeDensityProfile[s, sys, #, rules, x], {x, 0, 1},
        PlotLabel -> Row[{"density at ", #}], PlotRange -> densR] & /@ edges,
      Plot[edgeValueProfile[s, sys, #, rules, x], {x, 0, 1},
        PlotLabel -> Row[{"value function at ", #}], PlotRange -> valR] & /@ edges}];
    Export[FileNameJoin[{exeDir, label <> "result.pdf"}], grid];
    Print["Wrote ", label, "result.pdf  (solveScenario valid: ",
      isValidSystemSolution[sys, sol], ")"];
  ];

res["merge1", {{1,2},{3,3}}, {{4,5}}, {}, Id, mfAdj, {{1,2},{3,2},{2,4}}, {1,4},  {5,13}];
res["merge3", {{1,2},{3,3}}, {{4,5}}, {}, Lg, mfAdj, {{1,2},{3,2},{2,4}}, {0,25}, {5,13}];
res["fork1",  {{1,2}}, {{3,2},{4,5}}, {}, Id, mfAdj, {{1,2},{2,3},{2,4}}, {-0.5,2}, {2,7}];
res["fork2",  {{1,2}}, {{3,5},{4,2}}, {}, Id, mfAdj, {{1,2},{2,3},{2,4}}, {-0.5,2}, {2,7}];
res["fork3",  {{1,2}}, {{3,5},{4,5}}, {}, Id, mfAdj, {{1,2},{2,3},{2,4}}, {-0.5,2}, {2,9}];
res["corridor", {{1,100}}, {{4,0}}, {}, Id, corAdj, {{1,4},{1,2},{2,3},{3,4}}, {0,80}, {0,80}];
res["corridorpsi", {{1,100}}, {{4,0}}, {{1,2,3,40},{3,2,1,40}}, Id, corAdj,
    {{1,4},{1,2},{2,3},{3,4}}, {0,100}, {0,100}];

Print["Done."];
