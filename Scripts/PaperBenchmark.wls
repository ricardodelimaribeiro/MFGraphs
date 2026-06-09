#!/usr/bin/env wolframscript
(* Paper §5.4 benchmark (A.2) + C.9 smoke test.
   Usage: wolframscript -file Scripts/PaperBenchmark.wls [merge|fork|jamarat|all]
   Builds the exact §5 example networks, verifies topology (vertices/edges/q),
   smoke-tests the new edgeDensityProfile/edgeValueProfile helpers (C.9), and
   times the baseline (dnfReduceSystem) vs optimized (optimizedDNFReduceSystem)
   reducers directly on the built system (NOT solveScenario, which memoizes). *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

which = If[Length[$ScriptCommandLine] >= 2, ToLowerCase[$ScriptCommandLine[[2]]], "all"];

Print["Machine:     ", $Version];
Print["SystemID:    ", $SystemID];
Print["Date:        ", DateString[]];
Print["Selection:   ", which];
Print[StringRepeat["=", 60]];

Vfun = Function[x, (1/2) Sin[2 Pi (x + 1/4)]^2];
gId  = (# &);

(* --- the three §5 networks, built exactly as the EXE blocks do --- *)
buildMerge[] := amScenario[{1,2,3,4}, {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}},
    {{1,2},{3,3}}, {{4,5}}, {}, 1, Vfun, gId];
buildFork[]  := amScenario[{1,2,3,4}, {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}},
    {{1,2}}, {{3,2},{4,5}}, {}, 1, Vfun, gId];
(* Paper Jamarat: 7 interior vertices, entries v1,v2; exits v5,v6,v7.
   Interior edges e12,e24,e13,e34,e35,e46,e56,e57,e67 (Fig. jamarat). *)
buildJamarat[] := amScenario[Range[7],
    {{0,1,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,1,1,0,0},
     {0,0,0,0,0,1,0},{0,0,0,0,0,1,1},{0,0,0,0,0,0,1},{0,0,0,0,0,0,0}},
    {{1,100},{2,100}}, {{5,0},{6,0},{7,0}}, {}, 1, Vfun, gId];

countAlternatives[sys_] := Module[{cd},
    cd = systemData[sys, "ComplementarityData"];
    If[Head[cd] =!= mfgComplementarityData, Return[Missing["NoCompData"]]];
    Length[Lookup[First[cd], "AltFlows", {}]] +
    Length[Lookup[First[cd], "AltTransitionFlows", {}]] +
    Length[Lookup[First[cd], "AltOptCond", {}]]
];

structureReport[name_, s_, sys_] := Module[{aug},
    aug = augmentAuxiliaryGraph[sys];
    Print[name, ":"];
    Print["  augmented vertices: ", Length[aug["Vertices"]]];
    Print["  flow edges:         ", Length[aug["FlowEdges"]]];
    Print["  transition edges:   ", Length[aug["TransitionEdges"]]];
    Print["  alternatives (q):   ", countAlternatives[sys]];
];

(* --- C.9 smoke test: run merge1 end-to-end, confirm profiles are numeric --- *)
smokeC9[s_, sys_] := Module[{sol, rules, edges, dens, vals},
    sol = solveScenario[s];
    rules = If[AssociationQ[sol], sol["Rules"], sol];
    Print["  solveScenario valid: ", isValidSystemSolution[sys, sol]];
    edges = {{1,2},{3,2},{2,4}};
    dens = Table[edgeDensityProfile[s, sys, e, rules, xv],
        {e, edges}, {xv, {0, 0.25, 0.5, 0.75, 1.0}}];
    vals = Table[edgeValueProfile[s, sys, e, rules, xv],
        {e, edges}, {xv, {0, 0.25, 0.5, 0.75, 1.0}}];
    Print["  density samples all numeric: ", AllTrue[Flatten[dens], NumericQ]];
    Print["  value   samples all numeric: ", AllTrue[Flatten[vals], NumericQ]];
    Print["  density[e24] @ x=0.5: ", edgeDensityProfile[s, sys, {2,4}, rules, 0.5]];
    Print["  value  [e24] @ x=0.5: ", edgeValueProfile[s, sys, {2,4}, rules, 0.5]];
];

(* --- timing: reducer run directly on sys --- *)
timeFast[label_, sys_, fn_] := Module[{warm, rep, sol},
    {warm, sol} = AbsoluteTiming[fn[sys]];
    {rep, sol} = RepeatedTiming[fn[sys]];
    Print["  ", label, ": warmup ", Round[warm,0.001], " s | repeat ",
        Round[rep,0.001], " s | valid ", isValidSystemSolution[sys, sol]];
    <|"label" -> label, "warmup" -> warm, "repeat" -> rep, "sol" -> sol|>
];

results = <||>;

If[which === "merge" || which === "all",
    Module[{s, sys},
        s = buildMerge[]; sys = makeSystem[s];
        structureReport["ROAD MERGE", s, sys];
        Print["  -- C.9 smoke (merge1) --"];
        smokeC9[s, sys];
        Print["  -- timings --"];
        results["merge-baseline"]  = timeFast["baseline (dnf)", sys, dnfReduceSystem];
        results["merge-optimized"] = timeFast["optimized (optdnf)", sys, optimizedDNFReduceSystem];
        Print[];
    ];
];

If[which === "fork" || which === "all",
    Module[{s, sys},
        s = buildFork[]; sys = makeSystem[s];
        structureReport["ROAD FORK", s, sys];
        Print["  -- timings --"];
        results["fork-baseline"]  = timeFast["baseline (dnf)", sys, dnfReduceSystem];
        results["fork-optimized"] = timeFast["optimized (optdnf)", sys, optimizedDNFReduceSystem];
        Print[];
    ];
];

If[which === "jamarat" || which === "all",
    Module[{s, sys, t, sol},
        s = buildJamarat[]; sys = makeSystem[s];
        structureReport["JAMARAT", s, sys];
        Print["  -- timing (optimized only; baseline 2^95 is infeasible) --"];
        {t, sol} = AbsoluteTiming[optimizedDNFReduceSystem[sys]];
        Print["  optimized (optdnf): ", Round[t,0.01], " s | valid ",
            isValidSystemSolution[sys, sol]];
        results["jamarat-optimized"] = <|"warmup" -> t, "sol" -> sol|>;
        Print[];
    ];
];

(* persist solutions + timings *)
$out = FileNameJoin[{$RepoRoot, "Results", "paper_benchmark_" <> which <> ".wl"}];
Put[results, $out];
Print["Wrote: ", $out];
