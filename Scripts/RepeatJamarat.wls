#!/usr/bin/env wolframscript
(* Repeated timing of the Jamarat critical-congestion solve -- the "Optimized"
   Jamarat cell of the paper's tab:performance, which the paper's "under five
   minutes" claim rests on. A single run of that cell is close enough to the
   300 s threshold (~5% headroom) that one loaded run can cross it, so the cell
   is reported as a median over repeated runs rather than a single timing.

   Companion to PaperBenchmark.wls (which times the whole paper table once) and
   BaselinePruningBenchmark.wls (which contrasts the booleanReduceSystem
   baseline against the pruning solvers). This script does one thing: run the
   Jamarat cell N times and report the distribution.

   Scenario is built exactly as in BaselinePruningBenchmark.wls, so the numbers
   are comparable: V is numeric 0 because V does not enter the critical-case
   (alpha = 1) reduce -- only the post-hoc density recovery uses it.

   Solver is dnfReduceSystem, the default solveScenario reducer, timed directly
   on the built system. That bypasses the solveScenario memo cache; the solver's
   own $dnfSolveCache is Block-scoped and reset per call, so every run pays the
   full cost and the runs are independent.

   Usage:
     wolframscript -file Scripts/RepeatJamarat.wls [--runs N] [--timeout S] [--memcap GB]
   Defaults: --runs 5 --timeout 1800 --memcap 8
*)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

(* --- CLI --- *)
argValue[flag_String, default_] := Module[{pos},
    pos = FirstPosition[$ScriptCommandLine, flag];
    If[MissingQ[pos] || Length[$ScriptCommandLine] < First[pos] + 1,
        default,
        ToExpression[$ScriptCommandLine[[First[pos] + 1]]]]
];

nRuns   = argValue["--runs", 5];
tmo     = argValue["--timeout", 1800];
memCapGB = argValue["--memcap", 8];
memCap   = memCapGB * 10^9;

Print["Version:   ", $Version];
Print["SystemID:  ", $SystemID];
Print["Date:      ", DateString[]];
Print["Runs:      ", nRuns, "   (timeout ", tmo, " s, memory cap ",
      memCapGB, " GB)"];
Print[StringRepeat["=", 64]];

(* --- Scenario: identical to buildJamarat[] in BaselinePruningBenchmark.wls --- *)
gId = (# &);
buildJamarat[] := amScenario[Range[7],
    {{0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0},
     {0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0, 1},
     {0, 0, 0, 0, 0, 0, 0}},
    {{1, 100}, {2, 100}}, {{5, 0}, {6, 0}, {7, 0}}, {}, 1, 0, gId];

countAlternatives[sys_] := Module[{cd},
    cd = systemData[sys, "ComplementarityData"];
    If[Head[cd] =!= mfgComplementarityData, Return[Missing["NoCompData"]]];
    Length[Lookup[First[cd], "AltFlows", {}]] +
    Length[Lookup[First[cd], "AltTransitionFlows", {}]] +
    Length[Lookup[First[cd], "AltOptCond", {}]]
];

{tBuild, sys} = AbsoluteTiming[makeSystem[buildJamarat[]]];
Print["makeSystem: ", ToString[Round[tBuild, 0.001]], " s"];
Print["q (alternatives): ", countAlternatives[sys], "   (paper tab:performance: 36)"];
Print[StringRepeat["-", 64]];

(* --- Timed runs --- *)
times = {};
Do[
    Module[{t, sol, status, valid},
        {t, sol} = AbsoluteTiming[
            TimeConstrained[
                MemoryConstrained[dnfReduceSystem[sys], memCap, $MemOut],
                tmo, $TimedOut]];
        status = Which[sol === $TimedOut, "TIMEOUT", sol === $MemOut, "MEMOUT", True, "OK"];
        valid  = If[status === "OK", isValidSystemSolution[sys, sol], status];
        If[status === "OK", AppendTo[times, t]];
        Print["  run ", StringPadLeft[ToString[k], 2], ": ",
            StringPadLeft[If[status === "OK", ToString[Round[t, 0.01]], status], 10], " s",
            "  (", ToString[Round[t/60, 0.01]], " min)   valid=", valid];
    ],
    {k, nRuns}];

Print[StringRepeat["=", 64]];
If[times === {},
    Print["No run completed -- nothing to summarize."];
    ,
    Print["times (s):        ", Round[times, 0.01]];
    Print["median:           ", ToString[Round[Median[times], 0.01]], " s  = ",
          ToString[Round[Median[times]/60, 0.01]], " min"];
    Print["mean:             ", ToString[Round[Mean[times], 0.01]], " s"];
    Print["min / max:        ", ToString[Round[Min[times], 0.01]], " / ",
          ToString[Round[Max[times], 0.01]], " s"];
    If[Length[times] > 1,
        Print["median runs 2-n:  ", ToString[Round[Median[Rest[times]], 0.01]], " s",
              "   (kernel warm-up check: expect no material shift)"]];
    Print["under 300 s:      ", Count[times, t_ /; t < 300], " / ", Length[times], " runs"];
];
Print["\nDone at ", DateString[]];

(* :!CodeAnalysis::EndBlock:: *)
