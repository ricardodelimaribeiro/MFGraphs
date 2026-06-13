#!/usr/bin/env wolframscript
(* C1: measure the TRUE §5.4 baseline -- booleanReduceSystem (BooleanConvert the
   full constraint system to DNF, materializing all 2^q complementarity branches,
   then Reduce per disjunct) -- against the DNF-FIRST pruning solvers
   dnfReduceSystem / optimizedDNFReduceSystem, which never materialize 2^q.

   This is the comparison tab:performance actually claims. The merge/fork/Jamarat
   networks are built exactly as the §5 EXE blocks; a grid ladder finds the
   largest case where the boolean baseline still finishes (the pruning-gap row).

   V is numeric 0 here: this branch lacks V-as-function support, and V does not
   affect the critical-case (alpha=1) reduce timing (only the post-hoc density
   recovery uses it). Usage:
     wolframscript -file Scripts/BaselinePruningBenchmark.wls [small|grids|jamarat|all]
*)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

which = If[Length[$ScriptCommandLine] >= 2, ToLowerCase[$ScriptCommandLine[[2]]], "all"];

Print["Machine:   ", $Version];
Print["SystemID:  ", $SystemID];
Print["Date:      ", DateString[]];
Print["Selection: ", which];
Print[StringRepeat["=", 64]];

gId = (# &);

(* §5 networks, built exactly as the EXE blocks (V numeric 0; see header note) *)
buildMerge[] := amScenario[{1, 2, 3, 4}, {{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
    {{1, 2}, {3, 3}}, {{4, 5}}, {}, 1, 0, gId];
buildFork[]  := amScenario[{1, 2, 3, 4}, {{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
    {{1, 2}}, {{3, 2}, {4, 5}}, {}, 1, 0, gId];
buildJamarat[] := amScenario[Range[7],
    {{0, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0},
     {0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0}},
    {{1, 100}, {2, 100}}, {{5, 0}, {6, 0}, {7, 0}}, {}, 1, 0, gId];

countAlternatives[sys_] := Module[{cd},
    cd = systemData[sys, "ComplementarityData"];
    If[Head[cd] =!= mfgComplementarityData, Return[Missing["NoCompData"]]];
    Length[Lookup[First[cd], "AltFlows", {}]] +
    Length[Lookup[First[cd], "AltTransitionFlows", {}]] +
    Length[Lookup[First[cd], "AltOptCond", {}]]
];

(* Time one solver on sys with a hard time+memory cap. Returns assoc with seconds,
   validity, and status (OK | TIMEOUT | MEMOUT | FAILED). A first guarded warmup
   establishes status; cells whose warmup finishes under $RepeatBelow are then
   re-timed with RepeatedTiming so the reported time is a stable steady-state mean
   (matching the paper table, whose caption promises "stable mean over repeated
   runs"). Slower cells keep the single warmup time so we never repeat a
   multi-minute solve, and TIMEOUT/MEMOUT cells are never repeated at all. *)
$MemCap = 8 * 10^9;  (* 8 GB guard so a BooleanConvert blow-up can't wedge the box *)
$RepeatBelow = 30;   (* seconds: warmup-time threshold below which we RepeatedTiming *)
guardedRun[fn_, sys_, tmo_] :=
    TimeConstrained[MemoryConstrained[fn[sys], $MemCap, $MemOut], tmo, $TimedOut];
timeSolver[label_, sys_, fn_, tmo_] := Module[{warm, sol, status, valid, t, mode},
    {warm, sol} = AbsoluteTiming[guardedRun[fn, sys, tmo]];
    status = Which[
        sol === $TimedOut, "TIMEOUT",
        sol === $MemOut,   "MEMOUT",
        AssociationQ[sol] || MatchQ[sol, {__Rule}] || sol === True, "OK",
        True, "FAILED"];
    valid = If[status === "OK", isValidSystemSolution[sys, sol], "-"];
    (* For the steady-state mean we time the BARE call: the warmup already proved
       this cell completes safely and quickly, so the time+memory guards (which add
       measurable per-call overhead at the millisecond scale) would only distort it. *)
    {t, mode} = Which[
        status =!= "OK",     {Missing[status], "-"},
        warm < $RepeatBelow, {First @ RepeatedTiming[fn[sys]], "mean"},
        True,                {warm, "single"}];
    Print["  ", StringPadRight[label, 22], ": ",
        If[status === "OK", ToString[Round[t, 0.0001]] <> " s (" <> mode <> ")", status],
        "  (", tmo, "s cap)  valid=", valid];
    <|"label" -> label, "seconds" -> t, "mode" -> mode,
      "status" -> status, "valid" -> valid|>
];

results = <||>;

runCase[name_, build_, tmoBool_, tmoDnf_] := Module[{s, sys, q},
    s = build[]; sys = makeSystem[s];
    q = countAlternatives[sys];
    Print[name, "  (q = ", q, ")"];
    results[name <> "-q"] = q;
    results[name <> "-boolean"]   = timeSolver["baseline (boolean)",  sys, booleanReduceSystem,        tmoBool];
    results[name <> "-dnf"]       = timeSolver["dnf",                  sys, dnfReduceSystem,            tmoDnf];
    results[name <> "-optimized"] = timeSolver["optimized (optdnf)",   sys, optimizedDNFReduceSystem,   tmoDnf];
    Print[];
];

If[which === "small" || which === "all",
    runCase["ROAD-MERGE", buildMerge, 60, 60];
    runCase["ROAD-FORK",  buildFork,  60, 60];
];

If[which === "grids" || which === "all",
    runCase["grid-2x3", Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]],     60, 60];
    runCase["grid-2x4", Function[gridScenario[{2, 4}, {{1, 100}}, {{8, 0}}]],     90, 60];
    runCase["grid-3x3", Function[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]],    120, 60];
];

If[which === "jamarat" || which === "all",
    (* boolean baseline gets a long cap to confirm 2^36 materialization is infeasible;
       optimized is already known ~281 s and is re-timed with a matching cap. *)
    runCase["JAMARAT", buildJamarat, 300, 360];
];

Print[StringRepeat["=", 64]];
Print["SUMMARY (seconds; T/O = timeout, M/O = memory cap):"];
Print[StringPadRight["Case", 14], StringPadRight["q", 6],
    StringPadRight["baseline(boolean)", 20], StringPadRight["dnf", 14], "optimized"];
fmt[r_] := If[MissingQ[r["seconds"]], r["status"], ToString[Round[r["seconds"], 0.0001]]];
Do[
    If[KeyExistsQ[results, c <> "-boolean"],
        Print[StringPadRight[c, 14], StringPadRight[ToString[results[c <> "-q"]], 6],
            StringPadRight[fmt[results[c <> "-boolean"]], 20],
            StringPadRight[fmt[results[c <> "-dnf"]], 14],
            fmt[results[c <> "-optimized"]]]],
    {c, {"ROAD-MERGE", "ROAD-FORK", "grid-2x3", "grid-2x4", "grid-3x3", "JAMARAT"}}
];
Print["\nDone at ", DateString[]];

(* :!CodeAnalysis::EndBlock:: *)
