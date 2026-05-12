#!/usr/bin/env wolframscript
(* CompareSingleCaseBMR.wls — run only the booleanMinimizeReduceSystem stage
   on ONE case in a fresh kernel. Invoked by CompareSingleCase.wls so that
   BMR runs in its own process and cannot OOM-kill the kernel that holds the
   recursive / procedural / BooleanConvert / BooleanMinimize timings.

   Usage:
     wolframscript -file Scripts/CompareSingleCaseBMR.wls grid-4x5

   Output: prints exactly one machine-parseable line to stdout, prefixed
     BMR_RESULT|<verdict>|<seconds>
   where <verdict> is OK | TIMEOUT | FAIL and <seconds> is a Real. *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$BMRTimeout    = 600;    (* WL-level (often ignored by C kernel) *)
$WallclockKill = 800;    (* OS-level hard ceiling *)

$KnownCases = <|
    "grid-2x2" -> Function[gridScenario[{2, 2}, {{1, 100}}, {{4, 0}}]],
    "grid-3x2" -> Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]],
    "grid-2x3" -> Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]],
    "grid-2x4" -> Function[gridScenario[{2, 4}, {{1, 100}}, {{8, 0}}]],
    "grid-3x3" -> Function[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]],
    "grid-3x4" -> Function[gridScenario[{3, 4}, {{1, 100}}, {{12, 0}}]],
    "grid-4x4" -> Function[gridScenario[{4, 4}, {{1, 100}}, {{16, 0}}]],
    "grid-4x5" -> Function[gridScenario[{4, 5}, {{1, 100}}, {{20, 0}}]]
|>;

args = Rest[$ScriptCommandLine];
If[Length[args] < 1 || !KeyExistsQ[$KnownCases, First[args]],
    Print["BMR_RESULT|FAIL|0"];
    Exit[1]
];
caseName = First[args];

(* Watchdog: echo TIMEOUT just before SIGKILL — see CompareSingleCaseBC.wls. *)
Run["( sleep " <> ToString[$WallclockKill - 5] <> "; " <>
    "echo 'BMR_RESULT|TIMEOUT|" <> ToString[$WallclockKill] <> "'; " <>
    "sleep 5; " <>
    "kill -9 " <> ToString[$ProcessID] <>
    " ) 2>/dev/null &"];

scen = $KnownCases[caseName][];
sys  = makeSystem[scen];

{tBMR, resBMR} = AbsoluteTiming[
    TimeConstrained[
        Quiet[booleanMinimizeReduceSystem[sys,
              "ArmTimeout" -> 2, "DisjunctTimeout" -> 10],
              {booleanMinimizeReduceSystem::multisol}],
        $BMRTimeout, $Aborted
    ]
];
bmrVerdict = Which[
    resBMR === $Aborted,                   "TIMEOUT",
    ListQ[resBMR] || AssociationQ[resBMR], "OK",
    True,                                  "FAIL"
];

Print["BMR_RESULT|", bmrVerdict, "|", tBMR];
Exit[0];
