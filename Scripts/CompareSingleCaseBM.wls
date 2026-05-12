#!/usr/bin/env wolframscript
(* CompareSingleCaseBM.wls — run only the booleanMinimizeSystem stage on
   ONE case in a fresh kernel. Invoked by CompareSingleCase.wls so that
   booleanMinimizeSystem's known-bad TimeConstrained behavior (the inner
   per-disjunct loop, plus the C-level Boolean kernel's abort-deafness,
   can push wall-clock to many times the requested 60s limit) cannot
   block the main pipeline.

   A self-kill watchdog ensures hard wall-clock termination at
   $WallclockKill seconds even if TimeConstrained never fires inside
   the kernel.

   Usage:
     wolframscript -file Scripts/CompareSingleCaseBM.wls grid-4x5

   Output: prints exactly one machine-parseable line to stdout:
     BM_RESULT|<verdict>|<seconds>
   where <verdict> is OK | TIMEOUT | FAIL. *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength   = 0;

$BMTimeout     = 60;     (* WL-level (often ignored by C kernel) *)
$WallclockKill = 120;    (* OS-level hard ceiling *)

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
    Print["BM_RESULT|FAIL|0"];
    Exit[1]
];
caseName = First[args];

(* Spawn the OS-level watchdog. *)
Run["( sleep " <> ToString[$WallclockKill] <> "; kill -9 " <>
    ToString[$ProcessID] <> " ) >/dev/null 2>&1 &"];

scen = $KnownCases[caseName][];
sys  = makeSystem[scen];

{tBM, resBM} = AbsoluteTiming[
    TimeConstrained[
        Quiet[booleanMinimizeSystem[sys, "DisjunctTimeout" -> 10],
              {booleanMinimizeSystem::multisol}],
        $BMTimeout, $Aborted
    ]
];

bmVerdict = Which[
    resBM === $Aborted,                  "TIMEOUT",
    ListQ[resBM] || AssociationQ[resBM], "OK",
    True,                                "FAIL"
];

Print["BM_RESULT|", bmVerdict, "|", tBM];
Exit[0];
