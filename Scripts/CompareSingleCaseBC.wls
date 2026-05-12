#!/usr/bin/env wolframscript
(* CompareSingleCaseBC.wls — run only the BooleanConvert[..., "DNF"] stage
   on ONE case in a fresh kernel. Invoked by CompareSingleCase.wls so that
   BooleanConvert runs in its own process and its known-bad TimeConstrained
   behavior (C-level Boolean kernel ignores aborts for many minutes past
   the requested limit) cannot block the main pipeline.

   A self-kill watchdog ensures hard wall-clock termination:
       Run["( sleep <wallclock>; kill -9 <pid> ) >/dev/null 2>&1 &"]
   so even if TimeConstrained never fires inside the kernel, the process
   is SIGKILLed by the OS at $WallclockKill seconds.

   Usage:
     wolframscript -file Scripts/CompareSingleCaseBC.wls grid-4x5

   Output: prints exactly one machine-parseable line to stdout:
     BC_RESULT|<verdict>|<seconds>|<disjuncts>
   where <verdict> is OK | TIMEOUT and <disjuncts> is an integer or "?". *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength   = 0;

$BoolConvertTimeout = 60;     (* WL-level (often ignored by C kernel) *)
$WallclockKill      = 120;    (* OS-level hard ceiling *)

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
    Print["BC_RESULT|FAIL|0|?"];
    Exit[1]
];
caseName = First[args];

(* Watchdog: echo a TIMEOUT result line 5s before SIGKILL so the parent
   captures a meaningful verdict even when the C-level Boolean kernel
   ignores the WL TimeConstrained abort. The sibling shell inherits
   stdout from this wolframscript invocation, so its echo reaches the
   parent's RunProcess capture. If BooleanConvert finishes normally
   first, this process exits cleanly and the watchdog's later echo
   writes to a closed pipe (silently lost). *)
Run["( sleep " <> ToString[$WallclockKill - 5] <> "; " <>
    "echo 'BC_RESULT|TIMEOUT|" <> ToString[$WallclockKill] <> "|?'; " <>
    "sleep 5; " <>
    "kill -9 " <> ToString[$ProcessID] <>
    " ) 2>/dev/null &"];

scen        = $KnownCases[caseName][];
sys         = makeSystem[scen];
inputs      = solversTools`Private`buildSolverInputs[sys];
constraints = inputs[[1]];

{tBC, resBC} = AbsoluteTiming[
    TimeConstrained[
        BooleanConvert[constraints, "DNF"],
        $BoolConvertTimeout, $Aborted
    ]
];

bcVerdict = If[resBC === $Aborted, "TIMEOUT", "OK"];
bcDj = Which[
    resBC === $Aborted,   "?",
    Head[resBC] === Or,   ToString[Length[resBC]],
    True,                 "1"
];

Print["BC_RESULT|", bcVerdict, "|", tBC, "|", bcDj];
Exit[0];
