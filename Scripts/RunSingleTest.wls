#!/usr/bin/env wolframscript
(* Run a single MUnit .mt file and print TestReport summary.
   Usage: wolframscript -file Scripts/RunSingleTest.wls MFGraphs/Tests/<name>.mt *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

args = Select[Rest[$ScriptCommandLine], !StringEndsQ[#, ".wls"] &];
If[args === {}, Print["Missing test file argument"]; Exit[1]];
testFile = First[args];

Print["Running: ", testFile];
report = TestReport[testFile];

Print["Passed: ", report["TestsSucceededCount"]];
Print["Failed: ", report["TestsFailedCount"]];
If[report["TestsFailedCount"] > 0,
    Print["\nFailed tests:"];
    KeyValueMap[
        Print["  [", #2["Outcome"], "] ", #2["TestID"]] &,
        Select[report["TestResults"], !MatchQ[#["Outcome"], "Success"] &]
    ];
    Exit[1]
];
Exit[0];

(* :!CodeAnalysis::EndBlock:: *)
