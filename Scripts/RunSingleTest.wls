#!/usr/bin/env wolframscript
(* Run a single MUnit .mt file and print TestReport summary.
   Usage: wolframscript -file Scripts/RunSingleTest.wls MFGraphs/Tests/<name>.mt *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

args = Select[Rest[$ScriptCommandLine], !StringStartsQ[#, "-"] && !StringEndsQ[#, ".wls"] &];
If[args === {}, Print["Missing test file argument"]; Exit[1]];
testFile = First[args];
If[!FileExistsQ[testFile], Print["Test file not found: ", testFile]; Exit[1]];

Print["Running: ", testFile];
report = TestReport[testFile];

(* A failed TestReport or a file that ran zero tests must exit nonzero — a
   non-Boolean failure-count comparison would otherwise fall through to the
   unconditional Exit[0]. *)
If[!MatchQ[report, _TestReportObject],
    Print["TestReport failed to run (load error or bad file)."];
    Exit[1]
];
passed = report["TestsSucceededCount"];
failed = report["TestsFailedCount"];
If[!IntegerQ[passed] || !IntegerQ[failed] || passed + failed === 0,
    Print["No test results produced (file aborted or contains no tests)."];
    Exit[1]
];
Print["Passed: ", passed];
Print["Failed: ", failed];
If[failed > 0,
    Print["\nFailed tests:"];
    KeyValueMap[
        Print["  [", #2["Outcome"], "] ", #2["TestID"]] &,
        Select[report["TestResults"], !MatchQ[#["Outcome"], "Success"] &]
    ];
    Exit[1]
];
Exit[0];

(* :!CodeAnalysis::EndBlock:: *)
