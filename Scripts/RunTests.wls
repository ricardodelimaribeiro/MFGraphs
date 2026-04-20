#!/usr/bin/env wolframscript

Print["MFGraphs Test Runner"];
Print["==================="];

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

$TestsDir = FileNameJoin[{Directory[], "MFGraphs", "Tests"}];

$TestSuites = <|
    (* `fast` is the critical-only public gate.
       Keep archive files and unsupported legacy surfaces out of this list. *)
    "fast" -> {
        "package-loading.mt",
        "graphics-public-api.mt",
        "007-critical_congestion.mt",
        "007-critical_congestion_ok.mt",
        "012-critical_congestion.mt",
        "critical-numeric-backend.mt",
        "exit-pass-through.mt",
        "exit-pass-through-nonzero.mt",
        "exact-mode.mt",
        "infeasible-status.mt",
        "inconsistent-switching-critical.mt",
        "symbolic-underdetermined.mt",
        "scenario-kernel.mt",
        "solve-mfg.mt",
        "solver-contracts.mt"
    },
    "slow" -> {
        "J9F-critical_congestion.mt",
        "Jm9-critical_congestion.mt"
    }
|>;
$TestSuites["all"] = Join[$TestSuites["fast"], $TestSuites["slow"]];

args = Select[
    Rest[$ScriptCommandLine],
    !StringStartsQ[#, "-"] && !StringEndsQ[#, ".wls"] &
];
suiteName = ToLowerCase @ If[args === {}, "all", First[args]];
suiteFiles = Lookup[$TestSuites, suiteName, Missing["KeyAbsent", suiteName]];

If[MissingQ[suiteFiles],
    Print["Unknown suite: ", suiteName];
    Print["Available suites: ", StringRiffle[Keys[$TestSuites], ", "]];
    Exit[1];
];

Print["Suite: ", suiteName];
Print["Files: ", StringRiffle[suiteFiles, ", "]];
Print[];

report = TestReport[FileNameJoin[{$TestsDir, #}] & /@ suiteFiles];

Print["Tests run: ", report["TestsSucceededCount"] + report["TestsFailedCount"]];
Print["Passed: ", report["TestsSucceededCount"]];
Print["Failed: ", report["TestsFailedCount"]];

If[report["TestsFailedCount"] > 0,
    Print["\nFailed Tests:"];
    failed = Select[report["TestResults"], !MatchQ[#["Outcome"], "Success"] &];
    KeyValueMap[
        (Print["  [", #2["Outcome"], "] ", #2["TestID"]]) &,
        failed
    ];
    Exit[1]
,
    Exit[0]
];
